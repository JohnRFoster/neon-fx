## main transferability workflow script
## 1 - setup
## 3 - state data intake
## 4 - weather data intake
## 5 - create informative priors
## 6 - get initial conditions
## 7 - forecast step
## 7a - mice
## 7b - ticks
## 8 - analysis step
## 9 - save

# =========================================== #
#       1 - setup
# =========================================== #

library(tidyverse)
library(lubridate)
library(nimble)
library(parallel)
# nimbleOptions('MCMCjointlySamplePredictiveBranches' = FALSE)

update <- TRUE

dir.top <- "/projectnb/dietzelab/fosterj"
dir.analysis <- file.path(
	dir.top,
	"FinalOut/A_Correct/Analysis/dormantStages/nimbleModels"
)
dir.out <- "out"
# dir.out <- file.path(dir.top, "FinalOut/Chapter3/outConstraintForest")
if (update) {
	dir.out <- paste0(dir.out, "Update")
}

models <- c("Static", "WithMNAMice", "Weather", "WithWeatherAndMiceGlobal")
species <- c("Ixodes scapularis", "Amblyomma americanum")
neon.sites <- c(
	"BLAN",
	"HARV",
	"KONZ",
	"LENO",
	"OSBS",
	"SCBI",
	"SERC",
	"TALL",
	"TREE",
	"UKFS"
)

jobs <- expand_grid(
	model = models,
	species = species,
	site = neon.sites
)

jobs <- jobs %>%
	filter(
		!(site == "HARV" & species == "Amblyomma americanum"),
		!(site == "TREE" & species == "Amblyomma americanum"),
		!(site == "KONZ" & species == "Ixodes scapularis"),
		!(site == "OSBS" & species == "Ixodes scapularis"),
		!(site == "TALL" & species == "Ixodes scapularis"),
		!(site == "UKFS" & species == "Ixodes scapularis")
	)

# jobs <- jobs %>% filter(grepl("Weather", model))

job.num <- as.numeric(Sys.getenv("SGE_TASK_ID"))
if (is.na(job.num)) {
	job.num <- 47
}

site.job <- jobs$site[job.num]
species.job <- jobs$species[job.num]
model.job <- jobs$model[job.num]

message("====== Running simulation at ", site.job, " ======")
message(species.job)
message(model.job)


ua.cal <-
	if_else(
		model.job == "WithWeatherAndMiceGlobal",
		"mice_ic_parameter_process",
		"ic_parameter_process"
	)

n.slots <- Sys.getenv("NSLOTS") %>% as.numeric()
production <- TRUE
n.iter <- 50000
Nmc <- 2000
horizon <- 365

# =========================================== #
#       tick data intake
# =========================================== #
source("Functions/neon_tick_data.R")
neon.data <- neon_tick_data(species.job) %>% suppressMessages()

neon.job <- neon.data %>%
	filter(siteID == site.job, grepl("Forest", nlcd), time >= "2018-01-01") %>%
	arrange(time)

drag.dates <- neon.job$time %>% unique()
start.date <- first(drag.dates)
n.drags <- length(drag.dates)

# =========================================== #
#       get initial conditions
# =========================================== #

df.latent <- read_csv(file.path(dir.analysis, "dormantNymphTimeSeries.csv"))
month.get <- if_else(month(start.date) < 5, 4, month(start.date))
data.latent <- df.latent %>%
	mutate(model = gsub("DormantNymph", "", model)) %>%
	filter(
		model == model.job,
		type == "latent",
		statistic == "conf_50",
		ua == ua.cal,
		month(DATE) == month.get
	) %>%
	group_by(lifeStage) %>%
	summarise(mu = mean(value), prec = 1 / var(value)) %>%
	pivot_wider(names_from = lifeStage, values_from = c(mu, prec))

IC <- tibble(
	mu = c(
		pull(data.latent, mu_larvae),
		pull(data.latent, mu_dormant),
		pull(data.latent, mu_nymphs),
		pull(data.latent, mu_adults)
	),
	prec = c(
		pull(data.latent, prec_larvae),
		pull(data.latent, prec_dormant),
		pull(data.latent, prec_nymphs),
		pull(data.latent, prec_adults)
	)
) %>%
	as.matrix()

# =========================================== #
#       mouse data intake
# =========================================== #

source("../neon-tick-smam/Functions/capture_matrix.R")
neon.smam <- read_csv("../neon-tick-smam/Data/allSmallMammals.csv")
ch.ls <- capture_matrix(site.job, neon.smam)
ch <- ch.ls$ch
alive <- ch %in% 1:3
ch[alive] <- 1
ch[!alive] <- 0
ncaps <- rowSums(ch)
ch <- ch[ncaps > 0, ]
source("Functions/known_states.R")
ks <- known_states(ch)
mna <- colSums(ks)
mice.obs <- ymd(colnames(ch)) # unique sampling days: mice

# every day in mouse sequence
mice.seq <- seq.Date(mice.obs[1], mice.obs[length(mice.obs)], by = 1)

mna.all.days <- rep(NA, length(mice.seq))
mna.count <- 1
for (i in seq_along(mice.seq)) {
	if (mice.seq[i] %in% mice.obs) {
		mna.all.days[i] <- mna[mna.count]
		mna.count <- mna.count + 1
	} else {
		mna.all.days[i] <- mna[mna.count]
	}
}

# historical mna
source("Functions/mna_jags.R")
mna.hist <- mna_jags("Green Control", return.mean = TRUE)

# center and scale
mna.scaled <- tibble(
	mna.scaled = (mna.all.days - mna.hist$mean) / mna.hist$sd,
	Date = mice.seq
)

# =========================================== #
#       daymet intake and correction
# =========================================== #
source("Functions/daymet_downscale.R")
cgdd <- daymet_cumGDD(site.job) %>% suppressMessages()
maxTemp <- daymet_temp(site.job, minimum = FALSE) %>%
	select(Date, maxTempCorrect) %>%
	suppressMessages()
rh <- daymet_rh(site.job) %>%
	select(Date, maxRHCorrect, minRHCorrect) %>%
	suppressMessages()
precip <- daymet_precip(site.job) %>%
	select(Date, precipitation) %>%
	suppressMessages()

source("Functions/scale_met_forecast.R")
hist.means <- scale_met_forecast()

join1 <- left_join(maxTemp, rh, by = "Date")
join2 <- left_join(join1, precip, by = "Date")

df.daymet <- join2 %>%
	mutate(
		maxTempScale = (maxTempCorrect - hist.means$means["MAX_TEMP"]) /
			hist.means$sds["MAX_TEMP"],
		maxRHScale = (maxRHCorrect - hist.means$means["MAX_RH"]) /
			hist.means$sds["MAX_RH"],
		minRHScale = (minRHCorrect - hist.means$means["MIN_RH"]) /
			hist.means$sds["MIN_RH"],
		precipScale = (precipitation - hist.means$means["TOT_PREC"]) /
			hist.means$sds["TOT_PREC"]
	) %>%
	ungroup() %>%
	select(Date, contains("Scale"))

# =========================================== #
#       get informative priors
# =========================================== #
df.params <- read_csv(file.path(dir.analysis, "dormantNymphParams.csv"))
params.stats <- df.params %>%
	filter(model == model.job) %>%
	select(parameter, value) %>%
	group_by(parameter) %>%
	summarise(mu = mean(value), tau = 1 / var(value))

get_prior <- function(name) {
	pr <- numeric(2)
	xx <- params.stats %>%
		filter(parameter == name)
	pr[1] <- xx %>% pull(mu)
	pr[2] <- xx %>% pull(tau)
	pr
}

phi.l <- get_prior("phi.l.mu")
phi.n <- get_prior("phi.n.mu")
phi.a <- get_prior("phi.a.mu")
theta.l2n <- get_prior("theta.ln")
theta.n2a <- get_prior("theta.na")
repro <- get_prior("repro.mu")
repro.mu <- repro[1]

n.beta <- params.stats %>%
	filter(grepl("beta", parameter)) %>%
	nrow()

get_beta <- function(model.job) {
	if (model.job == "Static") {
		return(NA)
	} else {
		pr.beta <- matrix(NA, n.beta, 2)
		n.beta <- params.stats %>%
			filter(grepl("beta", parameter)) %>%
			nrow()

		pr.beta <- matrix(NA, n.beta, 2)
		for (i in seq_len(n.beta)) {
			if (model.job == "WithMNAMice") {
				node <- paste0("beta.m[", i, "]")
			} else {
				node <- paste0("beta[", i, "]")
			}
			pr.beta[i, ] <- get_prior(node)
		}
		return(pr.beta)
	}
}

pr.beta <- get_beta(model.job)

# function to approximate moment the inverse gamma
inv_gamma_mm <- function(x) {
	mu <- mean(x)
	v <- var(x)
	alpha <- (mu^2 / v) + 2
	beta <- mu * ((mu^2 / v) + 1)
	return(c("alpha" = alpha, "beta" = beta))
}

# get invgamma parameters
pr.sig <- df.params %>%
	filter(model == model.job, grepl("sig", parameter)) %>%
	select(parameter, value) %>%
	group_by(parameter) %>%
	summarise(alpha = inv_gamma_mm(value)[1], beta = inv_gamma_mm(value)[2])


# iterate ======================================================================================

ua.type <- c(
	"ic",
	"ic_parameter",
	"ic_parameter_driver",
	"ic_parameter_driver_process"
)

t = 1
k = 2

for (t in seq_len(n.drags)) {
	fx.start.date <- drag.dates[t]
	message("---------------------------------------------------")
	mm <- paste(fx.start.date, " (", round(t / n.drags * 100, 2), "%)")
	message(mm)

	for (k in seq_along(ua.type)) {
		ua.job <- ua.type[k]
		message("Simulating ", ua.job)

		# flags
		miceAndWeather <- model.job == "WithWeatherAndMiceGlobal"
		miceMNA <- model.job == "WithMNAMice"
		notStatic <- model.job != "Static"
		use.daymet <- grepl("Weather", model.job)
		ic <- grepl("ic", ua.job)
		parameter <- grepl("parameter", ua.job)
		driver <- grepl("driver", ua.job)
		process <- grepl("process", ua.job)

		dir.base <- file.path(
			dir.out,
			site.job,
			model.job,
			gsub(" ", "", species.job)
		)
		dir.save <- file.path(dir.base, ua.job)

		# initialize nimble lists
		constants <- data <- list()

		if (t == 1) {
			fx.sequence <- seq.Date(fx.start.date, by = 1, length.out = horizon)
			# n.days <- horizon
			y <- matrix(NA, 4, horizon)
		} else {
			# read last forecast parameters and state
			readDest <- file.path(
				dir.base,
				"ic_parameter_driver_process",
				drag.dates[t - 1]
			)

			last.params <- read_csv(file.path(readDest, "parameterSamples.csv")) %>%
				suppressMessages()

			# get parameter posterior summary
			params.stats <- last.params %>%
				rename(parameter = node) %>%
				select(parameter, value) %>%
				group_by(parameter) %>%
				summarise(mu = mean(value), tau = 1 / var(value))

			if (update) {
				phi.l <- get_prior("phi.l.mu")
				phi.n <- get_prior("phi.n.mu")
				phi.a <- get_prior("phi.a.mu")
				theta.l2n <- get_prior("theta.ln")
				theta.n2a <- get_prior("theta.na")
				# repro <- get_prior("repro.mu")

				if (model.job != "Static") {
					pr.beta <- matrix(NA, n.beta, 2)
					for (i in seq_len(n.beta)) {
						pr.beta[i, ] <- get_prior(paste0("beta[", i, "]"))
					}
				}

				# get invgamma parameters
				pr.sig <- last.params %>%
					rename(parameter = node) %>%
					filter(grepl("sig", parameter)) %>%
					select(parameter, value) %>%
					group_by(parameter) %>%
					summarise(
						alpha = inv_gamma_mm(value)[1],
						beta = inv_gamma_mm(value)[2]
					)
			}

			last.fx <- read_csv(file.path(readDest, "stateSamples.csv")) %>%
				suppressMessages()

			tick.stats <- last.fx %>%
				filter(time == fx.start.date) %>%
				group_by(lifeStage, time) %>%
				summarise(mu = mean(value), tau = 1 / var(value))

			IC <- matrix(NA, 4, 2)
			IC[1, 1] <- tick.stats %>% filter(lifeStage == "Larva") %>% pull(mu)
			IC[1, 2] <- tick.stats %>% filter(lifeStage == "Larva") %>% pull(tau)
			IC[2, 1] <- tick.stats %>% filter(lifeStage == "Dormant") %>% pull(mu)
			IC[2, 2] <- tick.stats %>% filter(lifeStage == "Dormant") %>% pull(tau)
			IC[3, 1] <- tick.stats %>% filter(lifeStage == "Nymph") %>% pull(mu)
			IC[3, 2] <- tick.stats %>% filter(lifeStage == "Nymph") %>% pull(tau)
			IC[4, 1] <- tick.stats %>% filter(lifeStage == "Adult") %>% pull(mu)
			IC[4, 2] <- tick.stats %>% filter(lifeStage == "Adult") %>% pull(tau)

			fx.sequence <- seq.Date(fx.start.date, by = 1, length.out = horizon)
			n.days <- length(fx.sequence)
		}

		if (use.daymet) {
			daymet.sub <- df.daymet %>%
				filter(Date %in% fx.sequence)
			data$maxtemp <- daymet.sub %>% pull(maxTempScale) %>% as.vector()
			data$maxrh <- daymet.sub %>% pull(maxRHScale) %>% as.vector()
			data$minrh <- daymet.sub %>% pull(minRHScale) %>% as.vector()
			data$precip <- daymet.sub %>% pull(precipScale) %>% as.vector()
		}

		obs <- neon.job %>%
			filter(time == fx.start.date)

		plots <- unique(obs$plotID)
		n.plots <- length(plots)

		y <- array(NA, dim = c(4, horizon, n.plots))
		area <- matrix(NA, horizon, n.plots)
		for (p in 1:n.plots) {
			obs.plot <- obs %>% filter(plotID == plots[p])
			y[1, 1, p] <- obs.plot %>% pull(Larva)
			y[3, 1, p] <- obs.plot %>% pull(Nymph)
			y[4, 1, p] <- obs.plot %>% pull(Adult)
			area[1, p] <- obs.plot %>% pull(totalSampledArea)
		}

		if (model.job == "Static") {
			n.beta <- 2
			pr.beta <- matrix(1, n.beta, 2)
		}

		# finalize data
		data$y <- y
		data$area <- area
		data$IC <- IC
		data$pr.phi.l <- phi.l
		data$pr.phi.n <- phi.n
		data$pr.phi.a <- phi.a
		data$pr.theta.l2n <- theta.l2n
		data$pr.theta.n2a <- theta.n2a
		# data$pr.repro <- repro
		data$repro.mu <- repro.mu
		data$pr.beta <- pr.beta
		data$pr.sig <- pr.sig %>% select(-parameter) %>% as.matrix()
		data$cgdd <- cgdd %>%
			filter(Date %in% ymd(fx.sequence)) %>%
			pull(cumGDD)
		data$max.cgdd <- max(data$cgdd) * 1.2
		data$xind <- matrix(1, 4, horizon)

		if (miceAndWeather || miceMNA) {
			data$mice <- mna.scaled %>%
				filter(Date %in% fx.sequence) %>%
				pull(mna.scaled)

			if (length(data$mice) < length(fx.sequence)) {
				horizon <- min(length(data$cgdd), length(data$mice))
				data$y <- y[, 1:horizon, ]
			}
		}

		if (year(fx.start.date) == 2021) {
			if (model.job == "Static" || model.job == "Weather") {
				horizon <- length(data$cgdd)
				data$y <- y[, 1:horizon, ]
			} else {
				horizon <- min(length(data$cgdd), length(data$mice))
				data$y <- y[, 1:horizon, ]
			}
		}

		# finalize constants
		constants$n.beta <- n.beta
		constants$n.plots <- n.plots
		constants$horizon <- horizon
		constants$ns <- 4

		area.init <- area
		nai <- which(is.na(area))
		area.init[nai] <- 160
		area.init[-nai] <- NA

		# build inits
		if (ua.job == "ic_parameter_driver_process") {
			inits <- function() {
				list(
					area = area.init,
					phi.l.mu = rnorm(1, phi.l[1], 1 / sqrt(phi.l[2])),
					phi.n.mu = rnorm(1, phi.n[1], 1 / sqrt(phi.n[2])),
					phi.a.mu = rnorm(1, phi.a[1], 1 / sqrt(phi.a[2])),
					theta.ln = rnorm(1, theta.l2n[1], 1 / sqrt(theta.l2n[2])),
					theta.na = rnorm(1, theta.n2a[1], 1 / sqrt(theta.n2a[2])),
					beta = rnorm(n.beta, pr.beta[, 1], 1 / sqrt(pr.beta[, 2])),
					sig = rinvgamma(4, pr.sig$alpha, pr.sig$beta),
					x = matrix(rpois(4 * horizon, 2) / 160 * 450, 4, horizon),
					Ex = matrix(rpois(4 * horizon, 2) / 160 * 450, 4, horizon),
					y = array(rpois(4 * horizon * n.plots, 5), dim = dim(data$y)),
					tau.temp = rexp(1),
					tau.maxrh = rexp(1),
					tau.minrh = rexp(1),
					tau.precip = rexp(1),
					tau.cgdd = rexp(1),
					x1 = jitter(data$maxtemp),
					x2 = jitter(data$maxrh),
					x3 = jitter(data$minrh),
					x4 = jitter(data$precip),
					gdd = jitter(data$cgdd)
				)
			}
		} else if (ua.job == "ic_parameter_driver") {
			inits <- function() {
				x.init <- matrix(rpois(4 * horizon, 2) / 160 * 450, 4, horizon)
				list(
					area = area.init,
					phi.l.mu = rnorm(1, phi.l[1], 1 / sqrt(phi.l[2])),
					phi.n.mu = rnorm(1, phi.n[1], 1 / sqrt(phi.n[2])),
					phi.a.mu = rnorm(1, phi.a[1], 1 / sqrt(phi.a[2])),
					theta.ln = rnorm(1, theta.l2n[1], 1 / sqrt(theta.l2n[2])),
					theta.na = rnorm(1, theta.n2a[1], 1 / sqrt(theta.n2a[2])),
					beta = rnorm(n.beta, pr.beta[, 1], 1 / sqrt(pr.beta[, 2])),
					x = x.init,
					Ex = x.init,
					y = array(rpois(4 * horizon * n.plots, 5), dim = dim(data$y)),
					tau.temp = rexp(1),
					tau.maxrh = rexp(1),
					tau.minrh = rexp(1),
					tau.precip = rexp(1),
					tau.cgdd = rexp(1),
					x1 = jitter(data$maxtemp),
					x2 = jitter(data$maxrh),
					x3 = jitter(data$minrh),
					x4 = jitter(data$precip),
					gdd = jitter(data$cgdd)
				)
			}
		} else if (ua.job == "ic_parameter") {
			inits <- function() {
				x.init <- matrix(rpois(4 * horizon, 2) / 160 * 450, 4, horizon)
				list(
					area = area.init,
					phi.l.mu = rnorm(1, phi.l[1], 1 / sqrt(phi.l[2])),
					phi.n.mu = rnorm(1, phi.n[1], 1 / sqrt(phi.n[2])),
					phi.a.mu = rnorm(1, phi.a[1], 1 / sqrt(phi.a[2])),
					theta.ln = rnorm(1, theta.l2n[1], 1 / sqrt(theta.l2n[2])),
					theta.na = rnorm(1, theta.n2a[1], 1 / sqrt(theta.n2a[2])),
					beta = rnorm(n.beta, pr.beta[, 1], 1 / sqrt(pr.beta[, 2])),
					x = x.init,
					Ex = x.init,
					y = array(rpois(4 * horizon * n.plots, 5), dim = dim(data$y))
				)
			}
		} else if (ua.job == "ic") {
			inits <- function() {
				x.init <- matrix(rpois(4 * horizon, 2) / 160 * 450, 4, horizon)
				list(
					area = area.init,
					x = x.init,
					Ex = x.init,
					y = array(rpois(4 * horizon * n.plots, 5), dim = dim(data$y))
				)
			}
		}

		source("Scripts/DA_neon.R")
		source("Functions/run_transfer_nimble.R")
		cl <- makeCluster(n.slots)
		out.nchains <- run_transfer_nimble(
			cl = cl,
			model = model.code,
			data = data,
			constants = constants,
			inits = inits,
			n.iter = n.iter,
			notStatic = notStatic,
			ic = ic,
			parameter = parameter,
			driver = driver,
			process = process,
			miceMNA = miceMNA,
			miceAndWeather = miceAndWeather,
			use.daymet = use.daymet
		)
		stopCluster(cl)

		# nimbleOptions('MCMCjointlySamplePredictiveBranches' = FALSE)
		# model <- nimbleModel(model.code,
		#                      constants = constants,
		#                      data = data,
		#                      inits = inits())
		# model$initializeInfo()
		# cModel <- compileNimble(model)
		# mcmcConf <- configureMCMC(cModel, onlyRW = TRUE)
		# mcmcBuild <- buildMCMC(mcmcConf)
		# compMCMC <- compileNimble(mcmcBuild)
		#
		# out.nchains <- list()
		# for(i in 1:3){
		#   compMCMC$run(niter = 3000, nburnin = 100)
		#   out.nchains[[i]] <- as.matrix(compMCMC$mvSamples)
		# }

		dat.hindcast <- do.call(rbind, out.nchains)

		if (ua.job == "ic_parameter_driver_process") {
			message("Checking convergence...")
			nodes <- colnames(out.nchains[[1]])
			gelman.keep <- numeric(length(nodes))
			for (ff in seq_along(nodes)) {
				mcmc.check <- list()
				col <- nodes[ff]
				for (c in seq_along(out.nchains)) {
					mcmc.check[[c]] <- coda::mcmc(out.nchains[[c]][, col])
				}
				gelman.keep[ff] <- try(coda::gelman.diag(
					mcmc.check,
					transform = TRUE
				)$psrf[1])
			}

			if (any(gelman.keep > 1.2)) {
				message("WARNING: Convergence not reached!")
				bad.nodes <- which(gelman.keep > 1.2)
				bad.params <- tibble(
					node = nodes[bad.nodes],
					psrf = as.numeric(gelman.keep[bad.nodes])
				) %>%
					arrange(psrf)
				print(tail(bad.params))
			} else {
				message("Convergence = TRUE")
			}
		}

		if (nrow(dat.hindcast) > 5000) {
			draws <- round(seq.int(1, nrow(dat.hindcast), length.out = 5000))
		} else {
			draws <- seq_len(nrow(dat.hindcast))
		}

		dat.draws <- dat.hindcast[draws, ]

		fileDest <- file.path(dir.save, fx.start.date)
		source("Functions/transferAnalysis.R")
		message("Running analysis...")
		transfer_analysis(
			fx.df = dat.draws,
			observations = neon.job,
			fx.dates = fx.sequence,
			model = model.job,
			ua = ua.job,
			spp = species.job,
			out.dir = fileDest
		)
		message(ua.job, " complete")
	}
}
