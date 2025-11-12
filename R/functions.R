#' get the neon tick data for a given species
#'
#' @param species the species to extract; "Ixodes scapularis" or "Amblyomma americanum"
#' @export

neon_tick_data <- function(species) {
	get_data <- function(i) {
		sub.df <- data |>
			filter(plotID == plots[i]) |>
			arrange(time) |>
			ungroup()

		if (all(sub.df$processedCount == 0)) {
			return(NULL)
		}

		drags <- unique(sub.df$time) |> sort()
		n.drags <- length(drags)
		days.sequence <- seq.Date(drags[1], drags[n.drags], by = 1)
		n.days <- length(days.sequence)

		df.l <- sub.df |>
			filter(lifeStage == "Larva") |>
			group_by(time, totalSampledArea) |>
			summarise(processedCount = sum(processedCount)) |>
			mutate(lifeStage = "Larva")

		df.na <- sub.df |>
			filter(lifeStage != "Larva") |>
			select(time, processedCount, totalSampledArea, lifeStage) |>
			distinct()

		counts <- bind_rows(df.l, df.na) |>
			pivot_wider(names_from = lifeStage, values_from = processedCount) |>
			mutate(
				plotID = plots[i],
				siteID = sub.df$siteID[1],
				species = species,
				nlcd = sub.df$nlcdClass[1],
				n.drags = n.drags,
				n.days = n.days,
				count.flag = nrow(df.l) * 2 == nrow(df.na)
			)
	}

	df <- read_csv(
		"data/tickTargets.csv"
	) |>
		suppressMessages()
	data <- df |>
		filter(scientificName %in% c(species, "AAorIX"))
	# time >= "2018-01-01")

	plots <- unique(data$plotID)
	n.plots <- length(plots)

	plot.info <- tibble()
	for (i in 1:n.plots) {
		plot.info <- bind_rows(plot.info, get_data(i))
	}

	return(plot.info)
}


#' this function reads the small mammal csv and creates a
#' capture history matrix for given scale for the NEON small mammal data
#' 1 = alive - ticks present
#' 2 = alive - ticks absent
#' 3 = alive - ticks unknown
#' 4 = dead - ticks present
#' 5 = dead - ticks absent
#' 6 = dead - ticks unknown
#' 7 = not seen
#'
#' @param site the unit to extract; "HARV"
#' @param neon.smam the small mammal neon.smam frame from the csv in /Data

capture_matrix <- function(site, neon.smam) {
	neon.df <- neon.smam %>%
		filter(siteID == site)

	all.days <- neon.df %>% pull(collectDate) %>% unique()

	df <- neon.df %>%
		filter(genusName == "Peromyscus")

	mice.days <- df %>% pull(collectDate) %>% unique()

	# the states that we need are:
	alive.p <- 1 # observed mouse with tick attached
	alive.a <- 2 # observed mouse without tick attached
	alive.u <- 3 # observed mouse with unknown tick status
	dead.p <- 4 # observed mouse with tick attached
	dead.u <- 5 # observed mouse with unknown tick status
	dead.a <- 6 # observed mouse without tick attached
	unobserved <- 7 # unobserved

	# if Peromyscus was not seen on all days, need to fill days back in
	fill.days <- tibble()
	if (length(mice.days) < length(all.days)) {
		missing.days <- all.days[which(!(all.days %in% mice.days))]
		fill.days <- neon.smam %>%
			filter(collectDate %in% missing.days) %>%
			distinct(collectDate, .keep_all = TRUE) %>%
			select(collectDate, tagID) %>%
			mutate(tagID = "noCapture", state = unobserved)
	}

	# the total number of unique tags
	total.ind <- df %>%
		filter(!is.na(tagID)) %>%
		pull(tagID) %>%
		unique() %>%
		length()

	# get the state of each mouse at each trap night
	df.state <- df %>%
		mutate(
			tickOn = if_else(
				adultTicksAttached == "U" | # any life stage unknown
					nymphalTicksAttached == "U" |
					larvalTicksAttached == "U",
				3,
				0
			),
			tickOn = if_else(
				adultTicksAttached == "Y" | # any life stage observed
					nymphalTicksAttached == "Y" |
					larvalTicksAttached == "Y",
				1,
				tickOn
			),
			tickOn = if_else(
				adultTicksAttached == "N" & # all life stages not observed
					nymphalTicksAttached == "N" &
					larvalTicksAttached == "N",
				2,
				tickOn
			),
			tickOn = if_else(
				adultTicksAttached == "U" & # all life stages unknown
					nymphalTicksAttached == "U" &
					larvalTicksAttached == "U",
				3,
				tickOn
			),
			tickOn = if_else(is.na(tickOn), 3, tickOn), # NAs get unknown status
			state = if_else(animalInTrap == 1 & tickOn == 3, alive.u, 0), # observed animal with unknown tick status
			state = if_else(animalInTrap == 1 & tickOn == 2, alive.a, state), # observed animal without tick attached
			state = if_else(animalInTrap == 1 & tickOn == 1, alive.p, state), # observed animal with tick attached
			state = if_else(animalInTrap == 2 & tickOn == 3, dead.u, state), # dead animal with unknown tick status
			state = if_else(animalInTrap == 2 & tickOn == 2, dead.a, state), # dead animal without tick attached
			state = if_else(animalInTrap == 2 & tickOn == 1, dead.p, state), # dead animal with tick attached
			state = if_else(animalInTrap == 0, unobserved, state), # unobserved
			tagID = if_else(is.na(tagID) & animalInTrap == 0, "noCapture", tagID)
		) %>% # placeholder for trap nights without any captures
		select(tagID, collectDate, state) %>%
		filter(!is.na(tagID))

	# states should all have >0 designation
	if (any(df.state$state == 0)) {
		stop("All possible states unaccounted for")
	}

	# add missing days
	df.all.days <- bind_rows(df.state, fill.days)

	ch <- df.all.days %>%
		group_by(collectDate, tagID) %>%
		distinct() %>%
		summarise(state = min(state)) %>% # conflicting states get unknown designation
		ungroup() %>%
		arrange(collectDate) %>%
		pivot_wider(
			names_from = collectDate,
			values_from = state,
			values_fill = unobserved
		) %>%
		filter(tagID != "noCapture") %>%
		select(-tagID)

	# check dimensions are what they should be
	if (ncol(ch) != length(all.days)) {
		stop("Conflicting days in data vs capture matrix")
	}
	if (nrow(ch) != total.ind) {
		stop("Conflicting individuals in data vs capture matrix")
	}

	# make sure mice that are recorded dead stay dead!
	# all days after found dead should be unobserved
	# assume that if there is a case of a zombie mouse,
	# the dead recording is true and the following capture
	# is an error in reading the tag number
	for (z in 1:nrow(ch)) {
		if (any(ch[z, ] %in% c(dead.u, dead.a, dead.p))) {
			# print(z)
			day.dead <- which(ch[z, ] %in% c(dead.u, dead.a, dead.p))
			if (day.dead < ncol(ch)) {
				if (!all(ch[z, (day.dead + 1):ncol(ch)] == unobserved)) {
					ch[z, (day.dead + 1):ncol(ch)] <- unobserved
				}
			}
		}
	}

	# dates and deltas
	# all days trapped
	capture.dates <- df.all.days %>%
		pull(collectDate) %>%
		unique() %>%
		sort()

	# every day in time series
	every.day <- seq.Date(first(capture.dates), last(capture.dates), by = 1)

	# the index where each capture occasion happens in time series
	capture.index <- which(every.day %in% capture.dates)

	# the number of days between capture occasions
	delta.days <- diff.Date(capture.dates) %>% as.numeric()

	return(list(
		ch = as.matrix(ch),
		alive.p = alive.p,
		alive.a = alive.a,
		alive.u = alive.u,
		dead.p = dead.p,
		dead.u = dead.u,
		dead.a = dead.a,
		unobserved = unobserved,
		capture.dates = capture.dates,
		every.day = every.day,
		capture.index = capture.index,
		delta.days = delta.days
	))
}

#' Known states function
#'
#' This function fills in the capture history matrix where individuals are known but not observed
#' @param ch Capture history matrix
#' @export
#' @examples known_states(ch)

known_states <- function(ch) {
	state <- ch
	for (i in 1:dim(ch)[1]) {
		n1 <- min(which(ch[i, ] != 0))
		n2 <- max(which(ch[i, ] != 0))
		if (n2 > n1) {
			state[i, n1:n2] <- 1
		}
		#state[i, n1:n2] <- 1
	}
	return(state)
}

#' Capture History Matrix Function for the Cary Mouse Data
#'
#' This function reads the raw csv file containing white footef mice capture data
#' @param path File path to csv, Default = "" and assumes csv is in working directory
#' @param grid What grid do you want to create matrix for? One of Green Control, Henry Control, Tea Control, Green Experimental, Henry Experimental, Tea Experimental
#' @examples ch.cary("Green Control")
#' @export

ch_cary <- function(grid, path = "") {
	file <- "Data/Cary_mouse.csv"

	dat <- read.csv(
		paste(path, file, sep = ""),
		na.strings = c(c("", " ", "      "), "NA")
	) # read in data

	smam <- dat[, c(
		"Grid",
		"Full.Date.1",
		"Full.Date.2",
		"Day.1",
		"Day.2",
		"Tag..",
		"Fate"
	)]

	alive <- c(1, 2) # codes for alive individuals
	#smam <- smam %>% filter(Fate %in% alive)            # extract only live individuals
	smam <- subset(smam, Fate == 1 | Fate == 2)

	smam[4:5] <- as.integer(!is.na(smam[4:5])) # converts trap histories to 1 or 0

	smam$Full.Date.1 <- as.Date(smam$Full.Date.1) # convert factors to dates
	smam$Full.Date.2 <- as.Date(smam$Full.Date.2) # convert factors to dates

	m <- subset(smam, Grid == grid)
	m <- m[, c("Tag..", "Day.1", "Day.2", "Full.Date.1", "Full.Date.2")]
	m <- subset(m, !is.na(Tag..))

	m$Tag.. <- as.character(m$Tag..)
	m.ls <- split(m, m$Full.Date.1) # split on first capture date for sampling occasion
	for (i in 1:length(m.ls)) {
		m.ls[[i]] <- m.ls[[i]][c(-4, -5)]
	}

	day.1 <- as.character(unique(m$Full.Date.1)) # 1st capture date of sampling occasion
	day.2 <- as.character(unique(m$Full.Date.2)) # 2nd capture date of sampling occasion
	days <- c(rbind(day.1, day.2)) # vector of unique trapping days (for colnames)

	ch.base <- merge(m.ls[[1]], m.ls[[2]], by = "Tag..", all = TRUE)
	l.m <- length(m.ls) * 1
	for (i in 3:l.m) {
		# loop through the rest
		g <- as.data.frame(m.ls[[i]])
		ch.base <- merge(ch.base, g, by = "Tag..", all = TRUE)
	}
	ch.base <- as.matrix(ch.base[, -1]) # convert all NAs to 0
	for (i in 1:nrow(ch.base)) {
		for (t in 1:ncol(ch.base)) {
			if (is.na(ch.base[i, t])) {
				ch.base[i, t] <- 0
			}
			as.numeric(ch.base[i, t])
		}
	}
	colnames(ch.base) <- days
	return(ch.base)
}

#' get the minimum number alive (mna) for jags input for a specific site

mna_jags <- function(site.run, return.mean = FALSE) {
	if (!grepl("Control", site.run)) {
		site.run <- paste(site.run, "Control")
	}
	ch <- suppressWarnings(ch_cary(site.run))
	ks <- known_states(ch)
	mna <- apply(ks, 2, sum)

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

	# tick data
	dat <- read.csv("data/tick_cleaned") # tick data
	tick <- dat[, c("Grid", "DATE", "n_larvae", "n_nymphs", "n_adults")]
	tick <- subset(tick, Grid == site.run)
	tick$DATE <- as.Date(tick$DATE)

	# match tick dates to mice dates
	start.tick <- which(tick$DATE[1] == mice.seq)
	end.tick <- which(tick$DATE[length(tick$DATE)] == mice.seq)

	# index mice estimates - convert sd to prec
	mna.for.jags <- mna.all.days[start.tick:end.tick]

	# center and scale
	mna.scaled <- scale(mna.for.jags)

	if (return.mean) {
		return(
			list(
				mna = mna.scaled,
				mean = mean(mna.for.jags),
				sd = sd(mna.for.jags)
			)
		)
	} else {
		return(mna.scaled)
	}
}


scale_met_forecast <- function() {
	met <- read.csv("data/Cary_Met_Data_Daily.csv")
	met <- met %>%
		filter(DATE >= "05-02-1995") %>%
		filter(DATE <= "08-17-2005") %>%
		select(c("MAX_TEMP", "MAX_RH", "MIN_RH", "TOT_PREC"))

	met.means <- apply(met, 2, mean, na.rm = TRUE)
	met.sd <- apply(met, 2, sd, na.rm = TRUE)

	met.scale <- apply(met, 2, scale)
	scale.max <- c(
		10000,
		(100 - met.means["MAX_RH"]) / met.sd["MAX_RH"],
		(100 - met.means["MIN_RH"]) / met.sd["MIN_RH"],
		10000
	)

	scale.min <- c(
		-10000,
		(0 - met.means["MAX_RH"]) / met.sd["MAX_RH"],
		(0 - met.means["MIN_RH"]) / met.sd["MIN_RH"],
		(0 - met.means["TOT_PREC"]) / met.sd["TOT_PREC"]
	)

	names(scale.max) <- c("MAX_TEMP", "MAX_RH", "MIN_RH", "TOT_PREC")
	names(scale.min) <- c("MAX_TEMP", "MAX_RH", "MIN_RH", "TOT_PREC")

	return(list(
		means = met.means,
		sds = met.sd,
		scale.max = scale.max,
		scale.min = scale.min
	))
}

#' Score forecasts
#'
#' CRPS is used in the main analysis, but other metrics are calculated as well
#'
#' @param df.pred.obs dataframe with forecast samples and observations
#' @param nmcmc number of MCMC samples

score <- function(df.pred.obs, nmcmc) {
	pred.obs <- df.pred.obs
	ls.vec <- unique(pred.obs$lifeStage)

	# crps scores and other metrics
	start.date <- min(pred.obs$time)
	obs.plots <- unique(pred.obs$plotID)
	all.scores <- tibble()
	for (p in seq_along(obs.plots)) {
		plot.subset <- pred.obs %>% filter(plotID == obs.plots[p])
		obs.dates <- unique(plot.subset$time)
		for (i in seq_along(ls.vec)) {
			ls.subset <- plot.subset %>%
				arrange(time) %>%
				filter(lifeStage == ls.vec[i])
			ls.observed <- ls.subset %>%
				select(time, observed) %>%
				distinct() %>%
				pull(observed)
			ls.predicted <- ls.subset %>%
				select(time, forecast) %>%
				group_by(time) %>%
				mutate(iter = 1:nmcmc) %>%
				pivot_wider(names_from = time, values_from = forecast) %>%
				select(-iter) %>%
				as.matrix() %>%
				t()

			rmse <- function(obs, pred) {
				if (length(obs) == 1) {
					n <- length(pred)
					rmse <- sqrt((1 / n) * sum((pred - obs)^2))
				} else {
					n <- ncol(pred)
					rmse <- obs
					for (i in 1:length(obs)) {
						rmse[i] <- sqrt((1 / n) * sum((pred[i, ] - obs[i])^2))
					}
				}
				return(rmse)
			}

			bayes_p <- function(obs, pred) {
				bayes.p.val <- rep(NA, nrow(pred))
				for (gg in 1:nrow(pred)) {
					cdf <- ecdf(pred[gg, ])
					bayes.p.val[gg] <- cdf(obs[gg])
				}
				return(bayes.p.val)
			}

			scores <- tibble(
				lifeStage = ls.vec[i],
				time = unique(ls.subset$time),
				plotID = obs.plots[p],
				horizon = as.numeric(time - start.date),
				percentBias = scoringutils::bias_sample(ls.observed, ls.predicted),
				crps = scoringutils::crps_sample(ls.observed, ls.predicted),
				rmse = rmse(ls.observed, ls.predicted),
				bayesP = bayes_p(ls.observed, ls.predicted)
			)

			all.scores <- bind_rows(all.scores, scores)
		}
	}
	return(all.scores)
}


#' Run the analysis on tick forecasts at NEON
#'
#' @param fx.df the forecast dataframe from run_transfer_nimble.R
#' @param observations dataframe of observations from the NEON site
#' @param fx.dates vector of dates for the forecast
#' @param model the model used for the forecast
#' @param ua the uncertainty type(s) used in the forecast
#' @param spp the species being forecasted
#' @param out.dir the output directory to save results

transfer_analysis <- function(
	fx.df,
	observations,
	fx.dates,
	model,
	ua,
	spp,
	out.dir
) {
	nmcmc <- nrow(fx.df)
	ls.tb <- tibble(
		lifeStage = c("Larva", "Dormant", "Nymph", "Adult"),
		ls.index = 1:4
	)

	time.tb <- tibble(
		time = fx.dates,
		time.index = 1:length(fx.dates)
	)

	fx.tb <- fx.df %>%
		as_tibble() %>%
		pivot_longer(cols = everything(), names_to = "node")

	state.nodes <- fx.tb %>%
		filter(grepl("x[", node, fixed = TRUE)) %>%
		pull(node) %>%
		unique()

	weather.nodes <- fx.tb %>%
		filter(
			grepl("x1", node, fixed = TRUE) |
				grepl("x2", node, fixed = TRUE) |
				grepl("x3", node, fixed = TRUE) |
				grepl("x4", node, fixed = TRUE)
		) %>%
		pull(node) %>%
		unique()

	weather <- length(weather.nodes) > 1

	states <- fx.tb %>%
		filter(node %in% state.nodes) %>%
		mutate(
			time.index = as.numeric(str_extract(node, "(?<=\\s)\\d*(?=\\])")),
			ls.index = as.numeric(str_extract(node, "(?<=x\\[)\\d(?=\\,)"))
		)

	join1 <- left_join(states, ls.tb, by = "ls.index")
	join2 <- left_join(join1, time.tb, by = "time.index")

	site <- observations$siteID %>% unique()

	fx.samples <- join2 %>%
		select(value, lifeStage, time) %>%
		mutate(siteID = site)

	data.site <- observations %>%
		filter(time %in% fx.sequence) %>%
		select(-n.drags, -n.days, -count.flag) %>%
		pivot_longer(
			cols = c("Larva", "Nymph", "Adult"),
			names_to = "lifeStage",
			values_to = "observed"
		)

	pred.obs <- tibble() # pred.obs will only have forecasts for the dates tick drags occured
	plots <- data.site$plotID %>% unique()
	for (p in seq_along(plots)) {
		data.sub <- data.site %>%
			filter(plotID == plots[p])
		plot.time <- unique(data.sub$time)
		fx.sub <- fx.samples %>%
			filter(time %in% plot.time, lifeStage != "Dormant")
		pred.obs.plot <- left_join(
			fx.sub,
			data.sub,
			by = c("time", "lifeStage", "siteID")
		)
		pred.obs <- bind_rows(pred.obs, pred.obs.plot)
	}

	fx.data <- pred.obs %>%
		mutate(
			forecast = value / 450 * totalSampledArea,
			start.date = start.date,
			ua = ua,
			model = model
		)

	fx.quantiles <- fx.data %>%
		group_by(
			lifeStage,
			time,
			siteID,
			start.date,
			ua,
			model,
			species,
			totalSampledArea,
			plotID,
			nlcd,
			observed
		) %>%
		summarise(
			lower95 = quantile(forecast, 0.025),
			lower75 = quantile(forecast, 0.125),
			median = median(forecast),
			mean = mean(forecast),
			upper75 = quantile(forecast, 0.875),
			upper95 = quantile(forecast, 0.975),
			variance = var(forecast)
		)

	# get scores
	scores <- score(fx.data, nmcmc) %>%
		mutate(siteID = site, ua = ua, species = spp, model = model)

	fx.out <- left_join(
		fx.quantiles,
		scores,
		by = c("lifeStage", "time", "species", "plotID", "model", "siteID", "ua")
	)

	# parameters
	param.samples <- fx.tb %>%
		filter(!node %in% state.nodes, !node %in% weather.nodes) %>%
		mutate(siteID = site, ua = ua, model = model)

	param.quants <- param.samples %>%
		group_by(node) %>%
		summarise(
			lower95 = quantile(value, 0.025),
			lower75 = quantile(value, 0.125),
			median = median(value),
			mean = mean(value),
			upper75 = quantile(value, 0.875),
			upper95 = quantile(value, 0.975),
			variance = var(value)
		) %>%
		ungroup() %>%
		mutate(
			siteID = site,
			ua = ua,
			species = spp,
			start.date == start.date,
			model = model
		)

	if (!dir.exists(out.dir)) {
		dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
	}

	# save
	message("  Writing files to ", out.dir)
	write_csv(ungroup(fx.out), file.path(out.dir, "fxQuantScore.csv"))
	write_csv(ungroup(fx.samples), file.path(out.dir, "stateSamples.csv"))
	write_csv(ungroup(param.samples), file.path(out.dir, "parameterSamples.csv"))
	write_csv(ungroup(param.quants), file.path(out.dir, "parameterSummary.csv"))

	# weather - if used
	if (weather) {
		daymet <- fx.tb %>%
			filter(node %in% weather.nodes) %>%
			group_by(node) %>%
			summarise(
				lower95 = quantile(value, 0.025),
				lower75 = quantile(value, 0.125),
				median = median(value),
				mean = mean(value),
				upper75 = quantile(value, 0.875),
				upper95 = quantile(value, 0.975),
				variance = var(value)
			) %>%
			ungroup() %>%
			mutate(time.index = as.numeric(str_extract(node, "(?<=\\[)\\d*(?=\\])")))

		daymet.time <- left_join(time.tb, daymet, by = "time.index") %>%
			select(-time.index) %>%
			mutate(
				siteID = site,
				ua = ua,
				species = spp,
				start.date == start.date,
				model = model
			)

		write_csv(daymet.time, file.path(out.dir, "weatherSummary.csv"))
	}
}
