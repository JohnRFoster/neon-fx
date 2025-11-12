run_transfer_nimble <- function(
	cl,
	model,
	data,
	constants,
	inits,
	n.iter,
	ic,
	parameter,
	driver,
	process,
	notStatic,
	miceMNA,
	miceAndWeather,
	use.daymet
) {
	library(parallel)
	library(nimble)
	library(coda)

	source("Functions/nimble_functions.R")
	# assign("dwtnorm", dwtnorm, envir = .GlobalEnv)
	# assign("rwtnorm", rwtnorm, envir = .GlobalEnv)

	n.cores <- length(cl) # number of cores used

	export.vec <- c(
		"model",
		"constants",
		"data",
		"n.iter",
		"if_else_nimble",
		"parameter",
		"notStatic",
		"driver",
		"process",
		"ic",
		"miceAndWeather",
		"miceMNA",
		"use.daymet"
	)

	clusterExport(cl, export.vec, envir = environment())

	# export inits to clusters
	for (j in seq_along(cl)) {
		set.seed(j)
		init <- inits()
		clusterExport(cl[j], "init", envir = environment())
	}

	message("Running mcmc...")
	out <- clusterEvalQ(cl, {
		# sample on each cluster
		library(nimble)
		library(coda)

		nimbleOptions('MCMCjointlySamplePredictiveBranches' = FALSE)

		model <- nimbleModel(
			model,
			constants = constants,
			data = data,
			inits = init
		)

		cModel <- compileNimble(model)
		mcmcConf <- configureMCMC(cModel, onlyRW = TRUE)
		# if(process){
		#   for(i in 1:nrow(data$y)){
		#     for(j in 1:ncol(data$y)){
		#       if(j > 1){
		#         node <- paste0('x[', i, ', ', j, ']')
		#         mcmcConf$addSampler(node, "RW")
		#       }
		#     }
		#   }
		# }
		Rmcmc <- buildMCMC(mcmcConf)
		Cmcmc <- compileNimble(Rmcmc)
		Cmcmc$run(niter = n.iter, nburnin = n.iter / 2)
		return(as.mcmc(as.matrix(Cmcmc$mvSamples)))
	})

	out.mcmc <- as.mcmc.list(out)
	return(out.mcmc)
}
