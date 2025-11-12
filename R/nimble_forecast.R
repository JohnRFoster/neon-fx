library(nimble)


if_else_nimble <- nimbleFunction(
	run = function(
		condition = integer(0),
		valueIf = double(0),
		valueElse = double(0)
	) {
		returnType(double(0))
		if (condition == TRUE) {
			return(valueIf)
		} else {
			return(valueElse)
		}
	}
)

model.code <- nimbleCode({
	### priors
	if (parameter) {
		phi.l.mu ~ dnorm(pr.phi.l[1], tau = pr.phi.l[2])
		phi.n.mu ~ dnorm(pr.phi.n[1], tau = pr.phi.n[2])
		phi.a.mu ~ dnorm(pr.phi.a[1], tau = pr.phi.a[2])
		theta.ln ~ dnorm(pr.theta.l2n[1], tau = pr.theta.l2n[2])
		theta.na ~ dnorm(pr.theta.n2a[1], tau = pr.theta.n2a[2])

		if (notStatic) {
			for (j in 1:n.beta) {
				beta[j] ~ dnorm(pr.beta[j, 1], tau = pr.beta[j, 2])
			}
		}
	} else {
		phi.l.mu <- pr.phi.l[1]
		phi.n.mu <- pr.phi.n[1]
		phi.a.mu <- pr.phi.a[1]
		theta.ln <- pr.theta.l2n[1]
		theta.na <- pr.theta.n2a[1]

		if (notStatic) {
			for (j in 1:n.beta) {
				beta[j] <- pr.beta[j, 1]
			}
		}
	}

	if (driver) {
		tau.temp ~ dexp(1)
		tau.maxrh ~ dexp(1)
		tau.minrh ~ dexp(1)
		tau.precip ~ dexp(1)
		tau.cgdd ~ dexp(1)
	}

	if (process) {
		for (i in 1:ns) {
			sig[i] ~ dinvgamma(pr.sig[i, 1], pr.sig[i, 2])
		}

		### precision priors with process error
		OMEGA[1, 1] <- sig[1]
		OMEGA[2, 2] <- sig[2]
		OMEGA[3, 3] <- sig[3]
		OMEGA[4, 4] <- sig[4]
	} else {
		### precision priors without process error
		OMEGA[1, 1] <- 1E-10
		OMEGA[2, 2] <- 1E-10
		OMEGA[3, 3] <- 1E-10
		OMEGA[4, 4] <- 1E-10
	}

	OMEGA[1, 2] <- 0
	OMEGA[1, 3] <- 0
	OMEGA[1, 4] <- 0
	OMEGA[2, 1] <- 0
	OMEGA[2, 3] <- 0
	OMEGA[2, 4] <- 0
	OMEGA[3, 1] <- 0
	OMEGA[3, 2] <- 0
	OMEGA[3, 4] <- 0
	OMEGA[4, 1] <- 0
	OMEGA[4, 2] <- 0
	OMEGA[4, 3] <- 0

	# Cholesky decomposition
	Ochol[1:ns, 1:ns] <- chol(OMEGA[1:ns, 1:ns])

	### first latent process
	if (ic) {
		for (i in 1:4) {
			x[i, 1] ~ T(dnorm(IC[i, 1], tau = IC[i, 2]), 0, Inf)
		}
	} else {
		for (i in 1:4) {
			x[i, 1] <- IC[i, 1]
		}
	}

	### define parameters
	for (t in 1:horizon) {
		# loop over every day in time series

		if (miceAndWeather) {
			# mice[t] ~ dnorm(mice.mean[t], tau = mice.prec[t])
			logit(l2n[t]) <- theta.ln + beta[13] * mice[t]
			logit(n2a[t]) <- theta.na + beta[14] * mice[t]
		} else if (miceMNA) {
			logit(l2n[t]) <- theta.ln + beta[1] * mice[t]
			logit(n2a[t]) <- theta.na + beta[2] * mice[t]
		} else {
			logit(l2n[t]) <- theta.ln
			logit(n2a[t]) <- theta.na
		}

		if (driver) {
			cgdd[t] ~ dnorm(gdd[t], tau = tau.cgdd)
			gdd[t] ~ dunif(0, max.cgdd)
		} else {
			gdd[t] <- cgdd[t]
		}

		theta.n2a[t] <- if_else_nimble(
			(gdd[t] <= 1000) | (gdd[t] >= 2500),
			n2a[t],
			0
		)
		lambda[t] <- if_else_nimble(
			(gdd[t] >= 1400) & (gdd[t] <= 2500),
			repro.mu,
			0
		)
		l2n.quest[t] <- if_else_nimble((gdd[t] >= 400) & (gdd[t] <= 2500), 1, 0)

		if (use.daymet) {
			if (driver) {
				maxtemp[t] ~ dnorm(x1[t], tau = tau.temp)
				maxrh[t] ~ dnorm(x2[t], tau = tau.maxrh)
				minrh[t] ~ dnorm(x3[t], tau = tau.minrh)
				precip[t] ~ dnorm(x4[t], tau = tau.precip)

				x1[t] ~ dnorm(0, 1)
				x2[t] ~ dnorm(0, 1)
				x3[t] ~ dnorm(0, 1)
				x4[t] ~ dnorm(0, 1)
			} else {
				x1[t] <- maxtemp[t]
				x2[t] <- maxrh[t]
				x3[t] <- minrh[t]
				x4[t] <- precip[t]
			}

			logit(phi.l[t]) <- phi.l.mu +
				beta[1] * x1[t] +
				beta[2] * x2[t] +
				beta[3] * x3[t] +
				beta[4] * x4[t]

			logit(phi.n[t]) <- phi.n.mu +
				beta[5] * x1[t] +
				beta[6] * x2[t] +
				beta[7] * x3[t] +
				beta[8] * x4[t]

			logit(phi.a[t]) <- phi.a.mu +
				beta[9] * x1[t] +
				beta[10] * x2[t] +
				beta[11] * x3[t] +
				beta[12] * x4[t]
		} else {
			logit(phi.l[t]) <- phi.l.mu
			logit(phi.n[t]) <- phi.n.mu
			logit(phi.a[t]) <- phi.a.mu
		}

		A[1, 1, t] <- phi.l[t] * (1 - l2n[t])
		A[1, 2, t] <- 0
		A[1, 3, t] <- 0
		A[1, 4, t] <- lambda[t]
		A[2, 1, t] <- phi.l[t] * l2n[t]
		A[2, 2, t] <- 1 - l2n.quest[t]
		A[2, 3, t] <- 0
		A[2, 4, t] <- 0
		A[3, 1, t] <- 0
		A[3, 2, t] <- l2n.quest[t]
		A[3, 3, t] <- phi.n[t] * (1 - theta.n2a[t])
		A[3, 4, t] <- 0
		A[4, 1, t] <- 0
		A[4, 2, t] <- 0
		A[4, 3, t] <- phi.n[t] * theta.n2a[t]
		A[4, 4, t] <- phi.a[t]

		### Data Model ###
		for (p in 1:n.plots) {
			dx[1, t, p] <- x[1, t] / 450 * area[t, p]
			dx[2, t, p] <- x[3, t] / 450 * area[t, p]
			dx[3, t, p] <- x[4, t] / 450 * area[t, p]
			y[1, t, p] ~ dpois(dx[1, t, p])
			y[3, t, p] ~ dpois(dx[2, t, p])
			y[4, t, p] ~ dpois(dx[3, t, p])
		}
	}

	for (t in 2:horizon) {
		# process error

		# expected number questing
		Ex[1:ns, t] <- A[1:ns, 1:ns, t - 1] %*% x[1:ns, t - 1]

		x[1:ns, t] ~
			dmnorm(mean = Ex[1:ns, t], cholesky = Ochol[1:ns, 1:ns], prec_param = 0)

		for (c in 1:ns) {
			xind[c, t] ~ dconstraint(x[c, t] >= 0)
		}

		# x[1,t] <- max(px[1,t], 0)
		# x[2,t] <- max(px[2,t], 0)
		# x[3,t] <- max(px[3,t], 0)
		# x[4,t] <- max(px[4,t], 0)
	}
})
