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
#' capture history matrix for given scale
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
