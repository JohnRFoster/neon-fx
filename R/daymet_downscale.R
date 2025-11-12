# =========================================================================== #
# functions for extracting and working with daymet
# daymet data has already been downloaded with R/0_intakeDayMet.R
# =========================================================================== #

library(tidyverse)

#' function that calculates cumulative growing degree days for each plot
#' @param site the site being modeled
#' @param org either "tick" or "smam"
daymet_cumGDD <- function(site) {
	df.all <- read_csv(
		"data/daymetSite_maxTemperature.csv"
	)
	df <- df.all %>%
		filter(siteID == site) %>%
		group_by(year) %>%
		mutate(
			growingDegree = if_else(maxTemperature > 10, maxTemperature - 10, 0),
			cumGDD = cumsum(growingDegree)
		) %>%
		select(Date, siteID, cumGDD, year)

	return(df)
}

## max temperature ==================================================================
daymet_temp <- function(site, minimum) {
	if (minimum) {
		df <- read_csv(
			"data/daymetSite_minTemperature.csv"
		)
		neon.col <- "tempTripleMinimum"
		daymet.col <- "minTemperature"
	} else {
		df <- read_csv(
			"data/daymetSite_maxTemperature.csv"
		)
		neon.col <- "tempTripleMaximum"
		daymet.col <- "maxTemperature"
	}

	df.temp <- df %>%
		filter(siteID == site) %>%
		group_by(yday) %>%
		select(-tile)

	neon.temp <- read_csv(
		"data/airTempDaily.csv"
	)
	neon.sub <- neon.temp %>%
		filter(siteID == site) %>%
		mutate(yday = yday(Date))

	neon.doy <- neon.sub %>%
		group_by(yday) %>%
		summarise(muNeon = mean(.data[[neon.col]])) %>%
		ungroup()

	daymet.doy <- df.temp %>%
		group_by(yday) %>%
		summarise(muDaymet = mean(.data[[daymet.col]])) %>%
		ungroup()

	tempbias <- left_join(neon.doy, daymet.doy, by = "yday") %>%
		mutate(tempBias = muNeon - muDaymet) %>%
		select(yday, tempBias)

	daymet.temp.bias <- left_join(df.temp, tempbias, by = "yday") %>%
		mutate(TempCorrect = .data[[daymet.col]] + tempBias)

	if (minimum) {
		daymet.temp.bias <- daymet.temp.bias %>%
			rename(minTempCorrect = TempCorrect)
	} else {
		daymet.temp.bias <- daymet.temp.bias %>%
			rename(maxTempCorrect = TempCorrect)
	}
	return(daymet.temp.bias)
}

# ndf <- neon.sub %>%
#   filter(Date >= "2018-01-01",
#          Date < "2021-01-01") %>%
#   select(Date, siteID, tempTripleMaximum)
# ddf <- daymet.temp.bias %>%
#   ungroup() %>%
#   filter(Date >= "2018-01-01",
#          Date < "2021-01-01") %>%
#   select(Date, siteID, maxTempCorrect)
#
# gg.max.temp <- left_join(ndf, ddf, by = c("siteID", "Date")) %>%
#   ggplot()+
#   aes(x = tempTripleMaximum, y = maxTempCorrect) +
#   geom_point() +
#   geom_abline() +
#   labs(title = "Daily Minimum Temperature (Deg. C)",
#        x = "NEON",
#        y = "Daymet Corrected") +
#   theme_pubr()

## relative humidity ==========================================================================

daymet_rh <- function(site) {
	df <- read_csv(
		"data/daymetSite_vaporPressure.csv"
	)
	df.vpd <- df %>%
		filter(siteID == site) %>%
		select(-tile)

	df <- read_csv(
		"data/daymetSite_maxTemperature.csv"
	)
	df.temp <- df %>%
		filter(siteID == site) %>%
		select(-tile)

	df.join <- left_join(
		df.vpd,
		df.temp,
		by = c(
			"siteID",
			"latitude",
			"longitude",
			"altitude",
			"year",
			"yday",
			"Date"
		)
	)
	df.dew <- df.join %>%
		ungroup() %>%
		mutate(rh = plantecophys::VPDtoRH(vaporPressure / 1000, maxTemperature))

	neon.temp <- read_csv(
		"data/RelativeHumidityDaily.csv"
	)
	neon.sub <- neon.temp %>%
		filter(siteID == site) %>%
		mutate(yday = yday(Date)) %>%
		select(Date, yday, RHMaximum, RHMinimum)

	neon.doy <- neon.sub %>%
		group_by(yday) %>%
		summarise(muRHmax = mean(RHMaximum), muRHmin = mean(RHMinimum))
	daymet.doy <- df.dew %>%
		group_by(yday) %>%
		summarise(muDaymet = mean(rh))

	df.join <- left_join(neon.doy, daymet.doy, by = "yday") %>%
		mutate(biasMax = muRHmax - muDaymet, biasMin = muRHmin - muDaymet)

	daymet.temp.bias <- left_join(df.dew, df.join, by = "yday") %>%
		mutate(
			maxRHCorrect = pmin(rh + biasMax, 100),
			minRHCorrect = pmin(rh + biasMin, 100)
		) %>%
		select(-maxTemperature, -vaporPressure)
	return(daymet.temp.bias)
}

# ndf <- neon.sub %>%
#   filter(Date >= "2018-01-01",
#          Date < "2021-01-01") %>%
#   select(-yday)
# ddf <- daymet.temp.bias %>%
#   ungroup() %>%
#   filter(Date >= "2018-01-01",
#          Date < "2021-01-01") %>%
#   select(Date, maxRHCorrect, minRHCorrect)
#
# rh.df <- left_join(ndf, ddf, by = c("Date"))
#
# rh.df %>%
#   ggplot() +
#   aes(x = `Date`) +
#   geom_line(aes(y = RHMaximum, linetype = "NEON")) +
#   geom_line(aes(y = maxRHCorrect, linetype = "Downscale")) +
#   theme_pubr()
#
# gg.max.rh <- rh.df %>%
#   ggplot() +
#   # aes(x = Date) +
#   aes(x = RHMaximum, y = maxRHCorrect) +
#   geom_point() +
#   # geom_abline() +
#   labs(title = "Daily Maximim RH (%)",
#        x = "NEON",
#        y = "Daymet Corrected") +
#   theme_pubr()
# gg.max.rh
# gg.min.rh <- rh.df %>%
#   ggplot()+
#   aes(x = RHMinimum, y = minRHCorrect) +
#   geom_point() +
#   geom_abline() +
#   labs(title = "Daily Maximim RH (%)",
#        x = "NEON",
#        y = "Daymet Corrected") +
#   theme_pubr()
# gg.min.rh

## Precipitation ==========================================================================

daymet_precip <- function(site) {
	neon.precip <- read_csv(
		"data/precipDaily.csv"
	)
	neon.sub <- neon.precip %>%
		filter(siteID == site) %>%
		mutate(year = year(Date)) %>%
		group_by(year) %>%
		summarise(sum.precip = sum(priPrecipTotal)) %>%
		pull(sum.precip) %>%
		mean()

	df <- read_csv(
		"data/daymetSite_precipitation.csv"
	)

	df.p <- df %>%
		filter(siteID == site)
	return(df.p)
}

# df %>%
#   filter(siteID == site) %>%
#   mutate(year = year(Date)) %>%
#   group_by(year) %>%
#   summarise(sum.precip = sum(precipitation)) %>%
#   pull(sum.precip) %>%
#   mean()
