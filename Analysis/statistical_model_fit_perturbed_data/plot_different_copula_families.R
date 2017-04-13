source('get_Rimage_data_univariate_distributions.R', local=TRUE)

# Read the summary statistics -- saved earlier for speed
store_var_list = readRDS('univariate_runs_summary_statistics.RDS')

# Read the original fit (based on un-perturbed data)
original_var_list = get_Rimage_data(
    '../statistical_model_fit/Rimages/session_univariate_distributions_FALSE_0.Rdata')

library(VineCopula)

#
# Show how changing the copula alters the relationship
#

duration_copula_type = unlist(sapply(store_var_list, f<-function(x) x$duration_season_copula$familyname))

series_simulation_files = Sys.glob("Rimages/session_series_simulation_TRUE_*.Rdata")

# Gaussian version
i1 = min(which(duration_copula_type == 'Gaussian'))
session_gauss = new.env()
load(series_simulation_files[i1], envir=session_gauss)

# Frank version
i2 = min(which(duration_copula_type == 'Frank'))
session_frank = new.env()
load(series_simulation_files[i2], envir=session_frank)

# Gumbel version
i3 = min(which(duration_copula_type == 'Rotated Gumbel 90 degrees'))
session_gum = new.env()
load(series_simulation_files[i3], envir=session_gum)

year2hours = 365.25 * 24

png('copula_family_and_seasonality.png', width=10, height=6, units='in', res=200)
par(mfrow=c(2,3))

contour_levels = seq(0.1, 1.5, by=0.05)
contour(session_gauss$duration_fit_conditional$var_season_copula, margins='unif',
    levels=contour_levels, col=rainbow(length(contour_levels+3)),
    xlab='Duration', ylab='Season', main='Gaussian copula (70 %)',
    cex.main=2)
points(session_gauss$duration_fit_conditional$copula_data, pch=19, cex=0.3)

contour_levels = seq(0.1, 1.5, by=0.05)
contour(session_frank$duration_fit_conditional$var_season_copula, margins='unif',
    levels=contour_levels, col=rainbow(length(contour_levels+3)),
    xlab='Duration', ylab='Season', main='Frank copula (29 %)',
    cex.main=2)
points(session_frank$duration_fit_conditional$copula_data, pch=19, cex=0.3)

contour_levels = seq(0.1, 1.5, by=0.05)
contour(session_gum$duration_fit_conditional$var_season_copula, margins='unif',
    levels=contour_levels, col=rainbow(length(contour_levels+3)),
    xlab='Duration', ylab='Season', main='Gumbel copula rotated (1%)',
    cex.main=2)
points(session_gum$duration_fit_conditional$copula_data, pch=19, cex=0.3)

# Boxplots of duration/season
tmn = session_gauss$synthetic_attr$startyear
month = ceiling((tmn-floor(tmn))*12)
boxplot((session_gauss$synthetic_attr$duration*year2hours) ~ month, names=month.abb, las=2,
    main='Duration Seasonality \n Gaussian copula', ylab='Hours',
    cex.main=1.5, ylim=c(0, 200))

tmn = session_frank$synthetic_attr$startyear
month = ceiling((tmn-floor(tmn))*12)
boxplot((session_frank$synthetic_attr$duration*year2hours) ~ month, names=month.abb, las=2,
    main='Duration Seasonality \n Frank copula', ylab='Hours',
    cex.main=1.5, ylim=c(0,200))

tmn = session_gum$synthetic_attr$startyear
month = ceiling((tmn-floor(tmn))*12)
boxplot((session_gum$synthetic_attr$duration*year2hours) ~ month, names=month.abb, las=2,
    main='Duration Seasonality \n Gumbel copula (rotated)', ylab='Hours',
    cex.main=1.5, ylim=c(0,200))
dev.off()


contour_levels = seq(0.1, 1.5, by=0.05)
contour(session_gauss$duration_fit_conditional$var_season_copula, margins='unif',
    levels=contour_levels, col=rainbow(length(contour_levels+3)),
    xlab='Duration', ylab='Season', main='Gaussian copula (70 %)',
    cex.main=2)
points(session_gauss$duration_fit_conditional$copula_data, pch=19, cex=0.3)

