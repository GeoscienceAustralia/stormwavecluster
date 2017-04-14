
## ##########################################
## ## Figure 1 -- regional setting ##
## 
## local_env = new.env()
## source('../DATA/Shand_et_al_database/figure/', chdir=TRUE, local=local_env)
## rm(local_env)
## 
## #######################################




##########################################
## Figure with seasonal trends ##

#load('Session_end.Rdata')
load('Rimages/session_series_simulation_FALSE_0.Rdata')

# Count events per month each year
month_num = as.numeric(format(event_statistics$time, '%m'))
years = floor(event_statistics$startyear)

# Strip off first and last years for comparison
month_year_combo = expand.grid(1:12, 1986:2015)
month_year_eventcount = month_year_combo[,1] * NA
for(i in 1:length(month_year_combo[,1])){
    month_year_eventcount[i] = sum( 
        (month_num == month_year_combo[i,1]) & 
        (years == month_year_combo[i,2]))
}

png('monthly_statistics.png', width=10, height=10, units='in', res=200)
par(mfrow=c(3,2))
par(mar=c(4,3,2,1))

boxplot(month_year_eventcount ~ month_year_combo[,1], xlab='', 
    ylab='Number of events', main='', names=month.abb,
    col='grey', cex.main=2.0, cex.axis=1.3, las=2)
#mtext('Month', side=1, line=2)
    #grid(col='brown')
title('Number of Storms', cex.main=2.0)

titles = c('Duration (hours)', 'Hsig (m)', 'Period (s)', 'Direction (degrees)', 'Tidal Residual (m)')

for(i in 1:5){
    boxplot(event_statistics[,i] ~ month_num, xlab='', 
        ylab=names(event_statistics)[i], col='grey', 
        names=month.abb, cex.axis=1.3, las=2)
    title(main = titles[i], cex.main=2.0)
}
dev.off()

###################################################################
#
# Pairwise scatterplots
#
###################################################################
nice_pairs_paper<-function( mydata, labels, extra_data=NULL, mydata_col='blue', add_signif_stars=TRUE, copula_type=NULL){
    # Prettier pairs plot
    pairs(mydata, labels,
        pch='.', cex=3, 
        upper.panel=function(x,y,...){ 
            points(x,y, col=mydata_col, ...); 
            grid(col='orange')
            if(!is.null(extra_data)){
                # Hack to identify the data columns
                icol = get('i', pos=parent.frame(n=2))
                jcol = get('j', pos=parent.frame(n=2))
                points(extra_data[,jcol], extra_data[,icol], col='black', pch='.', cex=2)
            }

        }, 
        diag.panel=DU$panel.hist,
        lower.panel=function(x,y,...){ 
            points(x,y, col=0); 

            spearman_cortest = try(cor.test(x, y, method='s'))
            if(class(spearman_cortest) != 'try-error'){
                spearman_cor = spearman_cortest$estimate
                title_word = as.character(round(spearman_cor, 3))
                font_size = 2.0

                if(add_signif_stars){
                    if(spearman_cortest$p.value < 0.05){
                        title_word = paste0(title_word, '*')
                        font_size = font_size + 0.5
                    }
                    if(spearman_cortest$p.value < 0.005){
                        title_word = paste0(title_word, '*')
                        font_size = font_size + 0.5
                    }
                }else{
                    font_size = font_size + 0.5 
                }

                if(!is.null(extra_data)){
                    # Hack to identify the data columns
                    icol = get('i', pos=parent.frame(n=2))
                    jcol = get('j', pos=parent.frame(n=2))
                    #points(extra_data[,jcol], extra_data[,icol], col='black', pch='.', cex=2)
                    extra_cortest = try(cor.test(extra_data[,jcol], extra_data[,icol], method='s'))
                    
                    title_word = paste0(title_word, '\n(', as.character(round(extra_cortest$estimate,3)), 
                        ') \n', copula_type[icol, jcol])
                }

                title(main=title_word, line=-8, col.main='black', cex.main=font_size) 
            }
        })
}

plotvar = c('hsig', 'duration', 'tideResid', 'tp1', 'dir') #c('duration', 'hsig', 'tp1', 'dir', 'tideResid')
plotvar_names = c('Hsig (m)', 'D (hours)', 'R (m)', 'T (s)', expression(theta ~ (degrees) )) #c('D', 'Hsig', 'T', 'theta', 'h')
event_stats_tmp = event_statistics[,plotvar]
names(event_stats_tmp) = plotvar_names
pdf('data_pairs.pdf', width=10, height=10)
nice_pairs_paper(event_stats_tmp, labels=plotvar_names)
dev.off()

### Plot the model results
#fit_env = new.env()
#load('Rimages/session_series_simulation_FALSE_0.Rdata', envir=fit_env)
#library(VineCopula)
#
## Make matrix recording copula types for plot
#copula_type = matrix(
#    BiCopName(fit_env$copula_model$copula_fit_mle$RVM$family, short=FALSE), 
#    ncol=5, nrow=5)
#copula_type = t(copula_type[5:1,5:1])
#copula_type[4,] = paste(copula_type[4,], '*', sep="")
#copula_type[5,4] = paste(copula_type[5,4], '*', sep="")
#copula_type[upper.tri(copula_type, diag=TRUE)] = ""
#copula_type = gsub('Independence', 'Indep', copula_type)
#copula_type = gsub('Survival Gumbel', '180_Gumbel', copula_type)
#
#pdf('model_data_pairs.pdf', width=10, height=10)
#synthetic_attr_tmp = fit_env$synthetic_attr[1:5000, plotvar]
#synthetic_attr_tmp$duration = synthetic_attr_tmp$duration * 365.25 * 24
#names(synthetic_attr_tmp) = plotvar_names
#nice_pairs_paper(synthetic_attr_tmp, labels=plotvar_names, #names(synthetic_attr_tmp), 
#    extra_data = event_stats_tmp, mydata_col='green', add_signif_stars=FALSE,
#    copula_type = copula_type)
#dev.off()
#
#
