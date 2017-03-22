# Sensitivity of the multivariate vine copula fit to random perturbations of the data
---------------------------------------------------------------------------------------

The code below investigates how the automatically selected copula model varies
due to the random data perturbation.

**Read the model results, and do some basic checks**

```r
library(VineCopula)
source('get_Rimage_data_vine_copula.R', local=TRUE)

# Read the summary statistics -- saved earlier for speed
store_var_list = readRDS('vine_copula_runs_summary_statistics.RDS')

# Read the original fit (based on un-perturbed data)
original_var_list = get_Rimage_data_vine_copula(
    '../statistical_model_fit/Rimages/session_vine_copula_FALSE_0.Rdata')

# Check that all the perturbed data sessions do jittering
stopifnot(all(sapply(store_var_list, f<-function(x) x$break_ties_with_jitter)))
# Check the original fit does not do jittering
stopifnot(original_var_list$break_ties_with_jitter == FALSE)

# Check that the event_statistics is unique in every session [i.e. the perturbed
# sessions really do randomly perturb event_statistics. We perturbed hsig, duration,
# tp1, and dir [tideResid was already unique].
#
# To do the check, compute the column sums of all event statistics. They should
# all be unique
#
max_es_vals = sapply(store_var_list, 
    f<-function(x) colSums(x$event_statistics[,1:4], na.rm=TRUE))
stopifnot(length(unique(max_es_vals)) == length(max_es_vals))


#
# Function to summarise parameters from perturbed runs
#
perturbed_summary<-function(variable_name){

    variable_vals = sapply(store_var_list,
        f<-function(x){
            num = x[[variable_name]] 
            denom = original_var_list[[variable_name]]
            # Handle errors gracefully
            if( (length(num) != length(denom)) || any(is.na(num)) || 
                (class(num) != class(denom))){
                num = denom*NA
            }
            return(num)
        }
    )

    variable_vals = t(variable_vals)
    print(summary(variable_vals))
    return(invisible())

}

#
# Function to summarise "relative errors" in perturbed runs
# i.e. (perturbed - original)/abs(original)
#
relative_error_summary<-function(variable_name){
    variable_differences = sapply(store_var_list,
        f<-function(x){
            num = x[[variable_name]] - original_var_list[[variable_name]]
            denom = original_var_list[[variable_name]]
            # Handle errors gracefully
            if( (length(num) != length(denom)) || any(is.na(num)) || 
                (class(num) != class(denom))){
                num = denom*NA
            }
            return(num/abs(denom))
        }
    )

    variable_differences = t(variable_differences)
    print(summary(variable_differences))
    return(invisible())
}
```

**Report on variations in the selected family for each pair in the C-Vine copula**

```r
# Print the original copula info
print(original_var_list$copula_model$copula_fit_mle)
```

```
## $value
## [1] 476.933
## 
## $convergence
## [1] 0
## 
## $message
## [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
## 
## $counts
## function gradient 
##       15       15 
## 
## $RVM
## C-vine copula with the following pair-copulas:
## Tree 1:
## 1,5  Independence 
## 1,4  Frank (par = 2.18, tau = 0.23) 
## 1,3  Gaussian (par = 0.54, tau = 0.36) 
## 1,2  Survival Gumbel (par = 2.51, tau = 0.6) 
## 
## Tree 2:
## 2,5;1  Frank (par = -1.02, tau = -0.11) 
## 2,4;1  Independence 
## 2,3;1  Gaussian (par = 0.23, tau = 0.15) 
## 
## Tree 3:
## 3,5;2,1  Frank (par = 1.45, tau = 0.16) 
## 3,4;2,1  Frank (par = -0.6, tau = -0.07) 
## 
## Tree 4:
## 4,5;3,2,1  Independence
```

```r
#
# Get matrix defining vine structure
#
original_fitmat = original_var_list$copula_model$copula_fit_mle$RVM$Matrix
fitmat = lapply(store_var_list, f<-function(x) x$copula_model$copula_fit_mle$RVM$Matrix)
#
# For the CVine copula, these should all be the same
#
fitmat_equal_to_original = sapply(fitmat, f<-function(x) all(x-original_fitmat == 0))
print(summary(fitmat_equal_to_original))
```

```
##    Mode    TRUE    NA's 
## logical     100       0
```

```r
stopifnot(all(fitmat_equal_to_original))

#
# Get matrix defining chosen copula families [with integer encoding, see ?BiCopName]
#
original_copula_family = original_var_list$copula_model$copula_fit_mle$RVM$family
copula_family = lapply(store_var_list, f<-function(x) x$copula_model$copula_fit_mle$RVM$family)
# Convert copula_family to array for easier summaries
copula_family_array = array(NA, dim=c(dim(copula_family[[1]]), length(copula_family)))
for(i in 1:length(copula_family)) copula_family_array[,,i] = copula_family[[i]]
#
# Print the copula families as integers. The mapping between these and the family names
# can be seen from ?BiCopName
#
nc = ncol(copula_family[[1]])
nr = nrow(copula_family[[1]])
var_names = names(original_var_list$es_cop_reorder)

# Loop over the family matrix, and report on the pairs
for(i in nr:2){

    print('')
    print('#')
    level_i = paste0(' Level ', (nr + 1 - i))
    var_i = var_names[nr+1 - i]
    print(paste0('# ' , level_i))
    print('#')

    for(j in (i-1):1){
        var_j = var_names[nr + 1 -j]
        cat('\n')
        print('-----')
        cat(paste0(var_i,  ' vs ', var_j, ' ; ', level_i, '\n'))
        copula_table = table(copula_family_array[i,j,])
        names_int = as.numeric(names(copula_table))
        names(copula_table) = BiCopName(names_int, short=FALSE)
        cat(' Perturbed models chose these copulas: \n')
        print(copula_table)
        cat(paste0('Raw data fit was ', BiCopName(original_copula_family[i,j], short=FALSE), '\n'))
    }

}
```

```
## [1] ""
## [1] "#"
## [1] "#  Level 1"
## [1] "#"
## 
## [1] "-----"
## hsig vs duration ;  Level 1
##  Perturbed models chose these copulas: 
##        Gaussian Survival Gumbel 
##              76              24 
## Raw data fit was Survival Gumbel
## 
## [1] "-----"
## hsig vs tideResid ;  Level 1
##  Perturbed models chose these copulas: 
## Gaussian 
##      100 
## Raw data fit was Gaussian
## 
## [1] "-----"
## hsig vs steepness ;  Level 1
##  Perturbed models chose these copulas: 
## Gaussian    Frank 
##       79       21 
## Raw data fit was Frank
## 
## [1] "-----"
## hsig vs dir ;  Level 1
##  Perturbed models chose these copulas: 
## Independence 
##          100 
## Raw data fit was Independence
## [1] ""
## [1] "#"
## [1] "#  Level 2"
## [1] "#"
## 
## [1] "-----"
## duration vs tideResid ;  Level 2
##  Perturbed models chose these copulas: 
##        Gaussian           Frank Survival Gumbel 
##              16              37              47 
## Raw data fit was Gaussian
## 
## [1] "-----"
## duration vs steepness ;  Level 2
##  Perturbed models chose these copulas: 
## Independence 
##          100 
## Raw data fit was Independence
## 
## [1] "-----"
## duration vs dir ;  Level 2
##  Perturbed models chose these copulas: 
## Frank 
##   100 
## Raw data fit was Frank
## [1] ""
## [1] "#"
## [1] "#  Level 3"
## [1] "#"
## 
## [1] "-----"
## tideResid vs steepness ;  Level 3
##  Perturbed models chose these copulas: 
## Frank 
##   100 
## Raw data fit was Frank
## 
## [1] "-----"
## tideResid vs dir ;  Level 3
##  Perturbed models chose these copulas: 
## Frank 
##   100 
## Raw data fit was Frank
## [1] ""
## [1] "#"
## [1] "#  Level 4"
## [1] "#"
## 
## [1] "-----"
## steepness vs dir ;  Level 4
##  Perturbed models chose these copulas: 
## Independence 
##          100 
## Raw data fit was Independence
```
