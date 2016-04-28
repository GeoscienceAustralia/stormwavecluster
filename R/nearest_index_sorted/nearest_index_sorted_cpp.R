# Rcpp code to find the index of a sorted vector x_c  (with x_c[i] <= x_c[i+1])
# nearest to a real number season_c.
#
# Uses bisection.

library(Rcpp)

cpp_nearest_index_sorted<-'
IntegerVector nearest_index_sorted_cpp(NumericVector x_c, NumericVector season_c, int check_is_sorted){
    //Rcpp::NumericVector x_c(x);        
    //Rcpp::NumericVector season_c(season);

    int n = x_c.size();
    int n2 = season_c.size();
    int i, j, upper, lower;
    IntegerVector out(n2);

    if(check_is_sorted == 1){
        for(j = 0; j < (n - 1); j++ ){
            if( x_c[j] > x_c[j+1] ){
                stop("First argument must be a monotonic non-decreasing vector"); 
            }
        }
    }


    for(j = 0; j < n2; j++){

        // Use bisection to find the nearest index

        if( season_c[j] < x_c[0]){ 
            // Quick exit
            out[j] = 1;
        }else{
            if( season_c[j] > x_c[n-1]){
                // Quick exit
                out[j] = n ;
            }else{
                // Bisect
                lower = 1;
                upper = n;
                i = floor(0.5*(lower + upper));

                while((x_c[i-1] > season_c[j])||(x_c[i] < season_c[j])){
                   if(x_c[i-1] > season_c[j]){
                       upper = i;
                   }else{
                       lower = i;
                   }
                   i = floor(0.5*(lower + upper));
                }

                if ( season_c[j] - x_c[i-1] > x_c[i] - season_c[j]){
                    out[j] = i+1;
                }else{
                    out[j] = i;
                }

            }
        }
    }
    return(out);
}

    '

# Convert it to an R function
cppFunction(cpp_nearest_index_sorted)

