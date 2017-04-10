#' Clustered event counting function -- consider recoding in Rcpp for speed
#'
#' This function should count the number of separate times that there are
#' 'nevents' events with hsig > hsig_threshold within a particular time-window.
#'
#' @param synthetic_attr data.frame or list with the storm statistics,
#' including duration (hours), startyear (years), hsig (m)
#' @param hsig_threshold threshold for all events in cluster
#' @param nevent Number of events in cluster
#' @param time_window numeric. maximum time (years) between the start-times of
#' the first/last events in the cluster
#' @param return_detailed logical. If TRUE, give the cluster start/end times
#' @param max_hsig_threshold numeric. We can require that the largest event in the cluster exceeds max_hsig_threshold
event_cluster_counter<-function(
    synthetic_attr, hsig_threshold = 5, 
    nevents = 3, time_window=6*7/365.25,
    return_detailed=FALSE, 
    max_hsig_threshold=hsig_threshold){

    stopifnot(max_hsig_threshold >= hsig_threshold)

    event_starts = synthetic_attr$startyear
    # Assume 'duration' is in hours
    event_duration = synthetic_attr$duration/(24*365.25)

    # Check for bad inputs
    stopifnot(min(diff(event_starts)) > 0)
    stopifnot(min(event_duration) > 0)

    event_ends = event_starts + event_duration
    event_hsig = synthetic_attr$hsig

    # First filter on events
    keepers = which(event_hsig >= hsig_threshold)

    if(length(keepers) < nevents){
        # We do not have enough events
        cluster_count = 0
    }else{

        event_starts = event_starts[keepers]
        event_ends = event_ends[keepers]
        event_hsig = event_hsig[keepers]
        lk = length(event_starts)

        i = 1
        cluster_starts = event_starts*NA
        cluster_ends = event_ends*NA
        cluster_count = 0
        les = length(event_starts)
        while(i < les){
            # Find potential cluster
            cluster_inds = which(event_starts[i:les] <= event_starts[i] + time_window) + i-1
            lci = length(cluster_inds)

            if( (lci < nevents) || 
                (max(event_hsig[cluster_inds]) < max_hsig_threshold) ){
                # This point is not the start of a cluster, move on
                i = i+1
            }else{
                # We have a new cluster
                cluster_count = cluster_count+1
                cluster_starts[cluster_count] = event_starts[cluster_inds[1]]
                cluster_ends[cluster_count] = event_ends[cluster_inds[lci]]
                i = i + lci
            }
        }

        ### Get the max hsig for each potential 'cluster' (taken as a sequence of
        ### 'nevents' events above hsig_threshold)
        ##cluster_max_hsig = event_starts[1:(lk-nevents+1)]*0 - 999999
        ##for(i in 1:nevents){
        ##    cluster_max_hsig = pmax(cluster_max_hsig, event_hsig[(1:(lk-nevents+1)) + (i-1)])
        ##}

        ### Find 
        ##cluster_start_index = which(event_starts[nevents:lk] - event_starts[1:(lk-nevents+1)] <= time_window)


        ##if(length(cluster_start_index) == 0 || 
        ##   max(cluster_max_hsig[cluster_start_index]) < max_hsig_threshold){
        ##    # We do not have enough events sufficiently closely spaced, or they do
        ##    # not achieve a high enough 'max hsig' within the cluster
        ##    cluster_count = 0
        ##}else{
        ##    # Typical case where we have events matching the definition

        ##    cluster_starts = event_starts[cluster_start_index]
        ##    cluster_ends = event_ends[cluster_start_index+(nevents-1)]
        ##    cluster_max_hsig = cluster_max_hsig[cluster_start_index]


        ##    # Ensure the first event has max_hsig > max_hsig_threshold
        ##    n = which(cluster_max_hsig > max_hsig_threshold)[1]
        ##    lc = length(cluster_starts)
        ##    cluster_starts = cluster_starts[n:lc]
        ##    cluster_ends = cluster_ends[n:lc]
        ##    cluster_max_hsig = cluster_max_hsig[n:lc]

        ##    # Remove non-independent events with a simple approach
        ##    # If there are many, this is slow
        ##    i = 1
        ##    while(i < length(cluster_starts)){
        ##        i = i+1
        ##        if(cluster_starts[i] < cluster_ends[i-1] | 
        ##           cluster_max_hsig[i] < max_hsig_threshold){
        ##            cluster_starts = cluster_starts[-i]
        ##            cluster_ends = cluster_ends[-i]
        ##            cluster_max_hsig = cluster_max_hsig[-i]
        ##            i = i-1
        ##        }
        ##    }
        ##
        ##    cluster_count = length(cluster_starts)
        ##}
    }
   
    # Might be useful to get the detailed information sometimes 
    if(return_detailed){
        if(cluster_count == 0){ 
            return(cbind(0,0))
        }else{
            return(cbind(cluster_starts, cluster_ends))
        }
    }

    return(cluster_count)
}



test_event_cluster_counter<-function(){

    # Make a synthetic data series for testing
    synthetic_attr = list(
        ########            --Cluster1---                --Cluster2--       ---Cluster3---
        hsig =      c(1,    3,   4,   3,    4, 0.5, 0.9,    3,    3,  0.8,   4,    5,    6,  0.2, 0.2),
        startyear = c(0, 0.38, 0.4, 0.5, 0.52, 0.7, 1.0, 1.01, 1.02, 1.22, 1.3, 1.31, 1.32, 1.33, 1.34),
        # Note duration is in hours, so is smaller than the time between events
        duration = rep(24, length=15))

    # Find the 3 clusters
    test1 = event_cluster_counter(synthetic_attr, hsig_threshold = 2.9, nevents=2, time_window=0.2)
    stopifnot(test1 == 3)
    print('PASS')

    # Only cluster1 matches
    test2 = event_cluster_counter(synthetic_attr, hsig_threshold = 2.9, nevents=4, time_window=0.2)
    stopifnot(test2 == 1)
    print('PASS')

    # Cluster1 is split in 2
    test3 = event_cluster_counter(synthetic_attr, hsig_threshold = 2.9, nevents=2, time_window=0.03)
    stopifnot(test3 == 4)
    print('PASS')

    # Hsig not large enough
    test4 = event_cluster_counter(synthetic_attr, hsig_threshold = 7, nevents=2, time_window=0.03)
    stopifnot(test4 == 0)
    print('PASS')

    # Max hsig not large enough
    test5 = event_cluster_counter(synthetic_attr, hsig_threshold = 2.9, nevents=2, time_window=0.03, max_hsig_threshold=7)
    stopifnot(test5 == 0)
    print('PASS')
   
    #Cluster 1 is split in 2, cluster2 is ignored 
    test6 = event_cluster_counter(synthetic_attr, hsig_threshold = 2.9, nevents=2, time_window=0.03, max_hsig_threshold=3.5)
    stopifnot(test6 == 3)
    print('PASS')

    # Cluster2 is ignored
    test7 = event_cluster_counter(synthetic_attr, hsig_threshold = 2.9, nevents=2, time_window=0.2, max_hsig_threshold=3.5)
    stopifnot(test7 == 2)
    print('PASS')

    # All are ignored
    test8 = event_cluster_counter(synthetic_attr, hsig_threshold = 2.9, nevents=4, time_window=0.2, max_hsig_threshold=5)
    stopifnot(test8 == 0)
    print('PASS')
    
}
