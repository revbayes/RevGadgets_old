################################################################################
#
# @brief Plotting the output of a population size inference in an Occurrence Birth Death Process analysis.
#
# @date Last modified: 2020-04-03
# @author Jeremy Andreoletti
#
# @param    Kt_mean                   data.frame        The processed output for plotting.
# @param    xlab                      character         The label of the x-axis.
# @param    ylab                      character         The label of the y-axis.
# @param    line.size                 numeric           The width of the lineage plot line.
# @param    interval.line.size        numeric           The width of the credence interval.
# @param    col.Hidden                character         The color of the hidden lineages plot line.
# @param    col.Observed              character         The color of the observed lineages plot line.
# @param    col.Total                 character         The color of the total lineages plot line.
# @param    col.Hidden.interval       character         The color of the credence interval lines around the hidden lineages plot.
# @param    col.Total.interval        character         The color of the credence interval lines around the total lineages plot.
# @param    palette.Hidden            character         The palette of the hidden lineages plot distribution.
# @param    palette.Total             character         The palette of the total lineages plot distribution.
# @param    show.Hidden               boolean           Wether to show the plot for hidden lineages.
# @param    show.Observed             boolean           Wether to show the plot for observed lineages.
# @param    show.Total                boolean           Wether to show the plot for total lineages.
# @param    show.intervals            boolean           Wether to show the credence intervals.
# @param    show.densities            boolean           Wether to show the diversity densities.
# @param    show.expectations         boolean           Wether to show the diversity expectations.
# @param    use.interpolate           boolean           Wether to interpolate densities.
#
#
################################################################################

rev.plot.nbLineages = function( Kt_mean,
                                xlab="Time",
                                ylab="Number of lineages",
                                col.Hidden = "dodgerblue3",
                                col.Observed = "gray25",
                                col.Total = "forestgreen",
                                col.Hidden.interval = "dodgerblue2",
                                col.Total.interval = "darkolivegreen4",
                                palette.Hidden = c("transparent", "dodgerblue2", "dodgerblue3", "dodgerblue4", "black"),
                                palette.Total = c("transparent", "green4", "forestgreen", "black"),
                                line.size=0.7,
                                interval.line.size=0.5,
                                show.Hidden=TRUE,
                                show.Observed=TRUE,
                                show.Total=TRUE,
                                show.intervals=TRUE,
                                show.densities=TRUE,
                                show.expectations=TRUE,
                                use.interpolate=TRUE ){
  
  N <- length(Kt_mean)-8                                                                # Maximal number of hidden lineages
  
  ## Format Kt_mean for plotting
  Kt_mean_plot <- Kt_mean %>% pivot_longer(-c("TimePoints", "NbObservedLin", "aggregNbHiddenLin", "aggregNbTotalLin", "NbHiddenLin0.025", "NbHiddenLin0.5", "NbHiddenLin0.975", "NbTotalLin0.025", "NbTotalLin0.5", "NbTotalLin0.975"), names_to="NbHiddenLin", values_to="ProbabilityDensity")
  Kt_mean_plot$NbHiddenLin <- as.integer(Kt_mean_plot$NbHiddenLin)
  # Kt_mean_plot <- Kt_mean_plot[Kt_mean_plot$ProbabilityDensity!=0,]                   # Remove 0-probability rows
  
  ## Get the distribution of the total number of lineages
  Kt_mean_plot$NbTotalLin <- Kt_mean_plot$NbObservedLin + Kt_mean_plot$NbHiddenLin
  Kt_mean_plot <- Kt_mean_plot[Kt_mean_plot$NbTotalLin < max(Kt_mean$NbTotalLin0.975)*1.1,]   # Remove lowest probability rows
  if (!show.Total){
    Kt_mean_plot <- Kt_mean_plot[Kt_mean_plot$NbHiddenLin < max(Kt_mean$NbHiddenLin0.975)*1.1,]   # Remove lowest probability rows
  }
  
  ## Plot densities and LTTs
  cols    <- c( "c1" = col.Hidden, "c2" = col.Observed, "c3" = col.Total, "H.95%" = col.Hidden.interval, "T.95%" = col.Total.interval )
  
  p <- ggplot(Kt_mean_plot, aes(x=TimePoints, y=NbTotalLin, z = ProbabilityDensity))
  
  ### Plot the number of hidden lineages : full distribution + aggregated expectation + credence interval
  if (show.Hidden){
    if (show.densities){
      p <- p + annotate(geom="raster", x=Kt_mean_plot$TimePoints, y=Kt_mean_plot$NbHiddenLin, interpolate = use.interpolate,
                        fill = colour_ramp(c(colorRampPalette(palette.Hidden)(N)))(Kt_mean_plot$ProbabilityDensity))
    }
    if (show.expectations){
      p <- p + geom_line(aes(y=aggregNbHiddenLin, color="c1"), size=line.size)
    }
    if (show.intervals){
      p <- p + geom_line(aes(y=NbHiddenLin, color="H.95%"), alpha=0) +  # Fake plot for the legend
        annotate(geom="line", x=Kt_mean$TimePoints, y=Kt_mean$NbHiddenLin0.025, color=cols["H.95%"], linetype="twodash", size=interval.line.size) +
        annotate(geom="line", x=Kt_mean$TimePoints, y=Kt_mean$NbHiddenLin0.975, color=cols["H.95%"], linetype="twodash", size=interval.line.size)
    }
  }
  
  ### Plot the total number lineages : full distribution + aggregated expectation + credence interval
  if (show.Total){
    if (show.densities){
      p <- p + annotate(geom="raster", x=Kt_mean_plot$TimePoints, y=Kt_mean_plot$NbTotalLin, interpolate = use.interpolate,
                        fill = colour_ramp(colorRampPalette(palette.Total)(N))(Kt_mean_plot$ProbabilityDensity))
      }
    if (show.expectations){
      p <- p + geom_line(aes(y=aggregNbTotalLin, color="c3"), size=line.size)
    }
    if (show.intervals){
      p <- p + geom_line(aes(y=NbTotalLin0.025, color="T.95%"), alpha=0) +  # Fake plot for the legend
        annotate(geom="line", x=Kt_mean$TimePoints, y=Kt_mean$NbTotalLin0.025, color=cols["T.95%"], linetype="twodash", size=interval.line.size) +
        annotate(geom="line", x=Kt_mean$TimePoints, y=Kt_mean$NbTotalLin0.975, color=cols["T.95%"], linetype="twodash", size=interval.line.size)
      }
  }
  
  ### Plot the number of observed lineages : mean
  if (show.Observed){
    p <- p +
      geom_line(aes(y=NbObservedLin, color="c2"), size=line.size)
  }
  
  ### Legend and axes
  p <- p + 
    scale_color_manual(name = "Lineages", breaks = c("c1", "H.95%", "c2", "c3", "T.95%"), 
                       values = cols, labels = c("Hidden", "95% credence interval", "Observed", "Total", "95% credence interval")) +
    scale_x_continuous(name = xlab, expand = c(0.01,0.01)) +   
    scale_y_continuous(name = ylab, expand = c(0.01,0.01)) + 
    theme(panel.background=element_rect(fill="white"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ggtitle("Probability density of the number of lineages through time")
  
  return (p)
}



################################################################################
#
# @brief Processing the output of the MCMC analysis + population size (diversity) inference in the Occurrence Birth Death Process.
#
# @date Last modified: 2020-04-03
# @author Jeremy Andreoletti
#
# @param    popSize_distribution_matrices_file      character      The number of expected diversification-rate changes.
# @param    trees_trace_file                        character      The number of expected diversification-rate changes.
# @param    weight_trees_posterior                  bool           Wether to combine trees uniformly or weighted according to their posterior probabilities.
#
#
################################################################################

rev.process.nbLineages = function( popSize_distribution_matrices_file,
                                   trees_trace_file, 
                                   weight_trees_posterior=T ){

  ## Import Kt : probability distribution of the number of hidden lineages through time
  Kt_trace_lines <- readLines(popSize_distribution_matrices_file)
  S <- diff(which(lapply(Kt_trace_lines, function(x){grep(";", x)})==1))[1]         # Get the number of time points (nb of lines between two ";" separators)
  Kt_trace_lines <- gsub("\\[| |\\]| ,|\t|;", "", Kt_trace_lines)[-1]               # Remove unwanted characters
  Kt_trace <- read.csv(text = Kt_trace_lines, header = FALSE, na.strings = "nan")
  Kt_trace[is.na(Kt_trace)] <- 0                                                    # Set NA values to 0
  N <- length(Kt_trace)-1                                                           # Maximal number of hidden lineages
  names(Kt_trace) <- 0:N                                                            # Set names to the number of hidden lineages

  ## Import the corresponding tree : get the number of observed lineages through time (LTT)
  trees <- read.table(trees_trace_file, header = T)
  trees$obd_tree <- sapply(trees$obd_tree, function(tree){read.tree(text=as.character(tree))})
  nb_trees <- length(trees$Iteration)                                               # Total number of trees
  burnin <- nb_trees-length(Kt_trace[,1])/S                                         # Number of trees in the burnin

  ## Add the iterations to Kt_trace
  iterations <- trees$Iteration[(burnin+1):nb_trees]
  Kt_trace$Iteration <- rep(iterations, each=S)                                     # Add an iteration number column

  ## Browse all iterations
  Kt_mean <- matrix(0, nrow=S, ncol=N+1)
  observedLin_mean <- rep(0, S)
  start_age <- function(tree){ltt.plot.coords(tree)[1]}
  timePoints <- seq(0, min(sapply(trees$obd_tree, start_age)), length.out = S)    # Get S time points between 0 and the oldest starting time
  posteriors <- exp(trees$Posterior-max(trees$Posterior))
  posteriors <- posteriors/sum(posteriors)
  for (i in (burnin+1):nb_trees){
    ### Increment the distribution of number of hidden lineades
    it <- trees$Iteration[i]
    print(paste("Tree", it, "out of", trees$Iteration[nb_trees]))
    Kt <- Kt_trace[Kt_trace$Iteration==it,-which(names(Kt_trace)=="Iteration")]
    lines_sum <- apply(Kt, 1, sum)
    for (j in 1:S){
      if (lines_sum[j]==0){
        Kt[j,1] <- 1.0                                      # Empty lines are considered to have 0 hidden lineages
      }
      else if (lines_sum[j]!=0){
        Kt[j,] <- Kt[j,]/lines_sum[j]                       # Normalise lines to 1 (several lines at 0.9999 or 1.0001)
      }
    }
    if (weight_trees_posterior){
      Kt_mean <- Kt_mean + Kt*posteriors[i]
    }
    else{ Kt_mean <- Kt_mean + Kt/length(iterations) }

    ### Get the LTT coordinates
    obd_tree <- collapse.singles(trees$obd_tree[[i]])       # Remove single nodes (ie. sampled ancestors)
    LTT <- data.frame(ltt.plot.coords(obd_tree))            # Extract LTT coordinates
    LTT$time <- round(LTT$time, 4)                          # Reduce precision (extant tips wrongly at time -0.000001)

    ### Increment number of observed lineages
    getNbObservedLin <- function(t){
      if (t < LTT$time[1]) return (0)
      first_consecutive_time_point <- which(LTT$time >= t)[1]
      return (LTT$N[first_consecutive_time_point])
    }
    observedLin <- sapply(timePoints, getNbObservedLin)
    if (weight_trees_posterior){
      observedLin_mean <- observedLin_mean + observedLin*posteriors[i]
    }
    else{ observedLin_mean <- observedLin_mean + observedLin/length(iterations) }
  }

  ## Get the most probable number of hidden lineages (weighted mean according to their respective probabilities)
  hiddenLin <- data.frame(weightedMean=matrix(apply(Kt_mean, 1, function(x){weighted.mean(as.integer(names(Kt_mean)), x)})),
                          maxProb=matrix(apply(Kt_mean, 1, function(x){as.integer(names(Kt_mean)[which(x==max(x))])})))
  Kt_mean$aggregNbHiddenLin <- hiddenLin$weightedMean
  # Kt_mean$aggregNbHiddenLin[is.na(Kt_mean$aggregNbHiddenLin)] <- 0           # Put NA values at 0 (`weighted.mean` artifact at t0)

  ## Get the aggregated number of observed and total lineages
  Kt_mean$aggregNbTotalLin <- observedLin_mean + Kt_mean$aggregNbHiddenLin
  Kt_mean$NbObservedLin <- round(observedLin_mean)
  
  ## Get the 95% credence interval around the number of hidden lineages
  Kt_mean[1:(N+1)] <- t(apply(Kt_mean[1:(N+1)], 1, function(row){row/sum(row)}))   # Force lines to sum to 1 (correct numerical uncertainties)
  Kt_mean_cumsum <- apply(Kt_mean[1:(N+1)], 1, cumsum)
  Kt_mean$NbHiddenLin0.025 <- apply(Kt_mean_cumsum, 2, function(col){which(col>0.025)[1]-1})
  Kt_mean$NbHiddenLin0.5 <- apply(Kt_mean_cumsum, 2, function(col){which(col>0.5)[1]})
  Kt_mean$NbHiddenLin0.975 <- apply(Kt_mean_cumsum, 2, function(col){which(col>0.975)[1]})

  ## Get the 95% credence interval around the total number of lineages
  Kt_mean$NbTotalLin0.025 <- observedLin_mean + Kt_mean$NbHiddenLin0.025
  Kt_mean$NbTotalLin0.5 <- observedLin_mean + Kt_mean$NbHiddenLin0.5
  Kt_mean$NbTotalLin0.975 <- observedLin_mean + Kt_mean$NbHiddenLin0.975
  
  ## Get the aggregated number of observed and total lineages
  Kt_mean$TimePoints <- timePoints

  return(Kt_mean)
}

