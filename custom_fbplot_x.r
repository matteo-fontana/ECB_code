fbplot.custom_x=function (Data, Depths = "MBD", Fvalue = 1.5, adjust = FALSE, 
                          display = TRUE, xlab = NULL, ylab = NULL, main = NULL, ylim=NULL, ygrid=NULL, x_range = NULL,...) 
{ 
  time_grid = seq(Data$t0, Data$tP, length.out = Data$P)
  if (is.character(Depths)) {
    Depths_spec = Depths
    if (Depths_spec == "MBD") {
      Depths = MBD(Data$values, manage_ties = TRUE)
    }
    else {
      Depths = eval(parse(text = paste(Depths, "( Data$values )", 
                                       sep = "")))
    }
  }
  else {
    stopifnot(length(Depths) == Data$N)
  }
  if (!is.list(adjust)) {
    out = .fbplot_fData(time_grid, Data$values, Depths, Fvalue)
  }
  else {
    nodenames = c("N_trials", "trial_size", "TPR", 
                  "F_min", "F_max", "tol", "maxiter", 
                  "VERBOSE")
    unused = setdiff(names(adjust), nodenames)
    if (length(unused) > 0) {
      for (i in unused) warning("Warning: unused parameter ", 
                                i, " in adjust argument of fbplot")
    }
    N_trials = ifelse(is.null(adjust$N_trials), 20, adjust$N_trials)
    trial_size = ifelse(is.null(adjust$trial_size), 8 * Data$N, 
                        adjust$trial_size)
    TPR = ifelse(is.null(adjust$TPR), 2 * pnorm(4 * qnorm(0.25)), 
                 adjust$TPR)
    F_min = ifelse(is.null(adjust$F_min), 0.5, adjust$F_min)
    F_max = ifelse(is.null(adjust$F_max), 5, adjust$F_max)
    tol = ifelse(is.null(adjust$tol), 0.001, adjust$tol)
    maxiter = ifelse(is.null(adjust$maxiter), 100, adjust$maxiter)
    VERBOSE = ifelse(is.null(adjust$VERBOSE), FALSE, adjust$VERBOSE)
    Cov = robustbase::covOGK(Data$values, sigmamu = robustbase::s_Qn)$cov
    CholCov <- chol(Cov)
    centerline = Data$values[which.max(Depths), ]
    Fvalues = rep(0, N_trials)
    cost_functional = function(F_curr) (length(.fbplot_fData(time_grid, 
                                                             Data_gauss, Depths = Depths_spec, Fvalue = F_curr)$ID_out)/trial_size - 
                                          TPR)
    for (iTrial in 1:N_trials) {
      if (VERBOSE > 0) {
        cat(" * * * Iteration ", iTrial, " / ", 
            N_trials, "\n")
      }
      Data_gauss = generate_gauss_fdata(trial_size, centerline, 
                                        CholCov = CholCov)
      if (VERBOSE > 0) {
        cat(" * * * * beginning optimisation\n")
      }
      opt = uniroot(cost_functional, interval = c(F_min, 
                                                  F_max), tol = tol, maxiter = maxiter)
      if (VERBOSE > 0) {
        cat(" * * * * optimisation finished.\n")
      }
      Fvalues[iTrial] = opt$root
    }
    Fvalue = mean(Fvalues)
    out = .fbplot_fData(time_grid, Data$values, Depths, Fvalue = 1.5)
  }
  
  ID_out = out$ID_out
  median_candidates=Data$values[Depths==max(Depths),]
  med=colMeans(matrix(median_candidates,ncol=dim(Data$values)[2]))
  check_down=t(t(Data$values[ID_out,])<med)
  check_down_vec=rowSums(check_down)>=ncol(Data$values)/2
  if(length(ID_out)<=1){check_down_vec=sum(check_down)>=ncol(Data$values)/2}
  check_down_vec=check_down_vec+1
  index_out_down=ID_out[check_down_vec==1]
  index_out_up=ID_out[check_down_vec==2]
  
  
  if (is.null(ylim)) ylim = range(Data$values)
  if (is.numeric(display)) {
    dev.set(display)
  }
  if (!display == FALSE) {
    #col_non_outlying = (scales::hue_pal(h = c(90, 180), l = 60))(Data$N - length(ID_out))
    #col_non_outlying = set_alpha(col_non_outlying, 0.5)
    col_non_outlying = 'grey90'
    
    if (length(ID_out) > 0) {
      
      col_outlying_ns_down = colorRampPalette(c('red',"#FFC6C6"))(length(index_out_down))
      col_outlying_ns_up = colorRampPalette(c('navyblue',"#C6C6FF"))(length(index_out_up))
      Depths_temp=Depths[ID_out]
      index_up_sorted=index_out_up[order(Depths[index_out_up])]
      index_down_sorted=index_out_down[order(Depths[index_out_down])]
      ID_out=rbind(index_down_sorted,index_up_sorted)
      col_outlying=rbind(col_outlying_ns_down,col_outlying_ns_up)
      
    }
    else {
      col_outlying = (scales::hue_pal(h = c(-90, 180), 
                                      c = 150))(1)
    }
    col_envelope = set_alpha("grey", alpha = 0.4)
    col_center = set_alpha("gray10", alpha = .8)
    col_fence_structure = set_alpha("gray10", alpha = .9)
    if(!is.null(ygrid)) time_grid=ygrid
   
     xlab = ifelse(is.null(xlab), "", xlab)
    ylab = ifelse(is.null(ylab), "", ylab)
    main = ifelse(is.null(main), "", main)
    if (length(ID_out) > 0) {
      matplot(time_grid, t(Data$values[-ID_out, ]), lty = 1, 
              type = "l", col = col_non_outlying, ylim = ylim, 
              xlab = xlab, ylab = ylab, main = main,xlim=x_range)
      max_envelope_limit = apply(Data$values[-ID_out, ], 
                                 2, max)
      min_envelope_limit = apply(Data$values[-ID_out, ], 
                                 2, min)
    }
    else {
      matplot(time_grid, t(Data$values), lty = 1, type = "l", 
              col = col_non_outlying, ylim = ylim, 
              xlab = xlab, ylab = ylab, main = main,xlim=x_range)
      max_envelope_limit = apply(Data$values, 2, max)
      min_envelope_limit = apply(Data$values, 2, min)
    }
    polygon(c(time_grid, rev(time_grid)), c(out$min_envelope_central, 
                                            rev(out$max_envelope_central)), col = col_envelope, 
            border = NA)
    lines(time_grid, out$max_envelope_central, lty = 1, col = col_center, 
          lwd = 3)
    lines(time_grid, out$min_envelope_central, lty = 1, col = col_center, 
          lwd = 3)
    lines(time_grid, Data$values[which.max(Depths), ], lty = 1, 
          type = "l", col = col_center, lwd = 3)
    lines(time_grid, max_envelope_limit, lty = 1, col = col_fence_structure, 
          lwd = 3)
    lines(time_grid, min_envelope_limit, lty = 1, col = col_fence_structure, 
          lwd = 3)
    half.time_grid = which.min(time_grid - 0.5)
    lines(c(time_grid[half.time_grid], time_grid[half.time_grid]), 
          c(out$max_envelope_central[half.time_grid], max_envelope_limit[half.time_grid]), 
          lty = 1, col = col_fence_structure, lwd = 3)
    lines(c(time_grid[half.time_grid], time_grid[half.time_grid]), 
          c(out$min_envelope_central[half.time_grid], min_envelope_limit[half.time_grid]), 
          lty = 1, col = col_fence_structure, lwd = 3)
    if (length(ID_out) > 0) {
      matplot(time_grid, t(toRowMatrixForm(Data$values[ID_out, 
      ])), lty = 1, type = "l", col = col_outlying, 
      lwd = 3, add = T)
    }
  }
  return(list(Depth = Depths, Fvalue = Fvalue, ID_outliers = ID_out,Col_outlying=col_outlying))
}

