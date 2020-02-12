#' Interpolation and extrapolation of phylogenetic diversity
#'
#' \code{inextPD_yhc}: Interpolation and extrapolation of phylogenetic diversities based on a framework of Hill numbers with order q.
#' @param data a matrix/data.frame of species abundances/incidences data.\cr Type (1) abundance data: When there are N assemblages, the
#' observed species abundances should be arranged as a species (in rows) by assemblage (in columns) matrix. The first row
#' (including N entries) lists the assemblage labels or site names for the N assemblages.\cr Type (2) incidence data:
#' The data input format for incidence data must be raw detection/non-detection data. That is, data for each community/assemblage
#' consist of a species-by-sampling-unit matrix. Users must first merge multiple-community data by species identity to obtain a pooled
#' list of species; then the rows of the input data refer to this pooled list. \cr
#' @param tree  a phylogenetic tree of class phylo among the observed species in the pooled assemblage.\cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}), default is "abundance". \cr
#' @param t_ a nammed vector of nonnegative integers specifying the sampling units in each community. Ignored if \code{datatype = "abundance"}.\cr
#' @param q a nonnegative integer value or integer vector specifying the diversity order of Hill numbers, default is 0. \cr
#' @param endpoint an interger specifying the endpoint for rarefaction and extrapolation range. If \code{NULL}, \code{endpoint} = double of the maximum
#' reference sample size. It will be ignored if \code{size} is given. \cr
#' @param knots an integer specifying the number of knot between 1 and the \code{endpoint}, default is 40.\cr
#' @param size a vector of nonnegative integers specifying the sample sizes for which diversity estimates will be calculated. If \code{NULL}, the diversity estimates will
#' be calculated for those sample sizes determined by the specified/default \code{endpoint} and \code{knot}. \cr
#' @param plot.type an integer or an integer vector specifying types of curves. Three types of plots: sample-size-based rarefaction and extrapolation curve (\code{plot.type = 1}); coverage-based rarefaction
#' and extrapolation curve (\code{plot.type = 2}); sample completeness curve (\code{plot.type = 3}), default is \code{plot.type = 1}. \cr
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
#' @param nboot an integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap; in this case,
#' the computation of s.e. and confidence intervals will be skipped, default is 50.
#' @param reftime a positive number specifying the reference time.
#' @import chaoUtility
#' @import ggplot2
#' @import dplyr
#' @importFrom stats rbinom
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom ape drop.tip
#' @importFrom phyclust get.rooted.tree.height
#' @return a list of four objects: \cr\cr
#' \code{$summary} summary of data. \cr\cr
#' \code{$refernce.point} root of phylogeny tree for the species in the pooled assemblage. \cr\cr
#' \code{$inext} table of various estimates for interpolated or extrapolated samples and their confidence intervals. \cr\cr
#' \code{$figure} showing chosen \code{plot.type} of sampling curves. \cr\cr
#' @examples
#' \donttest{
#' # Type (1) abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- inextPD_yhc(data = data, tree = tree,datatype = "abundance",q = c(0,1,2), plot.type = 1:3)
#' # Type (2) incidence data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' t_ <- data.inc$t
#' out <- inextPD_yhc(data = data, tree = tree,datatype = "incidence_raw",t_ = t_,q = c(0,1,2), plot.type = 1:3)
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609.\cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015) Rarefaction and extrapolation of phylogenetic diversity. Methods in Ecology and Evolution, 6, 380-388.\cr\cr
#' @export
inextPD_yhc <- function(data, tree, datatype = "abundance", t_, q = c(0,1,2), endpoint = NULL, knots = 40, size = NULL, plot.type = 1:3, conf = 0.95, nboot = 50,reftime=NULL) {
  DATATYPE <- c("abundance", "incidence_raw")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value", call. = FALSE)
  # qq <- 0:2
  # if(is.na(pmatch(q, qq) == T) == T) stop("invalid order of q, we only compute q = 0, 1 or 2", call. = FALSE)

  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf (confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)

  dat = list()
  name <- rownames(data)
  if(datatype=="incidence_raw"){
    if(ncol(data) != sum(t_)) stop("Number of columns does not euqal to the sum of key in sampling units", call. = FALSE)
    n <- 0
    for(i in 1:length(t_)){
      dat[[i]] <- data[,(n+1):(n+t_[i])]
      n <- n+t_[i]
    }
    if(is.null(names(t_))) {
      names(dat) <- paste0("site",1:length(t_))
    }else{
      names(dat) = names(t_)
    }
  }else{
    if(is.null(colnames(data))) {colnames(data) <- paste0("site",1:ncol(data))}
    dat <- lapply(1:ncol(data), function(i)  {x <- data[,i];names(x) <- name;x})
    names(dat) = colnames(data)
  }
  ###
  if (datatype=="incidence_raw") {
    pool.name = unique(unlist(sapply(dat, function(x)rownames(x)[rowSums(x)>0])))
    mydata = lapply(dat, function(X) X[pool.name, ])
    tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
    mytree <- drop.tip(tree,tip)
    H_max <- get.rooted.tree.height(mytree)$treeH
  }
  if (datatype=="abundance") {
    pool.name = unique(unlist(lapply(dat, function(x) names(x)[x>0] )))
    mydata = lapply(dat, function(X) X[pool.name])
    tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
    mytree <- drop.tip(tree,tip)
    H_max <- get.rooted.tree.height(mytree)
  }
  if(is.null(reftime)) {reft <- H_max
  }else if(reftime<=0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }else {reft <- reftime}

  ###########output1
  #atime <- Sys.time()
  if(class(mydata) == "list"){
    infos <- sapply(mydata, function(x){
      datainf_yhc(data = x, datatype, phylotr = mytree,reft = reft)})
  }else{
    return(NULL)
  }
  #btime <- Sys.time()
  #print(paste0('Info time:',btime-atime))
  ############output2
  if(datatype == "abundance") {
    up <- lapply(mydata, function(j) 2*sum(j))
  }else {
    up <- lapply(mydata, function(j) 2*ncol(j))
  }
  end <- max(unlist(up))

  FUN <- function(e){
    if(class(mydata)=="list"){
      Q = as.numeric(q)
      if(is.null(endpoint) == T) endpoint = end
      temp = iNEXTPD_yhc(datalist=mydata,phylotr=mytree,datatype,Q,nboot,conf,size,knots,endpoint,reft)
      temp$qPD.LCL[temp$qPD.LCL<0] <- 0;temp$SC.LCL[temp$SC.LCL<0] <- 0
      temp$SC.UCL[temp$SC.UCL>1] <- 1
      return(temp)
    }else{
      return(NULL)
    }
  }
  #atime <- Sys.time()
  # RE.table <- FUN(3)
  RE.table <- tryCatch(FUN(e), error = function(e){return()})
  #btime <- Sys.time()
  #print(paste0('R/E time:',btime-atime))
  ###############output3
  TYPE <- 1:3
  FUN2 <- function(e){
    if(is.na(sum(pmatch(plot.type, TYPE))) == F){
      temp2 <- lapply(plot.type, function(j) RE_plot_yhc(RE.table, datatype, j))
      allname <- c("RE.plot.size", "RE.plot.C", "RE.plot.sizeC")
      names(temp2) <- allname[plot.type]
      temp2
    }else{
      return("invalid plot type", call. = FALSE)
    }
  }
  #atime <- Sys.time()
  RE.plot <- tryCatch(FUN2(e), error = function(e){return()})
  #btime <- Sys.time()
  #print(paste0('plot time:',btime-atime))
  ans <- list(summary = infos,refernce_time = reft, inext = RE.table, figure = RE.plot)
  #class(ans) <- c("inextPD")
  ans

}

#' Asymptotic phylogenetic diversity profile for each community
#'
#' \code{PhdAsy}: An estimated asymptotic phylogenetic diversity profile for each community based on a framework of Hill numbers with order q and statistical method proposed in Chao et al. (2015) and Hsieh and Chao (2017).\cr
#' @param data a matrix/data.frame of species abundances/incidences data.\cr
#' See \code{\link{inextPD_yhc}} for data details.
#' @param tree a phylogenetic tree of class phylo among the observed species in the pooled assemblage. \cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}), default is "abundance". \cr
#' @param t_ a vector of nonnegative integers specifying the sampling units in each community.
#' Ignored if \code{datatype = "abundance"}.
#' @param q a nonnegative vector specifying the diversity order of Hill numbers. For phylogenetic diversity profiles,
#' \code{q} need to be a vector of length more than one, default is seq(0, 2, by = 0.25) \cr
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
#' @param nboot an integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap; in this case,
#' the computation of s.e. and confidence intervals will be skipped, default is 50.
#' @param reftime a positive number specifying the reference time.
#' @import chaoUtility
#' @import ggplot2
#' @import dplyr
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom ape drop.tip
#' @importFrom phyclust get.rooted.tree.height
#' @return  a list of four objects: \cr\cr
#' \code{$summary} summary of data. \cr\cr
#' \code{$refernce.point} root of the species in the pooled assemblage. \cr\cr
#' \code{$asy} a table showing all numerical values of diversity estimates and their s.e and confidence intervals. \cr\cr
#' \code{$figure} phylogenetic diversity profiles. \cr\cr
#' @examples
#' \donttest{
#' # Type (1) abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- PhdAsy_yhc(data = data, tree = tree, q = seq(0, 2, by = 0.25))
#' # Type (2) incidence data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' t_ <- data.inc$t
#' out <- PhdAsy_yhc(data = data, tree = tree,datatype = "incidence_raw", t_ = t_, q = seq(0, 2, by = 0.25))
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609. \cr\cr
#' @export
PhdAsy_yhc <- function(data, tree, datatype = "abundance", t_, q = seq(0, 2, by = 0.25), conf = 0.95, nboot = 50, reftime = NULL){
  dat = list()
  if(is.null(rownames(data))) stop("The rownames (species names) of data can't be empty. Species names in data must match those in the phylogenetic tree")
  pool.name <- rownames(data[rowSums(data)>0,])
  if (length(q) == 1) stop("length of q should be greater than one", call. = FALSE)
  if (sum(q<0)>=1) stop("q must be a positive number", call. = FALSE)
  if ((datatype != "incidence_raw") & (datatype != "abundance")) stop("invalid datatype", call. = FALSE)
  if ((conf < 0) | (conf < 0) | (is.numeric(conf)==F)) stop('conf"(confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(datatype=="incidence_raw"){
    if(ncol(data) != sum(t_)) stop("Number of columns does not euqal to the sum of key(t_) in sampling units", call. = FALSE)
    ntmp <- 0
    for(i in 1:length(t_)){
      dat[[i]] <- data[,(ntmp+1):(ntmp+t_[i])]
      ntmp <- ntmp+t_[i]
    }
    if(is.null(names(t_))) {
      names(dat) <- paste0("site",1:length(t_))
    }else{
      names(dat) = names(t_)
    }
  }else if (datatype == "abundance"){
    if(is.null(colnames(data))) {colnames(data) <- paste0("site",1:ncol(data))}
    dat <- lapply(1:ncol(data), function(i)  {x <- data[pool.name,i];names(x) <- pool.name;x})
    names(dat) = colnames(data)
  }
  ###
  if (datatype=="incidence_raw") {
    mydata = lapply(dat, function(X) X[pool.name, ])
  }
  if (datatype=="abundance") {
    mydata = lapply(dat, function(X) X[pool.name])
  }
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)$treeH
  if(is.null(reftime)) {reft <- H_max
  }else if(reftime<=0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }else {reft <- reftime}

  FUN = function(a){
    temp = AsyPD_yhc(datalist = mydata, datatype, phylotr = mytree, q, nboot, conf, reft)# mytree is pooled tree of class phylo
    ans <- list(summary = temp[[2]], refernce.point = reft, asy = temp[[1]], figure = Asy_plot_yhc(temp[[1]], 1))
    #class(ans) <- c("PhdAsy")
    return(ans)
  }
  #out <- FUN(3)
  out <- tryCatch(FUN(e), error = function(e){return()})
  return(out)
}

#' Observed q profile and time profile of phylogenetic diversity or phylogenetic Hill numbers for each community
#'
#' \code{PhdObs}: Observed phylogenetic diversity or phylogenetic Hill numbers profile for each community based on a framework of Hill numbers with order q and reference time.\cr
#' @param data a matrix/data.frame of species abundances/incidences data.\cr
#' See \code{\link{inextPD_yhc}} for data details.
#' @param tree a phylogenetic tree in the Newick format among the observed species in the pooled assemblage. \cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}), default is "abundance". \cr
#' @param t_ a named vector of nonnegative integers specifying the sampling units in each community. If the names are absent, each community will be
#' automatically assigned siteX as community names.
#' Ignored if \code{datatype = "abundance"}.
#' @param type two types for calculating phylogenetic profiles in Chiu, C.-H., Jost, L. and Chao(2014): Phylogenetic Diversity(\code{type = "PD"}) and
#' Phylogenetic Hill numbers(\code{type = "D"}).
#' @param profile two kinds of profiles to show phylogenetic profiles based on chosen \code{type}: q profile(\code{profile = "q"}) and time profile(\code{profile = "time"}). \cr
#' @param q a nonnegative vector specifying the diversity order of Hill numbers. For phylogenetic diversity profiles,
#' \code{q} need to be a vector of length more than one, default is seq(0, 2, by = 0.25) \cr
#' @param tprofile_times a nonnegative vector specifying the reference time points for which phylogenetic profiles will be computed.\cr
#' @param nboot an integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap; in this case,
#' the computation of s.e. and confidence intervals will be skipped, default is 50.
#' @param knots an integer specifying the number of points of \code{time} used to calculate AUC (Area Under Curve). Use 0 to skip calculation, default is 0\cr
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' @param reftime a positive number specifying the reference time.
#' @import chaoUtility
#' @import ggplot2
#' @import dplyr
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom ape drop.tip
#' @importFrom phyclust get.rooted.tree.height
#' @return Based on chosen \code{knots}, there are two kinds of output:\cr\cr
##' \itemize{
##'  \item{For \code{knots} an integer, a list of fourth objects:} \cr
#' \code{$summary} summary of data. \cr\cr
#' \code{$fortime_table} or \code{$forq_table} a table of observed phylogenetic diversity for each \code{time} or order \code{q} based on \code{type} and \code{profile}. \cr\cr
#' \code{$fortime_figure} or \code{$forq_figure} showing the curves for \code{$fortime_table} or \code{$forq_table} \cr\cr
#' \code{$AUC_table} table of Area Under Curve (AUC) of time profile for each site from 0.01 to root. \cr\cr
##'  \item{For \code{knots = NULL}, a list of three objects:} \cr
#' \code{$summary} for summarizing data information. \cr\cr
#' \code{$fortime_table} or \code{$forq_table} a table of observed estimates for each \code{time} or order \code{q} based on \code{type} and \code{profile}. \cr\cr
#' \code{$fortime_figure} or \code{$forq_figure} showing the curves for \code{$fortime_table} or \code{$forq_table} \cr\cr
##' }
#' @examples
#' \donttest{
#' # Type (1) abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- PhdObs_yhc(data = data, tree = tree, datatype = "abundance", type = "PD", profile = "q", q = seq(0, 2, by = 0.25))
#' # Type (2) incidence data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' t_ <- data.inc$t
#' out <- PhdObs_yhc(data = data, tree = tree, datatype = "incidence_raw", type = "PD", profile = "q", q = seq(0, 2, by = 0.25))
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609. \cr\cr
#' @export
PhdObs_yhc <- function(data, tree, datatype = "abundance", t_, type = "PD", profile = "q", q = seq(0, 2, by = 0.25), tprofile_times = NULL,
                       nboot = 50,conf = 0.95,knots = 0,reftime = NULL){
  dat = list()
  if (length(q) == 1) stop("length of q should be greater than one", call. = FALSE)
  if (sum(q<0)>=1) stop("q must be a positive number", call. = FALSE)
  if ((datatype != "incidence_raw") & (datatype != "abundance")) stop("invalid datatype", call. = FALSE)
  if ((type != "PD") & (type != "D")) stop("invalid type", call. = FALSE)
  if ((profile != "q") & (profile != "time")) stop("invalid profile", call. = FALSE)
  if (length(tprofile_times) == 1 & is.null(tprofile_times)==F) stop("length of time should be greater than one", call. = FALSE)
  if (sum(tprofile_times<0)>=1 & is.null(tprofile_times)==F) stop("time must be a positive number", call. = FALSE)
  if (is.null(knots) ==F) {
    if ((knots < 0) | (is.numeric(knots)==F) | (knots%%1>0)) {
      stop('knot must be a nonnegative integer, We use "knots" = 50 to calculate!', call. = FALSE)
    }
  }
  if(is.null(rownames(data))) stop("The rownames (species names) of data can't be empty. Species names in data must match those in the phylogenetic tree")
  pool.name <- rownames(data[rowSums(data)>0,])
  #if(is.null(tprofile_times)){ tprofile_times <- "unspecified" }
  if(datatype=="incidence_raw"){
    if(ncol(data) != sum(t_)) stop("Number of columns does not euqal to the sum of key(t_) in sampling units", call. = FALSE)
    ntmp <- 0
    for(i in 1:length(t_)){
      dat[[i]] <- data[,(ntmp+1):(ntmp+t_[i])]
      ntmp <- ntmp+t_[i]
    }
    if(is.null(names(t_))) {
      names(dat) <- paste0("site",1:length(t_))
    }else{
      names(dat) = names(t_)
    }
  }else if (datatype == "abundance"){
    if(is.null(colnames(data))) {colnames(data) <- paste0("site",1:ncol(data))}
    dat <- lapply(1:ncol(data), function(i)  {x <- data[pool.name,i];names(x) <- pool.name;x})
    names(dat) = colnames(data)
  }

  ###
  if (datatype=="incidence_raw") {
    mydata = lapply(dat, function(X) X[pool.name, ])
  }
  if (datatype=="abundance") {
    mydata = lapply(dat, function(X) X[pool.name])
  }
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)$treeH
  if(is.null(reftime)) {reft <- H_max
  }else if(reftime<=0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }else {reft <- reftime}

  FUN = function(e){
    ###########data information
    if(class(mydata) == "list"){
      infos <- sapply(mydata, function(x){
        datainf_yhc(data = x, datatype, phylotr = mytree,reft = reft)})
    }else{
      return(NULL)
    }
    ###########profiles
    if(profile == "q") {
      temp <- Phdqtable_yhc(datalist = mydata, phylotr = mytree, q, cal = type, datatype, nboot, conf, reft)
      ans <- list(summary = infos, forq_table = temp, forq_figure = Plotq_yhc(temp, type))
      class(ans) <- c("PhdObs")
      return(ans)
    }
    if(profile == "time") {
      if (is.null(tprofile_times)) {
        tprofile_times <- seq(0.01, H_max, length.out = 15) %>% unique() %>% sort
      } else {
        tprofile_times <- c(tprofile_times, 0.01, H_max) %>% unique() %>% sort
      }
      temp <- Phdttable_yhc(datalist = mydata, phylotr = mytree, times = tprofile_times,cal = type, datatype, nboot, conf)
      if (knots==0) {
        ans <- list(summary = infos, fortime_table = temp[[1]], fortime_figure = Plott_yhc(temp[[1]], type, temp[[2]]))
      } else {
        AUC <- AUC_one_table_yhc(datalist = mydata,phylotr = mytree,knot = knots,cal = type,datatype = datatype,nboot, conf,reft_max = max(tprofile_times))
        ans <- list(summary = infos, fortime_table = temp[[1]], fortime_figure = Plott_yhc(temp[[1]], type, temp[[2]]), AUC_table = AUC)
      }
      #class(ans) <- c("PhdObs")
      return(ans)
    }
  }
  temp <- tryCatch(FUN(e), error = function(e){return()})
  return(temp)
}

#' @useDynLib PhDyhc
#' @importFrom Rcpp sourceCpp
NULL
