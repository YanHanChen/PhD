#' Interpolation and extrapolation of phylogenetic diversity
#'
#' \code{iNEXTPD}: the seamless rarefaction and extrapolation sampling curves of phylogenetic diversity(PD) for q = 0, 1 and 2.
#' See Chao et al. (2010, 2015) and Hsieh and Chao (2017) for pertinent background and methods.
#' @param data a matrix/data.frame of species abundances/incidences data.\cr Type (1) abundance data: a S by N matrix/data.frame
#' where N is the number of assemblages. The element in i-th row and k-th is the abundance of species i in assemblage k. Please note
#' that the rownames of data must be the species names matching the species names in phylogeny tree and thus can't be empty.\cr
#' Type (2) incidence data: the sampling unit is quadrat or transect, the observed species was only recorded as presence(detection)/absence(non-detection)
#' data in each sampling unit. For single assemblage consisting of t_ sampling units, data is a S by t_ incidence matrix,
#' where S is the number of species. When there are N assemblages, users must first merge N data matrices by species identity to
#' obtain a large incidence matrix, where the rows of the matrix refer to all species presented in the pooled data. Likewise,
#' the rownames of data must be the species names matching the species names in phylogeny tree and thus can't be empty. \cr
#' @param tree a phylo object describing the Newick phylogeny tree for all observed species in the pooled assemblage.\cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}), default is \code{"abundance"}. \cr
#' @param t_ needed only when datatype = "incidence_raw", a sequence of named nonnegative integers specifying the sampling units in each community. Ignored if \code{datatype = "abundance"}.\cr
#' @param q a sequence of nonnegative integers specifying the diversity orders of PD. Default is \code{c(0,1,2)}. \cr
#' @param endpoint an positive interger specifying the endpoint for rarefaction and
#' extrapolation range. If \code{NULL}, \code{endpoint} = double of the maximum reference sample size. It will be ignored if \code{size} is given. \cr
#' @param knots a positive integer specifying the number of knots between 1 and the \code{endpoint}. Default is 40.\cr
#' @param size a sequence of positive integers specifying the sample sizes for which PD estimates will be calculated. If \code{NULL}, then PD estimates will be
#' calculated for those sample sizes determined by the specified/default \code{endpoint} and \code{knots}. \cr
#' @param plot.type a positive integer vector specifying types of curves. Three types of plots: sample-size-based rarefaction and extrapolation curve (\code{plot.type = 1});
#' coverage-based rarefaction and extrapolation curve (\code{plot.type = 2}); sample completeness curve (\code{plot.type = 3}). Default is \code{c(1,2,3)}. \cr
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
#' @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
#' in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
#' @param reftime a positive value specifying the reference time for tree. If \code{NULL}, \code{reftime} = the tree depth of pooled assemblage. Default is \code{NULL}.
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
#' \code{$summary} data information summary including reference sample size (\code{n}), number of observed species (\code{S.obs}),
#' observed branch length (Faith’s PD) (\code{PD.obs}), the first two species frequency counts (\code{f1*}, \code{f2*}) and their branch length sums,
#' respectively (\code{g1}, \code{g2}). \cr\cr
#' \code{$reference_time} the reference time for the caculation. If \code{reftime} is not specified by users or \code{NULL}, \code{$reference_time}
#' is the \code{tree} depth of pooled assemblage by default. \cr\cr
#' \code{$inext} a table of PD estimates and sample completeness for interpolated or extrapolated sample sizes along with their confidence intervals (if \code{nboot > 0}). \cr\cr
#' \code{$figure} plots of the estimated curves based on \code{$inext} using \code{ggplot2} package. Different choice of \code{plot.type} will yied different types of plot:
#' \itemize{
#'  \item{Sample-size-based R/E curve (\code{plot.type = 1}): the curve of PD estimates with confidence intervals (if \code{nboot > 0}) as a function of sample size.} \cr
#'  \item{Coverage-based R/E curve (\code{plot.type = 2}): the curve of the PD estimates with confidence intervals (if \code{nboot > 0}) as a function of sample coverage.} \cr
#'  \item{Sample completeness curve (\code{plot.type = 3}): the curve of sample coverage with respect to sample size.}
#' }
#' @examples
#' \donttest{
#' # Type (1) abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- iNEXTPD(data = data, tree = tree,datatype = "abundance",q = c(0,1,2), plot.type = 1:3)
#' # Type (2) incidence data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' t_ <- data.inc$t
#' out <- iNEXTPD(data = data, tree = tree,datatype = "incidence_raw",
#' t_ = t_,q = c(0,1,2), plot.type = 1:3)
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609.\cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015) Rarefaction and extrapolation of phylogenetic diversity. Methods in Ecology and Evolution, 6, 380-388.\cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. Systematic Biology 66, 100-111.
#' @export
iNEXTPD <- function(data, tree, datatype = "abundance", t_, q = c(0,1,2), endpoint = NULL, knots = 40, size = NULL, plot.type = 1:3, conf = 0.95, nboot = 50,reftime=NULL) {
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  DATATYPE <- c("abundance", "incidence_raw")
  if(is.na(pmatch(datatype, DATATYPE)) == T)
    stop("invalid datatype", call. = FALSE)
  if(class(q) != "numeric")
    stop("invlid class of order q, q should be a nonnegative integer value", call. = FALSE)
  # qq <- 0:2
  # if(is.na(pmatch(q, qq) == T) == T) stop("invalid order of q, we only compute q = 0, 1 or 2", call. = FALSE)

  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf (confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- as.matrix(data)
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can't be empty.", call. = FALSE)
  dat = list()
  data <- data[rowSums(data)>0,,drop=FALSE]
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
    H_max <- get.rooted.tree.height(mytree)
  }
  if (datatype=="abundance") {
    pool.name = unique(unlist(lapply(dat, function(x) names(x)[x>0] )))#200304 bugs coccurs when input contains rowsums==0.
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
      datainf(data = x, datatype, phylotr = mytree,reft = reft)})
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
      temp = inextPD(datalist=mydata,phylotr=mytree,datatype,Q,nboot,conf,size,knots,endpoint,reft)
      temp$qPD.LCL[temp$qPD.LCL<0] <- 0;temp$SC.LCL[temp$SC.LCL<0] <- 0
      temp$SC.UCL[temp$SC.UCL>1] <- 1
      return(temp)
    }else{
      return(NULL)
    }
  }
  #atime <- Sys.time()
  RE.table <- FUN(3)
  #RE.table <- tryCatch(FUN(e), error = function(e){return()})
  #btime <- Sys.time()
  #print(paste0('R/E time:',btime-atime))
  ###############output3
  TYPE <- 1:3
  FUN2 <- function(e){
    if(is.na(sum(pmatch(plot.type, TYPE))) == F){
      temp2 <- lapply(plot.type, function(j) RE_plot(RE.table, datatype, j))
      allname <- c("RE.plot.size", "RE.plot.C", "RE.plot.sizeC")
      names(temp2) <- allname[plot.type]
      temp2
    }else{
      return("invalid plot type", call. = FALSE)
    }
  }
  #atime <- Sys.time()
  RE.plot <- FUN2(3)
  #RE.plot <- tryCatch(FUN2(e), error = function(e){return()})
  #btime <- Sys.time()
  #print(paste0('plot time:',btime-atime))
  ans <- list(summary = infos,reference_time = reft, inext = RE.table, figure = RE.plot)
  #class(ans) <- c("iNEXTPD")
  ans

}

#' Asymptotic phylogenetic diversity q profile
#'
#' \code{PhdAsy}: estimates asymptotic diversity profile with respect to order \code{q} to infer true phylogenetic diversity(PD) profile
#' at fixed reference time (tree depth by default). It is based on statistical estimation of the unknown true phylogenetic diversity;
#' see Chao et al. (2015) and Hsieh and Chao (2017) for the statistical estimation detail.\cr
#' @param data a matrix/data.frame of species abundances/incidences data.\cr
#' See \code{\link{iNEXTPD}} for data details.
#' @param tree a phylo object describing the Newick phylogeny tree for all observed species in the pooled assemblage. \cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}). Default is "abundance". \cr
#' @param t_ needed only when datatype = "incidence_raw", a sequence of named nonnegative integers specifying the sampling units in each community. Ignored if \code{datatype = "abundance"}.\cr
#' Ignored if \code{datatype = "abundance"}.
#' @param q a nonnegative sequence specifying the diversity orders of PD. Default is seq(0, 2, by = 0.25). \cr
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
#' @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
#' in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
#' @param reftime a positive value specifying the reference time for tree. If \code{NULL}, \code{reftime} = the tree depth of pooled assemblage. Default is \code{NULL}.
#' @import chaoUtility
#' @import ggplot2
#' @import dplyr
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom ape drop.tip
#' @importFrom phyclust get.rooted.tree.height
#' @return  a list of four objects: \cr\cr
#' \code{$summary} summary of data; see \code{\link{iNEXTPD}} for summary details. \cr\cr
#' \code{$reference_time} the reference time for the caculation; see \code{\link{iNEXTPD}} for details.  \cr\cr
#' \code{$asy} table of estimated asymptotic PD q profile at fixed \code{reftime} along with their
#' confidence intervals (if \code{nboot > 0}). It is based on statistical estimation of the unknown true phylogenetic diversity \cr\cr
#' \code{$figure} plot of the estimated q profile based on \code{$asy} using \code{ggplot2} package. \cr\cr
#' @examples
#' \donttest{
#' # Type (1) abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- PhdAsy(data = data, tree = tree, q = seq(0, 2, by = 0.25))
#' # Type (2) incidence data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' t_ <- data.inc$t
#' out <- PhdAsy(data = data, tree = tree, datatype = "incidence_raw",
#' t_ = t_, q = seq(0, 2, by = 0.25))
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609. \cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015) Rarefaction and extrapolation of phylogenetic diversity. Methods in Ecology and Evolution, 6, 380-388.\cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. Systematic Biology 66, 100-111.
#' @export
PhdAsy <- function(data, tree, datatype = "abundance", t_, q = seq(0, 2, by = 0.25), conf = 0.95, nboot = 50, reftime = NULL){
  dat = list()
  if(is.null(rownames(data))) stop("The rownames (species names) of data can't be empty. Species names in data must match those in the phylogenetic tree")
  pool.name <- rownames(data[rowSums(data)>0,,drop=FALSE])
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
  H_max <- get.rooted.tree.height(mytree)
  if(is.null(reftime)) {reft <- H_max
  }else if(reftime<=0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }else {reft <- reftime}

  FUN = function(e){
    temp = AsyPD(datalist = mydata, datatype, phylotr = mytree, q, nboot, conf, reft)# mytree is pooled tree of class phylo
    ans <- list(summary = temp[[2]], reference_time = reft, asy = temp[[1]], figure = Asy_plot(temp[[1]], 1))
    #class(ans) <- c("PhdAsy")
    return(ans)
  }
  #out <- FUN(3)
  out <- tryCatch(FUN(e), error = function(e){return()})
  return(out)
}

#' Observed phylogenetic diversity or phylogenetic Hill numbers
#'
#' \code{PhdObs}: computes observed q profile and reference time profile of phylogenetic diversity(PD) or phylogenetic Hill numbers(D);
#' see Chao et al. (2010) for details of PD and D.\cr
#' @param data a matrix/data.frame of species abundances/incidences data.\cr
#' See \code{\link{iNEXTPD}} for data details.
#' @param tree a phylo object describing the Newick phylogeny tree for all observed species in the pooled assemblage. \cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}). Default is "abundance". \cr
#' @param t_ needed only when datatype = "incidence_raw", a sequence of named nonnegative integers specifying the sampling units in each community. Ignored if \code{datatype = "abundance"}.\cr
#' automatically assigned siteX as community names.
#' Ignored if \code{datatype = "abundance"}.
#' @param type desired diversity type: \code{type = "PD"} for phylogenetic diversity and \code{type = "D"} for phylogenetic Hill numbers. See Chao et al. (2010) for details. Default is \code{"PD"}.
#' @param profile desired profile type for chosen type, \code{profile = "q"} for order q profile and \code{profile = "time"} for reference time profile. Default is \code{"q"}. \cr
#' @param q a nonnegative sequence specifying the diversity orders for q profile if \code{profile == "q"}. Will be ignored if \code{profile = "time"}. Default is \code{seq(0, 2, by = 0.25)}.
#' @param tprofile_times a nonnegative sequence specifying the reference time points for time profile
#' if \code{profile = "time"}. Will be ignored if \code{profile = "q"}. If \code{NULL}, appropriate time pionts will be used based on the \code{tree} depth. Default is \code{NULL}.\cr
#' @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
#' in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' @param knots a positive integer specifying the number of reference time points used to calculate AUC (Area Under Curve) of time profile if \code{profile = "time"}.
#' Larger knots would yield more pricise integration but will take more time in calculation. Enter 0 to skip calculation. Will be ignored if \code{profile = "q"}. Default is 0.\cr
#' @param reftime a positive value specifying the reference time for tree. If \code{NULL}, \code{reftime} = the tree depth of pooled assemblage. Default is \code{NULL}.
#' @import chaoUtility
#' @import ggplot2
#' @import dplyr
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom ape drop.tip
#' @importFrom phyclust get.rooted.tree.height
#' @return a list containing three or four objects based on chosen \code{profile}: \cr\cr
#' \code{$summary} summary of data. \cr\cr
#' \code{$fortime_table} or \code{$forq_table} a table of observed phylogenetic diversity(\code{type = "PD"}) or phylogenetic Hill numbers(\code{type = "D"}) for each \code{time} or order \code{q} based on \code{profile}. \cr\cr
#' \code{$fortime_figure} or \code{$forq_figure} the curves for \code{$fortime_table} or \code{$forq_table} \cr\cr
#' \code{$AUC_table} table of Area Under Curve (AUC) of time profile for each site from 0.01 to the tree depth. Won't be provided if \code{knots = 0} or \code{profile = "q"}.  \cr\cr
#' @examples
#' \donttest{
#' # Type (1) abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- PhdObs(data = data, tree = tree, datatype = "abundance", type = "PD",
#'  profile = "q", q = seq(0, 2, by = 0.25))
#' # Type (2) incidence data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' t_ <- data.inc$t
#' out <- PhdObs(data = data, tree = tree, datatype = "incidence_raw", type = "PD",
#' profile = "q", q = seq(0, 2, by = 0.25))
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609. \cr\cr
#' @export
PhdObs <- function(data, tree, datatype = "abundance", t_, type = "PD", profile = "q", q = seq(0, 2, by = 0.25), tprofile_times = NULL,
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
  pool.name <- rownames(data[rowSums(data)>0,,drop=FALSE])
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
  H_max <- get.rooted.tree.height(mytree)
  if(is.null(reftime)) {reft <- H_max
  }else if(reftime<=0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }else {reft <- reftime}

  FUN = function(e){
    ###########data information
    if(class(mydata) == "list"){
      infos <- sapply(mydata, function(x){
        datainf(data = x, datatype, phylotr = mytree,reft = reft)})
    }else{
      return(NULL)
    }
    ###########profiles
    if(profile == "q") {
      temp <- Phdqtable(datalist = mydata, phylotr = mytree, q, cal = type, datatype, nboot, conf, reft)
      ans <- list(summary = infos, forq_table = temp, forq_figure = Plotq(temp, type))
      class(ans) <- c("PhdObs")
      return(ans)
    }
    if(profile == "time") {
      if (is.null(tprofile_times)) {
        tprofile_times <- seq(0.01, H_max, length.out = 15) %>% unique() %>% sort
      } else {
        tprofile_times <- c(tprofile_times, 0.01, H_max) %>% unique() %>% sort
      }
      temp <- Phdttable(datalist = mydata, phylotr = mytree, times = tprofile_times,cal = type, datatype, nboot, conf)
      if (knots==0) {
        ans <- list(summary = infos, fortime_table = temp[[1]], fortime_figure = Plott(temp[[1]], type, temp[[2]]))
      } else {
        AUC <- AUC_one_table(datalist = mydata,phylotr = mytree,knot = knots,cal = type,datatype = datatype,nboot, conf,reft_max = max(tprofile_times))
        ans <- list(summary = infos, fortime_table = temp[[1]], fortime_figure = Plott(temp[[1]], type, temp[[2]]), AUC_table = AUC)
      }
      #class(ans) <- c("PhdObs")
      return(ans)
    }
  }
  #temp <- FUN(3)
  temp <- tryCatch(FUN(e), error = function(e){return()})
  return(temp)
}

#' @useDynLib PhD
#' @importFrom Rcpp sourceCpp
NULL
