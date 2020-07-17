#' Interpolation and extrapolation of phylogenetic diversity/phylogenetic Hill numbers
#'
#' \code{iNEXTPD}: the seamless rarefaction and extrapolation sampling curves of phylogenetic diversity(PD) or phylogenetic Hill numbers(D).
#' See Chao et al. (2010, 2015) and Hsieh and Chao (2017) for pertinent background and methods.
#' @param data a matrix/data.frame of species abundances/incidences data.\cr Type (1) abundance data: a S by N matrix/data.frame
#' where N is the number of assemblages. The element in i-th row and k-th is the abundance of species i in assemblage k. Please note
#' that the rownames of data must be the species names matching the species names in phylogeny tree and thus can not be empty.\cr
#' Type (2) incidence data: the sampling unit is quadrat or transect, the observed species was only recorded as presence(detection)/absence(non-detection)
#' data in each sampling unit. For single assemblage consisting of t_ sampling units, data is a S by t_ incidence matrix,
#' where S is the number of species. When there are N assemblages, users must first merge N data matrices by species identity to
#' obtain a large incidence matrix, where the rows of the matrix refer to all species presented in the pooled data. Likewise,
#' the rownames of data must be the species names matching the species names in phylogeny tree and thus can not be empty. \cr
#' @param t_ needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the sampling units in each assemblage. Ignored if \code{datatype = "abundance"}.\cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}), default is \code{"abundance"}. \cr
#' @param tree a phylo object describing the Newick phylogeny tree for all observed species in the pooled assemblage.\cr
#' @param q a nonnegative value or sequence specifying the diversity order. Default is \code{c(0,1,2)}.
#' @param reftime a positive value or sequence specifying the reference time for tree. If \code{NULL}, \code{reftime} = the tree depth of pooled assemblage. Default is \code{NULL}.
#' @param type desired diversity type: \code{type = "PD"} for phylogenetic diversity and \code{type = "D"} for phylogenetic Hill number. See Chao et al. (2010) for details. Default is \code{"PD"}.
#' @param endpoint an positive interger specifying the endpoint for rarefaction and
#' extrapolation range. If \code{NULL}, \code{endpoint} = double of the maximum reference sample size. It will be ignored if \code{size} is given. \cr
#' @param knots a positive integer specifying the number of knots between 1 and the \code{endpoint}. Default is 40.\cr
#' @param size a sequence of positive integers specifying the sample sizes for which \code{PD} or \code{D} estimates will be calculated. If \code{NULL}, then estimates will be
#' calculated for those sample sizes determined by the specified/default \code{endpoint} and \code{knots}. \cr
#' @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
#' in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
#' @import chaoUtility
#' @import ggplot2
#' @import dplyr
#' @importFrom stats rbinom
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom ape drop.tip
#' @importFrom phyclust get.rooted.tree.height
#' @return a tibble of PD or D estimates and sample coverage for interpolated or extrapolated sample sizes along with their confidence intervals (if \code{nboot > 0}). \cr\cr
#' @examples
#' \donttest{
#' # Type (1) abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- iNEXTPD(data = data, tree = tree,datatype = "abundance",q = c(0,1,2),reftime = c(162.5,325))
#' # Type (2) incidence data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' t_ <- data.inc$t
#' out <- iNEXTPD(data = data, t_ = t_,datatype = "incidence_raw", tree = tree,q = c(0,1,2))
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609.\cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015) Rarefaction and extrapolation of phylogenetic diversity. Methods in Ecology and Evolution, 6, 380-388.\cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. Systematic Biology 66, 100-111.
#' @export
iNEXTPD <- function(data, t_, datatype = "abundance",tree,q = c(0,1,2),reftime=NULL,type = 'PD', endpoint = NULL, knots = 40, size = NULL, nboot = 50, conf = 0.95) {
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
  if(class(data)[1]=="numeric"|class(data)[1]=="integer"|class(data)[1]=="double" ) data <- as.matrix(data)
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)

  data <- data[rowSums(data)>0,,drop=FALSE]
  pool.name <- rownames(data)
  mydata = list()
  if(datatype=="incidence_raw"){
    if(ncol(data) != sum(t_)) stop("Number of columns does not euqal to the sum of key in sampling units", call. = FALSE)
    n <- 0
    for(i in 1:length(t_)){
      mydata[[i]] <- data[,(n+1):(n+t_[i])]
      n <- n+t_[i]
    }
    if(is.null(names(t_))) {
      names(mydata) <- paste0("assemblage",1:length(t_))
    }else{
      names(mydata) = names(t_)
    }
  }else{
    if(is.null(colnames(data))) {colnames(data) <- paste0("assemblage",1:ncol(data))}
    mydata <- lapply(1:ncol(data), function(i)  {x <- data[,i];names(x) <- pool.name;x})
    names(mydata) = colnames(data)
  }
  ###
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)

  # reft <- reftime
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
  #reftime <- ifelse(is.null(reftime),H_max,reftime)
  reftime <- sort(unique(reftime))
  if(sum(reftime<=0)>0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }

  ###########output1
  #atime <- Sys.time()
  # if(class(mydata) == "list"){
  #   infos <- lapply(mydata, function(x){
  #     datainf(data = x, datatype, phylotr = mytree,reft = reft) %>% mutate(Reference.time = reft)
  #     }) %>% do.call(rbind,.) %>% mutate(Assemblage = rep(names(mydata),each = length(reft))) %>%
  #     select(Assemblage,n,S.obs,PD.obs,`f1*`,`f2*`,g1,g2,Reference.time)
  # }else{
  #   return(NULL)
  # }
  #btime <- Sys.time()
  #print(paste0('Info time:',btime-atime))
  ############output2
  if(length(knots)!=length(mydata)) knots <- rep(knots,length(mydata))
  if(is.null(size)){
    if(is.null(endpoint)){
      if(datatype == "abundance") {
        endpoint <- sapply(mydata, function(x) 2*sum(x))
      }else if(datatype == "incidence_raw"){
        endpoint <- sapply(mydata, function(x) 2*ncol(x))
      }
    }else{
      if(length(endpoint)!=length(mydata)){
        endpoint <- rep(endpoint,length(mydata))
      }
    }
    size <- lapply(1:length(mydata),function(i){
      if(datatype == "abundance") {
        ni <- sum(mydata[[i]])
      }else if(datatype == "incidence_raw"){
        ni <- ncol(mydata[[i]])
      }

      if(endpoint[i] <= ni){
        mi <- floor(seq(1,endpoint[i],length.out = knots[i]))
      }else{
        mi <- floor(c(seq(1,ni,length.out = floor(knots[i]/2)),
                      seq(ni+1,endpoint[i],length.out = knots[i]-floor(knots[i]/2))))
      }
      unique(mi)
    })
  }else{
    if(class(size)=="numeric"|class(size)=="integer"|class(size)=="double"){
      size <- list(size = size)
    }
    if(length(size)!=length(mydata)) size <- lapply(1:length(mydata), function(x) size[[1]])
    size <- lapply(1:length(mydata),function(i){
      if(datatype == "abundance") {
        ni <- sum(mydata[[i]])
      }else if(datatype == "incidence_raw"){
        ni <- ncol(mydata[[i]])
      }
      if( sum(size[[i]] == ni) == 0 ) mi <- sort(c(ni,size[[i]]))
      else mi <- size[[i]]
      unique(mi)
    })
  }

  FUN <- function(e){
    if(class(mydata)=="list"){
      # temp = inextPD(datalist=mydata,phylotr=mytree,datatype,q,nboot,conf,m = size,reftime,cal = type)
      inextPD(datalist = mydata,datatype = datatype,phylotr = mytree,q = q,reft = reftime,m=size,
                     cal = type,nboot=nboot,conf = conf,unconditional_var = TRUE)
    }else{
      return(NULL)
    }
  }
  #atime <- Sys.time()
  #RE.table <- FUN(3)
  ans <- tryCatch(FUN(e), error = function(e){return()})
  #btime <- Sys.time()
  #print(paste0('R/E time:',btime-atime))
  ###############output3
  # TYPE <- 1:3
  # FUN2 <- function(e){
  #   if(is.na(sum(pmatch(plot.type, TYPE))) == F){
  #     temp2 <- lapply(plot.type, function(j) RE_plot(RE.table, datatype, j))
  #     allname <- c("RE.plot.size", "RE.plot.C", "RE.plot.sizeC")
  #     names(temp2) <- allname[plot.type]
  #     temp2
  #   }else{
  #     return("invalid plot type", call. = FALSE)
  #   }
  # }
  # #atime <- Sys.time()
  # RE.plot <- tryCatch(FUN2(e), error = function(e){return()})
  #btime <- Sys.time()
  #print(paste0('plot time:',btime-atime))
  # ans <- list(summary = infos, inext = RE.table, figure = RE.plot)
  #class(ans) <- c("iNEXTPD")
  ans

}

#' Asymptotic phylogenetic diversity/phylogenetic Hill numbers
#'
#' \code{PhdAsy}: estimates asymptotic diversity with respect to specified/default order \code{q} and reference time \code{reftime} to infer true phylogenetic diversity(PD) or phylogenetic Hill numbers(D).
#' It is based on statistical estimation of the unknown true PD or D; see Chao et al. (2015) and Hsieh and Chao (2017) for the statistical estimation detail.\cr
#' @param data a matrix/data.frame of species abundances/incidences data.\cr
#' See \code{\link{iNEXTPD}} for data details.
#' @param t_ needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the sampling units in each assemblage. Ignored if \code{datatype = "abundance"}.\cr
#' Ignored if \code{datatype = "abundance"}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}). Default is \code{"abundance"}. \cr
#' @param tree a phylo object describing the Newick phylogeny tree for all observed species in the pooled assemblage. \cr
#' @param q a nonnegative value or sequence specifying the diversity order. Default is \code{seq(0, 2, by = 0.25)}.
#' @param reftime a positive value or sequence specifying the reference time for tree. If \code{NULL}, \code{reftime} = the tree depth of pooled assemblage. Default is \code{NULL}.
#' @param type desired diversity type: \code{type = "PD"} for phylogenetic diversity and \code{type = "D"} for phylogenetic Hill number. See Chao et al. (2010) for details. Default is \code{"PD"}.
#' @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
#' in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' @import chaoUtility
#' @import ggplot2
#' @import dplyr
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom ape drop.tip
#' @importFrom phyclust get.rooted.tree.height
#' @return a tibble of estimated(asymptotic) phylogenetic diversity(\code{type = "PD"}) or phylogenetic Hill number(\code{type = "D"})
#' with respect to specified/default order \code{q} and reference time \code{reftime}.\cr\cr
#' @examples
#' \donttest{
#' # Type (1) abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- PhdAsy(data = data, datatype = "abundance", tree = tree, q = seq(0, 2, by = 0.25))
#' # Type (2) incidence data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' t_ <- data.inc$t
#' out <- PhdAsy(data = data, t_ = t_,datatype = "incidence_raw", tree = tree,q = seq(0, 2, by = 0.25))
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609. \cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015) Rarefaction and extrapolation of phylogenetic diversity. Methods in Ecology and Evolution, 6, 380-388.\cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. Systematic Biology 66, 100-111.
#' @export
PhdAsy <- function(data, t_, datatype = "abundance", tree, q = seq(0, 2, by = 0.25), reftime = NULL, type = 'PD', nboot = 50, conf = 0.95){
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  #if (length(q) == 1) stop("length of q should be greater than one", call. = FALSE)
  if (sum(q<0)>=1) stop("q must be a positive number", call. = FALSE)
  if ((datatype != "incidence_raw") & (datatype != "abundance")) stop("invalid datatype", call. = FALSE)
  if ((conf < 0) | (conf < 0) | (is.numeric(conf)==F)) stop('conf"(confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- as.matrix(data)
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)

  data <- data[rowSums(data)>0,,drop=FALSE]
  pool.name <- rownames(data)
  mydata = list()
  if(datatype=="incidence_raw"){
    if(ncol(data) != sum(t_)) stop("Number of columns does not euqal to the sum of key(t_) in sampling units", call. = FALSE)
    ntmp <- 0
    for(i in 1:length(t_)){
      mydata[[i]] <- data[,(ntmp+1):(ntmp+t_[i])]
      ntmp <- ntmp+t_[i]
    }
    if(is.null(names(t_))) {
      names(mydata) <- paste0("assemblage",1:length(t_))
    }else{
      names(mydata) = names(t_)
    }
  }else if (datatype == "abundance"){
    if(is.null(colnames(data))) {colnames(data) <- paste0("assemblage",1:ncol(data))}
    mydata <- lapply(1:ncol(data), function(i)  {x <- data[,i];names(x) <- pool.name;x})
    names(mydata) = colnames(data)
  }
  ###
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)

  # reft <- reftime
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
  #reftime <- ifelse(is.null(reftime),H_max,reftime)
  reftime <- sort(unique(reftime))
  if(sum(reftime<=0)>0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }

  # if(class(mydata) == "list"){
  #   infos <- sapply(mydata, function(x){
  #     datainf(data = x, datatype, phylotr = mytree,reft = reft)})
  # }else{
  #   return(NULL)
  # }

  FUN = function(e){
    AsyPD(datalist = mydata,datatype = datatype,phylotr = mytree,q = q,reft = reftime,cal = type,nboot,conf)# mytree is pooled tree of class phylo
  }
  #out <- FUN(3)
  ans <- tryCatch(FUN(e), error = function(e){return()})
  return(ans)
}

#' Observed phylogenetic diversity/phylogenetic Hill numbers
#'
#' \code{PhdObs}: computes emprircal phylogenetic diversity(PD) or phylogenetic Hill numbers(D) for specified/default order \code{q} and reference time \code{reftime};
#' see Chao et al. (2010) for details of PD and D.\cr
#' @param data a matrix/data.frame of species abundances/incidences data.\cr
#' See \code{\link{iNEXTPD}} for data details.
#' @param t_ needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the numbers of sampling units in each assemblage. Ignored if \code{datatype = "abundance"}.\cr
#' If \code{names(t_) = NULL}, automatically assign siteX as assemblage names. Ignored if \code{datatype = "abundance"}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}). Default is \code{"abundance"}. \cr
#' @param tree a phylo object describing the Newick phylogeny tree for all observed species in the pooled assemblage. \cr
#' @param q a nonnegative value or sequence specifying the diversity order. Default is \code{seq(0, 2, by = 0.25)}.
#' @param reftime a positive value or sequence specifying the reference time for tree. If \code{NULL}, \code{reftime} = the tree depth of pooled assemblage. Default is \code{NULL}.
#' @param type desired diversity type: \code{type = "PD"} for phylogenetic diversity and \code{type = "D"} for phylogenetic Hill number. See Chao et al. (2010) for details. Default is \code{"PD"}.
#' @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
#' in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95.
#' @import chaoUtility
#' @import ggplot2
#' @import dplyr
#' @importFrom stats rbinom
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom ape drop.tip
#' @importFrom phyclust get.rooted.tree.height
#' @return a tibble of empirical(observed) phylogenetic diversity(\code{type = "PD"}) or
#' phylogenetic Hill number(\code{type = "D"}) with respect to each \code{reftime} and order \code{q}.\cr\cr
#' @examples
#' \donttest{
#' # Type (1) abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- PhdObs(data = data,datatype = "abundance",tree = tree,
#' q = seq(0, 2, by = 0.25),reftime = c(162.5,325))
#' # Type (2) incidence data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' t_ <- data.inc$t
#' out <- PhdObs(data = data,t_ = t_, datatype = "incidence_raw",
#' tree = tree, q = seq(0, 2, by = 0.25))
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609. \cr\cr
#' @export
PhdObs <- function(data, t_, datatype = "abundance", tree, q = seq(0, 2, by = 0.25), reftime = NULL, type = "PD",
                       nboot = 50,conf = 0.95){
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  #if (length(q) == 1) stop("length of q should be greater than one", call. = FALSE)
  if (sum(q<0)>0) stop("q must be a positive number", call. = FALSE)
  if ((datatype != "incidence_raw") & (datatype != "abundance")) stop("invalid datatype", call. = FALSE)
  if ((type != "PD") & (type != "D")) stop("invalid type", call. = FALSE)
  # if ((profile != "q") & (profile != "time")) stop("invalid profile", call. = FALSE)
  # if (length(tprofile_times) == 1 & is.null(tprofile_times)==F) stop("length of time should be greater than one", call. = FALSE)
  # if (sum(tprofile_times<0)>=1 & is.null(tprofile_times)==F) stop("time must be a positive number", call. = FALSE)
  # if (is.null(knots) ==F) {
  #   if ((knots < 0) | (is.numeric(knots)==F) | (knots%%1>0)) {
  #     stop('knot must be a nonnegative integer, We use "knots" = 50 to calculate!', call. = FALSE)
  #   }
  # }
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- as.matrix(data)
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)

  data <- data[rowSums(data)>0,,drop=FALSE]
  pool.name <- rownames(data)
  mydata = list()

  if(datatype=="incidence_raw"){
    if(ncol(data) != sum(t_)) stop("Number of columns does not euqal to the sum of key(t_) in sampling units", call. = FALSE)
    ntmp <- 0
    for(i in 1:length(t_)){
      mydata[[i]] <- data[,(ntmp+1):(ntmp+t_[i])]
      ntmp <- ntmp+t_[i]
    }
    if(is.null(names(t_))) {
      names(mydata) <- paste0("assemblage",1:length(t_))
    }else{
      names(mydata) = names(t_)
    }
  }else if (datatype == "abundance"){
    if(is.null(colnames(data))) {colnames(data) <- paste0("assemblage",1:ncol(data))}
    mydata <- lapply(1:ncol(data), function(i)  {x <- data[,i];names(x) <- pool.name;x})
    names(mydata) = colnames(data)
  }

  ###
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)

  # reft <- reftime
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
  #reftime <- ifelse(is.null(reftime),H_max,reftime)
  reftime <- sort(unique(reftime))
  if(sum(reftime<=0)>0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }
  #=====old version=====
  # FUN = function(e){
  #   ###########data information
  #   if(class(mydata) == "list"){
  #     infos <- sapply(mydata, function(x){
  #       datainf(data = x, datatype, phylotr = mytree,reft = reftime)})
  #   }else{
  #     return(NULL)
  #   }
  #
  #   if(profile == "q") {
  #
  #     if(is.null(reft)){
  #       if (datatype=="incidence_raw") {
  #         da <- lapply(mydata, rowSums) %>% do.call(cbind, .) %>% rowSums()
  #       }else if (datatype=="abundance") {
  #         da <- do.call(cbind, mydata) %>% rowSums()
  #       }
  #       aL <- phyBranchAL_Abu(phylo = mytree,data = da,"abundance",
  #                             refT = reftime)$treeNabu %>%
  #         select(branch.abun,branch.length,tgroup)
  #       PD2 <- PD.qprofile(aL,q = 2, cal =  "PD",nt = sum(da))
  #       Q <- reftime-(reftime^2)/PD2
  #       reft = sort(c('Q'= Q, 'reftime' = reftime))
  #     }else{
  #       names(reft) <- NULL
  #       reft <- sort(reft)
  #     }
  #
  #     temp <- Phdqtable(datalist = mydata, phylotr = mytree, q, cal = type, datatype, nboot, conf, reft)
  #     ans <- list(summary = infos, forq_table = temp, forq_figure = Plotq(temp, type))
  #     #class(ans) <- c("PhdObs")
  #     return(ans)
  #   }
  #   if(profile == "time") {
  #     if (is.null(tprofile_times)) {
  #       tprofile_times <- seq(0.01, reftime, length.out = 15) %>% unique() %>% sort
  #     } else {
  #       tprofile_times <- c(tprofile_times, 0.01, reftime) %>% unique() %>% sort
  #     }
  #     temp <- Phdttable(datalist = mydata, phylotr = mytree, times = tprofile_times,cal = type, datatype, nboot, conf)
  #     if (knots==0) {
  #       ans <- list(summary = infos, fortime_table = temp[[1]], fortime_figure = Plott(temp[[1]], type, temp[[2]]))
  #     } else {
  #       AUC <- AUC_one_table(datalist = mydata,phylotr = mytree,knot = knots,cal = type,datatype = datatype,nboot, conf,reft_max = max(tprofile_times))
  #       ans <- list(summary = infos, fortime_table = temp[[1]], fortime_figure = Plott(temp[[1]], type, temp[[2]]), AUC_table = AUC)
  #     }
  #     #class(ans) <- c("PhdObs")
  #     return(ans)
  #   }
  # }
  #=====new version=====
  FUN <- function(e){
    EmpPD(datalist = mydata,datatype = datatype,phylotr = mytree,q = q,reft = reftime,cal = type,nboot,conf)
  }
  #temp <- FUN(3)
  ans <- tryCatch(FUN(e), error = function(e){return()})
  return(ans)
}


#' Compute phylogenetic diversity with a particular sample coverages.
#'
#' \code{EstimatePD}: estimates phylogenetic diversity(PD) or phylogenetic Hill numbers(D) with particular user-specified levels of sample coverage with respect to specified/default order \code{q} and reference time \code{reftime}.
#' @param data a matrix/data.frame of species abundances/incidences data.\cr
#' See \code{\link{iNEXTPD}} for data details.
#' @param t_ needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the sampling units in each assemblage. Ignored if \code{datatype = "abundance"}.\cr
#' Ignored if \code{datatype = "abundance"}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}). Default is "abundance". \cr
#' @param tree a phylo object describing the Newick phylogeny tree for all observed species in the pooled assemblage. \cr
#' @param q a nonnegative value or sequence specifying the diversity order. Default is \code{c(0,1,2)}.
#' @param reftime a positive value or sequence specifying the reference time for tree. If \code{NULL}, \code{reftime} = the tree depth of pooled assemblage. Default is \code{NULL}.
#' @param type desired diversity type: \code{type = "PD"} for phylogenetic diversity and \code{type = "D"} for phylogenetic Hill number. See Chao et al. (2010) for details. Default is \code{"PD"}.
#' @param level a positive value or sequence < 1 specifying a particular value of sample coverage.
#' If \code{NULL},then \code{level} will be chosen as the minimum coverage of all sites after extrapolating each assemblage to its double sample sizes. Default is \code{NULL}.
#' @param nboot a positive integer specifying the number of bootstrap replications. Enter 0 to skip bootstrap;
#' in this case, the caculation of standard errors and confidence intervals will be skipped. Default is 50.
#' @param conf a positive number < 1 specifying the level of confidence interval, default is 0.95. \cr
#' @import chaoUtility
#' @import dplyr
#' @importFrom stats qnorm
#' @importFrom stats sd
#' @importFrom stats optimize
#' @importFrom ape drop.tip
#' @importFrom phyclust get.rooted.tree.height
#' @return a tibble of the PD or D estimates with respect to specified/default order \code{q} and reference time \code{reftime} for the user-specified sample coverages.
#' In addition,the corresponding sample sizes and sample coverages are also provided. \cr\cr
#' @examples
#' \donttest{
#' # Type (1) abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- estimatePD(data = data, datatype = "abundance",tree = tree)
#' # Type (2) incidence data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' t_ <- data.inc$t
#' out <- estimatePD(data = data,t_ = t_, datatype = "incidence_raw", tree = tree)
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609. \cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015) Rarefaction and extrapolation of phylogenetic diversity. Methods in Ecology and Evolution, 6, 380-388.\cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. Systematic Biology 66, 100-111.
#' @export
estimatePD <- function(data,t_,datatype = "abundance", tree, q = c(0,1,2), reftime=NULL,type = 'PD',level = NULL, nboot = 50,conf = 0.95){
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  if (sum(q<0)>=1) stop("q must be a positive number", call. = FALSE)
  if ((datatype != "incidence_raw") & (datatype != "abundance")) stop("invalid datatype", call. = FALSE)
  if ((conf < 0) | (conf > 1) | (is.numeric(conf)==F)) stop('conf"(confidence level) must be a numerical value between 0 and 1, We use "conf" = 0.95 to calculate!', call. = FALSE)
  if ((nboot < 0) | (is.numeric(nboot)==F)) stop('nboot must be a nonnegative integer, We use "nboot" = 50 to calculate!', call. = FALSE)
  #if (length(level)>1) stop('Currently, we only accept one fixed level of coverage.')
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- as.matrix(data)
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)

  data <- data[rowSums(data)>0,,drop=FALSE]
  pool.name <- rownames(data)
  mydata = list()
  if(datatype=="incidence_raw"){
    if(ncol(data) != sum(t_)) stop("Number of columns does not euqal to the sum of t_ (number of sampling units for each assemblage).", call. = FALSE)
    ntmp <- 0
    for(i in 1:length(t_)){
      mydata[[i]] <- data[,(ntmp+1):(ntmp+t_[i])]
      ntmp <- ntmp+t_[i]
    }
    if(is.null(names(t_))) {
      names(mydata) <- paste0("assemblage",1:length(t_))
    }else{
      names(mydata) = names(t_)
    }
  }else if (datatype == "abundance"){
    if(is.null(colnames(data))) {colnames(data) <- paste0("assemblage",1:ncol(data))}
    mydata <- lapply(1:ncol(data), function(i)  {x <- data[pool.name,i];names(x) <- pool.name;x})
    names(mydata) = colnames(data)
  }
  ###

  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)

  # reft <- reftime
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
  #reftime <- ifelse(is.null(reftime),H_max,reftime)
  reftime <- sort(unique(reftime))
  if(sum(reftime<=0)>0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }

  if(is.null(level)){
    if(datatype=='abundance'){
      level <- sapply(mydata,function(x){
        ni <- sum(x)
        Coverage(data = x,datatype = datatype,m = 2*ni,nt = ni)
      })

    }else if(datatype=='incidence_raw'){
      level <- sapply(mydata,function(x){
        ni <- ncol(x)
        Coverage(data = x,datatype = datatype,m = 2*ni,nt = ni)
      })
    }
    level <- min(level)
  }

  out <- invChatPD(datalist = mydata, datatype = datatype,phylotr = mytree, q = q,
                   reft = reftime, cal = type,level = level, nboot, conf)
  return(out)
}


#' Summarize phylogenetic data infomation.
#'
#' \code{PDInfo}: summarize phylogenetic data statistics for specified/default reference time \code{reftime}.
#' @param data a matrix/data.frame of species abundances/incidences data.\cr
#' See \code{\link{iNEXTPD}} for data details.
#' @param t_ needed only when \code{datatype = "incidence_raw"}, a sequence of named nonnegative integers specifying the sampling units in each assemblage. Ignored if \code{datatype = "abundance"}.\cr
#' Ignored if \code{datatype = "abundance"}.
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}). Default is "abundance". \cr
#' @param tree a phylo object describing the Newick phylogeny tree for all observed species in the pooled assemblage. \cr
#' @param reftime a positive value or sequence specifying the reference time for tree. If \code{NULL}, \code{reftime} = the tree depth of pooled assemblage. Default is \code{NULL}.
#' @return a tibble of phylogenetic data statistics, including reference sample size \code{n}/ number of sampling units \code{T},
#' number of observed species \code{S.obs}, observed branch length (Faith's PD) \code{PD.obs}, the first two species frequency counts
#' (\code{f1*/Q1*}, \code{f2*/Q2*}) and their branch length sums, respectively (\code{g1/R1}, \code{g2/R2}) at specified/default reference time \code{reftime}. \cr\cr
#' @examples
#' \donttest{
#' # Type (1) abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- PDInfo(data = data,datatype = "abundance", tree = tree)
#' # Type (2) incidence data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' t_ <- data.inc$t
#' out <- PDInfo(data = data, t_ = t_, datatype = "incidence_raw", tree = tree)
#' }
#' @references
#' Chao, A., Chiu C.-H. and Jost, L. (2010). Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B., 365, 3599-3609. \cr\cr
#' Chao, A., Chiu, C.-H., Hsieh, T. C., Davis, T., Nipperess, D., and Faith, D. (2015) Rarefaction and extrapolation of phylogenetic diversity. Methods in Ecology and Evolution, 6, 380-388.\cr\cr
#' Hsieh, T. C. and Chao, A. (2017). Rarefaction and extrapolation: making fair comparison of abundance-sensitive phylogenetic diversity among multiple assemblages. Systematic Biology 66, 100-111.
#' @export
PDInfo <- function(data,t_,datatype = "abundance", tree, reftime=NULL){
  if(sum(c(duplicated(tree$tip.label),duplicated(tree$node.label[tree$node.label!=""])))>0)
    stop("The phylo tree should not contains duplicated tip or node labels, please remove them.", call. = FALSE)
  if ((datatype != "incidence_raw") & (datatype != "abundance")) stop("invalid datatype", call. = FALSE)
  if(class(data)=="numeric"|class(data)=="integer"|class(data)=="double" ) data <- as.matrix(data)
  if(is.null(rownames(data) ))
    stop("Row names of data must be the species names that match tip names in tree and thus can not be empty.", call. = FALSE)

  data <- data[rowSums(data)>0,,drop=FALSE]
  pool.name <- rownames(data)
  mydata = list()
  if(datatype=="incidence_raw"){
    if(ncol(data) != sum(t_)) stop("Number of columns does not euqal to the sum of key(t_) in sampling units", call. = FALSE)
    ntmp <- 0
    for(i in 1:length(t_)){
      mydata[[i]] <- data[,(ntmp+1):(ntmp+t_[i])]
      ntmp <- ntmp+t_[i]
    }
    if(is.null(names(t_))) {
      names(mydata) <- paste0("assemblage",1:length(t_))
    }else{
      names(mydata) = names(t_)
    }
  }else if (datatype == "abundance"){
    if(is.null(colnames(data))) {colnames(data) <- paste0("assemblage",1:ncol(data))}
    mydata <- lapply(1:ncol(data), function(i)  {x <- data[,i];names(x) <- pool.name;x})
    names(mydata) = colnames(data)
  }
  ###
  tip <- tree$tip.label[-match(pool.name,tree$tip.label)]
  mytree <- drop.tip(tree,tip)
  H_max <- get.rooted.tree.height(mytree)

  # reft <- reftime
  if(is.null(reftime)) reftime <- H_max else reftime <- reftime
  #reftime <- ifelse(is.null(reftime),H_max,reftime)
  reftime <- sort(unique(reftime))
  if(sum(reftime<=0)>0) {stop("Reference time must be greater than 0. Use NULL to set it to pooled tree height.",call. = FALSE)
  }

  if(datatype=='abundance'){
    infos <- lapply(mydata, function(x){
      datainf(data = x, datatype, phylotr = mytree,reft = reftime) %>% mutate(Reference.time = reftime)
    }) %>% do.call(rbind,.) %>% mutate(Assemblage = rep(names(mydata),each = length(reftime))) %>%
      select(Assemblage,n,S.obs,PD.obs,`f1*`,`f2*`,g1,g2,Reference.time)
  }else if (datatype=='incidence_raw'){
    infos <- lapply(mydata, function(x){
      datainf(data = x, datatype, phylotr = mytree,reft = reftime) %>% mutate(Reference.time = reftime)
    }) %>% do.call(rbind,.) %>% mutate(Assemblage = rep(names(mydata),each = length(reftime))) %>%
      select(Assemblage,`T`,S.obs,PD.obs,`Q1*`,`Q2*`,R1,R2,Reference.time)
  }

  return(infos)

}


#' Plot the outcome of iNEXTPD using \code{ggplot2} package.
#'
#' \code{ggiNEXTPD}: plot the outcome of iNEXTPD using ggplot2 package.
#' @param outcome the outcome of the function \code{iNEXTPD}
#' @param plot.type a positive integer vector specifying types of curves. Three types of plots: sample-size-based rarefaction and extrapolation curve (\code{plot.type = 1});
#' coverage-based rarefaction and extrapolation curve (\code{plot.type = 2}); sample coverage curve (\code{plot.type = 3}). Default is \code{c(1,2,3)}. \cr
#' @param datatype data type of input data: individual-based abundance data (\code{datatype = "abundance"}),
#' or species by sampling-units incidence matrix (\code{datatype = "incidence_raw"}). Default is "abundance". \cr
#' @return plots of the estimated curves based on \code{ggplot2} package. Different choice of \code{plot.type} will yied different types of plot:
#' \itemize{
#'  \item{Sample-size-based R/E curve (\code{plot.type = 1}): the curve of estimates as a function of sample size.} \cr
#'  \item{Sample completeness curve (\code{plot.type = 2}): the curve of sample coverage with respect to sample size.} \cr
#'  \item{Coverage-based R/E curve (\code{plot.type = 3}): the curve of the estimates as a function of sample coverage.} \cr
#'  }
#' @examples
#' \donttest{
#' # Type (1) abundance data
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- iNEXTPD(data = data, tree = tree,datatype = "abundance",q = c(0,1,2))
#' ggiNEXTPD(out)
#' # Type (2) incidence data
#' data(data.inc)
#' data <- data.inc$data
#' tree <- data.inc$tree
#' t_ <- data.inc$t
#' out <- iNEXTPD(data = data, t_ = t_,datatype = "incidence_raw", tree = tree,q = c(0,1,2))
#' ggiNEXTPD(out,datatype = "incidence_raw")
#' }
#' @export
ggiNEXTPD <- function(outcome,plot.type = 1:3,datatype = 'abundance'){
  TYPE <- c(1,2,3)
  if(is.na(sum(pmatch(plot.type, TYPE))) == F){
    temp2 <- lapply(plot.type, function(j) RE_plot(outcome, datatype, j))
    allname <- c('RE.plot.size', 'RE.plot.sizeC','RE.plot.C')
    names(temp2) <- allname[plot.type]
    temp2
  }
}

#' Plot the outcome of \code{PhdObs} or \code{PhdAsy} based on \code{ggplot2} package.
#'
#' \code{ggtqplot}: plot the outcome of \code{PhdObs} or \code{PhdAsy} using \code{ggplot2} package.
#' @param outcome the outcome of the function \code{PhdObs} or \code{PhdAsy}.
#' @param profile specifying the type of profile: \code{profile = 'q'} for q profile and \code{profile = 'time'} for time profile.
#' @return plot of the PD/D empirical or estimated curves based on \code{ggplot2} package.
#' @examples
#' \donttest{
#' # Type (1) q profile
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- PhdObs(data = data,datatype = "abundance",tree = tree,q = seq(0, 2, by = 0.25))
#' ggtqplot(out,profile = 'q')
#' # Type (2) time profile
#' data(data.abu)
#' data <- data.abu$data
#' tree <- data.abu$tree
#' out <- PhdObs(data = data,datatype = "abundance",tree = tree,q = c(0,1,2),
#' reftime = seq(0.1,325,length.out = 40))
#' ggtqplot(out,profile = 'time')
#' }
#' @export
ggtqplot <- function(outcome,profile = 'q'){
  if(profile=='q'){
    Plotq(out = outcome)
  }else if (profile=='time'){
    Plott(out = outcome)
  }
}

#' @useDynLib PhD
#' @importFrom Rcpp sourceCpp
NULL
