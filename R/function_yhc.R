#===============PhDObs==================
datainf_yhc <- function(data, datatype, phylotr,reft){
  new <- phyBranchAL_Abu(phylotr,data,datatype,reft)
  new$treeNabu$branch.length <- new$BLbyT[,1]
  data <- data[data>0]
  PD_obs <- sum(new$treeNabu$branch.length)
  fg1 <- new$treeNabu[1:length(data),] %>% filter(branch.abun==1)
  fg2 <- new$treeNabu %>% filter(branch.abun==2)
  f1 <- nrow(fg1);g1 <- sum(fg1$branch.length)
  f2 <- nrow(fg2);g2 <- sum(fg2$branch.length)
  a1 <- c(sum(data),length(data),PD_obs,f1,f2,g1,g2)
  names(a1) <- c("n", "S.obs", "PD.obs", "f1*", "f2*", "g1", "g2")
  return(a1)
}

PD.qprofile_yhc=function(Abun, aL, q, reft, cal, datatype,forboot=FALSE, nforboot = NULL) {
  #input:
  #phylotr : a class phylo tree
  #Abun: species*community matrix dataframe,where row is species and colmun is community or site.
  # **NOTE: species names in Abun should be exactly the same as the species names in newick format.
  # reft: a value specifying reference time ,where the time interval [0,-t] of phylo-tree would be consider.
  n <- ifelse(is.null(nforboot), ifelse(datatype=='incidence_raw', ncol(Abun), sum(Abun)), nforboot)
  #n <- ifelse(datatype=='incidence_raw', ncol(Abun), sum(Abun) )
  t_bar <- ifelse(forboot==FALSE,reft,sum(aL[,1] * aL[,2]) / n)
  pAbun <- aL[,1]
  bL <- aL[,2]

  Sub = function(q) {
    PD=ifelse(q==1, exp(-sum(bL*pAbun/t_bar/n*log(pAbun/t_bar/n)) ),
              sum(bL*(pAbun/t_bar/n)^q)^(1/(1-q)))
    ifelse(cal=="PD",PD,PD/t_bar)
  }
  sapply(q, Sub)
}
Phdqtable_yhc <- function(datalist, phylotr, q, cal, datatype, nboot, conf, reft){
  # Note 200117: currently, the reference time is automatically fixed at tree height of pooled species.
  qtile <- qnorm(1-(1-conf)/2)
  if (datatype=="incidence_raw") {
    chaotree <- phylo2phytree(phylotr)
    H_max <- chaotree$treeH
    da <- lapply(datalist, rowSums) %>% do.call(cbind, .) %>% rowSums()
    aLTable <- phyExpandData(x = da,labels = names(da),phy = chaotree,datatype = datatype)[,c(4,3,7)]
    PD2 <- PD.qprofile_yhc(Abun = da, chaotree = chaotree, aL = aLTable,q = 2, reft = H_max, cal =  "PD", datatype)
    Q <- H_max-(H_max^2)/PD2
  }
  if (datatype=="abundance") {
    H_max <- get.rooted.tree.height(phylotr)
    da <- do.call(cbind, datalist) %>% rowSums()
    aL <- phyBranchAL_Abu(phylo = phylotr,data = da,datatype,refT = H_max)$treeNabu %>%
      select(branch.abun,branch.length)
    PD2 <- PD.qprofile_yhc(Abun = da, aL,q = 2, reft = H_max, cal =  "PD", datatype)
    Q <- H_max-(H_max^2)/PD2
  }
  nms <- names(datalist)
  if(cal=="PD") refts <- unique(sort(c(Q,H_max,reft))) else refts <- unique(sort(c(0.01,Q,H_max,reft)))
  out <- lapply(1:length(datalist), function(i){
    aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = refts)
    x_no0 <- datalist[[i]]
    x_no0 <- x_no0[x_no0>0]
    emp <- sapply(1:ncol(aL$BLbyT), function(j){
      aL_table <- matrix(c(unlist(aL$treeNabu[,ncol(aL$treeNabu)]),aL$BLbyT[,j]),ncol = 2)
      PD.qprofile_yhc(Abun = x_no0, aL = aL_table, q,reft = refts[j], cal, datatype)
    }) %>% c()
    if(nboot!=0){
      Boots <- Boots.one_yhc(data=x_no0,phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = reft)
      Li_b <- Boots$Li
      f0 <- Boots$f0
      ses <- sapply(1:nboot, function(B){
        x_b <- Boots$boot_data[,B]
        x_b_tmp <- x_b[1:(length(x_no0)+f0)]
        out_b <- sapply(1:ncol(aL$BLbyT), function(j){
          if(f0==0) Li_refb <- aL$BLbyT[,j] else Li_refb <- c(aL$BLbyT[1:length(x_no0),j],
                                                              Boots$Li[(length(x_no0)+1):(length(x_no0)+f0)],
                                                              aL$BLbyT[-(1:length(x_no0)),j])
          Li_refb[Li_refb>refts[j]] <- refts[j]
          aL_table_b <- matrix(c(x_b,Li_refb),ncol = 2)
          aL_table_b <- aL_table_b[aL_table_b[,1]>0,]
          PD.qprofile_yhc(Abun = x_b_tmp[x_b_tmp>0], aL = aL_table_b, q,reft = refts[j], cal, datatype,
                          forboot=TRUE, nforboot = sum(x_no0))
        }) %>% c()
        out_b
      }) %>% apply(., 1, sd)
    }else{
      ses <- rep(0,length(emp))
    }
    output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
    output[output[,2]<0,2] <- 0
    output
  }) %>% do.call(rbind,.)

  if (cal=="PD") {
    RefTime <- rep(rep(c(paste0("Q = ", round(Q, 2)), paste0("root = ", round(H_max, 2))), each=length(q)),length(nms))
  } else {
    RefTime <- rep(rep(c("Taxonomic",paste0("Q = ",round(Q, 2)), paste0("root = ",round(H_max, 2))),
                       each=length(q)),length(nms))
  }
  nms_tmp <- rep(nms,each = length(q)*length(refts))
  Outputforq <- tibble(q = rep(q,length(refts)*length(nms)), Empirical = out[,1],LCL = out[,2], UCL = out[,3],
                       RefTime = RefTime,Community = nms_tmp)
  Outputforq <- Outputforq %>% mutate (method = ifelse(cal=="PD", "Phylogenetic Diversity", "Phylogenetic Hill numbers"))
  return(Outputforq)
}
color_nogreen <- function(n) {
  all <- c("red", "blue", "orange", "purple", "pink", "cyan", "brown", "yellow")
  all[1:n]
}
PD.Tprofile_yhc=function(Abun, ai, Lis, q, times, cal, datatype,forboot=FALSE, nforboot = NULL) {
  #q can be a vector
  n <- ifelse(is.null(nforboot), ifelse(datatype=='incidence_raw', ncol(Abun), sum(Abun)), nforboot)
  if(forboot==TRUE) times <- ai %*% Lis/n
  bL <- sapply(1:length(times), function(i) Lis[,i]/times[i])
  pAbun <- ai/n

  out <- sapply(q, function(j){
    if(j==1) as.numeric(-(pAbun*log(pAbun)) %*% bL) %>% exp()
    else as.numeric((pAbun)^j %*% bL) %>% .^(1/(1-j))
  }) %>% as.matrix()
  if(cal=="PD"){
    out <- sapply(1:length(times), function(i){
      out[i,]*times[i]
    }) %>% t()
  }
  out
}
Phdttable_yhc <- function(datalist, phylotr, times, cal, datatype, nboot, conf){
  # Note 200117: currently, the reference time is automatically fixed at tree height of pooled species.
  qtile <- qnorm(1-(1-conf)/2)
  if (datatype=="incidence_raw") {
    chaotree <- phylo2phytree(phylotr)
    H_max <- chaotree$treeH
    da <- lapply(datalist, rowSums) %>% do.call(cbind, .) %>% rowSums()
    aLTable <- phyExpandData(x = da,labels = names(da),phy = chaotree,datatype = datatype)[,c(4,3,7)]
    PD2 <- PD.qprofile_yhc(Abun = da, chaotree = chaotree, aL = aLTable,q = 2, reft = H_max, cal =  "PD", datatype)
    Q <- H_max-(H_max^2)/PD2
  }
  if (datatype=="abundance") {
    H_max <- get.rooted.tree.height(phylotr)
    da <- do.call(cbind, datalist) %>% rowSums()
    aL <- phyBranchAL_Abu(phylo = phylotr,data = da,datatype,refT = H_max)$treeNabu %>%
      select(branch.abun,branch.length)
    PD2 <- PD.qprofile_yhc(Abun = da, aL,q = 2, reft = H_max, cal =  "PD", datatype)
    Q <- H_max-(H_max^2)/PD2
  }
  nms <- names(datalist)
  q_int <- c(0, 1, 2)
  out <- lapply(1:length(datalist), function(i){
    aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = times)
    x_no0 <- datalist[[i]]
    x_no0 <- x_no0[x_no0>0]
    emp <- PD.Tprofile_yhc(Abun = x_no0,ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT,q=q_int,times,
                    cal,datatype) %>% c()
    if(nboot!=0){
      Boots <- Boots.one_yhc(data=x_no0,phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = H_max)

      f0 <- Boots$f0
      if(f0==0){
        Lis_b <- aL$BLbyT
      }else{
        Lis_b <- sapply(1:ncol(aL$BLbyT), function(l){
          tmp <- c(aL$BLbyT[1:length(x_no0),l],Boots$Li[(length(x_no0)+1):(length(x_no0)+f0)],
                   aL$BLbyT[-(1:length(x_no0)),l])
          tmp[tmp>times[l]] <- times[l]
          tmp
        })
      }
      ses <- sapply(1:nboot, function(B){
        x_b <- Boots$boot_data[,B]
        x_b_tmp <- x_b[1:(length(x_no0)+f0)]
        out_b <- PD.Tprofile_yhc(Abun = x_b_tmp[x_b_tmp>0],ai = x_b[x_b>0],
                                 Lis=Lis_b[x_b>0,],q=q_int,times,cal,
                                 datatype,forboot = TRUE,nforboot = sum(x_no0)) %>% c()
        out_b
      }) %>% apply(., 1, sd)
    }else{
      ses <- rep(0,length(emp))
    }
    output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
    output[output[,2]<0,2] <- 0
    output
  }) %>% do.call(rbind,.)


  Outputfort <- tibble(time = rep(times,length(q_int)*length(datalist)),
                       Empirical = out[,1],LCL = out[,2], UCL = out[,3],
                       q = rep(rep(q_int, each=length(times)),length(datalist)),
                       Community = rep(nms, each=length(times)*length(q_int)),
                       method=ifelse(cal=="PD", "Phylogenetic Diversity", "Phylogenetic Hill numbers"))
  out = list(fort = Outputfort, Q_Height=c(Q, H_max))
  return(out)
}
AUC_one_table_yhc <- function(datalist, phylotr, knot, cal, datatype, nboot, conf, reft_max) {
  qtile <- qnorm(1-(1-conf)/2)
  times_AUC <- seq(0.01, reft_max, length.out = knot)
  nms <- names(datalist)
  q_int <- c(0, 1, 2)
  AUC <- lapply(1:length(datalist),function(i){
    aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = times_AUC)
    x_no0 <- datalist[[i]]
    x_no0 <- x_no0[x_no0>0]
    emp <- PD.Tprofile_yhc(Abun = x_no0,ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT,q=q_int,times = times_AUC,
                           cal,datatype)
    LA <-  c(diff(times_AUC),0) %*% emp %>% as.numeric()
    RA <-   c(0,diff(times_AUC)) %*% emp %>% as.numeric()
    auc <- colMeans(rbind(LA,RA))
    if(nboot!=0){
      Boots <- Boots.one_yhc(data=x_no0,phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = reft_max)

      f0 <- Boots$f0
      if(f0==0){
        Lis_b <- aL$BLbyT
      }else{
        Lis_b <- sapply(1:ncol(aL$BLbyT), function(l){
          tmp <- c(aL$BLbyT[1:length(x_no0),l],Boots$Li[(length(x_no0)+1):(length(x_no0)+f0)],
                   aL$BLbyT[-(1:length(x_no0)),l])
          tmp[tmp>times[l]] <- times[l]
          tmp
        })
      }
      ses <- sapply(1:nboot, function(B){
        x_b <- Boots$boot_data[,B]
        x_b_tmp <- x_b[1:(length(x_no0)+f0)]
        out_b <- PD.Tprofile_yhc(Abun = x_b_tmp[x_b_tmp>0],ai = x_b[x_b>0],
                                 Lis=Lis_b[x_b>0,],q=q_int,times,cal,
                                 datatype,forboot = TRUE,nforboot = sum(x_no0))
        LA_b <-  c(diff(times_AUC),0) %*% out_b %>% as.numeric()
        RA_b <-   c(0,diff(times_AUC)) %*% out_b %>% as.numeric()
        auc_b <- colMeans(rbind(LA_b,RA_b))
        auc_b
      }) %>% apply(., 1, sd)
    }else{
      ses <- rep(0,length(auc))
    }
    output <- cbind(auc,auc-qtile*ses,auc+qtile*ses)
    output[output[,2]<0,2] <- 0
    output
  }) %>% do.call(rbind,.)

  AUC <- tibble(q = rep(q_int,length(nms)), Empirical = AUC[,1],LCL = AUC[,2], UCL = AUC[,3],
                       Community = rep(nms,each = length(q_int)))
  AUC
}
Plott_yhc <- function(out, cal, Q_Height){
  fort <- out
  fort$q = paste0("q = ", fort$q)
  fort$q <- factor(fort$q)
  Community <- unique(fort$Community)
  Q = Q_Height[1]
  root = Q_Height[2]
  # print(fort)
  # print(Q)
  if (cal=="PD") {
    if(length(Community)==1){
      p2 <- ggplot(fort, aes(x=time, y=Empirical)) + geom_line(size=1.5,aes(color=q))+
        geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=q),linetype = 0,alpha=0.3)
      p2 <-  p2 +xlab("time")+ylab("Phylogenetic Diversity")+theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
        geom_point(size=3, data=subset(fort, time%in%c(Q, root)), aes(x=time, y=Empirical, color=q))+
        annotate('text',x=Q, y=0.1,label="Q" ,parse = TRUE,size=5, color = "gray") +
        annotate('text',x=root, y=0.1, label="root",parse = TRUE,size=5, color = "gray") +
        geom_vline(xintercept = c(Q,root), linetype = "longdash",size=0.5, color = "gray")
    }else{
      p2 <- ggplot(fort, aes(x=time, y=Empirical, color=Community, linetype=Community)) + geom_line(size=1.5)  +
        geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=Community),linetype = 0,alpha=0.3)+
        scale_color_manual(values = color_nogreen(length(unique(fort$Community))))+
        scale_fill_manual(values = color_nogreen(length(unique(fort$Community))))+
        theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
        geom_point(size=3, data=subset(fort, time%in%c(Q, root)), aes(x=time, y=Empirical, color=Community))+
        annotate('text',x=Q, y=0.1,label="Q" ,parse = TRUE,size=5, color = "gray") +
        annotate('text',x=root, y=0.1, label="root",parse = TRUE,size=5, color = "gray") +
        geom_vline(xintercept = c(Q,root), linetype = "longdash",size=0.5, color = "gray") +
        facet_wrap(~q, scales = "free")
      p2 <-  p2 +xlab("time")+ylab("Phylogenetic Diversity")
    }
  } else {
    if(length(Community)==1){
      p2 <- ggplot(fort, aes(x=time, y=Empirical)) + geom_line(size=1.5,aes(color=q))+
        geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=q),linetype = 0,alpha=0.3)
      p2 <-  p2 +xlab("time")+ylab("Phylogenetic Hill numbers")+theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
        geom_point(size=3, data=subset(fort, time%in%c(0.01, Q, root)), aes(x=time, y=Empirical, color=q))+
        annotate('text',x=Q, y=0.1,label="Q" ,parse = TRUE,size=5, color = "gray") +
        annotate('text',x=0.01, y=0.1,label=0.01 ,parse = TRUE,size=5, color = "gray") +
        annotate('text',x=root, y=0.1, label="root",parse = TRUE,size=5, color = "gray") +
        geom_vline(xintercept = c(0.01, Q,root), linetype = "longdash",size=0.5, color = "gray")
    }else{
      p2 <- ggplot(fort, aes(x=time, y=Empirical, color=Community, linetype=Community)) + geom_line(size=1.5)  +
        geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=Community),linetype = 0,alpha=0.3)+
        scale_color_manual(values = color_nogreen(length(unique(fort$Community))))+
        scale_fill_manual(values = color_nogreen(length(unique(fort$Community))))+
        theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
        geom_point(size=3, data=subset(fort, time%in%c(0.01, Q, root)), aes(x=time, y=Empirical, color=Community))+
        annotate('text',x=Q, y=0.1,label="Q" ,parse = TRUE,size=5, color = "gray") +
        annotate('text',x=0.01, y=0.1,label=0.01 ,parse = TRUE,size=5, color = "gray") +
        annotate('text',x=root, y=0.1, label="root",parse = TRUE,size=5, color = "gray") +
        geom_vline(xintercept = c(0.01, Q,root), linetype = "longdash",size=0.5, color = "gray") +
        facet_wrap(~q, scales = "free")
      p2 <-  p2 +xlab("time")+ylab("Phylogenetic Hill numbers")
    }
  }
  return(p2)
}
Plotq_yhc <- function(out, cal){
  forq <- out
  forq$RefTime <- as.character(forq$RefTime)
  forq$RefTime <- factor(forq$RefTime, levels = unique(forq$RefTime <- factor(forq$RefTime)))
  Community <- unique(forq$Community)
  q1 <- unique(forq$q[(forq$q %% 1)==0])
  if (cal=="PD") {
    haha=sapply(strsplit(as.character(forq$RefTime), " = "), function(i) i) %>% as.vector() %>% unique()
    Q=haha[2] %>% as.numeric()
    root=haha[4] %>% as.numeric()
    if(length(Community)==1){
      p1 <- ggplot(forq, aes(x=q, y=Empirical, color=RefTime)) + geom_line(size=1.5)+
        geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=RefTime),linetype = 0,alpha=0.3)
      #lai 1006
      p1 <-  p1 +xlab("Order q")+ylab("Phylogenetic Diversity") +theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
        geom_point(size=3, data=subset(forq, q%in%q1), aes(x=q, y=Empirical, color=RefTime))
    }else{
      p1 <- ggplot(forq, aes(x=q, y=Empirical, color=Community, linetype=Community)) + geom_line(size=1.5)  +
        geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=Community),linetype = 0,alpha=0.3)+
        scale_color_manual(values = color_nogreen(length(unique(forq$Community))))+
        scale_fill_manual(values = color_nogreen(length(unique(forq$Community))))+
        theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
        geom_point(size=3, data=subset(forq, q%in%q1), aes(x=q, y=Empirical, color=Community))+
        facet_wrap(~RefTime, scales = "free")
      p1 <-  p1 +xlab("Order q")+ylab("Phylogenetic Diversity")
    }
  } else {
    haha=sapply(strsplit(as.character(forq$RefTime), " = "), function(i) i) %>% as.vector() %>% unique()
    Q=haha[[2]][2] %>% as.numeric()
    root=haha[[3]][2] %>% as.numeric()
    if(length(Community)==1){
      p1 <- ggplot(forq, aes(x=q, y=Empirical, color=RefTime)) + geom_line(size=1.5)+
        geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=RefTime),linetype = 0,alpha=0.3)
      p1 <-  p1 +xlab("Order q")+ylab("Phylogenetic Hill numbers") +theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
        geom_point(size=3, data=subset(forq, q%in%q1), aes(x=q, y=Empirical, color=RefTime))
    }else{
      p1 <- ggplot(forq, aes(x=q, y=Empirical, color=Community, linetype=Community)) + geom_line(size=1.5)  +
        geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=Community),linetype = 0,alpha=0.3)+
        scale_color_manual(values = color_nogreen(length(unique(forq$Community))))+
        scale_fill_manual(values = color_nogreen(length(unique(forq$Community))))+
        theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
        geom_point(size=3, data=subset(forq, q%in%q1), aes(x=q, y=Empirical, color=Community))+
        facet_wrap(~RefTime, scales = "free")
      p1 <-  p1 +xlab("Order q")+ylab("Phylogenetic Hill numbers")
    }
  }
  return(p1)
}
#===============PhDAsy==================
AsyPD_yhc <- function(datalist, datatype, phylotr, q, nboot, conf, reft){#change final list name
  nms <- names(datalist)
  qtile <- qnorm(1-(1-conf)/2)
  if(nboot!=0){
    Estoutput <- lapply(datalist,function(x){
      #atime <- Sys.time()
      aL <- phyBranchAL_Abu(phylo = phylotr,data = x,datatype,refT = reft)
      aL$treeNabu$branch.length <- aL$BLbyT[,1]
      aL <- aL$treeNabu %>% select(branch.abun,branch.length)
      x_no0 <- x[x>0]
      #btime <- Sys.time()
      #print(paste0("aL sample: ",btime-atime))
      #atime <- Sys.time()
      out <- PhD.q.est_yhc(Abun = x_no0,aL,q,reft)
      #btime <- Sys.time()
      #print(paste0("Est sample: ",btime-atime))
      est <- out[[1]]
      # atime <- Sys.time()
      Boots <- Boots.one_yhc(data=x_no0,phylo = phylotr,aL,datatype,nboot,reft = reft)
      Li_b <- Boots$Li
      f0 <- Boots$f0
      # btime <- Sys.time()
      # print(paste0("Creat boot sample: ",btime-atime))
      #atime_b_all <- Sys.time()
      ses <- sapply(1:nboot, function(B){
        #atime <- Sys.time()
        x_b <- Boots$boot_data[,B]
        x_b_tmp <- x_b[1:(length(x_no0)+f0)]
        aL_b <- tibble(branch.abun = x_b[x_b>0], branch.length=Li_b[x_b>0])
        outb <- PhD.q.est_yhc(Abun = x_b_tmp[x_b_tmp>0],aL = aL_b,q,reft,forboot = T,nforboot = sum(x))[[1]]
        #btime <- Sys.time()
        #print(paste0("Est boot sample",B,": ",btime-atime))
        return(outb)
      }) %>% apply(., 1, sd)
      #btime_b_all <- Sys.time()
      #print(paste0("Bootstrap computation time: ",btime_b_all-atime_b_all))
      est <- tibble(Order = q, Estimate = est, se = ses, LCL = est - qtile*ses, UCL = est + qtile*ses)
      info <- out[[2]]
      list(est = est,info =info)
    })
    info <- sapply(Estoutput, function(x) x[[2]])
    Estoutput <- Estoutput %>% lapply(., function(x) x[[1]]) %>%
      do.call(rbind,.) %>% mutate(Community = rep(names(datalist),each = length(q)))
    Estoutput$LCL[Estoutput$LCL<0] = 0
  }else{
    Estoutput <- lapply(datalist,function(x){
      # a <- Sys.time()
      aL <- phy_BranchAL_Abu(phylo = phylotr,data = x,datatype,refT = reft)
      aL$treeNabu$branch.length <- aL$BLbyT[[1]]
      aL <- aL$treeNabu %>% select(branch.abun,branch.length)
      x_no0 <- x[x>0]
      # b <- Sys.time()
      # print(b-a)
      # a <- Sys.time()
      out <- PhD.q.est_yhc(Abun = x_no0,aL,q,reft)
      # b<-Sys.time()
      # print(b-a)
      out
    })
    info <- sapply(Estoutput, function(x) x[[2]])
    Estoutput <- sapply(Estoutput, function(x) x[[1]]) %>% as.numeric() %>%
      tibble(Order = rep(q,length(datalist)), Estimate = ., se = 0, LCL = ., UCL = ., Community = rep(nms,each = length(q)))
  }
  return(list(Estoutput = Estoutput, info = info))
}
Asy_plot_yhc = function(output, type, method=NULL){##add title
  if(is.null(method) == T){
    ylab_ <- "Phylogenetic diversity"
  }else{
    if( substr(method,1,1) != "1" ){
      ylab_ <- paste("Phylogenetic", method ,"diversity")
    }else{
      ylab_ <- method
    }
  }
  if(is.null(method) == T & type == 1) {
    title_ <- "Phylogenetic diversity_estimated"
  } else if(is.null(method) == T & type == 2) {
    title_ <- "Phylogenetic diversity_empirical"
  } else if( substr(method,1,1) != "1"  & type == 1) {
    title_ <- paste("Phylogenetic", method ,"diversity", "_estimated")
  } else if( substr(method,1,1) != "1"  & type == 2) {
    title_ <- paste("Phylogenetic", method ,"diversity", "_empirical")
  } else if( substr(method,1,1) == "1"  & type == 1) {
    title_ <- paste(method, "_estimated")
  } else if( substr(method,1,1) == "1"  & type == 2) {
    title_ <- paste(method, "_empirical")
  } else if(type == 3) {
    title_ <- paste(method, "_estimated")
  } else {
    title_ <- paste(method, "_empirical")
  }
  if(ncol(output) %in% c(2, 5)) output = cbind(output, "Beta")
  if(ncol(output) == 6){
    colnames(output) = c("x", "y", "se", "LCL", "UCL","site")
    site <- unique(output[,6]) %>% unlist
  }else{
    colnames(output) = c("x", "y","site")
    site <- unique(output[,3]) %>% unlist
  }
  q <- unlist(output$x)
  q1<-q[(round(q) - q) ==0]
  if(length(site) == 1){
    p <- ggplot(output,aes(x=x,y=y))+geom_line(size=1.5,color="#F8766D")+xlab("Order q")+
      geom_point(size=3, data=subset(output, x%in%q1),color="#F8766D")+
      ylab(ylab_)+theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))
    if(ncol(output) == 6) p <- p + geom_ribbon(aes(ymin=LCL,ymax=UCL),alpha=0.3,fill="#F8766D")
  }else{
    p <- ggplot(output,aes(x=x,y=y,color=site,linetype=site))+geom_line(size=1.5)+xlab("Order q")+
      scale_color_manual(values = color_nogreen(length(unique(output$site))))+
      geom_point(size=3, data=subset(output, x%in%q1))+
      ylab(ylab_)+theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))
    if(ncol(output) == 6) p <- p + geom_ribbon(aes(ymin=LCL,ymax=UCL,fill=site),alpha=0.3, colour=NA)+scale_fill_manual(values = color_nogreen(length(unique(output$site))))
  }
  p <- p+ggtitle(title_)
  return(p)
}
PhD.q.est_yhc = function(Abun, aL, q, reft, forboot=FALSE, nforboot = NULL){ # reft is a single value
  n <- ifelse(is.null(nforboot),  sum(Abun), nforboot)
  PD_obs <- sum(aL$branch.length)
  fg1 <- aL[1:length(Abun),] %>% filter(branch.abun==1)
  fg2 <- aL %>% filter(branch.abun==2)
  f1 <- nrow(fg1);g1 <- sum(fg1$branch.length)
  f2 <- nrow(fg2);g2 <- sum(fg2$branch.length)
  if(f2 > 0){
    A = 2*f2/((n-1)*f1+2*f2)
  }else if(f2 == 0 & f1 > 0){
    A = 2/((n-1)*(f1-1)+2)
  }else{
    A = 1
  }
  if(forboot == TRUE){
    t_bar <- sum(aL[,1]/n*aL[,2])
  }else{
    t_bar <- reft
  }
  #print(paste0("t_bar:",t_bar))
  #t_bar <- sum(aL[,1]*aL[,2]/n)
  tmpaL <- aL %>% group_by(branch.abun, branch.length) %>% summarise(n_node = n()) %>% as.matrix()
  deltas <- sapply(0:(n-1), function(k){
    del_tmp <- tmpaL[tmpaL[,1]<=(n-k),,drop=FALSE]
    delta(del_tmpaL = del_tmp,k,n)
  })

  Sub <- function(q){
    if(q==0){
      ans <- PD_obs+Dq0(n,f1,f2,g1,g2,A)
    }else if(q==1){
      h2 <- Dq1_1(n,g1,A)
      h1 <- aL %>% filter(branch.abun<=(n-1)) %>% mutate(diga = digamma(n)-digamma(branch.abun)) %>%
        apply(., 1, prod) %>% sum(.)/n
      #print(paste0("A:",A," g1:",g1," h1:", h1, " h2:",h2))
      h <- h1+h2
      ans <- t_bar*exp(h/t_bar)
    }else if(q==2){
      ans <- Dq2(as.matrix(tmpaL),n,t_bar)
    }else{
      # timea <- Sys.time()
      k <- 0:(n-1)
      a <- (choose(q-1,k)*(-1)^k*deltas) %>% sum
      b <- ifelse(g1==0|A==1,0,(g1*((1-A)^(1-n))/n)*(A^(q-1)-sum(choose(q-1,k)*(A-1)^k)))
      ans <- ((a+b)/(t_bar^q))^(1/(1-q))
      # timeb <- Sys.time()
      # print(timeb-timea)
    }
    return(ans)
  }
  est <- sapply(q, Sub)
  info <- c("n" = n, "S.obs" = length(Abun), "PD.obs" = PD_obs, "f1*" = f1,
            "f2*" = f2, "g1" = g1, "g2" = g2)
  list(est = est, info = info)
}
Boots.one_yhc = function(data, phylo, aL, datatype, nboot,reft){
  if(datatype=='abundance'){
    n <- sum(data)
    f1 <- sum(data==1)
    f2 <- sum(data==2)
    f0 <- ceiling( ifelse( f2>0 , ((n-1) / n) * (((f1)^2) / (2*f2) ), ((n-1) / n) * (f1)*(f1-1) / 2 ) )
    c <- ifelse(f2>0, 1 - (f1/n)*((n-1)*f1/((n-1)*f1+2*f2)),
                1 - (f1/n)*((n-1)*(f1-1)/((n-1)*(f1-1)+2)))
    lambda <- (1-c) / sum((data/n)*(1- (data/n) )^n)
    p_hat <- (data/n) * (1-lambda*(1- (data/n) )^n)
    p_hat0 <- rep( (1-c) / f0 , f0 );names(p_hat0) <- paste0("notob",1:length(p_hat0))
    g1 <- aL$branch.length[aL$branch.abun==1] %>% sum
    g2 <- aL$branch.length[aL$branch.abun==2] %>% sum
    g0_hat <- ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
    #l_hat <- g0_hat/f0
    pL_b <- phyBranchAL_Abu(phylo, p_hat, datatype,refT = reft)
    pL_b$treeNabu$branch.length <- pL_b$BLbyT[,1]
    pL_b <- pL_b$treeNabu
    pL_b[length(p_hat)+1,"branch.abun"] <- 1
    Li <- c(unlist(pL_b[1:length(data),"branch.length"]), rep(g0_hat/f0,f0),
            unlist(pL_b[-(1:length(data)),"branch.length"]))
    p_hat <- c(p_hat,p_hat0,unlist(pL_b[-(1:length(data)),"branch.abun"]))
    # print(c(f0,c,g0_hat/f0))
    #set.seed(i)
    boot_data <- sapply(p_hat,function(p) rbinom(n = nboot,size = n,prob = p)) %>% t()

  }else if(datatype=='incidence_raw'){
    # n <- ncol(data)
    # y <- rowSums(data)
    # u <- sum(y)
    # f1 <- sum(y==1)
    # f2 <- sum(y==2)
    # f0 <- ceiling( ifelse( f2>0 , ((n-1) / n) * (((f1)^2) / (2*f2) ), ((n-1) / n) * (f1)*(f1-1) / 2 ) )
    # c <- ifelse(f2>0, 1 - (f1/u)*((n-1)*f1/((n-1)*f1+2*f2)),
    #             1 - (f1/u)*((n-1)*(f1-1)/((n-1)*(f1-1)+2)))
    # lambda <- (1-c) / sum((y/u)*(1- (y/n) )^n)
    # p_hat <- (y/n) * (1-lambda*(1- (y/n) )^n)
    # p_hat0 <- rep( (1-c) / f0 , f0 );names(p_hat0) <- paste0("notob",1:length(p_hat0))
    # phat <- c(p_hat, p_hat0)
    # g1 <- aL$branch.length[aL$branch.abun==1] %>% sum
    # g2 <- aL$branch.length[aL$branch.abun==2] %>% sum
    # g0_hat <- ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
    # l_hat <- g0_hat/f0
    # boot_data <- rmultinom(n = nboot,size = n,prob = phat)

  }
  return(list(boot_data=boot_data,Li = Li, f0 = f0))
}
#===============inextPD==================
iNEXTPD_yhc = function(datalist, phylotr, datatype, Q, nboot, conf=0.95, size=NULL, knots=40, endpoint, reft){
  nms <- names(datalist)
  qtile <- qnorm(1-(1-conf)/2)
  if(datatype!='incidence_raw') n <- unlist(lapply(datalist, sum))
  if(datatype=='incidence_raw') n <- unlist(lapply(datalist, ncol))
  m <- lapply(1:length(n),function(i){
    if(is.null(size) == T){
      if(endpoint <= n[i]){
        mi <- floor(seq(1,endpoint,length.out = knots))
      }else{
        mi <- floor(c(seq(1,n[i],length.out = floor(knots/2)),seq(n[i]+1,endpoint,length.out = knots-floor(knots/2))))
      }
      if(2*n[i] < knots & 2*n[i] > endpoint) m[[i]] <- 1:endpoint
    }else{
      mi <- size
      if( sum(size==n[i]) == 0 ) mi <- sort(c(n[i],mi))
    }
    unique(mi)
  })
  if(nboot!=0){
    Estoutput <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      aL$treeNabu$branch.length <- aL$BLbyT[,1]
      aL <- aL$treeNabu %>% select(branch.abun,branch.length)
      x_no0 <- datalist[[i]]
      x_no0 <- x_no0[x_no0>0]
      qPDm = sapply(Q,function(Qq){
        PhD.m.est_yhc(Abun = x_no0,aL = aL,m = m[[i]],Q = Qq,reft = reft, datatype)
      }) %>% as.numeric() #%>% tibble(order=rep(Q,each=length(m[[i]])),qPD=.)
      covm = Coverage_yhc(x_no0, datatype, m[[i]])
      Boots <- Boots.one_yhc(data=x_no0,phylo = phylotr,aL,datatype,nboot,reft = reft)
      Li_b <- Boots$Li
      f0 <- Boots$f0
      ses <- sapply(1:nboot, function(B){
        # atime <- Sys.time()
        x_b <- Boots$boot_data[,B]
        x_b_tmp <- x_b[1:(length(x_no0)+f0)]
        aL_b <- tibble(branch.abun = x_b[x_b>0], branch.length=Li_b[x_b>0])
        qPDm_b <- sapply(Q,function(Qq){
          PhD.m.est_yhc(Abun = x_b_tmp[x_b_tmp>0],aL=aL_b,m=m[[i]],Q=Qq,reft,datatype,forboot = T,nforboot = sum(datalist[[i]]))
        }) %>% as.numeric()
        covm_b = Coverage_yhc(x_b_tmp[x_b_tmp>0], datatype, m[[i]])
        # btime <- Sys.time()
        # print(paste0("Est boot sample",B,": ",btime-atime))
        return(c(qPDm_b,covm_b))
      }) %>% apply(., 1, sd)
      method <- ifelse(m[[i]]>n[i],'Extrapolation',ifelse(m[[i]]<n[i],'Rarefaction','Observed'))
      orderq <- rep(Q,each=length(m[[i]]))
      ses_cov <- ses[(length(ses)-length(m[[i]])+1):length(ses)]
      ses_pd <- ses[-((length(ses)-length(m[[i]])+1):length(ses))]
      tibble(m=rep(m[[i]],length(Q)),method=rep(method,length(Q)),order=orderq,
             qPD=qPDm,qPD.LCL=qPDm-qtile*ses_pd,qPD.UCL=qPDm+qtile*ses_pd,
             SC=rep(covm,length(Q)),SC.LCL=rep(covm-qtile*ses_cov,length(Q)),SC.UCL=rep(covm+qtile*ses_cov,length(Q)),
             site = nms[i])
    }) %>% do.call(rbind, .)
    if(datatype=='incidence_raw') colnames(Estoutput)[1] <-  "t"
  }else{
    Estoutput <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      aL$treeNabu$branch.length <- aL$BLbyT[,1]
      aL <- aL$treeNabu %>% select(branch.abun,branch.length)
      x_no0 <- datalist[[i]]
      x_no0 <- x_no0[x_no0>0]
      qPDm = sapply(Q,function(Qq){
        PhD.m.est_yhc(Abun = x_no0,aL = aL,m = m[[i]],Q = Qq,reft = reft, datatype, forboot=TRUE)
      }) %>% as.numeric() #%>% tibble(order=rep(Q,each=length(m[[i]])),qPD=.)
      covm = Coverage_yhc(x_no0, datatype, m[[i]])
      # Boots <- Boots.one_yhc(data=x_no0,phylo = phylotr,aL,datatype,nboot,reft = reft)
      # Li_b <- Boots$Li
      # f0 <- Boots$f0
      ses <- rep(0,length(qPDm)+length(covm))
      method <- ifelse(m[[i]]>n[i],'Extrapolation',ifelse(m[[i]]<n[i],'Rarefaction','Observed'))
      orderq <- rep(Q,each=length(m[[i]]))
      ses_cov <- ses[(length(ses)-length(m[[i]])+1):length(ses)]
      ses_pd <- ses[-((length(ses)-length(m[[i]])+1):length(ses))]
      tibble(m=rep(m[[i]],length(Q)),method=rep(method,length(Q)),order=orderq,
             qPD=qPDm,qPD.LCL=qPDm-qtile*ses_pd,qPD.UCL=qPDm+qtile*ses_pd,
             SC=rep(covm,length(Q)),SC.LCL=rep(covm-qtile*ses_cov,length(Q)),SC.UCL=rep(covm+qtile*ses_cov,length(Q)),
             site = nms[i])
    }) %>% do.call(rbind, .)
    if(datatype=='incidence_raw') colnames(Estoutput)[1] <-  "t"
  }
  return(Estoutput)
}
PhD.m.est_yhc = function(Abun, aL, m, Q, reft, datatype, forboot=FALSE, nforboot = NULL){
  n <- ifelse(is.null(nforboot), ifelse(datatype=='incidence_raw', ncol(Abun), sum(Abun)), nforboot)
  #n <- ifelse(datatype=='incidence_raw', ncol(Abun), sum(Abun) )
  t_bar <- ifelse(forboot==FALSE,reft,sum(aL[,1] * aL[,2]) / n)
  aL_matrix = as.matrix(aL)
  RPD_m = RPD(aL_matrix, n, n-1, Q)
  # print(RPD_m)
  obs = RPD(aL_matrix, n, n, Q)
  #asymptotic value
  asy <- PhD.q.est_yhc(Abun,aL,Q,reft,forboot=forboot,nforboot = n)$est
  #beta
  if(asy == obs) beta = 0
  if(asy != obs) beta =(obs-RPD_m)/(asy-RPD_m)
  #Extrapolation
  EPD = function(m,Q){
    m = m-n
    if( Q == 0 | Q == 1 ) EPD = obs+(asy-obs)*(1-(1-beta)^m)
    if( Q == 2 ) EPD = 1/sum( (aL_matrix[,2]/(t_bar)^2)*((1/(n+m))*(aL_matrix[,1]/n)+((n+m-1)/(n+m))*(aL_matrix[,1]*(aL_matrix[,1]-1)/(n*(n-1)))) )
    return(EPD)
  }
  Sub = function(m){
    if(m<n){
      RPD(aL_matrix,n,m,Q)
    }else if(m==n){
      obs
    }else{
      EPD(m,Q)
    }
  }
  sapply(m, Sub)
}
Coverage_yhc = function(data, datatype, m){
  n <- ifelse(datatype=='incidence_raw', ncol(data), sum(data) )
  if(datatype == "incidence_raw") datatype = "incidence"
  ifelse(class(data)=="matrix" || class(data)=="data.frame", type <- "raw", type <- "numeric")
  ifelse(type == "raw", x <- rowSums(data), x <- data )
  if(type=="raw" || datatype=='incidence') u<-sum(x)
  f1 = sum(x == 1)
  f2 = sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  Sub <- function(m){
    #if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
    if(m < n) {
      xx <- x[(n-x)>=m]
      out <- 1-sum(xx / n * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
    }
    if(m == n) out <- 1-f1/n*A
    if(m > n) out <- 1-f1/n*A^(m-n+1)
    out
  }
  Sub2 <- function(m){
    #if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
    if(m < n) {
      xx <- x[(n-x)>=m]
      out <- 1-sum(xx / u * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
    }
    if(m == n) out <- 1-f1/u*A
    if(m > n) out <- 1-f1/u*A^(m-n+1)
    out
  }
  sapply(m, FUN = function(i){
    ifelse(datatype=='incidence', Sub2(i), Sub(i) )
  })
}
RE_plot_yhc = function(data, datatype, type, method=NULL){
  data <- as.data.frame(data)
  if(is.null(method) == T){
    ylab_ <- "Phylogenetic diversity"
  }else{
    if( substr(method,1,1) == "q" ){
      ylab_ <- paste("Phylogenetic",substr(method,2,nchar(method)),"diversity")
    }else{
      ylab_ <- method
    }
  }
  title <- c("Sample-size-based sampling curve", "Coverage-based sampling curve", "Sample completeness curve")
  title <- title[type]
  if(ncol(data) %in% c(9, 10)){
    x <- ifelse(datatype=='incidence_raw', 'plots', "individuals")
    if(ncol(data) == 9) data = cbind(data, "Beta")
    site <- unique(data[,10])
    if(type == 1){
      output <- data[, c(1, 2, 4, 5, 6, 10,3)]
      xlab_ <- paste0("Number of ", x)
    }else if(type == 2){
      output <- data[, c(7, 2, 4, 5, 6, 10,3)]
      xlab_ <- "Sample Coverage"
    }else if(type == 3){
      output <- data[, c(1, 2, 7, 8, 9, 10,3)]
      xlab_ <- paste0("Number of ", x)
      ylab_ <- "Sample Coverage"
    }
    colnames(output) <- c("x", "Method", "y", "LCL", "UCL", "site","order")
    output[, 2] <- as.character(output[, 2])
    output[, 6] <- as.character(output[, 6])
    output2 <- output
    output2[output2[,2]=="Observed",2] <- "Rarefaction"
    output3 <- output
    output3[output3[,2]=="Observed",2] <- "Extrapolation"
    output2=rbind(output2, output3)
    if(length(site) == 1){
      ggplot(output2, aes(x = x, y = y))+
        geom_ribbon(aes(ymin = LCL, ymax = UCL),fill="#F8766D",alpha=0.3)+geom_line(size=1.5, aes(x = x, y = y, linetype=Method),color="#F8766D")+
        geom_point(size=3, data=subset(output, Method=="Observed"),color="#F8766D")+xlab(xlab_)+ylab(ylab_)+
        scale_linetype_manual(values = c("dashed", "solid"), name="Method",breaks=c("Interpolated", "Extrapolated"), labels=c("Interpolation", "Extrapolation"))+
        theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
        ggtitle(title)+facet_wrap(~order,scales = "free_y")+
        guides(linetype=guide_legend(keywidth=2.5))
    }else{
      ggplot(output2, aes(x = x, y = y, color=site))+geom_line(size=1.5, aes(x = x, y = y, color=site, linetype=Method))+
        scale_color_manual(values = color_nogreen(length(unique(output2$site))))+
        geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = site), alpha=0.3, colour=NA)+
        scale_fill_manual(values = color_nogreen(length(unique(output2$site))))+
        geom_point(size=3, data=subset(output, Method=="Observed"))+xlab(xlab_)+ylab(ylab_)+
        scale_linetype_manual(values = c("dashed", "solid"), name="Method",breaks=c("Interpolated", "Extrapolated"), labels=c("Interpolation", "Extrapolation"))+
        theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
        ggtitle(title)+facet_wrap(~order,scales = "free_y")+
        guides(linetype=guide_legend(keywidth=2.5))
    }
  }else if(ncol(data) %in% c(5, 6)){
    x <- ifelse(datatype=='incidence_raw', 'plots', "individuals")
    if(ncol(data) == 5) data = cbind(data, "Beta")
    site <- unique(data[,6])
    if(type == 1){
      output <- data[, c(1, 2, 4, 6)]
      xlab_ <- paste0("Number of ", x)
    }else if(type == 2){
      output <- data[, c(5, 2, 4, 6)]
      xlab_ <- "Sample Coverage"
    }else if(type == 3){
      output <- data[, c(1, 2, 5, 6)]
      xlab_ <- paste0("Number of ", x)
      ylab_ <- "Sample Coverage"
    }
    colnames(output) <- c("x", "Method", "y", "site")
    output[, 2] <- as.character(output[, 2])
    output2 <- output
    output2[output2[,2]=="Observed",2] <- "Interpolated"
    output3 <- output
    output3[output3[,2]=="Observed",2] <- "Extrapolated"
    output2=rbind(output2, output3)
    if(length(site) == 1){
      ggplot(output2, aes(x = x, y = y))+geom_line(size=1.5, aes(x = x, y = y, linetype=Method),color="#F8766D")+
        geom_point(size=3, data=subset(output, Method=="Observed"),color="#F8766D")+xlab(xlab_)+ylab(ylab_)+
        scale_linetype_manual(values = c("dashed", "solid"), name="Method",breaks=c("Interpolated", "Extrapolated"), labels=c("Interpolation", "Extrapolation"))+
        theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
        ggtitle(title)+
        guides(linetype=guide_legend(keywidth=2.5))
    }else{
      ggplot(output2, aes(x = x, y = y, color=site))+geom_line(size=1.5, aes(x = x, y = y, color=site, linetype=Method))+
        scale_color_manual(values = color_nogreen(length(unique(output2$site))))+
        geom_point(size=3, data=subset(output, Method=="Observed"))+xlab(xlab_)+ylab(ylab_)+
        scale_linetype_manual(values = c("dashed", "solid"), name="Method",breaks=c("Interpolated", "Extrapolated"), labels=c("Interpolation", "Extrapolation"))+
        theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
        ggtitle(title)+
        guides(linetype=guide_legend(keywidth=2.5))
    }
  }
}
