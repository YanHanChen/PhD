#===============PhDObs==================
datainf_yhc <- function(data, datatype, phylotr,reft){
  if(datatype == "abundance"){
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
  }else if(datatype == 'incidence_raw'){
    new <- phyBranchAL_Inc(phylotr,data,datatype,reft)
    new$treeNabu$branch.length <- new$BLbyT[,1]
    data <- data[rowSums(data)>0,]
    PD_obs <- sum(new$treeNabu$branch.length)
    Qr1 <- new$treeNabu %>% filter(branch.abun==1)
    Qr2 <- new$treeNabu %>% filter(branch.abun==2)
    Q1 <- nrow(Qr1);r1 <- sum(Qr1$branch.length)
    Q2 <- nrow(Qr2);r2 <- sum(Qr2$branch.length)
    a1 <- c(ncol(data),nrow(data),PD_obs,Q1,Q2,r1,r2)
    names(a1) <- c("T", "S.obs", "PD.obs", "Q1*", "Q2*", "R1", "R2")
  }
  return(a1)

}

PD.qprofile_yhc=function(aL, q, reft, splunits = NULL, cal ,datatype, nforboot = NULL) {
  #aL is a table of 3 columns: abundance, branch lengths and characters specifying the group of each node.
  pAbun <- unlist(aL[,1])
  bL <- unlist(aL[,2])
  if(datatype=="abundance"){
    Abun <- unlist(aL$branch.abun[aL$tgroup=="Tip"])
    n <- ifelse(is.null(nforboot),sum(Abun), nforboot)
    t_bar <- ifelse(is.null(nforboot),reft, sum(pAbun / n * bL)) #just to save computation time
  }else if (datatype=="incidence_raw"){
    n <- splunits
    t_bar <- sum(pAbun / n * bL)
  }

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
  # all assemblages.
  if (datatype=="incidence_raw") {
    H_max <- get.rooted.tree.height(phylotr)
    da <- lapply(datalist, rowSums) %>% do.call(cbind, .) %>% rowSums()
  }
  if (datatype=="abundance") {
    H_max <- get.rooted.tree.height(phylotr)
    da <- do.call(cbind, datalist) %>% rowSums()
  }
  if(abs(H_max - reft)<=1e-4) H_max <- reft
  # Note 200204: currently, to compute Q, we treat incidence data as abundance and absolutely pool
  aL <- phyBranchAL_Abu(phylo = phylotr,data = da,"abundance",refT = H_max)$treeNabu %>%
    select(branch.abun,branch.length,tgroup)
  PD2 <- PD.qprofile_yhc(aL,q = 2, reft = H_max, cal =  "PD",datatype =  "abundance")
  Q <- H_max-(H_max^2)/PD2
  nms <- names(datalist)
  if(cal=="PD") refts <- unique(sort(c(Q,H_max,reft))) else refts <- unique(sort(c(1,Q,H_max,reft)))
  if(datatype=="abundance"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = refts)
      x <- datalist[[i]] %>% .[.>0]
      n <- sum(x)
      emp <- sapply(1:ncol(aL$BLbyT), function(j){
        aL_table <- tibble(branch.abun = unlist(aL$treeNabu[,"branch.abun"]),branch.length = aL$BLbyT[,j],
                           tgroup = unlist(aL$treeNabu[,"tgroup"]))
        PD.qprofile_yhc(aL = aL_table, q = q,reft = refts[j], cal = cal, datatype = datatype)
      }) %>% c()
      if(nboot!=0){
        Boots <- Boots.one_yhc(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = refts, BLs = aL$BLbyT )
        Li_b <- Boots$Li
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",length(x)+f0),rep("Inode",nrow(Li_b)-length(x)-f0))
        aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        ses <- sapply(1:nboot, function(B){
          aL_table_b[,1] <- Boots$boot_data[,B]
          isn0 <- aL_table_b[,1]>0
          Li_b_tmp <- Li_b[isn0,]
          aL_table_b <- aL_table_b[isn0,]
          out_b <- sapply(1:ncol(aL$BLbyT), function(j){
            Li_refb <- Li_b_tmp[,j]
            Li_refb[Li_refb>refts[j]] <- refts[j]
            aL_table_b[,2] <- Li_refb
            PD.qprofile_yhc(aL = aL_table_b, q = q,reft = refts[j],cal = cal, datatype = datatype, nforboot = n)
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
  }else if(datatype=="incidence_raw"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Inc(phylo = phylotr,data = datalist[[i]],datatype,refT = refts)
      x <- datalist[[i]] %>% .[rowSums(.)>0,colSums(.)>0]
      n <- ncol(x)
      # For incidence type, we use occurence frequencies instead of raw data since we already have aL table.
      emp <- sapply(1:ncol(aL$BLbyT), function(j){
        aL_table <- tibble(branch.abun = unlist(aL$treeNabu[,"branch.abun"]),branch.length = aL$BLbyT[,j],
                           tgroup = unlist(aL$treeNabu[,"tgroup"]))
        PD.qprofile_yhc(aL = aL_table, q,reft = refts[j], cal, datatype,splunits=n)
      }) %>% c()
      if(nboot!=0){
        Boots <- Boots.one_yhc(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = refts, BLs = aL$BLbyT,splunits = n)
        Li_b <- Boots$Li
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",nrow(x)+f0),rep("Inode",nrow(Li_b)-nrow(x)-f0))
        aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        ses <- sapply(1:nboot, function(B){
          aL_table_b[,1] <- Boots$boot_data[,B]
          isn0 <- aL_table_b[,1]>0
          Li_b_tmp <- Li_b[isn0,]
          aL_table_b <- aL_table_b[isn0,]
          out_b <- sapply(1:ncol(aL$BLbyT), function(j){
            Li_refb <- Li_b_tmp[,j]
            Li_refb[Li_refb>refts[j]] <- refts[j]
            aL_table_b[,2] <- Li_refb
            PD.qprofile_yhc(aL = aL_table_b, q,reft = refts[j], splunits = n, cal, datatype)
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
  }
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
PD.Tprofile_yhc=function(ai,Lis, tgroup, q, times,splunits = NULL, cal, datatype,nforboot = NULL) {
  #q can be a vector
  if(datatype=="abundance"){
    Abun <- ai[tgroup=="Tip"]
    if(is.null(nforboot)){
      n <- sum(Abun)
      t_bars <- times
    }else{
      n <- nforboot
      t_bars <- as.numeric(ai %*% Lis/n)
    }
  }else if (datatype=="incidence_raw"){
    n <- splunits
    t_bars <- as.numeric(ai %*% Lis/n)
  }

  bL <- sapply(1:length(t_bars), function(i) Lis[,i]/t_bars[i])
  pAbun <- ai/n

  out <- sapply(q, function(j){
    if(j==1) as.numeric(-(pAbun*log(pAbun)) %*% bL) %>% exp()
    else as.numeric((pAbun)^j %*% bL) %>% .^(1/(1-j))
  }) %>% as.matrix()
  if(cal=="PD"){
    out <- sapply(1:length(t_bars), function(i){
      out[i,]*t_bars[i]
    }) %>% t()
  }
  out
}
Phdttable_yhc <- function(datalist, phylotr, times, cal, datatype, nboot, conf){
  # Note 200117: currently, the reference time is automatically fixed at tree height of pooled species.
  qtile <- qnorm(1-(1-conf)/2)
  # all assemblages.
  if (datatype=="incidence_raw") {
    H_max <- get.rooted.tree.height(phylotr)
    da <- lapply(datalist, rowSums) %>% do.call(cbind, .) %>% rowSums()
  }
  if (datatype=="abundance") {
    H_max <- get.rooted.tree.height(phylotr)
    da <- do.call(cbind, datalist) %>% rowSums()
  }
  if(abs(H_max - reft)<=1e-4) H_max <- reft
  # Note 200204: currently, to compute Q, we treat incidence data as abundance and absolutely pool
  aL <- phyBranchAL_Abu(phylo = phylotr,data = da,"abundance",refT = H_max)$treeNabu %>%
    select(branch.abun,branch.length,tgroup)
  PD2 <- PD.qprofile_yhc(aL,q = 2, reft = H_max, cal =  "PD",datatype =  "abundance")
  Q <- H_max-(H_max^2)/PD2
  nms <- names(datalist)

  q_int <- c(0, 1, 2)
  if(datatype=="abundance"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = times)
      x <- datalist[[i]] %>% .[.>0]
      n <- sum(x)
      emp <- PD.Tprofile_yhc(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT, tgroup = aL$treeNabu$tgroup,
                             q=q_int,times = times,cal = cal,datatype = datatype) %>% c()

      if(nboot!=0){
        Boots <- Boots.one_yhc(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = times, BLs = aL$BLbyT )
        Lis_b <- Boots$Li
        Lis_b <- sapply(1:length(times),function(l){
          tmp <- Lis_b[,l]
          tmp[tmp>times[l]] <- times[l]
          tmp
        })
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",length(x)+f0),rep("Inode",nrow(Lis_b)-length(x)-f0))
        ses <- sapply(1:nboot, function(B){
          x_b <- Boots$boot_data[,B]
          isn0 <- x_b>0
          Lis_b_tmp <- Lis_b[isn0,]
          tgroup_B_tmp <- tgroup_B[isn0]
          x_b <- x_b[isn0]
          out_b <- PD.Tprofile_yhc(ai = x_b,Lis=Lis_b_tmp,tgroup = tgroup_B_tmp,q=q_int,times = times,cal = cal,datatype =
                                     datatype,nforboot = n) %>% c()
          out_b
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(0,length(emp))
      }
      output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }else if(datatype=="incidence_raw"){
    out <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Inc(phylo = phylotr,data = datalist[[i]],datatype,refT = times)
      x <- datalist[[i]] %>% .[rowSums(.)>0,]
      n <- ncol(x)
      emp <- PD.Tprofile_yhc(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT, tgroup = aL$treeNabu$tgroup,
                             q=q_int,times = times,splunits = n,cal = cal,datatype = datatype) %>% c()

      if(nboot!=0){
        Boots <- Boots.one_yhc(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = times,
                               BLs = aL$BLbyT,splunits = n)
        Lis_b <- Boots$Li
        Lis_b <- sapply(1:length(times),function(l){
          tmp <- Lis_b[,l]
          tmp[tmp>times[l]] <- times[l]
          tmp
        })
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",nrow(x)+f0),rep("Inode",nrow(Lis_b)-nrow(x)-f0))
        ses <- sapply(1:nboot, function(B){
          x_b <- Boots$boot_data[,B]
          isn0 <- x_b>0
          Lis_b_tmp <- Lis_b[isn0,]
          tgroup_B_tmp <- tgroup_B[isn0]
          x_b <- x_b[isn0]
          out_b <- PD.Tprofile_yhc(ai = x_b,Lis=Lis_b_tmp,tgroup = tgroup_B_tmp,q=q_int,times = times,
                                   splunits = n,cal = cal,datatype = datatype) %>% c()
          out_b
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(0,length(emp))
      }
      output <- cbind(emp,emp-qtile*ses,emp+qtile*ses)
      output[output[,2]<0,2] <- 0
      output
    }) %>% do.call(rbind,.)
  }

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
  if(datatype=="abundance"){
    AUC <- lapply(1:length(datalist),function(i){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = times_AUC)
      x <- datalist[[i]] %>% .[.>0]
      n <- sum(x)
      emp <- PD.Tprofile_yhc(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT, tgroup = aL$treeNabu$tgroup,
                             q=q_int,times = times_AUC,cal = cal,datatype = datatype)
      # print(paste0("emp:",dim(emp)," AUC diff:",length(c(diff(times_AUC),0))))
      LA <-  c(diff(times_AUC),0) %*% emp %>% as.numeric()
      RA <-   c(0,diff(times_AUC)) %*% emp %>% as.numeric()
      auc <- colMeans(rbind(LA,RA))
      if(nboot!=0){
        Boots <- Boots.one_yhc(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = times_AUC, BLs = aL$BLbyT )
        Lis_b <- Boots$Li
        Lis_b <- sapply(1:length(times_AUC),function(l){
          tmp <- Lis_b[,l]
          tmp[tmp>times_AUC[l]] <- times_AUC[l]
          tmp
        })
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",length(x)+f0),rep("Inode",nrow(Lis_b)-length(x)-f0))

        ses <- sapply(1:nboot, function(B){
          x_b <- Boots$boot_data[,B]
          isn0 <- x_b>0
          Lis_b_tmp <- Lis_b[isn0,]
          tgroup_B_tmp <- tgroup_B[isn0]
          x_b <- x_b[isn0]
          out_b <- PD.Tprofile_yhc(ai = x_b,Lis=Lis_b_tmp,tgroup = tgroup_B_tmp,q=q_int,times = times_AUC,cal = cal,datatype =
                                     datatype,nforboot = n)
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
  }else if (datatype=="incidence_raw"){
    AUC <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Inc(phylo = phylotr,data = datalist[[i]],datatype,refT = times_AUC)
      x <- datalist[[i]] %>% .[rowSums(.)>0,]
      n <- ncol(x)
      emp <- PD.Tprofile_yhc(ai = aL$treeNabu$branch.abun,Lis=aL$BLbyT, tgroup = aL$treeNabu$tgroup,
                             q=q_int,times = times_AUC,splunits = n,cal = cal,datatype = datatype)
      LA <-  c(diff(times_AUC),0) %*% emp %>% as.numeric()
      RA <-   c(0,diff(times_AUC)) %*% emp %>% as.numeric()
      auc <- colMeans(rbind(LA,RA))
      if(nboot!=0){
        Boots <- Boots.one_yhc(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = times_AUC,
                               BLs = aL$BLbyT,splunits = n)
        Lis_b <- Boots$Li
        Lis_b <- sapply(1:length(times_AUC),function(l){
          tmp <- Lis_b[,l]
          tmp[tmp>times_AUC[l]] <- times_AUC[l]
          tmp
        })
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",nrow(x)+f0),rep("Inode",nrow(Lis_b)-nrow(x)-f0))
        ses <- sapply(1:nboot, function(B){
          x_b <- Boots$boot_data[,B]
          isn0 <- x_b>0
          Lis_b_tmp <- Lis_b[isn0,]
          tgroup_B_tmp <- tgroup_B[isn0]
          x_b <- x_b[isn0]
          out_b <- PD.Tprofile_yhc(ai = x_b,Lis=Lis_b_tmp,tgroup = tgroup_B_tmp,q=q_int,times = times_AUC,
                                   splunits = n,cal = cal,datatype = datatype)
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
  }


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
  if(datatype=="abundance"){
    Estoutput <- lapply(datalist,function(x){
      #atime <- Sys.time()
      x <- x[x>0]
      n <- sum(x)
      aL <- phyBranchAL_Abu(phylo = phylotr,data = x,datatype,refT = reft)
      aL$treeNabu$branch.length <- aL$BLbyT[,1]
      aL_table <- aL$treeNabu %>% select(branch.abun,branch.length,tgroup)
      out <- PhD.q.est_yhc(aL = aL_table,q = q,reft = reft,datatype = datatype)
      est <- out[[1]]
      if(nboot!=0){
        Boots <- Boots.one_yhc(phylo = phylotr,aL = aL$treeNabu,datatype,nboot,reft = reft, BLs = aL$BLbyT )
        Li_b <- Boots$Li
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",length(x)+f0),rep("Inode",nrow(Li_b)-length(x)-f0))
        aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        ses <- sapply(1:nboot, function(B){
          aL_table_b[,1] <- Boots$boot_data[,B]
          isn0 <- aL_table_b[,1]>0
          Li_b_tmp <- Li_b[isn0,]
          aL_table_b <- aL_table_b[isn0,]
          outb <- PhD.q.est_yhc(aL = aL_table_b,q = q,reft = reft,datatype = datatype,nforboot = n)[[1]]
          return(outb)
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(0,length(est))
      }
      est <- tibble(Order = q, Estimate = est, se = ses, LCL = est - qtile*ses, UCL = est + qtile*ses)
      info <- out[[2]]
      list(est = est,info =info)
    })
    info <- sapply(Estoutput, function(x) x[[2]])
  }else if(datatype=="incidence_raw"){
    Estoutput <- lapply(datalist,function(x){
      #atime <- Sys.time()
      x <- x[rowSums(x)>0,colSums(x)>0]
      n <- ncol(x)
      aL <- phyBranchAL_Inc(phylo = phylotr,data = x,datatype,refT = reft)
      aL$treeNabu$branch.length <- aL$BLbyT[,1]
      aL_table <- aL$treeNabu %>% select(branch.abun,branch.length,tgroup)
      out <- PhD.q.est_yhc(aL = aL_table,q = q,reft = reft,splunits = n,datatype = datatype)
      est <- out[[1]]
      if(nboot!=0){
        Boots <- Boots.one_yhc(phylo = phylotr,aL = aL$treeNabu,datatype = datatype,nboot = nboot,
                               splunits = n,reft = reft, BLs = aL$BLbyT )
        Li_b <- Boots$Li
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",nrow(x)+f0),rep("Inode",nrow(Li_b)-nrow(x)-f0))
        aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        ses <- sapply(1:nboot, function(B){
          aL_table_b[,1] <- Boots$boot_data[,B]
          isn0 <- aL_table_b[,1]>0
          Li_b_tmp <- Li_b[isn0,]
          aL_table_b <- aL_table_b[isn0,]
          outb <- PhD.q.est_yhc(aL = aL_table_b,q = q,reft = reft,splunits = n,
                                datatype = datatype,nforboot = n)[[1]]
          return(outb)
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(0,length(est))
      }
      est <- tibble(Order = q, Estimate = est, se = ses, LCL = est - qtile*ses, UCL = est + qtile*ses)
      info <- out[[2]]
      list(est = est,info =info)
    })
  }
  info <- sapply(Estoutput, function(x) x[[2]])
  Estoutput <- Estoutput %>% lapply(., function(x) x[[1]]) %>%
    do.call(rbind,.) %>% mutate(Community = rep(names(datalist),each = length(q)))
  Estoutput$LCL[Estoutput$LCL<0] = 0
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
PhD.q.est_yhc = function(aL, q, reft, splunits=NULL, datatype, nforboot = NULL){ # reft is a single value
  if(datatype=="abundance"){
    Abun <- unlist(aL$branch.abun[aL$tgroup=="Tip"])
    n <- ifelse(is.null(nforboot),sum(Abun), nforboot)
    t_bar <- ifelse(is.null(nforboot),reft, sum(aL[,1] / n * aL[,2])) #just to save computation time
  }else if (datatype=="incidence_raw"){
    n <- splunits
    inci_freq <- unlist(aL$branch.abun[aL$tgroup=="Tip"])
    t_bar <- sum(aL[,1] / n * aL[,2])
  }
  aL <- aL %>% select(branch.abun, branch.length)
  PD_obs <- sum(aL$branch.length)
  fg1 <- aL %>% filter(branch.abun==1)
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
      h1 <- aL %>% filter(branch.abun<=(n-1)) %>%
        mutate(diga = digamma(n)-digamma(branch.abun)) %>% apply(., 1, prod) %>% sum(.)/n
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
  if(datatype=="abundance")info <- c("n" = n, "S.obs" = length(Abun), "PD.obs" = PD_obs, "f1*" = f1,
                                     "f2*" = f2, "g1" = g1, "g2" = g2)
  else if (datatype == "incidence_raw") info <- c("T" = n, "S.obs" = length(inci_freq), "PD.obs" = PD_obs, "Q1*" = f1,
                                                  "Q2*" = f2, "R1" = g1, "R2" = g2)

  list(est = est, info = info)
}
Boots.one_yhc = function(phylo, aL, datatype, nboot,reft, BLs, splunits = NULL){
  if(datatype=='abundance'){
    data <- unlist(aL$branch.abun[aL$tgroup=="Tip"])
    names(data) <- rownames(BLs)[1:length(data)]
    n <- sum(data)
    f1 <- sum(data==1)
    f2 <- sum(data==2)
    f0 <- ceiling( ifelse( f2>0 , ((n-1) / n) * (((f1)^2) / (2*f2) ), ((n-1) / n) * (f1)*(f1-1) / 2 ) )
    c <- ifelse(f2>0, 1 - (f1/n)*((n-1)*f1/((n-1)*f1+2*f2)),
                1 - (f1/n)*((n-1)*(f1-1)/((n-1)*(f1-1)+2)))
    lambda <- (1-c) / sum((data/n)*(1- (data/n) )^n)
    ####Notice that the species order of data is different from that of aL####
    p_hat <- (data/n) * (1-lambda*(1- (data/n) )^n)
    p_hat0 <- rep( (1-c) / f0 , f0 );names(p_hat0) <- paste0("notob",1:length(p_hat0))
    g1 <- aL$branch.length[aL$branch.abun==1] %>% sum
    g2 <- aL$branch.length[aL$branch.abun==2] %>% sum
    g0_hat <- ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
    ###Notice that the species order of pL_b doesn't change even that of data changes. (property of phyBranchAL_Abu)
    pL_b <- phyBranchAL_Abu(phylo, p_hat, datatype,refT = reft[1])
    pL_b$treeNabu$branch.length <- pL_b$BLbyT[,1]
    pL_b <- pL_b$treeNabu
    pL_b[length(p_hat)+1,"branch.abun"] <- 1
    Li <- BLs %>% as.matrix()
    Li <- rbind(Li[1:length(data),,drop=F],matrix(g0_hat/f0,nrow = f0,ncol = ncol(Li)), Li[-(1:length(data)),,drop=F])
    p_hat <- c(p_hat,p_hat0,unlist(pL_b[-(1:length(data)),"branch.abun"]))
    boot_data <- sapply(p_hat,function(p) rbinom(n = nboot,size = n,prob = p)) %>% t()
  }else if(datatype=='incidence_raw'){
    n <- splunits
    data <- unlist(aL$branch.abun[aL$tgroup=="Tip"])
    u <- sum(data)
    f1 <- sum(data==1)
    f2 <- sum(data==2)
    f0 <- ceiling( ifelse( f2>0 , ((n-1) / n) * (((f1)^2) / (2*f2) ), ((n-1) / n) * (f1)*(f1-1) / 2 ) )
    c <- ifelse(f2>0, 1 - (f1/u)*((n-1)*f1/((n-1)*f1+2*f2)),
                1 - (f1/u)*((n-1)*(f1-1)/((n-1)*(f1-1)+2)))
    lambda <- u/n*(1-c) / sum((data/n)*(1- (data/n) )^n)
    p_hat0 <- rep( (u/n) * (1-c) / f0 , f0 );names(p_hat0) <- paste0("notob",1:length(p_hat0))
    g1 <- aL$branch.length[aL$branch.abun==1] %>% sum
    g2 <- aL$branch.length[aL$branch.abun==2] %>% sum
    g0_hat <- ifelse( g2>((g1*f2)/(2*f1)) , ((n-1)/n)*(g1^2/(2*g2)) , ((n-1)/n)*(g1*(f1-1)/(2*(f2+1))) )
    data_iB <- unlist(aL$branch.abun)
    pL_b <- (data_iB/n) * (1-lambda*(1- (data_iB/n) )^n)
    Li <- BLs %>% as.matrix()
    Li <- rbind(Li[1:length(data),,drop=F],matrix(g0_hat/f0,nrow = f0,ncol = ncol(Li)),
                Li[-(1:length(data)),,drop=F])
    p_hat <- c(pL_b[1:length(data)],p_hat0,pL_b[-(1:length(data))])
    boot_data <- sapply(p_hat,function(p) rbinom(n = nboot,size = n,prob = p)) %>% t()
  }
  return(list(boot_data=boot_data,Li = Li, f0 = f0))
}
#===============inextPD==================
iNEXTPD_yhc = function(datalist, phylotr, datatype, Q, nboot, conf=0.95, size=NULL, knots=40, endpoint, reft){
  nms <- names(datalist)
  qtile <- qnorm(1-(1-conf)/2)
  if(datatype=="abundance") ns <- unlist(lapply(datalist, sum))
  if(datatype=='incidence_raw') ns <- unlist(lapply(datalist, ncol))
  m <- lapply(1:length(ns),function(i){
    if(is.null(size) == T){
      if(endpoint <= ns[i]){
        mi <- floor(seq(1,endpoint,length.out = knots))
      }else{
        mi <- floor(c(seq(1,ns[i],length.out = floor(knots/2)),seq(ns[i]+1,endpoint,length.out = knots-floor(knots/2))))
      }
      if(2*ns[i] < knots & 2*ns[i] > endpoint) m[[i]] <- 1:endpoint
    }else{
      mi <- size
      if( sum(size==ns[i]) == 0 ) mi <- sort(c(ns[i],mi))
    }
    unique(mi)
  })
  if(datatype=="abundance"){
    Estoutput <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Abu(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      aL$treeNabu$branch.length <- aL$BLbyT[,1]
      aL_table <- aL$treeNabu %>% select(branch.abun,branch.length,tgroup)
      x <- datalist[[i]] %>% .[.>0]
      n <- sum(x)
      qPDm <- PhD.m.est_yhc(aL = aL_table,m = m[[i]],Q = Q,reft = reft, datatype = datatype,nforboot = NULL) %>%
        t() %>% as.numeric()
      covm = Coverage_yhc(x, datatype, m[[i]])
      if(nboot!=0){
        Boots <- Boots.one_yhc(phylo = phylotr,aL$treeNabu,datatype,nboot,reft,aL$BLbyT)
        Li_b <- Boots$Li
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",length(x)+f0),rep("Inode",nrow(Li_b)-length(x)-f0))
        aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        ses <- sapply(1:nboot, function(B){
          # atime <- Sys.time()
          aL_table_b[,1] <- Boots$boot_data[,B]
          isn0 <- aL_table_b[,1]>0
          Li_b_tmp <- Li_b[isn0,]
          aL_table_b <- aL_table_b[isn0,]
          qPDm_b <-  PhD.m.est_yhc(aL=aL_table_b,m=m[[i]],Q=Q,reft,datatype,nforboot = n) %>%
            t() %>% as.numeric()
          covm_b = Coverage_yhc(unlist(aL_table_b$branch.abun[aL_table_b$tgroup=="Tip"]), datatype, m[[i]])
          # btime <- Sys.time()
          # print(paste0("Est boot sample",B,": ",btime-atime))
          return(c(qPDm_b,covm_b))
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(0,length(c(qPDm,covm)))
      }
      method <- ifelse(m[[i]]>n,'Extrapolation',ifelse(m[[i]]<n,'Rarefaction','Observed'))
      orderq <- rep(Q,each=length(m[[i]]))
      ses_cov <- ses[(length(ses)-length(m[[i]])+1):length(ses)]
      ses_pd <- ses[-((length(ses)-length(m[[i]])+1):length(ses))]
      tibble(m=rep(m[[i]],length(Q)),method=rep(method,length(Q)),order=orderq,
             qPD=qPDm,qPD.LCL=qPDm-qtile*ses_pd,qPD.UCL=qPDm+qtile*ses_pd,
             SC=rep(covm,length(Q)),SC.LCL=rep(covm-qtile*ses_cov,length(Q)),SC.UCL=rep(covm+qtile*ses_cov,length(Q)),
             site = nms[i])
    }) %>% do.call(rbind, .)
  }else if(datatype=="incidence_raw"){
    Estoutput <- lapply(1:length(datalist), function(i){
      aL <- phyBranchAL_Inc(phylo = phylotr,data = datalist[[i]],datatype,refT = reft)
      aL$treeNabu$branch.length <- aL$BLbyT[,1]
      aL_table <- aL$treeNabu %>% select(branch.abun,branch.length,tgroup)
      x <- datalist[[i]] %>% .[rowSums(.)>0,colSums(.)>0]
      n <- ncol(x)
      qPDm <- PhD.m.est_yhc(aL = aL_table,m = m[[i]],Q = Q,reft = reft, datatype = datatype,nforboot = NULL,splunits = n) %>%
        t() %>% as.numeric()
      covm = Coverage_yhc(x, datatype, m[[i]])
      if(nboot!=0){
        Boots <- Boots.one_yhc(phylo = phylotr,aL$treeNabu,datatype,nboot,reft,aL$BLbyT,n)
        Li_b <- Boots$Li
        f0 <- Boots$f0
        tgroup_B <- c(rep("Tip",nrow(x)+f0),rep("Inode",nrow(Li_b)-nrow(x)-f0))
        aL_table_b <- tibble(branch.abun = 0, branch.length= Li_b[,1],tgroup = tgroup_B)
        ses <- sapply(1:nboot, function(B){
          # atime <- Sys.time()
          aL_table_b[,1] <- Boots$boot_data[,B]
          isn0 <- aL_table_b[,1]>0
          Li_b_tmp <- Li_b[isn0,]
          aL_table_b <- aL_table_b[isn0,]
          qPDm_b <-  PhD.m.est_yhc(aL=aL_table_b,m=m[[i]],Q=Q,reft,datatype,splunits = n) %>%
            t() %>% as.numeric()
          covm_b = Coverage_yhc(unlist(aL_table_b$branch.abun[aL_table_b$tgroup=="Tip"]), datatype, m[[i]],n)
          # btime <- Sys.time()
          # print(paste0("Est boot sample",B,": ",btime-atime))
          return(c(qPDm_b,covm_b))
        }) %>% apply(., 1, sd)
      }else{
        ses <- rep(0,length(c(qPDm,covm)))
      }
      method <- ifelse(m[[i]]>n,'Extrapolation',ifelse(m[[i]]<n,'Rarefaction','Observed'))
      orderq <- rep(Q,each=length(m[[i]]))
      ses_cov <- ses[(length(ses)-length(m[[i]])+1):length(ses)]
      ses_pd <- ses[-((length(ses)-length(m[[i]])+1):length(ses))]
      tibble(t=rep(m[[i]],length(Q)),method=rep(method,length(Q)),order=orderq,
             qPD=qPDm,qPD.LCL=qPDm-qtile*ses_pd,qPD.UCL=qPDm+qtile*ses_pd,
             SC=rep(covm,length(Q)),SC.LCL=rep(covm-qtile*ses_cov,length(Q)),SC.UCL=rep(covm+qtile*ses_cov,length(Q)),
             site = nms[i])
    }) %>% do.call(rbind, .)
  }
  return(Estoutput)
}
PhD.m.est_yhc = function(aL, m, Q, reft,datatype, nforboot = NULL, splunits = NULL){
  if(datatype=="abundance"){
    Abun <- unlist(aL$branch.abun[aL$tgroup=="Tip"])
    n <- ifelse(is.null(nforboot),sum(Abun), nforboot)
    t_bar <- ifelse(is.null(nforboot),reft, sum(aL[,1] / n * aL[,2])) #just to save computation time
  }else if (datatype=="incidence_raw"){
    n <- splunits
    inci_freq <- unlist(aL$branch.abun[aL$tgroup=="Tip"])
    t_bar <- sum(aL[,1] / n * aL[,2])
  }
  aL_matrix = as.matrix(aL[,c(1,2)])
  RPD_m = RPD(aL_matrix, n, n-1, Q)
  obs <- PD.qprofile_yhc(aL = aL, q = Q, reft=reft, cal="PD" ,datatype = datatype , nforboot = nforboot, splunits = splunits)
  #obs = RPD(aL_matrix, n, n, Q)
  #asymptotic value
  asy <- PhD.q.est_yhc(aL = aL,q = Q,reft = reft,datatype = datatype,nforboot = nforboot,splunits = splunits)$est
  #beta
  beta <- rep(0,length(Q))
  beta0plus <- which(asy != obs)
  beta[beta0plus] <- (obs[beta0plus]-RPD_m[beta0plus])/(asy[beta0plus]-RPD_m[beta0plus])
  # if(asy == obs) beta = 0
  # if(asy != obs) beta =(obs-RPD_m)/(asy-RPD_m)
  #Extrapolation
  EPD = function(m,Q){
    m = m-n
    out <- sapply(1:length(Q), function(i){
      if( Q[i] == 0 | Q[i] == 1 ) {
        obs[i]+(asy[i]-obs[i])*(1-(1-beta[i])^m)
      }else if( Q[i] == 2 ){
        1/sum( (aL_matrix[,2]/(t_bar)^2)*((1/(n+m))*(aL_matrix[,1]/n)+((n+m-1)/(n+m))*(aL_matrix[,1]*(aL_matrix[,1]-1)/(n*(n-1)))) )
      }
    })
    # if( Q == 0 | Q == 1 ) EPD = obs+(asy-obs)*(1-(1-beta)^m)
    # if( Q == 2 ) EPD = 1/sum( (aL_matrix[,2]/(t_bar)^2)*((1/(n+m))*(aL_matrix[,1]/n)+((n+m-1)/(n+m))*(aL_matrix[,1]*(aL_matrix[,1]-1)/(n*(n-1)))) )
    return(out)
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
Coverage_yhc = function(data, datatype, m, splunits = NULL){
  n <- ifelse(datatype=='incidence_raw', ifelse(is.null(splunits),ncol(data),splunits), sum(data) )
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
      output <- data[, c(1, 2, 7, 8, 9, 10,3)] %>% filter(order==1)
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
      outp <- ggplot(output2, aes(x = x, y = y))+
        geom_ribbon(aes(ymin = LCL, ymax = UCL),fill="#F8766D",alpha=0.3)+geom_line(size=1.5, aes(x = x, y = y, linetype=Method),color="#F8766D")+
        geom_point(size=3, data=subset(output, Method=="Observed"),color="#F8766D")+xlab(xlab_)+ylab(ylab_)+
        scale_linetype_manual(values = c("dashed", "solid"), name="Method",breaks=c("Interpolated", "Extrapolated"), labels=c("Interpolation", "Extrapolation"))+
        theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
        ggtitle(title)+guides(linetype=guide_legend(keywidth=2.5))
      if(type!=3) outp <- outp + facet_wrap(~order,scales = "free_y")
    }else{
      outp <- ggplot(output2, aes(x = x, y = y, color=site))+geom_line(size=1.5, aes(x = x, y = y, color=site, linetype=Method))+
        scale_color_manual(values = color_nogreen(length(unique(output2$site))))+
        geom_ribbon(aes(ymin = LCL, ymax = UCL, fill = site), alpha=0.3, colour=NA)+
        scale_fill_manual(values = color_nogreen(length(unique(output2$site))))+
        geom_point(size=3, data=subset(output, Method=="Observed"))+xlab(xlab_)+ylab(ylab_)+
        scale_linetype_manual(values = c("dashed", "solid"), name="Method",breaks=c("Interpolated", "Extrapolated"), labels=c("Interpolation", "Extrapolation"))+
        theme(text=element_text(size=20),legend.position="bottom",legend.key.width = unit(2,"cm"))+
        ggtitle(title)+guides(linetype=guide_legend(keywidth=2.5))
      if(type!=3) outp <- outp + facet_wrap(~order,scales = "free_y")
    }
  }
  return(outp)
}
