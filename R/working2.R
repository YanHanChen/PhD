atime <- Sys.time()
ans_old <- AsyPD_yhc1_old(datalist = mydata, datatype, phylotr = mytree, q, 100, conf, reft)
btime <- Sys.time()
btime - atime
atime <- Sys.time()
ans_new <- AsyPD_yhc1_new(datalist = mydata, datatype, phylotr = mytree, q, 5, conf, reft)
btime <- Sys.time()
btime - atime

library(ade4)
library(ape)
library(phytools)
library(dplyr)
library(ggplot2)
library(chaoUtility)
data("data.abu")
test_tr <- data.abu$tree
test_data <- data.abu$data
S <- 500
set.seed(123)
# tr <- rtree(n = S,rooted = T)
test_tr <- force.ultrametric(rtree(n = S,rooted = T), method=c("extend"))
#test_tr <- chronoMPL(tr)
test_data <- matrix(0,nrow = S,ncol = 3)
set.seed(123)
test_data[,1] <- rmultinom(n = 1,size = S*4,prob = rexp(S))
set.seed(456)
test_data[,2] <- rmultinom(n = 1,size = S*4,prob = rlnorm(n = S))
set.seed(789)
test_data[,3] <- rmultinom(n = 1,size = S*4,prob = rlnorm(n = S))
apply(test_data, 2, function(x) sum(x>0))
sum(rowSums(test_data)>0)
rownames(test_data) <- test_tr$tip.label
colnames(test_data) <- paste0("Site",1:ncol(test_data))

system.time(
  test <- PhdObs_yhc(data = test_data,tree = test_tr,datatype = "abundance",type = "PD",
                     profile = "time",q = seq(0, 2, by = 0.25))
)
system.time(
  test <- PhdAsy_yhc(data = test_data, tree = test_tr, datatype = "abundance",
                    q = seq(0, 2, by = 0.25), conf = 0.95, nboot = 0, reftime = NULL)
)

a<-Sys.time()
test <- PhdAsy_yhc1(data = test_data, tree = test_tr, datatype = "abundance",
                    q = seq(0, 2, by = 0.25), conf = 0.95, nboot = 3, reftime = NULL)
b<-Sys.time()
b-a

a<-Sys.time()
infos <- sapply(mydata, function(x){
  datainf_yhc1(data = x, datatype, phylotr = mytree,reft = reft)
})
b<-Sys.time()
b-a

# Test C.I accuracy
S <- 100
set.seed(123)
test_tr <- force.ultrametric(rtree(n = S,rooted = T), method=c("extend"))
set.seed(123)
test_data <- rmultinom(n = 1,size = S*4,prob = runif(S))
test_data <- rmultinom(n = 1,size = S*4,prob = rexp(S))
rownames(test_data) <- test_tr$tip.label
colnames(test_data) <- paste0("Site",1:ncol(test_data))


q <- c(seq(0,2,0.1))
tmp <- sapply(1:30, function(i){
  print(i)
  set.seed(i*2)
  test_data <- rmultinom(n = 1,size = S*4,prob = runif(S))
  rownames(test_data) <- test_tr$tip.label
  colnames(test_data) <- paste0("Site",1:ncol(test_data))
  use_mul <- PhdAsy_yhc(data = test_data,tree = test_tr,datatype = "abundance",q = q,nboot = 30)
  c(use_mul$asy$se,use_mul$asy$Estimate)
})
saveRDS(tmp,"tmp.rds")
tmp[c(1:length(q)),] %>% rowMeans()
tmp[-c(1:length(q)),] %>% apply(.,1,function(x) sd(x))

comp <- data.frame(q = q,se_true = apply(tmp[-c(1:length(q)),], 1, sd),se_mult = rowMeans(tmp[c(1:length(q)),]),
       se_bino = rowMeans(use_bin[c(1:length(q)),])) %>% melt(.,id.vars = "q",value.name = "se",
                                                       variable.name = "se_type") %>%
  as_tibble()
ggplot(comp)+geom_line(aes(x = q, y = se,col = se_type),size = 1.5) + theme_bw()+
  scale_color_manual(values=c(2,3,4))+
  theme(legend.position = "bottom",axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


se_true = apply(tmp[-c(1:11),], 1, sd)

bias <- data.frame(q = q, mult = rowMeans(apply(tmp[c(1:11),],2,function(x) x-se_true)),
                   bino = rowMeans(apply(use_bin[c(1:11),],2,function(x) x-se_true)))
RMSE_mult <-  sqrt(rowSums(apply(tmp[c(1:11),],2,function(x) (x-se_true)^2/50)))
RMSE_bino <-  sqrt(rowSums(apply(use_bin[c(1:11),],2,function(x) (x-se_true)^2/50)))
RMSE <- data.frame(q = q,mult = RMSE_mult,bino = RMSE_bino)
bias <- bias %>% melt(.,id.vars = "q",value.name = "Value",
                      variable.name = "method") %>% mutate(Error_type = "bias")
RMSE <- RMSE %>% melt(.,id.vars = "q",value.name = "Value",
                      variable.name = "method") %>% mutate(Error_type = "RMSE")
err <- rbind(bias,RMSE)

ggplot(err)+facet_wrap(~Error_type)+
  geom_line(aes(x = q, y = Value ,col = method),size = 1.5) + theme_bw()+
  scale_color_manual(values=c(3,4))+
    theme(legend.position = "bottom")
# Test iNEXTPD
S <- 100
set.seed(123)
test_tr <- force.ultrametric(rtree(n = S,rooted = T), method=c("extend"))
set.seed(456)
test_data <- rmultinom(n = 2,size = S*4,prob = rexp(S))
rownames(test_data) <- test_tr$tip.label
colnames(test_data) <- paste0("Site",1:ncol(test_data))
data = test_data;tree = test_tr;datatype = 'abundance';q = c(0,1,2);nboot = 50
endpoint = NULL; knots = 40; size = NULL; plot.type = 1; conf = 0.95; nboot = 50;reftime=NULL
a<-Sys.time()
new <- inextPD(data = test_data,tree = test_tr,datatype = 'abundance',q = c(1),nboot = 50)
b<-Sys.time()
b-a



