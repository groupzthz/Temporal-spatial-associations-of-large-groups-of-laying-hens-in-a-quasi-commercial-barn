# the SN were calculated separately for each week on the ubelix,
# which means there is an .Rdata file for each week. We need to put these together. 
# The script puts all the results together to then be able to work on it.
# This script is specific for the WINTEGARDEN
# NB: number of the pen must be changed

# LIBRARIES ####
library(asnipe)      # package that implements the Gaussian Mixture models
library(igraph)
library(qgraph)
library(data.table)
library(vegan)
library(corrplot)
library(readr)
library(dplyr)
library(ggplot2)
library(lme4)
library(scales)

# LOAD VARIABLES ####
pen='pen16'
plot.SN=FALSE

setwd(paste('data/', pen, '/WG/', sep=''))
# list of files
files = list.files(full.names = TRUE, pattern = "\\.RData$")
# load hens' IDS
hens <- readRDS(paste('data/other/hens_list_', pen, '_WG.rds', sep=''))
hens <- hens[!hens == 16162]

# weeks' number
weeks.collection <- c(44, 45, 46, 47, 48, 49, 50, 51, 52, 1, 2, 3, 4, 5, 7,
                      8, 9, 10, 11, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23,
                      24, 25, 26, 27, 28, 29)

age = seq(20,57,1)      # list with the age of the hens in weeks
age = age[! age %in% c(34,43)]  # these weeks are excluded as there was no data collection
# for pen 19 we have less weeks
if (pen == 'pen19'){
  weeks.collection = weeks.collection[! weeks.collection %in% c(seq(44,51,1))]
  age = age[! age %in% c(seq(20,27,1))]
}

`%notin%`<-Negate('%in%')

# in some weeks, some hens have little data. These individuals are excluded. 
# This file includes this information.
keep = read.csv(paste("./other/keep_weeks_SN_", pen, ".csv", sep=''))


in.WG = data.frame(matrix(nrow = 0, ncol = 2))
colnames(in.WG) <- c('date', 'numb hens')


# MERGE ASSOCIATION MATRICES ####################################
all.p.values.WG = c()
networks.WG <- array(0, c(length(files), length(hens), length(hens)))
networks.zeros.WG <- array(0, c(length(files), length(hens), length(hens)))
for (w in seq_along(weeks.collection)){
  print(weeks.collection[w])
  load(paste('association_matix_',pen,'_week', weeks.collection[w], '_WG.RData', sep=''))
  if (pen == 'pen16'){
  week_network <- week_network[-c(162), -c(162)]
  }
  keep.week = keep[keep$week== weeks.collection[w], ]   
  hens.keep = unique(keep.week$hen)
  difs <- setdiff(as.vector(hens), as.vector(hens.keep))
  if (pen == 'pen16'){
    if (length(difs) == 1){
      if (difs == 16162){
        difs = c()
        }
    } else {
      difs <- difs[!difs == 16162]
    }
  }
  all.p.values.WG[paste('week', weeks.collection[w], sep='_')] <- tail(p.v, n=1)
  discard.hen =unique(c(difs, elim.hens))
  diag(week_network) <- 0     # set the diagonal as 0
  week_network[as.character(discard.hen), ]<- 0  # data of flagged hens are 
  week_network[, as.character(discard.hen)]<- 0   # set to 0s
  networks.zeros.WG[w,,] = (week_network)
  diag(week_network) <- NA     # set the diagonal as NA
  week_network[as.character(discard.hen), ]<- NA   # data of flagged hens are 
  week_network[, as.character(discard.hen)]<- NA   # set to NAs
  networks.WG[w,,] = week_network
}

mean(all.p.values.WG)
min(all.p.values.WG)
max(all.p.values.WG)

nets.WG <- apply(networks.WG, 1, function(x)
{graph.adjacency(x, mode='undirected', diag=FALSE, weighted = TRUE)})

# PLOT SOCIAL NETWORKS ##############################
## here we plot and save all the SN if plot was set to True
if (plot.SN) {
  for (n in seq_along(nets.WG)){
    jpeg(filename = paste("SocialNetwork_week", weeks.collection[n],'_WG.jpeg', sep=''),
         width = 950, height = 500)
    par(mfrow = c(1, 2))
    hist(E(nets.WG[[n]])$weight, xlab='Weight', main='', breaks=20) # plot weight
    qgraph(networks.WG[n,,], minimum=0.10, cut=0.15, labels=TRUE)  #plot SN
    mtext(paste("Week", weeks.collection[n], 'in pen 16', sep = ' '), side = 3,
          line = -2, outer = TRUE)
    dev.off()
  }
}

# CORRELATIONS ACROSS WEEKS #########################
# We want to discover how similar are the association matrices, i.e., the SN
# we use the Mantel test. We contrast all the SN against each other.

all.comb =t(combn(length(weeks.collection), 2))   # t() is to transpose the table
all.corr.WG <- array(0, c(length(weeks.collection), length(weeks.collection)))
all.corr.pvalues.WG = c()

# do the correlations 
for (row in 1:nrow(all.comb)){
  print(row)
  first = all.comb[row, 1]
  second = all.comb[row, 2]
  st.week = networks.zeros.WG[first,,]
  nd.week = networks.zeros.WG[second,,]
  corr.SN = mantel(st.week, nd.week, method = "pearson", permutations = 1000)
  all.corr.WG[first, second] = corr.SN$statistic
  all.corr.pvalues.WG <- append(all.corr.pvalues.WG, corr.SN$signif)
}

#mirror the matrix
all.corr.WG[lower.tri(all.corr.WG)] = t(all.corr.WG)[lower.tri(all.corr.WG)]
rownames(all.corr.WG) <- age
colnames(all.corr.WG) <- age

# plot the correlation matrix
corrplot.mixed(all.corr.WG, is.corr=FALSE, col.lim = c(0, 1), lower='color',
               upper='pie', tl.col = 'black', tl.cex=1.1, cl.ratio=0.1, cl.cex=1,
               lower.col = colorRampPalette(c("khaki","yellowgreen", "springgreen4"))(10),
               upper.col = colorRampPalette(c("khaki", "yellowgreen", "springgreen4"))(10))

min((all.corr.WG)[lower.tri(all.corr.WG)])
max((all.corr.WG)[lower.tri(all.corr.WG)])
mean((all.corr.WG)[lower.tri(all.corr.WG)])

# adjusting p.value with FDR correction
p.values.WG.adj = p.adjust(all.corr.pvalues.WG, method = 'fdr', n=length(all.corr.pvalues.WG))
hist(p.values.WG.adj, breaks = 100)
max(p.values.WG.adj)

# ASSOCIATION INDEX EXTRACTION #######################
# extract ties across weeks
# the df "ties.WG" has all the AI of all couples along the weeks

ties.WG =  data.frame() 
for (n in seq_along(weeks.collection)){
  dur = networks.WG[n,,]
  dur=replace(dur, is.na(dur), 100) # replace NA with 100 to then recognise them
  diag(dur) <- NA
  dur[lower.tri(dur)] <- NA
  up.dur = c(t(dur))
  up.dur <- na.omit(up.dur)
  ties.WG <- rbind(ties.WG, up.dur)
}

colnames(ties.WG) = 
  unlist(apply(t(combn(hens, 2)) , 1, paste , collapse = "-" ))

ties.WG[ties.WG==100] <-NA      # values of 100 are substituted with Nas

model.ai.WG = ties.WG
model.ai.WG$hen.age = age
model.ai.WG = reshape2::melt(model.ai.WG, id.vars='hen.age', variable.name='couple', 
                             value.name = 'Ass.index')
couples = unlist(strsplit(as.character(model.ai.WG$couple), '-'))
model.ai.WG$hen1 = factor(couples[c(TRUE, FALSE)])
model.ai.WG$hen2 = factor(couples[c(FALSE, TRUE)])

sum(is.na(model.ai.WG))

model.ai.WG= model.ai.WG[complete.cases(model.ai.WG),] # rows with NAs are dropped

# statistics about the association indices
mean(model.ai.WG$Ass.index)
sd(model.ai.WG$Ass.index)
max(model.ai.WG$Ass.index)
min(model.ai.WG$Ass.index)

# distribution of AI
ggplot(model.ai.WG, aes(x=Ass.index)) + 
  geom_histogram(fill='azure3', bins=40) + 
  xlab('Distribution association index in WG') + 
  ylab('Frequency') + 
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  scale_y_continuous(labels = scales::label_number(scale_cut = cut_short_scale()))

# COMMUNITIES ####
# We want to extract communities 
# we try out different algorithms, to see which one fits best for our data.
# the best will be the one with the highest modularity
# function that applies a certain community function to the networks
# and extracts different variables -- returns a list
community <- function(nets, fun){
  comm = lapply(nets, function(x){fun(x)})  #applies the function "fun" to the network net
  modularity.comm = lapply(comm, function(x){modularity(x)})
  n.communities = lapply(comm, function(x){length(x)})
  size.community = lapply(comm, function(x){sizes(x)})
  return(list(modularity = list(modularity.comm), 
              n.communities = list(n.communities), 
              size = list(size.community)))
}
# different algorithms
communities.greedy.WG = community(nets.WG, fastgreedy.community)
communities.labelprop.WG = community(nets.WG, cluster_label_prop)
communities.eigen.WG = community(nets.WG, cluster_leading_eigen)
communities.walktrap.WG = community(nets.WG, cluster_walktrap)
communities.infomap.WG = community(nets.WG, cluster_infomap)

# put the results together in a dataframe
communities = list(communities.greedy.WG, communities.labelprop.WG, communities.eigen.WG, 
                   communities.walktrap.WG, communities.infomap.WG)

name.communities = c('Fast greedy', 'Propagating labels', 'Leading eigenvector', 
                     'Random walks', 'Infomap')

mod <- n.comm <- sizes <- name.comm <- list()   #initialize several empty lists

# we want to create a data frame out of all the parameters calculated for the communities
for (comm in seq_along(communities)) {
  mod = append(mod, unlist(communities[[comm]]$modularity))
  n.comm = append(n.comm, unlist(communities[[comm]]$n.communities))
  all.sizes = communities[[comm]]$size
  for (s in all.sizes[[1]]) {
    sizes = append(sizes, mean(s))
  }
  name.comm = append(name.comm, rep(name.communities[comm], 36))
}

# create the df
data.communities.WG = do.call(rbind, Map(data.frame, modularity=mod,
                                         n.communities=n.comm, mean.size=sizes,
                                         name.alg=name.comm))

boxplot(data.communities.WG$modularity ~ data.communities.WG$name.alg,
        ylab="Modularity" , xlab="Algorithm")

boxplot(data.communities.WG$n.communities ~ data.communities.WG$name.alg,
        ylab="# communities" , xlab="Algorithm")

boxplot(data.communities.WG$mean.size ~ data.communities.WG$name.alg,
        ylab="mean size communities" , xlab="Algorithm")


