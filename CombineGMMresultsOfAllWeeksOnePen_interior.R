# the SN were calculated separately for each week on the cluster,
# which means there is an .Rdata file for each week. We need to put these together. 
# The script puts all the results together to then be able to work on it.
# This script is specific for the INTERIOR part
# NB: number of the pen must be changed

# LOAD LIBRARIES ####
library(vegan)
library(corrplot)
library(asnipe)      # package that implements the Gaussian Mixture models
library(igraph)
library(qgraph)
library(data.table)
library(reshape2)
library(ggplot2)
library(lmerTest)
library(plyr)
library(dplyr)

# LOAD VARIABLES ####

# LOAD VARIABLES ####
pen='pen16'
plot.SN=FALSE

setwd(paste('data/', pen, '/interior/', sep=''))
# list of files
files = list.files(full.names = TRUE, pattern = "\\.RData$")
# load hens' IDS
hens <- readRDS(paste('data/other/hens_list_', pen, '_interior.rds', sep=''))
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

# MERGE ASSOCIATION MATRICES ####################################
all.p.values = c()        # storing p-value of randomness of SN
# stacked matrix to store the association matrices (AM)
networks <- array(0, c(length(weeks.collection), length(hens), length(hens)))  # this will be with NAs
networks.zeros <- array(0, c(length(weeks.collection), length(hens), length(hens))) # this will stay with the zeros

# loop to extract the AM
for (w in seq_along(weeks.collection)){
  print(weeks.collection[w])
  load(paste('association_matix_',pen,'_week', weeks.collection[w], 
             '_interior.RData', sep=''))
  keep.week = keep[keep$week==week, ]
  hens.keep = unique(keep.week$hen)
  difs <- setdiff(as.vector(hens), as.vector(hens.keep))
  if (pen == 'pen16'){
    if (16052 %notin% difs){
      difs<-append(difs, 16052)
    }
  }
  discard.hen =unique(c(difs, elim.hens))
  print(discard.hen)
  all.p.values[paste('week', week, sep='_')] <- tail(p.v, n=1)
  diag(week_network) <- 0     # set the diagonal as 0
  week_network[as.character(discard.hen), ]<- 0  # data of flagged hens are 
  week_network[, as.character(discard.hen)]<- 0   # set to 0s
  networks.zeros[w,,] = week_network
  diag(week_network) <- NA
  week_network[as.character(discard.hen), ]<- NA
  week_network[, as.character(discard.hen)]<- NA
  networks[w,,] = week_network
}

# control the p-values of the SN
mean(all.p.values)
max(all.p.values)
min(all.p.values)
plot(all.p.values~age)

nets <- apply(networks.zeros, 1, function(x)
{graph_from_adjacency_matrix(x, mode='undirected', diag=FALSE, weighted = TRUE)}) # SN

# ASSOCIATION INDEX INTERIOR PART ####
# extract ties across weeks
# the df "ties" has all the AI of all couples along the weeks

ties =  data.frame() 
for (n in seq_along(weeks.collection)){
  dur = networks[n,,]
  dur=replace(dur, is.na(dur), 100)
  diag(dur) <- NA
  dur[lower.tri(dur)] <- NA
  up.dur = c(t(dur))
  up.dur <- na.omit(up.dur)
  ties <- rbind(ties, up.dur)
}

colnames(ties) = 
  unlist(apply(t(combn(hens, 2)) , 1, paste , collapse = "-" ))
sum(is.na(ties))

ties[ties==100] <-NA

# model association index along weeks
model.ai.df = ties
model.ai.df$hen.age = age
model.ai.df = melt(model.ai.df, id.vars='hen.age', variable.name='couple', 
                   value.name = 'Ass.index')
couples = unlist(strsplit(as.character(model.ai.df$couple), '-'))
model.ai.df$hen1 = factor(couples[c(TRUE, FALSE)])
model.ai.df$hen2 = factor(couples[c(FALSE, TRUE)])
model.ai.df= model.ai.df[complete.cases(model.ai.df),]

mean(model.ai.df$Ass.index)
sd(model.ai.df$Ass.index)
max(model.ai.df$Ass.index)
min(model.ai.df$Ass.index)

# CORRELATIONS ####
# correlate the SN of the interior part across all weeks
all.comb =t(combn(length(weeks.collection), 2))   # t() is to transpose the table
all.corr.int <- array(0, c(length(weeks.collection), length(weeks.collection)))
all.corr.pvalues.int = c()

# compute the correlations
for (row in 1:nrow(all.comb)){
  print(row)
  first = all.comb[row, 1]
  second = all.comb[row, 2]
  st.week = networks.zeros[first,,]
  nd.week = networks.zeros[second,,]
  corr.SN = mantel(st.week, nd.week, method = "pearson", permutations = 1000)
  all.corr.int[first, second] = corr.SN$statistic
  all.corr.pvalues.int <- append(all.corr.pvalues.int, corr.SN$signif)
}

all.corr.int[lower.tri(all.corr.int)] = t(all.corr.int)[lower.tri(all.corr.int)]
rownames(all.corr.int) <- age
colnames(all.corr.int) <- age

# plot the correlation matrix for the interior part
corrplot.mixed(all.corr.int, is.corr=FALSE, col.lim = c(0, 1), lower='color',
               upper='pie', tl.col = 'black', tl.cex=1.1, cl.ratio=0.1, cl.cex=1,
               lower.col = colorRampPalette(c("khaki","yellowgreen", "springgreen4"))(10),
               upper.col = colorRampPalette(c("khaki", "yellowgreen", "springgreen4"))(10))


min((all.corr.int)[lower.tri(all.corr.int)])
max((all.corr.int)[lower.tri(all.corr.int)])
mean((all.corr.int)[lower.tri(all.corr.int)])

# adjusting p.value with FDR correction
p.values.int.adj = p.adjust(all.corr.pvalues.int, method = 'fdr', n=length(all.corr.pvalues.int))
hist(p.values.int.adj, breaks = 100)
max(p.values.int.adj)

# COMPARISON WINTERGARDEN-INTERNAL AREA ####
hens.int = hens
setwd()
# load environment WG - this is the result from the other script
load(paste("combined_results_WG_",pen, ".RData", sep=''))

# chech again the association indices
mean(model.ai.WG$Ass.index)
sd(model.ai.WG$Ass.index)
max(model.ai.WG$Ass.index)
min(model.ai.WG$Ass.index)

# bind together the ass. indices from the 2 different analsis
model.ai.df$data = rep('Interior', nrow(model.ai.df))
model.ai.WG$data = rep('WG', nrow(model.ai.WG))
ass.index = rbind(model.ai.df, model.ai.WG)
# plot
ggplot(data=ass.index, aes(x=Ass.index, fill=data)) + 
  geom_histogram(binwidth = 0.01, alpha=0.9, colour='black', size=0.2)+ 
  xlab('Distribution association index') + 
  ylab('Frequency') + 
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  scale_y_continuous(labels = scales::label_number(scale_cut = cut_short_scale())) + 
  scale_fill_manual(values=c('#3e424b', 'azure3')) + 
  guides(fill=guide_legend(title='Pen zone', title.theme=element_text(size=12),
                           label.theme = element_text(size=12), override.aes = list(size = 5)))

# compute correlation ass. matrices between WG and interior part within the SAME week
mantel.p = mantel.corr = c()
for (i in seq_along(weeks.collection)){
  print(i)
  SN = networks.zeros[i,,]
  rownames(SN) = hens.int   # set row and columns names
  colnames(SN) = hens.int
  SN = SN[order(row.names(x=SN)), order(colnames(x=SN))]
  SN.wg = networks.zeros.WG[i,,]
  rownames(SN.wg) = hens   # set row and columns names
  colnames(SN.wg) = hens
  SN.wg = SN.wg[order(row.names(x=SN.wg)), order(colnames(x=SN.wg))]
  #SN.wg <- SN.wg[-c(162), -c(162)]
  corr.SN = mantel(SN, SN.wg, method = "pearson", permutations = 1000)
  mantel.corr <- append(mantel.corr, corr.SN$statistic)
  mantel.p <- append(mantel.p, corr.SN$signif)
}

plot(mantel.corr)
min(mantel.corr)
max(mantel.corr)
mean(mantel.corr)

min(mantel.p)
hist(mantel.p, breaks = 100)
plot(mantel.corr~mantel.p)
abline(v=0.05, col='red')

# SAVE####
save.image(paste('result_',pen, '_corrIntWG.RData', sep=''))
