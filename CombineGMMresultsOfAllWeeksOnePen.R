# the SN were calculated separately for each week using the cluster,
# which means there is an .Rdata file for each week. 
# We need to put these results together. The script does this.
# Then it calculates the social networks and does several analyses on them.
# NB: number of the pen must be changed

# LIBRARIES & FUNCTIONS ####
library(asnipe)      # package that implements the Gaussian Mixture models
library(igraph)
library(qgraph)
library(data.table)
library(vegan)        # for the mantel test
library(corrplot)     # to plot correlation matrices & co.
library(reshape2)
library(ggplot2)
library(lmerTest)
library(plyr)
library(dplyr)
library(scales)
library(ggeffects)
library(sjPlot)


pen = 'pen17'

# LOAD VARIABLES ####
setwd(paste('E:/chickenrun/new_data2308/Solutions/data_GMM/', pen, sep=''))
plot.SN = FALSE  # if to plot the social network
save.AM <- FALSE # if to save the association matrices

# list of files - one for each week
files = list.files(full.names = TRUE, pattern = "\\.RData$")
# extract the hens identity
load(files[2])
#rm(list=setdiff(ls(), c("hens", "files", "pen")))

# extract the weeks' number
gmm.16.week = fread(paste('data_GMM_',pen,'.csv', sep=''), select = "week")
gmm.16.week = unlist(unique(gmm.16.week$week))
age = seq(20,57,1)      # list with the age of the hens in weeks
age = age[! age %in% c(34,43)]  # these weeks are excluded as there was no data collection
# for pen 19 we have less weeks
if (pen == 'pen19'){
  gmm.16.week = gmm.16.week[! gmm.16.week %in% c(seq(44,51,1))]
  age = age[! age %in% c(seq(20,27,1))]
}


# in some weeks, some hens have little data. These individuals are excluded. 
# This file includes this inforation.
keep = read.csv(paste("E:/chickenrun/new_data2308/Solutions/keep_weeks_SN_", 
                      pen, ".csv", sep=''))

# LOAD ENVIRONMENT ####
pen = 'pen20'
setwd(paste('E:/chickenrun/new_data2308/Solutions/data_GMM/', pen, '/Renv', sep=''))
load(paste('result_analysis_',pen, '.RData', sep=''))

# MERGE ASSOCIATION MATRICES ####################################
`%notin%`<-Negate('%in%')
all.p.values = c()        # storing p-value of randomness of SN
# stacked matrix to store the association matrices (AM)
# we have 2 because there is a function that is not accepting NAs in the script.
networks <- array(0, c(length(gmm.16.week), length(hens), length(hens)))  #this will be with NAs
networks.zeros <- array(0, c(length(gmm.16.week), length(hens), length(hens))) #this is with the zeros

# loop to extract the AM
for (w in seq_along(gmm.16.week)){
  print(gmm.16.week[w])
  load(paste('./result_GMM',pen,'_week', gmm.16.week[w], '.RData', sep=''))
  # we flag hens we do not have much data 
  keep.week = keep[keep$week==week, ]   
  hens.keep = unique(keep.week$hen)
  difs <- setdiff(as.vector(hens), as.vector(hens.keep))
  # in pen 16 we exclude this individual
  if (pen == 'pen16'){
    if (16052 %notin% difs){
      difs<-append(difs, 16052)
    }
  }
  discard.hen =unique(c(difs, elim.hens)) # hens that need to be excluded
  print(hens[0:10])
  all.p.values[paste('week', week, sep='_')] <- tail(p.v, n=1)
  diag(week_network) <- 0     # set the diagonal as 0
  week_network[as.character(discard.hen), ]<- 0  # data of flagged hens are 
  week_network[, as.character(discard.hen)]<- 0   # set to 0s
  networks.zeros[w,,] = week_network   # save the AM in the array with zeros
  diag(week_network) <- NA     # set the diagonal as NA
  week_network[as.character(discard.hen), ]<- NA   # data of flagged hens are 
  week_network[, as.character(discard.hen)]<- NA   # set to NAs
  networks[w,,] = week_network
  if (save.AM) {
    rownames(week_network) = hens   # set row and columns names
    colnames(week_network) = hens
    write.csv(week_network, paste("Association_matrix_week", gmm.16.week[w], 
                                  pen, '.csv', sep='_') , row.names=TRUE)
  }
}  

# remove from environment everything we do not need
rm(list=setdiff(ls(), c("hens", "files", "all.p.values","age", "networks.zeros",
                        "all.weeks", "networks", "gmm.16.week", "pen")))

# check all the p-values for random SN 
mean(all.p.values)
max(all.p.values)
min(all.p.values)
plot(all.p.values~age)

# extract SN
nets <- apply(networks.zeros, 1, function(x){
  graph_from_adjacency_matrix(x, mode='undirected', 
                              diag=FALSE, weighted = TRUE)}) # SN

# PLOT SOCIAL NETWORKS ####
# here we plot and save all the SN if plot was set to True
if (plot.SN) {
  for (n in seq_along(nets)){
    jpeg(filename = paste("SocialNetwork_week", gmm.16.week[n],'.jpeg', sep=''),
         width = 950, height = 500)
    par(mfrow = c(1, 2))
    hist(E(nets[[n]])$weight, xlab='Weight', main='', breaks=20) # plot weight
    qgraph(networks[n,,], minimum=0.10, cut=0.15, labels=TRUE)  #plot SN
    mtext(paste("Week", gmm.16.week[n], 'in pen 17', sep = ' '), side = 3,
          line = -2, outer = TRUE)
    perc.high = which(networks[n,,] > 0.1)
    print(length(perc.high)/length(networks[n,,]))
    dev.off()
  }}

# CORRELATIONS ACROSS WEEKS ####
# We want to discover how similar are the association matrices, i.e., the SN
# we use the Mantel test. We contrast all the SN against each other. 

all.comb=t(combn(length(age), 2))   # table of contrasts to do, t() is to transpose the table 
all.corr <- array(0, c(length(gmm.16.week), length(gmm.16.week))) # table to save correlation indices
all.corr.pvalues = c()  # p-values of correlations

# do the correlations
for (row in 1:nrow(all.comb)){
  print(row)
  first = all.comb[row, 1]
  second = all.comb[row, 2]
  st.week = networks.zeros[first,,]     # select the AMs from the stacked matrix
  nd.week = networks.zeros[second,,]
  corr.SN = mantel(st.week, nd.week, method = "pearson", 
                   permutations = 1000, na.rm=TRUE) # correlation with Mantel test
  all.corr[first, second] = corr.SN$statistic    # save the results
  all.corr.pvalues <- append(all.corr.pvalues, corr.SN$signif)
}

# mirror the matrix
all.corr[lower.tri(all.corr)] = t(all.corr)[lower.tri(all.corr)]
rownames(all.corr) <- age    # row and columns names
colnames(all.corr) <- age

# plotting the correlation matrix - double plot
corrplot.mixed(all.corr, is.corr=FALSE, col.lim = c(0, 1), lower='color',
               upper='pie', tl.col = 'black', tl.cex=1.1, cl.ratio=0.1, cl.cex=1,
               lower.col = colorRampPalette(c("khaki","yellowgreen", "springgreen4"))(10),
               upper.col = colorRampPalette(c("khaki", "yellowgreen", "springgreen4"))(10))

min((all.corr)[lower.tri(all.corr)])
max((all.corr)[lower.tri(all.corr)])
mean((all.corr)[lower.tri(all.corr)])

# adjusting p.value with FDR correction --> we did many comparisons
p.values.adj = p.adjust(all.corr.pvalues, method = 'fdr', 
                        n=length(all.corr.pvalues))
hist(p.values.adj, breaks = 10)

# ASSOCIATION INDEX EXTRACTION #######################
# extract ties across weeks
# the df "ties" has all the AI of all couples along the weeks

ties =  data.frame() 
for (n in seq_along(gmm.16.week)){
  dur = networks[n,,]
  dur=replace(dur, is.na(dur), 100) # replace NA with 100 to then recognise them
  diag(dur) <- NA
  dur[lower.tri(dur)] <- NA
  up.dur = c(t(dur)) 
  up.dur <- na.omit(up.dur)
  ties <- rbind(ties, up.dur)
}

# set the columns names of couples
colnames(ties) = 
  unlist(apply(t(combn(hens, 2)) , 1, paste , collapse = "-" ))

ties %>%
  mutate(across(colnames(ties), na_if, 100)) -> ties

### Model association index along weeks ####
# first we prepare the data
model.ai.df = ties
model.ai.df$hen.age = age
model.ai.df = melt(model.ai.df, id.vars='hen.age', variable.name='couple', 
                    value.name = 'Ass.index')
# we extract the nodes that are connected by a certain tie
couples = unlist(strsplit(as.character(model.ai.df$couple), '-')) 
model.ai.df$hen1 = factor(couples[c(TRUE, FALSE)])  #first hen of the couple
model.ai.df$hen2 = factor(couples[c(FALSE, TRUE)])  # second hen of the couple
model.ai.df= model.ai.df[complete.cases(model.ai.df),] # drop NA

# statistics about the association indices
mean(model.ai.df$Ass.index)
sd(model.ai.df$Ass.index)
max(model.ai.df$Ass.index)
min(model.ai.df$Ass.index)

# plot distribution of all the association indices
ggplot(model.ai.df, aes(x=Ass.index)) + 
  geom_histogram(fill='azure3', bins = 40) + 
  xlab('Distribution association index') + 
  ylab('Frequency') + 
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  scale_y_continuous(labels = scales::label_number(scale_cut = cut_short_scale()))

# model fitted with lmerTest - for pvalues
model.ai.reml = lmer(Ass.index~hen.age + (1|hen1) + (1|hen2), model.ai.df, REML=T)
detach(name="package:lmerTest")

# model fitted with lme4
model.ai = lmer(Ass.index~hen.age + (1|hen1) + (1|hen2), model.ai.df)

# check assumptions
# homogeneity & hetereogenity of residuals
hist(residuals(model.ai), probability = T, breaks=20)
x=seq(from=min(residuals(model.ai)),
      to=max(residuals(model.ai)),length.out=100)
lines(col = "steelblue",x=x, y=dnorm(x, mean=0, sd=sd(residuals(model.ai))))

qqnorm(residuals(model.ai))
qqline(residuals(model.ai), col = "steelblue", lwd = 2)  

plot(x=fitted(model.ai), y=residuals(model.ai), pch=19)
abline(h=0, lty=3)

# check model stability
# functions supplied by Roger Mundry for calculating 
# model stability and confidence intervals
source("E:/Elisas PhD stuff/linear_model_course_roger/R_course/code/glmm_stability.r")
source("E:/Elisas PhD stuff/linear_model_course_roger/R_course/code/boot_glmm.r")
model.stab=glmm.model.stab(model.res=model.ai, contr=NULL, para=F,
                          data=NULL)
table(model.stab$detailed$lme4.warnings)
table(model.stab$detailed$opt.warnings)
table.model.stab = round(model.stab$summary[, -1], 5)

# bootstrap confidence intervals
boot.model=boot.glmm.pred(model.res=model.ai, excl.warnings=F,
                         nboots=1000, para=T, n.cores=3, 
                         resol=1000, level=0.95)
round(boot.model$ci.estimates, 5)


# see the results of the model
summary(model.ai.reml)

# tab the results of the model
sjPlot::tab_model(model.ai, show.df = TRUE, show.re.var = FALSE, 
                  show.icc=FALSE, dv.labels = '', digits = 5)

# plot the model
to_plot_model = ggpredict(model.ai, terms ='hen.age')
ggplot() + 
  geom_boxplot(aes(x=hen.age, y=Ass.index, group=hen.age),
               fill="azure3", 
               outlier.colour = 'grey', 
               data = model.ai.df, 
               outlier.size =0.5) + 
  geom_ribbon(data = to_plot_model, 
              aes(x=x, ymin=conf.low, ymax=conf.high), 
              alpha=0.3, 
              fill='black') +
  geom_line(data=to_plot_model, 
            aes(x=x, y=predicted), 
            col='black', 
            lwd=1.3, 
            lineend = "round", 
            alpha=0.7)+
  xlab('Hens age in weeks') + 
  ylab('Association index') + 
  theme_bw() + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold")) +
  scale_y_continuous(limits = c(0, 0.12))

# STRONG AND DURABLE TIES (SDT) ####
# mean association index per individual per week
mean.ai = apply(networks, 1, function(x) {colMeans(x, na.rm=TRUE)})    
rownames(mean.ai) <- hens
colnames(mean.ai) <- age
mean.ai = t(mean.ai) # transpose the data
# calculate mean per individual 
tot.mean.ai = t(as.data.frame(colMeans(mean.ai, na.rm = TRUE)))   

durable.mean.wk = rowMeans(ties, na.rm = TRUE)  # mean AI within weeks
durable.sd.wk <- apply(ties, 1, sd, na.rm=TRUE)  # stand.dev AI within weeks
# 3rd quantile of AI within week, we'll use this to calculate SDT
durable.high.lim = durable.mean.wk + durable.sd.wk 

# we mark with 1 the AIs that are greater or equal than the 3rd quantile
most.durable = ifelse(ties >= durable.high.lim, 1, 0)
# count the 1 by columns (pair)
tot1 = unlist(apply(most.durable, 2, function(x){sum(x == 1)})) 

## SDT for 36 weeks ####
# select only the couples that have an AI higher than 3rd quantile for 36 weeks
weeks.friends = 36
friends = subset(tot1,tot1>= weeks.friends) 
length(friends)/length(tot1)   # proportion of these ties

# figure out who are the individuals involved in those ties
friends.ind36= unlist(strsplit(names(friends), '-'))
friendlier.hen = unique(friends.ind36) # hens involved
friends.ind36 = table(friends.ind36)   # count how much each individual appears   
friends.ind36 = friends.ind36[order(friends.ind36, decreasing = TRUE)]
print(friends.ind36)
# proportion hens have at least one strong tie for all of time
length(friends.ind36)/length(hens)
# individuals with highest number of ties, to be checked manually
numb.friends36 = names(friends.ind36[1:4])   

# plotting the social network of these limited amount of ties
# first need to prepare the df for plotting
plotting.net36 = ties[, names(friends)]
plotting.net36$week = gmm.16.week
plotting.net36$age = age
plotting.net36 = melt(plotting.net36, id.vars=c('week', 'age'), 
                    variable.name='couple', 
                    value.name = 'Ass.index')

plotting.net36$couple = as.character(plotting.net36$couple)
plotting.net36 %>%
  group_by(couple) %>% summarise(avg = mean(Ass.index)) -> plot.df 

# dataframe
couples = unlist(strsplit(plot.df$couple, '-'))
couples1 = couples[c(TRUE, FALSE)]
couples2 = couples[c(FALSE, TRUE)]
plot.df$hen1 = unlist(lapply(couples1, function(x){as.integer(substring(x, 3, 5))}))
plot.df$hen2 = unlist(lapply(couples2, function(x){as.integer(substring(x, 3, 5))}))
plot.df = plot.df[, 2:4]
plot.df$hen1 = as.character(plot.df$hen1)
plot.df$hen2 = as.character(plot.df$hen2)
colnames(plot.df) <- c('weight', 'from', 'to')

# plot the social network for 36 weeks 
qgraph(plot.df[, c(2,3,1)], directed=FALSE, layout='spring', 
       vsize=7, shape='circle', color='azure3', label.cex=1.3,
       edge.color='black', alpha=1, cut=0.11)

# calculate the number of ties for each node in this reduced SN
# CENTRALITY
reduced.SN36 <- plot.df[, c(2, 3, 1)]
reduced.SN36 = graph_from_data_frame(reduced.SN36)
degree36 = tibble::enframe(degree(reduced.SN36))

mean(degree36$value)
min(degree36$value)
max(degree36$value)

# look association indices over weeks for durable ties
mean.plot.net <- plotting.net36 %>% 
  group_by(age) %>% 
  summarise(Ass.index=mean(Ass.index))

# plotting
ggplot() +
  geom_point(plotting.net36, mapping=aes(x=age, y=Ass.index, col=couple)) + 
  geom_line(mean.plot.net, mapping=aes(x=age, y=Ass.index), color='black', linewidth=1) +
  geom_point(mean.plot.net, mapping=aes(x=age, y=Ass.index), color='black', size=2) +
  theme(legend.position = 'none') 

### SDT 36 weeks trend with hens age ####

# model fitted with lmerTest - for pvalues
library(lmerTest)
model.sdt36.reml = lmer(Ass.index~age + (1|couple), plotting.net36, REML=T)
detach(name="package:lmerTest")

# model fitted with lme4
model.sdt36 = lmer(Ass.index~age + (1|couple), plotting.net36)

# check assumptions
# homogeneity & hetereogenity of residuals
hist(residuals(model.sdt36), probability = T, breaks=20)
x=seq(from=min(residuals(model.sdt36)),
      to=max(residuals(model.sdt36)),length.out=100)
lines(col = "steelblue",x=x, y=dnorm(x, mean=0, sd=sd(residuals(model.sdt36))))

qqnorm(residuals(model.sdt36))
qqline(residuals(model.sdt36), col = "steelblue", lwd = 2)  

plot(x=fitted(model.sdt36), y=residuals(model.sdt36), pch=19)
abline(h=0, lty=3, col = "steelblue")

# check model stability
# functions supplied by Roger Mundry for calculating 
# model stability and confidence intervals
source("glmm_stability.r")
source("boot_glmm.r")

model.stab36=glmm.model.stab(model.res=model.sdt36, contr=NULL, para=F,
                           data=NULL)

table(model.stab36$detailed$lme4.warnings)
table(model.stab36$detailed$opt.warnings)
table.model.stab = round(model.stab36$summary[, -1], 5)

# bootstrap confidence intervals
boot.model36=boot.glmm.pred(model.res=model.sdt36, excl.warnings=F,
                          nboots=1000, para=T, n.cores=3, 
                          resol=1000, level=0.95)
round(boot.model36$ci.estimates, 5)

# see the results of the model
summary(model.sdt36.reml)

# plot the model
to_plot_model36 = ggpredict(model.sdt36, terms ='age')
ggplot() + 
  geom_boxplot(aes(x=age, y=Ass.index, group=age),
               fill="azure3", 
               outlier.colour = 'grey', 
               data = plotting.net36, 
               outlier.size =0.5) + 
  geom_ribbon(data = to_plot_model36, 
              aes(x=x, ymin=conf.low, ymax=conf.high), 
              alpha=0.3, 
              fill='steelblue') +
  geom_line(data=to_plot_model36, 
            aes(x=x, y=predicted), 
            col='steelblue', 
            lwd=1.3, 
            lineend = "round", 
            alpha=0.7)+
  xlab('Hens age in weeks') + 
  ylab('Association index for SDT') + 
  theme_bw() + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold")) #+
  #scale_y_continuous(limits = c(0, 0.12))

## SDT for 27 weeks ####
# select only the couples that have an AI higher 
# than 3rd quantile for 27 weeks or more

weeks.friends = 36 * 0.75
friends = subset(tot1,tot1>= weeks.friends) 
length(friends)/length(tot1)   # proportion of these ties

# figure out who are the individuals involved in those ties
friends.ind= unlist(strsplit(names(friends), '-'))
friendlier.hen = unique(friends.ind) # hens involved
friends.ind = table(friends.ind)   # count how much each individual appears   
friends.ind = friends.ind[order(friends.ind, decreasing = TRUE)]
print(friends.ind)
# prop hens have at least one strong tie for all of time
length(friends.ind)/length(hens)
# proportion of ties in which the hens selected above are participating
p = sum(friends.ind[numb.friends36])/sum(friends.ind)
print(p)

plotting.net = ties[, names(friends)]
plotting.net$week = gmm.16.week
plotting.net$age = age
plotting.net = melt(plotting.net, id.vars=c('week', 'age'), variable.name='couple', 
                    value.name = 'Ass.index')
plotting.net$couple = as.character(plotting.net$couple)

plotting.net %>%
  group_by(couple) %>% summarise(avg = mean(Ass.index)) -> plot.df 
couples = unlist(strsplit(plot.df$couple, '-'))
couples1 = couples[c(TRUE, FALSE)]
couples2 = couples[c(FALSE, TRUE)]
plot.df$hen1 = unlist(lapply(couples1, function(x){as.integer(substring(x, 3, 5))}))
plot.df$hen2 = unlist(lapply(couples2, function(x){as.integer(substring(x, 3, 5))}))
plot.df = plot.df[, 2:4]
plot.df$hen1 = as.character(plot.df$hen1)
plot.df$hen2 = as.character(plot.df$hen2)
colnames(plot.df) <- c('weight', 'from', 'to')

# social network for "27-weeks SDT
qgraph(plot.df[, c(2,3,1)], directed=FALSE, layout='spring', 
       vsize=4, shape='circle', color='#f9b208', label.cex=1.5,
       edge.color='black', cut=0.12)

reduced.SN27 <- plot.df[, c(2, 3, 1)]
reduced.SN27 = graph_from_data_frame(reduced.SN27)
degree27 = tibble::enframe(degree(reduced.SN27))

mean(degree27$value)
min(degree27$value)
max(degree27$value)

### SDT 27 weeks trend with hens age ####

# model fitted with lmerTest - for pvalues
library(lmerTest)
model.sdt27.reml = lmer(Ass.index~age + (1|couple), plotting.net, REML=T)
detach(name="package:lmerTest")

# model fitted with lme4
model.sdt27= lmer(Ass.index~age + (1|couple), plotting.net)

# check assumptions
# homogeneity & hetereogenity of residuals
hist(residuals(model.sdt27), probability = T, breaks=20)
x=seq(from=min(residuals(model.sdt27)),
      to=max(residuals(model.sdt27)),length.out=100)
lines(col = "steelblue",x=x, y=dnorm(x, mean=0, sd=sd(residuals(model.sdt27))))

qqnorm(residuals(model.sdt27))
qqline(residuals(model.sdt27), col = "steelblue", lwd = 2)  

plot(x=fitted(model.sdt27), y=residuals(model.sdt27), pch=19)
abline(h=0, lty=3, col = "steelblue")

# check model stability
# functions supplied by Roger Mundry for calculating 
# model stability and confidence intervals
model.stab27=glmm.model.stab(model.res=model.sdt27, contr=NULL, para=F,
                             data=NULL)

table(model.stab36$detailed$lme4.warnings)
table(model.stab36$detailed$opt.warnings)
table.model.stab = round(model.stab36$summary[, -1], 5)

# bootstrap confidence intervals
boot.model27=boot.glmm.pred(model.res=model.sdt27, excl.warnings=F,
                            nboots=1000, para=T, n.cores=3, 
                            resol=1000, level=0.95)
round(boot.model27$ci.estimates, 5)

# see the results of the model
summary(model.sdt27.reml)

# plot the model
to_plot_model27 = ggpredict(model.sdt27, terms ='age')
ggplot() + 
  geom_boxplot(aes(x=age, y=Ass.index, group=age),
               fill="azure3", 
               outlier.colour = 'grey', 
               data = plotting.net, 
               outlier.size =0.5) + 
  geom_ribbon(data = to_plot_model27, 
              aes(x=x, ymin=conf.low, ymax=conf.high), 
              alpha=0.3, 
              fill='steelblue') +
  geom_line(data=to_plot_model27, 
            aes(x=x, y=predicted), 
            col='steelblue', 
            lwd=1.3, 
            lineend = "round", 
            alpha=0.7)+
  xlab('Hens age in weeks') + 
  ylab('Association index for SDT') + 
  theme_bw() + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16,face="bold")) #+
  #scale_y_continuous(limits = c(0, 0.12))


# merge the dataframe with centrality 
degree.centrality = merge(degree36, degree27, by='name', 
                          all=TRUE, suffixes = c('36', '27'))

degree.centrality = melt(degree.centrality, id.vars='name', 
                         variable.name='weeks', value.name = 'centrality')

# select 50 random individual to plot
hen.plot = sample(degree.centrality$name, 50) 
degree.centrality.plot = degree.centrality[degree.centrality$name %in% hen.plot, ]

# plot degree centrality
ggplot(degree.centrality.plot, aes(x=name, y=centrality, colour = weeks)) + 
  geom_segment(aes(x=name, xend=name, y=0, 
                   yend=centrality), color = "grey50", size = 0.75) + 
  geom_point(size=4) + xlab('Hens ID') + ylab('Centrality degree') +
  scale_color_manual(values=c('black', 'springgreen4') ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), legend.position = "none") 


# NODE STRENGTH ####
# extract the strength
strength.nodes = as.data.frame(do.call(
  cbind, lapply(nets, function(x) {strength(x)})))    

rownames(strength.nodes) <- hens
colnames(strength.nodes) <- age
strength.nodes = t(strength.nodes)  # transpose
tot.mean.str = t(as.data.frame(colMeans(strength.nodes))) # mean strength
rownames(tot.mean.str) <- 'mean'
strength.nodes[strength.nodes==0] <- NA

# calculate mean, sd and range of strength values
mean(strength.nodes, na.rm=TRUE)
sd(strength.nodes, na.rm=TRUE)
min(strength.nodes, na.rm=TRUE)
max(strength.nodes, na.rm=TRUE)

# model strength along weeks
strength.df.plot = as.data.frame(strength.nodes)
strength.df.plot$hen.age = age
strength.df.plot = melt(strength.df.plot, id.vars='hen.age', 
                        variable.name='hen', value.name = 'strength')

# plot distribution of values
ggplot(strength.df.plot, aes(x=strength)) + 
  geom_histogram(fill='azure3', bins = 40) + 
  xlab('Distribution strength values') + 
  ylab('Frequency') + 
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  scale_y_continuous(labels = scales::label_number(scale_cut = cut_short_scale()))

# model strength over weeks
model.strength = lmer(strength ~ hen.age + (1|hen), data=strength.df.plot)

# check assumptions
# homogeneity & hetereogenity of residuals
hist(residuals(model.strength), probability = T, breaks=40)
x=seq(from=min(residuals(model.strength)),
      to=max(residuals(model.strength)),length.out=100)
lines(x=x, y=dnorm(x, mean=0, sd=sd(residuals(model.strength))))

qqnorm(residuals(model.strength))
qqline(residuals(model.strength), col = "steelblue", lwd = 2)  

plot(x=fitted(model.strength), y=residuals(model.strength), pch=19)
abline(h=0, lty=3, col = "steelblue",  lwd = 3)

# check model stability
strength.stab=glmm.model.stab(model.res= model.strength, contr=NULL, para=F,
                           data=NULL)
table(strength.stab$detailed$lme4.warnings)
table(strength.stab$detailed$opt.warnings)
table.model.stab = round(strength.stab$summary[, -1], 5)

# bootstrap confidence intervals
boot.strength=boot.glmm.pred(model.res=model.strength, excl.warnings=F,
                          nboots=1000, para=T, n.cores=3, 
                          resol=1000, level=0.95)
round(boot.strength$ci.estimates, 5)
m.stab.plot(boot.strength$ci.estimates)

# see the results of the model
# fit with lmerTest for df and pvalues
library(lmerTest)
model.strength.reml = lmer(strength ~ hen.age + (1|hen), 
                           REML = T, data=strength.df.plot)
summary(model.strength.reml)

# plot results of the model
plot_model_s = ggpredict(model.strength, terms ='hen.age')
ggplot() + 
  geom_boxplot(aes(x=hen.age, y=strength, group=hen.age),
               fill="azure3", 
               outlier.colour = 'grey', 
               data = strength.df.plot, 
               outlier.size =1) + 
  geom_ribbon(data = plot_model_s, 
              aes(x=x, ymin=conf.low, ymax=conf.high), 
              alpha=0.3, fill='black') +
  geom_line(data=plot_model_s, 
            aes(x=x, y=predicted), 
            col='black', 
            lwd=1.3,
            alpha=0.7) +
  xlab('Hens age in weeks') + 
  ylab('Strength') + 
  theme_bw() + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  scale_y_continuous(limits = c(0, 20))

# nodes with higher strength across weeks 
strength.df = as.data.frame(strength.nodes)
strength.df$mean.wk = rowMeans(strength.df, na.rm=TRUE)
strength.df$sd.wk <- apply(strength.df, 1, sd, na.rm=TRUE)
strength.df$high.lim = strength.df$mean.wk + strength.df$sd.wk
strength.df$low.lim = strength.df$mean.wk - strength.df$sd.wk

strongest = ifelse(strength.df[, 1:224] > strength.df$high.lim, 1, 0)
weakest = ifelse(strength.df[, 1:224] < strength.df$low.lim, 1, 0)

par(mfrow = c(2, 1))
#myPalette <- colorRampPalette(c("#05493C","#3C725B","#3C9F69","#A6D6C8"))
col =gray.colors(2)
corrplot(strongest, is.corr=FALSE, method='color',
         tl.col = 'black', tl.cex=0.5, cl.pos=FALSE, col=col,
         na.label='square', na.label.col = 'white')
corrplot(weakest, is.corr=FALSE, method='color',
         tl.col = 'black', tl.cex=0.5, cl.pos=FALSE, col=col,
         na.label='square', na.label.col = 'black')

strongest.hens = colMeans(strongest, na.rm=TRUE)
strongest.hens.36 = subset(strongest.hens, strongest.hens>=0.99)
strongest.hens.27 = subset(strongest.hens, strongest.hens>=0.75)
names(strongest.hens.27)
names(strongest.hens.36)

# DEGREE ####
degree.nodes = as.data.frame(do.call(cbind, lapply(nets, function(x) {degree(x)})))
rownames(degree.nodes) <- hens
colnames(degree.nodes) <- age
degree.nodes = t(degree.nodes)
degree.melt = melt(degree.nodes)
colnames(degree.melt) <- c('week', 'hen', 'degree')
min(apply(degree.nodes, 1, FUN = function(x) {min(x[x > 0])})) # min excluding zeros
max(degree.nodes)
median(degree.nodes)

par(mfrow = c(2, 1))
corrplot(mean.ai, is.corr=FALSE, col.lim = c(0, 0.1), method='color',
         tl.col = 'black', tl.cex=0.5)

corrplot(strength.nodes, is.corr = FALSE, method='color', 
         tl.col = 'black', tl.cex=0.5)

names(strongest.hens.sub) %in% friendlier.hen

# ASSORTATIVITY ####
# Assortativity depending on the use of the wintergarden
assortativity = c()
prop.WG = read.csv(paste(
  'D:/chickenrun/new_data/Solutions/WeeklyProportionTimeWG_',pen, '.csv', sep=''))
prop.WG$hen.age <- mapvalues(prop.WG$week, from=gmm.16.week, to=age)
prop.WG$hen.age = as.numeric(prop.WG$hen.age)

for (i in seq_along(gmm.16.week)){
  prop.WG.week = prop.WG[prop.WG$week == gmm.16.week[i], ]
  V(nets[[i]])$name <- hens
  V(nets[[i]])$prop.WG <- prop.WG.week$propWG
  assortativity <- append(assortativity, assortativity(nets[[i]], 
                                                       V(nets[[i]])$prop.WG))
}

inWG2 = read.csv(paste('D:/chickenrun/new_data/Solutions/Average_hen_inWG_perweek_', 
                       pen,'.csv', sep=''))
inWG2$hen.age <- mapvalues(inWG2$week, from=gmm.16.week, to=age)
inWG2$hen.age = as.numeric(inWG2$hen.age)
scaleFactor <- max(prop.WG$propWG) / max(inWG2$hen)

# plotting
ggplot() + 
  geom_boxplot(data = prop.WG, 
               aes(x=ordered(hen.age), y=propWG),
               fill="#2E8B57",
               outlier.colour="grey", 
               outlier.size=1) +
  geom_line(data=inWG2, 
            aes(x=ordered(hen.age), y=hen*scaleFactor, group=1), 
            col='#f9b208', 
            linewidth=1.2) +
  scale_y_continuous(name="Prop. day spent in WG", 
                     sec.axis=sec_axis(~./scaleFactor, name="Mean # hens in WG")) +
  xlab('Hens age in weeks')  + 
  theme_bw() +  
  theme(axis.text=element_text(size=13),
        axis.title=element_text(size=14,face="bold"),
        axis.title.y.left=element_text(color="#2E8B57"),
        axis.text.y.left=element_text(color="#2E8B57"),
        axis.title.y.right=element_text(color="#f9b208"),
        axis.text.y.right=element_text(color="#f9b208"))

# COMMUNITIES ###########################################
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
communities.greedy = community(nets, fastgreedy.community)
communities.labelprop = community(nets, cluster_label_prop)
communities.eigen = community(nets, cluster_leading_eigen)
communities.walktrap = community(nets, cluster_walktrap)
communities.infomap = community(nets, cluster_infomap)

# put the results together in a dataframe
communities = list(communities.greedy, communities.labelprop, communities.eigen, 
                   communities.walktrap, communities.infomap)

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

# create the df and save it
data.communities = do.call(rbind, Map(data.frame, modularity=mod,
                                      n.communities=n.comm, mean.size=sizes,
                                      name.alg=name.comm))

saveRDS(data.communities, 
        paste('D:/chickenrun/new_data2308/Solutions/communities_',pen,'.RDS', sep=''))

# SAVE ENVIRONMENT ####
setwd(paste('E:/chickenrun/new_data2308/Solutions/data_GMM/', pen, '/Renv', sep=''))
save.image(paste('result_analysis_',pen,'.RData', sep=''))
