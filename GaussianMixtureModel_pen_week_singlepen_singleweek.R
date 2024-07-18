# Gaussian Mixture Model 
# the script uses the Gaussian Mixture Model to calculate social network for each pen. 
# This script only considers "start time", i.e. when the hens are detected entering a certain area of the pen. 
# Pen is divided into 2 sides for the internal part (8 areas) + wintergarden
# this is the resolution we are using
# calculates also if social network we are observing are random or not
# Applies also different community algorithms to then choose which is the best one. 

library(asnipe)      # package that implements the Gaussian Mixture models
library(igraph)
library(qgraph)
library(plotrix)
library(dplyr) 
library(ggplot2)
library(tibble)

# define the function to plot the log-log graph
loglog_frequency <- function(df) {
  # pipeline with dplyr to group by date
  df %>%
    group_by(date) %>%
    arrange(df['timestamp'], .by_group = TRUE) %>%    # order the df by time
    # calculate time difference between subsequent rows within days
    mutate(Diff = timestamp - lag(timestamp)) -> df_diff    
  sum(df_diff$Diff < 0) # neg values
  differences  = as.list(df_diff['Diff'])  # column 'diff' as list
  differences = unlist(differences)         # we remove the character attribute
  differences = differences[!is.na(differences)]   # we remove NA
  
  count_diff = as.data.frame(table(differences))               # count the values 
  count_diff$differences = as.numeric(count_diff$differences) 
  
  # plot the log-log frequency distribution
  plot(count_diff$differences, count_diff$Freq, log='xy',
       xlab = "Time difference between successive novel detections (s)",
       ylab = "Frequency (log)",)
  return(df_diff)
}


# permutation of the SN, calculation of CV, plotting histogram, 
# to check if SN is random
# permutation of data is done within day and location
permutation_network <- function(big_gbi, days.gbi, week_network, location.gbi, cv.week){ 
  network.perm = network_permutation(big_gbi, data_format = 'GBI', 
                                     association_matrix = week_network,
                                     permutations = 500000,locations = location.gbi,
                                     days = days.gbi, within_day = TRUE, within_location = TRUE)
  cv.rand <- apply(network.perm,1,function(x) {sd(x)/mean(x)})
  p.value = (500000-length(cv.rand[cv.rand < cv.week]))/500000
  #hist(cv.rand, breaks = 100, col = 'black', main = p.value)
  #abline(v=cv.week, col='red', lwd=2)
  return.fun <- c(cv.rand, p.value)
  return(return.fun)
}


# function to delete little data for a certain week of life
# (can influence the SN)
hens.to.discard <- function(hen.week, keep.df){
  keep.w = keep.df[keep.df$week == hen.week,]
  keep.w0 = keep.w[keep.w$keep == 0, ]
  hens.discard = unique(keep.w0$hen)
  return(hens.discard)
}

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


pen = 'pen16'   # these number must be changed for the different pens/weeks
week = 44

# set the working directory
setwd('')
df <- read.csv(file = paste('data_GMM_', pen, '.csv', sep=''))    # open df hen data
df = df[, c('timestamp', 'hen', 'pen', 'date', 'week', 'location')]
keep <- read.csv(file=paste('keep_weeks_SN_', pen, '.csv', sep=''))
hens = unique(df$hen)

df = df[df$week == week, ]

# set the right data types
df$date = as.Date(df$timestamp)
df[['timestamp']] = strptime(df[['timestamp']], format = "%Y-%m-%d %H:%M:%S")


## define stacked matrix W x H x H to hold two association matrices
# where W is the # of weeks of life of the hens (we will build a SN for each week)
# and H is the # of hens
#networks <- array(0, c(length(hens), length(hens)))

p.values.weeks = c()       # to store all the p.values for non-random analysis
all.cv = c()               # to store all the coeff. of variation of true networks
all.cv.rand = list()       # to store all the coeff. of variation of random netwoks

# iteration along the weeks
days = as.character(unique(df$date))
c = 0
location.gbi = c()
days.gbi = c()
elim.hens = hens.to.discard(hen.week = week, keep.df = keep)    
df = subset(df, !(hen %in% elim.hens))     # we eliminate data from hens that have few data
for (day in days){
  print(day)
  reference = as.POSIXct(paste(c(day, "2:00:00"), collapse=' '))  # reference time to calculate how much time passed 
  df.day = df[df$date == day ,]     #subset the df to the day 
  df.day$passed = as.numeric(difftime(df.day$timestamp , reference, units='secs'))
  df.gmm = df.day[, c('hen', 'location', 'passed')]    # select columns for GMMevents
  df.gmm$hen = as.numeric(df.gmm$hen)
  # Gaussian Mixture Models
  gmm_model = gmmevents(time=df.gmm$passed, identity = df.gmm$hen, 
                        location = df.gmm$location, global_ids = hens, verbose = FALSE)
  gbi = gmm_model$gbi
  location.gbi = c(location.gbi, gmm_model$metadata$Location)   #vector of locations
  days.gbi = c(days.gbi, rep(day, nrow(gbi)))    # vector of days
  if (c == 0){
    big_gbi <- gbi
    c = c + 1
  } else {
    big_gbi = rbind(big_gbi, gbi)    
  }
}
week_network = get_network(big_gbi, data_format = 'GBI')  #get the SN
rownames(week_network) = hens
colnames(week_network) = hens
cv.week = sd(week_network)/mean(week_network) # calculate coefficient variation of association indices
start_time <- Sys.time()
p.v = permutation_network(big_gbi, days.gbi, week_network, location.gbi, cv.week)  # permutate SN
end_time <- Sys.time()

save.image(paste('result_GMM',pen , '_week', week, '.RData', sep=''))


