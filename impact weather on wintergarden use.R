# visits of WG compared to temperature - FOR ALL PENS TOGETHER ####

setwd()
library(glmmTMB)
library(performance)
library(DHARMa)
library(ggeffects)
library(car)
library(gt)
library(tibble)
source("glmmTMB_stability.r")
source("diagnostic_fcns.r")

# import weather data
weather = readr::read_csv2("data/weather_data.txt", col_names = TRUE)
colnames(weather) <- c('location', 'date', 'max_temp', 'min_temp', 'humidity', 'sunshine')
weather$date = as.Date(as.character(weather$date), format="%Y%m%d")
weather$max_temp = as.numeric(weather$max_temp)
weather$min_temp = as.numeric(weather$min_temp)

# import wintergarden data
prop.WG = read.csv('data/DailyProportionTimeWG.csv')
prop.WG$date = as.Date(prop.WG$date)
# merge two datasets together
data = merge(prop.WG, weather, by='date')

n.days = aggregate(data=data, date ~ hen, function(date) length(unique(date)))
prop.zero = prop.WG[prop.WG$prop_WG == 0, ]
prop.zero = as.data.frame(table(prop.zero$hen))
colnames(prop.zero) = c('hen', 'freq')
prop.zero = merge(n.days, prop.zero, by='hen')
prop.zero$noWG = ifelse(prop.zero$date == prop.zero$freq, 'yes', 'no')
table(prop.zero$noWG)

# visualize a sample of the data
sample.data = data[sample(nrow(data), 10000), ]
ggplot(data=sample.data) + geom_point(aes(x=max_temp, y=prop))
ggplot(data=sample.data) + geom_point(aes(x=sunshine, y=prop))

hist(data$prop_WG)
min(data$prop_WG)
max(data$prop_WG)
hist(data$max_temp)
hist(data$sunshine)

# scale continuous variables
data$prop_WG_tr = (data$prop*(nrow(data)-1) + 0.5)/nrow(data)
data$max_temp_z = scale(data$max_temp)
data$sunshine_z = scale(data$sunshine)

# model with beta family and logit link
model.temp = glmmTMB(prop_WG_tr~max_temp_z  * sunshine_z, 
                     family=beta_family(link='logit'), data=data)
# check dispersion
simulationOutput = simulateResiduals(fittedModel = model.temp)
testDispersion(simulationOutput)
# check collinearity
performance::check_collinearity(glmmTMB(prop_WG_tr~max_temp_z  + sunshine_z, 
                                        family=beta_family(link='logit'), data=data))
# results
summary(model.temp)
# table results
tabled.coefs <- summary(model.temp)[6]$coefficients$cond
tabled.coefs <- as.data.frame(round(tabled.coefs, 3))
row.header = c('(Intercept)', 'Max temp', 'Sunshine', 'Max temp:Sunshine' )
tabled.coefs = tibble::add_column(tabled.coefs, ' ' = row.header, .before = 'Estimate')

results.table2 <- gt(tabled.coefs) |>
  tab_header(
    title = "Weather influence on WG use")|> 
  tab_options( heading.background.color = "grey70")
print(results.table2)

# model stability
model.stab=glmmTMB.stab(model.res=model.temp, para=T, data=data)

# plotting the results
# max temperature
sample.data = data[sample(nrow(data), 5000), ]
to_plot_model_temp = ggpredict(model.temp, terms = 'max_temp_z[all]')
plot(to_plot_model_temp)

to_plot_model_temp$max_temp <- to_plot_model_temp$x *
  attr(data$max_temp_z,  'scaled:scale') + 
         attr(data$max_temp_z, 'scaled:center')  # scale back for plotting
# plot
ggplot() +
  geom_jitter(data=sample.data, aes(x=max_temp, y=prop_WG_tr), 
              alpha=0.5, col='grey50', width = 0.5) +
  geom_line(data=to_plot_model_temp, aes(x=max_temp, y=predicted),  
            lwd=1.3, col='blue') +
  theme_bw() + 
  ylim(c(0,0.15)) +
  labs(x='Maximum temperature °C', y='Daily proportion WG use') +
  theme(text = element_text(size = 15))

# sunshine %
to_plot_model_sun <- ggpredict(model.temp, terms = "sunshine_z[all]")
plot(to_plot_model_sun)

to_plot_model_sun$sun <- to_plot_model_sun$x *
  attr(data$sunshine_z,  'scaled:scale') + 
  attr(data$sunshine_z, 'scaled:center')

# plot
ggplot() +
  geom_jitter(data=sample.data, aes(x=sunshine, y=prop_WG_tr), 
              alpha=0.5, col='grey50', width = 0.7) +
  geom_line(data=to_plot_model_sun, aes(x=sun, y=predicted),  
            lwd=1.3, col='blue') +
  ylim(c(0,0.15)) +
  theme_bw() + 
  labs(x='Sunshine (%)', y='Daily proportion WG use') +
  theme(text = element_text(size = 15))

# SAVE ENVIRONMENT####
save.image(paste("effect_weather_WG_",pen,".RData", sep = ''))
