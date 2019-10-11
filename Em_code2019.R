# Load packages
library(ggplot2)
library(lmerTest)
library(lme4)
library(dplyr)
library(tidyverse)
library(lamW)
library(LambertW)
library(frair)


# Effluent exclusion experiment ---------------------------------

# Import data
feed_rate <- read.csv("feed_rate.csv")

# Clean up data frame a bit
feed_rate$crab_id<- as.factor(feed_rate$crab_id)
feed_rate$salmon_id<- as.factor(feed_rate$salmon_id)

# make the -'s 0
for (i in 1:length(feed_rate$crab_id)){
  if (feed_rate$dry_consumed[i] < 0) { feed_rate$dry_consumed[i] <- 0}} 

# Log transformation
feed_rate$log_dry <- log(feed_rate$dry_consumed + 1)

#Double check everything worked
str(feed_rate)


####### construct models
dry_model <- lm(log_dry ~ treatment*carapace , data = feed_rate)
anova(dry_model)
summary(dry_model)

dry_model_1 <- lm(dry_consumed ~ treatment*carapace , data = feed_rate)
anova(dry_model_1) # confirm log model and untransformed models have same result
summary(dry_model_1)

# confirm log transformed model is better
AIC(dry_model, dry_model_1)


# look at normality 
plot(dry_model)
shapiro.test(resid(dry_model))

kruskal.test(dry_consumed ~ carapace , data = feed_rate)
hist(resid(dry_model))

qqnorm(resid(dry_model))
qqline(resid(dry_model))

#### final graph
ggplot(data = feed_rate, aes(x = carapace, y = dry_consumed)) + 
  geom_point(aes(shape = treatment), size = 2) +  theme_classic() +
  geom_smooth(method = lm, colour = "black", aes(lty = treatment)) +
  xlab("Carapace Width (mm)") + ylab("Feeding rate (grams/hour)") +
  theme(axis.text=element_text(colour="black", 
                               size=13), text=element_text(colour="black", size=13)) +
  labs(lty = "Treatment", shape = "Treatment") +
  scale_linetype_discrete(name = "Treatmet", labels = c("Control", "Effluent")) +
  scale_shape_discrete(name = "Treatmet", labels = c("Control", "Effluent"))


# Flow rate calculations ---------------------------------
flow <- read.csv("flow_rate.csv")
str(flow)

mean(flow$flow_rate_ml_s[flow$treatment=="control"])
mean(flow$flow_rate_ml_s[flow$treatment=="effluent"])

mean(flow$flow_rate_ml_s[flow$tank_size=="small"])
sd(flow$flow_rate_ml_s[flow$tank_size=="small"])

mean(flow$flow_rate_ml_s[flow$tank_size=="large"])
sd(flow$flow_rate_ml_s[flow$tank_size=="large"])


# Functional response experiment ---------------------------------
oysters <- read.csv("functional_response.csv") #data file with green and red crab data for functional response
str(oysters)

# Take a quick look at the data
ggplot(data = oysters) + geom_jitter(aes(x = density, y = eaten, colour = species)) +
  scale_color_manual(values=c("green", "red")) + geom_smooth(method = lm, aes(x = density, y = eaten, colour = species))

ggplot(data = oysters) + geom_jitter(aes(x = cheliped, y = eaten, colour = species))

# GLM model ---------------------------------

#turn data into a matrix of successes and failures
success <- matrix(oysters$eaten, ncol = 1)
failures <- matrix(oysters$alive, ncol = 1)
oystmat <- cbind(success, failures)

###### Which Link is Best???
logit <- glm(oystmat ~ density + cheliped + species, data = oysters, 
             family = binomial(link = logit))
probit <- glm(oystmat ~ density + cheliped + species, data = oysters, 
              family = binomial(link = probit))
cauchit <- glm(oystmat ~ density + cheliped + species, data = oysters, 
               family = binomial(link = cauchit))
cloglog <- glm(oystmat ~ density + cheliped + species, data = oysters, 
               family = binomial(link = cloglog))
AIC(logit, probit, cauchit, cloglog)
#### We chose logit because it is the simplest and only has a slightly higher AIC than the alternatives

### which model is best ????
full <- glm(oystmat ~ density + cheliped + species, data = oysters, 
            family = binomial(link = logit))
den_che <- glm(oystmat ~ density + cheliped, data = oysters, 
               family = binomial(link = logit))
den_sp <- glm(oystmat ~ density + species, data = oysters, 
              family = binomial(link = logit))
den <- glm(oystmat ~ density, data = oysters, 
           family = binomial(link = logit))
int <- glm(oystmat ~ 1, data = oysters, 
           family = binomial(link = logit))
cara <-  glm(oystmat ~ density + carapace + species, data = oysters, 
             family = binomial(link = logit))
AIC(full, den_che, den_sp, den, int, cara)

### The best model is the one that includes density, carapace and species (cara)

# summary of all models
summary(full)
summary(den_che)
summary(den_sp)
summary(den)
summary(int)
summary(cara)

par(mfrow = c(2,2))

plot(cara)
qqplot(resid(cara))
hist(resid(cara))
shapiro.test(resid(cara)) 

par(mfrow = c(1,1))
plot(resid(cara) ~ factor(oysters$tank))

# Green functional response ---------------------------------
green <- filter(oysters, species == "green")
str(green)

# Test for type II
frair_test(eaten ~ density, data = green)

### Frair fit
frair_responses()

outII_g <- frair_fit(eaten ~ density, data = green, response = 'rogersII',
                     start = list(a = 0.2, h = 0.2), fixed = list(T=1))

# A linear fit
outI_g <- frair_fit(eaten ~ density, data = green, response = 'typeI',
                    start = list(a = 0.2), fixed=list(T=1))
# Visualise fits
plot(outII_g, pch=20, col=rgb(0,0,0,0.2), xlim=c(0,65))
lines(outII_g)
lines(outI_g, lty=3)

#####Greens resid
a <- outII_g$coefficients[1] # Get coeffs
h <- outII_g$coefficients[2]

fitsgrn <- data.frame(x = green$density) # Calculate 'a' for each value of x, where x is density of oysters

fitsgrn$Ne <- fitsgrn$x - lambertW0(a * h * fitsgrn$x * exp(-a * (1 - h * fitsgrn$x)))/(a * h) # calculate expected number of oysters eaten, based on frair flexnr equation, using lambert function

fitsgrn$actual <- green$eaten
fitsgrn$resid <- fitsgrn$Ne - fitsgrn$actual

plot(x = fitsgrn$x, y = fitsgrn$Ne)
plot(x = fitsgrn$Ne, y = fitsgrn$resid)
abline(h = 0, lty = 'dotted')

# Have a look at original fits returned by mle2 (*highly* recommended)
summary(outII_g$fit)
summary(outI_g$fit)
# Compare models using AIC
AIC(outI_g$fit,outII_g$fit)


#### BOOTSTRAP
set.seed(98765)
outII_g_boot <- frair_boot(outII_g, start = NULL, strata=NULL, nboot=2000,
                           para=TRUE, ncores=NaN, WARN.ONLY=FALSE)

# Illustrate bootlines
plot(outII_g_boot, xlim=c(0,70), ylim = c(0, 40), type='n', main='All bootstrapped lines')
lines(outII_g_boot, all_lines=TRUE)
points(outII_g_boot, pch=20)

# Illustrate bootpolys
plot(outII_g_boot, xlim=c(0,70), ylim = c(0, 40), type='n', main='Empirical 95 percent CI')
drawpoly(outII_g_boot, col=hsv(2/6,0.2, 0.8))
points(outII_g_boot, pch=20)
lines(outII_g_boot, all_lines=FALSE)


##### nlm to look at asymptote (analysis not in manuscript)
green_asymp <- nls(eaten ~ SSasymp(density, Asym, R0, lrc), 
                   data = green) # asymptotic curve with initial guesses
summary(green_asymp)



# Red functional response ---------------------------------
red <- filter(oysters, species == "red")
str(red)

# Test for type II or III
frair_test(eaten ~ density, data = red) 

### Frair fit
frair_responses()

outII_r <- frair_fit(eaten ~ density, data = red, response = 'rogersII',
                     start = list(a = 0.2, h = 0.2), fixed = list(T=1))

outIII_r <- frair_fit(eaten ~ density, data = red, response = 'flexpnr',
                      start = list(b=0.0029242, q = 1.9391532, h = 0.0312941), fixed = list(T=1))

outIII_r_2 <- frair_fit(eaten ~ density, data = red, response = 'hassIIInr',
                        start = list(b=0.0029242, c = 1.9391532, h = 0.0312941), fixed = list(T=1))

# A linear fit
outI_r <- frair_fit(eaten ~ density, data = red, response = 'typeI',
                    start = list(a = 0.2), fixed=list(T=1))
# Visualise fits
plot(outIII_r, pch=20, col=rgb(0,0,0,0.2), xlim=c(0,70))
lines(outIII_r)
lines(outI_r, lty=3)
lines(outII_r, lty = 4)
lines(outIII_r_2, lty = 5)


#### calculate predicted values?
b <- outIII_r$coefficients[1] # Get coeffs 
q <- outIII_r$coefficients[2]
h <- outIII_r$coefficients[3]

fits <- data.frame(a = b*red$density^q, x = red$density) # Calculate 'a' for each value of x, where x is density of oysters
fits$Ne <- fits$x - lambertW0(fits$a * h * fits$x * exp(-fits$a * (1 - h * fits$x)))/(fits$a * h) 
# calculate expected number of oysters eaten, based on frair flexnr equation, using lambert function

plot(x = fits$x, y = fits$Ne) # Plot 

fits$actual <- red$eaten
fits$resid <- fits$actual - fits$Ne
plot(y = fits$resid, x = fits$Ne)
abline(h = 0, lty = 'dotted')
##


# Have a look at original fits returned by mle2 (*highly* recommended)
summary(outIII_r$fit)
summary(outII_r$fit)
summary(outI_r$fit)
# Compare models using AIC
AIC(outI_r$fit,outII_r$fit, outIII_r$fit, outIII_r_2$fit)


#### BOOTSTRAP 
set.seed(604250)
outIII_r_boot <- frair_boot(outIII_r, start = NULL, strata=NULL, nboot=2000,
                            para=TRUE, ncores=NaN, WARN.ONLY=FALSE)

# Illustrate bootlines
plot(outIII_r_boot, xlim=c(0,70), ylim = c(0, 40), type='n', main='All bootstrapped lines')
lines(outIII_r_boot, all_lines=TRUE)
points(outIII_r_boot, pch=20)

# Illustrate bootpolys
plot(outIII_r_boot, xlim=c(0, 65), ylim = c(0, 70), type='n', main='Empirical 95 percent CI')
drawpoly(outIII_r_boot, border = NA, col=hsv(0/6,0.4, 0.8, alpha = 0.3))
points(outIII_r_boot, pch=20, col=hsv(0/6,0.4, 1))
lines(outIII_r_boot, lwd = 3, all_lines=FALSE, col=hsv(0/6,0.4, 1, alpha = 1))

##### nlm to look at asymptote (analysis not in manuscript)
red_asymp <- nls(eaten ~ SSasymp(density, Asym, R0, lrc), 
                 data = red) # asymptotic curve with initial guesses
summary(red_asymp)



# Functional response graph ---------------------------------


# Illustrate bootpolys
plot(outIII_r_boot, xlim=c(0, 65), ylim = c(0, 70), type='n',
     main='Functional Response Curves')
drawpoly(outIII_r_boot, border = NA, col=hsv(0/6,0.4, 0.8, alpha = 0.3))
points(outIII_r_boot, pch=20, col=hsv(0/6,0.4, 1))
lines(outIII_r_boot, lwd = 3, all_lines=FALSE, col=hsv(0/6,0.4, 1, alpha = 1))

#add green
drawpoly(outII_g_boot, border = NA, col=hsv(2/6,0.4, 0.8, alpha= 0.3))
points(outII_g_boot, pch=20,  col=hsv(2/6,0.4, 0.8, alpha= 1))
lines(outII_g_boot, lwd = 3, all_lines=FALSE, col=hsv(2/6,0.4, 0.8, alpha= 1))

# add legend
legend(x = "topleft", legend = c("Green Crabs", "Red Rock Crabs"), col = c("green", "red"), 
       lty = 1, cex = 1)