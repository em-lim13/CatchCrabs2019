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

# Effluent exclusion graph ---------------------------------

# Colour
### LEGEND SAYS TREATMET FOR SOME REASON KILL ME
ggplot(data = feed_rate, aes(x = carapace, y = dry_consumed, colour = treatment, 
                             shape = treatment, linetype = treatment)) + 
  geom_point(aes(shape = treatment, colour = treatment), size = 2.9) +  theme_classic() +
  geom_smooth(method = lm, size = 1.5, aes(lty = treatment, colour = treatment)) +
  scale_color_manual(values = c("blue", "red"), name = "Treatmet", labels = c("Control", "Effluent")) +
  xlab("Carapace Width (mm)") + ylab("Feeding rate (grams/hour)") +
  theme(axis.text=element_text(colour="black", face = "plain",
                               size = 14), text=element_text(colour="black", size=15, face = "bold")) +
  labs(lty = "Treatment", shape = "Treatment") +
  scale_linetype_discrete(name = "Treatmet", labels = c("Control", "Effluent")) +
  scale_shape_discrete(name = "Treatmet", labels = c("Control", "Effluent")) +
  theme(legend.position = c(.15, .85), legend.box.background = element_rect(colour = "black"), 
  legend.background = element_blank())

# BLACK AND WHITE
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
oysters2 <- read.csv("functional_response_2019.csv") #data file with 2019 data
oysters <- add_column(oysters, year = 2018) #add year column to 2018 data
oysters <- bind_rows(oysters, oysters2) #bind the two data frames
oysters$species <- as.factor(oysters$species) #fix variables
oysters$time <- as.factor(oysters$time) #fix variable
oysters$year <- as.factor(oysters$year)


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
green <- green[-c(20),] 
str(green)

# Make sure 2018 and 2019 are the same
#cheliped aren't signif diff
years <- lm(cheliped ~ year, data = green)
summary(years)
ggplot(green) + geom_boxplot(aes(x = year, y = cheliped))

# but 2019 crabs were bigger
years2 <- lm(carapace ~ year, data = green)
anova(years2)
summary(years2)
ggplot(green) + geom_boxplot(aes(x = year, y = carapace))

year3 <- lm(cheliped ~ carapace*year, data = green)
summary(year3)

# does year affect eaten?
ggplot(green) + geom_jitter(aes(x = density, y = proportion_eaten, colour = year))
year4 <- lm(proportion_eaten ~ density + cheliped*year, data = green)
summary(year4)

#cheli x cara graph
ggplot(green) + geom_point(aes(x = carapace, y = cheliped, colour = year)) + 
  geom_smooth(aes(x = carapace, y = cheliped, colour = year), method = lm)


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


####  OLD BOOTSTRAP
#set.seed(98765)
#outII_g_boot <- frair_boot(outII_g, start = NULL, strata=NULL, nboot=2000, para=TRUE, ncores=NaN, WARN.ONLY=FALSE)


# NEW bootstrap
set.seed(309331)
outII_g_boot <- frair_boot(outII_g, start = NULL, strata=green[,6], nboot=2000,
                           para=TRUE, ncores=NaN, WARN.ONLY=FALSE)

#outII_emd_g_boot <- frair_boot(emdII_g, start = NULL, strata=green[,6], nboot=2000, para=TRUE, ncores=NaN, WARN.ONLY=FALSE)

outII_g_boot
confint(outII_g_boot)
#outII_emd_g_boot

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
red <- red[-c(34),] # 32 no eat
red <- red[-c(24),] # 64 no eat
str(red)

# Make sure 2018 and 2019 are the same
#cheliped aren't signif diff
years_r <- lm(cheliped ~ year, data = red)
summary(years_r)
ggplot(red) + geom_boxplot(aes(x = year, y = cheliped))

# but 2019 crabs were bigger
years2_r <- lm(carapace ~ year, data = red)
anova(years2_r)
summary(years2_r)
ggplot(red) + geom_boxplot(aes(x = year, y = carapace))

year3_r <- lm(cheliped ~ carapace*year, data = red)
summary(year3_r)

# does year affect eaten?
ggplot(red) + geom_jitter(aes(x = density, y = proportion_eaten, colour = year))
year4_r <- lm(proportion_eaten ~ density + cheliped*year, data = red)
summary(year4_r)

#cheli x cara graph
ggplot(red) + geom_point(aes(x = carapace, y = cheliped, colour = year)) + 
  geom_smooth(aes(x = carapace, y = cheliped, colour = year), method = lm)



# Test for type II or III
frair_test(eaten ~ density, data = red) 

### Frair fit
frair_responses()

outII_r <- frair_fit(eaten ~ density, data = red, response = 'rogersII',
                     start = list(a = 0.2, h = 0.2), fixed = list(T=1))

outIII_r <- frair_fit(eaten ~ density, data = red, response = 'flexpnr',
                      start = list(b=0.0029242, h = 0.0312941), fixed = list(T=1, q = 1.9391532))

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
outIII_r_boot <- frair_boot(outIII_r, start = NULL, strata=red[,6], nboot=2000,
                            para=TRUE, ncores=NaN, WARN.ONLY=FALSE)


# is the asymptote ok
confint(outIII_r_boot)
outIII_r_boot

# Illustrate bootlines
plot(outIII_r_boot, xlim=c(0,70), ylim = c(0, 40), type='n', main='All bootstrapped lines')
lines(outIII_r_boot, all_lines=TRUE)
points(outIII_r_boot, pch=20)

# Illustrate bootpolys
plot(outIII_r_boot, xlim=c(0, 65), ylim = c(0, 40), type='n', main='Empirical 95 percent CI',
     xlab = 'Prey Density', ylab = 'Prey Eaten',
     cex.lab = 1.5,
     font.lab = 2,
     cex.axis = 1.2,
     cex.main = 1.5)
drawpoly(outIII_r_boot, border = NA, col=hsv(0/6,0.4, 0.8))
points(outIII_r_boot, pch=20, col=hsv(0/6,0, 0))
lines(outIII_r_boot, lwd = 2, all_lines=FALSE, col=hsv(0/6,0, 0, alpha = 1))


##### nlm to look at asymptote (analysis not in manuscript)
red_asymp <- nls(eaten ~ SSasymp(density, Asym, R0, lrc), 
                 data = red) # asymptotic curve with initial guesses
summary(red_asymp)



# Functional response graph ---------------------------------

plot(outIII_r_boot, xlim=c(0, 65), ylim = c(0, 40), type='n',
     xlab = "Initial Oyster Density",
     ylab="Oysters Consumed", 
     cex.lab = 1.5,
     font.lab = 2,
     cex.axis = 1.2,
     cex.main = 1.5)
lines(outII_g_boot, lwd = 3, all_lines=FALSE, col=hsv(2.5/6,0.9, 0.5, alpha= 1), lty = 2)
lines(outIII_r_boot, lwd = 3, all_lines=FALSE, col=hsv(6/6,0.9, 0.5, alpha = 1), lty = 1)
drawpoly(outIII_r_boot, border = NA, col=hsv(6/6,0.6, 0.8, alpha = 0.4))
drawpoly(outII_g_boot, border = NA, col=hsv(2.5/6,0.6, 0.8, alpha= 0.4))
points(outIII_r_boot, pch=17, col=hsv(6/6,0.8, 0.5), cex = 1.4)
points(outII_g_boot, pch=20,  col=hsv(2.5/6,0.8, 0.5, alpha= 1), cex = 1.4)
legend(x = "topleft", legend = c("Green Crabs", "Red Rock Crabs"), 
       col = c(hsv(2.5/6,0.9, 0.5, alpha= 1), c(hsv(6/6,0.9, 0.5))), 
       lty = c(2, 1), cex = 1.3, pch = c(20, 17))

# Comparing the distributions of the asymptotes of CP to CM ---------------
green.densities <- c(64)

boottest <- outII_g_boot$bootcoefs
greenboot <- green.densities - lambertW0(boottest[, 1] * boottest[, 2] * green.densities * exp(-boottest[, 1] * (1 - boottest[, 2] * green.densities)))/(boottest[, 1] * boottest[, 2]) 

red.densities <- c(64)
boottest <- outIII_r_boot$bootcoefs
redboota <- boottest[, 1] * red.densities^1.939
redboot <- red.densities - lambertW0(redboota * boottest[, 3] * red.densities * exp(-redboota * (1 - boottest[, 3] * red.densities)))/(redboota * boottest[, 3])

bootdisttest <- data.frame(green = greenboot, 
                           red = redboot)

ggplot(data = bootdisttest)+
  geom_histogram(aes(x = green), alpha = 0.3, fill = 'green', colour = 'black', binwidth = 2)+
  geom_histogram(aes(x = red), alpha = 0.3, fill = 'red', colour = 'black', binwidth = 2)+
  labs(x = 'Expected oysters eaten with 64 starting density')

download.packages("fBasics")
library(fBasics)
# Testing if there is a difference between the two asymptotes with kolmogorov-smirnov and t.tests
ks2Test(bootdisttest$green, bootdisttest$red)
t.test(bootdisttest$green, bootdisttest$red)

# Everything I've tried so far results in different means of the two
# species/different distributions that the two species are drawn from. t-test
# suggests difference in mean maximum consumptive rate of ~7-9. I don't think I
# can argue a non-statistically-significant difference without more data

# This is just a test of how the effect of sample size on type 1 error when bootstrapping
simfun <- function(n=6) {
  x <- rnorm(n)
  m.x <- mean(x)
  s.x <- sd(x)
  z <- m.x/(1/sqrt(n))
  t <- m.x/(s.x/sqrt(n))
  b <- replicate(2000, mean(sample(x, replace=TRUE)))
  c( t=abs(t) > qt(0.975,n-1), z=abs(z) > qnorm(0.975),
     z2 = abs(t) > qnorm(0.975), 
     b= (0 < quantile(b, 0.025)) | (0 > quantile(b, 0.975))
  )
}

out <- replicate(2000, simfun())
rowMeans(out)

# n = 3, nboot = 2000
# t      z     z2 b.2.5% 
# 0.0485 0.0500 0.1925 0.2545

# n = 3, nboot = 200
# t      z     z2 b.2.5% 
# 0.060  0.040  0.175  0.270 

# n = 6, nboot = 2000
# t      z     z2 b.2.5% 
# 0.0525 0.0510 0.1125 0.1430 



# Bootstrapping for studentized and BCa intervals -------------------------

CI_green <- confint(outII_g_boot)
CI_red <- confint(outIII_r_boot)
redframe <- data.frame(density = 0:64)
redframe$a_up <- CI_red[['b']]$stud$upper * redframe$density^1.939
redframe$a_lw <- CI_red[['b']]$stud$lower * redframe$density^1.939

0:64 - lambertW0(redboota_up * CI_red[['h']]$stud$upper * 0:64 * exp(-redboota_up * (1 - CI_red[['b']]$stud$upper * 0:64)))/(redboota_up * CI_red[['b']]$stud$upper)


# Illustrate bootpolys
plot(outIII_r_boot, xlim=c(0, 65), ylim = c(0, 60), type='n', xlab = expression(paste(italic('Crassostrea gigas'), ' Density')), ylab = expression(paste('Number of ', italic('Crassostrea gigas'),' Consumed')))
drawpoly(x = 0:64, 
         upper = 0:64 - lambertW0(CI_red * CI_red[['h']]$stud$lower * 0:64 * exp(-CI_red * (1 - CI_red[['h']]$stud$lower * 0:64)))/(CI_red * CI_red[['h']]$stud$lower),
         lower = 0:64 - lambertW0(redboota_lw * CI_red[['h']]$stud$upper * 0:64 * exp(-redboota_lw * (1 - CI_red[['h']]$stud$upper * 0:64)))/(redboota_lw * CI_red[['h']]$stud$upper),
         border = NA, col=hsv(0/6, 0.7, 0, alpha = 0.4))
points(outIII_r_boot, pch=19, col=hsv(0/6, 0.7, 0, alpha = 0.5))
lines(outIII_r_boot, lwd = 1, all_lines=FALSE)

#add green
drawpoly(x = 0:64, 
         upper = 0:64 - lambertW0(CI_green[['a']]$stud$upper * CI_green[['h']]$stud$upper * 0:64 * exp(-CI_green[['a']]$stud$upper * (1 - CI_green[['h']]$stud$upper * 0:64)))/(CI_green[['a']]$stud$upper * CI_green[['h']]$stud$upper),
         lower = 0:64 - lambertW0(CI_green[['a']]$stud$lower * CI_green[['h']]$stud$lower * 0:64 * exp(-CI_green[['a']]$stud$lower * (1 - CI_green[['h']]$stud$lower * 0:64)))/(CI_green[['a']]$stud$lower * CI_green[['h']]$stud$lower),
         border = NA, col=hsv(2/6, 1, 0, alpha = 0.2))
points(outII_g_boot, pch=17,  col=hsv(2/6,1, 0, alpha= 0.5))
lines(outII_g_boot, lwd = 1, all_lines=FALSE, lty = 2)

# add legend
legend(x = "topleft", legend = c(expression(italic("Carcinus maenas")), expression(italic("Cancer productus"))), pch = c(17, 19), 
       lty = c(2, 1), cex = 1, fill = c(hsv(2/6, 1, 0, 0.2), hsv(0/6, 1, 0, 0.4)))


# Number eaten per cm cheliped --------------------------------------------

cheliped <- read.csv('functional_response_cheliped.csv', header = T)
greenchel <- cheliped[cheliped$species == 'green',]
redchel <- cheliped[cheliped$species == 'red',]

par(mfrow = c(1, 1))
plot(cheliped$density_per_cm_round[cheliped$species == 'green'], cheliped$eaten_per_cm[cheliped$species == 'green'], pch = 19, col = 'forestgreen',
     xlab = 'Initial prey density', ylab = 'Number of prey eaten per cm cheliped height')
points(cheliped$density_per_cm_round[cheliped$species == 'red'], cheliped$eaten_per_cm[cheliped$species == 'red'], pch = 17, col = 'firebrick')
legend('topleft', col = c('forestgreen', 'firebrick'), pch = c(19, 17), legend = c('Carcinus maenas', 'Cancer productus'))
