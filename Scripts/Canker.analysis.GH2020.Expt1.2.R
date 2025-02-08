
cankers<-read.csv("../../WA-IN-Transplant2019-20.Necrosis.csv" )
mrphtrt<-read.csv("../../MorphoGH2020Expt1.csv")[,1:6]

names(cankers)
names(mrphtrt)[3] <- "Plant"

library(reshape2)
library(car)
library(lme4)
library(RColorBrewer)
library(ggpubr)
library(tidyverse)
library(agricolae)
library(gridExtra)
library(emmeans)
library(nlme)
library(asbio)
library(MuMIn)
library(dplyr)
library(DHARMa)
#source("Canker.functions.randomcorrected.R")

cank<-plyr::join(cankers,mrphtrt, by=c("State","Block","Plant"))

names(cank)[4]<-"mm2"

cank

# Analyze necrosis
#cank <- cank[cank$mm2 < 200,]
cank$mm2[cank$mm2 == 0] <- 0.00000001

# Make new random variables

cank$Plant.Rand <- with(cank, paste(State,Block,Plant,sep="")) %>% factor()
cank$Treat <- with(cank, paste(State,Block,sep="")) %>% factor()
cank$root.shoot <- cank$Displacement/cank$Caliper

## Analyze without covariate for root.shoot or caliper

	# store data in dataframe to global environment because 
	# boxCox is poorly written and relies heavily on global
	# environment, causing an errror; see the github page
	# https://stackoverflow.com/questions/31921425/apply-and-boxcox-function-error-in-r
	d <<- data.frame(ne=cank$mm2, tr=cank$Treat)

# start by fitting model for box cox
cank.model1 <- lm(ne ~ tr, data=d)

boxCox(cank.model1)
Anova(cank.model1)

## refit with box cox
lambda<- with(boxCox(cank.model1, plotit=FALSE), x[which.max(y)])

# remove d from global environment
remove(d,envir=.GlobalEnv)
lambda

cank.model.bc <- with(cank, lmer(mm2^lambda ~ Treat + (1|Treat:Plant.Rand)))

plot(cank.model.bc)

# perform Type III Anova
Anova(cank.model.bc, 3)

plot(cank.model.bc)

lcpairs<-as.data.frame(pairs(emmeans(cank.model.bc, "Treat"), adjust="tukey"))

lcpairs

plot(cank$root.shoot, cank$mm2)
plot(cank$Displacement, cank$mm2)
plot(cank$Caliper, cank$mm2)

d <<- data.frame(ne=cank$mm2, rsr=cank$root.shoot)
rsr.model1 <- lm(ne ~ rsr, data=d)
lambda<-with(boxCox(rsr.model1, plotit=FALSE), x[which.max(y)])
remove(d,envir=.GlobalEnv)

lambda

cank$Displacement
cank$Caliper

summary(with(cank, lm(mm2^.3 ~ Displacement)))
summary(with(cank, lm(mm2^.3 ~ Caliper)))
summary(with(cank, lm(mm2^.3 ~ Caliper+Displacement+root.shoot)))
summary(with(cank, lm(mm2^.3 ~ Displacement+root.shoot)))
summary(with(cank, lm(mm2^.3 ~ Caliper+root.shoot)))

with(cank, lmer(Displacement ~ State + (1|State:Block))) %>% summary
with(cank, lmer(Displacement ~ State+ (1|State:Block))) %>% Anova(type=3)
with(cank, lmer(Displacement ~ State+ (1|State:Block))) %>% simulateResiduals() %>% plot
with(cank, lm(Displacement ~ State)) %>% leveneTest()

with(cank, lmer(Caliper ~ State+ (1|State:Block))) %>% summary
with(cank, lmer(Caliper ~ State+ (1|State:Block))) %>% Anova(type=3)
with(cank, lmer(Caliper ~ State+ (1|State:Block))) %>% simulateResiduals() %>% plot

with(cank, lmer(root.shoot ~ State+ (1|State:Block))) %>% summary
with(cank, lmer(root.shoot ~ State+ (1|State:Block))) %>% Anova(type=3)
with(cank, lmer(root.shoot ~ State+ (1|State:Block))) %>% simulateResiduals() %>% plot

rsr.model.bc <- with(cank, lm(mm2^.3 ~ root.shoot))
str(rsr.model.bc)
summary(rsr.model.bc)

## Analyze with covariate for root.shoot
cank.model.bc.cov <- with(cank, lmer(mm2^lambda ~ Treat + root.shoot + (1|Treat:Plant.Rand)))

Anova(cank.model.bc.cov, 3)

lcpairs.cov<-as.data.frame(pairs(emmeans(cank.model.bc.cov, "Treat"), adjust="tukey"))

lcpairs.cov

levels(cank$Treat) <-
  c(
    "IN P1",
    "IN P2",
    "WA P1",
    "WA P2",
    "WA P3")

# look at the data
cank %>%
  filter(Inoculation =='Gm') %>%
  #with(.,boxplot(mm2 ~ Treat))
  ggplot()+
  geom_boxplot(mapping=aes(x = mm2, y = Treat))+
  geom_jitter(mapping=aes(x = mm2, y = Treat))

pn0<-
  cank %>%
    filter(Inoculation =='Gm') %>%
    #filter(mm2 < 200) %>%
    #filter(mm2 > 1) %>%
    with(plot.necrosis(treat=Treat,
                       cov=root.shoot,
                       random=Plant,
                       necr=mm2,
                       cook.denom=4,
                       #adjust="none",
                       return.analysis = T))

# don't get rid of any outliers
pn1 <-cank %>%
  filter(Inoculation =='Gm') %>%
  #.[-pn0$outliers,]%>%
  #filter(mm2 < 200) %>%
  #filter(mm2 > 1) %>%
  with(plot.necrosis(treat=Treat,
                     cov=root.shoot,
                     random=Plant,
                     necr=mm2,
                     cook.denom=NULL,
                     adjust="none",
                     return.analysis = T))

shapiro.test(residuals(pn1$Regression))

for (i in c(30:40)/100) {
  print(i)
  cank %>%
    filter(Inoculation =='Gm') %>%
    mutate(trans_mm2 = mm2^i) %>%
    with(.,
         lmer(trans_mm2 ~ Treat + root.shoot + (1|Treat:Plant))) %>%
    residuals %>%
    shapiro.test %>% print
  
}

# .36 was the max

cank %>%
  filter(Inoculation =='Gm') %>%
  mutate(trans_mm2 = mm2^.36) %>%
  with(.,
       lmer(trans_mm2 ~ State + root.shoot + (1|State:Treat:Plant))) %>% Anova

windows(); cank %>%
  filter(Inoculation =='Gm') %>%
  mutate(trans_mm2 = mm2^.36) %>%
  with(.,
       lmer(trans_mm2 ~ State + State:Treat + root.shoot + (1|Treat:Plant)))  %>% Anova(type=3)#simulateResiduals %>% plot

windows(); cank %>%
  filter(Inoculation =='Gm') %>%
  mutate(trans_mm2 = mm2^.36) %>%
  with(.,
       lmer(trans_mm2 ~ State + State:Treat + root.shoot + (1|Treat:Plant)))  %>% emmeans(~State | root.shoot, adjust='none')

windows(); cank %>%
  filter(Inoculation =='Gm') %>%
  mutate(trans_mm2 = mm2^.36) %>%
  with(.,
       lmer(trans_mm2 ~ State + State:Treat + root.shoot + (1|Treat:Plant)))  %>% emmeans(~State | root.shoot, adjust='none') %>% pairs(adjust='none')

(4.05 ^(1/.36)-3.55^(1/.36))/4.05 ^(1/.36)

((4.05 +0.196)^(1/.36)-(3.55-0.196)^(1/.36))/4.05 ^(1/.36)
((4.05 -0.196)^(1/.36)-(3.55+0.196)^(1/.36))/4.05 ^(1/.36)

(4.05 ^(1/.36)-3.55^(1/.36))/4.05 ^(1/.36)-((4.05 +0.196)^(1/.36)-(3.55-0.196)^(1/.36))/4.05 ^(1/.36)
(4.05 ^(1/.36)-3.55^(1/.36))/4.05 ^(1/.36)-((4.05 -0.196)^(1/.36)-(3.55+0.196)^(1/.36))/4.05 ^(1/.36)

anova(
  cank %>%
    filter(Inoculation =='Gm') %>%
    mutate(trans_mm2 = mm2^.36) %>%
    with(.,
         lmer(trans_mm2 ~ State:Treat + root.shoot + (1|Treat:Plant))) ,
  cank %>%
    filter(Inoculation =='Gm') %>%
    mutate(trans_mm2 = mm2^.36) %>%
    with(.,
         lmer(trans_mm2 ~ State + State:Treat + root.shoot + (1|Treat:Plant)))
)


# test State
# cank %>%
#  filter(Inoculation =='Gm') %>%
#  mutate(trans_mm2 = mm2^.36) %>%
#  with(
#    anova(
#      lmer(trans_mm2 ~ State + State:Treat + root.shoot + (1|Treat:Plant)),
#      lmer(trans_mm2 ~ State:Treat + root.shoot + (1|Treat:Plant))
#    )
#  )


pn_interaction <-
  cank %>%
  filter(Inoculation =='Gm') %>%
  #.[-pn0$outliers,]%>%
  #filter(mm2 < 200) %>%
  #filter(mm2 > 1) %>%
  mutate(treat = Treat,
         cov = root.shoot,
         random = Plant,
         n = mm2^.36) %>%
  with(., 
       lmer(n ~
              treat+ # main effect for treatment
              cov +
              cov:treat +    # effect of covariate depends on treatment
              (1|treat:random)   # random effect of plant
       )) %>%
  (function (x) list(Regression=x,
                     Summary=summary(x),
                     Anova=Anova(x, type=3),
                     Tukey=emmeans(x, ~ treat | cov)%>%
                       pairs()%>%
                       #pairs(adjust='none')%>%
                       as.data.frame)) %>%
  (function (x) {x[['Letters']]<- cldList(p.value ~ contrast, data = x[[4]]);
  x})

r.squaredLR(pn_interaction$Regression) # conditional R2
performance::r2(pn_interaction$Regression)
performance::r2_nakagawa(pn_interaction$Regression)


pn_interaction$Summary
pn_interaction$Anova


windows();plot(simulateResiduals(fittedModel = pn_interaction$Regression))

#cank %>%
#  filter(Inoculation =='Gm') %>%
##       lm(mm2 ~ State * root.shoot, data=.)%>%Anova
#  ggplot()+geom_point(aes(y=mm2, x=root.shoot))+facet_wrap('Treat')

pn <- cank %>%
  filter(Inoculation =='Gm') %>%
  #filter(mm2 < 200) %>%
  #filter(mm2 > 10) %>%
  mutate(treat = Treat,
         cov = root.shoot,
         random = Plant,
         n = mm2^.36) %>%
  with(.,
       lmer(n ~ treat + cov + (1|treat:random)))%>%
  (function (x) list(Regression=x,
                     Summary=summary(x),
                     Anova=Anova(x, type=3),
                     Tukey=emmeans(x, 'treat')%>%
                       pairs()%>%
                       #pairs(adjust='none')%>%
                       as.data.frame)) %>%
  (function (x) {x[['Letters']]<- cldList(p.value ~ contrast, data = x[[4]]);
  x})

performance::r2(pn$Regression)

windows();plot(simulateResiduals(fittedModel = pn$Regression))

#  cank %>%
#  filter(Inoculation =='Gm') %>%
##  .[-pn0$outliers,]%>%
#  filter(mm2 < 200) %>%
#  filter(mm2 > 1) %>%
#  with(.,
#       plot.necrosis(treat=Treat,
#                     cov=root.shoot,
#                     random=Plant,
#                     necr=mm2,
#                     cook.denom=NULL,
#                     adjust="tukey",
#                     xlab="Soil Plot",
#                     return.analysis = T))
pn
#pn$lambda
#pn$Regression
#pn$Regression@beta
#pn@beta

pn_interaction$Anova
pn$Anova
pn_interaction$Tukey
emmeans(pn_interaction$Regression, 'treat')

(4.64^(1/.36)-3.06^(1/.36))/4.64^(1/.36)

((4.64+0.316)^(1/.36)-(3.06-.274)^(1/.36))/4.64^(1/.36)
((4.64-0.316)^(1/.36)-(3.06+.274)^(1/.36))/4.64^(1/.36)

(4.64^(1/.36)-3.06^(1/.36))/4.64^(1/.36)-((4.64+0.316)^(1/.36)-(3.06-.274)^(1/.36))/4.64^(1/.36)
(4.64^(1/.36)-3.06^(1/.36))/4.64^(1/.36)-((4.64-0.316)^(1/.36)-(3.06+.274)^(1/.36))/4.64^(1/.36)

(4.18^(1/.36)-3.06^(1/.36))/4.64^(1/.36)

(4.18^(1/.36)-3.06^(1/.36))/4.18^(1/.36)-((4.18+0.310)^(1/.36)-(3.06-.274)^(1/.36))/4.18^(1/.36)
(4.18^(1/.36)-3.06^(1/.36))/4.18^(1/.36)-((4.18-0.310)^(1/.36)-(3.06+.274)^(1/.36))/4.18^(1/.36)


pn$Tukey
pn_interaction$Letters
pn$Letters
AIC(pn_interaction$Regression)
AIC(pn$Regression)

anova(pn_interaction$Regression,
      pn$Regression)

parameters::model_parameters(pn$Regression)
cank_rootshoot_int_params <- parameters::model_parameters(pn_interaction$Regression)
cank_rootshoot_int_params %>%
  as.data.frame %>%
  write.csv("cank_rootshoot_int_params.csv", row.names=F)

# windows()
# cank %>%
#  filter(Inoculation =='Gm') %>%
#  mutate(plant_size = scale(Displacement) + scale(Caliper))%>%
#  with(summary(lm(mm2 ~ plant_size*State)))
#  #with(plot(simulateResiduals(lm(mm2 ~ plant_size*State))))
#  #ggplot()+geom_point(aes(y=mm2, x=plant_size))+facet_wrap("State")
  
gp <- ggeffects::predict_response(
  #ggeffects::ggpredict(
  #pn$Regression,
  pn_interaction$Regression,
#  newdata =  cank %>%
#    filter(Inoculation =='Gm') %>%
#    mutate(treat = Treat,
#           cov = root.shoot,
#           random = Plant,
#           n = mm2^.36) %>% na.omit,
  terms = c('cov','treat'),
  # margin='empirical',
  # typical='weighted.mean',
  ci_level = .48)

#plot(gp)+facet_wrap('group')
attr(gp,'terms')[2] <- 'Treat'
attr(gp,'original.terms')[2] <- 'Treat'
names(attr(gp,'at.list'))[2] <- 'Treat'
attr(gp,'legend.title') <- 'State & Block'

levels(gp$group) <-
  levels(cank$Treat) <-
  c(
    "IN P1",
    "IN P2",
    "WA P1",
    "WA P2",
    "WA P3")

#windows(5,4);

gppoints <- cank %>%
  filter(Inoculation =='Gm') %>%
  #filter(mm2 < 200) %>%
  #filter(mm2 > 1) %>% 
  mutate(group_col = Treat) %>%
  mutate(group = Treat) %>%
  (function (y) inner_join(
    y, summarise(y,
                 mm2.mean = mean(mm2),
                 mm2.se   = sd(mm2)/sqrt(2),
                 .by='Plant.Rand')
  ) %>% dplyr::select(-mm2) %>% distinct %>% na.omit)

gp.df <- as.data.frame(gp) %>% left_join(gppoints %>% summarise(minx=min(root.shoot), maxx=max(root.shoot), .by=Treat)%>% mutate(group = Treat))


#gp$conf.low[gp$conf.high<=0] <- 0
gp$conf.high[gp$conf.high>=130^.36] <- 130^.36

gp$predicted[which(gp.df$x>=gp.df$maxx+5)] <- NA#130^.3
gp$conf.low[which(gp.df$x>=gp.df$maxx+5)] <- NA#130^.3
gp$conf.high[which(gp.df$x>=gp.df$maxx+5)] <- NA#130^.3

gp$predicted[which(gp.df$x<=gp.df$minx-5)] <- NA#130^.3
gp$conf.low[which(gp.df$x<=gp.df$minx-5)] <- NA#130^.3
gp$conf.high[which(gp.df$x<=gp.df$minx-5)] <- NA#130^.3

interactive_model_effects_plot <-
  plot(gp) +
  #plot(gp %>% left_join(dplyr::select(cank, State, Block, group = Treat)))+
  #ggplot()+
    geom_point(
      data = gppoints,
      mapping = aes(x = root.shoot,
                    y = mm2.mean^.36,
                    color=group_col))+
    geom_errorbar(
      data = gppoints,
      mapping = aes(x = root.shoot,
                    ymin = (mm2.mean-mm2.se)^.36,
                    ymax = (mm2.mean+mm2.se)^.36,
                    color=group_col), inherit.aes='false')+
    scale_y_continuous(#name = 'n',
                       transform = scales::transform_boxcox(1/.36),
                       breaks = seq(from=0, to=150, by=25)^.36,
                       labels = seq(from=0, to=150, by=25),
                       limits = c(0,130^.36)
                       )+
    labs(x="\nRoot to shoot ratio",
         y="Canker necrotic area (mm\u00b2)\n",
         title='')+
  facet_wrap('group')+
#    facet_grid(State ~ Block, scales='free', axes='all',margins=FALSE)+
    theme_minimal() +
    theme(legend.position = 'none',
          plot.title = element_blank())

#svg(Fig)
windows(5,4)

svg("Revision_Jan_2025/Figures_ordered/FigS4_interactive_model_effects.svg", 5, 4)
interactive_model_effects_plot; dev.off()

devEMF::emf("Revision_Jan_2025/Figures_ordered/FigS4_interactive_model_effects.emf", 5, 4)
interactive_model_effects_plot; dev.off()

windows(3.5,3.5);cank %>%
  filter(Inoculation =='Gm') %>%
  #filter(mm2 < 200) %>%
  #filter(mm2 > 1) %>%
#  mutate(mm2adj = pn$Regression@beta[1]+
#           as.integer(State=='IN' & Block=="SW")*pn$Regression@beta[2]+
#           as.integer(State=='WA' & Block=="N")*pn$Regression@beta[3]+
#           as.integer(State=='WA' & Block=="SE")*pn$Regression@beta[4]+
#           as.integer(State=='IN' & Block=="SW")*pn$Regression@beta[5]) %>%
#  mutate(mm2lambda=mm2^.3) %>%
#  dplyr::select(root.shoot, Displacement, Caliper,mm2lambda)%>%na.omit%>%
#  pca %>% plot(scaling=1,type='text')
  na.omit %>%
  GGally::ggpairs(
    columns=c('root.shoot', 'Displacement', 'Caliper'))#,'mm2lambda'))

  
######

plot(simulateResiduals(fittedModel = pn$Regression))

testDispersion(pn$Regression)

testUniformity(pn$Regression)
testOutliers(pn$Regression)


plot(pn$Regression, which=1)
qqnorm(residuals(pn$Regression))
shapiro.test(residuals(pn$Regression))
hist(residuals(pn$Regression))
