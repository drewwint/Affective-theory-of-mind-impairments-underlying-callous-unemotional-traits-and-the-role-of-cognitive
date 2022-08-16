
title: "Affective theory of mind impairments underlying callous-unemotional traits and the role of cognitive control: A pilot study"
author: "Drew E. Winters"


# Loading packages


library(lmerTest) ## for MLM with p values
library(ggplot2)  ## plotting
library(forcats)  ## to adjust for plot values
library(dplyr)    ## data manipulation 
library(lmtest)   ## for likelihood ratio tests
library(car)      ## re-coding and Durban Watson test
library(lavaan)   ## ML estimation 
library(olsrr)    ## Multicollinearity
library(MKinfer)  ## Bootstrapped T tests



# Demographics


## Age
psych::describe(T2_dat$age)


## Sex
a<-rbind(table(T2_dat$sex_male),round((table(T2_dat$sex_male)/sum(table(T2_dat$sex_male))*100),0))
colnames(a)<-c('female','male')
rownames(a)<-c('n','%')
a


## High to normative CU

a<-rbind(table(T2_dat$cu_group_9split),round((table(T2_dat$cu_group_9split)/sum(table(T2_dat$cu_group_9split))*100),0))
colnames(a)<-c('low_CU','hi_CU')
rownames(a)<-c('n','%')
a

## calculating clinical cutoffs

### CU traits 
clin<-ifelse(ifelse(T2_dat$sex_male==0 & T2_dat$age<=14 & T2_dat$icu_total >= 31, 1,0)+ifelse(T2_dat$sex_male==0 & T2_dat$age>=15 & T2_dat$icu_total >= 34, 1,0) + ifelse(T2_dat$sex_male==1 & T2_dat$age<=14 & T2_dat$icu_total >= 36, 1,0)+ifelse(T2_dat$sex_male==1 & T2_dat$age>=15 & T2_dat$icu_total >= 40, 1,0),1,0)


clindat<-data.frame(rbind(as.character(table(clin)),round(((table(clin)/NROW(clin))*100),2)))
colnames(clindat)<- c("non-clinical", "clinical")
rownames(clindat)<- c("n", "%")
clindat



### SDQ CD

sum(T2_dat$sdq_cd)
sum(T2_dat$sdq_cd)/nrow(T2_dat)*100



## Race & Ethnicity
### Multiple races
over_1_race<-which((T2_dat$race_Asian)+(T2_dat$race_Indian)+(T2_dat$race_Black)+
  (T2_dat$race_PacificIsland)+(T2_dat$race_White)+
  (T2_dat$race_other)>1)
    ### # of participants who selected more than one race 
length(over_1_race)
round(length(over_1_race)/NROW(T2_dat)*100,2)


### Describing all races
race<-(1*T2_dat$race_Asian)+(2*T2_dat$race_Indian)+(3*T2_dat$race_Black)+
  (5*T2_dat$race_PacificIsland)+(6*T2_dat$race_White)+
  (7*T2_dat$race_other)

    ## race # and proportion
a<-rbind(as.character(table(race)),round((table(race)/sum(table(race))*100),2))
rownames(a)<-c("n","%")
a


### Ethnicity 
a<-rbind(as.character(table(T2_dat$race_Latino)),round((table(T2_dat$race_Latino)/sum(table(T2_dat$race_Latino))*100),3))
colnames(a) <-c("non-latinx", "latinx")
rownames(a) <- c("n","%")
a

## Descriptives table 
### continuous descriptives
psych::describe(T2_dat %>% select("icu_total",
                                  "sdq_total_cd",
                                  "load2_max", 
                                  "mie1_c", 
                                  "mie2_c","age"))



### dichotomous descriptives
data.frame(rbind(as.character(table(T2_dat$bad_resp_2)), round(table(T2_dat$bad_resp_2)/(sum(table(T2_dat$bad_resp_2))/100),3)), row.names = c("bad_resp n","bad_resp %"))

data.frame(rbind(as.character(table(T2_dat$sex_male)), round(table(T2_dat$sex_male)/(sum(table(T2_dat$sex_male))/100),3)), row.names = c("sex_male n","sex_male %"))

data.frame(rbind(as.character(table(T2_dat$race_White)), round(table(T2_dat$race_White)/(sum(table(T2_dat$race_White))/100),3)), row.names = c("race_White n","race_White %"))


### correlations
data.frame(cor(T2_dat %>% select("icu_total",
                                 "sdq_total_cd",
                                 "load2_max", 
                                 "mie1_c", 
                                 "mie2_c",
                                 "age", 
                                 "sex_male",
                                 "race_White")))

PerformanceAnalytics::chart.Correlation(T2_dat %>% select("icu_total",
                                                          "sdq_total_cd",
                                                          "load2_max", 
                                                          "mie1_c", 
                                                          "mie2_c",
                                                          "age", 
                                                          "sex_male",
                                                          "race_White"))

data.frame(cor(T2_dat %>% select("icu_total",
                                 "sdq_total_cd",
                                 "load2_max", 
                                 "mie1_c", 
                                 "mie2_c",
                                 "age", 
                                 "sex_male",
                                 "race_White"), method= "spearman"))



# Preliminary analyses
## Assessing bias for those that are removed - T-tests
set.seed(1982)
lapply(T2_dat[,c("mie2_complex", 
                 "mie1_complex",
                 "icu_total" ,
                 "sex_male", 
                 "age", 
                 "race_White", 
                 "anny2_diff",
                 "WAI_total",
                 "sdq_total_cd")], 
       function(x) MKinfer::boot.t.test(x ~ T2_dat$bad_load_2, var.equal = FALSE))




## Assumption testing 

### Setting up dataframe for assumption testing 
T2_dat_test1 <- T2_dat[,c("record_id",
                          "mie2_complex", 
                          "mie1_complex",
                          "icu_total" ,
                          "sex_male", 
                          "age", 
                          "race_White", 
                          "anny2_diff",
                          "WAI_total",
                          "sdq_total_cd")]

### *Additivity - Multicollinerity*


corrt <- cor(T2_dat_test1[,-c(1)])
symnum(corrt) 

  # THere are no '+', '*', or 'B' suggesting a high correlation betwen variables so these are fine. 

 
### *Normality of residuals*
random <- rchisq(nrow(T2_dat_test1[,-c(1)]), 10)
  # creating a fake regression 
    # remove the ID variable since you are regression the entire data set
fake <- lm(random~.,data=T2_dat_test1[,-c(1)])
  # creating studentized residuals
studentized <- rstudent(fake)


### *Distribution of residuals*
hist(studentized,breaks = 8,
     main="Histogram of Residuals", 
     xlab="Studentized Residuals", 
     border="black", 
     col="blue")

  ## Mostly evenly distributed just a alight positive tail -- let's check the skewnes and kurtosis



psych::describe(studentized)

    # This looks okay - only a slight positive tail but the distribution is relativley intact. No significant concerns.

### *Mean of residuals close to 0*


mean(fake$residuals) ## mean is close to 0

#   Mean is very close to 0


### *Residual plot*
  # resigual plot
    # creating a z-scored fitted variable vector for interpretation 
fitted <- scale(fake$fitted.values)
  # residual plot - 
plot(fitted, studentized,
     main="Residual Plot of Standardized Fitted Values", 
     xlab="Fitted", 
     ylab = "Studentized",
     col="blue")
abline(0,0, col="red")
abline(v=0, col="red") 

    # ata look evenly dispursed across the figure in each sector

### *Autocorrelation* 
durbinWatsonTest(fake) # autocorrelation of residuals
  # Test is not statistically significant (p>.05) and the statistic is close to 2

### *Multicollinearity of predictors and controls*
library(olsrr)
ols_vif_tol(assumpt<-lm(random ~ mie2_complex + mie1_complex + icu_total + sex_male + age + race_White + anny2_diff + WAI_total, data = T2_dat))


  #  After including all variables that will be used as controls and predictors a only a small portion of each variable is filtered away meaning we do not have multicollinearity assumption met.- all VIF < 2.5 all tolerance is >=.5



# Measures Reliability Alpha
## ICU
psych::alpha(T2_dat[,c(grepl( "ICU." , names( T2_dat ) ))])


psych::describe(T2_dat$icu_total)


## SDQ - CD

psych::alpha(T2_dat[,c(grepl( "sdq\\d+" , names( T2_dat ) ))])



psych::describe(T2_dat$sdq_total_cd)




b<-psych::omega(mie_it1[-c(which(T2_dat$bad_load_2==1)),], 
                poly=TRUE,plot=FALSE,check.keys=TRUE)$omega.tot
b1<-psych::omega(mie_it1[-c(which(T2_dat$bad_load_2==1)),
                         c(1,2,4,6,7,17,19,21,25,26,28)], 
                 poly=TRUE,plot=FALSE,check.keys=TRUE)$omega.tot
b2<-psych::omega(mie_it1[-c(which(T2_dat$bad_load_2==1)),
                         c(3,5,8,9,10,11,12,13,14,15,16,18,20,22,23,24,27)], 
                 poly=TRUE,plot=FALSE,check.keys=TRUE)$omega.tot

b3<-psych::omega(mie_it2[-c(which(T2_dat$bad_load_2==1)),], 
                 poly=TRUE,plot=FALSE,check.keys=TRUE)$omega.tot
b4<-psych::omega(mie_it2[-c(which(T2_dat$bad_load_2==1)),
                         c(1,2,4,6,7,17,19,21,25,26,28)], 
                 poly=TRUE,plot=FALSE,check.keys=TRUE)$omega.tot
b5<-psych::omega(mie_it2[-c(which(T2_dat$bad_load_2==1)),
                         c(3,5,8,9,10,11,12,13,14,15,16,18,20,22,23,24,27)], 
                 poly=TRUE,plot=FALSE,check.keys=TRUE)$omega.tot


cbind(rbind(data.frame(T1_Total=b,row.names = "Omega - ToM"),
            data.frame(T1_Total=round(b,2),row.names = "Omega - ToM rounded")),
      rbind(data.frame(T1_Basic=b1,row.names = "Omega - ToM"),
            data.frame(T1_Basic=round(b1,2),row.names = "Omega - ToM rounded")),
      rbind(data.frame(T1_Complex=b2,row.names = "Omega - ToM"),
            data.frame(T1_Complex=round(b2,2),row.names = "Omega - ToM rounded")),
      rbind(data.frame(T2_Total=b3,row.names = "Omega - ToM"),
            data.frame(T2_Total=round(b3,2),row.names = "Omega - ToM rounded")),
      rbind(data.frame(T2_Basic=b4,row.names = "Omega - ToM"),
            data.frame(T2_Basic=round(b4,2),row.names = "Omega - ToM rounded")),
      rbind(data.frame(T2_Complex=b5,row.names = "Omega - ToM"),
            data.frame(T2_Complex=round(b5,2),row.names = "Omega - ToM rounded")))




# Mind in the eyes performance

## For overall ToM

library(lavaan)
mie<-'
mie1_c~icu_total+ 
        sex_male + 
        age + 
        race_White + 
        sdq_total_cd + 
        bad_resp2
'
mie_sem<-sem(mie,estimator='ML',data=T2_dat) 
summary(mie_sem, rsquare=TRUE,standardized=TRUE, nd=4)


### Bootstrapped CIs
set.seed(1982)
mie_semb<-sem(mie,estimator='ML',data=T2_dat[-c(T2_dat$bad_resp2),], se="bootstrap", bootstrap=5000)
parameterEstimates(mie_semb, boot.ci.type = "bca.simple")


## Basic and Complex

library(lavaan)
mie<-'
mie1_complex + mie1_basic ~ icu_total + 
                            sex_male + 
                            age + 
                            race_White  + 
                            sdq_total_cd + 
                            bad_resp2
'
mie_sem<-sem(mie,estimator='ML',data=T2_dat) 
parameterestimates(mie_sem, rsquare=TRUE,standardized=TRUE)[37:38,]



### Bootstrapped CIs
set.seed(1982)
mie_semb<-sem(mie,estimator='ML',data=T2_dat, se="bootstrap", bootstrap=5000)
parameterEstimates(mie_semb,rsquare=TRUE,standardized=TRUE, boot.ci.type = "bca.simple")[1:12,]


## Subscale analyses
> here we examine if any single subscale is driving the association observed.


library(lavaan)
mie2<-'
mie1_complex + mie1_basic ~ icu_callousness + 
                            icu_unemotional + 
                            icu_uncaring + 
                            sex_male + 
                            age + 
                            race_White  + 
                            sdq_total_cd + 
                            bad_resp2
'
mie_sem2<-sem(mie2,estimator='ML',data=T2_dat) 
parameterestimates(mie_sem2, rsquare=TRUE,standardized=TRUE)[c(1:16,56:57),]



## Figure



## still editing this one - need to get the labels right 

icu<-resid(lm(icu_total ~ 
            sex_male + 
            age + 
            race_White + 
            sdq_total_cd + 
            bad_resp2, data=T2_dat))

comp<-resid(lm(mie1_complex ~ 
            sex_male + 
            age + 
            race_White + 
            sdq_total_cd + 
            bad_resp2, data=T2_dat))

basic<-resid(lm(mie1_basic ~ 
            sex_male + 
            age + 
            race_White + 
            sdq_total_cd + 
            bad_resp2, data=T2_dat))

dat<-data.frame(icu,comp,basic)

library(latex2exp)
library(ggpubr)
my_y_titlea <- expression(paste(bold(underline("Complex")), " Affective Theory of Mind - T1"))
my_y_titleb <- expression(paste(bold(underline("Basic")), " Affective Theory of Mind - T1"))
a<-ggplot(dat,aes(x=icu,y=comp)) +
  geom_smooth(method="lm",col="black") +
  geom_point(color="dark grey")+ 
  theme_minimal()+
  scale_x_continuous("Callous-Unemotional Traits",n.breaks = 5, labels = c("20","25","30","35","40", "45", "50")) + #,n.breaks = 5, labels = c("20","25","30","35","40", "45", "50")
  scale_y_continuous(my_y_titlea, limits = c(-4,3), n.breaks = 4, labels = c("7.5","10","12.5","15","17.5")) + #, n.breaks = 4, labels = c("7.5","10","12.5","15","17.5")
  labs(title = TeX("$\\;$  $\\textit{R^2} = $\\textit{0.166^*}"))+
  theme(axis.text=element_text(size=rel(1.1)),axis.title=element_text(size=17),title =element_text(size=17))
#+ylim(-4,3)
# + scale_y_continuous(my_y_titlea,n.breaks =  5, labels = c("200","600","1000","1400","1800")

b<-ggplot(dat,aes(x=icu,y=basic)) +
  geom_smooth(method="lm",col="black") +
  geom_point(color="dark grey")+ 
  theme_minimal()+ 
  scale_x_continuous("Callous-Unemotional Traits",n.breaks = 5, labels = c("20","25","30","35","40", "45", "50")) + # n.breaks = 5, labels = c("20","25","30","35","40", "45", "50")
  scale_y_continuous(my_y_titleb,limits = c(-4,3),n.breaks = 4, labels = c("1","3","5","7","9")) + #,n.breaks = 4, labels = c("1","3","5","7","9")
  labs(title = TeX("$\\;$ $\\textit{R^2} = $\\textit{0.048^n^s}"))+
  theme(axis.text=element_text(size=rel(1.1)),axis.title=element_text(size=17),title =element_text(size=17))



d<-ggarrange(b,a, ncol=2,nrow=1)

d
## 550 x 550 or 5.5 in "in" units
#ggsave(path = "C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\DPRG main project\\Analyses and Data\\full_data_87\\1_Paper-1", plot=d, width = 5.5, height = 5.5, filename = "basic-comp_aff.tiff", device='tiff', dpi=700, limitsize = FALSE)





# Performance on Stroop signal task


sst<-'
load1_max ~ icu_total + 
            sex_male + 
            age + 
            race_White + 
            sdq_total_cd + 
            bad_resp2
'
sst_sem<-sem(sst,data=T2_dat[-c(which(T2_dat$bad_load_2==1)),])
parameterestimates(sst_sem,rsquare=TRUE,standardized=TRUE)[29,]



## Bootstrapped CIs
set.seed(198295)
sst_semb<-sem(sst,data=T2_dat[-c(which(T2_dat$bad_load_2==1)),], se="bootstrap", bootstrap=5000)
parameterEstimates(sst_semb,rsquare=TRUE,standardized=TRUE, boot.ci.type = "bca.simple")[c(1:6,29),]



## Subscale analyses 
> here we examine if any single subscale is driving the association observed.


sst<-'
load1_max ~ icu_unemotional + 
            icu_callousness + 
            icu_uncaring + 
            sex_male + 
            age + 
            race_White + 
            sdq_total_cd + 
            bad_resp2
'
sst_sem<-sem(sst,data=T2_dat[-c(which(T2_dat$bad_load_2==1)),])
parameterestimates(sst_sem,rsquare=TRUE,standardized=TRUE)



## Figures
icu<-resid(lm(icu_total ~ 
            sex_male + 
            age + 
            race_White + 
            sdq_total_cd + 
            bad_resp2, data=T2_dat))

load<-resid(lm(load1_max ~ 
            sex_male + 
            age + 
            race_White + 
            sdq_total_cd + 
            bad_resp2, data=T2_dat))
dat<-data.frame(icu,load)


library(latex2exp)
ggplot(dat,aes(x=icu,y=load)) +
  geom_smooth(method="lm",col="black",formula ='y ~ x') +
  geom_point(color="dark grey")+ 
  theme_minimal()+
  scale_y_continuous("Max Cognitive Load",n.breaks =  6, labels = c("200","600","1000","1400","1800")) +# ,n.breaks =  6, labels = c("200","600","1000","1400","1800")
  scale_x_continuous("Callous-Unemotional Traits",n.breaks = 5, labels = c("20", "25", "30", "35","40", "45", "4")) + #,n.breaks = 5, labels = c("20", "25", "30", "35","40", "45", "4")
  labs(title = TeX("$\\;$ $\\;$ $\\;$ $\\;$ $\\;$ $\\;$ $\\;$ $\\;$ $\\textit{R^2} = $\\textit{0.162^*}"))+
  theme(axis.text=element_text(size=rel(1.1)),axis.title=element_text(size=17),title =element_text(size=17))

## 450x450 or 4.5 x 4.5 in "in" units
#ggsave(path = "C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\DPRG main project\\Analyses and Data\\full_data_87\\1_Paper-1", width = 4.5, height = 4.5, filename = "max_load2.tiff", device='tiff', dpi=700, limitsize = FALSE)



>  Suggests those higher in CU traits had a more difficult time sustaining their attention for stop signal and pressed button too soon


### Counts of max load by CU
ggplot(T2_dat,aes(x=icu_total, y= load2_max))+ 
  geom_col(width = 2) + 
  scale_y_continuous("Max Cognitive Load", 
                     n.breaks = 8, 
                     labels = c("0","300","600","900","1200","1500","1800")) + 
  xlab("Callous-Unemotional Traits") + 
  xlim(18,49) + 
  theme_minimal()

#ggsave(path = "C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\DPRG main project\\Analyses and Data\\full_data_87\\1_Paper-1", width = 3, height = 3, filename = "max_load_dist_CU.tiff", device='tiff', dpi=700, limitsize = FALSE)


#### count stats

max_count <- T2_dat %>% group_by( load2_max) %>% summarize(count = n())

max_count[1:21,]

# maximum and min count for duration of CL
data.frame(rbind(max_count=max_count[which(max_count$count==max(max_count$count)),],min_count=max_count[which(max_count$count==min(max_count$count)),]))







# Change in mind in the eyes task - Identifying complex emotions

## Multilevel model comparing mean differneces form T1-T2 
> doesnt matter if I control for CD or not


### for ICU total




T2_dat_long2 <- tidyr::gather(T2_dat[-c(which(T2_dat$bad_load_2==1)),][,which(colnames(T2_dat) %in% c("record_id", "mie2_complex", "mie1_complex","icu_unemotional" ,"icu_total" , "bad_resp2", "sex_male", "age", "race_White", "anny2_diff","WAI_total","sdq_total_cd"))],time,value,mie1_complex,mie2_complex,mie1_complex:mie2_complex, factor_key=TRUE)

T2_dat_long2$time_rc <- as.numeric(car::recode(T2_dat_long2$time, "'mie2_complex'=1;'mie1_complex'=0"))


T2_dat_long2<-semTools::orthogonalize(T2_dat_long2, var1=c("icu_total", "time_rc"), var2=c("time_rc","icu_total"), match = TRUE)

llmod3<-lmer(value~ icu_total + 
               time_rc + 
               icu_total.time_rc + 
               anny2_diff + 
               sex_male + 
               age + 
               sdq_total_cd + 
               race_White +
               bad_resp2 +
               (1|record_id),
             data=T2_dat_long2)


reg2<- lm(value~icu_total* time + 
            anny2_diff + 
            sex_male + 
            age + 
            race_White+
            sdq_total_cd +
            bad_resp2,
          data=T2_dat_long2)


lrtest(reg2,llmod3)



> The addition of a random effect did improve model fit to the data - plausibly due to the random effect accounting for cluster confounding


anova(llmod3)
summary(llmod3)



#### R2 for model
round(MuMIn::r.squaredGLMM(llmod3),3)



#### R2 for each IV
r2glmm::r2beta(llmod3)



#### Confidence intervals

confint(llmod3)



##### Bootstrapped CIs
set.seed(101)
gg<- lmeresampler::bootstrap(llmod3, .f = fixef, type = "parametric", B = 5000)
bb<-confint(gg, type = "norm")
bb



### unemotional subscale

T2_dat_long2 <- tidyr::gather(T2_dat[-c(which(T2_dat$bad_load_2==1)),][,which(colnames(T2_dat) %in% c("record_id", "mie2_complex", "mie1_complex","icu_unemotional" , "bad_resp2","icu_total" ,"sex_male", "age", "race_White", "anny2_diff","WAI_total","sdq_total_cd"))],time,value,mie1_complex,mie2_complex,mie1_complex:mie2_complex, factor_key=TRUE)

T2_dat_long2$time_rc <- as.numeric(car::recode(T2_dat_long2$time, "'mie2_complex'=1;'mie1_complex'=0"))

T2_dat_long2<-semTools::orthogonalize(T2_dat_long2, var1=c("icu_unemotional", "time_rc"), var2=c("time_rc","icu_unemotional"), match = TRUE)


llmod2<-lmer(value~ icu_unemotional + 
               time_rc + 
               icu_unemotional.time_rc +
               anny2_diff + 
               sex_male + 
               age + 
               race_White+
               sdq_total_cd+
               bad_resp2+
               (1|record_id),
             data=T2_dat_long2)



reg1<- lm(value~icu_unemotional.time_rc + 
            icu_unemotional + 
            time_rc + 
            anny2_diff + 
            sex_male + 
            age + 
            race_White+
            sdq_total_cd+
            bad_resp2,
          data=T2_dat_long2)



#### LR test to see if random effects improves model estimation
lrtest(reg1,llmod2)


> the addition of a random effect did improve model fit to the data - plausibly due to the random effect accounting for cluster confounding

#### Mean differences for random effects model
anova(llmod2)



> we see a ICU/time interaction suggesting that yes there is a differnece from T1-T2


#### Model estimates
summary(llmod2)

#### R2 for entire model
round(MuMIn::r.squaredGLMM(llmod2),3)



#### R2 by each IV
r2glmm::r2beta(llmod2)



#### Bootstrapped confidence intervals
set.seed(101)
ff<- lmeresampler::bootstrap(llmod2, .f = fixef, type = "parametric", B = 5000)
aa<-confint(ff, type = "norm")
aa






### figures 

#### ICU unemotional change T1 and T2
T2_dat_long22<-T2_dat_long2[complete.cases(T2_dat_long2),]
T2_dat_long22$random.intercpet.preds <- predict(llmod2)
T2_dat_long22$time_rc<-car::recode(T2_dat_long22$time,"'mie2_complex'=2;'mie1_complex'=1", as.factor = FALSE)
T2_dat_long22$hi_lo<-ifelse(T2_dat_long22$icu_unemotional>mean(T2_dat_long22$icu_unemotional),"ICU Unemotional > Mean","ICU Unemotional < Mean")

ggplot(data=T2_dat_long22, aes(x=icu_unemotional, y=random.intercpet.preds, group = factor(time_rc), colour = factor(time_rc))) +
  geom_smooth(method='lm', formula = 'y~x',se=FALSE, size=2) + geom_jitter()+
  labs(x="Unemotional Traits", y="Complex Affective Theory of Mind") +
  scale_colour_grey("Time", start = .7, end= .2)  + theme_minimal() +
  theme(axis.text=element_text(size=rel(1.2)),axis.title=element_text(size=rel(1.2)), legend.title = element_text(colour = "black",size=rel(1.2))) 
## 500x500 or 5 x 5 in "in" units
#ggsave(path = "C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\DPRG main project\\Analyses and Data\\full_data_87\\1_Paper-1", width = 6, height = 5, filename = "unemotional_time_lines.tiff", device='tiff', dpi=700, limitsize = FALSE)



> A sharper decrease in identifying complex emotions at higher unemotional traits









#### ICU total change T1 and T2
T2_dat_long22$hi_lo2<-ifelse(T2_dat_long22$icu_total>mean(T2_dat_long22$icu_total),"ICU total > Mean","ICU total < Mean")

ggplot(data=T2_dat_long22, aes(x=icu_total, y=random.intercpet.preds, group = factor(time_rc), colour = factor(time_rc))) +
  geom_smooth(method='lm', formula = 'y~x',se=FALSE, size=2) + geom_point()+
  labs(x="Callous-Unemotional Traits", y="Complex Affective Theory of Mind") +
  scale_colour_grey("Time", start = .7, end= .2)  + theme_minimal() +
  theme(axis.text=element_text(size=rel(1.2)),axis.title=element_text(size=rel(1.2)), legend.title = element_text(colour = "black",size=rel(1.2)))

## 500x500 or 5 x 5 in "in" units
#ggsave(path = "C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\DPRG main project\\Analyses and Data\\full_data_87\\1_Paper-1", width = 6, height = 5, filename = "CU_time_lines.tiff", device='tiff', dpi=700, limitsize = FALSE)


#### Unemotional change random eff
T2_dat_long22$hi_lo2<-ifelse(T2_dat_long22$icu_total>(mean(T2_dat_long22$icu_total)+1.5),"Unemotion > Mean","Unemotion < Mean")

T1 <- expression(paste(bold("T1")))
T2 <- expression(paste(bold("T2")))

T2_dat_long22 %>%
  mutate(hi_lo = fct_relevel(hi_lo2, "Unemotional < Mean","Unemotional > Mean")) %>%
  ggplot(aes(x=factor(time), y=random.intercpet.preds, group = record_id, colour = icu_total)) +
  geom_smooth(method='lm',formula ='y ~ x') + geom_point()+
  labs(x=NULL, y="Complex Affective Theory of Mind") +
  scale_colour_gradient("Unemotional \nTraits", low="light grey",high="black") + facet_wrap(~hi_lo)  + scale_x_discrete(labels=c(T1, T2)) + theme_minimal() +
  theme(axis.text=element_text(size=rel(1.1)),axis.title=element_text(size=rel(1.6)), strip.text.x = element_text(colour = "black", face = "bold",size=rel(2)), legend.title = element_text(size=rel(1.3)), legend.text = element_text(size=rel(1.1)))

## 500x500 or 5 x 5 in "in" units
#ggsave(path = "C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\DPRG main project\\Analyses and Data\\full_data_87\\1_Paper-1", width = 6, height = 5, filename = "aff_load_unemotional.tiff", device='tiff', dpi=700, limitsize = FALSE)


#### ICU total change random eff
T2_dat_long2 <- tidyr::gather(T2_dat[-c(which(T2_dat$bad_load_2==1),c(41,7,5,29,40)),][,which(colnames(T2_dat) %in% c("record_id", "mie2_complex", "mie1_complex","icu_unemotional" ,"icu_total" , "sex_male", "age", "race_White", "anny2_diff","WAI_total","sdq_total_cd"))],time,value,mie1_complex,mie2_complex,mie1_complex:mie2_complex, factor_key=TRUE)

llmod3<-lmer(value~icu_total * time + anny2_diff + sex_male + age + race_White +(1|record_id),data=T2_dat_long2)
T2_dat_long22$random.intercpet.preds3 <- predict(llmod3)

T2_dat_long22$hi_lo2<-ifelse(T2_dat_long22$icu_total>(mean(T2_dat_long22$icu_total)+1.5),"CU Traits > Mean","CU Traits < Mean")

T1 <- expression(paste(bold("T1")))
T2 <- expression(paste(bold("T2")))

T2_dat_long22 %>%
  mutate(hi_lo = fct_relevel(hi_lo2, "CU Traits < Mean","CU Traits > Mean")) %>%
  ggplot(aes(x=factor(time), y=random.intercpet.preds3, group = record_id, colour = icu_total)) +
  geom_smooth(method='lm',formula ='y ~ x') + geom_point()+
  labs(x=NULL, y="Complex Affective Theory of Mind") +
  scale_colour_gradient("CU Traits", low="grey",high="black") + facet_wrap(~hi_lo)  + scale_x_discrete(labels=c(T1, T2)) + theme_minimal() +
  theme(axis.text=element_text(size=rel(1.1)),axis.title=element_text(size=rel(1.6)), strip.text.x = element_text(colour = "black", face = "bold",size=rel(2)), legend.title = element_text(size=rel(1.3)), legend.text = element_text(size=rel(1.1)))

## 500x500 or 5 x 5 in "in" units
#ggsave(path = "C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\DPRG main project\\Analyses and Data\\full_data_87\\1_Paper-1", width = 6, height = 5, filename = "aff_load_totalCU.tiff", device='tiff', dpi=700, limitsize = FALSE)




## Examining level of cognitive load
T2_dat_long4 <- tidyr::gather(T2_dat[-c(which(T2_dat$bad_load_2==1)),][,which(colnames(T2_dat) %in% c("record_id", "mie2_complex", "mie1_complex","icu_unemotional" ,"icu_total" , "bad_resp2", "sex_male", "age", "race_White", "anny2_diff","WAI_total","sdq_total_cd","load2_max"))],time,value,mie1_complex,mie2_complex,mie1_complex:mie2_complex, factor_key=TRUE)

T2_dat_long4$time_rc <- as.numeric(car::recode(T2_dat_long4$time, "'mie2_complex'=1;'mie1_complex'=0"))


llmod4<-lmer(value~ 
               factor(load2_max)*icu_total*time_rc + 
               (1|record_id),
             data=T2_dat_long4)
summary(llmod4)

z<-data.frame(r2glmm::r2beta(llmod4))[,c(1,6)]
z$Rsq<-round(z$Rsq,3)*100

z[grepl("\\d+:icu_total:time_rc" , z$Effect  ),]




eff<-llmod4@beta[58:70]
std.er<-c(0.59981,0.30210,0.63118,0.77534,0.52850,0.33528,0.30210,0.39964,0.23858,3.00819,0.25891,0.24698,0.20890)
hi<-eff+std.er+.2
lo<-eff-std.er-.2
index<-rep(1:13)
labels<-(as.numeric(levels(unique(llmod4@frame$`factor(load2_max)`)))[c(2,4,5,6,9,10,13,16:21)])
dd<-data.frame(index,labels,eff,std.er,hi,lo)

ggplot(dd,aes(y=index, x=eff, xmin=lo, xmax=hi))+ 
  geom_point()+ #this adds the effect sizes to the plot
  #this changes the features of the overall effects
  #one could resize, reshape, or recolor these points if desired
  geom_point(data=dd, color="Black", size=3)+ 
  geom_errorbarh(height=.3, size=1.25)+ #adds the CIs
  #sets the scales
  #note that I reverse the y axis to correctly order the effect #sizes based on my index variable
  scale_x_continuous(limits=c(-4,6), name="Complex Affective Theory of Mind")+
  scale_y_continuous(name = "Cognitive Load (ms)", breaks=1:13, labels = dd$labels, trans="reverse")+
  
  #adding a vertical line at the effect = 0 mark
  geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5)+
  geom_text(aes(label=paste0("B = ", round(eff,2),c("**","","","*","","","","","","","","",""),", R2 = ",z[grepl("\\d+:icu_total:time_rc" , z$Effect  ),]$Rsq[c(1,13,8,2,3,11,9,6,10,7,5,12,4)],"%")), vjust=-.7) + 
  # adding a line at the overall effect size finding
  # geom_vline(xintercept=mean(dd$eff), color="black", linetype="solid", alpha=.5)+
  
  #faceting based on my subgroups
  #facet_grid(design~., scales= "free", space="free")+
  
  #thematic stuff
  #ggtitle("CU Traits and Alcohol Use Frequency")+
  theme_classic()+
  theme(text=element_text(size=14, color="black"))+
  theme(panel.spacing = unit(1, "lines"))

## 500x500 or 7 x 4.5 in "in" units
#ggsave(path = "C:\\Users\\wintersd\\OneDrive - The University of Colorado Denver\\DPRG main project\\Analyses and Data\\full_data_87\\1_Paper-1", width = 6.5, height = 4.5, filename = "max_load_beta.tiff", device='tiff', dpi=700, limitsize = FALSE)









## Showing analysis for basic
### total CU

T2_dat_long2 <- tidyr::gather(T2_dat[-c(which(T2_dat$bad_load_2==1)),][,which(colnames(T2_dat) %in% c("record_id", "mie2_basic", "mie1_basic","icu_unemotional" ,"icu_total" , "bad_resp2", "sex_male", "age", "race_White", "anny2_diff","WAI_total","sdq_total_cd"))],time,value,mie1_basic,mie2_basic,mie1_basic:mie2_basic, factor_key=TRUE)

T2_dat_long2$time_rc <- as.numeric(car::recode(T2_dat_long2$time, "'mie2_basic'=1;'mie1_basic'=0"))


T2_dat_long2<-semTools::orthogonalize(T2_dat_long2, var1=c("icu_total", "time_rc"), var2=c("time_rc","icu_total"), match = TRUE)

llmod3<-lmer(value~ icu_total + 
               time_rc + 
               icu_total.time_rc + 
               anny2_diff + 
               sex_male + 
               age + 
               sdq_total_cd + 
               race_White +
               bad_resp2 +
               (1|record_id),
             data=T2_dat_long2)

summary(llmod3)


### unemotional 

T2_dat_long2 <- tidyr::gather(T2_dat[-c(which(T2_dat$bad_load_2==1)),][,which(colnames(T2_dat) %in% c("record_id", "mie2_basic", "mie1_basic","icu_unemotional" ,"icu_total" , "bad_resp2", "sex_male", "age", "race_White", "anny2_diff","WAI_total","sdq_total_cd"))],time,value,mie1_basic,mie2_basic,mie1_basic:mie2_basic, factor_key=TRUE)

T2_dat_long2$time_rc <- as.numeric(car::recode(T2_dat_long2$time, "'mie2_basic'=1;'mie1_basic'=0"))


T2_dat_long2<-semTools::orthogonalize(T2_dat_long2, var1=c("icu_unemotional", "time_rc"), var2=c("time_rc","icu_unemotional"), match = TRUE)

llmod3<-lmer(value~ icu_unemotional + 
               time_rc + 
               icu_unemotional.time_rc + 
               anny2_diff + 
               sex_male + 
               age + 
               sdq_total_cd + 
               race_White +
               bad_resp2 +
               (1|record_id),
             data=T2_dat_long2)

summary(llmod3)





















