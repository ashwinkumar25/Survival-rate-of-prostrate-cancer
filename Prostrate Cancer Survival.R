#*******************************************
# Enova Data Challenge
#*******************************************
# Ashwin Kumar Mathi (amathi6@uic.edu)
#======================================

library(ROCR)
library(caret)
library(mice)
library(DataExplorer)
library(e1071)
library(prediction)

#read files
getwd()
enova_trainset<-read.csv("training_data.csv", header=T, stringsAsFactors = T) 
enova_testset<-read.csv("(name)_score.csv", header=T, stringsAsFactors = T) 
backup <- enova_trainset
summary(enova_trainset)
sapply(enova_trainset,function(x) sum(is.na(x)))

# psa_6_months and tumor_6_months have more nulls and looks like the data was not missed
# at random hence dropping the columns. Also dropping ID and dates

View(enova_trainset)
enova_trainset <- enova_trainset[,-c(1,2,18,21)]
# Removing blank values in symptom
enova_trainset <- enova_trainset[!(enova_trainset$symptoms==""),]
View(enova_trainset)
str(enova_trainset)


enova_trainsetnew <- mice(enova_trainset, m = 3,method="pmm")
enova_trainsetnew <- complete(enova_trainsetnew,3)
backupenova_trainset <- enova_trainset

######################### UNIVARIATE ANALYSIS #################################
str(enova_trainset)
cols <- c("t_score","n_score","m_score","stage","race","family_history","first_degree_history",	"previous_cancer",	"smoker","side","symptoms",
          "rd_thrpy","h_thrpy","chm_thrpy","cry_thrpy","brch_thrpy","rad_rem","multi_thrpy","survival_1_year","survival_7_years")

enova_trainset[,cols]<-lapply(enova_trainset[,cols],factor)
enova_trainset <- enova_trainsetnew

write.csv(enova_trainset, file = "enova_trainset_clean.csv",row.names=FALSE)

#gleason_score
plot(density(enova_trainset$gleason_score,na.rm = T), col="red", lwd=2.5, main="gleason_score") 
hist(enova_trainset$gleason_score, breaks=10, col=c("steelblue", "red"))
#looks like a multimodal distribution, has few outliers


enova_trainset$gleason_score <- ifelse(enova_trainset$gleason_score<=6,"Critical-Low",ifelse(enova_trainset$gleason_score==7,"Critical-Medium","Critical-High"))
enova_trainset$gleason_score<- as.factor(enova_trainset$gleason_score)
table(enova_trainset$gleason_score)

t<-table(enova_trainset$gleason_score)
pt<-prop.table(t)*100
barplot(pt,col = c("green","orange"))

# t_score
t<-table(enova_trainset$t_score)
pt<-prop.table(t)*100
barplot(pt,col = c("steelblue","orange"))
#T1a  T1b  T1c  T2a  T2b  T2c  T3a  T3b  T3c   T4 
#907  864  922 1260 1220 1223 1544 1558 1460 4017

# n_score
t<-table(enova_trainset$n_score)
pt<-prop.table(t)*100
barplot(pt,col = c("steelblue","orange"))
#N0   N1   NX 
#9320 4126 1529 

# m_score
t<-table(enova_trainset$m_score)
pt<-prop.table(t)*100
barplot(pt,col = c("steelblue","orange"))
#M0   M1a   M1b   M1c 
#13764   435   295   481 

# stage
t<-table(enova_trainset$stage)
pt<-prop.table(t)*100
barplot(pt,col = c("steelblue","orange"))
#I  IIA  IIB  III   IV 
#516 1983 3371 2565 6540

# age
plot(density(enova_trainset$age,na.rm = T), col="red", lwd=2.5, main="Distribution of age") 
hist(enova_trainset$age, breaks=10, col=c("steelblue", "red"))
# looks like a multimodal distribution
# many patients are between 70 and 85, highest being 107

enova_trainset$age = cut(enova_trainset$age, c(1,70,80,108))
levels(enova_trainset$age)[levels(enova_trainset$age)== "(1,70]"] = "Age_Upto70"
levels(enova_trainset$age)[levels(enova_trainset$age)== "(70,80]"] = "Age_70to80"
levels(enova_trainset$age)[levels(enova_trainset$age)== "(80,108]"] = "Age_80to108"
table(enova_trainset$age)

# race
enova_trainset$race <- as.factor(enova_trainset$race)
t<-table(enova_trainset$race)
pt<-prop.table(t)*100
barplot(pt,col = c("steelblue","orange"))
# more patients belonging to race 4

# height
plot(density(enova_trainset$height,na.rm = T), col="red", lwd=2.5, main="Height") # multimodal

# weight
plot(density(enova_trainset$weight,na.rm = T), col="red", lwd=2.5, main="Weight") #almost normal

# create BMI 

enova_trainset$BMI <- ((enova_trainset$weight)/((enova_trainset$height)^2)*703)
View(enova_trainset)
summary(enova_trainset$BMI)
hist(enova_trainset$BMI, col=c("blue", "red"))

# family_history
t<-table(enova_trainset$family_history)
pt<-prop.table(t)*100
barplot(pt,col = c("steelblue","orange"))
#many patients have no history

# first_degree_history
t<-table(enova_trainset$first_degree_history)
pt<-prop.table(t)*100
barplot(pt,col = c("steelblue","orange"))
#many patients have zero first degree history

# previous_cancer
t<-table(enova_trainset$previous_cancer)
pt<-prop.table(t)*100
barplot(pt,col = c("steelblue","orange"))
# 0     1 
#13939  1036

# smoker
t<-table(enova_trainset$smoker)
pt<-prop.table(t)*100
barplot(pt,col = c("steelblue","orange"))

# side
t<-table(enova_trainset$side)
pt<-prop.table(t)*100
barplot(pt,col = c("steelblue","orange"))

# tumor_diagnosis
summary(enova_trainset$tumor_diagnosis)
plot(density(enova_trainset$tumor_diagnosis,na.rm = T), col="red", lwd=2.5, main="Tumor_diagnosis")
# almost normal 

# tumor_6_months removed

# tumor_1_year
summary(enova_trainset$tumor_1_year)
plot(density(enova_trainset$tumor_1_year,na.rm = T), col="red", lwd=2.5, main="Tumor_1_year")
skewness(enova_trainset$tumor_1_year)

#tumor_change
enova_trainset$tumor_change <- enova_trainset$tumor_1_year - enova_trainset$tumor_diagnosis

# psa_diagnosis
summary(enova_trainset$psa_diagnosis)
plot(density(enova_trainset$psa_diagnosis,na.rm = T), col="red", lwd=2.5, main="psa_diagnosis")
skewness(enova_trainset$psa_diagnosis) # normal 

# psa_6_months   removed
# psa_1_year
summary(enova_trainset$psa_1_year)
plot(density(enova_trainset$psa_1_year,na.rm = T), col="red", lwd=2.5, main="PSA_1_year")
skewness(enova_trainset$psa_1_year) # normal 

#psa_change
enova_trainset$psa_change <- enova_trainset$psa_1_year - enova_trainset$psa_diagnosis

# tea
summary(enova_trainset$tea)
plot(density(enova_trainset$tea,na.rm = T), col="red", lwd=2.5, main="tea")
skewness(enova_trainset$tea) # normal 


#symptoms

sym_symbols <- c("U03","U06","S07","U01","U02","S10","O11","U05","S04","P02","P01","O01","O09","O08","O10","P03")
symdataframe <- data.frame(matrix(ncol=16,nrow=0))
colnames(symdataframe) <- sym_symbols


for(i in 1:dim(enova_trainset)[1]){
  a <- as.list(strsplit(as.character(enova_trainset$symptoms[i]),",")[[1]])
  compare <- NULL
  for(j in 1:length(a)){
    compare <- rbind(compare,sym_symbols == a[j])
  }
  symdataframe <- rbind(symdataframe,apply(compare,2,sum))
} 

colnames(symdataframe) <- sym_symbols
View(symdataframe)

enova_trainset <- cbind(enova_trainset,symdataframe)
str(enova_trainset)
cols <- c("U03","U06","S07","U01","U02","S10","O11","U05","S04","P02","P01","O01","O09","O08","O10","P03")
enova_trainset[,cols]<-lapply(enova_trainset[,cols],factor)

# rd_thrpy
table(enova_trainset$rd_thrpy)
# 0    1 
# 6921 8054 

# h_thrpy
table(enova_trainset$h_thrpy)
# 0    1 
# 9833 5142 

# chm_thrpy
table(enova_trainset$chm_thrpy)
# 0    1 
# 5030 9945 

# cry_thrpy
table(enova_trainset$cry_thrpy)
# 0     1 
# 11396  3579 

# brch_thrpy
table(enova_trainset$brch_thrpy)
# 0     1 
# 11304  3671 

# rad_rem
table(enova_trainset$rad_rem)
# 0     1 
# 12363  2612 

# multi_thrpy
table(enova_trainset$multi_thrpy)
# 0     1 
# 3305 11670 

# survival_1_year
table(enova_trainset$survival_1_year)
# 0     1 
# 1549 13426 

# survival_7_years
table(enova_trainset$survival_7_years)
# 0    1 
# 8531 6444 

nrow(enova_trainset)
View(enova_trainset[,c(20)])

enova_trainset<-read.csv("enova_trainset_clean.csv", header=T, stringsAsFactors = T) 
str(enova_trainset)

View(enova_trainset)
# removing symptoms
enova_trainset <- enova_trainset[,-c(20)]

View(enova_trainset[,c(8,9)])

######################### BIVARIATE ANALYSIS #################################
# 
# gleason_score
chi1 <- table(enova_trainset$gleason_score,enova_trainset$survival_7_years)
chisq.test(chi1)  #significant
# t_score
chi2 <- table(enova_trainset$t_score,enova_trainset$survival_7_years)
chisq.test(chi2) #significant
# n_score
chi3 <- table(enova_trainset$n_score,enova_trainset$survival_7_years)
chisq.test(chi3)  # significant
# m_score
chi4 <- table(enova_trainset$m_score,enova_trainset$survival_7_years)
chisq.test(chi4)  # significant
# stage
chi5 <- table(enova_trainset$stage,enova_trainset$survival_7_years)
chisq.test(chi5)  # significant
# age
chi6 <- table(enova_trainset$age,enova_trainset$survival_7_years)
chisq.test(chi6)  # not significant

# race
chi7 <- table(enova_trainset$race,enova_trainset$survival_7_years)
chisq.test(chi7)   #not significant

# height
# not considered

# weight
# not considered

# family_history
chi8 <- table(enova_trainset$family_history,enova_trainset$survival_7_years)
chisq.test(chi8)   #not significant 

# first_degree_history
chi9 <- table(enova_trainset$first_degree_history,enova_trainset$survival_7_years)
chisq.test(chi9)   #not significant

# previous_cancer
chi10 <- table(enova_trainset$previous_cancer,enova_trainset$survival_7_years)
chisq.test(chi10)   #not significant

# smoker
chi11 <- table(enova_trainset$smoker,enova_trainset$survival_7_years)
chisq.test(chi11)   #not significant

# side
chi12 <- table(enova_trainset$side,enova_trainset$survival_7_years)
chisq.test(chi12)   #not significant

# tumor_diagnosis
a.aov <- aov(tumor_diagnosis~survival_7_years, data=enova_trainset)
a.aov
summary(a.aov)
# highly significant

# tumor_1_year
a.aov <- aov(tumor_1_year~survival_7_years, data=enova_trainset)
a.aov
summary(a.aov)
# highly significant

# psa_diagnosis
a.aov <- aov(psa_diagnosis~survival_7_years, data=enova_trainset)
a.aov
summary(a.aov)
# highly significant

# psa_1_year
a.aov <- aov(psa_1_year~survival_7_years, data=enova_trainset)
a.aov
summary(a.aov)
# highly significant

# rd_thrpy
chisq.test(table(enova_trainset$rd_thrpy,enova_trainset$survival_7_years))
# significant 

# h_thrpy
chisq.test(table(enova_trainset$h_thrpy,enova_trainset$survival_7_years))
# significant 

# chm_thrpy
chisq.test(table(enova_trainset$chm_thrpy,enova_trainset$survival_7_years))
# significant 

# cry_thrpy
chisq.test(table(enova_trainset$cry_thrpy,enova_trainset$survival_7_years))
# significant 

# brch_thrpy
chisq.test(table(enova_trainset$brch_thrpy,enova_trainset$survival_7_years))
# significant 

# rad_rem
chisq.test(table(enova_trainset$rad_rem,enova_trainset$survival_7_years))
#  almost significant 

# multi_thrpy
chisq.test(table(enova_trainset$multi_thrpy,enova_trainset$survival_7_years))
# significant 

# survival_1_year
chisq.test(table(enova_trainset$rd_thrpy,enova_trainset$survival_7_years))
# significant 

# survival_7_years

# U03
chisq.test(table(enova_trainset$U03,enova_trainset$survival_7_years))
# not significant 

# U06
chisq.test(table(enova_trainset$U06,enova_trainset$survival_7_years))
# not significant 

# S07
chisq.test(table(enova_trainset$S07,enova_trainset$survival_7_years))
#  significant 

# U01
chisq.test(table(enova_trainset$U01,enova_trainset$survival_7_years))
# not significant 

# U02
chisq.test(table(enova_trainset$U02,enova_trainset$survival_7_years))
# not significant 

# S10
chisq.test(table(enova_trainset$S10,enova_trainset$survival_7_years))
# not significant 

# O11
chisq.test(table(enova_trainset$O11,enova_trainset$survival_7_years))
# not significant 

# U05
chisq.test(table(enova_trainset$U05,enova_trainset$survival_7_years))
#  significant 

# S04
chisq.test(table(enova_trainset$S04,enova_trainset$survival_7_years))
# not significant 

# P02
chisq.test(table(enova_trainset$P02,enova_trainset$survival_7_years))
#  significant 

# P01
chisq.test(table(enova_trainset$P01,enova_trainset$survival_7_years))
#  significant 

# O01
chisq.test(table(enova_trainset$O01,enova_trainset$survival_7_years))
#  significant 

# O09
chisq.test(table(enova_trainset$O09,enova_trainset$survival_7_years))
#  significant 

# O08
chisq.test(table(enova_trainset$O08,enova_trainset$survival_7_years))
#  significant 

# O10
chisq.test(table(enova_trainset$O10,enova_trainset$survival_7_years))
#  significant 

# P03
chisq.test(table(enova_trainset$P03,enova_trainset$survival_7_years))
#  significant 

# BMI
a.aov <- aov(BMI~survival_7_years, data=enova_trainset)
a.aov
summary(a.aov)
#significant

# tumor_change
a.aov <- aov(tumor_change~survival_7_years, data=enova_trainset)
a.aov
summary(a.aov)
#significant

# psa_change
a.aov <- aov(psa_change~survival_7_years, data=enova_trainset)
a.aov
summary(a.aov)
#significant

View(enova_trainset)

# data split for train and validation
index<-sample(2,nrow(enova_trainset),replace=TRUE,prob=c(0.7,0.3))

data_train<-enova_trainset[index==1,]
data_test<-enova_trainset[index==2,]

# correlation of numeric variables

library(corrplot)
View(enova_trainset[,c(15,16,17,18,19,46)])
enova_cor <- enova_trainset[,c(15:19,45)]
cormat <- cor(enova_cor) 
round(cormat, 2) # Rounded to 2 decimals 
var1 <- corrplot(cormat, method="circle", addCoef.col="black",type="upper") 
corrplot(cormat, method="circle", addCoef.col="black") 
# tumor diagnosis and tumor 1 year are correlated

View(enova_trainset)

View(enova_trainset[,c(15:18)])

############################# LOGISTIC REGRESSION #################################

# model1
model1 <- glm(survival_7_years ~ ., data = data_train[,-c(8,9,15:18)], family = "binomial")
summary(model1)

#on training
predict <- predict(model1, type = 'response')
table(data_train$survival_7_years, predict > 0.5)
# 67% on training set 

# on testing
pred <- predict(model1,data_test[,-c(8,9,15:18)],type = "response")
pred <- predict(model1,data_test,type = "response")
p <- ifelse(pred>0.5,1,0)
p<-as.factor(p)
confusionMatrix(p,data_test$survival_7_years,positive = "1")


# model 2

View(data_train[,c(7,8,9,13,14,15:18,29,30,32,33,37)])

# removing insignificant variables from the bivariate analysis and previous results
model2 <- glm(survival_7_years ~ ., data = data_train[,-c(5,7:18,29,30,32,33,35,37)], family = "binomial")
summary(model2)

# performance 
# on train
pred <- predict(model2, type = 'response')
table(data_train$survival_7_years, pred > 0.5)

#on test
pred <- predict(model2,data_test,type = "response")
p <- ifelse(pred>0.5,1,0)
p<-as.factor(p)
confusionMatrix(p,data_test$survival_7_years,positive = "1")

# Model 3
model3 <- glm(survival_7_years ~ ., data = data_train[,-c(2,5,6,7:18,29,30,32,33,35,37)], family = "binomial")
summary(model3)

exp(cbind(Odd_Ratio=coef(model3),confint(model3)))

# performance 
# on train
pred <- predict(model3, type = 'response')
table(data_train$survival_7_years, pred > 0.5)

#on test
pred <- predict(model3,data_test,type = "response")
p <- ifelse(pred>0.5,1,0)
p<-as.factor(p)
confusionMatrix(p,data_test$survival_7_years,positive = "1")

#Print ROC Curve and AUC value

pr <- prediction(pred,data_test$survival_7_years)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
auc <- performance(pr,measure = "auc")
auc@y.values[[1]]

# step wise

null = glm(survival_7_years~1, data= data_train[,-c(8:15)], family = "binomial") 
full = glm(survival_7_years~., data= data_train[,-c(8:15)], family = "binomial")
step(null, scope=list(lower=null, upper=full), direction="both")


# model from stepwise
glm(formula = survival_7_years ~ survival_1_year + n_score + 
      tumor_change + gleason_score + rd_thrpy + U05 + S10 + multi_thrpy + 
      h_thrpy + tumor_1_year + rad_rem + brch_thrpy + age + BMI + 
      O09 + stage + P03 + race + S07 + psa_change + chm_thrpy + 
      O08 + P02 + P01, family = "binomial", data = data_train[,-c(8:15)])

model4  <- glm(survival_7_years ~ survival_1_year + n_score + 
                 tumor_change + gleason_score + rd_thrpy + U05 + S10 + multi_thrpy + 
                 h_thrpy + tumor_1_year + rad_rem + brch_thrpy + age + BMI + 
                 O09 + stage + P03 + race + S07 + psa_change + chm_thrpy + 
                 O08 + P02 + P01, data = data_train, family = "binomial")
summary(model4)

# performance 
# on train
pred <- predict(model4, type = 'response')
table(data_train$survival_7_years, pred > 0.5)

#on test
pred <- predict(model4,data_test,type = "response")
p <- ifelse(pred>0.5,1,0)
p<-as.factor(p)
confusionMatrix(p,data_test$survival_7_years,positive = "1")


                          # Preparing the score data set 
View(enova_testset)
nrow(enova_testset)
backenova_testset <- enova_testset


View(enova_testset[,c(1,2,18,21)])
enova_testset <- enova_testset[,-c(1,2,18,21)]
nrow(enova_testset)

# symptoms
sym_symbols <- c("U03","U06","S07","U01","U02","S10","O11","U05","S04","P02","P01","O01","O09","O08","O10","P03")
symdataframe <- data.frame(matrix(ncol=16,nrow=0))
colnames(symdataframe) <- sym_symbols


for(i in 1:dim(enova_testset)[1]){
  a <- as.list(strsplit(as.character(enova_testset$symptoms[i]),",")[[1]])
  compare <- NULL
  for(j in 1:length(a)){
    compare <- rbind(compare,sym_symbols == a[j])
  }
  symdataframe <- rbind(symdataframe,apply(compare,2,sum))
} 

colnames(symdataframe) <- sym_symbols
View(symdataframe)

enova_testset <- cbind(enova_testset,symdataframe)
str(enova_testset)
View(enova_testset)

write.csv(enova_testset, file = "enova_testset_clean.csv",row.names=FALSE)
write.csv(enova_testsetnew, file = "enova_testset_cleanNew.csv",row.names=FALSE)

enova_testsetnew <- enova_testset[,-c(20,29)] 
colsleft <- enova_testset[,c(20,29)] 
View(colsleft)
View(enova_testset)
View(enova_testsetnew)

# impute test set
View(enova_testsetnew)
enova_testsetnew <- mice(enova_testsetnew, m = 3,method="pmm")
enova_testsetnew <- complete(enova_testsetnew,3)
#
backupenova_testset <- enova_trainset

# bind the other columns
enova_testsetnew <- cbind(enova_testsetnew,colsleft)
View(enova_testsetnew)

sapply(enova_testsetnew,function(x) sum(is.na(x)))

# other variable changes

cols <- c("t_score","n_score","m_score","stage","race","family_history","first_degree_history",	"previous_cancer",	"smoker","side","symptoms",
          "rd_thrpy","h_thrpy","chm_thrpy","cry_thrpy","brch_thrpy","rad_rem","multi_thrpy","survival_1_year","survival_7_years")

enova_testsetnew[,cols]<-lapply(enova_testsetnew[,cols],factor)

cols <- c("U03","U06","S07","U01","U02","S10","O11","U05","S04","P02","P01","O01","O09","O08","O10","P03")
enova_testsetnew[,cols]<-lapply(enova_testsetnew[,cols],factor)

enova_testsetnew$gleason_score <- ifelse(enova_testsetnew$gleason_score<=6,"Critical-Low",ifelse(enova_trainset$gleason_score==7,"Critical-Medium","Critical-High"))
enova_testsetnew$gleason_score<- as.factor(enova_testsetnew$gleason_score)
table(enova_testsetnew$gleason_score)

t<-table(enova_testsetnew$gleason_score)
pt<-prop.table(t)*100
barplot(pt,col = c("green","orange"))

#age
enova_testsetnew$age = cut(enova_testsetnew$age, c(1,70,80,108))
levels(enova_testsetnew$age)[levels(enova_testsetnew$age)== "(1,70]"] = "Age_Upto70"
levels(enova_testsetnew$age)[levels(enova_testsetnew$age)== "(70,80]"] = "Age_70to80"
levels(enova_testsetnew$age)[levels(enova_testsetnew$age)== "(80,108]"] = "Age_80to108"
table(enova_testsetnew$age)

# create BMI 
enova_testsetnew$BMI <- ((enova_testsetnew$weight)/((enova_testsetnew$height)^2)*703)

#tumor_change
enova_testsetnew$tumor_change <- enova_testsetnew$tumor_1_year - enova_testsetnew$tumor_diagnosis

#psa_change
enova_testsetnew$psa_change <- enova_testsetnew$psa_1_year - enova_testsetnew$psa_diagnosis


# run on the models for score data

pred <- predict(model3,enova_testsetnew,type = "response")
p <- ifelse(pred>0.5,1,0)
p<-as.factor(p)
enova_testsetnew$survival_7_years <- p

table(enova_testsetnew$survival_7_years)

View(p)
count(p)

