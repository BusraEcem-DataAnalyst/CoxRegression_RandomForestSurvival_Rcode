library(ggplot2)
library(survminer)
library(survival)
library(gbm)
library(caret)
library(pROC)
library(tree)
library(ISLR)
library(vip)
library(e1071)
library(rminer)
library(tidyverse)
library(reshape2)
library(ggfortify)
library(rpart)
library(skimr)
library(kableExtra)
library(patchwork)
library(directlabels)
library(randomForest)
library(randomForestSRC)
library(pec)
library(prodlim)
library(ranger)
BreastCancer <- read.csv("C:/Users/info/Desktop/tez/Breastcancer_data.csv")
names(BreastCancer)<- c('zaman','Yas','olay','MenarsYas','OrtEmzirmeYılı','DoğumKontrolHapıKull','TümörTanı','NeoadjuvanKull')



#tanımlayıcı istatistik değerleri
summary(BreastCancer)
write.csv(summary(BreastCancer),"output.csv")
#frekans tabloları
library(epiDisplay)
tab1(BreastCancer$DoğumKontrolHapıKull, sort.group = "decreasing", cum.percent = TRUE)
tab1(BreastCancer$TümörTanı, sort.group = "decreasing", cum.percent = TRUE)
tab1(BreastCancer$NeoadjuvanKull, sort.group = "decreasing", cum.percent = TRUE)
tab1(BreastCancer$olay, sort.group = "decreasing", cum.percent = TRUE)

data.frame(zaman = BreastCancer$zaman, Yas = BreastCancer$Yas, olay=BreastCancer$olay,MenarsYas = BreastCancer$MenarsYas,
           OrtEmzirmeYili = BreastCancer$OrtEmzirmeYili, DogumKontrolHapiKull = BreastCancer$DogumKontrolHapiKull,
           TümörTani = BreastCancer$TümörTani, NeoadjuvanKull = BreastCancer$NeoadjuvanKull)

BreastCancer$DogumKontrolHapiKull <- factor(BreastCancer$DogumKontrolHapiKull,
                                            levels = c("1", "2"), 
                                            labels = c("Evet", "Hayır"))
BreastCancer$TümörTani <- factor(BreastCancer$TümörTani,
                                 levels = c("1", "2"), 
                                 labels = c("Evet", "Hayır"))
BreastCancer$NeoadjuvanKull <- factor(BreastCancer$NeoadjuvanKull,
                                      levels = c("1", "2"), 
                                      labels = c("Evet", "Hayır"))
#kaplan Meier
surv_object <- Surv(time = BreastCancer$zaman, event = BreastCancer$olay)
kaplan.meier<-survfit(surv_object~1,BreastCancer)
kaplan.meier
summary(kaplan.meier)
plot(kaplan.meier, col = c("Red"), xlab = "Gün", ylab =
       "Sağkalım Oranı", main =
       "Kaplan Meier Eğrisi", mark.time = TRUE)
legend(x = 500, y = 1, lty = 1:2, cex = .95, bty = "n")
log.rank1 <- survfit(surv_object~ TümörTani, data = BreastCancer)
summary(log.rank1)
print(log.rank1)
ggsurvplot(log.rank1, data = BreastCancer, pval = TRUE, xlab="zaman,gün", ylab="Sağkalım Eğrisi",
           legend.title="Tümör Tanı", legend.labs=c("Erken Evre","Geç Evre"))

log.rank2 <- survfit(surv_object~ DogumKontrolHapiKull, data = BreastCancer)
summary(log.rank2)
ggsurvplot(log.rank2, data = BreastCancer, pval = TRUE, xlab="zaman,gün", ylab="Sağkalım Eğrisi",
           legend.title="Doğum Kontrol Hapı Kullanımı", legend.labs=c("Evet","Hayır"))

log.rank3 <- survfit(Surv(zaman,olay)~ NeoadjuvanKull, data = BreastCancer)
summary(log.rank3)
ggsurvplot(log.rank3, data = BreastCancer, pval = TRUE, xlab="zaman,gün", ylab="Sağkalım Eğrisi",
           legend.title="Neoadjuvan Uygulanımı", legend.labs=c("Evet","Hayır"))
#Cox regresyon
Cox.Model<-coxph(Surv(zaman,olay)~.,data=BreastCancer,x = TRUE)
Cox.Model
ggforest(Cox.Model,data=BreastCancer)
summary(Cox.Model)
sum.surv <- summary(Cox.Model)
c_indexCox <- sum.surv$concordance
c_indexCox
#Kaplan Meier ve Cox regresyon karşılaştırma#

plot(kaplan.meier, conf.int = F, col = "black", main = "Model Uyumlarının Karşılaştırılması",
     xlab = "Zaman,gün", ylab = "Sağkalım oranı")

lines(survfit(Cox.Model, conf.int = F), col = "red")
legend(x = 700, y = 1, legend = c("Kaplan-Meier", "Cox-PH"), lty = 1,
       col = c("#238b45", "red"),
       cex = 1, bty = "n")
#randomforest
set.seed(800)
train <- sample(nrow(BreastCancer), 0.7*nrow(BreastCancer), replace = FALSE)
TrainSet <- BreastCancer[train,]
ValidSet <- BreastCancer[-train,]
summary(TrainSet)
summary(ValidSet)
fitform1 <- Surv(zaman,olay)~Yas+MenarsYas+OrtEmzirmeYili+DogumKontrolHapiKull+TümörTani+NeoadjuvanKull
set.seed(123)
fit<-rfsrc(fitform1, data = TrainSet, ntree = 80,splitrule = "logrank",importance = TRUE)
plot(fit)
get.cindex(time = TrainSet$zaman, censoring = TrainSet$olay, predicted = fit$predicted.oob)
plot.survival.rfsrc(fit,plots.one.page = FALSE,cens.model = "rfsrc")
#Prediction error curve
extends <- function(...) TRUE
library("doMC")
library("pec")
library("survival")
library(Rcpp)
registerDoMC()
set.seed(0692)
fitpec1 <- pec(list("CPH" = Cox.Model, "RSF" = fit), data = BreastCancer,
               formula = fitform1, splitMethod = "cv10", B = 6,
               keep.index = TRUE, keep.matrix = TRUE)
plot(fitpec1, what = "crossvalErr", xlim = c(0, 1022), legend = F)
legend(x = 890, y = 0.30, legend = c("KM", "CPH", "RSF"), lty = 1,
       col = c("black", "red", "green"), bty = "n", cex = 1, lwd = 2)
title("Comparison of Prediction Error Curves", line = 1, cex = 6)
#C-Index
startTime <- Sys.time()
set.seed(0692)
ApparrentCindex1 <- cindex(list("Cox" = Cox.Model, "RSF" = fit),
                           formula = fitform1, data = BreastCancer,
                           eval.times = seq(1, 1022, 1))
endTime <- Sys.time()
(totalRunTime <- endTime - startTime)
plot(ApparrentCindex1, legend = F, xlim=c(0,1000),ylim=c(0,1.0),col = c("red", "green"))
legend(x = 800, y = 1, legend = c("CPH", "RSF"), lty = 1, col = c("red", "green"), bty = "n", cex = 1, lwd = 2)
title("Comparison of Concordance", line = 1, cex = 6)
