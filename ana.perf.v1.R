# PERFORMANCE

# COMPARE C-indices ----

library(compareC)

compareC(dset.ana$TEMPO_FUP_GG, dset.ana$Decesso0no1si, fit16$linear.predictors, fit17$linear.predictors)
# ROC COX ----
# Si utilizza il rischio exp(lp) 

library(pROC)
library(verification)

ddist<-datadist(dset.cox)
options(datadist=ddist)

units(dset.cox$gg_voce) <- 'Day'
S<-Surv(dset.cox$gg_voce, as.numeric(dset.cox$voce))

fit1<-coxph(S ~ 
              hyp +
              prior_mi +
              tsg.bin +
              cluster(id),
            data = dset.cox
)

summary(fit1)

fit1.lp <- predict(fit1, type = "lp")

robj<-pROC::roc(dset.cox$voce, exp(fit1.lp),
                plot = F, smooth = F, auc = T, ci = T, algorithm = 1)

ra<-roc.area(as.numeric(as.character(dset.cox$voce)), exp(fit1.lp))

roc.out<-rbind(AUC=round(robj$auc,3),
               ci=round(robj$ci,3),
               p=round(ra$p.value,3))

write.csv2(roc.out, file=paste0(outdir, "/perf.csv"))

jpeg(file=paste0(outdir, "/roc.plot.jpg"), height=4, width=4, unit="in", res=300)

plot(roc(dset.cox$voce, exp(fit1.lp)), print.auc=F, col = "red")

text(0.4,0.2,paste0("AUC=", roc.out[1,1], " (", roc.out[2,1], " - ", roc.out[2,3], ")"))

#text(0.65,0.1,ifelse(roc.out[3,1]>0, paste0("p=", roc.out[3,1]), "p<0.001"))

dev.off()


# ROC COX TEST MODEL COMPARISON ----
# Si utilizza il rischio exp(lp)
library(pROC)


fit1<-coxph(S ~ POLMONITEMILD_0_NONMILD1,
            data = dset.ana,
            x=T, y=T
)

summary(fit1)

fit1.lp <- predict(fit1, type = "lp")

fit2<-coxph(S ~ RVSTRIANSI1NO0,
            data = dset.ana,
            x=T, y=T#, surv=TRUE
)


summary(fit2)

fit2.lp <- predict(fit2, type = "lp")

fit3<-coxph(S ~ POLMONITEMILD_0_NONMILD1 + RVSTRIANSI1NO0,
            data = dset.ana,
            x=T, y=T#, surv=TRUE
)

summary(fit3)

fit3.lp <- predict(fit3, type = "lp")

fit5<-coxph(S ~ POLMONITEMILD_0_NONMILD1 +
              RVSTRIANSI1NO0 +
              Calcificazioni,
            data = dset.ana,
            x=T, y=T#, surv=TRUE
)

fit5.lp <- predict(fit5, type = "lp")

fit6<-coxph(S ~ POLMONITEMILD_0_NONMILD1 +
              RVSTRIANSI1NO0 +
              Calcificazioni +
              AGATSTONSCORE +
              VERSAMENTOPERICARDICOmm,
            data = dset.ana,
            x=T, y=T#, surv=TRUE
)

fit6.lp <- predict(fit6, type = "lp")

fit7<-coxph(S ~ 
              AGATSTONSCORE +
              VERSAMENTOPERICARDICOmm,
            data = dset.ana,
            x=T, y=T#, surv=TRUE
)

fit7.lp <- predict(fit7, type = "lp")


fit8<-coxph(S ~ 
              AGATSTONSCORE +
              VERSAMENTOPERICARDICOmm +
              RVSTRIANSI1NO0,
            data = dset.ana,
            x=T, y=T#, surv=TRUE
)

fit8.lp <- predict(fit8, type = "lp")

fit9<-cph(S ~ 
            Età +
            AGATSTONSCORE +
            VERSAMENTOPERICARDICOmm +
            RVSTRIANSI1NO0,
          data = dset.ana,
          x=T, y=T#, surv=TRUE
)

fit9.lp <- predict(fit9, type = "lp")

fit12<-coxph(S ~ 
               Età +
               #SatO2 +
               #Sesso2F1M +
               #COMORB +
               POLMONITEMILD_0_NONMILD1,
             #RVSTRIANSI1NO0,
             data = dset.ana,
             x=T, y=T#, surv=TRUE
)

fit12.lp <- predict(fit12, type = "lp")

fit13<-coxph(S ~ 
               Età +
               #SatO2 +
               #Sesso2F1M +
               #COMORB +
               POLMONITEMILD_0_NONMILD1 +
               RVSTRIANSI1NO0,
             data = dset.ana,
             x=T, y=T#, surv=TRUE
)

fit13.lp <- predict(fit13, type = "lp")

fit14<-coxph(S ~ 
               Età +
               #LDH +
               #WBC +
               #PCR +
               #SatO2 +
               #Sesso2F1M +
               #COMORB +
               POLMONITE,
             #RVSTRIANSI1NO0,
             data = dset.ana,
             x=T, y=T#, surv=TRUE
)

fit14.lp <- predict(fit14, type = "lp")

fit15<-coxph(S ~ 
               Età +
               #LDH +
               #WBC +
               #PCR +
               #SatO2 +
               #Sesso2F1M +
               #COMORB +
               POLMONITE +
               RVSTRIANSI1NO0,
             data = dset.ana,
             x=T, y=T#, surv=TRUE
)

fit15.lp <- predict(fit15, type = "lp")

fit16<-coxph(S ~ 
               Età +
               #LDH +
               #WBC +
               PCR +
               #SatO2 +
               #Ddimero.gLpicco,# +
               #Sesso2F1M +
               #COMORB +
               POLMONITE,
             #RVSTRIANSI1NO0,
             data = dset.ana,
             x=T, y=T#, surv=TRUE
)

fit16.lp <- predict(fit16, type = "lp")

fit17<-coxph(S ~ 
               Età +
               #LDH +
               #WBC +
               PCR +
               #SatO2 +
               #Ddimero.gLpicco +
               #Sesso2F1M +
               #COMORB +
               POLMONITE +
               RVSTRIANSI1NO0,
             data = dset.ana,
             x=T, y=T#, surv=TRUE
)

fit17.lp <- predict(fit17, type = "lp")

plot(roc(dset.ana$Decesso0no1si, fit16.lp), print.auc=T, print.auc.x=0.4, print.auc.y=0.7, col = "blue")
plot(roc(dset.ana$Decesso0no1si, fit17.lp), print.auc=T, print.auc.x=0.9, print.auc.y=0.8, col = "red", add=T)
#plot(roc(dset.ana$Decesso0no1si, fit18.lp), print.auc=T, print.auc.x=0.6, print.auc.y=0.9, col = "black", add=T)
legend("topleft", legend=c("FULL", "REDUCED"),
       col=c("red", "blue"), lty=1, cex=0.8)
# plot(roc(dset.ana$Decesso0no1si, fit5.lp), print.auc=T, print.auc.x=0.5, print.auc.y=0.8, col = "orange", add=T)
# plot(roc(dset.ana$Decesso0no1si, fit6.lp), print.auc=T, print.auc.x=0.5, print.auc.y=0.8, col = "green", add=T)
rc<-roc.test(roc(dset.ana$Decesso0no1si, fit16.lp), roc(dset.ana$Decesso0no1si, fit17.lp))
text(0.4,0.2,paste0("p=",round(rc$p.value, 3)))

roc.test(roc(dset.ana$Decesso0no1si, fit16.lp), roc(dset.ana$Decesso0no1si, fit16.lp))

# ROC LOGISTC SINGLE PREDICTOR ---- 
library(pROC)
library(verification)

robj<-pROC::roc(dset.ana$Infezione, dset.ana$AAT.CONC,
                plot = T, smooth = T, auc = T, ci = T, algorithm = 1)

ra<-roc.area(dset.ana$Infezione, predict(logi.fit))

roc.out<-rbind(AUC=robj$auc,
               ci=robj$ci,
               p=ra$p.value)

write.csv2(roc.out, file=paste0(outdir, "/perf.csv"))

# ROC LOGISTC MULTIPLE PREDICTORS ---- 
library(pROC)
library(verification)

robj<-pROC::roc(dset.ana$Infezione, predict(fit, type = "lp"),
                plot = F, smooth = F, auc = T, ci = T, algorithm = 1)

ra<-roc.area(dset.ana$Infezione, predict(fit, type = "lp"))

roc.out<-rbind(AUC=round(robj$auc,3),
               ci=round(robj$ci,3),
               p=round(ra$p.value,3))

write.csv2(roc.out, file=paste0(outdir, "/perf.csv"))

jpeg(file=paste0(outdir, "/roc.plot.jpg"), height=4, width=4, unit="in", res=300)

plot(roc(dset.cox$voce, predict(fit, type = "lp")), print.auc=F, col = "red")

text(0.4,0.2,paste0("AUC=", roc.out[1,1], " (", roc.out[2,1], " - ", roc.out[2,3], ")"))

#text(0.65,0.1,ifelse(roc.out[3,1]>0, paste0("p=", roc.out[3,1]), "p<0.001"))

dev.off()

# CUT POINT DEP=DEATH, VAR=sVCAM.1..ng.mL..t2 ----
library(cutpointr)

set.seed(123)

opt_cut <- cutpointr(dset.ana, 
                     tsg, 
                     voce, 
                     direction = ">=", 
                     pos_class = 1,
                     neg_class = 0, 
                     method = maximize_metric, 
                     boot_run = 100, 
                     summary_func = mean,
                     metric = youden, silent = TRUE,
                     #metric = sum_sens_spec, silent = TRUE,
                     na.rm = T)

s<-summary(opt_cut)$cutpointr[[1]][1:13]
s

write.csv2(s, file=paste0(outdir, "/svcam1.cutpoint.t3.csv"))

plot(opt_cut)

dset.cox$tsg.bin<-ifelse(dset.cox$tsg<=0.01,0,1)

# NESTED MODEL COMPARISON ----

write.csv2(anova(fit17, fit16), file=paste0(outdir, "/model.comp.csv"))

write.csv2(BIC(fit17, fit16), file=paste0(outdir, "/AIC.comp.csv"))