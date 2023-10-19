# DCOX.TEMPLATE v1

source("/Users/mmanfrini/Analisi/template/codice/funbox.v2.R")

# PATH ----

wd<-setwd("~/Analisi/Stellarex/data.extraction/ESREFO26")
indir=paste0(wd, "/preprocess")
outdir=paste0(wd,"/out/reg")

# GLOBAL VARS ----

RAW=F

colrange=c()

# LOAD DATASET PRE PROCESSED ----

if(RAW==F){
  dset.ana <- read.csv2(paste0(indir,"/dset.ana.regclus.csv"),
                    stringsAsFactors=FALSE) #, na.strings = ""
}

# Load dataset descriptor file
nm<-read.csv2(paste0(indir,"/nm.dset.ana.filtered.regclus.csv"), 
              stringsAsFactors=T)


# SET VAR TYPE ----

setVarType(paste0(indir, "/nm.dset.ana.filtered.regclus.csv"))

if(length(numericList)>0){
  dset.ana <- applyNumeric(dset.ana, numericList)  
} else {
  print("No numeric vars")
}
if(length(integerList)>0){
  dset.ana <- applyInteger(dset.ana, integerList)
} else {
  print("No integer vars")
}
if(length(factorList)>0){
  dset.ana <- applyFactors(dset.ana, factorList)
} else {
  print("No factor vars")
}
if(length(orderedList)>0){
  dset.ana <- applyOrdered(dset.ana, orderedList)
} else {
  print("No ordered vars")
}
if(length(logicalList)>0){
  dset.ana <- applyLogical(dset.ana, logicalList)
} else {
  print("No logical vars")
}
if(length(dates)>0){
  dset.ana <- applyDate(dset.ana, dates, "%d/%m/%Y")
} else {
  print("No date vars")
}

# VISUALIZE ANALYSIS DATASET
vis_dat_wrapper(dset.ana)

# COX UNADJUSTED REGRESSION TEMPO_FUP_GG DEP=DEATH ----

selected<-nm[which(nm$rg==1 | nm$oc==1), "id"]
dset.cox<-dset.ana[, selected]

ddist<-datadist(dset.cox)
options(datadist=ddist)

confounders<-nm[which(nm$cf==1), "id"]

#expositionList<-nm[which(nm$exp==1), "id"]

cox.univar.rms(dset.cox, dset.cox$FUPTLR, dset.cox$TLR.1Y, "Day", confounders, outdir, "cox.uni.confounders.TLR1Y")

#cox.univar.rms(dset.cox, dset.cox$GIORNIFUP, dset.cox$MORTE_totale, "Day", expositionList, outdir, "cox.uni.exposition.MORTE_totale")


# MODC1 COX MODEL TEMPO_FUP_GG DEP=DEATH REG=POLMONITEMILD_0_NONMILD1 ----

ddist=datadist(dset.ana)
options(datadist='ddist')
options(prType = latex)

units(dset.ana$TEMPO_FUP_GG) <- 'Day'
S<-Surv(dset.ana$FUPTLR, dset.ana$TLR.1Y)
fit1<-coxph(S ~ Target.Limb,
          data = dset.ana,
          x=T, y=T
)

summary(fit1)

print(fit1, coef = TRUE)

cox.zph(fit1)

print(fastbw(fit1, type = c("individual"), force=(1)), k.aic=log(dim(dset.std)[1]))

validation<-validate(fit1, bw=FALSE, B=1000)
print(validation)
calibrate(fit1)
plot(Predict(fit1))
plot(anova(fit1), cex=0.8)

#AUC=0.5*validation[1,5]+0.5
AUC=CalculateAucFromDxy(validation)[8,5]
# AUC
print(AUC)
# VIF
print(vif(fit1))

# EXPORT RESULTS

cox.out=data.frame()

s<-summary(fit1, Sesso=0, Ipertensione.arteriosa=0)
pval<-get_model_stats(fit1)$coefs[dim(get_model_stats(fit1)$coefs)[2]]
cox.out<-rbind(cox.out,
               cbind(names(s[seq(1, dim(s)[1], 2), 1]),
                     round(s[seq(2, dim(s)[1], 2), 4],3),
                     paste0(round(s[seq(2, dim(s)[1], 2),6],3), " - ", round(s[seq(2, dim(s)[1], 2),7],3)),
                     pval)
)

colnames(cox.out)<-c("Variable", "HR", "95ci", "p")

cox.out

write.csv2(cox.out, file=paste0(outdir, "/cox.multi.morte.no-stent.CSV.VOLUME.TOT.BIN.csv"), row.names = F)
# COMPARE C-indices ----

library(compareC)

compareC(dset.ana$TEMPO_FUP_GG, dset.ana$Decesso0no1si, fit16$linear.predictors, fit17$linear.predictors)

# ROC ANALYSIS ----

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


#

write.csv2(anova(fit17, fit16), file=paste0(outdir, "/model.comp.csv"))

write.csv2(BIC(fit17, fit16), file=paste0(outdir, "/AIC.comp.csv"))



# KM MODc16 ----

jpeg(file=paste0(outdir, "/KM.MODC16.jpg"), height=6, width=4, unit="in", res=600)


ddist=datadist(dset.ana)
options(datadist='ddist')

units(dset.ana) <- 'Day'
S<-Surv(dset.ana$TEMPO_FUP_GG, dset.ana$Decesso0no1si==1)

fit16<-coxph(S ~ 
               Età +
               #LDH +
               #WBC +
               PCR +
               #SatO2 +
               Ddimero.gLpicco +
               #Sesso2F1M +
               #COMORB +
               strata(POLMONITE),
             #RVSTRIANSI1NO0,
             data = dset.ana,
             x=T, y=T#, surv=TRUE
)

pvalue<-round(summary(fit16)$sctest[3],3)


fit16<-coxph(S ~ 
               Età +
               #LDH +
               #WBC +
               PCR +
               #SatO2 +
               Ddimero.gLpicco +
               #Sesso2F1M +
               #COMORB +
               (POLMONITE),
             #RVSTRIANSI1NO0,
             data = dset.ana,
             x=T, y=T#, surv=TRUE
)


nd=data.frame(Età=rep(median(dset.ana$Età),2), 
              POLMONITE=factor(c(1,3)), 
              #RVSTRIANSI1NO0=factor(c(0,0,1)),
              PCR=rep(median(dset.ana$PCR),2),
              Ddimero.gLpicco=rep(median(dset.ana$Ddimero.gLpicco),2)
)


pfit <- survfit(fit16, 
                newdata = nd,
                conf.type="none"
                #data=dset.ana
)

p4<-ggsurvplot(pfit,
               conf.int = F,
               #fun="cumhaz",
               legend ="bottom",
               legend.title = "", 
               palette = c("Black", "gray50", "Red"),
               #palette = "uchicago",
               xlab="Time (days)",
               ylab="Survival rate",
               xlim=c(0,30), ylim=c(0,1),
               break.time.by = 10,
               #xscale = "d_m",
               censor = F,
               cumevents = F,
               pval = F,
               ggtheme = theme_minimal(),
               tables.theme = theme_cleantable(),
               risk.table = F,
               fontsize = 3,
               risk.table.title = "No. at risk",
               tables.height = 0.25,
               tables.y.text = FALSE,
               legend.labs = c("POLMONITE=1", "POLMONITE=3"),
               title = paste0("Adjusted survival curves for the Cox model. \nPOLMONITE"),
               data=nd
)

p4$table <- ggpubr::ggpar(p4$table,
                          font.title = list(size = 10, face = "bold"))

p4$table <- p4$table + theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  axis.text = element_text(size = 0)
)

p4$plot <- p4$plot + theme(
  panel.border = element_rect(colour = "black", fill=NA, size=0.5),
  # panel.background = element_rect(fill = "lightcyan",
  #                                 colour = "lightcyan",
  #                                 size = 0.5, linetype = "solid"),
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank()
)

print(p4)

p4$plot<-p4$plot+
  ggplot2::annotate("text", x = 2, y = 0.4, label = paste0("p=", pvalue), size = 4) #+
#   ggplot2::annotate("text", x = 40, y = 0.87, label = paste0("<",coff, " ng/L"), size = 4) +
#   ggplot2::annotate("text", x = 36, y = 0.45, label = paste0("\u2265",coff, " ng/L"), size = 4)



p4
# 
dev.off()  



