# LONGITUDINAL DATA ANALYSIS GLS No STRATA ----

library(ggeffects)

## GLS DATASET

#dset.gls<-dplyr::select(dset.ana, ... )

dset.gls<-dset.ana

dset.gls<-dset.gls[, -c(59, 74)]

# WIDE 2 LONG

dset.gls<-reshape(dset.gls, 
                  varying=352:384, 
                  direction = "long", 
                  idvar = "id", 
                  timevar = "TIME")

# TREAT FIRST TIME POINT AS COVAR

baseline<-subset(data.frame(dset.gls), TIME == 1, -TIME )
baseline<-upData(baseline, rename=c(sESelectin = "sESelectin.0", Thrombomodulin ="Thrombomodulin.0", IL6="IL6.0",
                                    TNF="TNF.0", Angiopoietin="Angiopoietin.0", Endoglin = "Endoglin.0", 
                                    Endothelin = "Endothelin.0", vWF="vWF.0", sICAM1 = "sICAM1.0", sVCAM1 = "sVCAM1.0",
                                    PSELECTIN="PSELECTIN.0"), print =FALSE )

followup<-subset(data.frame(dset.gls), TIME > 1, c(id, age, sex, TIME, 
                                                   sESelectin, Thrombomodulin, IL6, TNF, Angiopoietin, Endoglin, Endothelin, vWF, sICAM1, sVCAM1, PSELECTIN))

dset.gls.complete<-merge(baseline, followup, by="id" )

dset.gls$TIME<-as.factor(dset.gls$TIME)

dset.gls.complete$TIME<-as.factor(dset.gls.complete$TIME)

dset.gls$DEATH<-as.factor(dset.gls$DEATH)

dset.gls.complete$DEATH<-as.factor(dset.gls.complete$DEATH)


# MODELS

ddist<-datadist(dset.gls.complete)
options(datadist='ddist')

require(nlme)

# sVCAM1

cp<-list(corCAR1, corExp, corCompSymm, corLin, corGaus, corSpher)

z<-vector('list', length (cp))

for(k in 1:length(cp)){
  
  z[[k]]<-gls(sVCAM1 ~ TIME,
              data=dset.gls,
              correlation=cp[[k]](form=~TIME|id),
              na.action = na.exclude
  )
}

anova(z[[1]], z[[2]], z[[3]], z[[4]], z[[5]], z[[6]])

fit<-gls(sVCAM1 ~ TIME,
         data=dset.gls,
         correlation=corCompSymm(form=~TIME|id),
         na.action = na.omit)

summary(fit)

write.csv2(summary(fit)$tTable, file=paste0(outdir, "/sVCAM1.gls.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/sVCAM1.confint.csv"))

jpeg(file=paste0(outdir, "/sVCAM1.EFFECTS.jpg"), height=4, width=4, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME"))
write.csv2(dat, file=paste0(outdir, "/sVCAM1.marginal.csv"))
p<-plot(dat,connect.lines = TRUE)
p<-p+theme_cleveland()
p<-p+ggtitle("sVCAM1")
p<-p+theme(plot.title = element_text(hjust = 0.5))
p<-p+ylab("sVCAM1")
p<-p+annotate("text", x = 2, y = 1420, label = paste0("p=", round(summary(fit)$tTable[2, 4], 3)), size = 5)
p<-p+annotate("text", x = 2.9, y = 1380, label = paste0("p=", round(summary(fit)$tTable[3, 4], 3)), size = 5)
print(p)
dev.off()

emm <- emmeans(fit, ~ TIME, mode = "boot")
ct<-as.data.frame(contrast(emm, interaction = "pairwise"))
ct
write.csv2(ct, file=paste0(outdir, "/sVCAM1.ct.csv"))

# LONGITUDINAL DATA ANALYSIS GLS DEATH ----

library(ggeffects)

## GLS DATASET

#dset.gls<-dplyr::select(dset.ana, ... )

dset.gls<-dset.ana[which(dset.ana$GRP==1), ]

dset.gls<-dset.gls[, -c(59, 74)]

# WIDE 2 LONG

dset.gls<-reshape(dset.gls, 
                  varying=352:384, 
                  direction = "long", 
                  idvar = "id", 
                  timevar = "TIME")

# TREAT FIRST TIME POINT AS COVAR

baseline<-subset(data.frame(dset.gls), TIME == 1, -TIME )
baseline<-upData(baseline, rename=c(sESelectin = "sESelectin.0", Thrombomodulin ="Thrombomodulin.0", IL6="IL6.0",
                                    TNF="TNF.0", Angiopoietin="Angiopoietin.0", Endoglin = "Endoglin.0", 
                                    Endothelin = "Endothelin.0", vWF="vWF.0", sICAM1 = "sICAM1.0", sVCAM1 = "sVCAM1.0",
                                    PSELECTIN="PSELECTIN.0"), print =FALSE )

followup<-subset(data.frame(dset.gls), TIME > 1, c(id, age, sex, TIME, 
                                                   sESelectin, Thrombomodulin, IL6, TNF, Angiopoietin, Endoglin, Endothelin, vWF, sICAM1, sVCAM1, PSELECTIN))

dset.gls.complete<-merge(baseline, followup, by="id" )

dset.gls$TIME<-as.factor(dset.gls$TIME)

dset.gls.complete$TIME<-as.factor(dset.gls.complete$TIME)

dset.gls$DEATH<-as.factor(dset.gls$DEATH)

dset.gls.complete$DEATH<-as.factor(dset.gls.complete$DEATH)

# CV

cv.sVCAM1<-sd(dset.gls$sVCAM1, na.rm = T)/abs(mean(dset.gls$sVCAM1, na.rm= T))

cv.thrombomodulin<-sd(dset.gls$Thrombomodulin, na.rm = T)/abs(mean(dset.gls$Thrombomodulin, na.rm= T))

cv.endoglin<-sd(dset.gls$Endoglin, na.rm = T)/abs(mean(dset.gls$Endoglin, na.rm= T))

cv.endothelin<-sd(dset.gls$Endothelin, na.rm = T)/abs(mean(dset.gls$Endothelin, na.rm= T))

cv.sICAM1<-sd(dset.gls$sICAM1, na.rm = T)/abs(mean(dset.gls$sICAM1, na.rm= T))

cv.IL6<-sd(dset.gls$IL6, na.rm = T)/abs(mean(dset.gls$IL6, na.rm= T))

cv.TNF<-sd(dset.gls$TNF, na.rm = T)/abs(mean(dset.gls$TNF, na.rm= T))

cv.Angiopoietin<-sd(dset.gls$Angiopoietin, na.rm = T)/abs(mean(dset.gls$Angiopoietin, na.rm= T))

cv.vWF<-sd(dset.gls$vWF, na.rm = T)/abs(mean(dset.gls$vWF, na.rm= T))

cv.PSELECTIN<-sd(dset.gls$PSELECTIN, na.rm = T)/abs(mean(dset.gls$PSELECTIN, na.rm= T))

cv.ESELECTIN<-sd(dset.gls$sESelectin , na.rm = T)/abs(mean(dset.gls$PSELECTIN, na.rm= T))

cv<-cbind(cv.sVCAM1,
          cv.thrombomodulin,
          cv.endoglin,
          cv.endothelin,
          cv.sICAM1,
          cv.IL6,
          cv.TNF,
          cv.Angiopoietin,
          cv.vWF,
          cv.PSELECTIN,
          cv.ESELECTIN)

cv<-round(cv, 3)

write.csv2(cv, file=paste0(outdir, "/cv.csv"))

# MODELS

ddist<-datadist(dset.gls)
options(datadist='ddist')

require(nlme)

# sVCAM1



cp<-list(corCAR1, corExp, corCompSymm, corLin, corGaus, corSpher)

z<-vector('list', length (cp))

for(k in 1:length(cp)){
  
  z[[k]]<-gls(sVCAM1 ~ DEATH*TIME,
              data=dset.gls,
              correlation=cp[[k]](form=~TIME|id),
              na.action = na.exclude
              )
}

anova(z[[1]], z[[2]], z[[3]], z[[4]], z[[5]], z[[6]])

fit<-gls(sVCAM1 ~ DEATH*TIME,
         data=dset.gls,
         correlation=corCompSymm(form=~TIME|id),
         na.action = na.omit)

summary(fit)

fit.lp<-predict(fit, type="lp")

write.csv2(summary(fit)$tTable, file=paste0(outdir, "/sVCAM1.gls.DEATH.csv"))

write.csv2(confint(fit), file=paste0(outdir, "/sVCAM1.confint.DEATH.csv"))

jpeg(file=paste0(outdir, "/sVCAM1.DEATH.EFFECTS.jpg"), height=4, width=4, unit="in", res=600)
dat <- ggpredict(fit, terms = c("TIME", "DEATH"))
write.csv2(dat, file=paste0(outdir, "/sVCAM1.marginal.DEATH.csv"))
p<-plot(dat,connect.lines = TRUE)
p<-p+theme_cleveland()
p<-p+ggtitle("sVCAM1")
p<-p+theme(plot.title = element_text(hjust = 0.5))
p<-p+ylab("sVCAM1")
p<-p+annotate("text", x = 2, y = 1850, label = paste0("p=", round(summary(fit)$tTable[5, 4], 3)), size = 5)
p<-p+annotate("text", x = 2.9, y = 1950, label = paste0("p=", round(summary(fit)$tTable[6, 4], 3)), size = 5)
print(p)
dev.off()

emm <- emmeans(fit, ~ DEATH*TIME, mode = "boot")
ct<-as.data.frame(contrast(emm, interaction = "pairwise"))
ct
write.csv2(ct, file=paste0(outdir, "/sVCAM1.ct.DEATH.csv"))
