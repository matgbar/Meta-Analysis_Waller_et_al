#################################################################################
#Meta Analysis - Callous and unemotional traits

#Title: Opposite sides of the same coin: A meta-analysis... 

#Authors: Waller et al. 

#Analyses by Matthew Barstead (contact: barstead@umd.edu)

#################################################################################
library(metafor)
library(stargazer)
library(ggplot2)
#################################################################################
user<-Sys.getenv('USERPROFILE')
data.folder<-paste0(user, '\\Box Sync\\CU meta-analysis\\Meta_Raw_Data\\')
model.folder<-paste0(user, '\\Box Sync\\CU meta-analysis\\Output\\')
graphics.folder<-paste0(user, '\\Box Sync\\CU meta-analysis\\Graphics_Folder\\' )
#################################################################################
#Notes:   1. Correlations were transformed to Fisher's z
#         2. Standardized differences were transformed to r then to Fisher's z
#         3. Will require all results are exponentiated to return to -1:1 scale 
#################################################################################
dat<-read.csv(paste0(data.folder, 'Meta dataset_2-15-18.csv'), stringsAsFactors = F)
colnames(dat)[1]<-'id'
colnames(dat)[23]<-'Out_resp'

dat.Emp_tot<-dat[dat$Outcome=='empathy_tot',]
dat.Emp_aff<-dat[dat$Outcome=='empathy_aff',]
dat.Emp_cog<-dat[dat$Outcome=='empathy_cog',]
dat.glt<-dat[dat$Outcome=='guilt',]   #note probably have too few for guilt
dat.prosoc<-dat[dat$Outcome=='prosocial',]

#################################################################################
#some cleanup and recoding
dat$CU_resp.R[dat$CU_resp==0]<-'Self'
dat$CU_resp.R[dat$CU_resp==1]<-'Other'

dat$Out_resp.R[dat$Out_resp==0]<-'Self'
dat$Out_resp.R[dat$Out_resp==1]<-'Other'
dat$CU_resp.R[dat$citation=='Kimonis et al. (2016)']<-'Other'
dat$Out_resp.R[dat$citation=='Kimonis et al. (2016)']<-'Other'

#################################################################################
#Model for Total Empathy: 
fit.emp_tot<-rma(yi=Eff, vi=Eff_var, data=dat.Emp_tot, ni=N)
summary(fit.emp_tot)

sink(paste0(model.folder, 'Total Empathy Overall Model - no moderators.txt'))
summary(fit.emp_tot)
sink()

jpeg(paste0(graphics.folder, 'CU and Total Empathy - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_tot, order = 'obs', slab=dat.Emp_tot$citation)
title("Relation between CU Traits and Total Empathy (Fisher's z)")
dev.off()

jpeg(paste0(graphics.folder, 'CU and Total Empathy - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.emp_tot)
title("Funnel Plot of CU Relation with Total Empathy (Fisher's z)")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
REG.emp_tot<-regtest(fit.emp_tot, model = 'lm')
REG.emp_tot

#Trim and fill model
fit.emp_tot.TF_R<-trimfill(fit.emp_tot, estimator = 'R0')
fit.emp_tot.TF_R
funnel(fit.emp_tot.TF_R)

fit.emp_tot.TF_L<-trimfill(fit.emp_tot, estimator = 'L0')
fit.emp_tot.TF_L
funnel(fit.emp_tot.TF_L)


#file drawer analysis
fsn(yi=Eff, vi=Eff_var, data=dat.Emp_tot, type='Rosenthal', alpha = .05)

#################################################################################
#Model for Affective Empathy: 
fit.emp_aff<-rma(yi=Eff, vi=Eff_var, data=dat.Emp_aff)
summary(fit.emp_aff)

sink(paste0(model.folder, 'Affective Empathy Overall Model - no moderators.txt'))
summary(fit.emp_aff)
sink()

jpeg(paste0(graphics.folder, 'CU and Affective Empathy - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_aff, order = 'obs', slab=dat.Emp_aff$citation)
title("Relation between CU Traits and Affective Empathy (Fisher's z)")
dev.off()

jpeg(paste0(graphics.folder, 'CU and Affective Empathy - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.emp_aff)
title("Funnel Plot of CU Relation with Affective Empathy (Fisher's z)")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
REG.emp_aff<-regtest(fit.emp_aff, model = 'lm')
REG.emp_aff

#Trim and fill model - affective had no "imputed" studies
fit.emp_aff.TF_R<-trimfill(fit.emp_aff, estimator = 'R0')
fit.emp_aff.TF_R
funnel(fit.emp_aff.TF_R)

fit.emp_aff.TF_L<-trimfill(fit.emp_aff, estimator = 'L0')
fit.emp_aff.TF_L
funnel(fit.emp_aff.TF_L)

#file drawer analysis
fsn(yi=Eff, vi=Eff_var, data=dat.Emp_aff, type='Rosenthal', alpha = .05)

#################################################################################
#Model for Cognitive Empathy: 
fit.emp_cog<-rma(yi=Eff, vi=Eff_var, data=dat.Emp_cog)
summary(fit.emp_cog)

sink(paste0(model.folder, 'Cognitive Empathy Overall Model - no moderators.txt'))
summary(fit.emp_cog)
sink()

jpeg(paste0(graphics.folder, 'CU and Cognitive Empathy - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_cog, order = 'obs', slab=dat.Emp_cog$citation)
title("Relation between CU Traits and Cognitive Empathy (Fisher's z)")
dev.off()

jpeg(paste0(graphics.folder, 'CU and Cognitive Empathy - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.emp_cog)
title("Funnel Plot of CU Relation with Cognitive Empathy (Fisher's z)")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
REG.emp_cog<-regtest(fit.emp_cog, model = 'lm')
REG.emp_cog

#Trim and fill model
fit.emp_cog.TF_R<-trimfill(fit.emp_cog, estimator = 'R0')
fit.emp_cog.TF_R
funnel(fit.emp_cog.TF_R)

fit.emp_cog.TF_L<-trimfill(fit.emp_cog, estimator = 'L0')
fit.emp_cog.TF_L
funnel(fit.emp_cog.TF_L)

#file drawer analysis
fsn(yi=Eff, vi=Eff_var, data=dat.Emp_cog, type='Rosenthal', alpha = .05)

#################################################################################
#Model for Prosociality: 
fit.prosoc<-rma(yi=Eff, vi=Eff_var, data=dat.prosoc)
summary(fit.prosoc)

sink(paste0(model.folder, 'Prosociality Overall Model - no moderators.txt'))
summary(fit.prosoc)
sink()

jpeg(paste0(graphics.folder, 'CU and Prosociality - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.prosoc, order = 'obs', slab=dat.prosoc$citation)
title("Relation between CU Traits and Prosociality (Fisher's z)")
dev.off()

jpeg(paste0(graphics.folder, 'CU and Prosociality - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.prosoc)
title("Funnel Plot of CU Relation with Prosociality (Fisher's z)")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
REG.prosoc<-regtest(fit.prosoc, model = 'lm')
REG.prosoc

#Trim and fill model
fit.prosoc.TF_R<-trimfill(fit.prosoc, estimator = 'R0')
fit.prosoc.TF_R
funnel(fit.prosoc.TF_R)

fit.prosoc.TF_L<-trimfill(fit.prosoc, estimator = 'L0')
fit.prosoc.TF_L
funnel(fit.prosoc.TF_L)

#file drawer analysis
fsn(yi=Eff, vi=Eff_var, data=dat.prosoc, type='Rosenthal', alpha = .05)

#################################################################################
#Model for Guilt: 
fit.glt<-rma(yi=Eff, vi=Eff_var, data=dat.glt, method = 'FE')
summary(fit.glt)

sink(paste0(model.folder, 'Guilt Overall Model - no moderators.txt'))
summary(fit.glt)
sink()

jpeg(paste0(graphics.folder, 'CU and Guilt - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.glt, order = 'obs', slab=dat.glt$citation)
title("Relation between CU Traits and Guilt (Fisher's z)")
dev.off()

jpeg(paste0(graphics.folder, 'CU and Guilt - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.glt)
title("Funnel Plot of CU Relation with Guilt (Fisher's z)")
dev.off()

#############################################################################
#Combining and transforming back to Pearson's r
#Function for extracting vector of results (for use in creating a summary table)
rtoz_summary<-function(DF, ci, mod, outcome){
  #Extracting info about the data
  k<-length(DF[,1])
  N<-sum(DF$N)
  #Extracting Model information
  z<-mod$beta[,1]
  z.se<-mod$se
  z.LB<-z-z.se*qnorm(1-(1-ci)/2)
  z.UB<-z+z.se*qnorm(1-(1-ci)/2)
  rho<-(exp(z)-1)/(exp(z)+1)
  rho.LB<-(exp(z.LB)-1)/(exp(z.LB)+1)
  rho.UB<-(exp(z.UB)-1)/(exp(z.UB)+1)
  DF.temp<-data.frame(Outcome=outcome,
                      N=N, 
                      k=k, 
                      rho=round(rho, digits = 3),
                      rho.LB=round(rho.LB, digits = 3),
                      rho.UB=round(rho.UB, digits = 3)
                      )
  return(DF.temp)
}

DF.summary<-rbind(rtoz_summary(dat.Emp_tot, .95, fit.emp_tot, 'Total Empathy'),
                  rtoz_summary(dat.Emp_aff, .95, fit.emp_aff, 'Affective Empathy'),
                  rtoz_summary(dat.Emp_cog, .95, fit.emp_cog, 'Cognitive Empathy'),
                  rtoz_summary(dat.prosoc, .95, fit.prosoc, 'Prosociality'), 
                  rtoz_summary(dat.glt, .95, fit.glt, 'Guilt'),
                  rtoz_summary(dat.emp_comp, .95, fit.emp_comp, 'Affective-Cognitive')
                  )
#Vectors of Q values & I^2 vals
Q_vals<-c(fit.emp_tot$QE, fit.emp_aff$QE, fit.emp_cog$QE, fit.prosoc$QE, fit.glt$QE, fit.emp_comp$QE)
Q_p<-c(fit.emp_tot$QEp, fit.emp_aff$QEp, fit.emp_cog$QEp, fit.prosoc$QEp, fit.glt$QEp, fit.emp_comp$QEp)
I_2<-c(round(fit.emp_tot$I2, digits = 2), 
       round(fit.emp_aff$I2, digits = 2),
       round(fit.emp_cog$I2, digits = 2), 
       round(fit.prosoc$I2, digits =2),
       round(fit.glt$I2, digits = 2),
       round(fit.emp_comp$I2, digits = 2))

DF.summary<-data.frame(DF.summary, 
                       Q = Q_vals,
                       p = Q_p, 
                       I = I_2) 

stargazer(DF.summary, summary = F, out=paste0(model.folder,'Model Summary - no moderators.txt'), rownames = F, 
          notes = 'LB and UB based on 95% CI')

#Plotting summary effect sizes 
g1<-ggplot(aes(x=rho, y=Outcome), data=DF.summary)
g2<-g1+geom_errorbarh(aes(xmin=rho.LB, xmax=rho.UB), height=.25)
g3<-g2+geom_point(aes(size=N), pch=18)
g4<-g3+scale_size_continuous(range = c(1,12))
g5<-g4+xlab(expression(rho))+ylab('')+
  ggtitle(paste0('Correlation between Outcome Measures and CU Traits'))+
  theme(plot.title = element_text(hjust=.5))
g6<-g5+scale_x_continuous(limits = c(-.35, .05))+geom_vline(color='red', lty='dashed',
                                                          xintercept = 0)
g6

jpeg(paste0(graphics.folder, 'Model Summary Graphic.jpeg'), res=300, units='in', height=7, width=7)
g6
dev.off()

###########################################################################################
#Model comparing two forms of empathy - within study differences in effects
#merging and cleaning data set
dat.emp_comp<-merge(dat.Emp_aff, dat.Emp_cog, by='id')
cols<-c(1:5, 12,14,15,22,23,24,37,39,40)
dat.emp_comp<-dat.emp_comp[,cols]
colnames(dat.emp_comp)<-c('id', 'citation', 'female', 'age', 'N', 'R.affective',
                          'Eff_affective', 'Eff_var_affective', 'CU_resp', 'Out_resp',
                          'Sample', 'R_cognitive', 'Eff_cognitive', 'Eff_var_cognitive')

dat.emp_comp$CU_resp.R[dat.emp_comp$CU_resp==0]<-'Self'
dat.emp_comp$CU_resp.R[dat.emp_comp$CU_resp==1]<-'Other'

dat.emp_comp$Out_resp.R[dat.emp_comp$Out_resp==0]<-'Self'
dat.emp_comp$Out_resp.R[dat.emp_comp$Out_resp==1]<-'Other'

#Correlations pulled from studies: 
Aff_cog_cor<-c(.063, .063, .351, 0, 0, .49, 0, .5, 0, .76, .32,0,-.01,.08,0,.45,0,0)

#Calculating difference in effects and variance terms adjusted for outcome correlation
Eff<-dat.emp_comp$Eff_affective-dat.emp_comp$Eff_cognitive
Eff_var<-dat.emp_comp$Eff_var_affective+dat.emp_comp$Eff_var_cognitive-2*Aff_cog_cor*sqrt(dat.emp_comp$Eff_var_affective*dat.emp_comp$Eff_var_cognitive)

dat.emp_comp$Eff<-Eff
dat.emp_comp$Eff_var<-Eff_var
#------------------------------------------------------------------------------------------
#Random effects model for difference in effects:
fit.emp_comp<-rma(yi=Eff, vi=Eff_var, data=dat.emp_comp)
summary(fit.emp_comp)

sink(paste0(model.folder, 'Difference in Empathy Dimension Model - no moderators.txt'))
summary(fit.emp_comp)
sink()

jpeg(paste0(graphics.folder, 'Difference in Empathy Model - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_comp, order = 'obs', slab=dat.emp_comp$citation)
title("Difference in Relation between CU and Empathy Dimensions (Fisher's z)")
dev.off()

jpeg(paste0(graphics.folder, 'Difference in Empathy Dimension Model - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.emp_comp)
title("Funnel Plot of Difference in Relation between CU and Empathy Dimensions")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
REG.emp_comp<-regtest(fit.emp_comp, model = 'lm')
REG.emp_comp

#Trim and fill model
fit.emp_comp.TF_R<-trimfill(fit.emp_comp, estimator = 'R0')
fit.emp_comp.TF_R
funnel(fit.emp_comp.TF_R)

fit.emp_comp.TF_L<-trimfill(fit.emp_comp, estimator = 'L0')
fit.emp_comp.TF_L
funnel(fit.emp_comp.TF_L)

#file drawer analysis
fsn(yi=Eff, vi=Eff_var, data=dat.emp_comp, type='Rosenthal', alpha = .05)
###########################################################################################
#Unconditional Models (no moderators - publication graphics)
#Fitting two summary graphics (one for affective vs. cognitive) & one for summary vals
dat.graph1<-dat[dat$Outcome=='empathy_aff' | dat$Outcome=='empathy_cog',]
dat.graph1$cite.fac<-as.factor(dat.graph1$citation)
fit.graph1<-rma(yi=Eff, vi=Eff_var, data=dat.graph1)
par(mar=c(4,4,1,2))

jpeg(paste0(graphics.folder,'Affective vs. Cognitive Empathy.jpeg'), res=300, units = 'in', height = 8.5, width=11)
par(cex=.75, font=1)
forest(fit.graph1, xlim=c(-10, 2), 
       order = order(dat.graph1$Outcome), 
       ilab = cbind(dat.graph1$N, 
                    dat.graph1$female, 
                    dat.graph1$age, 
                    dat.graph1$CU_resp.R, 
                    dat.graph1$Out_resp.R),
       ilab.xpos = c(-7, -6.5, -5.5, -4.5, -3), 
       rows = c(3:20, 25:43),
       ylim = c(-1,47), 
       cex=.75, 
       xlab="Fisher's z", mlab="", 
       addfit = F,
       slab = dat.graph1$cite.fac)

par(cex=.75, font=4)
text(-10, c(21, 44), pos=4, c('Affective Empathy', 'Cogntive Empathy'))

par(font=4)
text(-7, 46, 'N')
par(font=2)
text(c(-6.5, -5.5, -4.5, -3), 46, c('%Female', 'Mean Age', 'CU Respondent', 'Outcome Respondent'))
text(-10, 46, 'Citation', pos=4)

addpoly(fit.emp_aff, row=1.5, cex=.7, mlab = "")
addpoly(fit.emp_cog, row=23.5, cex=.7, mlab = "")

text(-9, 1.5, 
     pos=4, 
     cex=0.75,
     bquote(
       paste("RE Model for Affective Empathy (Q = ", 
             .(formatC(fit.emp_aff$QE, digits=2, format="f")), 
             ", df = ", 
             .(fit.emp_aff$k - fit.emp_aff$p),
             ", p = ", 
             .(formatC(fit.emp_aff$QEp, digits=2, format="f")), 
             "; ", I^2, " = ",
             .(formatC(fit.emp_aff$I2, digits=1, format="f")), "%)")))

text(-9, 23.5, 
     pos=4, 
     cex=0.75,
     bquote(
       paste("RE Model for Cognitive Empathy (Q = ", 
             .(formatC(fit.emp_cog$QE, digits=2, format="f")), 
             ", df = ", 
             .(fit.emp_cog$k - fit.emp_cog$p),
             ", p = ", 
             .(formatC(fit.emp_cog$QEp, digits=2, format="f")), 
             "; ", I^2, " = ",
             .(formatC(fit.emp_cog$I2, digits=1, format="f")), "%)")))
dev.off()

#Unconditional Models (no moderators - publication graphics)
#Fitting two summary graphics (one for affective vs. cognitive) & one for summary vals
dat.graph2<-dat[dat$Outcome!='empathy_aff',]
dat.graph2<-dat.graph2[dat.graph2$Outcome!='empathy_cog',]

dat.graph2$cite.fac<-as.factor(dat.graph2$citation)
fit.graph2<-rma(yi=Eff, vi=Eff_var, data=dat.graph2)

#round down age and % female for purposes of plotting
dat.graph2$female<-round(dat.graph2$female, digits = 2)
dat.graph2$age<-round(dat.graph2$age, digits = 2)

par(mar=c(4,4,1,2))

jpeg(paste0(graphics.folder,'Main Effects summary.jpeg'), res=300, units = 'in', height = 8.5, width=11)
par(cex=.75, font=1)
forest(fit.graph2, xlim=c(-10, 2), 
       order = order(dat.graph2$Outcome), 
       ilab = cbind(dat.graph2$N, 
                    dat.graph2$female, 
                    dat.graph2$age, 
                    dat.graph2$CU_resp.R, 
                    dat.graph2$Out_resp.R),
       ilab.xpos = c(-7, -6.5, -5.5, -4.5, -3), 
       rows = c(3:27, 32:34, 39:54),
       ylim = c(-1,58), 
       cex=.75, 
       xlab="Fisher's z", mlab="", 
       addfit = F,
       slab = dat.graph2$cite.fac)

par(cex=.75, font=4)
text(-10, c(28, 35, 55), pos=4, c('Total Empathy', 'Guilt', 'Prosociality'))

par(font=4)
text(-7, 57, 'N')
par(font=2)
text(c(-6.5, -5.5, -4.5, -3), 57, c('%Female', 'Mean Age', 'CU Respondent', 'Outcome Respondent'))
text(-10, 57, 'Citation', pos=4)

addpoly(fit.prosoc, row=1.5, cex=.7, mlab = "")
addpoly(fit.glt, row=30.5, cex=.7, mlab = "")
addpoly(fit.prosoc, row=37.5, cex=.7, mlab = "")

text(-9, 1.5, 
     pos=4, 
     cex=0.75,
     bquote(
       paste("RE Model for Total Empathy (Q = ", 
             .(formatC(fit.emp_tot$QE, digits=2, format="f")), 
             ", df = ", 
             .(fit.emp_tot$k - fit.emp_tot$p),
             ", p = ", 
             .(formatC(fit.emp_tot$QEp, digits=2, format="f")), 
             "; ", I^2, " = ",
             .(formatC(fit.emp_tot$I2, digits=1, format="f")), "%)")))

text(-9, 30.5, 
     pos=4, 
     cex=0.75,
     bquote(
       paste("FE Model for Guilt (Q = ", 
             .(formatC(fit.glt$QE, digits=2, format="f")), 
             ", df = ", 
             .(fit.glt$k - fit.glt$p),
             ", p = ", 
             .(formatC(fit.glt$QEp, digits=2, format="f")), 
             "; ", I^2, " = ",
             .(formatC(fit.glt$I2, digits=1, format="f")), "%)")))

text(-9, 37.5, 
     pos=4, 
     cex=0.75,
     bquote(
       paste("RE Model for Prosociality (Q = ", 
             .(formatC(fit.prosoc$QE, digits=2, format="f")), 
             ", df = ", 
             .(fit.prosoc$k - fit.prosoc$p),
             ", p = ", 
             .(formatC(fit.prosoc$QEp, digits=2, format="f")), 
             "; ", I^2, " = ",
             .(formatC(fit.prosoc$I2, digits=1, format="f")), "%)")))

dev.off()

#######################################################################################
#Moderation analyses
fit.emp_tot.mod_age<-rma(yi=Eff, 
                          vi=Eff_var, 
                          mods = ~age, 
                          data=dat.Emp_tot, 
                          ni=N, 
                          knha=T)
summary(fit.emp_tot.mod_age)

fit.emp_tot.mod_female<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~female, 
                         data=dat.Emp_tot, 
                         ni=N, 
                         knha=T)
summary(fit.emp_tot.mod_female)

fit.emp_tot.mod_demo<-rma(yi=Eff, 
                            vi=Eff_var, 
                            mods = ~female+age, 
                            data=dat.Emp_tot, 
                            ni=N, 
                            knha=T)
summary(fit.emp_tot.mod_demo)

fit.emp_tot.mod_CU<-rma(yi=Eff, 
                        vi=Eff_var, 
                        mods = ~CU_resp, 
                        data=dat.Emp_tot, 
                        ni=N, 
                        knha=T)
summary(fit.emp_tot.mod_CU)

fit.emp_tot.mod_CU.all<-rma(yi=Eff, 
                        vi=Eff_var, 
                        mods = ~female+age+CU_resp, 
                        data=dat.Emp_tot, 
                        ni=N, 
                        knha=T)
summary(fit.emp_tot.mod_CU.all)

fit.emp_tot.mod_Out<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~Out_resp, 
                         data=dat.Emp_tot, 
                         ni=N, 
                         knha=T)
summary(fit.emp_tot.mod_Out)

fit.emp_tot.mod_Out.all<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~female+age+Out_resp, 
                         data=dat.Emp_tot, 
                         ni=N, 
                         knha=T)
summary(fit.emp_tot.mod_Out.all)

#None of the moderators were significant predictors... 

#------------------------------------------------------------------------------------------
#Prosocial Model:
fit.prosoc.mod_age<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~age, 
                         data=dat.prosoc, 
                         ni=N, 
                         knha=T)
summary(fit.prosoc.mod_age)

fit.prosoc.mod_female<-rma(yi=Eff, 
                            vi=Eff_var, 
                            mods = ~female, 
                            data=dat.prosoc, 
                            ni=N, 
                            knha=T)
summary(fit.prosoc.mod_female)

fit.prosoc.mod_demo<-rma(yi=Eff, 
                          vi=Eff_var, 
                          mods = ~female+age, 
                          data=dat.prosoc, 
                          ni=N, 
                          knha=T)
summary(fit.prosoc.mod_demo)

fit.prosoc.mod_CU<-rma(yi=Eff, 
                        vi=Eff_var, 
                        mods = ~CU_resp, 
                        data=dat.prosoc, 
                        ni=N, 
                        knha=T)
summary(fit.prosoc.mod_CU)

fit.prosoc.mod_CU.all<-rma(yi=Eff, 
                            vi=Eff_var, 
                            mods = ~female+age+CU_resp, 
                            data=dat.prosoc, 
                            ni=N, 
                            knha=T)
summary(fit.prosoc.mod_CU.all)

fit.prosoc.mod_Out<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~Out_resp, 
                         data=dat.prosoc, 
                         ni=N, 
                         knha=T)
summary(fit.prosoc.mod_Out)

fit.prosoc.mod_Out.all<-rma(yi=Eff, 
                             vi=Eff_var, 
                             mods = ~female+age+Out_resp, 
                             data=dat.prosoc, 
                             ni=N, 
                             knha=T)
summary(fit.prosoc.mod_Out.all)
#------------------------------------------------------------------------------------------
#Affective Empathy Model:
fit.emp_aff.mod_age<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~age, 
                         data=dat.Emp_aff, 
                         ni=N, 
                         knha=T)
summary(fit.emp_aff.mod_age)

fit.emp_aff.mod_female<-rma(yi=Eff, 
                            vi=Eff_var, 
                            mods = ~female, 
                            data=dat.Emp_aff, 
                            ni=N, 
                            knha=T)
summary(fit.emp_aff.mod_female)

fit.emp_aff.mod_demo<-rma(yi=Eff, 
                          vi=Eff_var, 
                          mods = ~female+age, 
                          data=dat.Emp_aff, 
                          ni=N, 
                          knha=T)
summary(fit.emp_aff.mod_demo)

fit.emp_aff.mod_CU<-rma(yi=Eff, 
                        vi=Eff_var, 
                        mods = ~CU_resp, 
                        data=dat.Emp_aff, 
                        ni=N, 
                        knha=T)
summary(fit.emp_aff.mod_CU)

fit.emp_aff.mod_CU.all<-rma(yi=Eff, 
                            vi=Eff_var, 
                            mods = ~female+age+CU_resp, 
                            data=dat.Emp_aff, 
                            ni=N, 
                            knha=T)
summary(fit.emp_aff.mod_CU.all)

fit.emp_aff.mod_Out<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~Out_resp, 
                         data=dat.Emp_aff, 
                         ni=N, 
                         knha=T)
summary(fit.emp_aff.mod_Out)

fit.emp_aff.mod_Out.all<-rma(yi=Eff, 
                             vi=Eff_var, 
                             mods = ~female+age+Out_resp, 
                             data=dat.Emp_aff, 
                             ni=N, 
                             knha=T)
summary(fit.emp_aff.mod_Out.all)

#------------------------------------------------------------------------------------------
#Cognitive Empathy Model:
fit.emp_cog.mod_age<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~age, 
                         data=dat.Emp_cog, 
                         ni=N, 
                         knha=T)
summary(fit.emp_cog.mod_age)

fit.emp_cog.mod_female<-rma(yi=Eff, 
                            vi=Eff_var, 
                            mods = ~female, 
                            data=dat.Emp_cog, 
                            ni=N, 
                            knha=T)
summary(fit.emp_cog.mod_female)

fit.emp_cog.mod_demo<-rma(yi=Eff, 
                          vi=Eff_var, 
                          mods = ~female+age, 
                          data=dat.Emp_cog, 
                          ni=N, 
                          knha=T)
summary(fit.emp_cog.mod_demo)

fit.emp_cog.mod_CU<-rma(yi=Eff, 
                        vi=Eff_var, 
                        mods = ~CU_resp, 
                        data=dat.Emp_cog, 
                        ni=N, 
                        knha=T)
summary(fit.emp_cog.mod_CU)

fit.emp_cog.mod_CU.all<-rma(yi=Eff, 
                            vi=Eff_var, 
                            mods = ~female+age+CU_resp, 
                            data=dat.Emp_cog, 
                            ni=N, 
                            knha=T)
summary(fit.emp_cog.mod_CU.all)

fit.emp_cog.mod_Out<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~Out_resp, 
                         data=dat.Emp_cog, 
                         ni=N, 
                         knha=T)
summary(fit.emp_cog.mod_Out)

fit.emp_cog.mod_Out.all<-rma(yi=Eff, 
                             vi=Eff_var, 
                             mods = ~female+age+Out_resp, 
                             data=dat.Emp_cog, 
                             ni=N, 
                             knha=T)
summary(fit.emp_cog.mod_Out.all)

#------------------------------------------------------------------------------------------
#Empathy Difference Model:
fit.emp_comp.mod_age<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~age, 
                         data=dat.emp_comp, 
                         ni=N, 
                         knha=T)
summary(fit.emp_comp.mod_age)

fit.emp_comp.mod_female<-rma(yi=Eff, 
                            vi=Eff_var, 
                            mods = ~female, 
                            data=dat.emp_comp, 
                            ni=N, 
                            knha=T)
summary(fit.emp_comp.mod_female)

fit.emp_comp.mod_demo<-rma(yi=Eff, 
                          vi=Eff_var, 
                          mods = ~female+age, 
                          data=dat.emp_comp, 
                          ni=N, 
                          knha=T)
summary(fit.emp_comp.mod_demo)

fit.emp_comp.mod_CU<-rma(yi=Eff, 
                        vi=Eff_var, 
                        mods = ~CU_resp, 
                        data=dat.emp_comp, 
                        ni=N, 
                        knha=T)
summary(fit.emp_comp.mod_CU)

fit.emp_comp.mod_CU.all<-rma(yi=Eff, 
                            vi=Eff_var, 
                            mods = ~female+age+CU_resp, 
                            data=dat.emp_comp, 
                            ni=N, 
                            knha=T)
summary(fit.emp_comp.mod_CU.all)

fit.emp_comp.mod_Out<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~Out_resp, 
                         data=dat.emp_comp, 
                         ni=N, 
                         knha=T)
summary(fit.emp_comp.mod_Out)

fit.emp_comp.mod_Out.all<-rma(yi=Eff, 
                             vi=Eff_var, 
                             mods = ~female+age+Out_resp, 
                             data=dat.emp_comp, 
                             ni=N, 
                             knha=T)
summary(fit.emp_comp.mod_Out.all)
