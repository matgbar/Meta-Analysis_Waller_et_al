#################################################################################
#Meta Analysis - Callous and unemotional traits

#Title: Opposite sides of the same coin: A meta-analysis... 

#Authors: Waller et al. 

#Analyses by Matthew Barstead (contact: barstead@umd.edu)

#################################################################################
library(metafor)
library(stargazer)
library(ggplot2)
library(ggridges)
library(mice)
library(bayesplot)
#################################################################################
#For Windows Laptop
user<-Sys.getenv('USERPROFILE')
data.folder<-paste0(user, '/Box Sync/CU meta-analysis/Meta_Raw_Data/')
model.folder<-paste0(user, '/Box Sync/CU meta-analysis/Output/')
graphics.folder<-paste0(user, '/Box Sync/CU meta-analysis/Graphics_Folder/' )

#################################################################################
#For Linux workstation
wd<-'/home/mbarsted/GitHub/Meta-Analysis_Waller_et_al/'
data.folder<-'/home/mbarsted/GitHub/Meta-Analysis_Waller_et_al/Meta_Raw_Data/'
dir.create(data.folder)
model.folder<-'/home/mbarsted/GitHub/Meta-Analysis_Waller_et_al/Output/'
dir.create(model.folder)
graphics.folder<-'/home/mbarsted/GitHub/Meta-Analysis_Waller_et_al/Graphics_Folder/'
dir.create(graphics.folder)

#################################################################################
#Notes:   1. Correlations were transformed to Fisher's z
#         2. Standardized differences were transformed to r then to Fisher's z
#         3. Will require all results are exponentiated to return to -1:1 scale
#################################################################################

#################################################################################
#Analytic Notes: 
#   Based on feedback from JH, decided to take psychometric meta-anlatyic approach
#   Invovles correcting effect sizes & variances for attenuation 
#   See for more details: 
#     Schmidt, F. L., & Hunter, J. E. (2015). Methods of meta-analysis: Correcting 
#       error and bias in research findings (3rd ed.). Washington, DC: SAGE.  
#################################################################################

#################################################################################
#Or can load latest data from GitHub
load(paste0(wd, '/Meta_Raw_Data/Final_Data_110918.RData'))

#Imputing missing reliabilities prior to recoding. 
dat$R_var<-((1-dat$R^2)^2)/(dat$N-1)
psych::describe(dat)

#Need to create a series of binary variables - cannot handle categorical variables
dat$Emp_tot<-ifelse(dat$Outcome=='empathy_tot', 1, 0)
dat$Emp_aff<-ifelse(dat$Outcome=='empathy_aff', 1, 0)
dat$Emp_cog<-ifelse(dat$Outcome=='empathy_cog', 1, 0)
dat$prosocial<-ifelse(dat$Outcome=='prosocial', 1, 0)
dat$guilt<-ifelse(dat$Outcome=='guilt', 1, 0)

#selecting total empathy as reference value - largest N
dat.imp<-dat[c('id', 
               'female', 
               'age', 
               'N', 
               'Emp_aff', 
               'Emp_cog', 
               'prosocial',
               'guilt',
               'R', 
               'R_var', 
               'CU_resp', 
               'Out_resp', 
               'Samp_typ', 
               'CU_Rel', 
               'Out_Rel')]

#Imputation Model - imputing missing reliabilities before the calculation of attenuation corrected scores
imp<-mice(dat.imp, maxit=100, m=5)

png(paste0(graphics.folder, 'MICE_Imputation_convergence.png'), 
    res=300, 
    units='in', 
    height = 4, 
    width = 6)
plot(imp)
dev.off()

#Extracting and obtaining estimates for missing values (averaging across the 5 imputed datasets)
dat.imp<-complete(imp)
dat.imp$citation<-dat$citation
dat.imp<-merge(dat.imp, dat.descrip[,c(1,10)], by='id', all=T)

#Note missing two values from Pechorro et al., regarding ICU use
#Adding those back in below (don't need this step when re-running analyses - just load .RData file)
dat.imp$ICU[85]<-0
dat.imp$ICU[86]<-0

#some cleanup and recoding
dat.imp$CU_resp.R[dat.imp$CU_resp==0]<-'Self'
dat.imp$CU_resp.R[dat.imp$CU_resp==1]<-'Other'
dat.imp$Out_resp.R[dat.imp$Out_resp==0]<-'Self'
dat.imp$Out_resp.R[dat.imp$Out_resp==1]<-'Other'
dat.imp$Same_diff[dat.imp$CU_resp.R==dat.imp$Out_resp.R]<-0
dat.imp$Same_diff[dat.imp$CU_resp.R!=dat.imp$Out_resp.R]<-1

dat.imp$Same_diff.R[dat.imp$Same_diff==0]<-'Same'
dat.imp$Same_diff.R[dat.imp$Same_diff==1]<-'Different'

dat.imp$Samp_typ.R[dat.imp$Samp_typ==1]<-'Referred/Clinical'
dat.imp$Samp_typ.R[dat.imp$Samp_typ==0]<-'Community/Non-clincal'

#need to correct d's first, then I can convert them to r's that can be analyzed
#r's will have their own conversion below (were initially correlations in the report)
dat<-dat[order(dat$id),]
dat.imp<-dat.imp[order(dat.imp$id),]

dat.imp$d_cor[!is.na(dat$d)]<-dat$d[!is.na(dat$d)]*sqrt(dat.imp$Out_Rel[!is.na(dat$d)])

#Calculating variances specific to d's (will not assume equality of population variances)
#Formula from Schmidt & Hunter (2015) p. 293
dat.imp$Eff_var_d[!is.na(dat$d)]<-
  (dat$N_CU[!is.na(dat$d)]+dat$N_ctrl[!is.na(dat$d)])/(dat$N_CU[!is.na(dat$d)]*dat$N_ctrl[!is.na(dat$d)])+
  dat.imp$d_cor[!is.na(dat$d)]^2/(2*(dat$N_CU[!is.na(dat$d)]+dat$N_ctrl[!is.na(dat$d)]))
  
dat.imp$Eff_var_d[!is.na(dat$d)]<-dat.imp$Eff_var_d[!is.na(dat$d)]/dat.imp$Out_Rel[!is.na(dat$d)]

#The conversion of corrected d's to corrected r's 
#Formulas taken from Borenstein et al. (2009) pp. 48-49
a<-(dat$N_CU+dat$N_ctrl)^2/(dat$N_CU*dat$N_ctrl)
dat.imp$Eff[!is.na(dat$d)]<-dat.imp$d_cor[!is.na(dat$d)]/sqrt(dat.imp$d_cor[!is.na(dat$d)]^2+a[!is.na(dat$d)])

#Converting these variance estimates to r scale
dat.imp$Eff_var[!is.na(dat$d)]<-
  (a[!is.na(dat$d)]^2*dat.imp$Eff_var_d[!is.na(dat$d)])/(dat.imp$d_cor[!is.na(dat$d)]^2+a[!is.na(dat$d)])^3

#Creating Attenuation-corrected values for raw correlations
dat.imp$Eff[is.na(dat$d)]<-dat.imp$R[is.na(dat$d)]/sqrt(dat.imp$CU_Rel[is.na(dat$d)]*dat.imp$Out_Rel[is.na(dat$d)])
hist(dat.imp$Eff)

#Adding variances for what were originally correlations. 
dat.imp$Eff_var[is.na(dat$d)]<-dat.imp$R_var[is.na(dat$d)]/sqrt(dat.imp$CU_Rel[is.na(dat$d)]*dat.imp$Out_Rel[is.na(dat$d)])

#This is the new analytic data set
sink(paste0(model.folder, 'Imputed_Descriptives.txt'))
psych::describe(dat.imp)
sink()

dat.Emp_tot<-dat.imp[rowSums(dat.imp[c('Emp_aff',
                                       'Emp_cog', 
                                       'guilt', 
                                       'prosocial')]
                             )==0,]

dat.Emp_aff<-dat.imp[dat.imp$Emp_aff==1,]
dat.Emp_cog<-dat.imp[dat.imp$Emp_cog==1,] 
dat.glt<-dat.imp[dat.imp$guilt==1,]   #note probably have too few for guilt
dat.prosoc<-dat.imp[dat.imp$prosocial==1,]

#################################################################################
#Applying Attenuation Correction to Effects and Variances - based on JH's reccomendations
#Removing antiquated use of fail-safe N 
#Using metafor to achieve Hunter-Schmidt approach: 
#see: http://www.metafor-project.org/doku.php/tips:hunter_schmidt_method
#################################################################################

#################################################################################
#Model for Total Empathy: 
fit.emp_tot<-rma(yi=Eff, 
                 vi=Eff_var, 
                 weights = 1/Eff_var, 
                 data=dat.Emp_tot, 
                 ni=N,
                 method = 'HS')

summary(fit.emp_tot)

sink(paste0(model.folder, 'Total_Empathy_AC.txt'))
summary(fit.emp_tot)
sink()

tiff(paste0(graphics.folder, 'CU and Total Empathy - Forest.tiff'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_tot, order = 'obs', slab=dat.Emp_tot$citation)
title("Relation between CU Traits and Total Empathy")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias - needs to be tied to original effects 
#(makes no sense to perform this test on corrected values)
fit.emp_tot_bias<-rma(yi=R, 
                      vi=R_var,
                      data=dat.Emp_tot, 
                      ni=N)

REG.emp_tot<-regtest(fit.emp_tot_bias, model = 'lm')
sink(paste0(model.folder, 'Reg_test_Total_Empathy.txt'))
print(REG.emp_tot)
sink()

fit.emp_tot.TF_L<-trimfill(fit.emp_tot_bias, estimator = 'L0')
fit.emp_tot.TF_L

sink(paste0(model.folder, 'Total_Empathy_L0.txt'))
summary(fit.emp_tot.TF_L)
sink()

#Funnel plot for uncorrected correlations: 
tiff(paste0(graphics.folder, 'CU and Total Empathy - Funnels.tiff'), res = 300, width = 11, height = 8, units = 'in')
par(mfrow=c(1,2))
funnel(fit.emp_tot_bias, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Raw Total Empathy Correlations'))
funnel(fit.emp_tot.TF_L, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Trim and Fill Results for Total Empathy (L'[0]*' Estimator)'))

dev.off()

#################################################################################
#Model for Affective Empathy: 
fit.emp_aff<-rma(yi=Eff, 
                 vi=Eff_var, 
                 weights = 1/Eff_var, 
                 data=dat.Emp_aff, 
                 ni=N,
                 method = 'HS')

summary(fit.emp_aff)

sink(paste0(model.folder, 'Affective_Empathy_AC.txt'))
summary(fit.emp_aff)
sink()

tiff(paste0(graphics.folder, 'CU and Affective Empathy - Forest.tiff'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_aff, order = 'obs', slab=dat.Emp_aff$citation)
title("Relation between CU Traits and Affective Empathy")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
fit.emp_aff_bias<-rma(yi=R, 
                      vi=R_var,
                      data=dat.Emp_aff, 
                      ni=N)

REG.emp_aff<-regtest(fit.emp_aff_bias, model = 'lm')
sink(paste0(model.folder, 'Reg_test_Affective_Empathy.txt'))
print(REG.emp_aff)
sink()

fit.emp_aff.TF_L<-trimfill(fit.emp_aff_bias, estimator = 'L0')
fit.emp_aff.TF_L

sink(paste0(model.folder, 'Affective_Empathy_L0.txt'))
summary(fit.emp_aff.TF_L)
sink()

#Funnel plot for uncorrected correlations: 
tiff(paste0(graphics.folder, 'CU and Affective Empathy - Funnels.tiff'), res = 300, width = 11, height = 8, units = 'in')
par(mfrow=c(1,2))
funnel(fit.emp_aff_bias, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Trim and Fill Results for Affective Empathy Correlations'))
funnel(fit.emp_aff.TF_L, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Trim and Fill Results for Affective Empathy (L'[0]*' Estimator)'))

dev.off()

#################################################################################
#Model for Cognitive Empathy: 
fit.emp_cog<-rma(yi=Eff, 
                 vi=Eff_var, 
                 weights = 1/Eff_var, 
                 data=dat.Emp_cog, 
                 ni=N,
                 method = 'HS')

summary(fit.emp_cog)

sink(paste0(model.folder, 'Cognitive_Empathy_AC.txt'))
summary(fit.emp_cog)
sink()

tiff(paste0(graphics.folder, 'CU and Cognitive Empathy - Forest.tiff'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_cog, order = 'obs', slab=dat.Emp_cog$citation)
title("Relation between CU Traits and Cognitive Empathy")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
fit.emp_cog_bias<-rma(yi=R, 
                      vi=R_var,
                      data=dat.Emp_cog, 
                      ni=N)

REG.emp_cog<-regtest(fit.emp_cog_bias, model = 'lm')
sink(paste0(model.folder, 'Reg_test_Cognitive_Empathy.txt'))
print(REG.emp_cog)
sink()

fit.emp_cog.TF_L<-trimfill(fit.emp_cog_bias, estimator = 'L0')
fit.emp_cog.TF_L

sink(paste0(model.folder, 'Cognitive_Empathy_L0.txt'))
summary(fit.emp_cog.TF_L)
sink()

#Funnel plot for uncorrected correlations: 
tiff(paste0(graphics.folder, 'CU and Cognitive Empathy - Funnels.tiff'), res = 300, width = 11, height = 8, units = 'in')
par(mfrow=c(1,2))
funnel(fit.emp_cog_bias, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Raw Cognitive Empathy Correlations'))
funnel(fit.emp_cog.TF_L, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Trim and Fill Results for Cognitive Empathy (L'[0]*' Estimator)'))

dev.off()

#------------------------------------------------------------------------------------------
#Model for Prosocialty: 
fit.prosoc<-rma(yi=Eff, 
                vi=Eff_var, 
                weights = 1/Eff_var, 
                data=dat.prosoc, 
                ni=N,
                method = 'HS')

summary(fit.prosoc)

sink(paste0(model.folder, 'Prosocial_AC.txt'))
summary(fit.prosoc)
sink()

tiff(paste0(graphics.folder, 'CU and Prosocialty - Forest.tiff'), res=300, width = 7, height=7, units='in')
forest.rma(fit.prosoc, order = 'obs', slab=dat.prosoc$citation)
title("Relation between CU Traits and Prosocialty")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
fit.prosoc_bias<-rma(yi=R, 
                      vi=R_var,
                      data=dat.prosoc, 
                      ni=N)

REG.prosoc<-regtest(fit.prosoc_bias, model = 'lm')
sink(paste0(model.folder, 'Reg_test_Prosociality.txt'))
print(REG.prosoc)
sink()

fit.prosoc.TF_L<-trimfill(fit.prosoc_bias, estimator = 'L0')
fit.prosoc.TF_L

sink(paste0(model.folder, 'Prosociality_L0.txt'))
summary(fit.prosoc.TF_L)
sink()

#Funnel plot for uncorrected correlations: 
tiff(paste0(graphics.folder, 'CU and Prosociality - Funnels.tiff'), res = 300, width = 11, height = 8, units = 'in')
par(mfrow=c(1,2))
funnel(fit.prosoc_bias, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Raw Prosociality Correlations'))
funnel(fit.prosoc.TF_L, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Trim and Fill Results for Prosociality (L'[0]*' Estimator)'))

dev.off()

#################################################################################
#Model for Guilt: 
fit.glt<-rma(yi=Eff, 
             vi=Eff_var, 
             weights = 1/Eff_var, 
             data=dat.glt, 
             ni=N,
             method = 'FE')
summary(fit.glt)

sink(paste0(model.folder, 'Guilt Overall Model - no moderators.txt'))
summary(fit.glt)
sink()

tiff(paste0(graphics.folder, 'CU and Guilt - Forest.tiff'), res=300, width = 7, height=7, units='in')
forest.rma(fit.glt, order = 'obs', slab=dat.glt$citation)
title("Relation between CU Traits and Guilt")
dev.off()

###########################################################################################
#Model comparing two forms of empathy - within study differences in effects
#merging and cleaning data set
dat.emp_comp<-merge(dat.Emp_aff, dat.Emp_cog, by='id')

cols<-c(1,16,2:4,9,25,26,18,19,13,17,34,50,51)
dat.emp_comp<-dat.emp_comp[,cols]
colnames(dat.emp_comp)<-c('id', 'citation', 'female', 'age', 'N', 'R_affective',
                          'Eff_affective', 'Eff_var_affective', 'CU_resp', 'Out_resp',
                          'Sample','ICU', 'R_cognitive', 'Eff_cognitive', 'Eff_var_cognitive')

#Correlations pulled from studies when available, ordered by citation: 
Aff_cog_cor<-c(.063, .063, .35, 0, 0, .49, 0, .5, 0, .76, .32, 0, -.01, .08,0,.45,0,0)

#Calculating difference in effects and variance terms adjusted for outcome correlation
Eff<-dat.emp_comp$Eff_affective-dat.emp_comp$Eff_cognitive
Eff_var<-dat.emp_comp$Eff_var_affective+dat.emp_comp$Eff_var_cognitive-2*Aff_cog_cor*sqrt(dat.emp_comp$Eff_var_affective*dat.emp_comp$Eff_var_cognitive)

#note that I stupidly changed naming conventions here - sorry
#the Eff variable is corrected for attentuation, censored when an effect was < -1, and adjusted for the dependent nature of the correlations where possible
dat.emp_comp$Eff<-Eff
dat.emp_comp$Eff_var<-Eff_var

#------------------------------------------------------------------------------------------
#Random effects model for difference in effects:
fit.emp_comp<-rma(yi=Eff, 
                  vi=Eff_var, 
                  weights = 1/Eff_var, 
                  data=dat.emp_comp, 
                  ni=N,
                  method = 'HS'
)
summary(fit.emp_comp)

sink(paste0(model.folder, 'Difference in Empathy Dimension Model - no moderators.txt'))
summary(fit.emp_comp)
sink()

tiff(paste0(graphics.folder, 'Difference in Empathy Model - Forest.tiff'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_comp, order = 'obs', slab=dat.emp_comp$citation)
title("Difference in Relation between CU and Empathy Dimensions")
dev.off()

#############################################################################
#Extracting model summary for Each constrast
Ints<-c(fit.emp_tot$b, fit.emp_aff$b, fit.emp_cog$b, fit.prosoc$b, fit.glt$b, fit.emp_comp$b)
ses<-c(fit.emp_tot$se, fit.emp_aff$se, fit.emp_cog$se, fit.prosoc$se, fit.glt$se, fit.emp_comp$se)
zvals<-c(fit.emp_tot$zval, fit.emp_aff$zval, fit.emp_cog$zval, fit.prosoc$zval, fit.glt$zval, fit.emp_comp$zval)
pvals<-c(fit.emp_tot$pval, fit.emp_aff$pval, fit.emp_cog$pval, fit.prosoc$pval, fit.glt$pval, fit.emp_comp$pval)
ci.lbs<-c(fit.emp_tot$ci.lb, fit.emp_aff$ci.lb, fit.emp_cog$ci.lb, fit.prosoc$ci.lb, fit.glt$ci.lb, fit.emp_comp$ci.lb)
ci.ubs<-c(fit.emp_tot$ci.ub, fit.emp_aff$ci.ub, fit.emp_cog$ci.ub, fit.prosoc$ci.ub, fit.glt$ci.ub, fit.emp_comp$ci.ub)

#Vectors of Q values & I^2 vals
Q_vals<-c(fit.emp_tot$QE, fit.emp_aff$QE, fit.emp_cog$QE, fit.prosoc$QE, fit.glt$QE, fit.emp_comp$QE)
Q_p<-c(fit.emp_tot$QEp, fit.emp_aff$QEp, fit.emp_cog$QEp, fit.prosoc$QEp, fit.glt$QEp, fit.emp_comp$QEp)
I_2<-c(round(fit.emp_tot$I2, digits = 2), 
       round(fit.emp_aff$I2, digits = 2),
       round(fit.emp_cog$I2, digits = 2), 
       round(fit.prosoc$I2, digits =2),
       round(fit.glt$I2, digits = 2),
       round(fit.emp_comp$I2, digits = 2))

DF.summary<-data.frame(Outcome=c('Total Empathy', 
                                 'Affective Empathy', 
                                 'Cognitive Empathy', 
                                 'Prosocialty', 
                                 'Guilt', 
                                 'Aff-Cog Empathy'),
                       k=c(27, 18, 19, 19, 3, 19),
                       rho = Ints, 
                       se =ses, 
                       zval=zvals, 
                       pval=pvals, 
                       ci.lb=ci.lbs, 
                       ci.ub=ci.ubs,
                       Q_val=Q_vals, 
                       Q_p, 
                       I_2
)

stargazer(DF.summary, 
          out = paste0(model.folder,'Unconditional Models.txt'), 
          summary = F)

#Plotting summary effect sizes 
Emp_tot.dist<-rnorm(100000, 
                    mean=DF.summary$rho[DF.summary$Outcome=='Total Empathy'], 
                    sd=DF.summary$se[DF.summary$Outcome=='Total Empathy'])

Emp_aff.dist<-rnorm(100000, 
                    mean=DF.summary$rho[DF.summary$Outcome=='Affective Empathy'], 
                    sd=DF.summary$se[DF.summary$Outcome=='Affective Empathy'])

Emp_cog.dist<-rnorm(100000, 
                    mean=DF.summary$rho[DF.summary$Outcome=='Cognitive Empathy'], 
                    sd=DF.summary$se[DF.summary$Outcome=='Cognitive Empathy'])

Prosoc.dist<-rnorm(100000, 
                   mean=DF.summary$rho[DF.summary$Outcome=='Prosocialty'], 
                   sd=DF.summary$se[DF.summary$Outcome=='Prosocialty'])

Glt.dist<-rnorm(100000, 
                mean=DF.summary$rho[DF.summary$Outcome=='Guilt'], 
                sd=DF.summary$se[DF.summary$Outcome=='Guilt'])

Emp_comp.dist<-rnorm(100000, 
                     mean=DF.summary$rho[DF.summary$Outcome=='Aff-Cog Empathy'], 
                     DF.summary$se[DF.summary$Outcome=='Aff-Cog Empathy'])

DF.plot.dist2<-cbind(Emp_tot.dist, Emp_aff.dist, Emp_cog.dist, Prosoc.dist, Glt.dist, Emp_comp.dist)

colnames(DF.plot.dist2)<-c('Total Empathy', 'Affective Empathy', 'Cognitive Empathy', 'Prosociality', 'Guilt', 'Aff-Cog')

library(bayesplot)

tiff(paste0(graphics.folder, 'Unconditional Model Summary Graphic.tiff'), res=300, units='in', height=7, width=7)
mcmc_areas(DF.plot.dist2, prob=.95)+
  xlab(expression('Population Esimates of'~rho))+
  ggtitle('Simulated Distributions of Effect Sizes')
dev.off()

###########################################################################################
#Unconditional Models (no moderators - publication graphics)
#Fitting two summary graphics (one for affective vs. cognitive) & one for summary vals
#Code for graphic models was adapted from: http://www.metafor-project.org/doku.php/plots:forest_plot_with_subgroups
dat.graph1<-dat.imp
dat.graph1<-dat.graph1[rowSums(dat.graph1[,5:8])==0,]

dat.graph1$citation[dat.graph1$citation=='Anastassiou-Hadjicharalambous & Warden (2008)']<-'Anastassiou-H. & Warden (2008)'
dat.graph1$cite.fac<-as.factor(dat.graph1$citation)
fit.graph1<-rma(yi=Eff, vi=Eff_var, data=dat.graph1)

#Adding back in Outcome as a variable
dat.graph1$female<-round(dat.graph1$female, digits = 2)
dat.graph1$age<-round(dat.graph1$age, digits = 2)
par(mar=c(4,4,1,2))

tiff(paste0(graphics.folder,'Total Empathy.tiff'), res=1200, units = 'in', height = 7, width=10)
par(cex=1, font=1)
forest(fit.graph1, xlim=c(-10, 2), 
       order = order(dat.graph1$Eff), 
       ilab = cbind(dat.graph1$N, 
                    dat.graph1$female, 
                    dat.graph1$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:29),
       ylim = c(-1,32), 
       cex=.75, 
       xlab=expression(~rho), 
       mlab="", 
       addfit = F,
       slab = dat.graph1$cite.fac)

par(cex=1, font=4)
par(font=4)
text(-7, 31, 'N')
par(font=2)
text(c(-5, -3), 31, c('%Female', 'Mean Age'))
text(-10, 31, 'Citation', pos=4)
text(1.25, 31, expression(paste(bolditalic('r'), '[LB,UB]')))

addpoly(fit.emp_tot, row=1.5, cex=.7, mlab = "")

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

title('Corrected Effects for Total Empathy')
dev.off()

#-------------------------------------------------------------------------------
dat.graph2<-dat.imp
dat.graph2<-dat.graph2[dat.graph2$Emp_aff==1,]

dat.graph2$cite.fac<-as.factor(dat.graph2$citation)
fit.graph2<-rma(yi=Eff, vi=Eff_var, data=dat.graph2)

#Adding back in Outcome as a variable
dat.graph2$female<-round(dat.graph2$female, digits = 2)
dat.graph2$age<-round(dat.graph2$age, digits = 2)
par(mar=c(4,4,1,2))

tiff(paste0(graphics.folder,'Affective Empathy.tiff'), res=1200, units = 'in', height = 7, width=10)
par(cex=1, font=1)
forest(fit.graph2, xlim=c(-10, 2), 
       order = order(dat.graph2$Eff), 
       ilab = cbind(dat.graph2$N, 
                    dat.graph2$female, 
                    dat.graph2$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:20),
       ylim = c(-1,23), 
       cex=.75, 
       xlab=expression(~rho), 
       mlab="", 
       addfit = F,
       slab = dat.graph2$cite.fac)

par(cex=1, font=4)
par(font=4)
text(-7, 22, 'N')
par(font=2)
text(c(-5, -3), 22, c('%Female', 'Mean Age'))
text(-10, 22, 'Citation', pos=4)
text(1.25, 22, expression(paste(bolditalic('r'), '[LB,UB]')))

addpoly(fit.emp_aff, row=1.5, cex=.7, mlab = "")

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

title('Corrected Effects for Affective Empathy')
dev.off()

#-------------------------------------------------------------------------------
dat.graph3<-dat.imp
dat.graph3<-dat.graph3[dat.graph3$Emp_cog==1,]

dat.graph3$cite.fac<-as.factor(dat.graph3$citation)
fit.graph3<-rma(yi=Eff, vi=Eff_var, data=dat.graph3)

#Adding back in Outcome as a variable
dat.graph3$female<-round(dat.graph3$female, digits = 2)
dat.graph3$age<-round(dat.graph3$age, digits = 2)
par(mar=c(4,4,1,2))

tiff(paste0(graphics.folder,'Cognitive Empathy.tiff'), res=1200, units = 'in', height = 7, width=10)
par(cex=1, font=1)
forest(fit.graph3, xlim=c(-10, 2), 
       order = order(dat.graph3$Eff), 
       ilab = cbind(dat.graph3$N, 
                    dat.graph3$female, 
                    dat.graph3$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:21),
       ylim = c(-1,24), 
       cex=.75, 
       xlab=expression(~rho), 
       mlab="", 
       addfit = F,
       slab = dat.graph3$cite.fac)

par(cex=1, font=4)
par(font=4)
text(-7, 23, 'N')
par(font=2)
text(c(-5, -3), 23, c('%Female', 'Mean Age'))
text(-10, 23, 'Citation', pos=4)
text(1.25, 23, expression(paste(bolditalic('r'), '[LB,UB]')))

addpoly(fit.emp_cog, row=1.5, cex=.7, mlab = "")

text(-9, 1.5, 
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

title('Corrected Effects for Cognitive Empathy')
dev.off()

#-------------------------------------------------------------------------------
dat.graph4<-dat.imp
dat.graph4<-dat.graph4[dat.graph4$prosocial==1,]

dat.graph4$cite.fac<-as.factor(dat.graph4$citation)
fit.graph4<-rma(yi=Eff, vi=Eff_var, data=dat.graph4)

#Adding back in Outcome as a variable
dat.graph4$female<-round(dat.graph4$female, digits = 2)
dat.graph4$age<-round(dat.graph4$age, digits = 2)
par(mar=c(4,4,1,2))

tiff(paste0(graphics.folder,'Prosociality.tiff'), res=1200, units = 'in', height = 7, width=10)
par(cex=1, font=1)
forest(fit.graph4, xlim=c(-10, 2), 
       order = order(dat.graph4$Eff), 
       ilab = cbind(dat.graph4$N, 
                    dat.graph4$female, 
                    dat.graph4$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:21),
       ylim = c(-1,24), 
       cex=.75, 
       xlab=expression(~rho), 
       mlab="", 
       addfit = F,
       slab = dat.graph4$cite.fac)

par(cex=1, font=4)
par(font=4)
text(-7, 23, 'N')
par(font=2)
text(c(-5, -3), 23, c('%Female', 'Mean Age'))
text(-10, 23, 'Citation', pos=4)
text(1.25, 23, expression(paste(bolditalic('r'), '[LB,UB]')))

addpoly(fit.prosoc, row=1.5, cex=.7, mlab = "")

text(-9, 1.5, 
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

title('Corrected Effects for Prosociality')
dev.off()

#-------------------------------------------------------------------------------
dat.graph5<-dat.imp
dat.graph5<-dat.graph5[dat.graph5$guilt==1,]

dat.graph5$cite.fac<-as.factor(dat.graph5$citation)
fit.graph5<-rma(yi=Eff, vi=Eff_var, data=dat.graph5)

#Adding back in Outcome as a variable
dat.graph5$female<-round(dat.graph5$female, digits = 2)
dat.graph5$age<-round(dat.graph5$age, digits = 2)
par(mar=c(4,4,1,2))

tiff(paste0(graphics.folder,'Guilt.tiff'), res=1200, units = 'in', height = 3.5, width=10)
par(cex=1, font=1)
forest(fit.graph5, xlim=c(-10, 2), 
       order = order(dat.graph5$Eff), 
       ilab = cbind(dat.graph5$N, 
                    dat.graph5$female, 
                    dat.graph5$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:5),
       ylim = c(-1,8), 
       cex=.75, 
       xlab=expression(~rho), 
       mlab="", 
       addfit = F,
       slab = dat.graph5$cite.fac)

par(cex=1, font=4)
par(font=4)
text(-7, 7, 'N')
par(font=2)
text(c(-5, -3), 7, c('%Female', 'Mean Age'))
text(-10, 7, 'Citation', pos=4)
text(1.25, 7, expression(paste(bolditalic('r'), '[LB,UB]')))

addpoly(fit.glt, row=1.5, cex=.7, mlab = "")

text(-9, 1.5, 
     pos=4, 
     cex=0.75,
     bquote(
       paste("RE Model for Guilt (Q = ", 
             .(formatC(fit.glt$QE, digits=2, format="f")), 
             ", df = ", 
             .(fit.glt$k - fit.glt$p),
             ", p = ", 
             .(formatC(fit.glt$QEp, digits=2, format="f")), 
             "; ", I^2, " = ",
             .(formatC(fit.glt$I2, digits=1, format="f")), "%)")))

title('Corrected Effects for Guilt')
dev.off()

#######################################################################################
#Moderation analyses
#######################################################################################
#Total Empathy
fit.emp_tot.mod_all<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var, 
                         mods = ~female+age+Samp_typ+CU_resp+Out_resp+ICU, 
                         data=dat.Emp_tot, 
                         ni=N, 
                         knha=T, 
                         method = 'HS')
summary(fit.emp_tot.mod_all)

#Prosocial Model:
table(dat.prosoc$CU_resp, dat.prosoc$Out_resp)
#Perfectly equal - just pick one
fit.prosoc.mod_all<-rma(yi=Eff, 
                        vi=Eff_var,
                        weights = 1/Eff_var,
                        mods = ~female+age+Samp_typ+CU_resp+ICU, 
                        data=dat.prosoc, 
                        ni=N, 
                        knha=T, 
                        method = 'HS')
summary(fit.prosoc.mod_all)  

#Affective Empathy:
table(dat.Emp_aff$CU_resp, dat.Emp_aff$Out_resp)
fit.Emp_aff.mod_all<-rma(yi=Eff, 
                        vi=Eff_var,
                        weights = 1/Eff_var,
                        mods = ~female+age+Samp_typ+CU_resp+Out_resp+ICU, 
                        data=dat.Emp_aff, 
                        ni=N, 
                        knha=T, 
                        method = 'HS')
summary(fit.Emp_aff.mod_all)  

#Cognitive Empathy Model:
table(dat.Emp_cog$CU_resp, dat.Emp_cog$Out_resp)
fit.Emp_cog.mod_all<-rma(yi=Eff, 
                        vi=Eff_var,
                        weights = 1/Eff_var,
                        mods = ~female+age+Samp_typ+CU_resp+Out_resp+ICU, 
                        data=dat.Emp_cog, 
                        ni=N, 
                        knha=T, 
                        method = 'HS')
summary(fit.Emp_cog.mod_all)  

##############################################################################################
#Exploring single moderators now - Female
fit.emp_tot.mod_female<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var, 
                         mods = ~female, 
                         data=dat.Emp_tot, 
                         ni=N, 
                         knha=T, 
                         method = 'HS')
summary(fit.emp_tot.mod_female)

#Prosocial Model:
fit.prosoc.mod_female<-rma(yi=Eff, 
                        vi=Eff_var,
                        weights = 1/Eff_var,
                        mods = ~female, 
                        data=dat.prosoc, 
                        ni=N, 
                        knha=T, 
                        method = 'HS')
summary(fit.prosoc.mod_female)  

#Affective Empathy:
fit.Emp_aff.mod_female<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~female, 
                         data=dat.Emp_aff, 
                         ni=N, 
                         knha=T, 
                         method = 'HS')
summary(fit.Emp_aff.mod_female)  

#Cognitive Empathy Model:
fit.Emp_cog.mod_female<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~female, 
                         data=dat.Emp_cog, 
                         ni=N, 
                         knha=T, 
                         method = 'HS')
summary(fit.Emp_cog.mod_female)  

#############################################################################################
#Exploring Moderators - age only
fit.emp_tot.mod_age<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var, 
                         mods = ~age, 
                         data=dat.Emp_tot, 
                         ni=N, 
                         knha=T, 
                         method = 'HS')
summary(fit.emp_tot.mod_age)

#Prosocial Model:
fit.prosoc.mod_age<-rma(yi=Eff, 
                        vi=Eff_var,
                        weights = 1/Eff_var,
                        mods = ~age, 
                        data=dat.prosoc, 
                        ni=N, 
                        knha=T, 
                        method = 'HS')
summary(fit.prosoc.mod_age)  

#Affective Empathy:
fit.Emp_aff.mod_age<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~age, 
                         data=dat.Emp_aff, 
                         ni=N, 
                         knha=T, 
                         method = 'HS')
summary(fit.Emp_aff.mod_age)  

#Cognitive Empathy Model:
fit.Emp_cog.mod_age<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~age, 
                         data=dat.Emp_cog, 
                         ni=N, 
                         knha=T, 
                         method = 'HS')
summary(fit.Emp_cog.mod_age)  

#######################################################################################
#Exploring moderators - Sample Type
fit.emp_tot.mod_samp<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var, 
                         mods = ~Samp_typ, 
                         data=dat.Emp_tot, 
                         ni=N, 
                         knha=T, 
                         method = 'HS')
summary(fit.emp_tot.mod_samp)

#Prosocial Model:
fit.prosoc.mod_samp<-rma(yi=Eff, 
                        vi=Eff_var,
                        weights = 1/Eff_var,
                        mods = ~Samp_typ, 
                        data=dat.prosoc, 
                        ni=N, 
                        knha=T, 
                        method = 'HS')
summary(fit.prosoc.mod_samp)  

#Affective Empathy:
fit.Emp_aff.mod_samp<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~Samp_typ, 
                         data=dat.Emp_aff, 
                         ni=N, 
                         knha=T, 
                         method = 'HS')
summary(fit.Emp_aff.mod_samp)  

#Cognitive Empathy Model:
fit.Emp_cog.mod_samp<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~Samp_typ, 
                         data=dat.Emp_cog, 
                         ni=N, 
                         knha=T, 
                         method = 'HS')
summary(fit.Emp_cog.mod_samp)

#######################################################################################
#Exploring moderators - CU_respondent
fit.emp_tot.mod_CU_resp<-rma(yi=Eff, 
                          vi=Eff_var,
                          weights = 1/Eff_var, 
                          mods = ~CU_resp, 
                          data=dat.Emp_tot, 
                          ni=N, 
                          knha=T, 
                          method = 'HS')
summary(fit.emp_tot.mod_CU_resp)

#Prosocial Model:
fit.prosoc.mod_CU_resp<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~CU_resp, 
                         data=dat.prosoc, 
                         ni=N, 
                         knha=T, 
                         method = 'HS')
summary(fit.prosoc.mod_CU_resp)  

#Affective Empathy:
fit.Emp_aff.mod_CU_resp<-rma(yi=Eff, 
                          vi=Eff_var,
                          weights = 1/Eff_var,
                          mods = ~CU_resp, 
                          data=dat.Emp_aff, 
                          ni=N, 
                          knha=T, 
                          method = 'HS')
summary(fit.Emp_aff.mod_CU_resp)  

#Cognitive Empathy Model:
fit.Emp_cog.mod_CU_resp<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~CU_resp, 
                         data=dat.Emp_cog, 
                         ni=N, 
                         knha=T, 
                         method = 'HS')
summary(fit.Emp_cog.mod_CU_resp)  

#######################################################################################
#Exploring moderators - Outcome respondent
fit.emp_tot.mod_Out_resp<-rma(yi=Eff, 
                             vi=Eff_var,
                             weights = 1/Eff_var, 
                             mods = ~Out_resp, 
                             data=dat.Emp_tot, 
                             ni=N, 
                             knha=T, 
                             method = 'HS')
summary(fit.emp_tot.mod_Out_resp)

#Prosocial Model: - exactly the same as the CU respondent model
fit.prosoc.mod_Out_resp<-rma(yi=Eff, 
                            vi=Eff_var,
                            weights = 1/Eff_var,
                            mods = ~Out_resp, 
                            data=dat.prosoc, 
                            ni=N, 
                            knha=T, 
                            method = 'HS')
summary(fit.prosoc.mod_Out_resp)  

#Affective Empathy:
fit.Emp_aff.mod_Out_resp<-rma(yi=Eff, 
                             vi=Eff_var,
                             weights = 1/Eff_var,
                             mods = ~Out_resp, 
                             data=dat.Emp_aff, 
                             ni=N, 
                             knha=T, 
                             method = 'HS')
summary(fit.Emp_aff.mod_Out_resp)  

#Cognitive Empathy Model:
fit.Emp_cog.mod_Out_resp<-rma(yi=Eff, 
                             vi=Eff_var,
                             weights = 1/Eff_var,
                             mods = ~Out_resp, 
                             data=dat.Emp_cog, 
                             ni=N, 
                             knha=T, 
                             method = 'HS')
summary(fit.Emp_cog.mod_Out_resp)  

#######################################################################################
#Exploring moderators - ICU respondent
fit.emp_tot.mod_ICU<-rma(yi=Eff, 
                              vi=Eff_var,
                              weights = 1/Eff_var, 
                              mods = ~ICU, 
                              data=dat.Emp_tot, 
                              ni=N, 
                              knha=T, 
                              method = 'HS')
summary(fit.emp_tot.mod_ICU)

#Prosocial Model: - exactly the same as the CU respondent model
fit.prosoc.mod_ICU<-rma(yi=Eff, 
                             vi=Eff_var,
                             weights = 1/Eff_var,
                             mods = ~ICU, 
                             data=dat.prosoc, 
                             ni=N, 
                             knha=T, 
                             method = 'HS')
summary(fit.prosoc.mod_ICU)  

#Affective Empathy:
fit.Emp_aff.mod_ICU<-rma(yi=Eff, 
                              vi=Eff_var,
                              weights = 1/Eff_var,
                              mods = ~ICU, 
                              data=dat.Emp_aff, 
                              ni=N, 
                              knha=T, 
                              method = 'HS')
summary(fit.Emp_aff.mod_ICU)  

#Cognitive Empathy Model:
fit.Emp_cog.mod_ICU<-rma(yi=Eff, 
                              vi=Eff_var,
                              weights = 1/Eff_var,
                              mods = ~ICU, 
                              data=dat.Emp_cog, 
                              ni=N, 
                              knha=T, 
                              method = 'HS')
summary(fit.Emp_cog.mod_ICU)  

#############################################################################################################################
#Graphing moderation
#Exploring moderation graphically involving dichotomous ICU variable (Cognitive Empathy Model)
ICU_resp<-c(0,1)
emp_aff.ICU.pred<-predict(fit.Emp_aff.mod_ICU, newmods = ICU_resp, level = 95)

ICU_self<-emp_aff.ICU.pred$pred[1]
ICU_self.se<-emp_aff.ICU.pred$se[1]

ICU_other<-emp_aff.ICU.pred$pred[2]
ICU_other.se<-emp_aff.ICU.pred$se[2]

ICU_self.dist<-rnorm(100000, mean=ICU_self, sd=ICU_self.se)
ICU_other.dist<-rnorm(100000, mean=ICU_other, sd=ICU_other.se)

emp_aff.ICU.pred.DF<-cbind(ICU_self.dist, ICU_other.dist)
colnames(emp_aff.ICU.pred.DF)<-c('CU Measure: Other', 'CU Measure: ICU')


tiff(paste0(graphics.folder, 'Aff_Moderated_by_ICU_Resp.tiff'), res=300, units='in', height=8, width=11)
mcmc_areas(emp_aff.ICU.pred.DF, prob = .95)+
  xlab(expression(~rho))+
  ggtitle('Association between Affective Empathy and CU Traits as a Function of CU Measure')
dev.off()

#Exploring moderation graphically involving dichotomous Out_resp variable (Cognitive Empathy Model)
Out_resp<-c(0,1)
emp_cog.Out.pred<-predict(fit.Emp_cog.mod_Out_resp, newmods = Out_resp, level = 95)

Out_self<-emp_cog.Out.pred$pred[1]
Out_self.se<-emp_cog.Out.pred$se[1]

Out_other<-emp_cog.Out.pred$pred[2]
Out_other.se<-emp_cog.Out.pred$se[2]

Out_self.dist<-rnorm(100000, mean=Out_self, sd=Out_self.se)
Out_other.dist<-rnorm(100000, mean=Out_other, sd=Out_other.se)

emp_cog.Out.pred.DF<-cbind(Out_self.dist, Out_other.dist)
colnames(emp_cog.Out.pred.DF)<-c('Respondent: Self', 'Respondent: Other')

tiff(paste0(graphics.folder, 'Cog_Moderated_by_Out_Resp.tiff'), res=300, units='in', height=8, width=11)
mcmc_areas(emp_cog.Out.pred.DF, prob = .95)+
  xlab(expression(~rho))+
  ggtitle('Association between Cognitive Empathy and CU Traits as a Function of Cognitive Empathy Respondent')
dev.off()

#############################################################################################################################
#############################################################################################################################
length(unique(dat$citation))#48 studies

library(dplyr)
#Total unique N included
temp.df.N <- dat %>% 
  distinct(citation, N) %>% 
  select(N)
sum(temp.df.N$N)

temp.df.age <- dat %>% 
  distinct(citation, age, N)

sum(temp.df.age$age*temp.df.age$N)/sum(temp.df.age$N)

temp.df.female <- dat %>% 
  distinct(citation, female, N)

sum(temp.df.female$female*temp.df.female$N)/sum(temp.df.female$N)/100

#############################################################################################################################
#creating a series of graphics that summarize changes in attentuation corrected vs. original estimates.
#############################################################################################################################

dat.Emp_tot$citation[dat.Emp_tot$citation=='Anastassiou-Hadjicharalambous & Warden (2008)']<-'Anastassiou-H. & Warden (2008)'

#Separating out rows for plotting mutliple effects 
#Note there are some studies that reported male and female separately
#We chose to treat these as effects related to different populations\

table(dat.Emp_tot$citation)
dat.Emp_tot$citation[1]<-'Dadds et al. (2009m)'   #male sample
dat.Emp_tot$citation[2]<-'Dadds et al. (2009f)'   #female sample

#Need separate effects on different lines (treated distinct )
g.Emp_tot<-ggplot()+
  geom_point(data=dat.Emp_tot, 
             aes(y=citation, 
                 x=R),
             fill = "#f3b810", 
             color = "#f3b810", 
             alpha=.5)+
  geom_errorbarh(data=dat.Emp_tot, 
                 aes(y=citation, 
                     xmin = R+qnorm(.025)*sqrt(R_var), 
                     xmax = R+qnorm(.975)*sqrt(R_var) 
                     ), 
                 color = "#f3b810", 
                 alpha=.5)+
  geom_point(data=dat.Emp_tot, 
             aes(y=citation, 
                 x=Eff),
             fill = "#8fb6ce", 
             color = "#8fb6ce", 
             alpha=.5)+
  geom_errorbarh(data=dat.Emp_tot, 
                 aes(y=citation, 
                     xmin = Eff+qnorm(.025)*sqrt(Eff_var), 
                     xmax = Eff+qnorm(.975)*sqrt(Eff_var) 
                 ), 
                 color = "#8fb6ce")+
  geom_vline(aes(xintercept =0), color='red', lty='dashed')+
  ylab('Citation')+
  xlab('CU - Total Empathy Effects - Corrected and Uncorrected')+
  ggtitle('Comparing Attenuation-Corrected (Blue) and Raw (Orange) Correlations: Total Empathy ~ CU')+
  scale_x_continuous(limits = c(-1.2, .6))
  
g.Emp_tot

######################################
table(dat.Emp_aff$citation)
dat.Emp_aff$citation[1]<-'Dadds et al. (2009m)'
dat.Emp_aff$citation[2]<-'Dadds et al. (2009f)'
dat.Emp_aff$citation[17]<-'Brouns et al. (2013m)'
dat.Emp_aff$citation[18]<-'Brouns et al. (2013f)'

g.Emp_aff<-ggplot()+
  geom_point(data=dat.Emp_aff, 
             aes(y=citation, 
                 x=R),
             fill = "#f3b810", 
             color = "#f3b810", 
             alpha=.5)+
  geom_errorbarh(data=dat.Emp_aff, 
                 aes(y=citation, 
                     xmin = R+qnorm(.025)*sqrt(R_var), 
                     xmax = R+qnorm(.975)*sqrt(R_var) 
                 ), 
                 color = "#f3b810", 
                 alpha=.5)+
  geom_point(data=dat.Emp_aff, 
             aes(y=citation, 
                 x=Eff),
             fill = "#8fb6ce", 
             color = "#8fb6ce", 
             alpha=.5)+
  geom_errorbarh(data=dat.Emp_aff, 
                 aes(y=citation, 
                     xmin = Eff+qnorm(.025)*sqrt(Eff_var), 
                     xmax = Eff+qnorm(.975)*sqrt(Eff_var) 
                 ), 
                 color = "#8fb6ce")+
  geom_vline(aes(xintercept =0), color='red', lty='dashed')+
  ylab('Citation')+
  xlab('CU - Affective Empathy Effects - Corrected and Uncorrected')+
  ggtitle('Comparing Attenuation-Corrected (Blue) and Raw (Orange) Correlations: Affective Empathy ~ CU')+
  scale_x_continuous(limits = c(-1.2, .6))

g.Emp_aff

######################################
table(dat.Emp_cog$citation)
dat.Emp_cog$citation[1]<-'Dadds et al. (2009m)'
dat.Emp_cog$citation[2]<-'Dadds et al. (2009f)'
dat.Emp_cog$citation[18]<-'Brouns et al. (2013m)'
dat.Emp_cog$citation[19]<-'Brouns et al. (2013f)'

g.Emp_cog<-ggplot()+
  geom_point(data=dat.Emp_cog, 
             aes(y=citation, 
                 x=R),
             fill = "#f3b810", 
             color = "#f3b810", 
             alpha=.5)+
  geom_errorbarh(data=dat.Emp_cog, 
                 aes(y=citation, 
                     xmin = R+qnorm(.025)*sqrt(R_var), 
                     xmax = R+qnorm(.975)*sqrt(R_var) 
                 ), 
                 color = "#f3b810", 
                 alpha=.5)+
  geom_point(data=dat.Emp_cog, 
             aes(y=citation, 
                 x=Eff),
             fill = "#8fb6ce", 
             color = "#8fb6ce", 
             alpha=.5)+
  geom_errorbarh(data=dat.Emp_cog, 
                 aes(y=citation, 
                     xmin = Eff+qnorm(.025)*sqrt(Eff_var), 
                     xmax = Eff+qnorm(.975)*sqrt(Eff_var) 
                 ), 
                 color = "#8fb6ce")+
  geom_vline(aes(xintercept =0), color='red', lty='dashed')+
  ylab('Citation')+
  xlab('CU - Cognitive Empathy Effects - Corrected and Uncorrected')+
  ggtitle('Comparing Attenuation-Corrected (Blue) and Raw (Orange) Correlations: Cognitive Empathy ~ CU')+
  scale_x_continuous(limits = c(-1.2, .6))

g.Emp_cog

######################################
table(dat.prosoc$citation)
dat.prosoc$citation[13]<-'Pechorro et al. (2013m)'
dat.prosoc$citation[14]<-'Pechorro et al. (2013f)'
dat.prosoc$citation[18]<-'Pechorro et al. (2015am)'
dat.prosoc$citation[19]<-'Pechorro et al. (2015af)'

g.prosoc<-ggplot()+
  geom_point(data=dat.prosoc, 
             aes(y=citation, 
                 x=R),
             fill = "#f3b810", 
             color = "#f3b810", 
             alpha=.5)+
  geom_errorbarh(data=dat.prosoc, 
                 aes(y=citation, 
                     xmin = R+qnorm(.025)*sqrt(R_var), 
                     xmax = R+qnorm(.975)*sqrt(R_var) 
                 ), 
                 color = "#f3b810", 
                 alpha=.5)+
  geom_point(data=dat.prosoc, 
             aes(y=citation, 
                 x=Eff),
             fill = "#8fb6ce", 
             color = "#8fb6ce", 
             alpha=.5)+
  geom_errorbarh(data=dat.prosoc, 
                 aes(y=citation, 
                     xmin = Eff+qnorm(.025)*sqrt(Eff_var), 
                     xmax = Eff+qnorm(.975)*sqrt(Eff_var) 
                 ), 
                 color = "#8fb6ce")+
  geom_vline(aes(xintercept =0), color='red', lty='dashed')+
  ylab('Citation')+
  xlab('CU - Prosociality Effects - Corrected and Uncorrected')+
  ggtitle('Comparing Attenuation-Corrected (Blue) and Raw (Orange) Correlations: Prosociality ~ CU')+
  scale_x_continuous(limits = c(-1.2, .6))

g.prosoc

######################################
table(dat.glt$citation)
g.glt<-ggplot()+
  geom_point(data=dat.glt, 
             aes(y=citation, 
                 x=R),
             fill = "#f3b810", 
             color = "#f3b810", 
             alpha=.5)+
  geom_errorbarh(data=dat.glt, 
                 aes(y=citation, 
                     xmin = R+qnorm(.025)*sqrt(R_var), 
                     xmax = R+qnorm(.975)*sqrt(R_var) 
                 ), 
                 color = "#f3b810", 
                 alpha=.5)+
  geom_point(data=dat.glt, 
             aes(y=citation, 
                 x=Eff),
             fill = "#8fb6ce", 
             color = "#8fb6ce", 
             alpha=.5)+
  geom_errorbarh(data=dat.glt, 
                 aes(y=citation, 
                     xmin = Eff+qnorm(.025)*sqrt(Eff_var), 
                     xmax = Eff+qnorm(.975)*sqrt(Eff_var) 
                 ), 
                 color = "#8fb6ce")+
  geom_vline(aes(xintercept =0), color='red', lty='dashed')+
  ylab('Citation')+
  xlab('CU - Guilt Effects - Corrected and Uncorrected')+
  ggtitle('Comparing Attenuation-Corrected (Blue) and Raw (Orange) Correlations: Guilt ~ CU')+
  scale_x_continuous(limits = c(-1.2, .6))

g.glt

#Note that some studies have multiple effects - have not partialed separted these (thought it would be possible)
#Will keep graphic for GitHub - but will not include in supplement without separating out within study effects 

tiff(paste0(graphics.folder,'Attenuation Correction vs Raw Effects.tiff'), res=300, units = 'in', height = 8, width=20)
cowplot::plot_grid(g.Emp_tot, g.Emp_aff, g.Emp_cog, g.prosoc, g.glt, 
                        labels = c('A', 'B', 'C', 'D', 'E'), nrow = 3)
dev.off()
#############################################################################################################################
library(metaSEM)
library(RColorBrewer)
#############################################################################################################################

#With Total Empathy Scores (this seems totally crazy to me...):
dat.meta3<-dat.imp[dat.imp$Emp_cog!=1,]
dat.meta3<-dat.meta3[dat.meta3$Emp_aff!=1,]

#Total of 49 effects left..
y<-dat.meta3$Eff
v<-dat.meta3$Eff_var
ID<-dat.meta3$id    #Note that this currently still has some within study break down by gender

#Prepping covariate matrix - for effects (at level 2)
x<-matrix(ncol=3, nrow = 49)
x[,1]<-dat.meta3$prosocial  #prosocial as the target category.
x[,2]<-dat.meta3$guilt    #guilt as the target category
x[,3]<-ifelse(dat.meta3$prosocial+dat.meta3$guilt>=1, 0, 1)

colnames(x)<-c('Prosocial', 
               'Guilt', 
               'Total Empathy')

#fitting null model
fit.meta_null<-meta3(y=y,
                     v=v,
                     cluster = ID
                     )

sink(paste0(model.folder, 'Three-level_Null_Model.txt'))
summary(fit.meta_null)  #provides a sense of overall effect across all three constructs
sink()

#Modeling Intercept as Total Empathy
fit.meta_Emp_tot<-meta3(y=y,
                        v=v,
                        cluster = ID,
                        x = x[,1:2]
)

sink(paste0(model.folder, 'Three-level_Emp_tot_Model.txt'))
summary(fit.meta_Emp_tot)
sink()

#Variance Decomposition
dat.meta3[,16:18]<-x
colnames(dat.meta3)[16:18]<-c('Prosocial', 
                              'Guilt', 
                              'Total Empathy')

Tau<-fit.meta_Emp_tot$mx.fit$Tau$values[2,2]
gamma_b<-fit.meta_Emp_tot$mx.fit$Inter$values[1,1]

gamma_w<-c(fit.meta_Emp_tot$mx.fit$Beta$values[1,1],
           fit.meta_Emp_tot$mx.fit$Beta$values[1,2])

within_cov<-c(16,17)
between_cov<-NULL

sigma2<-fit.meta_Emp_tot$mx.fit$Tau$values[1,1]

var.decomp<-r2MLM(data=dat.meta3, 
                  within_covs = within_cov,
                  between_covs = between_cov,
                  random_covs = NULL,
                  gamma_w = gamma_w, 
                  gamma_b = gamma_b, 
                  Tau = Tau,
                  sigma2 = sigma2, 
                  clustermeancentered = F)

sink(paste0(model.folder, 'Three-level_VarDecomp.txt'))
print(var.decomp)
sink()

#------------------------------------------------------------
#Getting all pairwise comparisons...
fit.meta_prosoc<-meta3(y=y,
                    v=v,
                    cluster = ID,
                    x = x[,c(2,3)]
)
summary(fit.meta_prosoc)

#Quick summary - in the joint model: 
#     1. All three significantly predict CU scores
#     2. Prosociality is more strongly correlated with CU scores that either other covariate
#     3. There is no difference between guilt and total empathy in their relations with CU 

############################################################################################################################
#Uncorrected Tables for Supplemental
dat.graphS1<-dat.imp
dat.graphS1<-dat.graphS1[dat.graphS1$Emp_aff==1 | dat.graphS1$Emp_cog==1,]

dat.graphS1$cite.fac<-as.factor(dat.graphS1$citation)
fit.graph1<-rma(yi=R, vi=R_var, data=dat.graphS1)

fit.emp_cogR<-rma(yi=R, 
                  vi=R_var,
                  weights = 1/R_var,
                  data=dat.Emp_cog, 
                  ni=N)

summary(fit.emp_cogR)

fit.emp_affR<-rma(yi=R, 
                  vi=R_var,
                  weights = 1/R_var,
                  data=dat.Emp_aff, 
                  ni=N)

summary(fit.emp_affR)

#Adding back in Outcome as a variable
dat.graphS1$Outcome[dat.graphS1$Emp_aff==1]<-'Affective Empathy'
dat.graphS1$Outcome[dat.graphS1$Emp_cog==1]<-'Cognitive Empathy'
dat.graphS1$female<-round(dat.graphS1$female, digits = 2)
dat.graphS1$age<-round(dat.graphS1$age, digits = 2)

par(mar=c(4,4,1,2))

tiff(paste0(graphics.folder,'Affective vs Cognitive Empathy_Uncorrected.tiff'), res=300, units = 'in', height = 8.5, width=11)
par(cex=.80, font=1)
forest(fit.graph1, xlim=c(-10, 2), 
       order = order(dat.graphS1$Outcome, dat.graphS1$Eff), 
       ilab = cbind(dat.graphS1$N, 
                    dat.graphS1$female, 
                    dat.graphS1$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:20, 25:43),
       ylim = c(-1,47), 
       cex=.75, 
       xlab=expression(~rho), 
       mlab="", 
       addfit = F,
       slab = dat.graphS1$cite.fac)

par(cex=.80, font=4)
text(-10, c(21, 44), pos=4, c('Affective Empathy', 'Cogntive Empathy'))

par(font=4)
text(-7, 46, 'N')
par(font=2)
text(c(-5, -3), 46, c('%Female', 'Mean Age'))
text(-10, 46, 'Citation', pos=4)
par(font=4)
text(1.5,46, "r")

addpoly(fit.emp_affR, row=1.5, cex=.7, mlab = "")
addpoly(fit.emp_cogR, row=23.5, cex=.7, mlab = "")

text(-9, 1.5, 
     pos=4, 
     cex=0.75,
     bquote(
       paste("RE Model for Affective Empathy (Q = ", 
             .(formatC(fit.emp_affR$QE, digits=2, format="f")), 
             ", df = ", 
             .(fit.emp_affR$k - fit.emp_affR$p),
             ", p = ", 
             .(formatC(fit.emp_affR$QEp, digits=2, format="f")), 
             "; ", I^2, " = ",
             .(formatC(fit.emp_affR$I2, digits=1, format="f")), "%)")))

text(-9, 23.5, 
     pos=4, 
     cex=0.75,
     bquote(
       paste("RE Model for Cognitive Empathy (Q = ", 
             .(formatC(fit.emp_cogR$QE, digits=2, format="f")), 
             ", df = ", 
             .(fit.emp_cogR$k - fit.emp_cogR$p),
             ", p = ", 
             .(formatC(fit.emp_cogR$QEp, digits=2, format="f")), 
             "; ", I^2, " = ",
             .(formatC(fit.emp_cogR$I2, digits=1, format="f")), "%)")))
title('Uncorrected Effects for Cognitive and Affective Empathy')
dev.off()

#Unconditional Models (no moderators - publication graphics)
#Fitting two summary graphics (one for affective vs. cognitive) & one for summary vals
dat.graphS2<-dat.imp

dat.graphS2<-dat.graphS2[dat.graphS2$Emp_aff!=1,]
dat.graphS2<-dat.graphS2[dat.graphS2$Emp_cog!=1,]
dat.graphS2$citation[dat.graphS2$citation=='Anastassiou-Hadjicharalambous & Warden (2008)']<-'Anastassiou-H. & Warden (2008)'

dat.graphS2$cite.fac<-as.factor(dat.graphS2$citation)
fit.graph2<-rma(yi=R, vi=R_var, data=dat.graphS2)

fit.emp_totR<-rma(yi=R, 
                  vi=R_var, 
                  weights = 1/R_var,
                  data=dat.Emp_tot, 
                  ni=N)

summary(fit.emp_totR)

fit.prosocR<-rma(yi=R, 
                 vi=R_var, 
                 weights = 1/R_var,
                 data=dat.prosoc, 
                 ni=N)

summary(fit.prosocR)

fit.gltR<-rma(yi=R, 
              vi=R_var, 
              weights = 1/R_var,
              data=dat.glt, 
              ni=N, 
              method = 'FE')

summary(fit.gltR)

#round down age and % female for purposes of plotting
dat.graphS2$female<-round(dat.graphS2$female, digits = 2)
dat.graphS2$age<-round(dat.graphS2$age, digits = 2)
dat.graphS2$Outcome[dat.graphS2$prosocial==1]<-'prosoc'
dat.graphS2$Outcome[dat.graphS2$guilt==1]<-'guilt'
dat.graphS2$Outcome<-ifelse(!is.na(dat.graphS2$Outcome), dat.graphS2$Outcome, 'emp_tot')

#Note to self - adding spaces into variables - even string variables makes ordering of effects difficult on 
#accompanying table/graphic

par(mar=c(4,4,1,2))

tiff(paste0(graphics.folder,'Main Effects summary_Uncorrected.tiff'), res=300, units = 'in', height = 8.5, width=11)
par(cex=.80, font=1)
forest(fit.graph2, xlim=c(-10, 2), 
       order = order(dat.graphS2$Outcome, dat.graphS2$R), 
       ilab = cbind(dat.graphS2$N, 
                    dat.graphS2$female, 
                    dat.graphS2$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:29, 34:36, 41:59),
       ylim = c(-1,63), 
       cex=.75, 
       xlab = expression(~rho), 
       mlab="", 
       addfit = F,
       slab = dat.graphS2$cite.fac)

par(cex=.80, font=4)
text(-10, c(30, 37, 60), pos=4, c('Total Empathy', 'Guilt', 'Prosociality'))

par(font=4)
text(-7, 62, 'N')
par(font=2)
text(c(-5, -3), 62, c('%Female', 'Mean Age'))
text(-10, 62, 'Citation', pos=4)
par(font=4)
text(1.5,62, "r")

addpoly(fit.emp_totR, row=1.5, cex=.7, mlab = "")
addpoly(fit.gltR, row=32.5, cex=.7, mlab = "")
addpoly(fit.prosocR, row=39.5, cex=.7, mlab = "")

text(-9, 1.5, 
     pos=4, 
     cex=0.75,
     bquote(
       paste("RE Model for Total Empathy (Q = ", 
             .(formatC(fit.emp_totR$QE, digits=2, format="f")), 
             ", df = ", 
             .(fit.emp_totR$k - fit.emp_totR$p),
             ", p = ", 
             .(formatC(fit.emp_totR$QEp, digits=2, format="f")), 
             "; ", I^2, " = ",
             .(formatC(fit.emp_totR$I2, digits=1, format="f")), "%)")))

text(-9, 32.5, 
     pos=4, 
     cex=0.75,
     bquote(
       paste("FE Model for Guilt (Q = ", 
             .(formatC(fit.gltR$QE, digits=2, format="f")), 
             ", df = ", 
             .(fit.gltR$k - fit.gltR$p),
             ", p = ", 
             .(formatC(fit.gltR$QEp, digits=2, format="f")), 
             "; ", I^2, " = ",
             .(formatC(fit.gltR$I2, digits=1, format="f")), "%)")))

text(-9, 39.5, 
     pos=4, 
     cex=0.75,
     bquote(
       paste("RE Model for Prosociality (Q = ", 
             .(formatC(fit.prosocR$QE, digits=2, format="f")), 
             ", df = ", 
             .(fit.prosocR$k - fit.prosocR$p),
             ", p = ", 
             .(formatC(fit.prosocR$QEp, digits=2, format="f")), 
             "; ", I^2, " = ",
             .(formatC(fit.prosocR$I2, digits=1, format="f")), "%)")))
title('Uncorrected Effects for Prosociality, Guilty, and Total Empathy')
dev.off()

#Creating separate tables for each construct from main effects table: 
dat.graph2a<-dat.graph2[dat.graph2$Outcome=='emp_tot',]
fit.graph2a<-rma(yi=Eff, vi=Eff_var, data=dat.graph2a)

tiff(paste0(graphics.folder,'Total Empathy Summary.tiff'), res=600, units = 'in', height = 7.5, width=11)
par(cex=.80, font=1)
forest(fit.graph2a, xlim=c(-10, 2), 
       order = order(dat.graph2a$Eff), 
       ilab = cbind(dat.graph2a$N, 
                    dat.graph2a$female, 
                    dat.graph2a$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:29),
       ylim = c(0,32), 
       cex=.75, 
       xlab = expression(~rho), 
       mlab="", 
       addfit = F,
       slab = dat.graph2a$cite.fac)

par(font=4)
text(-7, 31, 'N')
par(font=2)
text(c(-5, -3), 31, c('%Female', 'Mean Age'))
text(-10, 31, 'Citation', pos=4)
par(font=4)
text(1.5, 31, 'r')

addpoly(fit.emp_tot, row=1.5, cex=.7, mlab = "")

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
title('Corrected Effects for Total Empathy')
dev.off()

##
dat.graph2b<-dat.graph2[dat.graph2$Outcome=='prosoc',]
fit.graph2b<-rma(yi=Eff, vi=Eff_var, data=dat.graph2b)

tiff(paste0(graphics.folder,'Prosociality Summary.tiff'), res=600, units = 'in', height = 6.5, width=11)
par(cex=.80, font=1)
forest(fit.graph2b, xlim=c(-10, 2), 
       order = order(dat.graph2b$Eff), 
       ilab = cbind(dat.graph2b$N, 
                    dat.graph2b$female, 
                    dat.graph2b$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:21),
       ylim = c(0,24), 
       cex=.75, 
       xlab = expression(~rho), 
       mlab="", 
       addfit = F,
       slab = dat.graph2b$cite.fac)

par(font=4)
text(-7, 23, 'N')
par(font=2)
text(c(-5, -3), 23, c('%Female', 'Mean Age'))
text(-10, 23, 'Citation', pos=4)
par(font=4)
text(1.5, 23, 'r')

addpoly(fit.prosoc, row=1.5, cex=.7, mlab = "")

text(-9, 1.5, 
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
title('Corrected Effects for Prosociality')
dev.off()

##
dat.graph2c<-dat.graph2[dat.graph2$Outcome=='guilt',]
fit.graph2c<-rma(yi=Eff, vi=Eff_var, data=dat.graph2c)

tiff(paste0(graphics.folder,'Guilt Summary.tiff'), res=600, units = 'in', height = 3, width=11)
par(cex=.80, font=1)
forest(fit.graph2c, xlim=c(-10, 2), 
       order = order(dat.graph2c$Eff), 
       ilab = cbind(dat.graph2c$N, 
                    dat.graph2c$female, 
                    dat.graph2c$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:5),
       ylim = c(0,9), 
       cex=.75, 
       xlab = expression(~rho), 
       mlab="", 
       addfit = F,
       slab = dat.graph2c$cite.fac)

par(font=4)
text(-7, 8, 'N')
par(font=2)
text(c(-5, -3), 8, c('%Female', 'Mean Age'))
text(-10, 8, 'Citation', pos=4)
par(font=4)
text(1.5, 8, 'r')

addpoly(fit.glt, row=1.5, cex=.7, mlab = "")

text(-9, 1.5, 
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
title('Corrected Effects for Guilt')
dev.off()

#Empathy Comparison Model
dat.emp_comp_uc<-merge(dat.Emp_aff, dat.Emp_cog, by='id')

cols<-c(1,16,2:4,9, 10, 25,26,18,19,13,17,34,35,50,51)
dat.emp_comp_uc<-dat.emp_comp_uc[,cols]
colnames(dat.emp_comp_uc)<-c('id', 'citation', 'female', 'age', 'N', 'R_affective', 'R_var_affective',
                          'Eff_affective', 'Eff_var_affective', 'CU_resp', 'Out_resp',
                          'Sample','ICU', 'R_cognitive', 'R_var_cognitive', 'Eff_cognitive', 'Eff_var_cognitive')

#Correlations pulled from studies when available, ordered by citation: 
Aff_cog_cor<-c(.063, .063, .35, 0, 0, .49, 0, .5, 0, .76, .32, 0, -.01, .08,0,.45,0,0)
#Calculating difference in effects and variance terms adjusted for outcome correlation
Eff<-dat.emp_comp_uc$R_affective-dat.emp_comp_uc$R_cognitive
Eff_var<-dat.emp_comp_uc$R_var_affective+dat.emp_comp_uc$R_var_cognitive-2*Aff_cog_cor*sqrt(dat.emp_comp_uc$R_var_affective*dat.emp_comp_uc$R_var_cognitive)

#note that I stupidly changed naming conventions here - sorry
#the Eff variable is corrected for attentuation, censored when an effect was < -1, and adjusted for the dependent nature of the correlations where possible
dat.emp_comp_uc$Eff<-Eff
dat.emp_comp_uc$Eff_var<-Eff_var

fit.emp_comp_uc<-rma(yi=Eff, 
                  vi=Eff_var, 
                  weights = 1/Eff_var, 
                  data=dat.emp_comp_uc, 
                  ni=N
)

summary(fit.emp_comp_uc)
