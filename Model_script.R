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

dat<-read.csv(paste0(wd, 'Model_Data_w-Reliability_081218.csv'), stringsAsFactors = F)
colnames(dat)[1]<-'id'

#Imputing missing reliabilities prior to recoding. 
psych::describe(dat)

#Need to create a series of binary variables - cannot handle categorical variables
dat$Emp_tot<-ifelse(dat$Outcome=='empathy_tot', 1, 0)
dat$Emp_aff<-ifelse(dat$Outcome=='empathy_aff', 1, 0)
dat$Emp_cog<-ifelse(dat$Outcome=='empathy_cog', 1, 0)
dat$prosocial<-ifelse(dat$Outcome=='prosocial', 1, 0)
dat$guilt<-ifelse(dat$Outcome=='guilt', 1, 0)

dat$R_var<-((1-dat$R^2)^2)/(dat$N-1)

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

#Imputation Model
imp<-mice(dat.imp, maxit=100, m=5)

png(paste0(graphics.folder, 'MICE_Imputation_convergence.png'), 
    res=300, 
    units='in', 
    height = 4, 
    width = 6)
plot(imp)
dev.off()

#Extracting and obtaining estimates for missing values (averaging across the 20 imputed datasets)
dat.imp<-complete(imp)
dat.imp$citation<-dat$citation

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

#Creating Attenuation-corrected values
dat.imp$Eff<-dat.imp$R/sqrt(dat.imp$CU_Rel*dat.imp$Out_Rel)

#There are some values that exceed -1 threshold. 
#For now randomly selecting betwen -.9 and -.999 
#Better alternatives may exist
#3/4 of effects with issues involved small samples, smallish reliabilities, and d's instead of r's
#(note the d's were converted to r's)

dat.imp$Eff.r<-ifelse(dat.imp$Eff<=-1, 
                      runif(n=1, max = -.9, min=-.99),
                      dat.imp$Eff)
dat.imp$Eff_var<-dat.imp$R_var/sqrt(dat.imp$CU_Rel*dat.imp$Out_Rel)

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
fit.emp_tot<-rma(yi=Eff.r, 
                 vi=Eff_var, 
                 weights = 1/Eff_var, 
                 data=dat.Emp_tot, 
                 ni=N,
                 method = 'HS')

summary(fit.emp_tot)

sink(paste0(model.folder, 'Total_Empathy_AC.txt'))
summary(fit.emp_tot)
sink()

jpeg(paste0(graphics.folder, 'CU and Total Empathy - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_tot, order = 'obs', slab=dat.Emp_tot$citation)
title("Relation between CU Traits and Total Empathy")
dev.off()

jpeg(paste0(graphics.folder, 'CU and Total Empathy - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.emp_tot)
title("Funnel Plot of CU Relation with Total Empathy")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
REG.emp_tot<-regtest(fit.emp_tot, model = 'lm')
REG.emp_tot

#Trim and fill model
fit.emp_tot.TF_R<-trimfill(fit.emp_tot, estimator = 'R0')
fit.emp_tot.TF_R

sink(paste0(model.folder, 'Total_Empathy_R0.txt'))
summary(fit.emp_tot.TF_R)
sink()

jpeg(paste0(graphics.folder, 'R0_estimator_CU-tot_emp_Sup1A.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.emp_tot.TF_R)
title(expression('Trim and Fill Results for Total Empathy (R'[0]*' Estimator)'))
dev.off()

#If correct, interpretations are completely different now
#Bias against publishing negative effects (may be that it is not novel enough?)
#Should report the corrected values using trim-and-fill - graph both, will stick with R0 in MS 
fit.emp_tot.TF_L<-trimfill(fit.emp_tot, estimator = 'L0')
fit.emp_tot.TF_L

sink(paste0(model.folder, 'Total_Empathy_L0.txt'))
summary(fit.emp_tot.TF_L)
sink()

jpeg(paste0(graphics.folder, 'L0_estimator_CU-tot_emp_Sup1B.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.emp_tot.TF_L)
title(expression('Trim and Fill Results for Total Empathy (L'[0]*' Estimator)'))
dev.off()

#################################################################################
#Model for Affective Empathy: 
fit.emp_aff<-rma(yi=Eff.r, 
                 vi=Eff_var, 
                 weights = 1/Eff_var, 
                 data=dat.Emp_aff, 
                 ni=N,
                 method = 'HS')

summary(fit.emp_aff)

sink(paste0(model.folder, 'Affective_Empathy_AC.txt'))
summary(fit.emp_aff)
sink()

jpeg(paste0(graphics.folder, 'CU and Affective Empathy - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_aff, order = 'obs', slab=dat.Emp_aff$citation)
title("Relation between CU Traits and Affective Empathy")
dev.off()

jpeg(paste0(graphics.folder, 'CU and Affective Empathy - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.emp_aff)
title("Funnel Plot of CU Relation with Affective Empathy")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
REG.emp_aff<-regtest(fit.emp_aff, model = 'lm')
REG.emp_aff

#Trim and fill model
fit.emp_aff.TF_R<-trimfill(fit.emp_aff, estimator = 'R0')
fit.emp_aff.TF_R

sink(paste0(model.folder, 'Affective_Empathy_R0.txt'))
summary(fit.emp_aff.TF_R)
sink()

jpeg(paste0(graphics.folder, 'R0_estimator_CU-aff_emp_Sup2A.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.emp_aff.TF_R)
title(expression('Trim and Fill Results for Affective Empathy (R'[0]*' Estimator)'))
dev.off()

#If correct, interpretations are completely different now
#Bias against publishing negative effects (may be that it is not novel enough?)
#Should report the corrected values using trim-and-fill - graph both, will stick with R0 in MS 
fit.emp_aff.TF_L<-trimfill(fit.emp_aff, estimator = 'L0')
fit.emp_aff.TF_L

sink(paste0(model.folder, 'Affective_Empathy_L0.txt'))
summary(fit.emp_aff.TF_L)
sink()

jpeg(paste0(graphics.folder, 'L0_estimator_CU-aff_emp_Sup2B.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.emp_aff.TF_L)
title(expression('Trim and Fill Results for Affective Empathy (L'[0]*' Estimator)'))
dev.off()

#################################################################################
#Model for Cognitive Empathy: 
fit.emp_cog<-rma(yi=Eff.r, 
                 vi=Eff_var, 
                 weights = 1/Eff_var, 
                 data=dat.Emp_cog, 
                 ni=N,
                 method = 'HS')

summary(fit.emp_cog)

sink(paste0(model.folder, 'Cognitive_Empathy_AC.txt'))
summary(fit.emp_cog)
sink()

jpeg(paste0(graphics.folder, 'CU and Cognitive Empathy - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_cog, order = 'obs', slab=dat.Emp_cog$citation)
title("Relation between CU Traits and Cognitive Empathy")
dev.off()

jpeg(paste0(graphics.folder, 'CU and Cognitive Empathy - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.emp_cog)
title("Funnel Plot of CU Relation with Cognitive Empathy")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
REG.emp_cog<-regtest(fit.emp_cog, model = 'lm')
REG.emp_cog

#Trim and fill model
fit.emp_cog.TF_R<-trimfill(fit.emp_cog, estimator = 'R0')
fit.emp_cog.TF_R

sink(paste0(model.folder, 'Cognitive_Empathy_R0.txt'))
summary(fit.emp_cog.TF_R)
sink()

jpeg(paste0(graphics.folder, 'R0_estimator_CU-cog_emp_Sup3A.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.emp_cog.TF_R)
title(expression('Trim and Fill Results for Cognitive Empathy (R'[0]*' Estimator)'))
dev.off()

#If correct, interpretations are completely different now
#Bias against publishing negative effects (may be that it is not novel enough?)
#Should report the corrected values using trim-and-fill - graph both, will stick with R0 in MS 
fit.emp_cog.TF_L<-trimfill(fit.emp_cog, estimator = 'L0')
fit.emp_cog.TF_L

sink(paste0(model.folder, 'Cognitive_Empathy_L0.txt'))
summary(fit.emp_cog.TF_L)
sink()

jpeg(paste0(graphics.folder, 'L0_estimator_CU-cog_emp_Sup3B.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.emp_cog.TF_L)
title(expression('Trim and Fill Results for Cognitive Empathy (L'[0]*' Estimator)'))
dev.off()

#################################################################################
#Model for Guilt: 
fit.glt<-rma(yi=Eff.r, 
             vi=Eff_var, 
             weights = 1/Eff_var, 
             data=dat.glt, 
             ni=N,
             method = 'FE')
summary(fit.glt)

sink(paste0(model.folder, 'Guilt Overall Model - no moderators.txt'))
summary(fit.glt)
sink()

jpeg(paste0(graphics.folder, 'CU and Guilt - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.glt, order = 'obs', slab=dat.glt$citation)
title("Relation between CU Traits and Guilt")
dev.off()

jpeg(paste0(graphics.folder, 'CU and Guilt - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.glt)
title("Funnel Plot of CU Relation with Guilt")
dev.off()

###########################################################################################
#Model comparing two forms of empathy - within study differences in effects
#merging and cleaning data set
dat.emp_comp<-merge(dat.Emp_aff, dat.Emp_cog, by='id')

cols<-c(1,16,3:5,10,23,24,17,18,13,32,46,47)
dat.emp_comp<-dat.emp_comp[,cols]
colnames(dat.emp_comp)<-c('id', 'citation', 'female', 'age', 'N', 'R_affective',
                          'Eff_affective', 'Eff_var_affective', 'CU_resp', 'Out_resp',
                          'Sample', 'R_cognitive', 'Eff_cognitive', 'Eff_var_cognitive')

#Correlations pulled from studies, ordered by citation: 
Aff_cog_cor<-c(.063, .063, .35, 0, 0, .49, 0, .5, 0, .76, .32, 0, -.01, .08,0,.45,0,0)

#Calculating difference in effects and variance terms adjusted for outcome correlation
Eff<-dat.emp_comp$Eff_affective-dat.emp_comp$Eff_cognitive
Eff_var<-dat.emp_comp$Eff_var_affective+dat.emp_comp$Eff_var_cognitive-2*Aff_cog_cor*sqrt(dat.emp_comp$Eff_var_affective*dat.emp_comp$Eff_var_cognitive)

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

jpeg(paste0(graphics.folder, 'Difference in Empathy Model - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_comp, order = 'obs', slab=dat.emp_comp$citation)
title("Difference in Relation between CU and Empathy Dimensions")
dev.off()

jpeg(paste0(graphics.folder, 'Difference in Empathy Dimension Model - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.emp_comp)
title("Funnel Plot of Difference in Relation between CU and Empathy Dimensions")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
#Trim and fill model
fit.emp_comp.TF_R<-trimfill(fit.emp_comp, estimator = 'R0')
fit.emp_comp.TF_R

sink(paste0(model.folder, 'Cognitive_Empathy_R0.txt'))
summary(fit.emp_comp.TF_R)
sink()

jpeg(paste0(graphics.folder, 'R0_estimator_CU-comp_emp_Sup4A.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.emp_comp.TF_R)
title(expression('Trim and Fill Results for Cognitive Empathy (R'[0]*' Estimator)'))
dev.off()

#If correct, interpretations are completely different now
#Bias against publishing negative effects (may be that it is not novel enough?)

fit.emp_comp.TF_L<-trimfill(fit.emp_comp, estimator = 'L0')
fit.emp_comp.TF_L

sink(paste0(model.folder, 'Cognitive_Empathy_L0.txt'))
summary(fit.emp_comp.TF_L)
sink()

jpeg(paste0(graphics.folder, 'L0_estimator_CU-comp_emp_Sup4B.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.emp_comp.TF_L)
title(expression('Trim and Fill Results for Cognitive Empathy (L'[0]*' Estimator)'))
dev.off()

#------------------------------------------------------------------------------------------
#Model for Cognitive Empathy: 
fit.prosoc<-rma(yi=Eff.r, 
                vi=Eff_var, 
                weights = 1/Eff_var, 
                data=dat.prosoc, 
                ni=N,
                method = 'HS')

summary(fit.prosoc)

sink(paste0(model.folder, 'Prosocial_AC.txt'))
summary(fit.prosoc)
sink()

jpeg(paste0(graphics.folder, 'CU and Prosocialty - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.prosoc, order = 'obs', slab=dat.prosoc$citation)
title("Relation between CU Traits and Prosocialty")
dev.off()

jpeg(paste0(graphics.folder, 'CU and Prosocialty - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.prosoc)
title("Funnel Plot of CU Relation with Prosocialty")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
REG.prosoc<-regtest(fit.prosoc, model = 'lm')
REG.prosoc

#Trim and fill model
fit.prosoc.TF_R<-trimfill(fit.prosoc, estimator = 'R0')
fit.prosoc.TF_R

sink(paste0(model.folder, 'Prosocial_R0.txt'))
summary(fit.prosoc.TF_R)
sink()

jpeg(paste0(graphics.folder, 'R0_estimator_CU-prosoc_Sup5A.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.prosoc.TF_R)
title(expression('Trim and Fill Results for Prosocialty (R'[0]*' Estimator)'))
dev.off()

#If correct, interpretations are completely different now
#Bias against publishing negative effects (may be that it is not novel enough?)
#Should report the corrected values using trim-and-fill - graph both, will stick with R0 in MS 
fit.prosoc.TF_L<-trimfill(fit.prosoc, estimator = 'L0')
fit.prosoc.TF_L

sink(paste0(model.folder, 'Prosocial_L0.txt'))
summary(fit.prosoc.TF_L)
sink()

jpeg(paste0(graphics.folder, 'L0_estimator_CU-prosoc_Sup5B.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.prosoc.TF_L)
title(expression('Trim and Fill Results for Prosocialty (L'[0]*' Estimator)'))
dev.off()

#############################################################################
#Function below not necessary if effects not converted r-to-z transfomrmation. 

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
                       k=c(27, 18, 19, 19, 3, 18),
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

jpeg(paste0(graphics.folder, 'Unconditional Model Summary Graphic.jpeg'), res=300, units='in', height=7, width=7)
mcmc_areas(DF.plot.dist2, prob=.95)+
  xlab(expression('Population Esimates of'~rho))
dev.off()

###########################################################################################
#Unconditional Models (no moderators - publication graphics)
#Fitting two summary graphics (one for affective vs. cognitive) & one for summary vals
#Code for graphic models was adapted from: http://www.metafor-project.org/doku.php/plots:forest_plot_with_subgroups
dat.graph1<-dat.AC
dat.graph1<-data.frame(citation = dat$citation,
                       Outcome = dat$Outcome,
                       dat.graph1)

dat.graph1<-dat.graph1[dat.graph1$Emp_aff==1 | dat.graph1$Emp_cog==1,]

dat.graph1$cite.fac<-as.factor(dat.graph1$citation)
fit.graph1<-rma(yi=Eff, vi=Eff_var, data=dat.graph1)
par(mar=c(4,4,1,2))

jpeg(paste0(graphics.folder,'Affective vs. Cognitive Empathy.jpeg'), res=300, units = 'in', height = 8.5, width=11)
par(cex=.80, font=1)
forest(fit.graph1, xlim=c(-10, 2), 
       order = order(dat.graph1$Outcome, dat.graph1$Eff), 
       ilab = cbind(dat.graph1$N, 
                    dat.graph1$female, 
                    dat.graph1$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:20, 25:43),
       ylim = c(-1,47), 
       cex=.75, 
       xlab=expression(~rho), mlab="", 
       addfit = F,
       slab = dat.graph1$cite.fac)

par(cex=.80, font=4)
text(-10, c(21, 44), pos=4, c('Affective Empathy', 'Cogntive Empathy'))

par(font=4)
text(-7, 46, 'N')
par(font=2)
text(c(-5, -3), 46, c('%Female', 'Mean Age'))
text(-10, 46, 'Citation', pos=4)
par(font=4)
text(0,46, "r")

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
dat.graph2<-dat.AC
dat.graph2<-data.frame(citation = dat$citation,
                       Outcome = dat$Outcome,
                       dat.graph2)

dat.graph2<-dat[dat$Outcome!='empathy_aff',]
dat.graph2<-dat.graph2[dat.graph2$Outcome!='empathy_cog',]
dat.graph2$citation[dat.graph2$citation=='Anastassiou-Hadjicharalambous & Warden (2008)']<-'Anastassiou-H. & Warden (2008)'

dat.graph2$cite.fac<-as.factor(dat.graph2$citation)
fit.graph2<-rma(yi=Eff, vi=Eff_var, data=dat.graph2)

#round down age and % female for purposes of plotting
dat.graph2$female<-round(dat.graph2$female, digits = 2)
dat.graph2$age<-round(dat.graph2$age, digits = 2)
par(mar=c(4,4,1,2))

jpeg(paste0(graphics.folder,'Main Effects summary.jpeg'), res=300, units = 'in', height = 8.5, width=11)
par(cex=.80, font=1)
forest(fit.graph2, xlim=c(-10, 2), 
       order = order(dat.graph2$Outcome, dat.graph2$Eff), 
       ilab = cbind(dat.graph2$N, 
                    dat.graph2$female, 
                    dat.graph2$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:29, 34:36, 41:59),
       ylim = c(-1,63), 
       cex=.75, 
       xlab="Fisher's z", mlab="", 
       addfit = F,
       slab = dat.graph2$cite.fac)

par(cex=.80, font=4)
text(-10, c(30, 37, 60), pos=4, c('Total Empathy', 'Guilt', 'Prosociality'))

par(font=4)
text(-7, 62, 'N')
par(font=2)
text(c(-5, -3), 62, c('%Female', 'Mean Age'))
text(-10, 62, 'Citation', pos=4)
par(font=4)
text(0,62, "r")

addpoly(fit.emp_tot, row=1.5, cex=.7, mlab = "")
addpoly(fit.glt, row=32.5, cex=.7, mlab = "")
addpoly(fit.prosoc, row=39.5, cex=.7, mlab = "")

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

text(-9, 32.5, 
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

text(-9, 39.5, 
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
#######################################################################################
#Total Empathy
fit.emp_tot.mod_all<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var, 
                         mods = ~female+age+Samp_typ+CU_resp+Out_resp, 
                         data=dat.Emp_tot, 
                         ni=N, 
                         knha=T)
summary(fit.emp_tot.mod_all)

#Prosocial Model:
table(dat.prosoc$CU_resp, dat.prosoc$Out_resp)
#Perfectly equal - just pick one
fit.prosoc.mod_all<-rma(yi=Eff, 
                        vi=Eff_var,
                        weights = 1/Eff_var,
                        mods = ~female+age+Samp_typ+CU_resp, 
                        data=dat.prosoc, 
                        ni=N, 
                        knha=T)
summary(fit.prosoc.mod_all)  

#Affective Empathy:
table(dat.Emp_aff$CU_resp, dat.Emp_aff$Out_resp)

fit.Emp_aff.mod_all<-rma(yi=Eff, 
                        vi=Eff_var,
                        weights = 1/Eff_var,
                        mods = ~female+age+Samp_typ+CU_resp+Out_resp, 
                        data=dat.Emp_aff, 
                        ni=N, 
                        knha=T)
summary(fit.Emp_aff.mod_all)  

#Cognitive Empathy Model:
table(dat.Emp_cog$CU_resp, dat.Emp_cog$Out_resp)
#Perfectly equal - just pick one
fit.Emp_cog.mod_all<-rma(yi=Eff, 
                        vi=Eff_var,
                        weights = 1/Eff_var,
                        mods = ~female+age+Samp_typ+CU_resp+Out_resp, 
                        data=dat.Emp_cog, 
                        ni=N, 
                        knha=T)
summary(fit.Emp_cog.mod_all)  

##############################################################################################
#Exploring single moderators now - Female
fit.emp_tot.mod_female<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var, 
                         mods = ~female, 
                         data=dat.Emp_tot, 
                         ni=N, 
                         knha=T)
summary(fit.emp_tot.mod_female)

#Prosocial Model:
fit.prosoc.mod_female<-rma(yi=Eff, 
                        vi=Eff_var,
                        weights = 1/Eff_var,
                        mods = ~female, 
                        data=dat.prosoc, 
                        ni=N, 
                        knha=T)
summary(fit.prosoc.mod_female)  

#Affective Empathy:
fit.Emp_aff.mod_female<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~female, 
                         data=dat.Emp_aff, 
                         ni=N, 
                         knha=T)
summary(fit.Emp_aff.mod_female)  

#Cognitive Empathy Model:
fit.Emp_cog.mod_female<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~female, 
                         data=dat.Emp_cog, 
                         ni=N, 
                         knha=T)
summary(fit.Emp_cog.mod_female)  

#############################################################################################
#Exploring Moderators - age only
fit.emp_tot.mod_age<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var, 
                         mods = ~age, 
                         data=dat.Emp_tot, 
                         ni=N, 
                         knha=T)
summary(fit.emp_tot.mod_age)

#Prosocial Model:
fit.prosoc.mod_age<-rma(yi=Eff, 
                        vi=Eff_var,
                        weights = 1/Eff_var,
                        mods = ~age, 
                        data=dat.prosoc, 
                        ni=N, 
                        knha=T)
summary(fit.prosoc.mod_age)  

fit.prosoc.mod_age2<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~age+I(age^2), 
                         data=dat.prosoc, 
                         ni=N, 
                         knha=T)
summary(fit.prosoc.mod_age2)

#Affective Empathy:
fit.Emp_aff.mod_age<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~age, 
                         data=dat.Emp_aff, 
                         ni=N, 
                         knha=T)
summary(fit.Emp_aff.mod_age)  

#Cognitive Empathy Model:
fit.Emp_cog.mod_age<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~age, 
                         data=dat.Emp_cog, 
                         ni=N, 
                         knha=T)
summary(fit.Emp_cog.mod_age)  

#Exploring moderation graphically. 
age<-seq(3, 18, by=.05)
prosoc.age.pred<-predict(fit.prosoc.mod_age, newmods = age, level = 95)
prosoc.age.pred.DF<-data.frame(Age=age, Pred_vals=prosoc.age.pred$pred, CI_UB=prosoc.age.pred$ci.ub, CI_LB=prosoc.age.pred$ci.lb)

g1<-ggplot()+
  geom_line(data=prosoc.age.pred.DF, aes(x=Age, y=Pred_vals),color="#03396c")+
  geom_ribbon(data=prosoc.age.pred.DF, aes(x=Age,ymin=CI_LB, ymax=CI_UB), alpha=.50, color="#03396c", fill="#d1e1ec")+
  geom_point(data=dat.prosoc, aes(x=age, y=Eff))+
  geom_errorbar(data=dat.prosoc, aes(x=age, ymin=Eff-sqrt(Eff_var), ymax=Eff+sqrt(Eff_var)))+
  ylab(expression(paste("Effect Size: Pearson's ", italic('r'))))+
  xlab("Sample Mean Age (yrs)")+
  geom_hline(yintercept = 0, lty='dashed')

jpeg(paste0(graphics.folder, 'Prosoc_Moderated_by_Age.jpeg'), res=300, units='in', height = 8, width = 8)
g1
dev.off()

fit.prosoc.mod_female<-rma(yi=Eff, 
                           vi=Eff_var, 
                           mods = ~female, 
                           data=dat.prosoc, 
                           ni=N, 
                           knha=T)
summary(fit.prosoc.mod_female)

fit.prosoc.mod_smpl<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~Samp_typ, 
                         data=dat.prosoc, 
                         ni=N, 
                         knha=T)
summary(fit.prosoc.mod_smpl)

fit.prosoc.mod_CU<-rma(yi=Eff, 
                       vi=Eff_var, 
                       mods = ~CU_resp, 
                       data=dat.prosoc, 
                       ni=N, 
                       knha=T)
summary(fit.prosoc.mod_CU)#Significant (also collinear with outcome rater - only need one model)

CU_resp<-c(0,1)
prosoc.CU.pred<-predict(fit.prosoc.mod_CU, newmods = CU_resp, level = 95)

CU_self<-prosoc.CU.pred$pred[1]
CU_self.se<-prosoc.CU.pred$se[1]

CU_other<-prosoc.CU.pred$pred[2]
CU_other.se<-prosoc.CU.pred$se[2]

CU_self.dist<-rnorm(100000, mean=CU_self, sd=CU_self.se)
CU_other.dist<-rnorm(100000, mean=CU_other, sd=CU_other.se)

prosoc.CU.pred.DF<-cbind(CU_self.dist, CU_other.dist)
colnames(prosoc.CU.pred.DF)<-c('Respondent: Self', 'Respondent: Other')

jpeg(paste0(graphics.folder, 'Prosoc_Moderated_by_Resp.jpeg'), res=300, units='in', height=8, width=8)
mcmc_areas(prosoc.CU.pred.DF, prob = .95)+
  xlab("Effect Size (Fisher's z)")
dev.off()

fit.prosoc.mod_Out<-rma(yi=Eff, 
                        vi=Eff_var, 
                        mods = ~Out_resp, 
                        data=dat.prosoc, 
                        ni=N, 
                        knha=T)
summary(fit.prosoc.mod_Out)

#------------------------------------------------------------------------------------------
#Affective Empathy Model:
table(dat.)
fit.emp_aff.mod_age<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~female+age+CU_resp+Out_resp, 
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

fit.emp_aff.mod_smpl<-rma(yi=Eff, 
                          vi=Eff_var, 
                          mods = ~Samp_typ, 
                          data=dat.Emp_aff, 
                          ni=N, 
                          knha=T)
summary(fit.emp_aff.mod_smpl)

fit.emp_aff.mod_CU<-rma(yi=Eff, 
                        vi=Eff_var, 
                        mods = ~CU_resp, 
                        data=dat.Emp_aff, 
                        ni=N, 
                        knha=T)
summary(fit.emp_aff.mod_CU)

fit.emp_aff.mod_Out<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~Out_resp, 
                         data=dat.Emp_aff, 
                         ni=N, 
                         knha=T)
summary(fit.emp_aff.mod_Out)

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

fit.emp_cog.mod_smpl<-rma(yi=Eff, 
                          vi=Eff_var, 
                          mods = ~Samp_typ, 
                          data=dat.Emp_cog, 
                          ni=N, 
                          knha=T)
summary(fit.emp_cog.mod_smpl)

fit.emp_cog.mod_CU<-rma(yi=Eff, 
                        vi=Eff_var, 
                        mods = ~CU_resp, 
                        data=dat.Emp_cog, 
                        ni=N, 
                        knha=T)
summary(fit.emp_cog.mod_CU)

fit.emp_cog.mod_Out<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~Out_resp, 
                         data=dat.Emp_cog, 
                         ni=N, 
                         knha=T)
summary(fit.emp_cog.mod_Out)#Significant

Out_resp<-c(0,1)
emp_cog.Out.pred<-predict(fit.emp_cog.mod_Out, newmods = Out_resp, level = 95)

Out_self<-emp_cog.Out.pred$pred[1]
Out_self.se<-emp_cog.Out.pred$se[1]

Out_other<-emp_cog.Out.pred$pred[2]
Out_other.se<-emp_cog.Out.pred$se[2]

Out_self.dist<-rnorm(100000, mean=Out_self, sd=Out_self.se)
Out_other.dist<-rnorm(100000, mean=Out_other, sd=Out_other.se)

emp_cog.Out.pred.DF<-cbind(Out_self.dist, Out_other.dist)
colnames(emp_cog.Out.pred.DF)<-c('Respondent: Self', 'Respondent: Other')

jpeg(paste0(graphics.folder, 'CogEmp_Moderated_by_OutResp.jpeg'), res=300, units='in', height=8, width=8)
mcmc_areas(emp_cog.Out.pred.DF, prob = .95)+
  xlab("Effect Size (Fisher's z)")
dev.off()
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

fit.emp_comp.mod_smpl<-rma(yi=Eff, 
                           vi=Eff_var, 
                           mods = ~Sample, 
                           data=dat.emp_comp, 
                           ni=N, 
                           knha=T)
summary(fit.emp_comp.mod_smpl)

fit.emp_comp.mod_CU<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~CU_resp, 
                         data=dat.emp_comp, 
                         ni=N, 
                         knha=T)
summary(fit.emp_comp.mod_CU)

#Exploring marginally significant distributions. 
CU_resp<-c(0,1)
emp_comp.CU.pred<-predict(fit.emp_comp.mod_CU, newmods = CU_resp, level = 95)

CU_self<-emp_comp.CU.pred$pred[1]
CU_self.se<-emp_comp.CU.pred$se[1]

CU_other<-emp_comp.CU.pred$pred[2]
CU_other.se<-emp_comp.CU.pred$se[2]

CU_self.dist<-rnorm(100000, mean=CU_self, sd=CU_self.se)
CU_other.dist<-rnorm(100000, mean=CU_other, sd=CU_other.se)

emp_comp.CU.pred.DF<-cbind(CU_self.dist, CU_other.dist)
colnames(emp_comp.CU.pred.DF)<-c('Respondent: Other', 'Respondent: Self')

jpeg(paste0(graphics.folder, 'CompEmp_Moderated_by_CUResp.jpeg'), res=300, units='in', height=8, width=8)
mcmc_areas(emp_comp.CU.pred.DF, prob = .95)+
  xlab("Effect Size (Fisher's z)")
dev.off()


fit.emp_comp.mod_Out<-rma(yi=Eff, 
                          vi=Eff_var, 
                          mods = ~Out_resp, 
                          data=dat.emp_comp, 
                          ni=N, 
                          knha=T)
summary(fit.emp_comp.mod_Out)

###########################################################################################################
#Testing whether ICU mesure used is a moderator... 
dat.icu<-read.csv(paste0(data.folder, 'descriptive table.csv'))
colnames(dat.icu)[1]<-'id'

##
dat.Emp_tot<-merge(dat.Emp_tot, dat.icu[,c(1,10)], by='id')
fit.emp_tot.mod_icu<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~ICU, 
                         data=dat.Emp_tot, 
                         ni=N, 
                         knha=T)
summary(fit.emp_tot.mod_icu)

##
dat.Emp_aff<-merge(dat.Emp_aff, dat.icu[,c(1,10)], by='id')
fit.emp_aff.mod_icu<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~ICU, 
                         data=dat.Emp_aff, 
                         ni=N, 
                         knha=T)
summary(fit.emp_aff.mod_icu)

##
dat.Emp_cog<-merge(dat.Emp_cog, dat.icu[,c(1,10)], by='id')
fit.emp_cog.mod_icu<-rma(yi=Eff, 
                         vi=Eff_var, 
                         mods = ~ICU, 
                         data=dat.Emp_cog, 
                         ni=N, 
                         knha=T)
summary(fit.emp_cog.mod_icu)

##
dat.prosoc<-merge(dat.prosoc, dat.icu[,c(1,10)], by='id')
fit.prosoc.mod_icu<-rma(yi=Eff, 
                        vi=Eff_var, 
                        mods = ~ICU, 
                        data=dat.prosoc, 
                        ni=N, 
                        knha=T)
summary(fit.prosoc.mod_icu)

##
dat.emp_comp<-merge(dat.emp_comp, dat.icu[,c(1,10)], by='id')
fit.emp_comp.mod_icu<-rma(yi=Eff, 
                          vi=Eff_var, 
                          mods = ~ICU, 
                          data=dat.emp_comp, 
                          ni=N, 
                          knha=T)
summary(fit.emp_comp.mod_icu)

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
library(metaSEM)
library(RColorBrewer)

#With Total Empathy Scores (this seems totally crazy to me...):
dat.meta3<-dat.AC[dat.AC$Emp_cog!=1,]
dat.meta3<-dat.meta3[dat.meta3$Emp_aff!=1,]

#Total of 49 effects left..

y<-dat.meta3$Eff
v<-dat.meta3$Eff_var
ID<-dat.meta3$id    #Note that this currently still has some within study break down by gender... (grrrrrrrrrr)

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
summary(fit.meta_null)
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
                  sigma2 = sigma2)

sink(paste0(model.folder, 'Three-level_VarDecomp.txt'))
print(var.decomp)
sink()

meta_Emp_tot.var<-as.data.frame(var.decomp$Decompositions)
meta_Emp_tot.var<-meta_Emp_tot.var[1:5,]

meta_Emp_tot.var1<-data.frame(Component = rep(rownames(meta_Emp_tot.var), 3), 
                              Variance = c(as.numeric(meta_Emp_tot.var[,1]), 
                                           as.numeric(meta_Emp_tot.var[,2]), 
                                           as.numeric(meta_Emp_tot.var[,3])
                              ), 
                              Denom = c(rep('Total', 5), 
                                        rep('Within', 5), 
                                        rep('Between', 5)
                              )
                            Variance = c(as.numeric(meta_Emp_tot.var[,1]), 
                                         as.numeric(meta_Emp_tot.var[,2]), 
                                         as.numeric(meta_Emp_tot.var[,3])
                            ), 
                            Denom = c(rep('Total', 5), 
                                      rep('Within', 5), 
                                      rep('Between', 5)
                            )
)

rownames(meta_Emp_tot.var1)<-NULL

#This took a lot of messing around - not sure I understand the logic of the reordering

meta_Emp_tot.var1$Component<-factor(meta_Emp_tot.var1$Component, levels = c("sigma2",
                                                                            "mean variation",
                                                                            "fixed, between",
                                                                            "fixed, within",
                                                                            "slope variation")
)

meta_Emp_tot.var1$Denom<-factor(meta_Emp_tot.var1$Denom, levels = c('Total', 
                                                                    'Within', 
                                                                    'Between')
                                                                        "mean variation",
                                                                        "fixed, between",
                                                                        "fixed, within",
                                                                        "slope variation")
)

meta_Emp_tot.var1$Denom<-factor(meta_Emp_tot.var1$Denom, levels = c('Total', 
                                                                'Within', 
                                                                'Between')
)

myColors <- brewer.pal(5,"Set3")
names(myColors) <- levels(meta_Emp_tot.var1$Component)

myColors[1]<-"#D3D3D3"
myColors[2]<-"#808080"

colScale <- scale_colour_manual(name = "Component",values = myColors)
fillScale <- scale_fill_manual(name = "Component",values = myColors)

g1<-ggplot(meta_Emp_tot.var1)+
  geom_bar(aes(y=Variance, 
               fill=Component,
               color=Component, 
               x=Denom
  ), 
  stat = 'identity', 
  position = 'stack')+
  colScale+
  fillScale+
  ggtitle('Three-level Variance Decomposition')+
  xlab('Source of Variability')+
  ylab('Proportion of Variance Explained')

png(paste0(study2.graphics, '/S1_JOYMood_Decomposition.png'), 
    res=300,
    height = 8, 
    width = 8, 
    units = 'in')
g1
dev.off()

#------------------------------------------------------------
#Modeling prooscialty
fit.meta_prosoc<-meta3(y=y,
                       v=v,
                       cluster = ID,
                       x = x[,2:3]
)

sink(paste0(model.folder, 'Three-level_Prosocial_Model.txt'))
summary(fit.meta_prosoc)
sink()

#Modeling glt
fit.meta_glt<-meta3(y=y,
                    v=v,
                    cluster = ID,
                    x = x[,c(1,3)]
)

sink(paste0(model.folder, 'Three-level_Guilt_Model.txt'))
summary(fit.meta_glt)
sink()

fit.meta_0int<-meta3(y=y,
                     v=v,
                     cluster = ID,
                     x = x, 
                     intercept.constraints = 0
)

summary(fit.meta_0int)  

#Adding additional covariates: 
x2<-cbind(x, dat$age, dat$female, dat$CU_resp, dat$Out_resp, dat$Samp_typ)
colnames(x2)[6:10]<-c('age', 'female', 'CU_resp', 'Out_resp', 'Samp_typ')

fit.table2<-cbind(as.matrix(coef(fit.meta2)[1:10]), confint(fit.meta2)[1:10,])
row.names(fit.table2)<-colnames(x2)
colnames(fit.table2)<-c("Estimate", "Lower Bound", "Upper Bound")
fit.table2

#Separate Empathy Total from Empathy Affect/Cog [Empathy overall bigger focus]
