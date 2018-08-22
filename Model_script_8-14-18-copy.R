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

#Tables 2-5 (created for the supplement)
#Add something about the join model - MetaSEM approach. 
#Results - 
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
dat<-read.csv(paste0(wd, 'Model_Data_w-Reliability_082218.csv'), stringsAsFactors = F)
colnames(dat)[1]<-'id'

#################################################################################
#Or can load latest data from GitHub
load(paste0(wd, '/Meta_Raw_Data/Updated_Data_081318.RData'))

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

#Will stick with the original variance estimates using this approach however
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
sink(paste0(model.folder, 'Reg_test_Total_Empathy.txt'))
print(REG.emp_tot)
sink()

#Trim and fill models
fit.emp_tot.TF_R<-trimfill(fit.emp_tot, estimator = 'R0')
fit.emp_tot.TF_R

sink(paste0(model.folder, 'Total_Empathy_R0.txt'))
summary(fit.emp_tot.TF_R)
sink()

jpeg(paste0(graphics.folder, 'R0_estimator_CU-tot_emp_Sup1A.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.emp_tot.TF_R)
title(expression('Trim and Fill Results for Total Empathy (R'[0]*' Estimator)'))
dev.off()

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
sink(paste0(model.folder, 'Reg_test_Affective_Empathy.txt'))
print(REG.emp_aff)
sink()

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
sink(paste0(model.folder, 'Reg_test_Cognitive_Empathy.txt'))
print(REG.emp_cog)
sink()

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

fit.emp_cog.TF_L<-trimfill(fit.emp_cog, estimator = 'L0')
fit.emp_cog.TF_L

sink(paste0(model.folder, 'Cognitive_Empathy_L0.txt'))
summary(fit.emp_cog.TF_L)
sink()

jpeg(paste0(graphics.folder, 'L0_estimator_CU-cog_emp_Sup3B.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.emp_cog.TF_L)
title(expression('Trim and Fill Results for Cognitive Empathy (L'[0]*' Estimator)'))
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
sink(paste0(model.folder, 'Reg_test_Prosociality.txt'))
print(REG.prosoc)
sink()

#Trim and fill model
fit.prosoc.TF_R<-trimfill(fit.prosoc, estimator = 'R0')
fit.prosoc.TF_R

sink(paste0(model.folder, 'Prosocial_R0.txt'))
summary(fit.prosoc.TF_R)
sink()

#Some of these "missing" effects likely due more to LB of correlation 
#Not obvious that the asymmetry in these plots can really be interpreted 
#(at least when missing effects are to the left)
jpeg(paste0(graphics.folder, 'R0_estimator_CU-prosoc_Sup5A.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.prosoc.TF_R)
title(expression('Trim and Fill Results for Prosocialty (R'[0]*' Estimator)'))
dev.off()

fit.prosoc.TF_L<-trimfill(fit.prosoc, estimator = 'L0')
fit.prosoc.TF_L

sink(paste0(model.folder, 'Prosocial_L0.txt'))
summary(fit.prosoc.TF_L)
sink()

jpeg(paste0(graphics.folder, 'L0_estimator_CU-prosoc_Sup5B.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.prosoc.TF_L)
title(expression('Trim and Fill Results for Prosocialty (L'[0]*' Estimator)'))
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

cols<-c(1,16,2:4,9,24,25,17,18,13,33,48,49)
dat.emp_comp<-dat.emp_comp[,cols]
colnames(dat.emp_comp)<-c('id', 'citation', 'female', 'age', 'N', 'R_affective',
                          'Eff_affective', 'Eff_var_affective', 'CU_resp', 'Out_resp',
                          'Sample', 'R_cognitive', 'Eff_cognitive', 'Eff_var_cognitive')

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

sink(paste0(model.folder, 'Cognitive_Affective_Empathy_R0.txt'))
summary(fit.emp_comp.TF_R)
sink()

jpeg(paste0(graphics.folder, 'R0_estimator_CU-comp_emp_Sup4A.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.emp_comp.TF_R)
title(expression('Trim and Fill Results for Cognitive vs Affective Empathy (R'[0]*' Estimator)'))
dev.off()

fit.emp_comp.TF_L<-trimfill(fit.emp_comp, estimator = 'L0')
fit.emp_comp.TF_L

sink(paste0(model.folder, 'Cognitive_Affective_Empathy_L0.txt'))
summary(fit.emp_comp.TF_L)
sink()

jpeg(paste0(graphics.folder, 'L0_estimator_CU-comp_emp_Sup4B.jpeg'), res = 300, width = 7, height = 7, units = 'in')
funnel(fit.emp_comp.TF_L)
title(expression('Trim and Fill Results for Cognitive vs Affective Empathy (L'[0]*' Estimator)'))
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

jpeg(paste0(graphics.folder, 'Unconditional Model Summary Graphic.jpeg'), res=300, units='in', height=7, width=7)
mcmc_areas(DF.plot.dist2, prob=.95)+
  xlab(expression('Population Esimates of'~rho))+
  ggtitle('Simulated Distributions of Effect Sizes')
dev.off()

###########################################################################################
#Unconditional Models (no moderators - publication graphics)
#Fitting two summary graphics (one for affective vs. cognitive) & one for summary vals
#Code for graphic models was adapted from: http://www.metafor-project.org/doku.php/plots:forest_plot_with_subgroups
dat.graph1<-dat.imp
dat.graph1<-dat.graph1[dat.graph1$Emp_aff==1 | dat.graph1$Emp_cog==1,]

dat.graph1$cite.fac<-as.factor(dat.graph1$citation)
fit.graph1<-rma(yi=Eff, vi=Eff_var, data=dat.graph1)

#Adding back in Outcome as a variable
dat.graph1$Outcome[dat.graph1$Emp_aff==1]<-'Affective Empathy'
dat.graph1$Outcome[dat.graph1$Emp_cog==1]<-'Cognitive Empathy'
dat.graph1$female<-round(dat.graph1$female, digits = 2)
dat.graph1$age<-round(dat.graph1$age, digits = 2)

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
       xlab=expression(~rho), 
       mlab="", 
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
dat.graph2<-dat.imp

dat.graph2<-dat.graph2[dat.graph2$Emp_aff!=1,]
dat.graph2<-dat.graph2[dat.graph2$Emp_cog!=1,]
dat.graph2$citation[dat.graph2$citation=='Anastassiou-Hadjicharalambous & Warden (2008)']<-'Anastassiou-H. & Warden (2008)'

dat.graph2$cite.fac<-as.factor(dat.graph2$citation)
fit.graph2<-rma(yi=Eff, vi=Eff_var, data=dat.graph2)

#round down age and % female for purposes of plotting
dat.graph2$female<-round(dat.graph2$female, digits = 2)
dat.graph2$age<-round(dat.graph2$age, digits = 2)
dat.graph2$Outcome[dat.graph2$prosocial==1]<-'prosoc'
dat.graph2$Outcome[dat.graph2$guilt==1]<-'guilt'
dat.graph2$Outcome<-ifelse(!is.na(dat.graph2$Outcome), dat.graph2$Outcome, 'emp_tot')

#Note to self - adding spaces into variables - even string variables makes ordering of effects difficult on 
#accompanying table/graphic

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
       xlab = expression(~rho), 
       mlab="", 
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

#######################################################################################
#Exploring moderators - Sample Type
fit.emp_tot.mod_samp<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var, 
                         mods = ~Samp_typ, 
                         data=dat.Emp_tot, 
                         ni=N, 
                         knha=T)
summary(fit.emp_tot.mod_samp)

#Prosocial Model:
fit.prosoc.mod_samp<-rma(yi=Eff, 
                        vi=Eff_var,
                        weights = 1/Eff_var,
                        mods = ~Samp_typ, 
                        data=dat.prosoc, 
                        ni=N, 
                        knha=T)
summary(fit.prosoc.mod_samp)  

#Affective Empathy:
fit.Emp_aff.mod_samp<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~Samp_typ, 
                         data=dat.Emp_aff, 
                         ni=N, 
                         knha=T)
summary(fit.Emp_aff.mod_samp)  

#######################################################################################
#Exploring moderators - CU_respondent
fit.emp_tot.mod_CU_resp<-rma(yi=Eff, 
                          vi=Eff_var,
                          weights = 1/Eff_var, 
                          mods = ~CU_resp, 
                          data=dat.Emp_tot, 
                          ni=N, 
                          knha=T)
summary(fit.emp_tot.mod_CU_resp)

#Prosocial Model:
fit.prosoc.mod_CU_resp<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~CU_resp, 
                         data=dat.prosoc, 
                         ni=N, 
                         knha=T)
summary(fit.prosoc.mod_CU_resp)  

#Affective Empathy:
fit.Emp_aff.mod_CU_resp<-rma(yi=Eff, 
                          vi=Eff_var,
                          weights = 1/Eff_var,
                          mods = ~CU_resp, 
                          data=dat.Emp_aff, 
                          ni=N, 
                          knha=T)
summary(fit.Emp_aff.mod_CU_resp)  

#Cognitive Empathy Model:
fit.Emp_cog.mod_CU_resp<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~CU_resp, 
                         data=dat.Emp_cog, 
                         ni=N, 
                         knha=T)
summary(fit.Emp_cog.mod_CU_resp)  

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
######################################
dat.Emp_tot$citation[dat.Emp_tot$citation=='Anastassiou-Hadjicharalambous & Warden (2008)']<-'Anastassiou-H. & Warden (2008)'
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

#Note that some studies have multiple effects - have not partialed separted these (thought it would be possible)
#Will keep graphic for GitHub - but will not include in supplement without separating out within study effects 

jpeg(paste0(graphics.folder,'Attenuation Correction vs Raw Effects.jpeg'), res=300, units = 'in', height = 8, width=20)
cowplot::plot_grid(g.Emp_tot, g.Emp_aff, g.Emp_cog, g.prosoc, 
                        labels = c('A', 'B', 'C', 'D'))
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
fit.meta_0int<-meta3(y=y,
                     v=v,
                     cluster = ID,
                     x = x, 
                     intercept.constraints = 0
)

summary(fit.meta_0int)  

#Getting all pairwise comparisons...
fit.meta_glt<-meta3(y=y,
                    v=v,
                    cluster = ID,
                    x = x[,c(1,3)]
)
summary(fit.meta_glt)

#Quick summary - in the joint model: 
#     1. All three significantly predict CU scores
#     2. Prosociality is more strongly correlated with CU scores that either other covariate
#     3. There is no difference between guilt and total empathy in their relations with CU 

