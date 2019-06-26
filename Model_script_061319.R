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
wd<-'~/GitHub/Meta-Analysis_Waller_et_al'
data.folder<-paste0(wd, '/Meta_Raw_Data')
model.folder<-paste0(wd, '/Output')
graphics.folder<-paste0(wd, '/Graphics_Folder' )

#################################################################################
#Notes:   1. Correlations were transformed to Fisher's z
#         2. Standardized mean differences were transformed to r then to Fisher's z
#         3. Will require that all results are exponentiated to return to -1:1 scale
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
dat_raw<-read.csv(paste0(data.folder, '/Updated_Effect_Data_061119.csv'))
head(dat_raw) #Column 34 looks mainly empty (NAs)

dat_raw<-dat_raw[,-34]  #removing empty column

#Imputing missing reliabilities prior to recoding. 
psych::describe(dat_raw)

#Also extra rows in the data
dat_raw<-dat_raw[1:106,]

#Working data set: 
dat <- dat_raw

#Need to create a series of binary variables - cannot handle categorical variables
dat$Emp_tot<-ifelse(dat$Outcome=='empathy_tot', 1, 0)
table(dat$Emp_tot)

dat$Emp_aff<-ifelse(dat$Outcome=='empathy_aff', 1, 0)
table(dat$Emp_aff)

dat$Emp_cog<-ifelse(dat$Outcome=='empathy_cog', 1, 0)
table(dat$Emp_cog)

dat$prosocial<-ifelse(dat$Outcome=='prosocial', 1, 0)
table(dat$prosocial)

dat$guilt<-ifelse(dat$Outcome=='guilt', 1, 0)
table(dat$guilt)

#selecting total empathy as reference value - largest N
dat.imp<-dat[c('Id', 
               'female', 
               'age', 
               'N', 
               'Emp_aff', 
               'Emp_cog', 
               'prosocial',
               'guilt',
               'R', 
               'CU_resp', 
               'Out_resp', 
               'Samp_typ', 
               'CU_Rel', 
               'Out_Rel')]

#Imputation Model - imputing missing reliabilities before the calculation of attenuation corrected scores
ini <- mice(dat.imp, maxit = 0)
ini$nmis  #Missing 12 CU reliabilities and 13 outcome reliabilities

imp<-mice(dat.imp[,2:ncol(dat.imp)], maxit=100, m=5)

png(paste0(graphics.folder, '/MICE_Imputation_convergence.png'), 
    res=300, 
    units='in', 
    height = 4, 
    width = 6)
plot(imp)
dev.off()

#Extracting and obtaining estimates for missing values (averaging across the 5 imputed datasets)
dat.imp<-complete(imp)
dat.imp$citation<-dat$citation
dat.imp$Id <- dat$Id

dat.imp$R_var<-((1-dat$R^2)^2)/(dat$N-1)  #Getting correlation variance 
                                          #Needed for attenuation correction later on

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

#Correcting d's for attentuation
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
sink(paste0(model.folder, '/Imputed_Descriptives.txt'))
psych::describe(dat.imp)
sink()

dat_descrip<-read.csv(paste0(data.folder, '/descriptive_060519.csv'))
colnames(dat_descrip)[1]<-"Id"
dat_descrip$ICU<-as.numeric(str_detect(dat_descrip$CU.traits.measure, 'ICU'))
table(dat_descrip$ICU)

dat.imp<-merge(dat.imp, dat_descrip[,c("Id", "ICU")], by = "Id")

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

sink(paste0(model.folder, '/Total_Empathy_AC.txt'))
summary(fit.emp_tot)
sink()

tiff(paste0(graphics.folder, '/CU and Total Empathy - Forest.tiff'), res=900, width = 7, height=7, units='in')
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

summary(fit.emp_tot_bias)

REG.emp_tot<-regtest(fit.emp_tot_bias, model = 'lm')
sink(paste0(model.folder, '/Reg_test_Total_Empathy.txt'))
print(REG.emp_tot)
sink()

fit.emp_tot.TF_L<-trimfill(fit.emp_tot_bias, estimator = 'L0', ilim = c(-1,1))
fit.emp_tot.TF_L

sink(paste0(model.folder, '/Total_Empathy_L0.txt'))
summary(fit.emp_tot.TF_L)
sink()

#Funnel plot for uncorrected correlations: 
tiff(paste0(graphics.folder, '/CU and Total Empathy - Funnels.tiff'), res = 300, width = 11, height = 8, units = 'in')
par(mfrow=c(1,2))
funnel(fit.emp_tot_bias, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Raw Total Empathy Correlations'))
funnel(fit.emp_tot.TF_L, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Trim and Fill Results for Total Empathy (L'[0]*' Estimator)'))

dev.off()
par(mfrow=c(1,1))

#################################################################################
#Model for Affective Empathy: 
fit.emp_aff<-rma(yi=Eff, 
                 vi=Eff_var, 
                 weights = 1/Eff_var, 
                 data=dat.Emp_aff, 
                 ni=N,
                 method = 'HS')

summary(fit.emp_aff)

sink(paste0(model.folder, '/Affective_Empathy_AC.txt'))
summary(fit.emp_aff)
sink()

tiff(paste0(graphics.folder, '/CU and Affective Empathy - Forest.tiff'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_aff, order = 'obs', slab=dat.Emp_aff$citation)
title("Relation between CU Traits and Affective Empathy")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
fit.emp_aff_bias<-rma(yi=R, 
                      vi=R_var,
                      data=dat.Emp_aff, 
                      ni=N)

summary(fit.emp_aff_bias)

REG.emp_aff<-regtest(fit.emp_aff_bias, model = 'lm')
sink(paste0(model.folder, '/Reg_test_Affective_Empathy.txt'))
print(REG.emp_aff)
sink()

fit.emp_aff.TF_L<-trimfill(fit.emp_aff_bias, estimator = 'L0', ilim = c(-1,1))
fit.emp_aff.TF_L

sink(paste0(model.folder, '/Affective_Empathy_L0.txt'))
summary(fit.emp_aff.TF_L)
sink()

#Funnel plot for uncorrected correlations: 
tiff(paste0(graphics.folder, '/CU and Affective Empathy - Funnels.tiff'), res = 300, width = 11, height = 8, units = 'in')
par(mfrow=c(1,2))
funnel(fit.emp_aff_bias, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Trim and Fill Results for Affective Empathy Correlations'))
funnel(fit.emp_aff.TF_L, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Trim and Fill Results for Affective Empathy (L'[0]*' Estimator)'))

dev.off()
par(mfrow=c(1,1))

#################################################################################
#Model for Cognitive Empathy: 
fit.emp_cog<-rma(yi=Eff, 
                 vi=Eff_var, 
                 weights = 1/Eff_var, 
                 data=dat.Emp_cog, 
                 ni=N,
                 method = 'HS')

summary(fit.emp_cog)

sink(paste0(model.folder, '/Cognitive_Empathy_AC.txt'))
summary(fit.emp_cog)
sink()

tiff(paste0(graphics.folder, '/CU and Cognitive Empathy - Forest.tiff'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_cog, order = 'obs', slab=dat.Emp_cog$citation)
title("Relation between CU Traits and Cognitive Empathy")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
fit.emp_cog_bias<-rma(yi=R, 
                      vi=R_var,
                      data=dat.Emp_cog, 
                      ni=N)

summary(fit.emp_cog_bias)

REG.emp_cog<-regtest(fit.emp_cog_bias, model = 'lm')
sink(paste0(model.folder, '/Reg_test_Cognitive_Empathy.txt'))
print(REG.emp_cog)
sink()

fit.emp_cog.TF_L<-trimfill(fit.emp_cog_bias, estimator = 'L0', ilim = c(-1,1))
fit.emp_cog.TF_L

sink(paste0(model.folder, '/Cognitive_Empathy_L0.txt'))
summary(fit.emp_cog.TF_L)
sink()

#Funnel plot for uncorrected correlations: 
tiff(paste0(graphics.folder, '/CU and Cognitive Empathy - Funnels.tiff'), res = 300, width = 11, height = 8, units = 'in')
par(mfrow=c(1,2))
funnel(fit.emp_cog_bias, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Raw Cognitive Empathy Correlations'))
funnel(fit.emp_cog.TF_L, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Trim and Fill Results for Cognitive Empathy (L'[0]*' Estimator)'))

dev.off()
par(mfrow=c(1,1))

#------------------------------------------------------------------------------------------
#Model for Prosocialty: 
fit.prosoc<-rma(yi=Eff, 
                vi=Eff_var, 
                weights = 1/Eff_var, 
                data=dat.prosoc, 
                ni=N,
                method = 'HS')

summary(fit.prosoc)

sink(paste0(model.folder, '/Prosocial_AC.txt'))
summary(fit.prosoc)
sink()

tiff(paste0(graphics.folder, '/CU and Prosocialty - Forest.tiff'), res=300, width = 7, height=7, units='in')
forest.rma(fit.prosoc, order = 'obs', slab=dat.prosoc$citation)
title("Relation between CU Traits and Prosocialty")
dev.off()

#---------------------------------------------------------------------------------
#testing for publication bias 
fit.prosoc_bias<-rma(yi=R, 
                      vi=R_var,
                      data=dat.prosoc, 
                      ni=N)
summary(fit.prosoc_bias)

REG.prosoc<-regtest(fit.prosoc_bias, model = 'lm')
sink(paste0(model.folder, '/Reg_test_Prosociality.txt'))
print(REG.prosoc)
sink()

fit.prosoc.TF_L<-trimfill(fit.prosoc_bias, estimator = 'L0', ilim = c(-1,1))
fit.prosoc.TF_L

sink(paste0(model.folder, '/Prosociality_L0.txt'))
summary(fit.prosoc.TF_L)
sink()

#Funnel plot for uncorrected correlations: 
tiff(paste0(graphics.folder, '/CU and Prosociality - Funnels.tiff'), res = 300, width = 11, height = 8, units = 'in')
par(mfrow=c(1,2))
funnel(fit.prosoc_bias, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Raw Prosociality Correlations'))
funnel(fit.prosoc.TF_L, 
       xlab = expression(paste('Uncorrected Correlations (', italic('r'), ')')), 
       main =  expression('Trim and Fill Results for Prosociality (L'[0]*' Estimator)'))

dev.off()
par(mfrow=c(1,1))

#################################################################################
#Model for Guilt: 
fit.glt<-rma(yi=Eff, 
             vi=Eff_var, 
             weights = 1/Eff_var, 
             data=dat.glt, 
             ni=N,
             method = 'FE')
summary(fit.glt)

sink(paste0(model.folder, '/Guilt Overall Model - no moderators.txt'))
summary(fit.glt)
sink()

tiff(paste0(graphics.folder, '/CU and Guilt - Forest.tiff'), res=300, width = 7, height=7, units='in')
forest.rma(fit.glt, order = 'obs', slab=dat.glt$citation)
title("Relation between CU Traits and Guilt")
dev.off()

fit.glt2<-rma(yi=R, 
             vi=R_var, 
             data=dat.glt, 
             ni=N,
             method = 'FE')
summary(fit.glt2)

###########################################################################################
#Model comparing two forms of empathy - within study differences in effects
#merging and cleaning data set
dat.emp_comp<-merge(dat.Emp_aff, dat.Emp_cog, by=c('Id', 
                                                   'citation', 
                                                   'female', 
                                                   'age', 
                                                   'N', 
                                                   'ICU'))

cols<-c(1:6, 11:22, 25:26, 31, 36:37, 45:46)
colnames(dat.emp_comp)[cols]
dat.emp_comp<-dat.emp_comp[,cols]
colnames(dat.emp_comp)<-c('Id', 'citation', 'female', 'age', 'N', 'ICU', 'R_Aff', 
                          'CU_resp', 'Out_resp', 'Samp_typ', 'CU_Rel', 'Aff_Rel', 'R_var_Aff',
                          'CU_resp.R', 'Out_resp.R', 'Same_diff', 'Same_diff.R', 
                          'Samp_typ.R', 'Eff_Aff', 'Eff_var_Aff', 'R_cog', 'Cog_Rel', 'R_var_Cog', 
                          'Eff_Cog', 'Eff_var_Cog')

#Correlations pulled from studies when available: 
#Set to 0 when unavailable or unreported  
dep_cor<-matrix(c("Dadds et al. (2012)", 0,       
                  "Georgiou et al. (2018)", .13,    
                  "Gillen et al. (2018)", 0,       
                  "Kahn et al. (2017)", -0.01,        
                  "Kimonis et al. (2016)", 0,      
                  "Liu et al. (2018)", .70,          
                  "Lui et al. (2016)", .08,         
                  "McDonald et al. (2018)", -.52,    
                  "Munoz et al. (2011)", .50,        
                  "O'Kearney et al. (2017)", 0, 
                  "Antoniadou et al. (2016)", .76,  
                  "Pasalich et al. (2014)", .32,    
                  "Pechorro et al. (2015b)", .49,   
                  "Pechorro et al. (2016)", 0,     
                  "Pechorro et al. (2017)", 0,     
                  "Pijper et al. (2016)", .45,      
                  "Raine and Chen (2018)", 0,     
                  "Brouns et al. (2013)", 0,       
                  "Brouns et al. (2013)", 0,       
                  "van Vugt et al. (2012)", 0,     
                  "Ciucci & Baroncelli (2014)", .35,
                  "Dadds et al. (2009)", .063,        
                  "Dadds et al. (2009)", .063), 
                byrow = TRUE, ncol = 2)

colnames(dep_cor)<-c("citation", "R_Aff_Cog")
dep_cor<-as.data.frame(dep_cor)
dep_cor$R_Aff_Cog<-as.numeric(as.character(dep_cor$R_Aff_Cog))

#Calculating difference in effects and variance terms adjusted for outcome correlation
Eff<-dat.emp_comp$Eff_Aff-dat.emp_comp$Eff_Cog
Eff_var<-dat.emp_comp$Eff_var_Aff+dat.emp_comp$Eff_var_Cog-2*dep_cor$R_Aff_Cog*sqrt(dat.emp_comp$Eff_var_Aff*dat.emp_comp$Eff_var_Cog)

R<-dat.emp_comp$Eff_Aff-dat.emp_comp$R_Aff-dat.emp_comp$R_cog
R_var<-dat.emp_comp$R_var_Aff+dat.emp_comp$R_var_Cog-2*dep_cor$R_Aff_Cog*sqrt(dat.emp_comp$R_var_Aff*dat.emp_comp$R_var_Cog)

dat.emp_comp$Eff<-Eff
dat.emp_comp$Eff_var<-Eff_var

dat.emp_comp$R<-R
dat.emp_comp$R_var<-R_var

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

sink(paste0(model.folder, '/Difference in Empathy Dimension Model - no moderators.txt'))
summary(fit.emp_comp)
sink()

tiff(paste0(graphics.folder, '/Difference in Empathy Model - Forest.tiff'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_comp, order = 'obs', slab=dat.emp_comp$citation)
title("Difference in Relation between CU and Empathy Dimensions")
dev.off()

fit.emp_comp_bias<-rma(yi=R, 
                       vi=R_var,
                       data=dat.emp_comp, 
                       ni=N)

summary(fit.emp_comp_bias)

#forest.rma(fit.emp_comp_bias, order = 'obs', slab=dat.emp_comp$citation)
#title("Difference in Relation between CU and Empathy Dimensions")

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
                       k=c(nrow(dat.Emp_tot), 
                           nrow(dat.Emp_aff), 
                           nrow(dat.Emp_cog), 
                           nrow(dat.prosoc), 
                           nrow(dat.glt), 
                           nrow(dat.emp_comp)),
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
          out = paste0(model.folder,'/Unconditional Models.txt'), 
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

tiff(paste0(graphics.folder, '/Unconditional Model Summary Graphic.tiff'), res=300, units='in', height=7, width=7)
mcmc_areas(DF.plot.dist2, prob=.95)+
  xlab(expression('Population Esimates of'~rho))+
  ggtitle('Simulated Distributions of Effect Sizes')+theme_bw()+
  geom_vline(xintercept = 0, lty = 'dashed', color = 'red')
dev.off()

###########################################################################################
#Unconditional Models (no moderators - publication graphics)
#Fitting two summary graphics (one for affective vs. cognitive) & one for summary vals
#Code for graphic models was adapted from: http://www.metafor-project.org/doku.php/plots:forest_plot_with_subgroups
dat.graph1<-dat.imp
dat.graph1<-dat.graph1[rowSums(dat.graph1[,5:8])==0,]
dat.graph1$citation<-as.character(dat.graph1$citation)
dat.graph1$citation[dat.graph1$citation=="Raine and Chen (2018)"]<-"Raine & Chen (2018)"  

dat.graph1$cite.fac<-as.factor(dat.graph1$citation)

#Fixing formatting of citation for consistency
fit.graph1<-rma(yi=Eff, vi=Eff_var, data=dat.graph1)

#Adding back in Outcome as a variable
dat.graph1$female<-round(dat.graph1$female, digits = 2)
dat.graph1$age<-round(dat.graph1$age, digits = 2)
par(mar=c(4,4,1,2))

tiff(paste0(graphics.folder,'/Total Empathy.tiff'), res=1200, units = 'in', height = 7, width=10)
par(cex=1, font=1)
forest(fit.graph1, xlim=c(-10, 2), 
       order = order(dat.graph1$Eff), 
       ilab = cbind(dat.graph1$N, 
                    dat.graph1$female, 
                    dat.graph1$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:(nrow(dat.graph1)+2)),
       ylim = c(-1, nrow(dat.graph1)+5), 
       cex=.75, 
       xlab=expression(~rho), 
       mlab="", 
       addfit = F,
       slab = dat.graph1$cite.fac)

par(cex=1, font=4)
par(font=4)
text(-7, nrow(dat.graph1)+4, 'N')
par(font=2)
text(c(-5, -3), nrow(dat.graph1)+4, c('%Female', 'Mean Age'))
text(-10, nrow(dat.graph1)+4, 'Citation', pos=4)
text(1.25, nrow(dat.graph1)+4, expression(paste(bolditalic('r'), '[LB,UB]')))

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
dat.graph2$citation<-as.character(dat.graph2$citation)
dat.graph2$citation[dat.graph2$citation=="Raine and Chen (2018)"]<-"Raine & Chen (2018)"  

dat.graph2$cite.fac<-as.factor(dat.graph2$citation)
fit.graph2<-rma(yi=Eff, vi=Eff_var, data=dat.graph2)

#Adding back in Outcome as a variable
dat.graph2$female<-round(dat.graph2$female, digits = 2)
dat.graph2$age<-round(dat.graph2$age, digits = 2)
par(mar=c(4,4,1,2))

tiff(paste0(graphics.folder,'/Affective Empathy.tiff'), res=1200, units = 'in', height = 7, width=10)
par(cex=1, font=1)
forest(fit.graph2, xlim=c(-10, 2), 
       order = order(dat.graph2$Eff), 
       ilab = cbind(dat.graph2$N, 
                    dat.graph2$female, 
                    dat.graph2$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:(nrow(dat.graph2)+2)),
       ylim = c(-1,nrow(dat.graph2)+5), 
       cex=.75, 
       xlab=expression(~rho), 
       mlab="", 
       addfit = F,
       slab = dat.graph2$cite.fac)

par(cex=1, font=4)
par(font=4)
text(-7, nrow(dat.graph2)+4, 'N')
par(font=2)
text(c(-5, -3), nrow(dat.graph2)+4, c('%Female', 'Mean Age'))
text(-10, nrow(dat.graph2)+4, 'Citation', pos=4)
text(1.25, nrow(dat.graph2)+4, expression(paste(bolditalic('r'), '[LB,UB]')))

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
dat.graph3$citation<-as.character(dat.graph3$citation)
dat.graph3$citation[dat.graph3$citation=="Raine and Chen (2018)"]<-"Raine & Chen (2018)"  

dat.graph3$cite.fac<-as.factor(dat.graph3$citation)
fit.graph3<-rma(yi=Eff, vi=Eff_var, data=dat.graph3)

#Adding back in Outcome as a variable
dat.graph3$female<-round(dat.graph3$female, digits = 2)
dat.graph3$age<-round(dat.graph3$age, digits = 2)
par(mar=c(4,4,1,2))

tiff(paste0(graphics.folder,'/Cognitive Empathy.tiff'), res=1200, units = 'in', height = 7, width=10)
par(cex=1, font=1)
forest(fit.graph3, xlim=c(-10, 2), 
       order = order(dat.graph3$Eff), 
       ilab = cbind(dat.graph3$N, 
                    dat.graph3$female, 
                    dat.graph3$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:(nrow(dat.graph3)+2)),
       ylim = c(-1, nrow(dat.graph3)+5), 
       cex=.75, 
       xlab=expression(~rho), 
       mlab="", 
       addfit = F,
       slab = dat.graph3$cite.fac)

par(cex=1, font=4)
par(font=4)
text(-7, nrow(dat.graph3)+4, 'N')
par(font=2)
text(c(-5, -3), nrow(dat.graph3)+4, c('%Female', 'Mean Age'))
text(-10, nrow(dat.graph3)+4, 'Citation', pos=4)
text(1.25, nrow(dat.graph3)+4, expression(paste(bolditalic('r'), '[LB,UB]')))

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
dat.graph4$citation<-as.character(dat.graph4$citation)

dat.graph4$cite.fac<-as.factor(dat.graph4$citation)
fit.graph4<-rma(yi=Eff, vi=Eff_var, data=dat.graph4)

#Adding back in Outcome as a variable
dat.graph4$female<-round(dat.graph4$female, digits = 2)
dat.graph4$age<-round(dat.graph4$age, digits = 2)
par(mar=c(4,4,1,2))

tiff(paste0(graphics.folder,'/Prosociality.tiff'), res=1200, units = 'in', height = 7, width=10)
par(cex=1, font=1)
forest(fit.graph4, xlim=c(-10, 2), 
       order = order(dat.graph4$Eff), 
       ilab = cbind(dat.graph4$N, 
                    dat.graph4$female, 
                    dat.graph4$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:(nrow(dat.graph4)+2)),
       ylim = c(-1,nrow(dat.graph4)+5), 
       cex=.75, 
       xlab=expression(~rho), 
       mlab="", 
       addfit = F,
       slab = dat.graph4$cite.fac)

par(cex=1, font=4)
par(font=4)
text(-7, nrow(dat.graph4)+4, 'N')
par(font=2)
text(c(-5, -3), nrow(dat.graph4)+4, c('%Female', 'Mean Age'))
text(-10, nrow(dat.graph4)+4, 'Citation', pos=4)
text(1.25, nrow(dat.graph4)+4, expression(paste(bolditalic('r'), '[LB,UB]')))

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

tiff(paste0(graphics.folder,'/Guilt.tiff'), res=1200, units = 'in', height = 3.5, width=10)
par(cex=1, font=1)
forest(fit.graph5, xlim=c(-10, 2), 
       order = order(dat.graph5$Eff), 
       ilab = cbind(dat.graph5$N, 
                    dat.graph5$female, 
                    dat.graph5$age),
       ilab.xpos = c(-7, -5, -3), 
       rows = c(3:(nrow(dat.graph5)+2)),
       ylim = c(-1, nrow(dat.graph5)+5), 
       cex=.75, 
       xlab=expression(~rho), 
       mlab="", 
       addfit = F,
       slab = dat.graph5$cite.fac)

par(cex=1, font=4)
par(font=4)
text(-7, nrow(dat.graph5)+4, 'N')
par(font=2)
text(c(-5, -3), nrow(dat.graph5)+4, c('%Female', 'Mean Age'))
text(-10, nrow(dat.graph5)+4, 'Citation', pos=4)
text(1.25, nrow(dat.graph5)+4, expression(paste(bolditalic('r'), '[LB,UB]')))

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
summary(fit.Emp_cog.mod_age)  #So this is significant again (have to remake that figure)  

#Empathy Difference Model:
fit.emp_comp.mod_age<-rma(yi=Eff, 
                         vi=Eff_var,
                         weights = 1/Eff_var,
                         mods = ~age, 
                         data=dat.emp_comp, 
                         ni=N, 
                         knha=T, 
                         method = 'HS')
summary(fit.emp_comp.mod_age)  #So this is significant again (have to remake that figure)  



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

#Interesting that we have opposing directions here (worth exploring visually)
fit.emp_comp.mod_CU_resp<-rma(yi=Eff, 
                             vi=Eff_var,
                             weights = 1/Eff_var,
                             mods = ~CU_resp, 
                             data=dat.emp_comp, 
                             ni=N, 
                             knha=T, 
                             method = 'HS')
summary(fit.emp_comp.mod_CU_resp)  

#Empathy Difference Model Robust test of respondent/age:
fit.emp_comp.mod_age_CU_resp<-rma(yi=Eff, 
                                  vi=Eff_var,
                                  weights = 1/Eff_var,
                                  mods = ~age + CU_resp, 
                                  data=dat.emp_comp, 
                                  ni=N, 
                                  knha=T, 
                                  method = 'HS')
summary(fit.emp_comp.mod_age_CU_resp)  #So this is significant again (have to remake that figure)  

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
summary(fit.Emp_cog.mod_Out_resp) #Significant - need this plotted for sure  

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
#First plot is contrast of marginal interaction effects for CU respondent & association w/ affective vs. cognitive empathy
CU_resp<-c(0,1)
emp_aff.CU_resp.pred<-predict(fit.Emp_aff.mod_CU_resp, newmods = CU_resp, level = 95)

#Graphing moderation
#Exploring moderation graphically involving dichotomous ICU variable (Cognitive Empathy Model)

CU_self<-emp_aff.CU_resp.pred$pred[1]
CU_self.se<-emp_aff.CU_resp.pred$se[1]

CU_other<-emp_aff.CU_resp.pred$pred[2]
CU_other.se<-emp_aff.CU_resp.pred$se[2]

CU_self.dist<-rnorm(100000, mean=CU_self, sd=CU_self.se)
CU_other.dist<-rnorm(100000, mean=CU_other, sd=CU_other.se)

emp_aff.CU.pred.DF<-cbind(CU_self.dist, CU_other.dist)
colnames(emp_aff.CU.pred.DF)<-c('CU Respondent: Self', 'CU Respondent: Other')

#------------------------------------------------------------------------------
emp_cog.CU_resp.pred<-predict(fit.Emp_cog.mod_CU_resp, newmods = CU_resp, level = 95)

#Graphing moderation
#Exploring moderation graphically involving dichotomous ICU variable (Cognitive Empathy Model)

CU_self<-emp_cog.CU_resp.pred$pred[1]
CU_self.se<-emp_cog.CU_resp.pred$se[1]

CU_other<-emp_cog.CU_resp.pred$pred[2]
CU_other.se<-emp_cog.CU_resp.pred$se[2]

CU_self.dist<-rnorm(100000, mean=CU_self, sd=CU_self.se)
CU_other.dist<-rnorm(100000, mean=CU_other, sd=CU_other.se)

emp_cog.CU.pred.DF<-cbind(CU_self.dist, CU_other.dist)
colnames(emp_cog.CU.pred.DF)<-c('CU Respondent: Self', 'CU Respondent: Other')

#------------------------------------------------------------------------------
g.aff<-mcmc_areas(emp_aff.CU.pred.DF, prob = .95)+
  xlab(expression(~rho))+
  ggtitle('Association between Affective Empathy and CU Traits as a Function of CU Respondent')+
  theme_bw()+
  geom_vline(xintercept = 0, lty='dashed', color='red')+
  scale_x_continuous(limits = c(-1, 0.75))

g.cog<-mcmc_areas(emp_cog.CU.pred.DF, prob = .95)+
  xlab(expression(~rho))+
  ggtitle('Association between Cognitive Empathy and CU Traits as a Function of CU Respondent')+
  theme_bw()+
  geom_vline(xintercept = 0, lty='dashed', color='red')+
  scale_x_continuous(limits = c(-1, 0.75))
###############################################################################
Out_resp<-c(0,1)
emp_cog.Out_resp.pred<-predict(fit.Emp_cog.mod_Out_resp, newmods = Out_resp, level = 95)

#Graphing moderation
#Exploring moderation graphically involving dichotomous ICU variable (Cognitive Empathy Model)

Out_self<-emp_cog.Out_resp.pred$pred[1]
Out_self.se<-emp_cog.Out_resp.pred$se[1]

Out_other<-emp_cog.Out_resp.pred$pred[2]
Out_other.se<-emp_cog.Out_resp.pred$se[2]

Out_self.dist<-rnorm(100000, mean=Out_self, sd=Out_self.se)
Out_other.dist<-rnorm(100000, mean=Out_other, sd=Out_other.se)

emp_cog.Out.pred.DF<-cbind(Out_self.dist, Out_other.dist)
colnames(emp_cog.Out.pred.DF)<-c('Out Respondent: Self', 'Out Respondent: Other')

#------------------------------------------------------------------------------
g.cog2<-mcmc_areas(emp_cog.Out.pred.DF, prob = .95)+
  xlab(expression(~rho))+
  ggtitle('Association between Cognitive Empathy and CU Traits as a Function of Outcome Respondent')+
  theme_bw()+
  geom_vline(xintercept = 0, lty='dashed', color='red')+
  scale_x_continuous(limits = c(-1, 0.1))

tiff(paste0(graphics.folder, '/Cog_Moderated_by_Out_Resp.tiff'), res=300, units='in', height=8, width=11)
g.cog2
dev.off()

###############################################################################
CU_resp<-c(0,1)
emp_comp.CU_resp.pred<-predict(fit.emp_comp.mod_CU_resp, newmods = CU_resp, level = 95)

#Graphing moderation
#Exploring moderation graphically involving dichotomous ICU variable (Cognitive Empathy Model)

CU_self<-emp_comp.CU_resp.pred$pred[1]
CU_self.se<-emp_comp.CU_resp.pred$se[1]

CU_other<-emp_comp.CU_resp.pred$pred[2]
CU_other.se<-emp_comp.CU_resp.pred$se[2]

CU_self.dist<-rnorm(100000, mean=CU_self, sd=CU_self.se)
CU_other.dist<-rnorm(100000, mean=CU_other, sd=CU_other.se)

emp_comp.CU.pred.DF<-cbind(CU_self.dist, CU_other.dist)
colnames(emp_comp.CU.pred.DF)<-c('CU Respondent: Self', 'CU Respondent: Other')

#------------------------------------------------------------------------------
g.diff<-mcmc_areas(emp_comp.CU.pred.DF, prob = .95)+
  xlab(expression(~Delta~rho))+
  ggtitle('Association between Affective-Cognitive Empathy Difference and CU Traits as a Function of CU Respondent')+
  theme_bw()+
  geom_vline(xintercept = 0, lty='dashed', color='red')+
  scale_x_continuous(limits = c(-1, 0.75))

tiff(paste0(graphics.folder, '/Emp_Diff_Moderated_by_CU_Resp.tiff'), res=300, units='in', height=8, width=11)
g.diff
dev.off()

tiff(paste0(graphics.folder, '/AffvCog_Moderated_by_CU_Resp.tiff'), res=300, units='in', height=8, width=11)
cowplot::plot_grid(g.aff, g.cog, g.diff,
                   labels=c("A", "B", "C"), 
                   align = "v", 
                   nrow = 3)
dev.off()

###############################################################################
age<-seq(min(dat.Emp_cog$age), max(dat.Emp_cog$age), by = .01)

emp_cog.age.pred<-predict(fit.Emp_cog.mod_age, newmods = age, level = 95)
emp_cog.age.orig<-predict(fit.Emp_cog.mod_age, newmods = dat.Emp_cog$age, level = 95)

emp_cog.age.DF<-data.frame(pred_y = emp_cog.age.pred$pred, 
                           CI_LB = emp_cog.age.pred$ci.lb, 
                           CI_UB = emp_cog.age.pred$ci.ub, 
                           Age = age)
dat.Emp_cog$CI_age_LB<-emp_cog.age.orig$ci.lb-emp_cog.age.orig$pred
dat.Emp_cog$CI_age_UB<-emp_cog.age.orig$ci.ub-emp_cog.age.orig$pred


g.cog3<-ggplot()+
  geom_ribbon(data = emp_cog.age.DF, 
              aes(ymin = CI_LB, 
                  ymax = CI_UB, 
                  x = Age),
              fill = color_scheme_get(scheme = "blue")[[2]], 
              alpha = .75)+
  geom_line(data = emp_cog.age.DF, 
            aes(x = Age, 
                y = pred_y),
            color = color_scheme_get(scheme = "blue")[[6]], 
            lwd = 2, 
            inherit.aes = FALSE)+
  geom_point(data = dat.Emp_cog, 
             aes(x = age, 
                 y = Eff, 
                 size = N), 
             inherit.aes = FALSE, 
             alpha = .5)+
  geom_errorbar(data = dat.Emp_cog,
                aes(ymin = Eff + CI_age_LB, 
                    ymax = Eff + CI_age_UB, 
                    x = age), 
                inherit.aes = FALSE, 
                alpha = .5, 
                width = 0)+
  geom_hline(yintercept = 0, lty='dashed', color = 'red')+
  scale_size_continuous(expression(italic('N')), 
                        range = c(1, 10), 
                        breaks = seq(250, 1250, by = 250))+
  ylab(expression(~rho))+
  xlab('Age (years)')+
  ggtitle("Moderation of Association between CU and Cognitive Empathy by Sample Mean Age")+
  theme_bw()


tiff(paste0(graphics.folder, '/Cog_Moderated_by_Age.tiff'), res=300, units='in', height=8, width=11)
g.cog3
dev.off()

#############################################################################################################################
#############################################################################################################################
length(unique(dat$citation))

library(dplyr)
#Total unique N included
temp.df.N <- dat %>% 
  distinct(citation, N) 

sum(temp.df.N$N)

temp.df.age <- dat %>% 
  distinct(citation, age, N)

sum(temp.df.age$age*temp.df.age$N)/sum(temp.df.age$N)

temp.df.female <- dat %>% 
  distinct(citation, female, N)

sum(temp.df.female$female*temp.df.female$N)/sum(temp.df.female$N)/100

###############################################################################
#Putting together Comparison of corrected vs. uncorrected correlations...
dat.Emp_tot$citation<-as.character(dat.Emp_tot$citation)
table(dat.Emp_tot$citation)

dat.Emp_tot$citation[dat.Emp_tot$Id==9 & dat.Emp_tot$female==0]<-'Dadds et al. (2009m)'   #male sample
dat.Emp_tot$citation[dat.Emp_tot$Id==9 & dat.Emp_tot$female==100]<-'Dadds et al. (2009f)'   #female sample
dat.Emp_tot$citation[dat.Emp_tot$citation=="Raine and Chen (2018)"]<-"Raine & Chen (2018)"

corr_Emp_tot_R<-reshape2::melt(dat.Emp_tot[,c("citation", 
                                              "R", 
                                              "Eff")], 
                               id.var = "citation")

corr_Emp_tot_var<-reshape2::melt(dat.Emp_tot[,c("citation", 
                                              "R_var", 
                                              "Eff_var")], 
                               id.var = "citation")

corr_Emp_tot_R$citation==corr_Emp_tot_var$citation

corr_Emp_tot<-data.frame(citation = corr_Emp_tot_R$citation, 
                         Eff = corr_Emp_tot_R$value,
                         Type = ifelse(corr_Emp_tot_R$variable=="R", "Raw", "Corrected"), 
                         LB = corr_Emp_tot_R$value + qnorm(.025)*sqrt(corr_Emp_tot_var$value), 
                         UB = corr_Emp_tot_R$value + qnorm(.975)*sqrt(corr_Emp_tot_var$value))

######################################
#Need separate effects on different lines (treated distinct )
g.Emp_tot<-ggplot(data = corr_Emp_tot, 
                  aes(x = Eff, 
                      y = fct_rev(as_factor(citation)),
                      color = Type, 
                      shape = Type)
                  )+
  geom_point(size = 5)+
  geom_errorbarh(aes(xmin = LB, 
                     xmax = UB), 
                 alpha = .5, 
                 lwd = 1.5)+
geom_vline(aes(xintercept =0), color='red', lty='dashed')+
  ylab('Citation')+
  xlab(expression(italic('r')))+
  ggtitle('Comparing Attenuation-Corrected and Raw Correlations: Total Empathy & CU')+
  scale_x_continuous(limits = c(-1.2, .6))+
  labs(color = element_blank(), 
       shape = element_blank())+
  theme_bw()+
  theme(axis.text.y = element_text(size = 14), 
        legend.text = element_text(size = 18))

g.Emp_tot

######################################
dat.Emp_aff$citation<-as.character(dat.Emp_aff$citation)
table(dat.Emp_aff$citation)

dat.Emp_aff$citation[dat.Emp_aff$Id==5 & dat.Emp_aff$female==0]<-'Brouns et al. (2013m)'   #male sample
dat.Emp_aff$citation[dat.Emp_aff$Id==5 & dat.Emp_aff$female==100]<-'Brouns et al. (2013f)'   #female sample
dat.Emp_aff$citation[dat.Emp_aff$Id==9 & dat.Emp_aff$female==0]<-'Dadds et al. (2009m)'   #male sample
dat.Emp_aff$citation[dat.Emp_aff$Id==9 & dat.Emp_aff$female==100]<-'Dadds et al. (2009f)'   #female sample

dat.Emp_aff$citation[dat.Emp_aff$citation=="Raine and Chen (2018)"]<-"Raine & Chen (2018)"

corr_Emp_aff_R<-reshape2::melt(dat.Emp_aff[,c("citation", 
                                              "R", 
                                              "Eff")], 
                               id.var = "citation")

corr_Emp_aff_var<-reshape2::melt(dat.Emp_aff[,c("citation", 
                                                "R_var", 
                                                "Eff_var")], 
                                 id.var = "citation")

corr_Emp_aff_R$citation==corr_Emp_aff_var$citation

corr_Emp_aff<-data.frame(citation = corr_Emp_aff_R$citation, 
                         Eff = corr_Emp_aff_R$value,
                         Type = ifelse(corr_Emp_aff_R$variable=="R", "Raw", "Corrected"), 
                         LB = corr_Emp_aff_R$value + qnorm(.025)*sqrt(corr_Emp_aff_var$value), 
                         UB = corr_Emp_aff_R$value + qnorm(.975)*sqrt(corr_Emp_aff_var$value))

######################################
#Need separate effects on different lines (treated distinct )
g.Emp_aff<-ggplot(data = corr_Emp_aff, 
                  aes(x = Eff, 
                      y = fct_rev(as_factor(citation)),
                      color = Type, 
                      shape = Type)
)+
  geom_point(size = 5)+
  geom_errorbarh(aes(xmin = LB, 
                     xmax = UB), 
                 alpha = .5, 
                 lwd = 1.5)+
  geom_vline(aes(xintercept =0), color='red', lty='dashed')+
  ylab('Citation')+
  xlab(expression(italic('r')))+
  ggtitle('Comparing Attenuation-Corrected and Raw Correlations: Affective Empathy & CU')+
  scale_x_continuous(limits = c(-1.2, .6))+
  labs(color = element_blank(), 
       shape = element_blank())+
  theme_bw()+
  theme(axis.text.y = element_text(size = 14), 
        legend.text = element_text(size = 18))

g.Emp_aff

######################################
dat.Emp_cog$citation<-as.character(dat.Emp_cog$citation)
table(dat.Emp_cog$citation)

dat.Emp_cog$citation[dat.Emp_cog$Id==5 & dat.Emp_cog$female==0]<-'Brouns et al. (2013m)'   #male sample
dat.Emp_cog$citation[dat.Emp_cog$Id==5 & dat.Emp_cog$female==100]<-'Brouns et al. (2013f)'   #female sample
dat.Emp_cog$citation[dat.Emp_cog$Id==9 & dat.Emp_cog$female==0]<-'Dadds et al. (2009m)'   #male sample
dat.Emp_cog$citation[dat.Emp_cog$Id==9 & dat.Emp_cog$female==100]<-'Dadds et al. (2009f)'   #female sample

dat.Emp_cog$citation[dat.Emp_cog$citation=="Raine and Chen (2018)"]<-"Raine & Chen (2018)"

corr_Emp_cog_R<-reshape2::melt(dat.Emp_cog[,c("citation", 
                                              "R", 
                                              "Eff")], 
                               id.var = "citation")

corr_Emp_cog_var<-reshape2::melt(dat.Emp_cog[,c("citation", 
                                                "R_var", 
                                                "Eff_var")], 
                                 id.var = "citation")

corr_Emp_cog_R$citation==corr_Emp_cog_var$citation

corr_Emp_cog<-data.frame(citation = corr_Emp_cog_R$citation, 
                         Eff = corr_Emp_cog_R$value,
                         Type = ifelse(corr_Emp_cog_R$variable=="R", "Raw", "Corrected"), 
                         LB = corr_Emp_cog_R$value + qnorm(.025)*sqrt(corr_Emp_cog_var$value), 
                         UB = corr_Emp_cog_R$value + qnorm(.975)*sqrt(corr_Emp_cog_var$value))

######################################
#Need separate effects on different lines (treated distinct )
g.Emp_cog<-ggplot(data = corr_Emp_cog, 
                  aes(x = Eff, 
                      y = fct_rev(as_factor(citation)),
                      color = Type, 
                      shape = Type)
)+
  geom_point(size = 5)+
  geom_errorbarh(aes(xmin = LB, 
                     xmax = UB), 
                 alpha = .5, 
                 lwd = 1.5)+
  geom_vline(aes(xintercept =0), color='red', lty='dashed')+
  ylab('Citation')+
  xlab(expression(italic('r')))+
  ggtitle('Comparing Attenuation-Corrected and Raw Correlations: Cognitive Empathy & CU')+
  scale_x_continuous(limits = c(-1.2, .6))+
  labs(color = element_blank(), 
       shape = element_blank())+
  theme_bw()+
  theme(axis.text.y = element_text(size = 14), 
        legend.text = element_text(size = 18))

g.Emp_cog

######################################
dat.prosoc$citation<-as.character(dat.prosoc$citation)
table(dat.prosoc$citation)

dat.prosoc$citation[dat.prosoc$Id==42 & dat.prosoc$female==0]<-'Pechorro et al. (2013m)'   #male sample
dat.prosoc$citation[dat.prosoc$Id==42 & dat.prosoc$female==100]<-'Pechorro et al. (2013f)'   #female sample
dat.prosoc$citation[dat.prosoc$Id==43 & dat.prosoc$female==0]<-'Pechorro et al. (2015am)'   #male sample
dat.prosoc$citation[dat.prosoc$Id==43 & dat.prosoc$female==100]<-'Pechorro et al. (2015af)'   #female sample

corr_prosoc_R<-reshape2::melt(dat.prosoc[,c("citation", 
                                              "R", 
                                              "Eff")], 
                               id.var = "citation")

corr_prosoc_var<-reshape2::melt(dat.prosoc[,c("citation", 
                                                "R_var", 
                                                "Eff_var")], 
                                 id.var = "citation")

corr_prosoc_R$citation==corr_prosoc_var$citation

corr_prosoc<-data.frame(citation = corr_prosoc_R$citation, 
                         Eff = corr_prosoc_R$value,
                         Type = ifelse(corr_prosoc_R$variable=="R", "Raw", "Corrected"), 
                         LB = corr_prosoc_R$value + qnorm(.025)*sqrt(corr_prosoc_var$value), 
                         UB = corr_prosoc_R$value + qnorm(.975)*sqrt(corr_prosoc_var$value))

######################################
#Need separate effects on different lines (treated distinct )
g.prosoc<-ggplot(data = corr_prosoc, 
                  aes(x = Eff, 
                      y = fct_rev(as_factor(citation)),
                      color = Type, 
                      shape = Type)
)+
  geom_point(size = 5)+
  geom_errorbarh(aes(xmin = LB, 
                     xmax = UB), 
                 alpha = .5, 
                 lwd = 1.5)+
  geom_vline(aes(xintercept =0), color='red', lty='dashed')+
  ylab('Citation')+
  xlab(expression(italic('r')))+
  ggtitle('Comparing Attenuation-Corrected and Raw Correlations: Prosociality & CU')+
  scale_x_continuous(limits = c(-1.2, .6))+
  labs(color = element_blank(), 
       shape = element_blank())+
  theme_bw()+
  theme(axis.text.y = element_text(size = 14), 
        legend.text = element_text(size = 18))

g.prosoc

######################################
dat.glt$citation<-as.character(dat.glt$citation)
table(dat.glt$citation)

corr_glt_R<-reshape2::melt(dat.glt[,c("citation", 
                                              "R", 
                                              "Eff")], 
                               id.var = "citation")

corr_glt_var<-reshape2::melt(dat.glt[,c("citation", 
                                                "R_var", 
                                                "Eff_var")], 
                                 id.var = "citation")

corr_glt_R$citation==corr_glt_var$citation

corr_glt<-data.frame(citation = corr_glt_R$citation, 
                         Eff = corr_glt_R$value,
                         Type = ifelse(corr_glt_R$variable=="R", "Raw", "Corrected"), 
                         LB = corr_glt_R$value + qnorm(.025)*sqrt(corr_glt_var$value), 
                         UB = corr_glt_R$value + qnorm(.975)*sqrt(corr_glt_var$value))

######################################
#Need separate effects on different lines (treated distinct )
g.glt<-ggplot(data = corr_glt, 
                  aes(x = Eff, 
                      y = fct_rev(as_factor(citation)),
                      color = Type, 
                      shape = Type)
)+
  geom_point(size = 5)+
  geom_errorbarh(aes(xmin = LB, 
                     xmax = UB), 
                 alpha = .5, 
                 lwd = 1.5)+
  geom_vline(aes(xintercept =0), color='red', lty='dashed')+
  ylab('Citation')+
  xlab(expression(italic('r')))+
  ggtitle('Comparing Attenuation-Corrected and Raw Correlations: Guilt & CU')+
  scale_x_continuous(limits = c(-1.2, .6))+
  labs(color = element_blank(), 
       shape = element_blank())+
  theme_bw()+
  theme(axis.text.y = element_text(size = 14), 
        legend.text = element_text(size = 18))

g.glt


tiff(paste0(graphics.folder,'/Attenuation Correction vs Raw Effects.tiff'), res=300, units = 'in', height = 24, width=16)
cowplot::plot_grid(g.Emp_tot, g.Emp_aff, g.Emp_cog, g.prosoc, g.glt, 
                   labels = c('A', 'B', 'C', 'D', 'E'), 
                   nrow = 5, 
                   rel_heights = c(5, 3, 3, 3, 1),
                   align = "v")
dev.off()

#############################################################################################################################
library(metaSEM)
library(RColorBrewer)
#############################################################################################################################

#IMPORTANT - Section is exploratory and not updated with current report

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
