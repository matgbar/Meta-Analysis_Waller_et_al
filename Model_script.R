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
dat<-read.csv(paste0(data.folder, 'Meta_raw.csv'), stringsAsFactors = F)
colnames(dat)[1]<-'id'
dat.Emp_tot<-dat[dat$Outcome=='empathy_tot',]
dat.Emp_aff<-dat[dat$Outcome=='empathy_aff',]
dat.Emp_cog<-dat[dat$Outcome=='empathy_cog',]
dat.glt<-dat[dat$Outcome=='guilt',]   #note probably have too few for guilt
dat.prosoc<-dat[dat$Outcome=='prosocial',]

#################################################################################
#Model for Total Empathy: 
fit.emp_tot<-rma(yi=Eff, vi=Eff_var, data=dat.Emp_tot)
summary(fit.emp_tot)

sink(paste0(model.folder, 'Total Empathy Overall Model - no moderators.txt'))
summary(fit.emp_tot)
sink()

jpeg(paste0(graphics.folder, 'CU and Total Empathy - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_tot)
title("Relation between CU Traits and Total Empathy (Fisher's z)")
dev.off()

jpeg(paste0(graphics.folder, 'CU and Total Empathy - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.emp_tot)
title("Funnel Plot of CU Relation with Total Empathy (Fisher's z)")
dev.off()

#################################################################################
#Model for Affective Empathy: 
fit.emp_aff<-rma(yi=Eff, vi=Eff_var, data=dat.Emp_aff)
summary(fit.emp_aff)

sink(paste0(model.folder, 'Affective Empathy Overall Model - no moderators.txt'))
summary(fit.emp_aff)
sink()

jpeg(paste0(graphics.folder, 'CU and Affective Empathy - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_aff)
title("Relation between CU Traits and Affective Empathy (Fisher's z)")
dev.off()

jpeg(paste0(graphics.folder, 'CU and Affective Empathy - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.emp_aff)
title("Funnel Plot of CU Relation with Affective Empathy (Fisher's z)")
dev.off()

#################################################################################
#Model for Cognitive Empathy: 
fit.emp_cog<-rma(yi=Eff, vi=Eff_var, data=dat.Emp_cog)
summary(fit.emp_cog)

sink(paste0(model.folder, 'Cognitive Empathy Overall Model - no moderators.txt'))
summary(fit.emp_cog)
sink()

jpeg(paste0(graphics.folder, 'CU and Cognitive Empathy - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.emp_cog)
title("Relation between CU Traits and Cognitive Empathy (Fisher's z)")
dev.off()

jpeg(paste0(graphics.folder, 'CU and Cognitive Empathy - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.emp_cog)
title("Funnel Plot of CU Relation with Cognitive Empathy (Fisher's z)")
dev.off()

#################################################################################
#Model for Prosociality: 
fit.prosoc<-rma(yi=Eff, vi=Eff_var, data=dat.prosoc)
summary(fit.prosoc)

sink(paste0(model.folder, 'Prosociality Overall Model - no moderators.txt'))
summary(fit.prosoc)
sink()

jpeg(paste0(graphics.folder, 'CU and Prosociality - Forest.jpeg'), res=300, width = 7, height=7, units='in')
forest.rma(fit.prosoc)
title("Relation between CU Traits and Prosociality (Fisher's z)")
dev.off()

jpeg(paste0(graphics.folder, 'CU and Prosociality - Funnel.jpeg'), res=300, width = 7, height=7, units='in')
funnel(fit.prosoc)
title("Funnel Plot of CU Relation with Prosociality (Fisher's z)")
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
                  rtoz_summary(dat.prosoc, .95, fit.prosoc, 'Prosociality')
                  )
stargazer(DF.summary, summary = F, out=paste0(model.folder,'Model Summary.txt'), rownames = F, 
          notes = 'LB and UB based on 95% CI')

#Plotting summary effect sizes 
g1<-ggplot(aes(x=rho, y=Outcome), data=DF.summary)
g2<-g1+geom_errorbarh(aes(xmin=rho.LB, xmax=rho.UB), height=.25)
g3<-g2+geom_point(aes(size=N), pch=18)
g4<-g3+scale_size_continuous(range = c(8,12))
g5<-g4+xlab(expression(rho))+ylab('')+
  ggtitle(paste0('Correlation between Outcome Measures and CU Traits'))+
  theme(plot.title = element_text(hjust=.5))
g5

jpeg(paste0(graphics.folder, 'Model Summary Graphic.jpeg'), res=300, units='in', height=7, width=7)
g5
dev.off()
