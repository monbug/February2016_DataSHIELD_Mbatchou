################################################################
#
#             DATASHIELD DEMO MONBUG: FEDERATED ANALYSIS
#
###############################################################
options(warn = -1) #remove annoying warnings

# Install all DataSHIELD packages
#install.packages('datashieldclient', repos=c(getOption('repos'), 'http://cran.datashield.org'), dependencies=TRUE)

# Source the utils environnment
source('utils_analysis.R',echo=F,print.eval=F)

#load libraries
library(datashieldclient)


#load login credentials
load('monbug-login.rda')

#login to datashield and assign data to 'D' as default
opals <- datashield.login(logins=logindata,assign=TRUE)
#take a look at the opals servers
print(opals)

#look at the objects in the sever(s) side
datashield.symbols(opals)

#check the variables of dataframe D in servers side
ds.colnames('D')

#compute the mean and sd of HDL
run.desc.stats('LAB_HDL',data='D')

#compute the stats of BMI_CATEGORIES
run.desc.stats('PM_BMI_CATEGORICAL',data='D')

#compute chi-square BMI_CAT x GENDER
run.table2d('GENDER','PM_BMI_CATEGORICAL',data='D')

###-------DIABETES
run.desc.stats('DIS_DIAB',data='D')
##DISCLOSIVE CASE when a cell count is < 5 particpants
run.table2d('DIAB','PM_BMI_CATEGORICAL',ref= 'study1',data='D')


###-------DIS_CVA (cardio vascular disease). one study is empty for this variable
run.desc.stats('DIS_CVA',data='D')
#remove the empty study in your analysis for DIS_CVA : study2 is empty
run.desc.stats('DIS_CVA',data='D',datasources = opals[-2])

######------- GLM DIAB ~ GENDER + LAB_GLUC_ADJUSTED + LAB_HDL
run.dummy.study('D') #create dummy study variable to account for heterogeneity in your statiscal model
formula <-  'D$DIS_DIAB~D$GENDER + D$LAB_GLUC_ADJUSTED + D$PM_BMI_CATEGORICAL + D$LAB_HDL'

# run the glm adding heterogeneity effect from each study using study1 as reference
glm.res <- run.meta.glm(formula=formula,family='binomial',ref = 'study1',check = T)

#glm.res <- ds.glm(formula,family='binomial') #without heterogeneity( this model does not fit the reality because studies are all different)

#extract glm statistics
run.extract.glm.stats(glm.res)

