################################################################
#
#         POOLED ANALYSIS DEMO MONBUG: ILMA
#
###############################################################
# This script will run a glm using the physically pooled data.
# compare the result with the one obtained using the federated approach via DataSHIELD

# Source the utils environnment
source('utils_analysis.R',echo=F,print.eval=F)

library(opal) #will be installed as a dependency the datashield R packages
library(data.table) # install the package first

#connect to opal demo
url<-"https://opal-demo.obiba.org:8443"
o<-opal.login('administrator', 'password', url)

##get the tables from opal server to your client R [because you act as an administrator you have access to all the data] the pool everything 
study1 <- var.assign(opal = o,table = 'CNSIM1',datasource = 'ds')
idx1 <- row.names(study1)
study2 <- var.assign(opal = o,table = 'CNSIM2',datasource = 'ds')
idx2 <- row.names(study2)
study3 <- var.assign(opal = o,table = 'CNSIM3',datasource = 'ds')
idx3 <- row.names(study3)
#-------- Pool all studies (ILMA)
study.pool <- rbind(study1,study2,study3)
study.pool <- data.table(study.pool,keep.rownames = T)
setkey(study.pool,rn)

#--create s dummy variable to account for heterogeneity of studies
study.pool[idx1,study:=1,nomatch = NA]
study.pool[idx2,study:=2,nomatch = NA]
study.pool[idx3,study:=3,nomatch = NA]
#put s as a factorial variable
study.pool[,study:=as.factor(study)]

####-------GLM POOLED via ilma: DIS_DIAB ~ GENDER + LAB_GLUC_ADJUSTED + LAB_HDL
#with studies heterogeneity effects using study1 as reference
res <- study.pool[,glm(DIS_DIAB~study+GENDER+LAB_GLUC_ADJUSTED+PM_BMI_CATEGORICAL+LAB_HDL,family = 'binomial')] 

#res <- study.pool[,glm(DIS_DIAB~GENDER+LAB_GLUC_ADJUSTED+PM_BMI_CATEGORICAL+LAB_HDL,family = 'binomial')] #without heterogeneity effects (does not fit the reality)

#see coefficient AND compare with the result from the federated DataSHIELD version: Results might be slightly differents because on decimal rounding  
g.extract(res)


