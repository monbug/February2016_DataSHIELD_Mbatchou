bioshare.env <- new.env()

############################################################################
# utils functions
###########################################################################
# check if an object is defined
bioshare.env$isDefined <- dsBaseClient:::isDefined
#get opals login object(s) in the env
bioshare.env$findLoginObjects <- dsBaseClient:::findLoginObjects
#check the class of an object
bioshare.env$checkClass <- dsBaseClient:::checkClass

#check if an object is assigned
bioshare.env$isAssigned <- dsBaseClient:::isAssigned
#extract elements and object form server side vector (e.g. D$AGE_YRS)
bioshare.env$extract <- dsBaseClient:::extract
#pooled mean
bioshare.env$run.pooled.mean<- dsBaseClient:::getPooledMean

#vectorize list like structure
bioshare.env$.vectorize <- function(x,subscript,simplify = T) 
{ 
  if(simplify) sapply(x,'[[',subscript)
  else lapply(x,'[[',subscript)
}

#function get vars from opal
bioshare.env$var.assign<-function(opal,datasource,table,variables=NULL)
{
  datafile<-paste0(datasource, ".", table)
  if(length(variables)==1){
    variable<-paste0(datafile,':',variables)
    data.err <- try(opal.assign(opal,'D', variable,missings=T),silent = T)
  }else{
    data.err <- try(opal.assign(opal,'D', datafile,variables,missings=T),silent = T)
  }
  
  if(inherits(data.err,'try-error')) { 
    message(paste0('-- The required data is not assigned\n'),(data.err),'-- Check the error message and try again!')
  }else { cat(paste0('The required data from ', datafile,' was correctly assigned.'))}
  
  var.value<-opal.execute(opal,'D')
  return (var.value)
}

#extract glm result from R not datashield
bioshare.env$g.extract <- function(g)
{
  f <- g$family
  if(f$family == 'gaussian') g.coef <- coef(summary(g))[,1]
  else if (grepl('poisson',f$family)) g.coef <- exp(coef(summary(g)))[,1]
  else if (grepl('binomial',f$family)) {
    g.coef <- coef(summary(g))
    g.coef[-1,1] <- exp(g.coef[-1,1])
    g.coef[1,1]<- f$linkinv(g.coef[1,1])
  }
  g.pvalue <- coef(summary(g))[,4]
  pval<-format(g.pvalue,digits = 4); pval <- sapply(pval,function(x){if(as.numeric(x)<2.2e-16) {x <-'<2.2e-16'}; x});
  est <-format(g.coef[,1],digits=3) 
  g.coef[,1]<- paste0(est,'(',pval,')'); 
  d.f(effect_and_pvalue = g.coef[,1])
}


#########################################################################################
#    analytics functions here below 
########################################################################################

################ categorical (table2d) for multiple variables ######################
bioshare.env$run.cat<-function(subset,vars.list,type = NULL,save=F, print= F)
{
  if(missing(subset)) stop('subset variable is required ...',call.=F)
  if(missing(vars.list)) stop('var(s) is(are) required ...',call.=F)
  type <- match.arg(type,c('combine','split'))
  vars.list <-as.list(vars.list)
    
  #starting message
  message(paste0('--',toupper(type),' analysis'))
  message('--Categorical variables list: \n',paste0(vars.list,collapse='; '))
  
  ###preparing progress bar
  total <- length(vars.list)
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  i<-1
  
  ####start looping now
  result<-NULL
  for(var in vars.list){
    #compute table2D
    message(paste0('=> Computing chi-square ',var,' X ',subset,'\nDo not interrupt!...'))
    chi_square<-ds.table2D(x=paste0('D$',var), y=paste0('D$',subset), type=type,warningMessage=F)   ##### <----- TO CHANGE use my functions
    
    #arranging final result
    w<-structure(list(chi_square),.Names=var)
    result<-c(result,w)
    setTxtProgressBar(pb, i)
    i<-i+1
  }
  #close progress bar
  close(pb)
  
  if(as.logical(save)){
    ##Saving Result object in file versioned
    date<-format(Sys.Date(),'%d%b%y')
    resultname<-paste0('Cat_split_',subset,'_',date,'.rda')
    message(paste0('\n\n***\tSaving Results for Categorical Variables in <',resultname,'> file.'))
    save(result,file=resultname)
  }
  
  #print result
  if(as.logical(print)){
    return(result)
  }else{
    return(invisible(result))
  }
}


########INTERNAL FUNCTION
bioshare.env$run.get.subclass<-function(subclass = NULL,vars.list=NULL,data = NULL, datasources=NULL){
  
  #subset by sub_by first to generate the first subset
  if(is.null(datasources)) datasources <- findLoginObjects()
  ds <- datasources
  
  #verify data
  if(is.null(data)) {stop('Please specify the data(dataframe) to subset',call.=F)}
  
  if(is.null(subclass)) {
    stop('Nothing to do, sublass var is missing required...',call.=F)
  }else {
    #Sanity check : subclass should always be a factor
    subset.class <- checkClass(ds,paste0(data,'$',subclass))
    if(subset.class != 'factor') stop('subclass must be a categorical or factor type variable ...',call.=F)
    
    if (!is.null(vars.list)){
      vars<-unique(c(subclass,unlist(vars.list)))
      message(paste0('=> Creating sub classes of  selected vars.list in ',data, ' by ',subclass,'\nWait please do not interrupt!...'))
      
      cally <- call('subsetDS',dt = data, complt=F, rs=NULL, cs=vars)
      datashield.assign(ds,'newdt',cally)   
      #sublclassing 
      cally <- call('subsetByClassDS','newdt', subclass)
      datashield.assign(ds,'subobj',cally)
      to.rm <- c('newdt','subobj')
      
    } else { #sublcass only
      message(paste0('=> Creating sub classes of all vars in ',data, ' by ',subclass,'\nWait please do not interrupt!...'))
  
      cally <- call('subsetByClassDS',data, subclass)
      datashield.assign(ds,'subobj',cally)
      to.rm <- 'subobj'
    }
  }
  
  #### arranging subclass result in server side
  
  #define infoname
  message('=> Assigning the subsetted dataframe(s) on server side...\nWait please do not interrupt!...')
  cally <- as.name('namesDS(subobj)')
  subinfo<-datashield.aggregate(ds,cally)[[1]] #get names of subsetted object
    
  #define name of object to be assigned
  subinfobj<-sapply(subinfo,function(x){  #assign new subsetted object in server and return their name
    newname<-paste0(data,'.',sub('\\.level_(\\d+)$','\\1',x))  #transform names
    toassign<-paste0('subobj$',x)
    datashield.assign(ds,newname,as.name(toassign))    ###assigned objs are in server
    return(newname)
  },USE.NAMES = F)
  
  
  
  to.rm <- c("complt","cs","dt","rs",to.rm)
  for(x in to.rm) {datashield.rm(ds,x)} #Clean space
  
  message(paste0('-- dataframe ',subinfobj ,' is assigned',collapse='\n'))
  info <- paste0("ds.colnames('",subinfobj,"')",collapse=' OR ')
  cat(paste0("You may check the assigned subsetted dataframes with the following datashield command: ",info,""))
  return(invisible(subinfobj))
}

####################################################################################################
# this function recode the levels of a factor variable ex: bin a bmi 1,2,3 to bmi 0, 1 
bioshare.env$run.changelevel <- function(var,new.levels,new.varname=NULL,data=NULL,datasources) {
  
  if(missing(var)) stop('var is required')
  if(missing(new.levels)) stop('new.levels is required')
  if(missing(datasources)) datasources <- findLoginObjects()
  
  ds <- datasources
  
  # if no name was provided for the new variable give it a default name
  if(is.null(new.varname)) { 
    new.varname <- paste0(var,'.new')
    warning(paste0('new.varname is not specified ... the new var will be assigned to ',new.varname),immediate.=T,call.=F)
  }
  new.var <- new.varname
  
  if(!is.null(data)) var.dollar <- paste0(data,'$',var)
  else var.dollar <- var 
  
  # call the internal function that checks the input object is of the same class in all studies.
  typ <- checkClass(ds, var.dollar)
  
  # if input vector is not a factor stop
  if(typ != 'factor') stop("The input var must be a factor!", call.=FALSE)
  
  # get the current number of levels
  cally <- paste0("levels(", var.dollar, ")")
  xx <- datashield.aggregate(ds[1], as.symbol(cally))  ###check info only for one server reduce time
  if((ll<-length(unique(unlist(xx)))) != length(new.levels)){
    stop(gettextf("Please new.levels must be a vector of size %d",ll), call.=FALSE)
  }
  message(paste0('->recoding levels of ',var,' ...'))
  new.levels.char <- paste0('c(',paste0(new.levels,collapse=','),')') #ex:  "c(0,0,1)"
  cally <- paste0('recodeLevelsDS(',var.dollar,',',new.levels.char,')')
  datashield.assign(ds,new.var,as.name(cally))
  
  if(!is.null(data)) {
    # put in the dataframe data
    callycbind <- paste0('cbind(',data,',',new.var,')')
    datashield.assign(ds,data,as.name(callycbind))
    #clean server side
    datashield.rm(ds,new.var)   
    cat(paste0("You may check the level-recoded variable with the following commands: run.desc.stats('",new.var,"',data='",data,"')"))  
  }else{
    cat(paste0("You may check the level-recoded variable with the following commands: run.desc.stats('",new.var,"')"))
  }
  
}



#########################################################################################################
########################################### MODELING UTILS ##############################################


#######################################################################
# this function remove the sugroup term in original model: useful for subgroup analysis

bioshare.env$run.adjust.model <- function(model,term){
  x <- unlist(strsplit(model,'\\+'))
  x.ind <- which(sapply(x,function(k) grepl(term,k,ignore.case=T)))
  if(length(x.ind)) paste0(x[-x.ind],collapse='+') 
  else paste0(x,collapse='+')
} 

#################### create dummy study effect vars #####

bioshare.env$run.dummy.study <- function (data,datasources)
{
  if(missing(data)) stop('data is mandatory ...\nplease add data name as character!')
  if(missing(datasources)) datasources <- findLoginObjects()
  ds <- datasources
  # call length
  cally <- paste0('NROW(',data,')')
  ln <- datashield.aggregate(ds,as.symbol(cally))
  #make zeros vector in each server 
  for (i in 1:length(ds)) { datashield.assign(ds[i],'zeros.dummy',call('rep',0,ln[[i]])) }
  
  message(paste0(paste0(names(ds),collapse=', '), ' dummy variable(s) will be stored in dataframe ',data))
  
  for (i in 1:length(ds)){
    
    effect_name<-names(ds[i])
    #assign 1 to study and 0 to others
    message(paste0('---processing ',effect_name,' dummy ...'))
    datashield.assign(ds[i],effect_name,as.name('zeros.dummy+1'))
    datashield.assign(ds[-i],effect_name,as.name('zeros.dummy'))
    
    callcbind<-paste0('cbind(',data,',',effect_name,')')
    datashield.assign(ds[i],data,as.symbol(callcbind))
    datashield.assign(ds[-i],data,as.symbol(callcbind))
    
    datashield.rm(ds[i],effect_name)
    datashield.rm(ds[-i],effect_name)
    #cat(capture.output(datashield.aggregate(ds,as.name(paste0('meanDS(',effect_name,')'))))) 
  }
  cat(paste0("You may check the stored dummy variables with the following datashield command: ds.colnames('",data,"')"))
}


####################################   GLM   #######################
bioshare.env$run.meta.glm<-function(formula, family, ref, datasources,save = F, print = T,...)
{
  if(missing(formula)){
    stop('formula is required ...',call.=F)
  }
  if(missing(family)){
    stop('family is required ...',call.=F)
  }
  if(missing(datasources)){
    datasources<-findLoginObjects()
  }
  ds <- datasources
  
  if(length(ds) > 1){
    if(missing(ref) || !ref %in% names(ds)) stop('ref study is required when running a glm with more than one studies!\nPlease check that reference study is correctly spelled...',call.=F)
    
    formulasplit<-unlist(strsplit(formula,'~'))
    outcomevar<-formulasplit[1]
    explanvars<-formulasplit[2]
    
    data <- extract(outcomevar)$holders
    effect_name <- paste0(data,'$',names(ds))
    effect_name <- effect_name[-(which(grepl(ref,effect_name)))]
    effect.vars.in.formula <- paste(effect_name,collapse='+')
    
    if(!is.na(explanvars)) explanvars <- paste(effect.vars.in.formula,explanvars,sep = '+')
    else explanvars<-effect.vars.in.formula
    
    #update formula with study effect dummies vars
    formula<-paste0(outcomevar,'~',explanvars)  
  }
  
  message('->running glm...\n wait do not interrupt!\n')
  message(paste0('formula for glm: ',formula))  
  message(paste0('family for glm: ',family))
  
  #run glm now
  glm.result<-ds.glm(formula=formula,family=family,datasources=ds,...)
  
  
  ##Saving Result object in file versioned
  if(as.logical(save)){
    
    date<-format(Sys.Date(),'%d%b%y')
    save_file<-paste0('glm_result','_',date,'.rda')
    message(paste0('saving result in ',save_file,'...'))
    save(glm.result,file=save_file)
  }
  
  #print result
  if(as.logical(print)){
    return(glm.result)
  }else{
    return(invisible(glm.result))
  }
  
}


####################################################################################
#this function create a formula according to the model, outcome and exposition vars

bioshare.env$run.update.formula<-function(outcome,expo,model,data)
{
  mf <- match.call(expand.dots = FALSE)
  arg.call <- names(mf)[-1]
  arg.names<-c('outcome','expo','model','data')
  missing.arg<-which(!arg.names %in% arg.call)
  missing.call <- paste(arg.names[missing.arg],collapse=' and ')
  if(length(missing.arg)>1) stop(paste0(missing.call,' are required'),call.=F)
  else if (length(missing.arg)==1) stop(paste0(missing.call,' is required'),call.=F)
  
  fm <-  paste0(data,'$',outcome,'~',data,'$',model,'+',expo)
  fm <- gsub('+',paste0('+',data,'$'),fm,fixed=T)
  fm <- gsub('*',paste0('*',data,'$'),fm,fixed=T)  
  fm
}



################################################################################################
#this function run a glm based on the model
#param outcome = the outcome variable in character (ex: 'SYM_SBREATH')
#param expo =  the exposition variable in character (ex: 'NO2_ESCAPE')
#param model = the model to apply in glm, it's a truncated formula made of -->
#-->counfounding variables (ex: '~D$AGE_YRS+D$EDU_HIGHEST_2+D$SMK_STATUS') without outcome and expo
#param family = glm family (ex: 'gaussian' or 'binomial',...) 
#param ...= any params that go to glm (except formula) (i.e: offset, data, weights,wiewIter, datasource)
#param Ncases = a boolean (TRUE or FALSE(default)) wether to compute N cases or not in the final result  

bioshare.env$run.model<-function(outcome,expo,model,family,data,Ncases=FALSE,pval=FALSE, ...)
{
  if(missing(outcome)) stop('outcome is required...',call.=F)
  if(missing(expo)) stop('exposition variable is required...',call.=F)
  if(missing(model)) stop('model formula is required...',call.=F)
  #verify data
  if(missing(data)) stop('data is mandatory',call.=F)
  
  #update formula
  formula <- run.update.formula(outcome,expo,model,data)
  
  #run glm 
  glm.res <- try(run.meta.glm(formula,family,print=T,...),silent=T)
  glm.err <- inherits(glm.res,'try-error')
  glm.ok <-  !( glm.err || is.null(glm.res)) 
  
  #extract glm stats and process result
  if(glm.ok) glm.stats<-run.extract.glm.stats(glm.res,Ncases=Ncases,pval=pval)
  else if (glm.err) return(message(glm.res))
  else return(glm.res)
  
  glm.stats$results <- glm.stats$stats[expo,,F]   #<--display variable names and statistics
  
  
  #print glm stats
  cat(glm.stats$formula,'\n\n')
  print(glm.stats$results)
  return(invisible(glm.stats$results))
}


##############################################
#function to extract glm result :P_OR(p.value) 
bioshare.env$run.extract.glm.stats <- function(glm.result,pval=FALSE,Ncases=FALSE,rdigit =3)
{
  if(missing(glm.result)) stop('Please provide a valid glm result...',call.=F)
  glm.family <- glm.result$family$family
  glm.coef <- coef(glm.result)
  if(grepl("poisson|binomial", glm.family)){
    stats <- data.frame(OR_CI = apply(glm.coef,1,function(x) {
        OR <- round(x['P_OR'],rdigit)
        pvalue<- format(x['p-value'],digits=4) ; pvalue <- if(as.numeric(pvalue)< 2.2e-16) {"<2.2e-16"}else{pvalue}
        low <- round(x[length(x)-1],rdigit)
        high <- round(x[length(x)],rdigit)
        
        res.extract <- paste0(OR,' [',low,' - ',high,']')
        #add N valid (complete cases)
        if(Ncases) {res.extract <- paste0(res.extract,' (n = ',glm.result$nsubs,')')}
        #add pvalue in result
        if(pval) {res.extract <- paste0(res.extract, '(p=',pvalue,')')}
        return (res.extract)
      
      }),stringsAsFactors = F
    )
  }else if (glm.family == 'gaussian'){
    stats <- data.frame(Estimate_CI = apply(glm.coef,1,function(x) {
        estimate <- round(x['Estimate'],rdigit)
        pvalue<- format(x['p-value'],digits=4) ; pvalue <- if(as.numeric(pvalue)< 2.2e-16) {"<2.2e-16"}else{pvalue}
        low <- round(x[length(x)-1],rdigit)
        high <- round(x[length(x)],rdigit)
        
        res.extract <- paste0(estimate,' [',low,' - ',high,']')
        #add N valid (complete cases)
        if(Ncases) {res.extract <- paste0(res.extract,' (n = ',glm.result$nsubs,')')}
        #add pvalue in result
        if(pval) {res.extract <- paste0(res.extract, '(p=',pvalue,')')}
        return (res.extract)
      }),stringsAsFactors = F
    )
  }
  glm.coef[2,dim(glm.coef)[2]-1]
  formula <-  gsub('\\s+',' ',formula(glm.result)) #fix spaces in formula
  list(formula = formula, stats = stats)
}


#################################################################
#DOC TO DO 
#
bioshare.env$run.NA.glm.subset<-function(formula,glm.result,NAsubset=NULL,datasources=NULL)
{
  mf <- match.call(expand.dots = FALSE)
  arg.names <- names(mf)
  
  if ((! 'glm.result' %in% arg.names) && (!'formula' %in% arg.names)) {stop('Either a glm run or a formula is required ...',call.=F)}
  if(!'formula' %in% arg.names){
    glm.formula <- gsub('\\s+','',glm.result$formula)
  }else {
    glm.formula <- formula
  }
  
  #parsing var.list from glm formula
  vars.list<-strsplit(glm.formula,'\\+|~|\\*|\\|+')
  
  if(is.null(NAsubset)){ 
    data <- unique(extract(vars.list)$holders)
    NAsubset <- paste0('NA.',data)
    warning(paste0("NAsubset is not specified ...the subset will be saved in ", NAsubset, " dataframe on server side"),call.=F,immediate.=T)
  }
  
  if(is.null(datasources)) { datasources = findLoginObjects()}
  ds <- datasources
  
  #define vars.2df and vars.names from var.list 
  vars.2df<-unlist(vars.list)
  vars.names<-extract(vars.list)$elements
  
  # create a df for the selected variables
  cally<-paste0('dataframeDS(list(',paste0(vars.2df,collapse=','),')',',NULL,FALSE,TRUE,','c(',paste0("'",vars.names,"'",collapse=','),')',',TRUE,FALSE)')
  datashield.assign(ds,'RD',as.symbol(cally))
  
  # define complete cases boolean var
  callcc<-'complete.cases(RD)'
  datashield.assign(ds,'cc',as.symbol(callcc))
  
  #define new df RD with cc variable
  callcbind<-'cbind(RD,cc)'
  datashield.assign(ds,'RD',as.symbol(callcbind))
  
  #define subset of RD --> RDF with values == NA in cc
  callSUBSET <- call('subsetDS', dt='RD', complt=FALSE, rs=NULL, cs=NULL, lg=5, th=FALSE, varname='cc')
  datashield.assign(ds, NAsubset, callSUBSET)
  
  #clean server workspace
  to_rm <- c("cc","complt","cs","dt","lg" , "rs", "th","varname","RD")
  invisible(sapply(to_rm,function(x) datashield.rm(ds,x)))
  #info to user
  cat(paste0("You may check the assigned ",NAsubset," dataframe with the following datashield commands: 
   datashield.symbols(opals), ds.colnames('",NAsubset,"') AND ds.dim('",NAsubset,"')")) 
  
  return(NAsubset)
}

#####################-----------STATS --------------------####################################

#this function compute desc statistic: table1d for a factor variable or meansd for a continuous variable
#across all studies and pool the results
bioshare.env$run.desc.stats<-function(var,data = NULL,datasources = NULL)
{
  if(is.null(datasources)) { datasources = findLoginObjects()}
  ds <- datasources
  
  #update var to na.data$var
  if(!is.null(data)) var<- paste0(data,'$',var)
  
  class.check <- checkClass(ds[1],var) 
  is.num <- class.check %in% c('integer','numeric')
  is.factor <- class.check == 'factor'
  is.class.null <- class.check == "NULL"
  
  if(is.factor) {
    tocall <- paste0('table1dDS(',var,')')
    rs <- datashield.aggregate(ds,as.name(tocall))
    
    #XXXXXXXX 
    rs.tab<-.vectorize(rs,'table')
    rs.mess <- .vectorize(rs,'message')
    #get the message validity
    rs.ok <- !any(grepl('invalid',rs.mess,ignore.case=T))
    #compute pooled tab
    Pooled <- apply(rs.tab,1,function(x) sum(as.numeric(x)))
    
    #put everything together
    rs.tab.with.pooled <- cbind(rs.tab,Pooled)
    k <- apply(rs.tab.with.pooled,2,function(x) {x<-as.numeric(x); pct <- round(x/x[length(x)]*100,2); paste0(pct,'(n = ',x,')')})
    res <- data.frame(k,row.names=row.names(rs.tab),stringsAsFactors=F)
    
  }else if (is.num){ #is numeric
    
    #tocall <- paste0('quantileMeanDS(',var,')')
    tocallmean <- paste0('meanDS(',var,')')
    tocallength <- paste0('length(',var,')')
    tocallnumna<-paste0('numNaDS(',var,')')
    tocallvar <- paste0('varDS(',var,')')
    rs.mean <- .vectorize(datashield.aggregate(ds,as.name(tocallmean)),1); 
    rs.numna <- .vectorize(datashield.aggregate(ds,as.name(tocallnumna)),1)
    rs.length <- .vectorize(datashield.aggregate(ds,as.name(tocallength)),1)
    rs.var <- .vectorize(datashield.aggregate(ds,as.name(tocallvar)),1)
    
    validN <- rs.length - rs.numna 
    
    .var.pooled <- function(variances,means,weights )
    {
      v<- variances
      m <- means
      l <- weights
      l.minus1 <- l-1
      m.total <- weighted.mean(m,l)
      err.ss <- sum(l.minus1*v) # overall error sum of squares
      tg.ss <- sum(l*(m-m.total)^2) # total (overall) group sum of squares
      l.total <- sum(l)
      v.total <- (err.ss + tg.ss)/(l.total-1)
      #s.ag <- sqrt(v.ag)
      return (v.total)
    }
    
  
    rs.mean.agg <- weighted.mean(rs.mean,validN,na.rm=T)
    rs.var.agg <- .var.pooled(rs.var, rs.mean, validN)
    validN.agg <- sum(validN)
    names.w.p <- c(names(ds),'Pooled')
    
    rs.mean.with.pooled <- c(rs.mean,rs.mean.agg)
    names(rs.mean.with.pooled) <- names.w.p
    rs.var.with.pooled <- c(rs.var,rs.var.agg)
    names(rs.var.with.pooled) <- names.w.p
    validN.w.p <- c(validN,validN.agg) 
    names(validN.w.p) <- names.w.p
     
    sd.w.p<-round(sqrt(rs.var.with.pooled),2)
    m.w.p <- round(rs.mean.with.pooled,2) 
    meanSd <- paste0(m.w.p,'(',sd.w.p,') (n = ',validN.w.p,')')
  
    res <- data.frame(meanSd,row.names=names.w.p,stringsAsFactors=F)
  } else if (is.class.null){
    stop(paste0('No such variable in ',names(ds[1])),call.=F)
  }  
  res
}


#this function computes 2 x 2 table and chi2 
bioshare.env$run.table2d <- function(x,y, data = NULL, col.percent = F,row.percent = F, chisq.test = T,split = F,datasources = NULL) 
{
  if(missing(x)) stop ('x variable (character) is required ...',call.=F)
  if(missing(y)) stop ('y variable (character) is required ...',call.=F)
  if(is.null(datasources)) datasources = findLoginObjects()
  ds <- datasources
  
  if(is.null(data)) callt2 <- paste0('table2dDS(',x,',',y,')')
  else callt2 <- paste0('table2dDS(',data,'$',x,',',data,'$',y,')')
  
  t2.res <-  datashield.aggregate(ds,as.symbol(callt2))
  
  process.result <- function (result)
  {
    with.pooled <- length(result) > 1
    t2.res <- result  #start here
    
    t2.mess <- .vectorize(t2.res,'message')
    t2.bad.idx <- which(grepl('invalid',t2.mess,ignore.case=T))
    t2.ok <- !(length(t2.bad.idx))
    
    t2d.i <- .vectorize(t2.res,subscript='table',simplify=F)
    
    #get basic infos
    cnames <- colnames(t2d.i[[1]])
    rnames <- row.names(t2d.i[[1]])
    #combine table
    if(with.pooled){
      Pooled <- Reduce('+',t2d.i)
      t2d.i.and.pool <- c(t2d.i,list(Pooled = Pooled))
      res <- t2d.i.and.pool
    }else{
      res <- t2d.i
    }
    
    # rm margin
    .rm.margin <- function(x) { as.matrix(x[-nrow(x),-ncol(x)])}
    t2.nomarg <- lapply(res,.rm.margin)
    
    if(row.percent && !col.percent){
      res <- lapply(t2.nomarg, function(x){
        p <- prop.table(x,1)*100
        p.marg<- margin.table(p,1)
        x.marg <- margin.table(x,1)
        p <- round(cbind(p,p.marg),2) #prop
        x <- round(cbind(x,x.marg),2) #count
        p.x <- paste0(p,'(n = ',x,')')
        p.x <- data.frame(matrix(p.x,nrow(p),ncol(p)))
        colnames(p.x) <- cnames
        row.names(p.x) <- row.names(x)
        p.x
      } )
    }
    if(col.percent){
      res <- lapply(t2.nomarg, function(x){
        p <- prop.table(x,2)*100
        p.marg <- margin.table(p,2)
        x.marg <- margin.table(x,2)
        p <- round(rbind(p,p.marg),2) #prop
        x <- round(rbind(x,x.marg),2) #count
        p.x <- paste0(p,'(n = ',x,')')
        p.x <- data.frame(matrix(p.x,nrow(p),ncol(p)))
        row.names(p.x) <- rnames
        colnames(p.x) <- colnames(x)
        p.x
      } )
    }
    
    #correct result: invalid result display total margin only
    if(!t2.ok) {
      res[t2.bad.idx] <- t2d.i[t2.bad.idx]
      if(with.pooled) res$Pooled <- t2d.i.and.pool$Pooled
    }
    
    #now compute chisquare test 
    if(chisq.test) {
      chi2 <-if(!t2.ok) {
        'Invalid table(s): all entries of the table must be nonnegative and finite'
      }else{
        if(with.pooled) {
          Pooled.data <- t2.nomarg$Pooled
          try(chisq.test(Pooled.data,correct=F),silent=T)
        }else{
          Study.data <- t2.nomarg[[1]]
          try(chisq.test(Study.data,correct=F),silent=T)
        }
      }
      res <- c(res,list(chi2 = chi2))
    }
    return (res)
  }
  
  message <- paste0('\n-----------  ',x,'(row) X ',y,'(col)  ------------')
  cat(message,'\n\n')
  if(split) {
    final <- c()
    for(i in 1:length(t2.res)){final <- c(final,process.result(t2.res[i])) }
    final
  }
  else process.result(t2.res)
}


##########################################

#close everything
bioshare.env$run.close<-function(all=F)
{
  #detect opals 
  objs <- ls(name=.GlobalEnv)
  if(length(objs)){
    for(obj in objs){
      obj.val<- eval(parse(text=obj))
      if (is.list(obj.val) && (class(obj.val[[1]]) == 'opal')){
        message('Closing opal(s) server connection(s)')
        obj.opal <- obj.val
        datashield.logout(obj.opal)
        rm(list= obj,pos=search())
        cat(paste0( names(obj.opal),' server is disconnected...'),sep='\n')
      }
    }
  }

 if(as.logical(all)) rm(list=objs,envir=.GlobalEnv)
 else rm(bioshare.env,pos=search())  
 detach(bioshare.env,pos=search())
 cat('bioshare environnment is now detached from memory...')  
}


# attach bioshare env
attach(bioshare.env)






