################################################################################################
################################################################################################
####
#### Script: rhinovirus.R
#### Description: pathway-pathway interactions binomialRF / HRV
#### author: Samir Rachid Zaim
#### date: 11/18/2020
####
################################################################################################
################################################################################################
rm(list=ls())

### load libraries
require(data.table)
require(randomForest)
require(foreign)
require(binomialRF)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

################################################################################################
################################################################################################
################################################################################################


#### Load data 
rhv.pheno <- foreign::read.arff('~/Dropbox/Virogram/Java/ClassificationVirogram/HRV/summary.binary.signed.txt.arff')
exvivo.pheno <- foreign::read.arff('~/Dropbox/Virogram/Java/ClassificationVirogram/Ex Vivo/summary.binary.signed.txt.arff')

#### load pathway ontology
go.bp.description <- fread("~/Dropbox/Virogram/Pathways/go.bp.description.txt")
go.bp.description$Pathway <- go.bp.description$path_id

################################################################################################
################################################################################################
################################################################################################


#### Data preprocessing 
#### format 

create_training_testing_sets <- function(rhv.pheno,exvivo.pheno, greaterThanPatient=greaterThanPatients){
  rhv_df = rhv.pheno[,! names(rhv.pheno) %in% c('Class', 'Name')]
  rhv_y = rhv.pheno[, names(rhv.pheno) %in% 'Class']
  
  exvivo_df<- exvivo.pheno[, !names(exvivo.pheno) %in% c("Class",'Name')]
  exvivo_y <-  exvivo.pheno[, names(exvivo.pheno) %in% c("Class")]
  
  
  full.intersect = intersect(colnames(rhv_df),  colnames(exvivo_df))
  rhv_df = rhv_df[, full.intersect]
  exvivo_df = exvivo_df[, full.intersect]
  
  rhv_df <- data.frame(do.call(cbind,lapply(1:ncol(rhv_df), function(zz) as.numeric(as.character(rhv_df[,zz]))))); names(rhv_df) <- full.intersect
  exvivo_df <- data.frame(do.call(cbind,lapply(1:ncol(exvivo_df), function(zz) as.numeric(as.character(exvivo_df[,zz]))))); names(exvivo_df) <- full.intersect
  
  remove_cols= which(sapply(1:ncol(rhv_df), function(x) all(rhv_df[,x]==0))==T)
  remove_col3= which(sapply(1:ncol(exvivo_df), function(x) all(exvivo_df[,x]==0))==T)
  removeCols = union(remove_cols, remove_col3)
  
  rhv_df = rhv_df[,-removeCols]
  exvivo_df=exvivo_df[,-removeCols]
  
  limited_columns =sapply(1:ncol(rhv_df), function(x) sum(rhv_df[,x] !=0))
  cols_to_keep = which(limited_columns > greaterThanPatients)
  
  rhv_df = rhv_df[, cols_to_keep]
  exvivo_df=exvivo_df[, cols_to_keep]
  
  return(list(rhv_X = rhv_df, 
              rhv_y = rhv_y, 
              exvivo_df= exvivo_df,
              exvivo_y= exvivo_y ))
}




#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### 
#### INTERACTION STUDY
####  
#### look at Pathway-Pathway Interaction 
#### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

cbinom500 = read.csv('~/Dropbox/Pathway-Pathway Interaction/data/cbinom_500trials.csv')
cbinom1000 = read.csv('~/Dropbox/Pathway-Pathway Interaction/data/cbinom_1000trials.csv')
cbinom2000 = read.csv('~/Dropbox/Pathway-Pathway Interaction/data/cbinom_2000trials.csv')

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

rho=0.4
ntrees = 2000

if(ntrees==500){
  cbinom=cbinom500$x
} else if (ntrees==1000){
  cbinom = cbinom1000$x
} else if (ntrees==2000){
  cbinom=cbinom2000$x
}

require(reshape)
require(binomialRF)

sample_sim <- function(greaterThanPatients, ntrees, training_datasets){
  
  rhv_df = training_datasets$rhv_X
  rhv_y  = as.factor(ifelse(training_datasets$rhv_y=='symptomatic','Asymptomatic', 'Symptomatic'))
  
  exv_df = training_datasets$exvivo_df
  exv_y = as.factor(ifelse(training_datasets$exvivo_y == 'NE', 'Asymptomatic', 'Symptomatic'))
  
  
  bin.rf <- binomialRF(rhv_df,rhv_y, ntrees =ntrees, percent_features = .8, fdr.threshold = 0.05, 
                       user_cbinom_dist = cbinom,sampsize =  nrow(rhv_df)*rho)
  print(head(bin.rf))
  
  bin.rf <- bin.rf[bin.rf$significance < .2,]
  bin.rf <- data.frame(bin.rf, stringsAsFactors = F)
  pathways = bin.rf$variable
  
  if(nrow(bin.rf) >= 2){
    naive.model = randomForest(rhv_df, rhv_y,ntree = 2000)
    final.model = randomForest(rhv_df[,pathways], rhv_y,ntree = 2000)
    
    result_df = data.frame(
      PatientFilter = greaterThanPatients,
      NumPathways=nrow(bin.rf),
      BinomialRF = mean(predict(final.model, exv_df)==exv_y),
      Naive = mean(predict(naive.model, exv_df)==exv_y)
    )
    
  } else {
    
    result_df = data.frame(
      PatientFilter = greaterThanPatients,
      NumPathways=nrow(bin.rf),
      BinomialRF = NA,
      Naive = NA
    )
    
  }
  print(result_df)
  return(result_df)
}


require(parallel)

training_datasets = create_training_testing_sets(rhv.pheno, exvivo.pheno, greaterThanPatient=1)
filter1 = do.call(rbind, mclapply(1:10, function(x) sample_sim(greaterThanPatients=1, ntrees = ntrees, training_datasets)))

training_datasets = create_training_testing_sets(rhv.pheno, exvivo.pheno, greaterThanPatient=3)
filter3 = do.call(rbind, mclapply(1:10, function(x) sample_sim(greaterThanPatients=3, ntrees = ntrees,training_datasets)))

training_datasets = create_training_testing_sets(rhv.pheno, exvivo.pheno, greaterThanPatient=5)
filter5 = do.call(rbind, mclapply(1:10, function(x) sample_sim(greaterThanPatients=5, ntrees = ntrees,training_datasets)))

training_datasets = create_training_testing_sets(rhv.pheno, exvivo.pheno, greaterThanPatient=7)
filter7 = do.call(rbind, mclapply(1:10, function(x) sample_sim(greaterThanPatients=7, ntrees = ntrees,training_datasets)))

full_results = data.table(rbind(filter1, filter3, filter5, filter7 ))
msd = function(val){ paste(round(mean(val),2),' (',round(sd(val),1),')',sep='')}

final_results = full_results[, list(NaiveRF= msd(Naive), BinomialRF =msd(BinomialRF), Pathways=msd(NumPathways)), by=PatientFilter]
final_results
