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
# require(foreign)
require(binomialRF)
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

################################################################################################
################################################################################################
################################################################################################


# #### Load data 
# rhv.pheno <- foreign::read.arff('~/Dropbox/Virogram/Java/ClassificationVirogram/HRV/summary.binary.signed.txt.arff')
# exvivo.pheno <- foreign::read.arff('~/Dropbox/Virogram/Java/ClassificationVirogram/Ex Vivo/summary.binary.signed.txt.arff')

#### load pathway ontology
go.bp.description <- fread("~/Dropbox/Virogram/Pathways/go.bp.description.txt")
go.bp.description$Pathway <- go.bp.description$path_id

# save(rhv.pheno, exvivo.pheno, file='~/Dropbox/Pathway-Pathway Interaction/data/hrv.RData')
load('~/Dropbox/Pathway-Pathway Interaction/data/hrv.RData')
################################################################################################
################################################################################################
################################################################################################

greaterThanPatients = 1

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

cbinom500 = read.csv('~/Dropbox/Pathway-Pathway Interaction/data/cbinom_500_interaction_trials.csv')
cbinom1000 = read.csv('~/Dropbox/Pathway-Pathway Interaction/data/cbinom_1000_interaction_trials.csv')
cbinom2000 = read.csv('~/Dropbox/Pathway-Pathway Interaction/data/cbinom_2000_interaction_trials.csv')

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

set.seed(1992)
rho=0.4
ntrees = 500
percent_features = 0.4

if(ntrees==500){
  cbinom=cbinom500$x
} else if (ntrees==1000){
  cbinom = cbinom1000$x
} else if (ntrees==2000){
  cbinom=cbinom2000$x
}

require(reshape)
source('~/Dropbox/Samir/binomialRF/R/k_binomialRF.R')

# teraction.matrix <- k.bin.rf$Interaction
# interaction.matrix$pathway1 <- as.character(interaction.matrix$pathway1)
# interaction.matrix$pathway2 <- as.character(interaction.matrix$pathway2)
# colnames(interaction.matrix)[1] <- 'path_id'
# 
# go.bp.description <- data.frame(go.bp.description)
# 
# interaction.matrix <- plyr::join(interaction.matrix, go.bp.description[,c('path_id','description')])
# interaction.matrix <- dplyr::left_join(interaction.matrix, go.bp.description[,c('path_id','description')], c('pathway2'='path_id'))
# 
# final.interaction.matrix <- cbind(k.bin.rf, interaction.matrix[, c('description.x','description.y')])
# fname = paste('~/Dropbox/Pathway-Pathway Interaction/results/pathway_pathway_interactions_patients',greaterThanPatients,'.csv',sep='')
# write.csv(final.interaction.matrix,fname , row.names = F)
# 
# pathways = unique(c(as.character(final.interaction.matrix$Interaction$pathway1),
#                     as.character(final.interaction.matrix$Interaction$pathway2)))
# 
# final.model = randomForest(rhv_df[,pathways], rhv_y,ntree = 500)

require(reshape)
require(binomialRF)

sample_sim <- function(greaterThanPatients, ntrees, training_datasets, rho=0.4){
  
  rhv_df = training_datasets$rhv_X
  rhv_y  = as.factor(ifelse(training_datasets$rhv_y=='symptomatic','Asymptomatic', 'Symptomatic'))
  
  exv_df = training_datasets$exvivo_df
  exv_y = as.factor(ifelse(training_datasets$exvivo_y == 'NE', 'Asymptomatic', 'Symptomatic'))
  
  ## train binomialRF interaction
  k.bin.rf <- k_binomialRF(rhv_df,rhv_y, ntrees =ntrees, percent_features = .9, fdr.threshold = 0.05, K=2, cbinom_dist = cbinom,
                           sampsize = nrow(rhv_df)*rho)
  # print(head(k.bin.rf))
  
  
  k.bin.rf <- k.bin.rf[k.bin.rf$significance < .05,]
  k.bin.rf <- data.frame(k.bin.rf, stringsAsFactors = F)
  
  if(nrow(k.bin.rf)>=2){
    k.bin.rf <- transform(k.bin.rf, Interaction = colsplit(Interaction, split = " \\| ", names = c('pathway1', 'pathway2')))
    pathways = unique(c(as.character(k.bin.rf$Interaction$pathway1), as.character(k.bin.rf$Interaction$pathway2)))
    # print(pathways)
    
    # combos = CJ(pathways, pathways, unique = TRUE)
    # crossproduct = function(XX,df){
    #   rhv_df = df
    #   b = rhv_df[unlist(combos[2,])]
    #   b$inter = b[,1]* b[,2]
    #   names(b)[3] = paste(names(b)[1],names(b)[2], sep='_')
    #   mat =  as.data.frame(b[,3]); names(mat) = names(b)[3]
    #   return(mat)
    # }
    # 
    # newrhv_DF = do.call(cbind, lapply(1:nrow(combos), function(xx) crossproduct(xx, rhv_df)))
    # new_exv_df= do.call(cbind, lapply(1:nrow(combos), function(xx) crossproduct(xx, exv_df)))
    #
    # idx = sample( nrow(rhv_df), size=13)
    
    baseRFmodel = randomForest(rhv_df, rhv_y,ntree = ntrees)
    implicitBinRFmodel = randomForest(rhv_df[,pathways], rhv_y,ntree = ntrees)
    # explicitBinRFmodel = randomForest(newrhv_DF, rhv_y,ntree = ntrees)
    
    result_df = data.frame(
      PatientFilter = greaterThanPatients,
      NumPathways=nrow(k.bin.rf),
      implicitBinomialRF = mean(predict(implicitBinRFmodel, exv_df)==exv_y),
      # explicitBinomialRF = mean(predict(explicitBinRFmodel, new_exv_df)==exv_y),
      baseRF = mean(predict(baseRFmodel, exv_df)==exv_y)
    
    )

  } else {
    result_df = data.frame(
      PatientFilter = greaterThanPatients,
      NumPathways=nrow(k.bin.rf),
      baseRF = NA, 
      implicitBinomialRF = NA
      # explicitBinomialRF= NA
    )
    
  }
  return(result_df)
}

nsim =50
print(paste('Conducting', nsim, 'simulations'))
cat('\n\n')
require(parallel)
require(data.table)

training_datasets = create_training_testing_sets(rhv.pheno, exvivo.pheno, greaterThanPatient=1)
filter1 = do.call(rbind, mclapply(1:nsim, function(x) sample_sim(greaterThanPatients=1, ntrees = ntrees, training_datasets=training_datasets)))
filter1 

training_datasets = create_training_testing_sets(rhv.pheno, exvivo.pheno, greaterThanPatient=3)
filter3 = do.call(rbind, mclapply(1:nsim, function(x) sample_sim(greaterThanPatients=3, ntrees = ntrees,training_datasets=training_datasets)))

training_datasets = create_training_testing_sets(rhv.pheno, exvivo.pheno, greaterThanPatient=5)
filter5 = do.call(rbind, mclapply(1:nsim, function(x) sample_sim(greaterThanPatients=5, ntrees = ntrees,training_datasets=training_datasets)))

training_datasets = create_training_testing_sets(rhv.pheno, exvivo.pheno, greaterThanPatient=7)
filter7 = do.call(rbind, mclapply(1:nsim, function(x) sample_sim(greaterThanPatients=7, ntrees = ntrees,training_datasets=training_datasets)))

full_results = data.table(rbind(filter1, filter3, filter5, filter7 ))
msd = function(val){ paste(round(mean(val, na.rm = T),3),' (',round(sd(val, na.rm=T),1),')',sep='')}

final_results = full_results[, list(baseRF= msd(baseRF), implicitRF =msd(implicitBinomialRF),
                                    InteractionsPredicted=msd(NumPathways)), by=PatientFilter]
final_results

write.csv(final_results,'~/Dropbox/Pathway-Pathway Interaction/results/HRVsignal_degradation500.csv')
save(final_results, full_results,file= '~/Dropbox/Pathway-Pathway Interaction/results/interaction_res500.RData')

# require(ggplot2)
# ggplot(full2, aes(x=PatientFilter, y=jitter(value), color=variable)) + geom_point()
