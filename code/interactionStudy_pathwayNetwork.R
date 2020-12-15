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
ntrees = 1000
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
  
  
  k.bin.rf <- k.bin.rf[k.bin.rf$adjSignificance < .2,]
  k.bin.rf <- data.frame(k.bin.rf, stringsAsFactors = F)
  return(k.bin.rf)
}

nsim =500
print(paste('Conducting', nsim, 'simulations'))
cat('\n\n')
require(parallel)
require(data.table)
require(reshape)

training_datasets = create_training_testing_sets(rhv.pheno, exvivo.pheno, greaterThanPatient=5)
filter5 = do.call(rbind, mclapply(1:nsim, function(x) sample_sim(greaterThanPatients=5, ntrees = ntrees,training_datasets=training_datasets)))
filter5 = data.table(filter5)
filter5$Interaction = filter5$Interaction

resfilter5 <- transform(filter5$Interaction, Interaction = colsplit(filter5$Interaction, split = "\\|", names = c('pathway1', 'pathway2')))
resfilter5 = data.table(resfilter5[,c('Interaction.pathway1','Interaction.pathway2')])

final_res = resfilter5[, .N, by=list(Interaction.pathway1,Interaction.pathway2)]
write.csv(final_res,'~/Dropbox/Pathway-Pathway Interaction/results/HRV_final_network.csv')

################################################################################################
################################################################################################
################################################################################################
### sample network 
# require(ggplot2)
# ggplot(full2, aes(x=PatientFilter, y=jitter(value), color=variable)) + geom_point()


set.seed(13)
rhv_df = training_datasets$rhv_X
rhv_y  = as.factor(ifelse(training_datasets$rhv_y=='symptomatic','Asymptomatic', 'Symptomatic'))

exv_df = training_datasets$exvivo_df
exv_y = as.factor(ifelse(training_datasets$exvivo_y == 'NE', 'Asymptomatic', 'Symptomatic'))

## train binomialRF interaction
k.bin.rf <- k_binomialRF(rhv_df,rhv_y, ntrees =ntrees, percent_features = .9, fdr.threshold = 0.05, K=2, cbinom_dist = cbinom,
                         sampsize = nrow(rhv_df)*rho)
print(head(k.bin.rf))


k.bin.rf <- k.bin.rf[k.bin.rf$adjSignificance < .2,]
k.bin.rf <- data.frame(k.bin.rf, stringsAsFactors = F)

resfilter5 <- transform(k.bin.rf$Interaction, Interaction = colsplit(k.bin.rf$Interaction, split = "\\|", names = c('pathway1', 'pathway2')))
resfilter5 = data.table(resfilter5[,c('Interaction.pathway1','Interaction.pathway2')])
write.csv(resfilter5, file='~/Dropbox/Pathway-Pathway Interaction/results/sampleInteraction.csv')
################################################################################################
################################################################################################

