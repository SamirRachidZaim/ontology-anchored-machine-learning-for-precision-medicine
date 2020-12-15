# interpret results
rm(list=ls())
require(data.table)
internetwork = fread('~/Dropbox/Pathway-Pathway Interaction/results/figure_sampleInteraction_GO-MODULE.csv',stringsAsFactors = F)
internetwork

#### load pathway ontology
go.bp.description <- fread("~/Dropbox/Virogram/Pathways/go.bp.description.txt",stringsAsFactors = F)

interaction.matrix <- dplyr::left_join(internetwork, go.bp.description[,c('path_id','description')], c('Interaction.pathway1'='path_id'))
colnames(interaction.matrix)[3] = 'Description_pathway1'
interaction.matrix <- dplyr::left_join(interaction.matrix, go.bp.description[,c('path_id','description')], c('Interaction.pathway2'='path_id'))
colnames(interaction.matrix) = c('Pathway1','Pathway2','Description_pathway1','Description_pathway2')
interaction.matrix
write.csv(interaction.matrix, file='interaction_matrix_GO-BP.csv')
