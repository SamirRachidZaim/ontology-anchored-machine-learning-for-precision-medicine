# train cbinom distributions 
require(binomialRF)
set.seed(1)
rho=0.4
percent_features = 0.4
nfeatures= 2444
interaction_level = 2
interaction_prob = binomialRF::calculateBinomialP_Interaction(nfeatures, percent_features, K=interaction_level)


cbinom500 = correlbinom::correlbinom(rho, successprob = interaction_prob, trials=500, model='kuk')
cbinom1000 = correlbinom::correlbinom(rho, successprob = interaction_prob, trials=1000, model='kuk')
cbinom2000 = correlbinom::correlbinom(rho, successprob = interaction_prob, trials=2000, model='kuk')

write.csv(cbinom500, file='~/Dropbox/Pathway-Pathway Interaction/data/cbinom_500_interaction_trials.csv')
write.csv(cbinom1000, file='~/Dropbox/Pathway-Pathway Interaction/data/cbinom_1000_interaction_trials.csv')
write.csv(cbinom2000, file='~/Dropbox/Pathway-Pathway Interaction/data/cbinom_2000_interaction_trials.csv')
print('done')