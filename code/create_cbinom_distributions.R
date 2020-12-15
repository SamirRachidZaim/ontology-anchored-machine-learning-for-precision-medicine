# train cbinom distributions 
set.seed(1)
rho=0.4
percent_features = 0.4


cbinom500 = correlbinom::correlbinom(rho, successprob = 1/2444, trials=500, model='kuk')
cbinom1000 = correlbinom::correlbinom(rho, successprob = 1/2444, trials=1000, model='kuk')
cbinom2000 = correlbinom::correlbinom(rho, successprob = 1/2444, trials=2000, model='kuk')

write.csv(cbinom500, file='~/Dropbox/Pathway-Pathway Interaction/data/cbinom_500trials.csv')
write.csv(cbinom1000, file='~/Dropbox/Pathway-Pathway Interaction/data/cbinom_1000trials.csv')
write.csv(cbinom2000, file='~/Dropbox/Pathway-Pathway Interaction/data/cbinom_2000trials.csv')
print('done')