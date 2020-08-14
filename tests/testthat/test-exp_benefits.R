# test functions in distributions.R that are related to the calculation of cost and benefits
setwd('~/GIT/vaccineEarlyInvest')
?socialCost

test_that('priceTakerCost() works only for sensible inputs',{
  expect_gte(priceTakerCost(c(10,12,140),10),0)
  expect_equal(class(priceTakerCost(124,10)),'numeric')
  expect_length(priceTakerCost(c(10,12,140),10),1)
  expect_error(priceTakerCost(c(-10,12,140),10))
  expect_error(priceTakerCost(c(10,12,140),-10))
})

test_that('socialCost works only for sensible inputs',{
  targets = c("Spike","Recombinant","Other")
  probs = c(0.1, 0.1, 0.3)
  targetPermutations = getTargetPermutations(targets,probs)
  par = Parameters$new()
  d = loadData(par=par)
  d$Target <- "Other"
  d$Target[1:5]<-"Spike"
  d$Target[10:15]<-"Recombinant"
  dordered <- candidatesFung(d, par)$dordered
  dordered <- dordered[,1:11]
  dcandidate = copy(dordered)
  dplatforms <- unique(dordered[, .(Platform, pplat)])
  capacities = sample(seq(0,10),nrow(dcandidate),replace = T)
  dcandidate[,capacity:=capacities]
  dist = overallDistribution(dcandidate, targetPermutations, dplatforms, grid=1,
                             poverall=par$poverall, psubcat=par$psubcat)
  cost = socialCost(capacities,dist,par)
  expect_gte(cost,0)
  expect_equal(class(cost),'numeric')
  expect_length(cost,1)
  capacities[1] = -1
  expect_error(socialCost(capacities,dist,par))
  expect_error(socialCost(capacities,dist,par))
})

test_that('countryExpectedBenefits works only for sensible inputs',{
  targets = c("Spike","Recombinant","Other")
  probs = c(0.1, 0.1, 0.3)
  targetPermutations = getTargetPermutations(targets,probs)
  par = Parameters$new()
  d = loadData(par=par)
  d$Target <- "Other"
  d$Target[1:5]<-"Spike"
  d$Target[10:15]<-"Recombinant"
  dordered <- candidatesFung(d, par)$dordered
  dordered <- dordered[,1:11]
  dcandidate = copy(dordered)
  dplatforms <- unique(dordered[, .(Platform, pplat)])
  capacities = sample(seq(0,10),nrow(dcandidate),replace = T)
  benefits = countryExpectedBenefits(capacities,dcandidate,targetPermutations,dplatforms,par)
  expect_gte(benefits,0)
  expect_equal(class(benefits),'numeric')
  expect_length(benefits,1)
  expect_error(countryExpectedBenefits(capacities,dcandidate,targetPermutations,dplatforms,par,grid = -1))
})

?countryNetBenefits
test_that('countryNetBenefits() works only for sensible input',{
  targets = c("Spike","Recombinant","Other")
  probs = c(0.1, 0.1, 0.3)
  targetPermutations = getTargetPermutations(targets,probs)
  par = Parameters$new()
  d = loadData(par=par)
  d$Target <- "Other"
  d$Target[1:5]<-"Spike"
  d$Target[10:15]<-"Recombinant"
  dordered <- candidatesFung(d, par)$dordered
  dordered <- dordered[,1:11]
  dcandidate = copy(dordered)
  dplatforms <- unique(dordered[, .(Platform, pplat)])
  countryData = loadCountryData('inst/extdata/countryData.xlsx')
  benefitsTable = getBenefitsTable(countryData)
  capacities = sample(seq(0,10),nrow(dcandidate),replace = T)
  netben = countryNetBenefits(capacities,dcandidate,targetPermutations,dplatforms,grid = 1, price = 10,par=par)
  expect_length(netben,1)
  expect_equal(class(netben),'numeric')
  expect_error(countryNetBenefits(capacities,dcandidate,targetPermutations,dplatforms,grid = 1, price = -10,par=par))
  expect_error(countryNetBenefits(capacities,dcandidate,targetPermutations,dplatforms,grid = -1, price = 10,par=par))
})
