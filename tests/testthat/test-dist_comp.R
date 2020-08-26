# test a set of functions related to the computation of distributions, including candidateDraws(), 
# gettargetPermutationss(), permutationDistribution(), platformDistribution, subcatDistribution(),
# countryDistribution()


test_that("candidateDraws() works", {
  par = Parameters$new()
  d = loadData(par=par)
  d$Target <- "Other"
  d$Target[1:5]<-"Spike"
  d$Target[10:15]<-"Recombinant"
  cand = candidatesFung(d,par,seed = 10)
  draw1 = candidateDraws(cand$dordered[,1:11], par,seed=10)
  draw2 = candidateDraws(cand$dordered[,1:11], par,seed=20)
  expect_equal(class(draw1),c('data.table','data.frame'))
  expect_gt(sum(draw1$ysubcand!=draw2$ysubcand),0)
})


test_that("getTargetPermutations() works",{
  targets = c("Spike","Recombinant","Other")
  par = Parameters$new()
  probs = c(0.1, 0.5, 0.3)
  expect_equal(class(getTargetPermutations(targets,probs)),c('data.table','data.frame'))
  probs = c(1,1,1)
  expect_true(all(getTargetPermutations(targets,probs)$success==1))
})

test_that("getTargetPermutations() stops when probs are wrong",{
  targets = c("Spike","Recombinant","Other")
  par = Parameters$new()
  probs = c('a', 0.5, 0.3)
  expect_error(getTargetPermutations(targets,probs))
  probs = c(-1,1,1)
  expect_error(getTargetPermutations(targets,probs))
  probs = c(100,1,1)
  expect_error(getTargetPermutations(targets,probs))
  targets = c(1,2,4)
  expect_error(getTargetPermutations(targets,probs))
})

test_that("overallDistribution() works",{
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
  capacity = sample(seq(0,100,1),nrow(dcandidate),replace = T)
  dcandidate[,capacity:=capacity]
  expect_silent(overallDistribution(dcandidate = dcandidate,targetPermutations = targetPermutations, dplatforms = dplatforms,0.5,0.9))
  expect_equal(class(overallDistribution(dcandidate = dcandidate,targetPermutations = targetPermutations, dplatforms = dplatforms,0.5,0.9)),
               c('data.table','data.frame'))
  expect_silent(overallDistribution(dcandidate = dcandidate, dplatforms = dplatforms,poverall = 0.5,psubcat = 0.9,target=F))
  expect_equal(class(overallDistribution(dcandidate = dcandidate, dplatforms = dplatforms,poverall = 0.5,psubcat = 0.9,target=F)),
               c('data.table','data.frame'))
})

test_that("overallDistribution() reports error for incorrect input",{
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
  capacity = sample(seq(0,100,1),nrow(dcandidate),replace = T)
  dcandidate[,capacity:=capacity]
  expect_error(overallDistribution(dcandidate = dcandidate,targetPermutations = targetPermutations, dplatforms = dplatforms,0.5,100))
  expect_error(overallDistribution(dcandidate = dcandidate,targetPermutations = targetPermutations, dplatforms = dplatforms,-1,0.9))
  expect_error(overallDistribution(dcandidate = dcandidate,targetPermutations = targetPermutations, dplatforms = dplatforms,-1,0.9,grid = -1))
})

test_that("permutationDistribution works",{
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
  capacity = sample(seq(0,100,1),nrow(dcandidate),replace = T)
  dcandidate[,capacity:=capacity]
  dcandidate[, pcandperm := pcand]
  expect_silent(permutationDistribution(dcandidate = dcandidate, dplatforms = dplatforms,0.5,0.4))
  expect_equal(class(permutationDistribution(dcandidate = dcandidate, dplatforms = dplatforms,0.5,0.4)),c('data.table','data.frame'))
})


test_that("permutationDistribution() reports error for incorrect input",{
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
  capacity = sample(seq(0,100,1),nrow(dcandidate),replace = T)
  dcandidate[,capacity:=capacity]
  dcandidate[, pcandperm := pcand]
  expect_error(permutationDistribution(dcandidate = dcandidate, dplatforms = dplatforms,-0.5,0.4))
  expect_error(permutationDistribution(dcandidate = dcandidate, dplatforms = dplatforms,0.5,100.4))
})

test_that('subcatDistribution() works if sub is in dcandidate%sub',{
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
  capacity = sample(seq(0,100,1),nrow(dcandidate),replace = T)
  dcandidate[,capacity:=capacity]
  dcandidate[, pcandperm := pcand]
  expect_error(subcatDistribution('ds',dcandidate))
  expect_silent(subcatDistribution('S protein',dcandidate))
})

test_that('subcatDistribution() only works for character',{
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
  capacity = sample(seq(0,100,1),nrow(dcandidate),replace = T)
  dcandidate[,capacity:=capacity]
  dcandidate[, pcandperm := pcand]
  expect_error(subcatDistribution(c('S protein','Other DNA'),dcandidate))
})

test_that('platformDistribution() works if plat is in dcandidate%Platform',{
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
  capacity = sample(seq(0,100,1),nrow(dcandidate),replace = T)
  dcandidate[,capacity:=capacity]
  dcandidate[, pcandperm := pcand]
  subcatDists <- rbindlist(lapply(unique(dcandidate$Subcategory), subcatDistribution, dcandidate))
  setkey(subcatDists, Platform, Subcategory, capacity)
  expect_error(platformDistribution('asd',subcatDists,0.8))
  expect_silent(platformDistribution('RNA',subcatDists,0.8))
  expect_error(platformDistribution('RNA',subcatDists,-0.8))
})

test_that('platformDistribution() only works for character',{
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
  capacity = sample(seq(0,100,1),nrow(dcandidate),replace = T)
  dcandidate[,capacity:=capacity]
  dcandidate[, pcandperm := pcand]
  subcatDists <- rbindlist(lapply(unique(dcandidate$Subcategory), subcatDistribution, dcandidate))
  setkey(subcatDists, Platform, Subcategory, capacity)
  expect_error(platformDistribution(c('RNA','DNA'),subcatDists,0.8))
})

test_that('countryDistribution works',{
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
  capacity = sample(seq(0,10),nrow(dcandidate),replace = T)
  expect_silent(countryDistribution(capacity, dcandidate, targetPermutations,dplatforms,par=Parameters$new()))
  country_dist = countryDistribution(capacity, dcandidate, targetPermutations,dplatforms,par=Parameters$new())
  expect_equal(class(country_dist),c('data.table','data.frame'))
})

