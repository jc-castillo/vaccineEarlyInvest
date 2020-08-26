# test functions in countries.R
test_that('loadCountryData works',{
  expect_silent(loadCountryData(system.file('extdata/countryData.xlsx', package = 'vaccineEarlyInvest')))
  expect_equal(class(loadCountryData(system.file('extdata/countryData.xlsx', package = 'vaccineEarlyInvest'))),c('data.table','data.frame'))
  expect_error(loadCountryData(countrydata))
})

test_that('getBenefitsTable works',{
  countryData = loadCountryData(system.file('extdata/countryData.xlsx', package = 'vaccineEarlyInvest'))
  expect_silent(getBenefitsTable(countryData))
  benefitsTable = getBenefitsTable(countryData)
  expect_equal(class(benefitsTable),c('data.table','data.frame'))
  expect_silent(getBenefitsTable(countryData, max = 90,grid=0.1))
  expect_error(getBenefitsTable(countryData, max = -90,grid=0.1))
  expect_error(getBenefitsTable(countryData, max = -90,grid=-0.1))
  expect_error(getBenefitsTable(countryData, max = 90,grid=-0.1))
  expect_error(getBenefitsTable(countryData, max = 'a',grid=-0.1))
  expect_true(all(getBenefitsTable(countryData)>=0,na.rm=T))
})

test_that('countryParameters works',{
  countryData = loadCountryData(system.file('extdata/countryData.xlsx', package = 'vaccineEarlyInvest'))
  expect_silent(countryParameters(countryData))
  countryPar = countryParameters(countryData)
  expect_equal(class(countryPar),'list')
})

test_that('globalNetBenefits() works',{
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
  countryData = loadCountryData(system.file('extdata/countryData.xlsx', package = 'vaccineEarlyInvest'))
  benefitsTable = getBenefitsTable(countryData)
  expect_silent(globalNetBenefits(capacities = sample(seq(1,10),nrow(dcandidate),replace = T),dcandidate = dcandidate, targetPermutations = targetPermutations, 
                                  dplatforms=dplatforms,benefitsTable = benefitsTable,par = Parameters$new()))
  netben = globalNetBenefits(capacities = sample(seq(1,10),nrow(dcandidate),replace = T),dcandidate = dcandidate, targetPermutations = targetPermutations, 
                             dplatforms=dplatforms,benefitsTable = benefitsTable,par = Parameters$new())
  expect_equal(class(netben),'numeric')
  expect_length(netben,1)
  expect_error(globalNetBenefits(capacities = sample(seq(1,10),nrow(dcandidate)-1,replace = T),dcandidate = dcandidate, targetPermutations = targetPermutations, 
                                 dplatforms=dplatforms,benefitsTable = benefitsTable,par = Parameters$new()))
})

test_that('expectedBenefitsTable() works',{
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
  countryData = loadCountryData(system.file('extdata/countryData.xlsx', package = 'vaccineEarlyInvest'))
  benefitsTable = getBenefitsTable(countryData)
  dcandidate[,capacity:=sample(seq(1,10),nrow(dcandidate),replace = T)]
  dist = overallDistribution(dcandidate, targetPermutations, dplatforms, grid=1,
                             poverall=par$poverall, psubcat=par$psubcat)
  exp_ben = expectedBenefitsTable(dist,benefitsTable)
  expect_silent(expectedBenefitsTable(dist,benefitsTable))
  expect_equal(class(exp_ben),'numeric')
  expect_length(exp_ben,1)
  expect_error(expectedBenefitsTable(dist,benefitsTable,grid = -1))
})

test_that('globalNetBenefits() works',{
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
  countryData = loadCountryData(system.file('extdata/countryData.xlsx', package = 'vaccineEarlyInvest'))
  benefitsTable = getBenefitsTable(countryData)
  capacities = sample(seq(0,10),nrow(dcandidate),replace = T)
  expect_silent(globalNetBenefits(capacities,dcandidate,targetPermutations,dplatforms,benefitsTable,par))
  global_netben = globalNetBenefits(capacities,dcandidate,targetPermutations,dplatforms,benefitsTable,par)
  expect_length(global_netben,1)
  expect_equal(class(global_netben),'numeric')
})

test_that('globalNetBenefits() reports error for incorrect input',{
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
  countryData = loadCountryData(system.file('extdata/countryData.xlsx', package = 'vaccineEarlyInvest'))
  benefitsTable = getBenefitsTable(countryData)
  capacities = sample(seq(0,10),nrow(dcandidate),replace = T)
  expect_error(globalNetBenefits(capacities,dcandidate,targetPermutations,dplatforms,benefitsTable,par,grid=-1))
  expect_error(globalNetBenefits(capacities,dcandidate,targetPermutations,dplatforms,benefitsTable,grid=-1))
})