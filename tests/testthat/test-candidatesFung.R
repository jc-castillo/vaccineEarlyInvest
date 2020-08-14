#Test files for function candidatesFung() in models.R
setwd('~/GIT/vaccineEarlyInvest')
test_that('returns a list with two named data.tables',{
  par = Parameters$new()
  d0 = loadData(par = par)
  d0$Target = 'Others'
  list1 = candidatesFung(d0,par,computeExpComp = F)
  list2 = candidatesFung(d0,par,computeExpComp = T)
  list3 = candidatesFung(d0,par,computeExpComp = F,seed = 5)
  expect_equal(class(list1),'list')
  expect_equal(names(list1),c('dordered','dcanddraws'))
  expect_equal(class(list2),'list')
  expect_equal(names(list2),c('dordered','dcanddraws'))
  expect_equal(class(list3),'list')
  expect_equal(names(list3),c('dordered','dcanddraws'))
})

test_that('seed needs to be an integer',{
  par = Parameters$new()
  d0 = loadData(par = par)
  expect_error(candidatesFung(d0,par,computeExpComp = F,seed = 'a'))
})
