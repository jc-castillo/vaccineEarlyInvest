#Test files for function loadData in models.R
setwd('~/GIT/vaccineEarlyInvest')
test_that('loadData report error if the includeVaccines contains vaccines that are not prescribed',{
  par = Parameters$new()
  expect_error(loadData(par = par, candidateFile = 'inst/extdata/vaccinesSummary.csv',includeVaccines = 1))
  expect_error(loadData(par = par, candidateFile = 'inst/extdata/vaccinesSummary.csv',includeVaccines = c('BCG',123)))
  expect_error(loadData(par = par, candidateFile = 'inst/extdata/vaccinesSummary.csv',includeVaccines = 'asd'))
})


test_that('includeVaccines works',{
  par = Parameters$new()
  d0 = loadData(par = par, candidateFile = 'inst/extdata/vaccinesSummary.csv')
  d1 = loadData(par = par, candidateFile = 'inst/extdata/vaccinesSummary.csv',includeVaccines = c('BCG'))
  expect_true('Live attenuated bacteria' %in% d1$Subcategory)
  expect_equal(setdiff(d1$Subcategory,d0$Subcategory),'Live attenuated bacteria')
})


test_that('loadData returns a data.table',{
  par = Parameters$new()
  d0 = loadData(par = par, candidateFile = 'inst/extdata/vaccinesSummary.csv')
  d1 = loadData(par = par, candidateFile = 'inst/extdata/vaccinesSummary.csv',includeVaccines = c('BCG'))
  expect_equal(class(d0),class(d1))
  expect_equal(class(d0),c('data.table','data.frame'))
})

test_that('loadData works with no candidateFile entered',{
  par = Parameters$new()
  d0 = loadData(par = par)
  d1 = loadData(par = par,includeVaccines = c('BCG'))
  expect_equal(class(d0),c('data.table','data.frame'))
  expect_equal(class(d1),c('data.table','data.frame'))
})

