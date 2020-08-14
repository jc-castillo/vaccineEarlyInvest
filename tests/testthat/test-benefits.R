#Test files for functions benefitIntegral(), benefitIntegralDisc(), and benefits() in models.R

test_that('benefitIntegral() takes in only frac between 0 and 1',{
  par = Parameters$new()
  expect_error(benefitIntegral(-1,par = par))
  expect_error(benefitIntegral(100, par = par))
  expect_error(benefitIntegral(c(0.1,200), par=par))
})

test_that('benefitIntegral() reports error for certain inputs',{
  par = Parameters$new()
  expect_error(benefitIntegral('a',par = par))
  expect_error(benefitIntegral(c('a',0.2),par = par))
})

test_that('benefitIntegral() works for differnet benefit distributions',{
  par = Parameters$new(benefitdist = "pnorm",alpha=0.5)
  benefitIntegral(0.1,par)
  par = Parameters$new(benefitdist = "piecewiseLinear")
  expect_equal(class(benefitIntegral(0.1,par)),'numeric')
  
})

test_that('benefitIntegral() works when frac is a vector',{
  par = Parameters$new()
  expect_vector(benefitIntegral(c(0.1,0.2,0.3),par = par))
  frac = c(0.1,0.2,0.3)
  expect_equal(length(frac),length(benefitIntegral(frac,par = par)))
})

test_that('benefitIntegral() returns numric',{
  par = Parameters$new()
  frac = c(0.1,0.2,0.3)
  expect_equal(class(benefitIntegral(0.1, par=par)),'numeric')
  expect_equal(class(benefitIntegral(frac = frac, par=par)),'numeric')
})

test_that('benefitIntegralDisc() work for vectors',{
  par = Parameters$new()
  expect_vector(benefitIntegralDisc(frac1 = rep(0,4),frac2 = seq(from=0.1,to=0.4,by=0.1), par = par))
  expect_vector(benefitIntegralDisc(frac1 = rep(0,4),frac2 = 1, par = par))
})

test_that('benefitIntegralDisc() does not work for non-numeric input',{
  par = Parameters$new()
  expect_error(benefitIntegralDisc(frac1 = 'a',frac2 = 0.3, par = par))
  expect_error(benefitIntegralDisc(frac1 = 0,frac2 = 'a', par = par))
  expect_error(benefitIntegralDisc(frac1 = 0,frac2 = 0.3, par = par,share1 = 'a'))
  expect_error(benefitIntegralDisc(frac1 = 0,frac2 = 0.3, par = par,share2 = 'a'))
})

test_that('benefitIntegralDisc() work if benefitdist is pnorm',{
  par1 = Parameters$new(benefitdist = 'pnorm',alpha=1)
  par0 = Parameters$new(benefitdist = 'pnorm',alpha=0.2)
  expect_equal(class(benefitIntegralDisc(frac1 = 0,frac2 = 0.3, par = par1)),'numeric')
  expect_vector(benefitIntegralDisc(frac1 = rep(0,4),frac2 = seq(from=0.1,to=0.4,by=0.1), par = par1))
  expect_equal(class(benefitIntegralDisc(frac1 = 0,frac2 = 0.3, par = par0)),'numeric')
  expect_vector(benefitIntegralDisc(frac1 = rep(0,4),frac2 = seq(from=0.1,to=0.4,by=0.1), par = par0))
})

test_that('benefitIntegralDisc() reports error unimplemented method',{
  par = Parameters$new()
  expect_error(benefitIntegral(frac1=0.1,frac2 = 0.3,par=par))
  expect_error(benefitIntegral(frac1=0,frac2 = 0.3,par=par,share1 = 0.5))
  expect_error(benefitIntegral(frac1=0,frac2 = 0.3,par=par,share2 = 0.5))
  par = Parameters$new(benefitdist = "piecewiseLinear")
  expect_error(benefitIntegral(frac1=0,frac2 = 0.3,par=par))
})

test_that('benefitIntegralDisc() reports error when frac and share are not between 0 and 1',{
  par = Parameters$new()
  expect_error(benefitIntegral(frac=-0.1,frac2 = 0.3,par=par))
  expect_error(benefitIntegral(frac=0,frac2 = 1.3,par=par))
  expect_error(benefitIntegral(frac=0,frac2 = 0.3,par=par,share1 = 1.5))
  expect_error(benefitIntegral(frac=0,frac2 = 0.3,par=par,share2 = 1.5))
  expect_error(benefitIntegral(frac=0,frac2 = 0.3,par=par,share1 = -1.5))
  expect_error(benefitIntegral(frac=0,frac2 = 0.3,par=par,share2 = -1.5))
})

test_that('benefits() reports error for incorrect input type/range',{
  par = Parameters$new()
  expect_error(benefits(100,list(c(0.5,1)),c(6,20),par))
  expect_error(benefits('a',list(c(0.5,1)),6,par))
  expect_error(benefits(100,list(c(0.5,1),c('a',0.7)),c(6,20),par))
  expect_error(benefits(-100,list(c(0.5,1),c(1,0.7)),c(6,20),par))
  expect_error(benefits(100,list(c(0.5,1),c(1,-0.7)),c(6,20),par))
  expect_error(benefits(-100,list(c(0.5,1),c(1,0.7)),c(-6,20),par))
  expect_error(benefits(100, 10,-10,par))
})
