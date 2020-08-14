#Test files for methods of Parameters
test_that("clone works", {
  par = Parameters$new()
  par1 = par$clone()
  par2 = par$clone(deep = T)
  expect_true(all.equal.environment(par,par1))
  expect_true(all.equal.environment(par,par2))
})

test_that("setDerivedParameters() works", {
  par = Parameters$new(mortality = NULL)
  expect_silent(par$setDerivedParameters())
  par = Parameters$new(global = T)
  expect_silent(par$setDerivedParameters())
  par = Parameters$new(global = T,benefitdist = 'pnorm')
  expect_silent(par$setDerivedParameters())
  par = Parameters$new(global = T,benefitdist = 'piecewiseLinear')
  expect_silent(par$setDerivedParameters())
})

test_that("Parameters$new reports error if only one of pop, gdp, and frac is provided",{
  expect_error(Parameters$new(population = 1))
  expect_error(Parameters$new(gdp_pc = 342))
  expect_error(Parameters$new(frac_high_risk = 0.1))
})
