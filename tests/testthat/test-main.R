test_that('demandPriceTaker() works',{
  skip_on_cran()
  population <- 31.99 # In millions
  gdp_pc <- 6.71 # In thousand dollars
  frac_high_risk <- 0.131
  loss2yr <- 0.269
  parameters <- Parameters$new(population=population, gdp_pc=gdp_pc, 
                               frac_high_risk=frac_high_risk, loss2yr=loss2yr)
  demand = demandPriceTaker(parameters)
  expect_equal(class(demand),'list')
  expect_length(demand,2)
  opt = demand$optimizations
  expect_true(all(opt>=0))
  expect_is(demand$allCapacities,c("matrix", "array" ))
  expect_true(all(demand$allCapacities>=0))
  demand = demandPriceTaker(population = population, gdp_pc = gdp_pc, frac_high_risk = frac_high_risk, loss2yr = loss2yr)
  expect_equal(class(demand),'list')
  expect_length(demand,2)
  opt = demand$optimizations
  expect_true(all(opt>=0))
  expect_is(demand$allCapacities,c("matrix", "array" ))
  expect_true(all(demand$allCapacities>=0))
})

test_that('portfolioPriceTaker() works',{
  skip_on_cran()
  population <- 31.99 # In millions
  gdp_pc <- 6.71 # In thousand dollars
  frac_high_risk <- 0.131
  loss2yr <- 0.269
  parameters <- Parameters$new(population=population, gdp_pc=gdp_pc, 
                               frac_high_risk=frac_high_risk, loss2yr=loss2yr)
  port = portfolioPriceTaker(parameters,price = 10)
  expect_equal(class(port),'list')
  expect_length(port,6)
  expect_true(all(port$capacities>=0))
  expect_gte(port$totCapacity,0)
  expect_gte(port$cost,0)
  expect_gte(port$expBenefits,0)
  expect_gte(port$expCapacity,0)
  expect_true(all(port$distribution>=0,na.rm = T))
  port = portfolioPriceTaker(population=population, gdp_pc=gdp_pc, 
                             frac_high_risk=frac_high_risk, loss2yr=loss2yr,price = 10)
  expect_equal(class(port),'list')
  expect_length(port,6)
  expect_true(all(port$capacities>=0))
  expect_gte(port$totCapacity,0)
  expect_gte(port$cost,0)
  expect_gte(port$expBenefits,0)
  expect_gte(port$expCapacity,0)
  expect_true(all(port$distribution>=0,na.rm = T))
})

