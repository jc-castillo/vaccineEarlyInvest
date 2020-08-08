
test_that("return numeric", {
  x = c(0.1,0.9)
  y = c(0.2,0.2,0.6)
  expect_is(sumdist(x,y),'numeric')
})

