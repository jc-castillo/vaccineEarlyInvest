#Test files for functions sumdist() and sumdistSelf() in distributions.R

test_that('Report errors when the input is not numeric',{
  dist1 = c(0.1,0.2,0.7)
  dist2 = c(0.1,'a')
  dist3 = c('a','b')
  expect_error(sumdist(dist1,dist2))
  expect_error(sumdist(dist1,dist3))
  expect_error(sumdistSelf(dist2))
  expect_error(sumdistSelf(dist3))
})

test_that('Report error if the input is not a distribution',{
  dist_neg = c(-0.5,1.5)
  dist1 = c(0.1,0.9)
  dist_no1 = c(1,2,4)
  expect_error(sumdist(dist_neg,dist1))
  expect_error(sumdist(dist1,dist_neg))
  expect_error(sumdist(dist1,dist_no1))
  expect_error(sumdist(dist_no1,dist1))
  expect_error(sumdistSelf(dist_neg))
  expect_error(sumdistSelf(dist_no1))
})

test_that('Return a discrete distribution',{
  dist1 = c(0.1,0.2,0.7)
  dist2 = c(0.2,0.8)
  expect_equal(sum(sumdist(dist1,dist2)),1)
  expect_equal(sum(sumdist(dist1,dist1)),1)
  expect_equal(sum(sumdist(dist2,dist2)),1)
  expect_equal(sum(sumdistSelf(dist1)),1)
})
