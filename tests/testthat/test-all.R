rm(list=ls())
library(MNP)
library(testthat)
context("tests MNP")

# set random seed
set.seed(12345)

test_that("tests MNP on the detergent data", {
  # load the detergent data
  data(detergent)
  # run the standard multinomial probit model with intercepts and the price
  res1 <- mnp(choice ~ 1, choiceX = list(Surf=SurfPrice, Tide=TidePrice,
              Wisk=WiskPrice, EraPlus=EraPlusPrice, Solo=SoloPrice, All=AllPrice),
              cXnames = "price", data = detergent, n.draws = 500, burnin = 100,
              thin = 3, verbose = TRUE)
  # summarize the results
  x <- summary(res1)
  expect_that(length(x), is_equivalent_to(8))
  expect_true("coef.table" %in% names(x))
  expect_that(round(x$coef.table[4, 1], 0), equals(2))
  expect_that(round(x$coef.table["(Intercept):Solo", "2.5%"], 0), equals(1))
  expect_that(round(x$cov.table[10, 1], 0), equals(1))
  expect_that(round(x$cov.table["Tide:Tide", "std.dev."], 0), equals(0))
  
  # calculate the quantities of interest for the first 3 observations
  x <- predict(res1, newdata = detergent[1:3,])
  expect_that(length(x), is_equivalent_to(4))
  expect_true("p" %in% names(x))
  expect_that(dim(x$o), is_equivalent_to(c(3, 6, 100)))
  expect_that(as.numeric(round(x$p[1, "Tide"], 1)), equals(0.1))
  expect_that(as.numeric(round(x$p[3, "Wisk"], 1)), equals(0.3))
})  


# set random seed
set.seed(12345)

test_that("tests MNP on the Japanese election census", {
  # load the Japanese election data
  data(japan)
  # run the multinomial probit model with ordered preferences
  res2 <- mnp(cbind(LDP, NFP, SKG, JCP) ~ gender + education + age, data = japan, verbose = TRUE)
  # summarize the results
  x <- summary(res2)
  expect_that(length(x), is_equivalent_to(8))
  expect_true("coef.table" %in% names(x))
  expect_that(round(x$coef.table[12,2], 1), is_equivalent_to(0.0))
  expect_that(round(x$coef.table["education:SKG", "mean"], 0), is_equivalent_to(0))
  expect_that(round(x$cov.table[2,1], 0), is_equivalent_to(1))
  expect_that(round(x$cov.table["LDP:LDP", "mean"], 0), is_equivalent_to(1))
  
  # calculate the predicted probabilities for the 10th observation
  # averaging over 100 additional Monte Carlo draws given each of MCMC draw.
  x <- predict(res2, newdata = japan[10,], type = "prob", n.draws = 100, verbose = TRUE)
  expect_that(length(x), is_equivalent_to(2))
  expect_true("p" %in% names(x))
  expect_that(dim(x$p), is_equivalent_to(c(1, 4, 5000)))
  expect_that(x$x[1, "age:LDP"], is_equivalent_to(50))
  expect_that(x$x[1, 1], is_equivalent_to(1))
})  

# set random seed
set.seed(123456)

test_that("tests MNP to discover the difference between local and travis-ci", {
  # load the detergent data
  data(detergent)
  # run the standard multinomial probit model with intercepts and the price
  res1 <- mnp(choice ~ 1, choiceX = list(Surf=SurfPrice, Tide=TidePrice, Wisk=WiskPrice, 
                                         EraPlus=EraPlusPrice, Solo=SoloPrice, All=AllPrice),
              cXnames = "price", data = detergent, n.draws = 100, burnin = 10,thin = 3, 
              verbose = TRUE)
  # summarize the results
  x <- summary(res1)
  expect_that(length(x), is_equivalent_to(8))
  expect_true("coef.table" %in% names(x))
  
  ############################################################
  # this only works for travis-ci
  # expect_that(round(x$coef.table[4, 1], 5), equals(2.00358))
  # this only works for local "R CMD check --as-cran" (random seed 12345)
  # expect_that(round(x$coef.table[4, 1], 5), equals(2.01363))
  # with random seed 123456, the above should be
  expect_that(round(x$coef.table[4, 1], 5), equals(1.9492))
  ############################################################
  
  # this happen to works for both local "R CMD check --as-cran" and travis-ci
  # expect_that(round(x$coef.table["(Intercept):Solo", "2.5%"], 3), equals(1.077))
  # the previous works for random seed 12345, with random seed 123456, it should be
  expect_that(round(x$coef.table["(Intercept):Solo", "2.5%"], 3), equals(1.033))
})  

