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
                                         Wisk=WiskPrice, EraPlus=EraPlusPrice, 
                                         Solo=SoloPrice, All=AllPrice), 
              cXnames = "price", data = detergent, n.draws = 500, burnin = 100, 
              thin = 3, verbose = TRUE)
  # summarize the results
  # x <- summary(res1)
  # expect_that(length(x), is_equivalent_to(8))
  # expect_true("coef.table" %in% names(x))
  # expect_that(round(x$coef.table[2,2], 3), is_equivalent_to(0.115))
  # expect_that(round(x$coef.table["price", "mean"], 1), is_equivalent_to(-63.4))
  # expect_that(round(x$cov.table[2,4], 3), is_equivalent_to(0.721))
  # expect_that(round(x$cov.table["Tide:Tide", "mean"], 3), is_equivalent_to(0.658))
})

if (0) {
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
  expect_that(round(x$coef.table[2,2], 3), is_equivalent_to(0.115))
  expect_that(round(x$coef.table["price", "mean"], 1), is_equivalent_to(-63.4))
  expect_that(round(x$cov.table[2,4], 3), is_equivalent_to(0.721))
  expect_that(round(x$cov.table["Tide:Tide", "mean"], 3), is_equivalent_to(0.658))
  
  # calculate the quantities of interest for the first 3 observations
  x <- predict(res1, newdata = detergent[1:3,])
  expect_that(length(x), is_equivalent_to(4))
  expect_true("p" %in% names(x))
  expect_that(dim(x$o), is_equivalent_to(c(3, 6, 100)))
  expect_that(x$o[1,2,3], is_equivalent_to(4))
  expect_that(round(x$p[2, "Tide"], 2), is_equivalent_to(0.33))
})  
}

# set random seed
set.seed(12345)

if (0) {
test_that("tests MNP on the Japanese election census", {
  # load the Japanese election data
  data(japan)
  # run the multinomial probit model with ordered preferences
  res2 <- mnp(cbind(LDP, NFP, SKG, JCP) ~ gender + education + age, data = japan, verbose = TRUE)
  # summarize the results
  x <- summary(res2)
  expect_that(length(x), is_equivalent_to(8))
  expect_true("coef.table" %in% names(x))
  expect_that(round(x$coef.table[2,1], 3), is_equivalent_to(1.129))
  expect_that(round(x$coef.table["education:LDP", "mean"], 3), is_equivalent_to(-0.109))
  expect_that(round(x$cov.table[2,3], 3), is_equivalent_to(0.937))
  expect_that(round(x$cov.table["LDP:NFP", "mean"], 3), is_equivalent_to(1.017))
  
  # calculate the predicted probabilities for the 10th observation
  # averaging over 100 additional Monte Carlo draws given each of MCMC draw.
  x <- predict(res2, newdata = japan[10,], type = "prob", n.draws = 100, verbose = TRUE)
  expect_that(length(x), is_equivalent_to(2))
  expect_true("p" %in% names(x))
  expect_that(dim(x$p), is_equivalent_to(c(1, 4, 5000)))
  expect_that(x$p[1,2,3], is_equivalent_to(0.33))
  expect_that(round(x$p[1, "JCP", 5000], 2), is_equivalent_to(0.13))
})  
}