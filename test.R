## load the package
library(MNP)

## load the detergent data
data(detergent)
## run the standard multinomial probit model with intercepts and the price
res1 <- mnp(choice ~ 1, choiceX = list(Surf=SurfPrice, Tide=TidePrice,
                                       Wisk=WiskPrice, EraPlus=EraPlusPrice,
                                       Solo=SoloPrice, All=AllPrice),
            cXnames = "price", data = detergent, n.draws = 500, burnin = 100,
            thin = 3, verbose = TRUE)
## summarize the results
summary(res1)
## calculate the quantities of interest for the first 3 observations
pre1 <- predict(res1, newdata = detergent[1:3,])

## load the Japanese election data
data(japan)
## run the multinomial probit model with ordered preferences
res2 <- mnp(cbind(LDP, NFP, SKG, JCP) ~ gender + education + age, data = japan,
            verbose = TRUE)
## summarize the results
summary(res2)
## calculate the predicted probabilities for the 10th observation
## averaging over 100 additional Monte Carlo draws given each of MCMC draw.
pre2 <- predict(res2, newdata = japan[10,], type = "prob", n.draws = 100,
                verbose = TRUE)
}

