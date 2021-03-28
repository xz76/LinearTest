
library(matrixStats)
library(truncdist)
library(survival)
library(extraDistr)

source("patp_test.R")
source("get_data.R")

## Example
tmatrix <- trans(state_names = c("health", "illness", "death"),
                 from = c(1, 1, 1, 2, 2),
                 to = c(2, 2, 3, 3, 1))
## Data
tmp <- get_data(n = 10, M0 = 10, M1 = 20)

## Point Estimator
point_res <- sop_t(tmp, S = 1:3, T_c = 1:2, ipw = 0, trans = tmatrix, times = NULL)

## Linear Test
test_res <- patp_test(tmp, tmat = tmatrix, j = 2, B = 1000)
