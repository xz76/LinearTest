README
================

### Description

R code for linear test.

### Main Functions

The main functions is `patp_test()` and `get_data()`. `get_data()` is a
function to create an example dataset.

### Dependencies

`patp_test()` requires the R packages `survival` and `matrixStats` to be
installed and loaded.`get_data()` requires the R packages `truncdist`
and `extraDistr` to be installed and loaded.

### Input data

The input data need to be a data frame in the long format generated by
`get_data()`. The data frame should contain the variables

-   `id`: variable name that identifies the individual observations.
-   `cid`: variable name that identifies the clusters.
-   `from`: the state of the process at Tstart. The possible values are
    1,…,k.
-   `Tstart`: starting time of the interval in the record.
-   `Tstop`: ending time of the interval in record.
-   `trans`: an integer that uniquely identifies the transition.
-   `status`: indicator variable. If status=1, the corresponding
    transition has been observed.
-   `group`: variable name of the binary grouping variable.

### Function `get_data()`

The function `get_data()` helps to create example data.

-   `n`: number of clusters.
-   `M0`: the minimum number of sample size within one cluster.
-   `M1`: the maximun number of sample size within one cluster.
-   `tmat`: a matrix of indicator transitions between states of the
    process where different transitions are identified by TRUE or FALSE.

### Function `patp_test()`

The function `patp_test()` calculates the p-value for the comparison of
the population-averaged transition probability
*P**r*(*X*(*t*) = *j*\|*X*(*s*) = *h*) between two groups, using a
linear test. The function performs has following arguments:

-   `data`: a data.frame in the long format follows `get_data()`
    requirements.
-   `tmat`: a matrix of indicator transitions between states of the
    process where different transitions are identified by TRUE or FALSE.
-   `id`: variable name that identifies the individual observations.
-   `cid`: variable name that identifies the clusters.
-   `group`: variable name of the binary grouping variable.
-   `h`: the state h in *P**r*(*X*(*t*) = *j*\|*X*(*s*) = *h*).
-   `j`: the state j in *P**r*(*X*(*t*) = *j*\|*X*(*s*) = *h*).
-   `s`: the time s in *P**r*(*X*(*t*) = *j*\|*X*(*s*) = *h*). The
    default value is 0.
-   `weighted`: logical value. If TRUE, the estimators are weighted by
    the inverse of cluster sizes. This is useful when cluster size is
    random and expected to be informative. The default value is FALSE.
-   `LMAJ`: logical value. If TRUE, the landmark version of the
    estimator is used in the test. This is useful when s&gt;0 and the
    Markov assumption is not plausible. The default value is FALSE.
-   `B`: number of nonparametric cluster bootstrap replications. The
    default value is 1000.

### Example

The artificial dataset generated by `get_data()` contains clustered
observations from an illness-death process without recovery . The matrix
`tmatrix` of possible transition looks as follows.

``` r
tmatrix <- trans(state_names = c("health", "illness", "death"),from = c(1, 1, 1, 2, 2),
                 to = c(2, 2, 3, 3, 1))
tmatrix
```

    ##         health illness death
    ## health   FALSE    TRUE  TRUE
    ## illness   TRUE   FALSE  TRUE
    ## death    FALSE   FALSE FALSE

Create an example dataset with 10 clusters called `tmp`. The dataset can
be obtained as follows

``` r
tmp <- get_data(n = 10, M0 = 10, M1 = 20, tmat = tmatrix)
head(tmp)
```

    ##     id cid from to trans    Tstart     Tstop      time status group
    ## 1.1  1   1    1  2     1 0.0000000 0.6449366 0.6449366      1     0
    ## 1.2  1   1    2  1     3 0.6449366 0.9683062 0.3233695      0     0
    ## 1.3  1   1    2  3     4 0.6449366 0.9683062 0.3233695      0     0
    ## 1.4  1   1    1  2     1 0.0000000 0.6449366 0.6449366      1     0
    ## 1.5  1   1    2  1     3 0.6449366 0.9683062 0.3233695      0     0
    ## 1.6  1   1    2  3     4 0.6449366 0.9683062 0.3233695      0     0

Two-sample comparison of the transition probability
*P*(*X*(*t*) = 2\|*X*(0) = 1) between the groups defined by the variable
group can be performed as follows

``` r
set.seed(1234)
patp_test(data = tmp, tmat = tmatrix, cid = "cid", id = "id", group = "group",
          h = 1, j = 2, s = 0, weighted = FALSE, LMAJ = FALSE, B = 1000)
```

    ## p-value at State2 
    ##         0.6167089

It is recommended to use at least 1000 cluster bootstrap replications
when performing two-sample hypothesis testing.
