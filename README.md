
<!-- README.md is generated from README.rmd. Please edit that file -->

    ## ℹ Loading plt

# plt

Persistence landscapes are a vectorization of persistence data (also
called persistence diagrams), [originated by Peter
Bubenik](https://jmlr.csail.mit.edu/papers/v16/bubenik15a.html), that
have useful statistical properties including linearity and an inner
product. This is an R package interface to [Paweł Dłotko’s Persistence
Landscapes
Toolbox](https://www2.math.upenn.edu/~dlotko/persistenceLandscape.html),
[developed with
Bubenik](https://www.sciencedirect.com/science/article/pii/S0747717116300104)
to efficienctly compute and calculate with persistence landscapes. It
was adapted from [Jose Bouza’s **tda-tools**
package](https://github.com/jjbouza/tda-tools), a TDA pipeline used by
[Bubenik’s
lab](https://people.clas.ufl.edu/peterbubenik/researchgroup/).

## Installation

Since the package is not on CRAN, we suggest that you install
**remotes** on your system. To do this just open an R session and run:

``` r
install.packages("remotes")
```

Then install the package from the GitHub repository as follows:

``` r
remotes::install_github("corybrunson/plt")
```

Alternatively, you can clone or download the code repository and install
the package from source (from within the directory):

``` r
devtools::install()
```

You should now be able to load the package normally from an R session:

``` r
library(plt)
```

## Quickstart Guide

The **plt** package supports various operations involving persistence
landscapes:

-   Compute persistence landscapes from persistence data
-   Perform Hilbert space operations (vector space operations and an
    inner product) on persistence landscapes
-   Plot persistence landscapes

Examples and tests in **plt** rely on other packages to simulate data
and to compute persistence diagrams from data:

-   **tdaunif** provides functions to sample uniformly from various
    immersed manifolds.
-   **ripserr** and **TDA** provide functions to compute persistence
    data from point clouds and distance matrices.

### Calculation

**plt** introduces the ‘Rcpp_PersistenceLandscape’ S4 class, which is
exposed using **Rcpp** from the underlying ‘PersistenceLandscape’ C++
class. Instances of this class can be created using `new()` but the
recommended way is to use `landscape()`. This function accepts either a
single matrix of persistence data or a persistence diagram object, which
is a specially formatted list containing the first and most important
element `$pairs`. The `$pairs` entry is itself a list, of a 2-column
matrix of persistence pairs for each homological degree from 0
(`$pairs[[1]]`) to the maximum degree calculated. The generic conversion
function `as_persistence()` includes methods for outputs from
`ripserr::vietoris_rips()` and from `TDA::*Diag()`.

To begin an illustration, we noisily sample 60 points from a figure
eight and compute the persistence diagram of the point cloud:

``` r
set.seed(513611L)
pc <- tdaunif::sample_lemniscate_gerono(60, sd = .1)
pd <- ripserr::vietoris_rips(pc, dim = 1, threshold = 2, p = 2)
print(pd)
#>    dimension     birth      death
#> 1          0 0.0000000 0.01918952
#> 2          0 0.0000000 0.01947548
#> 3          0 0.0000000 0.02604350
#> 4          0 0.0000000 0.04218479
#> 5          0 0.0000000 0.04542467
#> 6          0 0.0000000 0.05941691
#> 7          0 0.0000000 0.06030423
#> 8          0 0.0000000 0.06260854
#> 9          0 0.0000000 0.06478082
#> 10         0 0.0000000 0.06766925
#> 11         0 0.0000000 0.07158685
#> 12         0 0.0000000 0.07398253
#> 13         0 0.0000000 0.07623591
#> 14         0 0.0000000 0.07662517
#> 15         0 0.0000000 0.07896342
#> 16         0 0.0000000 0.08043306
#> 17         0 0.0000000 0.08642298
#> 18         0 0.0000000 0.08698163
#> 19         0 0.0000000 0.08850341
#> 20         0 0.0000000 0.08850349
#> 21         0 0.0000000 0.09024332
#> 22         0 0.0000000 0.09502309
#> 23         0 0.0000000 0.09673399
#> 24         0 0.0000000 0.10001714
#> 25         0 0.0000000 0.10007503
#> 26         0 0.0000000 0.10147098
#> 27         0 0.0000000 0.10290514
#> 28         0 0.0000000 0.11388282
#> 29         0 0.0000000 0.11527285
#> 30         0 0.0000000 0.12074671
#> 31         0 0.0000000 0.12093951
#> 32         0 0.0000000 0.12170977
#> 33         0 0.0000000 0.13435379
#> 34         0 0.0000000 0.13492348
#> 35         0 0.0000000 0.13563400
#> 36         0 0.0000000 0.13846410
#> 37         0 0.0000000 0.14438541
#> 38         0 0.0000000 0.14806602
#> 39         0 0.0000000 0.16101999
#> 40         0 0.0000000 0.16433348
#> 41         0 0.0000000 0.16487642
#> 42         0 0.0000000 0.17046509
#> 43         0 0.0000000 0.18270262
#> 44         0 0.0000000 0.18502783
#> 45         0 0.0000000 0.18551972
#> 46         0 0.0000000 0.19104680
#> 47         0 0.0000000 0.19117990
#> 48         0 0.0000000 0.19144013
#> 49         0 0.0000000 0.19311162
#> 50         0 0.0000000 0.19403676
#> 51         0 0.0000000 0.20145942
#> 52         0 0.0000000 0.20234674
#> 53         0 0.0000000 0.20856429
#> 54         0 0.0000000 0.21921402
#> 55         0 0.0000000 0.24334202
#> 56         0 0.0000000 0.24700342
#> 57         0 0.0000000 0.24971202
#> 58         0 0.0000000 0.25881722
#> 59         0 0.0000000 0.37692215
#> 60         1 0.4809292 0.63582254
#> 61         1 0.3016234 0.60751718
#> 62         1 0.2504500 0.27279150
#> 63         1 0.2251884 0.23008714
```

We the convert the persistence data to the preferred persistence diagram
format and inspect some of its features:

``` r
pd <- as_persistence(pd)
print(head(pd$pairs[[1]]))
#>      [,1]       [,2]
#> [1,]    0 0.01918952
#> [2,]    0 0.01947548
#> [3,]    0 0.02604350
#> [4,]    0 0.04218479
#> [5,]    0 0.04542467
#> [6,]    0 0.05941691
print(head(pd$pairs[[2]]))
#>           [,1]      [,2]
#> [1,] 0.4809292 0.6358225
#> [2,] 0.3016234 0.6075172
#> [3,] 0.2504500 0.2727915
#> [4,] 0.2251884 0.2300871
```

This allows us to compute a persistence landscape—in this case, for the
1-dimensional features. Here we compute the landscape exactly, which can
be cost-prohibitive for larger persistence data:

``` r
pl1 <- landscape(pd, degree = 1, exact = TRUE)
print(pl1)
#> Persistence landscape (exact format) of 2 envelopes over (0.2,0.6)
```

### Class

The object `pl1` is not an array, but rather an object that encapsulates
both the data that encode a landscape and several basic operations that
can be performed on it. This allows us to work with persistence
landscapes without worrying about pre-processing their representations.
At any point, the underlying landscape can be extracted using
`$getInternal()`, which in the case of an exactly calculated landscape
returns a list of 2-column matrices, each matrix containing coordinates
that define one level of the landscape as a piecewise linear function:

``` r
print(length(pl1$getInternal()))
#> [1] 2
print(pl1$getInternal())
#> [[1]]
#>            [,1]        [,2]
#>  [1,]      -Inf 0.000000000
#>  [2,] 0.2251884 0.000000000
#>  [3,] 0.2276378 0.002449358
#>  [4,] 0.2300871 0.000000000
#>  [5,] 0.2504500 0.000000000
#>  [6,] 0.2616207 0.011170764
#>  [7,] 0.2727915 0.000000000
#>  [8,] 0.3016234 0.000000000
#>  [9,] 0.4545703 0.152946885
#> [10,] 0.5442232 0.063293964
#> [11,] 0.5583759 0.077446647
#> [12,] 0.6358225 0.000000000
#> [13,]       Inf 0.000000000
#> 
#> [[2]]
#>           [,1]       [,2]
#> [1,]      -Inf 0.00000000
#> [2,] 0.4809292 0.00000000
#> [3,] 0.5442232 0.06329396
#> [4,] 0.6075172 0.00000000
#> [5,]       Inf 0.00000000
```

The length of this list is the number of levels of the landscape.

An alternative, approximate construction computes the value of each
level of the landscape at each point on a 1-dimensional grid, ranging
from `min_x` to `max_x` at increments of `by`. A landscape constructed
using a discrete approximation is stored as a 3-dimensional array of
dimensions (levels, values, 2), with one level per feature (some of
which may be trivial) and one value per grid point, stored as $x,y$
pairs along the third dimension.

``` r
b_ran <- pl_support(pl1)
pl1d <- landscape(pd, degree = 1,
                  min_x = b_ran[[1L]], max_x = b_ran[[2L]], by = 0.02)
print(dim(pl1d$getInternal()))
#> [1]  4 21  2
print(pl1d$getInternal())
#> , , 1
#> 
#>           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
#> [1,] 0.2251884 0.2451884 0.2651884 0.2851884 0.3051884 0.3251884 0.3451884
#> [2,] 0.2251884 0.2451884 0.2651884 0.2851884 0.3051884 0.3251884 0.3451884
#> [3,] 0.2251884 0.2451884 0.2651884 0.2851884 0.3051884 0.3251884 0.3451884
#> [4,] 0.2251884 0.2451884 0.2651884 0.2851884 0.3051884 0.3251884 0.3451884
#>           [,8]      [,9]     [,10]     [,11]     [,12]     [,13]     [,14]
#> [1,] 0.3651884 0.3851884 0.4051884 0.4251884 0.4451884 0.4651884 0.4851884
#> [2,] 0.3651884 0.3851884 0.4051884 0.4251884 0.4451884 0.4651884 0.4851884
#> [3,] 0.3651884 0.3851884 0.4051884 0.4251884 0.4451884 0.4651884 0.4851884
#> [4,] 0.3651884 0.3851884 0.4051884 0.4251884 0.4451884 0.4651884 0.4851884
#>          [,15]     [,16]     [,17]     [,18]     [,19]     [,20]     [,21]
#> [1,] 0.5051884 0.5251884 0.5451884 0.5651884 0.5851884 0.6051884 0.6251884
#> [2,] 0.5051884 0.5251884 0.5451884 0.5651884 0.5851884 0.6051884 0.6251884
#> [3,] 0.5051884 0.5251884 0.5451884 0.5651884 0.5851884 0.6051884 0.6251884
#> [4,] 0.5051884 0.5251884 0.5451884 0.5651884 0.5851884 0.6051884 0.6251884
#> 
#> , , 2
#> 
#>      [,1] [,2]        [,3] [,4]        [,5]       [,6]       [,7]       [,8]
#> [1,]    0    0 0.007603077    0 0.003565016 0.02356502 0.04356502 0.06356502
#> [2,]    0    0 0.000000000    0 0.000000000 0.00000000 0.00000000 0.00000000
#> [3,]    0    0 0.000000000    0 0.000000000 0.00000000 0.00000000 0.00000000
#> [4,]    0    0 0.000000000    0 0.000000000 0.00000000 0.00000000 0.00000000
#>            [,9]    [,10]    [,11]    [,12]     [,13]       [,14]      [,15]
#> [1,] 0.08356502 0.103565 0.123565 0.143565 0.1423288 0.122328755 0.10232875
#> [2,] 0.00000000 0.000000 0.000000 0.000000 0.0000000 0.004259174 0.02425917
#> [3,] 0.00000000 0.000000 0.000000 0.000000 0.0000000 0.000000000 0.00000000
#> [4,] 0.00000000 0.000000 0.000000 0.000000 0.0000000 0.000000000 0.00000000
#>           [,16]      [,17]      [,18]      [,19]       [,20]      [,21]
#> [1,] 0.08232875 0.06425917 0.07063412 0.05063412 0.030634119 0.01063412
#> [2,] 0.04425917 0.06232875 0.04232875 0.02232875 0.002328755 0.00000000
#> [3,] 0.00000000 0.00000000 0.00000000 0.00000000 0.000000000 0.00000000
#> [4,] 0.00000000 0.00000000 0.00000000 0.00000000 0.000000000 0.00000000
```

Exactly computed landscapes can be converted to discrete landscape
objects, but the other direction is not well-defined. Below, we view a
portion of the discretized exact landscape and make a failed attempt to
exactify the discrete one:

``` r
# default conversion to discrete uses `by = 0.001`
# TODO: Export method with `min_x,max_x,by` parameters.
print(dim(pl1$getDiscrete()))
#> [1]   2 636   2
# print first 12 x-coordinates
pl1$getDiscrete()[, seq(220L, 280L), , drop = FALSE]
#> , , 1
#> 
#>      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#> [1,] 0.22 0.221 0.222 0.223 0.224 0.225 0.226 0.227 0.228 0.229  0.23 0.231
#> [2,] 0.22 0.221 0.222 0.223 0.224 0.225 0.226 0.227 0.228 0.229  0.23 0.231
#>      [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
#> [1,] 0.232 0.233 0.234 0.235 0.236 0.237 0.238 0.239  0.24 0.241 0.242 0.243
#> [2,] 0.232 0.233 0.234 0.235 0.236 0.237 0.238 0.239  0.24 0.241 0.242 0.243
#>      [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32] [,33] [,34] [,35] [,36]
#> [1,] 0.244 0.245 0.246 0.247 0.248 0.249  0.25 0.251 0.252 0.253 0.254 0.255
#> [2,] 0.244 0.245 0.246 0.247 0.248 0.249  0.25 0.251 0.252 0.253 0.254 0.255
#>      [,37] [,38] [,39] [,40] [,41] [,42] [,43] [,44] [,45] [,46] [,47] [,48]
#> [1,] 0.256 0.257 0.258 0.259  0.26 0.261 0.262 0.263 0.264 0.265 0.266 0.267
#> [2,] 0.256 0.257 0.258 0.259  0.26 0.261 0.262 0.263 0.264 0.265 0.266 0.267
#>      [,49] [,50] [,51] [,52] [,53] [,54] [,55] [,56] [,57] [,58] [,59] [,60]
#> [1,] 0.268 0.269  0.27 0.271 0.272 0.273 0.274 0.275 0.276 0.277 0.278 0.279
#> [2,] 0.268 0.269  0.27 0.271 0.272 0.273 0.274 0.275 0.276 0.277 0.278 0.279
#>      [,61]
#> [1,]  0.28
#> [2,]  0.28
#> 
#> , , 2
#> 
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7]  [,8]  [,9] [,10] [,11]  [,12]  [,13]
#> [1,]    0    0    0    0    0    0    0 0.001 0.002 0.001     0 -0.001 -0.001
#> [2,]    0    0    0    0    0    0    0 0.000 0.000 0.000     0  0.000  0.000
#>       [,14]  [,15]  [,16]  [,17]  [,18]  [,19]  [,20]  [,21]  [,22]  [,23]
#> [1,] -0.001 -0.001 -0.001 -0.001 -0.001 -0.001 -0.001 -0.001 -0.001 -0.001
#> [2,]  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000
#>       [,24]  [,25]  [,26]  [,27]  [,28]  [,29]  [,30]  [,31]  [,32]
#> [1,] -0.001 -0.001 -0.001 -0.001 -0.001 -0.001 -0.001 -0.001 -0.001
#> [2,]  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000  0.000
#>              [,33] [,34] [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42] [,43]
#> [1,] -2.385245e-18 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009  0.01
#> [2,]  0.000000e+00 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000  0.00
#>      [,44] [,45] [,46] [,47] [,48] [,49] [,50] [,51] [,52]         [,53]  [,54]
#> [1,] 0.009 0.008 0.007 0.006 0.005 0.004 0.003 0.002 0.001 -2.168404e-18 -0.001
#> [2,] 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000  0.000000e+00  0.000
#>       [,55]  [,56]  [,57]  [,58]  [,59]  [,60]  [,61]
#> [1,] -0.001 -0.001 -0.001 -0.001 -0.001 -0.001 -0.001
#> [2,]  0.000  0.000  0.000  0.000  0.000  0.000  0.000
try(pl1d$getExact())
#> Error in pl1d$getExact() : 
#>   Error: Can not convert a discrete PL to an exact PL.
```

### Visualization

**plt** provides a `plot()` method for the ‘Rcpp_PersistenceLandscape’
class. It uses **grDevices** to build color palettes, and as such its
default palette is viridis; but the user may supply the name of a
recognized palette or a sequence of colors between which to interpolate:

``` r
n_env <- max(pl_num_envelopes(pl1), pl_num_envelopes(pl1d))
par(mfrow = c(2L, 1L), mar = c(2, 2, 0, 2))
plot(pl1, palette = "terrain", n_levels = n_env, asp = 1)
plot(pl1d, palette = "terrain", n_levels = n_env, asp = 1)
```

![](man/figures/README-unnamed-chunk-12-1.png)<!-- -->

### Vector Operations

### Inner Product
