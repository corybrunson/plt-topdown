
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
pd <- ripserr::vietoris_rips(pc, max_dim = 1, threshold = 2, p = 2)
print(pd)
#> PHom object containing persistence data for 63 features.
#> 
#> Contains:
#> * 59 0-dim features
#> * 4 1-dim features
#> 
#> Radius/diameter: min = 0; max = 0.63582.
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
objects, but the other direction is not well-defined:

``` r
# default conversion to discrete uses `by = 0.001`
# TODO: Export method with `min_x,max_x,by` parameters.
print(dim(pl1$getDiscrete()))
#> [1]   2 636   2
# print first 12 x-coordinates
pl1$getDiscrete()[, seq(12L), , drop = FALSE]
#> , , 1
#> 
#>       [,1]  [,2]  [,3]  [,4]  [,5]  [,6]  [,7]  [,8]  [,9] [,10] [,11] [,12]
#> [1,] 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009  0.01 0.011 0.012
#> [2,] 0.001 0.002 0.003 0.004 0.005 0.006 0.007 0.008 0.009  0.01 0.011 0.012
#> 
#> , , 2
#> 
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
#> [1,]    0    0    0    0    0    0    0    0    0     0     0     0
#> [2,]    0    0    0    0    0    0    0    0    0     0     0     0
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
