# Modified IPOD function

This repository contains a modified version of the `IPOD()` function from the R package `leapp`.

The original function has a `length.out` argument, but this argument is not used when constructing the sequence of lambda values. Instead, the sequence is generated using a fixed step size of 0.1.

For data on a large numerical scale, this can result in a very large number of lambda values and consequently a long computation time.

In this modified version, the lambda sequence is generated using

```r
length.out = length.out
```

so that the number of lambda values can be controlled by the `length.out` argument.

## Usage

The function can be loaded directly from GitHub using `source()`:

```r
source("https://raw.githubusercontent.com/fsakaori-ds/IPOD-modified/main/IPOD_modified.R")
```

The default value is `length.out = 50`. For example,

```r
result <- IPOD(X, Y, H, length.out = 50)
```

## Modification

The original code constructs the lambda sequence using

```r
lambdas = seq(..., 0, by = -0.1)
```

This version changes it to

```r
lambdas = seq(..., 0, length.out = length.out)
```

in both branches of the function.

No other changes have been made to the original `IPOD()` function.

## Reference

She, Y. and Owen, A. B. (2011).  
Outlier Detection Using Nonconvex Penalized Regression.  
*Journal of the American Statistical Association*, 106, 626–639.

## License

This repository contains modified code from the R package `leapp`, which is distributed under GPL (>= 2).

This modified version is distributed under the GNU General Public License version 3 or later.
