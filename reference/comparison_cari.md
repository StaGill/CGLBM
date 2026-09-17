# Co-clustering adjusted Rand index

Compare two row/column partitions without constructing cellwise
character labels.

## Usage

``` r
comparison_cari(z_a, w_a, z_b, w_b)
```

## Arguments

- z_a, w_a:

  Row and column labels for the first co-clustering.

- z_b, w_b:

  Row and column labels for the second co-clustering.

## Value

A numeric co-clustering adjusted Rand index.

## Examples

``` r
comparison_cari(c(1, 1, 2), c(1, 2), c(2, 2, 1), c(2, 1))
#> [1] 1
```
