# Reference manual

```@meta
CurrentModule = SparseElementaryDivisors
```

```@docs
SparseElementaryDivisors
```

```@index
Pages = ["reference_manual.md"]
```

## Elementary divisors

```@docs
sparse_elementary_divisors
```

## Settings and statistics

```@docs
ReductionConfig
ReductionStats
```

## Integer types and overflow

```@docs
SparseElementaryDivisors.DEFAULT_VALUE_TYPES
CoefficientOverflowError
```

## Matrix Market files

```@docs
read_matrix_market
write_matrix_market
```

## Comparison with Magma

```@docs
magma_elementary_divisors
timed_magma_elementary_divisors
```

## Developer tools

Public names that are not exported: they are reached as `SparseElementaryDivisors.ReductionTrace` and `SparseElementaryDivisors.RowsData`.

```@docs
SparseElementaryDivisors.ReductionTrace
SparseElementaryDivisors.RowsData
```
