# Predict from a connectivity model

Projects new feature values (for example future climate) onto the PCA
loadings fitted to present-day data, so the components keep their
meaning, then predicts with the Cubist model.

## Usage

``` r
# S3 method for class 'corridor_model'
predict(object, newdata, ...)
```

## Arguments

- object:

  A `corridor_model` from
  [`fit_connectivity()`](https://paulmaier.github.io/corridoR/reference/fit_connectivity.md).

- newdata:

  A data.frame with the same feature columns as the data the model was
  fit to.

- ...:

  Not used.

## Value

Numeric predictions.
