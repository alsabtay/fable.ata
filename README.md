# fable.ata
Fable Modelling Wrappers for ATAforecasting Package

This package provides a tidy R interface to the Ata method
procedure using [fable](https://github.com/tidyverts/fable). This
package makes use of the [ATAforecasting
package](https://cran.r-project.org/package=ATAforecasting) for R.

## Installation


You can install the **development** version from
[Github](https://github.com/alsabtay/fable.ata) with:

``` r
# install.packages("remotes")
remotes::install_github("alsabtay/fable.ata")
```

## Example


USAccDeaths: Accidental Deaths in the US 1973--1978

``` r
library(fable.ata)
as_tsibble(USAccDeaths) %>% model(ata = ATAM(value ~ trend("A") + season("A", method = "stl"))) %>% forecast(h=24)
``` 
