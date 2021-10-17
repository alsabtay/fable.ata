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
as_tsibble(USAccDeaths) %>% model(ata = AutoATA(value ~ trend("A") + season("A", method = "stl"))) %>% forecast(h=24)
``` 

## Links

[Github page](https://github.com/alsabtay/fable.ata)

[Github.io page](https://alsabtay.github.io/fable.ata/index.html)

[Github - ATAforecasting](https://github.com/alsabtay/ATAforecasting)

[Github.io - ATAforecasting](https://alsabtay.github.io/ATAforecasting/)

[Project team website](https://atamethod.wordpress.com/)


## License
This package is free and open source software, licensed under GPL-3.
