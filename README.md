# coffee

## Chronological Ordering for Fossils and Environmental Events

It's quite fortunate that the above title can be acronymised to *coffee*. While individual calibrated radiocarbon dates can span several centuries, combining multiple dates together with any chronological constraints can make a chronology much more robust and precise.

This package uses Bayesian methods to enforce the chronological ordering of radiocarbon and other dates, 
for example for trees with multiple radiocarbon dates spaced at exactly known intervals (e.g., 10 annual rings). 
Another example is sites where the relative chronological position of the dates is taken into account - 
the ages of dates further down a site must be older than those of dates further up.

Please check out the vignettes folder for a tutorial. In short, install the package from github, load it and run the two main functions:

```{r, eval=FALSE}
require(devtools)
install_github("Maarten14C/coffee")
library(coffee)
rings()
strat()
```
