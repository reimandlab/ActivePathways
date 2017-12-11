# activeDriverPW 0.0.0.9002

## Major changes
* Rename to activeDriverPW

# mpea 0.0.0.9001

## Major Changes
* Change column contribution method. Column contribution is now reported as the
log-fold-change when the column is excluded. ie, -log10(p_val_with_column / p_val_without_column)

## Minor Changes
* No longer raises an error if no significant terms are found and `return.all` 
is false. Instead returns an empty `data.table` and issues a warning.
* Issue a warning if genes are filtered out for not being found in the background
* Change default p-value adjustment method to "holm"
* If return.all==FALSE, calculates columnContribution only for terms that will
returned. Speeds up runtime by roughly a factor of 2
* Add "none" option to p-value adjustment methods
* New implementation of orderedHypergeometric function which is several times faster.

## Bug fixes
* Fix Brown's method when all p-values in a column are the same
* Fix in orderedHypergeometric which added an extra NA to the `complement`, leading
to small errors in the calculated p-value
* Fixed another bug in orderedHypergeometric which added an extra NA to the 
`complement` in some cases, leading to small errors in the calculated p-value

# mpea 0.0.0.9000

* Initial Build

* Added a `NEWS.md` file to track changes to the package.
