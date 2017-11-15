# mpea 0.0.0.90001

## Minor Changes
* No longer raises an error if no significant terms are found and `return.all` 
is false. Instead returns an empty `data.table` and issues a warning.
* Issue a warning if genes are filtered out for not being found in the background
* Change default p-value adjustment method to "holm"
* If return.all==FALSE, calculates columnContribution only for terms that will
returned. Speeds up runtime by roughly a factor of 2
* Add "none" option to p-value adjustment methods


## Bug fixes
* Fix Brown's method when all p-values in a column are the same
* Fix in orderedHypergeometric which added an extra NA to the `complement`, leading
to small errors in the calculated p-value

# mpea 0.0.0.9000

* Initial Build

* Added a `NEWS.md` file to track changes to the package.



