## Test environments
* local OS X install, R 4.0.0
* Win-builder
* Travis-CI / Github: 
	* oldrel: passed
	* release: passed
	* devel: failed, because a dependency was not available on R-devel (tcltk). But I think this is rather a problem of Travis.

## R CMD check results
There were no ERRORs, WARNINGs, or relevant NOTEs on oldrel and release.

## Downstream dependencies
None.