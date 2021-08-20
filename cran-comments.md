## Test environments
* local Fedora install, R 4.0.5
* win-builder (devel and release)
* rhub

## R CMD check results
There were no ERRORs, NOTEs or WARNINGs.

This is a resubmission, with the following corrections as kindly suggested by Julia Haider:

* references to the methods have been added to the description field of the DESCRIPTION file
* all instances of cat() were replaced by message()
* the original par settings are now saved if the function is exited, using oldpar <- par(no.readonly=TRUE); on.exit(par(oldpar))
