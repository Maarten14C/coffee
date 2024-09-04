# coffee 0.4.0
* now linking to both the 'rice' and 'rintcal' packages
* updated how gaps in strats are plotted

# coffee 0.3.0
* updated the vignette
* added functions scissors and thinner for manipulation of MCMC output
* enhanced strat code so it runs around 40% faster
* previous strat runs can now be reloaded (using run=FALSE)
* files are now read and written to faster (using fread, fwrite)
* strat iterations can either be loaded into memory or written to files
* added fit measure (and also offset for rings)

# coffee 0.2.0
* now links to the rintcal package (renamed from IntCal)
* the strat function now handles undated levels, blocks of dates, and different types of gaps (of exactly known length, or according to a normal or gamma distribution)

# coffee 0.1.1

* added references to the methods to the DESCRIPTION file
* replaced all instances of cat() with message()
* the original par settings are now saved if the function is exited

# coffee  0.1.0

* this is the first CRAN iteration of the coffee package

