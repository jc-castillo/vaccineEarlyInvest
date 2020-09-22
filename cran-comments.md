## cran-comments for vaccineEarlyInvest v0.1.0

Test environments:

* Ubuntu 16.04.6 LTS x86_64-pc-linux-gnu (64-bit) R 4.0.2 (2020-06-22) via Travis
* local OS X install, R 4.0.2
* builder.r-hub.io (platforms = NULL)
* win-builder (devel and release)

There are no WARNINGs or ERRORs and 4 NOTEs:

N  checking CRAN incoming feasibility (932ms)      New submission   Maintainer: 'Juan Camilo Castillo <jccast@upenn.edu>'

N  checking for future file timestamps   unable to verify current time

N  checking files in 'vignettes'     'figure'   The following directory looks like a leftover from 'knitr':   Please remove from your package.   

N checking examples ... NOTE
Examples with CPU or elapsed time > 5s
                 user system elapsed
candidatesFung 13.237  0.113  28.535


* The vignettes/figure directory has precompiled figures for the vignettes. This is necessary because a full run of the code in the vignettes takes around 10 minutes to run locally.