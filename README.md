# ESCUDDO

This library contains the functions for performing the ESCUDDO analysis  

This library can be downloaded using devtools::install_github("sampsonj74/ESCUDDO")  

The vignette is buried in vignettes/my-vignette_cache/html/ESCUDDOv.html  

##########################################################################

I need to remember what I need to do when building future copies:

library("devtools","testthat","roxygen2")

devtools::document("/Users/sampsonjn/Documents/WORK/ESCUDDO/ESCUDDO")  

devtools::build_vignettes("/Users/sampsonjn/Documents/WORK/ESCUDDO/ESCUDDO")  

file.copy("/Users/sampsonjn/Documents/WORK/ESCUDDO/ESCUDDO/doc/my-vignette.html",
          "/Users/sampsonjn/Documents/WORK/ESCUDDO/ESCUDDO/vignettes/my-vignette_cache/html/ESCUDDOv.html", overwrite=TRUE)  
          
git add -f /Users/sampsonjn/Documents/WORK/ESCUDDO/ESCUDDO/vignettes/my-vignette_cache/html/ESCUDDOv.html


