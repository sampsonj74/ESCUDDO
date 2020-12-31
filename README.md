# ESCUDDO

This library contains the functions for performing the ESCUDDO analysis  

This library can be downloaded using devtools::install_github("sampsonj74/ESCUDDO")  

The vignette can be found by clicking on http://htmlpreview.github.io/?https://github.com/sampsonj74/Vignettes/blob/main/escuddo.html

##########################################################################

I need to remember what I need to do when building future copies:

library("devtools","testthat","roxygen2")

devtools::document("/Users/sampsonjn/Documents/WORK/ESCUDDO/ESCUDDO")  

devtools::build_vignettes("/Users/sampsonjn/Documents/WORK/ESCUDDO/ESCUDDO")  

##############################################################################

I originally handled the vignette using the following code:

file.copy("/Users/sampsonjn/Documents/WORK/ESCUDDO/ESCUDDO/doc/my-vignette.html",
          "/Users/sampsonjn/Documents/WORK/ESCUDDO/ESCUDDO/vignettes/my-vignette_cache/html/ESCUDDOv.html", overwrite=TRUE)  
          
git add -f /Users/sampsonjn/Documents/WORK/ESCUDDO/ESCUDDO/vignettes/my-vignette_cache/html/ESCUDDOv.html

I now just save the file to the Vignettes repositiory:

All I need to do is rename it as ESCUDDO
