# Enhancement
Simulation model of genetic fitness impacts of enhancement for Chinook salmon

Author: Carrie Holt
Developed: 2016/2017 
Updated: Feb 2023
Results published: Withler et al (2018) CSAS Res. Doc.

The main function that runs and plots the code is ContourSA_fitness.R. The underlying population model is described in the function run.lever.model is within the file "run.level.model.r". 

The akima package v0.6-2 needs to be installed for contour plots. The new version does not provide smooth edges on contours

I installed the older package this way:
detach("package:akima", unload=TRUE)
install.packages("C:/Users/HoltC/AppData/Local/Programs/R/R-4.2.1/library/akima_0.6-2.tar.gz", repos = NULL, type="source")