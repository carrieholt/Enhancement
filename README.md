# Enhancement
Simulation model of genetic fitness impacts of enhancement for Chinook salmon

Author: Carrie Holt
Developed: 2016/2017 
Updated: Feb 2023
Results published: Withler et al (2018) CSAS Res. Doc.

The main function that runs and plots the code is ContourSA_fitness.R. The underlying population model is described in the function run.lever.model is within the file "run.level.model.r". 

In Feb 2023, the code was updated to account for two corrections:
(1) The equation for numbers of hatchery smolts (Hatch.sm in FuncDefs.R) was changed to include the proportion of the brood stock that is female. Previously this proportion was not included and all brood was assumed to be female.
(2) The equation for the optimum phenotype value in the natural environment (Pnat in run.lever.model.r) was changed so that it used pHOSeff instead of pHOS.


The akima package v0.6-2 needs to be installed for contour plots*. The updated versions of this package do not provide smooth edges on contour plots. 


*I installed the older package this way:
detach("package:akima", unload=TRUE)
install.packages("C:/Users/HoltC/AppData/Local/Programs/R/R-4.2.1/library/akima_0.6-2.tar.gz", repos = NULL, type="source")