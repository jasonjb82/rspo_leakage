## RSPO Leakage

#### This repository contains the `R` project set up with the `packrat` package for code used for the analysis in the paper *'Deforestation spillovers from oil palm sustainability certification' (DOI)*. This project has a private project library that installs version of libraries and any dependencies used in the code for this project.


**Instructions on how to download and use the code**
1. Download or clone repository to local machine

2. Open `R` project file `rspo_leakage.Rproj` in RStudio. Once the project opens, `packrat` will automatically be installed into the project private library. The `R` session will restart when this is completed.

3. Run `packrat::status()` in the console to see list of packages used in the code that need to be installed for you to run the code.

4. Run `packrat::restore()` in the console to all required packages install packages and dependencies. If required, a window will pop up requesting `RTools` to be installed. Install RTools and then re-run `packrat::restore()` to finishing installing any other packages and dependencies.

5. Once all packages have been installed, open `rspo_leakage.R` and run the code.




#### Additional information and resources on `packrat` at the links below

[1] https://rstudio.github.io/packrat/

[2] http://rstudio-pubs-static.s3.amazonaws.com/221948_fb7215fecb0d49ac903f701fd8d45132.html

[3] https://blog.methodsconsultants.com/posts/using-packrat-with-git-submodules/




 
