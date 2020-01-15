## RSPO Leakage
#### This repository contains the RStudio project set up with `packrat` for code used on 'Deforestation spillovers from oil palm sustainability certification' (DOI). This project has a private project library that installs version of libraries used in setting up the code for this project and any dependencies.

**Instructions on how to download and use the code**
1. Download or clone repository to local machine
2. Open RStudio project file rspo_leakage.rnw
3. Run `packrat::status()` in the console to see packages that need to be installed
4. Run `packrat::restore()` to install packages and dependencies. If required, a window will popup requesting `RTools` to be installed. Install RTools and then re-run `packrat::restore()`
5. Once all packages have been installed, you should be able to run the code with all the required libraries.

Additional information and resources on `packrat` at the links below

[1] https://rstudio.github.io/packrat/

[2] http://rstudio-pubs-static.s3.amazonaws.com/221948_fb7215fecb0d49ac903f701fd8d45132.html

[3] https://blog.methodsconsultants.com/posts/using-packrat-with-git-submodules/




 
