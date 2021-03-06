# RSPO Leakage

#### This repository contains the `R` project set up with the `packrat` package for code used for the analysis in the paper *'Deforestation spillovers from oil palm sustainability certification' (https://doi.org/10.1088/1748-9326/ab7f0c)*. This project has a private project library that installs version of libraries and any dependencies used in the code for this project. In addition, the repo contains all data necessary to replicate the paper's analyses.


Instructions on how to download and use the code
----------------------------------------------------
1. Download or clone repository to local machine. If the location of the project folder is in a location that is being backed up by a cloud storage service like `Dropbox` or `Microsoft OneDrive`, pause the syncing process to prevent any issues with the package installation process.

2. Extract files in `data.zip`(`long_kali.csv`, `master_mill_data.csv` and `long_kali_metadata.csv`) in the `input` folder into the same folder. 

3. Open `R` project file `rspo_leakage.Rproj` in RStudio. Once the project opens, `packrat` will automatically install into the project private library. The `R` session will restart when this is completed.

4. Run `packrat::status()` in the console to see list of packages used in the code that need to be installed for you to run the code.

5. Run `packrat::restore()` in the console to install all required packages and dependencies. If required, a window will pop up requesting `RTools` to be installed. Install RTools and then re-run `packrat::restore()` to finishing installing any other packages and dependencies.

6. Once all packages have been installed, open `leakage.R` and run the code. Figures, plot and tables will be saved into the `output` folder.




Additional information and resources on `packrat` 
-------------------------------------------------

[1] https://rstudio.github.io/packrat/

[2] http://rstudio-pubs-static.s3.amazonaws.com/221948_fb7215fecb0d49ac903f701fd8d45132.html

[3] https://blog.methodsconsultants.com/posts/using-packrat-with-git-submodules/




 
