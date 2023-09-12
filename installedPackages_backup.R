
# Save a list of packages installed in my current Rstudio/ computer device
installedPackages <- as.data.frame(installed.packages())
write.csv(installedPackages, "202330912_previously_installed_packages.csv")
getwd()

# Create a list of libraries from my old list that are not already installed when I fresshly downloaded R 
# from my new device 

installedPreviously <- read.csv("/Library/Frameworks/R.framework/Versions/4.1/Resources/bin/exec/Rscripts/202330912_previously_installed_packages.csv")

baseR <- as.data.frame(installed.packages())

toInstall <- setdiff(installedPreviously, baseR)
 
# OR ALTERNATIVELY if previous line doesn't work properly

toInstall <- setdiff(installedPreviously$X, baseR$Package)

# Download the list of libraries

install.packages(toInstall)