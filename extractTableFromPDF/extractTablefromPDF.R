
---
  title: "extract tables"
author: "Ming Tang"
date: "6/9/2023"
output: html_document
---
  
  ### how to extract tables from PDF file using Tabulizer https://github.com/ropensci/tabulizer
  
  
  #install.packages("tabulizer") not working, have to install it from github
  if (!require("remotes")) {
    install.packages("remotes")
  }
remotes::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"))
library(tabulizer)

out <- extract_tables("/Users/nikolasgiannakis/Desktop/Papers/Patsalos_et_al.pdf", pages = 11, guess = TRUE, 
                      output = "data.frame")



out[[1]]
View(out[[1]])
