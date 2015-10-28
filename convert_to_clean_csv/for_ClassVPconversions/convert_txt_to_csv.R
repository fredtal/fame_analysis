# ------------------ #
# Convert GC txt files to csv format
# Talia Young
# April 2015
# ------------------ #

#rm(list=ls())
library(stringr)

# text files ----

importdir <- "~/Documents/Rutgers/fatty acids/Jen's sample data/txt files (original)"
exportdir <- "~/Documents/Rutgers/fatty acids/Jen's sample data/csv files (converted, messy)"

#make the list of files you want to import
filelist <- list.files(importdir)

#for each file
for(i in 1:length(filelist)) {
  print(paste("Converting '", filelist[i], "'...", sep=""))
  
  #import file
  file <- readLines(paste(importdir, filelist[i], sep="/"))

  #parse filename
  pieces <- strsplit(as.character(filelist[i]), ".txt") #split on period
  
  #export the file as a csv
  writeLines(text=file, 
             con=paste(exportdir, "/", sapply(pieces, "[", 1), ".csv", 
                       sep="")) 
  
  rm(file, pieces)
}

rm(importdir, exportdir, filelist, i)


