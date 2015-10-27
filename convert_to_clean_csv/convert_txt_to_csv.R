# ------------------ #
# Convert GC txt files to csv format
# Talia Young
# April 2015
# ------------------ #

#rm(list=ls())
require(stringr)

# text files ----

importdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay shared/Magdalena Bay data/GC data and reports/GC data-cleaned/reports-all, organized by file type/txt files"
exportdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay shared/Magdalena Bay data/GC data and reports/GC data-cleaned/reports-all, organized by file type/txt-to-csv"

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

rm(importdir, exportdir, filelist)


