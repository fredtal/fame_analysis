# ------------------ #
# Convert GC are files to csv format
# Talia Young
# May 2015
# ------------------ #

#rm(list=ls())
#require(stringr)

# text files ----

importdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay shared/Magdalena Bay data/GC data and reports/GC data-cleaned/reports-all, organized by file type/are files"

# 1. convert ares to txt (they emerge as tab-delimited)

filelist <- list.files(importdir)

for(i in 1:length(filelist)) {
  print(paste("Converting '", filelist[i], "'...", sep=""))
  
  if(grepl(".are", filelist[i])) {
    newfilename <- sub(".are", ".txt", filelist[i])
    file.rename(from=paste(importdir, filelist[i], sep="/"), 
                to=paste(importdir, newfilename, sep="/"))
  }
}  
rm(newfilename, filelist)

