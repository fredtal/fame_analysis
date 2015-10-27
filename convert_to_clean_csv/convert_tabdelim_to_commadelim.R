# ------------------ #
# Convert GC tab-delimited text files to comma-delimited text files
# Talia Young
# April 2015
# ------------------ #

# (I accidentally saved some of the Lab Solutions files as tab-delimited
# instead of comma-delimited.)
# (I converted the .are files to tab-delimited .txt files.)

#rm(list=ls())

importdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay shared/Magdalena Bay data/GC data and reports/GC data-cleaned/reports-all, organized by file type/are files"
exportdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay shared/Magdalena Bay data/GC data and reports/GC data-cleaned/reports-all, organized by file type/are-to-csv"

#make the list of files you want to import
filelist <- list.files(importdir)

#for each file
for(i in 1:length(filelist)) {
  print(paste("Importing '", filelist[i], "'...", sep=""))
  #import file
  file <- readLines(paste(importdir, filelist[i], sep="/"))

  #replace tabs with commas
  newfile <- gsub(pattern="\\t", replacement=",", x=file)
  
  #reparse the filename
  pieces <- strsplit(as.character(filelist[i]), ".txt") #split on period
    
  #export the file as a comma-delimited text file
  writeLines(text=newfile, con=paste(exportdir, "/", sapply(pieces, "[", 1), ".csv", sep="")) 
  #writeLines(text=newfile, con=paste(importdir, "/", filelist[i], "_commas.txt", sep="")) 
  
  rm(file, newfile)
}
rm(i)
