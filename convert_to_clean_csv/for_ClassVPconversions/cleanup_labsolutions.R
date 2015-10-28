# ------------------ #
# Clean up LabSolutions .csv files
# (remove large header)
# Talia Young
# May 2015
# ------------------ #

#rm(list=ls())
#require(stringr)

# text files ----

importdir <- "~/Documents/Rutgers/fatty acids/Jen's sample data/original data/csv files (converted, messy)"
exportdir <- "~/Documents/Rutgers/fatty acids/Jen's sample data/original data/csv files (clean)"

#make the list of files you want to import
filelist <- list.files(importdir)

#for each file
for(i in 1:length(filelist)) {
  print(paste("Converting '", filelist[i], "'...", sep=""))
  
  #import file
  file <- readLines(paste(importdir, filelist[i], sep="/"))
  
  #remove header
  #startrow <- grep("\\[Original", file) #the useful part starts with the line that has "[Original]" on it
  startrow <- grep("Peak#", file) #the useful part starts with the line that has "[Original]" on it
  endrow <- grep("\\[Compound", file)
  newfile <- file[startrow:(endrow-2)]
  
  #parse filename
  pieces <- strsplit(as.character(filelist[i]), ".csv") #split on period
  #export the file as a csv
  writeLines(text=newfile, 
             con=paste(exportdir, "/", sapply(pieces, "[", 1), "_noheader.csv", 
                       sep=""))   
  rm(file, pieces)
}

rm(importdir, exportdir, filelist)


