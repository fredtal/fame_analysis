# ------------------ #
# Convert LabSolutions text files to clean csvs
# Talia Young
# October 2015
# ------------------ #
# This code does three things: 
# 1. Converts Lab Solutions text files to csv files (tab-delimited)
# 2. Converts those tab-delimited files to comma-delimited files
# 3. Removes header junk from those comma-delimited files
# ------------------ #
#rm(list=ls())
library(stringr)

## ** NOTE: each filename must be in the following format: 
# uniqueidnumber_species_otherinfo
#  * id number must be unique in the data set and cannot have spaces or underscores
#  * species can have spaces or dashes but no underscores
#  * otherinfo can have spaces or dashes but no underscores

#  * there must be 3 pieces of data separated between 2 underscores
#  * if there is no species, it should say "id_unknown".  
#  * blanks should be: blank15_blank_otherinfo

# ------------ #
# 1. Set up directories ----
# ------------ #
# You need to set up 5 directories, one with the original files, and 4 empty ones: 
# importdir: with the original text files
importdir <- "~/Documents/Rutgers/fatty acids/Jen's sample data/1.txt files (original)"

# tempdir1: empty (this will hold the tab-delimited files converted from text)
tempdir1 <- "~/Documents/Rutgers/fatty acids/Jen's sample data/2.csv files (tab-delimited, with header)"

# tempdir2: empty (this will hold the comma-delimited files converted from 
#                  tempdir1 with all the header junk still at the top)
tempdir2 <- "~/Documents/Rutgers/fatty acids/Jen's sample data/3.csv files (comma-delimited, with header)"

#tempdir3: empty (this will hold csv files with no header but all columns
tempdir3 <- "~/Documents/Rutgers/fatty acids/Jen's sample data/4.csv files (no header, all columns)"

#exportdir: empty (this will hold csv files with no header but all columns
exportdir <- "~/Documents/Rutgers/fatty acids/Jen's sample data/5.csv files (rt, ap only)"

# ------------ #
# 2. Do the conversions ----
# ------------ #

# ------------ #
# a. Make the list of files you want to import
# ------------ #
filelist <- list.files(importdir)

# ------------ #
# b. Convert txt files to csv files (tab-delimited)
# ------------ #
for(i in 1:length(filelist)) {
  print(paste("Converting '", filelist[i], "'...", sep=""))
  
  #import file
  file <- readLines(paste(importdir, filelist[i], sep="/"))
  
  #parse filename
  pieces <- strsplit(as.character(filelist[i]), ".txt") #split on period
  
  #export the file as a csv
  writeLines(text=file, 
             con=paste(tempdir1, "/", sapply(pieces, "[", 1), ".csv", 
                       sep="")) 
  
  rm(file, pieces)
}
rm(i, filelist)

# ------------ #
# c. Convert csv files (tab-delimited) to csv files (comma-delimited)
# ------------ #
filelist <- list.files(tempdir1)

#for each file
for(i in 1:length(filelist)) {
  print(paste("Importing '", filelist[i], "'...", sep=""))
  #import file
  file <- readLines(paste(tempdir1, filelist[i], sep="/"))
  
  #replace tabs with commas
  newfile <- gsub(pattern="\\t", replacement=",", x=file)
  
  #reparse the filename
  pieces <- strsplit(as.character(filelist[i]), ".csv") #split on period
  
  #export the file as a comma-delimited text file
  writeLines(text=newfile, con=paste(tempdir2, "/", sapply(pieces, "[", 1), ".csv", sep="")) 
  
  rm(file, newfile, pieces)
}
rm(i, filelist)

# ------------ #
# c. Remove header junk from the top of the csv file
# ------------ #

filelist <- list.files(tempdir2)

#for each file
for(i in 1:length(filelist)) {
  print(paste("Converting '", filelist[i], "'...", sep=""))
  
  #import file
  file <- readLines(paste(tempdir2, filelist[i], sep="/"))
  
  #remove header
  #startrow <- grep("\\[Original", file) #the useful part starts with the line that has "[Original]" on it
  startrow <- grep("Peak#", file) #the useful part starts with the line that has "[Original]" on it
  #endrow <- grep("\\[Compound", file)
  #newfile <- file[startrow:(endrow-2)]
  newfile <- file[startrow:(length(file)-1)]
  
  #parse filename
  pieces <- strsplit(as.character(filelist[i]), ".csv") #split on period
  #export the file as a csv
  writeLines(text=newfile, 
             con=paste(tempdir3, "/", sapply(pieces, "[", 1), "_noheader.csv", 
                       sep=""))   
  rm(file, newfile, startrow, pieces)
}

rm(i, filelist)

# ------------ #
# d. Extract retention time and area percent columns only
# ------------ #
filelist <- list.files(tempdir3)

for(i in 1:length(filelist)) {
  print(paste("Processing", filelist[i]))
  
  thisfilewhole <- read.csv(paste(tempdir3, filelist[i], sep="/"), 
                            row.names=NULL, header=TRUE,
                            stringsAsFactors=FALSE)
  
  # b. Identify which columns to extract ----
  # area column
  areacolindices <- grep("Area", colnames(thisfilewhole)) #all the columns with area
  #columns with area in them but we don't want
  otherareacolindices <- c(grep("Area%", colnames(thisfilewhole)), 
                           grep("Area %", colnames(thisfilewhole)), 
                           grep("Area Percent", colnames(thisfilewhole)), 
                           grep("Area Ratio", colnames(thisfilewhole)),
                           grep("Area.Ratio", colnames(thisfilewhole))) 
  #all columns with area minus other columns with area = area column that I want
  areacolindex <- setdiff(areacolindices, otherareacolindices) 
  rm(areacolindices, otherareacolindices)
  
  # rt column
  rtcolindex <- grep("R.Time", colnames(thisfilewhole))
  
  # c. Pull out just desired columns (RT and AP) and rows (first peak to cutoff)----
  thisfile <- thisfilewhole[,c(rtcolindex, areacolindex)]
  rm(rtcolindex, areacolindex)
  
  # d. Clean up a little bit----
  colnames(thisfile) <- c("retention_time", "area")  #rename columns!
  thisfile <- apply(thisfile, c(1,2), as.numeric)  #convert all to numeric
  thisfile <- as.data.frame(thisfile) #convert to data frame in order to be able to call column names 
  
  #remove any NAs at the bottom, if there are any
  if(sum(is.na(thisfile$retention_time))>0) {
    nastartrow <- min(which(is.na(thisfile$retention_time)))
    thisfile <- thisfile[1:(nastartrow-1),]
  }  
  # e. Name and export the file----
  # get the fish id and species
  # first number = id.  text between underscores = species
  splitname <- strsplit(filelist[i], "[_.]+")[[1]] #split on underscore OR period
  #thank you: http://stackoverflow.com/questions/10738729/r-strsplit-with-multiple-unordered-split-arguments
  #also: http://stat.ethz.ch/R-manual/R-patched/library/base/html/regex.html
  
  thisfishid <- splitname[1] 
  #convert blank to all lowercase
  if(grepl("Blank", thisfishid) | grepl("blank", thisfishid)) thisfishid <- tolower(thisfishid)
  
  thisfishsp <- tolower(splitname[2])
  if(thisfishsp=="striped marlin") thisfishsp="stm"
    
  if(thisfishsp=="csv") { #if blank, there's no species
    newfilename <- thisfishid 
  } else {
    newfilename <- paste(thisfishid, thisfishsp, sep="_")
  }
  rm(splitname)
  
  write.csv(thisfile, 
            paste(exportdir, "/", newfilename, ".csv", sep=""), 
            row.names=FALSE)
  
  rm(thisfilewhole, thisfile, thisfishid, thisfishsp, newfilename)
}
rm(i, filelist)


rm(importdir, tempdir1, tempdir2, tempdir3, exportdir)
