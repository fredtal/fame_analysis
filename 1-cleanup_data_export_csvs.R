# ------------------ #
# 0. Import GC data in various forms, clean it up, and then export all as CSV files
# Talia Young
# 5/2015, based on "1-import GC data" from 2012
#
# This script cleans up the data files as they come direct from the GC.  
# You only need to run it once.
#
# Input: it imports everything in one directory
# Output: it spits out correctly formatted csv files in another directory. 
# ------------------ #

# ------------------ #
# 1. Set up
# ------------------ #

#rm(list=ls())

#this setting prevents warnings when searching for the cutoff text
Sys.setlocale('LC_ALL','C') 
#http://r.789695.n4.nabble.com/Strings-from-different-locale-td3023176.html

# set directories
importdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay shared/Magdalena Bay data/GC data and reports/GC data-cleaned/reports-all, organized by file type/csv files-ALL-original"
exportdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay shared/Magdalena Bay data/GC data and reports/GC data-cleaned/reports-all, organized by file type/csv files-ALL-trimmed"

# generate list of files
filelist <- list.files(importdir)

# ------------------ #
# 2. For each file ----
# ------------------ #
for(i in 1:length(filelist)) {
  print(paste("Processing", filelist[i]))

  # a. Import and clean up the file ----  
  thisfilewhole <- read.csv(paste(importdir, filelist[i], sep="/"), 
                            row.names=NULL, header=FALSE,
                            stringsAsFactors=FALSE)
    
  # b. Identify which columns and rows to extract ----
  # area column
  areacolindices <- grep("Area", thisfilewhole) #all the columns with area
  #columns with area in them but we don't want
  otherareacolindices <- c(grep("Area%", thisfilewhole), 
                           grep("Area %", thisfilewhole), 
                           grep("Area Percent", thisfilewhole), 
                           grep("Area Ratio", thisfilewhole)) 
  #all columns with area minus other columns with area = area column that I want
  areacolindex <- setdiff(areacolindices, otherareacolindices) 
  
  # rt column
  # check to see if any of the RT search strings are in the file
  rttext <- c("Retention Time", "R.Time", "Ret. Time")  
  for(j in 1:length(rttext)) { #for each of the possible rt column headers
    if(length(grep(rttext[j], thisfilewhole))>0) { #if you find it the file
      rtcolindex <- grep(rttext[j], thisfilewhole) #find the column index
      break #then break the for loop
    }
  }

  # data start row
  firstpkrowindex <- grep(rttext[j], thisfilewhole[,rtcolindex]) 
  #rm(rttext, j)
  
  # bottom cutoff row
  cutofftext <- c("Totals", "\\[Compound")  
  for(k in 1:length(cutofftext)) { #for each of the possible cutoff text options
    if(length(grep(cutofftext[k], thisfilewhole[,1]))>0) { #if you find it the file
      #the cutoff row is one above "Totals", so the grep row minus 1
      cutoffrow <- grep(cutofftext[k], thisfilewhole[,1]) - 1 #
      break #then break the for loop
    } else {
      cutoffrow <- nrow(thisfilewhole) #if the search string doesn't come up, it's the whole file
    }
  } #rm(k, cutofftext)
  
  # c. Pull out just desired columns (RT and AP) and rows (first peak to cutoff)----
  thisfile <- thisfilewhole[(firstpkrowindex+1):cutoffrow,
                            c(rtcolindex, areacolindex)]
  #rm(rtcolindex, areacolindex, firstpkrowindex, cutoffrow)
  
  # d. Clean up a little bit----
  colnames(thisfile) <- c("retention_time", "area")  #rename columns!
  #  thisfile <- thisfile[which(!is.na(thisfile[,1])),]  #remove the NAs
  #  thisfile <- thisfile[which(!thisfile[,1]==""),]  #remove the blanks  
  thisfile <- apply(thisfile, c(1,2), as.numeric)  #convert all to numeric
  thisfile <- as.data.frame(thisfile) #convert to data frame in order to be able to call column names 
  
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
  #rm(splitname)
 
  write.csv(thisfile, 
            paste(exportdir, "/", newfilename, ".csv", sep=""), 
            row.names=FALSE)
  
  #rm(thisfilewhole, thisfile, thisfishid, thisfishsp, newfilename)
}

#rm(i)
#rm(importdir, exportdir)
  
  