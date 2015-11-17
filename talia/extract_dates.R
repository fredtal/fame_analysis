# extract dates from gc files
# nov 2015
#
# this code only works for class-vp files right now because i already removed the headers from the LS files.  :/

#rm(list=ls())

# set directories
importdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay shared/Magdalena Bay data/GC data/GC data-clean/csv files-various/csv files-ALL-original"
exportdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay shared/Magdalena Bay data/GC data/GC data-clean/csv files-various/csv files-nov15-trimmed"

# generate list of files
filelist <- list.files(importdir)

library(stringr) #for str_trim
datetable <- data.frame(matrix(nrow=length(filelist), ncol=2))
colnames(datetable) <- c("file_name", "date")
datetable$date <- as.Date(datetable$date)

for(i in 1:length(filelist)) {
  print(paste("Processing", filelist[i]))
  if(!grepl("Blank", filelist[i]) & !grepl("blank", filelist[i])) { #skip the blanks
    
    # a. Import and clean up the file ----  
    thisfilewhole <- read.csv(paste(importdir, filelist[i], sep="/"), 
                              row.names=NULL, header=FALSE,
                              stringsAsFactors=FALSE)
    
    #in the class VP files, the date is next to the "Run time"
    if(sum(grepl("Class-VP", thisfilewhole))>0) { #if it's a Class-VP file...
      if(sum(grepl("Run time", thisfilewhole))>0) { #if it has Run Time
        colindex <- grep("Run time", thisfilewhole) 
        rowindex <- grep("Run time", thisfilewhole[,colindex]) 
        datemessy <- thisfilewhole[rowindex, colindex+1]
        datemessy <- str_trim(datemessy, side="left")
        #http://stackoverflow.com/questions/5992082/how-to-remove-all-whitespace-from-a-string
        #http://stackoverflow.com/questions/2261079/how-to-trim-leading-and-trailing-whitespace-in-r
      } else { #else on run time
        colindex <- grep("Date", thisfilewhole)
        rowindex <- grep("Date", thisfilewhole[,colindex])
        datemessy <- thisfilewhole[rowindex+1, colindex]
      } #close else
      date <- as.Date(datemessy, format="%m/%d/%Y")
      datetable$date[i] <- date
      rm(date, datemessy, colindex, rowindex)
      
      
      splitname <- strsplit(filelist[i], "[_.]+")[[1]] #split on underscore OR period
      #thank you: http://stackoverflow.com/questions/10738729/r-strsplit-with-multiple-unordered-split-arguments
      #also: http://stat.ethz.ch/R-manual/R-patched/library/base/html/regex.html
      
      thisfishid <- splitname[1] 
      thisfishsp <- tolower(splitname[2])
      if(thisfishsp=="striped marlin") thisfishsp="stm"
      if(thisfishsp=="csv") { #if blank, there's no species
        newfilename <- thisfishid 
      } else {
        newfilename <- paste(thisfishid, thisfishsp, sep="_")
        newfilename <- paste(newfilename, "csv", sep=".")
      }
      #rm(splitname)
      
      datetable$file_name[i] <- newfilename
      
      
    } #close class VP
    rm(thisfilewhole)
  } #close blank check
} # close for

write.csv(datetable, "~/Dropbox/Rutgers-Dropbox/Magdalena Bay not shared/MagBay analysis/data-exported/date_table.csv")

