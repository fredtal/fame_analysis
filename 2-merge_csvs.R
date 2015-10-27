# ------------------ #
# 1: Lining up GC data
# Talia Young
# started 12/2012, modified 9/2013, 5/2015
# 
# This script imports csvs of GC data (cleaned up in 0-cleanup_data_export_csvs)
# and syncs all the peaks!  
# (Then you have to finish the syncing by hand.)

# Input: directory with csvs
# Output: single csv file with semi-synced peaks
# ------------------ #

#rm(list=ls())

#for testing
importdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay shared/Magdalena Bay data/GC data and reports/GC data-cleaned/reports-all, organized by file type/testing"
dataexportdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay not shared/Talia code (graphs and stats)/MagBay analysis/data-exported"

importdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay shared/Magdalena Bay data/GC data and reports/GC data-cleaned/reports-all, organized by file type/csv files-ALL-trimmed"
dataexportdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay not shared/Talia code (graphs and stats)/MagBay analysis/data-exported"
codedir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay not shared/Talia code (graphs and stats)/MagBay analysis/code"

peakthreshold <- 0.5

# ------------------ #
#  1. Import files, ids, and species -----
# ------------------ #

#make the list of files you want to import
filenames <- list.files(importdir)

#get list of ids and sp
#pieces <- strsplit(filenames, "[_.]+")
#idvec <- sapply(pieces, "[", 1)
#spvec <- sapply(pieces, "[", 2)
#rm(pieces)

#make list of files
filevec <- vector(mode="list", length(filenames))
filevec2 <- vector(mode="list", length(filenames))

for(i in 1:(length(filenames))) {
  filevec[[i]] <- read.csv(paste(importdir, filenames[i], sep="/"), 
                            row.names=NULL, header=TRUE,
                            stringsAsFactors=FALSE)
} 
rm(i)
#spuniq <- sort(unique(spvec))

#do together: squid, sharks, stm, dolphinfish, flying fish

#check species
#calculate 

# ------------------ #
# Algorithm! ----
# I. Generate list of all the files.  
# II. For each file, 
#   1. Import and clean up file.
#   2. Calibrate file based on palmitic acid.
#   3. Now, for the peaks.
#     >> if this is the first file...
#     a. Turn all the peaks into a peak list.  
#     >> if this is the not the first file...
#     a. Match the peaks to the peak list.
#     b. Deal with the unmatched peaks
#  III. For each file again,
#   1. Rematch each file to the final peak list.
#   2. Merge thisfile with alldata.  
#   3. Assign DFP base values to any new FAMEs.  
#  IV. Export data.
# ------------------ #


# ------------------ #
#  II. GENERATE LIST OF ALL PEAKS ----
# ------------------ #

for(i in 1:length(filenames)) {
  print(paste("Processing", filenames[i], "..."))
  thisfile <- filevec[[i]]
  colnames(thisfile) <- c("rt", "a") #shorten column names
  
  #remove all rows with RT < 15 min
  lessthan15inds <- thisfile$rt < 15
  thisfile <- thisfile[!lessthan15inds,]
  rm(lessthan15inds)

  thisfile$ap <- (thisfile$a/sum(thisfile$a))*100 #add area percent column
  
  #renumber the rows!  so you can call the indices later.  
  rownames(thisfile) <- 1:nrow(thisfile)
  #thank you: http://stackoverflow.com/questions/12505712/renumbering-rows-after-ordering-in-r-programme
  
  # ------------------ #
  #  2. Calibrate the file based on palmitic acid  ----
  # ------------------ #
  
  #rt for palmitic acid will be less than 24
  lessthan24inds <- thisfile$rt < 24 
  #just look at the retention times less than 24 min because 
  #palmitic acid is probably the biggest peak less than 24
  palmiticind <- which.max(thisfile$a[lessthan24inds])
  rm(lessthan24inds)
  
  #check to make sure it's reasonable that this peak is actually palmitic acid
  palmiticap <- thisfile$ap[palmiticind]
  palmiticrt <- thisfile$rt[palmiticind]
  if(palmiticap < 10) print (paste("WARNING!: The palmitic peak is only", 
                                   round(palmiticap, digits=2), "% of area."))
  if(palmiticrt < 21 | palmiticrt > 23) 
    print (paste("WARNING!: The palmitic retention time is", 
                 palmiticrt, "and not between 21 and 23 minutes."))
  rm(palmiticind, palmiticap)
           
  #create a column that is distance from palmitic acid
  thisfile$dfp <- thisfile$rt - palmiticrt
  rm(palmiticrt)
  
  #remove late peaks
  timeendthreshold <- 22.4 # timeendthreshhold ----
  thisfile <- thisfile[thisfile$dfp < timeendthreshold,] #remove peaks after 22.4 minutes DFP
  rm(timeendthreshold)
  
  #renumber the rows!  so you can call the indices later.  
  rownames(thisfile) <- 1:nrow(thisfile)
  #thank you: http://stackoverflow.com/questions/12505712/renumbering-rows-after-ordering-in-r-programme
  
  # ------------------ #
  #  3. Remove the C19 peak  ----
  # ------------------ #
  
  #The C19 peak usually comes out between 7.5 and 8.5 DFP and is bigger than 4%
  c19 <- subset(thisfile, dfp > 7.5 & dfp < 8.5 & ap > 4)
  
  if(nrow(c19)==0) { 
    #no c19 spike.  no worries.
  } else if (nrow(c19)>1) {
    #If there's more than one, take the bigger one
    print("Note: there's more than one peak bigger than 4% in the C19 range between 7.5 and 8.5.")
    print(c19)
    c19 <- c19[which.max(c19$Areapercent),]
  }   
  #then, remove the row.  
  if(nrow(c19)==1) {
    c19index <- as.numeric(rownames(c19))
    thisfile <- thisfile[-c19index,] #remove the line
    rownames(thisfile) <- 1:nrow(thisfile) #renumber
    #recalculate the area percents without the c19 peak.  
    thisfile$ap <- thisfile$ap / sum(thisfile$ap) * 100
  }  
  rm(c19, c19index)
  
  # ------------------ #
  #  4. Remove small peaks  ----
  # ------------------ #
  #remove all peaks < 0.2%
  thisfile <- thisfile[thisfile$ap > peakthreshold,]
  
  #renumber the rows!  so you can call the indices later.  
  rownames(thisfile) <- 1:nrow(thisfile)
  #thank you: http://stackoverflow.com/questions/12505712/renumbering-rows-after-ordering-in-r-programme
  
  
  # ------------------ #
  #  5. Now, generate list of all peaks  ----
  # ------------------ #
  
  #add FAMEid and FAMEname columns to thisfile
  thisfile["fameid"] <- ""
  thisfile["famename"] <- ""
  thisfile["flag"] <- ""
  
  noname <- "unnamed"  
  #  >> If this is the first file...  ----
  if(i==1) {  
    
    #  a. Turn all the peaks into a peak list for future reference ----    
    #give an ID and name to all the fames
    for(j in 1:(nrow(thisfile))) {
      thisfile$fameid[j] <- paste("FAME @", 
                                  round(thisfile$dfp[j], digits=2))
      thisfile$famename[j] <- noname
    } #j in (nrow(thisfile))
    rm(j)
    
    #the peak list is all same as thisfile minus the flag
    peaklist <- subset(thisfile, select=-flag)
    print(paste("Note:", nrow(peaklist), "peaks created for initial peak list."))
    
    #add the lowerwindow and upperwindow columns
    peaklist["lowerwindow"] <- NA
    peaklist["upperwindow"] <- NA    
    
    #sort by dist_from_palmitic
    peaklist <- peaklist[with(peaklist, order(dfp)),]
    rownames(peaklist) <- 1:nrow(peaklist)
    
    #run the calculate.RT.windows function on the peak list ----
    source(paste(codedir, "calculate.RT.windows.R", sep="/"))
    defaultwindow <- 0.8 #defaultwindow ----
    
    peaklist <- calculate.RT.windows(datatable=peaklist, 
                                     defaultwindow=defaultwindow)
    
    #  >> If this is NOT the first file...  ----
    
  } else { #if(i>1) { 
    
    #  a. Match the peaks to the peak list ----
    source(paste(codedir, "match.peaks.R", sep="/"))
    thisfile <- match.peaks(thisfile, peaklist)
    
  #  b. Deal with the unmatched peaks ----
    unmatched <- subset(thisfile, fameid=="")
    unmatchedindices <- as.numeric(rownames(unmatched))
    
    if(nrow(unmatched)!=0) { #if there are any unmatched       
      
      #  ii. Name the unmatched peaks based on retention time ----
      thisfile[unmatchedindices,]$fameid <- paste("FAME @", 
                                                  round(thisfile[unmatchedindices,]$dfp,digits=2))
      thisfile[unmatchedindices,]$famename <- noname
      
      #  iii. Then add the unmatched peaks to the peak list ----
      #I had to call the unmatched again in order to get the name with it this time
      #also, when appending it to the peak list, we don't want Pkno or flag
      unmatchedwithname <- thisfile[unmatchedindices, !colnames(thisfile) %in% "flag"]
      #thank you: http://stackoverflow.com/questions/4605206/drop-columns-r-data-frame
      unmatchedwithname["lowerwindow"] <- NA
      unmatchedwithname["upperwindow"] <- NA
      
      peaklist <- rbind(peaklist, unmatchedwithname)
      peaklist <- peaklist[with(peaklist, order(dfp)),]
      rownames(peaklist) <- 1:nrow(peaklist) #reassign row numbers
      peaklist <- calculate.RT.windows(datatable=peaklist, 
                                       defaultwindow=defaultwindow)
      rm(unmatchedwithname)
    } #if(nrow(unmatched)!=0)
  rm(unmatched, unmatchedindices)
  } #else if filecount > 1
  rm(noname)
  
  #d. Save thisfile so we can look at it later.
  filevec2[[i]] <- thisfile
  
  rm(thisfile)
  
} #i in 1:filenames
rm(defaultwindow)
rm(i)



# ------------------ #
# III. Now match each file to the master peak list ----
# ------------------ #

#for(i in 1:3) {
for(i in 1:length(filenames)) {
  thisfile2<- filevec2[[i]]
  
  #thisfishid ----
  #splitname <- strsplit(filenames[i], "[[:punct:][:space:]]")[[1]]
  splitname <- strsplit(filenames[i], "[_.]+")[[1]]  
  #must use double brackets to vector which is first item in the list which is result of strsplit.  ugh
  #thank you: http://www.dummies.com/how-to/content/how-to-split-strings-in-r.html
  thisfishid <- splitname[1]
  rm(splitname)
  
  #1. rematch each file with the final peak list!----
  print(paste("Matching peaks for", filenames[i]))
  thisfile2 <- match.peaks(thisfile2, peaklist, lastmatch=TRUE)
  
  unmatched <- subset(thisfile2, fameid=="")
  if(nrow(unmatched) > 0) print(paste("Error: unmatched peaks remain in", 
                                      filenames[i], "after final match with peaklist."))
  rm(unmatched)
  
  # >> IF THIS IS THE FIRST FILE...
  if(i==1) {
    #  2. Create big file  ----   
    thisfishdata <- thisfile2[c("fameid", "dfp", "dfp", "rt", "ap")] #dfp twice once for base, once for regular
    
    alldata <- data.frame(thisfishdata)
    #use data.frame instead of cbind to get a data frame (which can hold multiple data types)
    #instead of a matrix.
    #thank you: http://stackoverflow.com/questions/11151339/r-numeric-vector-becoming-non-numeric-after-cbind-of-dates
    
    #rename the columns
    colnames(alldata) <- c("fameid", 
                           "dfp.base",
                           paste("dfp.", thisfishid, sep=""), 
                           paste("rt.", thisfishid, sep=""), 
                           paste("ap.", thisfishid, sep=""))
    rm(thisfishdata)
  } else if(i > 1) {
    
    # 2. merge thisfile's data with alldata, by fameid ----
    thisfishdata <- thisfile2[c("fameid", "dfp", "rt", "ap")]
    
    #rename the columns
    colnames(thisfishdata) <- c("fameid", 
                                paste("dfp.", thisfishid, sep=""),
                                paste("rt.", thisfishid, sep=""), 
                                paste("ap.", thisfishid, sep=""))
    rm(thisfishid)
    
    mergedbit <- merge(alldata, thisfishdata, 
                       by.x="fameid", by.y="fameid", all.x=TRUE, all.y=TRUE)
    rm(thisfishdata)
    
    # 3. assign DFP base values to any new FAMEs ----
    
    #pull out all the rows that don't have an RT.base value (except for the Totals line)
    emptyDFPbase <- subset(mergedbit, is.na(dfp.base))
    emptyDFPbasenameindices <- rownames(emptyDFPbase)
    
    #add DFP.base values
    if(nrow(emptyDFPbase) != 0) { #if there are any new ones
      #for each row that doesn't have an RT.base value...
      for(j in 1:length(emptyDFPbasenameindices)) {
        nonNAs <- which(!is.na(emptyDFPbase[j,])) #returns col indices of non-NA vals of this row        baseRTindex <- min(nonNAs[2:length(nonNAs)]) 
        #the base RT index should be the index of the first value after the first column that is not NA in that row
        #thank you: http://stackoverflow.com/questions/6808621/find-the-index-position-of-the-first-non-na-value-in-an-r-vector
        DFPbaseindex <- min(nonNAs[2:length(nonNAs)])
        #set the RT.base value to the RTbaseindex that you found
        mergedbit[emptyDFPbasenameindices[j],]$dfp.base <- 
          mergedbit[emptyDFPbasenameindices[j],DFPbaseindex]
        rm(nonNAs, DFPbaseindex)
      } #j
      rm(j)
    } #if
    rm(emptyDFPbase, emptyDFPbasenameindices)
    
    #sort merged bit by dist_from_palmitic
    mergedbit <- mergedbit[with(mergedbit, order(dfp.base)),]
    
    # 4. change all NA values in areapercent columns after first 2 columns to 0
    #index vector: which colnames contain "Areapercent"
    areapercentindices <- grep("ap", colnames(mergedbit), value=FALSE)
    
    #for all the columns with areapercent, 
    #find the values that are na, and assign them to 0 (yay!)
    mergedbit[,areapercentindices][is.na(mergedbit[,areapercentindices])] <- 0
    alldata <- mergedbit
    rownames(alldata) <- 1:nrow(alldata)  
    rm(areapercentindices)
    rm(mergedbit)
  } #else if filecount > 1 
  rm(thisfile2)
} #for i in 1:length(filenames)
rm(i)

write.csv(alldata, paste(dataexportdir, "alldata.csv", sep="/"), row.names=F)
write.csv(peaklist, paste(dataexportdir, "peaklist.csv", sep="/"), row.names=F)

# ------------------ #
# IV. EXPORT DATA ----
# ------------------ #

# extract ap and dfp columns
apdfpcolindices <- sort(union(grep("ap", colnames(alldata)), 
                              grep("dfp", colnames(alldata))))
#basecolindex <- grep("base", colnames(alldata))
#apdfpcolindices <- apdfpcolindices[apdfpcolindices!=basecolind] #remove the base column
APDFP <- alldata[, apdfpcolindices] 
rm(apdfpcolindices)

#get ids for species matching
splitnames <- strsplit(colnames(APDFP), ".", fixed=TRUE) #split on period
idsonly <- sapply(splitnames, "[[", 2)
idsonly <- idsonly[2:length(idsonly)] #remove the "base" column name for species matching below
#thank you: http://stackoverflow.com/questions/2803460/how-to-get-the-second-sub-element-of-every-element-in-a-list-in-r#then rename the columns
rm(splitnames)

id_sp_batch <- read.csv("~/Dropbox/Rutgers-Dropbox/Magdalena Bay not shared/Talia code (graphs and stats)/MagBay analysis/data-to import/id_species_batchnum.csv", 
                        stringsAsFactor=F)
id_sp <- id_sp_batch[c("fish_id", "common_name")]

# match ids to species
idsonlydf <- as.data.frame(idsonly, stringsAsFactors=F) #from the APs
id_sp_match <- merge(idsonlydf, id_sp, by.x="idsonly", by.y="fish_id", all.x=T)
#id_sp_match[is.na(id_sp_match$common_name),] #check to see if they all matched
#rm(idsonlydf)
sp <- id_sp_match$common_name#make a column of just the species, so the col name is right later.  
sp <- c("", id_sp_match$common_name) #add a blank at the front for the "base" column

#bind the species to the APs
#it's just easier for me to do the sort vertically, so i am transposing and and then retransposing it
APDFPt <- t(APDFP)
APDFPt <- cbind.data.frame(sp, APDFPt)
APDFPtsp <- APDFPt[with(APDFPt, order(sp)), ]
APDFP <- t(APDFPtsp)
write.csv(APDFP, paste(dataexportdir, "APDFP.csv", sep="/"), row.names=F)




