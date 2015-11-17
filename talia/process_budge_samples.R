# ------------------ #
# Process Budge samples so I can see what they look like
# ------------------ #

#rm(list=ls())

# importdir = where you have the files you want to sync
importdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay shared/Magdalena Bay data/GC data/GC data-clean/budge samples/original"

# dataexportdir = where you want to dump the csv end product
dataexportdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay shared/Magdalena Bay data/GC data/GC data-clean/budge samples/processed"

# codedir = where you have the 2 functions below stored.
#codedir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay not shared/fame_analysis/functions"

# Import function to calculate retention time windows
# (This function sets the window around a given peak to be the default window.
#  If that window runs into another peak window, the function resets the window
#  to equal half the difference.)
#source(paste(codedir, "calculate.RT.windows.R", sep="/"))
#defaultwindow <- 0.8 #defaultwindow ----    

# Import function to match peaks
#source(paste(codedir, "match.peaks.R", sep="/"))

# set target num peaks (all files should be pared to ~this target number of peaks, to get comparable resolution)
targetpeaknum <- 50

# ------------------ #
#  1. Import files, ids, and species -----
# ------------------ #

#make the list of files you want to import
filenames <- list.files(importdir)

#get list of ids and sp
pieces <- strsplit(filenames, "[_.]+")
#idvec <- sapply(pieces, "[", 1)
filenamevec <- paste(sapply(pieces, "[", 1), sapply(pieces, "[", 2), sep="_")
spvec <- sapply(pieces, "[", 2)
rm(pieces)

#make list of files
filevec <- vector(mode="list", length(filenames))
#filevec2 <- vector(mode="list", length(filenames))

for(i in 1:(length(filenames))) {
  filevec[[i]] <- read.csv(paste(importdir, filenames[i], sep="/"), 
                           row.names=NULL, header=TRUE,
                           stringsAsFactors=FALSE)
} 
rm(i)

# make table for peakcounts
peakcounts <- data.frame(matrix(nrow=length(filenames), ncol=9))
colnames(peakcounts) <- c("file_name", "num_peaks", "peak_threshold", "palmitic_rt", "palmitic_ap", "c19_rt", "c19_dfp", "c19_ap", "notes")


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
  peakcounts$palmitic_ap[i] <- palmiticap
  peakcounts$palmitic_rt[i] <- palmiticrt
  #print(paste0("  Palmitic acid is ", round(palmiticap, digits=1), "%, at ", round(palmiticrt, digits=2), " minutes."))
  if(palmiticap < 10) {
    #print(paste("  WARNING!: The palmitic peak is only", round(palmiticap, digits=2), "% of area."))
    peakcounts$notes[i]  <- paste("WARNING!: The palmitic peak is only", round(palmiticap, digits=2), "% of area.")
  } #if pap <10
  if(palmiticrt < 21 | palmiticrt > 23) {
    #print(paste("WARNING!: The palmitic retention time is", palmiticrt, "and not between 21 and 23 minutes."))
    peakcounts$notes[i]  <- paste(peakcounts$notes[i], "WARNING!: The palmitic peak is only", round(palmiticap, digits=2), "% of area.")
  }
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
    c19 <- c19[which.max(c19$Areapercent),] #set c19 to the biggest peak in that range
    peakcounts$notes[i] <- paste(peakcounts$notes[i], "Note: there's more than one peak bigger than 4% in the C19 range between 7.5 and 8.5.")
  }   
  #then, remove the row.  
  if(nrow(c19)==1) {
    c19index <- as.numeric(rownames(c19))
    peakcounts$c19_rt[i] <- c19$rt
    peakcounts$c19_ap[i] <- c19$ap
    peakcounts$c19_dfp[i] <- c19$dfp
    thisfile <- thisfile[-c19index,] #remove the line
    rownames(thisfile) <- 1:nrow(thisfile) #renumber
    #recalculate the area percents without the c19 peak.  
    thisfile$ap <- thisfile$ap / sum(thisfile$ap) * 100
    rm(c19index)
  }  
  rm(c19)
  
  # ------------------ #
  #  4. Remove small peaks  ----
  # ------------------ #
  #remove all peaks < peakthreshold, which is set at top
  #  thisfile <- thisfile[thisfile$ap > peakthreshold,]
  
  #renumber the rows!  so you can call the indices later.  
  #  rownames(thisfile) <- 1:nrow(thisfile)
  #thank you: http://stackoverflow.com/questions/12505712/renumbering-rows-after-ordering-in-r-programme
  
  # ------------------ #
  #  4. Remove extra peaks in files with a lot of them  ----
  # ------------------ #
  # remove peaks with 0 ap
  thisfile <- thisfile[thisfile$ap!=0,]  
  rownames(thisfile) <- 1:nrow(thisfile)
  
  peakcounts$file_name[i] <- filenames[i]
  peakcounts$num_peaks[i] <- nrow(thisfile)
  
  # we want to get all the files to a similar resolution, to have about the same number of peaks
  # the target number of peaks is set above
  if(nrow(thisfile)>(1.5*targetpeaknum)){ #only pare peaks if the file has > 1.5*targetnumpeaks
    #find a threshold that gives us the right number
    peakthreshold <- sort(thisfile$ap, decreasing=T)[1.5*targetpeaknum] #cutoff is the ap where the peaknum = 1.5*targetnumpeaks
    numpeaksremoved <- length(thisfile$ap[thisfile$ap<=peakthreshold])
    thisfile <- thisfile[thisfile$ap >= peakthreshold,] #remove all peaks below that threshold
    thisfile$ap <- thisfile$ap/sum(thisfile$ap)*100
    peakcounts$peak_threshold[i] <- peakthreshold #save the peakthreshold so you can look at it later
    #renumber the rows!  so you can call the indices later.  
    rownames(thisfile) <- 1:nrow(thisfile)
    #thank you: http://stackoverflow.com/questions/12505712/renumbering-rows-after-ordering-in-r-programme    
    rm(peakthreshold)
  }
  
  write.csv(thisfile, paste0(dataexportdir, "/", filenamevec[i], "_processed.csv"))
  rm(thisfile, numpeaksremoved)
}
rm(i)
rm(peakcounts, filenames, filenamevec, filevec, spvec, targetpeaknum)
rm(importdir, dataexportdir)