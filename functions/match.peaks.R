# -------------- #
# Function: match.peaks
# Input: a file of GC data to match, master list of peaks, and a flag to indicate if this is the last round of matching or not
# Output: that same file of GC data with each peak matched to the master list
# -------------- #
# for every peak on the master list, 
# see if it can be found GC file, and if so, add it in
# -------------- #

match.peaks <- function(file, masterlist, lastmatch=FALSE) {
  
  if(!"fameid" %in% colnames(file)) { #if the fameid column doesn't exist yet
    file$fameid <- ""
    file$famename <- ""
  }
  
  for(k in 1:(nrow(masterlist))) { #for each of the peaks in the masterlist
      #see if there are any matches in the list of peaks to match
      #we have to do it this way because 
      #otherwise, sometimes two or more peaks will match to a peak in the masterlist
      #but you won't know it.  
      #it's inefficient - it requires going through the entire master list
      #even if you only have a few peaks to match
      #but i think it has to be this way.
      
      thismasterpeakDFP <- masterlist$dfp[k]
      lower <- thismasterpeakDFP - masterlist$lowerwindow[k]
      upper <- thismasterpeakDFP + masterlist$upperwindow[k]
      
      matches <- subset(file, dfp > lower & dfp < upper)
      matchindex <- as.numeric(rownames(matches))
      
      if(nrow(matches)==0) {   #no match!
#       print("0 matches to this peak!")
        #file$flag[matchindices] <- "this peak in the master doesn't match any in the new list"
      } else if(nrow(matches)==1) {   #1 match!
#        print("this peak in the master has exactly one match in the new list")
          file$fameid[matchindex] <- masterlist$fameid[k]
          file$famename[matchindex] <- masterlist$famename[k]
      } else {#(nrow(matches)>1) #more than 1 match!
#        print(paste("Note! this peak in the master matches", nrow(matches), "FAMEs in the new list."))
        file$flag[matchindex] <- "This peak is one of more than one that matched a single peak in the master - the first time."
        if(lastmatch==TRUE) {
          file$flag[matchindex] <- "This peak is one of more than one that matched a single peak in the master - the last time."
          print("Error! On the last match, a master peak matched more than one FAME!")
          print(paste("The peak in question is: DFP:", thismasterpeakDFP))
          print(matches)
          break
        } 
        
        #otherwise, if lastmatch==FALSE, 
        #match the one with the closest DFP and 
        #let the other one get created in the list as a separate peak
        
        #this gets you the index within the matches list
        closestindexwithinmatches <- which(abs(matches$dfp-thismasterpeakDFP)==
                                min(abs(matches$dfp-thismasterpeakDFP)))
        #thank you: https://stat.ethz.ch/pipermail/r-help/2008-July/167216.html
        if(length(closestindexwithinmatches) > 1) { 
          #if there are two peaks equidistant from the target peak
          #choose the smaller one
          closestindexwithinmatches <- 1
        }
        
        #but you actually need the index from the file
        closestindex <- as.numeric(rownames(matches[closestindexwithinmatches,]))
        
        file$fameid[closestindex] <- masterlist$fameid[k]
        file$famename[closestindex] <- masterlist$famename[k]
        
       #file$flag[thispeakindex] <- paste("peak matched to", nrow(matches), "FAMEs in secondary list")
      } #else          
    } #k in nrow(masterlist)
  
  return(file)
  
} #function
