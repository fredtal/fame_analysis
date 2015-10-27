#create retention time windows

calculate.RT.windows <- function(datatable, defaultwindow) {
  
  #set all window values that are NA equal to the default window instead
  datatable$lowerwindow[is.na(datatable$lowerwindow)] <- defaultwindow
  datatable$upperwindow[is.na(datatable$upperwindow)] <- defaultwindow
  
  for(k in 1:(nrow(datatable)-1)) { #only go to the second-to-last peak
    thisDFP <- datatable$dfp[k]
    thisDFPplus <- thisDFP + datatable$upperwindow[k]
    nextDFP <- datatable$dfp[k+1]     
    nextDFPminus <- nextDFP - datatable$lowerwindow[k+1]
    if(thisDFPplus > nextDFPminus) { #if the windows run into each other 
      #split the difference between the peaks for the window ----
      splitdiff <- (nextDFP - thisDFP) / 2
      
      #set the window above for this one, and the window below for the next one
      #to be the difference between the two
      datatable$upperwindow[k] <- splitdiff
      datatable$lowerwindow[k+1] <- splitdiff
    } #if thisDFPplus > nextDFPminus
  } #k in datatable  
  
#  rownames(datatable) <- 1:nrow(datatable)
  return(datatable)
}