#assemble list of minimum retention times

minrts <- data.frame(matrix(nrow=length(filenames), ncol=2))
colnames(minrts) <- c("file_name", "min_rt")

for(i in 1:length(filenames)) {
  minrts$file_name[i] <- filenames[i]
  
  thisfile <- filevec[[i]]
  colnames(thisfile) <- c("rt", "a") #shorten column names
  
  minrts$min_rt[i] <- min(thisfile$rt, na.rm=T)
}

minrts

minrts[with(minrts, order(min_rt)),]

