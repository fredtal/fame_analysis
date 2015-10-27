# ------------------ #
# 2: Calculate lipid content
# Talia Young
# June 2015
# 
# This script calculates lipid content based on the C19 addition and the BHT addition.
# It also merges the C:N ratios in to allow comparisons.  

# The calculated lipid contents right now are waaaay off.  
# Many are more than 1 (not possible!).
# Next steps remain unclear.  -ty 6/29/15

# Inputs: 
#   + folder with csv files of GC data
#   + csv with id and FA sample weight
#   + csv with id and c_to_n ratio
# Output: table with columns:
#           + id
#           + sp
#           + min retention time (because we didn't record peaks before 13 min for many samples, which will affect calculations)
#           + start (all or late, based on min retention time)
#           + lipid content based on C19
#           + lipid content based on BHT
#           + C:N ratio
# ------------------ #

#rm(list=ls())

# ------------------ #
# 0. Prep ----
# ------------------ #

# a. Get all the GC data into a list
importdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay shared/Magdalena Bay data/GC data/GC data-clean/csv files-ALL-trimmed"
filenames <- list.files(importdir)
filevec <- vector(mode="list", length(filenames)) #put each file into a list
for(i in 1:(length(filenames))) {
  filevec[[i]] <- read.csv(paste(importdir, filenames[i], sep="/"), 
                           row.names=NULL, header=TRUE,
                           stringsAsFactors=FALSE)
  colnames(filevec[[i]]) <- c("rt", "a")
  
  #add some columns: distance from palmitic, area > 15, ap > 15
  
  #a. add area and ap > 15
  morethan15inds <- which(filevec[[i]]$rt > 15)  
  filevec[[i]]$ap <- (filevec[[i]]$a / sum(filevec[[i]]$a)) * 100
  filevec[[i]]$morethan15a <- ifelse(test=filevec[[i]]$rt > 15, 
                        yes=filevec[[i]]$a, 
                        no=0)
  filevec[[i]]$morethan15ap <- (filevec[[i]]$morethan15a / sum(filevec[[i]]$morethan15a)) * 100
  
  #b. add dfp
  #rt for palmitic acid will be less than 24
  lessthan24inds <- filevec[[i]]$rt < 24 
  #just look at the retention times less than 24 min because 
  #palmitic acid is probably the biggest peak less than 24
  palmiticind <- which.max(filevec[[i]]$morethan15a[lessthan24inds])
  rm(lessthan24inds)
  
  #check to make sure it's reasonable that this peak is actually palmitic acid
  palmiticap <- filevec[[i]]$morethan15ap[palmiticind]
  palmiticrt <- filevec[[i]]$rt[palmiticind]
  if(palmiticap < 10) print (paste(filenames[i], "WARNING!: The palmitic peak is only", 
                                   round(palmiticap, digits=2), "% of area."))
  if(palmiticrt < 21 | palmiticrt > 23) 
    print (paste(filenames[i], "WARNING!: The palmitic retention time is", 
                 palmiticrt, "and not between 21 and 23 minutes."))
  rm(palmiticind, palmiticap)
  
  #create a column that is distance from palmitic acid
  filevec[[i]]$dfp <- filevec[[i]]$rt - palmiticrt
  rm(palmiticrt)
} 
rm(i)

# b. make a table to hold all the lipid content calcs
lipids <- as.data.frame(matrix(NA, nrow=length(filenames), ncol=6))
colnames(lipids) <- c("id", "species", "minrt", "start", "c19", "bht")

pieces <- strsplit(filenames, "[_.]+")
lipids$id <- sapply(pieces, "[", 1)
lipids$species <- sapply(pieces, "[", 2)
rm(pieces)

#designate the digestive glands separately
lipids$species[grep("DG", lipids$id)] <- "squid-DG"

# c. add sample weights
wts <- read.csv("~/Dropbox/Rutgers-Dropbox/Magdalena Bay not shared/Talia code (graphs and stats)/MagBay analysis/data-to import/FA_sample_weights.csv")
colnames(wts)[colnames(wts)=="FA_sample_dry_wt_g"] <- "wt"
wts <- subset(wts, select=c("fish_id", "wt"))
lipids <- merge(lipids, wts, by.x="id", by.y="fish_id", all.x=T, all.y=F)
rm(wts)

# ------------------ #
# 1. Calculate lipid content based on C19 or BHT ----
# ------------------ #

for(i in 1:length(filenames)) {
  print(paste("Processing", filenames[i], "..."))
  thisfile <- filevec[[i]]
  
  # where do peaks start?
  lipids$minrt[i] <- min(thisfile$rt, na.rm=T)
  if(lipids$minrt[i] > 13) lipids$start[i] <- "late" else lipids$start[i] <- "all"

  # ------------------ #
  # a. Calculate lipid content based on C19 (for those samples with C19 additions)
  # ------------------ #
  # C19 solution: 0.5g / 1000 ml = 0.0005 / 1 ml
  # We added 0.5 ml of this solution = 0.00025 g in each sample
  # Add up the area percents of everything except C19, 
  # divide it by area percent of C19
  # then multiply by 0.00025.  
  #
  # Because: area percent of everything:area percent of C19 ~  weight of total lipids:weight of c19
  #
  # lipid percent: [ (area percent everything / area percent of c19) * 0.00025 ] / weight of sample
  # this approach assumes that we took all the lipid out of the sample.  
  # if we didn’t remove all lipids, then the actual numerator would be higher, so we are underestimating lipid content
  
  # find C19 peak
  #The C19 peak usually comes out between 7.5 and 8.5 DFP and is bigger than 4%
  c19 <- subset(thisfile, dfp > 7.5 & dfp < 8.5 & morethan15ap > 4)
  c19ind <- as.numeric(rownames(c19))
 
  if(nrow(c19)>1) { #If there's more than one, take the bigger one
    print("Note: there's more than one peak bigger than 4% in the C19 range between 7.5 and 8.5.")
    print(c19)
    c19 <- c19[which.max(c19$morethan15ap),]
  } 
  
  if(nrow(c19)==1) { #if there is a c19 peak, use it to calculate lipid content
    # add up all the other peaks, divide it by the c19 ap, then multiply by ?0.00025
    lipids$c19[i] <- ((sum(thisfile$a[-c19ind]) / c19$a) * (0.5/1000/2)) / lipids$wt[i]     
  } else if (nrow(c19)==0) { #no c19 spike.
    # if there isn't one, skip this file
    lipids$c19[i] <- NA  
  }
    
  # ------------------ #
  # b. Calculate lipid content based on BHT (for those samples whose plots start before 15)
  # ------------------ #
  # bht usually comes out around 13 minutes and will be bigger than 0.03%
  bhtthreshold <- 0.03
  bht <- subset(thisfile,  rt > 12.8 & rt < 13.5 & ap > bhtthreshold)
  bhtind <- as.numeric(rownames(bht))
  
  if(nrow(bht)>1) { #If there's more than one, take the bigger one
    print(paste("Note: there's more than one peak bigger than", bhtthreshold*100, 
                "% in the BHT range between 12.8 and 13.5."))
    print(bht)
    bht <- bht[which.max(bht$ap),]
  } #if nrow > 1
  
  if(nrow(bht)==1) { #if there is a bht peak, use it to calculate lipid content
    # add up all the other peaks, divide it by the bht ap, then multiply by 0.00025
    lipids$bht[i] <- ((sum(thisfile$a[-bhtind]) / bht$a) * 0.001) / lipids$wt[i]    
  } else if (nrow(bht)==0) { #no bht spike.
    # if there isn't one, skip this file
    lipids$bht[i] <- NA  
  }#else if
} #for i in 1:length(filenames)



# ------------------ #
# 2. Add C:N ratios ----
# ------------------ #
# c. add sample weights
isotopes <- read.csv("~/Dropbox/Rutgers-Dropbox/Magdalena Bay not shared/Talia code (graphs and stats)/MagBay analysis/data-to import/all_isotope_data.csv")
colnames(isotopes) <- tolower(colnames(isotopes))
cton <- subset(isotopes, select=c("fish_id", "c_to_n"))
colnames(cton)[colnames(cton)=="c_to_n"] <- "cton"
lipids <- merge(lipids, cton, by.x="id", by.y="fish_id", all.x=T, all.y=F)

# ------------------ #
# 3. Have a look at some stuff ----
# ------------------ #
# have a look
pdf(paste(plotdir, "lipid_content_looking2.pdf", sep="/"), width=10, height=4)
par(mfrow=c(1,3))
plot(lipids$c19, lipids$bht)
plot(lipids$c19, lipids$cton)
plot(lipids$bht, lipids$cton)
par(mfrow=c(1,1))
dev.off()

lipidstrimmed <- subset(lipids, !(is.na(c19) & is.na(bht)))
lipidsc19 <- subset(lipids, !(is.na(c19)))
lipidsbht <- subset(lipids, !(is.na(bht)))

#require(ggplot2)
plotdir <- ("~/Dropbox/Rutgers-Dropbox/Magdalena Bay not shared/Talia code (graphs and stats)/MagBay analysis/plots")

pdf(paste(plotdir, "lipid_content_looking_c19.pdf", sep="/"), height=20, width=20)
par(mfrow=c(6,6))
for(i in 1:length(unique(lipidsc19$species))) {
  thissp <- unique(lipidsc19$species)[i]
  this <- subset(lipidsc19, species==thissp)
  plot(x=rep(1, nrow(this)), y=this$c19, main=thissp, xlab="", xaxt="n")
  #rm(thissp, this)
}
par(mfrow=c(1,1))
dev.off()

pdf(paste(plotdir, "lipid_content_looking_bht.pdf", sep="/"), height=20, width=20)
par(mfrow=c(6,6))
for(i in 1:length(unique(lipidsc19$species))) {
  thissp <- unique(lipidsc19$species)[i]
  this <- subset(lipidsc19, species==thissp)
  plot(x=rep(1, nrow(this)), y=this$c19, main=thissp, xlab="", xaxt="n")
  #rm(thissp, this)
}
par(mfrow=c(1,1))
dev.off()



