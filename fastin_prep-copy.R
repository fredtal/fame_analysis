# ------------------ #
# Prep data for FASTIN
# Talia Young
# May 2015
#
# This script assumes that you've sync-ed up all the FA rows.
#  + It imports a file of SI data.  
#  + It imports a file of FA data. 
#    (usually with one AP and one RT row for each fish) 
#    (transposed from the hand-synced file)
#  + It will separate out the predators and prey and 
#    export the data for fastin. 
#  + It will also make supplemental files (cc.mean.var, etc.)
#  + There's plot code at the bottom.  
# ------------------ #

#rm(list=ls())
fastindir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay not shared/Talia code (graphs and stats)/MagBay analysis/data-to import/data for fastin"
plotdir <- "~/Dropbox/Rutgers-Dropbox/Magdalena Bay not shared/Talia code (graphs and stats)/MagBay analysis/plots/fastin"

#import SI and FA data
source('~/Dropbox/Rutgers-Dropbox/Magdalena Bay not shared/Talia code (graphs and stats)/MagBay analysis/code/import_si_and_fa_data.R')

# ------------------ #
# 1. Export SI data ----
# ------------------ #
# predators SI ----
stmsif <- subset(stmsi, select=c("id", "d13c", "d15n"))
write.csv(stmsif, paste0(fastindir, "/striped marlin_si.csv"),
          row.names=F)

# prey SI ----
#remove end members
preysif <- subset(preysi, !sp=="zooplankton")
preysif <- subset(preysi, !cat=="mollusk")

preysif <- subset(preysi, !sp %in% prednames[prednames != "dolphinfish"], 
                 select=c("cat", "d13c", "d15n"))
preysif <- preysif[preysif$cat !="other",] #remove prey not in one of the categories

# export all prey in one file, with species as the first column
write.csv(preysif, paste0(fastindir, "/prey_si.csv"), row.names=F)


# ------------------ #
# 2. Generate SI supplemental files ----
# ------------------ #
prey.list.si <- unique(as.character(preysif$cat))
if(length(which(prey.list.si=="other")) > 0) {
  prey.list.si <- prey.list.si[-which(prey.list.si=="other")] #remove other
}
n.prey.si <- length(prey.list.si)
n.si <- 2

# SI fractionation coefficient means
fc.mean.si <- as.data.frame(matrix(nrow=n.prey.si, ncol=n.si,
                     dimnames=list(prey.list.si, c("d13c", "d15n"))))
fc.mean.si$d13c <- 1.8 
fc.mean.si$d15n <- 1.9
# per Madigan 2012

# SI fractionation coefficient var
fc.var.si <- as.data.frame(matrix(nrow=n.prey.si, ncol=n.si,
                     dimnames=list(prey.list.si, c("d13c", "d15n"))))
fc.var.si[is.na(fc.var.si)] <- 0.2 #set all equal to 0.2 for now

# export ----
write.csv(fc.mean.si, paste(fastindir, "fc_mean_si.csv", sep="/"))
write.csv(fc.var.si, paste(fastindir, "fc_var_si.csv", sep="/"))

rm(prey.list.si, n.prey.si, n.si)
rm(fc.mean.si, fc.var.si)
rm(preysif, stmsif)


# ------------------ #
# 3. Separate out predators and prey FA data and export ----
# ------------------ #
# predators ----
stmfaf <- stmfa
rownames(stmfaf) <- stmfaf$id
stmfaf <- subset(stmfa, select=-c(sp, predator_cat))
write.csv(stmfaf, paste0(fastindir, "/striped marlin_fa.csv", sep=""), 
          row.names=F)
rm(stmfaf)

preyfaf <- preyfa
preyfaf$sp <- preyfaf$cat
preyfaf <- subset(preyfaf, cat!="other", select=-c(id, cat))
write.csv(preyfaf, paste0(fastindir, "/", "prey_fa.csv"), row.names=F)

# ------------------ #
# 4. Generate FA supplemental files ----
# ------------------ #
prey.list.fa <- unique(as.character(preyfa$cat))
prey.list.fa <- prey.list.fa[-which(prey.list.fa=="other")] #remove other
n.prey.fa <- length(prey.list.fa)
n.fa <- length(colnames(preyfa)[grep("X", colnames(preyfa))])

# FA correlation coefficient means and variances ----
# FA CC means
cc.mean.fa <- matrix(nrow=n.prey.fa, ncol=n.fa,
                     dimnames=list(prey.list.fa, 
                                   colnames(preyfa)[grep("X", colnames(preyfa))]))
cc.mean.fa[is.na(cc.mean.fa)] <- 1 #set all equal to 1 for now
# FA CC vars
cc.var.fa <- matrix(nrow=n.prey.fa, ncol=n.fa,
                    dimnames=list(prey.list.fa, 
                                  colnames(preyfa)[grep("X", colnames(preyfa))]))
cc.var.fa[is.na(cc.var.fa)] <- 0.05 #set all equal to 0.05 for now

# fat content: ----
# default: mean=5, var=0.5
fat.cont <- data.frame(row.names=prey.list.fa,
                       rep(5, times=n.prey.fa), rep(0.5, times=n.prey.fa))
colnames(fat.cont) <- c("mean", "var")

# set some specific species based on literature (for now): 
# use lipid content calculations later

# piscivore
# greenjack 3.74% and sierra 3.09% (Murillo 2014)
#fat.cont$mean[rownames(fat.cont)=="piscivore"] <- 3.5

# squid 1.1-13.2% (Iverson 2002)
fat.cont$mean[grep("squid", rownames(fat.cont))] <- mean(c(1.1, 13.2))

# octopus 0.9-1.4% (Iverson 2002)
fat.cont$mean[grep("octopus", rownames(fat.cont))] <- mean(c(0.9, 1.4))

# cephalopod
fat.cont$mean[grep("cephalopod", rownames(fat.cont))] <- mean(c(0.9, 1.1, 1.4, 13.2))

# pelagic schooling fish
# atlantic mackerel 3-19% (Wallace 1991), king mackerel 1.96% (Fernandes 2014)
# mediterranean sardine 1.2-18.2% (Bandarra 2006), 
# herring 0.5-10.4% (Iverson 2002) #set herring default
# Atlantic thread herring (Fernandes 2014) 9.03, herring 0.5-10.4% (Iverson 2002)
fat.cont$mean[rownames(fat.cont)=="pelagic fish"] <- mean(c(3, 19, 1.96, 0.5, 10.4, 1.1, 13.2))
fat.cont$mean[rownames(fat.cont)=="fish"] <- mean(c(3, 19, 1.96, 0.5, 10.4, 1.1, 13.2))

# PRC 5-14% (Goytortua-Bores 2006)
fat.cont$mean[rownames(fat.cont)=="pelagic red crab"] <- mean(c(5,14))

# mediterranean sardine 1.2-18.2% (Bandarra 2006), 
fat.cont$mean[grep("sardine", rownames(fat.cont))] <- mean(c(1.2, 18.2))

# export ----
write.csv(cc.mean.fa, paste(fastindir, "cc_mean_fa.csv", sep="/"))
write.csv(cc.var.fa, paste(fastindir, "cc_var_fa.csv", sep="/"))

write.table(fat.cont, paste(fastindir, "fat_content.csv", sep="/"), sep=",",
            row.names=T, col.names=F) #must be write table to get row.names=T and col.names=F

rm(prey.list.fa, n.prey.fa, n.fa)
rm(cc.mean.fa, cc.var.fa, fat.cont)
rm(preyfaf)

print("fastinR data prepped!")
