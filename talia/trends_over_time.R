datetablebackup <- datetable
#datetable <- datetable[!duplicated(datetable),]

datetable$file_name2 <- ifelse(test=grepl("-A_", datetable$file_name), yes=datetable$file_name, 
                               no=ifelse(test=grepl("A_", datetable$filename), 
                                      yes=gsub("A_", "-A_", datetable$file_name), 
                                      no=datetable$file_name))

datetable$file_name2 <- gsub("A_", "-A_", datetable$file_name)
datetable$file_name2 <- ifelse(test=grepl("-A_", datetable$file_name), yes=datetable$file_name, 
                               no=datetable$file_name2)
datetable$file_name2 <- ifelse(test=grepl("_sardine", datetable$file_name), yes=datetable$file_name, 
                               no=datetable$file_name2)

datetable2 <- datetable[c("date", "file_name2")]
colnames(datetable2)  <- c("date", "file_name")
datetable2 <- datetable2[!is.na(datetable2),]

pd <- merge(peakcounts, datetable2)


plot(pd$palmitic_rt, pd$date)
plot(pd$date, pd$palmitic_rt)

plot(pd$c19_rt, pd$date)
plot(pd$date, pd$palmitic_ap)
