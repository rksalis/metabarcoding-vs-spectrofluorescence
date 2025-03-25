# data visualisation and analysis of species Nov 26 2021
# R. Salis - prop

library(data.table)
library(ggplot2)
library(ggpubr)
library(rlist)#
library(tidyverse)
library(dplyr)
library(tidyr)
library(tidytext)#
library(RColorBrewer)
library(plyr)
library(scales)
library(patchwork)


asv_table1 <- fread("asv_table_eDNA1_18S_pr2_mesocosms_vs.csv")
tax_table <- fread("tax_table_eDNA1_18S_pr2_mesocosms_vs.csv",header=TRUE)
sample_data <- fread("sample_data_eDNA1_18S_pr2_mesocosms_vs.csv")

asv_table <- dcast(melt(asv_table1, id.vars = "V1"), variable ~ V1)
dt_ASVsamp <- sample_data[asv_table, on = .(V1 = variable)]
#convert to long format and retain only the variables want to keep
dt_long <- melt(dt_ASVsamp,
                id.vars = c("Sample.date","sampling.week","sampling.occasion","treatment","Mesocosm","Sample.number","rep","Heated","Invasion"),
                measure.vars = patterns("^ASV"),
                variable.name = "ASV",
                value.name = c("normalised_reads"))

#add taxonomy information
dt_long_sp <- dt_long[tax_table, on = .(ASV = V1)]

#convert back to wide
dt_wide_sp <- dcast(dt_long_sp, Kingdom+Supergroup+Division+Class+Order+Family+Genus+Species+treatment+Mesocosm+rep~Sample.date,
                    value.var =  c("normalised_reads"))
dt_wide_sp

####### look at algal groups
dt_long_sp_Chlorophyta <- dt_long_sp[Division == "Chlorophyta", ]
dt_long_sp_Chlorophyta[,Algae:= "Green"]
dt_long_sp_Mesostigmatophyceae <- dt_long_sp[Class == "Mesostigmatophyceae", ]
dt_long_sp_Mesostigmatophyceae[,Algae:= "Green"]
dt_long_sp_Zygnemophyceae <- dt_long_sp[Class == "Zygnemophyceae", ]
dt_long_sp_Zygnemophyceae[,Algae:= "Green"]
dt_long_sp_Coleochaetophyceae <- dt_long_sp[Class == "Coleochaetophyceae", ]
dt_long_sp_Coleochaetophyceae[,Algae:= "Green"]
dt_long_sp_Chlorophytes <- rbind(dt_long_sp_Chlorophyta,dt_long_sp_Mesostigmatophyceae, 
                          dt_long_sp_Zygnemophyceae, dt_long_sp_Coleochaetophyceae)
write.csv(dt_long_sp_Chlorophytes, file= "dt_long_sp_Chlorophytes.csv")
#the green group (Chlorophytes)

dt_long_sp_Bacillariophyta <- dt_long_sp[Class == "Bacillariophyta", ]
dt_long_sp_Bacillariophyta[,Algae:= "Chromophytes"]
dt_long_sp_Dinoflagellata <- dt_long_sp[Class == "Dinoflagellata", ]
dt_long_sp_Dinoflagellata[,Algae:= "Chromophytes"]
dt_long_sp_Chrysophyceae <- dt_long_sp[Class == "Chrysophyceae", ]
dt_long_sp_Chrysophyceae[,Algae:= "Chromophytes"]
dt_long_sp_Chromophytes <- rbind(dt_long_sp_Bacillariophyta,dt_long_sp_Dinoflagellata,dt_long_sp_Chrysophyceae)
write.csv(dt_long_sp_Chromophytes, file= "dt_long_sp_Chromophytes.csv")
#the brown group (Chromophytes, which includes Chrysophytes, Diatoms and Dinoflagellates) 

dt_long_sp_Cryptophyta <- dt_long_sp[Division == "Cryptophyta", ]
dt_long_sp_Cryptophyta[,Algae:= "Cryptophytes"]
dt_long_sp_Rhodophyta <- dt_long_sp[Division == "Rhodophyta", ]
dt_long_sp_Rhodophyta[,Algae:= "Cryptophytes"]
dt_long_sp_Cryptophytes <- rbind(dt_long_sp_Cryptophyta,dt_long_sp_Rhodophyta)
write.csv(dt_long_sp_Cryptophytes, file= "dt_long_sp_Cryptophytes.csv")
#mixed- group (Cryptophytes and phycoerythrin-containing algae, Beutler et al., 2002).

dt_long_sp_algae <- rbind(dt_long_sp_Chlorophytes,dt_long_sp_Chromophytes,dt_long_sp_Cryptophytes)
write.csv(dt_long_sp_algae, file= "table_normcount_algae_raw.csv")
dt_long_sp_algae_sum <- dt_long_sp_algae[, sum(normalised_reads), by=list(Sample.date, sampling.week, sampling.occasion, Mesocosm, treatment, Heated,Algae)]
setnames(dt_long_sp_algae_sum, "V1", "Norm_count")

#remove ASVs with zero reads
dt_long_sp_algae_noz <- dt_long_sp_algae[(normalised_reads !=0),]
#calculate the number of ASVs in each sample
dt_long_sp_algae_count <- dt_long_sp_algae_noz[, .(count = length(unique(ASV))), by=list(Sample.date, sampling.week, sampling.occasion, Mesocosm, treatment, Heated, Algae)]
#calculate the number of ASVs in each mesocosm (sum sampling dates)
dt_long_sp_algae_count_mesocosm <- dt_long_sp_algae_noz[, .(count = length(unique(ASV))), by=list(Mesocosm, Heated, Algae)]
#calculate the number of ASVs in each algal group (sum sampling dates and mesocosm) 
dt_long_sp_algae_count_algae <- dt_long_sp_algae_noz[, .(count = length(unique(ASV))), by=list(Algae)]

#add 16S data
asv_table16S1 <- fread("~/Library/CloudStorage/GoogleDrive-romana.salis@gmail.com/My Drive/LUND/Limnoscenes/Salis16S/asv_table_eDNA1_16S_noNEG_vs.csv")
tax_table16S <- fread("~/Library/CloudStorage/GoogleDrive-romana.salis@gmail.com/My Drive/LUND/Limnoscenes/Salis16S/tax_table_eDNA1_16S_noNEG_vs.csv",header=TRUE)
sample_data16S <- fread("~/Library/CloudStorage/GoogleDrive-romana.salis@gmail.com/My Drive/LUND/Limnoscenes/Salis16S/sample_data_eDNA1_16S_noNEG_vs.csv")
asv_table16S <- dcast(melt(asv_table16S1, id.vars = "V1"), variable ~ V1)
dt_ASVsamp16S <- sample_data16S[asv_table16S, on = .(V1 = variable)]
#convert to long format and retain only the variables want to keep
dt_long16S <- data.table::melt(dt_ASVsamp16S,
                   id.vars = c("Sample.date","sampling.week","sampling.occasion","treatment","Mesocosm","Sample.number","rep","Heated","Invasion"),
                   measure.vars = patterns("^ASV"),
                   variable.name = "ASV",
                   value.name = c("normalised_reads"))
#add taxonomy information
dt_long_sp_16S <- dt_long16S[tax_table16S, on = .(ASV = V1)]
dt_long_sp_16S_Cyanobacteria <- dt_long_sp_16S[Phylum == "Cyanobacteria", ]
dt_long_sp_16S_Cyanobacteria[,Algae:= "Cyanobacteria"]
write.csv(dt_long_sp_16S_Cyanobacteria, file= "dt_long_sp_16S_Cyanobacteria.csv")

dt_long_sp_16S_Cyanobacteria_sum <- dt_long_sp_16S_Cyanobacteria[, sum(normalised_reads), by=list(Sample.date,sampling.week, sampling.occasion, Mesocosm, treatment,Heated,Algae)]
setnames(dt_long_sp_16S_Cyanobacteria_sum, "V1", "Norm_count")

#remove ASVs with zero reads
dt_long_sp_16S_Cyanobacteria_noz <- dt_long_sp_16S_Cyanobacteria[(normalised_reads !=0),]
#calculate the number of ASVs in each sample
dt_long_sp_16S_Cyanobacteria_count <- dt_long_sp_16S_Cyanobacteria_noz[, .(count=length(unique(ASV))), by=list(Sample.date, sampling.week, sampling.occasion, Mesocosm, treatment, Heated, Algae)]
#calculate the number of ASVs in each mesocosm (sum sampling dates)
dt_long_sp_16S_Cyanobacteria_count_mesocosm <- dt_long_sp_16S_Cyanobacteria_noz[, .(count=length(unique(ASV))), by=list(Mesocosm, Heated, Algae)]
#calculate the number of ASVs in each algal group (sum sampling dates and mesocosm)
dt_long_sp_16S_Cyanobacteria_count_algae <- dt_long_sp_16S_Cyanobacteria_noz[, .(count=length(unique(ASV))), by=list(Algae)]


dt_Algae_16S_18S_sum <- rbind(dt_long_sp_16S_Cyanobacteria_sum, dt_long_sp_algae_sum)
write.csv(dt_Algae_16S_18S_sum, file= "table_normcount_algae_sum.csv")

dt_Algae_16S_18S_count <- rbind(dt_long_sp_16S_Cyanobacteria_count, dt_long_sp_algae_count)
write.csv(dt_Algae_16S_18S_count, file= "table_normcount_algae_count.csv")

dt_Algae_16S_18S_count_mesocosm <- rbind(dt_long_sp_16S_Cyanobacteria_count_mesocosm, dt_long_sp_algae_count_mesocosm)


dt_long_sp_algae[,Supergroup:=NULL]
setnames(dt_long_sp_algae, "Division", "Phylum")
dt_Algae_16S_18S <- rbind(dt_long_sp_16S_Cyanobacteria, dt_long_sp_algae)
write.csv(dt_Algae_16S_18S, file= "table_normcount_algae.csv")


#import ALA data to compare
dt_ALA_averages <- fread("~/Library/CloudStorage/GoogleDrive-romana.salis@gmail.com/My Drive/LUND/Limnoscenes/methods_ms/ALA1_averages.csv")
dt_temps <- fread("~/Library/CloudStorage/GoogleDrive-romana.salis@gmail.com/My Drive/LUND/Limnoscenes/methods_ms/temp_weeks_39.csv")
dt_temps_long <- data.table::melt(dt_temps,
                                  id.vars = c("sampling.week","AveTempC","AveTempH","Ave_tempdiff","DegreeWeeks","Ave_tempdiff_R"),
                                  measure.vars = patterns("^M"),
                                  variable.name = "Mesocosm",
                                  value.name = c("Average_temperature"))
dt_temps_long$Mesocosm<-gsub("M","",as.character(dt_temps_long$Mesocosm))
dt_temps_long$Mesocosm<-as.integer(dt_temps_long$Mesocosm)
dt_ALA_averagest <- dt_ALA_averages[dt_temps_long, on = .(sampling.week = sampling.week, Mesocosm =Mesocosm)]
dt_ave_long <- data.table::melt(dt_ALA_averagest[treatment!="NEG"],
                    id.vars = c("Sample.date","sampling.week","sampling.week.end","treatment","Heated","Invasion","Mesocosm","rep", "AveTempH","AveTempC","Ave_tempdiff","DegreeWeeks","Average_temperature"),
                    measure.vars = patterns("^A_"),
                    variable.name = "Algae",
                    value.name = c("concentration"))
dt_ave_long$Algae<-gsub("A_","",as.character(dt_ave_long$Algae))
dt_Algae_all_matching1 <- dt_ave_long[dt_Algae_16S_18S_sum , on = c(Algae = "Algae", 
                                                               Mesocosm ="Mesocosm", Sample.date = "Sample.date", treatment="treatment", Heated="Heated")]
write.csv(dt_Algae_all_matching1, file= "dt_Algae_all_matching.csv")

dt_Algae_all_matching_count1 <- dt_Algae_all_matching1[dt_Algae_16S_18S_count , on = c(Algae = "Algae", 
                                                                   Mesocosm ="Mesocosm", Sample.date = "Sample.date", treatment="treatment", Heated="Heated")]
write.csv(dt_Algae_all_matching_count1, file= "dt_Algae_all_matching_count.csv")


dt_Algae_all_matching1$Norm_count[dt_Algae_all_matching1$Norm_count == 0] <- NA
dt_Algae_all_matching1$concentration[dt_Algae_all_matching1$concentration == 0] <- NA
dt_Algae_all_matching_count1$count[dt_Algae_all_matching_count1$count == 0] <- NA
dt_Algae_all_matching_count1$concentration[dt_Algae_all_matching_count1$concentration == 0] <- NA

dt_Algae_all_matching<- na.omit(dt_Algae_all_matching1)
dt_Algae_all_matching_count<- na.omit(dt_Algae_all_matching_count1)
#n=376

dt_Algae_all_matching$Norm_count_round <- round(dt_Algae_all_matching$Norm_count)
dt_Algae_all_matching$log_conc <- log(dt_Algae_all_matching$concentration)
dt_Algae_all_matching$log_norm <- log(dt_Algae_all_matching$Norm_count)
dt_Algae_all_matching_count$log_conc <- log(dt_Algae_all_matching_count$concentration)
dt_Algae_all_matching_count$log_ASVcount <- log(dt_Algae_all_matching_count$count)

dt_Algae_all_matching = dt_Algae_all_matching[sampling.week != "0", ]
dt_Algae_all_matching_count = dt_Algae_all_matching_count1[sampling.week != "0", ]

dt_Algae_all_matching_cyano = dt_Algae_all_matching[Algae == "Cyanobacteria", ]
dt_Algae_all_matching_chrom = dt_Algae_all_matching[Algae == "Chromophytes", ]
dt_Algae_all_matching_green = dt_Algae_all_matching[Algae == "Green", ]
dt_Algae_all_matching_crypt = dt_Algae_all_matching[Algae == "Cryptophytes", ]

dt_Algae_all_matching_count_cyano = dt_Algae_all_matching_count0[Algae == "Cyanobacteria", ]
dt_Algae_all_matching_count_chrom = dt_Algae_all_matching_count0[Algae == "Chromophytes", ]
dt_Algae_all_matching_count_green = dt_Algae_all_matching_count0[Algae == "Green", ]
dt_Algae_all_matching_count_crypt = dt_Algae_all_matching_count0[Algae == "Cryptophytes", ]

pdf("hist_Norm_count_round.pdf", width = 7, height = 5)
hist(dt_Algae_all_matching$Norm_count_round, main = "",
     xlab = "Normalised read count", col = "lightgray", border = "black")
dev.off()
pdf("hist_Norm_count_cyano.pdf", width = 7, height = 5)
hist(dt_Algae_all_matching_cyano$Norm_count_round, main = "",
     xlab = "Normalised read count", col = "lightgray", border = "black")
dev.off()
pdf("hist_Norm_count_chrom.pdf", width = 7, height = 5)
hist(dt_Algae_all_matching_chrom$Norm_count_round, main = "",
     xlab = "Normalised read count", col = "lightgray", border = "black")
dev.off()
pdf("hist_Norm_count_cyano.pdf", width = 7, height = 5)
hist(dt_Algae_all_matching_cyano$Norm_count_round, main = "",
     xlab = "Normalised read count", col = "lightgray", border = "black")
dev.off()
hist(dt_Algae_all_matching_chrom$Norm_count)
#hist(dt_Algae_all_matching_chrom$log_norm)
hist(dt_Algae_all_matching_crypt$Norm_count)
#hist(dt_Algae_all_matching_crypt$log_norm)
hist(dt_Algae_all_matching_green$Norm_count)
#hist(dt_Algae_all_matching_green$log_norm)
hist(dt_Algae_all_matching$concentration)
#hist(dt_Algae_all_matching$log_conc)
hist(dt_Algae_all_matching_cyano$concentration)
#hist(dt_Algae_all_matching_cyano$log_conc)
hist(dt_Algae_all_matching_chrom$concentration)
#hist(dt_Algae_all_matching_chrom$log_conc)
hist(dt_Algae_all_matching_crypt$concentration)
#hist(dt_Algae_all_matching_crypt$log_conc)
hist(dt_Algae_all_matching_green$concentration)
#hist(dt_Algae_all_matching_green$log_conc)

#now load full data - NOT ONLY samples where the algal groups were detected by both methods
#dt_Algae_all <- fread("dt_Algae_all_matching.csv")
#dt_Algae_count <- fread("dt_Algae_all_matching_count.csv")
dt_Algae_all <- dt_Algae_all_matching_count1[sampling.occasion!="0"]
#transform factors
dt_Algae_all$sampling.occasion<- as.factor(dt_Algae_all$sampling.occasion)
dt_Algae_all$Mesocosm<- as.factor(dt_Algae_all$Mesocosm)
dt_Algae_all$Heated<- factor(dt_Algae_all$Heated, levels=c("N","Y"))
dt_Algae_all$Mesocosm<- factor(dt_Algae_all$Mesocosm,levels=c("2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","19","20","21","23","24","25","26","27"))
dt_Algae_all$sampling.occasion<- factor(dt_Algae_all$sampling.occasion,levels=c("1","2","3","4","5"))

#calculate relative abundance
dt_Algae_all[is.na(dt_Algae_all)] <- 0
RelAbundConc <- unlist(lapply(split(dt_Algae_all, list(dt_Algae_all$sampling.week,dt_Algae_all$Mesocosm)), function(x){x$concentration / sum(x$concentration)}))
RelAbundConc<- as.data.frame(RelAbundConc)
RelAbundConc$MesoDate <- rownames(RelAbundConc)
RelAbundConc<-as.data.table(RelAbundConc %>% separate(MesoDate, c( 'sampling.week','Mesocosm')))
RelAbundConc$group = as.character(lapply(strsplit(as.character(RelAbundConc$Mesocosm), split=""),
                                         tail, n=1))
RelAbundConc$Mesocosm <- gsub('.{1}$', '', RelAbundConc$Mesocosm)
RelAbundConc$Mesocosm<- factor(RelAbundConc$Mesocosm,levels=c("2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","19","20","21","23","24","25","26","27"))
RelAbundConc$sampling.week<- factor(RelAbundConc$sampling.week,levels=c("10","12","14","18","22"))

RelAbundReads <- unlist(lapply(split(dt_Algae_all, list(dt_Algae_all$sampling.week,dt_Algae_all$Mesocosm)), function(x){x$Norm_count / sum(x$Norm_count)}))
RelAbundReads<- as.data.frame(RelAbundReads)
RelAbundReads[is.na(RelAbundReads)] <- 0
RelAbundReads$MesoDate <- rownames(RelAbundReads)
RelAbundReads<-as.data.table(RelAbundReads %>% separate(MesoDate, c( 'sampling.week','Mesocosm')))
RelAbundReads$group = as.character(lapply(strsplit(as.character(RelAbundReads$Mesocosm), split=""),
                                          tail, n=1))
RelAbundReads$Mesocosm <- gsub('.{1}$', '', RelAbundReads$Mesocosm)
RelAbundReads$Mesocosm<- factor(RelAbundReads$Mesocosm,levels=c("2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","19","20","21","23","24","25","26","27"))
RelAbundReads$sampling.week<- factor(RelAbundReads$sampling.week,levels=c("10","12","14","18","22"))

RelAbundASVs <- unlist(lapply(split(dt_Algae_count, list(dt_Algae_count$sampling.week,dt_Algae_count$Mesocosm)), function(x){x$count / sum(x$count)}),use.names=TRUE)
RelAbundASVs<- as.data.frame(RelAbundASVs)
RelAbundASVs[is.na(RelAbundASVs)] <- 0
RelAbundASVs$MesoDate <- rownames(RelAbundASVs)
RelAbundASVs<-as.data.table(RelAbundASVs %>% separate(MesoDate, c( 'sampling.week','Mesocosm')))
RelAbundASVs$group = as.character(lapply(strsplit(as.character(RelAbundASVs$Mesocosm), split=""),
                                         tail, n=1))
RelAbundASVs$Mesocosm <- gsub('.{1}$', '', RelAbundASVs$Mesocosm)
RelAbundASVs$Mesocosm<- factor(RelAbundASVs$Mesocosm,levels=c("2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","19","20","21","23","24","25","26","27"))
RelAbundASVs$sampling.week<- factor(RelAbundASVs$sampling.week,levels=c("10","12","14","18","22"))

keycols<-c("sampling.week","Mesocosm")
setkeyv(RelAbundConc,keycols)
setkeyv(RelAbundReads,keycols)
setkeyv(RelAbundASVs,keycols)
setkeyv(dt_Algae_all,keycols)

dt_Algae_all_relabund <- cbind(dt_Algae_all, RelAbundConc$RelAbundConc, RelAbundReads$RelAbundReads, RelAbundASVs$RelAbundASVs)   
setnames(dt_Algae_all_relabund, "V2", "RelAbundConc")
setnames(dt_Algae_all_relabund, "V3", "RelAbundReads")
setnames(dt_Algae_all_relabund, "V4", "RelAbundASVs")

write.csv(dt_Algae_all_relabund, file= "dt_Algae_all2.csv")

colnames(dt_Algae_all_relabund)[colnames(dt_Algae_all_relabund) == "RelAbundConc"] <- "ALA"
colnames(dt_Algae_all_relabund)[colnames(dt_Algae_all_relabund) == "RelAbundReads"] <- "Reads"
colnames(dt_Algae_all_relabund)[colnames(dt_Algae_all_relabund) == "RelAbundASVs"] <- "ASVs"

sum_ALA_chrom <- sum(dt_Algae_all_relabund[Algae=="Chromophytes"]$ALA)
sum_ALA_cyano <- sum(dt_Algae_all_relabund[Algae=="Cyanobacteria"]$ALA)
sum_ALA_green <- sum(dt_Algae_all_relabund[Algae=="Green"]$ALA)
sum_ALA_crypt <- sum(dt_Algae_all_relabund[Algae=="Cryptophytes"]$ALA) 
ALA_total <- sum_ALA_chrom + sum_ALA_crypt + sum_ALA_cyano + sum_ALA_green
sum_ALA_chrom/ALA_total*100
sum_ALA_cyano/ALA_total*100
sum_ALA_green/ALA_total*100
sum_ALA_crypt/ALA_total*100

sum_Reads_chrom <- sum(dt_Algae_all_relabund[Algae=="Chromophytes"]$Reads)
sum_Reads_cyano <- sum(dt_Algae_all_relabund[Algae=="Cyanobacteria"]$Reads)
sum_Reads_green <- sum(dt_Algae_all_relabund[Algae=="Green"]$Reads)
sum_Reads_crypt <- sum(dt_Algae_all_relabund[Algae=="Cryptophytes"]$Reads) 
Reads_total <- sum_Reads_chrom + sum_Reads_crypt + sum_Reads_cyano + sum_Reads_green
sum_Reads_chrom/Reads_total*100
sum_Reads_cyano/Reads_total*100
sum_Reads_green/Reads_total*100
sum_Reads_crypt/Reads_total*100

sum_ALA_chrom <- sum(dt_Algae_all_relabund[Algae=="Chromophytes"]$ALA)
sum_ALA_cyano <- sum(dt_Algae_all_relabund[Algae=="Cyanobacteria"]$ALA)
sum_ALA_green <- sum(dt_Algae_all_relabund[Algae=="Green"]$ALA)
sum_ALA_crypt <- sum(dt_Algae_all_relabund[Algae=="Cryptophytes"]$ALA) 
ALA_total <- sum_ALA_chrom + sum_ALA_crypt + sum_ALA_cyano + sum_ALA_green
sum_ALA_chrom/ALA_total*100
sum_ALA_cyano/ALA_total*100
sum_ALA_green/ALA_total*100
sum_ALA_crypt/ALA_total*100
n <- nrow(dt_Algae_all_relabund)
dt_Algae_all_relabund$ALA_adj <- (dt_Algae_all_relabund$ALA * (n - 1) + 0.5) / n

#dt_Algae_all_relabund1 <- dt_Algae_all_relabund
#dt_Algae_all_relabund1$Reads[dt_Algae_all_relabund1$Reads == 0] <- NA
#dt_Algae_all_relabund1$ALA[dt_Algae_all_relabund1$ALA == 0] <- NA
#dt_Algae_all_relabund1$ASVs[dt_Algae_all_relabund1$ASVs == 0] <- NA
#dt_Algae_all_relabund_nz<- na.omit(dt_Algae_all_relabund1)
#dt_Algae_all_relabund_cyano = dt_Algae_all_relabund_nz[Algae == "Cyanobacteria", ]
#dt_Algae_all_relabund_chrom = dt_Algae_all_relabund_nz[Algae == "Chromophytes", ]
#dt_Algae_all_relabund_green = dt_Algae_all_relabund_nz[Algae == "Green", ]
#dt_Algae_all_relabund_crypt = dt_Algae_all_relabund_nz[Algae == "Cryptophytes", ]
dt_Algae_all_relabund_cyano = dt_Algae_all_relabund[Algae == "Cyanobacteria", ]
dt_Algae_all_relabund_chrom = dt_Algae_all_relabund[Algae == "Chromophytes", ]
dt_Algae_all_relabund_green = dt_Algae_all_relabund[Algae == "Green", ]
dt_Algae_all_relabund_crypt = dt_Algae_all_relabund[Algae == "Cryptophytes", ]



## GLMMs
##histograms of data distributionz
library(gridExtra)  # For arranging multiple plots

# Function to create histograms with density overlay
plot_histogram <- function(data, var, xlabel, log_x = FALSE) {
  p <- ggplot(data, aes(x = .data[[var]])) +
    geom_histogram(aes(y = ..density..), fill = "grey70", color = "black", bins = 30) +
    geom_density(colour = "red", linewidth = 1, adjust = 1.5) +
    labs(x = xlabel, y = "Density") +
    theme_classic(base_size = 14) +
    theme(axis.title.y = element_text(size = 14), 
          axis.title.x = element_text(size = 14),
          plot.margin = margin(5, 5, 5, 5))  # Adjust margins for consistency
  
  # Log-transform x-axis for log-normal distribution
  if (log_x) {
    p <- p + scale_x_log10()
  }
  return(p)
}

# Create individual plots
p1 <- plot_histogram(dt_Algae_all_matching, "concentration", "ALA biomass concentration (µg/L)", log_x = TRUE) + ggtitle("A")
p2 <- plot_histogram(dt_Algae_all_matching, "Norm_count_round", "Normalised read count") + ggtitle("C")
p3 <- plot_histogram(dt_Algae_all_matching_count0, "count", "Number of ASVs") + ggtitle("E")
p4 <- plot_histogram(dt_Algae_all_relabund, "ALA", "Proportion of biomass") + ggtitle("B")
p5 <- plot_histogram(dt_Algae_all_relabund, "Reads", "Proportion of reads") + ggtitle("D")
p6 <- plot_histogram(dt_Algae_all_relabund, "ASVs", "Proportion of ASVs") + ggtitle("F")

# Arrange plots in a grid
grid.arrange(p1, p4, p2, p5, p3, p6, ncol = 2)

ggsave("p_histogram.jpeg", width = 5, height = 8)
ggsave("p_histogram.pdf", width = 5, height = 8)

mean(dt_Algae_all_matching$Norm_count_round)
var(dt_Algae_all_matching$Norm_count_round) # var >> mean -> overdispersion
mean(dt_Algae_all_matching_count0$count)
var(dt_Algae_all_matching_count0$count) # var >> mean -> overdispersion

# Load required package
library(glmmTMB)
library(car)
library(DHARMa)
library(performance)

## test treatment effects - normalised counts
glmm_poisson_treat <- glmmTMB(Norm_count_round ~ Heated*Invasion + (1|Mesocosm) + (1|sampling.week) + (1|Algae), 
                              data = dt_Algae_all_matching,
                              family = poisson)
performance::check_singularity(glmm_poisson_treat)
plot(simulateResiduals(fittedModel = glmm_poisson_treat, plot = F))# best fitting
summary(glmm_poisson_treat)
Anova(glmm_poisson_treat, type = "III")

glmm_nb_treat <- glmmTMB(Norm_count_round ~ Heated*Invasion + (1|Mesocosm) + (1|sampling.week)+ (1|Algae), 
                         data = dt_Algae_all_matching,
                         family = nbinom2)
performance::check_singularity(glmm_nb_treat)
plot(simulateResiduals(fittedModel = glmm_nb_treat, plot = F))# best fitting
summary(glmm_nb_treat)
Anova(glmm_nb_treat, type = "III")

glmm_nb_treat_crypt <- glmmTMB(Norm_count_round ~ Heated*Invasion+ (1|Mesocosm) + (1|sampling.week) , 
                         data = dt_Algae_all_matching_crypt,
                         family = nbinom2)
performance::check_singularity(glmm_nb_treat_crypt)# + (1|Mesocosm) + (1|sampling.week) causing singularity
plot(simulateResiduals(fittedModel = glmm_nb_treat_crypt, plot = F))# best fitting
summary(glmm_nb_treat_crypt)
Anova(glmm_nb_treat_crypt, type = "III")# significant invasion
r2_nakagawa(glmm_nb_treat_crypt)

glmm_nb_treat_chrom <- glmmTMB(Norm_count_round ~ Heated*Invasion + (1|Mesocosm) + (1|sampling.week), 
                               data = dt_Algae_all_matching_chrom,
                               family = nbinom2)
performance::check_singularity(glmm_nb_treat_chrom)
plot(simulateResiduals(fittedModel = glmm_nb_treat_chrom, plot = F))# best fitting
summary(glmm_nb_treat_chrom)
Anova(glmm_nb_treat_chrom, type = "III")

glmm_nb_treat_cyano <- glmmTMB(Norm_count_round ~ Heated*Invasion + (1|Mesocosm) + (1|sampling.week), 
                               data = dt_Algae_all_matching_cyano,
                               family = nbinom2)
performance::check_singularity(glmm_nb_treat_cyano)
plot(simulateResiduals(fittedModel = glmm_nb_treat_cyano, plot = F))# best fitting
summary(glmm_nb_treat_cyano)
Anova(glmm_nb_treat_cyano, type = "III")

glmm_nb_treat_green <- glmmTMB(Norm_count_round ~ Heated*Invasion + (1|Mesocosm)  , 
                               data = dt_Algae_all_matching_green,
                               family = nbinom2)
performance::check_singularity(glmm_nb_treat_green)# + (1|sampling.week) causing singularity
plot(simulateResiduals(fittedModel = glmm_nb_treat_green, plot = F))# best fitting
summary(glmm_nb_treat_green)
Anova(glmm_nb_treat_green, type = "III")

# treatment number ASVs
glmm_poisson_treat_ASV <- glmmTMB(count ~ Heated*Invasion + (1|Mesocosm) + (1|sampling.week) + (1|Algae), 
                              data = dt_Algae_all_matching_count0,
                              family = poisson)
performance::check_singularity(glmm_poisson_treat_ASV)
plot(simulateResiduals(fittedModel = glmm_poisson_treat_ASV, plot = F))# best fitting
summary(glmm_poisson_treat_ASV)
Anova(glmm_poisson_treat_ASV, type = "III")

glmm_nb_treat_ASV <- glmmTMB(count ~ Heated*Invasion + (1|Mesocosm) + (1|sampling.week)+ (1|Algae), 
                         data = dt_Algae_all_matching_count0,
                         family = nbinom2)
performance::check_singularity(glmm_nb_treat_ASV)
plot(simulateResiduals(fittedModel = glmm_nb_treat_ASV, plot = F))# best fitting
summary(glmm_nb_treat_ASV)
Anova(glmm_nb_treat_ASV, type = "III")

glmm_nb_treat_crypt_ASV <- glmmTMB(count ~ Heated*Invasion + (1|Mesocosm) + (1|sampling.week), 
                               data = dt_Algae_all_matching_count_crypt,
                               family = nbinom2)
performance::check_singularity(glmm_nb_treat_crypt_ASV)
plot(simulateResiduals(fittedModel = glmm_nb_treat_crypt_ASV, plot = F))# best fitting
summary(glmm_nb_treat_crypt_ASV)
Anova(glmm_nb_treat_crypt_ASV, type = "III")# significant invasion
r2_nakagawa(glmm_nb_treat_crypt_ASV)

glmm_nb_treat_chrom_ASV <- glmmTMB(count ~ Heated*Invasion + (1|Mesocosm) + (1|sampling.week), 
                               data = dt_Algae_all_matching_count_chrom,
                               family = nbinom2)
performance::check_singularity(glmm_nb_treat_chrom_ASV)
plot(simulateResiduals(fittedModel = glmm_nb_treat_chrom_ASV, plot = F))# best fitting
summary(glmm_nb_treat_chrom_ASV)
Anova(glmm_nb_treat_chrom_ASV, type = "III")

glmm_nb_treat_cyano_ASV <- glmmTMB(count ~ Heated*Invasion + (1|Mesocosm) + (1|sampling.week), 
                               data = dt_Algae_all_matching_count_cyano,
                               family = nbinom2)
performance::check_singularity(glmm_nb_treat_cyano_ASV)
plot(simulateResiduals(fittedModel = glmm_nb_treat_cyano_ASV, plot = F))# best fitting
summary(glmm_nb_treat_cyano_ASV)
Anova(glmm_nb_treat_cyano_ASV, type = "III")

glmm_nb_treat_green_ASV <- glmmTMB(count ~ Heated*Invasion + (1|Mesocosm) , 
                               data = dt_Algae_all_matching_count_green,
                               family = nbinom2)
performance::check_singularity(glmm_nb_treat_green_ASV)# + (1|sampling.week) causing singularity
plot(simulateResiduals(fittedModel = glmm_nb_treat_green_ASV, plot = F))# best fitting
summary(glmm_nb_treat_green_ASV)
Anova(glmm_nb_treat_green_ASV, type = "III")

# treatment ALA concentration
glmm_lognormal_treat_ALA <- glmmTMB(concentration ~ Heated*Invasion + (1|Mesocosm) + (1|sampling.week)+ (1|Algae), 
                             data = dt_Algae_all_matching_count0,
                             family = lognormal)
performance::check_singularity(glmm_lognormal_treat_ALA)
plot(simulateResiduals(fittedModel = glmm_lognormal_treat_ALA, plot = F))# best fitting
summary(glmm_lognormal_treat_ALA)
Anova(glmm_lognormal_treat_ALA, type = "III")

glmm_lognormal_treat_crypt_ALA <- glmmTMB(concentration ~ Heated*Invasion + (1|Mesocosm) + (1|sampling.week), 
                                   data = dt_Algae_all_matching_count_crypt,
                                   family = lognormal)
performance::check_singularity(glmm_lognormal_treat_crypt_ALA)
plot(simulateResiduals(fittedModel = glmm_lognormal_treat_crypt_ALA, plot = F))# best fitting
summary(glmm_lognormal_treat_crypt_ALA)
Anova(glmm_lognormal_treat_crypt_ALA, type = "III")# significant invasion

glmm_lognormal_treat_chrom_ALA <- glmmTMB(concentration ~ Heated*Invasion , 
                                   data = dt_Algae_all_matching_count_chrom,
                                   family = lognormal)
performance::check_singularity(glmm_lognormal_treat_chrom_ALA)#  + (1|Mesocosm) + (1|sampling.week) causing singularity
plot(simulateResiduals(fittedModel = glmm_lognormal_treat_chrom_ALA, plot = F))# best fitting
summary(glmm_lognormal_treat_chrom_ALA)
Anova(glmm_lognormal_treat_chrom_ALA, type = "III")

glmm_lognormal_treat_cyano_ALA <- glmmTMB(concentration ~ Heated*Invasion + (1|Mesocosm) + (1|sampling.week), 
                                   data = dt_Algae_all_matching_count_cyano,
                                   family = lognormal)
performance::check_singularity(glmm_lognormal_treat_cyano_ALA)
plot(simulateResiduals(fittedModel = glmm_lognormal_treat_cyano_ALA, plot = F))# best fitting
summary(glmm_lognormal_treat_cyano_ALA)
Anova(glmm_lognormal_treat_cyano_ALA, type = "III")

glmm_lognormal_treat_green_ALA <- glmmTMB(concentration ~ Heated*Invasion + (1|Mesocosm)+ (1|sampling.week) , 
                                   data = dt_Algae_all_matching_count_green,
                                   family = lognormal)
performance::check_singularity(glmm_lognormal_treat_green_ALA)
plot(simulateResiduals(fittedModel = glmm_lognormal_treat_green_ALA, plot = F))# best fitting
summary(glmm_lognormal_treat_green_ALA)
Anova(glmm_lognormal_treat_green_ALA, type = "III")

# beta reg rel abundance ASV counts
glmm_relab_treat_ASVs_all <- glmmTMB(ASVs ~ Heated*Invasion + (1 | Mesocosm) + (1 | sampling.week)+ (1 | Algae),
                                     data = dt_Algae_all_relabund,
                                     family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_treat_ASVs_all)# + (1 | Mesocosm) + (1 | sampling.week) causing singularity
glmm_relab_treat_ASVs_all <- glmmTMB(ASVs ~ Heated*Invasion + (1 | Algae),
                               data = dt_Algae_all_relabund,
                               family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_treat_ASVs_all)
summary(glmm_relab_treat_ASVs_all)
Anova(glmm_relab_treat_ASVs_all, type = "III")
plot(simulateResiduals(glmm_relab_treat_ASVs_all, plot = F))

glmm_relab_treat_ASVs_crypt <- glmmTMB(ASVs ~ Heated*Invasion   + (1 | Mesocosm) + (1 | sampling.week),
                                 data = dt_Algae_all_relabund_crypt,
                                 family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_treat_ASVs_crypt)
summary(glmm_relab_treat_ASVs_crypt)
Anova(glmm_relab_treat_ASVs_crypt, type = "III")
plot(simulateResiduals(glmm_relab_treat_ASVs_crypt, plot = F))

glmm_relab_treat_ASVs_cyano <- glmmTMB(ASVs ~ Heated*Invasion   + (1 | Mesocosm) + (1 | sampling.week),
                                 data = dt_Algae_all_relabund_cyano,
                                 family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_treat_ASVs_cyano)
summary(glmm_relab_treat_ASVs_cyano)
Anova(glmm_relab_treat_ASVs_cyano, type = "III")
plot(simulateResiduals(glmm_relab_treat_ASVs_cyano, plot = F))

glmm_relab_treat_ASVs_chrom <- glmmTMB(ASVs ~ Heated*Invasion   + (1 | Mesocosm) + (1 | sampling.week),
                                 data = dt_Algae_all_relabund_chrom,
                                 family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_treat_ASVs_chrom)
summary(glmm_relab_treat_ASVs_chrom)
Anova(glmm_relab_treat_ASVs_chrom, type = "III")
r2_nakagawa(glmm_relab_treat_ASVs_chrom)
plot(simulateResiduals(glmm_relab_treat_ASVs_chrom, plot = F))

glmm_relab_treat_ASVs_green <- glmmTMB(ASVs ~ Heated*Invasion   + (1 | Mesocosm) + (1 | sampling.week),
                                 data = dt_Algae_all_relabund_green,
                                 family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_treat_ASVs_green)
summary(glmm_relab_treat_ASVs_green)
Anova(glmm_relab_treat_ASVs_green, type = "III")
r2_nakagawa(glmm_relab_treat_ASVs_green)
plot(simulateResiduals(glmm_relab_treat_ASVs_green, plot = F))

# beta reg rel abundance reads
glmm_relab_reads_all <- glmmTMB(Reads ~ Heated*Invasion + (1 | Algae) + (1 | Mesocosm) + (1 | sampling.week),
                                data = dt_Algae_all_relabund,
                                family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_reads_all)# + (1 | Mesocosm) + (1 | sampling.week) causing singularity
# is singular
glmm_relab_reads_all <- glmmTMB(Reads ~ Heated*Invasion + (1 | Algae),
                                data = dt_Algae_all_relabund,
                                family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_reads_all)
summary(glmm_relab_reads_all)
Anova(glmm_relab_reads_all, type = "III")
r2_nakagawa(glmm_relab_reads_all)
plot(simulateResiduals(glmm_relab_reads_all, plot = F))

hist(dt_Algae_all_relabund_crypt$Reads)
glmm_relab_treat_ASVs_crypt <- glmmTMB(Reads ~ Heated*Invasion + (1 | Mesocosm) + (1 | sampling.week),
                                 data = dt_Algae_all_relabund_crypt,
                                 family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_treat_ASVs_crypt)
summary(glmm_relab_treat_ASVs_crypt)
Anova(glmm_relab_treat_ASVs_crypt, type = "III")
r2_nakagawa(glmm_relab_treat_ASVs_crypt)
plot(simulateResiduals(glmm_relab_treat_ASVs_crypt, plot = F))

glmm_relab_treat_ASVs_cyano <- glmmTMB(Reads ~ Heated*Invasion + (1 | Mesocosm) + (1 | sampling.week),
                                 data = dt_Algae_all_relabund_cyano,
                                 family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_treat_ASVs_cyano)
summary(glmm_relab_treat_ASVs_cyano)
Anova(glmm_relab_treat_ASVs_cyano, type = "III")
r2_nakagawa(glmm_relab_treat_ASVs_cyano)
plot(simulateResiduals(glmm_relab_treat_ASVs_cyano, plot = F))

glmm_relab_treat_ASVs_chrom <- glmmTMB(Reads ~ Heated*Invasion + (1 | Mesocosm) + (1 | sampling.week),
                                 data = dt_Algae_all_relabund_chrom,
                                 family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_treat_ASVs_chrom)
summary(glmm_relab_treat_ASVs_chrom)
Anova(glmm_relab_treat_ASVs_chrom, type = "III")
r2_nakagawa(glmm_relab_treat_ASVs_chrom)
plot(simulateResiduals(glmm_relab_treat_ASVs_chrom, plot = F))

glmm_relab_treat_ASVs_green <- glmmTMB(Reads ~ Heated*Invasion + (1 | Mesocosm) + (1 | sampling.week),
                                 data = dt_Algae_all_relabund_green,
                                 family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_treat_ASVs_green)
summary(glmm_relab_treat_ASVs_green)
Anova(glmm_relab_treat_ASVs_green, type = "III")
r2_nakagawa(glmm_relab_treat_ASVs_green)
plot(simulateResiduals(glmm_relab_treat_ASVs_green, plot = F))

# beta reg rel abundance conc
glmm_relab_ALA_all <- glmmTMB(ALA_adj ~ Heated*Invasion + (1 | Algae) + (1 | Mesocosm) + (1 | sampling.week),
                                data = dt_Algae_all_relabund,
                                family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_ALA_all)# + (1 | Mesocosm) + (1 | sampling.week) causing singularity
# is singular
glmm_relab_ALA_all <- glmmTMB(ALA_adj ~ Heated*Invasion + (1 | Algae),
                                data = dt_Algae_all_relabund,
                                family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_ALA_all)
summary(glmm_relab_ALA_all)
Anova(glmm_relab_ALA_all, type = "III")
r2_nakagawa(glmm_relab_ALA_all)
plot(simulateResiduals(glmm_relab_ALA_all, plot = F))

hist(dt_Algae_all_relabund_crypt$ALA_adj)
glmm_relab_treat_ALA_crypt <- glmmTMB(ALA_adj ~ Heated*Invasion ,
                                       data = dt_Algae_all_relabund_crypt,
                                       family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_treat_ALA_crypt)#+ (1 | Mesocosm) + (1 | sampling.week)
summary(glmm_relab_treat_ALA_crypt)
Anova(glmm_relab_treat_ALA_crypt, type = "III")
plot(simulateResiduals(glmm_relab_treat_ALA_crypt, plot = F))

glmm_relab_treat_ALA_cyano <- glmmTMB(ALA_adj ~ Heated*Invasion + (1 | Mesocosm) ,
                                       data = dt_Algae_all_relabund_cyano,
                                       family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_treat_ALA_cyano)# + (1 | sampling.week) causing singularity
summary(glmm_relab_treat_ALA_cyano)
Anova(glmm_relab_treat_ALA_cyano, type = "III")
r2_nakagawa(glmm_relab_treat_ALA_cyano)
plot(simulateResiduals(glmm_relab_treat_ALA_cyano, plot = F))

glmm_relab_treat_ALA_chrom <- glmmTMB(ALA_adj ~ Heated*Invasion + (1 | Mesocosm) + (1 | sampling.week),
                                       data = dt_Algae_all_relabund_chrom,
                                       family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_treat_ALA_chrom)
summary(glmm_relab_treat_ALA_chrom)
Anova(glmm_relab_treat_ALA_chrom, type = "III")
r2_nakagawa(glmm_relab_treat_ALA_chrom)
plot(simulateResiduals(glmm_relab_treat_ALA_chrom, plot = F))

glmm_relab_treat_ALA_green <- glmmTMB(ALA_adj ~ Heated*Invasion + (1 | Mesocosm) ,
                                       data = dt_Algae_all_relabund_green,
                                       family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_treat_ALA_green)# + (1 | sampling.week) causing singularity
summary(glmm_relab_treat_ALA_green)
Anova(glmm_relab_treat_ALA_green, type = "III")
r2_nakagawa(glmm_relab_treat_ALA_green)
plot(simulateResiduals(glmm_relab_treat_ALA_green, plot = F))
glmm_relab_treat_ALA_green

######
## models with out testing treatments
# Poisson GLMM
glmm_poisson <- glmmTMB(Norm_count_round ~ log_conc + (1 | Mesocosm) + (1 | sampling.week)+ (1 | Algae),
                        data = dt_Algae_all_matching,
                        family = poisson)
performance::check_singularity(glmm_poisson)
summary(glmm_poisson)
Anova(glmm_poisson, type = "II")
r2_nakagawa(glmm_poisson)
plot(simulateResiduals(glmm_poisson))

# negbinom normalised read counts
glmm_nb_all <- glmmTMB(Norm_count_round ~ log_conc + (1 | Mesocosm) + (1 | sampling.week) + (1 | Algae),
                         data = dt_Algae_all_matching,
                         family = nbinom2)
performance::check_singularity(glmm_nb_all)
summary(glmm_nb_all)
Anova(glmm_nb_all, type = "II")
r2_nakagawa(glmm_nb_all)
plot(simulateResiduals(glmm_nb_all))
# Compare model fits using AIC
AIC(glmm_poisson, glmm_nb_all) # use binomial
glmm_nb_all

#To estimate the variation explained by algal group identity, we fit a simplified model with Algae as the sole random effect and compared marginal and conditional R² values.
glmm_nb_all_algae <- glmmTMB(Norm_count_round ~ log_conc + (1 | Algae),
                       data = dt_Algae_all_matching,
                       family = nbinom2)
performance::check_singularity(glmm_nb_all_algae)
summary(glmm_nb_all_algae)
Anova(glmm_nb_all_algae, type = "II")
r2_nakagawa(glmm_nb_all_algae)
plot(simulateResiduals(glmm_nb_all_algae))

glmm_nb_chrom <- glmmTMB(Norm_count_round ~ log_conc + (1 | Mesocosm) + (1 | sampling.week),
                         data = dt_Algae_all_matching_chrom,
                         family = nbinom2)
performance::check_singularity(glmm_nb_chrom)
summary(glmm_nb_chrom)
Anova(glmm_nb_chrom, type = "II")
r2_nakagawa(glmm_nb_chrom)
plot(simulateResiduals(glmm_nb_chrom))
glmm_nb_chrom

glmm_nb_crypt <- glmmTMB(Norm_count_round ~ log_conc + (1 | Mesocosm),
                         data = dt_Algae_all_matching_crypt,
                         family = nbinom2)
performance::check_singularity(glmm_nb_crypt)# + (1 | sampling.week) causing singularity
summary(glmm_nb_crypt)
Anova(glmm_nb_crypt, type = "II")
r2_nakagawa(glmm_nb_crypt)
plot(simulateResiduals(glmm_nb_crypt))
glmm_nb_crypt

glmm_nb_cyano <- glmmTMB(Norm_count_round ~ log_conc + (1 | Mesocosm) + (1 | sampling.week),
                         data = dt_Algae_all_matching_cyano,
                         family = nbinom2)
performance::check_singularity(glmm_nb_cyano)
summary(glmm_nb_cyano)
Anova(glmm_nb_cyano, type = "II")
r2_nakagawa(glmm_nb_cyano)
plot(simulateResiduals(glmm_nb_cyano))
glmm_nb_cyano

glmm_nb_green <- glmmTMB(Norm_count_round ~ log_conc + (1 | Mesocosm) +(1 | sampling.week),
                         data = dt_Algae_all_matching_green,
                         family = nbinom2)
performance::check_singularity(glmm_nb_green)
summary(glmm_nb_green)
Anova(glmm_nb_green, type = "II")
r2_nakagawa(glmm_nb_green)
plot(simulateResiduals(glmm_nb_green))
glmm_nb_green

# negbinom ASV counts
glmm_nb_ASVs_all <- glmmTMB(count ~ log_conc + (1 | Mesocosm) + (1 | sampling.week)+ (1 | Algae),
                       data = dt_Algae_all_matching_count0,
                       family = nbinom2)
performance::check_singularity(glmm_nb_ASVs_all) 
summary(glmm_nb_ASVs_all)
Anova(glmm_nb_ASVs_all, type = "II")
r2_nakagawa(glmm_nb_ASVs_all)
plot(simulateResiduals(glmm_nb_ASVs_all))
glmm_nb_ASVs_all

#To estimate the variation explained by algal group identity, we fit a simplified model with Algae as the sole random effect and compared marginal and conditional R² values.
glmm_nb_ASVs_all_algae <- glmmTMB(count ~ log_conc + (1 | Algae),
                            data = dt_Algae_all_matching_count0,
                            family = nbinom2)
performance::check_singularity(glmm_nb_ASVs_all_algae) 
summary(glmm_nb_ASVs_all_algae)
Anova(glmm_nb_ASVs_all_algae, type = "II")
r2_nakagawa(glmm_nb_ASVs_all_algae)
plot(simulateResiduals(glmm_nb_ASVs_all_algae))
glmm_nb_ASVs_all_algae

glmm_nb_ASVs_chrom <- glmmTMB(count ~ log_conc + (1 | Mesocosm) + (1 | sampling.week),
                              data = dt_Algae_all_matching_count_chrom,
                              family = nbinom2)
performance::check_singularity(glmm_nb_ASVs_chrom) 
summary(glmm_nb_ASVs_chrom)
Anova(glmm_nb_ASVs_chrom, type = "II")
r2_nakagawa(glmm_nb_ASVs_chrom)
plot(simulateResiduals(glmm_nb_ASVs_chrom))
glmm_nb_ASVs_chrom

glmm_nb_ASVs_crypt <- glmmTMB(count ~ log_conc + (1 | Mesocosm),
                         data = dt_Algae_all_matching_count_crypt,
                         family = nbinom2)
performance::check_singularity(glmm_nb_ASVs_crypt) #+ (1 | sampling.week), led to singularity
summary(glmm_nb_ASVs_crypt)
Anova(glmm_nb_ASVs_crypt, type = "II")
r2_nakagawa(glmm_nb_ASVs_crypt)
plot(simulateResiduals(glmm_nb_ASVs_crypt))
glmm_nb_ASVs_crypt

glmm_nb_ASVs_cyano <- glmmTMB(count ~ log_conc + (1 | Mesocosm) + (1 | sampling.week),
                         data = dt_Algae_all_matching_count_cyano,
                         family = nbinom2)
performance::check_singularity(glmm_nb_ASVs_all_algae) 
summary(glmm_nb_ASVs_cyano)
Anova(glmm_nb_ASVs_cyano, type = "II")
r2_nakagawa(glmm_nb_ASVs_cyano)
plot(simulateResiduals(glmm_nb_ASVs_cyano))
glmm_nb_ASVs_cyano

glmm_nb_ASVs_green <- glmmTMB(count ~ log_conc + (1 | Mesocosm) + (1 | sampling.week),
                         data = dt_Algae_all_matching_count_green,
                         family = nbinom2)
performance::check_singularity(glmm_nb_ASVs_green) 
summary(glmm_nb_ASVs_green)
Anova(glmm_nb_ASVs_green, type = "II")
r2_nakagawa(glmm_nb_ASVs_green)
plot(simulateResiduals(glmm_nb_ASVs_green))
glmm_nb_ASVs_green

#  relative abundance normalised reads - check beta regression model for proportional responses - relative abundance of ASVs and read numbers
# general trend across algae groups
# We modeled the relative abundance of reads using a beta regression with a logit link. 
# Due to near-zero variance estimates for both Mesocosm and sampling.week, 
# we fitted a simplified model with ALA as a fixed effect and only Algae group as a random. 
# The model explained 10.7% of the variance due to fixed effects (marginal R² = 0.107) 
# and 43.8% when including the random effect of Algae (conditional R² = 0.438)
hist(dt_Algae_all_relabund$Reads)
glmm_relab_reads_all <- glmmTMB(Reads ~ ALA + (1 | Algae) + (1 | Mesocosm) + (1 | sampling.week),
                                data = dt_Algae_all_relabund,
                                family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_reads_all)#+ (1 | Mesocosm) + (1 | sampling.week) cause singularity
# is singular
glmm_relab_reads_all <- glmmTMB(Reads ~ ALA + (1 | Algae),
                            data = dt_Algae_all_relabund,
                            family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_reads_all)
summary(glmm_relab_reads_all)
Anova(glmm_relab_reads_all, type = "II")
r2_nakagawa(glmm_relab_reads_all)
plot(simulateResiduals(glmm_relab_reads_all))
glmm_relab_reads_all

hist(dt_Algae_all_relabund_crypt$Reads)
glmm_relab_reads_crypt <- glmmTMB(Reads ~ ALA + (1 | Mesocosm) + (1 | sampling.week),
                              data = dt_Algae_all_relabund_crypt,
                              family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_reads_crypt)
summary(glmm_relab_reads_crypt)
Anova(glmm_relab_reads_crypt, type = "II")
r2_nakagawa(glmm_relab_reads_crypt)
plot(simulateResiduals(glmm_relab_reads_crypt))
glmm_relab_reads_crypt

glmm_relab_reads_cyano <- glmmTMB(Reads ~ ALA + (1 | Mesocosm) + (1 | sampling.week),
                              data = dt_Algae_all_relabund_cyano,
                              family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_reads_cyano)
summary(glmm_relab_reads_cyano)
Anova(glmm_relab_reads_cyano, type = "II")
r2_nakagawa(glmm_relab_reads_cyano)
plot(simulateResiduals(glmm_relab_reads_cyano))
glmm_relab_reads_cyano

glmm_relab_reads_chrom <- glmmTMB(Reads ~ ALA + (1 | Mesocosm) + (1 | sampling.week),
                              data = dt_Algae_all_relabund_chrom,
                              family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_reads_chrom)
summary(glmm_relab_reads_chrom)
Anova(glmm_relab_reads_chrom, type = "II")
r2_nakagawa(glmm_relab_reads_chrom)
plot(simulateResiduals(glmm_relab_reads_chrom))
glmm_relab_reads_chrom

glmm_relab_reads_green <- glmmTMB(Reads ~ ALA + (1 | Mesocosm) + (1 | sampling.week),
                              data = dt_Algae_all_relabund_green,
                              family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_reads_green)
summary(glmm_relab_reads_green)
Anova(glmm_relab_reads_green, type = "II")
r2_nakagawa(glmm_relab_reads_green)
plot(simulateResiduals(glmm_relab_reads_green))
glmm_relab_reads_green

# beta reg rel abundance ASV counts
glmm_relab_ASVs_all <- glmmTMB(ASVs ~ ALA + (1 | Algae),
                            data = dt_Algae_all_relabund,
                            family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_ASVs_all)#+ (1 | Mesocosm) + (1 | sampling.week) cause singularity
summary(glmm_relab_ASVs_all)
Anova(glmm_relab_ASVs_all, type = "II")
r2_nakagawa(glmm_relab_ASVs_all)
plot(simulateResiduals(glmm_relab_ASVs_all))
glmm_relab_ASVs_all

glmm_relab_ASVs_crypt <- glmmTMB(ASVs ~ ALA  + (1 | Mesocosm) + (1 | sampling.week),
                              data = dt_Algae_all_relabund_crypt,
                              family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_ASVs_crypt)
summary(glmm_relab_ASVs_crypt)
Anova(glmm_relab_ASVs_crypt, type = "II")
r2_nakagawa(glmm_relab_ASVs_crypt)
plot(simulateResiduals(glmm_relab_ASVs_crypt))
glmm_relab_ASVs_crypt

glmm_relab_ASVs_cyano <- glmmTMB(ASVs ~ ALA  + (1 | Mesocosm) + (1 | sampling.week),
                              data = dt_Algae_all_relabund_cyano,
                              family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_ASVs_cyano)
summary(glmm_relab_ASVs_cyano)
Anova(glmm_relab_ASVs_cyano, type = "II")
r2_nakagawa(glmm_relab_ASVs_cyano)
plot(simulateResiduals(glmm_relab_ASVs_cyano))
glmm_relab_ASVs_cyano

glmm_relab_ASVs_chrom <- glmmTMB(ASVs ~ ALA  + (1 | Mesocosm) + (1 | sampling.week),
                              data = dt_Algae_all_relabund_chrom,
                              family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_ASVs_chrom)
summary(glmm_relab_ASVs_chrom)
Anova(glmm_relab_ASVs_chrom, type = "II")
r2_nakagawa(glmm_relab_ASVs_chrom)
plot(simulateResiduals(glmm_relab_ASVs_chrom))
glmm_relab_ASVs_chrom

glmm_relab_ASVs_green <- glmmTMB(ASVs ~ ALA  + (1 | Mesocosm) + (1 | sampling.week),
                              data = dt_Algae_all_relabund_green,
                              family = beta_family(link = "logit"))
performance::check_singularity(glmm_relab_ASVs_green)
summary(glmm_relab_ASVs_green)
Anova(glmm_relab_ASVs_green, type = "II")
r2_nakagawa(glmm_relab_ASVs_green)
plot(simulateResiduals(glmm_relab_ASVs_green))
glmm_relab_ASVs_green


#### Plots

# correlation plots with predicted regression lines from glmmTMB models added

# create plot function
plot_model <- function(model, predictor, response_var,
                       xlab = NULL, ylab = NULL,
                       title = NULL, ylims = NULL,
                       raw_data = NULL, log_y = FALSE, log_x = FALSE,
                       offset_y = 0.5) {
  
  # Get raw data
  if (is.null(raw_data)) raw_data <- model.frame(model)
  if (is.null(xlab)) xlab <- predictor
  if (is.null(ylab)) ylab <- response_var
  
  # Get predictions
  # Define custom prediction range to match data
  x_min <- min(raw_data[[predictor]], na.rm = TRUE)
  x_max <- max(raw_data[[predictor]], na.rm = TRUE)
  pred <- ggpredict(model, terms = paste0(predictor, " [", x_min, ":", x_max,  " by=0.005]"))
  pred_df <- as.data.frame(pred)

  # Plot
  p <- ggplot() +
    geom_jitter(data = raw_data,
                aes(x = .data[[predictor]], y = .data[[response_var]]),
                alpha = 0.4, color = "black",, width = 0.01, height = 0.05) +
    geom_line(data = pred_df, aes(x = x, y = predicted),
              color = "black", linewidth = 1) +
    geom_ribbon(data = pred_df, aes(x = x, ymin = conf.low, ymax = conf.high),
                fill = "black", alpha = 0.2) +
    theme_bw() +
    xlab(xlab) +
    ylab(ylab) +
    ggtitle(title) +
    theme(plot.title = element_text(size = 10))
  
  # Y-limits
  if (!is.null(ylims)) {
    p <- p + coord_cartesian(ylim = ylims)
  }
  
  return(p)
}



#Normalised count plots
# Create plots
library(patchwork)
plot(ggeffect(glmm_nb_all, terms = c("log_conc")))

pl_all <- plot_model(glmm_nb_all, "log_conc", "Norm_count_round",
                     title = "All groups", ylab = "Normalised read count", xlab = "log Biomass (µg/L)", 
                     ylims = c(1, 800)) + scale_y_log10()
pl_chrom <- plot_model(glmm_nb_chrom, "log_conc", "Norm_count_round",
                       title = "Chromophytes", ylab = "Normalised read count", xlab = "log Biomass (µg/L)", 
                       ylims = c(1, 800)) + scale_y_log10()
pl_crypt <- plot_model(glmm_nb_crypt, "log_conc", "Norm_count_round",
                       title = "Cryptophytes", ylab = "Normalised read count", xlab = "log Biomass (µg/L)", 
                       ylims = c(1, 800)) + scale_y_log10()
pl_cyano <- plot_model(glmm_nb_cyano, "log_conc", "Norm_count_round",
                       title = "Cyanobacteria", ylab = "Normalised read count", xlab = "log Biomass (µg/L)", 
                       ylims = c(1, 800)) + scale_y_log10()
pl_green <- plot_model(glmm_nb_green, "log_conc", "Norm_count_round",
                       title = "Green algae", ylab = "Normalised read count", xlab = "log Biomass (µg/L)", 
                       ylims = c(1, 800)) + scale_y_log10()
Fig_NormReads <- pl_chrom + pl_crypt + pl_cyano + pl_green  +
  plot_layout(nrow = 1) +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect") & theme(legend.position = "bottom")

Fig_NormReads
ggsave("Fig_NormReads.pdf", Fig_NormReads, width = 10, height = 3.2)



pl_ASVs_all <- plot_model(glmm_nb_ASVs_all, "log_conc", "count",
                     title = "All groups", ylab = "Number of ASVs", xlab = "log Biomass (µg/L)", 
                     ylims = c(1, 250)) + scale_y_log10()
pl_ASVs_chrom <- plot_model(glmm_nb_ASVs_chrom, "log_conc", "count",
                       title = "Chromophytes", ylab = "Number of ASVs", xlab = "log Biomass (µg/L)",
                       ylims = c(1, 250)) + scale_y_log10()
pl_ASVs_crypt <- plot_model(glmm_nb_ASVs_crypt, "log_conc", "count",
                       title = "Cryptophytes", ylab = "Number of ASVs", xlab = "log Biomass (µg/L)",
                       ylims = c(1, 250)) + scale_y_log10()
pl_ASVs_cyano <- plot_model(glmm_nb_ASVs_cyano, "log_conc", "count",
                       title = "Cyanobacteria", ylab = "Number of ASVs", xlab = "log Biomass (µg/L)",
                       ylims = c(1, 250)) + scale_y_log10()
pl_ASVs_green <- plot_model(glmm_nb_ASVs_green, "log_conc", "count",
                       title = "Green algae", ylab = "Number of ASVs", xlab = "log Biomass (µg/L)",
                       ylims = c(1, 250)) + scale_y_log10()
Fig_ASVs <- pl_ASVs_chrom + pl_ASVs_crypt + pl_ASVs_cyano + pl_ASVs_green +
  plot_layout(nrow = 1) +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect") & theme(legend.position = "bottom")

Fig_ASVs
ggsave("Fig_ASVs.pdf", Fig_ASVs, width = 10, height = 3.2)


pl_relab_ASVs_all <- plot_model(glmm_relab_ASVs_all, "ALA", "ASVs",
                          title = "All groups", ylab = "Relative abundance of ASVs", xlab = "Relative abundance of total Chlorophyll a",
                          ylims=c(0,1))
pl_relab_ASVs_chrom <- plot_model(glmm_relab_ASVs_chrom, "ALA", "ASVs",
                            title = "Chromophytes", ylab = "Relative abundance of ASVs", xlab = "Relative abundance of total Chlorophyll a",
                            ylims=c(0,1))
pl_relab_ASVs_crypt <- plot_model(glmm_relab_ASVs_crypt, "ALA", "ASVs",
                            title = "Cryptophytes", ylab = "Relative abundance of ASVs", xlab = "Relative abundance of total Chlorophyll a",
                            ylims=c(0,1))
pl_relab_ASVs_cyano <- plot_model(glmm_relab_ASVs_cyano, "ALA", "ASVs",
                            title = "Cyanobacteria", ylab = "Relative abundance of ASVs", xlab = "Relative abundance of total Chlorophyll a",
                            ylims=c(0,1))
pl_relab_ASVs_green <- plot_model(glmm_relab_ASVs_green, "ALA", "ASVs",
                            title = "Green algae", ylab = "Relative abundance of ASVs", xlab = "Relative abundance of total Chlorophyll a",
                            ylims=c(0,1))
  
Fig_relab_ASVs <- pl_relab_ASVs_chrom + pl_relab_ASVs_crypt + pl_relab_ASVs_cyano + pl_relab_ASVs_green +
  plot_layout(nrow = 1) +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect") & theme(legend.position = "bottom")
Fig_relab_ASVs
ggsave("Fig_relab_ASVs.pdf", Fig_relab_ASVs, width = 10, height = 3.2)

pl_relab_Reads_all <- plot_model(glmm_relab_reads_all, "ALA", "Reads",
                                title = "All groups", ylab = "Relative abundance of reads", xlab = "Relative abundance of total Chlorophyll a",
                                ylims=c(0,1))
pl_relab_Reads_chrom <- plot_model(glmm_relab_reads_chrom, "ALA", "Reads",
                                  title = "Chromophytes", ylab = "Relative abundance of reads", xlab = "Relative abundance of total Chlorophyll a",
                                  ylims=c(0,1))
pl_relab_Reads_crypt <- plot_model(glmm_relab_reads_crypt, "ALA", "Reads",
                                  title = "Cryptophytes", ylab = "Relative abundance of reads", xlab = "Relative abundance of total Chlorophyll a",
                                  ylims=c(0,1))
pl_relab_Reads_cyano <- plot_model(glmm_relab_reads_cyano, "ALA", "Reads",
                                  title = "Cyanobacteria", ylab = "Relative abundance of reads", xlab = "Relative abundance of total Chlorophyll a",
                                  ylims=c(0,1))
pl_relab_Reads_green <- plot_model(glmm_relab_reads_green, "ALA", "Reads",
                                  title = "Green algae", ylab = "Relative abundance of reads", xlab = "Relative abundance of total Chlorophyll a",
                                  ylims=c(0,1))

Fig_relab_Reads <- pl_relab_Reads_chrom + pl_relab_Reads_crypt + pl_relab_Reads_cyano + pl_relab_Reads_green +
  plot_layout(nrow = 1) +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect") & theme(legend.position = "bottom")
Fig_relab_Reads
ggsave("Fig_relab_Reads.pdf", Fig_relab_Reads, width = 10, height = 3.2)

Fig_Totals <- pl_all + pl_ASVs_all + pl_relab_ASVs_all + pl_relab_Reads_all +
  plot_layout(nrow = 1) +
  plot_layout(guides = "collect", axis_titles = "collect", axes = "collect") & theme(legend.position = "bottom")
Fig_Totals
ggsave("Fig_Totals.pdf", Fig_Totals, width = 10, height = 3.2)




ggplot(dt_Algae_all) +
  geom_bar(aes(x=factor(sampling.occasion), y=concentration, fill=Algae),position = "stack", stat="identity")+
  theme(legend.position = "right") +  theme_bw() +
  scale_fill_brewer(palette = "Paired") + xlab("Sample") +
  facet_wrap(~Mesocosm , nrow=4, ncol=6)
ggplot(dt_Algae_all) +
  geom_bar(aes(x=factor(sampling.occasion), y=Norm_count, fill=Algae),position = "stack", stat="identity")+
  theme(legend.position = "right") +  theme_bw() +
  scale_fill_brewer(palette = "Paired") + ylab("Normalised reads") + xlab("Sample")+
  facet_wrap(~Mesocosm , nrow=4, ncol=6)
ggplot(dt_Algae_all) +
  geom_bar(aes(x=factor(sampling.occasion), y=count, fill=Algae),position = "stack", stat="identity")+
  theme(legend.position = "right") +  theme_bw() +
  scale_fill_brewer(palette = "Paired") + ylab("ASVs") + xlab("Sample")+
  facet_wrap(~Mesocosm , nrow=4, ncol=6)

ggplot(dt_Algae_all) +
  geom_bar(aes(x=factor(sampling.occasion), y=concentration, fill=Algae),position = "fill", stat="identity")+
  theme(legend.position = "right") +  theme_bw() +scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Paired") +  ylab("Relative Abundance (% of total chl a)") + xlab("Sample")+
  facet_wrap(~Mesocosm , nrow=4, ncol=6)
ggplot(dt_Algae_all) +
  geom_bar(aes(x=factor(sampling.occasion), y=Norm_count, fill=Algae),position = "fill", stat="identity")+
  theme(legend.position = "right") +  theme_bw() +scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Paired") + ylab("Relative Abundance (% of total reads)") + xlab("Sample")+
  facet_wrap(~Mesocosm , nrow=4, ncol=6)
ggplot(dt_Algae_all) +
  geom_bar(aes(x=factor(sampling.occasion), y=count, fill=Algae),position = "fill", stat="identity")+
  theme(legend.position = "right") +  theme_bw() +scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Paired") + ylab("Relative Abundance (% of total ASVs)") + xlab("Sample")+
  facet_wrap(~Mesocosm , nrow=4, ncol=6)

# plots showing sums for mesocosms - for ASV counts use dt_Algae_16S_18S_count_mesocosm which contains the number of ASVs per mesocosm, for read numbers and concentration sums are calculated automatically
p_count <- ggplot(dt_Algae_all) +
  geom_bar(aes(x=factor(Mesocosm), y=count, fill=Algae, color=Algae),position = "fill", stat="identity")+
  theme_bw()  +scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Paired") + scale_color_brewer(palette = "Paired") +  ylab("Number of ASVs") + xlab("Mesocosm")+
  theme(axis.text.y   = element_text(size=12), axis.text.x   = element_text(size=12),
        axis.title.y  = element_text(size=12), axis.title.x  = element_text(size=12),
        panel.background = element_blank(), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),
        legend.position="", axis.line = element_line(colour = "black"))

p_ala <- ggplot(dt_Algae_all) +
  geom_bar(aes(x=factor(Mesocosm), y=concentration, fill=Algae, color=Algae),position = "fill", stat="identity")+
  theme_bw() +scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Paired") + scale_color_brewer(palette = "Paired") +  ylab("Relative Abundance (% of total chl a)") + xlab("")+
  theme(axis.text.y   = element_text(size=12), axis.text.x   = element_blank(),
        axis.title.y  = element_text(size=12), axis.title.x  = element_blank(),
        panel.background = element_blank(), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),
        legend.position="right", axis.line = element_line(colour = "black"))

p_reads <- ggplot(dt_Algae_all) +
  geom_bar(aes(x=factor(Mesocosm), y=Norm_count, fill=Algae, color=Algae),position = "fill", stat="identity") +  theme_bw() +scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Paired") + scale_color_brewer(palette = "Paired") + ylab("Relative Abundance (% of total reads)") + xlab("")+
  theme(axis.text.y   = element_text(size=12), axis.text.x   = element_blank(),
        axis.title.y  = element_text(size=12), axis.title.x  = element_blank(),
        panel.background = element_blank(), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),
        legend.position="none", axis.line = element_line(colour = "black"))

Fig_relativeabund<- p_ala+p_reads+p_count+plot_layout(nrow=3)+ plot_annotation(tag_levels = 'A')& 
  theme(plot.tag = element_text(size = 14, face="bold"), plot.tag.position="topleft")
Fig_relativeabund


DT_full_long <- melt(dt_Algae_all_relabund,
                     id.vars = c("Sample.date","sampling.week","treatment","Mesocosm","sampling.occasion","rep","Heated","Invasion","Average_temperature","Algae"),
                     measure.vars = c("ALA","ASVs","Reads"),
                     variable.name = "method",
                     value.name = c("value"))
DT_full_long


DT_full_long=DT_full_long[sampling.week!="0"]

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- plyr::ddply(data, groupvars, .drop=.drop,
                       .fun = function(xx, col) {
                         c(N    = length2(xx[[col]], na.rm=na.rm),
                           mean = mean   (xx[[col]], na.rm=na.rm),
                           sd   = sd     (xx[[col]], na.rm=na.rm)
                         )
                       },
                       measurevar
  )
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

#calculate number of ASVs for each algal group across the mesocosms
sum_dt_Algae_16S_18S_count_group <- as.data.table(summarySE(dt_Algae_16S_18S_count_mesocosm, measurevar="count", 
                                                            groupvars=c("Algae")))
sum_dt_Algae_16S_18S_count_group

p_all <- ggplot(DT_full_long) +
  geom_bar(aes(x=factor(method), y=value, fill=Algae, color=Algae),position = "fill", stat="identity")+
  theme(legend.position = "right") +  theme_bw() +scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Paired")  + scale_color_brewer(palette = "Paired") +  ylab("Relative Abundance") + xlab("")+
  theme(axis.text.y   = element_text(size=12), axis.text.x   = element_text(size=12,angle = 90, vjust = 1, hjust=1),
        axis.title.y  = element_text(size=12), axis.title.x  = element_text(size=12),
        panel.background = element_blank(), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),
        legend.position="right", axis.line = element_line(colour = "black"),panel.border = element_blank())+
  facet_grid(~factor(sampling.occasion))
p_all
ggsave("p_all.jpeg", width = 8, height = 5)
ggsave("p_all.pdf", width = 8, height = 5)





#Venn Diagrams comparing detection of groups in Samples via eDNA and ALA
install.packages("VennDiagram")   # Install & load VennDiagram package
library("VennDiagram")

#Chromo
grid.newpage()                    
draw.pairwise.venn(area1 = 100,   
                   area2 = 91,
                   cross.area = 73)

#Crypto
grid.newpage()                    
draw.pairwise.venn(area1 = 22,   
                   area2 = 134,
                   cross.area = 22)

#Cyano
grid.newpage()                    
draw.pairwise.venn(area1 = 132,   
                   area2 = 120,
                   cross.area = 111)

#Green
grid.newpage()                    
draw.pairwise.venn(area1 = 143,   
                   area2 = 143,
                   cross.area = 142)