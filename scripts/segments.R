library(data.table)
library(splitstackshape)
help(package='splitstackshape')
library(diverdt)

seg_l <- 
  unlist(
    lapply(
      readLines('../data/segments/EUR/chro1.segs'),
      function(line_) strsplit(line_, " ")
      )
    )

seg_dt <- data.table(CHR = 1, WD_POS_ID = seg_l)
seg_dt[, WD_ID := as.integer(WD_POS_ID)]
seg_dt[, NREP := shift(WD_ID, type='lead') -  WD_ID]
wd_dt <- expandRows(seg_dt[!is.na(NREP)], 'NREP')
           

eur_dt <- load_bim_frq('../data/reich_data/EUR/EUR', 
                       one_allele_perline = FALSE)

recomb_data <- 
  fread('../data/genetic_map_b36/genetic_map_chr1_combined_b37.txt')
names(recomb_data) <- c('POS', 'RATE', 'CM_R')

recomb_data

eur_dt[CHR == 1, WD_ID := wd_dt$WD_ID]
eur_dt1 <- eur_dt[CHR == 1]

eur_dt1[, HTZ := MAF * ( 1 - MAF)]

reur1 <- recomb_data[eur_dt1, on = .(POS)]
reur1[!is.na(RATE)]

eur_dt1[, mean(HTZ), by = WD_ID]
