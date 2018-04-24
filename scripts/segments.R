library(data.table)
library(splitstackshape)
library(diverdt)
library(ggplot)

seg_l <- 
  unlist(
    lapply(
      readLines('../data/segments/EUR/chro1.segs'),
      function(line_) lapply(strsplit(line_, " "), as.integer)
      )
    )

# Using plus one since R vectors start in one
# Ask Florencia why her index array starts in 0
seg_dt <- data.table(CHR = 1, WD_ID = seg_l + 1) 
seg_dt[, NREP := shift(WD_ID, type='lead') -  WD_ID]
wd_dt <- expandRows(seg_dt[!is.na(NREP)], 'NREP')

eur_dt <- load_bim_frq('../data/reich_data/EUR/EUR', 
                       one_allele_perline = FALSE)

recomb_data <- 
  fread('../data/genetic_map_b36/genetic_map_chr1_b36.txt')
names(recomb_data) <- c('POS', 'RATE', 'CM_R')

eur_dt1 <- eur_dt[CHR == 1]

eur_dt1[, `:=`(POS_ID = .I, 
               WD_ID = wd_dt$WD_ID,
               HTZ = MAF * ( 1 - MAF))]
eur_dt1[, HTZ_WD := mean(HTZ), by = WD_ID]

reur1 <- recomb_data[eur_dt1, on = .(POS)
                     ][!is.na(RATE)]

reur1[WD_ID == POS_ID]
