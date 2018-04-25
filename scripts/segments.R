library(data.table)
library(splitstackshape)
library(diverdt)
library(ggplot)

files_path <- 
  '../data/Results-01-phased-filtered'

files_list <- 
  grep('1e', 
       list.files(files_path, full.names=TRUE), 
       value = TRUE)
files_list

fp_re <- 
  paste0("^(", 
         files_path, 
         ")(/EUR1-1e-)", 
         "([0-9]{2})$")

f_sufix <- 
  gsub(fp_re, "\\1\\2", files_list[[1]])

exponents <- 
  gsub(fp_re, "\\3", files_list)

seg_dt <-
  rbindlist(
    lapply(exponents, 
           function(e) fread(paste0(f_sufix, e))[,EXP := e]
           ))

setnames(seg_dt, 
         paste0("V", 1:4), 
         c("IL", "IR", "WD_SIZE", "AVG_HOM"))




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
