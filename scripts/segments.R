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


eur_segs_ref <-
  fread(paste0(files_path, '/EURmap'), 
        select = 1:3,
        col.names = c("CHR", "POS", "SNP"))

eur_dt <- 
  load_bim_frq('../data/reich_data/EUR/EUR', 
               one_allele_perline = FALSE)

eur_dt1 <- 
  eur_dt[eur_segs_ref, on = .(CHR, POS, SNP)
         ][,`:=`(POSD_ID = .I, HTZ = MAF * ( 1 - MAF))]

recomb_data <- 
  fread('../data/genetic_map_b36/genetic_map_chr1_b36.txt',
        col.names = c('POS', 'RATE', 'CM_R'))

wd_dt <- expandRows(seg_dt, 'WD_SIZE')

eur_wd <- 
  eur_dt1[rep(1:.N, length(exponents)), 
          `:=`(EXP = wd_dt$EXP, 
           WD_ID = wd_dt$IR)]

eur_wd[, HTZ_WD := mean(HTZ), by = .(WD_ID, EXP)]

reur1 <- recomb_data[eur_dt1, on = .(POS)
                     ][!is.na(RATE)]

reur1[WD_ID == POS_ID]
