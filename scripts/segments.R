library(data.table)
library(splitstackshape)
library(diverdt)
library(ggplot2)

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
         ][,`:=`(POS_ID = .I, HTZ = MAF * ( 1 - MAF))]


wd_dt <- expandRows(seg_dt, 'WD_SIZE')

eur_wd <- 
  eur_dt1[rep(1:.N, length(exponents))
          ][, `:=`(EXP = wd_dt$EXP,
                   IL = wd_dt$IL,
                   IR = wd_dt$IR)]

eur_wd[, HTZ_WD := mean(HTZ), by = .(IR, EXP)]

recomb_data <- 
  fread('../data/genetic_map_b36/genetic_map_chr1_b36.txt',
        col.names = c('POS', 'RATE', 'CM_R'))


reur <- recomb_data[eur_wd, on = .(POS)
                    ][!is.na(RATE)]

reur[,`:=`(WD_AVG_RATE = mean(RATE),
            WD_SD_RATE = sd(RATE)),
      by = .(EXP, IR)]

reur[POS_ID == IL | POS_ID == IR , 
     .(mean(RATE), var(RATE)), 
     by = EXP]

reur[, .(mean(RATE), var(RATE)), by = EXP]






p_reur <- reur[POS_ID %in% 1200:1500]
segs_p <- unique(p_reur[EXP %in% c('01', '05', '10', '15'), 
                 .(IR, IL, WD_AVG_RATE, WD_SD_RATE, EXP)])
ggplot(p_reur[EXP == '01']) + 
  geom_line(aes(x=POS_ID, y = RATE)) + 
  theme_bw() + 
  scale_color_brewer(palette = 'Set1') + 
  geom_segment(data = segs_p,
               aes(x = IL, xend = IR, 
                   y = WD_AVG_RATE, yend = WD_AVG_RATE, 
                   size = WD_SD_RATE,
                   color = EXP 
                   ),
               alpha = 0.7
            )
