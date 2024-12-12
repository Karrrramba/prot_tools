# prot_tools
library()
## hydropathy
# hyd <- "Ala:  1.800  Arg: -4.500  Asn: -3.500  Asp: -3.500  Cys:  2.500  Gln: -3.500  Glu: -3.500  Gly: -0.400  His: -3.200  Ile:  4.500  Leu:  3.800  Lys: -3.900  Met:  1.900  Phe:  2.800  Pro: -1.600  Ser: -0.800  Thr: -0.700  Trp: -0.900  Tyr: -1.300  Val:  4.200"
# hyd_vals <- as.numeric(unlist(stringr::str_extract_all(hyd, "[\\d+|\\.|\\-]+")))
# hyd_names <- unlist(stringr::str_extract_all(hyd, "([a-zA-Z]+)"))
# 
# hydropath_ref <- setNames(hyd_vals, aa_code_1[1:20])
# hydropath_ref
# save(hydropath_ref, file = "data/hdrpth_ref.rda")

s <- "FYPYPC"
load('data/hdrpth_ref.rda')

calculate_hydropathy_index <- function(aa_seq) {
  s <- unlist(strsplit(aa_seq, split = ""))
  h <- 0
  for (i in s) {
     if (i %in% names(hydropath_ref)) {
       h <- h + unname(hydropath_ref[i])
     }
  }
  return(h)
}

# identify modifications
unlist(stringr::str_extract_all(s, "([A-Z]\\([a-z]+\\))"))
# extract mod
# check with mod_list
# extract mod_mass and add to aa_mass

## mass calculation
# accepts U/Sel for selenocysteine and O/Pyl for pyrrolysine
# Three-letter code is accepted but generally the output will be in one-letter code.
# aa_masses <- "
# monoisotopic	average
# Alanine (A)	71.03711	71.0788
# Arginine (R)	156.10111	156.1875
# Asparagine (N)	114.04293	114.1038
# Aspartic acid (D)	115.02694	115.0886
# Cysteine (C)	103.00919	103.1388
# Glutamic acid (E)	129.04259	129.1155
# Glutamine (Q)	128.05858	128.1307
# Glycine (G)	57.02146	57.0519
# Histidine (H)	137.05891	137.1411
# Isoleucine (I)	113.08406	113.1594
# Leucine (L)	113.08406	113.1594
# Lysine (K)	128.09496	128.1741
# Methionine (M)	131.04049	131.1926
# Phenylalanine (F)	147.06841	147.1766
# Proline (P)	97.05276	97.1167
# Serine (S)	87.03203	87.0782
# Threonine (T)	101.04768	101.1051
# Tryptophan (W)	186.07931	186.2132
# Tyrosine (Y)	163.06333	163.1760
# Valine (V)	99.06841	99.1326
# 
# Selenocysteine (U)	150.953636	150.0388
# Pyrrolysine (O)	237.147727	237.3018
# "
# 
# aa_masses_mono <- gsub("\t", "", unlist(stringr::str_extract_all(aa_masses, "\\t(\\d+\\.\\d+)\\t")))
# aa_masses_avg <- gsub("\n", "", unlist(stringr::str_extract_all(aa_masses, "(\\d+\\.\\d+)\n")))
# aa_names <- unlist(stringr::str_remove_all(gsub("\n", "", unlist(stringr::str_extract_all(aa_masses, "\n.+\t"))), "\t.+"))[-1]
# 
# aa_codes <- unlist(stringr::str_extract_all(aa_names, "\\b([a-zA-Z]{1,3})"))
# aa_codes <- na.omit(gsub("^[a-z]{3}", NA, aa_codes))
# aa_code_3 <- na.omit(gsub("[A-Z]$", NA, aa_codes))
# aa_code_3 <- gsub("Sel", "Sec", aa_code_3)
# aa_code_3 <- gsub("Pyr", "Pyl", aa_code_3)
# aa_code_1 <- unlist(stringr::str_extract_all(aa_codes, "[A-Z]$"))
# 
# aa_dict <- setNames(aa_code_1, aa_code_3)
# save(aa_dict, file = "data/aa_dict.rda")
# aa_mw_mono <- setNames(aa_masses_mono, aa_code_1)
# save(aa_mw_mono, file = "data/aa_mw_mono.rda")
# aa_mw_avg <- setNames(aa_masses_avg, aa_code_1)
# save(aa_mw_avg, file = "data/aa_mw_avg.rda")

load("data/aa_dict.rda")
load("data/aa_mw_mono.rda")
load("data/aa_mw_avg.rda")

calculate_aa_n <- function(aa_seq, n, m = NA){
  
}
