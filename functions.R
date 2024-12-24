# prot_tools

###load data
load('data-raw/hydropath_ref.rda')
load("data-raw/aa_dict.rda")
load("data-raw/aa_mw_mono.rda")
load("data-raw/aa_mw_avg.rda")

#### example data
s2 <- list(
  h1_a = "PYVALFEKCCLIGCTKRSLAKYC", 
  h1_b = "VAAKWKDDVIKLCGRELVRAQIAICGMSTWS",
  h2_a = "LYSALANKCCHVGCTKRSLARFC",
  h2_b = "DSWMEEVIKLCGRELVRAQIAICGMSTWS",
  h3_a = "DVLAGLSSSCCKWGCSKSEISS",
  h3_b = "RAAPYGVRLCGREFIRAVIFTCGGSRW"
)
## loading and formatting
### format protein sequence to upper case and single letter code
transform_sequence <- function(aa_seq) {
  s <- unlist(strsplit(toupper(aa_seq), NULL))
  return(s)
}

### read fasta
read_fasta <- function(data){
  d <- return(Biostrings::readAAStringSet(data))
  
  data.frame(
    "Sequence" = d@ranges@seq,
    "UniprotID" = stringr::str_extract(d@ranges@NAMES, "\\|(.+)\\|"),\
    "GeneName" = stringr::str_extract(d@ranges@NAMES, "")
  )
  # results <- setNames(list(stringr::str_extract(data$names, "\\|(.+)\\|"), ), stringr::str_extract(data, "\\|(.+)$"))
}
str(read_fasta("../relaxin/relaxin_superfamily_uniprotkb_2024_11_26.fasta"))
read_fasta("../relaxin/relaxin_superfamily_uniprotkb_2024_11_26.fasta")@ranges@NAMES

read_fasta2 <- function(data){
  data.frame(
    "GeneName" = NULL,
    "UniprotID" = NULL,
    "Sequence" = NULL
    )
  
  current_entry <- NULL
  
  
  results <- rbind(results)
}


## helper functions
### hydropathy
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

### identify modifications
unlist(stringr::str_extract_all(s, "([A-Z]\\([a-z]+\\))"))
# extract mod
# check with mod_list
# extract mod_mass and add to aa_mass

### mass calculations
# accepts U/Sel for selenocysteine and O/Pyl for pyrrolysine
# Three-letter code is accepted but generally the output will be in one-letter code.

calculate_mass_avg <- function(aa_seq){
  mass <- 0
  for (a in unlist(strsplit(aa_seq, ""))) {
    if (a %in% names(aa_mw_avg)){
      mass <- mass + unname(aa_mw_avg[a])
    }
  }
  return(mass)
}

calculate_mass_mono <- function(aa_seq){
  mass <- 0
  for (a in unlist(strsplit(aa_seq, ""))) {
    if (a %in% names(aa_mw_mono)){
      mass <- mass + unname(aa_mw_mono[a])
    }
  }
  return(mass)
}

# calculate n 
calculate_n <- function(aa_seq, mass = NA, avg = TRUE){
  if (!is.numeric(m)){
    stop("Mass must be numeric.")
  }
  
  if (!is.character(aa_seq) | !is.list(aa_seq)) {
    stop("Aa_seq must be a character vector or list")
  } else if (all(sapply(aa_seq, is.character))) {
    stop("All entries must be character strings.")
  }
  
  if (avg == TRUE) {
    pep_mw <- calculate_mass_avg(aa_seq)
  } else {
    pep_mw <- calcualte_mass_mono(aa_seq)
  }
  n <- m / pep_mw
  return(n)
}


#splits protein sequences into peptides by index 
split_sequence_by_index <- function(aa_seq, index = NA){
  
  if (!is.character(index)) {
    stop("Indices must be character, use '10-35' instead of 10-35.")
  }
  
  s <- transform_sequence(aa_seq)
  
  peptides <- NULL
  current_pep <- NULL
  for (i in index) {
    start_i <- as.numeric(strsplit(i, "-")[[1]][1])
    print(start_i)
    stop_i <- as.numeric(strsplit(i, "-")[[1]][2])
    print(stop_i)
    # throw error if start>stop
    if (start_i > stop_i) {
      stop("Start index is smaller than stop index.")
    }
    if (start_i < 1 | start_i > length(s) | stop_i > length(s)) {
      stop("At least index is outside of protein sequence.")
    }
    print(s[start_i:stop_i])
    current_pep <- paste0(s[start_i:stop_i], collapse = "")
    peptides <- append(peptides, current_pep)
  }
  return(peptides)
}

calculate_aa_mass <- function(aa_seq, aa, avg = TRUE, prot_m = 1) {
  
  if (!is.character(aa_seq) || !is.character(aa)) {
    stop("Both aa_seq and aa should be character strings.")
  }
  
  s <- transform_sequence(aa_seq)
  aa <- transform_sequence(aa)
  
  if (!all(unique(aa) %in% names(aa_mw_avg))) {
    stop("Non-canonical amino acid detected. Please use single letter codes as found in aa_dict.")
  }
  
  if (!all(unique(aa) %in% unique(s))) {
    warning("One or more specified amino acids not found in sequence.")
  }
  
  st <- table(s)
  #  if "use_all" use all amino acids in protein sequence
  if (any(aa == "use_all")) {
    aoi <- st
  } else {
    # extract only specified amino acids
    aoi <- st[aa]
  }
  # calculate total protein mass
  protein_aa_m <- setNames(rep(0, length(aoi)), names(aoi))
  for (a in names(st)) {
    protein_aa_m[a] <- st[[a]] * aa_mw_avg[[a]]
  }
  protein_m <- sum(protein_aa_m)
  
  # calculate aa masses 
  aa_percent_m <- setNames(rep(0, length(aoi)), names(aoi))
  aa_m_total <- setNames(rep(0, length(aoi)), names(aoi))
  
  for (a in names(aoi)) {
   aa_percent_m[a] <- ((aoi[[a]] * (aa_mw_avg[[a]] + 18)) / protein_m) * 100
   aa_m_total[a] <- prot_m * (aa_percent_m[[a]] / 100)
  }
  
  return(aa_m_total)
}

calculate_aa_mass("acafg", c("a", "c"), prot_m = 1)

## protein summary
protein_summary <- function(aa_seq){
  s <- transform_sequence(aa_seq)
  r <- table(s)/length(s) * 100
  
  return(data.frame(
    AminoAcid = names(s),
    OccurrencesTotal = as.vector(s),
    OccurrencesRelative = as.vector(r)
    ))
  return(data.frame(
    AcidicAATotal = names()
  ))
# #acidic and basic aa
# theoretical #peptides and avg pep length?
  
}



# digest proteins
digest_protein <- function(aa_seq, 
                           missed = 0, 
                           specificity = c("P", "A"), 
                           combined = TRUE,
                           min_length = 1, 
                           tabulate = TRUE) {
  # allow multi-protease digestion
  
  # loop specificities
  
  # if (combined) {
  #   digest all peptides sequentially
  # }
  
  
  spec_start <- NULL
  spec_stop <- NULL
  
  results <- NULL
  
  for (i in seq_along(aa_seq)){
    
    aa_seq <- stringr::str_extract(aa_seq[i][[1]], "(?<=\\.)[A-Z]+(?=\\.|$)")
  
    peptides <- list()
    current_pep <- ""
    aa_index <- ""
  
    for (i in unlist(strsplit(aa_seq, ""))) {
      aa <- substring(aa_seq, i, i)
      current_pep <- paste0(current_pep, aa)
      
      if (aa %in% specificity){
        peptides <- append(peptides, current_pep)
        if(missed == 0){
          current_pep <- ""
        }
      }
      # if max(missed): add to peptides and reset current_pep
      if (sum(stringr::str_count(current_pep, specificity) == missed + 1)){
        peptides <- append(peptides, current_pep)
        current_pep <- ""
      }
    }  
      
    if (nchar(current_pep) > 0){
      peptides <- append(peptides, current_pep)
    }
    
    peptides <- Filter(function(i) nchar(i) >= min_length, peptides)
    
    prot_pep <- list(peptides = peptides)
    results[[names(aa_seqs[i])]] <- prot_pep
      
      # count basic aas in peptide
      # basic_aa <- c("K", "R", "H")
      # acidic_aa <- c("D", "E")
      # n_base <- NA
      # n_acid <- NA
    
      # for (p in peptides) {
      #   for (aa in p)
      #   n_base <- sum(stringr::str_count(aa, basic_aa))
      #   n_acid <- sum(stringr::str_count(aa, acidic_aa))
      # }
      
      # return(results)
    
    if (tabulate == TRUE) {
      results <- data.frame(
        "protein" = NA,
        "peptide" = NA,
        "missed_cleavages" = NA,
        "basic_aas" = n_base,
        "acidic_aas" = n_acid
      )
    }
  }
  return(results)
}

digest_protein(s)

basic_aa <- c("K", "R", "H")
acidic_aa <- c("D", "E")
n_base <- NA
n_acid <- NA
for (p in peptides) {
  for (aa in p)
    n_base <- sum(stringr::str_count(aa, basic_aa))
    n_acid <- sum(stringr::str_count(aa, acidic_aa))
}
