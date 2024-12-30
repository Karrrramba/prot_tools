# prot_tools

###load data
load('data-raw/hydropath_ref.rda')
load("data-raw/aa_dict.rda")
load("data-raw/aa_mw_mono.rda")
load("data-raw/aa_mw_avg.rda")

#### example data
s1 <- "../relaxin/relaxin_superfamily_uniprotkb_2024_11_26.fasta"

s2 <- list(
  h1_a = "PYVALFEKCCLIGCTKRSLAKYC", 
  h1_b = "VAAKWKDDVIKLCGRELVRAQIAICGMSTWS",
  h2_a = "LYSALANKCCHVGCTKRSLARFC",
  h2_b = "DSWMEEVIKLCGRELVRAQIAICGMSTWS",
  h3_a = "DVLAGLSSSCCKWGCSKSEISS",
  h3_b = "RAAPYGVRLCGREFIRAVIFTCGGSRW"
)

bsa_seq <- "DTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEPERNECFLSHKDDSPDLPKLKPDPNTLCDEFKADEKKFWGKYLYEIARRHPYFYAPELLYYANKYNGVFQECCQAEDKGACLLPKIETMREKVLASSARQRLRCASIQKFGERALKAWSVARLSQKFPKAEFVEVTKLVTDLTKVHKECCHGDLLECADDRADLAKYICDNQDTISSKLKECCDKPLLEKSHCIAEVEKDAIPENLPPLTADFAEDKDVCKNYQEAKDAFLGSFLYEYSRRHPEYAVSVLLRLAKEYEATLEECCAKDDPHACYSTVFDKLKHLVDEPQNLIKQNCDQFEKLGEYGFQNALIVRYTRKVPQVSTPTLVEVSRSLGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVPKAFDEKLFTFHADICTLPDTEKQIKKQTALVELLKHKPKATEEQLKTVMENFVAFVDKCCAADDKEACFAVEGPKLVVSTQTALA"
## loading and formatting
### format protein sequence to upper case and single letter code
transform_sequence <- function(aa_seq) {
  s <- unlist(strsplit(toupper(aa_seq), NULL))
  return(s)
}

### read fasta
read_fasta <- function(file, same_org = TRUE) {
  # Read the file line by line
  fasta <- readLines(file)
  # Identify header lines
  ind <- grep(">", fasta)
  # Identify the sequence lines
  s <- data.frame(ind = ind, from = ind + 1, to = c((ind - 1)[-1], length(fasta)))
  # Process sequence lines
  seqs <- rep(NA, length(ind))
  for(i in 1:length(ind)) {
    seqs[i] <- paste(fasta[s$from[i]:s$to[i]], collapse = "")
  }
  # Create a data frame 
  d <- data.frame(
    # Extract UniprotID
    UniprotID = stringr::str_match(string = fasta[ind],pattern = "\\|(.+?)\\|")[, 2],
    # Extract gene name
    GeneName = stringr::str_match(string = fasta[ind],pattern = "[A-Z0-9]\\|([A-Z0-9_]+)")[, 2],
    Sequence = seqs
    )
  # Remove organism from gene name
  if (same_org) {
    d$GeneName <- gsub("_.*", "", d$GeneName)
  }
  
  # Return the data frame as a result object from the function
  return(d)
}

data <- read_fasta(s1)

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
  aa_seq <- transform_sequence(aa_seq)
  mass <- 0
  for (a in aa_seq) {
    if (a %in% names(aa_mw_avg)){
      mass <- mass + unname(aa_mw_avg[a])
    }
  }
  return(mass)
}

calculate_mass_mono <- function(aa_seq){
  aa_seq <- transform_sequence(aa_seq)
  mass <- 0
  for (a in aa_seq) {
    if (a %in% names(aa_mw_mono)){
      mass <- mass + unname(aa_mw_mono[a])
    }
  }
  return(mass)
}

#splits protein sequences into peptides by index ---- 
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

# calculate molar amount of protein ----
# aa_seq = Amino acid sequence.
# mass = protein mass in g
# 
calculate_prot_n <- function(aa_seq, mass = NA, m_avg = TRUE){
  if (!is.numeric(mass)){
    stop("Mass must be numeric.")
  }
  
  if (!is.character(aa_seq) | !is.list(aa_seq)) {
    stop("Aa_seq must be a character vector or list")
  } else if (all(sapply(aa_seq, is.character))) {
    stop("All entries must be character strings.")
  }
  
  aa_seq <- transform_sequence(aa_seq)
  
  if (m_avg == TRUE) {
    pep_mw <- calculate_mass_avg(aa_seq)
  } else {
    pep_mw <- calcualte_mass_mono(aa_seq)
  }
  n <- mass / pep_mw
  return(n)
}

# aa mass----
calculate_aa_m <- function(aa_seq, aa, avg = TRUE, prot_m = 1) {
  
  if (!is.character(aa_seq) || !is.character(aa)) {
    stop("Both aa_seq and aa should be character strings.")
  }
  
  s <- transform_sequence(aa_seq)
  aa <- transform_sequence(aa)
  
  if (!all(unique(aa) %in% names(aa_mw_avg))) {
    stop("Non-canonical amino acid detected.")
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
  aa_n <- setNames(rep(0, length(aoi)), names(aoi))
  for (a in names(aoi)) {
   aa_percent_m[a] <- ((aoi[[a]] * (aa_mw_avg[[a]] + 18)) / protein_m) * 100
   aa_m_total[a] <- prot_m * (aa_percent_m[[a]] / 100)
   aa_n[a] <- aa_m_total[a] / aa_mw_avg[[a]]
  }
  
  return(aa_m_total)
}
calculate_aa_m(bsa_seq, "y", prot_m = 1E-3)

## protein summary -----
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

# specificty: 
digest_protein <- function(aa_seq, 
                           missed = 0, 
                           specificity = c("PA"),  #given as string of one-letter code divided by a pipe ("|") for exceptions, e.g. "KR|P" for trypsin
                           combined = TRUE,
                           min_length = 1, 
                           min_charge = 1,
                           ommit = NA) {
  
  # allow multi-protease digestion
  aa_seq <- transform_sequence(aa_seq)
  
  # Store created peptides
  peptides <- character()
  
  # Function to perform digestion for a given protease specificity
  digest_for_protease <- function(aa_seq, specificity, missed) {
    
    current_pep <- ""
    temp_peptides <- character()
    
    for (aa in aa_seq) {
      current_pep <- paste0(current_pep, aa)
      
      # If aa is protease-specific cut peptide
      if (aa %in% specificity) {
        temp_peptides <- c(temp_peptides, current_pep)
        if (missed == 0) {
          current_pep <- ""
        }
      }
      
      # If max(missed): add to peptides and reset current_pep
      if (sum(stringr::str_count(current_pep, specificity) == missed + 1)) {
        temp_peptides <- c(temp_peptides, current_pep)
        current_pep <- ""
      }
    }
    
    if (nchar(current_pep) > 0) {
      temp_peptides <- c(temp_peptides, current_pep)
    }
    
    return(temp_peptides)
  }
  
  # Loop through proteases if combined == FALSE, or digest them together if combined == TRUE
  if (combined) {
    # Combine all specificity rules into a single digestion pass
    spec <- unlist(strsplit(toupper(gsub("\\|", "", specificity)), NULL))
    peptides <- digest_for_protease(aa_seq, spec, missed)
  } else {
    # Perform separate digestions for each specificity
    for (spec in specificity) {
      spec <- unlist(strsplit(toupper(gsub("\\|", "", spec)), NULL))
      peptides <- c(peptides, digest_for_protease(aa_seq, spec, missed))
    }
  }
  
  # Filter peptides based on minimum length
  peptides <- peptides[nchar(peptides) >= min_length]
  
  # Filter peptides for specified minimum charge
  peptides <- peptides[sapply(peptides, function(pep) stringr::str_count(pep, "K|H|R") >= min_charge - 1)]
  
  # Filter peptides for unwanted amino acids
  peptides <- peptides[!stringr::str_detect(peptides, paste(ommit, collapse = "|"))]
  
  return(peptides)
}

digest_protein(bsa_seq, specificity = c("KR"), min_length = 3, min_charge = 3)
