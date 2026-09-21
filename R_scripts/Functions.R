##########Editing Rate #####################
for (n in 1:nrow(filename)){
  if (filename$REF[n] == "A"){
    filename$Editing.rate[n] = filename$G[n]/filename$TOTAL[n]
  }
  else{
    filename$Editing.rate[n] = filename$C[n]/filename$TOTAL[n]
  }
} 


###############Calculate Variance####################
########assuming per sample editing rates are split by"_" in a column

Var <- function(x) {
  values <- as.numeric(unlist(strsplit(x,"_")))
  var(values)           
}

#########Extract FASTA from ref genome#####################
#Load ref_genome
ref_genome <- readDNAStringSet("mm10.fa", format = "fasta")

#Define window size (+- 100bp from edited position)
window_size <- 100
###Create Function###################
########chromosome names in mm10 assembly start with "chr" prefix. Modify your chromosome names to include chr
extract_FASTA <- function(chr, position,ref_genome) {  ###looks for three vectors chr, position and ref_fasta seq, 
  start_pos <- position - window_size
  end_pos <- position + window_size
  start_pos <- max(start_pos, 1)####return at least 1 value if length is not equal to start_pos
  end_pos <- min(end_pos, length(ref_genome[[chr]])) ####seq names in ref genome are as chr1, chr2 etc choose the appropriate name for a given fasta file
  seq <- ref_genome[[chr]][start_pos:end_pos] ####extract from ref_genome using "chr" names for length start_pos:end_pos
  return(as.character(seq))
}

########################END########################################





