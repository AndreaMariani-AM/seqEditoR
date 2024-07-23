sub_CGI <- function(seqs, orientation, seq_to_add, offset, direction=c("upstream", "downstream")) {
  # convert the DNA strings to a char vector
  dna_chars <- as.character(seqs)
  chars_to_add <- as.character(seq_to_add)

  # Offset controls how far from the TSS/center the sequence can be added.
  # If offset < length(chars_to_add) the sequence can cover the TSS and go downstream of it
  middle_point <- nchar(dna_chars)/2

  # Check whether to add upstream and downstream. Consider the direction of the transcription
  if(match.arg(direction) == "upstream"){
    if(orientation == "+"){
      start_point <- middle_point - offset
      substr(dna_chars, start=start_point, stop=start_point+nchar(seq_to_add)) <- seq_to_add
    } else {
      start_point <- middle_point + offset
      substr(dna_chars, start=start_point-nchar(seq_to_add), stop=start_point) <- seq_to_add
    }
  } else {
    if(orientation == "+"){
      start_point <- middle_point + offset
      substr(dna_chars, start=start_point-nchar(seq_to_add), stop=start_point) <- seq_to_add
    } else {
      start_point <- middle_point - offset
      substr(dna_chars, start=start_point, stop=start_point+nchar(seq_to_add)) <- seq_to_add
    }
  }
  #convert to DNAString
  return(Biostrings::DNAString(dna_chars))
}


substitute_CGI <- function(seqs, orientation, seq_to_add, offset, direction=c("upstream", "downstream")) {
  Biostrings::DNAstringSet(
    mapply(sub_CGI, seqs, orientation, as.character(seq_to_add), offset, direction))
}

