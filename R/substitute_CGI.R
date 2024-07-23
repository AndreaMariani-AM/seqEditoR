#' Add a user provided CG Island anywhere in the sequence. It operates on a sequence and its corresponding orientation
#'
#' @param seqs A list of Biostrings sequences to perform operations on
#' @param orientation A vector containing the orientation, "+", or "-" of the sequence.
#' @param seq_to_add The string to be added
#' @param offset How many bps from the center of the sequence should the sequence be added
#' @param direction Which direction the substitutions needs to happen. Valid inputs are "upstream" or "downstream".
#' This is affected by the orientation
#'
#' @return A DNAStringSet object
#' @import Biostrings
#' @export
#'
#' @examples
#' random_seq <- Biostrings::DNAStringSet("ATCGTATCGATCGATCGATCATCGATCGATCGATCAG")
#' OR <- c("+")
#' substitute_CGI(seqs=random_seq, orientation=OR, seq_to_add="CGCG", offset=10, direction="upstream")

substitute_CGI <- function(seqs, orientation, seq_to_add, offset, direction=c("upstream", "downstream")) {
  Biostrings::DNAStringSet(
    mapply(sub_CGI, seqs, orientation, as.character(seq_to_add), offset, direction))
}


#' Utilily for substitute_CGI
#'
#' @param seqs A Biostrings sequence to perform operations on
#' @param orientation The orientation of the sequence, "+", or "-".
#' @param seq_to_add The string to be added
#' @param offset How many bps from the center of the sequence should the sequence be substituted.
#' @param direction Which direction the substitutions needs to happen. Valid inputs are "upstream" or "downstream".
#' This is affected by the orientation
#'
#' @return A DNAString object
#' @import Biostrings
#' @export
#'
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
