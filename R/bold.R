#' Get statistics on marker sequences in BOLD by taxon
#'
#' @param taxon taxon name
#' @return a tibble
#' @export
bold_statistics <- function(taxon) {
  res <- bold::bold_seqspec(taxon)
  if (!is.data.frame(res)) {
    return(tibble())
  }
  df <- res %>%
    bind_rows() %>%
    select(sampleid, institution_storing, sequenceID, markercode, genbank_accession, nucleotides, seq_primers)
  marker_stats <- df %>%
    group_by(markercode) %>%
    summarize(sequences = n(), .groups = "drop") %>%
    filter(!is.na(markercode) & markercode != "")
  if (nrow(marker_stats) == 0) {
    return(tibble())
  }
  return(as_tibble(marker_stats))
}
