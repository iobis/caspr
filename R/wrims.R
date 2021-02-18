#' Get the complete WRiMS checklist
#'
#' @export
wrims_checklist <- function() {
    limit <- 1000
    results_list <- list()
    counter <- 0
    while (TRUE) {
      res <- rgbif::name_usage(
        datasetKey = "0a2eaf0c-5504-4f48-a47f-c94229029dc8",
        limit = limit,
        start = limit * counter
      )
      results_list[[counter + 1]] <- res$data
      if (res$meta$endOfRecords) {

        result <- bind_rows(results_list) %>%
          select(taxonID, canonicalName, taxonomicStatus, rank) %>%
          mutate(taxonID = as.numeric(stringr::str_remove(taxonID, "urn:lsid:marinespecies.org:taxname:"))) %>%
          filter(rank == "SPECIES" & taxonomicStatus == "ACCEPTED") %>%
          distinct() %>%
          as_tibble()

        taxa <- taxon(result$taxonID)
        result <- result %>%
          select(taxonID, canonicalName) %>%
          left_join(taxa, by = "taxonID") %>%
          mutate(
            scientificName = ifelse(is.na(scientificName), canonicalName, scientificName),
            taxonRank = "Species",
            wrims = TRUE
          ) %>%
          select(-canonicalName)

        return(result)
      }
      counter <- counter + 1
    }
}
