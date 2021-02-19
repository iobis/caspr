#' WKT string for the South Pacific
#' @export
wkt_south_pacific <- "POLYGON ((-218 -62, -218 14, -86 14, -67 -23, -74 -62, -218 -62))"

#' Fetch expert priority checklists from the PacMAN GitHub repository.
#'
#' @param summarize group by taxon and concatenate references and remarks.
#' @return The expert priority checklist.
#' @export
expert_checklists <- function(summarize = TRUE) {
  df <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1s45jM7oR4lbURG5JN1r6XamoWjGMZt7hnexAgc8yTYY/edit?usp=sharing") %>%
    select(taxonID = AphiaID_accepted, scientificName = scientificName_accepted, references, remarks = taxonRemarks, taxonRank) %>%
    mutate(taxonID = suppressWarnings(as.numeric(taxonID))) %>%
    filter(!is.na(taxonID))
  records <- bind_rows(worrms::wm_record_(df$taxonID)) %>%
    select(AphiaID, phylum, class, family)
  df <- df %>%
    left_join(records, by = c("taxonID" = "AphiaID"))
  if (summarize == TRUE) {
    df <- df %>%
      group_by(taxonID, scientificName, taxonRank, phylum, class, family) %>%
      summarize(references = paste0(stri_remove_empty_na(unique(references)), collapse = "; "), remarks = paste0(stri_remove_empty_na(unique(remarks)), collapse = "; ")) %>%
      ungroup()
  }
  df
}

#' Fetch a checklist from OBIS.
#'
#' @param wrims limit checklist to WRiMS taxa.
#' @param area limit checklist to area, options are south_pacific and fiji.
#' @return An OBIS checklist.
#' @export
obis_checklist <- function(wrims = TRUE, area = "south_pacific") {
  geometry <- NULL
  areaid <- NULL
  if (area == "south_pacific") {
    geometry <- wkt_south_pacific
  } else if (area == "fiji") {
    areaid <- 68
  }
  robis::checklist(wrims = wrims, geometry = geometry, areaid = areaid)
}
