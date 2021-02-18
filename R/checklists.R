#' WKT string for the South Pacific
#' @export
wkt_south_pacific <- "POLYGON ((-218 -62, -218 14, -86 14, -67 -23, -74 -62, -218 -62))"

#' Fetch expert priority checklists from the PacMAN GitHub repository.
#'
#' @param summarize group by taxon and concatenate references and remarks.
#' @param enhance add complete taxonomy from OBIS.
#' @return The expert priority checklist.
#' @export
expert_checklists <- function(summarize = TRUE, enhance = FALSE) {
  df <- read.csv("https://raw.githubusercontent.com/iobis/pacman/main/SpeciesLists/SpeciesList.tsv?token=AADXUOMJQ24AIINJURWKUYDAFLISW", sep = "\t") %>%
    select(taxonID = AphiaID_accepted, scientificName = scientificName_accepted, references, remarks = taxonRemarks, taxonRank) %>%
    mutate(taxonID = suppressWarnings(as.numeric(taxonID))) %>%
    filter(!is.na(taxonID))
  if (summarize == TRUE) {
    df <- df %>%
      group_by(taxonID, scientificName, taxonRank) %>%
      summarize(references = paste0(stri_remove_empty_na(references), collapse = "; "), remarks = paste0(stri_remove_empty_na(remarks), collapse = "; ")) %>%
      ungroup()
  }
  if (enhance == TRUE) {
    taxa <- taxon(df$taxonID)
    df <- df %>%
      select(-c(taxonRank)) %>%
      left_join(taxa, by = "taxonID") %>%
      mutate(scientificName.y = ifelse(is.na(scientificName.y), scientificName.x, scientificName.y)) %>%
      select(-scientificName.x) %>%
      rename(scientificName = scientificName.y)
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
