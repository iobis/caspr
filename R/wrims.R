#' Get the complete WRiMS checklist. This may not be complete due to missing macroalgae.
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
          select(taxonID, canonicalName, taxonomicStatus, rank, phylum, class, family) %>%
          mutate(taxonID = as.numeric(stringr::str_remove(taxonID, "urn:lsid:marinespecies.org:taxname:"))) %>%
          rename(scientificName = canonicalName, taxonRank = rank) %>%
          filter(taxonRank == "SPECIES" & taxonomicStatus == "ACCEPTED") %>%
          select(-taxonomicStatus) %>%
          mutate(taxonRank = "Species") %>%
          distinct() %>%
          as_tibble()

        #taxa <- taxon(result$taxonID)
        #result <- result %>%
        #  select(taxonID, canonicalName) %>%
        #  left_join(taxa, by = "taxonID") %>%
        #  mutate(
        #    scientificName = ifelse(is.na(scientificName), canonicalName, scientificName),
        #    taxonRank = "Species",
        #    wrims = TRUE
        #  )

        return(result)
      }
      counter <- counter + 1
    }
}

mr_geometries <- function(mrgid) {
  url <- glue::glue("https://marineregions.org/rest/getGazetteerGeometries.jsonld/{mrgid}/")
  message(url)
  geom <- jsonlite::fromJSON(url)$"mr:hasGeometry"$"gsp:asWKT"
  if (is.null(geom)) {
    return(st_sfc(st_polygon()))
  } else {
    geom_fixed <- stringr::str_replace(geom, "<.*>\\s+", "")
    geom_union <- st_union(st_as_sfc(geom_fixed))
    return(geom_union)
  }
}

wrims_distribution <- function(aphiaid, simplify = 0.5) {
  dist <- worrms::wm_distribution(aphiaid) %>%
    mutate(mrgid = stringr::str_replace_all(locationID, "http://marineregions.org/mrgid/", ""))
  mrgids <- unique(dist$mrgid)
  geom <- sapply(mrgids, mr_geometries)
  spatial <- tibble(geom = geom, mrgid = names(geom)) %>%
    as.data.frame()
  dist <- dist %>%
    left_join(spatial, by = "mrgid") %>%
    st_as_sf(crs = 4326)
  if (!is.na(simplify) && !is.null(simplify)) {
    sf_use_s2(FALSE)
    dist <- dist %>%
      st_simplify(dTolerance = simplify)
  }
  return(dist)
}

world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")

plot_wrims_dist <- function(dist) {
  ggplot() +
    geom_sf(data = world, fill = NA, size = 0.1, color = "#000000") +
    geom_sf(data = dist %>% filter(!is.na(establishmentMeans)), aes(fill = establishmentMeans, color = establishmentMeans), alpha = 0.5, size = 0.1) +
    ggpattern::geom_sf_pattern(data = dist %>% filter(is.na(establishmentMeans)), alpha = 1, size = 0.2, fill = NA, color = "#02a4cc", pattern_size = 0.1, pattern_color = "#02a4cc", pattern_spacing = 0.02, pattern_density = 0.5, pattern_fill = NA) +
    scale_fill_manual(values = c("Alien" = "#9c491f", "Native" = "#649c1f"), na.value = NA) +
    scale_color_manual(values = c("Alien" = "#9c491f", "Native" = "#649c1f"), na.value = NA) +
    theme_void() +
    theme(legend.position = "bottom")
}

plot_obis_dist <- function(aphiaid) {
  url <- glue::glue("https://api.obis.org/occurrence/grid/3?taxonid={aphiaid}")
  json <- httr::GET(URLencode(url)) %>% httr::content(as = "text")
  dist <- geojsonsf::geojson_sf(json)
  ggplot() +
    geom_sf(data = world, fill = NA, size = 0.1, color = "#000000") +
    geom_sf(data = dist, alpha = 0.7, size = 0.5, color = "#fc5603", fill = "#fc5603") +
    theme_void() +
    theme(legend.position = "bottom")
}
