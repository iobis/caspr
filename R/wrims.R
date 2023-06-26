sf_use_s2(FALSE)

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

#' Get geometries by MRGID.
#'
#' @export
mr_geometries <- function(mrgids) {
  message("Fetching MarineRegions geometries")
  sapply(mrgids, function(mrgid) {
    url <- glue::glue("https://marineregions.org/rest/getGazetteerGeometries.jsonld/{mrgid}/")
    geom <- jsonlite::fromJSON(url)$"mr:hasGeometry"$"gsp:asWKT"
    if (is.null(geom)) {
      return(st_sfc(st_polygon()))
    } else {
      geom_fixed <- stringr::str_replace(geom, "<.*>\\s+", "")
      geom_fixed <- geom_fixed[geom_fixed != ""]
      geom_union <- st_union(st_as_sfc(geom_fixed))
      return(geom_union)
    }
  })
}

#' Get WRiMS distribution records.
#'
#' @export
wrims_distribution <- function(aphiaid, simplify = FALSE, cleanup = FALSE) {
  message("Fetching WRiMS distribution")
  sf_use_s2(FALSE)

  dist <- worrms::wm_distribution(aphiaid) %>%
    mutate(mrgid = stringr::str_replace_all(locationID, "http://marineregions.org/mrgid/", ""))
  mrgids <- unique(dist$mrgid)
  geom <- mr_geometries(mrgids)
  spatial <- tibble(geom = geom, mrgid = names(geom)) %>%
    as.data.frame()
  dist_spatial <- dist %>%
    left_join(spatial, by = "mrgid") %>%
    st_as_sf(crs = 4326)
  if (simplify) {
    sf_use_s2(FALSE)
    dist_spatial <- dist_spatial %>%
      st_simplify(dTolerance = simplify)
  }
  if (cleanup) {
    dist_spatial <- dist_spatial %>%
      filter(!st_is(., "LINESTRING") & !st_is(., "MULTILINESTRING") & !st_is_empty(.))
  }

  return(dist_spatial)
}

#' Create summary of WRiMS distribution.
#'
#' @export
wrims_summary <- function(wrims_dist, unknown = 0.75, add_background = FALSE, simplify = 0.1, hole_max_area = 100000000, union_buffer = 0.01, buffer = 1) {
  message("Summarizing WRiMS distribution")
  sf_use_s2(FALSE)

  dist <- wrims_dist %>%
    st_buffer(dist = union_buffer) %>%
    group_by(establishmentMeans) %>%
    summarize() %>%
    nngeo::st_remove_holes(max_area = hole_max_area)

  if (simplify) {
    sf_use_s2(FALSE)
    dist <- dist %>%
      st_simplify(dTolerance = simplify)
  }

  dist_native <- dist %>%
    filter(establishmentMeans == "Native") %>%
    st_union()
  dist_alien <- dist %>%
    filter(establishmentMeans == "Alien") %>%
    st_union() %>%
    st_difference(dist_native)

  summary <- data.frame(
    establishmentMeans = c("native", "alien"),
    geometry = c(dist_native, dist_alien)
  ) %>% st_as_sf()

  if (buffer) {
    summary_buffer <- summary %>%
      st_collection_extract() %>%
      st_set_crs(NA) %>%
      regional_seas(
        group = "establishmentMeans",
        dist = buffer,
        density = 4
      )
    summary <- summary %>%
      bind_rows(summary_buffer) %>%
      group_by(establishmentMeans) %>%
      summarize()
  }

  summary <- summary %>%
    mutate(risk_score = ifelse(establishmentMeans == "native", 0, 1))

  if (add_background) {
    background <- st_as_sfc(st_bbox(c(xmin = -180, xmax = 180, ymax = 90, ymin = -90), crs = st_crs(4326))) %>%
      st_difference(summary %>% st_union())
    summary <- bind_rows(
      summary,
      data.frame(risk_score = unknown, establishmentMeans = NA, geometry = background) %>% st_as_sf()
    )
  }

  summary <- summary %>%
    st_crop(xmin = -180, xmax = 180, ymax = 90, ymin = -90)

  return(summary)
}

#' Create summary of WRiMS distribution.
#'
#' @export
wrims_summary_for_aphiaid <- function(aphiaid, unknown = 0.75, add_background = FALSE, simplify = 0.1, hole_max_area = 100000000, union_buffer = 0.01, buffer = 1) {
  wrims_dist <- wrims_distribution(aphiaid)
  wrims_summary(wrims_dist, unknown = unknown, add_background = add_background, simplify = simplify, hole_max_area = hole_max_area, union_buffer = union_buffer, buffer = buffer)
}

#' Plot summary of WRiMS distribution.
#'
#' @export
plot_wrims_summary <- function(summary) {
  ggplot() +
    geom_sf(data = summary, aes(fill = as.factor(risk_score)), alpha = 1) +
    scale_fill_brewer(name = "risk score", palette = "BuPu", direction = -1) +
    theme_void() +
    theme(legend.position = "bottom")
}

#' Plot WRiMS distribution records.
#'
#' @export
plot_wrims_dist <- function(dist) {
  p <- ggplot() +
    theme_void() +
    theme(legend.position = "bottom")
  p %>%
    add_wrims_dist(dist)
}

#' Add WRiMS distribution records to a map.
#'
#' @export
add_wrims_dist <- function(p, dist) {
  p +
    geom_sf(data = get_world(), fill = NA, size = 0.1, color = "#000000") +
    geom_sf(data = dist %>% filter(!is.na(establishmentMeans)), aes(fill = establishmentMeans, color = establishmentMeans), alpha = 0.5, size = 0.1) +
    ggpattern::geom_sf_pattern(data = dist %>% filter(is.na(establishmentMeans)), alpha = 1, size = 0.2, fill = NA, color = "#02a4cc", pattern_size = 0.1, pattern_color = "#02a4cc", pattern_spacing = 0.02, pattern_density = 0.5, pattern_fill = NA) +
    scale_fill_manual(values = c("Alien" = "#9c491f", "Native" = "#649c1f"), na.value = NA) +
    scale_color_manual(values = c("Alien" = "#9c491f", "Native" = "#649c1f"), na.value = NA)
}

#' Get OBIS distribution
#'
#' @export
obis_distribution <- function(aphiaid, res = 3) {
  url <- glue::glue("https://api.obis.org/occurrence/grid/{res}?taxonid={aphiaid}")
  json <- httr::GET(URLencode(url)) %>% httr::content(as = "text")
  geojsonsf::geojson_sf(json)
}

#' Plot OBIS distribution.
#'
#' @export
plot_obis_dist <- function(dist) {
  ggplot() +
    geom_sf(data = get_world(), fill = NA, size = 0.1, color = "#000000") +
    theme_void() +
    theme(legend.position = "bottom") %>%
  add_obis_dist(dist)
}

#' Plot OBIS distribution.
#'
#' @export
add_obis_dist <- function(p, dist) {
  p +
    geom_sf(data = get_world(), fill = NA, size = 0.1, color = "#000000") +
    geom_sf(data = dist, alpha = 0.7, size = 0.5, color = "#fc5603", fill = "#fc5603") +
    theme_void() +
    theme(legend.position = "bottom")
}
