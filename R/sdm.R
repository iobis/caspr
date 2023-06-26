sf_use_s2(FALSE)

get_layer_polygon <- function(layer, tmpdir = tempdir()) {
  layer[!is.na(layer)] <- 1
  fr <- tempfile(pattern = "file", fileext = ".tif", tmpdir = tmpdir)
  fp <- tempfile(pattern = "file", fileext = ".gpkg", tmpdir = tmpdir)
  unlink(fr)
  unlink(fp)
  raster::writeRaster(layer, fr, "GTiff", overwrite = TRUE)
  cmd <- glue("/Users/pieter/miniconda3/envs/pygdal/bin/gdal_polygonize.py {fr} -f GPKG {fp}")
  message(cmd)
  system(cmd)
  pol <- read_sf(fp) %>% st_set_crs(NA) %>% st_union() %>% st_set_crs(4326)
  return(pol)
}

extract_layer_values <- function(layer, points) {
  pol <- get_layer_polygon(layer)
  pol_buffered <- pol %>%
    st_set_crs(NA) %>%
    st_buffer(-0.01) %>% # TODO: adapt to raster resolution
    st_set_crs(4326)
  sf::sf_use_s2(FALSE)
  closest <- st_nearest_points(points, pol_buffered)
  sf::sf_use_s2(TRUE)
  pts = st_cast(closest, "POINT")
  snapped <- pts[seq(2, length(pts), 2)] %>% st_as_sf()
  env_extracted <- raster::extract(layer, snapped)
  return(env_extracted)
}

add_env_data <- function(points, env) {
  for (varname in names(env)) {
    env_extracted <- extract_layer_values(env[[varname]], points)
    points <- points %>%
      mutate({{varname}} := env_extracted)
  }
  return(points)
}

#' Get OBIS and GBIF occurrences.
#'
#' @export
get_all_occurrences <- function(aphiaid) {
  cols <- c("id" = NA_character_, "scientificName" = NA_character_, "year" = NA_integer_, "decimalLongitude" = NA_real_, "decimalLatitude" = NA_real_, "coordinateUncertaintyInMeters" = NA_character_, "establishmentMeans" = NA_character_, "degreeOfEstablishment" = NA_character_, "pathway" = NA_character_, "occurrenceStatus" = NA_character_, "basisOfRecord" = NA_character_)

  scientificname <- worrms::wm_record(aphiaid)$scientificname

  message("Fetching occurrences from OBIS")
  occ_obis <- robis::occurrence(scientificname) %>%
    mutate(source = "obis", year = date_year) %>%
    tibble::add_column(!!!cols[!names(cols) %in% names(.)]) %>%
    mutate(coordinateUncertaintyInMeters = as.numeric(coordinateUncertaintyInMeters)) %>%
    select(id, scientificName, year, decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters, establishmentMeans, degreeOfEstablishment, pathway, occurrenceStatus, basisOfRecord)

  message("Fetching occurrences from GBIF")
  occ_gbif <- rgbif::occ_search(scientificName = scientificname, limit = 1000000)$data %>%
    mutate(source = "gbif", id = key) %>%
    tibble::add_column(!!!cols[!names(cols) %in% names(.)]) %>%
    select(id, scientificName, year, decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters, establishmentMeans, degreeOfEstablishment, pathway, occurrenceStatus, basisOfRecord)

  occ <- bind_rows(occ_obis, occ_gbif) %>%
    filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
    filter(!(decimalLongitude == 0 & decimalLatitude == 0)) %>%
    filter(!stringr::str_detect(tidyr::replace_na(occurrenceStatus, ""), stringr::regex("absent", ignore_case = TRUE))) %>%
    filter(!stringr::str_detect(tidyr::replace_na(basisOfRecord, ""), stringr::regex("fossil", ignore_case = TRUE))) %>%
    group_by(scientificName, year, decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters, establishmentMeans, degreeOfEstablishment, pathway, occurrenceStatus, basisOfRecord) %>%
    summarize(ids = paste0(id, collapse = ";")) %>%
    ungroup() %>%
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = FALSE)

    return(occ)
}

#' Plot occurrences.
#'
#' @export
plot_occurrences <- function(occ) {
  if ("presence" %in% names(occ)) {
    ggplot() +
      geom_sf(data = get_world(), fill = NA, size = 0.1, color = alpha("#000000", 0.3)) +
      geom_sf(data = occ %>% st_as_sf() %>% filter(!as.logical(presence)), alpha = 1, color = "#9dc5cc", pch = 3, stroke = 1.2, size = 2) +
      geom_sf(data = occ %>% st_as_sf() %>% filter(as.logical(presence)), alpha = 1, color = "#fc5603", pch = 3, stroke = 1.2, size = 2) +
      theme_void()
  } else {
    ggplot() +
      geom_sf(data = get_world(), fill = NA, size = 0.1, color = alpha("#000000", 0.3)) +
      geom_sf(data = occ, alpha = 1, color = "#fc5603", pch = 3, stroke = 1.5, size = 3) +
      theme_void()
  }
}

#' Run the SDM.
#'
#' @export
run_model <- function(occ, bg = NULL, predictors = c("BO_sstmean", "BO_salinity")) {
  message("Running distribution model")

  # data preparation

  env <- sdmpredictors::load_layers(predictors, equalarea = FALSE, rasterstack = FALSE) %>%
    raster::stack()
  if (is.null(bg)) {
    bg <- generate_background_points()
  }
  occ <- occ |>
    add_env_data(env)
  bg <- bg |>
    add_env_data(env)
  df <- bind_rows(
      occ %>% mutate(presence = 1),
      bg %>% mutate(presence = 0)
    ) %>%
      as.data.frame()

  # glm model

  formula <- as.formula(paste0("presence ~ ", paste0(paste0("poly(", predictors, ", 2)"), collapse = " + ")))
  m_glm <- glm(formula, data = df, family = "binomial")

  # global prediction

  env_pred_glm <- raster::predict(env, m_glm, type = "response")

  return(env_pred_glm)

}

#' Plot prediction.
#'
#' @export
plot_prediction <- function(pred) {
  ggplot() +
    stars::geom_stars(data = stars::st_as_stars(pred), aes(x, y, fill = layer), alpha = 0.85) +
    geom_sf(data = get_world(), fill = alpha("#ffffff", 0.8), color = NA) +
    scale_fill_distiller(name = "risk score", palette = "BuPu", direction = 1, na.value = "#ffffff") +
    theme_void() +
    theme(legend.position = "bottom")
}

#' Plot risk score.
#'
#' @export
plot_risk <- function(risk, mask = TRUE) {
  ggplot() +
    geom_sf(data = risk, aes(fill = as.factor(risk_score)), alpha = 1, color = NA) +
    scale_fill_brewer(name = "risk score", palette = "Spectral", direction = -1) +
    theme_void() +
    theme(legend.position = "bottom")
}

#' Generate background points.
#'
#' @export
generate_background_points <- function() {
  message("Generating background points")
  land <- landr::get_land_polygons(simplified = TRUE)
  bounds <- st_bbox(c(xmin = -20026376.39, xmax = 20026376.39, ymax = 20048966.10 * 0.9, ymin = -20048966.10 * 0.9), crs = st_crs(3857)) %>%
    st_as_sfc()
  buffered <- land %>%
    sf::st_buffer(dist = units::set_units(50, km)) %>%
    st_crop(bounds) %>%
    st_union() %>%
    st_as_sf()
  buffer <- buffered %>%
    qgisprocess::qgis_run_algorithm("native:difference", INPUT = ., OVERLAY = land %>% select(-FID), .quiet = FALSE) %>%
    qgisprocess::qgis_output("OUTPUT") %>%
    read_sf()
  background <- st_sample(buffer %>% st_transform("ESRI:54030"), 1000) %>%
    st_as_sf() %>%
    st_transform(4326) %>%
    rename(geometry = x)
  return(background)
}
