#' Plot distances from a point location on a map.
#'
#' @export
plot_location_distances <- function(coord, distances) {
  ggplot() +
    theme_void() %>%
    add_location_distances(coord, distances)
}

#' Add distances from a point location to a map.
#'
#' @export
add_location_distances <- function(p, coord, distances) {
  point <- st_sfc(st_point(coord), crs = 4326)
  p <- p +
    geom_sf(data = point, fill = NA, size = 3, color = "#93e312", pch = 4, stroke = 2)

  #for (distance in distances) {
  #  buffer <- st_buffer(point, distance, max_cells = 10000) %>%
  #    st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
  #  p <- p + geom_sf(data = buffer, fill = NA, color = "#f54260", pch = 21, stroke = 1)
  #}

  # override until GDAL update?
  # https://github.com/r-spatial/sf/issues/1999
  wkt <- readr::read_file("wkt.txt")
  buffer <- st_as_sfc(wkt, crs = 4326)
  p <- p + geom_sf(data = buffer, fill = NA, color = "#f54260", pch = 21, stroke = 1)

  p
}

#' @export
risk_analysis_map <- function(aphiaid, coord) {
  obis_dist <- obis_distribution(aphiaid, 3)
  wrims_dist <- wrims_distribution(aphiaid)
  plot_wrims_dist(wrims_dist) %>%
    add_obis_dist(obis_dist) %>%
    add_location_distances(coord, c(1000000, 5000000))
}

#' Plot summary of WRiMS distribution.
#'
#' @export
plot_risk_score <- function(risk) {
  ggplot() +
    geom_sf(data = risk, aes(fill = risk_score), alpha = 1, color = NA) +
    geom_sf(data = get_world(), alpha = 1, color = "#ffffff", fill = NA, size = 0.1) +
    scale_fill_distiller(name = "risk score", palette = "BuPu", direction = 1) +
    theme_void() +
    theme(legend.position = "bottom")
}

#' Calculate risk from WRiMS summary and reclassified predictions.
#'
#' @export
calculate_risk <- function(wrims_summ, pred_reclassified) {
  sf_use_s2(FALSE)
  wrims_summ %>%
    select(risk_score) %>%
    st_intersection(pred_reclassified) %>%
    mutate(risk_score = risk_score * risk_score.1) %>%
    group_by(risk_score) %>%
    summarize()
}

#' @export
reclassify_prediction <- function(pred) {
  pred %>%
    raster::reclassify(
      matrix(c(
        0, 0.2, 0,
        0.2, 0.5, 0.75,
        0.5, Inf, 1
      )),
      include.lowest = TRUE
    ) %>%
    stars::st_as_stars() %>%
    st_as_sf(as_points = FALSE, merge = TRUE) %>%
    rename(c(risk_score = layer))
}

#' Run risk analysis model.
#'
#' @export
risk_analysis <- function(aphiaid, coord) {
  sf_use_s2(FALSE)

  occ <- get_all_occurrences(aphiaid)
  bg <- generate_background_points()
  pred <- run_model(occ, bg)
  pred_reclassified <- reclassify_prediction(pred)

  wrims_dist <- wrims_distribution(aphiaid, simplify = 0.1)
  wrims_summ <- wrims_summary_for_aphiaid(aphiaid, add_background = TRUE)
  risk <- calculate_risk(wrims_summ, pred_reclassified)

  risk_score <- risk$risk_score[st_nearest_feature(st_sfc(st_point(coord), crs = 4326), risk)]
  points <- bind_rows(occ %>% mutate(presence = 1), bg %>% mutate(presence = 0)) %>% as.data.frame()

  return(list(
    caspr_version = as.character(packageVersion("caspr")),
    time = Sys.time(),
    aphiaid = aphiaid,
    location = coord,
    risk = risk,
    risk_score = risk_score,
    points = points,
    pred = pred,
    pred_reclassified = pred_reclassified,
    wrims_dist = wrims_dist,
    wrims_summ = wrims_summ
  ))
}

#' @export
image_meta <- function(filename, name, s3_root) {
  list(
    "@id" = unbox(filename),
    "@type" = unbox("File"),
    "name" = unbox(name),
    "encodingFormat" = unbox("image/png"),
    "url" = unbox(glue("{s3_root}{filename}"))
  )
}

#' Upload risk analysis to S3
#'
#' @export
upload_risk_analysis <- function(ra) {
  tmp <- tempdir()

  time <- format(ra$time, format = "%Y-%m-%dT%H:%M:%S")
  location <- paste0(round(ra$location, digits = 1), collapse = "_")
  aphiaid <- ra$aphiaid
  prefix <- glue("{aphiaid}_{time}_{location}/")
  s3_root <- glue("https://pacman-risk-analysis.s3.eu-central-1.amazonaws.com/{prefix}")
  caspr_version <- as.character(ra$caspr_version)
  metadata_file <- glue("{tmp}/ro-crate-metadata.json")
  risk_analysis_file <- glue("{tmp}/risk_analysis.json")

  # generate metadata

  meta <- create_rocrate_meta(title = glue("Risk analysis {aphiaid}"), date = time, caspr_version = caspr_version, s3_root = s3_root)

  # add images

  filename <- "wrims_distribution.png"
  meta[["@graph"]] <- append(meta[["@graph"]], list(image_meta(filename, "WRiMS distribution", s3_root)))
  meta[["@graph"]][[2]][["hasPart"]] <- append(meta[["@graph"]][[2]][["hasPart"]], list(list("@id" = unbox(filename))))
  plot_wrims_dist(ra$wrims_dist) %>% ggsave(glue("{tmp}/{filename}"), plot = ., width = 10, height = 6, dpi = 300, scale = 1.2, bg = "white")

  filename <- "wrims_summary.png"
  meta[["@graph"]] <- append(meta[["@graph"]], list(image_meta(filename, "WRiMS distribution summary", s3_root)))
  meta[["@graph"]][[2]][["hasPart"]] <- append(meta[["@graph"]][[2]][["hasPart"]], list(list("@id" = unbox(filename))))
  plot_wrims_summary(ra$wrims_summ) %>% ggsave(glue("{tmp}/{filename}"), plot = ., width = 10, height = 6, dpi = 300, scale = 1.2, bg = "white")

  filename <- "occurrences.png"
  meta[["@graph"]] <- append(meta[["@graph"]], list(image_meta(filename, "Occurrences from OBIS and GBIF, background points", s3_root)))
  meta[["@graph"]][[2]][["hasPart"]] <- append(meta[["@graph"]][[2]][["hasPart"]], list(list("@id" = unbox(filename))))
  plot_occurrences(ra$points) %>% ggsave(glue("{tmp}/{filename}"), plot = ., width = 10, height = 6, dpi = 300, scale = 1.2, bg = "white")

  filename <- "prediction.png"
  meta[["@graph"]] <- append(meta[["@graph"]], list(image_meta(filename, "Distribution model prediction", s3_root)))
  meta[["@graph"]][[2]][["hasPart"]] <- append(meta[["@graph"]][[2]][["hasPart"]], list(list("@id" = unbox(filename))))
  plot_prediction(ra$pred) %>% ggsave(glue("{tmp}/{filename}"), plot = ., width = 10, height = 6, dpi = 300, scale = 1.2, bg = "white")

  filename <- "prediction_reclassified.png"
  meta[["@graph"]] <- append(meta[["@graph"]], list(image_meta(filename, "Distribution model prediction (reclassified)", s3_root)))
  meta[["@graph"]][[2]][["hasPart"]] <- append(meta[["@graph"]][[2]][["hasPart"]], list(list("@id" = unbox(filename))))
  plot_wrims_summary(ra$pred_reclassified) %>% ggsave(glue("{tmp}/{filename}"), plot = ., width = 10, height = 6, dpi = 300, scale = 1.2, bg = "white")

  filename <- "risk_score.png"
  meta[["@graph"]] <- append(meta[["@graph"]], list(image_meta(filename, "Risk score", s3_root)))
  meta[["@graph"]][[2]][["hasPart"]] <- append(meta[["@graph"]][[2]][["hasPart"]], list(list("@id" = unbox(filename))))
  plot_risk_score(ra$risk) %>% ggsave(glue("{tmp}/{filename}"), plot = ., width = 10, height = 6, dpi = 300, scale = 1.2, bg = "white")

  # write JSON files

  meta_json <- jsonlite::prettify(jsonlite::toJSON(meta))
  ra_json <- jsonlite::toJSON(ra[c("aphiaid", "risk_score", "caspr_version", "time", "location")], auto_unbox = TRUE)
  write(meta_json, metadata_file)
  write(ra_json, risk_analysis_file)

  # clean up

  unlink(list.files(tmp, pattern = "rs-graphics", full.names = TRUE), recursive = TRUE)

  # sync

  s3sync(tmp, prefix = prefix, bucket = "pacman-risk-analysis", direction = "upload", region = Sys.getenv("AWS_S3_REGION_NAME"))

  return(glue("{s3_root}ro-crate-metadata.json"))
}

#' Create RO-Crate metadata JSON.
#'
#' @export
create_rocrate_meta <- function(title, date, caspr_version, s3_root) {
  list(
    "@context" = unbox("https://w3id.org/ro/crate/1.1/context"),
    "@graph" = list(
      list(
        "@type" = unbox("CreativeWork"),
        "@id" = unbox("ro-crate-metadata.json"),
        "conformsTo" = list("@id" = unbox("https://w3id.org/ro/crate/1.1")),
        "about" = list("@id" = unbox("./"))
      ),
      list(
        "@id" = unbox("./"),
        "@type" = unbox("Dataset"),
        "name" = unbox(title),
        "description" = unbox(title),
        "datePublished" = unbox(date),
        "license" = unbox("https://creativecommons.org/licenses/by/4.0/"),
        "hasPart" = list(
          list("@id" = unbox("risk_analysis.json"))
        )
      ),
      list(
        "@id" = unbox("risk_analysis.json"),
        "@type" = unbox("File"),
        "name" = unbox("Risk analysis results"),
        "encodingFormat" = unbox("application/json"),
        "url" = unbox(glue("{s3_root}risk_analysis.json"))
      ),
      list(
        "@id" = unbox("https://github.com/iobis/caspr"),
        "@type" = unbox("SoftwareApplication"),
        "url" = unbox("https://github.com/iobis/caspr"),
        "name" = unbox("caspr"),
        "version" = unbox(caspr_version)
      ),
      list(
        "@id" = unbox("risk_analysis"),
        "@type" = unbox("CreateAction"),
        "@agent" = list("@id" = unbox("https://github.com/iobis/caspr")),
        "name" = unbox("Risk analysis"),
        "result" = list("@id" = unbox("risk_analysis.json"))
      )
    )
  )
}
