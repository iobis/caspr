#' Get unioned world geometry.
#'
#' @export
get_world <- function() {
  rnaturalearth::ne_countries(scale = "large", returnclass = "sf") %>% st_union()
}
