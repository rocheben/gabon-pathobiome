# Helper used by all figure scripts: writes each figure twice, as a vector
# file (SVG, used in the manuscript) and as a 600 dpi PNG (for journals that
# require raster files). Width is set to the text width of the manuscript
# (6.4 in) so that fonts appear at their nominal size.
save_fig <- function(name, draw, width = 6.4, height = 4) {
  png(file.path("figures", paste0(name, ".png")), width = width, height = height, units = "in", res = 600, bg = "white")
  draw(); invisible(dev.off())
  svg(file.path("figures", paste0(name, ".svg")), width = width, height = height, bg = "white")
  draw(); invisible(dev.off())
}
