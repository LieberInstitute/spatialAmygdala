library(here)
library(ggplot2)
library(patchwork)
library(jsonlite)
library(png)
library(SpatialExperiment)

plots_dir  <- here("plots", "VisiumHD", "07_img_cropping")
assets_dir <- here("processed-data", "VisiumHD", "07_img_cropping")
dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)

PAD  <- 20
GENE <- "SLC17A7"

# ---- load finalized object --------------------------------------------------
spe <- readRDS(here("code", "VisiumHD", "Samui_app",
                    "spe_cells_with_domain.rds"))
samples <- unique(spe$sample_id)

# ---- restore full-res coords ------------------------------------------------
# The finalized object has per-sample median-centred coords (ranges ~ +/-13000),
# which put every cell near the image origin. Pull the original pixel coords
# back in from the build object.
#
# Cell ids are NOT unique across samples (every sample has its own
# cellid_...-1 series), so matching on colnames alone silently returns the
# first sample's coordinates for every sample. Match on sample_id + cell id.
raw <- readRDS(here("processed-data", "VisiumHD", "02_build_spe",
                    "spe_cells_combined.rds"))

stopifnot(min(spatialCoords(raw)) >= 0)      # build object must be uncentred

key_spe <- paste(spe$sample_id, colnames(spe), sep = "|")
key_raw <- paste(raw$sample_id, colnames(raw), sep = "|")
stopifnot(!anyDuplicated(key_raw))

idx <- match(key_spe, key_raw)
stopifnot(!anyNA(idx))

spatialCoords(spe) <- spatialCoords(raw)[idx, , drop = FALSE]
rm(raw)

# per-sample ranges should differ between samples; all-identical means the
# match collapsed onto one sample again
print(sapply(samples, function(s)
  apply(spatialCoords(spe)[spe$sample_id == s, ], 2, range)))

# ---- add hires images -------------------------------------------------------
spatial_dirs <- c(
  Br9280_CeA = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-W369TJK_A1/outs/segmented_outputs/spatial",
  Br8325_MeA = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-937FVHX_A1/outs/segmented_outputs/spatial",
  Br8325_CeA = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-937FVHX_D1/outs/segmented_outputs/spatial",
  Br9280_ITC = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-HW9VGBW_A1/outs/segmented_outputs/spatial",
  Br9280_MeA = "/dcs04/lieber/marmaypag/spatialAMY_LIBD4125/spatialAmygdala/processed-data/VisiumHD/01_spaceranger/H1-HW9VGBW_D1/outs/segmented_outputs/spatial"
)
stopifnot(setequal(names(spatial_dirs), samples))

for (s in samples) {
  d   <- spatial_dirs[[s]]
  sfj <- fromJSON(txt = readLines(file.path(d, "scalefactors_json.json"), warn = FALSE))
  if (!any(imgData(spe)$sample_id == s & imgData(spe)$image_id == "hires")) {
    spe <- addImg(spe, sample_id = s, image_id = "hires",
                  imageSource = file.path(d, "tissue_hires_image.png"),
                  scaleFactor = sfj$tissue_hires_scalef, load = TRUE)
  }
}

sf_of <- function(s) {
  idt <- imgData(spe)
  sf  <- idt$scaleFactor[idt$sample_id == s & idt$image_id == "hires"]
  stopifnot(length(sf) == 1)
  sf
}

# ---- QC before cropping -----------------------------------------------------
plot_full <- function(s) {
  sf <- sf_of(s)
  img <- imgRaster(spe, sample_id = s, image_id = "hires")
  H <- nrow(img); W <- ncol(img)
  keep <- spe$sample_id == s
  xy <- spatialCoords(spe)[keep, , drop = FALSE]
  df <- data.frame(x = xy[, 1] * sf, y = H - xy[, 2] * sf,
                   expr = counts(spe)[GENE, keep])
  df <- df[order(df$expr), ]
  ggplot(df, aes(x, y, colour = expr)) +
    annotation_raster(img, xmin = 0, xmax = W, ymin = 0, ymax = H) +
    geom_point(size = 0.05, alpha = 0.6) +
    scale_colour_viridis_c(option = "magma", name = GENE) +
    coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE) +
    ggtitle(paste(s, "full")) + theme_void()
}

# ---- crop, write PNGs, record offsets ---------------------------------------
crops <- do.call(rbind, lapply(samples, function(s) {

  sf  <- sf_of(s)
  m   <- as.matrix(imgRaster(spe, sample_id = s, image_id = "hires"))
  rgb <- col2rgb(m) / 255
  arr <- array(0, dim = c(nrow(m), ncol(m), 3))
  for (k in 1:3) arr[, , k] <- matrix(rgb[k, ], nrow(m), ncol(m))

  xy <- spatialCoords(spe)[spe$sample_id == s, , drop = FALSE] * sf
  x0 <- max(floor(min(xy[, 1])) - PAD, 0); x1 <- min(ceiling(max(xy[, 1])) + PAD, ncol(m))
  y0 <- max(floor(min(xy[, 2])) - PAD, 0); y1 <- min(ceiling(max(xy[, 2])) + PAD, nrow(m))

  writePNG(arr[(y0 + 1):y1, (x0 + 1):x1, , drop = FALSE],
           file.path(assets_dir, paste0("sample_", s, "_cropped_hires.png")))

  data.frame(sample_id = s, scaleFactor = sf, x0 = x0, y0 = y0,
             W = x1 - x0, H = y1 - y0)
}))
print(crops)

full_plots <- lapply(samples, plot_full)
names(full_plots) <- samples

# ---- write back so the object is self-contained ------------------------------
# imgData has scaleFactor (a multiply) but no offset slot, so the crop origin
# comes off spatialCoords, in full-res units.
for (s in samples) {
  cr   <- crops[crops$sample_id == s, ]
  keep <- spe$sample_id == s

  spatialCoords(spe)[keep, 1] <- spatialCoords(spe)[keep, 1] - cr$x0 / cr$scaleFactor
  spatialCoords(spe)[keep, 2] <- spatialCoords(spe)[keep, 2] - cr$y0 / cr$scaleFactor

  spe <- rmvImg(spe, sample_id = s, image_id = "hires")
  spe <- addImg(spe, sample_id = s, image_id = "hires",
                imageSource = file.path(assets_dir,
                                        paste0("sample_", s, "_cropped_hires.png")),
                scaleFactor = cr$scaleFactor, load = TRUE)

  # lowres scaleFactor still refers to uncropped coords -> would plot wrong
  if (any(imgData(spe)$sample_id == s & imgData(spe)$image_id == "lowres")) {
    spe <- rmvImg(spe, sample_id = s, image_id = "lowres")
  }
}

saveRDS(spe, file.path(assets_dir, "spe_cropped.rds"))

# ---- QC after ---------------------------------------------------------------
plot_crop <- function(s) {
  sf <- sf_of(s)
  img <- imgRaster(spe, sample_id = s, image_id = "hires")
  H <- nrow(img); W <- ncol(img)
  keep <- spe$sample_id == s
  xy <- spatialCoords(spe)[keep, , drop = FALSE]
  df <- data.frame(x = xy[, 1] * sf, y = H - xy[, 2] * sf,
                   expr = counts(spe)[GENE, keep])
  df <- df[order(df$expr), ]
  ggplot(df, aes(x, y, colour = expr)) +
    annotation_raster(img, xmin = 0, xmax = W, ymin = 0, ymax = H) +
    geom_point(size = 0.05, alpha = 0.6) +
    scale_colour_viridis_c(option = "magma", name = GENE) +
    coord_fixed(xlim = c(0, W), ylim = c(0, H), expand = FALSE) +
    ggtitle(paste(s, "cropped")) + theme_void()
}

pdf(file.path(plots_dir, "crop_check_hires.pdf"), width = 14, height = 7)
for (s in samples) print(full_plots[[s]] | plot_crop(s))
dev.off()

for (s in samples) {
  img <- imgRaster(spe, sample_id = s, image_id = "hires")
  xy  <- spatialCoords(spe)[spe$sample_id == s, ] * sf_of(s)
  cat(sprintf("%-12s img %d x %d | x [%.0f, %.0f] y [%.0f, %.0f]\n",
              s, ncol(img), nrow(img),
              min(xy[,1]), max(xy[,1]), min(xy[,2]), max(xy[,2])))
}


# ====== Export the cropped hires images for Samui app ==========================================

for (s in unique(spe$sample_id)) magick::image_write(magick::image_read(imgRaster(spe, sample_id = s, image_id = "hires")), file.path(plots_dir, paste0(s, "_hires.tif")), format = "tiff")