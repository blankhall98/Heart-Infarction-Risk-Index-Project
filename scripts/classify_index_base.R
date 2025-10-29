# scripts/classify_index_base.R
# ============================================================
# Loads data/index_base.(rds|csv), classifies rows into A/B/C/D
# based on SET1 (280) and TAPSE1_PSAP1 (0.32), overwrites the base
# in /data with the new column, and saves a scatter plot in /images.
#
# Classes:
#  A: SET1 > 280  & TAPSE1_PSAP1 > 0.32
#  B: SET1 < 280  & TAPSE1_PSAP1 > 0.32
#  C: SET1 > 280  & TAPSE1_PSAP1 < 0.32
#  D: SET1 < 280  & TAPSE1_PSAP1 < 0.32
#
# Note on WD pitfalls:
# - From project root: files live under data/ and images/
# - From scripts/: use ../data and ../images
# The resolver below finds paths automatically.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# ---------- Path helpers (auto-locate /data and /images) ----------
find_data_dir <- function() {
  # try common locations relative to WD
  for (p in c("data", "../data", "./")) {
    if (dir.exists(p) && any(tolower(list.files(p)) %in% c("index_base.rds","index_base.csv")))
      return(normalizePath(p, winslash = "/"))
  }
  stop("Could not find /data with index_base.(rds|csv). WD: ", getwd())
}

DATA_DIR <- find_data_dir()
IMG_DIR  <- normalizePath(file.path(DATA_DIR, "..", "images"), winslash = "/", mustWork = FALSE)
if (!dir.exists(IMG_DIR)) dir.create(IMG_DIR, recursive = TRUE, showWarnings = FALSE)

# Input/Output files
RDS_PATH <- file.path(DATA_DIR, "index_base.rds")
CSV_PATH <- file.path(DATA_DIR, "index_base.csv")
IMG_PATH <- file.path(IMG_DIR,  "index_base_quadrants.png")

# ---------- Load base (prefer RDS, else CSV) ----------
if (file.exists(RDS_PATH)) {
  index_base <- readRDS(RDS_PATH)
} else if (file.exists(CSV_PATH)) {
  index_base <- readr::read_csv(CSV_PATH, show_col_types = FALSE)
} else {
  stop("index_base.rds or index_base.csv not found in: ", DATA_DIR)
}

# Sanity: ensure needed columns exist & numeric
needed <- c("SET1", "TAPSE1_PSAP1")
missing <- setdiff(needed, names(index_base))
if (length(missing)) stop("Missing columns in index_base: ", paste(missing, collapse = ", "))

to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
index_base <- index_base %>%
  mutate(
    SET1          = to_num(SET1),
    TAPSE1_PSAP1  = to_num(TAPSE1_PSAP1)
  )

# ---------- Classify A/B/C/D ----------
set_thr   <- 280
ratio_thr <- 0.32
eps <- 1e-9  # tolerance for floating point

index_base <- index_base %>%
  mutate(
    ClassABCD = case_when(
      SET1 >= set_thr - eps & TAPSE1_PSAP1 >= ratio_thr - eps ~ "A",
      SET1 <  set_thr - eps & TAPSE1_PSAP1 >= ratio_thr - eps ~ "B",
      SET1 >= set_thr - eps & TAPSE1_PSAP1 <  ratio_thr - eps ~ "C",
      SET1 <  set_thr - eps & TAPSE1_PSAP1 <  ratio_thr - eps ~ "D",
      TRUE ~ NA_character_
    )
  )

# Optional: if you want to drop rows that landed in NA (exactly on a threshold)
# index_base <- index_base %>% filter(!is.na(ClassABCD))

# Quick console summary
print(table(index_base$ClassABCD, useNA = "ifany"))

# ---------- Save updated base (overwrite same files) ----------
saveRDS(index_base, RDS_PATH)
readr::write_csv(index_base, CSV_PATH)

cat("\n✅ Saved updated base with ClassABCD to:",
    "\n  -", RDS_PATH,
    "\n  -", CSV_PATH, "\n")

# ---------- Plot & save image ----------
p <- ggplot(index_base, aes(x = SET1, y = TAPSE1_PSAP1, color = ClassABCD)) +
  geom_point(alpha = 0.85, size = 2) +
  geom_vline(xintercept = set_thr, linetype = "dashed") +
  geom_hline(yintercept = ratio_thr, linetype = "dashed") +
  scale_color_manual(values = c(A = "#1f77b4", B = "#2ca02c", C = "#ff7f0e", D = "#d62728"),
                     na.translate = TRUE, na.value = "grey50") +
  labs(
    title = "Classification by SET1 and TAPSE1/PSAP1",
    x = "SET1",
    y = "TAPSE1/PSAP1",
    color = "Class"
  ) +
  theme_minimal(base_size = 12)

ggsave(filename = IMG_PATH, plot = p, width = 8, height = 6, dpi = 300)

cat("🖼️  Saved plot to: ", IMG_PATH, "\n")

        