# scripts/build_index_mace_matched.R
# ============================================================
# Build a cleaned, matched dataset from:
#   - Sheet "statadatasheet"  -> select Registro, TAPSE1, PSAP1, SET1, gender, Age, bmi
#   - Sheet "MACE"            -> select REGISTRO + 5 day columns
#
# Outputs (saved in /data):
#   - index_mace_matched.csv
#   - index_mace_matched.rds
#
# Note: Working Directory pitfalls!
#  - If you run from project ROOT, the data is at: data/ECO PHASE.xlsx
#  - If you run from /scripts, the data is at:    ../data/ECO PHASE.xlsx
#  - The resolver below auto-detects the correct path.
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(readr)      # write_csv
  library(stringr)
})

# ---------- Helper: robust path resolution for ECO PHASE.xlsx ----------
find_input <- function() {
  # Prefer here::here() if available (works great with RStudio Projects)
  if (requireNamespace("here", quietly = TRUE)) {
    cand <- here::here("data", "ECO PHASE.xlsx")
    if (file.exists(cand)) return(cand)
  }
  # If running inside RStudio, anchor on the open script location
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(p)) {
      cand <- file.path(dirname(p), "..", "data", "ECO PHASE.xlsx")
      if (file.exists(cand)) return(cand)
    }
  }
  # If called via Rscript (looks for --file=...)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(file_arg)) {
    cand <- file.path(dirname(normalizePath(file_arg)), "..", "data", "ECO PHASE.xlsx")
    if (file.exists(cand)) return(cand)
  }
  # Common fallbacks
  for (cand in c(
    file.path("data", "ECO PHASE.xlsx"),         # running from root
    file.path("..", "data", "ECO PHASE.xlsx"),   # running from /scripts
    "ECO PHASE.xlsx"                              # running inside /data
  )) if (file.exists(cand)) return(cand)
  
  stop(
    "Could not find 'ECO PHASE.xlsx'.\n",
    "wd: ", getwd(), "\n",
    "Tried: data/, ../data/, and current dir.\n",
    "Ensure the file (with the space in its name) exists in the /data folder."
  )
}

INPUT_XLSX <- find_input()
message("✅ Using Excel file at: ", normalizePath(INPUT_XLSX, winslash = "/"))

# Output dir = the 'data' folder that contains the Excel
OUT_DIR <- dirname(INPUT_XLSX)
OUT_CSV <- file.path(OUT_DIR, "index_mace_matched.csv")
OUT_RDS <- file.path(OUT_DIR, "index_mace_matched.rds")

# ---------- Helper: clean → numeric ----------
# - Treat "." and "" as NA
# - Convert commas to dots (e.g., "12,3" → "12.3")
# - Coerce quietly to numeric
to_numeric_clean <- function(x) {
  v <- as.character(x)
  v[v %in% c(".", "")] <- NA
  v <- gsub(",", ".", v, fixed = FALSE)
  suppressWarnings(as.numeric(v))
}

# ---------- 1) Read sheets ----------
# Note: we also specify na= to catch obvious missing marks at import time.
statadatasheet <- read_excel(
  INPUT_XLSX, sheet = "statadatasheet",
  na = c("", ".", "NA")
)
mace <- read_excel(
  INPUT_XLSX, sheet = "MACE",
  na = c("", ".", "NA")
)

# ---------- 2) Build index_df (statadatasheet) ----------
index_cols <- c("Registro", "TAPSE1", "PSAP1", "SET1", "gender", "Age", "bmi")
missing_index_cols <- setdiff(index_cols, names(statadatasheet))
if (length(missing_index_cols)) {
  warning("Missing columns in statadatasheet: ", paste(missing_index_cols, collapse = ", "))
}

index_df <- statadatasheet %>%
  select(any_of(index_cols)) %>%
  mutate(Registro = trimws(as.character(Registro))) %>%
  mutate(
    TAPSE1 = to_numeric_clean(TAPSE1),
    PSAP1  = to_numeric_clean(PSAP1),
    SET1   = to_numeric_clean(SET1)
  ) %>%
  # Keep only rows where ALL three echo vars are numeric (no NAs after coercion)
  filter(!is.na(TAPSE1) & !is.na(PSAP1) & !is.na(SET1)) %>%
  # If duplicated IDs exist, keep first occurrence
  distinct(Registro, .keep_all = TRUE)

# ---------- 3) Reduce MACE to requested columns and coerce to numeric ----------
mace_keep <- c(
  "REGISTRO",
  "MACE_días",
  "MUERTECV_días",
  "REINFARTO_días",
  "CHOQUE_días",
  "SANGRADOMAYOR_días"
)
missing_mace_cols <- setdiff(mace_keep, names(mace))
if (length(missing_mace_cols)) {
  warning("Missing columns in MACE: ", paste(missing_mace_cols, collapse = ", "))
}

mace_sub <- mace %>%
  select(any_of(mace_keep)) %>%
  rename(Registro = REGISTRO) %>%
  mutate(Registro = trimws(as.character(Registro)))

mace_day_cols <- c(
  "MACE_días", "MUERTECV_días", "REINFARTO_días", "CHOQUE_días", "SANGRADOMAYOR_días"
)

mace_sub_num <- mace_sub %>%
  mutate(across(all_of(mace_day_cols), to_numeric_clean)) %>%
  # Keep only rows where ALL 5 day columns are numeric
  filter(if_all(all_of(mace_day_cols), ~ !is.na(.))) %>%
  distinct(Registro, .keep_all = TRUE)

# ---------- 4) Match by Registro (inner join) ----------
index_mace_matched <- index_df %>%
  inner_join(mace_sub_num, by = "Registro")

# ---------- 4.5) Add ratio and rename final base ----------
# (Uses the object created in Step 4: index_mace_matched)
index_base <- index_mace_matched %>%
  mutate(TAPSE1_PSAP1 = TAPSE1 / PSAP1)

# ---------- 5) Save outputs (now using index_base) ----------
OUT_CSV <- file.path(OUT_DIR, "index_base.csv")
OUT_RDS <- file.path(OUT_DIR, "index_base.rds")

write_csv(index_base, OUT_CSV)
saveRDS(index_base, OUT_RDS)

cat(
  "\n✅ Saved outputs:",
  "\n  - CSV:", normalizePath(OUT_CSV, winslash = "/"),
  "\n  - RDS:", normalizePath(OUT_RDS, winslash = "/"),
  "\nRows:", nrow(index_base),
  " | Columns:", ncol(index_base), "\n\n"
)

