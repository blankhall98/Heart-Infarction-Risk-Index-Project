# scripts/km_sub_experiments.R
# ============================================================
# Sub-experiments for MACE only, given sparse C & D:
#   1) A vs B
#   2) A vs B vs CD (merge C & D into "CD")
#
# Inputs (auto-located): data/index_base.rds or data/index_base.csv
# Outputs (white background):
#   images/sub-experiment/km_MACE_A_vs_B.png
#   images/sub-experiment/km_MACE_A_vs_B_risk.png
#   images/sub-experiment/km_MACE_A_B_CD.png
#   images/sub-experiment/km_MACE_A_B_CD_risk.png
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(forcats)
})

# ---------- Locate /data and /images ----------
find_data_dir <- function() {
  for (p in c("data", "../data", "./")) {
    if (dir.exists(p) && any(tolower(list.files(p)) %in% c("index_base.rds","index_base.csv")))
      return(normalizePath(p, winslash = "/"))
  }
  stop("Could not find /data with index_base.(rds|csv). WD: ", getwd())
}
DATA_DIR <- find_data_dir()
IMG_ROOT <- normalizePath(file.path(DATA_DIR, "..", "images"), winslash = "/", mustWork = FALSE)
if (!dir.exists(IMG_ROOT)) dir.create(IMG_ROOT, recursive = TRUE, showWarnings = FALSE)
IMG_DIR  <- normalizePath(file.path(IMG_ROOT, "sub-experiment"), winslash = "/", mustWork = FALSE)
if (!dir.exists(IMG_DIR)) dir.create(IMG_DIR, recursive = TRUE, showWarnings = FALSE)

RDS_PATH <- file.path(DATA_DIR, "index_base.rds")
CSV_PATH <- file.path(DATA_DIR, "index_base.csv")

# ---------- Load base ----------
if (file.exists(RDS_PATH)) {
  index_base <- readRDS(RDS_PATH)
} else if (file.exists(CSV_PATH)) {
  index_base <- readr::read_csv(CSV_PATH, show_col_types = FALSE)
} else {
  stop("index_base.rds or index_base.csv not found in: ", DATA_DIR)
}

# ---------- Prep ----------
req <- c("ClassABCD", "MACE_días")
miss <- setdiff(req, names(index_base))
if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))

df <- index_base %>%
  mutate(
    ClassABCD = factor(ClassABCD, levels = c("A","B","C","D")),
    MACE_días = suppressWarnings(as.numeric(MACE_días))
  ) %>%
  filter(!is.na(ClassABCD), !is.na(MACE_días)) %>%
  mutate(
    time  = pmin(MACE_días, 30, na.rm = TRUE),
    event = ifelse(MACE_días < 30, 1L, 0L)
  )

# Palettes
pal_AB_named  <- c(A="#1f77b4", B="#2ca02c")
pal_AB_unnamed <- unname(pal_AB_named)    # for unnamed fallback
pal_ABCDm_named <- c(A="#1f77b4", B="#2ca02c", CD="#ff7f0e")

# -------- Helpers --------
logrank_pval <- function(dat, group_var) {
  g <- factor(dat[[group_var]])
  g <- droplevels(g)
  if (length(levels(g)) < 2) return(NA_real_)
  if (sum(dat$event, na.rm = TRUE) == 0)  return(NA_real_)
  res <- survdiff(Surv(time, event) ~ g, data = transform(dat, g = g))
  1 - pchisq(res$chisq, df = length(res$n) - 1)
}

resolve_palette <- function(palette, levs) {
  # If palette is named, align by names; else recycle in order
  if (!is.null(names(palette))) {
    return(unname(palette[levs]))
  } else {
    if (length(palette) < length(levs)) palette <- rep(palette, length.out = length(levs))
    return(unname(palette[seq_along(levs)]))
  }
}

save_km <- function(sf, dat, group_var, title, xlab, ylab, legend_title, palette, file_stub) {
  # Ensure group is a factor with dropped unused levels
  dat[[group_var]] <- droplevels(factor(dat[[group_var]]))
  levs <- levels(dat[[group_var]])
  pal  <- resolve_palette(palette, levs)
  
  # Compute p-value if valid
  pval_ok <- !is.na(logrank_pval(dat, group_var))
  
  plt <- ggsurvplot(
    sf, data = dat,
    conf.int = FALSE,
    pval = pval_ok,
    risk.table = TRUE,
    risk.table.height = 0.22,
    ggtheme = theme_bw(base_size = 12),  # white background
    palette = pal,
    title = title,
    xlab = xlab, ylab = ylab,
    legend.title = legend_title
  )
  if (!pval_ok) {
    plt$plot <- plt$plot + ggtitle(title, subtitle = "P-value not shown: requires ≥2 strata with events.")
  }
  
  ggsave(filename = file.path(IMG_DIR, paste0(file_stub, ".png")),
         plot = plt$plot, width = 9, height = 6, dpi = 300, bg = "white")
  ggsave(filename = file.path(IMG_DIR, paste0(file_stub, "_risk.png")),
         plot = plt$table + theme_bw(base_size = 10), width = 9, height = 3, dpi = 300, bg = "white")
}

# ---------- Experiment 1: A vs B only ----------
df_AB <- df %>%
  filter(ClassABCD %in% c("A","B")) %>%
  mutate(ClassABCD = droplevels(factor(ClassABCD, levels = c("A","B"))))

if (nrow(df_AB) > 0) {
  sf_AB <- survfit(Surv(time, event) ~ ClassABCD, data = df_AB)
  # palette can be named or unnamed; both work now
  save_km(
    sf = sf_AB, dat = df_AB, group_var = "ClassABCD",
    title = "KM: MACE (A vs B) — 30-day",
    xlab = "Days", ylab = "Survival probability",
    legend_title = "Class",
    palette = pal_AB_named,
    file_stub = "km_MACE_A_vs_B"
  )
} else {
  message("No data for A vs B.")
}

# ---------- Experiment 2: A, B, and CD (C+D merged) ----------
df_ABCDm <- df %>%
  mutate(ClassABCDm = fct_collapse(ClassABCD, CD = c("C","D"))) %>%
  filter(ClassABCDm %in% c("A","B","CD")) %>%
  mutate(ClassABCDm = droplevels(factor(ClassABCDm, levels = c("A","B","CD"))))

if (nrow(df_ABCDm) > 0) {
  sf_ABCDm <- survfit(Surv(time, event) ~ ClassABCDm, data = df_ABCDm)
  save_km(
    sf = sf_ABCDm, dat = df_ABCDm, group_var = "ClassABCDm",
    title = "KM: MACE (A, B, CD) — 30-day",
    xlab = "Days", ylab = "Survival probability",
    legend_title = "Class",
    palette = pal_ABCDm_named,
    file_stub = "km_MACE_A_B_CD"
  )
} else {
  message("No data for A, B, CD merged.")
}

cat("\n✅ Sub-experiments complete.\nImages written to:\n  ", IMG_DIR, "\n\n")

