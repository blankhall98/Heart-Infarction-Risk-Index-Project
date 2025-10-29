# scripts/km_analysis.R
# ============================================================
# Kaplan–Meier analyses for five endpoints, stratified by ClassABCD.
# Study duration: 30 days (value == 30 => censored, value < 30 => event).
#
# Inputs  (auto-located): data/index_base.rds or data/index_base.csv
# Outputs:
#   images/km_<endpoint>_by_class.png
#   images/km_<endpoint>_by_class_risk.png
#   data/km_overall_summary.csv
#   data/km_medians_by_class.csv
#
# Note: Most "file not found" errors come from working directory differences.
# This script auto-finds /data and /images relative to where you run it.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(broom)
  library(rlang)
  library(tidyr)
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
IMG_DIR  <- normalizePath(file.path(DATA_DIR, "..", "images"), winslash = "/", mustWork = FALSE)
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

# ---------- Basic checks & prep ----------
required_cols <- c("ClassABCD", "SET1", "TAPSE1_PSAP1",
                   "MACE_días","MUERTECV_días","REINFARTO_días","CHOQUE_días","SANGRADOMAYOR_días")
missing <- setdiff(required_cols, names(index_base))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

index_base <- index_base %>%
  mutate(
    ClassABCD = factor(ClassABCD, levels = c("A","B","C","D")),
    SET1 = suppressWarnings(as.numeric(SET1)),
    TAPSE1_PSAP1 = suppressWarnings(as.numeric(TAPSE1_PSAP1)),
    gender = if ("gender" %in% names(.)) as.factor(gender) else NULL,
    Age    = if ("Age" %in% names(.)) suppressWarnings(as.numeric(Age)) else NULL,
    bmi    = if ("bmi" %in% names(.)) suppressWarnings(as.numeric(bmi)) else NULL
  )

# Endpoints and labels
endpoints <- c("MACE_días","MUERTECV_días","REINFARTO_días","CHOQUE_días","SANGRADOMAYOR_días")
endpoint_labels <- c(
  MACE_días = "Any MACE",
  MUERTECV_días = "Cardiovascular Death",
  REINFARTO_días = "Reinfarction",
  CHOQUE_días = "Shock",
  SANGRADOMAYOR_días = "Major Bleeding"
)

# Base palette for A/B/C/D
base_class_cols <- c(A = "#1f77b4", B = "#2ca02c", C = "#ff7f0e", D = "#d62728")

# ---------- Helpers ----------
logrank_pval <- function(df_ep) {
  # Needs ≥2 strata AND at least one event overall
  if (length(levels(df_ep$ClassABCD)) < 2) return(NA_real_)
  if (sum(df_ep$event, na.rm = TRUE) == 0)  return(NA_real_)
  res <- survdiff(Surv(time, event) ~ ClassABCD, data = df_ep)
  1 - pchisq(res$chisq, df = length(res$n) - 1)
}

cox_tidy_or_null <- function(formula, data) {
  # Needs ≥2 strata and events in ≥2 strata for a stable fit
  if (length(levels(data$ClassABCD)) < 2) return(NULL)
  by_stratum_events <- with(data, tapply(event, ClassABCD, function(z) sum(z, na.rm = TRUE)))
  if (sum(by_stratum_events > 0, na.rm = TRUE) < 2) return(NULL)
  fit <- tryCatch(coxph(formula, data = data, ties = "efron", iter.max = 50),
                  warning = function(w) w, error = function(e) NULL)
  if (inherits(fit, "warning") || is.null(fit)) return(NULL)
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
}

median_by_class <- function(fit) {
  med <- tryCatch(surv_median(fit), error = function(e) NULL)
  if (is.null(med)) return(NULL)
  med %>% transmute(ClassABCD = gsub("^ClassABCD=", "", strata), median = median)
}

extract_hrs <- function(tt) {
  if (is.null(tt)) return(tibble(term = character(), estimate = numeric(), conf.low = numeric(), conf.high = numeric()))
  tt %>% dplyr::filter(grepl("^ClassABCD", term)) %>%
    transmute(term, HR = estimate, CI_low = conf.low, CI_high = conf.high)
}
to_wide <- function(df, prefix) {
  if (nrow(df) == 0) return(tibble())
  df %>%
    mutate(comp = sub("^ClassABCD", "", term)) %>%
    select(comp, HR, CI_low, CI_high) %>%
    tidyr::pivot_wider(
      names_from = comp,
      values_from = c(HR, CI_low, CI_high),
      names_glue = paste0(prefix, "{comp}_{.value}")
    )
}

# ---------- Main loop ----------
overall_summary <- list()
medians_all <- list()

for (ep in endpoints) {
  ep_sym <- rlang::sym(ep)
  
  df_ep <- index_base %>%
    filter(!is.na(!!ep_sym), !is.na(ClassABCD)) %>%
    mutate(
      time  = pmin(as.numeric(!!ep_sym), 30, na.rm = TRUE),
      event = ifelse(as.numeric(!!ep_sym) < 30, 1L, 0L)
    ) %>%
    droplevels()   # drop unused ClassABCD levels for this endpoint
  
  if (nrow(df_ep) == 0) next
  
  # Palette for present levels only (unnamed to avoid level-name mismatch warnings)
  levs <- levels(df_ep$ClassABCD)
  pal  <- unname(base_class_cols[levs])
  
  # Survfit
  sf <- survfit(Surv(time, event) ~ ClassABCD, data = df_ep)
  
  # Log-rank
  p_logrank <- logrank_pval(df_ep)
  show_pval <- !is.na(p_logrank)
  
  # KM plot (white background)
  plt <- ggsurvplot(
    sf, data = df_ep,
    conf.int = FALSE,
    pval = show_pval,
    risk.table = TRUE,
    risk.table.height = 0.22,
    ggtheme = theme_bw(base_size = 12),   # white bg
    palette = pal,
    title = paste0("KM: ", endpoint_labels[[ep]], " (30-day)"),
    xlab = "Days",
    ylab = "Survival probability",
    legend.title = "Class"
  )
  
  if (!show_pval) {
    plt$plot <- plt$plot + ggtitle(
      paste0("KM: ", endpoint_labels[[ep]], " (30-day)"),
      subtitle = "P-value not shown: requires ≥2 strata with events."
    )
  }
  
  # Save images (white bg)
  img_name <- paste0("km_", sub("_días$", "", ep), "_by_class.png")
  ggsave(filename = file.path(IMG_DIR, img_name), plot = plt$plot,
         width = 9, height = 6, dpi = 300, bg = "white")
  ggsave(filename = file.path(IMG_DIR, sub(".png$", "_risk.png", img_name)),
         plot = plt$table + theme_bw(base_size = 10),
         width = 9, height = 3, dpi = 300, bg = "white")
  
  # Summaries
  n_total  <- nrow(df_ep)
  n_events <- sum(df_ep$event, na.rm = TRUE)
  n_cens   <- n_total - n_events
  
  # Cox models
  cox_unadj <- cox_tidy_or_null(Surv(time, event) ~ ClassABCD, df_ep)
  
  # Adjusted covariates usable in this subset
  cand_covs <- intersect(c("gender","Age","bmi","SET1","TAPSE1_PSAP1"), names(df_ep))
  usable <- cand_covs[sapply(cand_covs, function(v) {
    vv <- df_ep[[v]]
    if (all(is.na(vv))) return(FALSE)
    if (is.numeric(vv)) return(sd(vv, na.rm = TRUE) > 0)
    if (is.factor(vv) || is.character(vv)) return(length(unique(na.omit(vv))) > 1)
    TRUE
  })]
  cox_adj <- NULL
  if (length(usable) > 0) {
    f_adj <- as.formula(paste("Surv(time, event) ~ ClassABCD +", paste(usable, collapse = " + ")))
    cox_adj <- cox_tidy_or_null(f_adj, df_ep)
  }
  
  # Medians by class
  med_df <- median_by_class(sf)
  if (!is.null(med_df)) {
    med_df$endpoint <- ep
    med_df$endpoint_label <- endpoint_labels[[ep]]
    medians_all[[ep]] <- med_df
  }
  
  # HRs wide
  wide_u <- to_wide(extract_hrs(cox_unadj), "U_")
  wide_a <- to_wide(extract_hrs(cox_adj),   "A_")
  
  base_row <- tibble(
    endpoint = ep,
    endpoint_label = endpoint_labels[[ep]],
    n = n_total,
    events = n_events,
    censored = n_cens,
    logrank_p = p_logrank
  )
  overall_summary[[ep]] <- bind_cols(base_row, wide_u, wide_a)
}

# ---------- Write summaries ----------
overall_summary_df <- overall_summary %>% bind_rows()
medians_df <- medians_all %>% bind_rows() %>%
  tidyr::pivot_wider(names_from = ClassABCD, values_from = median, names_prefix = "median_")

readr::write_csv(overall_summary_df, file.path(DATA_DIR, "km_overall_summary.csv"))
readr::write_csv(medians_df,        file.path(DATA_DIR, "km_medians_by_class.csv"))

cat("\n✅ KM analysis complete.",
    "\nImages saved in: ", IMG_DIR,
    "\nSummary CSVs:",
    "\n - ", file.path(DATA_DIR, "km_overall_summary.csv"),
    "\n - ", file.path(DATA_DIR, "km_medians_by_class.csv"), "\n\n")

