# 01_build_analytic_dataset.R
rm(list = ls())

library(dplyr)
library(tidyr)
library(readr)

# =========================
# Paths
# =========================
raw_data_dir <- "data/raw"
derived_data_dir <- "data/derived"

if (!dir.exists(derived_data_dir)) dir.create(derived_data_dir, recursive = TRUE)

# =========================
# Read raw data
# =========================
data_participant <- read.csv(
  file = file.path(raw_data_dir, "data_participant.csv"),
  header = TRUE
)

data_death <- read.csv(
  file = file.path(raw_data_dir, "data_death.csv"),
  header = TRUE
)

data_death_cause <- read.csv(
  file = file.path(raw_data_dir, "data_death_cause.csv"),
  header = TRUE
)

# =========================
# Rename fields
# =========================
nm_map <- c(
  "eid"       = "eid",
  "p21022"    = "age_at_recruitment",
  "p31"       = "sex",
  "p53_i0"    = "date_at_recruitment",
  "p21001_i0" = "bmi",
  "p40023"    = "death_record",
  "p191"      = "date_lost_to_followup",
)

old_names <- names(data_participant)
names(data_participant) <- ifelse(
  old_names %in% names(nm_map),
  nm_map[old_names],
  old_names
)

# =========================
# Keep primary cause-of-death record if needed for bookkeeping
# (not used in the final outcome definition, but harmless to keep)
# =========================
death_cause_primary <- data_death_cause %>%
  mutate(cause_rank = readr::parse_number(sub(".*-(\\d+)$", "\\1", dnx_death_cause_id))) %>%
  arrange(eid, cause_rank) %>%
  group_by(eid) %>%
  slice(1) %>%
  ungroup() %>%
  select(eid, cause_icd10)

# =========================
# Merge core death information
# =========================
df_primary <- data_participant %>%
  left_join(data_death %>% select(eid, date_of_death), by = "eid") %>%
  left_join(death_cause_primary, by = "eid")

death_registry_end <- as.Date("2023-12-31")

# =========================
# Build full analytic dataset
# Outcome: all-cause mortality after baseline
# Censoring: lost to follow-up or registry end
# =========================
analytic_full <- df_primary %>%
  rename(
    id            = eid,
    age           = age_at_recruitment,
    baseline_date = date_at_recruitment,
    date_ltfu     = date_lost_to_followup
  ) %>%
  mutate(
    baseline_date = as.Date(baseline_date),
    date_ltfu     = as.Date(date_ltfu),
    date_of_death = as.Date(date_of_death)
  ) %>%
  filter(
    !is.na(id),
    !is.na(bmi),
    bmi > 15 & bmi < 40
  ) %>%
  mutate(
    death_record = as.integer(replace_na(death_record, 0) > 0),
    
    death_after = if_else(
      !is.na(date_of_death) & date_of_death > baseline_date,
      date_of_death,
      as.Date(NA)
    ),
    
    event_date = death_after,
    
    censor_date = if_else(
      !is.na(date_ltfu) & date_ltfu < death_registry_end,
      date_ltfu,
      death_registry_end
    ),
    
    old_group = as.integer(age > 65),
    
    outcome = if_else(
      !is.na(event_date) & event_date <= censor_date,
      1L,
      0L
    ),
    
    obs_end = if_else(
      outcome == 1L,
      event_date,
      censor_date
    ),
    
    time = as.numeric(obs_end - baseline_date) / 365.25
  ) %>%
  filter(
    time >= 0,
    !(death_record == 1 & outcome == 0)
  ) %>%
  transmute(
    id,
    bmi    = as.numeric(bmi),
    age    = as.numeric(age),
    sex,
    old_group = as.integer(old_group),
    time   = as.numeric(time),
    d      = as.integer(outcome)
  )

# =========================
# Final dataset for analysis
# =========================
SCALE <- 1

df_for_analysis <- analytic_full %>%
  transmute(
    bmi    = bmi,
    BMIc   = (bmi - 25) / SCALE,
    gender = if_else(sex == "Male", 1L, 0L),
    old    = if_else(old_group %in% c(0, "0"), 0L, 1L),
    time   = time,
    d      = d
  ) %>%
  filter(
    is.finite(bmi),
    is.finite(BMIc),
    is.finite(time),
    !is.na(d)
  )

write.csv(
  df_for_analysis,
  file.path(derived_data_dir, "df_for_analysis.csv"),
  row.names = FALSE
)

cat("Saved df_for_analysis with", nrow(df_for_analysis), "rows\n")