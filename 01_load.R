# ============================================================================
# VetReg Data Cleaning - 01.Load & prepare data
# ============================================================================

Sys.setlocale("LC_CTYPE", "en_US.UTF-8")

if (!require(pacman)) install.packages("pacman")
pacman::p_load(
  # core tidyverse (includes ggplot2, dplyr, stringr, lubridate, readr, tidyr, purrr, forcats)
  tidyverse, data.table, lubridate,
  
  # data I/O & formats
  rio, readxl, openxlsx, googlesheets4, DBI, duckdb, jsonlite,
  
  # string & ID utilities
  stringr, openssl, uuid, clipr,
  
  # cleaning & manipulation
  janitor,
  
  # plotting — core extensions
  ggplot2, scales, RColorBrewer, paletteer,
  ggh4x, ggrepel, ggridges, ggdist, ggalluvial, ggthemes,
  
  # plotting — layout & export
  gridExtra, patchwork, corrplot, magick, webshot2,
  
  # tables & reporting
  knitr, gt, networkD3,
  
  # statistical methods & outlier detection
  outliers, EnvStats, robustbase, dbscan, isotree,
  
  # parallel processing
  parallel,
  
  # project structure
  here, glue, rstudioapi
)

# Set project root to the folder containing this script
project_root <- tryCatch(
  dirname(rstudioapi::getSourceEditorContext()$path),
  error = function(e) here::here()
)
if (is.null(project_root) || project_root == "") project_root <- here::here()
setwd(project_root)


# ============================================================================
# Constants & file paths
# ============================================================================

pharmacy     <- "Utlevering til dyrehold fra apotek m.m."
veterinarian <- "Melding om dyrehelsepersonells bruk av legemidler"
# ============================================================================
# ATC code prefixes — INCLUDED in analysis
# ============================================================================

atc_prefixes <- c(
  # --- A: Alimentary tract — intestinal antiinfectives ----------------------
  "QA07AA", "A07AA",     # Intestinal antibiotics (e.g. neomycin, colistin)
  "QA07AB", "A07AB",     # Intestinal sulfonamides (e.g. phthalylsulfathiazole)
  "QA07AX03", "A07AX03", # Nifuroxazide (nitrofuran intestinal antiinfective)
  "QA07AX04", "A07AX04", # Nifurzide (nitrofuran intestinal antiinfective)
  
  # --- G: Genito-urinary — gynecological/intrauterine antiinfectives --------
  "QG01AA", "G01AA",     # Gynecological antibiotics (e.g. oxytetracycline)
  "QG01AE", "G01AE",     # Gynecological sulfonamides (e.g. sulfatolamide)
  "QG01BA", "G01BA",     # Gynecological antibiotics + corticosteroids
  "QG01BE", "G01BE",     # Gynecological sulfonamides + corticosteroids
  "QG51AA",              # Intrauterine antibiotics
  "QG51AG",              # Intrauterine antiinfective/antiseptic combos
  
  # --- J: Systemic antiinfectives — antibacterials --------------------------
  "QJ01", "J01",         # Antibacterials for systemic use
  "QJ51",                # Antibacterials for intramammary use
  
  # --- P: Antiparasitic — antiprotozoals with AB activity -------------------
  "QP51AG"               # Sulfonamides used as antiprotozoals (e.g. sulfadimidine)
)

# ============================================================================
# ATC code prefixes — EXCLUDED from analysis
# (dermatologicals, ophthalmics, otics — topical/sensory routes)
# ============================================================================

atc_ignored <- c(
  # --- D: Dermatologicals ---------------------------------------------------
  "D06AA",               # Antibiotics for dermatological use (e.g. tetracyclines)
  "D06AX",               # Other antibiotics for topical/dermatological use
  "D06BA",               # Sulfonamides for dermatological use
  "QD06AA",              # Vet: dermatological antibiotics
  "QD06BA",              # Vet: dermatological sulfonamides
  
  # --- J: Antimycobacterials ------------------------------------------------
  "J04AB",               # Antibiotics for mycobacterial infections (e.g. rifampicin)
  "QJ02AA",              # Vet: systemic antimycotics (amphotericin B)
  
  # --- S: Sensory organs (eye & ear) ----------------------------------------
  "QS01AA",              # Vet: ophthalmic antibiotics
  "QS02CA",              # Vet: otic corticosteroid + antiinfective combos
  "S01AA",               # Ophthalmic antibiotics
  "S01AX",               # Other ophthalmic antiinfectives
  "S01CA",               # Ophthalmic corticosteroid + antiinfective combos
  "S02AA",               # Otic antiinfectives
  "S03CA"                # Ophthalmological/otological corticosteroid + antiinfective combos
)

atc_pattern     <- paste0("^(", paste(atc_prefixes, collapse = "|"), ")")
atc_ign_pattern <- paste0("^(", paste(atc_ignored,  collapse = "|"), ")")

# file paths from year range
vetreg_years      <- 2018:2023
vetreg_files_paths <- file.path(project_root, "raw_data", paste0("vetreg-", vetreg_years, ".csv"))
files_sql          <- paste0("'", vetreg_files_paths, "'", collapse = ", ")

# ============================================================================
# Load all VetReg data into DuckDB, then query from it
# ============================================================================
con <- dbConnect(duckdb())

parquet_path <- "vetreg_all.parquet"

# Build parquet only if it doesn't exist yet; delete the file to rebuild
if (!file.exists(parquet_path)) {
  message("Parquet not found \u2014 building from CSVs...")
  dbExecute(con, glue("
    COPY (
      SELECT
        reportid, reseptid, registrertdato, utlevertdato, utleveringstype,
        kilde, helsepersonell, tilfoertav, mottakerpostnr, mottakerpoststed,
        aktivitet, mottakers_produsentnr, dyrekategori, merke, antalldyr,
        beskrivelse, diagnose, varenummer, varenavn, atckode, levert_mengde,
        enhet_mengde, antall_pakninger, planavsluttbehandling,
        tilbakeholdelsestid, feilkoder,
        EXTRACT(YEAR FROM TRY_STRPTIME(utlevertdato, '%Y-%m-%d')) AS year
      FROM read_csv(
        [{files_sql}],
        delim = ';', header = true, all_varchar = true,
        filename = true, quote = '\"',
        strict_mode = false, ignore_errors = true
      )
      WHERE EXTRACT(YEAR FROM TRY_STRPTIME(utlevertdato, '%Y-%m-%d')) >= 2018
    ) TO '{parquet_path}' (FORMAT PARQUET)
  "))
  message("Parquet saved.")
}

# Read cattle + AB subset from parquet
vetreg <- dbGetQuery(con, glue("
  SELECT *
  FROM '{parquet_path}'
  WHERE dyrekategori ILIKE '%storfe%'
    AND regexp_matches(atckode, '{atc_pattern}')
")) |> setDT()

dbDisconnect(con)

# ============================================================================
# Reference data: SPC, dose limits, unit fixes
# ============================================================================

# -- SPC product info (from Digivet) ----------------------------------------
load(file.path(project_root, "data", "Varenr_Virkestoff_unique.rds"))

spc_base <- Varenr_Virkestoff_unique |>
  mutate(varenummer = str_pad(varenummer, width = 6, side = "left", pad = "0"),
         lmp_antall = as.numeric(lmp_antall),
         lmp_antall = if_else(is.na(lmp_antall) | lmp_antall == 0, 1, lmp_antall)
  )

spc_info <- spc_base |>
  bind_rows(
    spc_base |> filter(varenavn == "Doxylin") |> slice(1) |>
      mutate(varenummer = "002971", lmp_pakningsstr = "50",
             lmp_mengde = "50", datakilde = "Manual"),
    spc_base |> filter(varenavn == "Carepen vet", varenummer == "001994") |> slice(1) |>
      mutate(varenummer = "482630"),
    spc_base |> filter(varenavn == "Carepen vet", varenummer == "019403") |> slice(1) |>
      mutate(varenummer = "486198")
  )

# Subset used for joins later
ref_all <- spc_info |>
  select(varenummer, lmp_enhet_pakning_v, lmp_mengde, lmp_antall, legemiddelform_kort_dn) |>
  distinct()

# -- Unit-mismatch fix rules ------------------------------------------------
# Products where reported unit doesn't match SPC → apply conversion
ref_except <- read_csv(file.path(project_root, "data", "refs_ap1 - medications_without_matching_units.csv")) |>
  select(utleveringstype, varenummer, enhet_mengde, fix_method, fix_value, comment) |>
  bind_rows(
    expand_grid(
      utleveringstype = c(veterinarian, pharmacy),
      varenummer      = "128131",
      enhet_mengde    = c("enpac", "ml", NA_character_)
    ) |>
      mutate(fix_method = "multiply_by", fix_value = 1, comment = "meant stk")
  ) |>
  distinct()

# -- VMP dose limits --------------------------------------------------------
vmp_limits <- read_csv(file.path(project_root, "data", "VMP_limits.csv"))


# -- DDDvet reference values from EMA ---------------------------------
dddvet <- read_csv(file.path(project_root, "data", "dddvet_values_ema.csv"))

# -- Diagnosis codes (in English) ----------------------------------------------
diag_codes_english <- read_csv2(file.path(project_root, "data", "diag_codes_english.csv")) |>
  select(diagnosis, disease_category = main_category)

# ============================================================================
# Population data
# ============================================================================

mimiro_herd   <- read_csv(file.path(project_root, "data", "mimiro_herd.csv"))
animalia_herd <- read_csv2(file.path(project_root, "data", "animalia_herd.txt"))
cattle_pop    <- read_tsv(file.path(project_root, "data", "cattle_pop_ssb.tsv"))