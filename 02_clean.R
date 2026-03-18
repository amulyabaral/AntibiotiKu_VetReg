
# ============================================================================
# VetReg — Data cleaning & preparation
# ============================================================================
# This script takes the cattle antibiotic VetReg extract (`vetreg`) produced by the loading
# script 01_load.R and applies sequential cleaning steps:
#
#   vetreg_1 — Date parsing, DOT (days on treatment) calculation
#   vetreg_2 — Animal-count imputation, treatment-unit classification &
#              reference-data joins (SPC product info, diagnosis codes)
#   vetreg_3 — Recode group → individual where ear-tags match animal count
#   vetreg_4 — Unit harmonisation, dose calculation & plausibility flagging
#
# Each step is kept as a separate object so intermediate results can be
# inspected during development.
# ============================================================================


# ============================================================================
# Step 1 — Parse dates & compute days-on-treatment (DOT)
# ============================================================================
# • Pad `varenummer` to a fixed 6-character width (leading zeros).
# • Convert date columns from character to Date.
# • Correct obvious data-entry errors in `planavsluttbehandling`:
#     - Year 9999 → replaced with 2023 (confirmed by manual inspection).
# • Derive `dot` = planned treatment duration in days (inclusive).
# • Coerce `tilbakeholdelsestid` (withdrawal time) to numeric.
# ============================================================================

vetreg_1 <- vetreg |>
  mutate(
    varenummer            = str_pad(varenummer, width = 6, side = "left", pad = "0"),
    registrertdato        = as_date(registrertdato),
    utlevertdato          = as_date(utlevertdato),
    year                  = year(utlevertdato),
    planavsluttbehandling = as_date(planavsluttbehandling),
    
    # --- Fix erroneous year = 9999 (confirmed via manual checking) ----------
    planavsluttbehandling = if_else(
      year(planavsluttbehandling) == 9999,
      make_date(2023, month(planavsluttbehandling), day(planavsluttbehandling)),
      planavsluttbehandling
    ),
    # --- Fix DOT > 365: re-anchor year to utlevertdato ---------------------
    # There are 2 records where they have reported 2021 when it should be 2020
    planavsluttbehandling = if_else(
      as.integer(planavsluttbehandling - utlevertdato) + 1 > 365,
      make_date(year(utlevertdato), month(planavsluttbehandling), day(planavsluttbehandling)),
      planavsluttbehandling
    ),
    # --- Derived columns ----------------------------------------------------
    dot                   = as.integer(planavsluttbehandling - utlevertdato) + 1,
    tilbakeholdelsestid   = as.numeric(tilbakeholdelsestid)
  )


# ============================================================================
# Step 2 — Impute animal counts, classify treatment unit & join reference data
# ============================================================================
# Goal: every row should have a usable `antalldyr` >= 1 and a classification
# of the treatment as "Individual", "Group", or "Unknown".
#
# Logic:
#   • If `merke` (ear-tag ID) is present but antalldyr == 0 → set to 1
#     (the animal is identified, so at least one animal was treated).
#   • If `merke` is missing and antalldyr == 0 → impute 1 (assumption;
#     discussed as a limitation).
#   • `treatment_unit` is assigned BEFORE the second imputation so that the
#     original zero-count rows are labelled "Unknown" rather than "Individual".
#
# After imputation, reference data is joined:
#   • `ref_all`             — SPC packaging unit & strength from Digivet.
#   • `diag_codes_english`  — Maps Norwegian diagnosis codes to English
#                             disease categories.
#   • Fix encoding artefact: "sprÃ¸yte" → "sprøyte".
#
# NOTE: When `antalldyr` > 1 AND unique ear-tags in the same report match
# `antalldyr`, the animals were individually identified.  These are recoded
# to antalldyr = 1 per merke in Step 3 (below). This is not done in Step 2
# because we need to make descriptive data table with this first.
# ============================================================================

vetreg_2 <- vetreg_1 |>
  mutate(
    antalldyr_old    = antalldyr,

    # --- Describe what happened to the count --------------------------------
    antalldyr_status = case_when(
      !is.na(merke) & antalldyr == 0 ~ "corrected to 1",
      is.na(merke) & antalldyr == 0 ~ "not given, imputed 1",
      TRUE                            ~ "used as is"
    ),
    merke_status     = if_else(!is.na(merke), "has_merke", "no_merke"),

    # --- First correction: merke present, count = 0 → 1 --------------------
    antalldyr        = if_else(!is.na(merke) & antalldyr == 0, "1", antalldyr),

    # --- Classify treatment unit -----------
    treatment_unit   = case_when(
      antalldyr >  1 ~ "Group",
      antalldyr == 1 ~ "Individual",
      antalldyr == 0 ~ "Unknown",
      TRUE           ~ NA_character_
    ),

    # --- Second correction: no merke, count still 0 → impute 1 -------------
    antalldyr        = if_else(is.na(merke) & antalldyr == 0, "1", antalldyr),

    # --- Final numeric coercion ---------------------------------------------
    antalldyr        = as.numeric(antalldyr)
  )

# --- Join reference data (SPC product info & diagnosis codes) ---------------

vetreg_2 <- vetreg_2 |>
  mutate(
    varenummer = str_pad(as.character(varenummer), width = 6, side = "left", pad = "0")
  ) |>

  # --- Join SPC packaging & strength info -----------------------------------
  left_join(ref_all, by = "varenummer") |>
  
    # --- Join English diagnosis categories ------------------------------------
  left_join(diag_codes_english, by = c("diagnose" = "diagnosis")) |>
  
    # --- Fix UTF-8 encoding artefact in packaging unit ------------------------
  mutate(
    lmp_enhet_pakning_v = if_else(
      lmp_enhet_pakning_v == "sprÃ¸yte", "sprøyte", lmp_enhet_pakning_v
    )
  )

# ============================================================================
# Step 3 — Recode group treatments where each animal is individually tagged
# ============================================================================
# Within a reportid, if the sum of antalldyr across rows (where antalldyr > 1)
# equals the number of distinct ear-tags (merke), then each row actually
# describes a single individually identified animal → set antalldyr to 1.
# e.g. a report with two rows (antalldyr = 5 and antalldyr = 2) and 7 distinct
# ear-tags means all 7 animals are tagged individually.
# Without this step, downstream dose calculations would incorrectly divide
# levert_mengde across multiple animals per row.
# ============================================================================
vetreg_3 <- vetreg_2 |>
  mutate(antalldyr_num = parse_number(as.character(antalldyr))) |>
  group_by(reportid) |>
  mutate(
    n_merke     = n_distinct(merke[antalldyr_num > 1], na.rm = TRUE),
    sum_dyr     = sum(antalldyr_num[antalldyr_num > 1], na.rm = TRUE),
    recode_flag = treatment_unit == "Group" & antalldyr_num > 1 & n_merke == sum_dyr
  ) |>
  ungroup() |>
  mutate(antalldyr = if_else(recode_flag, 1, antalldyr_num)) |>
  select(-antalldyr_num, -n_merke, -sum_dyr, -recode_flag)

# ============================================================================
# Step 4 — Unit harmonisation, dose calculation & plausibility flagging
# ============================================================================
# Only veterinarian-reported rows carry dosing information; pharmacy rows are
# marked "Not Applicable".
#
# Sub-steps:
#   a) Prepare VMP dose limits (min / max per product & packaging unit).
#   b) Join dose limits + unit-mismatch fix rules onto the data.
#   c) Coerce amount columns to numeric.
#   d) Handle the "dose → stk" reclassification for syringe products
#      (if reported unit is "dose", packaging is syringe, and amount is
#       plausible relative to DOT × animal count).
#   e) Determine unit-match category (`is_match`):
#        1 = exact match
#        2 = syringe reported as "stk"/"spr"
#        3 = bulk (g/ml) reported as "stk"/"spr" → multiply by pack strength
#        0 = no match (apply fix rule if available)
#   f) Calculate `calculated_dose` using the appropriate conversion.
#   g) Flag each row's dose plausibility:
#        "OK"                          — within [min_dose, max_dose × antalldyr]
#        "Flagged - Low"               — below minimum
#        "Flagged - High"              — above maximum (scaled by animal count)
#        "Not Checked - No Dose Info"  — reference limits unavailable
#        "Not Checked - Calculation Error" — conversion returned NA
#        "Not Applicable"              — pharmacy dispensing (no dosing data)
#   h) For rows that pass the check (dose_flag == "OK") with a unit mismatch,
#      overwrite `levert_mengde` and `enhet_mengde` with the corrected values.
# ============================================================================

# --- Prepare VMP dose-limit reference ----------------------------------

vmp_limits_clean <- vmp_limits |>
  mutate(
    varenummer = str_pad(as.character(varenummer), width = 6, side = "left", pad = "0"),
    min_dose   = readr::parse_number(as.character(`min dose`)),
    max_dose   = `Max dose*2`
  ) |>
  select(varenummer, lmp_enhet_pakning_v, min_dose, max_dose) |>
  distinct()

# --- Join, calculate, flag & correct ---------------------------------

vetreg_3 <- vetreg_3 |>
  
  # --- Joins ----------------------------------------------------------------
  left_join(vmp_limits_clean, by = c("varenummer", "lmp_enhet_pakning_v")) |>
  left_join(ref_except,       by = c("utleveringstype", "varenummer", "enhet_mengde")) |>
  
  # --- Coerce amounts to numeric --------------------------------------------
mutate(
  levert_mengde = readr::parse_number(as.character(levert_mengde)),
  fix_value     = readr::parse_number(as.character(fix_value)),
  lmp_mengde    = readr::parse_number(as.character(lmp_mengde)),
  lmp_antall    = readr::parse_number(as.character(lmp_antall)),
  antall_pakninger    = readr::parse_number(as.character(antall_pakninger))
  
)

vetreg_4 <- vetreg_3 |>
  mutate(
  # --- Determine unit-match category -----------------------------------
    is_match = case_when(
      enhet_mengde == lmp_enhet_pakning_v                              ~ 1,
      lmp_enhet_pakning_v == "sprøyte" & enhet_mengde %in% c("stk", "spr") ~ 2,
      lmp_enhet_pakning_v %in% c("g", "ml") & enhet_mengde %in% c("stk", "spr") ~ 3,
      TRUE                                                             ~ 0
    )) |>
  
  mutate(
    # --- Reclassify "dose" → "stk" for syringe products -----------------
    adjust_doser = utleveringstype == veterinarian &
      grepl("dose", enhet_mengde, ignore.case = TRUE) &
      lmp_enhet_pakning_v == "sprøyte" &
      levert_mengde <= 4 * dot * antalldyr,
    
    enhet_mengde = if_else(adjust_doser, "stk", enhet_mengde)) |>
  mutate(
    # --- Calculate dose (veterinarian rows only) -------------------------
    calculated_dose = if_else(
      utleveringstype == veterinarian,
      case_when(
        is_match == 0 & fix_method == "multiply_by" ~ levert_mengde * fix_value,
        is_match == 0 & fix_method == "divide_by"   ~ levert_mengde / fix_value,
        is_match %in% c(1, 2)                       ~ levert_mengde,
        is_match == 3                                ~ levert_mengde * lmp_mengde,
        TRUE                                         ~ NA_real_
      ),
      NA_real_
    ),
    
    # --- Flag dose plausibility ------------------------------------------
    no_dose_info = if_else(
      utleveringstype == veterinarian,
      is.na(min_dose) | is.na(max_dose),
      NA
    ),
    
    dose_flag = if_else(
      utleveringstype == veterinarian,
      case_when(
        no_dose_info == TRUE                                    ~ "Not Checked - No Dose Info",
        !is.na(calculated_dose) & calculated_dose < min_dose    ~ "Flagged - Low",
        !is.na(calculated_dose) & calculated_dose > (max_dose * antalldyr) ~ "Flagged - High",
        is.na(calculated_dose)                                  ~ "Not Checked - Calculation Error",
        TRUE                                                    ~ "OK"
      ),
      "Not Applicable"
    ),
    
    # --- Overwrite amount & unit for corrected rows ----------------------
    levert_mengde = if_else(
      utleveringstype == veterinarian & dose_flag == "OK" & is_match == 0,
      calculated_dose,
      levert_mengde
    ),
    
    enhet_mengde = if_else(
      utleveringstype == veterinarian & dose_flag == "OK" & is_match == 0,
      lmp_enhet_pakning_v,
      enhet_mengde
    )
  )

#### from pharmacies, 'g' enhet = they just mean sprøyte (manual verificaiton) 
## same with g spr they mean same unit 
##same with mlhgl 

vetreg_4 <- vetreg_4 |>
  mutate(
    calculated_dose_pharmacy = if_else(
      utleveringstype == pharmacy,
      antall_pakninger * lmp_antall * lmp_mengde,
      NA_real_
    ),
    
    no_dose_info_pharmacy = if_else(
      utleveringstype == pharmacy,
      is.na(min_dose) | is.na(max_dose),
      NA
    ),
    
    dose_flag_pharmacy = if_else(
      utleveringstype == pharmacy,
      case_when(
        no_dose_info_pharmacy == TRUE                                              ~ "Not Checked - No Dose Info",
        !is.na(calculated_dose_pharmacy) & calculated_dose_pharmacy < min_dose     ~ "Flagged - Low",
        !is.na(calculated_dose_pharmacy) & calculated_dose_pharmacy > (max_dose * antalldyr) ~ "Flagged - High",
        is.na(calculated_dose_pharmacy)                                            ~ "Not Checked - Calculation Error",
        TRUE                                                                       ~ "OK"
      ),
      "Not Applicable"
    )
  )
