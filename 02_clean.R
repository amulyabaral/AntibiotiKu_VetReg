# hello

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
# inspected during development.But they are not necessary and can be deleted..
# ============================================================================
vetreg <- vetreg |>
    mutate(
      varenavn = case_when(
        varenummer == "001994" ~ "Carepen vet",
        varenummer == "011015" ~ "Streptocillin Forte vet",
        varenummer == "361411" ~ "Synulox LC Plus vet",
        .default = varenavn
      )
    )

# ============================================================================
# Parse dates & compute days of treatment / days-on-treatment (DOT)
# ============================================================================
# Pad `varenummer` to a fixed 6-character width (leading zeros).
# Convert date columns from character to Date.
# Correct obvious data-entry errors in `planavsluttbehandling`:
#     - Year 9999 → replaced with 2023 (confirmed by manual inspection).
# Derive `dot` = planned treatment duration in days (inclusive).
# Coerce `tilbakeholdelsestid` (withdrawal time) to numeric.
# ============================================================================

vetreg_1 <- vetreg |>
  #filter out drugs that are just straight up not allowed to be used in cattle so most likely dogs/cats
  filter(!str_detect(varenavn, "Abboticin|Antirobe|Convenia|Therios|Clindabactin|Dalacin|Denagard|Keflex|Noroclav|Ronaxan|Oxolinsyre vet")) |>
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
    # There are 2 records where they have reported 2021 when it should be 2020
    planavsluttbehandling = if_else(
      as.integer(planavsluttbehandling - utlevertdato) + 1 > 365,
      make_date(year(utlevertdato), month(planavsluttbehandling), day(planavsluttbehandling)),
      planavsluttbehandling
    ),
    # --- Derived columns ----------------------------------------------------
    dot                   = as.integer(planavsluttbehandling - utlevertdato) + 1,
    tilbakeholdelsestid   = as.numeric(tilbakeholdelsestid)
  ) |> 
  mutate(varenavn = if_else(varenummer == "019403" & varenavn == "Benzylpenicillin Vetcare", 
                                 "Carepen vet", 
                                 varenavn))

# ============================================================================
#  Recode animal counts, classify treatment unit & join reference data
# ============================================================================
# Every row should have a usable `antalldyr` >= 1 and a classification
# of the treatment as "Individual", "Group", or "Unknown".
#
# Logic:
#   If 'merke' (animal ID) is present but antalldyr == 0 then set it to 1
#     (the animal is identified, so at least one animal was treated).
#   If 'merke' is missing and antalldyr == 0 → impute 1 (this is an assumption;
#     and is discussed as a limitation).
#   'treatment_unit' is assigned BEFORE the second imputation so that the
#     original zero-count rows are labelled "Unknown" rather than "Individual" treatments.
#
# After imputation, reference data is joined:
#    'ref_all'             — SPC packaging unit & strength from Digivet.
#    'diag_codes_english'  — Maps Norwegian diagnosis codes to English
#                             disease categories.
#    Fix encoding artefact: "sprÃ¸yte" → "sprøyte".
#
# NOTE: When `antalldyr` > 1 AND number of unique ear-tags in the same report match
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
#  Recode group treatments 
# ============================================================================
# Within a reportid, if the value of antalldyr across rows (where antalldyr > 1)
# equals the number of distinct ear-tags (merke), then each row actually
# describes a single individually identified animal → set antalldyr to 1.
# e.g. 7 reports with two values of antalldyr (5 reports with antalldyr = 5 and antalldyr = 2) and 7 distinct
# ear-tags means all 7 animals are tagged individually.
# Without this step, downstream dose calculations would incorrectly divide
# levert_mengde across multiple animals per row.
# ============================================================================

vetreg_3 <- vetreg_2 |>
  mutate(
    # Within each group, check if the condition is met:
    # 1) the row has a merke and antalldyr > 1
    # 2) the sum of unique antalldyr values equals the number of distinct merke
    # (e.g. antalldyr is 7 and 7 unique merkes, or antalldyr has "5, 2" and there are 7 distinct merke)
    condition_met = !is.na(merke) & 
      antalldyr > 1 &
      sum(unique(antalldyr[!is.na(antalldyr)])) == n_distinct(merke, na.rm = TRUE),
    
    # If the condition is met, recode antalldyr to 1, otherwise keep as-is
    antalldyr = if_else(condition_met, 1L, antalldyr),
    
    .by = c(reportid, utlevertdato, registrertdato, utleveringstype, kilde, mottakers_produsentnr)
  ) |>
  select(-condition_met) # drop the helper column

# ============================================================================
# Unit harmonisation, dose calculation & plausibility flagging
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
#        0 = no match
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
    max_dose   = coalesce(`Max dose*2`, readr::parse_number(as.character(`Max dose`)))
  ) |>
  select(varenummer, lmp_enhet_pakning_v, min_dose, max_dose) |>
  distinct()

# --- Join ---------------------------------

vetreg_4 <- vetreg_3 |>
  
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


# =============================================================================
# Dose calculation, validation, and unit normalisation
#
# vetreg_5 works through three steps:
#
# 1. UNIT MATCHING (is_match)
#    if the dispensed unit (enhet_mengde) matches the drug's
#    packaging unit (lmp_enhet_pakning_v):
#      1 = exact match
#      2 = close match (syringe product dispensed as "stk" or "spr")
#      3 = partial match (g/ml product dispensed as stk/spr so it requires
#          multiplying by pack size to get the right quantity)
#      0 = no match — a fix must be applied
#
# 2. DOSE CALCULATION (calculated_dose)
#    Reclassifies "dose" → "stk" for syringe products where the quantity is
#    plausibly small (≤ 4 × daily dose × number of animals), then calculates
#    the actual dose based on dispenser type and match score:
#      Vet + is_match 0      → apply correction factor (multiply or divide)
#      Vet + is_match 1 or 2 → use dispensed amount directly
#      Vet + is_match 3      → dispensed amount × pack size
#      Pharmacy              → number of packs × pack count × pack size
#
# 3. PLAUSIBILITY FLAG (dose_flag)
#    Checks whether the calculated dose falls within expected min/max bounds
#    (max scaled by number of animals). Possible values:
#      "OK" | "Flagged - Low" | "Flagged - High" | "No Dose Info" |
#      "Calculation Error"
#
# =============================================================================



vetreg_5 <- vetreg_4 |>
  mutate(
    # --- Unit-match category -------------------------
    is_match = case_when(
      enhet_mengde == lmp_enhet_pakning_v                                       ~ 1,
      lmp_enhet_pakning_v == "sprøyte" & enhet_mengde %in% c("stk", "spr")      ~ 2,
      lmp_enhet_pakning_v %in% c("g", "ml") & enhet_mengde %in% c("stk", "spr") ~ 3,
      TRUE                                                                      ~ 0
    ),
    
    # after this is done, we can now start correcting
    # we do a small reclassification: if a vet gave something in "doses" 
    # for a syringe product, and the quantity is plausibly small 
    # (≤ 4× days of treatment × number of animals), 
    # it relabels the unit as "stk" 
    
    # --- Reclassify "dose" or its variations → "stk" for syringe products -----------------
    adjust_doser = utleveringstype == veterinarian &
      grepl("dose", enhet_mengde, ignore.case = TRUE) &
      lmp_enhet_pakning_v == "sprøyte" &
      levert_mengde <= 4 * dot * antalldyr,
    enhet_mengde = if_else(adjust_doser, "stk", enhet_mengde),
    
    #With units sorted, we now calculate the actual dose depending on who dispensed it
    # --- Calculate dose --------------------------------------------------
    calculated_dose = case_when(
      utleveringstype == veterinarian & is_match == 0 & fix_method == "multiply_by" ~ levert_mengde * fix_value,
      utleveringstype == veterinarian & is_match == 0 & fix_method == "divide_by"   ~ levert_mengde / fix_value,
      utleveringstype == veterinarian & is_match %in% c(1, 2)                       ~ levert_mengde,
      utleveringstype == veterinarian & is_match == 3                               ~ levert_mengde * lmp_mengde,
      utleveringstype == pharmacy                                                   ~ antall_pakninger * lmp_antall * lmp_mengde,
      TRUE                                                                          ~ NA_real_
    ),
    # --- Flag dose plausibility ------------------------------------------
    no_dose_info = if_else(
      utleveringstype %in% c(veterinarian, pharmacy),
      is.na(min_dose) | is.na(max_dose),
      NA
    ),
    dose_flag = case_when(
      no_dose_info == TRUE                                      ~ "No Dose Info",
      !is.na(calculated_dose) & calculated_dose < min_dose      ~ "Flagged - Low",
      !is.na(calculated_dose) & calculated_dose > (max_dose * antalldyr) ~ "Flagged - High",
      is.na(calculated_dose)                                    ~ "Calculation Error",
      TRUE                                                      ~ "OK"
    )
  )

# =============================================================================

# Unit normalisation
#    For vet rows where dose was within limits but units didn't originally match
#    aka (is_match == 0), writes the fixed values (both amount and unit based on the correct unit) 
#   back into the source columns

# =============================================================================

vetreg_6 <- vetreg_5 |>
  mutate(
    # --- Overwrite amount & unit for vet rows that are within dose, 
    #if the units do not match  ------------------
    levert_mengde = if_else(
      utleveringstype == veterinarian & dose_flag == "OK" & is_match == 0,
      calculated_dose, levert_mengde
    ),
    enhet_mengde = if_else(
      dose_flag == "OK" & is_match == 0,
      lmp_enhet_pakning_v, enhet_mengde
    )
  )

##### Now we impute!

# varenummers with at least one OK record
has_ok <- vetreg_6 |>
  filter(dose_flag == "OK") |>
  distinct(varenummer)

no_ok <- vetreg_6 |>
  summarise(n = n(), .by = c(varenummer, varenavn)) |>
  anti_join(has_ok, by = "varenummer") |>
  arrange(desc(n))

nrow(no_ok)
no_ok

# Enrich no_ok with SPC info
no_ok_spc <- no_ok |>
  left_join(Varenr_Virkestoff_unique, by = c("varenummer", "varenavn"))

# Extract the actual records from vetreg_6
no_ok_records <- vetreg_6 |>
  semi_join(no_ok, by = "varenummer")

#### there are 11 records that do not have a single correct OK dose varenummer 
### all of htem are human medication except one, which is 36 kg of 424075 clamoxyl vet medisinpellet - those have therefore been removed
## therefore those have been removed

###bactrim is fixed manually because they most likely meant 10 instead of 1 

###manual fixes 

vetreg_7 <- vetreg_6 |>
  mutate(
    calculated_dose = if_else(
      utleveringstype == veterinarian & varenummer == "364639" & enhet_mengde == "ml" & levert_mengde == 1,
      10, calculated_dose
    ),
    levert_mengde = if_else(
      utleveringstype == veterinarian & varenummer == "364639" & enhet_mengde == "ml" & levert_mengde == 1,
      10, levert_mengde
    )) |>
  mutate(
    dose_flag = if_else(
      utleveringstype == veterinarian & varenummer == "364639" & enhet_mengde == "ml" & dose_flag !="OK",
      "OK", dose_flag
                        )
        )
  
########### test

# varenummers with at least one OK record
has_ok <- vetreg_7 |>
  filter(dose_flag == "OK") |>
  distinct(varenummer)

no_ok <- vetreg_7 |>
  summarise(n = n(), .by = c(varenummer, varenavn)) |>
  anti_join(has_ok, by = "varenummer") |>
  arrange(desc(n))

# Enrich no_ok with SPC info
no_ok_spc <- no_ok |>
  left_join(Varenr_Virkestoff_unique, by = c("varenummer", "varenavn"))

# Extract the actual records from vetreg_6
no_ok_records <- vetreg_7 |>
  semi_join(no_ok, by = "varenummer")
###########################################

## manually removing entries that have a varenummer that does not have a single 'OK' dosage use, and 
## are likely errors for some records for these varenummers 
# 170628 - Gentamicin B. Braun (human drug)
# 424075 - Clamoxyl vet powder (poultry, withdrawn)
# 037995 - Penomax (human drug, pivmecillinam)

# 222851 - Sulfametizol alternova (human drug)  ##last one 17950393 is for 3000 animals wrong entry 
manual_removal_reportid <- c("14267595", "16503674", "16973378", "17136395", "20547893")

vetreg_8 <- vetreg_7 |>
  filter(!reportid %in% manual_removal_reportid) |>
  filter(!(reportid == "17950393" & antalldyr > 3000))

#### calculating mean dose  ####


product_mean_dose <- vetreg_8 |>
  filter(dose_flag == "OK") |>
  summarise(
   mean_dose = mean(calculated_dose, na.rm = TRUE),
    lmp_enhet_pakning_v = first(lmp_enhet_pakning_v, na.rm = TRUE),
    .by = c(varenummer, varenavn)
  )

product_median_dose <- vetreg_8 |>
  filter(dose_flag == "OK") |>
  summarise(
    median_dose = median(calculated_dose/antalldyr, na.rm = TRUE),
    lmp_enhet_pakning_v = first(lmp_enhet_pakning_v, na.rm = TRUE),
    .by = c(varenummer, varenavn)
  )

product_mean_dose_all <- vetreg_8 |>
  summarise(
    mean_dose = mean(calculated_dose, na.rm = TRUE),
    lmp_enhet_pakning_v = first(lmp_enhet_pakning_v, na.rm = TRUE),
    .by = c(varenummer, varenavn)
  )

product_median_dose_all <- vetreg_8 |>
  summarise(
    median_dose = median(calculated_dose/antalldyr, na.rm = TRUE),
    lmp_enhet_pakning_v = first(lmp_enhet_pakning_v, na.rm = TRUE),
    .by = c(varenummer, varenavn)
  )



###imputing
vetreg_9 <- vetreg_8 |>
  left_join(
    product_mean_dose |> select(varenummer, mean_dose),
    by = c("varenummer"),
  ) |>
  mutate(
    levert_mengde = case_when(
      levert_mengde == 500 & varenummer == "287359" ~ 1,
      dose_flag != "OK" & utleveringstype == veterinarian ~ mean_dose,
      dose_flag != "OK" & utleveringstype == pharmacy ~ mean_dose, #this does not change the amounts, just for checking
      .default = levert_mengde
    ),
    enhet_mengde = case_when(
      dose_flag != "OK" & utleveringstype == veterinarian ~ lmp_enhet_pakning_v,
      .default = enhet_mengde
    ),
  #  antall_pakninger = case_when(
  #   dose_flag != "OK" & utleveringstype == pharmacy ~ mean_dose / (lmp_antall * lmp_mengde),
   #   .default = antall_pakninger )
    )

  