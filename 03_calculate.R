
# We remove NAs from the following columns and replace it with a blank. This is not necessary but
# helps avoid any problems with merging in the future.
vetreg_10 <- vetreg_9 %>%
  mutate(
    reportid = replace(reportid, is.na(reportid), "")
  ) %>%
  mutate(
    reseptid = replace(reseptid, is.na(reseptid), "")
  ) %>%
  mutate(
    merke = replace(merke, is.na(merke), "")
  )

### Merging VetReg with FEST----
# We finally merge VetReg with FEST, we use inner join to get all the data in the VetReg that shows up in FEST

#first we have to add some missing info

Varenr_Virkestoff_unique <- Varenr_Virkestoff_unique |>
  mutate(varenummer = str_pad(varenummer, width = 6, side = "left", pad = "0"),
         lmp_antall = as.numeric(lmp_antall),
         lmp_antall = if_else(is.na(lmp_antall) | lmp_antall == 0, 1, lmp_antall),
         across(where(is.character), \(x) str_replace_all(x, "sprÃ¸yte", "sprøyte"))
  ) 

Varenr_Virkestoff_unique <- Varenr_Virkestoff_unique |>
  bind_rows(
    Varenr_Virkestoff_unique |> filter(varenavn == "Doxylin") |> slice(1) |>
      mutate(varenummer = "002971", lmp_pakningsstr = "50",
             lmp_mengde = "50", datakilde = "Manual"),
    Varenr_Virkestoff_unique |> filter(varenavn == "Carepen vet", varenummer == "001994") |> slice(1) |>
      mutate(varenummer = "482630"),
    Varenr_Virkestoff_unique |> filter(varenavn == "Carepen vet", varenummer == "019403") |> slice(1) |>
      mutate(varenummer = "486198")
  ) 

#lets start like we have been before
vetreg_11 <- vetreg_10 |>
  select("reportid", "reseptid", "registrertdato","utlevertdato","utleveringstype","kilde","helsepersonell", "tilfoertav",
         "mottakerpostnr",  "mottakerpoststed", "aktivitet", "mottakers_produsentnr", "dyrekategori","merke",
         "antalldyr","beskrivelse", "diagnose", "varenummer","varenavn","atckode","levert_mengde","enhet_mengde",
         "antall_pakninger", "planavsluttbehandling", "tilbakeholdelsestid","feilkoder")


vetreg_11 <- vetreg_11 |>
  mutate(unique_row_id = UUIDgenerate(n = n()))

VetReg_merged <- vetreg_11 |>
  inner_join(Varenr_Virkestoff_unique, by = "varenummer", relationship = "many-to-many")
# Okay to overlook the warning message on unexpected many-to-many relationship (it is expected)

# The ATC codes are almost always consistent, but in case of any discrepancies, we trust the code in FEST
# So we remove atckode and varenavn column that we got from VetReg and other columns not needed from the final merged dataframe
# We also standardize the varenavn column.
VetReg_merged <- VetReg_merged |>
  select(-c(atckode,varenavn.x)) |>
  rename(varenavn = varenavn.y)

### Information----
# The VetReg data that has been standardized will be used here. The purpose of this script is to get the antibiotic
# specific records from VetReg and standardise it.

### Input data to work on----
# Read support registers:

# Antibiotics ATC codes
AB_atc <- import("data/ATC_codes_antimicrobials.csv") %>%
  clean_names()

# Antibiotics class
# Here we use encoding = "Latin-1" to read Norwegian letters
# Using header as FALSE because there are two columns with merged header here so we will assign header row later
# Also we only need the first 6 columns so we use [1:6]

AB_class <- import("data/ABvirkestoff tilegnet ABklasser.csv", encoding = "Latin-1", header = F)[,1:6]
AB_class <- AB_class %>%
  row_to_names(row_number = 1) %>%
  clean_names()
AB_class <- AB_class %>%
  rename(ab_klasse = a_bklasse,
         ab_rapport = a_brapport,
         ab_klasse_2 = abklasse_2,
         ab_klasse_engelsk = abklasse_engelsk)    #again, some manual column name cleaning is required
# Select unique combinations of (virkestoff_navn, ab_klasse, ab_rapport and remove the rows where virkestoff_navn is not given
AB_class <- AB_class %>%
  distinct(virkestoff_navn, ab_klasse, ab_rapport) %>%
  filter(virkestoff_navn !="")

# Antibiotics conversion factor
AB_conversion <- import("data/AB conversion factors.csv", encoding = "Latin-1") %>%
  clean_names()
# Only select the relevant columns and rename those columns to make the column names standard across all support registers
AB_conversion <- AB_conversion %>%
  select(virkestoff_navn = derivat_virkestoff_i_iu_norsk,
         active_moiety = aktivt_virkestoffstoff_engelsk,
         iu_to_mg = iu_to_mg_conversion_factors,
         product_to_active = conversion_factors_of_derivatives)
# The active moiety for Bacitracin is Bacitracin, we add this information:
AB_conversion <- AB_conversion %>%
  dplyr::mutate(active_moiety = replace(active_moiety, virkestoff_navn == "Bacitracin", "Bacitracin"))
# There are commas that need to be replaced with decimals in AB_conversion
AB_conversion <- AB_conversion %>%
  mutate(
    across(
      .cols = c(iu_to_mg, product_to_active),
      \(x) as.numeric(gsub(",",".",x))
    )
  )
# Some more cleaning on AB_conversion: remove rows where there is no data, use 1 instead of NA in column product_to_active
AB_conversion <- AB_conversion %>%
  filter(!is.na(iu_to_mg) & virkestoff_navn !="")
AB_conversion <- AB_conversion %>%
  mutate(product_to_active = case_when(is.na(product_to_active) ~ 1,
                                       TRUE ~ product_to_active)) %>%
  distinct()

# Antibiotics short names
form_short_name <- import("data/Legemiddelform_norsk_engelsk.csv", encoding = "Latin-1") %>%
  clean_names()
# Select unique combinations of legemiddelform_kort_dn, legemiddelform_vi and remove rows where legemiddelform_kort_dn is not given
form_short_name <- form_short_name %>%
  distinct(legemiddelform_kort_dn, legemiddelform_vi) %>%
  filter(legemiddelform_kort_dn !="")

VetReg <- VetReg_merged %>%
  mutate(
    across(                                                        
      .cols = c(planavsluttbehandling, registrertdato, utlevertdato),                      # the columns that need to be modified
      .fns = \(x) lubridate::ymd(x)                                                        # We use lubridate::ymd() to change character to Date class
    )
  )
# Select the columns where the delivery date is in the year required
#VetReg <- VetReg %>%
#  filter(format(utlevertdato, "%Y") == "2020") #change here

# Select antibiotic specific data. Select the rows where the ATC code in VetReg is in AB_atc support register
antibiotics <- VetReg %>%
  filter(substr(atc_kode,1,8) %in% AB_atc$atc_code | substr(atc_kode,1,7) %in% AB_atc$atc_code |
           substr(atc_kode,1,6) %in% AB_atc$atc_code | substr(atc_kode,1,5) %in% AB_atc$atc_code | 
           substr(atc_kode,1,4) %in% AB_atc$atc_code | substr(atc_kode,1,3) %in% AB_atc$atc_code)

# Join the antibiotic dataframe with form_short_name (by legemiddelform_kort_dn)
# and AB_class support register (by virkestoff_navn)
antibiotics <- left_join(antibiotics, form_short_name, by = "legemiddelform_kort_dn")
antibiotics <- left_join(antibiotics, AB_class, by = "virkestoff_navn")

# After manually checking against SPC, change the active substance name from Benzylpenicillin to Benzylpenicillinprokain for the following varenummers
antibiotics <- antibiotics %>%
  dplyr::mutate (
    virkestoff_navn = replace(
      virkestoff_navn, (varenummer == "001994" | varenummer == "019403" | varenummer == "027429" |
                          varenummer == "027442" | varenummer == "076460" | varenummer == "144907" |
                          varenummer == "402735" | varenummer == "454918" | varenummer == "511675" |
                          varenummer == "529597" | varenummer == "560289") 
      & virkestoff_navn == "Benzylpenicillin", "Benzylpenicillinprokain"
    )
  )

# After manually checking against SPC, change the active substance name from Benzylpenicillin to Benzylpenicillinbenzatin for the following varenummers
antibiotics <- antibiotics %>%
  dplyr::mutate (
    virkestoff_navn = replace(
      virkestoff_navn, (varenummer == "264994") 
      & virkestoff_navn == "Benzylpenicillin", "Benzylpenicillinbenzatin"
    )
  )

# After manually checking against SPC, change the active substance name from Benzylpenicillin to Penicillinbenethamin for the following varenummers
antibiotics <- antibiotics %>%
  dplyr::mutate (
    virkestoff_navn = replace(
      virkestoff_navn, (varenummer == "205876" | varenummer == "279687") 
      & virkestoff_navn == "Benzylpenicillin", "Penicillinbenethamin"
    )
  )

# Join the antibiotics dataframe with AB_conversion support register
# Virkestoff_navn is the name of the column for which the strength is given 
# This column is a mixture of derivative and active moiety
antibiotics <- left_join(antibiotics, AB_conversion[,c("virkestoff_navn","active_moiety")],  by = "virkestoff_navn")
antibiotics <- antibiotics %>%
  relocate(active_moiety, .after = virkestoff_navn)
antibiotics <- antibiotics %>%
  mutate(active_moiety = case_when(
    is.na(active_moiety) ~ virkestoff_navn,
    TRUE ~ active_moiety
  ))

antibiotics <- left_join(antibiotics, AB_conversion[,c("virkestoff_navn","iu_to_mg")],  by = "virkestoff_navn")
antibiotics <- left_join(antibiotics, AB_conversion[,c("virkestoff_navn","product_to_active")], by = "virkestoff_navn")
antibiotics <- antibiotics %>%
  mutate(
    product_to_active = replace(product_to_active, is.na(product_to_active), 1)
  )

### Formatting----

#removing all variables from the environment except antibiotics
#rm (list = setdiff(ls(),"antibiotics"))

# Change the delivery type where the vet has probably used the wrong form of service
antibiotics <- antibiotics %>%
  mutate(
    utleveringstype = replace(
      utleveringstype, utleveringstype == "Utlevering til dyrehold fra apotek m.m." &
        kilde == "Skjematjeneste" & 
        aktivitet != "Akvakultur fisk" &
        legemiddelform_vi != "Premiks", "Melding om dyrehelsepersonells bruk av legemidler"
    )
  )

#I like using replace_tag, this will be used in the next step
antibiotics <- antibiotics %>%
  mutate(
    replace_tag = NA
  )

# Correct the strength and strength denominator where package is in kg or L
# The following changes need to be made:
# When strength is in mg and strength denominator is in g and package in in kg, change strength to g and strength denominator to kg
# When strength is in mg and strength denominator is in ml and package in in L, change strength to g and strength denominator to L
# When strength denominator volume is not given and strength denominator unit is kg and package in in kg, change strength denominator volume to 1
# When strength denominator volume and unit are not given and package is in ml or not given and the form is injection, change strength denominatior volume to 1 and unit to ml
antibiotics <- antibiotics %>%
  mutate(
    replace_tag = 
      replace(replace_tag, styrke_u == "mg" & styrke_nevner_u == "g" & lmp_enhet_pakning_v == "kg", 1),
    replace_tag =
      replace(replace_tag, styrke_u == "mg" & styrke_nevner_u == "ml" & lmp_enhet_pakning_v == "L", 2),
    replace_tag =
      replace(replace_tag, is.na(styrke_nevner_v) & styrke_nevner_u == "kg" & lmp_enhet_pakning_v == "kg", 3),
    replace_tag =
      replace(replace_tag,is.na(styrke_nevner_v) & is.na(styrke_nevner_u) & lmp_enhet_pakning_v %in% c("","ml") & legemiddelform_vi == "Injeksjon", 4)
  ) 

antibiotics <- antibiotics %>%
  mutate(
    styrke_u = replace(styrke_u, replace_tag == 1, "g"),
    styrke_u = replace(styrke_u, replace_tag == 2, "g"),
    styrke_nevner_u = replace(styrke_nevner_u, replace_tag == 1, "kg"),
    styrke_nevner_u = replace(styrke_nevner_u, replace_tag == 2, "L"),
    styrke_nevner_v = replace(styrke_nevner_v, replace_tag == 3, 1),
    styrke_nevner_u = replace(styrke_nevner_u, replace_tag == 4, "ml"),
    styrke_nevner_v = replace(styrke_nevner_v, replace_tag == 4, 1)
  ) %>%
  mutate(
    replace_tag = NULL
  )


source(file = "format_antibiotic_data_for_digivet.R", encoding = "UTF-8")




#########################


#================================================================




ab_norwegian_to_english <- read_csv2("data/ABvirkestoff tilegnet ABklasser.csv", locale = locale(encoding = "latin1"))


antibiotics_report$activeingredient_clean <- trimws(antibiotics_report$ActiveIngredient)

antibiotics_report_1 <- antibiotics_report |>
  mutate(virkestoff_navn = str_trim(virkestoff_navn))
ab_class <- read_csv("data/ab_class.csv")

antibiotics_report_1 <- antibiotics_report_1 |>
  left_join(ab_class, by = "activeingredient_clean")
# "Others" category for any unmatched active ingredients

antibiotics_report_2 <- antibiotics_report_1 |>
  filter(include_in_amu == TRUE)


diag_codes_english <- read.csv2("data/diag_codes_english.csv")

antibiotics_report_3 <- antibiotics_report_2 %>%
  mutate(
    DateTreatment = as.Date(DateTreatment),
    Year = format(DateTreatment, "%Y"),
    Month = factor(format(DateTreatment, "%b"), 
                   levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")),
    ActiveIngredientCategory = factor(drug_class, 
                                      levels = unique(drug_class))
  ) %>%
  left_join(read.csv2("data/diag_codes_english.csv") %>%
              dplyr::select(diagnosis, main_category), 
            by = c("Diagnosis" = "diagnosis")) %>%
  mutate(
    main_category = str_replace_all(main_category, regex("conditions|disorders", ignore_case = TRUE), ""),
    AdministrationMethod = str_replace_all(AdministrationMethod, regex(" products", ignore_case = TRUE), "")
  ) 

antibiotics_report_4 <- antibiotics_report_3 %>%
  mutate(
    Route = case_when(
      grepl("Oral oppløsning|Tablett|Oralpasta|Oralt pulver|Oral oppløsning_pulver|Medisinpellet", legemiddelform_vi) ~ "Oral",
      grepl("Injeksjon|Infusjon|Injeksjon_pulver|Infusjon_pulver", legemiddelform_vi) ~ "Parenteral",
      grepl("Uteritorie", legemiddelform_vi) ~ "Intrauterine",
      legemiddelform_vi == "Intramammarie" ~ "Intramammary",
      legemiddelform_vi == "Øyepreparat" ~ "Ophthalmic",
      legemiddelform_vi == "Ørepreparat" ~ "Otic",
      grepl("Hudpreparat|Pasta", legemiddelform_vi) ~ "Topical",
      grepl("Tyggetablett|Medisinpellet|Oralpasta", legemiddelform_kort_dn) ~ "Oral",
      grepl("Injeksjonsvæske|Infusjonsvæske", legemiddelform_kort_dn) ~ "Parenteral",
      grepl("Intramammarie", legemiddelform_kort_dn) ~ "Intramammary",
      grepl("Øyedråper|Øyesalve", legemiddelform_kort_dn) ~ "Ophthalmic",
      grepl("Øredråper", legemiddelform_kort_dn) ~ "Otic",
      grepl("Salve|Krem|Pasta", legemiddelform_kort_dn) ~ "Topical",
      TRUE ~ NA_character_
    ),
    ActiveIngredient = english_name)


active_ingredient_summary <- antibiotics_report_4 %>%
  group_by(MedicalProductCode) %>%
  summarise(
    has_multiple_active_ingredients = n_distinct(ActiveIngredient) > 1,
    formatted_active_ingredients = ifelse(n_distinct(ActiveIngredient) > 1, paste(sort(unique(ActiveIngredient)), collapse = ", "), first(ActiveIngredient))
  )

antibiotics_report_5 <- antibiotics_report_4 %>%
  left_join(active_ingredient_summary, by = "MedicalProductCode")


antibiotics_report_7 <- antibiotics_report_5 %>%
  mutate(
    # Add "_TMP" if there are multiple active ingredients, it contains "trimethoprim", and active contains "sulfa"
    ActiveIngredient = ifelse(
      has_multiple_active_ingredients == TRUE & 
        str_detect(formatted_active_ingredients, "Trimethoprim") & 
        str_detect(ActiveIngredient, regex("sulfa", ignore_case = TRUE)),
      paste0(ActiveIngredient, "_TMP"),
      ActiveIngredient
    ),
    # Add "_sulfa" if there are any sulfa medications in formatted ingredients and active ingredient is "trimethoprim"
    ActiveIngredient = ifelse(
      str_detect(formatted_active_ingredients, regex("sulfa", ignore_case = TRUE)) & 
        ActiveIngredient == "Trimethoprim",
      paste0(ActiveIngredient, "_sulfa"),
      ActiveIngredient
    )
  )


antibiotics_report_8 <- antibiotics_report_7 %>%
  mutate(
    # Fix 1 & 2: Assign Parenteral route for injectable products identified by "SA" in
    # product name, OR for Vetmast (intramuscular injectable despite the mastitis indication)
    Route = ifelse(
      is.na(Route) &
        (str_detect(MedicalProductName, regex("\\bSA\\b")) |
           MedicalProductName == "Vetmast"),
      "Parenteral",
      Route
    ),
    AdministrationMethod = ifelse(
      is.na(AdministrationMethod) &
        (str_detect(MedicalProductName, regex("\\bSA\\b")) |
           MedicalProductName == "Vetmast"),
      "Injectable",
      AdministrationMethod
    )
  )


antibiotics_imputed_ddd <- antibiotics_report_8 %>%
  left_join(dddvet, by = c("Route", "ActiveIngredient" = "Substance"))

antibiotics_imputed_ddd <- antibiotics_imputed_ddd %>%
  mutate(
    DDDvet = ifelse(grepl("intramammarie", legemiddelform_kort_dn, ignore.case = TRUE), 1, DDDvet),
    Unit   = ifelse(grepl("intramammarie", legemiddelform_kort_dn, ignore.case = TRUE), "UD/teat", Unit),
    
    DDDvet = ifelse(grepl("uteritorie", legemiddelform_kort_dn, ignore.case = TRUE), 1, DDDvet),
    Unit   = ifelse(grepl("uteritorie", legemiddelform_kort_dn, ignore.case = TRUE), "UD/animal", Unit)
  )


missing_matches <- antibiotics_imputed_ddd %>%
  filter(is.na(DDDvet))


missing_matches_final <- antibiotics_imputed_ddd %>%
  filter(is.na(DDDvet)) %>%
  group_by(ActiveIngredient, MedicalProductName, Route, AdministrationMethod, lmp_pakningstype_dn) %>%
  summarise(
    count = n()
  )

matches_final <- antibiotics_imputed_ddd %>%
  group_by(ActiveIngredient, Route, AdministrationMethod) %>%
  summarise(
    count = n()
  )

#================================================================

antibiotics_fixed_ddd <- antibiotics_imputed_ddd |>
  filter(!(is.na(DDDvet)))

antibiotic_use_summary <- antibiotics_fixed_ddd %>%
  group_by(Year, Month, ActiveIngredientCategory) %>%
  summarise(TotalAntibioticKg = sum(ActiveSubstanceKg, na.rm = TRUE)) %>%
  ungroup()

yearly_totals <- antibiotic_use_summary %>%
  group_by(Year) %>%
  summarise(TotalKg = sum(TotalAntibioticKg, na.rm = TRUE)) %>%
  ungroup()

#================================================================

cattle_pop <- data.frame(
  Year = c(2018, 2019, 2020, 2021, 2022, 2023, 2024),
  Total_Cattle = c(877388, 862550, 876776, 882968, 893771, 881320, 851697),
  Total_Cows = c(309781, 307520, 306076, 317073, 319590, 311300, 309108),
  Dairy_Cows = c(218613, 215069, 207855, 212629, 211058, 202876, 202771),
  Beef_Cows = c(91168, 92451, 98221, 104444, 108532, 108424, 106337),
  Other_Cattle = c(567607, 555030, 570700, 565895, 574181, 570020, 542589)
)
cattle_pop$Year = as.character(cattle_pop$Year)



animalia_herd_en <- animalia_herd %>%
  rename(
    producer_id           = PRODNR10,
    year                  = AAR,
    yearly_cows           = ÅRSKYR,
    live_sold_bulls       = SOLGT_LIV_OKSE,
    live_sold_heifers     = SOLGT_LIV_KVIGE,
    purchased_bulls       = KJØPT_OKSE,
    purchased_heifers     = KJØPT_KVIGE,
    other_exit_bulls      = ANNEN_UTMELDING_OKSE,
    other_exit_heifers    = ANNEN_UTMELDING_KVIGE,
    slaughtered_bulls     = SLAKT_OKSE,
    slaughtered_heifers   = SLAKT_KVIGE,
    heifer_calves_born    = KALVING_ANT_KVIGEKALVER,
    bull_calves_born      = KALVING_ANT_OKSEKALVER,
    stillborn_count       = ANT_DODFODT,
    died_count            = ANT_KREPERT,
    tagged_died_count     = ANT_KREPERT_MERKET
  ) %>%
  mutate(producer_id = as.character(producer_id))

mimiro_herd_en <- mimiro_herd %>%
  rename(
    producer_id           = PRODUSENT,
    organization_number   = `Organization Number`,
    owner_id              = OWNER_ID,
    postal_code           = POSTALCODE,
    yearly_bulls          = `Yearly Bull`,
    yearly_cows           = `Yearly Cows`,
    yearly_heifers        = `Yearly Heifers`,
    yearly_milking_cows   = `Yearly Milking Cows`,
    amount                = AMOUNT,
    conventional_amount   = CONVENTIONAL,
    organic_amount        = ORGANIC,
    year                  = YEAR,
    milking_system        = `Milking System`,
    barn_type             = `Barn Type`
  ) %>%
  mutate(producer_id = as.character(producer_id))


antibiotics_pop <- antibiotics_fixed_ddd %>%
  left_join(
    cattle_pop, by = "Year"
  )

##################

### taotal kg by month ####

# Calculate yearly totals

library(paletteer)
tab10 <- as.character(paletteer_d("ggthemes::Tableau_10"))

# Order matches Figure 5's factor levels — each class keeps its Fig 5 color
class_palette <- c(
  "Beta-lactamase-sensitive penicillins" = tab10[1],
  "Tetracyclines"                        = tab10[2],
  "Broad-spectrum penicillins"           = tab10[3],
  "Sulfonamides"                         = tab10[4],
  "Diaminopyrimidines"                   = tab10[5],
  "Aminoglycosides"                      = tab10[6],
  "Amphenicols"                          = tab10[7],
  "Fluoroquinolones"                     = tab10[8],
  "Beta-lactamase-resistant penicillins" = tab10[9],
  "Macrolides"                           = tab10[10]
)

# Define the temperature palette
temperature_palette <- paletteer_d("ggthemes::Tableau_10")

# Aggregate data by Year, Month, and ActiveIngredientCategory
antibiotic_use_summary <- antibiotics_pop %>%
  group_by(Year, Month, ActiveIngredientCategory) %>%
  summarise(TotalAntibioticKg = sum(ActiveSubstanceKg, na.rm = TRUE)) %>%
  ungroup()

yearly_totals <- antibiotic_use_summary %>%
  group_by(Year) %>%
  summarise(TotalKg = sum(TotalAntibioticKg, na.rm = TRUE)) %>%
  ungroup()

yearly_monthly_totals <- antibiotic_use_summary %>%
  group_by(Year, Month) %>%
  summarise(TotalKg = sum(TotalAntibioticKg, na.rm = TRUE)) %>%
  ungroup()

# Merge yearly_totals with antibiotic_use_summary to include total per year
antibiotic_use_summary <- antibiotic_use_summary %>%
  left_join(yearly_totals, by = "Year")

# Create a named vector for the custom labels with yearly totals
year_labels <- yearly_totals %>%
  mutate(label = paste0(Year, "\n", round(TotalKg, 1), " kg")) %>%
  dplyr::select(Year, label) %>%
  deframe()

# Create the plot with custom facet labels
p <- ggplot(antibiotic_use_summary, aes(x = Month, y = TotalAntibioticKg, fill = ActiveIngredientCategory)) +
  geom_bar(stat = "identity", position = "stack") +
  
  # Apply custom facet labels using the labeller
  facet_grid(cols = vars(Year), scales = "free_y", 
             labeller = labeller(Year = year_labels)) +
  
  scale_fill_manual(values = class_palette, drop = TRUE) +
  labs(
    title = "",
    x = "Month",
    y = "kg active substance",
    fill = "Active Ingredient"
  ) +
  scale_y_continuous(breaks = pretty_breaks(n = 10), expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),
    text = element_text(family = "Arial", color = "black"),
    axis.title.x = element_text(size = 20, color = "black"),
    axis.title.y = element_text(size = 20, color = "black"),
    strip.text = element_text(color = "black", size = 20),  # Facet label text
    axis.text.x = element_text(color = "black", size = 12, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 20, color = "black"),
   # panel.grid = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
  # panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "right",
  )

p

ggsave("use_by_month_total.png", p, width = 20, height = 8)


yearly_biomass <- antibiotics_pop %>%
  group_by(Year) %>%
  summarise(biomass_kg = first(Total_Cattle) * 425) %>%   # 425 kg avg live weight
  ungroup()

antibiotic_use_disease <- antibiotics_pop %>%
  group_by(Year, main_category) %>%
  summarise(TotalAntibioticKg = sum(ActiveSubstanceKg, na.rm = TRUE), .groups = "drop") %>%
  left_join(yearly_biomass, by = "Year") %>%
  mutate(mg_per_kg_biomass = (TotalAntibioticKg * 1e6) / biomass_kg)   # kg → mg: ×1,000,000

antibiotic_use_route <- antibiotics_pop %>%
  group_by(Year, Route) %>%
  summarise(TotalAntibioticKg = sum(ActiveSubstanceKg, na.rm = TRUE), .groups = "drop") %>%
  left_join(yearly_biomass, by = "Year") %>%
  mutate(mg_per_kg_biomass = (TotalAntibioticKg * 1e6) / biomass_kg)

top3_diseases <- antibiotic_use_disease %>%
  group_by(main_category) %>%
  summarise(total = sum(mg_per_kg_biomass)) %>%
  slice_max(total, n = 3) %>%
  pull(main_category)

disease_plot_data <- antibiotic_use_disease %>%
  mutate(main_category = if_else(main_category %in% top3_diseases, main_category, "Other")) %>%
  group_by(Year, main_category) %>%
  summarise(mg_per_kg_biomass = sum(mg_per_kg_biomass), .groups = "drop")

shared_theme <- list(
  theme_minimal(),
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    axis.text.x = element_text(color = "black", size = 20, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 20, color = "black", hjust = 1),
    strip.text = element_text(color = "black", size = 20),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.position = "bottom",
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )
)

p1 <- ggplot(disease_plot_data, aes(x = Year, y = mg_per_kg_biomass, color = main_category, group = main_category)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_y_continuous(breaks = pretty_breaks(n = 10),
                     limits = c(0, NA)
                     ) +
  scale_x_discrete(expand = expansion(add = c(0.1, 0.1))) +    # <- discrete, not continuous
  scale_color_brewer(palette = "Set1") +
  labs(title = NULL, y = "mg/kg animal biomass", x = "Year", color = "Diagnosis category") +
  guides(color = guide_legend(ncol = 1)) +
  shared_theme

p2 <- ggplot(antibiotic_use_route, aes(x = Year, y = mg_per_kg_biomass, color = Route, group = Route)) +
  geom_line(linewidth = 1) + geom_point(size = 2) +
  scale_y_continuous(breaks = pretty_breaks(n = 10)) +
  scale_x_discrete(expand = expansion(add = c(0.1, 0.1))) +    # <- discrete, not continuous
  
  scale_color_brewer(palette = "Dark2") +
  labs(title = NULL, y = "mg/kg animal biomass", x = "Year", color = "Route") +
  guides(color = guide_legend(ncol = 1)) +
  shared_theme

p1 + p2
ggsave("use_mg_kg_biomass.png", p1 + p2, width = 10, height = 6)

common_theme <- theme(
  plot.tag = element_text(size = 20, face = "bold")
)

p_fixed  <- p  + common_theme
p1_fixed <- p1 + common_theme
p2_fixed <- p2 + common_theme

# free() releases the top plot from panel-width alignment with the bottom row
combined_plot_fig_3 <- (free(p_fixed) / (p1_fixed | p2_fixed)) +
  plot_layout(heights = c(3, 1)) +
  plot_annotation(
    tag_levels = list(c("A", "B", "C"))
    )

combined_plot_fig_3
ggsave("Figure_3_combined.png", combined_plot_fig_3, width = 16, height = 12, dpi = 300)
ggsave("Figure_3.tiff", combined_plot_fig_5, width = 16, height = 12, dpi = 300)


mg_per_kg_biomass <- antibiotics_pop %>%
  group_by(Year) %>%
  summarise(TotalAntibioticKg = sum(ActiveSubstanceKg, na.rm = TRUE), .groups = "drop") %>%
  left_join(yearly_biomass, by = "Year") %>%
  mutate(mg_per_kg_biomass = (TotalAntibioticKg * 1e6) / biomass_kg)

mg_per_kg_biomass_active_sub <- antibiotics_pop %>%
  group_by(Year, ActiveIngredientCategory) %>%
  summarise(TotalAntibioticKg = sum(ActiveSubstanceKg, na.rm = TRUE), .groups = "drop") %>%
  left_join(yearly_biomass, by = "Year") %>%
  mutate(mg_per_kg_biomass = (TotalAntibioticKg * 1e6) / biomass_kg)


mg_per_kg_biomass_indication <- antibiotics_pop %>%
  group_by(Year, main_category) %>%
  summarise(TotalAntibioticKg = sum(ActiveSubstanceKg, na.rm = TRUE), .groups = "drop") %>%
  left_join(yearly_biomass, by = "Year") %>%
  mutate(mg_per_kg_biomass = (TotalAntibioticKg * 1e6) / biomass_kg)


mg_per_kg_biomass_route <- antibiotics_pop %>%
  group_by(Year, Route) %>%
  summarise(TotalAntibioticKg = sum(ActiveSubstanceKg, na.rm = TRUE), .groups = "drop") %>%
  left_join(yearly_biomass, by = "Year") %>%
  mutate(mg_per_kg_biomass = (TotalAntibioticKg * 1e6) / biomass_kg)


mg_per_kg_biomass_active_sub |> write.table("clipboard", sep = "\t", row.names = FALSE)
mg_per_kg_biomass_indication |> write.table("clipboard", sep = "\t", row.names = FALSE)
mg_per_kg_biomass_route |> write.table("clipboard", sep = "\t", row.names = FALSE)

#================================================================

##############

# nDDDvet calculations

antibiotics_intram_records <- antibiotics_pop %>%
  filter(grepl("intramam", tolower(legemiddelform_kort_dn)))

intram_counts <- antibiotics_intram_records |>
  mutate(
    n_syringes_intramam = case_when(
      MedicalProductCode == "000002" ~ levert_mengde/0.004,
    
      ##if someone is reading this, important that this must be first — would otherwise incorrectly 
      ## match the veterinarian case below it
      # 4 mL/dose was assumed from concentration (250 mg/mL DHS, standard 1g dose)
      TypeOfDataSource == "Pharmacy sale (prescription)" & lmp_enhet_pakning_v == "sprøyte"
      ~ antall_pakninger * lmp_antall * lmp_mengde,
      TypeOfDataSource == "Pharmacy sale (prescription)" & lmp_enhet_pakning_v != "sprøyte"
      ~ antall_pakninger * lmp_antall,
      TypeOfDataSource == "Journal entry (use by veterinarian)" & lmp_enhet_pakning_v == "sprøyte" 
      & enhet_mengde %in% c("spr", "stk", "sprøyte") ~ levert_mengde,
      TypeOfDataSource == "Journal entry (use by veterinarian)" & lmp_enhet_pakning_v != "sprøyte" 
      & enhet_mengde %in% c("spr", "stk", "sprøyte") ~ levert_mengde,
      TypeOfDataSource == "Journal entry (use by veterinarian)" & lmp_enhet_pakning_v != "sprøyte" 
      & enhet_mengde == lmp_enhet_pakning_v ~ levert_mengde / lmp_mengde,
      TRUE ~ NA_real_
    )
  )

intram_counts_summary <- intram_counts |>
  distinct(unique_row_id, .keep_all = TRUE) |>
  group_by(Year, TypeOfDataSource, formatted_active_ingredients, MedicalProductName, MedicalProductCode,
           main_category, drug_class,Total_Cows,legemiddelform_kort_dn, lmp_enhet_pakning_v, lmp_mengde, lmp_antall ) |>
  summarize(total_syringes = sum(n_syringes_intramam, na.rm = TRUE))

dry_cow_names <- c(
  "Ubrostar dry cow vet", 
  "Orbenin extra vet",
  "Orbenin vet",
  "Siccalactin vet",
  "Benestermycin vet"
)

lactating_names <- c(
  "Synulox Comp vet",
  "Synulox LC Plus vet",
  "Carepen vet",
  "Mastipen vet",
  "Prokainpenicillin/dihydrostreptomycin",
  "Streptocillin Forte vet",
  "Streptocillin vet"
)


intram_counts_summary <- intram_counts_summary |>
  mutate(
    cow_treatment_type = case_when(
      MedicalProductName %in% dry_cow_names ~ "Dry cow",
      MedicalProductName %in% lactating_names ~ "Lactating cow",
      TRUE ~ NA_character_  #flag anything unexpected
    )
  )


intram_ddd_summary <- intram_counts_summary |>
  group_by(Year, formatted_active_ingredients, cow_treatment_type) |>
  summarize(
    total_syringes = sum(total_syringes, na.rm = TRUE),
    Total_Cows = first(Total_Cows),
    .groups = "drop"
  ) |>
  mutate(
    ddd_divisor = if_else(cow_treatment_type == "Dry cow", 4, 1),
    doses_per_1000_cows = ((total_syringes / ddd_divisor) / Total_Cows) * 1000
    )

all_ingredients <- intram_ddd_summary |>
  pull(formatted_active_ingredients) |>
  unique() |>
  sort()



antibiotics_iup <- antibiotics_pop |>
  filter(grepl("uteritorie", tolower(legemiddelform_kort_dn))) |>
  mutate(levert_mengde_adj = if_else(
    grepl("pharmacy", tolower(TypeOfDataSource)),
    antall_pakninger * lmp_antall * lmp_mengde,
    levert_mengde
  ))

iup_dddvet_total <- antibiotics_iup |>
  group_by(Year, formatted_active_ingredients) |>
  summarise(
    total_iup = sum(levert_mengde_adj),
    Total_Cows = first(Total_Cows),
    dddvet_per_1000_cows = total_iup * 1000 / Total_Cows
  )


antibiotics_non_intram_records <- antibiotics_pop %>%
  filter(!grepl("intramam", tolower(legemiddelform_kort_dn)),
         !grepl("uteritorie", tolower(legemiddelform_kort_dn))
  )


tk_calculated <- antibiotics_non_intram_records %>%
  mutate(tk_dddvet = (as.numeric(ActiveSubstanceKg) * 1e6) / as.numeric(DDDvet),
  )



# group by drug_class instead of ActiveIngredient — at the national level
# there are too many ingredients to read in a stacked bar
national_dddvet_by_class <- tk_calculated |>
  group_by(drug_class, Year, Route) |>
  summarise(
    total_weight_kg            = first(Total_Cattle) * 425,
    sum_tk_dddvet         = sum(tk_dddvet, na.rm = TRUE),
    dddvet_per_kg_biomass      = sum_tk_dddvet / total_weight_kg,
    .groups = "drop"
  )


#####################

###plotting #######

all_ingredients <- sort(unique(c(
  intram_ddd_summary$formatted_active_ingredients,
  iup_dddvet_total$formatted_active_ingredients
)))

# Map each ingredient to its drug class; combination products are assigned
# to the penicillin class by convention (all combos here are penicillin-based)
ingredient_class <- c(
  "Amoxicillin"                                          = "Broad-spectrum penicillins",
  "Benzylpenicillin"                                     = "Beta-lactamase-sensitive penicillins",
  "Benzylpenicillin, Dihydrostreptomycin"                = "Beta-lactamase-sensitive penicillins",
  "Benzylpenicillin, Framycetin, Penethamate hydriodide" = "Beta-lactamase-sensitive penicillins",
  "Cloxacillin"                                          = "Beta-lactamase-resistant penicillins",
  "Dihydrostreptomycin, Procaine benzylpenicillin"       = "Beta-lactamase-sensitive penicillins",
  "Oxytetracycline"                                      = "Tetracyclines"
)
if (!all(all_ingredients %in% names(ingredient_class))) {
  stop("Unmapped ingredient(s): ",
       paste(setdiff(all_ingredients, names(ingredient_class)), collapse = ", "))
}

# Lighten a colour towards white by `prop` (0 = unchanged, 1 = white)
mix_white <- function(col, prop) {
  v <- col2rgb(col) * (1 - prop) + 255 * prop
  rgb(v[1], v[2], v[3], maxColorValue = 255)
}

# Ingredients share their class colour from class_palette (matching the
# national/kg figures); siblings within a class get lightening steps so
# they remain distinguishable in stacked bars
shared_palette <- vapply(all_ingredients, function(ing) {
  cls      <- ingredient_class[[ing]]
  siblings <- sort(names(ingredient_class)[ingredient_class == cls])
  mix_white(class_palette[[cls]], (match(ing, siblings) - 1) * 0.25)
}, character(1)) |>
  setNames(all_ingredients)

plot_theme <- theme_minimal() +
  theme(
    text             = element_text(family = "Arial", color = "black"),
    axis.title       = element_text(size = 15, color = "black", face = "bold"),
    axis.text.x      = element_text(color = "black", size = 15,
                                    angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y      = element_text(size = 15, color = "black", hjust = 1),
    strip.text       = element_text(color = "black", size = 20, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.text      = element_text(size = 10),
    legend.title     = element_text(size = 10),
    legend.position  = "bottom",
    plot.title       = element_text(size = 15, face = "bold", hjust = 0.5), 
    panel.border     = element_rect(color = "black", fill = NA, size = 1),
    plot.tag         = element_text(size = 20, face = "bold")
  )

fill_scale_products <- scale_fill_manual(
  values = shared_palette,
  drop   = TRUE,
  name   = "Active substance",
  guide  = guide_legend(ncol = 1)
)

p_dry <- ggplot(
  intram_ddd_summary |> filter(cow_treatment_type == "Dry cow"),
  aes(x = factor(Year), y = doses_per_1000_cows, fill = formatted_active_ingredients)
) +
  geom_col() +
  fill_scale_products +
  scale_y_continuous(breaks = pretty_breaks(n = 6), limits = c(0, NA)) +
  labs(
    title = "Intramammary products for dry cows",
    x = "Year", y = "nDCDvet per 1000 cows",
  ) +
  plot_theme

p_lactating <- ggplot(
  intram_ddd_summary |> filter(cow_treatment_type == "Lactating cow"),
  aes(x = factor(Year), y = doses_per_1000_cows, fill = formatted_active_ingredients)
) +
  geom_col() +
  scale_y_continuous(breaks = pretty_breaks(n = 6), limits = c(0, NA)) +
  fill_scale_products +
  labs(
    title = "Intramammary products for lactating cows",
    x = "Year", y = "nDDDvet per 1000 cows",
  ) +
  plot_theme


p_iup <- ggplot(
  iup_dddvet_total,
  aes(x = factor(Year), y = dddvet_per_1000_cows, fill = formatted_active_ingredients)
) +
  geom_col() +
  scale_y_continuous(breaks = pretty_breaks(n = 6), limits = c(0, NA)) +
  fill_scale_products +
  labs(
    title = "Intrauterine products",
    x = "Year", y = "nDDDvet per 1000 cows",
  ) +
  plot_theme

p_national <- ggplot(
  national_dddvet_by_class,
  aes(x = factor(Year), y = dddvet_per_kg_biomass, fill = drug_class)
) +
  scale_fill_manual(values = class_palette, drop = TRUE) +
  geom_col() +
  facet_wrap(~ Route, scales = "free_y") +
  scale_y_continuous(breaks = pretty_breaks(n = 6), limits = c(0, NA)) +
  labs(
    x = "Year", y = "nDDDvet per kg animal biomass",
    fill = "Antimicrobial class"
  ) +
  guides(fill = guide_legend(nrow = 2)) +
  plot_theme

products_row <- (p_dry | p_lactating | p_iup) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# ── combine: national on top (full width), three product plots below ────
combined_plot <- p_national / products_row +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(tag_levels = "A")


combined_plot

ggsave("DDDvet_combined_new.png", combined_plot,
       height = 14, width = 14, dpi = 300)


ggsave("DDDvet_combined_new.tiff", combined_plot,
       height = 14, width = 14)


####################################################
