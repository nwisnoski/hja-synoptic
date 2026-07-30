# All input paths, raw-to-derived names, and filters used by the preparation
# scripts are declared here. Nothing in this file modifies a source workbook.

prep_config <- list(
  paths = list(
    master_environment = "data/hja-env_data_clean.csv",
    fticr_workbook = "data/hja-FTICRMS.xlsx",
    sample_manifest = "data/derived/hja_2016_sample_manifest.csv",
    soil_metadata = "data/hja-synoptic_env-data-soils.csv",
    dada2_output = "results/dada2_2016",
    output_root = "data/derived/analysis_inputs"
  ),
  fticr = list(
    metadata_columns = 1:14,
    site_intensity_columns = 15:74,
    lab_standard_columns = 75:76,
    expected_metadata_headers = c(
      "Measured Ionic Mass (m/z)", "C", "H", "O", "N", "C13", "S", "P",
      "Na", "El_comp", "Class", "NeutralMass", "Error_ppm", "Candidates"
    ),
    metadata_analysis_names = c(
      "measured_mz", "carbon_count", "hydrogen_count", "oxygen_count",
      "nitrogen_count", "c13_indicator", "sulfur_count", "phosphorus_count",
      "sodium_count", "elemental_composition_source",
      "molecular_class_source", "neutral_mass", "mass_error_ppm",
      "candidate_count"
    ),
    primary_filter = list(
      minimum_mass = 200,
      maximum_mass = 900,
      maximum_c13_indicator = 0,
      minimum_detected_site_count = 2,
      formula_assignment_rule = "carbon_count > 0"
    )
  )
)

environment_variable <- function(
    source_sample_type, source_column, analysis_column, block, units = "",
    recommended_transformation = "none", primary_44_site = TRUE, notes = "") {
  data.frame(
    source_file = prep_config$paths$master_environment,
    source_sample_type = source_sample_type,
    source_column = source_column,
    analysis_column = analysis_column,
    predictor_block = block,
    units = units,
    recommended_transformation = recommended_transformation,
    primary_44_site = primary_44_site,
    notes = notes,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}
ev <- environment_variable

environment_variable_map <- do.call(rbind, list(
  # Landscape and network position. Site-level values are read from Stream rows.
  ev("Stream", "Stream_Centerline_UTM_Zone_10N_X_m", "site_utm_x_m", "landscape", "m"),
  ev("Stream", "Stream_Centerline_UTM_Zone_10N_Y_m", "site_utm_y_m", "landscape", "m"),
  ev("Stream", "Segment_Number", "site_segment_number", "landscape"),
  ev("Stream", "Drainage_area_ha", "site_drainage_area_ha", "landscape", "ha", "log10"),
  ev("Stream", "Valley_Slope_m/m", "site_valley_slope_m_m", "landscape", "m/m"),
  ev("Stream", "Stream_Slope_m/m", "site_stream_slope_m_m", "landscape", "m/m"),
  ev("Stream", "Valley_Width_m", "site_valley_width_m", "landscape", "m"),
  ev("Stream", "Segment_Valley_Avg_Slope_m/m", "site_segment_valley_avg_slope_m_m", "landscape", "m/m"),
  ev("Stream", "Segment_Stream_Avg_Slope_m/m", "site_segment_stream_avg_slope_m_m", "landscape", "m/m"),
  ev("Stream", "Segment_Valley_Avg_Width_m", "site_segment_valley_avg_width_m", "landscape", "m"),
  ev("Stream", "Stream_Order", "site_stream_order", "landscape", "ordinal"),
  ev("Stream", "Segment_Sinuosity_m/m", "site_segment_sinuosity_m_m", "landscape", "m/m"),
  ev("Stream", "Valley_Distance_to_outlet_m", "site_distance_to_outlet_m", "landscape", "m", "log10"),

  # Hydrology and exchange.
  ev("Sediment", "geometric_Mean_Hydaulic_Conductivity_m/s", "sediment_hydraulic_conductivity_geomean_m_s", "hydrology", "m/s", "log10"),
  ev("Sediment", "Piezo_depth_below_streambed_(to_top)_m", "sediment_piezometer_depth_top_m", "hydrology", "m"),
  ev("Sediment", "Piezo_depth_below_streambed_(to_bottom)_m", "sediment_piezometer_depth_bottom_m", "hydrology", "m"),
  ev("Stream", "Channel_width", "site_channel_width_m", "hydrology", "m"),
  ev("Stream", "Channel_depth_m", "site_channel_depth_m", "hydrology", "m"),
  ev("Stream", "Tracer_Downstream_discahrge_m3/s", "site_discharge_downstream_m3_s", "hydrology", "m3/s", "log10", FALSE, "Reduced sensitivity subset."),
  ev("Stream", "Tracer_Upstream_discharge_m3/s", "site_discharge_upstream_m3_s", "hydrology", "m3/s", "log10", FALSE, "Reduced sensitivity subset."),
  ev("Stream", "Peak_travel_time_hr", "site_peak_travel_time_hr", "hydrology", "hr", "log10", FALSE),
  ev("Stream", "Velocity_m/s", "site_velocity_m_s", "hydrology", "m/s", "log10", FALSE),

  # Surface-water chemistry, absorbance, and fluorescence.
  ev("Stream", "Water_Temp_arrival_deg_C", "surface_water_temperature_c", "surface_water", "deg C"),
  ev("Stream", "DO_arrival_mg/L", "surface_do_mg_l", "surface_water", "mg/L"),
  ev("Stream", "DO_arrival_perc_sat", "surface_do_percent_saturation", "surface_water", "%"),
  ev("Stream", "NPOC_mg/L", "surface_npoc_mg_l", "surface_dom", "mg/L", "log10"),
  ev("Stream", "SUVA254", "surface_suva254", "surface_dom"),
  ev("Stream", "Spectral_Slope_Ratio", "surface_spectral_slope_ratio", "surface_dom"),
  ev("Stream", "TDN_mg/L", "surface_tdn_mg_l", "surface_nutrients", "mg/L", "log10"),
  ev("Stream", "Integrated_Absorbance_300-400_nm_1/m", "surface_integrated_absorbance_300_400_1_m", "surface_dom", "1/m", "log10"),
  ev("Stream", "Decadic_Absorbance_at_254_nm_1/m", "surface_absorbance_254_1_m", "surface_dom", "1/m", "log10"),
  ev("Stream", "Decadic_Absorbance_at_300_nm_1/m", "surface_absorbance_300_1_m", "surface_dom", "1/m", "log10"),
  ev("Stream", "Naperian_Absorbance_at_340_nm_1/m", "surface_absorbance_340_1_m", "surface_dom", "1/m", "log10"),
  ev("Stream", "Naperian_Absorbance_at_380_nm_1/m", "surface_absorbance_380_1_m", "surface_dom", "1/m", "log10"),
  ev("Stream", "Spectral_Slope_Ratio_(Slope_of_Abs275-295/300-350)", "surface_spectral_slope_ratio_275_295_to_300_350", "surface_dom"),
  ev("Stream", "Integrated_Fluorescence_Raman_units", "surface_integrated_fluorescence_ru", "surface_dom", "Raman units", "log10"),
  ev("Stream", "Intensity_of_EEM_peak_A_(250-450)_Raman_units", "surface_eem_peak_a_ru", "surface_dom", "Raman units", "log10"),
  ev("Stream", "Intensity_of_EEM_peak_C_(350-450)_Raman_units", "surface_eem_peak_c_ru", "surface_dom", "Raman units", "log10"),
  ev("Stream", "Intensity_of_EEM_peak_T_(275-340)_Raman_units", "surface_eem_peak_t_ru", "surface_dom", "Raman units", "log10"),
  ev("Stream", "Fluorescence_Index_(Cory_et_al._2010)", "surface_fluorescence_index", "surface_dom"),
  ev("Stream", "Chloride_mg/L", "surface_chloride_mg_l", "surface_ions", "mg/L", "log10"),
  ev("Stream", "Sulfate_mg/L", "surface_sulfate_mg_l", "surface_ions", "mg/L", "log10"),
  ev("Stream", "Sodium_mg/L", "surface_sodium_mg_l", "surface_ions", "mg/L", "log10"),
  ev("Stream", "Potassium_mg/L", "surface_potassium_mg_l", "surface_ions", "mg/L", "log10"),
  ev("Stream", "Magnesium_mg/L", "surface_magnesium_mg_l", "surface_ions", "mg/L", "log10"),
  ev("Stream", "Calcium_mg/L", "surface_calcium_mg_l", "surface_ions", "mg/L", "log10"),
  ev("Stream", "PO4_mg/L", "surface_po4_mg_l", "surface_nutrients", "mg/L", "log10"),
  ev("Stream", "Surface_NO2_NO3_mg/L_as_NO3", "surface_no2_no3_mg_l_as_no3", "surface_nutrients", "mg/L", "log10"),
  ev("Stream", "NH3_mg/L_as_NH3", "surface_nh3_mg_l_as_nh3", "surface_nutrients", "mg/L", "log10"),

  # Hyporheic-water chemistry, absorbance, and fluorescence.
  ev("Hyporheic", "Water_Temp_arrival_deg_C", "hyporheic_water_temperature_c", "hyporheic_water", "deg C"),
  ev("Hyporheic", "DO_arrival_mg/L", "hyporheic_do_mg_l", "hyporheic_water", "mg/L"),
  ev("Hyporheic", "DO_arrival_perc_sat", "hyporheic_do_percent_saturation", "hyporheic_water", "%"),
  ev("Hyporheic", "NPOC_mg/L", "hyporheic_npoc_mg_l", "hyporheic_dom", "mg/L", "log10"),
  ev("Hyporheic", "SUVA254", "hyporheic_suva254", "hyporheic_dom"),
  ev("Hyporheic", "Spectral_Slope_Ratio", "hyporheic_spectral_slope_ratio", "hyporheic_dom"),
  ev("Hyporheic", "TDN_mg/L", "hyporheic_tdn_mg_l", "hyporheic_nutrients", "mg/L", "log10"),
  ev("Hyporheic", "Integrated_Absorbance_300-400_nm_1/m", "hyporheic_integrated_absorbance_300_400_1_m", "hyporheic_dom", "1/m", "log10"),
  ev("Hyporheic", "Decadic_Absorbance_at_254_nm_1/m", "hyporheic_absorbance_254_1_m", "hyporheic_dom", "1/m", "log10"),
  ev("Hyporheic", "Intensity_of_EEM_peak_A_(250-450)_Raman_units", "hyporheic_eem_peak_a_ru", "hyporheic_dom", "Raman units", "log10"),
  ev("Hyporheic", "Intensity_of_EEM_peak_C_(350-450)_Raman_units", "hyporheic_eem_peak_c_ru", "hyporheic_dom", "Raman units", "log10"),
  ev("Hyporheic", "Intensity_of_EEM_peak_T_(275-340)_Raman_units", "hyporheic_eem_peak_t_ru", "hyporheic_dom", "Raman units", "log10"),
  ev("Hyporheic", "Fluorescence_Index_(Cory_et_al._2010)", "hyporheic_fluorescence_index", "hyporheic_dom"),
  ev("Hyporheic", "Chloride_mg/L", "hyporheic_chloride_mg_l", "hyporheic_ions", "mg/L", "log10"),
  ev("Hyporheic", "Sulfate_mg/L", "hyporheic_sulfate_mg_l", "hyporheic_ions", "mg/L", "log10"),
  ev("Hyporheic", "Sodium_mg/L", "hyporheic_sodium_mg_l", "hyporheic_ions", "mg/L", "log10"),
  ev("Hyporheic", "Potassium_mg/L", "hyporheic_potassium_mg_l", "hyporheic_ions", "mg/L", "log10"),
  ev("Hyporheic", "Magnesium_mg/L", "hyporheic_magnesium_mg_l", "hyporheic_ions", "mg/L", "log10"),
  ev("Hyporheic", "Calcium_mg/L", "hyporheic_calcium_mg_l", "hyporheic_ions", "mg/L", "log10"),
  ev("Hyporheic", "PO4_mg/L", "hyporheic_po4_mg_l", "hyporheic_nutrients", "mg/L", "log10"),
  ev("Hyporheic", "Surface_NO2_NO3_mg/L_as_NO3", "hyporheic_no2_no3_mg_l_as_no3", "hyporheic_nutrients", "mg/L", "log10"),
  ev("Hyporheic", "NH3_mg/L_as_NH3", "hyporheic_nh3_mg_l_as_nh3", "hyporheic_nutrients", "mg/L", "log10"),

  # Sediment material and enzyme-activity measurements.
  ev("Sediment", "sediment_Organic_Content_%", "sediment_organic_content_percent", "sediment", "%"),
  ev("Sediment", "sediment_Organic_Content_%_R2", "sediment_organic_content_duplicate_percent", "sediment_qc", "%", "none", FALSE),
  ev("Sediment", "sediment_NAG_µmol/(g*hr)", "sediment_nag_umol_g_hr", "sediment_enzyme", "umol/(g*hr)", "log10"),
  ev("Sediment", "sediment_LAP_µmol/(g*hr)", "sediment_lap_umol_g_hr", "sediment_enzyme", "umol/(g*hr)", "log10"),
  ev("Sediment", "sediment_GLU_µmol/(g*hr)", "sediment_glu_umol_g_hr", "sediment_enzyme", "umol/(g*hr)", "log10"),
  ev("Sediment", "sediment_AP_µmol/(g*hr)", "sediment_ap_umol_g_hr", "sediment_enzyme", "umol/(g*hr)", "log10"),
  ev("Sediment", "sediment_NAG_(DUP)_µmol/(g*hr)", "sediment_nag_duplicate_umol_g_hr", "sediment_qc", "umol/(g*hr)", "none", FALSE),
  ev("Sediment", "sediment_LAP_(DUP)_µmol/(g*hr)", "sediment_lap_duplicate_umol_g_hr", "sediment_qc", "umol/(g*hr)", "none", FALSE),
  ev("Sediment", "sediment_GLU_(DUP)_µmol/(g*hr)", "sediment_glu_duplicate_umol_g_hr", "sediment_qc", "umol/(g*hr)", "none", FALSE),
  ev("Sediment", "sediment_AP_(DUP)_µmol/(g*hr)", "sediment_ap_duplicate_umol_g_hr", "sediment_qc", "umol/(g*hr)", "none", FALSE),
  ev("Sediment", "D10_um", "sediment_d10_um", "sediment_texture", "um", "log10", FALSE, "Reduced sensitivity subset."),
  ev("Sediment", "D50_um", "sediment_d50_um", "sediment_texture", "um", "log10", FALSE, "Reduced sensitivity subset."),
  ev("Sediment", "D90_um", "sediment_d90_um", "sediment_texture", "um", "log10", FALSE, "Reduced sensitivity subset."),
  ev("Sediment", "perc_gravel", "sediment_gravel_percent", "sediment_texture", "%", "none", FALSE),
  ev("Sediment", "perc_sand", "sediment_sand_percent", "sediment_texture", "%", "none", FALSE),
  ev("Sediment", "perc_mud", "sediment_mud_percent", "sediment_texture", "%", "none", FALSE)
))
rm(ev)
stopifnot(!anyDuplicated(environment_variable_map$analysis_column))
