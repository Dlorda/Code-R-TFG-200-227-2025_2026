
###
#Para instalar el paquete de datos TwoSampleMR:

#install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))

#Librerias necesarias 
library(TwoSampleMR)
library(ieugwasr)
library(openxlsx)
library(dplyr)  
library(viridis)
library(ggplot2)

### 
# Seed: Para generar la misma aleatoriedad cada vez que se ejecute
set.seed(1234)  #Permite reproducibilidad en clumping
###

###
# Modificaciones generales del estudio
###

# pval
s = -6
pval <- as.numeric(paste0("5e", s))

# número métodos en un exposure con resultados significativos para realizar análisis posteriores
min_sig = 2

# Lista de siglas para las enfermedades exposure
names_list <- c("ClcD", "HD", "IBD", "CrD", "UC", "ME", "PV", "RA", "SLE", "T1D")

# Identificador del Outcome en IeuOpenGWAS project
out <- 'ieu-b-2'

# EXPOSURES: Identificador en IeuOpenGWAS project asignado a cada sigla de cada exposure
exp <- c(
  "ClcD" = "ebi-a-GCST000612",
  "HD"   = "ebi-a-GCST90018855",
  "IBD"  = "ieu-a-294",
  "CrD"  = "ieu-a-12",
  "UC"   = "ieu-a-970",
  "ME"   = "ieu-b-18",
  "PV"   = "ebi-a-GCST90019016",
  "RA"   = "ebi-a-GCST90013534",
  "SLE"  = "ebi-a-GCST003156",
  "T1D"  = "ebi-a-GCST90014023"
)

# Métodos MR 

#Acrónimos
method_names_display <- c(
  "MR Egger", 
  "Inverse variance weighted", 
  "Simple median", 
  "Weighted mode",  
  "Weighted median"
)
#Comandos 
method_names_commands <- c(
  "mr_egger_regression", 
  "mr_ivw", 
  "mr_simple_mode",
  "mr_weighted_mode", 
  "mr_weighted_median"
)

method_map <- setNames(method_names_commands, method_names_display)

# Etiquetas de rasgos (para corregir el formato de las columnas) (trait)
trait_names <- c(
  "ClcD" = "Celiac disease || id:ebi-a-GCST000612",
  "HD"   = "Hypothyroidism/Hashimoto’s disease || id:ebi-a-GCST90018855",
  "IBD"  = "Inflammatory bowel disease || id:ieu-a-294",
  "CrD"  = "Crohn's disease || id:ieu-a-12",
  "UC"   = "Ulcerative Colitis || id:ieu-a-970",
  "ME"   = "Multiple sclerosis || id:ieu-b-18",
  "PV"   = "Psoriasis vulgaris || id:ebi-a-GCST90019016",
  "RA"   = "Rheumatoid arthritis || id:ebi-a-GCST90013534",
  "SLE"  = "Systemic lupus erythematosus || id:ebi-a-GCST003156",
  "T1D"  = "Type 1 diabetes || id:ebi-a-GCST90014023"
)

###
# Aquí empieza el estudio 
###

###
# Extracción Exposure 
###
for (i in seq_along(names_list)) {
  if (!names_list[i] %in% names(exp)) {
    warning("No hay ID en 'exp' para: ", names_list[i]); next
  }
  if (!names_list[i] %in% names(trait_names)) {
    warning("No hay 'trait_names' para: ", names_list[i]); next
  }
  exp_id <- unname(exp[names_list[i]])
  
  result <- tryCatch(
    extract_instruments(
      outcomes = exp_id,
      p1 = pval,
      clump = TRUE,
      p2 = pval,
      r2 = 0.001,
      kb = 10000,
      opengwas_jwt = ieugwasr::get_opengwas_jwt(),
      force_server = TRUE
    ),
    error = function(e) { 
      message("Error en extract_instruments para ", names_list[i], " (", exp_id, "): ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(result) || nrow(result) == 0) {
    message("No se han encontrado SNPs significativos para ", names_list[i])
    next
  }
  
  var_name <- paste0(names_list[i], "_exp", i-1, "_dat")
  assign(var_name, result, envir = .GlobalEnv)
  
  # Añadir etiqueta 'trait'
  exposure_trait <- trait_names[names_list[i]]
  tmp <- get(var_name)
  tmp$trait_label <- rep(exposure_trait, nrow(tmp))
  assign(var_name, tmp, envir = .GlobalEnv)
}
message("Se han añadido trait_label a las exposiciones.")

###
# Clumping
###
for (i in seq_along(names_list)) {
  var_name <- paste0(names_list[i], "_exp", i-1, "_dat")
  clumped_var_name <- paste0(names_list[i], "_exp", i-1, "_cdat")
  if (!exists(var_name)) { message("No existe el objeto de exposición: ", var_name); next }
  
  cd <- try(clump_data(get(var_name)), silent = TRUE)
  if (!inherits(cd, "try-error") && !is.null(cd) && nrow(cd) > 0) {
    assign(clumped_var_name, cd, envir = .GlobalEnv)
  } else {
    message("!!Clumping falló o devolvió 0 SNPs para: ", var_name,"!!")
  }
}

###
# Extracción datos Outcome 
###
for (i in seq_along(names_list) ) {
  exp_var <- paste0(names_list[i], "_exp", i-1, "_cdat")  
  out_var <- paste0(names_list[i], i-1,"_out_dat")  
  if (!exists(exp_var)) next
  
  assign(out_var, extract_outcome_data(
    snps = get(exp_var)$SNP,
    outcomes = out),
    envir = .GlobalEnv)
}

###
# Armonización
###
for (i in seq_along(names_list) ) {
  exposure_var <- paste0(names_list[i], "_exp", i-1, "_cdat")
  outcome_var  <- paste0(names_list[i], i-1, "_out_dat")
  if (!exists(exposure_var) || !exists(outcome_var)) next
  
  exposure_dat <- get(exposure_var)
  outcome_dat  <- get(outcome_var)
  
  required_any <- list( 
    SNP = c("SNP"),
    beta = c("beta.outcome","beta"),
    se   = c("se.outcome","se"),
    ea   = c("effect_allele.outcome","effect_allele"),
    oa   = c("other_allele.outcome","other_allele"),
    id   = c("id.outcome","id"),
    out  = c("outcome","trait")
  )
  ok <- all(vapply(required_any, function(opts) any(opts %in% names(outcome_dat)), logical(1)))
  if (!ok) { 
    message("Outcome data para ", names_list[i], " no tiene todas las columnas necesarias. Saltando...")
    next 
  }
  
  harmonised_data <- harmonise_data(
    exposure_dat = exposure_dat,
    outcome_dat  = outcome_dat
  )
  assign(paste0("dat", i-1, "_", names_list[i]), harmonised_data, envir = .GlobalEnv)
}

###
# Realizar MR + conteo de métodos significativos
###
sig_method_counts <- setNames(integer(length(names_list)), names_list)
sig_method_detail  <- setNames(vector("list", length(names_list)), names_list)
names(sig_method_detail) <- names_list

for (i in seq_along(names_list)) {
  dat_name <- paste0("dat", i-1, "_", names_list[i])
  if (!exists(dat_name)) next
  
  result <- mr(get(dat_name), method_list = method_names_commands)
  
  # Contar métodos con p < 0.05 y anotar
  n_sig <- sum(result$pval < 0.05, na.rm = TRUE)
  result$n_sig_methods <- n_sig
  sig_method_counts[names_list[i]] <- n_sig
  sig_method_detail[[names_list[i]]] <- dplyr::mutate(
    result[, c("method", "pval")],
    significant = pval < 0.05
  )
  
  # Guardar
  assign(paste0("res", i-1, "_", s,"_", names_list[i]), result, envir = .GlobalEnv)
}

###
# Función para decidir si seguir con análisis posteriores
###
should_proceed <- function(exposure_key) {
  if (!exists("sig_method_counts", envir = .GlobalEnv)) return(FALSE)
  cnt <- sig_method_counts[exposure_key]
  isTRUE(!is.na(cnt) && cnt >= min_sig)
}

###
# Guardar todos los resultados MR en un Excel con hoja resumen
###
save_all_MR_results_to_excel <- function(filename_prefix = paste0("MR_results", "_", s)) {
  wb <- createWorkbook()
  significant_results <- list()
  greenStyle <- createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
  
  for (i in seq_along(names_list)) {
    res_name <- paste0("res", i-1, "_", s,"_", names_list[i])
    if (!exists(res_name)) next
    
    sheet_name <- names_list[i]
    data <- get(res_name)
    
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, data)
    
    if ("pval" %in% colnames(data)) {
      pval_col <- which(colnames(data) == "pval")
      conditionalFormatting(
        wb, sheet = sheet_name,
        cols = pval_col,
        rows = 2:(nrow(data)+1),
        rule = "<0.05",
        style = greenStyle,
        type = "expression"
      )
      sig_data <- data[data$pval < 0.05, ]
      if (nrow(sig_data) > 0) {
        sig_data$Fuente <- sheet_name
        significant_results[[sheet_name]] <- sig_data
      }
    }
  }
  
  # Añadir hoja de resultados significativos
  if (length(significant_results) > 0) {
    all_significant <- do.call(rbind, significant_results)
    if ("pval" %in% colnames(all_significant)) {
      all_significant <- all_significant[order(all_significant$pval), ]
    }
    addWorksheet(wb, "Resultados significativos")
    writeData(wb, "Resultados significativos", all_significant)
  }
  
  # Añadir hoja resumen con número de métodos significativos
  if (exists("sig_method_counts", envir = .GlobalEnv)) {
    addWorksheet(wb, "Resumen_significancia")
    resumen_df <- data.frame(
      Exposure = names(sig_method_counts),
      n_sig_methods = as.integer(sig_method_counts),
      Proceed = sig_method_counts >= min_sig
    )
    writeData(wb, "Resumen_significancia", resumen_df)
  }
  
  # Guardar archivo evitando sobrescribir
  base_name  <- paste0(filename_prefix, ".xlsx")
  final_name <- base_name
  count <- 1
  while (file.exists(final_name)) {
    final_name <- paste0(filename_prefix, "_", count, ".xlsx")
    count <- count + 1
  }
  saveWorkbook(wb, final_name, overwrite = FALSE)
  message("Resultados guardados en: ", normalizePath(final_name))
}

save_all_MR_results_to_excel()
if (interactive()) browseURL(getwd())

#####
# Análisis de sensibilidad
#####
p <- 0.05

###
## Heterogeneidad
###
for (i in seq_along(names_list) ) {
  if (i == 1){
    res_heterogeneity_egger_list <- list()
    res_heterogeneity_ivw_list   <- list()
    res_heterogeneity_egger_namelist <- list()
    res_heterogeneity_ivw_namelist   <- list()
  }
  
  # Comprobar número de métodos con resultados significativos
  this_exp <- names_list[i]
  if (!should_proceed(this_exp)) {
    message("Saltando heterogeneidad para ", this_exp, 
            " (n_sig_methods = ", sig_method_counts[this_exp], 
            "; umbral = ", min_sig, ").")
    next
  }
  
  res_name <- paste0("res", i-1, "_", s,"_", names_list[i])
  if (!exists(res_name)) { message("No existe: ", res_name); next }
  
  res_i <- get(res_name)
  egger_result <- res_i %>% dplyr::filter(grepl("^MR Egger", method))
  ivw_result   <- res_i %>% dplyr::filter(grepl("^Inverse variance weighted", method))
  
  dataset_name <- paste0("dat", i-1, "_", names_list[i])
  if (!exists(dataset_name)) next
  
  has_egger <- nrow(egger_result) > 0
  has_ivw   <- nrow(ivw_result)   > 0
  if (!has_egger && !has_ivw) {
    message("No hay filas de Egger ni IVW en: ", dataset_name)
    next
  }
  
  if (has_egger && any(egger_result$pval < p, na.rm = TRUE)) {
    heterogeneity_egger_result <- mr_heterogeneity(get(dataset_name), method_list = "mr_egger_regression")
    res_heterogeneity_egger_list[[dataset_name]] <- heterogeneity_egger_result
    res_heterogeneity_egger_namelist[[length(res_heterogeneity_egger_namelist) + 1]] <- dataset_name 
    assign(make.names(paste0("het_egger_", i-1, "_", s,"_", names_list[i])), heterogeneity_egger_result, envir = .GlobalEnv)
  }
  if (has_ivw && any(ivw_result$pval < p, na.rm = TRUE)) {
    heterogeneity_ivw_result <- mr_heterogeneity(get(dataset_name), method_list = "mr_ivw")
    res_heterogeneity_ivw_list[[dataset_name]] <- heterogeneity_ivw_result
    res_heterogeneity_ivw_namelist[[length(res_heterogeneity_ivw_namelist) + 1]] <- dataset_name 
    assign(make.names(paste0("het_ivw_", i-1, "_", s,"_", names_list[i])), heterogeneity_ivw_result, envir = .GlobalEnv)
  }
  
  if (i == length(names_list) ) {
    print(res_heterogeneity_egger_namelist)
    print(res_heterogeneity_ivw_namelist)
    message("Cantidad de elementos en res_heterogeneity_egger_list: ", length(res_heterogeneity_egger_list))
    message("Cantidad de elementos en res_heterogeneity_ivw_list: ", length(res_heterogeneity_ivw_list))
  }
}

###
## Pleiotropía horizontal
###
res_h_pleiotropy_list <- list()
res_h_pleiotropy_namelist <- character(0)

for (i in seq_along(names_list)) {
  # Comprobar número de métodos con resultados significativos
  this_exp <- names_list[i]
  if (!should_proceed(this_exp)) {
    message("Saltando pleiotropía para ", this_exp, 
            " (n_sig_methods = ", sig_method_counts[this_exp], 
            "; umbral = ", min_sig, ").")
    next
  }
  
  res_name <- paste0("res", i-1, "_", s, "_", names_list[i])
  if (!exists(res_name, inherits = TRUE)) { message("No existe objeto de resultados: ", res_name); next }
  
  egger_result <- get(res_name) %>% dplyr::filter(grepl("^MR Egger", method))
  dataset_name <- paste0("dat", i-1, "_", names_list[i])
  if (nrow(egger_result) == 0) { message("Sin filas MR Egger para: ", dataset_name); next }
  if (!any(egger_result$pval < p, na.rm = TRUE)) { message("MR Egger no significativo para: ", dataset_name); next }
  if (!exists(dataset_name, inherits = TRUE)) { message("No existe dataset armonizado: ", dataset_name); next }
  dat_h <- get(dataset_name)
  if (NROW(dat_h) < 3) { message("Menos de 3 SNPs; no pleiotropía: ", dataset_name); next }
  
  pleiotropy_result <- TwoSampleMR::mr_pleiotropy_test(dat_h)
  pleiotropy_result <- as.data.frame(pleiotropy_result)
  pleiotropy_result$dataset <- dataset_name
  
  assign(make.names(paste0("h_pleiotropy_", i-1, "_", s, "_", names_list[i])), pleiotropy_result, envir = .GlobalEnv)
  res_h_pleiotropy_list[[dataset_name]] <- pleiotropy_result
  res_h_pleiotropy_namelist <- c(res_h_pleiotropy_namelist, dataset_name)
}

if (length(res_h_pleiotropy_list)) {
  print(res_h_pleiotropy_namelist)
  message("Cantidad de elementos en res_h_pleiotropy_list: ", length(res_h_pleiotropy_list))
}

## 
# Análisis single SNP (mr_singlesnp) — un data.frame por dataset
##
if (!exists("res_singleSNP_list")) res_singleSNP_list <- list()

for (i in seq_along(names_list)) {
  # Comprobar número de métodos con resultados significativos
  this_exp <- names_list[i]
  if (!should_proceed(this_exp)) {
    message("Saltando heterogeneidad para ", this_exp, 
            " (n_sig_methods = ", sig_method_counts[this_exp], 
            "; umbral = ", min_sig, ").")
    next
  }
  
  res_name <- paste0("res", i-1, "_", s, "_", this_exp)
  if (!exists(res_name)) { message("No existe resultados MR para: ", res_name); next }
  res_mr <- get(res_name)
  if (!all(c("method","pval") %in% names(res_mr))) { message("Resultados MR sin method/pval para ", this_exp); next }
  
  sig_methods_display <- res_mr %>% dplyr::filter(pval < 0.05) %>% dplyr::pull(method)
  if (length(sig_methods_display) == 0) { message("Sin métodos significativos en ", this_exp, ". Saltando…"); next }
  
  sig_cmds <- unname(na.omit(method_map[sig_methods_display]))
  sig_cmds <- unique(sig_cmds)
  if (length(sig_cmds) == 0) { message("Sin comandos válidos entre los métodos significativos en ", this_exp, "."); next }
  
  dat_name <- paste0("dat", i-1, "_", this_exp)
  if (!exists(dat_name)) { message("No existe objeto armonizado: ", dat_name); next }
  dat <- get(dat_name)
  
  obj_single <- mr_singlesnp(
    dat,
    single_method = "mr_wald_ratio",
    all_method    = sig_cmds
  )
  res_singleSNP_list[[dat_name]] <- obj_single
  message("SingleSNP para ", this_exp, " con all_method = {", paste(sig_cmds, collapse = ", "), "}.")
}

message("Cantidad de datasets en res_singleSNP_list: ", length(res_singleSNP_list))

##
# Análisis Leave-one-out
##
res_leaveoneout_list <- list()
valid_loo <- c("mr_ivw", "mr_egger_regression")  # soporte estable
prefer_loo_safe <- c("Inverse_variance_weighted", "MR_Egger")  # nombres display saneados (para plots)

for (i in seq_along(names_list) ) {
  this_exp <- names_list[i]
  if (!should_proceed(this_exp)) {
    message("Saltando Leave-one-out para ", this_exp, 
            " (n_sig_methods = ", sig_method_counts[this_exp], 
            "; umbral = ", min_sig, ").")
    next
  }
  
  res_name <- paste0("res", i-1, "_", s,"_", names_list[i])
  if (!exists(res_name)) { message("No MR para: ", res_name); next }
  
  sig_methods <- get(res_name) %>% dplyr::filter(pval <= 0.05) %>% pull(method)
  if (length(sig_methods) == 0) { message("Sin métodos significativos en: ", names_list[i]); next }
  
  dataset_name <- paste0("dat", i-1, "_", names_list[i])
  if (!exists(dataset_name)) { message("No existe dataset armonizado: ", dataset_name); next }
  dat <- get(dataset_name)
  
  for (method_name in sig_methods) {
    method_command   <- method_map[[method_name]]
    method_name_safe <- gsub("[[:space:]()]", "_", method_name)
    
    if (!is.null(method_command) && method_command %in% valid_loo) {
      obj <- try(mr_leaveoneout(dat, method = get(method_command)), silent = TRUE)
      if (!inherits(obj, "try-error")) {
        if (is.null(res_leaveoneout_list[[dataset_name]])) res_leaveoneout_list[[dataset_name]] <- list()
        res_leaveoneout_list[[dataset_name]][[method_name_safe]] <- obj
        message("Leaveoneout ok: ", dataset_name, " con ", method_command)
      } else {
        message("mr_leaveoneout falló: ", dataset_name, " con ", method_command)
      }
    } else {
      message("Leave-one-out no soportado para: ", method_name, " (", method_command, ")")
    }
  }
}
message("Cantidad de datasets en res_leaveoneout_list: ", length(res_leaveoneout_list))

##
# GUARDAR RESULTADOS Análisis de sensibilidad
##
save_sensitivity_results <- function(filename_prefix = paste0("Sensitivity_results", "_", s)) {
  wb <- openxlsx::createWorkbook()
  p <- 0.05
  greenStyle <- openxlsx::createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
  
  ## 1) HETEROGENEIDAD
  openxlsx::addWorksheet(wb, "Heterogeneity")
  hetero_lists <- list()
  if (exists("res_heterogeneity_egger_list")) hetero_lists <- c(hetero_lists, res_heterogeneity_egger_list)
  if (exists("res_heterogeneity_ivw_list"))   hetero_lists <- c(hetero_lists, res_heterogeneity_ivw_list)
  if (length(hetero_lists) > 0) {
    hetero_results <- tryCatch(dplyr::bind_rows(hetero_lists, .id = "Dataset"), error = function(e) NULL)
  } else {
    hetero_results <- NULL
  }
  if (!is.null(hetero_results) && nrow(hetero_results) > 0) {
    hetero_results$Analysis <- "Heterogeneity"
    openxlsx::writeData(wb, "Heterogeneity", hetero_results)
    if (any(c("Q_pval", "p", "pval") %in% colnames(hetero_results))) {
      col_p <- if ("pval" %in% colnames(hetero_results)) "pval" else if ("p" %in% colnames(hetero_results)) "p" else "Q_pval"
      pval_col <- which(colnames(hetero_results) == col_p)
      openxlsx::conditionalFormatting(
        wb, "Heterogeneity",
        cols = pval_col, rows = 2:(nrow(hetero_results) + 1),
        rule = "<=0.05", style = greenStyle,
        type = "expression"
      )
    }
  } else {
    openxlsx::writeData(wb, "Heterogeneity", "No había datos significativos en los resultados de Egger/IVW para realizar el análisis")
  }
  
  ## 2) PLEIOTROPÍA
  openxlsx::addWorksheet(wb, "Horizontal Pleiotropy")
  if (exists("res_h_pleiotropy_list") && length(res_h_pleiotropy_list) > 0) {
    pleio_results <- tryCatch(dplyr::bind_rows(res_h_pleiotropy_list, .id = "Dataset"), error = function(e) NULL)
    if (!is.null(pleio_results) && nrow(pleio_results) > 0) {
      pleio_results$Analysis <- "Pleiotropy"
      openxlsx::writeData(wb, "Horizontal Pleiotropy", pleio_results)
      if ("pval" %in% colnames(pleio_results)) {
        pval_col <- which(colnames(pleio_results) == "pval")
        openxlsx::conditionalFormatting(
          wb, "Horizontal Pleiotropy",
          cols = pval_col, rows = 2:(nrow(pleio_results) + 1),
          rule = "<=0.05", style = greenStyle,
          type = "expression"
        )
      }
    } else {
      openxlsx::writeData(wb, "Horizontal Pleiotropy", "Se calcularon tests, pero no hay filas combinables.")
    }
  } else {
    openxlsx::writeData(wb, "Horizontal Pleiotropy", "No había datos significativos en MR Egger para realizar el análisis")
  }
  
  ## 3) SINGLE SNP
  openxlsx::addWorksheet(wb, "Single_SNP")
  flatten_single <- list()
  if (exists("res_singleSNP_list") && length(res_singleSNP_list) > 0) {
    for (ds in names(res_singleSNP_list)) {
      df <- res_singleSNP_list[[ds]]
      if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) next
      df$Dataset <- ds
      flatten_single[[ds]] <- df
    }
  }
  if (length(flatten_single) > 0) {
    snp_df <- tryCatch(dplyr::bind_rows(flatten_single, .id = "Key"), error = function(e) NULL)
    if (!is.null(snp_df) && nrow(snp_df) > 0) {
      snp_df$Analysis <- "Single_SNP"
      openxlsx::writeData(wb, "Single_SNP", snp_df)
      if (any(c("p", "pval") %in% colnames(snp_df))) {
        col_p <- if ("pval" %in% colnames(snp_df)) "pval" else "p"
        pval_col <- which(colnames(snp_df) == col_p)
        openxlsx::conditionalFormatting(
          wb, "Single_SNP",
          cols = pval_col, rows = 2:(nrow(snp_df) + 1),
          rule = "<=0.05", style = greenStyle,
          type = "expression"
        )
      }
    } else {
      openxlsx::writeData(wb, "Single_SNP", "No hay datos disponibles para Single_SNP analysis")
    }
  } else {
    openxlsx::writeData(wb, "Single_SNP", "No se encontraron resultados en la lista de Single_SNP")
  }
  
  ## 4) LEAVE-ONE-OUT
  openxlsx::addWorksheet(wb, "Leave-one-out")
  flatten_loo <- list()
  if (exists("res_leaveoneout_list") && length(res_leaveoneout_list) > 0) {
    for (ds in names(res_leaveoneout_list)) {
      for (mtd in names(res_leaveoneout_list[[ds]])) {
        df <- res_leaveoneout_list[[ds]][[mtd]]
        if (is.null(df)) next
        df$Dataset <- ds
        df$Method  <- mtd
        flatten_loo[[paste0(ds, "__", mtd)]] <- df
      }
    }
  }
  if (length(flatten_loo) > 0) {
    leaveoneout_df <- tryCatch(dplyr::bind_rows(flatten_loo, .id = "Key"), error = function(e) NULL)
    if (!is.null(leaveoneout_df) && nrow(leaveoneout_df) > 0) {
      leaveoneout_df$Analysis <- "Leave-one-out"
      openxlsx::writeData(wb, "Leave-one-out", leaveoneout_df)
      if (any(c("p", "pval") %in% colnames(leaveoneout_df))) {
        col_p <- if ("pval" %in% colnames(leaveoneout_df)) "pval" else "p"
        pval_col <- which(colnames(leaveoneout_df) == col_p)
        openxlsx::conditionalFormatting(
          wb, "Leave-one-out",
          cols = pval_col, rows = 2:(nrow(leaveoneout_df) + 1),
          rule = ">0.05", style = greenStyle,
          type = "expression"
        )
      }
    } else {
      openxlsx::writeData(wb, "Leave-one-out", "No hay datos disponibles para Leave-one-out analysis")
    }
  } else {
    openxlsx::writeData(wb, "Leave-one-out", "No se encontraron resultados en la lista de Leave-one-out")
  }
  
  ## 3.1) Sig_Single
  if (length(flatten_single) > 0) {
    snp_all <- tryCatch(dplyr::bind_rows(flatten_single, .id = "Key"), error = function(e) NULL)
    if (!is.null(snp_all) && nrow(snp_all) > 0) {
      p_col <- if ("pval" %in% names(snp_all)) "pval" else if ("p" %in% names(snp_all)) "p" else NA_character_
      snp_col <- if ("SNP" %in% names(snp_all)) "SNP" else if ("snp" %in% names(snp_all)) "snp" else NA_character_
      if (!is.na(p_col)) {
        snp_sig <- dplyr::filter(snp_all, .data[[p_col]] <= 0.05)
        if (!is.na(snp_col)) snp_sig <- dplyr::filter(snp_sig, !grepl("^All\\b", .data[[snp_col]]))
        if (nrow(snp_sig) > 0) {
          openxlsx::addWorksheet(wb, "Sig_Single")
          openxlsx::writeData(wb, "Sig_Single", snp_sig)
          p_idx <- which(names(snp_sig) == p_col)
          openxlsx::conditionalFormatting(
            wb, "Sig_Single",
            cols = p_idx, rows = 2:(nrow(snp_sig) + 1),
            rule = "<=0.05", style = greenStyle,
            type = "expression"
          )
        }
      }
    }
  }
  
  ## 4.1) Nonsig_Leave: p > 0.05 o NA
  if (length(flatten_loo) > 0) {
    loo_all <- tryCatch(dplyr::bind_rows(flatten_loo, .id = "Key"), error = function(e) NULL)
    if (!is.null(loo_all) && nrow(loo_all) > 0) {
      p_col <- if ("pval" %in% names(loo_all)) "pval" else if ("p" %in% names(loo_all)) "p" else NA_character_
      snp_col <- if ("SNP" %in% names(loo_all)) "SNP" else if ("snp" %in% names(loo_all)) "snp" else NA_character_
      if (!is.na(p_col)) {
        loo_nonsig <- dplyr::filter(loo_all, is.na(.data[[p_col]]) | .data[[p_col]] > 0.05)
        if (!is.na(snp_col)) loo_nonsig <- dplyr::filter(loo_nonsig, !grepl("^All\\b", .data[[snp_col]]))
        if (nrow(loo_nonsig) > 0) {
          openxlsx::addWorksheet(wb, "Nonsig_Leave")
          openxlsx::writeData(wb, "Nonsig_Leave", loo_nonsig)
          p_idx <- which(names(loo_nonsig) == p_col)
          openxlsx::conditionalFormatting(
            wb, "Nonsig_Leave",
            cols = p_idx, rows = 2:(nrow(loo_nonsig) + 1),
            rule = ">0.05", style = greenStyle,
            type = "expression"
          )
        }
      }
    }
  }
  
  ## Guardar Excel evitando sobrescribir
  count <- 0
  final_name <- paste0(filename_prefix, ".xlsx")
  while (file.exists(final_name)) {
    final_name <- paste0(filename_prefix, "_", count, ".xlsx")
    count <- count + 1
  }
  openxlsx::saveWorkbook(wb, final_name, overwrite = FALSE)
  message("Archivo guardado en: ", normalizePath(final_name))
}

#Ejecutar el comando personalizado
save_sensitivity_results()
if (interactive()) browseURL(getwd())

# ###
# # ---------------
# # Gráficas
# # ---------------
# ###

library(ggplot2)
library(stringr)
library(viridis)

# Función para guardar sin sobrescribir
guardar_plot_sin_sobrescribir <- function(plot, base_filename, width = 7, height = 5, dpi = 300) {
  count <- 0
  final_filename <- paste0(base_filename, ".png")
  while (file.exists(final_filename)) {
    count <- count + 1
    final_filename <- paste0(base_filename, "_", count, ".png")
  }
  ggsave(
    filename = final_filename, plot = plot,
    width = width, height = height, dpi = dpi,
    bg = "white", limitsize = FALSE
  )
  message("Gráfico guardado en: ", normalizePath(final_filename))
}

# Tema MR reutilizable
mr_theme <- function(base_size = 12,
                     axis_text = 11,
                     axis_title = 12,
                     legend_text = 9,
                     legend_title = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      
      axis.text  = element_text(size = axis_text),
      axis.title = element_text(size = axis_title),
      
      legend.position = "bottom",
      legend.title = element_text(size = legend_title),
      legend.text  = element_text(size = legend_text),
      
      legend.key.size  = unit(0.6, "lines"),
      legend.spacing.x = unit(0.3, "lines"),
      
      plot.margin = margin(8, 12, 12, 12)
    )
}

# Guías para evitar leyendas extra (linetype/shape/fill) y controlar tamaño/filas
mr_guides <- function(legend_title = "Método", nrow = 2) {
  guides(
    colour = guide_legend(
      title = legend_title,
      nrow = nrow,
      byrow = TRUE,
      override.aes = list(size = 3, linewidth = 1.1)
    ),
    linetype = "none",
    shape    = "none",
    fill     = "none"
  )
}

# Opcional: partir etiquetas largas de la leyenda
wrap_legend_labels <- function(plot, width = 18) {
  plot + scale_colour_discrete(labels = function(x) str_wrap(x, width = width))
}

# ---------------------
# Scatter plots 
# Gráficas de dispersión (uno por dataset) 
# ---------------------
modified_scatter_plots <- list()

for (i in seq_along(names_list)) {
  this_exp <- names_list[i]
  
  if (!should_proceed(this_exp)) {
    message("Saltando scatter para ", this_exp,
            " (n_sig_methods = ", sig_method_counts[this_exp],
            "; umbral = ", min_sig, ").")
    next
  }
  
  res_name <- paste0("res", i-1, "_", s, "_", this_exp)
  dat_name <- paste0("dat", i-1, "_", this_exp)
  
  if (!exists(res_name) || !exists(dat_name)) next
  
  res_all <- get(res_name)
  dat_obj <- get(dat_name)
  
  # Filtrar solo métodos significativos
  res_sig <- res_all %>%
    dplyr::filter(!is.na(pval) & pval < 0.05)
  
  # Si no hay métodos significativos, no hacer plot
  if (nrow(res_sig) == 0) {
    message("No hay métodos significativos para el scatter de ", this_exp)
    next
  }
  
  scatter <- mr_scatter_plot(res_sig, dat_obj)[[1]]
  
  plot_mod <- scatter +
    scale_color_viridis_d(option = "D", begin = 0, end = 0.85) +
    mr_theme(base_size = 12, axis_text = 11, axis_title = 12, legend_text = 9, legend_title = 10) +
    mr_guides("Método", nrow = 2)
  
  modified_scatter_plots[[this_exp]] <- plot_mod
  
  guardar_plot_sin_sobrescribir(
    plot = plot_mod,
    base_filename = paste0("mr_scatterplot_sig_", this_exp, s),
    width = 7, height = 5, dpi = 300
  )
}

# ---------------------
# Gráfica Leave-one-out (uno por dataset) 
# ---------------------
#modified_leaveoneout_plots <- list()
#
#if (length(res_leaveoneout_list) > 0) {
#  for (dataset_name in names(res_leaveoneout_list)) {
#    
#    exposure_key <- sub("^dat\\d+_", "", dataset_name)
#    
#    if (!should_proceed(exposure_key)) {
#      message("Saltando Leave-one-out para ", exposure_key,
#              " (n_sig_methods = ", sig_method_counts[exposure_key],
#              "; umbral = ", min_sig, ").")
#      next
#    }
#    
#    obj <- res_leaveoneout_list[[dataset_name]]
#    
#    if (is.list(obj) && !is.data.frame(obj)) {
#      avail <- intersect(prefer_loo_safe, names(obj))
#      if (length(avail) == 0) {
#        message("Sin método LOO válido en ", dataset_name)
#        next
#      }
#      method_label <- avail[1]
#      res_leave <- obj[[method_label]]
#    } else {
#      method_label <- "mr_ivw_or_egger"
#      res_leave <- obj
#    }
#    
#    if (is.null(res_leave) || !is.data.frame(res_leave) || nrow(res_leave) == 0) next
#    leave <- mr_leaveoneout_plot(res_leave)[[1]]
#    
#    n_snps <- length(unique(res_leave$SNP))
#    
#    plot_mod <- leave +
#      scale_color_viridis_d(option = "D", begin = 0, end = 0.85) +
#      theme_minimal(base_size = 9) +
#      theme(
#        legend.position = "none",
#        
#        # SNPs un poco más grandes
#        axis.text.y  = element_text(size = 10),
#        
#        # resto más pequeño
#        axis.text.x  = element_text(size = 8),
#        axis.title.x = element_text(size = 8),
#        plot.title   = element_text(size = 10, face = "bold"),
#        
#        plot.margin  = margin(8, 18, 12, 18)
#      ) +
#      guides(colour = "none", linetype = "none", shape = "none", fill = "none") +
#      labs(colour = NULL, linetype = NULL, shape = NULL, fill = NULL) +
#      ggtitle(paste0("Leave-one-out - ", exposure_key, " [", method_label, "]"))
#    
#    # Ajustamos solo la ALTURA según número de SNPs (suave)
#    ggsave(
#      filename = paste0("mr_leaveoneout_", exposure_key, "_", method_label, "_", s, ".png"),
#      plot = plot_mod,
#      width = 12,
#      height = 0.14 * n_snps + 4,   # más espacio vertical real
#      dpi = 300,
#      bg = "white",
#      limitsize = FALSE
#    )
#    
#  }
#} else {
#  message("No se generaron leave-one-out plots.")
#}
#

#
# ---------------------
# Funnel plots (uno por dataset)
# ---------------------
modified_funnel_plots <- list()

if (length(res_singleSNP_list) > 0) {
  for (dataset_name in names(res_singleSNP_list)) {
    
    exposure_key <- sub("^dat\\d+_", "", dataset_name)
    
    if (!should_proceed(exposure_key)) {
      message("Saltando Funnel para ", exposure_key,
              " (n_sig_methods = ", sig_method_counts[exposure_key],
              "; umbral = ", min_sig, ").")
      next
    }
    
    res_funnel <- res_singleSNP_list[[dataset_name]]
    if (is.null(res_funnel) || !is.data.frame(res_funnel) || nrow(res_funnel) == 0) next
    
    funnel <- mr_funnel_plot(res_funnel)[[1]]
    
    plot_mod <- funnel +
      scale_color_viridis_d(option = "D", begin = 0, end = 0.85) +
      mr_theme(base_size = 12, axis_text = 11, axis_title = 12, legend_text = 9, legend_title = 10) +
      mr_guides("Método", nrow = 2) +
      ggtitle(paste0("Funnel plot - ", exposure_key))
    
    modified_funnel_plots[[dataset_name]] <- plot_mod
    
    guardar_plot_sin_sobrescribir(
      plot = plot_mod,
      base_filename = paste0("mr_funnel_", exposure_key, "_", s),
      width = 7, height = 5, dpi = 300
    )
  }
} else {
  message("No se generaron funnel plots.")
}

