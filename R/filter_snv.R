filter_snv <- function(data=NULL, coverage=TRUE, source='WXS', exac=0.001, vaf=0.0) {
  if (is.null(data)) stop("Mandatory argument data is missing")

  # check for valid entries for the source argument
  if (!(source %in% c('WGS', 'WXS', 'CPANEL')) & coverage == TRUE) stop("Invalid source argument, must be one of WGS, WXS, or CPANEL")

  # use MuTect's internal filter and keep those listed as "PASS"
  data <- data %>%
    filter(gatk_filter == 'PASS') %>%
    filter(gatk_mutation_type == 'snv')

  if (coverage == TRUE) {
    if (source == 'WXS') {
      data <- data %>%
        filter(gatk_normal_depth >= 20) %>%
        filter(gatk_tumour_depth >= 30)
    } else if (source == 'WGS') {
      data <- data %>%
        filter(gatk_normal_depth >= 10) %>%
        filter(gatk_tumour_depth >= 10)
    } else if (source == 'CPANEL') {
      data <- data %>%
        filter(gatk_normal_depth >= 50) %>%
        filter(gatk_tumour_depth >= 50)
    } else {
      stop("Invalid source, must be either WGS, WXS, or CPANEL")
      }
  }

  # check if clinvar was used to filter the data (older data sets may not have this)
  if ('annovar_clinvar' %in% names(data)) {
    data <- data %>%
      filter(is.na(annovar_dbsnp) | annovar_dbsnp == '' | annovar_dbsnp == 0 | grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic')) %>%
      filter(is.na(annovar_complete_genomics) | annovar_complete_genomics != '' | annovar_complete_genomics == 0 | grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic')) %>%
      filter(is.na(annovar_1000g) | annovar_1000g == '' | annovar_1000g == 0 | grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic'))
  } else {
    data <- data %>%
      filter(dbsnp_site != 'DBSNP') %>%
      filter(is.na(annovar_dbsnp) | annovar_dbsnp == '' | annovar_dbsnp == 0) %>%
      filter(is.na(annovar_complete_genomics) | annovar_complete_genomics == '' | annovar_complete_genomics == 0) %>%
      filter(is.na(annovar_1000g) | annovar_1000g == '' | annovar_1000g == 0)
    }

  # filter variants based on their variant allele fraction
  data <- data %>% filter(gatk_tumour_allele_fraction >= vaf)

  # there is an exome specific filter that can be applied using the ExAC and ESP databases
  if ('annovar_exac' %in% names(data) & source == 'WXS') {
    data <- data %>%
      filter(annovar_exac <= exac | is.na(annovar_exac)) %>%
      filter(is.na(annovar_esp) | annovar_esp == '' | annovar_esp == 0 | grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic'))
    }
  return(data)
  }
