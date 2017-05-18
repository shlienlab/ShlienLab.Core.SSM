filter_indel <- function(data=NULL, source='WXS', coverage=TRUE, indels=c('ins', 'del'), exac=0.001) {
  if (is.null(data)) stop("Mandatory argument data is missing")

  # only process PASS indel data (i.e. gatk_mutation_type is "ins" or "del")
  data <- data %>%
    filter(gatk_filter == 'PASS') %>%
    filter(gatk_mutation_type %in% indels)

  # we haven't done a detailed analysis of depth of coverage for indels, regardless of the library type
  # we will use the same coverage cutoffs until further investigation
  if (coverage == TRUE) {
    if (source == 'WXS') {
      data <- data %>%
        filter(gatk_normal_depth >= 10) %>%
        filter(gatk_tumour_depth >= 20)
    } else if (source == 'WGS') {
      data <- data %>%
        filter(gatk_normal_depth >= 10) %>%
        filter(gatk_tumour_depth >= 20)
    } else if (source == 'CPANEL') {
      data <- data %>%
        filter(gatk_normal_depth >= 10) %>%
        filter(gatk_tumour_depth >= 20)
    } else {
      stop("Invalid source, must be either WGS, WXS, or CPANEL")
    }
  }
  #data <- data %>% filter(gatk_normal_depth >= normal.depth & gatk_tumour_depth >= tumour.depth)

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

  # there is an exome specific filter that can be applied using the ExAC and ESP databases
  if ('annovar_exac' %in% names(data) & source == 'WXS') {
    data <- data %>%
      filter(annovar_exac <= exac | is.na(annovar_exac)) %>%
      filter(is.na(annovar_esp) | annovar_esp == '' | annovar_esp == 0 | grepl(x=annovar_clinvar, pattern='CLINSIG=pathogenic'))
  }

  return(data)
  }
