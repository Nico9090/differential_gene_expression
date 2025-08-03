#pre-preparation 
raw_proteins<-readr::read_delim("proteinGroups.txt")
proteins<-raw_proteins %>%
  #<--1. don't want contaminated proteins
  dplyr::filter(is.na(`Potential contaminant`),
                is.na(`Reverse`)) %>%
  dplyr::select(`Protein IDs`,`Majority protein IDs`,
                `Peptide counts (razor+unique)`,`Gene names`,
                `Peptide sequences`,
                dplyr::starts_with("Intensity"))
#<---subsetting some replicates
dis_normal_cell_dc_nc_ws<-proteins%>%
  dplyr::select(`Protein IDs`,`Majority protein IDs`,
                `Peptide counts (razor+unique)`,
                dplyr::starts_with("Intensity 37"),
                dplyr::starts_with("Intensity 38"),
                dplyr::starts_with("Intensity 39"),
                dplyr::starts_with("Intensity 16"),
                dplyr::starts_with("Intensity 17"),
                dplyr::starts_with("Intensity 18"))

#<---limma analysis---->
intensities<-dis_normal_cell_dc_nc_ws %>%
  select(dplyr::starts_with("Intensity"))
intensities<-base::log2(intensities)
meta_data<-dplyr::tibble(samples=colnames(intensities),
                         condition=c(rep("normal",9),rep("disease",9))
)
meta_data[["condition"]]<-factor(meta_data[["condition"]],
                                 levels=c("normal","disease"))
design<-stats::model.matrix(~0 + meta_data[["condition"]]) 

dc_nc_ws_dep<-limma_fit_for_deg(counts = intensities,
                                design = design,
                                meta_data = meta_data,
                                sample_identifier_by = "condition"
)
dc_nc_ws_dep[["protein"]]<-dis_normal_cell_dc_nc_ws[["gene"]]

#<-----------actual organoid tissue, not cell culture----->
dis_normal_organ<-proteins%>%
  dplyr::select(`Protein IDs`,`Majority protein IDs`,
                `Peptide counts (razor+unique)`,
                `Peptide sequences`,
                dplyr::starts_with("Intensity 67"),
                dplyr::starts_with("Intensity 68"),
                dplyr::starts_with("Intensity 69"),
                dplyr::starts_with("Intensity 70"),
                dplyr::starts_with("Intensity 71"),
                dplyr::starts_with("Intensity 72"),
                dplyr::starts_with("Intensity 58"),
                dplyr::starts_with("Intensity 59"),
                dplyr::starts_with("Intensity 60"),
                dplyr::starts_with("Intensity 61"),
                dplyr::starts_with("Intensity 62"),
                dplyr::starts_with("Intensity 63"),
                dplyr::starts_with("Intensity 64"),
                dplyr::starts_with("Intensity 65"),
                dplyr::starts_with("Intensity 66")
  )
intensities_organ<-dis_normal_organ %>%
  select(dplyr::starts_with("Intensity"))
intensities_organ<-base::log2(intensities_organ+1)
meta_data2<-dplyr::tibble(samples=colnames(intensities_organ),
                          condition=c(rep("normal",18),rep("disease",27))
)
meta_data2[["condition"]]<-factor(meta_data2[["condition"]],
                                  levels=c("normal","disease"))
design2<-stats::model.matrix(~0 + meta_data2[["condition"]]) 
#<---limma analysis again
dis_normal_organ_dep<-limma_fit_for_deg(counts = intensities_organ,
                                        design = design2,
                                        meta_data = meta_data2,
                                        sample_identifier_by = "condition"
)

dis_normal_organ_dep[["protein"]]<-proteins$`Gene names`
