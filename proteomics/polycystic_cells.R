source("~/bioinformatics_projects/scripts/functions.R")
dis_normal_cell_dc_nc_ws<-readr::read_csv("dc9_nc1_protein_ints.csv")

#<---differential protein expression---->
dc_nc_ws_dep<-readr::read_csv("dc9_nc1_diff_expr_prot.csv")

dc_nc_ws_dep[["sequence"]]<-proteins$`Peptide sequences`

dc_nc_ws_gsea<-readRDS("dc_nc_ws_gsea.rds")

dc_nc_ws.go<-dc_nc_ws_gsea[["go_obj"]]

sig_dc_nc_ws_dep<-dc_nc_ws_dep%>%
  dplyr::filter(adj.P.Val<=0.05)%>%
  dplyr::filter(B>0)
#<----GO term mining---->
dc_nc_ws_go_terms<-tibble(GO_term=dc_nc_ws.go@result[["ID"]],
                          GO_desc=dc_nc_ws.go@result[["Description"]]
)
go_genes<-dc_nc_ws.go@geneSets[["GO:0048660"]]

androgen_metab<-sig_dc_nc_ws_dep[sig_dc_nc_ws_dep[["protein"]]%in% 
                                   go_genes,]

most_enriched_go<-dc_nc_ws.go@result[dc_nc_ws.go@result[["enrichmentScore"]]>=0.5,]

pos.reg.homotypic.celladh<-sig_dc_nc_ws_dep[sig_dc_nc_ws_dep[["protein"]]%in% 
                                              go_genes,]
endocardium_development<-sig_dc_nc_ws_dep[sig_dc_nc_ws_dep[["protein"]]%in% 
                                            go_genes,]
smooth_muscle_cell_prolif<-sig_dc_nc_ws_dep[sig_dc_nc_ws_dep[["protein"]]%in% 
                                              go_genes,]
