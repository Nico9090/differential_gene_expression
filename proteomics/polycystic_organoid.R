#<----pkd organoid----->
source("~/bioinformatics_projects/scripts/functions.R")
#Load data
dis_normal_organ<-readr::read_csv("dis_normal_organ_ints.csv") #intensities from mass spec
dis_normal_organ_dep<-readr::read_csv("dis_normal_organ_diff_expr_prot.csv") #DEP
dis_normal_organ_gsea<-readRDS("dis_normal_organ_gsea.rds") #gsea analysis data
dis_normal_organ.go<-dis_normal_organ_gsea[["go_obj"]] #gsea object
#<--------------------------------------------------------->
#want to use B>0 and p.val<0.05 for high certainty of DEP
sig_dis_normal_organ_dep<-dis_normal_organ_dep%>%
  dplyr::filter(adj.P.Val<=0.05)%>%
  dplyr::filter(B>0)
#<----GO term mining---->
#list of go terms and ids
dis_normal_organ_go_terms<-tibble(GO_term=dis_normal_organ.go@result[["ID"]],
                          GO_desc=dis_normal_organ.go@result[["Description"]])
#select go id/term of interest
go_genes<-dis_normal_organ.go@geneSets[["GO:0033365"]]

#DEPs filtered by belonging to selected go terms
prot_to_organelle<-sig_dis_normal_organ_dep[sig_dis_normal_organ_dep[["protein"]]%in% 
                                  go_genes,]

stereocilium<-dis_normal_organ_dep[dis_normal_organ_dep[["protein"]]%in% 
                                     go_genes,]
N_acet<-sig_dis_normal_organ_dep[sig_dis_normal_organ_dep[["protein"]]%in% 
                                   go_genes,]


dis_normal_organ.go <- enrichplot::pairwise_termsim(dis_normal_organ.go)

#<---go terms of interest
#<--reconstructing gsea based on interested terms
cyto_skel_categories<-c("microtubule"," actin ",
                        "axoneme","cytoskeleton",
              "actin ","microvillus",
              "stereocilium",
              "transmembrane transport"
              )
interested_cyto_skel <- dis_normal_organ.go@result %>%
  filter(
    grepl(paste(cyto_skel_categories, collapse = "|"),
          Description,
          ignore.case = TRUE)
  )

interested_cyto_skel_gsea_go <- new("gseaResult",
                                    result = interested_cyto_skel)

interested_cyto_skel_gsea_go<-enrichplot::pairwise_termsim(interested_cyto_skel_gsea_go)

most_enriched_go<-dis_normal_organ.go@result[dis_normal_organ.go@result[["enrichmentScore"]]>=0.5,]

most_enriched_go<-most_enriched_go %>% dplyr::filter(setSize>=20)
