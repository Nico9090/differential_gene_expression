#Gene set enrichment analysis
dc_nc_ws_gsea<-gseago_fn(dc_nc_ws_dep%>%
                           dplyr::filter(adj.P.Val<=0.05)%>%
                           dplyr::filter(B>0),
                         "protein",
                         "logFC",
                         org.Hs.eg.db)
dis_normal_organ_gsea<-gseago_fn(dis_normal_organ_dep%>%
                                   dplyr::filter(adj.P.Val<=0.05),
                                 "protein",
                                 "logFC",
                                 org.Hs.eg.db)
