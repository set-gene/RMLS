h1 = oncoPrint(MLS_mut_mat ,show_pct = FALSE,right_annotation =NULL,
               left_annotation = rowAnnotation("VAF" = anno_barplot(top20_MLS_vaf$VAF,axis_param = list(direction = "reverse",gp=gpar(fontsize = 4),at = c(0, 0.2, 0.4))), "q-value" = anno_barplot(top20_MLS_vaf$q_value,axis_param = list(direction = "reverse",gp=gpar(fontsize = 4),at = c(0, 1, 2))),
                                               width = unit(1.5, "cm"),
                                               annotation_name_gp = gpar(fontsize=5)),
               alter_fun = alter_fun, col = col,
               column_labels = colnames(MLS_mut_mat),
               show_column_names = TRUE,
               remove_empty_columns = FALSE,
               remove_empty_rows = FALSE,top_annotation = NULL,
               heatmap_legend_param = heatmap_legend_param,
               column_names_gp=gpar(fontsize = 5),
               row_names_gp=gpar(fontsize = 5)) 

h2 =  oncoPrint(RMLS_mut_mat ,show_row_names = FALSE, show_pct = FALSE, right_annotation = rowAnnotation("q-value" = anno_barplot(top20_RMLS_vaf$q_value,axis_param=list(gp=gpar(fontsize = 4),at = c(0, 1, 2))),"VAF" = anno_barplot(top20_RMLS_vaf$VAF,axis_param=list(gp=gpar(fontsize = 4),at = c(0, 0.2, 0.4))),width = unit(1.5, "cm"),
                                                                                                         annotation_name_gp = gpar(fontsize=5)),
                alter_fun = alter_fun2, col = col2, column_labels = colnames(RMLS_mut_mat),
                show_column_names = TRUE,
                remove_empty_columns = FALSE, top_annotation = NULL,
                remove_empty_rows = FALSE ,row_labels = NULL,show_heatmap_legend = FALSE,column_names_gp=gpar(fontsize = 5))


ht_np=oncoPrint(npmat,get_type = function(x) strsplit(x, ";")[[1]], alter_fun = alter_fun, col = col, row_names_side = "right", show_column_names=TRUE, show_pct = FALSE, heatmap_legend_param=list(title ="Mutations", title_gp= gpar(fontsize = 6),labels_gp = gpar(fontsize = 4)), cluster_row_slices=FALSE, row_names_gp = gpar(fontsize = 0), pct_gp = gpar(fontsize = 2), column_order=columnorder_np, row_split = factor(c(rep("Significantly\nDifferent", 7), rep("dN/dScv", 9), rep("Pathway", 44),rep("Other Significantly\nDifferent", 66), rep("Other drivers", 47)), levels=c("Significantly\nDifferent", "dN/dScv", "Pathway", "Other Significantly\nDifferent", "Other drivers")), remove_empty_rows = TRUE, column_gap = unit(8, "mm"), row_title_gp= gpar(fontsize = 5), column_title_gp = gpar(fontsize = 5), column_names_gp=gpar(fontsize = 4), show_heatmap_legend = FALSE)
