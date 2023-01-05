

def ActivePathways(scores, gmt, background = makeBackground(gmt),
                            geneset_filter = [5, 1000], cutoff = 0.1, significant = 0.05,
                            merge_method = c("Brown", "Fisher", "Stouffer","Strube"),
                            correction_method = c("holm", "fdr", "hochberg", "hommel",
                                                  "bonferroni", "BH", "BY", "none"),
                            cytoscape_file_tag = Null, color_palette = NULL, custom_colors = NULL, 
                            color_integrated_only = "#FFFFF0", scores_direction = NULL, 
                            expected_direction = NULL):

    return

