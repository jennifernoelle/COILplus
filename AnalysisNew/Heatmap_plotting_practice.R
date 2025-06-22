superheat(mtcars,
          # scale the matrix columns
          scale = TRUE,
          # add row dendrogram
          row.dendrogram = TRUE)
set.seed(2016113)

superheat(mtcars,
          # scale the matrix columns
          scale = TRUE,
          # generate three column clusters
          n.clusters.rows = 3)

gears <- paste(mtcars$gear, "gears")

set.seed(2016113)
superheat(mtcars,
          # scale the matrix columns
          scale = TRUE,
          # cluster by gears
          membership.rows = gears)

grid.lines()
# Using id (NOTE: locations are not in consecutive blocks)
superheat(mtcars,
          # scale the matrix columns
          scale = TRUE,
          # cluster by gears
          membership.rows = gears)
grid.lines(x = c(1,3,5), y = c(1,2), lwd = 5)

grid.polyline(x=c((0:4)/10, rep(.5, 5), (10:6)/10, rep(.5, 5)),
              y=c(rep(.5, 5), (10:6/10), rep(.5, 5), (0:4)/10),
              id=rep(1:5, 4),
              gp=gpar(col=1:5, lwd=3))


# row_membership <- data.frame(vfamily = v.taxa.ordered[keep_v_index, 4])
# rownames(row_membership) <- v.taxa.ordered[keep_v_index, 1]
# col_membership <- data.frame(pfamily = p.taxa[keep_p_index, 3])
# rownames(col_membership) <- p.taxa[keep_p_index, 1]
# 
# grid.lines()
# pheatmap(plot_pred[keep_v_index, keep_p_index],
#          annotation_row = row_membership,
#          annotation_col = col_membership,
#          cluster_rows = FALSE,
#          cluster_cols = FALSE,
#          show_rownames = FALSE,
#          show_colnames = FALSE,
#          color = colorRampPalette(c("white", "darkgrey"))(100), 
#          na_col = "black", 
#          fontfamily = "serif"
#          )
# 
# # Where do the families change?
# vfamilies <- unique(row_membership$vfamily)
# grds_v <- sapply(vfamilies, function(x) min(which(row_membership$vfamily == x)))
# pfamilies <- unique(col_membership$pfamily)
# grds_p <- sapply(pfamilies, function(x) min(which(col_membership$pfamily == x)))
# 
# nr <- length(keep_v_index)
# nc <- length(keep_p_index)
# for (k_v in grds_v) {
#   for(k_p in grds_p){
#   grid.lines(x=c(0,1), y=k_p/nc, gp=gpar(col="black", lwd=2))
#   grid.lines(x=k_v/nr, y=c(0,1), gp=gpar(col="black", lwd=2))
# }}



