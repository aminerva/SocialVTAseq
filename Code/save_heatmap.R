save_pheatmap_pdf <- function(x, filename, width=10, height=10) {
   
   # Function for saving the heatmap from pheatmap function

   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}