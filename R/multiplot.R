#' Multiplot function from R cookbook
#'
#' @param ...       comma separated list of plots to include
#' @param plotlist  list, naming plots to include
#' @param cols      number of columns for grid of plots
#' 
#' @export
#' 
multiplot <- function(..., plotlist=NULL, cols) {
    #require(grid)
    # function from R cookbook: 
    ##http://wiki.stdout.org/rcookbook/Graphs/Multiple%20
    ##          graphs%20on%20one%20page%20(ggplot2)/
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # Make the panel
    plotCols = cols # Number of columns of plots
    plotRows = ceiling(numPlots/plotCols) # Number of rows needed
    
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
    vplayout <- function(x, y)
        viewport(layout.pos.row = x, layout.pos.col = y)
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
        curRow = ceiling(i/plotCols)
        curCol = (i-1) %% plotCols + 1
        print(plots[[i]], vp = vplayout(curRow, curCol ))
    }
    
}