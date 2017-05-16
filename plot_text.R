plot.text <- function(t) {
    text = t
    ggplot() +
        annotate("text", x = 4, y = 25, size=8, label = text) +
        theme_bw() +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())
}
