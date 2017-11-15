plot.text <- function(t, size=8) {
    text = t
    ggplot() +
        annotate("text", x = 4, y = 25, size=size, label = text) +
        theme_bw() +
        theme(panel.grid.major=element_blank(),
              panel.grid.minor=element_blank())
}
