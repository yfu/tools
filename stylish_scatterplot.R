require(grid)
require(ggplot2)
require(scales)
require(ggthemes)

stylish_scatterplot <- function(df, cn1, cn2, col="group", name1, name2, main="my.title", lim) {
    ## cn1 and cn2: colname1 and colname2
    sample <- df
    gg=ggplot( sample, aes_string(x = cn1, y = cn2, color=col) ) +
        theme (
            plot.margin=unit(c(1,1,0,0),"lines"),
            axis.text=element_text (size=4),
            axis.title=element_text(size=6),
            legend.margin=unit(0,"lines"),
            panel.margin=unit(0, "lines"),
            axis.ticks.margin=unit(0,"lines"),
            legend.key.size=unit(0.5,"lines")
            ) +
                scale_x_log10 ( limits = c(1,lim), breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(10^.x) ) ) +
                                   scale_y_log10 ( limits = c(1,lim), breaks = trans_breaks("log10", function(x) 10^x),
                                                  labels = trans_format("log10", math_format(10^.x) ) ) +
                                                      annotation_logticks () +
                                                          scale_size_manual( values=c(0,8) ) +
                                                              geom_abline (intercept=0, slope=1, colour="darkgrey", linetype='dashed') +
                                                                  geom_abline (intercept=log10(2), slope=1, colour="darkgrey", linetype='twodash') +
                                                                      geom_abline (intercept=-log10(2), slope=1, colour="darkgrey", linetype='twodash') +

                                                                          theme_few () +
                                                                              scale_fill_continuous(guide = "legend") +
                                                                                  geom_point(size=5, alpha=0.75, na.rm=T) +
                                                                                      scale_colour_manual(values=c("black", "lightblue" ,"darkgreen", "red", "blue", "orange")) +
                                                                                          guides(colour = guide_legend (title=expression (paste ("type")), title.position = "top")) +
                                                                                              xlab ( substitute ( paste(italic(name1)), list(name1=name1, name2=name2))) +
                                                                                                  ylab ( substitute ( paste(italic(name2)), list(name1=name1, name2=name2))) +
                                                                                                      coord_fixed(ratio=1, xlim = c(1, lim), ylim = c(1, lim)) + 
                                                                                                          ggtitle(main)
    gg
}
