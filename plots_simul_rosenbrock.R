library(ggplot2)
library(patchwork)

#xfull<-x
x<-xfull[1:100,]
x<-xfull

df <- as.data.frame(x)
colnames(df) <- paste0("x", 1:ncol(df))
expr_labels <- list(expression(x[1]), expression(x[2]), expression(x[3]))
plot_grid <- matrix(list(), nrow = d, ncol = d)

legend_shown <- FALSE

for (i in 1:d) {
  for (j in 1:d) {
    p <- NULL
    xvar <- paste0("x", j)
    yvar <- paste0("x", i)
    
    if (i == j) {
      # Diagonal: marginal density
      p <- ggplot(df, aes(x = .data[[xvar]])) +
        geom_density(fill = "white", color = "black") +
        labs(x = expr_labels[[j]]) +
        theme_bw(base_size = 10)
    } else if (i > j) {
      # Lower triangle: scatterplot
      p <- ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
        geom_point(color = "black", alpha = 0.5, size = 0.5) +
        labs(x = expr_labels[[j]], y = expr_labels[[i]]) +
        theme_bw(base_size = 10)
    } else {
      # Upper triangle: log-scaled contour
      show_legend <- legend_shown
      p <- ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
        stat_density_2d(aes(color = after_stat(level)), contour = TRUE, bins=13, size = 0.3) +
        scale_color_gradient(low = "gray70", high = "black", trans="log", name = "log level") +
        labs(x = expr_labels[[j]], y = expr_labels[[i]]) +
        theme_bw(base_size = 10)
      
      if (!show_legend) {
        p <- p + theme(legend.position = "none")
      } else {
        legend_shown <- TRUE
      }
    }
    
    plot_grid[[i, j]] <- p
  }
}

final_plot <- wrap_plots(t(plot_grid), nrow = d, ncol = d, guides = "collect") 
final_plot

ggsave("density_ex2_bw.png", width = 10, height = 8, dpi = 300, units = "in")
