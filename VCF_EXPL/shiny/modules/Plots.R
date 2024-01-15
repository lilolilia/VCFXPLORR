### PLOTS FUNCTIONS ###
library(ggplot2)

plot_histogram <- function(data, column) {
  ggplot(data, aes(x = .data[[column]])) +
    geom_histogram() +
    labs(x = column, y = "Counts") +
    scale_x_log10() +
    theme_classic()
}

Bar_plot1 <- function(data) {
  ggplot(data, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat = 'identity') +
    labs(x = "", y = "Counts", fill = "") +
    theme_classic()
}

Bar_plot2 <- function(data, column) {
  # Create a bar plot for the specified column
  ggplot(data, aes(x = .data[[column]], fill = .data[[column]])) +
    geom_bar(position = "dodge", color = "black") +
    labs(title = paste("Genotype Counts -", column),
         x = "Genotypes",
         y = "Count") +
    theme_classic()
}

Bar_plot3 <- function(data) {
  ggplot(data, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat = 'identity') +
    labs(x = "", y = "Counts", fill = "") +
    theme(legend.position = "none")
}

Bar_plot4 <- function(data) {
  ggplot(data, aes(x = APOBEC_target, y = Freq, fill = APOBEC_target)) +
    geom_bar(stat = 'identity') +
    labs(x = "", y = "Variants", fill = "") +
    theme(legend.position = "none")
}

Bar_plot5<- function(data) {
  ggplot(data, aes(x = Var1, y = Freq, fill = Var1)) +
    geom_bar(stat = 'identity') +
    labs(x = "", y = "Counts", fill = "in_dbSNP") +
    theme(legend.position = "none") +
    labs(title = "Variations in dbSNP")
}

Bar_plot6 <- function(data) {
  ggplot(data, aes(x = reorder(hgnc_symbol, AF), y = AF, fill = hgnc_symbol)) +
  geom_bar(stat = "identity") +
  labs(title = "Allele Frequency Bar Plot by hgnc_symbol", x = "hgnc_symbol", 
       y = "Allele Frequency") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_lolliplot <- function(data, x_col, y_col, gene_col, pch = 20, col = "red", cex = 2, 
                           xlab = "Genomic Position", ylab = "Allele Frequency", 
                           main = "Lolliplot Based on Allele Frequency") {
  # Convert columns to numeric if not already
  x_values <- as.numeric(data[[x_col]])
  y_values <- as.numeric(data[[y_col]])
  
  # Scatter plot
  plot(
    x_values,
    y_values,
    pch = pch,
    col = col,
    cex = cex,
    xlab = xlab,
    ylab = ylab,
    main = main
  )
  
  # Add lollipops
  segments(
    x0 = x_values,
    y0 = rep(0, nrow(data)),  # Replicate 0 for the length of the data
    x1 = x_values,
    y1 = y_values,
    col = "blue"
  )
  
  # Add gene names vertically below each data point using text
  text(
    x = x_values,
    y = par("usr")[3] - 0.05,  # Adjust the position below the x-axis
    labels = data[[gene_col]],
    pos = 1,        # Position the text below the point
    cex = 0.8,      # Adjust the size of the text
    xpd = NA,       # Allow plotting outside the plotting region
    srt = 30       # Rotate text vertically
  )
}
