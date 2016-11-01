library(ggplot2)
library(reshape2)

file_list <- list.files(pattern = "*DoC$")

for (file in file_list)
{
  png(file = paste(file, 'png', sep = ".") , width = 1200 , height = 800)
  x <- read.table(file, header = TRUE , sep = "\t")
  df <- data.frame(x$Locus , x$Total_Depth)
  names(df)[names(df) == "x.Total_Depth"] <- file
  mdf <- melt(df)
  plot <- ggplot(mdf , aes(x = x.Locus , y = value , colour = value < 20))
  plot <-  plot + geom_point() + ggtitle(file)
  plot <-  plot + scale_colour_manual (name = 'Coverage' , values = setNames(c('red','green'),c(T,F)))
  plot <- plot + xlab("") + ylab ("Coverage")
  plot(plot)
  dev.off()
}
