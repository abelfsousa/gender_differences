suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
options(bitmapType = "cairo")
set.seed(123)


data1 <- data.frame(x = c(rep("XX", 10), rep("YY", 10)), y = c(c(1:10), c(11:20))) %>%
  as.tibble()

plot1 <- ggplot(data = data1, mapping=aes(x=x, y=y, fill=x)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size=15),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color="black", size=12),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    aspect.ratio=1) +
  scale_fill_brewer(type = "seq", palette = "Greys") +
  scale_x_discrete(labels=c("Tumour", "Normal")) +
  labs(x = "", y = "Gene expression")
ggsave(filename="boxplot1.png", plot=plot1, height=2, width=2, path = "./plots/other_plots/")
unlink("boxplot1.png")


plot2 <- ggplot(data = data1, mapping=aes(x=x, y=y, fill=x)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.title = element_text(size=15),
    axis.text.y = element_blank(),
    axis.text.x = element_text(color="black", size=12),
    axis.ticks.y = element_blank(),
    legend.position = "none",
    aspect.ratio=1) +
  scale_fill_brewer(type = "seq", palette = "Greys") +
  scale_x_discrete(labels=c("Male", "Female")) +
  labs(x = "", y = "Gene expression")
ggsave(filename="boxplot2.png", plot=plot2, height=2, width=2, path = "./plots/other_plots/")
unlink("boxplot2.png")
