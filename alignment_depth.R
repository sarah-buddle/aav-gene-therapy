library(ggplot2)
library(tidyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

input <- args[1]
output <- args[2]

depth <- read.delim(input, sep = "\t", header = FALSE, row.names = NULL, col.names = c("chromosome", "position", "depth"))

complete_sequence <- data.frame(position = seq(1, max(depth$position)))

chromosome <- unique(depth$chromosome)

depth_clean <- expand_grid(chromosome, complete_sequence) %>%
  left_join(depth, by = c("chromosome", "position")) %>%
  mutate(depth = ifelse(is.na(depth), 0, depth))

ggplot(depth_clean, aes(x = position, y = depth)) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank()) +
  ylab("Read depth") +
  scale_x_continuous(expand = c(0,0.1)) +
  scale_y_continuous(expand = c(0,0.1)) +
  geom_area(fill = "lightblue") +
  facet_grid(rows = vars(chromosome))

ggsave(output, height = 1.5, width = 6)
