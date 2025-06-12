library('ggplot2')
data=read.csv(file ="P2_genetable_collapsed.txt", sep='\t', header=TRUE, as.is=TRUE,stringsAsFactors=FALSE)

ranked_data=data[order(data[,14]),]


new=ranked_data[1:50,]
ggplot(new, aes(x=new[,14], y=-log10(new[,13]))) +
  geom_point()+ geom_text(
    data = new,
    aes(x = new[,14], y = -log10(new[,13]), label = new[,1]),
    hjust = 1,
    vjust = 2
  )

#load libraries

library('ggplot2')
library('ggrepel')
#import genetable_collapsed.txt from process_experiments.py as dataframe

genetable <- read.table('P2_genetable_collapsed.txt', header = FALSE, skip = 4, sep = '\t', na.strings = ' ')
genetable <- na.omit(genetable)

names(genetable) <- c('gene', 'transcript', 'pval', 'gamma', 'sgrna_count_pval', 'sgrna_count_gamma')

head(genetable)
#definte dataframes with hits and pseudogenes

is_pseudo <- data.frame(subset(genetable, grepl('pseudo', genetable$gene, fixed = TRUE)))

is_hit <- subset(genetable, ((genetable$gamma / sd(is_pseudo$gamma) * -1 * log(genetable$pval, 10) >= 6 |
                                genetable$gamma / sd(is_pseudo$gamma) * -1 * log(genetable$pval, 10) <= -6 &
                                !(genetable$gene %in% is_pseudo$gene))))


head(is_pseudo)
head(is_hit)
#define thresholds for plot

draw_gamma_threshold_neg <- function(x) {
  -6*sd(is_pseudo$gamma)/x
}

draw_gamma_threshold_pos <- function(x) {
  6*sd(is_pseudo$gamma)/x
}
#assign categories to genes for ggplot2 color and legend

genetable_color <- rbind(is_hit, is_pseudo)
genetable_color$cat <- ifelse(genetable_color$gene %in% is_hit$gene, 'gene_hit', 'negative_control')

head(genetable_color)
#plot gamma-value volcano plot with pseudo-genes and hits past empirically-defined threshold

p <- ggplot(genetable, aes(x = gamma, y = -1 * log(pval, 10), label = gene)) +
  geom_point(color = 'gray80', size = 1) +
  geom_point(data = genetable_color, aes(color = cat), size = 1) +
  theme_classic() +
  xlim(-6.0, 6.0) +
  ylim(0.0, 5) +
  xlab(expression('log'[2] * '(high AR / low AR expression)')) +
  ylab(expression('-log'[10] * '(Mann-Whitney p-value)')) +
  
  theme(axis.text.x = element_text(size = rel(2), color = 'black'),
        axis.text.y = element_text(size = rel(2), color = 'black'),
        axis.title.x = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        legend.title = element_blank(),
        legend.text = element_text(size = rel(1.5), color = 'black')
  ) +
  
  stat_function(fun = draw_gamma_threshold_neg, linetype = 'dashed', alpha = 0.4) +
  stat_function(fun = draw_gamma_threshold_pos, linetype = 'dashed', alpha = 0.4) +
  
  geom_text_repel(data = subset(is_hit, abs(is_hit$gamma / sd(is_pseudo$gamma) * -1 * log(is_hit$pval, 10)) > 22),
                  segment.alpha = 0.2, force = 8) +
  
  scale_color_manual(values = c('thistle', 'sandybrown'), labels = c('Hit', 'Negative control')) +
  guides(color = guide_legend(override.aes = list(size = 2)))

p

ggsave('volcano.pdf', width = 8, height = 6)
