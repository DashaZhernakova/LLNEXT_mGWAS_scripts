args <- commandArgs(trailingOnly = TRUE)

library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)

infile = args[1]
outfile = args[2]
d <- read.delim(infile, header = F, as.is = T, sep = "\t", check.names = F)
colnames(d) <- c("timepoint", "bac", "SNPChr", "SNP", "SNPChrPos", "EA", "OA", "N", "MAF", "beta", "se", "P")

d$SNPChrPos <- as.numeric(d$SNPChrPos)

for (tp in unique(d$timepoint))

don <- d %>% 
  
  # Compute chromosome size
  group_by(SNPChr) %>% 
  summarise(chr_len=max(SNPChrPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(d, ., by=c("SNPChr"="SNPChr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(SNPChr, SNPChrPos) %>%
  mutate( BPcum=SNPChrPos+tot)

# Add highlight 
#mutate( is_highlight=ifelse(SNP %in% snps_to_highlight, "yes", "no") %>%


axisdf <- don %>% group_by(SNPChr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png(outfile, width = 15, height = 5, units = 'in', res = 400)
ggplot(don, aes(x=BPcum, y=-log10(PValue))) +
  # Show all points
  geom_point( aes(color=as.factor(SNPChr)),  size=1) +
  scale_color_manual(values = rep(c("#495DA0", "#74AFDF"), 22 )) +
  geom_hline(yintercept=-1*log(1.872659e-10,10), color = "#EF3B2C") +
  geom_hline(yintercept=-1*log(5e-08,10), color = "orange") +
  # custom X axis:
  xlab("Chromosome") +
  scale_x_continuous( label = axisdf$SNPChr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) , limits = c(0,45)) +
  
  # Add annotation
  #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=ProbeName), size=2, max.overlaps = 20) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.y = element_line(color="lightgrey", size = 0.5),
    axis.ticks.x = element_blank()
    
  )
dev.off()