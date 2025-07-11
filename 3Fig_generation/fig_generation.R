packages <- c("tidyverse", "data.table", "sf", "viridis")
lapply(packages, suppressPackageStartupMessages(require), character.only=T)

setwd("C:/Users/gg9824/Dropbox/00ESRC Fellowship/Projects/UKB Geog Var/rev_results/")

exp_list <- c("Scz", "Dep", "Well")
out_list <- c("Greenness", "Greenspace") 

dodge <- position_dodge(width=0.5)

for(x in exp_list){
  for(y in out_list){
      df_raw <- fread(paste0(x, "/", y, "_", x, "_PRS_results_lm.csv")) %>% select(N, score, beta, Upper_CI, Lower_CI)
      
      df <- fread(paste0(x, "/", y, "_", x, "_PRS_MSOA_results_Mundlak_risk_score.csv")) %>% 
      select(N, score, beta, Upper_CI, Lower_CI)
      
      df_mean <- fread(paste0(x, "/", y, "_", x, "_PRS_MSOA_results_Mundlak_risk_score_mean.csv")) %>% 
      select(N, score, beta, Upper_CI, Lower_CI)
      
      df_raw$level <- "single-level"
      df_mean$level <- "contextual"
      df$level <- "within-area"
      
      df_master <- rbind(df_raw,df_mean, df)
      df_master$score <- factor(df_master$score, levels=c("polygenic_score_5e8", "polygenic_score_000001", "polygenic_score_00001", 
                                                          "polygenic_score_0001", "polygenic_score_001", "polygenic_score_01", 
                                                          "polygenic_score_05", "polygenic_score_1", "polygenic_score_2", 
                                                          "polygenic_score_3", "polygenic_score_4", "polygenic_score_5"))
      
      resplot <- df_master %>% ggplot(
        aes(x=beta, y=score, colour=level)) +
        geom_point(position=dodge)+ 
        geom_errorbar(aes(xmax=Upper_CI, xmin=Lower_CI), position=dodge, width=0.5)+
        geom_vline(xintercept = 0, linetype="dotted") +
        ggtitle(paste0(x, "_", y))
      
      nm <- paste0(x, "-", y, "_df") 
      assign(nm, df_master)
      
      ggsave(paste0(x, "_", y, ".png"), resplot)
  }
}

`Dep-Greenness_df`$expo <- "Depression"
`Scz-Greenness_df`$expo <- "Schizophrenia"
`Well-Greenness_df`$expo <- "Wellbeing"

### Bind all three exposures together for Greenness

greenness_master <- rbind(`Dep-Greenness_df`, `Scz-Greenness_df`, `Well-Greenness_df`)
greenness_master
greenness_master$level <- factor(greenness_master$level, levels=c("single-level", "within-area", "contextual"))

### 

levels(greenness_master$score)  <- gsub("polygenic_score", "PGI", levels(greenness_master$score))


### Gen manual colour palette from Viridis (omit bright yellow)

dodge <- position_dodge(width = 0.6)

graphic_pal <- c("#E69F00", "#56B4E9", "#009E73")

greenness_plot <- ggplot(greenness_master %>% arrange(level) %>% filter(level!="between"),
                         aes(x= beta,
                             y= score,
                             colour= level)) +
  scale_colour_manual(values=graphic_pal)+
  geom_point (size=1.5, position = dodge) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(size = 8, angle = 30, hjust = 0.9, vjust = 0.2),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  xlab("Beta")+
  ylab("PGI Threshold") +
  facet_grid(cols=vars(expo)) +
  xlim(-0.08, 0.08)+
  geom_errorbar(aes(xmax=Upper_CI, xmin=Lower_CI), position=dodge, width=0.5)

greenness_plot

ggsave("greenness_all.png", greenness_plot, dpi=500, width = 10, height = 6)

######## Exclude contextual

greenness_plot_nocontext <- ggplot(greenness_master %>% filter(level!= "contextual" & level!= "between") %>%
                           arrange(level),
                         aes(x= beta,
                             y= score,
                             colour= level)) +
  scale_colour_manual(values=graphic_pal)+
  geom_point (size=1.5, position = dodge) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(size = 8, angle = 30, hjust = 0.9, vjust = 0.2),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  xlab("Beta")+
  ylab("PGI Threshold") +
  facet_grid(cols=vars(expo)) +
  xlim(-0.02, 0.02)+
  geom_errorbar(aes(xmax=Upper_CI, xmin=Lower_CI), position=dodge, width=0.5)

greenness_plot_nocontext

ggsave("greenness_nocontext.png", greenness_plot_nocontext, dpi=500, width = 8, height = 6)
ggsave("greenness_nocontext_wide.png", greenness_plot_nocontext, dpi=500, width = 10, height = 6)


######## Repeat Fig for Greenspace outcome

`Dep-Greenspace_df`$expo <- "Depression"
`Scz-Greenspace_df`$expo <- "Schizophrenia"
`Well-Greenspace_df`$expo <- "Wellbeing"

### Bind all three exposures together for Greenspace

greenspace_master <- rbind(`Dep-Greenspace_df`, `Scz-Greenspace_df`, `Well-Greenspace_df`)


greenspace_master$level <- factor(greenspace_master$level, levels=c("single-level", "within-area", "contextual"))

#### Replace long pRS labels
levels(greenspace_master$score)  <- gsub("polygenic_score", "PGI", levels(greenspace_master$score))


greenspace_plot <- ggplot(greenspace_master %>% arrange(level) %>% filter(level!="between"),
                         aes(x= beta,
                             y= score,
                             colour= level)) +
  scale_colour_manual(values=graphic_pal)+
  geom_point (size=1.5, position = dodge) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(size = 8, angle = 30, hjust = 0.9, vjust = 0.2),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  xlab("Beta")+
  ylab("PGI Threshold") +
  facet_grid(cols=vars(expo)) +
  xlim(-0.08, 0.08)+
  geom_errorbar(aes(xmax=Upper_CI, xmin=Lower_CI), position=dodge, width=0.5)

greenspace_plot

ggsave("greenspace_all.png", greenspace_plot, dpi=500, width = 10, height = 6)

#### greenspace exclude contextual

######## Exclude contextual

greenspace_plot_nocontext <- ggplot(greenspace_master %>% filter(level!= "contextual" & level != "between") %>%
                                 arrange(level),
                               aes(x= beta,
                                   y= score,
                                   colour= level)) +
  scale_colour_manual(values=graphic_pal)+
  geom_point (size=1.5, position = dodge) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5) +
  theme_bw()+
  theme(axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(size = 8, angle = 30, hjust = 0.9, vjust = 0.2),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 10),
        axis.title = element_text(size = 10),
        legend.position = "bottom",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
  xlab("Beta")+
  ylab("PGI Threshold") +
  facet_grid(cols=vars(expo)) +
  xlim(-0.012, 0.012)+
  geom_errorbar(aes(xmax=Upper_CI, xmin=Lower_CI), position=dodge, width=0.5)


greenspace_plot_nocontext

ggsave("greenspace_nocontext.png", greenspace_plot_nocontext, dpi=500, width = 8, height = 6)
ggsave("greenspace_nocontext_wide.png", greenspace_plot_nocontext, dpi=500, width = 10, height = 6)

