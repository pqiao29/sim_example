library(ggplot2)
RDR_var_pool <- seq(0.2, 3, 0.4)
K_range <- 3:7

df_I <- df_H <- data.frame()
for(RDR_var in RDR_var_pool){
  
  for(K in K_range){
    
    res_fin <- readRDS(paste0("result/RDRvar", RDR_var*10, "_K", K, ".rds"))
    res_fin <- res_fin[unlist(lapply(res_fin, function(x) length(x) > 1))]
    N <- length(res_fin)
    
    I_hclust <- unlist(lapply(res_fin, function(x) x[3, 1]))
    I_RDR <- unlist(lapply(res_fin, function(x) x[1, 1]))
    I_both <- unlist(lapply(res_fin, function(x) x[2, 1]))
    H_RDR <- unlist(lapply(res_fin, function(x) x[1, 2]))
    H_both <- unlist(lapply(res_fin, function(x) x[2, 2]))
    
    I_accuracy <- c(I_hclust, I_RDR, I_both)
    H_accuracy <- c(H_RDR, H_both)
    
    methods_I <- c(rep(paste0("hclust", K), N), 
                   rep(paste0("RDR", K), N), 
                   rep(paste0("both", K), N))
    methods_H <- c(rep(paste0("RDR", K), N), rep(paste0("both", K), N))
    
    df_I <- rbind(df_I, data.frame("RDR_var" = RDR_var, "accuracy" = I_accuracy, "Method" = methods_I))
    methods_H <- c(rep(paste0("RDR", K), N), rep(paste0("both", K), N))
  }
  
}

##### Plot clustering
method_order <- paste0(rep(c("hclust", "RDR", "both"), each = length(K_range)), rep(K_range, 3))
method_labels <- paste0(rep(c("Hierarchical (K = ", "RDR (K = ", "RDR + BAF (K = "), each = length(K_range)), rep(K_range, 3), ")")
df_I$RDR_var <- factor(df_I$RDR_var, levels = RDR_var_pool, labels = 1:length(RDR_var_pool))
df_I$Method <- factor(df_I$Method, levels = method_order, labels = method_labels)
method_colors <- c("#F6E8C3", "#DFC27D", "#BF812D", "#8C510A", "#FDDBC7","#F4A582", "#D6604D", "#B2182B", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC")[c(1:4, 4, 5:8, 8, 9:12, 12
)]

xend <- length(RDR_var_pool) + 0.5
gg_I <- ggplot(df_I, aes(x = RDR_var, y = accuracy, fill = Method, color = Method))  +
  geom_boxplot(alpha = 0.8,outlier.size = 0.1)  + 
  scale_fill_manual(values = method_colors) +
  scale_color_manual(values = c(rep("#5B1A18", 2*length(K_range)), rep("#084594", length(K_range)))) +
  geom_segment(aes(x = 0.5, y = 90, xend = xend, yend = 90), color = rgb(160, 52, 48, maxColorValue = 255), linetype="dotted", linewidth=0.4) + 
  geom_segment(aes(x = 0.5, y = 80, xend = xend, yend = 80), color = rgb(160, 52, 48, maxColorValue = 255), linetype="dotted", linewidth=0.4) + 
  geom_segment(aes(x = 0.5, y = 70, xend = xend, yend = 70), color = rgb(160, 52, 48, maxColorValue = 255), linetype="dotted", linewidth=0.4) + 
  geom_segment(aes(x = 0.5, y = 60, xend = xend, yend = 60), color = rgb(160, 52, 48, maxColorValue = 255), linetype="dotted", linewidth=0.4) + 
  geom_segment(aes(x = 0.5, y = 50, xend = xend, yend = 50), color = rgb(160, 52, 48, maxColorValue = 255), linetype="dotted", linewidth=0.4) + 
  geom_segment(aes(x = 0.5, y = 40, xend = xend, yend = 40), color = rgb(160, 52, 48, maxColorValue = 255), linetype="dotted", linewidth=0.4) + 
  scale_y_continuous(breaks = seq(40, 100, 10), labels = paste0(seq(40, 100, 10), "%"), position = "right") +
  scale_x_discrete(breaks = 1:length(RDR_var_pool), labels = RDR_var_pool) +
  theme(legend.position = "none", 
        axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_text(size = 34, color = "#5B1A18", margin = margin(r = 0, unit = "cm")),
        axis.text.x = element_text(size = 34, color = "#5B1A18", margin = margin(r = 0, unit = "cm")),
        axis.title.y = element_blank(), 
        panel.background = element_blank())  +
  guides(fill = guide_legend(byrow = TRUE))

ggsave(gg_I, width=24, height=12,  file = paste0(home_dir, "/result/Clusters.png"))

##### Plot states
df <- df_H
method_order <- paste0(rep(c("RDR", "both"), each = length(K_range)), rep(K_range, 2))
method_label <- paste0(rep(c("RDR (K = ", "RDR + BAF (K = "), each = length(K_range)), rep(K_range, 2), ")")
df$RDR_var <- factor(df$RDR_var, levels = RDR_var_pool, labels = 1:length(RDR_var_pool))
df$Method <- factor(df$Method, levels = method_order, labels = method_label)
xend <- length(RDR_var_pool + 2)
method_colors <- c("#F6E8C3", "#DFC27D", "#BF812D", "#8C510A", "#FDDBC7","#F4A582", "#D6604D", "#B2182B", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC")[c(5:8, 8, 9:12, 12)]


gg_H <- ggplot(df, aes(x = RDR_var, y = accuracy, fill = Method))  +
  geom_boxplot(color = "#5B1A18", alpha = 0.8,outlier.size = 0.1)  +
  scale_fill_manual(values = method_colors) +
  geom_segment(aes(x = 0.5, y = 90, xend = xend, yend = 90), color = rgb(160, 52, 48, maxColorValue = 255), linetype="dotted", size=0.4) +
  geom_segment(aes(x = 0.5, y = 80, xend = xend, yend = 80), color = rgb(160, 52, 48, maxColorValue = 255), linetype="dotted", size=0.4) +
  geom_segment(aes(x = 0.5, y = 70, xend = xend, yend = 70), color = rgb(160, 52, 48, maxColorValue = 255), linetype="dotted", size=0.4) +
  geom_segment(aes(x = 0.5, y = 60, xend = xend, yend = 60), color = rgb(160, 52, 48, maxColorValue = 255), linetype="dotted", size=0.4) +
  geom_segment(aes(x = 0.5, y = 50, xend = xend, yend = 50), color = rgb(160, 52, 48, maxColorValue = 255), linetype="dotted", size=0.4) +
  geom_segment(aes(x = 0.5, y = 40, xend = xend, yend = 40), color = rgb(160, 52, 48, maxColorValue = 255), linetype="dotted", size=0.4) +
  scale_y_continuous(breaks = seq(40, 100, 10), labels = paste0(seq(40, 100, 10), "%"), position = "right") +
  scale_x_discrete(breaks = c(1:length(RDR_var_pool)), labels = RDR_var_pool) +
  theme(legend.position = "none",
        axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 34, color = "#5B1A18", margin = margin(r = 0, unit = "cm")),
        axis.text.x = element_text(size = 34, color = "#5B1A18", margin = margin(t = -1, unit = "cm")),
        axis.title.y = element_blank(),
        panel.background = element_blank()) +
  guides(fill = guide_legend(byrow = TRUE))
ggsave(gg_H, width=24, height=12,  file = paste0(home_dir, "/result/States.png"))
