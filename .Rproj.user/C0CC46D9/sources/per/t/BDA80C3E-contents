library(ggplot2)
library(dplyr)
library(scales)
library(patchwork)

df_seir <- df_all_totals |> 
  filter(Khoang != "5. Miễn dịch (V)")

df_v <- df_all_totals |> 
  filter(Khoang == "5. Miễn dịch (V)")

p_seir <- ggplot(df_seir, aes(x = time, y = So_Luong, color = Khoang)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ KichBan, ncol = 3) +
  theme_bw() +
  scale_y_continuous(labels = comma) +
  scale_color_manual(values = c("1. Cảm nhiễm (S)" = "#E41A1C",  
                                "2. Ủ bệnh (E)" = "#FF7F00",     
                                "3. Đang nhiễm (I)" = "#984EA3", 
                                "4. Khỏi bệnh (R)" = "#4DAF4A")) + 
  labs(title = "A. Diễn biến các khoang dịch tễ (S, E, I, R)",
       x = NULL, 
       y = "Số lượng người", 
       color = "Nhóm dịch tễ") +
  theme(legend.position = "bottom",
        text = element_text(size = 12),
        strip.text = element_text(face = "bold", size = 11),
        plot.title = element_text(face = "bold"))

p_v <- ggplot(df_v, aes(x = time, y = So_Luong)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ KichBan, ncol = 3) +
  theme_bw() +
  scale_y_continuous(labels = comma) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none") +
  labs(title = "B. Biến động nhóm có miễn dịch (Khoang V)",
       x = "Thời gian (Ngày)", 
       y = "Số lượng người") +
  theme(legend.position = "bottom",
        text = element_text(size = 12),
        plot.title = element_text(face = "bold"))


p <- p_seir / p_v + plot_layout(heights = c(1.5, 1))

print(p)

# ggsave("Seir_V.png", plot = p, width = 14, height = 9, dpi = 300)