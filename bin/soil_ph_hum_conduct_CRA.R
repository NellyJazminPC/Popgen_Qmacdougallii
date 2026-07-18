# -------------------- Paquetes --------------------
library(tidyverse)
library(readxl)
library(janitor)
library(stringr)
library(forcats)
library(patchwork)

# -------------------- Config común --------------------
pal_season   <- c("SECAS" = "#F0E442", "LLUVIAS" = "#0072B2") # color-blind
pal_texture  <- c("Limo" = "#999999", "Arena" = "#E69F00", "Arcilla" = "#D55E00")

# Función helper para limpiar/ordenar etiquetas de sitio "1_CZ" -> "CZ" y respetar orden
clean_site <- function(x) {
  num  <- as.numeric(str_extract(as.character(x), "^[0-9]+"))
  lbl  <- str_replace(as.character(x), "^[0-9]+_", "")
  fct_reorder(lbl, num, .fun = min, .na_rm = TRUE)
}

# -------------------- Carga y preparación --------------------
df0 <- read_excel("../data/soil/soil.xlsx", sheet = 1) %>%
  clean_names() %>%
  # renombres a llaves usadas abajo
  rename(
    ph            = p_h,
    humedad_pct   = percent_humedad
  ) %>%
  mutate(
    temporada = factor(str_to_upper(temporada), levels = c("SECAS","LLUVIAS")),
    sitio_lbl = clean_site(sitio),
    zona      = factor(str_to_title(zona), levels = c("Zona Norte","Zona Sur"))
  )

# -------------------- (1) BOXPLOTS por sitio/temporada --------------------
vars_plot <- c("ph","conductividad","humedad_pct","cra")

df_long <- df0 %>%
  pivot_longer(all_of(vars_plot), names_to = "variable", values_to = "valor") %>%
  mutate(
    variable = factor(variable,
                      levels = c("ph","conductividad","humedad_pct","cra"),
                      labels = c("pH","Conductivity (µS/cm)","Moisture (%)","WHC")
    )
  )

plot_var <- function(dat, y_limits = NULL, add_ph_ref = FALSE, title_lab = NULL) {
  ggplot(dat, aes(x = sitio_lbl, y = valor, fill = temporada)) +
    geom_boxplot(outlier.alpha = 0.35,
                 position = position_dodge(width = 0.7), width = 0.6) +
    stat_summary(fun = median, geom = "point", shape = 21, size = 2,
                 color = "black", fill = "white",
                 position = position_dodge(width = 0.7)) +
    scale_fill_manual(values = pal_season, labels = c("Dry","Rainy")) +
    labs(x = "Site", y = NULL, title = title_lab, fill = "Season") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top") +
    {if (!is.null(y_limits)) scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0,0.05))) else NULL} +
    {if (add_ph_ref) geom_hline(yintercept = 7, linetype = 2, color = "grey40") else NULL}
}

p_ph   <- df_long %>% filter(variable == "pH")                     %>% plot_var(c(3,5.5), TRUE,  "pH")
p_cond <- df_long %>% filter(variable == "Conductivity (µS/cm)")   %>% plot_var(NULL,  FALSE, "Conductivity (µS/cm)")
p_hum  <- df_long %>% filter(variable == "Moisture (%)")           %>% plot_var(NULL,  FALSE, "Moisture (%)")
p_cra  <- df_long %>% filter(variable == "WHC")                    %>% plot_var(NULL,  FALSE, "WHC")

p_final <- (p_ph | p_cond) / (p_hum | p_cra) +
  plot_annotation(
    title = "Soil properties by site and season",
    subtitle = "Boxplots (Dry vs. Rainy): pH (cap at 5.5), Conductivity, Moisture, and WHC",
    theme = theme(plot.title = element_text(face = "bold", size = 14),
                  plot.subtitle = element_text(size = 12))
  )

print(p_final)

ggsave("../data/soil/boxplots_suelo_por_sitio_colorblind_opt1.png",
       p_final, width = 12, height = 9, dpi = 300)

# -------------------- (2) TEXTURA por SITIO (promedio) --------------------
g_tex_site <- df0 %>%
  group_by(zona, sitio_lbl) %>%
  summarise(
    Limo    = mean(percent_limo,   na.rm = TRUE),
    Arena   = mean(percent_arena,  na.rm = TRUE),
    Arcilla = mean(percent_arcilla,na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(Limo, Arena, Arcilla), names_to = "fraccion", values_to = "porcentaje") %>%
  mutate(fraccion = factor(fraccion, levels = c("Limo","Arena","Arcilla"))) %>%
  ggplot(aes(x = sitio_lbl, y = porcentaje, fill = fraccion)) +
  geom_col(width = 0.8, color = "white", size = 0.2) +
  coord_flip() +
  facet_grid(zona ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = pal_texture, name = NULL) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20), expand = c(0,0)) +
  labs(x = NULL, y = "Porcentaje") +
  theme_minimal(base_size = 12) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    panel.spacing.y = unit(8,"pt"),
    legend.position = "top"
  )
print(g_tex_site)
ggsave("../data/soil/textura_por_sitio.png", g_tex_site, width = 5.5, height = 6.5, dpi = 300)

# -------------------- (3) TEXTURA por INDIVIDUO (lista continua) --------------------
g_tex_ind <- df0 %>%
  select(zona, nom_muestra, percent_limo, percent_arena, percent_arcilla) %>%
  pivot_longer(
    c(percent_limo, percent_arena, percent_arcilla),
    names_to = "fraccion", values_to = "porcentaje"
  ) %>%
  mutate(
    fraccion = recode(fraccion,
                      percent_limo = "Limo",
                      percent_arena = "Arena",
                      percent_arcilla = "Arcilla"),
    # ordenar individuos 1 -> 9 de arriba a abajo
    nom_muestra = fct_rev(factor(nom_muestra, levels = unique(nom_muestra)))
  ) %>%
  ggplot(aes(x = nom_muestra, y = porcentaje, fill = fraccion)) +
  geom_col(width = 0.8, color = "white", size = 0.2) +
  coord_flip() +
  facet_grid(zona ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = pal_texture, name = NULL) +
  scale_y_continuous(limits = c(0,100), breaks = seq(0,100,20), expand = c(0,0)) +
  labs(x = NULL, y = "Percentage") +
  theme_minimal(base_size = 11) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    panel.spacing.y = unit(6,"pt"),
    legend.position = "top",
    axis.text.y = element_text(size = 7)
  )
print(g_tex_ind)
ggsave("../data/soil/textura_por_individuo.png", g_tex_ind, width = 8, height = 10, dpi = 300)
