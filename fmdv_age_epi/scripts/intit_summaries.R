
rm(list=ls())
pacman::p_load(tidyverse,
               here,
               ggplot2,
               ComplexHeatmap,
               patchwork,
               sf,
               ggspatial,
               mgcv,
               lme4)



loc <- "C:/Users/sarac/Desktop/fmdv/buffalo_fmdv/fmdv_age_epi"

#############################
# Clean up, merge, and subset data
#############################
mass_cap <- read.csv(here(loc, "dat_from_bree/masscap_fielddata (1).csv"))
recap_1 <- read.csv(here(loc, "dat_from_bree/recap1_0923_1023_fielddata (1).csv"))
recap_2 <- read.csv(here(loc, "dat_from_bree/recap2_0424_0524_fielddata (1).csv"))
recap_3 <- read.csv(here(loc, "dat_from_bree/recap3_0924_1024_fielddata (1).csv"))
recap_4 <- read.csv(here(loc, "dat_from_bree/recap4_0425_0525_fielddata (1).csv"))



recaps <- dplyr::full_join(recap_1, recap_2) %>%
  dplyr::full_join(recap_3) 

recaps$biologger_no = as.character(recaps$biologger_no)

recaps <- recaps %>%
  dplyr::full_join(recap_4)


all_caps <- dplyr::full_join(mass_cap, recaps) %>%
  dplyr::select(animal_id, unique_id, recap_no, firstcap_age, sex, capture_date, capture_lat, capture_long)
all_caps$capture_date <- as.Date(all_caps$capture_date, format = "%m/%d/%Y")
all_caps$capture_num <- ifelse(all_caps$capture_date < '2023-08-01', 0, all_caps$recap_no)
summary(as.factor(all_caps$capture_num))
str(all_caps)
all_caps[all_caps$animal_id==132, "firstcap_age"] <- 1

swab_dat <- read.csv(here(loc, "dat_from_bree/tonsil_swab_dat.csv"))
test_dat <- read.csv(here(loc, "dat_from_bree/masscap_elisa.csv"))

all_caps$animal_id <- as.character(all_caps$animal_id)
#swab_dat$capture_num <- ifelse(is.na(swab_dat$recap_no) ==TRUE, 0, swab_dat$recap_n)
test_dat$animal_id <- as.character(test_dat$animal_id)
#swab_dat$vi[swab_dat$vi==""] <- "NO"

swab_dat <- dplyr::full_join(swab_dat, all_caps) %>%
  dplyr::full_join(test_dat) 
summary(swab_dat$capture_lat)
str(swab_dat)


swab_dat <- swab_dat %>%
  dplyr::mutate(ts_pir_vi = case_when(vi == "YES" | vi == "yes" ~ 1,
                                      vi == "" ~ 0),
                ts_pir_qpcr_interp = dplyr::case_when(rt_qpcr <= 30 ~ 'strong_pos',
                                              rt_qpcr > 30 & rt_qpcr <= 40 ~ 'weak_pos',
                                              rt_qpcr > 40 ~ 'neg',
                                              rt_qpcr == "Undetermined" ~ 'neg'),
                ts_pir_qpcr_pos = dplyr::case_when(rt_qpcr <= 40 ~ 1,
                                                   rt_qpcr > 40 ~ 0)) %>%
  dplyr::rename(ts_pir_qpcr_ct = rt_qpcr) %>%
  dplyr::select(-c(sample_type, vi)) %>%
  dplyr::filter(is.na(capture_date) == FALSE)
swab_dat[swab_dat==""] <- NA

summary(as.factor(swab_dat$ts_pir_qpcr_ct))
swab_dat$ts_pir_qpcr_ct = dplyr::replace_values(swab_dat$ts_pir_qpcr_ct, "Undetermined" ~ "45") %>%
  as.numeric()
swab_dat$ts_arc_qpcr_ct = dplyr::replace_values(swab_dat$ts_arc_qpcr_ct, 0 ~ 45) %>%
  as.numeric()
summary((swab_dat$ts_arc_qpcr_ct))
swab_dat$probang_arc_qpcr_ct = dplyr::replace_values(swab_dat$probang_arc_qpcr_ct, 0 ~ 45) %>%
  as.numeric()
swab_dat$ts_arc_qpcr_pos = dplyr::replace_values(swab_dat$ts_arc_qpcr_pos, c("pos","Pos") ~ "1", c("neg","Neg") ~ "0") %>%
  as.numeric()
swab_dat$probang_arc_qpcr_pos = dplyr::replace_values(swab_dat$probang_arc_qpcr_pos, c("pos","Pos") ~ "1", c("neg","Neg") ~ "0") %>%
  as.numeric()
swab_dat$sat1_elisa_pos = dplyr::replace_values(swab_dat$sat1_elisa_pos, "Pos" ~ "1", "Neg" ~ "0", "Incl" ~ "2") %>%
  as.numeric()
swab_dat$sat2_elisa_pos = dplyr::replace_values(swab_dat$sat2_elisa_pos, "Pos" ~ "1", "Neg" ~ "0", "Incl" ~ "2") %>%
  as.numeric()
swab_dat$sat3_elisa_pos = dplyr::replace_values(swab_dat$sat3_elisa_pos, "Pos" ~ "1", "Neg" ~ "0", "Incl" ~ "2") %>%
  as.numeric()
swab_dat$nsp_elisa_pos = dplyr::replace_values(swab_dat$nsp_elisa_pos, "Pos" ~ "1", "Neg" ~ "0") %>%
  as.numeric()
swab_dat$sex = stringr::str_to_lower(swab_dat$sex)

str(swab_dat)



series <- all_caps %>%
  dplyr::select(animal_id, capture_date, firstcap_age, capture_num) %>%
  dplyr::group_by(animal_id) %>%
  dplyr::distinct(capture_date, .keep_all = TRUE) %>%
  dplyr::mutate(capt_order = row_number()) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(id_cols = c("animal_id","firstcap_age"), names_from = "capt_order", names_prefix = "capture_date_", values_from = c("capture_date")) %>%
  dplyr::mutate(cap_age_1 = firstcap_age,
                cap_age_2 = as.numeric(difftime(capture_date_2, capture_date_1, units = "days"))/365 + firstcap_age,
                cap_age_3 = as.numeric(difftime(capture_date_3, capture_date_1, units = "days"))/365 + firstcap_age,
                cap_age_4 = as.numeric(difftime(capture_date_4, capture_date_1, units = "days"))/365  + firstcap_age,
                cap_age_5 = as.numeric(difftime(capture_date_5, capture_date_1, units = "days"))/365  + firstcap_age) %>%
  tidyr::pivot_longer(cols = matches("^(capture_date|cap_age)_\\d$"),
                      names_to = c(".value", "cap_num"),
                      names_pattern = "(.*)_(\\d+)") %>%
  dplyr::rename(capture_age = cap_age) %>%
  dplyr::select(-c(cap_num)) %>%
  na.omit()


swab_dat <- dplyr::full_join(swab_dat, series)
swab_dat <- swab_dat %>%
  dplyr::select(c(animal_id, firstcap_age, sex, capture_num, capture_age, unique_id, sample_id, capture_date, capture_lat, capture_long, 
                  ts_pir_qpcr_ct, ts_pir_qpcr_interp, ts_pir_qpcr_pos, ts_pir_vi,
                  ts_arc_qpcr_ct, ts_arc_qpcr_pos, probang_arc_qpcr_ct, probang_arc_qpcr_pos,
                  nsp_elisa_pi, nsp_elisa_pos, sat1_elisa_pi, sat1_elisa_pos, sat2_elisa_pi, sat2_elisa_pos, sat3_elisa_pi, sat3_elisa_pos))

str(swab_dat)

write.csv(swab_dat, here(loc, "clean_dat/fmdv_detection_dat.csv"), row.names = FALSE)

#############################
# read the data back in & make some wuick subsets
#############################
rm(list=ls())
loc <- "C:/Users/sarac/Desktop/fmdv/buffalo_fmdv/fmdv_age_epi"
dat <- read.csv(here(loc, "clean_dat/fmdv_detection_dat.csv"))
dat$capture_date <- as.Date(dat$capture_date)

dat_sub <- dat %>%
  dplyr::ungroup() %>%
  dplyr::filter(firstcap_age <= 3)

dat_first <- dat_sub %>%
  dplyr::arrange(capture_date) %>%
  dplyr::distinct(animal_id, .keep_all = TRUE)

#############################
# Some quick summaries for the readme--getting cutoff values
#############################
dat %>% 
  dplyr::group_by(ts_arc_qpcr_pos) %>%
  dplyr::summarise(min = min(ts_arc_qpcr_ct, na.rm = TRUE),
                   mean = mean(ts_arc_qpcr_ct, na.rm = TRUE),
                   max = max(ts_arc_qpcr_ct, na.rm = TRUE)) %>%
  flextable::flextable()

dat %>% 
  dplyr::group_by(probang_arc_qpcr_pos) %>%
  dplyr::summarise(min = min(probang_arc_qpcr_ct, na.rm = TRUE),
                   mean = mean(probang_arc_qpcr_ct, na.rm = TRUE),
                   max = max(probang_arc_qpcr_ct, na.rm = TRUE)) %>%
  flextable::flextable()

dat %>% 
  dplyr::group_by(nsp_elisa_pos) %>%
  dplyr::summarise(min = min(nsp_elisa_pi, na.rm = TRUE),
                   mean = mean(nsp_elisa_pi, na.rm = TRUE),
                   max = max(nsp_elisa_pi, na.rm = TRUE)) %>%
  flextable::flextable()

dat %>% 
  dplyr::group_by(sat1_elisa_pos) %>%
  dplyr::summarise(min = min(sat1_elisa_pi, na.rm = TRUE),
                   mean = mean(sat1_elisa_pi, na.rm = TRUE),
                   max = max(sat1_elisa_pi, na.rm = TRUE)) %>%
  flextable::flextable()

dat %>% 
  dplyr::group_by(sat2_elisa_pos) %>%
  dplyr::summarise(min = min(sat2_elisa_pi, na.rm = TRUE),
                   mean = mean(sat2_elisa_pi, na.rm = TRUE),
                   max = max(sat2_elisa_pi, na.rm = TRUE)) %>%
  flextable::flextable()

dat %>% 
  dplyr::group_by(sat3_elisa_pos) %>%
  dplyr::summarise(min = min(sat3_elisa_pi, na.rm = TRUE),
                   mean = mean(sat3_elisa_pi, na.rm = TRUE),
                   max = max(sat3_elisa_pi, na.rm = TRUE)) %>%
  flextable::flextable()


dat_first %>% 
  janitor::tabyl(ts_pir_qpcr_pos, ts_pir_vi)
dat_first %>% 
  janitor::tabyl(ts_pir_qpcr_interp, ts_pir_vi)

summary(as.factor(dat_first$sat1_elisa_pos))
summary(as.factor(dat_first$sat2_elisa_pos))
summary(as.factor(dat_first$sat3_elisa_pos))
summary(as.factor(dat_first$nsp_elisa_pos))

dat_first %>% 
  janitor::tabyl(nsp_elisa_pos, sat1_elisa_pos)
dat_first %>% 
  janitor::tabyl(nsp_elisa_pos, sat2_elisa_pos)
dat_first %>% 
  janitor::tabyl(nsp_elisa_pos, sat3_elisa_pos)

#############################
# Some quick figures just for fun
#############################

ggplot(dat_sub, aes(x = capture_date, y = ts_pir_qpcr_ct, color = ts_pir_qpcr_interp, group = animal_id)) +
  geom_line() +
  geom_point() +
  theme_bw()

ggplot(dat_sub, aes(x = capture_lat, y = ts_pir_qpcr_ct, color = ts_pir_qpcr_interp, group = animal_id)) +
  geom_point() +
  theme_bw()

ggplot(dat_first, aes(x = capture_age, y = ts_pir_qpcr_ct, color = as.factor(capture_age))) +
  geom_boxplot(width = 0.2) +
  geom_point(position = "jitter") +
  theme_bw()

ggplot(dat_first, aes(x = capture_lat, y = ts_pir_qpcr_ct, color = ts_pir_qpcr_interp, group = animal_id)) +
  geom_point() +
  theme_bw()
ggplot(dat_first, aes(x = capture_lat, y = ts_pir_qpcr_ct, color = capture_age, group = animal_id)) +
  geom_point() +
  theme_bw()

dat_first_sf <- dat_first %>%
  dplyr::filter(is.na(capture_long)==FALSE) %>%
  sf::st_as_sf(coords = c("capture_long", "capture_lat"), crs = 4326)

plot <- ggplot() +
    geom_sf(data = dat_first_sf, 
            aes(color=ts_pir_qpcr_pos), 
            size = 4, alpha = 0.8) +
    annotation_scale(location = "br", 
                     width_hint = 0.3, 
                     text_cex = 0.8,
                     pad_x = unit(0.2, "cm"), 
                     pad_y = unit(0.2, "cm")) +
    annotation_north_arrow(location = "br", 
                           which_north = "true",
                           style = north_arrow_fancy_orienteering,
                           height = unit(2.2, "cm"), 
                           width = unit(2.2, "cm"),
                           pad_x = unit(1.0, "cm"), 
                           pad_y = unit(0.5, "cm")) +
    coord_sf(crs = sf::st_crs(4326)) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      axis.title.x = element_text(size = 20, face = "bold"),
      axis.title.y = element_text(size = 20, face = "bold"),
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5)
    )
  
ggsave(file.path(here(loc, "outputs/figs/qpcr_plot.png")), plot = plot, width = 12, height = 18, units = "in", dpi = 700)


dat_series <- dat_sub %>%
  group_by(animal_id) %>%
  filter(n() > 1) %>%
  ungroup()


ggplot(dat_series, aes(x = capture_date, y = ts_pir_qpcr_ct, color = (capture_age), group = animal_id)) +
  geom_line() +
  geom_point() +
  theme_bw()

hist(dat_series$firstcap_age)
summary(as.factor(dat_series$firstcap_age))


ggplot(dat_series, aes(x = capture_age, y = ts_pir_qpcr_ct, color = capture_num, group = animal_id)) +
  geom_line() +
  geom_point() +
  theme_bw()


ggplot(dat_first, aes(x = capture_age, y = sat1_elisa_pi, color = as.factor(sat1_elisa_pos))) +
  geom_point(position = "jitter") +
  theme_bw()

ggplot(dat_first, aes(x = capture_age, y = sat2_elisa_pi, color = as.factor(sat2_elisa_pos))) +
  geom_point(position = "jitter") +
  theme_bw()

ggplot(dat_first, aes(x = capture_age, y = sat3_elisa_pi, color = as.factor(sat3_elisa_pos))) +
  geom_point(position = "jitter") +
  theme_bw()

ggplot(dat_first, aes(x = capture_age, y = nsp_elisa_pi, color = as.factor(nsp_elisa_pos))) +
  geom_point(position = "jitter") +
  theme_bw()


#############################
# Now some actual data informed models
#############################

# check qpcr values
hist(dat_first$ts_pir_qpcr_ct)
hist(log(dat_first$ts_pir_qpcr_ct))
shapiro.test(dat_first$ts_pir_qpcr_ct)
shapiro.test(log(dat_first$ts_pir_qpcr_ct))

str(dat_first)
fit1 <- glm(ts_pir_qpcr_ct ~  poly(firstcap_age, 1) + sex + capture_lat, 
            data = dat_first, family = Gamma)
summary(fit1)

fit2 <- glm(ts_pir_qpcr_ct ~  poly(firstcap_age, 2) + sex + capture_lat, 
            data = dat_first, family = Gamma)
summary(fit2)

fit3 <- glm(ts_pir_qpcr_ct ~  poly(firstcap_age, 1)*capture_lat + sex + capture_lat, 
            data = dat_first, family = Gamma)
summary(fit3)

fit4 <- glm(ts_pir_qpcr_ct ~  poly(firstcap_age, 2)*capture_lat + sex + capture_lat , 
            data = dat_first, family = Gamma)
summary(fit4)


fit5 <- glm(ts_pir_qpcr_pos ~  poly(firstcap_age, 1) + sex + capture_lat, 
            data = dat_first, family = binomial)
summary(fit5)

fit6 <- glm(ts_pir_qpcr_pos ~  poly(firstcap_age, 2) + sex + capture_lat, 
            data = dat_first, family = binomial)
summary(fit6)

fit7 <- glm(ts_pir_qpcr_pos ~  poly(firstcap_age, 1)*capture_lat + sex + capture_lat, 
            data = dat_first, family = binomial)
summary(fit7)

fit8 <- glm(ts_pir_qpcr_pos ~  poly(firstcap_age, 2)*capture_lat + sex + capture_lat, 
            data = dat_first, family = binomial)
summary(fit8)


fit9 <- mgcv::gam(ts_pir_qpcr_pos ~ poly(capture_age, 2), data = dat_first, family = binomial)
summary(fit9)
fit10 <- mgcv::gam(ts_pir_qpcr_pos ~ s(capture_age), data = dat_first, family = binomial)
summary(fit10)
fit11 <- mgcv::gam(ts_pir_qpcr_pos ~ s(capture_age) + s(capture_lat), data = dat_first, family = binomial)
summary(fit11)


par(mfrow = c(1,2))
plot(fit11, all.terms = TRUE)

AIC(fit10, fit11)

fit12 <- mgcv::gam(ts_pir_qpcr_pos ~ s(firstcap_age, capture_lat), data = dat_first, family = binomial)
summary(fit12)
par(mfrow = c(1,1))
plot(fit12, page = 1, scheme = 2)
vis.gam(fit12, view = c("firstcap_age","capture_lat"),
        theta = 50, n.grid = 50, lwd = 0.2)


fit13 <- mgcv::gam(ts_pir_qpcr_pos ~ s(capture_age, capture_lat), data = dat_sub, family = binomial)
summary(fit13)
par(mfrow = c(1,1))
plot(fit13, page = 1, scheme = 2)




fit20 <- mgcv::gam(ts_pir_vi ~ s(capture_age), data = dat_first, family = binomial)
summary(fit20)
fit21 <- mgcv::gam(ts_pir_vi ~ s(capture_age) + s(capture_lat), data = dat_first, family = binomial)
summary(fit21)


par(mfrow = c(1,2))
plot(fit21, all.terms = TRUE)

AIC(fit20, fit21)

fit22 <- mgcv::gam(ts_pir_vi ~ s(firstcap_age, capture_lat), data = dat_first, family = binomial)
summary(fit22)
AIC(fit21, fit22)
par(mfrow = c(1,1))
plot(fit22, page = 1, scheme = 2)
vis.gam(fit22, view = c("firstcap_age","capture_lat"),
        theta = 50, n.grid = 50, lwd = 0.2)



fit23 <- mgcv::gam(ts_pir_vi ~ s(capture_age, capture_lat), data = dat_sub, family = binomial)
summary(fit23)
par(mfrow = c(1,1))
plot(fit23, page = 1, scheme = 2)


dat_series <- dat_sub %>%
  group_by(animal_id) %>%
  filter(n() > 1) %>%
  ungroup()

str(dat_series)


ggplot(dat_series, aes(x = capture_age, y = ts_pir_qpcr_ct, color = as.factor(capture_num), group = animal_id)) +
  geom_line() +
  geom_point() +
  theme_bw()

ggplot(dat_series, aes(x = capture_num, y = ts_pir_qpcr_ct, color = ts_pir_qpcr_interp, group = animal_id)) +
  geom_line() +
  geom_point() +
  theme_bw()

ggplot(dat_series, aes(x = capture_age, y = ts_pir_qpcr_ct, color = ts_pir_qpcr_interp, group = animal_id)) +
  geom_line() +
  geom_point() +
  theme_bw()


pcr_mat <- dat_series %>%
  dplyr::group_by(animal_id) %>%
  dplyr::distinct(capture_num, .keep_all = TRUE) %>%
  dplyr::filter(is.na(ts_pir_qpcr_interp) == FALSE) %>%
  dplyr::filter(n() > 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(c(animal_id, ts_pir_qpcr_interp, firstcap_age, capture_num)) %>%
  dplyr::arrange(capture_num) %>%
  dplyr::arrange(firstcap_age) %>%
  tidyr::pivot_wider(id_cols = c("animal_id","firstcap_age"), names_from = "capture_num", names_prefix = "cap_", values_from = c("ts_pir_qpcr_interp"), names_sort = TRUE) %>%
  dplyr::arrange(firstcap_age) 

pcr_mat %>%
  dplyr::group_by(animal_id) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n > 1L)

#pcr_mat[pcr_mat==0] <- NA
age <- pcr_mat$firstcap_age
pcr_mat <- tibble::column_to_rownames(pcr_mat, "animal_id") %>%
  dplyr::select(-c(firstcap_age))
pcr_mat <- as.matrix(pcr_mat)

# make all the heatmaps
map_pcr <- ComplexHeatmap::Heatmap(pcr_mat, 
                                    cluster_rows = FALSE, 
                                    cluster_columns = FALSE,
                                    col = c("strong_pos" = "darkred", "weak_pos" = "red3", "neg" = "gray60"), na_col = "white",
                                    #rect_gp = gpar(col = "gray20", lwd = 1),
                                    show_heatmap_legend = FALSE,
                                    row_title_gp = gpar(fontsize = 10),
                                    row_names_side = "left", row_names_gp = gpar(fontsize = 8),
                                    show_row_names =  FALSE,
                                    left_annotation = rowAnnotation(Age = age), 
                                    #                               simple_anno_size = unit(2, "mm"), show_annotation_name = FALSE, border = TRUE),
                                    column_title = "qPCR status", column_title_side = "top",
                                    column_title_gp = gpar(fontsize = 16, hjust = 0),
                                    column_names_rot = 45,
                                    column_names_gp = gpar(fontsize = 12))

map_pcr

lgd_pos <- ComplexHeatmap::Legend(title = "Detection",
                                  title_gp = gpar(fontsize = 12),
                                  labels = c("Strong positive", "Weak positve", "No detection"),
                                  labels_gp = gpar(fontsize = 10),
                                  legend_gp = gpar(fill = c("darkred", "red3", "gray60")),
                                  grid_height = unit(0.25, "cm"),
                                  grid_width = unit(0.25, "cm"))
hm_comb <- ComplexHeatmap::draw(map_pcr,  annotation_legend_list = lgd_pos) 

tiff(here(loc,"outputs/figs/qpcr_heatmap.tiff"), res = 600, width = 14, height = 24, units = "cm" )
hm_comb
dev.off()

pcr_mat <- dat_series %>%
  dplyr::group_by(animal_id) %>%
  dplyr::distinct(capture_num, .keep_all = TRUE) %>%
  dplyr::filter(is.na(ts_pir_qpcr_interp) == FALSE) %>%
  dplyr::filter(n() > 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(c(animal_id, ts_pir_qpcr_interp, firstcap_age, capture_age)) %>%
  dplyr::arrange(capture_age) %>%
  dplyr::arrange(firstcap_age) %>%
  dplyr::mutate(capture_age_cat = dplyr::case_when(capture_age <1 ~ "< 1",
                                                   capture_age >=1 & capture_age <1.25 ~ "1-1.25",
                                                   capture_age >=1.25 & capture_age <1.5 ~ "1.25-1.5",
                                                   capture_age >=1.5 & capture_age <1.75 ~ "1.5-1.75",
                                                   capture_age >=1.75 & capture_age <2 ~ "1.75-2",
                                                   capture_age >=2 & capture_age <2.25 ~ "2-2.25",
                                                   capture_age >=2.25 & capture_age <2.5 ~ "2.25-2.5",
                                                   capture_age >=2.5 & capture_age <2.75 ~ "2.5-2.75",
                                                   capture_age >=2.75 & capture_age <3 ~ "2.75-3",
                                                   capture_age >=3 & capture_age <3.25 ~ "3-3.25",
                                                   capture_age >=3.25 & capture_age <3.5 ~ "3.25-3.5",
                                                   capture_age >=3.5 ~ "3.5+")) %>%
  tidyr::pivot_wider(id_cols = c("animal_id","firstcap_age"), names_from = "capture_age_cat", names_prefix = "cap_", values_from = c("ts_pir_qpcr_interp"), names_sort = TRUE) %>%
  dplyr::arrange(firstcap_age) 

pcr_mat %>%
  dplyr::group_by(animal_id) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n > 1L)

#pcr_mat[pcr_mat==0] <- NA
age <- pcr_mat$firstcap_age
pcr_mat <- tibble::column_to_rownames(pcr_mat, "animal_id") %>%
  dplyr::select(-c(firstcap_age))
pcr_mat <- as.matrix(pcr_mat)

# make all the heatmaps
map_pcr <- ComplexHeatmap::Heatmap(pcr_mat, 
                                   cluster_rows = FALSE, 
                                   cluster_columns = FALSE,
                                   col = c("strong_pos" = "darkred", "weak_pos" = "red3", "neg" = "gray60"), na_col = "white",
                                   #rect_gp = gpar(col = "gray20", lwd = 1),
                                   show_heatmap_legend = FALSE,
                                   row_title_gp = gpar(fontsize = 10),
                                   row_names_side = "left", row_names_gp = gpar(fontsize = 8),
                                   show_row_names =  FALSE,
                                   left_annotation = rowAnnotation(Age = age), 
                                   #                               simple_anno_size = unit(2, "mm"), show_annotation_name = FALSE, border = TRUE),
                                   column_title = "qPCR status", column_title_side = "top",
                                   column_title_gp = gpar(fontsize = 16, hjust = 0),
                                   column_names_rot = 45,
                                   column_names_gp = gpar(fontsize = 12))

map_pcr

lgd_pos <- ComplexHeatmap::Legend(title = "Detection",
                                  title_gp = gpar(fontsize = 12),
                                  labels = c("Strong positive", "Weak positve", "No detection"),
                                  labels_gp = gpar(fontsize = 10),
                                  legend_gp = gpar(fill = c("darkred", "red3", "gray60")),
                                  grid_height = unit(0.25, "cm"),
                                  grid_width = unit(0.25, "cm"))
hm_comb <- ComplexHeatmap::draw(map_pcr,  annotation_legend_list = lgd_pos) 





ggplot(dat_series, aes(x = capture_age, y = ts_pir_vi, color = as.factor(capture_num), group = animal_id)) +
  geom_line() +
  geom_point() +
  theme_bw()

ggplot(dat_series, aes(x = capture_date, y = ts_pir_vi, color = (ts_pir_vi), group = animal_id)) +
  geom_line() +
  geom_point() +
  theme_bw()


vi_mat <- dat_series %>%
  dplyr::group_by(animal_id) %>%
  dplyr::distinct(capture_num, .keep_all = TRUE) %>%
  dplyr::filter(is.na(ts_pir_qpcr_interp) == FALSE) %>%
  dplyr::filter(n() > 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(c(animal_id, ts_pir_vi, firstcap_age, capture_num)) %>%
  dplyr::arrange(capture_num) %>%
  dplyr::arrange(firstcap_age) %>%
  tidyr::pivot_wider(id_cols = c("animal_id","firstcap_age"), names_from = "capture_num", names_prefix = "cap_", values_from = c("ts_pir_vi"), names_sort = TRUE) %>%
  dplyr::arrange(firstcap_age) 

vi_sums <- vi_mat %>%
  tibble::column_to_rownames("animal_id") %>%
  dplyr::select(dplyr::starts_with("cap"))%>%
  dplyr::rowwise() %>%
  dplyr::mutate(row_sum  = sum(c_across(starts_with("cap")), na.rm = TRUE),
                row_n =  sum(c_across(starts_with("cap")) >= 0, na.rm = TRUE),
                prop_pos = row_sum / row_n)
head(vi_sums)

vi_mat %>%
  dplyr::group_by(animal_id) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n > 1L)

#vi_mat[vi_mat==0] <- NA
age <- vi_mat$firstcap_age
vi_mat <- tibble::column_to_rownames(vi_mat, "animal_id") %>%
  dplyr::select(-c(firstcap_age))
vi_mat <- as.matrix(vi_mat)

# make all the heatmaps
map_vi <- ComplexHeatmap::Heatmap(vi_mat, 
                                   cluster_rows = FALSE, 
                                   cluster_columns = FALSE,
                                   col = c("gray60", "darkred"), na_col = "white",
                                   #rect_gp = gpar(col = "gray20", lwd = 1),
                                   show_heatmap_legend = FALSE,
                                   row_title_gp = gpar(fontsize = 10),
                                   row_names_side = "left", row_names_gp = gpar(fontsize = 8),
                                   show_row_names =  FALSE,
                                   left_annotation = rowAnnotation(Age = age), 
                                   #                               simple_anno_size = unit(2, "mm"), show_annotation_name = FALSE, border = TRUE),
                                   column_title = "VI status", column_title_side = "top",
                                   column_title_gp = gpar(fontsize = 16, hjust = 0),
                                   column_names_rot = 45,
                                   column_names_gp = gpar(fontsize = 12))

map_vi

lgd_pos <- ComplexHeatmap::Legend(title = "Detection",
                                  title_gp = gpar(fontsize = 12),
                                  labels = c("Positive", "No detection"),
                                  labels_gp = gpar(fontsize = 10),
                                  legend_gp = gpar(fill = c("darkred", "gray60")),
                                  grid_height = unit(0.25, "cm"),
                                  grid_width = unit(0.25, "cm"))
hm_comb <- ComplexHeatmap::draw(map_vi,  annotation_legend_list = lgd_pos) 

tiff(here(loc,"outputs/figs/vi_heatmap.tiff"), res = 600, width = 14, height = 24, units = "cm" )
hm_comb
dev.off()



fit22 <-  mgcv::gam(pcr_pos ~ s(cap_age) + season, data = dat_first, family = binomial)
summary(fit12)
par(mfrow = c(1,1))
plot(fit12, page = 1, scheme = 2)
vis.gam(fit12, view = c("capture_lat","firstcap_age"),
        theta = 50, n.grid = 50, lwd = 0.2)






dat_first <- dat_first %>%
  dplyr::mutate(sat1_elisa_bin = dplyr::case_when(sat1_elisa_pos > 0 ~ 1,
                                                  sat1_elisa_pos == 0 ~ 0),
                sat2_elisa_bin = dplyr::case_when(sat2_elisa_pos > 0 ~ 1,
                                                  sat2_elisa_pos == 0 ~ 0),
                sat3_elisa_bin = dplyr::case_when(sat3_elisa_pos > 0 ~ 1,
                                                  sat3_elisa_pos == 0 ~ 0)) 

# SAT1 ELISA

fit30 <- mgcv::gam(sat1_elisa_bin ~ s(capture_age), data = dat_first, family = binomial)
summary(fit30)
fit31 <- mgcv::gam(sat1_elisa_bin ~ s(capture_age) + s(capture_lat), data = dat_first, family = binomial)
summary(fit31)
par(mfrow = c(1,2))
plot(fit31, all.terms = TRUE)

AIC(fit30, fit31)

fit32 <- mgcv::gam(sat1_elisa_bin ~ s(firstcap_age, capture_lat), data = dat_first, family = binomial)
summary(fit32)
AIC(fit31, fit32)
par(mfrow = c(1,1))
plot(fit32, page = 1, scheme = 2)



# SAT2 ELISA

fit30 <- mgcv::gam(sat2_elisa_bin ~ s(capture_age), data = dat_first, family = binomial)
summary(fit30)
fit31 <- mgcv::gam(sat2_elisa_bin ~ s(capture_age) + s(capture_lat), data = dat_first, family = binomial)
summary(fit31)
par(mfrow = c(1,2))
plot(fit31, all.terms = TRUE)

AIC(fit30, fit31)

fit32 <- mgcv::gam(sat2_elisa_bin ~ s(firstcap_age, capture_lat), data = dat_first, family = binomial)
summary(fit32)
AIC(fit31, fit32)
par(mfrow = c(1,1))
plot(fit32, page = 1, scheme = 2)



# SAT3 ELISA

fit30 <- mgcv::gam(sat3_elisa_bin ~ s(capture_age), data = dat_first, family = binomial)
summary(fit30)
fit31 <- mgcv::gam(sat3_elisa_bin ~ s(capture_age) + s(capture_lat), data = dat_first, family = binomial)
summary(fit31)
par(mfrow = c(1,2))
plot(fit31, all.terms = TRUE)

AIC(fit30, fit31)

fit32 <- mgcv::gam(sat3_elisa_bin ~ s(firstcap_age, capture_lat), data = dat_first, family = binomial)
summary(fit32)
AIC(fit31, fit32)
par(mfrow = c(1,1))
plot(fit32, page = 1, scheme = 2)



