rm(list = ls())

# install.packages("pacman")
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survival, survminer, skimr, lattice, reshape2, grDevices, ggrepel)

##Ready list_of_patients_WARD
list_of_patients_WARD <- read.csv2("Original_Data/EB1.0 43-plex cleaned.csv")
list_of_patients_WARD <- list_of_patients_WARD %>%
  filter(!grepl("_D28", ID)) %>%
  distinct()
list_of_patients_WARD <- list_of_patients_WARD[list_of_patients_WARD$ICU != "NoICU_HV", ]

##delete unnecessary columns
list_of_patients_WARD$Plaat <- NULL
list_of_patients_WARD$ICU <- NULL
list_of_patients_WARD$FerGr<- NULL
list_of_patients_WARD$Sepsis<- NULL
list_of_patients_WARD$Study<- NULL
list_of_patients_WARD$Pentraxin_3 <- NULL

##rename the colnames
colnames(list_of_patients_WARD)[colnames(list_of_patients_WARD) == "Lipocalin_2_NGAL"] <- "NGAL"
colnames(list_of_patients_WARD)[colnames(list_of_patients_WARD) == "Factor_XIV_protein_C"] <- "protein_C"

##CRP* 1000000
list_of_patients_WARD$CRP <- as.numeric(list_of_patients_WARD$CRP) * 1000000

# import studydata and cts information
studydata_cts <- read.csv("Original_Data/studydata_with_cts.csv")

list_of_patients_WARD <- list_of_patients_WARD[list_of_patients_WARD$ID %in% studydata_cts$EB_id, ]
list_of_patients_WARD <- merge(list_of_patients_WARD, 
                               studydata_cts[, c("EB_id", "CTS")],
                               by.x = "ID", by.y = "EB_id",
                               all.x = TRUE)

## ======= Load packages =======
pacman::p_load(tidyverse, rstatix, forcats, patchwork, ggplot2)

## ======= Load data =======
# list_of_patients_WARD should already be in env and contain CTS
# list_of_patients_WARD <- read.csv("Original_data/list_of_patients_WARD.csv", check.names = FALSE)

## ======= Define domain marker lists =======
domain_markers <- list(
  "Coagulation & Endothelium" = c(
    "D_dimer","Thrombomodulin","TFPI","Tie_2","Angiopoietin_1","Angiopoietin_2",
    "Endocan","Fractalkine","E_selectin","VCAM_1","vWF","ADAMTS13","Gas6"
  ),
  "Cytokines" = c(
    "IL_1RA","IL_6","IL_8","IL_10","IL_17","IL_23","IL_27","April","MFG_E8","PCSK9"
  ),
  "Neutrophil degranulation" = c(
    "NGAL","MPO","Proteinase_3","CD163","Resistin","TREM_1"
  ),
  "Systemic inflammation & organ damage" = c(
    "Tenascin_C","TFF3","Ferritin","NTproBNP","Procalcitonin",
    "CRP","Cystatin_C","Cardiac_Myoglobin","Syndecan"
  )
)

## ======= Data cleaning =======
df <- list_of_patients_WARD %>%
  mutate(across(-c(ID, CTS), as.numeric)) %>%
  mutate(CTS = as.factor(CTS)) %>%
  pivot_longer(-c(ID, CTS), names_to = "marker", values_to = "value") %>%
  mutate(value_num = log10(pmax(value, 1e-9)))

## ======= Hedges' g function (keep original calculation contrasts) =======
get_effsize <- function(data) {
  # DO NOT change calculation directions:
  # keep (1 vs 2), (2 vs 3), (1 vs 3)
  pairs <- tibble::tribble(
    ~first, ~second,
    "1", "2",
    "2", "3",
    "1", "3"
  )
  
  purrr::pmap_dfr(pairs, function(first, second) {
    dat <- data %>%
      filter(CTS %in% c(first, second)) %>%
      mutate(CTS = fct_relevel(CTS, second, first))
    
    n_groups <- length(unique(dat$CTS[!is.na(dat$value_num)]))
    if (n_groups < 2) {
      message("Skip ", first, " vs ", second, ": one group missing data")
      return(NULL)
    }
    
    # Compute effect size (auto fallback if hedges correction not supported)
    es <- tryCatch({
      dat %>%
        group_by(marker) %>%
        rstatix::cohens_d(value_num ~ CTS,
                          var.equal = FALSE,
                          hedges.correction = TRUE,
                          ci = TRUE) %>%
        ungroup()
    }, error = function(e) {
      dat %>%
        group_by(marker) %>%
        rstatix::cohens_d(value_num ~ CTS,
                          var.equal = FALSE,
                          hedges.correction = FALSE,
                          ci = TRUE) %>%
        ungroup()
    })
    
    # Auto-detect effect size column and rename to g
    colnames_es <- names(es)
    effect_col <- intersect(colnames_es, c("effsize","d","cohens_d","estimate","effect_size"))
    if (length(effect_col) == 0) stop("Effect size column not recognized; got: ", paste(colnames_es, collapse=", "))
    colnames(es)[match(effect_col[1], colnames_es)] <- "g"
    
    # t-test + BH correction
    tt <- dat %>%
      group_by(marker) %>%
      rstatix::t_test(value_num ~ CTS, var.equal = FALSE) %>%
      ungroup() %>%
      mutate(p_BH = p.adjust(p, "BH")) %>%
      select(marker, statistic, p_BH)
    
    # Combine
    es %>%
      left_join(tt, by = "marker") %>%
      mutate(
        contrast_lab = sprintf("CTS%s vs CTS%s", first, second),
        p_BH.signif = case_when(
          p_BH < 1e-4 ~ "****",
          p_BH < 1e-3 ~ "***",
          p_BH < 1e-2 ~ "**",
          p_BH < 0.05 ~ "*",
          TRUE ~ "ns"
        ),
        dir = ifelse(statistic > 0, "UP", "DOWN"),
        signif_dir = ifelse(p_BH.signif == "ns", "ns", paste(p_BH.signif, dir))
      )
  })
}

## ======= Color scheme for Hedges' g =======
color_mapping <- c(
  "ns"="#888888",
  "* UP"="#ffd4c4","** UP"="#FF9F94","*** UP"="#ff5c33","**** UP"="#cc2900",
  "* DOWN"="#c2d1ff","** DOWN"="#99b3ff","*** DOWN"="#668cff","**** DOWN"="blue"
)

## ======= Desired display order and label mapping (ONLY titles change) =======
# Keep calculation labels as: "CTS1 vs CTS2", "CTS2 vs CTS3", "CTS1 vs CTS3"
# But DISPLAY as (and order in facets): "CTS2 vs CTS1", "CTS3 vs CTS1", "CTS3 vs CTS2"
contrast_order_calc  <- c("CTS1 vs CTS2", "CTS1 vs CTS3", "CTS2 vs CTS3")
contrast_display_map <- c(
  "CTS1 vs CTS2" = "CTS2 vs CTS1",
  "CTS1 vs CTS3" = "CTS3 vs CTS1",
  "CTS2 vs CTS3" = "CTS3 vs CTS2"
)

## ======= Plot Hedges' g for each domain (titles/order changed only) =======
for (domain in names(domain_markers)) {
  markers_now <- domain_markers[[domain]]
  df_sub <- df %>% filter(marker %in% markers_now)
  
  eff_df <- get_effsize(df_sub) %>%
    mutate(
      # Order facets by desired display sequence while keeping underlying calc labels
      contrast_lab = factor(contrast_lab, levels = contrast_order_calc)
    )
  
  p <- ggplot(eff_df, aes(x = g, y = fct_reorder(marker, g))) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.3) +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.3, size = 0.4) +
    geom_point(aes(fill = signif_dir), shape = 22, size = 3) +
    scale_fill_manual(values = color_mapping, na.value = "#888888") +
    facet_wrap(
      ~contrast_lab, nrow = 1,
      labeller = labeller(contrast_lab = as_labeller(contrast_display_map))
    ) +
    labs(
      title = paste0(domain, " — Hedges’ g across CTS groups"),
      x = "Hedges’ g (positive = higher in first group)",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    )
  
  print(p)
  ggsave(paste0("Figure/Hedges_", gsub("[^A-Za-z]", "_", domain), ".svg"),
         p, width = 10, height = 6, dpi = 300)
}

## ======= Plot boxplots for each domain (unchanged) =======
colors_cts <- c("1"="#f28e8e","2"="#5bc0de","3"="#a1d99b")

for (domain in names(domain_markers)) {
  markers_now <- domain_markers[[domain]]
  df_sub <- df %>% filter(marker %in% markers_now)
  
  p <- ggplot(df_sub, aes(x = CTS, y = value_num, fill = CTS)) +
    geom_boxplot(width = 0.6, outlier.size = 0.8) +
    geom_jitter(width = 0.1, alpha = 0.4, size = 0.8) +
    facet_wrap(~marker, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = colors_cts) +
    labs(
      title = paste0(domain, " — Marker distributions across CTS"),
      x = NULL, y = "log10(value)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )
  
  print(p)
  ggsave(paste0("Figure/Boxplot_", gsub("[^A-Za-z]", "_", domain), ".svg"),
         p, width = 12, height = 6, dpi = 300)
}

