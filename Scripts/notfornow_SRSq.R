rm(list = ls())

## =========================
## 0) Load dependencies
## =========================
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(edgeR)
  library(SepstratifieR)
})

## =========================
## 1) File paths
## =========================
file_wbs  <- "Original_data/whole_blood_stimuli_long.csv"   # ID, group
file_key  <- "Original_data/SK_Amsterdam.csv"               # Sample.ID, RNAseq, Sample
file_cnts <- "Original_data/Combined_counts_Amsterdam.csv"  # raw counts (columns = samples)
stopifnot(file.exists(file_wbs), file.exists(file_key), file.exists(file_cnts))

## =========================
## 2) Read metadata & keep CAP (adjust if needed)
## =========================
wbs   <- read.csv(file_wbs)
dfwbs <- wbs %>% select(ID, group) %>% distinct()
target_ids <- dfwbs %>% filter(group %in% c("CAP")) %>% pull(ID)

key <- read.csv(file_key)
Elder <- key %>%
  filter(grepl("EB|ELDER", Sample.ID), RNAseq == 1)

Elder$ID <- stringr::str_match(Elder$Sample.ID, pattern = "\\d\\d\\d\\d_")[,1]
Elder$ID <- stringr::str_remove(Elder$ID, "_")

dfwb_sel <- Elder %>% filter(ID %in% target_ids)

## =========================
## 3) Read counts and match samples
## =========================
## 3) Read counts and set Ensembl IDs as rownames (strip version)
cnt <- read.csv(file_cnts, row.names = 1, check.names = FALSE)  # <- row 1 are IDs
cnt <- as.data.frame(cnt, check.names = FALSE)

ids <- sub("\\..*$", "", rownames(cnt))   # remove .16/.7 etc.
stopifnot(length(ids) == nrow(cnt))
rownames(cnt) <- make.unique(ids)         # ensure uniqueness (just in case)

cnt[[ens_col]] <- NULL

stopifnot("Sample" %in% colnames(dfwb_sel))
matched_samples <- intersect(dfwb_sel$Sample, colnames(cnt))
if (length(matched_samples) < 25)
  stop("Matched samples < 25. Check mapping between key$Sample and counts colnames.")

expr_counts <- cnt[, matched_samples, drop = FALSE]

## —— merge back group info for matched samples —— ##
meta <- dfwb_sel %>%
  filter(Sample %in% matched_samples) %>%
  mutate(ID = as.character(ID)) %>%
  left_join(dfwbs %>% mutate(ID = as.character(ID)), by = "ID") %>%
  distinct(Sample, .keep_all = TRUE)

## =========================
## 4) counts -> logCPM (genes x samples)
## =========================
dge <- edgeR::DGEList(counts = expr_counts)
expr_logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

## =========================
## 5) Run SRS classifiers (no manual gene list needed)
##    stratifyPatients expects samples x genes, so transpose.
## =========================
X <- t(expr_logcpm)  # samples x genes

pred7  <- SepstratifieR::stratifyPatients(X, gene_set = "davenport") # 7-gene
pred19 <- SepstratifieR::stratifyPatients(X, gene_set = "extended")  # 19-gene

## =========================
## 6) Extract SRSq + collapse SRS3 -> nearest SRS1/2
## =========================
## X is what you passed into stratifyPatients (samples x genes)
samp_names <- rownames(X)

grab <- function(pred, samp, suffix){
  slots <- slotNames(pred)
  
  # be robust across versions
  srs  <- if ("SRS"       %in% slots) slot(pred, "SRS")
  else if ("SRS_group" %in% slots) slot(pred, "SRS_group")
  else rep(NA_character_, length(samp))
  
  srsq <- if ("SRSq"      %in% slots) slot(pred, "SRSq")
  else if ("score"     %in% slots) slot(pred, "score")
  else rep(NA_real_, length(samp))
  
  out <- data.frame(
    Sample = samp,
    SRS_raw = as.character(srs),
    SRSq    = as.numeric(srsq),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  
  if ("SRS_probs" %in% slots && !is.null(slot(pred, "SRS_probs"))) {
    probs <- as.data.frame(slot(pred, "SRS_probs"), check.names = FALSE)
    names(probs) <- paste0(names(probs), "_", suffix)
    out <- cbind(out, probs)
  }
  
  names(out)[names(out) == "SRS_raw"] <- paste0("SRS_", suffix, "_raw")
  names(out)[names(out) == "SRSq"]    <- paste0("SRSq_", suffix)
  out
}

## Ensure X exists and has sample names
X <- t(expr_logcpm)                         # samples x genes
if (is.null(rownames(X)) || any(rownames(X) == "")) {
  rownames(X) <- colnames(expr_logcpm)      # fallback
}
samp_names <- rownames(X)

## run the models
pred7  <- SepstratifieR::stratifyPatients(X, gene_set = "davenport")
pred19 <- SepstratifieR::stratifyPatients(X, gene_set = "extended")

## now this works
res7  <- grab(pred7,  samp_names, "7")
res19 <- grab(pred19, samp_names, "19")


# collapse function: if SRS3, choose higher of P(SRS1) vs P(SRS2)
collapse_2class <- function(df, suffix){
  raw_col <- paste0("SRS_", suffix, "_raw")
  p1 <- paste0("SRS1_", suffix)
  p2 <- paste0("SRS2_", suffix)
  call2 <- ifelse(df[[raw_col]] == "SRS3",
                  ifelse(df[[p1]] >= df[[p2]], "SRS1", "SRS2"),
                  df[[raw_col]])
  factor(call2, levels = c("SRS1","SRS2"))
}

res7$SRS_7   <- collapse_2class(res7, "7")
res19$SRS_19 <- collapse_2class(res19, "19")

## =========================
## 7) Assemble final table (SRSq + SRS1/2 only)
## =========================
out <- meta %>%
  select(Sample, ID, group, everything()) %>%
  left_join(res7  %>% select(Sample, SRSq_7,  SRS_7),  by = "Sample") %>%
  left_join(res19 %>% select(Sample, SRSq_19, SRS_19), by = "Sample") %>%
  distinct(Sample, .keep_all = TRUE)

## =========================
## 8) Save
## =========================
dir.create("Results", showWarnings = FALSE)
write.csv(out, "Results/SRSq_and_SRS12_only.csv", row.names = FALSE)
cat("\n✅ Wrote Results/SRSq_and_SRS12_only.csv with columns: Sample, ID, group, SRSq_7, SRS_7, SRSq_19, SRS_19\n")

## (Optional) quick sanity checks
# table(out$SRS_7); table(out$SRS_19)
# cor(out$SRSq_7, out$SRSq_19, use="complete.obs")

dir.create("Results", showWarnings = FALSE)
write.csv(out, "Original_data/SRSq_and_SRS12_only.csv", row.names = FALSE)

## =========================####
## 9) Quick QC (optional)
## =========================
if (all(c("SRSq_7","SRSq_19") %in% names(out))) {
  plot(out$SRSq_7, out$SRSq_19,
       xlab="SRSq (7 genes)", ylab="SRSq (19 genes)",
       main="Concordance of SRSq: 7 vs 19")
  abline(0,1,lty=2)
  mtext(sprintf("Pearson r = %.3f", cor(out$SRSq_7, out$SRSq_19, use="complete.obs")),
        side=3, line=0.25)
}
if (all(c("group","SRSq_7","SRSq_19") %in% names(out))) {
  par(mfrow=c(1,2))
  boxplot(SRSq_7 ~ group,  data=out, main="SRSq (7) by group")
  boxplot(SRSq_19 ~ group, data=out, main="SRSq (19) by group")
  par(mfrow=c(1,1))
  cat("\nWilcoxon (7):\n");  print(wilcox.test(SRSq_7  ~ group, data=out))
  cat("\nWilcoxon (19):\n"); print(wilcox.test(SRSq_19 ~ group, data=out))
}
if (all(c("SRS_7","SRS_19") %in% names(out))) {
  cat("\nCross-tab SRS (7 vs 19):\n")
  print(table(out$SRS_7, out$SRS_19, useNA="ifany"))
}

cat("\n✅ Done: Results/SRS_SRSq_7vs19_merged.csv has been generated.\n")

library(ggplot2)

out$agree <- ifelse(out$SRS_7 == out$SRS_19, "Agree", "Disagree")

ggplot(out, aes(x = SRSq_7, y = SRSq_19, color = agree)) +
  geom_point(size=2) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(x="SRSq (7 genes)", y="SRSq (19 genes)", color="7 vs 19") +
  theme_classic()

library(ggplot2)
library(tidyr)

plot_df <- out %>%
  dplyr::select(Sample, SRSq_7, SRSq_19) %>%
  pivot_longer(cols = c(SRSq_7, SRSq_19), names_to = "Signature", values_to = "SRSq")

ggplot(plot_df, aes(x = Signature, y = SRSq, group = Sample)) +
  geom_line(alpha=0.3) +
  geom_point(aes(color = Signature)) +
  theme_classic() +
  labs(y="SRSq", x="Signature (7 vs 19)")

tab <- table(out$SRS_7, out$SRS_19)
tab

library(ggplot2)
tab_df <- as.data.frame(tab)
ggplot(tab_df, aes(x=Var2, y=Var1, fill=Freq)) +
  geom_tile(color="white") +
  geom_text(aes(label=Freq)) +
  scale_fill_gradient(low="white", high="steelblue") +
  labs(x="SRS (19 genes)", y="SRS (7 genes)", fill="Count") +
  theme_classic()
#####
