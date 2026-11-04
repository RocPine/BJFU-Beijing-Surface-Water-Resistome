# SEM: Structural equation model linking human activities, water quality,
#      microbial communities, and ARGs
# 数据: 07b_edfig4_sem_input.csv (由 07b_edfig4_sem_build_input.py 从 SI Tables 构建)
# 输出: Extended Fig. 4a/4b 的标准化路径系数表
# 运行: Rscript 07_ed_figures/07b_edfig4_sem.R
library(piecewiseSEM)
library(dplyr)

rm(list = ls())

if (!dir.exists("outputs")) dir.create("outputs")

df <- read.csv("07_ed_figures/07b_edfig4_sem_input.csv", check.names = FALSE)

Rt <- df %>%
  group_by(group) %>%
  summarise(
    across(where(is.numeric), \(x) sum(x, na.rm = TRUE) / n()),
    .groups = "drop"
  ) %>%
  as.data.frame()

wv <- c("W1","W2","W3","W4","W5","W6","W7","W9","W11","W12","W13","W14","W15")
hv <- c("H8","H9","H10","H11")

build_lv <- function(resp, vars) {
  m <- lm(as.formula(paste(resp, "~", paste(vars, collapse = "+"))), data = Rt)
  cf <- coef(m)
  as.numeric(as.matrix(Rt[, vars, drop = FALSE]) %*% cf[vars])
}

Rt$WQM <- build_lv("Microbe", wv)
Rt$HAM <- build_lv("Microbe", hv)
Rt$WQA <- build_lv("ARG", wv)
Rt$HAA <- build_lv("ARG", hv)
Rt$WQm <- build_lv("Microbe_m", wv)
Rt$HAm <- build_lv("Microbe_m", hv)
Rt$WQp <- build_lv("Microbe_p", wv)
Rt$HAp <- build_lv("Microbe_p", hv)

model <- piecewiseSEM::psem(
  lm(ARG ~ WQA + HAA + Microbe + Microbe_m + Microbe_p, data = Rt),
  lm(Microbe ~ WQM + HAM, data = Rt),
  lm(Microbe_m ~ WQm + HAm, data = Rt),
  lm(Microbe_p ~ WQp + HAp, data = Rt),
  lm(WQM ~ HAM, data = Rt)
)

s <- summary(model, .progressBar = FALSE)
print(s$coefficients)

coef_tab <- as.data.frame(s$coefficients)
if ("Response" %in% names(coef_tab) && "Predictor" %in% names(coef_tab)) {
  coef_tab <- coef_tab[, c("Response", "Predictor", "Std.Estimate", "P.Value")]
  names(coef_tab) <- c("Response", "Predictor", "Std_Estimate", "P_value")
} else {
  coef_tab$Predictor <- rownames(coef_tab)
  coef_tab <- coef_tab[, c("Predictor", "Std.Estimate", "P.Value")]
  names(coef_tab) <- c("Predictor", "Std_Estimate", "P_value")
}
write.csv(coef_tab, file.path("outputs", "07b_edfig4a_sem_coefficients.csv"),
          row.names = FALSE)
cat("\n已写入 outputs/07b_edfig4a_sem_coefficients.csv\n")

# ---- Extended Fig. 4b: 个体环境变量对各响应组的标准化系数 ----
# W 变量 (水质) 与 H 变量 (城市化) 分别回归 (与原 SEM 脚本 model1-8 口径一致)
wv_all <- c("W1","W2","W3","W4","W5","W6","W7","W9","W11","W12","W13","W14","W15")
hv_all <- c("H8","H9","H10","H11")
W_LABEL <- c(W1="pH", W2="DO", W3="Temperature", W4="TDS", W5="EC",
             W6="Turbidity", W7="COD", W9="TOC", W11="TN", W12="NO3-N",
             W13="NO2-N", W14="NH3-N", W15="TP")
H_LABEL <- c(H8="NTL", H9="PISA", H10="NDVI", H11="HSCI")
RESP_LABEL <- c(Microbe="M", Microbe_m="m", Microbe_p="p", ARG="A")

fig2d <- list()
for (resp in c("Microbe", "Microbe_m", "Microbe_p", "ARG")) {
  m_w <- lm(as.formula(paste(resp, "~", paste(wv_all, collapse = "+"))), data = Rt)
  m_h <- lm(as.formula(paste(resp, "~", paste(hv_all, collapse = "+"))), data = Rt)
  for (m in list(m_w, m_h)) {
    cf <- as.data.frame(piecewiseSEM::coefs(m, standardize = "scale"))
    for (i in seq_len(nrow(cf))) {
      v <- cf$Predictor[i]
      label <- if (v %in% names(W_LABEL)) W_LABEL[v] else if (v %in% names(H_LABEL)) H_LABEL[v] else NULL
      if (!is.null(label)) {
        fig2d[[length(fig2d) + 1]] <- data.frame(
          Variable = label, Response = RESP_LABEL[[resp]],
          Std_Estimate = cf$Std.Estimate[i], P_value = cf$P.Value[i]
        )
      }
    }
  }
}
fig2d_tab <- do.call(rbind, fig2d)
write.csv(fig2d_tab, file.path("outputs", "07b_edfig4b_sem_coefficients.csv"),
          row.names = FALSE)
print(fig2d_tab)
cat("\n已写入 outputs/07b_edfig4b_sem_coefficients.csv\n")
