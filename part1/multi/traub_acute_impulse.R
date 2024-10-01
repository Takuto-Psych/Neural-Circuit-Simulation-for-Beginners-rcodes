rm(list = ls())
# i_injの8番目の要素に値を設定
# EXE_VALUES$I_Inj[9] <- 1e-4 / COMPARTMENT_VALUES$area[9]
# 条件に基づいて注入電流を設定
# EXE_VALUES$I_Inj[9] <- if (t >= 1000 && t < 1001) { 1e-2 / COMPARTMENT_VALUES$area[9] } else { 0 }
# の処理の後者


library(tidyverse)


VALUES <- data.frame(
  NC = 19 ##コンパートメント
  ,DT = 0.01 ##時間刻み
  ,T = 1200 ##シミュレーション時間、2000ミリ秒
  
  ,Cm = 3.0  ##膜容量 (μF / cm^2)
  ,Ri = 0.1 ##内部抵抗 (KΩ-cm)
  ,Rm = 10.0 ##膜抵抗 (KΩ-cm^2) (未使用)
  ,Beta = 0.075 ##定数
  
  ,V_LEAK = -60.0 ##リーク電位
)

VALUES <- VALUES %>% mutate(
  .
  ,NT = T / DT
  
  ,V_Na = 115.0 + V_LEAK ##ナトリウムのリーク電位
  ,V_Ca = 140.0 + V_LEAK ##カルシウムのリーク電位
  ,V_K = -15.0 + V_LEAK ##カリウムのリーク電位
)


##ゲート変数生成
GATE_VALUES <- data.frame(
  t(rep(1:12))
)

names(GATE_VALUES) <- c(
  
  "M" ## ナトリウムチャネルの活性化ゲート。
  ,"S" ## カルシウムチャネルの活性化ゲート。
  ,"N" ## 遅延整流カリウムチャネルの活性化ゲート。
  ,"C" ## カルシウム依存性カリウムチャネルの活性化ゲート。
  ,"A" ## A型カリウムチャネルの活性化ゲート。
  ,"H" ## ナトリウムチャネルの不活性化ゲート。
  ,"R" ## カルシウムチャネルの不活性化ゲート。
  ,"B" ## A型カリウムチャネルの不活性化ゲート。
  ,"Q" ## カリウムAHPチャネルに関連する変数。
  ,"V" ##膜電位
  ,"XI" ## 細胞内のカルシウム濃度に関連する変数。
  ,"N_VARS" ##ゲート変数の数 + 1
)

# コンパートメントごとの定数の定義
COMPARTMENT_VALUES <- data.frame(
  g_Na = c(0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 0.0, 15.0, 30.0, 15.0, 0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
  ,g_K_DR = c(0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 0.0, 5.0, 15.0, 5.0, 0.0, 20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
  ,g_K_A = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
  ,g_K_C = c(0.0, 5.0, 5.0, 10.0, 10.0, 10.0, 5.0, 20.0, 10.0, 20.0, 5.0, 15.0, 15.0, 15.0, 15.0, 15.0, 5.0, 5.0, 0.0)
  ,g_K_AHP = c(0.0, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.0)
  ,g_Ca = c(0.0, 5.0, 5.0, 12.0, 12.0, 12.0, 5.0, 8.0, 4.0, 8.0, 5.0, 17.0, 17.0, 17.0, 10.0, 10.0, 5.0, 5.0, 0.0)
  ,g_leak = rep(0.1, VALUES$NC)
  ,phi = c(7769, 7769, 7769, 7769.0, 7769.0, 7769.0, 7769.0, 34530.0, 17402.0, 26404.0, 5941.0, 5941.0, 5941.0, 5941.0, 5941.0, 5941.0, 5941.0, 5941.0, 5941.0)
  ,rad = c(2.89e-4, 2.89e-4, 2.89e-4, 2.89e-4, 2.89e-4, 2.89e-4, 2.89e-4, 2.89e-4, 4.23e-4, 2.42e-4, 2.42e-4, 2.42e-4, 2.42e-4, 2.42e-4, 2.42e-4, 2.42e-4, 2.42e-4, 2.42e-4, 2.42e-4)
  ,len = c(1.20e-2, 1.20e-2, 1.20e-2, 1.20e-2, 1.20e-2, 1.20e-2, 1.20e-2, 1.20e-2, 1.25e-2, 1.10e-2, 1.10e-2, 1.10e-2, 1.10e-2, 1.10e-2, 1.10e-2, 1.10e-2, 1.10e-2, 1.10e-2, 1.10e-2)
  ,area = c(2.188e-5, 2.188e-5, 2.188e-5, 2.188e-5, 2.188e-5, 2.188e-5, 2.188e-5, 2.188e-5, 3.320e-5, 1.673e-5, 1.673e-5, 1.673e-5, 1.673e-5, 1.673e-5, 1.673e-5, 1.673e-5, 1.673e-5, 1.673e-5, 1.673e-5)
)


##円の面積求めたり、alphaの式を設定するなどの基礎的な関数群
FUNCTIONS_BASE <- list(
  pi_r_2 = function(rad) { ##円の面積を求める
    pi * rad * rad
  }
  ,alpha_m =  function(v){
    0.32 * (13.1 - (v - VALUES$V_LEAK)) / 
      (exp((13.1 - (v - VALUES$V_LEAK)) / 4.0) - 1)
  }
  ,alpha_s = function(v) {
    1.6 / 
      (1 + exp(-0.072 * ((v - VALUES$V_LEAK) - 65)))
  }
  ,alpha_n = function(v) {
    0.016 * (35.1 - (v - VALUES$V_LEAK)) / 
      (exp((35.1 - (v - VALUES$V_LEAK)) / 5.0) - 1)
  }
  ,alpha_c = function(v) {
    # ベクトル条件を定義
    condition <- v <= (50 + VALUES$V_LEAK)
    
    # 条件に応じた計算をベクトル処理で行う
    result <- ifelse(
      condition
      ,exp(
        ((v - VALUES$V_LEAK) - 10) / 11.0 -
          ((v - VALUES$V_LEAK) - 6.5) / 27.0
      ) /
        18.975
      ,2 * exp(-((v - VALUES$V_LEAK) - 6.5) / 27.0)
    )
    return(result)
  }
  ,alpha_a = function(v) {
    0.02 * (13.1 - (v - VALUES$V_LEAK)) /
      (exp((13.1 - (v - VALUES$V_LEAK)) / 10.0) - 1)
  }
  ,alpha_h = function(v) {
    0.128 * exp((17 - (v - VALUES$V_LEAK)) / 18.0)
  }
  ,alpha_r = function(v) {
    # ベクトル条件を定義
    condition <- v <= (0 + VALUES$V_LEAK)
    
    # 条件に応じた計算をベクトル処理で行う
    result <- ifelse(
      condition,
      0.005,
      exp(-(v - VALUES$V_LEAK) / 20.0) / 200.0
    )
    return(result)
  }
  ,alpha_b = function(v) {
    0.0016 * exp((-13 - (v - VALUES$V_LEAK)) / 18.0)
  }
  ,alpha_q = function(x) {
    pmin((0.2e-4) * x, 0.01)
  }
  
  ,beta_m = function(v) {
    0.28 * ((v - VALUES$V_LEAK) - 40.1) / 
      (exp(((v - VALUES$V_LEAK) - 40.1) / 5.0) - 1)
  }
  ,beta_s = function(v) {
    0.02 * ((v - VALUES$V_LEAK) - 51.1) /
      (exp(((v - VALUES$V_LEAK) - 51.1) / 5.0) - 1)
  }
  ,beta_n = function(v) {
    0.25 * exp((20 - (v - VALUES$V_LEAK)) / 40.0)
  }
  ,beta_c = function(v) {
    # ベクトル条件を定義
    condition <- v <= (50 + VALUES$V_LEAK)
    
    # 条件に応じた計算をベクトル処理で行う
    result <- ifelse(
      condition
      ,2 * exp(-((v - VALUES$V_LEAK) - 6.5) / 27.0) - FUNCTIONS_BASE$alpha_c(v)
      ,0
    )
    
    return(result)
  }
  ,beta_a = function(v) {
    0.0175 * ((v - VALUES$V_LEAK) - 40.1) / 
      (exp(((v - VALUES$V_LEAK) - 40.1) / 10.0) - 1)
  }
  ,beta_h = function(v) {
    4.0 /
      (1 + exp((40 - (v - VALUES$V_LEAK)) / 5.0))
  }
  ,beta_r = function(v) {
    # ベクトル条件を定義
    condition <- v <= (0 + VALUES$V_LEAK)
    
    # 条件に応じた計算をベクトル処理で行う
    result <- ifelse(
      condition,
      0,
      0.005 - FUNCTIONS_BASE$alpha_r(v)
    )
    
    return(result)
  }
  
  ,beta_b = function(v) {
    0.05 /
      (1 + exp((10.1 - (v - VALUES$V_LEAK)) / 5.0))
  }
  ,beta_q = function(x) {
    0.001
  }
)

FUNCTIONS_TAU_INF <- list(
  inf_m = function(v) {
    FUNCTIONS_BASE$alpha_m(v) / (FUNCTIONS_BASE$alpha_m(v) + FUNCTIONS_BASE$beta_m(v))
  }
  ,inf_s = function(v) {
    FUNCTIONS_BASE$alpha_s(v) / (FUNCTIONS_BASE$alpha_s(v) + FUNCTIONS_BASE$beta_s(v))
  }
  ,inf_n = function(v) {
    FUNCTIONS_BASE$alpha_n(v) / (FUNCTIONS_BASE$alpha_n(v) + FUNCTIONS_BASE$beta_n(v))
  }
  ,inf_c = function(v) {
    FUNCTIONS_BASE$alpha_c(v) / (FUNCTIONS_BASE$alpha_c(v) + FUNCTIONS_BASE$beta_c(v))
  }
  ,inf_a = function(v) {
    FUNCTIONS_BASE$alpha_a(v) / (FUNCTIONS_BASE$alpha_a(v) + FUNCTIONS_BASE$beta_a(v))
  }
  ,inf_h = function(v) {
    FUNCTIONS_BASE$alpha_h(v) / (FUNCTIONS_BASE$alpha_h(v) + FUNCTIONS_BASE$beta_h(v))
  }
  ,inf_r = function(v) {
    FUNCTIONS_BASE$alpha_r(v) / (FUNCTIONS_BASE$alpha_r(v) + FUNCTIONS_BASE$beta_r(v))
  }
  ,inf_b = function(v) {
    FUNCTIONS_BASE$alpha_b(v) / (FUNCTIONS_BASE$alpha_b(v) + FUNCTIONS_BASE$beta_b(v))
  }
  ,inf_q = function(v) {
    FUNCTIONS_BASE$alpha_q(v) / (FUNCTIONS_BASE$alpha_q(v) + FUNCTIONS_BASE$beta_q(v))
  }
  
  ,tau_m = function(v) {
    1 / (FUNCTIONS_BASE$alpha_m(v) + FUNCTIONS_BASE$beta_m(v))
  }
  ,tau_s = function(v) {
    1 / (FUNCTIONS_BASE$alpha_s(v) + FUNCTIONS_BASE$beta_s(v))
  }
  ,tau_n = function(v) {
    1 / (FUNCTIONS_BASE$alpha_n(v) + FUNCTIONS_BASE$beta_n(v))
  }
  ,tau_c = function(v) {
    1 / (FUNCTIONS_BASE$alpha_c(v) + FUNCTIONS_BASE$beta_c(v))
  }
  ,tau_a = function(v) {
    1 / (FUNCTIONS_BASE$alpha_a(v) + FUNCTIONS_BASE$beta_a(v))
  }
  ,tau_h = function(v) {
    1 / (FUNCTIONS_BASE$alpha_h(v) + FUNCTIONS_BASE$beta_h(v))
  }
  ,tau_r = function(v) {
    1 / (FUNCTIONS_BASE$alpha_r(v) + FUNCTIONS_BASE$beta_r(v))
  }
  ,tau_b = function(v) {
    1 / (FUNCTIONS_BASE$alpha_b(v) + FUNCTIONS_BASE$beta_b(v))
  }
  ,tau_q = function(v) {
    1 / (FUNCTIONS_BASE$alpha_q(v) + FUNCTIONS_BASE$beta_q(v))
  }
)




###初期値の設定をします。関数に入れてもいいけど、デバッグ大変なので

EXE_VALUES <- list(
  Vars = data.frame(
    matrix(0, nrow = VALUES$NC, ncol = GATE_VALUES$N_VARS)
  )
  ,I_Inj = rep(0, VALUES$NC)
  ,G_Comp = matrix(0, nrow = VALUES$NC, ncol = 2)  # G_COMPの初期化
  ,I_Comp = rep(0, VALUES$NC)
)


names(EXE_VALUES$Vars) <- names(GATE_VALUES)

EXE_VALUES$Vars$V <- VALUES$V_LEAK

EXE_VALUES$Vars$M <- FUNCTIONS_TAU_INF$inf_m(EXE_VALUES$Vars$V)
EXE_VALUES$Vars$S <- FUNCTIONS_TAU_INF$inf_s(EXE_VALUES$Vars$V)
EXE_VALUES$Vars$N <- FUNCTIONS_TAU_INF$inf_n(EXE_VALUES$Vars$V)
EXE_VALUES$Vars$C <- FUNCTIONS_TAU_INF$inf_c(EXE_VALUES$Vars$V)
EXE_VALUES$Vars$A <- FUNCTIONS_TAU_INF$inf_a(EXE_VALUES$Vars$V)
EXE_VALUES$Vars$H <- FUNCTIONS_TAU_INF$inf_h(EXE_VALUES$Vars$V)
EXE_VALUES$Vars$R <- FUNCTIONS_TAU_INF$inf_r(EXE_VALUES$Vars$V)
EXE_VALUES$Vars$B <- FUNCTIONS_TAU_INF$inf_b(EXE_VALUES$Vars$V)

i_Ca <- COMPARTMENT_VALUES$g_Ca * COMPARTMENT_VALUES$area *
  EXE_VALUES$Vars$S * EXE_VALUES$Vars$S * EXE_VALUES$Vars$R  *
  (EXE_VALUES$Vars$V - VALUES$V_Ca)

EXE_VALUES$Vars$XI <- - i_Ca * COMPARTMENT_VALUES$phi / VALUES$Beta
EXE_VALUES$Vars$Q <- FUNCTIONS_TAU_INF$inf_q(EXE_VALUES$Vars$XI)



for (i in 1:VALUES$NC) {
  EXE_VALUES$G_Comp[i, 1] <- if (i == 1) {
    0
  } else {
    2.0 / (
      (VALUES$Ri * COMPARTMENT_VALUES$len[i - 1]) / 
        FUNCTIONS_BASE$pi_r_2(COMPARTMENT_VALUES$rad[i - 1]) +
        (VALUES$Ri * COMPARTMENT_VALUES$len[i]) / 
        FUNCTIONS_BASE$pi_r_2(COMPARTMENT_VALUES$rad[i])
    )
  }
  EXE_VALUES$G_Comp[i, 2] <- if (i == VALUES$NC) {
    0
  } else {
    2.0 / (
      (VALUES$Ri * COMPARTMENT_VALUES$len[i]) / 
        FUNCTIONS_BASE$pi_r_2(COMPARTMENT_VALUES$rad[i]) +
        (VALUES$Ri * COMPARTMENT_VALUES$len[i + 1]) / 
        FUNCTIONS_BASE$pi_r_2(COMPARTMENT_VALUES$rad[i + 1])
    )
  }
}


###オイラー法で実行していく処理を関数にします。
####デルタのデータたちをリストでまとめます


solve_euler <- function(vars, i_inj, g_comp) {
  ##関数内部でのデバック用
  # vars = EXE_VALUES$Vars
  # i_inj = EXE_VALUES$I_Inj
  # g_comp = EXE_VALUES$G_Comp
  
  DELTA_VALUES <- list(
    Delta_Vars = data.frame(
      matrix(0, nrow = VALUES$NC, ncol = GATE_VALUES$N_VARS)
    )
    ,i_ion = rep(0, VALUES$NC)
    ,i_comp = rep(0, VALUES$NC)
  )
  
  names(DELTA_VALUES$Delta_Vars) <- names(GATE_VALUES)
  
  v <- vars$V
  
  ##各コンダクタンスを設定
  for (var_name in c("M", "S", "N", "C", "A", "H", "R", "B")) {
    DELTA_VALUES$Delta_Vars[[var_name]] <- (VALUES$DT / FUNCTIONS_TAU_INF[[paste0("tau_", tolower(var_name))]](v)) *
      (-vars[[var_name]] + FUNCTIONS_TAU_INF[[paste0("inf_", tolower(var_name))]](v))
  }
  
  DELTA_VALUES$Delta_Vars$Q <- (VALUES$DT / FUNCTIONS_TAU_INF$tau_q(vars$XI)) *
    (-vars$Q + FUNCTIONS_TAU_INF$inf_q(vars$XI))
  
  
  DELTA_VALUES$i_ion <-
    COMPARTMENT_VALUES$g_leak * (v - VALUES$V_LEAK) +
    COMPARTMENT_VALUES$g_Na * vars$M * vars$M * vars$H * (v - VALUES$V_Na) +
    COMPARTMENT_VALUES$g_Ca * vars$S * vars$S * vars$R * (v - VALUES$V_Ca) +
    COMPARTMENT_VALUES$g_K_DR * vars$N * (v - VALUES$V_K) +
    COMPARTMENT_VALUES$g_K_A * vars$A * vars$B * (v - VALUES$V_K) + 
    COMPARTMENT_VALUES$g_K_AHP * vars$Q * (v - VALUES$V_K) +
    COMPARTMENT_VALUES$g_K_C * vars$C * pmin(1, vars$XI / 250.0) * (v - VALUES$V_K)
  
  
  # i_compの計算
  i_comp1 <- rep(0, VALUES$NC)
  i_comp2 <- rep(0, VALUES$NC)
  # # 最初の要素の処理
  for(i in 1:VALUES$NC){
    i_comp1[i] <- ifelse(
      i == 1
      ,0
      ,(g_comp[i, 1] * (v[i - 1] - v[i]) / COMPARTMENT_VALUES$area[i])
    )
    i_comp2[i] <-ifelse(
      i == VALUES$NC
      ,0
      , (g_comp[i, 2] * (v[i + 1] - v[i]) / COMPARTMENT_VALUES$area[i])
    )
    DELTA_VALUES$i_comp[i] <- i_comp1[i] + i_comp2[i]
  }
  
  
  
  DELTA_VALUES$Delta_Vars$V <- (VALUES$DT / VALUES$Cm) * (- DELTA_VALUES$i_ion + DELTA_VALUES$i_comp + i_inj)
  
  i_Ca <-
    COMPARTMENT_VALUES$g_Ca * COMPARTMENT_VALUES$area * vars$S *　vars$S * vars$R *
    (vars$V - VALUES$V_Ca)
  
  DELTA_VALUES$Delta_Vars$XI <-
    VALUES$DT * (- COMPARTMENT_VALUES$phi * i_Ca - VALUES$Beta * vars$XI)
  
  vars <- vars + DELTA_VALUES$Delta_Vars
  
  
  result_list <- list(
    Vars = vars
    ,I_Inj = i_inj
    ,G_Comp = g_comp
    ,I_Comp = DELTA_VALUES$i_comp
    ,I_Ion = DELTA_VALUES$i_ion
  )
  
  return(result_list)
}


DF <- data.frame(
  time = numeric(VALUES$NT)
  ,matrix(0, nrow = VALUES$NT, ncol = VALUES$NC)
)
colnames(DF)[2:(VALUES$NC + 1)] <- paste0("V_", 1:VALUES$NC)

##表に出てこない変数の確認用に書いていたコード
# VAR_DF <- data.frame(
#   time = numeric(VALUES$NT)
#   ,matrix(0, nrow = VALUES$NT, ncol = VALUES$NC)
# )
# colnames(VAR_DF)[2:(VALUES$NC + 1)] <- paste0("V_", 1:VALUES$NC)

RESULT_LIST <- list()

for(nt in 1:VALUES$NT){
  t = VALUES$DT * (nt - 1)
  
  RESULT_LIST[[nt]] <- c(t, EXE_VALUES$Vars$V)
  
  ##表に出てこない変数の確認用に書いていたコード
  # VAR_DF[nt, 1] <- t
  # VAR_DF[nt, 2:(VALUES$NC + 1)] <- EXE_VALUES$I_Comp
  
  EXE_VALUES$I_Inj[9] <- if (t >= 1000 && t < 1001) { 1e-2 / COMPARTMENT_VALUES$area[9] } else { 0 }
  
  EXE_VALUES <- solve_euler(vars = EXE_VALUES$Vars, i_inj = EXE_VALUES$I_Inj, g_comp = EXE_VALUES$G_Comp)
  print(nt)
}

DF <- do.call(rbind, RESULT_LIST)
DF <- data.frame(DF)
names(DF) <- c("time", paste0("V_", rep(1:VALUES$NC)))


DF_LONG <- DF %>% pivot_longer(
  cols = starts_with("V_")
  ,names_to = "compartment"
  ,values_to = "V"
) %>% mutate(
  .
  ,compartment = gsub("V_", "", compartment)
) %>% mutate(
  .
  ,compartment = factor(
    compartment
    ,levels = sort(unique(as.numeric(compartment)))
  )
)


DF_LONG %>% filter(
  .
  ,as.numeric(compartment) >= 9
  ,time >= 980
  ,time <= 1050
  
) %>%  ggplot(
  .
  ,aes(x = time, y = V, colour = compartment)
) + geom_line() 




