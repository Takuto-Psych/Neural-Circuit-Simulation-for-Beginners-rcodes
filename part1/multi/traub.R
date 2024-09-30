rm(list = ls())

library(tidyverse)


VALUES <- data.frame(
  NC = 19 ##コンパートメント
  ,DT = 0.01 ##時間刻み
  ,T = 200 ##シミュレーション時間、200ミリ秒
  
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




##プログラムの実行に関数群。配列処理できるようにして軽めに作ることを意識。
FUNCTIONS_EXE <- list(
  initialize = function() {
    ##外部電流の初期値
    I_inj <- rep(0, VALUES$NC)
    
    ##各コンポートメント（19個）のゲート変数（12(実際は11)）VARSの初期化
    VARS <- data.frame(
      matrix(0, nrow = VALUES$NC, ncol = GATE_VALUES$N_VARS)
    )
    names(VARS) <- names(GATE_VALUES)
    
    VARS$V <- VALUES$V_LEAK
    
    VARS$M <- FUNCTIONS_TAU_INF$inf_m(VARS$V)
    VARS$S <- FUNCTIONS_TAU_INF$inf_a(VARS$V)
    VARS$N <- FUNCTIONS_TAU_INF$inf_n(VARS$V)
    VARS$C <- FUNCTIONS_TAU_INF$inf_c(VARS$V)
    VARS$A <- FUNCTIONS_TAU_INF$inf_a(VARS$V)
    VARS$H <- FUNCTIONS_TAU_INF$inf_h(VARS$V)
    VARS$R <- FUNCTIONS_TAU_INF$inf_r(VARS$V)
    VARS$B <- FUNCTIONS_TAU_INF$inf_b(VARS$V)
    
    
    i_Ca <- COMPARTMENT_VALUES$g_Ca *
      COMPARTMENT_VALUES$area * VARS$S * VARS$R *
      (VARS$V - VALUES$V_Ca)
    
    VARS$XI <- - i_Ca * COMPARTMENT_VALUES$phi / VALUES$Beta
    VARS$Q <- FUNCTIONS_TAU_INF$inf_q(VARS$XI)
    
    
  
    G_COMP <- matrix(0, nrow = VALUES$NC, ncol = 2)  # G_COMPの初期化
    
    # iが1のときは0、そうでなければ計算を行う
    G_COMP[1, 1] <- 0
    G_COMP[2:VALUES$NC, 1] <- sapply(2:VALUES$NC, function(i) {
      2.0 / (
        (VALUES$Ri * COMPARTMENT_VALUES$len[i - 1]) / (pi * COMPARTMENT_VALUES$rad[i - 1]^2) + 
          (VALUES$Ri * COMPARTMENT_VALUES$len[i]) / (pi * COMPARTMENT_VALUES$rad[i]^2)
      )
    })
    
    # iがNCのときは0、そうでなければ計算を行う
    G_COMP[VALUES$NC, 2] <- 0
    G_COMP[1:(VALUES$NC - 1), 2] <- sapply(1:(VALUES$NC - 1), function(i) {
      2.0 / (
        (VALUES$Ri * COMPARTMENT_VALUES$len[i]) / (pi * COMPARTMENT_VALUES$rad[i]^2) + 
          (VALUES$Ri * COMPARTMENT_VALUES$len[i + 1]) / (pi * COMPARTMENT_VALUES$rad[i + 1]^2)
      )
    })
    
    int_list <- list(
      VARS
      ,I_inj
      ,G_COMP  
    )
    
    names(int_list) <- c("VARS", "I_Inj", "G_COMP")
    
    return(int_list)
  }
  
  ,solve_euler = function(vars, i_inj, g_comp){
    vars = EXE_VALUES$VARS
    i_inj = EXE_VALUES$I_inj
    g_comp = EXE_VALUES$G_COMP

    DVARS <- data.frame(
      matrix(0, nrow = VALUES$NC, ncol = GATE_VALUES$N_VARS)
    )
    names(DVARS) <- names(GATE_VALUES)
    
    i_ion <- rep(0, VALUES$NC)
    i_comp <- rep(0, VALUES$NC)
    
    ##各コンポートメント（19個）のゲート変数（12(実際は11)）VARSの更新
    v <- vars$V
    
    DVARS$M <- (VALUES$DT / FUNCTIONS_TAU_INF$tau_m(DVARS$V)) * 
      (- vars$V + FUNCTIONS_TAU_INF$inf_m(v))
    DVARS$S <- (VALUES$DT / FUNCTIONS_TAU_INF$tau_s(DVARS$V)) * 
      (- vars$V + FUNCTIONS_TAU_INF$inf_s(DVARS$V))
    DVARS$N <- (VALUES$DT / FUNCTIONS_TAU_INF$tau_n(DVARS$V)) * 
      (- vars$V + FUNCTIONS_TAU_INF$inf_n(DVARS$V))
    DVARS$C <- (VALUES$DT / FUNCTIONS_TAU_INF$tau_c(DVARS$V)) * 
      (- vars$V + FUNCTIONS_TAU_INF$inf_c(DVARS$V))
    DVARS$A <- (VALUES$DT / FUNCTIONS_TAU_INF$tau_a(DVARS$V)) * 
      (- vars$V + FUNCTIONS_TAU_INF$inf_a(DVARS$V))
    DVARS$H <- (VALUES$DT / FUNCTIONS_TAU_INF$tau_h(DVARS$V)) * 
      (- vars$V + FUNCTIONS_TAU_INF$inf_h(DVARS$V))
    DVARS$R <- (VALUES$DT / FUNCTIONS_TAU_INF$tau_r(DVARS$V)) * 
      (- vars$V + FUNCTIONS_TAU_INF$inf_r(DVARS$V))
    DVARS$B <- (VALUES$DT / FUNCTIONS_TAU_INF$tau_b(DVARS$V)) * 
      (- vars$V + FUNCTIONS_TAU_INF$inf_b(DVARS$V))
    
    DVARS$Q <- (VALUES$DT / FUNCTIONS_TAU_INF$tau_q(vars$XI)) *
      (- vars$V + FUNCTIONS_TAU_INF$inf_q(vars$XI))
    
    
    ##イオンの更新
    i_ion <- 
      (COMPARTMENT_VALUES$g_leak + DVARS$V - VALUES$V_LEAK) +
      COMPARTMENT_VALUES$g_Na * vars$M * vars$M * vars$H *(DVARS$V - VALUES$V_Na) +
      COMPARTMENT_VALUES$g_Ca * vars$S * vars$S * vars$R * (DVARS$V - VALUES$V_Ca) +
      COMPARTMENT_VALUES$g_K_DR * vars$N * (DVARS$V - VALUES$V_K) + 
      COMPARTMENT_VALUES$g_K_A * vars$A * vars$B *(DVARS$V - VALUES$V_K) + 
      COMPARTMENT_VALUES$g_K_AHP * vars$Q * (DVARS$V - VALUES$V_K) + 
      COMPARTMENT_VALUES$g_K_A * vars$C * pmin(vars$XI / 250) * (DVARS$V - VALUES$V_K)
    
    
    # 条件に応じて計算
    i_comp[1] <- 0  # iが1のときは0
    
    # 2からNC-1までの計算
    i_comp[2:(VALUES$NC - 1)] <- 
      g_comp[2:(VALUES$NC - 1), 1] * (DVARS$V[1:(VALUES$NC - 2)] - DVARS$V[2:(VALUES$NC - 1)]) / COMPARTMENT_VALUES$area[2:(VALUES$NC - 1)] +
      g_comp[2:(VALUES$NC - 1), 2] * (DVARS$V[3:VALUES$NC] - DVARS$V[2:(VALUES$NC - 1)]) / COMPARTMENT_VALUES$area[2:(VALUES$NC - 1)]
    
    # NCのときは0
    i_comp[VALUES$NC] <- 0
    
    DVARS$V <- (VALUES$DT / VALUES$Cm) * (- i_ion + i_comp + i_inj) 
    
    i_Ca <- COMPARTMENT_VALUES$g_Ca * COMPARTMENT_VALUES$area * vars$S * vars$R * (vars$V - VALUES$V_Ca)
    
    DVARS$XI <- DVARS$XI + 
      VALUES$DT * (COMPARTMENT_VALUES$phi * i_Ca - VALUES$Beta * vars$XI)
    
    vars <- vars + DVARS
    
    list_result <- list(
      vars
      ,i_inj
      ,g_comp
    )
    names(list_result) <- c("VARS", "I_Inj", "G_COMP")
    return(list_result)
  }
)




EXE_VALUES <- FUNCTIONS_EXE$initialize()

EXE_VALUES <- FUNCTIONS_EXE$solve_euler(EXE_VALUES$VARS, EXE_VALUES$I_Inj, EXE_VALUES$G_COMP)
EXE_VALUES


DF <- data.frame(
  time = numeric(VALUES$NT)
  ,matrix(0, nrow = VALUES$NT, ncol = VALUES$NC)
)
colnames(DF)[2:(VALUES$NC + 1)] <- paste0("V_", 1:VALUES$NC)


EXE_VALUES <- FUNCTIONS_EXE$solve_euler(vars = EXE_VALUES$VARS, i_inj = EXE_VALUES$I_inj, g_comp = EXE_VALUES$G_COMP)

for(nt in 1:VALUES$NT){
  t = VALUES$DT * (nt - 1)
  
  DF[nt, 1] <- t
  DF[nt, 2:(VALUES$NC + 1)] <- EXE_VALUES$VARS$V
  
  # i_injの8番目の要素に値を設定
  EXE_VALUES$I_inj[8] <- 1e-4 / COMPARTMENT_VALUES$area[8]
  # 条件に基づいて注入電流を設定
  # EXE_VALUES$I_inj[8] <- if (t >= 1000 && t < 1001) { 1e-2 / COMPARTMENT_VALUES$area[8] } else { 0 }
  
  EXE_VALUES <- FUNCTIONS_EXE$solve_euler(vars = EXE_VALUES$VARS, i_inj = EXE_VALUES$I_Inj, g_comp = EXE_VALUES$G_COMP)
  
  print(nt)
}










