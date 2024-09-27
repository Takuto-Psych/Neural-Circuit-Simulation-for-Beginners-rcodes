rm(list = ls())
library(tidyverse)
set.seed(2525)

VALUES <- data.frame(
  n = 4000    ##総ニューロン数
  ,n_e = 0.8 ## 興奮性ニューロン率
  
  ,t = 1000 ##総シミュレーション時間 1000ms
  ,dt = 1.0 ##deltaの時間 1ms
  
  ,tau_m = 20.0 ###膜の時定数(ms)
  ,tau_e = 5.0 ###興奮性シナプスの時定数(ms)
  ,tau_i = 10.0 ###抑制性シナプスの時定数(ms)
  
  ,v_rest = -49.0 ##静止膜電位
  ,v_init = - 60.0 ##初期値
  ,v_reset = - 60.0 ##リセット電位
  ,theta = -50.0　##閾値
  
  ,p = 0.3 ##シナプス結合の形成確立
  
  ,d = 10 ###集団を出すための時間幅の設定
)


VALUES <- VALUES %>% mutate(
  .,
  ,n_e = n * n_e
  ,n_i = n - n_e ## 抑制性ニューロン数
  ,nt = t / dt ##ステップ数
  ,g_e = 1.62 / tau_e ##興奮性シナプスのスパイク一発当たりの後シナプス電位の変化量
  ,g_i = - 9.0 / tau_i ##抑制性シナプスのスパイク一発当たりの後シナプス電位の変化量
  ,nd = t / d ##時間幅の個数の設定
)

initialize <- function() { ##ネットワークの構造体の初期化を行うための関数
  list(
    v = VALUES$v_init + runif(VALUES$n, 0, 10)      # 膜電位
    ,ge = rep(0, VALUES$n)                    # 興奮性シナプス入力
    ,gi = rep(0, VALUES$n)                    # 抑制性シナプス入力
    ,w = matrix( #シナプス重み行列。結合があれば1、なければ0
      runif(VALUES$n) < VALUES$p, nrow = VALUES$n, ncol = VALUES$n
    ) 
    ,s = rep(FALSE, VALUES$n)                  # スパイク状態
  )
}




caluculate_synaptic_input <- function(network) {
  ##興奮と抑制それぞれの重みを分けて行列演算
  w_exc <- network$w[,1:VALUES$n_e] ##興奮性シナプスの行列取り出し
  w_inh <- network$w[,(VALUES$n_e + 1) : VALUES$n] ##抑制性シナプスの行列取り出し
  
  s_exc <- network$s[1:VALUES$n_e] ##興奮性シナプスのスパイク状態取り出し
  s_inh <- network$s[(VALUES$n_e + 1) : VALUES$n] ##抑制性ニューロンのスパイク状態取り出し
  
  ##行列積で各ニューロンの興奮性及び抑制性シナプス入力を計算
  ###スパイクがあっても、重みが0なら伝達されない
  re <- w_exc %*% s_exc
  ri <- w_inh %*% s_inh
  
  ##膜電位の更新
  network$ge <- exp(- VALUES$dt / VALUES$tau_e) * network$ge + re
  network$gi <- exp(- VALUES$dt / VALUES$tau_i) * network$gi + ri
  
  return(network)
}


update_cell_parameters <- function(network) {
  network$v <- network$v +
    VALUES$dt * (
      -(network$v - VALUES$v_rest) +
        VALUES$g_e * network$ge + 
        VALUES$g_i * network$gi
    ) / 
    VALUES$tau_m
  
  network$s <- network$v > VALUES$theta
  
  network$v <- ifelse(
    network$s
    ,VALUES$v_reset
    ,network$v
  )
  return(network)
}




output_spike <- function(nt, network, spike_df) {
  spike_neuron_id <- which(network$s)
  spike_neuron_num <- length(spike_neuron_id)
  
  if(spike_neuron_num > 0) {
    new_spike_df <- data.frame(
      t = rep((nt) * VALUES$dt, spike_neuron_num) 
      ,spike_neuron_id = spike_neuron_id
    )
    
    spike_df <- rbind(spike_df, new_spike_df)
  }
  return(spike_df)
}







###ここから主になるところを実行
network <- initialize()

##結果格納用データフレーム作成
DF <- data.frame(
  t = numeric()
  ,spike_neuron_id = numeric()
)

for(nt in 1:VALUES$nt) {
  network <- caluculate_synaptic_input(network) ##シナプス入力の作成
  network <- update_cell_parameters(network) ##膜電位の更新
  DF <- output_spike(nt, network, DF)
}



DF %>% ggplot(
  .
  ,aes(x = t, y = spike_neuron_id)
) + geom_point(
  alpha = 0.6
  ,color = 'black'
)



##追加分析

####毎時の発火シナプス数

DF_DELTA <- data.frame(
  delta = rep(1:VALUES$nd, each = VALUES$d)
  ,t = 1:VALUES$t
)
DF_DELTA <- DF_DELTA %>%  mutate(
  delta_time = (delta-1) * VALUES$d
)


DF_WITH_DELTA <- left_join(
  x = DF
  ,y = DF_DELTA
  ,by = "t"
)



SUM_DF_DELTA <- DF_WITH_DELTA %>% group_by(
  .
  ,delta_time
) %>% summarise(
  .
  ,spike_num = length(spike_neuron_id)
) %>% mutate(
  .
  ,rate = spike_num / (VALUES$n * VALUES$d) * 1000
)


SUM_DF_DELTA %>% ggplot(
  .
  ,aes(x = delta_time, y = rate)
) + geom_line(
) + geom_ribbon(
  aes(ymin = 0, ymax = rate)
)



SUM_DF_DELTA %>% ggplot(
  .
  ,aes(x = delta_time, y = rate)
) + geom_line(
) + geom_ribbon(
  aes(ymin = 0, ymax = rate)
)

