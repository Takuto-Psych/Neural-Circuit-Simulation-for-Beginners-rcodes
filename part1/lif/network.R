rm(list = ls())
library(tidyverse)
VALUES <- data.frame(
  tau = 20.0, ## ms
  v_rest = -65.0, ## mV
  v_reset = -65.0, ## mV
  theta = -55.0, ## mV
  r_m = 1.0, ##MOhm
  dt = 1.0, ## ms
  t = 1000.0, ##ms
  i_ext = 12.0, ##nA
  n = 2,
  tau_syn = 5.0, ##ms
  r_syn = 1.0, ##M0hm
  w = 2.0  ###結合強度
)

VALUES <-VALUES %>% mutate(
  .,
  nt = t / dt, ##times
)

DF <- data.frame(
  t = numeric(), 
  v1 = numeric(), v2 = numeric(), 
  s1 = numeric(), s2 = numeric() 
)


v <- c(VALUES$v_rest, VALUES$v_rest - 15 )
i_syn <- c(0.0, 0.0)
s <- c(0, 0)


### function for the array updating



for(i in 1:VALUES$nt){
  t <- VALUES$dt * i
  DF <- rbind(
    DF,
    data.frame(
      t = t, 
      v1 = v[1], v2 = v[2], 
      s1 = s[1], s2 = s[2]
    )
  )
  i_syn <- exp( - VALUES$dt / VALUES$tau_syn) * 
    i_syn + VALUES$w * s[rep(2:1, length.out = VALUES$n)]
  
  v <- v + 
    VALUES$dt *
    (
      -(v - VALUES$v_rest) +
        VALUES$r_syn * i_syn + 
        VALUES$r_m * VALUES$i_ext
    ) / 
    VALUES$tau
  
  s <- v > VALUES$theta
  
  
  if(s[1]){
    DF <- rbind(
      DF,
      data.frame(
        t = t + VALUES$dt,
        v1 = v[1], v2 = v[2],
        s1 = s[1], s2 = s[2]
      )
    )
    DF <- rbind(
      DF,
      data.frame(
        t = t + VALUES$dt,
        v1 = 0.0, v2 = v[2],
        s1 = s[1], s2 = s[2]
      )
    )
  }
  if(s[2]){
    DF <- rbind(
      DF,
      data.frame(
        t = t + VALUES$dt,
        v1 = v[1], v2 = v[2],
        s1 = s[1], s2 = s[2]
      )
    )
    DF <- rbind(
      DF,
      data.frame(
        t = t + VALUES$dt,
        v1 = v[1], v2 = 0,
        s1 = s[1], s2 = s[2]
      )
    )
  }
  
  v <- ifelse(s, VALUES$v_reset, v)
}

DF_LONG <- DF %>% pivot_longer(
  .,
  cols = c(v1, v2, s1, s2),
  names_to = c(".value", "neuron_label"),
  names_pattern = "([vs])(\\d)"
)

network_spike <- DF_LONG %>% ggplot(
  .,
  aes(x = t, y = v, colour = factor(neuron_label))
) + geom_line()
network_spike
