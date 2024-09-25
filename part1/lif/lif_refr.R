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
  t_refr = 5.0 ## ms
)
VALUES <- VALUES %>% mutate(
  .,
  nt = t / dt, ##times
  nt_refr = t_refr / dt
)

DF <- data.frame(t = numeric(), v = numeric(), s = numeric() )

v <- VALUES$v_rest
refr <- 0 ##couter for refractory period

for(i in 0:(VALUES$nt - 1)){
  t <- VALUES$dt * i
  DF <- rbind(
    DF,
    data.frame(t = t, v = v, s = 0)
  )
  v <- v + 
    VALUES$dt * (-(v - VALUES$v_rest) + VALUES$r_m * VALUES$i_ext) / VALUES$tau
  
  s <- v > VALUES$theta
  
  
  
  if(s){
    DF <- rbind(
      DF,
      data.frame(t = t + VALUES$dt, v = v, s = 1)
    )
    DF <- rbind(
      DF,
      data.frame(t = t + VALUES$dt, v = 0.0, s = 1)
    )
    DF <- rbind(
      DF,
      data.frame(t = t + VALUES$dt, v = VALUES$v_reset, s = 1)
    )
  }
  
  refr <- ifelse(s, VALUES$nt_refr, refr - 1)
  
  v <- ifelse(refr > 0, VALUES$v_reset, v)
  print(refr)
}


lif_spike <- DF %>% ggplot(
  .,
  aes(x = t, y = v)
) + geom_line()
lif_spike
