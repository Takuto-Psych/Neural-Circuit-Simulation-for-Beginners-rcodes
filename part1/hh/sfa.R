##環境リセット
rm(list = ls())
library(tidyverse)

###変数設定
e_rest <- -65.0 ##mV
c <- 1.0 ##micro F / cm^2
g_leak <- 0.3 ##ms /cm^2
e_leak <- 10.6 + e_rest ##mV
g_na <- 120.0 ##mS /cm^2
e_na <- 115.0 + e_rest ##mS /cm^2
g_k <- 36.0 ##mS/cm^2
e_k <- -12.0 + e_rest ##mV
i_ext <- 40.0 ##micro A/cm^2
tau_ahp <- 200 ##ms
g_ahp <- 1400 ##ms/cmm^2

dt <- 0.01 ##10 micro s
t <- 1000 ##100ms
nt <- t/dt ## t/dt

VALUES <- cbind(e_rest, c, g_leak, e_leak, g_na, e_na, g_k, e_k, i_ext, dt, t, nt, tau_ahp, g_ahp) %>% 
  data.frame()


##create functions

###define the formula
alpha_m <- function(v) {
  (2.5 - 0.1 * (v - VALUES$e_rest) ) /
  (exp(2.5 - 0.1 * (v - VALUES$e_rest)) - 1.0 ) %>%
    as.double()
}

beta_m <- function(v) {
  4.0 * exp(-(v - VALUES$e_rest) / 18.0) %>% 
    as.double()
}

alpha_h <- function(v) {
  0.07 * exp(-(v - VALUES$e_rest) / 20.0) %>% 
    as.double()
}


beta_h <- function(v) {
  1.0 /
    (exp(3.0 - 0.1 * (v - VALUES$e_rest)) + 1.0 ) %>%
    as.double()
}

alpha_n <- function(v) {
  (0.1 - 0.01 * (v - VALUES$e_rest)) /
    (exp(1.0 - 0.1 * (v - VALUES$e_rest)) - 1.0 ) %>%
    as.double()
}

beta_n <- function(v) {
  0.125 * exp(- (v - VALUES$e_rest) / 80.0) %>% 
    as.double()
}


## how to set the first value
m0 <- function(v) {
  alpha_m(v) / (alpha_m(v) + beta_m(v) ) %>% 
    as.double()
}
  
h0 <- function(v) {
  alpha_h(v) / (alpha_h(v) + beta_h(v) ) %>% 
    as.double()
}

n0 <- function(v) {
  alpha_n(v) / (alpha_n(v) + beta_n(v) ) %>% 
    as.double()
}

tau_m <- function(v) {
  1 / (alpha_m(v) + beta_m(v) ) %>% 
    as.double()
}
  
tau_h <- function(v) {
  1 / (alpha_h(v) + beta_h(v) ) %>% 
    as.double()
}

tau_n <- function(v) {
  1 / (alpha_n(v) + beta_n(v) ) %>% 
    as.double()
}


### how to update the opening ratio of channel

dmdt <- function(v, m) {
   (1.0 / tau_m(v)) * (-m + m0(v)) %>% 
    as.double()
}

dhdt <- function(v, h) {
  (1.0 / tau_h(v)) * (-h + h0(v)) %>% 
    as.double()
}

dndt <- function(v, n) {
  (1.0 / tau_n(v)) * (-n + n0(v)) %>% 
    as.double()
}

dadt <- function(s, a){
  (1.0 / VALUES$tau_ahp) * (-a + s) %>% 
    as.double()
}


dvdt <- function(v, m, h, n, a, i_ext) {
  (- VALUES$g_leak * (v - VALUES$e_leak) -
     VALUES$g_na * m * m * m * h * (v - VALUES$e_na) -
     VALUES$g_k * n * n * n * n * (v - VALUES$e_k) -
     VALUES$g_ahp * a * (v - VALUES$e_k) +  ##i_ahp
     VALUES$i_ext) /
    VALUES$c %>% as.double() 
}

###before for rooping, set the first value

v <- VALUES$e_rest
m <- m0(v)
h <- h0(v)
n <- n0(v)
a <- 0.0
v_old <- VALUES$e_rest
DF <- data.frame(t = numeric(VALUES$nt), v = numeric(VALUES$nt), m = numeric(VALUES$nt),
                 h = numeric(VALUES$nt), n = numeric(VALUES$nt), a = numeric(VALUES$nt))




###for roop the updating
i <- 0
for(i in 1:nt) {
  t = dt * i
  DF[i,] <- c(t, v, m, h, n, a)
  
  s <- ifelse(
    v > 0 & v_old < 0,
    1,
    0
  )
  v_old <- v
  
  dmdt1 <- dmdt(v, m)
  dhdt1 <- dhdt(v, h)
  dndt1 <- dndt(v, n)
  dadt1 <- dadt(s, a)
  dvdt1 <- dvdt(v, m, h, n, a, VALUES$i_ext)
  
  dmdt2 <- dmdt(v + 0.5 * VALUES$dt * dvdt1, m + 0.5 * VALUES$dt * dmdt1)
  dhdt2 <- dhdt(v + 0.5 * VALUES$dt * dvdt1, h + 0.5 * VALUES$dt * dhdt1)
  dndt2 <- dndt(v + 0.5 * VALUES$dt * dvdt1, n + 0.5 * VALUES$dt * dndt1)
  dadt2 <- dadt(s                          , a + 0.5 * VALUES$dt * dadt1)
  dvdt2 <- dvdt(v + 0.5 * VALUES$dt * dvdt1, m + 0.5 * VALUES$dt * dmdt1, 
                h + 0.5 * VALUES$dt * dhdt1, n + 0.5 * VALUES$dt * dndt1,
                a + 0.5 * VALUES$dt * dadt1, VALUES$i_ext)
  
  dmdt3 <- dmdt(v + 0.5 * VALUES$dt * dvdt2, m + 0.5 * VALUES$dt * dmdt2)
  dhdt3 <- dhdt(v + 0.5 * VALUES$dt * dvdt2, h + 0.5 * VALUES$dt * dhdt2)
  dndt3 <- dndt(v + 0.5 * VALUES$dt * dvdt2, n + 0.5 * VALUES$dt * dndt2)
  dadt3 <- dadt(s                          , a + 0.5 * VALUES$dt * dadt2)
  dvdt3 <- dvdt(v + 0.5 * VALUES$dt * dvdt2, m + 0.5 * VALUES$dt * dmdt2, 
                h + 0.5 * VALUES$dt * dhdt2, n + 0.5 * VALUES$dt * dndt2,
                a + 0.5 * VALUES$dt * dadt2, VALUES$i_ext)
  
  dmdt4 <- dmdt(v + VALUES$dt * dvdt3, m + VALUES$dt * dmdt3)
  dhdt4 <- dhdt(v + VALUES$dt * dvdt3, h + VALUES$dt * dhdt3)
  dndt4 <- dndt(v + VALUES$dt * dvdt3, n + VALUES$dt * dndt3)
  dadt4 <- dadt(s                    , a + VALUES$dt * dadt3)
  dvdt4 <- dvdt(v + VALUES$dt * dvdt3, m + VALUES$dt * dmdt3,
                h + VALUES$dt * dhdt3, n + VALUES$dt * dndt3,
                a + VALUES$dt * dadt3, VALUES$i_ext)  
  
  m <- m + VALUES$dt * (dmdt1 + 2 * dmdt2 + 2 * dmdt3 + dmdt4) / 6
  h <- h + VALUES$dt * (dhdt1 + 2 * dhdt2 + 2 * dhdt3 + dhdt4) / 6
  n <- n + VALUES$dt * (dndt1 + 2 * dndt2 + 2 * dndt3 + dndt4) / 6
  a <- a + VALUES$dt * (dadt1 + 2 * dadt2 + 2 * dadt3 + dadt4) / 6
  v <- v + VALUES$dt * (dvdt1 + 2 * dvdt2 + 2 * dvdt3 + dvdt4) / 6
  
  ### im waiting for 'for' you but tell me if you are not stacking
  print(i)
  
}

sfa_spike <- DF %>% ggplot(
  .,
  aes(x = t, y = v)
) + geom_line()
ggsave("part1/hh/sfa_spike.png", sfa_spike, dpi = 320)
sfa_spike




DF %>% filter(
  .,
  t <= 100
) %>% ggplot(
  .,
  aes(x = t, y = v)
) + geom_line()


DF %>% filter(
  .,
  t >= 900
) %>% ggplot(
  .,
  aes(x = t, y = v)
) + geom_line()
