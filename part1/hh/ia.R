##環境リセット
rm(list = ls())
library(tidyverse)

###変数設定
c <- 1.0 ##micro F / cm^2
##i_ext <- 9.0 ## mirco A /cm^2
i_ext <- 8.60941453 ## mirco A /cm^2
g_leak <- 0.3 ##ms /cm^2
e_leak <- -17.0
g_na <- 120.0 ##mS /cm^2
e_na <- 55.0 ##mS /cm^2
g_k <- 20.0 ##mS/cm^2
e_k <- -75.0 ##mV
g_a <- 47.7 ##mS/cm^2
e_a <- -75.0 ##mV
mshft <- -5.3
hshft <- -12.0
nshft <- -4.3


dt <- 0.01 ##10 micro s
t <- 1000 ##100ms
nt <- t/dt ## t/dt

VALUES <- cbind(c, g_leak, e_leak, g_na, e_na, g_k, e_k, g_a, e_a,
                i_ext, dt, t, nt,
                mshft, hshft, nshft) %>% data.frame()


##create functions

###define the formula
alpha_m <- function(v) {
  (- 0.1 * (v + 35.0 + VALUES$mshft) ) /
    (exp(-0.1 * (v + 35.0 + VALUES$mshft)) - 1.0 ) %>%
    as.double()
}

beta_m <- function(v) {
  4.0 * exp(-(v + 60.0 + VALUES$mshft) / 18.0) %>% 
    as.double()
}

alpha_h <- function(v) {
  0.07 * exp(-(v + 60.0 + VALUES$hshft) / 20.0) %>% 
    as.double()
}


beta_h <- function(v) {
  1.0 /
    (exp(-0.1 * (v + 30.0 + VALUES$hshft)) + 1.0 ) %>%
    as.double()
}


alpha_n <- function(v) {
  (- 0.01 * (v + 50.0 + VALUES$nshft)) /
    (exp(-0.1 * (v + 50.0 + VALUES$nshft)) - 1.0 ) %>%
    as.double()
}

beta_n <- function(v) {
  0.125 * exp(- (v + 60.0 + VALUES$nshft) / 80.0) %>% 
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

a0 <- function(v) {
  ((0.0761 * exp((v + 94.22) / 31.84) / 
             (1 + exp((v + 1.17) / 28.93)))
   ) ^ 0.33333 %>% as.double()
}

b0 <- function(v) {
  1.0 / ((1 + exp ((v + 53.3) / 14.54)) ^ 4) %>% as.double()
}

tau_a <- function(v) {
  0.3632 + 1.158 / ( 1. + exp (( v + 55.96 ) / 20.12))
}

tau_b <- function(v) {
  1.24   + 2.678 / ( 1. + exp ((v + 50) / 16.027)) 
}

tau_m <- function(v) {
  (1.0 / 3.8) * (1 / (alpha_m(v) + beta_m(v))) %>% 
    as.double()
}

tau_h <- function(v) {
  (1.0 / 3.8) * (1 / (alpha_h(v) + beta_h(v))) %>% 
    as.double()
}

tau_n <- function(v) {
  (2.0) * (1 / (alpha_n(v) + beta_n(v))) %>% 
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

dadt <- function(v, a){
  (1.0 / tau_a(v)) * (-a + a0(v)) %>% 
    as.double()
}

dbdt <- function(v, b){
  (1.0 / tau_b(v)) * (-b + b0(v)) %>% 
    as.double()
}


dvdt <- function(v, m, h, n, a, b, i_ext) {
  (- VALUES$g_leak * (v - VALUES$e_leak) -
     VALUES$g_na * m * m * m * h * (v - VALUES$e_na) -
     VALUES$g_k * n * n * n * n * (v - VALUES$e_k) -
     VALUES$g_a * a * a * a * b * (v - VALUES$e_a) +  ##I_A
     VALUES$i_ext) /
    VALUES$c %>% as.double() 
}

###before for rooping, set the first value

v <- VALUES$e_leak
m <- m0(v)
h <- h0(v)
n <- n0(v)
a <- a0(v)
b <- b0(v)
DF <- data.frame(t = numeric(VALUES$nt), v = numeric(VALUES$nt), m = numeric(VALUES$nt),
                 h = numeric(VALUES$nt), n = numeric(VALUES$nt), a = numeric(VALUES$nt),
                 b = numeric(VALUES$nt))




###for roop the updating

for(i in 1:nt) {
  t = VALUES$dt * i
  DF[i,] <- c(t, v, m, h, n, a, b)

  
  dmdt1 <- dmdt(v, m)
  dhdt1 <- dhdt(v, h)
  dndt1 <- dndt(v, n)
  dadt1 <- dadt(v, a)
  dbdt1 <- dbdt(v, b)
  dvdt1 <- dvdt(v, m, h, n, a, b, VALUES$i_ext)
  
  dmdt2 <- dmdt(v + 0.5 * VALUES$dt * dvdt1, m + 0.5 * VALUES$dt * dmdt1)
  dhdt2 <- dhdt(v + 0.5 * VALUES$dt * dvdt1, h + 0.5 * VALUES$dt * dhdt1)
  dndt2 <- dndt(v + 0.5 * VALUES$dt * dvdt1, n + 0.5 * VALUES$dt * dndt1)
  dadt2 <- dadt(v + 0.5 * VALUES$dt * dvdt1, a + 0.5 * VALUES$dt * dadt1)
  dbdt2 <- dbdt(v + 0.5 * VALUES$dt * dvdt1, b + 0.5 * VALUES$dt * dbdt1)
  dvdt2 <- dvdt(v + 0.5 * VALUES$dt * dvdt1, m + 0.5 * VALUES$dt * dmdt1, 
                h + 0.5 * VALUES$dt * dhdt1, n + 0.5 * VALUES$dt * dndt1,
                a + 0.5 * VALUES$dt * dadt1, b + 0.5 * VALUES$dt * dbdt1,
                VALUES$i_ext)
  
  dmdt3 <- dmdt(v + 0.5 * VALUES$dt * dvdt2, m + 0.5 * VALUES$dt * dmdt2)
  dhdt3 <- dhdt(v + 0.5 * VALUES$dt * dvdt2, h + 0.5 * VALUES$dt * dhdt2)
  dndt3 <- dndt(v + 0.5 * VALUES$dt * dvdt2, n + 0.5 * VALUES$dt * dndt2)
  dadt3 <- dadt(v + 0.5 * VALUES$dt * dvdt2, a + 0.5 * VALUES$dt * dadt2)
  dbdt3 <- dbdt(v + 0.5 * VALUES$dt * dvdt2, b + 0.5 * VALUES$dt * dbdt2)
  dvdt3 <- dvdt(v + 0.5 * VALUES$dt * dvdt2, m + 0.5 * VALUES$dt * dmdt2, 
                h + 0.5 * VALUES$dt * dhdt2, n + 0.5 * VALUES$dt * dndt2,
                a + 0.5 * VALUES$dt * dadt2, b + 0.5 * VALUES$dt * dbdt2,
                VALUES$i_ext)
  
  dmdt4 <- dmdt(v + VALUES$dt * dvdt3, m + VALUES$dt * dmdt3)
  dhdt4 <- dhdt(v + VALUES$dt * dvdt3, h + VALUES$dt * dhdt3)
  dndt4 <- dndt(v + VALUES$dt * dvdt3, n + VALUES$dt * dndt3)
  dadt4 <- dadt(v + VALUES$dt * dvdt3, a + VALUES$dt * dadt3)
  dbdt4 <- dbdt(v + VALUES$dt * dvdt3, b + VALUES$dt * dbdt3)
  dvdt4 <- dvdt(v + VALUES$dt * dvdt3, m + VALUES$dt * dmdt3,
                h + VALUES$dt * dhdt3, n + VALUES$dt * dndt3,
                a + VALUES$dt * dadt3, b + VALUES$dt * dbdt3,
                VALUES$i_ext)  
  
  m <- m + VALUES$dt * (dmdt1 + 2 * dmdt2 + 2 * dmdt3 + dmdt4) / 6
  h <- h + VALUES$dt * (dhdt1 + 2 * dhdt2 + 2 * dhdt3 + dhdt4) / 6
  n <- n + VALUES$dt * (dndt1 + 2 * dndt2 + 2 * dndt3 + dndt4) / 6
  a <- a + VALUES$dt * (dadt1 + 2 * dadt2 + 2 * dadt3 + dadt4) / 6
  b <- b + VALUES$dt * (dbdt1 + 2 * dbdt2 + 2 * dbdt3 + dbdt4) / 6
  v <- v + VALUES$dt * (dvdt1 + 2 * dvdt2 + 2 * dvdt3 + dvdt4) / 6
  
  ### im waiting for 'for' you but tell me if you are not stacking
  print(i)
  
}

ia_spike <- DF %>% ggplot(
  .,
  aes(x = t, y = v)
) + geom_line()
ggsave("part1/hh/ia_spike.png", ia_spike, dpi = 320)
ia_spike



