###環境クリーンアップ
rm(list = ls())
library(tidyverse)

### 初期値設定
t <- 0 ##時刻 最初は0秒から
x <- 0 ##位置 最初は原点（0m）から
v <- 0 ##速度 最初は停止（０m/sから）
dt <- 1.0 ##時間の刻み幅 1秒ずつ進める
a <- 1.0 ##加速度 1m/s^2

DF <- cbind(NA, NA) %>% 
  data.frame() %>% na.omit()
names(DF) <- c("x", "t")

while (t < 10.0) {
  TMP_DF <- cbind(x, t) %>% data.frame() 
  DF <- rbind(TMP_DF, DF) 
  
  x <- x + v * dt ##次の時刻の位置を計算
  v <- v + a * dt ##次の時刻の速度を計算
  t <- t + dt ##時間をdt秒進める
}


p <- DF %>% ggplot(
  aes(x = t, y = x)
) + geom_point() + stat_function( ## 解析解を線で引く
  fun = function(num) 0.5 * (num ^2)
)

p

