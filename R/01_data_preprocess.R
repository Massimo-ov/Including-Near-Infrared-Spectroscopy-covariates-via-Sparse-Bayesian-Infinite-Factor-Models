#BISCUIT NIRS DATASET, OSBORNE 1984
library(fds)

data(labp)
data(labc)
data(nirp)
data(nirc)
df_train    <- t(nirc$y)
df_train <- df_train[-23,]
df_test <- t(nirp$y)
#nirs p reduction like paper BROWN 1380-2400 (WEST 1202 -2400)
start_wl <- 1380
end_wl   <- 2400
i0 <- (start_wl - 1100) / 2 + 1
i1 <- (end_wl   - 1100) / 2 + 1
step_cols <- 4 / 2
keep_idx <- seq(from = i0, to = i1, by = step_cols)
df_train <- df_train[, keep_idx]
df_test <- df_test[, keep_idx]

labels <- nirc$x
labels <- labels[keep_idx] 

#preprocess, Dipendent <- 1=Fat 2=Sugar 3=Flour 4=Water
Dip <- scale(labc[1,])
df_train<- scale(df_train,center=TRUE, scale= TRUE)
df_train <- cbind(Dip[-23,],df_train)

y_test<-scale(labp[1,])
df_test <-  scale(df_test, center= TRUE, scale = TRUE)
df_test <- cbind(y_test,df_test)
