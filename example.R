####有delay####
####birth death####
# N数量
sample <- 1000
k <- c(10)
S_matrix <- c(1)
S_matrix <- matrix(S_matrix,nrow = 1) #按照列
S_matrix_delay <- -c(1)
S_matrix_delay <- matrix(S_matrix_delay,nrow = 1)
tmax <- 20
n_initial <- matrix(c(0),nrow = 1)
t_initial <- 0
delay_type <- matrix(c(2),nrow = 1)
delaytime_list <- list()
delaytime_list <- append(delaytime_list,5)
product_list <- matrix(c(0),nrow = 1)


####oscillation####
# X数量 Y数量
# 0->x ->y
# y->0
sample <- 1000
# k <- c(1/(1+n[2]),1/(1+n[2]))
S_matrix <- c(1,0,0,-1)
S_matrix <- matrix(S_matrix,nrow = 2) #按照列
S_matrix_delay <- c(-1,1,0,0)
S_matrix_delay <- matrix(S_matrix_delay,nrow = 2)
tmax <- 400
n_initial <- matrix(c(0,0),nrow = 2)
t_initial <- 0
delay_type <- matrix(c(2,0),nrow = 1)
delaytime_list <- list()
delaytime_list <- append(delaytime_list,20)
delaytime_list <- append(delaytime_list,0)
product_list <- matrix(c(0,0,0,1),nrow = 2)

k <- function(n){  k <- c(1/(1 + (n[2])^2), 1/(1 + n[2]))}


####bursty####
sample <- 10000
a <- 0.0282
b <- 3.46
j <- 30
k <- c(sapply(1:j, function(i) a * b^i / (1 + b)^(i + 1)))
S_matrix <- c(1:j)
S_matrix <- matrix(S_matrix,nrow = 1) #按照列
S_matrix_delay <- -c(1:j)
S_matrix_delay <- matrix(S_matrix_delay,nrow = 1)
tmax <- 200
n_initial <- matrix(c(0),nrow = 1)
t_initial <- 0
delay_type <- matrix(rep(c(2),times=j),nrow = 1)
delaytime_list <- list()
for (i in 1:j) {
  delaytime_list <- append(delaytime_list,120) 
}

product_list <- matrix(rep(0,j),nrow = 3)



####on off####
#off on p
sample <- 10000
k <- c(0.03,0.04,10)
S_matrix=c(-1,1,0,1,-1,0,0,0,1)
S_matrix <- matrix(S_matrix,nrow = 3)
S_matrix_delay <- c(0,0,0,0,0,0,0,0,-1)
S_matrix_delay <- matrix(S_matrix_delay,nrow = 3)
tmax <- 100
n_initial <- matrix(c(1,0,0),nrow = 3)
t_initial <- 0
delay_type <- matrix(c(0,0,2),nrow = 1)
delaytime_list <- list()
delaytime_list <- append(delaytime_list,0)
delaytime_list <- append(delaytime_list,0)
delaytime_list <- append(delaytime_list,5)

product_list <- matrix(c(1,0,0,0,1,0,0,1,0),nrow = 3)


####sir####
# S E I R
sample <- 5
k <- c(1e-4,1e-2)
S_matrix=c(-1,1,0,0,0,0,-1,1)
S_matrix <- matrix(S_matrix,nrow = 4)
S_matrix_delay <- c(0,-1,1,0,0,0,0,0)
S_matrix_delay <- matrix(S_matrix_delay,nrow = 4)
tmax <- 400
n_initial <- matrix(c(999,0,1,0),nrow = 4)
t_initial <- 0
delay_type <- matrix(c(2,0),nrow = 1)
delaytime_list <- list()
delaytime_list <- append(delaytime_list,20)
delaytime_list <- append(delaytime_list,0)
product_list <- matrix(c(1,0,1,0,0,0,1,0),nrow = 4)



####delay interrupt####
sample <- 5000
k <- c(10,1)
S_matrix=c(1,-1)
S_matrix <- matrix(S_matrix,nrow = 1)
S_matrix_delay <- c(-1,0)
S_matrix_delay <- matrix(S_matrix_delay,nrow = 1)
tmax <- 50
n_initial <- matrix(c(0),nrow = 1)
t_initial <- 0
delay_type <- matrix(c(2,0),nrow = 1)
delaytime_list <- list()
delaytime_list <- append(delaytime_list,5)
delaytime_list <- append(delaytime_list,0)
product_list <- matrix(c(0,1),nrow = 1)

####delay interrupt example2####
sample <- 10000
k <- c(2,0.1,0.5,0.1)
S_matrix=c(1,0,-1,0,-1,1,0,-1)
S_matrix <- matrix(S_matrix,nrow = 2)
S_matrix_delay <- c(0,0,0,0,0,-1,0,0)
S_matrix_delay <- matrix(S_matrix_delay,nrow = 2)
tmax <- 50
n_initial <- matrix(c(0,0),nrow = 2)
t_initial <- 0
delay_type <- matrix(c(0,0,2,0),nrow = 1)
delaytime_list <- list()
delaytime_list <- append(delaytime_list,0)
delaytime_list <- append(delaytime_list,0)
delaytime_list <- append(delaytime_list,15)
delaytime_list <- append(delaytime_list,0)
product_list <- matrix(c(0,0,1,0,1,0,0,1),nrow = 2)


####没有delay####
####test1####
sample <- 10000
k <- c(30,1)
S_matrix <- c(1, -1)
S_matrix <- matrix(S_matrix,nrow = 1)
tmax <- 10
n_initial <- 0
t_initial <- 0
fun_fr <- function(k,n){
  return(k*c(1,n))
}

####test2####
sample <- 10000
j=40
a=0.0282
b=3.46
k <- c(sapply(1:j, function(i) a * b^i / (1 + b)^(i + 1)),1/120)
S_matrix=c(1:j,-1)
S_matrix <- matrix(S_matrix,nrow = 1)
tmax <- 1000
t_initial <- 0
n_initial <- c(0)
n_initial <- matrix(n_initial,nrow = 1)
fun_fr <- function(k,n){
  return(k*c(rep(1,j),n))
}

####test3####
# g* g n 
sample <- 10000
k <- c(0.03,0.04,10,1)
S_matrix=c(-1,1,0,0,1,-1,0,0,0,0,1,-1)
S_matrix <- matrix(S_matrix,nrow = 3, byrow = TRUE)
tmax <- 100
t_initial <- 0
n_initial <- c(1,0,0)
n_initial <- matrix(n_initial,nrow = 3)
fun_fr <- function(k,n){
  return(k*c(n[1],n[2],n[2],n[3]))
}
