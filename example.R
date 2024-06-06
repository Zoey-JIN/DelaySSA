####with delay####
####birth death####
# 0->N=>0 k1=10 tao=5

# order: N
sample <- 1000
k <- c(10)
S_matrix <- c(1)
S_matrix <- matrix(S_matrix,nrow = 1) 
S_matrix_delay <- c(-1)
S_matrix_delay <- matrix(S_matrix_delay,nrow = 1)
tmax <- 20
n_initial <- matrix(c(0),nrow = 1)
t_initial <- 0
delay_type <- matrix(c(2),nrow = 1)
delaytime_list <- list()
delaytime_list <- append(delaytime_list,5)
product_matrix <- matrix(c(0),nrow = 1)


####oscillation####
# 0->X =>Y k1=1/(1+Y^2) tau=20
# Y->0 k2=1/(1+Y)

# order: X Y
sample <- 1000
S_matrix <- c(1,0,0,-1)
S_matrix <- matrix(S_matrix,nrow = 2) 
S_matrix_delay <- c(-1,1,0,0)
S_matrix_delay <- matrix(S_matrix_delay,nrow = 2)
tmax <- 400
n_initial <- matrix(c(0,0),nrow = 2)
t_initial <- 0
delay_type <- matrix(c(2,0),nrow = 1)
delaytime_list <- list()
delaytime_list <- append(delaytime_list,20)
delaytime_list <- append(delaytime_list,0)
product_matrix <- matrix(c(0,0,0,1),nrow = 2)

k <- function(n){  k <- c(1/(1 + (n[2])^2), 1/(1 + n[2]))}


####bursty####
# 0->iN=>0 k=ab^i/(1+b)^(i+1) tau=120
# a = 0.0282
# b = 3.46
# i = 1,...,30

# order: N
sample <- 10000
a <- 0.0282
b <- 3.46
j <- 30
k <- c(sapply(1:j, function(i) a * b^i / (1 + b)^(i + 1)))
S_matrix <- c(1:j)
S_matrix <- matrix(S_matrix,nrow = 1) 
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

product_matrix <- matrix(rep(0,j),nrow = 1)



####on off/ Telegraph ####
#off->on k1=0.03
#on->off k2=0.04
#on->on+N N=>0 k3=10 tau=5

#order: off on N
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

product_matrix <- matrix(c(1,0,0,0,1,0,0,1,0),nrow = 3)


####sir####
# S+I->E+I E=>I k1=0.0001 tau=20
# I->R k2=0.02

# order: S E I R
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
product_matrix <- matrix(c(1,0,1,0,0,0,1,0),nrow = 4)



####delay interrupt####
# 0->N=>0 k1=10 tau=5
# N->0 k2=1

# order: N
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
product_matrix <- matrix(c(0,1),nrow = 1)

####delay interrupt example2####
# 0->A k1=2
# A->0 k2=0.1
# A->I=>0 k3=0.5 tau=15
# I->0 k4=0.1

# order: A I
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
product_matrix <- matrix(c(0,0,1,0,1,0,0,1),nrow = 2)


####without delay####
####birth death####
# 0->N k1=30
# N->0 k2=1 

# order: N
sample <- 10000
k <- c(30,1)
S_matrix <- c(1, -1)
S_matrix <- matrix(S_matrix,nrow = 1)
tmax <- 10
n_initial <- matrix(c(0),nrow = 1)
t_initial <- 0
product_matrix <- matrix(c(0,1),nrow = 1)

####bursty####
# 0->iN k=ab^i/(1+b)^(i+1) 
# N->0 k2=1/120
# a = 0.0282
# b = 3.46
# i = 1,...,40

# order: N
sample <- 10000
j=40
a=0.0282
b=3.46
k <- c(sapply(1:j, function(i) a * b^i / (1 + b)^(i + 1)),1/120)
S_matrix=c(1:j,-1)
S_matrix <- matrix(S_matrix,nrow = 1)
tmax <- 1000
t_initial <- 0
n_initial <- matrix(0,nrow = 1)
n_initial <- matrix(n_initial,nrow = 1)
product_matrix <- matrix(c(rep(0,j),1),nrow = 1)

####on off / Telegraph####
# off->on k1=0.03
# on->off k2=0.04
# on->on+N k3=10 
# N=>0 k4=1

# order: off on N
sample <- 10000
k <- c(0.03,0.04,10,1)
S_matrix=c(-1,1,0,1,-1,0,0,0,1,0,0,-1)
S_matrix <- matrix(S_matrix,nrow = 3)
tmax <- 100
t_initial <- 0
n_initial <- c(1,0,0)
n_initial <- matrix(n_initial,nrow = 3)
product_matrix <- matrix(c(1,0,0,0,1,0,0,1,0,0,0,1),nrow = 3)