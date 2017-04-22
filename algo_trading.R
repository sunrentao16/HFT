# algorithm trading project
# based on Trading algorithm with learning in latent alpha models

# read data------------------------------
library(R.matlab)
data<-readMat("Data/FB_20141103.mat")
#data<-readMat("FB_20141103.mat")
# df<-as.data.frame(data$data)
dt<-data$data
listName=c('Event','SellVolume','SellPrice',
           'BuyVolume','BuyPrice')
#
sapply(seq_along(listName),
       function(i){
         eval(parse(text=paste( listName[i],'<<- dt[[',i,']]' ,sep="")))
       })
#
# time----------------------------------
#op<- options(digits.secs=3)
#init=strptime("2014-11-03 00:00:00.001", "%Y-%m-%d %H:%M:%OS")
market_start_time = 34200000 
market_end_time = 57600000

market_start_rowN = 0 # find row number
market_end_rowN = 0
for( i in 1:nrow(Event)){
  if(Event[i,1] <=  market_start_time){
    market_start_rowN = i
  }else if( Event[i,1] <=  market_end_time){
    market_end_rowN = i
  }else {break}
}

# merge data---------------------------
midprice = ( SellPrice[market_start_rowN : market_end_rowN,1] 
          + BuyPrice[market_start_rowN : market_end_rowN,1] ) /2/10000 # convert to dollar unit
data_merged = cbind(Event[market_start_rowN : market_end_rowN ,]
                    ,midprice)
data_market_order = data_merged[data_merged[,3] == 69 | data_merged[,3] == 70 ,]
# modify unit of data
# data_market_order[,4] = data_market_order[,4] / 100
# convert time to [0,1]
# data_market_order[,1] = (data_market_order[,1] - market_start_time) / (market_end_time - market_start_time)

# cumulated market order
N_p = rep(0, nrow(data_market_order))
N_m = rep(0, nrow(data_market_order))
for(i in 2:nrow(data_market_order)){
  if(data_market_order[i,7] == 0){
    N_p[i] = N_p[i-1] + data_market_order[i,4]
    N_m[i] = N_m[i-1]
  }else{
    N_p[i] = N_p[i-1]
    N_m[i] = N_m[i-1] + data_market_order[i,4]
  }
}
### model: S - S_0= b*(N_P - N_M)
df=data.frame(x=(N_p-N_m)/10000,y=data_market_order[,8])
model_midprice=glm(y~x, data=df, family=gaussian)
#
b = model_midprice$coefficients[2]
S_0 = model_midprice$coefficients[1]
#
plot(df$x,df$y)
lines(df$x, fitted(model_midprice), col="red")

### handle Poisson Process-------------------------
lambda_p = rep(0, length(N_p))
lambda_m = rep(0, length(N_m))
#
step_size = 100

for(i in step_size:length(N_p)){
    #lambda_p[i] = N_p[i] / (data_market_order[i,1] - market_start_time)
    lambda_p[i] = (N_p[i] - N_p[(i-step_size+1)])/ (data_market_order[i,1] - data_market_order[i-step_size+1,1])
    lambda_m[i] = (N_m[i] - N_m[(i-step_size+1)])/ (data_market_order[i,1] - data_market_order[i-step_size+1,1])
}
# handle infinte
lambda_p[is.infinite(lambda_p)] = max(lambda_p[is.finite(lambda_p)])
lambda_m[is.infinite(lambda_m)] = max(lambda_m[is.finite(lambda_m)])

sum(is.infinite(lambda_p))
sum(is.infinite(lambda_m))
# set theta values and x
theta_0 = mean(data_market_order[,8])
x = theta_0 - data_market_order[,8]
x_p=x; x_m= -1*x;
x_p[x_p < 0] = 0
x_m[x_m < 0] = 0
#plot(x_p, lambda_p)

# linear regress lambda; find mu and kappa
df_l_p = data.frame(x_p,y=lambda_p)
model_lambda_p = glm(y~x, data=df_l_p, family=gaussian)

df_l_m = data.frame(x_m,y=lambda_m)
model_lambda_m = glm(y~x, data=df_l_m, family=gaussian)

# got mu and kappa
mu = (model_lambda_p$coefficients[1] + model_lambda_m$coefficients[1])/2
kappa = (model_lambda_p$coefficients[2] + model_lambda_m$coefficients[2])/2

### Theta ----------------------------------------
hist(midprice) # interesting observation
# values of theta
theta_1 = mean(midprice[midprice < mean(midprice)])
theta_2 = mean(midprice[midprice > mean(midprice)])
theta=c(theta_1, theta_2)
# transition matrix C
C = matrix(c(.9,0.2, 0.1, 0.8), nrow=2, ncol=2)

### estimate Porbability of Theta-------------------------------------
p = matrix(0, nrow=2, ncol=length(N_p))
# initial
p[,1] = 0.5
# set time increment
delta_t = rep(1, nrow(data_market_order))
delta_t[2:length(delta_t)] = data_market_order[2:nrow(data_market_order),1]-data_market_order[1:nrow(data_market_order)-1,1]
delta_t[delta_t==0] = 1
# estimate unconditional probability
for(k in 2:length(N_p)){
  for(j in 1:2){
    p[j,k] = p[j,k-1] * exp(2*(1-mu-0.5*kappa*abs(theta[j] - midprice[k-1]))*delta_t[k] 
                            +  C[j, ,drop=FALSE] %*% p[,k-1,drop=FALSE]/p[j,k-1]*delta_t[k] 
                            + (N_p[k]-N_p[k-1] ) * log( mu+kappa*(theta[j]-midprice[k-1])*((theta[j]-midprice[k-1])>0) )  
                            + (N_m[k]-N_m[k-1] ) * log( mu+kappa*(theta[j]-midprice[k-1])*((theta[j]-midprice[k-1])<0)*(-1) )
                            )
  }
}

# log Probability
p_log = matrix(0, nrow=2, ncol=length(N_p))
# initial
p_log[,1] = log(0.5)
# estimate unconditional probability
for(k in 2:length(N_p)){
  for(j in 1:2){
    p_log[j,k] = p_log[j,k-1] + (2*(1-mu-0.5*kappa*abs(theta[j] - midprice[k-1]))*delta_t[k] 
                            + (N_p[k]-N_p[k-1] ) * log( mu+kappa*(theta[j]-midprice[k-1])*((theta[j]-midprice[k-1])>0) )  
                            + (N_m[k]-N_m[k-1] ) * log( mu+kappa*(theta[j]-midprice[k-1])*((theta[j]-midprice[k-1])<0)*(-1) )
    )
  }
}

### trading strategy------------------------------
v_ac = rep(0, length(N_p))
Q = rep(0, length(N_p))
#
phi = 2.5e-5
a = 1e-5
gamma = sqrt(phi/a)

for(t in 1:length(N_p)){
  v_ac[t] = -1 * gamma * tanh(gamma*(market_start_time - data_market_order[t,1])) * Q[t]
  Q[t+1] = Q[t] + v_ac[t]
}
