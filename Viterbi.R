library(MASS)
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(formattable)
#set pi transition vecotr
pi_tran =c(0.2,0.5,0.3)
#set parameters for the gaussian
emission_parm = matrix(c(-20,5,10,5,20,5),ncol=2,byrow = T)
#number of data points to sample
num_data = 100
#generate data
state = 1
observations<-0
for(i in 2:num_data){
uni<-runif(1)
if(runif(1)<pi_tran[1]){
state[i]<-1
observations[i]<-rnorm(n=1,mean = emission_parm[1,1],sd=emission_parm[1,2])
}else if(uni<pi_tran[2]){
state[i]<-2
observations[i]<-rnorm(n=1,mean = emission_parm[2,1],sd=emission_parm[2,2])
}else{
state[i]<-3
observations[i]<-rnorm(n=1,mean = emission_parm[3,1],sd=emission_parm[3,2])
}

}
#plot data
p1<-data.frame(point=1:num_data,observation=observations)%>%
ggplot(aes(x=point,y=observation))+geom_line()+theme_pubclean()

#initalize forward filter matrix
forward_filter = matrix(0,nrow = 3,ncol = num_data)
#intial conditions
forward_filter[1,1] = pnorm(observations[1],emission_parm[1,1],emission_parm[1,2])*.4
forward_filter[2,1] = pnorm(observations[1],emission_parm[2,1],emission_parm[2,2])*.3
forward_filter[3,1] = pnorm(observations[1],emission_parm[3,1],emission_parm[3,2])*.3
#calculate stable foriward filter
C_stable = sum(forward_filter[1:3,1])
forward_filter[1,1] = forward_filter[1,1]/C_stable
forward_filter[2,1] = forward_filter[2,1]/C_stable
forward_filter[3,1] = forward_filter[3,1]/C_stable
#calculate forward filter using stable method
for(i in 2:num_data){
forward_filter[1,i] = dnorm(observations[i],mean = emission_parm[1,1],sd = emission
_parm[1,2])*
sum(forward_filter[1:3,i-1]*pi_tran)
forward_filter[2,i] = dnorm(observations[i],emission_parm[2,1],emission_parm[2,2])*
sum(forward_filter[1:3,i-1]*pi_tran)
forward_filter[3,i] = dnorm(observations[i],emission_parm[3,1],emission_parm[3,2])*
sum(forward_filter[1:3,i-1]*pi_tran)

C_stable = sum(forward_filter[1:3,i])
forward_filter[1,i] = forward_filter[1,i]/C_stable
forward_filter[2,i] = forward_filter[2,i]/C_stable
forward_filter[3,i] = forward_filter[3,i]/C_stable
}
#get most probable state at each point
prob_state=apply(forward_filter,2,which.max)
p2<-data.frame(point=1:num_data,viterbi_path=prob_state,true=state)%>%
melt(id.var="point")%>%mutate(variable=factor(variable,levels = c("true","viterbi_p
ath")))%>%ggplot(aes(x=point,y=value))+geom_line()+facet_wrap(variable~.,ncol = 1)+sc
ale_y_continuous(breaks = 1:2)+ylab("state")+theme_pubclean()

ggarrange(p1,p2,ncol = 1)

table(viterbi_path=apply(forward_filter,2,which.max),true_state=state)

data.frame(gene=1:num_data,expression=c("low","med","high")[apply(forward_filter,2,wh
ich.max)])
