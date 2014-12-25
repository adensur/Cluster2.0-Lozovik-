##mean squared radial deviation

##initialization
N=27         #number of particles
data=data.frame()
velocity=function(aggr){   #function calculates velocity distribution over all the array. 
        M=dim(aggr)[3]
        N=dim(aggr)[2]
        vel=NULL
        for(m in 1:M){
                
                vel=c(vel,rad(aggr[4:6,,m]))
                
        }
        vel
}
aggr.rad=function(aggr){
        M=dim(aggr)[3]
        radial=NULL
        for(m in 1:M){
                radial=cbind(radial,rad(aggr[,,m]))
        }
        radial
}

##setting the desired standard deviation points, in which we calculate sd(rad)
t=seq(0.001,0.2,by=0.003)
K1=500
K2=15
dt=0.05
#

##calculate!
for(i in 1:length(t)){
r=reinit(N) #initial equilibrium state, T=0
r=temp(r, temperature = t[i])    ##setting temperature
arr=molecular2(r,K1,K2,dt=dt)    ##calculating an array of coordinates/velocities
rad1=aggr.rad(arr)               ##calculating mean radius
sd1=NULL
for(n in 1:N){
        sd1=c(sd1,sd(rad1[n,]))
}
#
sd10=mean(sd1)
vector=c(t[i],Temp,mean(sd1))
data=rbind(data,vector)
#print()
}
colnames(data)=c("sd","temp","rad")


