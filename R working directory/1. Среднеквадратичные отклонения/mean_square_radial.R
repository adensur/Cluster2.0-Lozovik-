##mean squared radial deviation
N=27         #number of particles
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
t=c(0.0001,0.001,0.01,0.1,0.2)

data=matrix(nrow=0,ncol=2)
colnames(data)=c("temperature","sd")
#
for(i in 1:length(t)){
r=reinit(27) #initial equilibrium state, T=0
r=temp(r, temperature = t[i])
arr=molecular2(r,1000,50,dt=0.01)
rad1=aggr.rad(arr)
hist(rad1,50)
sd1=NULL
for(n in 1:N){
        sd1=c(sd1,sd(rad1[n,]))
}
#
sd10=mean(sd1)
vector=c(Temp,mean(sd1))
data=rbind(data,vector)
}
