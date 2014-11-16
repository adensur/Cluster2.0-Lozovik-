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
t=0.01

#
r=reinit(27) #initial equilibrium state, T=0
r=temp(r, temperature = t)
arr=molecular2(r,1000,50,dt=0.01)

vel1=velocity(arr)
hist(vel1,50)

rad1=aggr.rad(arr)
hist(rad1,50)
for(n in 1:N){
        sd
}
#