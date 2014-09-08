##Data frame molecular dynamics
N<-27
T.vector<-seq(0.01,0.1,by=0.01)
K1=10
dt=0.1          ##K1 and dt determine precision of molecular dynamics. dt*K1 should be equal to 10 to provide optimal sampling times.
                ##For K1=100 and dt=0.1 we get average precision of modeling.
                ##For K1=1000 and dt=0.01 precision is better, but all the calculations take much more time. 
K2=10           ##index for number of sampled points for each particle
r<-reinit(N)

data<-matrix(ncol=6)
colnames(data)<-c("Temperature","particle number","potential energy U","total radius","inner shell?","outer shell?")
data<-as.data.frame(data)

i<-1
for(t in T.vector){
        r<-reinit(N)
        r<-temp(r,t)
        for(k in 1:K2){
                for (n in 1:N){
                        data[i,1]<-Temp
                        data[i,2]<-n
                        data[i,3]<-Un(r,n)
                        data[i,4]<-Un
                        r<-molecular(r,K=K1,dt=dt)
                }
        }
}
















##Initialization
r<-reinit(27)
r<-temp(r,0.03)


##Molecular dynamics
K=10
for(i in 1:K){
r<-molecular(r)
        print(E(r))
}


#Molecular dynamics with aggregation:
K=30
aggr<-array(dim=c(6,27,K))
for(i in 1:K){
        r<-molecular(r)
        aggr[,,i]<-r
}
aggr.plot(aggr)


##T=0.0475
##Calculating potential energy distribution and mean
K=100
EnergyU<-NULL
for(i in 1:K){
        r<-molecular(r)
        EnergyU<-c(EnergyU,U(r))
        print(i)
}
hist(EnergyU,breaks=50)
mean(EnergyU)

##T=0.0475
##radial distribution and mean, total and inner/outer shells:
K=100
radius<-matrix(nrow=27,ncol=K)
for(i in 1:K){
        print(i)
        r<-molecular(r)
        radius[,i]<-rad(r)
}
radv<-as.vector(radius)
hist(radv,breaks=50)
rad.inn<-radv[radv<1]
rad.out<-radv[radv>1]
hist(rad.inn,breaks=50)
hist(rad.out,breaks=50)
mean(rad.inn)
mean(rad.out)




##Cycle for many temperatures:
sd=c("0.04","0.05","0.06") ## and so on
for(s in sd){
        r<-reinit(27)
        r<-temp(r,sd)
        ##Molecular dynamics
        K=10
        for(i in 1:K){
                r<-molecular(r)
                print(E(r))
        }
}


