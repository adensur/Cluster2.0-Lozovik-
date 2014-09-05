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


