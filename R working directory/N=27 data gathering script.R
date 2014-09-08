##Data frame molecular dynamics
N<-27
T.vector<-seq(0.01,0.1,by=0.01)
K1=1000
dt=0.01          ##K1 and dt determine precision of molecular dynamics. dt*K1 should be equal to 10 to provide optimal sampling times.
                ##For K1=100 and dt=0.1 we get average precision of modeling.
                ##For K1=1000 and dt=0.01 precision is better, but all the calculations take much more time. 
K2=1000           ##index for number of sampled points for each particle

data<-matrix(ncol=5)
colnames(data)<-c("Temperature","Particle_number","Potential_energy","Total_radius","Inner/outer shell")
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
                        data[i,4]<-radn(r,n)
                        if(radn(r,n)<1)data[i,5]<-as.factor("inner shell")
                        else data[i,5]<-("outer shell")
                        write.table(data[i,],file="/results/Macroscopic values, N=27/data.txt", rownames=FALSE,
                                    append=TRUE)
                        r<-molecular(r,K=K1,dt=dt)
                        i=i+1
                }
                print(k)
        }
        print(c("T=",t)
}











