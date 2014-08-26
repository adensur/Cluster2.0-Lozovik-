gradient.descent<-function(N,r,alfa=1,K=10,print=FALSE){
        if(missing(r)){
                r<-init(N)
                print("missing r. Reinitializing randoms")
        }
        else{
                N<-ncol(r)
        }
        Break<-FALSE
        for(i in 1:K){
                string<-paste(c("i=",i,"of 100000", " U=",U(r)),collapse="")
                if(print){print(string)}
                deltaU<- 1    ##difference between current U and calculated on next iteration U
                t<-1          ##counter to exclude the possibility of infinite loop
                while(deltaU>=0 & !Break){
                        r2<-delta(r,alfa=alfa)
                        deltaU<-U(r2)-U(r)
                        ##print(c("delta=U",deltaU))
                        if(deltaU>0){
                                alfa<-alfa/10
                                ##print(c("alfa=",alfa))
                        }
                        t<=t+1
                        if(t>=100000){Break<-TRUE}##to avoid infinite loop
                        break
                }
                if(Break){
                        print("Too many iterations. exiting loop")
                        break
                }
                r<-r2
        }
        print(U(r))
        r
}
grad.U<-function(r,var,k){##takes numeric vectors x,y,z, and number of particles N
        ##var is the variable, for which the gradient is calculated.var=1 - x, var=2 - y, var=3 - z
        ##k is the index of the variable in vector. (i.e., number of particle)
        ##so for var="x" and k=3 the gradient U by x[3] will be calculated.
        N<-ncol(r)      ##number of particles
        if(k>N){print("grad.U: invalid k. k is more then number of particles")}
        sum<-0
        for(i in 1:N){
                if(i!=k){
                        sum<-sum+6*(r[var,k]-r[var,i])*(rki(r,k,i)^(-8))
                }
        }
        sum<-2*r[var,k]-sum
        sum
}
vgrad.U<-function(r){##based on grad.U. returns the 3*N vector of gradient
        N<-ncol(r)
        null<-rep(0,times=N)
        grad<-rbind(null,null,null)##the initial value for vector grad is all zeroes; 3*N matrix
        for(k in 1:N){
                for(var in 1:3){
                        grad[var,k]<-grad.U(r,var,k)
                }
        }
        grad
}
rki<-function(r,k,i){##calculates the distance between particle i and particle k; 
        ##x, y, z (vectors) specify coordinates of all particles
        rki<-NULL
        if(k==i){print("rki error! k=i!")}
        else{
                rki<-((r[1,k]-r[1,i])^2+(r[2,k]-r[2,i])^2+(r[3,k]-r[3,i])^2)^0.5
        }
        names(rki)<-"Distance"
        rki
}
U<-function(r){##returns a value of the potential energy of the system
        N<-ncol(r)
        U<-0
        for(i in 1:N){
                sum<-0
                if(i<N){
                        for(k in (i+1):N){
                                sum<-sum+rki(r,k,i)^(-6)
                        }
                }
                U<-U+sum+r[1,i]^2+r[2,i]^2+r[3,i]^2
        }
        names(U)<-"Potential energy"
        U
}
init<-function(N){##initializes matrix 3 times N with correct rownames and random values
        r<-rbind(rnorm(N),rnorm(N),rnorm(N))
        rownames(r)<-(c("x","y","z"))
        r
}
delta<-function(r,alfa=1){##calculate vector of difference (one gradient descent iteration)
        r2<-r
        N<-ncol(r)
        for(k in 1:N){
                for(var in 1:3){
                        r2[var,k]<-r[var,k]-alfa*grad.U(r,var,k)
                }
        }
        r2
}
reinit<-function(N){##loads r from file; returns it
        file<-paste(c("data_init/",N,".csv"),collapse="")
        r<-as.matrix(read.csv(file))
        r
}
myplot<-function(r,...){##rad is a vector of shell radiuses. Particles within radk, rad(k+1) will be drown in same color
        library(rgl)
        plot3d(r[1,],r[2,],r[3,],...)##дописать функцию, чтобы рисовала частицы разных оболочек разными цветами
}##extends possibility of plot3d to plot a matrix 3xN
descent<-function(N=1:100,print=FALSE,alfa=1){##calcs descent over vector of N's and write each to a file
        for(i in N){
                r<-gradient.descent(i,alfa=alfa,K=100000, print=print)
                file<-paste(c("data_init/",i,".csv"),collapse="")
                write.csv(r,file=file,row.names=FALSE)
                print(paste(c("printed file",file)),collapse="")
        }
        string<-paste(c("files for N=",paste(N,collapse=", ")," printed succesfully"), collapse="")
        print(string)
}
rad<-function(r){##calculates the distance from particle k to the beginning of the coordinates
        N<-ncol(r)
        rad<-rep(0,times=N)
        for(k in 1:N){
                rad[k]<-(r[1,k]^2+r[2,k]^2+r[3,k]^2)^(0.5)
        }
        rad
}
temp<-function(r,T=0.001){##adds a small, random increments to the coordinates (r) and velocities of the system
        N<-ncol(r)
        dr<-rbind(rnorm(N,sd=T),rnorm(N,sd=T),rnorm(N,sd=T)) ##for now T is just the standart
                                                                ##deviation of the new distribution
        dv<-rbind(rnorm(N,sd=T),rnorm(N,sd=T),rnorm(N,sd=T))
        r<-r+dr
        r<-rbind(r,dv)
}
rstep<-function(r, dt=1){##returns a next step of the "leap frog" iteration process
        N<-ncol(r)
        deltar<-dt*r[4:6,]##this is now the matrix 3*N of small random increments
        null<-rep(0,times=N)
        null<-rbind(null,null,null)##null is now a matrix 3*N of all zeroes
        deltar<-rbind(deltar,null)
        rnew<-r+deltar
        rnew
}
vstep<-function(r,dt=1){##same as rstep, but for velocities v
        N<-ncol(r)
        deltav<-dt*vgrad.U(r[1:3,])
        null<-rep(0,times=N)
        null<-rbind(null,null,null)##null is now a matrix 3*N of all zeroes
        deltav<-rbind(null,deltav)
        vnew<-r-deltav
        vnew
}
array.descent<-function(N=27,M=20,sd=1,alfa=0.5,K=20000,print=FALSE){
        ##calculates a series of grad. descents to compare local minimums
        arr<-array(dim=c(3,27,M))
        for(m in 1:M){
                r<-rbind(rnorm(N,sd=sd),rnorm(N,sd=sd),rnorm(N,sd=sd))
                arr[,,m]<-gradient.descent(N=N,r=r,alfa=alfa,K=K,print=print)
        }
        arr
}
array.U<-function(arr,M){
        ##calculates a vector of potential energies from the array of 3*N*M M different r vectors.
        vector.U<-NULL
        for(m in 1:M){
                vector.U<-c(vector.U,U(arr[,,m]))
        }
}
myplot2<-function(r,neightbours=5,...){
        t<-0
        plot3d(r[1,],r[2,],r[3,],...)
        library(rgl)
        N<-ncol(r)
        if(neightbours>N){
                neightbours<-N
        }
        for(i in 1:N){
                plotted<-0
                lim<-0.1
                while(plotted<=neightbours & t<10000){
                        for(k in i:N){
                                t<-t+1
                                if(i!=k){
                                        r2<-r[,c(i,k)]
                                        a<-r[,i]
                                        b<-r[,k]
                                        theta <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )
                                        ##this calculates angle between vectors
                                        if(theta<lim){
                                                plotted<-plotted+1
                                                plot3d(r2[1,],r2[2,],r2[3,],type="l",add=TRUE,...)
                                        }
                                }
                        }
                        lim<-lim+0.1
                }
        }
}
molecular<-function(r,K,dt=0.1,print=FALSE,plot=FALSE){
        add=FALSE
        for (i in 1:K){
                if(print)print(U(r))
                if(plot)myplot(r,add=add)
                add<-TRUE
                r<-rstep(r,dt)
                r<-vstep(r,dt)
        }
        if(print)print(U(r))
        r
}
##r<-gradient.descent(N=27,r=r,alfa=0.5,K=5000, print = TRUE)
##r<-reinit(N)
##ra<-rad(r)
##plot(sort(ra))
##myplot(r)
##U(r)