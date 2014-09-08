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
        grad<-matrix(nrow=3,ncol=N)
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
Un<-function(r,n){
        ##returns potential energy of the particle #n
        N<-ncol(r)
        U<-0
        for(i in 1:N){
                if(i!=n){
                        U<-U+rki(r,n,i)^(-6)
                }
        }
        U<-U/2+r[1,n]^2+r[2,n]^2+r[3,n]^2
        names(U)<-"Potential energy"
        U
}
Kinn<-function(r,n){
        ##calculates kinetic energy of the particle n in the system specified by r.
        sum<-r[4,n]^2+r[5,n]^2+r[6,n]^2
        names(sum)<-"Kinetic energy"
        sum
}
E<-function(r){
        N<-dim(r)[2]
        E<-0
        for(n in 1:N){
                E<-E+r[4,n]^2+r[5,n]^2+r[6,n]^2
        }
        E+U(r)
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
reinit<-function(N,type="N",K,dt,T){##loads r from file; returns it
        if(type=="N"){
                print("type=N, getting energy minimum configurations")
                file<-paste(c("data_init/",N,".csv"),collapse="")
                r<-as.matrix(read.csv(file))
                r
        }
        else if(type=="arr"){
                file<-paste(c("data_init/arr",N,".csv"),collapse="")
                arr2<-read.csv(file=file) ##a data frame 3 by (N*M)
                L<-dim(arr2)[2]
                M<-L/N
                arr<-array(dim=c(3,N,M))
                for(m in 1:M){
                        for(n in 1:N){
                                arr[,n,m]<-arr2[,(N*(m-1)+n)]
                        }
                }
                arr
        }
        else if(type=="r.aggr"){
                print("type = r.aggr specified. getting array data for r aggregations")
                file<-paste(c("data_molecular/r.aggr_N=",N,"_K=",K,"_dt=",dt,"_T=",T,".csv"),collapse="")
                aggr2<-read.csv(file=file) ##data frame of 6*N*K
                L<-dim(aggr2)[2]
                K2<-L/N
                if(K2!=K){stop("wrong K in the dimension of file. aborting")}
                else{
                        agrr<-array(dim=c(3,N,K))
                        for(k in 1:K){
                                for(n in 1:N){
                                        aggr[,n,k]<-aggr2[,(N*(k-1)+n)]
                                }
                        }
                        agrr
                }
        }
}
myplot<-function(r,...){##rad is a vector of shell radiuses. Particles within radk, rad(k+1) will be drown in same color
        library(rgl)
        plot3d(r[1,],r[2,],r[3,],...)##äîïèñàòü ôóíêöèþ, ÷òîáû ðèñîâàëà ÷àñòèöû ðàçíûõ îáîëî÷åê ðàçíûìè öâåòàìè
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
radn<-function(r,n){##calculates mean radius of particle n in the cluster specified by r
        (r[1,n]^2+r[2,n]^2+r[3,n]^2)^(0.5)
}
temp<-function(r,temperature=0.001){##adds a small, random increments to the coordinates (r)
        ##T is actually just a standart deviation
        ##true "T" should be retrieved as a difference in energy before and after adding
        ##increments
        N<-ncol(r)
        U1<-U(r)
        dr<-rbind(rnorm(N,sd=temperature),rnorm(N,sd=temperature),rnorm(N,sd=temperature)) 
        r<-r+dr
        r<-rbind(r,0,0,0)
        Temp<<-U(r)-U1
        r
}
rstep<-function(r, dt=1){##returns a next step of the "leap frog" iteration process
        r[1:3,]<-r[1:3,]+dt*r[4:6,]
        r
}
vstep<-function(r,dt=1){##same as rstep, but for velocities v
        r[4:6,]<-r[4:6,]-dt*vgrad.U(r[1:3,])
        r
}
array.descent<-function(N=27,M=20,sd=1,alfa=0.5,K=20000,print=FALSE){
        ##calculates a series of grad. descents to compare local minimums
        arr<-array(dim=c(3,N,M))
        for(m in 1:M){
                r<-rbind(rnorm(N,sd=sd),rnorm(N,sd=sd),rnorm(N,sd=sd))
                arr[,,m]<-gradient.descent(N=N,r=r,alfa=alfa,K=K,print=print)
        }
        file<-paste(c("data_init/arr",N,".csv"),collapse="")
        write.csv(arr,file=file,row.names=FALSE)
        arr
}
array.U<-function(arr,M){
        ##calculates a vector of potential energies from the array of 3*N*M M different r vectors.
        vector.U<-NULL
        for(m in 1:M){
                vector.U<-c(vector.U,U(arr[,,m]))
        }
        vector.U
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
molecular<-function(r=r,K=100,dt=0.1){
        ##simplified function to make K molecular dynamic steps;
        ##r is a standart object 6*N provided by the temp() function
        for(k in 1:K){
                r<-rstep(r,dt)
                r<-vstep(r,dt)
        }
        r
}
molecular2<-function(r,K1,K2,dt=0.1,fun="r.aggregate"){
        N=ncol(r)
        if(fun=="r.aggregate"){
                print("specified r=r.aggregate. aggregating array")
                r.aggr<-array(dim=c(6,N,K1))
                for(k1 in 1:K1){
                        for(k2 in 1:K2){
                                r<-rstep(r,dt)
                                r<-vstep(r,dt)
                        }
                        aggr[,,k1]<-r
                        print(paste(c("k1=",k1,"E=",E(r))))
                }
                return(r.aggr)
        }
}
mywrite<-function(object,type,N,K,dt,T){
        ##writes into a file an object, class of which is specified by "type"
        if(type=="r.aggr"){
                print("type=r.aggr specified. Writing into file array of r's, marking N, K, dt and T in the filename")
                file<-paste(c("data_molecular/r.aggr_N=",N,"_K=",K,"_dt=",dt,"_T=",T,".csv"),collapse="")
                write.csv(object,file=file,row.names=FALSE)
        }
}
ro.perp<-function(x,a){
        ##calculates the distance between point M(x,y,z) and line s(a,b,c), where a is a vector of length 3 (a,b,c), and x 
        ##is a vector of length 3 x,y,z
        sqrt((x[2]*a[3]-a[2]*x[3])^2+(x[3]*a[1]-x[1]*a[3])^2+(x[1]*a[2]-x[2]*a[1])^2)/sqrt(a[1]^2+a[2]^2+a[3]^2)
}
ro.perp.aggrn<-function(aggrn,a){
        ##calculates the ro.perp function between a and vectors
        ##x, specified by aggrn[,k]ths elements
        ##aggrn is an array of 3*K
        K<-dim(aggrn)[2]
        ro<-NULL
        for(i in 1:K){
                ro<-c(ro,ro.perp(aggrn[,i],a))
        }
        ro
}
ro.perp.NK<-function(aggr,line){
        ##calculates perpendicular distance from each point in an array aggr
        ##returns N*K matrix
        N=dim(aggr)[2]
        K=dim(aggr)[3]
        matrix.roperp<-matrix(nrow=N,ncol=K)
        for(n in 1:N){
                for(k in 1:K){
                        matrix.roperp[n,k]<-ro.perp(aggr[1:3,n,k],line)
                }
        }
        matrix.roperp
}
ro.par.NK<-function(aggr,line){
        ##calculates perpendicular distance from each point in an array aggr
        ##returns N*K matrix
        N=dim(aggr)[2]
        K=dim(aggr)[3]
        matrix.ropar<-matrix(nrow=N,ncol=K)
        for(n in 1:N){
                for(k in 1:K){
                        matrix.ropar[n,k]<-ro.par(aggr[1:3,n,k],line)
                }
        }
        matrix.ropar
}
ro.par<-function(x,a){
        ##calculates the distance between point of perpendicular from M(x,y,z) on the line s(a,b,c) to the beginning of coords
        ##, where a is a vector of length 3 (a,b,c), and x 
        ##is a vector of length 3 x,y,z
        (x[1]*a[1]+x[2]*a[2]+x[3]*a[3])/sqrt(a[1]^2+a[2]^2+a[3]^2)
}
line.search<-function(aggr,S, line=rnorm(3)){
        ##A function to search for the axis of symmetry of the system using random
        ##search. 
        
        ## aggr - specifies the array of r objects with aggr[,,k] for each k
        ##being a vector of 3*N of coordinates of N particles.
        ##line is specified by 3 parametres a, b, c with a^2+b^2+c^2=1.
        
        ## n - specifies a particle, for which the line of symmetry is found.
        N=dim(aggr)[2]
        alfa=1
        sd1=vector(length=N)
        sd2=vector(length=N)
        for(i in 1:N){
                sd1[i]=sd(ro.perp.aggr(aggr[,i,],line))
        }
        for(i in 1:S){
                a2<-line-rnorm(3)/alfa
                for(i in 1:N){
                        sd2[i]=sd(ro.perp.aggr(aggr[,i,],a2))
                        
                }
                if(sum(sd2)<sum(sd1)){
                        line<-a2
                        sd1<-sd2
                        alfa<-alfa*0.9
                }
        }
        line<-line/sqrt(line[1]^2+line[2]^2+line[3]^2)
        print(c("returning vector for sd=",sd1))
        line
}
aggr.plot<-function(aggr,line=c(0,0,0),skip=1){
        ##skip - an integer representing, how often to plot a vector
        add<-FALSE
        library(rgl)
        K<-dim(aggr)[3]
        for(i in 1:K){
                        if(i%%skip==0){
                        myplot(aggr[,,i],add=add)
                        add<-TRUE
                }
        }
        x<-c(0,line[1])
        y<-c(0,line[2])
        z<-c(0,line[3])
        plot3d(x,y,z,type="l",add=TRUE)
}
        
##r<-gradient.descent(N=27,r=r,alfa=0.5,K=5000, print = TRUE)
##r<-reinit(N)
##ra<-rad(r)
##plot(sort(ra))
##myplot(r)
##U(r)

##hist(fun) plots distribution
