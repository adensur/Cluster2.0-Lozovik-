init<-function(N){##initializes matrix 3 x N with correct rownames and random values
        r<-cbind(rnorm(N),rnorm(N),rnorm(N))
        colnames(r)<-(c("x","y","z"))
        r
}