#script to calculate macr values for many temperatures
T=c(0.01,0.02)        #list of pre-temperatures
N=27
for(t in T){
        r=reinit(27)
        r=temp(r,t)
        
}