description of functions and variables, used in cluster2.0.R

function U
1 U<-function(r){##returns a value of the potential energy of the system
2        N<-ncol(r)
3        U<-0
4        for(i in 1:N){
5               sum<-0
6                if(i<N){
7                        for(k in (i+1):N){
8                                sum<-sum+rki(r,k,i)^(-5)
9                        }
10                }
11                U<-U+sum+r[1,i]^2+r[2,i]^2+r[3,i]^2
12        }
13        names(U)<-"Potential energy"
14        U
15 }
input: vector "r" 3 to N of x,y,z coordinates of particles.
U is assumed to be equal to:

sum(i=1:N)(r[i]^2)+sum(i,k=1:N, k>i)(1/r[i,k]^5)	(1)

The first term in (1) corresponds to parabolic confinement potential; minimization of it cause particles
to be attracted to the center of the confinement. The second term corresponds to the repulsion of high order
of particles between each other; the minimisation of it should cause them all to just scatter into infinity;
in the presence of the confinement the particles will (supposedly) form a stable structure.

In order to calculate sums and double sums in (1), "for" loops are used.

The first term has sum using only 1 index; it can be calculated easily with just 1 "for" loop.
In my function U() there is the outer "for" loop, cycling through all indexes from 1 to N (where N is the number
of particles), and each iteration adds r[i]^2 to the whole potential energy. (see lines 4, 11,12)

The second term uses sum by 2 indexes 1 and k. Basically, it has all the possible pairs of (i,k) with i!=k,
where both i and k take values of 1 to N. It is possible to just cycle through all i's and all k's, excluding the
iterations where i equals to k. However, that way each pair will be counted twice, for example: (2,3) and (3,2) 
should be a single term in potential energy. In order to avoid that, I calculate sum for index i in 1:N,
and then for index k from i+1 to N. THat way, k never equals to i, and each pair is counted one and only one time.

This algorhytm uses 2 "for" cycles one inside another. (lines 4 and 7). The "if" condition on line 6 is put there,
because otherwise when i equals N, the vector (i+1):N has 2 elements: N+1 and N. We dont want anything to be calculated 
when i equals N for the 2nd term, because that "pair" was already counted when i was N-1, and k was N.

FUNCTION grad.U.
grad.U is a function that returns the gradient of energy U, or, to be exact, the derivative of U by the coordinate var
(var=1 means x, var = 2 means y and so on) of particle k. So grad.U(var=1, k=3) means derivative of U by x[3].
In order to calculate derivative, we have to derive every term in two sums in (1).
The terms of the first sum are all independent, and only one derivative is not equal to zero.
For the second sum we will have total of (N-1) non-zero derivatives, for each pair of indexes that contain certain k.
For example, for k=3 and i=5 there would be 4 terms (1,3), (2,3), (4,3) and (5,3), derivative of which is not equal to zero.
In order to calculate it all, I use one "for" loop with N iterations for calculating derivatives of sum#2, and "if"
condition to exclute iteration with i=k.










