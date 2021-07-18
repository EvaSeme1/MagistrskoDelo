#The program reads the file with the values of signals and draws a graph of the ratio of the area under the exponential tail and the area under the signal in dependence of the area under the signal.

manj=718
N=718
import matplotlib.pyplot as plt
import numpy as np
import random

M=[[0]]*manj
stetje=0
file = open('signali1.txt', 'r') 
w=[0]*(2504)
n=0
while n<manj:
    z=[]
    tz=[]
    j=0
    D=40
    Dd=1
    x=int(random.random()*D)
    i=0
    while i<2504:
        w[i]=file.readline()
        tz.append((-1)*float(w[i]))
        i=i+1 
    x=int(random.random()*D)
    i=0
    while i<2504-2*D:
        if i%(Dd*D)==0:
            S=tz[i+x]
            k=1
            while k<D:
                S=S+tz[i+k+x]
                k=k+1
            S=S/D
            z.append(S)
        i=i+1
    p=np.amax(tz)
    if p>0:
        M[stetje]=z
        stetje=stetje+1
    n=n+1
file.close()
manj=stetje
N=stetje

#The function gets the measurements in the form of matrix. 
#C is our constant
#The function returns two vectors of sequential points. The value C is between the first and the second point.
#K is the vector of locations of points before the signal reaches C
def choose1(matrix, C):
    w=[0]*N
    v=[0]*N
    K=[0]*N
    i=0
    while i<N:
       row=matrix[i]
       n=len(row)
       j=0
       while j<n:
           if row[j]>C:
              w[i]=row[j]
              v[i]=row[j-1]
              K[i]=j-1
              break
           j=j+1
       i=i+1
    return (v, w, K)
C=5
vN,wN,KN=choose1(M,C)

#C is our constant
#v is the vector of points before C
#w is the vector of points after C
#The function returns the vector t with random t*, obtained with linear interpolation. 
def tstar1(v, w, C):
    y=[0, 1]
    i=0
    N=len(v) 
    t=[0]*N
    import scipy.interpolate
    while i<N:
        p=[v[i], w[i]]
        f=scipy.interpolate.interp1d(p, y)
        t[i]=float(f(C))
        i=i+1
    return(t)    

#The function returns the vector with elements arranged from the smallest to the biggest.
def order(s):
    import numpy
    i=0
    n=len(s)
    t=[0]*n
    while i<n:
        a=numpy.amin(s)
        t[i]=a
        s.remove(a)
        i=i+1
    return t

tN=tstar1(vN, wN, C)
tsN=[0]*len(tN)
i=0
while i<len(tN):
    tsN[i]=tN[i]
    i=i+1
tN=order(tN)  

#t is the vector of t*-s
#a is the number of columns
#g is the vector of heigths of columns in the histogram
def histogram(t, a):
    import math
    g=[0]*a
    n=len(t)
    i=0
    while i<n:
        if t[i]!=1:
           j=math.floor(t[i]*a)
           g[j]=g[j]+1
           i=i+1 
        else:
            g[a-1]=g[a-1]+1
            i=i+1
    return g
gN=histogram(tN,10)

#N is the number of measurements
#a is the number of columns in the histogram
#g is the vector of heights of columns in the histogram    
#The function returns the vector of points on the cumulative function.
def cumulative(N, a, g):
    s=0
    k=[0]*(a+1)
    i=1
    while i<a+1:
        s=s+g[i-1]
        k[i]=s
        i=i+1
    i=1
    while i<a+1:
        k[i]=k[i]/s
        i=i+1
    return(k)    
kN=cumulative(N,10,gN)

#k is the vector of points on the cumulative function
#a is the number of columns
#The function gets the value obtained with interpolation in the form of x. It returns a suitable correct time.
def function(k, a, x):
    #a=len(k)
    i=1
    p=1/a
    while i<a+1:
        if x<i*p:
          return k[i-1]+(x-(i-1)*p)*a*(k[i]-k[i-1])
        i=i+1
        
PtN=[0]*N
i=0
while i<N:
    PtN[i]=function(kN,10,tsN[i])
    i=i+1
xN=[]
yN=[]

#The function gets matrix of values of signals and gives ratios of areas and the area under the signal.
def ploscine(mat,K,Pt):
    areas=[0]*N
    areas2=[0]*N
    ratio=[0]*N
    i=0
    while i<N:
        row=mat[i]
        j=int(800/(D*Dd))
        while j<int(2000/(D*Dd)):
            areas[i]=areas[i]+row[j]
            j=j+1
        L=int(1050/(D*Dd))
        areas2[i]=areas2[i]+row[L]*0.5*(1-Pt[i])**2
        areas2[i]=areas2[i]+row[L+1]*((2-Pt[i])/2+(1-Pt[i])*Pt[i]/2)
        j=L+2
        while j<int(2000/(D*Dd)):
            areas2[i]=areas2[i]+row[j]
            j=j+1
        ratio[i]=areas2[i]/areas[i]
        i=i+1
    return(ratio, areas)
ratioN, areasN=ploscine(M,KN,PtN)
i=0
while i<len(ratioN):
    if ratioN[i]<0.4 and ratioN[i]>0.05:
        xN.append(PtN[i])
        yN.append(ratioN[i])
    i=i+1
    
import matplotlib.pyplot as plt

#Graph of the ratio of areas in dependence of the area under the signal:
#plt.plot(areasN,ratioN,'ro')
#plt.ylabel('razmerje ploščin')
#plt.xlabel('ploščina pod signalom')

#Distribution of signals by their ratio of areas (histogram):
#plt.hist(yN,20)

#Gaph of the ratio of areas in dependence of phase of the signal:
#plt.plot(xN,yN,'ro')  