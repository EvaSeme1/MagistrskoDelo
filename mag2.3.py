# -*- coding: utf-8 -*-
#The program draws the graph of FWHM of errors in dependence of artificial signal's sigma. It compares linear interpolation, improved linear interpoltion, third order interpolation and improved third order interpolation. 

import math
#q=[sigma, amplitude]
#The function calculates the value of a function (Gaussian distribution or inverse Gaussian distribution) at the given point (x or C).
def Gauss(q, x):
    return math.exp(-(x/q[0])**2/2)*q[1]

def invGauss(q, C):
    return (math.sqrt(-2*math.log(C/q[1]))*q[0]) 

#q=[sigma, amplitude]
#a is time between two samples - in our case 1
#s is the starting point
#N is the number of signals
#The function returns the matrix of sets of measurements in form of mat.
#Z is vector with times between starting points and the first points  
def matrix(q, a, s, N):
    import random
    mat=[0]*N
    i=0
    Z=[0]*N
    while i<N:
        w=[]
        z=a*random.random()
        Z[i]=z
        x=s+z
        y=Gauss(q, x)
        w.append(y)
        while x<(-s):
            x=x+a
            y=Gauss(q, x)
            w.append(y)
        mat[i]=w
        i=i+1
    return (mat, Z)

#The function gets the measurements in form of matrix. 
#C is our constant
#The function return four vectors of sequential points. The value C is between the second and the third point.
#K is vector of locations of points before the signal reaches C
def choose3(matrix, C):
    N=len(matrix)
    w=[0]*N
    v=[0]*N
    K=[0]*N
    l=[0]*N
    u=[0]*N
    i=0
    while i<N:
       row=matrix[i]
       n=len(row)
       j=0
       while j<n:
           if row[j]>C:
              w[i]=row[j]
              v[i]=row[j-1]
              l[i]=row[j-2]
              u[i]=row[j+1]
              K[i]=j-1
              break
           j=j+1
       i=i+1
    return (l, v, w, u, K)

#C is our constant
#v is vector of points before C
#w is vector of points after C
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

#The function returns a vector with elements arranged from the smallest to the biggest.
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

#The function is the composition of function tstar1 and function order.
def ALL1(v, w, C):
    s=tstar1(v, w, C)
    t=order(s)
    return t

#s is starting point
#a is time between two samples 
#q=[sigma, amplitude]
#C is our constant
#Z is matrix with times between starting points and the first points 
#K is vector of locations of points before the signal reaches C
#N is number of signals
#The function returns the correct t-s in vector T.
def correct(s, a, q, C, Z, K, N):
    j=0
    S=[0]*N
    while j<N:
        y=invGauss(q, C)
        a=1
        t=(-s-Z[j]-a*K[j]-y)/a
        S[j]=t
        j=j+1
    T=order(S)   
    return(T)

#The function gets the vector of wrong values and the correct ones and returns the vector of absolute errors.
#only for unimproved interpolations
def absolute(t, c):
    n=len(t)
    A=[0]*n
    i=0
    while i<n:
        A[i]=c[i]-t[i]
        i=i+1
    return(A)
    
#The function gets vector of absolute errors.
#It returns FWHM.
def error(A):
    import numpy
    import math
    FWHM=2*math.sqrt(2*math.log(2))*numpy.std(A)
    return(FWHM)

#C is our constant
#l, v, wand u are the vectors of sequential points chosen so that the signal reaches C between the second and the third point.
#The function returns times goten with third order interpolation.    
def tstar3(l, v, w, u, C):
    y=[-1, 0, 1, 2]
    i=0
    N=len(v)
    t=[0]*N
    import scipy.interpolate
    import numpy
    while i<N:
        W=scipy.interpolate.lagrange(y, [l[i], v[i], w[i], u[i]])
        X=numpy.roots(W-C)
        j=0
        while j<3:
            if abs(X[j].imag)<1e-5:
                R=X[j].real
                if R<1:
                    if R>0:
                        t[i]=R     
            j=j+1
        i=i+1
    return(t)    

##The function is the composition of function tstar3 and function order.
def ALL3(l, v, w, u, C):
    s=tstar3(l, v, w, u, C)
    t=order(s)
    return t
    
#t is vector of t*-s
#a is number of columns
#g is vector of heigths of columns in the histogram
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

#N is number of measurements
#a is number of columns in histogram
#g is vector of heights of columns in histogram    
#The function returns vector of points on the cumulative function.
def cumulative(N, a, g):
    s=0
    k=[0]*a
    i=0
    while i<a:
        s=s+g[i]/N
        k[i]=s
        i=i+1
    return(k)    

#k is vector of points on the cumulative function
#a is number of columns
#The function gets the value obtained with interpolation in form of x. It returns the suitable correct time.
def function(k, a, x):
    i=1
    p=1/a
    while i<a+1:
        if x<i*p:
          return k[i-1]#+(x-(i-1)*p)*a*(k[i]-k[i-1]) #(linear interpolation of cumulative function for better results)  
        i=i+1

#k is vector of points on the cumulative function
#a is number of columns
#The function gets the vector of wrong values (vector t) and the correct ones (vector c) and returns the vector of absolute errors.
#only for improved interpolations
def myabsolute(k, a, t, c):
    n=len(t)
    A=[0]*n
    i=0
    p=[0]*n
    while i<n:
        p[i]=function(k, a, t[i])
        A[i]=c[i]-p[i]
        i=i+1
    return(A)

#Main part of the code:
import numpy
#Setting variables:
q=[5, 1]
C=0.2
s=-30
N=1000000
a=256

#Calculating FWHM (improved linear interpolation), FWHM1 (linear interpolation) and FWHM3 (cubic interpolation) for a particular q[0]:
mat, Z=matrix(q, 1, s, N)
l, v, w, u, K=choose3(mat, C)
t1=ALL1(v, w, C) 
t3=ALL3(l, v, w, u, C)
c=correct(s, 1, q, C, Z, K, N)
A1=absolute(t1, c)
A3=absolute(t3, c)
FWHM1=(error(A1))
FWHM3=(error(A3))
g1=histogram(t1, a)
k1=cumulative(N, a, g1)
NM=10000
mat, Z=matrix(q, 1, s, NM)
l, v, w, u, K=choose3(mat, C)
t1=ALL1(v, w, C) 
c=correct(s, 1, q, C, Z, K, NM)
A=myabsolute(k1, a, t1, c)
FWHM=error(A)