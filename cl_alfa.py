import numpy as np
import matplotlib.pyplot as plt
#given
rho=1.22
N=1500
r=0.355
sigma=0.0573 #(n*b/pi*r)
cl_alfa=5.73
cd=0.01
const=rho*3.141*r*r*np.power((2*3.141*N*r/60),2)

def result_ct(theta):
  def lmda(r):
    return sigma*cl_alfa/16*(np.power(1+(32*theta*r/(sigma*cl_alfa)),0.5)-1)
 
  a=0
  b=1
  n=100
  h=(b-a)/n
  

  #Trapezoidal
  def ct(r):
    return 0.5*sigma*cl_alfa*(theta*r*r-lmda(r)*r)
  m=0
  arr=[]
  for i in range (n+1):
    k=a+h*i
    ct(k)
    arr=np.append(arr,ct(k))

  for i in range(1,n):
    m=m+arr[i]
  y=(h/2)*(float(arr[0])+float(arr[n])+(2*float(m)))
  return y*const


 
yy=[]
xx=[]
 
