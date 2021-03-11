import numpy as np
import matplotlib.pyplot as plt
 
p=0.2 #payload kg
et=20*60 #endurance time sec
n=2
c=0.012  #in m   
rho=1.22    #kg/m3 
cl_alfa=5.73
v=2*3.7           #3 cell battery volt
#cp calculation
def cp(sigma):
  def lmda(r):
    theta=(22-17*r)*3.141/180
    return sigma*cl_alfa/16*(np.power(1+(32*theta*r/(sigma*cl_alfa)),0.5)-1)
 
  a=0
  b=1
  n=100
  h=(b-a)/n
  cd=0.01
 
  #Trapezoidal
  def cq(r):
    theta=(22-17*r)*3.141/180
    return 0.5*sigma*cl_alfa*(theta*r*r-lmda(r)*r)*lmda(r)
  m=0
  arr=[]
  for i in range (n+1):
    k=a+h*i
    cq(k)
    arr=np.append(arr,cq(k))
 
  for i in range(1,n):
    m=m+arr[i]
  y=(h/2)*(float(arr[0])+float(arr[n])+(2*float(m)))
  return (y*1.15+(0.5*sigma*cd*0.25))
 
 
 
#ct calculation
def ct(sigma):
  a=0
  b=1
  n=100
  h=(b-a)/n
  def lamda(r):
    theta=(22-17*r)*3.141/180
    return sigma*cl_alfa/16*((1+(32*theta*r/(sigma*cl_alfa)))**0.5-1)
 
  #Trapezoidal
  def ct(r):
    theta=(22-17*r)*3.141/180
    return 0.5*sigma*cl_alfa*(theta*r*r-lamda(r)*r)
  m=0
  arr=[]
  for i in range (n+1):
    k=a+h*i
    ct(k)
    arr=np.append(arr,ct(k))
 
  for i in range(1,n):
    m=m+arr[i]
  y=(h/2)*(float(arr[0])+float(arr[n])+(2*float(m)))
 
  return (y)
 
# Algorithm 
xx=[]
yy=[]
 
for dl in range(50,100):                     #disk loading N/m2
  
 
  gtow=p*3*1000                            # grams
  condition=True
  while (condition):
    
    T=gtow*0.001*9.81/4       # in Newton each rotor
    r=(T/(3.141*dl))**0.5      #r in m
    sigma=n*c/(3.141*r)
    cp(sigma)
    ct(sigma)
    omega=(T/(ct(sigma)*rho*3.141*r**4))**0.5
    N=omega*60/(2*3.141)                                #N rpm
    const_cp=rho*3.141*r*r*np.power((2*3.141*N*r/60),3)
    const_ct=rho*3.141*r*r*np.power((2*3.141*N*r/60),2)
    pow=cp(sigma)*const_cp*2.82                  #hover power
    thrust=ct(sigma)*const_ct                 #hover thrust
    cap=4*pow*et/(3600*0.001*v)                                #mAh
    i=pow/v                                #imax=2**1.5*pow/v
    kv=N/v
    l_bl=(4.891)*(i**0.1751)*(pow**0.2476)     #l_bl in mm
    d_bl=(41.45)*kv**(-0.1919)*(pow**0.1935)   #d_bl in mm
    m_esc=0.977*i**0.8483                        #m_esc in grams
    m_bat=0.0418*(cap**0.9327)*(3**1.0725)       #cell=3, mbat_in grams
    m_frame=1.3119*(r**1.2767)*(m_bat)**0.4587    # Airframe in grams
    m_bldc=0.0109*(kv**0.5122)*(pow**(-0.1902))*((np.log10(l_bl))**2.5582)*((np.log10(d_bl))**12.8502)    #m_bldc in grams
    m_rotor=0.0195*(r**2.0859)*(sigma**(-0.2038))*(n**0.5344)                                #m_rotor in grams
    total_weight=((4*m_esc)+m_bat+m_frame+4*(m_bldc)+(4*m_rotor)+(p*1000))           #m_total in grams
  
    condition=(total_weight-gtow)>1
    gtow=gtow+1                                                             #gtow increament 1  grams
 
 
  error=abs(total_weight-gtow)
  print(dl,error)
  
  xx.append(dl)
  yy.append(error)
 
  dl=dl+1
#print(yy)
plt.plot(xx,yy,'--o')
plt.xlabel('Disk loading N/m^2')
plt.ylabel('error')
print('minimum result ')
print(xx[np.argmin(yy)],min(yy))
 
print('mass of esc=',m_esc)
print('mass of bldc=',m_bldc)
print('mass of battery=',m_bat)
print('battery capacity (mAh)=',cap)
print('motor kv=',kv)
print('motor rpm', N)
print('hovering Power each motor=',pow)
print('hovering Thrust each motor=',thrust)
print('mass of total weight (grams)=',total_weight)
