import matplotlib.pyplot as plt
import numpy as np

"""
    parameters of the simulation
"""
r = 8.314
Ma=28.8e-3
ra=r/Ma
cp=1005.0
grav=9.81
jp=201
kp=600
dz=24.0
dy=50.0
dt=1.0
z=np.mgrid[0.0:dz*kp:dz]
zn=z+dz*0.5
y=np.mgrid[0.0:dy*jp:dy]
y=y-np.mean(y)
yn=y+dy*0.5
# just say first km, gamma=-ra/cp, then 6 degC/km
gamma1=9.81/cp
gamma2=6.3/1000.
cbase=1000.
psurf=1.e5
tsurf=290.
tn=np.zeros_like(zn)
pn=np.zeros_like(zn)
t=np.zeros_like(zn)
p=np.zeros_like(zn)
p[0]=psurf
t[0]=tsurf
tn[0]=t[0]-gamma1*dz*0.5


# reduce temperature based on lapse rate to cbase
k=1
while z[k] < cbase:
    t[k]=t[k-1]-gamma1*dz
    tn[k]=tn[k-1]-gamma1*dz
    k = k + 1
    

# take account of discontinuity
dz1=np.minimum(cbase-z[k-1],dz)
dz2=np.minimum(z[k]-cbase,dz)
dz1n=np.minimum(cbase-zn[k-1],dz)
dz2n=np.minimum(zn[k]-cbase,dz)
t[k]=t[k-1]-gamma1*dz1-gamma2*dz2
tn[k]=tn[k-1]-gamma1*dz1n-gamma2*dz2n
k=k+1

# do the rest based on new lapse rate
while k < kp:
    t[k]=t[k-1]-gamma2*dz1
    tn[k]=tn[k-1]-gamma2*dz1
    k = k + 1


# forward euler step for first one
pn[0]=p[0]-p[0]*grav*0.5*dz/(ra*t[0])
# now use centred -step
k=1
while k < kp:
    p[k]=p[k-1]-pn[k-1]*grav*dz/(ra*tn[k-1])
    pn[k]=pn[k-1]-p[k]*grav*dz/(ra*t[k])
    k = k + 1


rhon=pn/ra/tn
rho=p/ra/t


"""
    start here
"""
time=0.0

# page 60 of Vaughan's thesis - calculate fth,n, the chracteristic parameter for 
# each thermal
def fth_n(thermal_top,thermal_height,zbraking):
    if (thermal_top-zbraking) >= thermal_height:
        f=1.0
    elif ((thermal_top-zbraking) < thermal_height) & (thermal_height<=thermal_top):
        f=0.9*(thermal_top-thermal_height) / zbraking + 0.1
    else:
        f=0.
    return f

# page 60 of vaughan's thesis, calculate the reference profile
def reference_prof_w(thermal_top,thermal_height,zbraking,wmin,wmax,therm_base):
    zmax=thermal_top-zbraking
    w=0.0
    if((thermal_height >= therm_base) and (thermal_height<=zmax)):
        w=wmin+(thermal_height-therm_base)/(zmax-therm_base)*(wmax-wmin)
    elif((thermal_height > zmax) and (thermal_height<=thermal_top)):
        w=wmax*(1.0-(thermal_height-zmax)/(thermal_top-zmax))
        
    return w
    
    
# page 60 of vaughan's thesis, rate of vertical prop
def rate_of_vertical_prop(fth,thermal_top,thermal_height,zbraking,wmin,wmax,therm_base):
    wr=reference_prof_w(thermal_top,thermal_height,zbraking,wmin,wmax,therm_base)
    return fth*wr*0.5


    
# equation 4.2, vaughan's thesis
def peak_central_w(wbasic, fcn,fth,\
    thermal_top,thermal_height,zbraking,wmin,wmax,therm_base):
    
    Wth=rate_of_vertical_prop(fth,\
        thermal_top,thermal_height,zbraking,wmin,wmax,therm_base)

    Wthdash=Wth-wbasic
    
    wth=wbasic+Wthdash + 1.5*fcn*np.maximum(Wthdash,0.)
    
    return wth

# equation 4.4, vaughan's thesis
def lateral_w(wbasic, fcn,fth,\
    thermal_top,thermal_height,zbraking,wmin,wmax,therm_base):
    
    Wth=rate_of_vertical_prop(fth,\
        thermal_top,thermal_height,zbraking,wmin,wmax,therm_base)

    Wthdash=Wth-wbasic
    
    wth=wbasic+Wthdash + 1.5*fcn*np.maximum(Wthdash,0.)
    
    wl = wth*(1.0-1.5*(1.0-Wth/wth))
    return wl

if __name__=="__main__":
    print('Simulation started')
    # distance below clout-top below which thermal starts to brake
    z_braking=[500.,500.,500.,500.]
    # height each thermal reaches
    zmax=[2000.,3000.,5000.,6000.]
    fcn=[1.,1.,1.,1.]
    n=0

    times=np.mgrid[0:3600.:dt]
    thermal_height_store=np.zeros_like(times)
    prop_w_store=np.zeros_like(times)
    peak_w_store=np.zeros_like(times)
    lateral_w_store=np.zeros_like(times)
    wmin=2.0
    wmax=5.0
    wbasic=2.
    therm_base=500.0
    thermal_top=zmax[n]
    thermal_height=therm_base
    zbraking=z_braking[n]
    for nn in range(len(times)):
        # factor for each thermal
        fth=fth_n(thermal_top,thermal_height,zbraking)
        # calculate the speed of the nth thermal
        Wth=rate_of_vertical_prop(fth,thermal_top,\
            thermal_height,zbraking,wmin,wmax,therm_base)
            
        # calculate the lateral speed on edge of thermal
        Wlat=lateral_w(wbasic,fcn[n],fth,thermal_top,\
            thermal_height,zbraking,wmin,wmax,therm_base)
            
        # peak central
        Wcen=peak_central_w(wbasic,fcn[n],fth,thermal_top,\
            thermal_height,zbraking,wmin,wmax,therm_base)
            
            
        thermal_height=thermal_height+Wth*dt
        
        thermal_height_store[nn]=thermal_height
        prop_w_store[nn]=Wth
        peak_w_store[nn]=Wcen
        lateral_w_store[nn]=Wlat

    print('Simulation finished')
