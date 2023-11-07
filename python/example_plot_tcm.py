import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# fig, ax = plt.subplots()
# ax.set_xlim([0, 10])
# 
# scat = ax.scatter(1, 0)
# x = np.linspace(0, 10)

from netCDF4 import Dataset 

# open file
nc = Dataset('/tmp/output3.nc')

# extract needed variables
iivar=23
x = nc['x'][:]
z = nc['z'][:]
q =  nc['q'][:,:,:,iivar]
qmin=np.min(q)
qmax=np.max(q)
q = nc['q'][0,:,:,iivar]
# q = nc['w'][0,:,:]
time=nc['time'][:]


# plot out
# fig=plt.figure(figsize=[10,2.5])
fig, ax = plt.subplots()
ax=plt.pcolormesh(x,z,q[:-1,:-1].transpose(),vmin=qmin,vmax=qmax,shading='flat')
plt.title('time: ' + str(time[0]) + 's')
plt.xlabel('x (m)')
plt.ylabel('z (m)')

#plt.clim((-0.05,0.05))
plt.colorbar()

plt.gca().set_aspect('equal')
plt.gca().autoscale(tight=True)

def update(i):
    q = nc['q'][i*10,:,:,iivar]
    q = q[:-1, :-1]
#     q = np.copy(nc['w'][i*10,:,:])
    ax.set_array(q.transpose().flatten())
    ax.axes.set_title('time: ' + str(np.copy(time[i*10])) + 's')

ani = animation.FuncAnimation(fig, update, frames=int(len(time)/10),blit=False,repeat=False)


ani.save('/tmp/plot_tcm.gif')
nc.close()

plt.close('all')
