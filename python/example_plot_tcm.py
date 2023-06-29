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
x = nc['x'][:]
z = nc['z'][:]
q = nc['q'][0,:,:,15]
time=nc['time'][:]


# plot out
# fig=plt.figure(figsize=[10,2.5])
fig, ax = plt.subplots()
ax=plt.pcolormesh(x,z,q.transpose(),shading='gouraud')
plt.title('time: ' + str(time[0]) + 's')
plt.xlabel('x (m)')
plt.ylabel('z (m)')

#plt.clim((-0.05,0.05))
plt.colorbar()

plt.gca().set_aspect('equal')
plt.gca().autoscale(tight=True)

def update(i):
    q = np.copy(nc['q'][i*10,:,:,15])
    ax.set_array(q.transpose().ravel())
    ax.axes.set_title('time: ' + str(np.copy(time[i*10])) + 's')

ani = animation.FuncAnimation(fig, update, frames=int(len(time)/10))


ani.save('/tmp/plot_tcm.gif')
nc.close()

plt.close('all')
