import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import glob 


# Snapshots for different k
srcdir = './snaps_alpha=0.5/'
files = glob.glob(srcdir+'k=*.out')
data = np.array([np.loadtxt(f) for f in files])
k_list = [float(f.split('/')[-1][2:-4]) for f in files]
k_list = k_list[::2]
snaps = data.shape[2]-2

xx = data[0, :, 0]
colors = pl.cm.jet(np.linspace(0,1,len(k_list)))

for i_snap in range(1, snaps+1): 
	plt.figure()
	for i, k in enumerate(k_list):
		plt.plot(xx, data[i, :, i_snap+1], color=colors[i], label=r'$k=%s$' % k)
	plt.xlim(xx.min(), xx.max())
	plt.xlabel(r'$X$')
	plt.ylabel(r'$C_A$')
	plt.title(str(i_snap/snaps)+r'$T_{max}$')
	plt.legend()
	plt.savefig(srcdir+'snap'+str(i_snap)+'.png', dpi=300)


# Snapshots for different alpha 
srcdir = './snaps_k=2/'
files = glob.glob(srcdir+'alpha=*.out')
data = np.array([np.loadtxt(f) for f in files])
alpha_list = [float(f.split('/')[-1][6:-4]) for f in files]
snaps = data.shape[2]-2

xx = data[0, :, 0]
colors = pl.cm.jet(np.linspace(0,1,len(alpha_list)))

for i_snap in range(1, snaps+1): 
	plt.figure()
	for i, alpha in enumerate(alpha_list):
		plt.plot(xx, data[i, :, i_snap+1], color=colors[i], label=r'$\alpha=%s$' % alpha)
	plt.xlim(xx.min(), xx.max())
	plt.xlabel(r'$X$')
	plt.ylabel(r'$C_A$')
	plt.title(r'$%s T_{max}$' % (i_snap/snaps))
	plt.legend()
	plt.savefig(srcdir+'snap'+str(i_snap)+'.png', dpi=300)
# plt.show()
