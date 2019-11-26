import matplotlib
matplotlib.use('Agg')

import nbodypy as npy
import numpy as np
import matplotlib.pyplot as plt
import os
for i in range(0,1281):

    cluster=npy.load_snapshot(filename='%s.dat' % str(i).zfill(5),units='realkpc',origin='galaxy')

	npy.starplot(cluster,coords='xy',filename='%s.png' % str(i).zfill(5))
	plt.close()

try:
	os.system('ffmpeg -start_number 0  -i %05d.png movie.mp4')
except:
	print('FFMPEG NOT FOUND TO CONVERT PNG FILES INTO MP4')
