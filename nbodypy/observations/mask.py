from ..main.cluster import sub_cluster

import numpy as np
from scipy.spatial import ConvexHull

#Place a mask over the simulated data in order to replicate an observational dataset

class ObservationalMask(object):

	def __init__(self,name):

		self.name=name
		self.cluster=None

		if name=='m15' or name=='M15':
			#From Emmanulle Dalessandro - TODO reference GO Proposals
			self.ra_gc=322.4929716
			self.dec_gc=12.1670539
			self.dist_gc=0.
			self.pmra_gc=0.
			self.pmdec_gc=0.
			self.vlos_gc=0.

			self.rt=900.0/3600.0
			self.rm=78.0/3600.0

		elif name=='m30' or name=='M30':
			self.ra_gc=325.092194
			self.dec_gc=-23.1798325
			self.dist_gc=0.
			self.pmra_gc=0.
			self.pmdec_gc=0.
			self.vlos_gc=0.

			self.rt=1540.2/3600.0
			self.rm=61.8/3600.0		

		self.mbintype='num'
		self.rbintype='fix'


	def set_mcorr(self,mcorr):
		#Se mass correction for entire dataset
		self.mcorr=mcorr

	def set_cluster(self,cluster):
		#Set base cluster to use values of rm for
		self.cluster=cluster

	def apply_orbit(self):
		#Apply orbit to self.cluster based on mask
		self.cluster.add_orbit(self.ra_gc,self.dec_gc,self.dist_gc,self.pmra_gc,self.pmdec_gc,self.vlos_gc)

	def apply_key_params(self):
		#Force key parameters to be specific values - useful if only part of the data is provided
		self.cluster.rm=self.rm
		self.cluster.rt=self.rt

	def apply_mask(self):
		#Mask simulated data to reflect fields of view of observational

		self.cluster.to_centre()
		if self.cluster.rm is None:
			self.cluster.key_params(do_order=True)

		if self.name=='m15' or self.name=='M15':

			#Set limits in terms of x/rm and y/rm
			q1lims=np.array([-10.442896000073247,-5.881643789792685,0.6284492307691526,6.75208615384611])
			q2lims=np.array([-10.236322576156747,-5.877316149603139,-7.435661538461547,-1.5543600000000606])
			q3lims=np.array([-4.870638209602899,-0.2689031587450201,0.5989707692307068,6.760892307692282])
			q4lims=np.array([-4.889138688343503,-0.493742739092162,-7.426458461538488,-1.558675384615454])

			q1indx=((self.cluster.x/self.cluster.rm) >= q1lims[0])*((self.cluster.x/self.cluster.rm) <= q1lims[1]) * ((self.cluster.y/self.cluster.rm) >= q1lims[2]) * ((self.cluster.y/self.cluster.rm) <= q1lims[3])
			q2indx=((self.cluster.x/self.cluster.rm) >= q2lims[0])*((self.cluster.x/self.cluster.rm) <= q2lims[1]) * ((self.cluster.y/self.cluster.rm) >= q2lims[2]) * ((self.cluster.y/self.cluster.rm) <= q2lims[3])
			q3indx=((self.cluster.x/self.cluster.rm) >= q3lims[0])*((self.cluster.x/self.cluster.rm) <= q3lims[1]) * ((self.cluster.y/self.cluster.rm) >= q3lims[2]) * ((self.cluster.y/self.cluster.rm) <= q3lims[3])
			q4indx=((self.cluster.x/self.cluster.rm) >= q4lims[0])*((self.cluster.x/self.cluster.rm) <= q4lims[1]) * ((self.cluster.y/self.cluster.rm) >= q4lims[2]) * ((self.cluster.y/self.cluster.rm) <= q4lims[3])
			
			xhull=np.array([-0.025200275014222413 , -0.023248840101820767 , -0.022252171603265446 , -0.020876374968530564 , -0.019768232948539806 , -0.019516293671548123 , -0.018310827507562306 , -0.016141939740057017 , -0.01507099008868237 , 0.020622932660628257 , 0.030769270910864063 , 0.03402070578819945 , 0.034327150432339686 , 0.03464686440992152 , 0.034355000816019575 , 0.03374778542629317 , 0.02442140775679065 , 0.024403514928413285 , 0.023904902487531996 , 0.022511976486394615 , 0.020470072739238877 , 0.0188020846312626 , 0.017870995541351138 , 0.017591948439801587 , 0.008300117344556341 , 0.007085393874204251 , -0.009972763393508404 , -0.011261330698878028 , -0.02214748809652206 , -0.026122585612488373 , -0.03346958022086817 , -0.03399466578178241 , -0.03309639437024884 , -0.03100587078182435 , -0.02904970980645645 , -0.028262363464591378 , -0.027903603829392717 , -0.025745234255737832])/self.rm
			yhull=np.array([-0.03149380369138381 , -0.03327898128037248 , -0.03394116648375825 , -0.03479537796031218 , -0.03535176259000022 , -0.035427181195303775 , -0.03539306698617067 , -0.03488420765700969 , -0.034618570591543 , -0.025546899011511297 , -0.022961518217483424 , -0.021911221962715968 , -0.021612782576963474 , -0.02111804111942903 , -0.01924967910897742 , -0.01598325703837688 , 0.03095362077113272 , 0.03101861911798951 , 0.03242387358258812 , 0.03366825171526268 , 0.034802686363492465 , 0.03536156299431959 , 0.03552819870817594 , 0.035491980093389976 , 0.03296722781833852 , 0.03263449270445504 , 0.027920886042690843 , 0.02755543757181496 , 0.02444082224849501 , 0.023292683397820806 , 0.021116007388034132 , 0.02039697408638987 , 0.013997660928067528 , 0.0012260088299883426 , -0.010216512222578881 , -0.014745097038631718 , -0.016765934883835022 , -0.02861945184298367])/self.rm

			hindx=inhull(xhull,yhull,self.cluster.x/self.cluster.rm,self.cluster.y/self.cluster.rm)

			self.oindx=np.logical_or(np.logical_or(np.logical_or(q1indx,q2indx),np.logical_or(q3indx,q4indx)),hindx)
		elif self.name=='m30' or self.name=='M30':
			self.oindx=np.ones(self.cluster.ntot,bool)


		ocluster=sub_cluster(self.cluster,indx=self.oindx)

		if self.cluster is not None:
			ocluster.rm=self.cluster.rm
		else:
			ocluster.rm=self.rm

		return ocluster

	def load_rbins(self):
		#Load predefined radial bins. Scale by rm if self.cluster is provided
		if self.name=='m15' or self.name=='M15':
			self.r_lower=np.array([25.,40.,55.,70.,85.,250.,350.,450.,550.,650.])/3600.
			self.r_upper=np.array([40.,55.,70.,85.,100.,350.,450.,550.,650.,750.])/3600.
		elif self.name=='m30' or self.name=='M30':
			self.r_lower=np.array([10.,20.,40.,200.,250.,350.,650.])/3600.
			self.r_upper=np.array([20.,40.,100.,250.,350.,650.,1000.])/3600.
			
		self.r_mean=(self.r_lower+self.r_upper)/2.


		#Scale by half-mass radius if simlated cluster is given
		if self.cluster is not None:
			self.r_lower*=self.cluster.rm/self.rm
			self.r_upper*=self.cluster.rm/self.rm
			self.r_mean*=self.cluster.rm/self.rm

	def load_mbins(self):
		#Load predefined mass bins
		if self.name=='m15' or self.name=='M15':
			self.m_lower=np.array([0.4,0.45,0.5,0.55,0.6,0.65,0.7])
			self.m_upper=np.array([0.45,0.5,0.55,0.6,0.65,0.7,0.75])
		elif self.name=='m30' or self.name=='M30':
			#self.m_lower=np.array([0.45,0.5,0.55,0.6,0.65,0.7])
			#self.m_upper=np.array([0.5,0.55,0.6,0.65,0.7,0.75])
			step=(0.779999971-0.4)/10.0
			self.m_lower=0.4+step*np.linspace(0,9,10)
			self.m_upper=0.4+step*np.linspace(1,10,10)

		self.m_mean=(self.m_lower+self.m_upper)/2.

def inhull(xhull,yhull,x,y):

	H=np.column_stack([xhull,yhull])
	hull=ConvexHull(H)

	indx=np.zeros(len(x),bool)
	for i in range(0,len(x)):
		inside = True
		for ind in range(1, len(hull.vertices)):
			res = cross(H[hull.vertices[ind-1]], H[hull.vertices[ind]], (x[i],y[i]))

			if res < 0:
				inside = False

		res = cross(H[hull.vertices[-1]], H[hull.vertices[0]], (x[i],y[i]))
		if res < 0:
			inside = False

		indx[i]=inside

	return indx

def cross(o, a, b):
    """ 2D cross product of OA and OB vectors,
     i.e. z-component of their 3D cross product.
    :param o: point O
    :param a: point A
    :param b: point B
    :return cross product of vectors OA and OB (OA x OB),
     positive if OAB makes a counter-clockwise turn,
     negative for clockwise turn, and zero
     if the points are colinear.
    """
 
    return (a[0] - o[0]) * (b[1] - o[1]) -\
           (a[1] - o[1]) * (b[0] - o[0])
