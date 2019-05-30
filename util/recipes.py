#Generic recipes for making key calculations

import numpy as np

def nbinmaker(x,nbin=10,nsum=False):
    #Split an array into nbin bin's of equal number elements
   
    x=np.asarray(x)
 
    xorder=sorted(range(0,len(x)),key=lambda k:x[k])

    x_lower=np.array([])
    x_upper=np.array([])
    x_hist=np.zeros(nbin)
    x_sum=np.zeros(nbin)
    x_mid=np.array([])

    for i in range(0,nbin):
        indx=int(float(i)*float(len(x))/float(nbin))
        x_lower=np.append(x_lower,x[xorder[indx]])
        indx=int(float(i+1)*float(len(x))/float(nbin))-1
        x_upper=np.append(x_upper,x[xorder[indx]])

    for j in range(0,nbin):
        indx=(x>=x_lower[j]) * (x<=x_upper[j])
        x_hist[j]=len(x[indx])
        x_sum[j]=np.sum(x[indx])
        x_mid=np.append(x_mid,x_sum[j]/x_hist[j])

    if nsum:
        return x_lower,x_mid,x_upper,x_hist,x_sum
    else:
        return x_lower,x_mid,x_upper,x_hist

def binmaker(x,nbin=10,nsum=False,steptype='linear'):
    #Split an array into nbin bin's of equal size
    
    x_hist=np.zeros(nbin)
    x_sum=np.zeros(nbin)
    x=np.array(x)
 
    if steptype=='linear':
        steps=np.linspace(np.amin(x),np.amax(x),nbin+1)
    else:
        steps=np.logspace(np.log10(np.amin(x)),np.log10(np.amax(x)),nbin+1)

    x_lower=steps[:-1]
    x_upper=steps[1:]
    x_mid=(x_upper+x_lower)/2.

    for j in range(0,nbin):
        indx=(x>=x_lower[j]) * (x<=x_upper[j])
        x_hist[j]=len(x[indx])
        x_sum[j]=np.sum(x[indx])

    if nsum:
        return x_lower,x_mid,x_upper,x_hist,x_sum
    else:
        return x_lower,x_mid,x_upper,x_hist

def dx_function(x,nx=10):
    #Find distribution function using nx bins containing an equal number of points
    
    x_lower,x_mean,x_upper,x_hist=nbinmaker(x,nx)
   
    lx_mean=np.log10(x_mean)
    dx=x_hist/(x_upper-x_lower)
    ldx=np.log10(dx)

    (alpha,yalpha),V=np.polyfit(lx_mean,ldx,1,cov=True)
    ealpha=np.sqrt(V[0][0])
    eyalpha=np.sqrt(V[1][1])

    return x_mean,x_hist,dx,alpha,ealpha,yalpha,eyalpha

def mean_prof(x,y,nbin=10,bintype='fix',steptype='linear',nsum=False,median=False):
    #Calculate mean profile of parameter y that depends on x  
    if bintype=='num':
        x_lower,x_mid,x_upper,x_hist=nbinmaker(x,nbin,nsum)
    else:
        x_lower,x_mid,x_upper,x_hist=binmaker(x,nbin,nsum,steptype=steptype)

    y_bin=[]
    y_sig=[]
    x_bin=[]

    for i in range(0,nbin):
        indx=(x>=x_lower[i]) * (x<=x_upper[i])
        
        if True in indx:
            x_bin=np.append(x_bin,x_mid[i])
            if x_hist[i]>1:
                if median:
                    y_bin=np.append(y_bin,np.median(y[indx]))
                else:
                    y_bin=np.append(y_bin,np.mean(y[indx]))

                y_sig=np.append(y_sig,np.std(y[indx]))
            elif x_hist[i]==1:
                y_bin=np.append(y_bin,y[indx])
                y_sig=np.append(y_sig,0.0)
            else:
                y_bin=np.append(y_bin,0.0)
                y_sig=np.append(y_sig,0.0)

    return np.array(x_bin),np.array(y_bin),np.array(y_sig)

def smooth(x,y,nbin=10,bintype='fix',median=False,**kwargs):
    #Smooth a profile

    x=np.array(x)
    y=np.array(y)
    
    #Smooth by intervals in dx
    if 'dx' in kwargs:
        dx=float(kwargs.get('dx'))
        nbin=int((np.max(x)-np.min(x))/dx)
    
    if bintype=='num':
        x_lower,x_mid,x_upper,x_hist=nbinmaker(x,nbin)
    else:
        x_lower,x_mid,x_upper,x_hist=binmaker(x,nbin)

    y_bin=[]
    y_sig=[]
    x_bin=[]
    y_max=[]
    y_min=[]

    for i in range(0,nbin):
        indx=(x>=x_lower[i]) * (x<=x_upper[i])
        
        if True in indx:
            x_bin=np.append(x_bin,x_mid[i])
            if x_hist[i]>1:
                if median:
                    y_bin=np.append(y_bin,np.median(y[indx]))
                else:
                    y_bin=np.append(y_bin,np.mean(y[indx]))
                y_sig=np.append(y_sig,np.std(y[indx]))
                y_min=np.append(y_min,np.min(y[indx]))
                y_max=np.append(y_max,np.max(y[indx]))
            elif x_hist[i]==1:
                y_bin=np.append(y_bin,y[indx])
                y_sig=np.append(y_sig,0.0)
                y_min=np.append(y_min,y[indx])
                y_max=np.append(y_max,y[indx])
            else:
                y_bin=np.append(y_bin,0.0)
                y_sig=np.append(y_sig,0.0)
                y_min=np.append(y_min,0.0)
                y_max=np.append(y_max,0.0)

    return x_bin,y_bin,y_sig,y_min,y_max

def interpolate(r1,r2,x=None,y=None):

    x1,y1=r1
    x2,y2=r2

    m=(y2-y1)/(x2-x1)
    b=y1-m*x1

    if x!=None:
        return m*x+b
    elif y!=None:
        return (y-b)/m
    else:
        print('NO INTERMEDIATE COORDINATE GIVEN')
        return 0
        
def rotate(x,y,z,thetax=0.,thetay=0.,thetaz=0.):

    #rotate about x axis:
    y2=y*np.cos(thetax)-z*np.sin(thetax)
    z2=-y*np.sin(thetax)+z*np.cos(thetax)

    y=y2
    z=z2

    #rotate about y axis:
    x2=x*np.cos(thetay)-z*np.sin(thetay)
    z2=-x*np.sin(thetay)+z*np.cos(thetay)

    x=x2
    z=z2

    #rotate about z axis:
    x2=x*np.cos(thetaz)-y*np.sin(thetaz)
    y2=-x*np.sin(thetaz)+y*np.cos(thetaz)

    x=x2
    y=y2

    return np.array(x),np.array(y),np.array(z)






