import numpy as np
from collections import Counter

#Split an array into nbin bin's of equal number elements
def nbinmaker(x,nbin,nsum=False):
    
    xorder=sorted(range(0,len(x)),key=lambda k:x[k])

    x_lower=np.array([])
    x_upper=np.array([])
    x_hist=np.zeros(nbin)
    x_sum=np.zeros(nbin)

    for i in range(0,nbin):
        indx=int(float(i)*float(len(x))/float(nbin))
        x_lower=np.append(x_lower,x[xorder[indx]])
        indx=int(float(i+1)*float(len(x))/float(nbin))-1
        x_upper=np.append(x_upper,x[xorder[indx]])

    for j in range(0,nbin):
        indx=(x>=x_lower[j]) * (x<=x_upper[j])
        x_hist[j]=len(x[indx])
        x_sum[j]=np.sum(x[indx])


    x_mid=np.array([])
    for i in range(0,nbin):
        x_mid=np.append(x_mid,x_sum[i]/x_hist[i])

    if nsum:
        return x_lower,x_mid,x_upper,x_hist,x_sum
    else:
        return x_lower,x_mid,x_upper,x_hist

#Split an array into nbin bin's of equal size
def binmaker(x,nbin,nsum=False):
    

    x_lower=np.array([])
    x_upper=np.array([])
    x_mid=np.array([])
    x_hist=np.zeros(nbin)
    x_sum=np.zeros(nbin)

    step=abs(np.amax(x)-np.amin(x))/float(nbin)

    for i in range(0,nbin):
        x_lower=np.append(x_lower,np.amin(x)+step*float(i))
        x_upper=np.append(x_upper,np.amin(x)+step*float(i+1))
        x_mid=np.append(x_mid,(x_lower[-1]+x_upper[-1])/2.)

    for j in range(0,nbin):
        indx=(x>=x_lower[j]) * (x<=x_upper[j])
        x_hist[j]=len(x[indx])
        x_sum[j]=np.sum(x[indx])

    if nsum:
        return x_lower,x_mid,x_upper,x_hist,x_sum
    else:
        return x_lower,x_mid,x_upper,x_hist

def mean_prof(x,y,nbin=10,bintype='fix',nsum=False,median=False):
    
    if bintype=='num':
        x_lower,x_mid,x_upper,x_hist=nbinmaker(x,nbin,nsum)
    else:
        x_lower,x_mid,x_upper,x_hist=binmaker(x,nbin,nsum)

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
                y_bin=np.append(y_bin,y[indx])
                y_sig=np.append(y_sig,0.0)

    return x_bin,y_bin,y_sig

