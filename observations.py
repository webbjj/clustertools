#Routines for making Nbody Clusters look like Observations

#Make radial bins based on an observed dataset
def obsrbinmaker(r,rm,obs_cut):

    r_hist=[]
    r_sum=[]
    r_mean=[]

    if obs_cut=='M30':
        rh=61.800000000000004
        #In arcseconds:
        r_lower=[10.0,20.0,40.0,200.0,250.0,350.0,650.0]
        r_upper=[20.0,40.0,100.0,250.0,350.0,650.0,1000.0]

        for i in range(0,len(r_lower)):
            r_lower[i]=rm*r_lower[i]/rh
            r_upper[i]=rm*r_upper[i]/rh
            r_hist.append(0.0)
            r_sum.append(0.0)

        for i in range(0,len(r)):
            for j in range(0,len(r_lower)):
                if r[i]>=r_lower[j] and r[i]<=r_upper[j]:
                    r_hist[j]+=1.0
                    r_sum[j]+=r[i]


        for i in range(0,len(r_lower)):
            if r_hist[i]>0:
                r_mean.append(r_sum[i]/r_hist[i])
            else:
                r_mean.append((r_lower[i]+r_upper[i])/2.0)


    return r_lower,r_mean,r_upper,r_hist
