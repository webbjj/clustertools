#Extrct key cluster properties (Output modelled after extrct.f from Jarrod Hurley)
    def extrct(self):

        #Shift to clustercentric origin and to pc if not so already:
        origin0=self.origin
        units0=self.units
        if origin0 !='cluster':
            xvshift(self,self.xgc,self.ygc,self.zgc,self.vxgc,self.vygc,self.vzgc,'cluster')
        if units0 !='realpc':
            kpctopc(self)

        self.mgc=0.0
        self.mc=0.0
        self.mrt=0.0
        for i in range(0,self.ntot):
            if self.r[i]<self.rc:
                self.mc=self.mc+self.m[i]
            if self.r[i]<=self.rtide:
                self.mgc=self.mgc+self.m[i]
            else:
                self.mrt=self.mrt+self.m[i]

        rorder=sorted(range(0,self.ntot),key=lambda k:self.r[k])

        self.rmax=self.r[rorder[-1]]

        #Find 10% Lagrange radius
        msum=0.0
        for i in range(0,self.ntot):
            msum+=self.m[rorder[i]]
            if msum>=0.1*(self.mgc):
                self.r10=self.r[rorder[i]]
                self.n10=i
                break
        
        #Find 50% Lagrange radius
        msum=0.0
        self.r10=0.0
        self.r50=0.0
        self.num1=0
        self.numrc=0

        for i in range(0,self.ntot):
            msum+=self.m[rorder[i]]

            if msum>=0.1*(self.mgc) and self.r10==0.0:
                self.r10=self.r[rorder[i]]
                self.n10=i

            if msum>=0.5*(self.mgc) and self.r50==0.0:
                self.r50=self.r[rorder[i]]
                self.n50=i

        #Find number of stars within 1 pc
            if self.r[rorder[i]] <= 1.0:
                self.num1+=1

            if self.r[rorder[i]] <= self.rc:
                self.numrc=i


#self.trh=(float(self.ntot)/math.log(0.4*float(self.ntot),10))*0.858*math.sqrt(math.pow(self.r50,3.0)/self.mgc)

        self.trh=Half_Mass_Relaxation_Time(self.ntot,self.r50,self.mgc+self.mrt)

        print (self.ns,self.nb,self.tphys,self.trh,self.mgc+self.mrt,self.mc,self.mrt,self.rmax,self.r50,self.r10,self.rc,self.n50,self.num1,self.n10,self.numrc,self.rtide)
        #Shift back to original origin
        
        if origin0 == 'galaxy' and self.origin=='cluster':
            xvshift(self,self.xgc,self.ygc,self.zgc,self.vxgc,self.vygc,self.vzgc,'cluster')
        if units0 =='realkpc' and self.units=='realpc':
            pctokpc(self)
