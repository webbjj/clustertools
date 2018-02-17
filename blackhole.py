#Extrct black hole and its properties from a cluster
#INCOMPLETE
def get_bhtdf(base):

    filename=base+'bhtdf.dat'
    bhfile=open(filename,'r')
    data=bhfile.readline().split()
    
    t=[]
    x=[]
    y=[]
    z=[]
    r=[]
    v=[]
    e=[]
    A=[]
    kin=[]
    pot=[]
    lx=[]
    ly=[]
    lz=[]
    l=[]

    while len(data) > 1:

        t.append(float(data[0]))
        x.append(float(data[25]))
        y.append(float(data[26]))
        z.append(float(data[27]))
        r.append(float(data[1]))
        v.append(float(data[2]))
        
        e.append(float(data[3]))
        A.append(float(data[5]))
        kin.append(float(data[16]))
        pot.append(float(data[17]))


        lx.append(float(data[63]))
        ly.append(float(data[64]))
        lz.append(float(data[65]))
        l.append(float(data[66]))
                  
        data=bhfile.readline().split()

    bh=BlackHole(t,x,y,z,r,v)
    bh.energies(e,A,kin,pot)
    bh.amomenta(lx,ly,lz,l)

    return bh

class BlackHole(cluster):

    def __init__(self,t,x,y,z,r,v):
        self.t=t
        self.x=x
        self.y=y
        self.z=z
        self.r=r
        self.v=v

    def energies(self,e,kin,pot,A):
        self.e=e
        self.e=e
        self.e=e
        self.A=A

    def amomenta(self,lx,ly,lz,l):
        self.lx=lx
        self.ly=ly
        self.lz=lz
        self.l=l




