from math import *
from random import *

ParticleN = 32

fi = open (str(ParticleN),"w")

spacex = 0.9

x=[]
y=[]
z=[]


x.append(0.1)
y.append(0.0)
z.append(0.0)

i = 1
while (i<ParticleN):
    count =0
    xi = uniform(0,0.1)
    yi = uniform(0,0.01)
    x.append (x[0]+spacex*i+xi)
    y.append (y[0]+yi)
    z.append (0.0)
    i=i+1

for i in range (ParticleN):
    print (x[i],y[i],z[i])
    print (x[i],y[i],z[i], file = fi)


for i in range (ParticleN):
    for j in range (ParticleN):
        if (i != j):
            r = sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]))
            if (r < 0.8):
                print ("------------------------------------------------")
                print (i, j, r)
