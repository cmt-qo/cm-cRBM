import numpy as np


def jv(s):
  if s==0:
    return 3
  if s==1:
    return 2
  if s==2:
    return 1
  if s==3:
    return 0


def jp(s):
 if s==0:
   return 7
 if s==1:
   return 6
 if s==2:
   return 5
 if s==3:
   return 4

print "b:"
for ai in xrange(2):
  for s in xrange(4):
    print "ai,i",ai,s
    if ai==0:
      print "j",jv(s)
    else:
      print "j",jp(s)

print " "
print " "
print "W:"
sym=np.loadtxt("symmetries_L2.txt")
sym=np.transpose(sym)
for ai in xrange(2):
  for i in xrange(8):
    print "ai,i",ai,i,":"
    i0=int(sym[i][0])
    i1=int(sym[i][1])
    i2=int(sym[i][2])
    i3=int(sym[i][3])
    if ai==0:
      print "i,j",i0*8+jv(0),i1*8+jv(1),i2*8+jv(2),i3*8+jv(3)
    if ai==1:
      print "i,j",i0*8+jp(0),i1*8+jp(1),i2*8+jp(2),i3*8+jp(3)
    print " "
