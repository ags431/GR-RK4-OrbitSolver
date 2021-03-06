import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

M = 1.0 #62*1.99*(10**30)

G = 1.0 #6.67*(10**-11)

c = 1.0 #3*(10**8)

rinner = 6.0*G*M/(c**2) #548629.733

hVariation = 1 #set to 1 to automatically vary h with speed; set to 0 to use a fixed h ===!!!!!===
# FOR FIXED H
h = 1.0 #timestep to use when not varying h; proper time
# FOR VARYING H
hscale = 0.1 #the time step to use when v = c
maxH = 10000.0 #maximum H
minH = 0.1 #minimum H

maxOrbits = 2.0

def tPP(t, r, phi, ut, ur, uphi):
    return -2.0*G*M*ut*ur/((1.0-(2.0*G*M)/((c**2)*r))*(r**2)*(c**2))

def rPP(t, r, phi, ut, ur, uphi):
    tempRatio = (1.0-(2.0*G*M)/((c**2)*r))
    #print("temp ratio is", tempRatio)
    return G*M*(ur**2)/(tempRatio*(r**2)*(c**2)) + tempRatio*r*(uphi**2) - tempRatio*M*G*(ut**2)/((c**2)*(r**2))

def phiPP(t, r, phi, ut, ur, uphi):
    return -2*ur*uphi/r
    
def linearSpeed(r, ur, uphi):
    return np.sqrt(ur**2 + (uphi*r)**2)

properi = 0.0
ti = 0.0
ri = 1640450.646 #1.5*(10**11)
phii = 0.0
uti = 1.0/np.sqrt(1-(3.0*M)/(ri))
uri = 0.0
uphii = np.sqrt(G*M/((ri**3)))*uti

propers = [properi]
ts = [ti]
rs = [ri]
phis = [phii]
uts = [uti]
urs = [uri]
uphis = [uphii]

index = 0

def transitOrbit(rstop, urboost, uphiboost):
    
    global index
    global h
    
    startingRadius = rs[index]
    
    urs[index] += urboost
    uphis[index] += uphiboost
    
    #normalize
    gtt = -(1.0 - 2.0*G*M/(rs[index]*(c**2)))
    grr = 1.0/(-gtt)
    gphiphi = rs[index]**2
    uts[index] = (-1.0 - grr*(urs[index]**2) - gphiphi*(uphis[index]**2))/gtt
    
    print("time to loop!")
    
    while True:
        
        if rs[index] < rstop:
            return 0
        elif rs[index] > startingRadius:
            propers.pop()
            ts.pop()
            rs.pop()
            phis.pop()
            uts.pop()
            urs.pop()
            uphis.pop()
            return 1
        elif phis[index] > maxOrbits*2.0*np.pi:
            return 2
        
        if hVariation == 1:
            speed = linearSpeed(rs[index], urs[index], uphis[index])
            h = (hscale*c/speed)**2
            
            if h > maxH:
                h = maxH
            elif h < minH:
                h = minH
                
            print("h", h, "speed", speed)
        
        #using p to denote prime
        #using c to denote current
        #using f to denote final
        
        tc = ts[index]
        rc = rs[index]
        phic = phis[index]
        utc = uts[index]
        urc = urs[index]
        uphic = uphis[index]
        #first run
        utp1 = tPP(tc, rc, phic, utc, urc, uphic)
        urp1 = rPP(tc, rc, phic, utc, urc, uphic)
        uphip1 = phiPP(tc, rc, phic, utc, urc, uphic)
        tp1 = utc
        rp1 = urc
        phip1 = uphic
        #print(utp1, urp1, uphip1, tp1, rp1, phip1)
        #second run
        utp2 = tPP(tc + 0.5*h*(tp1), rc + 0.5*h*(rp1), phic + 0.5*h*(phip1), utc + 0.5*h*utp1, urc + 0.5*h*urp1, uphic + 0.5*h*uphip1)
        urp2 = rPP(tc + 0.5*h*(tp1), rc + 0.5*h*(rp1), phic + 0.5*h*(phip1), utc + 0.5*h*utp1, urc + 0.5*h*urp1, uphic + 0.5*h*uphip1)
        uphip2 = phiPP(tc + 0.5*h*(tp1), rc + 0.5*h*(rp1), phic + 0.5*h*(phip1), utc + 0.5*h*utp1, urc + 0.5*h*urp1, uphic + 0.5*h*uphip1)
        tp2 = utc + 0.5*h*utp1
        rp2 = urc + 0.5*h*urp1
        phip2 = uphic + 0.5*h*uphip1
        #third run
        utp3 = tPP(tc + 0.5*h*(tp2), rc + 0.5*h*(rp2), phic + 0.5*h*(phip2), utc + 0.5*h*utp2, urc + 0.5*h*urp2, uphic + 0.5*h*uphip2)
        urp3 = rPP(tc + 0.5*h*(tp2), rc + 0.5*h*(rp2), phic + 0.5*h*(phip2), utc + 0.5*h*utp2, urc + 0.5*h*urp2, uphic + 0.5*h*uphip2)
        uphip3 = phiPP(tc + 0.5*h*(tp2), rc + 0.5*h*(rp2), phic + 0.5*h*(phip2), utc + 0.5*h*utp2, urc + 0.5*h*urp2, uphic + 0.5*h*uphip2)
        tp3 = utc + 0.5*h*utp2
        rp3 = urc + 0.5*h*urp2
        phip3 = uphic + 0.5*h*uphip2
        #fourth run
        utp4 = tPP(tc + h*(tp3), rc + h*(rp3), phic + h*(phip3), utc + h*utp3, urc + h*urp3, uphic + h*uphip3)
        urp4 = rPP(tc + h*(tp3), rc + h*(rp3), phic + h*(phip3), utc + h*utp3, urc + h*urp3, uphic + h*uphip3)
        uphip4 = phiPP(tc + h*(tp3), rc + h*(rp3), phic + h*(phip3), utc + h*utp3, urc + h*urp3, uphic + h*uphip3)
        tp4 = utc + h*utp3
        rp4 = urc + h*urp3
        phip4 = uphic + h*uphip3
             
        utf = utc + (h/6.0)*(utp1 + 2.0*utp2 + 2.0*utp3 + utp4)
        urf = urc + (h/6.0)*(urp1 + 2.0*urp2 + 2.0*urp3 + urp4)
        uphif = uphic + (h/6.0)*(uphip1 + 2.0*uphip2 + 2.0*uphip3 + uphip4)
        tf = tc + (h/6.0)*(tp1 + 2.0*tp2 + 2.0*tp3 + tp4)
        rf = rc + (h/6.0)*(rp1 + 2.0*rp2 + 2.0*rp3 + rp4)
        phif = phic + (h/6.0)*(phip1 + 2.0*phip2 + 2.0*phip3 + phip4)
        
        properf = propers[index] + h
        
        print(properf, tf, rf, phif, utf, urf, uphif)
        
        propers.append(properf)
        ts.append(tf)
        rs.append(rf)
        phis.append(phif)
        uts.append(utf)
        urs.append(urf)
        uphis.append(uphif)
        index += 1

# polar plot

# ============ ADD TRANSIT ORBITS HERE ==================
# transitOrbit( rStop, urBoost, uphiBoost)
# rStop - the radius to "stop the transit orbit" at; this can be used to make successive transit orbits, in which this would be the radius at which to switch to another orbit
# urBoost - the boost in ur given to the pod to changed to the new transit orbit
# uphiBoost - the boost in uphi given to the pod to change to the new transit orbit
# function will return 0 if the orbit reached rStop normally, 1 if the orbit returned to the original radius again without reaching rStop

# determine ur and uphi impulse from router



# example orbit
ret = transitOrbit(rinner, 0, (-0.997)*uphii)
if ret == 0:
    # successive orbits here
    #ret = transitOrbit(rinner, 0, -uphis[len(uphis)-1]*0.999)
    if ret == 0:
        # each successive orbit is wrapped in another if
        print("placeholder")

# =======================================================

if ret == 0:
    print("hit inner orbit")
if ret == 1:
    print("returned to the starting radius!")
if ret == 2:
    print("reached max orbits")

xs = rs*np.cos(phis)
ys = rs*np.sin(phis)
fancy = plt.subplot(111, projection='3d')
fancy.plot(xs, ys, zs=ts, color='b', linewidth=2)

blackholexs = [0]*len(ts)
blackholeys = blackholexs
fancy.plot(blackholexs, blackholeys, zs=ts, color='black', linewidth=1)

fancy.set_xlabel("x")
fancy.set_ylabel("y")
fancy.set_zlabel("Coordinate Time")

plt.figure()
polar = plt.subplot(111, projection='polar')
polar.plot(phis, rs, color='r', linewidth=2)

plt.figure()
rect = plt.subplot(111)
rect.plot(phis, rs, color='g', linewidth=2)

rect.set_xlabel("phi")
rect.set_ylabel("r")

plt.figure()
timeDifference = range(0, len(ts))
for x in range(0,len(ts)):
    timeDifference[x] = ts[x]-propers[x]
rect2 = plt.subplot(111)
rect2.plot(propers, timeDifference, color='black', linewidth=2)

rect2.set_xlabel("Proper Time")
rect2.set_ylabel("Coordinate/Proper Time Difference")

plt.show()