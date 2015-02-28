# -*- coding: UTF-8 -*-
"""
  Python implementation of the ReFH  
  Copyright (C) 2013 Neil Nutt
  
  neilnutt [at] googlemail.com

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""
'''v0.0.1 early development, unchecked.

'''

from math import exp
from math import log

def calc_tp(propwet,dplbar,urbext,dpsbar):
    """Function calculates and returns time to peak"""
    tp = 1.56 * ( propwet ** (-1.09)) * (dplbar ** 0.6) * (( 1 + urbext)** (-3.34)) * ( dpsbar ** (-0.28))
    return tp

def calc_bl(bfihost,dplbar,propwet,urbext):
    """Function calculates and returns baseflow lag"""
    bl = 25.5 * (bfihost**0.47) * (dplbar**0.21) * (propwet**-0.53) * ((1 + urbext)**-3.01) 
    return bl
    
def calc_br(bfihost,propwet):
    """Function calculates and returns baseflow recharge"""
    br = 3.75 * (bfihost**1.08) * (propwet**0.36)
    return br
    
def calc_estCriticalDuration(tp,saar):
    """Function calculates and returns the estimated critical duration"""
    d = tp * ( 1.0 + (saar/1000.0))
    return d
    
def calc_cmax(bfihost,propwet):
    """Function calculates and returns maximum soil moisture"""
    cmax = 596.7 * (bfihost**0.95) * (propwet**-0.24)
    return cmax
    
def calc_designRProfile(steps,season,depth):
    """Function calculates and returns the rainfall profile"""
    profile = list()
    depths = list()
    
    if season == "summer":
        a = 0.1
        b = 0.815
    elif season == "winter":
        a = 0.060
        b = 1.026
    else:
        print "Invalid season"
        return

    def proportion(x,a,b):
        z = x * b
        y = (1 - a **z) / (1-a)
        return y
    
    stepCount = 1.0
    x = stepCount/steps
    profile.append(proportion(x,a,b))

    while True:
        stepCount +=  2.0
        x = stepCount/steps
        rainfallProportion = (proportion(x,a,b)-sum(profile))/2.0
        profile.append(rainfallProportion)
        profile.insert(0,rainfallProportion)
        x += 2
        if stepCount >= steps: break
    
    for i in profile:
        depths.append(depth*i)
    
    return profile,depths
    
def calc_cini(season,cmax,bfihost,propwet):
    """Function calculates and returns initial soil moisture"""
    if season == "summer":
        cini = 0.5*cmax* ( 0.09 - 0.82*bfihost + 0.43*propwet)
    elif season == "winter":
        cini = 0.5*cmax* ( 1.20 - 1.07*bfihost + 0.82*propwet)
    else:
        print "Invalid season"
        return
    
    if cini < 0.0:
        cini = 0.0
    
    return cini
    
def calc_alpha(season,T):
    """Function calculates and returns the seasonal adjustment factor for cini"""
    if T < 5:
        alpha = 1.0
    elif season == "summer":
        alpha = 1.444 * T**-0.182
    elif season == "winter":
        alpha = 1.166 * T**-0.073
    else:
        print "Invalid season"
        return
    
    return alpha
    
def calc_runoff(depths,cmax,cini,alpha,intervals):
    """Function calculates and returns the basic pr, netR, sm"""
    prs = list()
    netRs = list()
    sms = list()
    interval = list()

    sm = float(cini) * alpha

    sms.append(sm)
    
    for i in range(intervals):
        interval.append(i + 1)
        pr = ( sm / cmax ) + ( depths[i] / (2 * cmax))
        
        pr = min(pr,1.0)
        
        netR = pr * depths[i]
        sm = sm + depths[i]
        
        sm = min(sm,cmax)
        
        prs.append(pr)
        netRs.append(netR)
        sms.append(sm)
    
    results=dict()
    results['interval'] = interval
    results['percentageRunoff'] = prs
    results['netRainfall'] = netRs
    results['soilMoisture'] = sms
    
    return results
    
def calc_ibaseflow(season,cini,saar,area):
    """Function calculates and returns the initial baseflow"""
    if season == "summer":
        ibf = (33.9 * ( cini - 85.4) + 3.14 * saar) * area / 100000.0
    elif season == "winter":
        ibf = (63.8 * ( cini - 120.8) + 5.54 * saar) * area / 100000.0
    else:
        print "Invalid season"
        return
    
    if ibf < 0.0:
        ibf = 0.0
    
    return ibf
    
    
    
def calc_ks(dt,br,bl):
    """Function calculates and returns the k baseflow parameters"""
    k3 = exp((-dt/bl)*(1.0+br))
    k1 = (br / (1.0 +br))     *    ( ((bl/dt)*(1.0-k3) / (1.0+br)) -k3 )
    k2 = (br / (1.0 +br))     *    ( (1.0 - bl*(1.0-k3)/(dt*(1.0+br))) )
    
    # k3 = exp((-dt/bl)*(1.0+br))
    # k1 = br *( (bl/dt) * (1.0-k3) -k3)
    # k2 = br *( 1.0 - (1.0-k3) * (bl/dt))

    return k1,k2,k3
    
def calc_baseflow(k1,k2,k3,qt,previous_qt,previous_z):
    """Function calculates and returns the next baseflow"""

    # print previous_qt,qt,previous_z
    
    # z = k1 * previous_qt  +  k2*qt +  k3*previous_z # this is the equation in the refh report
    # z = ( 1.0/(1.0-k2)) * (k1*previous_qt  +k2*qt  + (k1 +k3 )*previous_z)

    
    z = ((k1+k3)*previous_z+k2*qt+k1*previous_qt)/(1.0-k2)  ### !!!!THIS IS DIFFERENT TO WHAT IS IN THE REFH REPORT!!!!  (eqn comes from sheet supplied by Konrad Adams
    

    # # block that sort of worked for NN (February 2013) when trying to figure out the error in the Refh report
    # # decay = exp(-dt/bl)
    # # recharge=dt*(br/bl)
    # # print previous_z,decay,qt,recharge,qt*recharge
    # # z = previous_z*decay+qt*recharge
    
    
    if z < 0.0:
        z = 0.0
     
    return z 
    
def calc_seasonCorrectFactor(saar,duration,season):
    """Function calculates and returns the seasonal correction factor"""
    
    if season == "summer":
        if duration <=1.0:
            alpha = -8.03E-5
            beta = 1.04
        elif duration >= 24.0:
            alpha = -10.26E-5
            beta = 1.05
        else:
            x = (1.0, 2.0, 6.0, 24.0)
            a = (-8.03E-5,-6.87E-5,-4.93E-5,-10.26E-5)
            b= (1.04,1.03,1.02,1.05)

            alpha = -3.89E-07*duration**2 + 8.68E-06 * duration - 8.68E-05
            beta = 2.34E-04 * duration **2 - 5.32E-03 * duration + 1.04E+00
        scf = alpha * saar + beta

    elif season == "winter":
        if duration <=1.0:
            alpha = 0.0004
            beta = 0.4
        elif duration >= 24.0:
            alpha =0.0011
            beta = 0.5333
        else:
            x = (1.0, 2.0, 6.0, 24.0)
            a = (0.0004,0.0006,0.0009,0.0011)
            b= (0.4,0.4454,0.4672,0.5333)

            alpha = 2.21E-04*log(duration) + 4.37E-04
            beta = 4.07E-01*duration**8.45E-02
        scf = ( 1.0 - exp(-alpha*saar) ) ** beta
    
    return scf
    
def calc_runoffRouting(area,tp,dt):
    """Function calculates and returns the Unit Hydrograph"""
    t=0
    up = 0.65
    uk = 0.8
    scale = area / (3.6*tp)
    iuh = list()
    uh = list()
    output = dict()
    
    tbt = 2.0 * tp /up
    uc=up *(tbt-2.0*tp) / (tbt-tp)
    tb= tp*(1.0+2.0*(1.0-up)/(uk*uc ))
    
    # slope and intercept of uh recession
    # line 2
    m2 = ( uk * uc - up ) / tp
    c2 = up - m2 * tp
    # print 'u = ',m2,'t +',c2
    # line 3
    m3 = (0.0 - uk * uc) / (tb - 2 *tp)
    c3 = -m3 * tb
    # print 'u = ',m3,'t +',c3
    
    st_previous = 0.0
    
    while True:
        if t <= tp:
            u = up * t /tp
            st = 0.5 * t * u
            marker = 1
        elif t<=2*tp:
            u = m2 * t + c2
            st = (0.5 * tp * up) + (0.5 * m2 * t**2 + c2 * t) - ( 0.5 * m2 * tp**2 + c2 * tp)
            constant = st
            marker = 2
        elif t>2*tp:
            u = m3 * t + c3
            st = (0.5 * tp * up) + (0.5 * (up +uk*uc)*tp) + (0.5 * m3 * t**2 + c3 * t)- ( 0.5 * m3 * (2.*tp)**2 + c3 * (2.*tp))
            marker = 3
        
        if u < 0.0:
            break        
        
        ds = (st - st_previous)/dt
        
        
        iuh.append(u)
        uh.append(scale*ds)
        
        # print t,u,scale*u,ds*scale,marker
        
        t = t+dt
        st_previous = st
        
    output['iuh']=iuh
    output['uh']=uh
    
    return output
    
def calc_hydrographs(runoffSeries,uhSeries,dt,br,bl,numberSteps,ibf,k1,k2,k3):
    """Function calculates and returns the hydrographs"""
    output=dict()
    q = list()
    z = list()
    totalflow = list()
    z.append(ibf)
    
    for i in range(len(2*runoffSeries)+1):    
        q.append(0)
    for r in runoffSeries:
        i =runoffSeries.index(r)
        for u in uhSeries:
            if i < len(q):
                q[i]=u*r+q[i]
            else:
                q.append(u*r)
            i = i+1
    output['rapidRuoff']=q
    
    i = 0
    for rapid in output['rapidRuoff']:
        #i = output['rapidRuoff'].index(rapid)
        if i == 0:
            previous_qt = 0.0
        else:
            previous_qt = output['rapidRuoff'][i-1]
        baseflow= calc_baseflow(k1,k2,k3,rapid,previous_qt,z[i-1])
        z.append(baseflow)
        i=i+1
        
    output['baseflow']=z
    
    for i in range(numberSteps):
        totalflow.append(output['rapidRuoff'][i]+output['baseflow'][i])
        #print i,output['rapidRuoff'][i],output['baseflow'][i],output['rapidRuoff'][i]+output['baseflow'][i]
    
    output['totalflow'] = totalflow
    
    return output
    
def optimiseDuration(cds,rp):
    # cds can be a dictionary of cds or a path to the cds.csv
    # load the catchment descriptors if a dictionary hasn't been passed
    if type(cds) != type(dict()):
        import cds_reader
        cds = cds_reader.csvCds(cds)
        
    durationWithPeak = list()
    

    tp = calc_tp(cds['PROPWET'],cds['DPLBAR'],cds['URBEXT1990'],cds['DPSBAR'])
    dt = tp/4
    d = calc_estCriticalDuration(tp,cds['SAAR'])
    # print "Recommended design storm duration: "+ str(d) +"hrs"
    
    for factor in range(5,75,2):
        duration = dt * factor
        # duration = d
        # print "Adopted design storm duration: "+ str(duration) +"hrs"
        
        
        # print "Adopted timestep: "+ str(dt) +"hrs"
        
        pointRainfall = pointRainfallDepth(cds,duration,rp)
        # print 'Point rainfall depth: '+str(pointRainfall)+'mm'
        
        arf = areaReductionFactor(cds['AREA'],duration)
        # print 'Area reduction factor: ' + str(arf)
        
        rainfall = arf * pointRainfall
        # print 'Catchment average rainfall depth: ' +str(rainfall)+'mm'
        
        if cds['URBEXT1990'] < 0.125:
            season = "winter"
        else:
            season = "summer"
        # print "Season is: "+ season
        
        scf = calc_seasonCorrectFactor(cds['SAAR'],duration,season)
        # print "Seasonal correction factor: "+ str(scf)
        
        designRainfall = pointRainfall*arf*scf
        # print "Design rainfall: "+ str(designRainfall)+"mm"
        
        cmax = calc_cmax(cds['BFIHOST'],cds['PROPWET'])
        # print "Cmax: "+ str(cmax)+"mm"
        
        alpha = calc_alpha(season,rp)
        # print "Alpha: "+ str(alpha)+"mm"

        cini = calc_cini(season,cmax,cds['BFIHOST'],cds['PROPWET'])*alpha
        # print "Cini: "+ str(cini)+"mm"  ### The spreadsheet is double counting alpha!!!!
        
        # print "Tp: "+ str(tp)+"hrs"
        
        bl = calc_bl(cds['BFIHOST'],cds['DPLBAR'],cds['PROPWET'],cds['URBEXT1990'])
        # print "BL: "+ str(bl)+"hrs"
        
        br = calc_br(cds['BFIHOST'],cds['PROPWET'])
        # print "BR: "+ str(br)+"hrs"
        
        ibf = calc_ibaseflow(season,cini,cds['SAAR'],cds['AREA'])
        # print "Initial baseflow: "+ str(ibf)+"cumecs"
        
        k1,k2,k3 = calc_ks(dt,br,bl)
        # print "k1: "+str(k1)
        # print "k2: "+str(k2)
        # print "k3: "+str(k3)
        

        
        steps = int(duration/dt)
        
        rainfallProfile = calc_designRProfile(steps,season,designRainfall)

        
        sm = calc_runoff(rainfallProfile[1],cmax,cini,alpha,steps)
        rainfallIntensity = max(sm['netRainfall'])/dt
            
        # print
        # print "Step, hr, rain (mm), pr (%), runoff(mm), soil moisture (mm)"
        # for i in range(len(rainfallProfile[1])):
            # try:
                # print i,i*dt,rainfallProfile[1][i],sm['percentageRunoff'][i],sm['netRainfall'][i],sm['soilMoisture'][i]
            # except:
                # pass
        
        # print
        # print
        uh = calc_runoffRouting(cds['AREA'],tp,dt)['uh']

        hydrographs = calc_hydrographs(sm['netRainfall'],uh,dt,br,bl,steps*2,ibf,k1,k2,k3)
        
        maxFlow = max(hydrographs['totalflow'])
        entry = [duration,designRainfall,rainfallIntensity,maxFlow]
        
        durationWithPeak.append(entry)
    
    
    # for i in range(steps*2):
        # print i,i*dt,hydrographs['rapidRuoff'][i],hydrographs['baseflow'][i],hydrographs['totalflow'][i]

    return durationWithPeak
    
if __name__ == '__main__':
    printAll = True
    optimiseDur = True
    from ddf import *
    
    if optimiseDur == True:
        durationWithPeaks = optimiseDuration('Ettrick.csv',200)
        for dur,depth,intensity,peak in durationWithPeaks:
            print dur,depth,intensity,peak
        
    if printAll == True:
        
        import cds_reader
        cds = cds_reader.csvCds('Ettrick.csv')
        
        from ddf import *
        tp = calc_tp(cds['PROPWET'],cds['DPLBAR'],cds['URBEXT1990'],cds['DPSBAR'])
        d = calc_estCriticalDuration(tp,cds['SAAR'])
        print "Recommended design storm duration: "+ str(d) +"hrs"
        duration = 20.5
        
        print "Adopted design storm duration: "+ str(duration) +"hrs"
        
        dt = 0.5
        print "Adopted timestep: "+ str(dt) +"hrs"
        
        rp = 200.0
        pointRainfall = pointRainfallDepth(cds,duration,rp)
        print 'Point rainfall depth: '+str(pointRainfall)+'mm'
        
        arf = areaReductionFactor(cds['AREA'],duration)
        print 'Area reduction factor: ' + str(arf)
        
        rainfall = arf * pointRainfall
        print 'Catchment average rainfall depth: ' +str(rainfall)+'mm'
        
        if cds['URBEXT1990'] < 0.125:
            season = "winter"
        else:
            season = "summer"
        print "Season is: "+ season
        
        scf = calc_seasonCorrectFactor(cds['SAAR'],duration,season)
        print "Seasonal correction factor: "+ str(scf)
        
        designRainfall = pointRainfall*arf*scf
        print "Design rainfall: "+ str(designRainfall)+"mm"
        
        cmax = calc_cmax(cds['BFIHOST'],cds['PROPWET'])
        print "Cmax: "+ str(cmax)+"mm"
        
        alpha = calc_alpha(season,rp)
        print "Alpha: "+ str(alpha)+"mm"

        cini = calc_cini(season,cmax,cds['BFIHOST'],cds['PROPWET'])*alpha
        print "Cini: "+ str(cini)+"mm"  ### The spreadsheet is double counting alpha!!!!
        
        print "Tp: "+ str(tp)+"hrs"
        
        bl = calc_bl(cds['BFIHOST'],cds['DPLBAR'],cds['PROPWET'],cds['URBEXT1990'])
        print "BL: "+ str(bl)+"hrs"
        
        br = calc_br(cds['BFIHOST'],cds['PROPWET'])
        print "BR: "+ str(br)+"hrs"
        
        ibf = calc_ibaseflow(season,cini,cds['SAAR'],cds['AREA'])
        print "Initial baseflow: "+ str(ibf)+"cumecs"
        
        k1,k2,k3 = calc_ks(dt,br,bl)
        print "k1: "+str(k1)
        print "k2: "+str(k2)
        print "k3: "+str(k3)
        

        
        steps = int(duration/dt)
        rainfallProfile = calc_designRProfile(steps,season,designRainfall)
        
        sm = calc_runoff(rainfallProfile[1],cmax,cini,alpha,steps)
        
        
        print
        print "Step, hr, rain (mm), pr (%), runoff(mm), soil moisture (mm)"
        for i in range(len(rainfallProfile[1])):
            try:
                print i,i*dt,rainfallProfile[1][i],sm['percentageRunoff'][i],sm['netRainfall'][i],sm['soilMoisture'][i]
            except:
                pass
        
        print
        print
        uh = calc_runoffRouting(cds['AREA'],tp,dt)['uh']
        
        hydrographs = calc_hydrographs(sm['netRainfall'],uh,dt,br,bl,steps*2,ibf,k1,k2,k3)

        
        for i in range(steps*2):
            print i,i*dt,hydrographs['rapidRuoff'][i],hydrographs['baseflow'][i],hydrographs['totalflow'][i]