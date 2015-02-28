"""
  Halcrow Group Limited - Python module to assist
  with running the FEH DDF.
  
  Copyright (C) 2013 Neil Nutt for Halcrow Group Limited
  
  neil [dot] nutt [at] ch2m.com

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

'''v0.0.1 early development, unchecked.'''

from math import log,exp
def pointRainfallDepth(cds,duration,returnPeriod):
    
    c=cds['C']
    d1=cds['D1']
    d2=cds['D2']
    d3=cds['D3']
    e=cds['E']
    f=cds['F']
    
    y = -log(-log(1-returnPeriod**(-1.0)))
    
    lnR12 = (c*y+d1)*log(12.0)+e*y+f
    lnR48 = lnR12 + (c*y+d2)*(log(48.0)-log(12.0))
    
    if duration <= 12.0:
        lnR = (c*y+d1)*log(duration)+e*y+f
    elif duration > 12.0 and duration <= 48.0:
        lnR = lnR12 + (c*y+d2)*(log(duration)-log(12.0))
    elif duration > 48.0:
        lnR = lnR48 + (c*y+d3)*(log(duration)-log(48.0))
        
    pointDepth = exp(lnR)
    
    return pointDepth
    
def areaReductionFactor(area,dur):
    if area < 0.0:
        return None
        
    elif area < 20.0:
        arf1 = 0.40 - 0.0208*log(4.6-log(area))
        arf2 = 0.0394 * area**0.354
    elif area < 100.0:
        arf1 = 0.40 - 0.00382*(4.6-log(area))**2
        arf2 = 0.0394 * area**0.354
    elif area < 500.0:
        arf1 = 0.40 - 0.00382*(4.6-log(area))**2
        arf2 = 0.0627 * area**0.254
    elif area < 1000.0:
        arf1 = 0.40 - 0.0208*log(log(area)-4.60)
        arf2 = 0.0627 * area**0.254
    elif area >= 1000.0:
        arf1 = 0.40 - 0.0208*log(log(area)-4.60)
        arf2 = 0.1050 * area**0.180
    
    arf = 1.0 - arf2*dur**(-arf1)
    
    return arf
    
if __name__ == '__main__':
    testDdf = 1
    testArf = 0
    
    if testDdf == 1:
        #duration = 5.0 #hrs
        returnPeriod = 20.0 #yrs
        durations = [0.5,1.0,11.9,12.0,12.1,47.9,48.0,48.1,500.0]
        returnPeriods = [1.5,2.0,5.0,25.0,100.0,1000.0]
        
        
        cds =dict()
        cds['C'] = -0.02293
        cds['D1'] = 0.41911
        cds['D2'] = 0.51466
        cds['D3'] = 0.34721
        cds['E'] = 0.28183
        cds['F'] = 2.3596
        
        for duration in durations:
            print
            print '--------------'
            print
            for returnPeriod in returnPeriods:
                print returnPeriod,"yr , ",duration,"hr : ",pointRainfallDepth(cds,duration,returnPeriod)
        
    if testArf == 1:
        areas = [1.0,2.0,19.0,21.0,99.0,101.0,499.0,501.0,999.0,1001.0,5000.0]
        durations = [0.5,1.0,3.0,6.0,12.0,24.0]
        for area in areas:
            print '--------'
            print
            for duration in durations :
                print duration,area,areaReductionFactor(area,duration)    