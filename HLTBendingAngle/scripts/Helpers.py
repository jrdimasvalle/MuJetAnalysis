from ROOT import *
from math import *
import array

def has_delta_y():
    return TCut("has_delta_y > 0")

def has_csc(stx = 1):
    if stx == 1:
        return TCut("nlayerscsc>=4")
    if stx == 2:
        return TCut("nlayers_st2>=4")
    if stx == 3:
        return TCut("nlayers_st3>=4")
    if stx == 4:
        return TCut("nlayers_st4>=4")


def dxy_cut(dxy_min, dxy_max):
    return TCut("%f <= abs(dxy) && abs(dxy)<%f"%(dxy_min, dxy_max))


def eta_sh_cut(eta_min, eta_max, stx = 2):
    if stx ==1:
        return TCut("%f <= abs(csc_gp_eta) && abs(csc_gp_eta)<= %f"%(eta_min, eta_max))
    if stx == 2:
        return TCut("%f <= abs(csc_gp_second_st2) && abs(csc_gp_second_st2)<= %f"%(eta_min, eta_max))
    if stx == 3:
            return TCut("%f <= abs(csc_gp_second_st3) && abs(csc_gp_second_st3)<= %f"%(eta_min, eta_max))
    if stx == 3:
            return TCut("%f <= abs(csc_gp_second_st4) && abs(csc_gp_second_st4)<= %f"%(eta_min, eta_max))
    

def has_csc_sh_pairs(var):
    return TCut("has_csc_%d>0"%var)


def SimTrack_pt_cut(minpt, maxpt):
    return TCut("pt_SimTrack > %d && pt_SimTrack < %d"%(minpt,maxpt))


def endcap_csc(value = 1):
    if value == 0:
        return TCut("endcap_st1 >0 && endcap_st2 >0 && endcap_st3 >0 && endcap_st4 >0")
    
    if value == 1:
        return TCut("endcap_st1 == 1 && endcap_st2 == 1 && endcap_st3 == 1 && endcap_st4 == 1")

    if value == 2:
        return TCut("endcap_st1 == 2 && endcap_st2 == 2 && endcap_st3 == 2 && endcap_st4 == 2")


def ME1X_only(X):
    if X == 1:
        return TCut("(csc_ring == 1  || csc_ring == 4)")

    if X ==2:

        return TCut("csc_ring == 2")

    if X ==3:
        
        return TCut("csc_ring == 3")

def parity_csc(case):
    if case ==0:
        return TCut("csc_chamber%2 == 1 && csc_chamber_st2%2 == 0 && csc_chamber_st3%2 == 0")
    if case ==1:
        return TCut("csc_chamber%2 == 1 && csc_chamber_st2%2 == 1 && csc_chamber_st3%2 == 1")
    if case == 2:
        return TCut("csc_chamber%2 == 0 && csc_chamber_st2%2 == 0 && csc_chamber_st3%2 == 0")
    if case ==3:
        return TCut("csc_chamber%2 == 0 && csc_chamber_st2%2 == 1 && csc_chamber_st3%2 == 1")


def pT_from_Direction(etamin = 1.6, stx = 1, parity_case=0, minptcut=10):

    slope = 1.0
    intercept = 1.0

    if stx ==1:
        if parity_case ==0 :
            if etamin == 1.6:
                slope = 4.53
                intercept = 9.422
            if etamin == 1.8:
                slope= 5.788
                intercept = 9.743
            if etamin == 2.0:
                slope = 8.367
                intercept = 10.22
            if etamin == 2.2:
                slope = 11.02
                intercept = 14.84

        if parity_case ==1 :
            if etamin == 1.6:
                slope = 4.779
                intercept = 9.954
            if etamin == 1.8:
                slope= 6.273
                intercept = 11.91
            if etamin == 2.0:
                slope = 9.315
                intercept =12.21
            if etamin == 2.2:
                slope = 10.34
                intercept = 11.02
                
        if parity_case ==2 :
            if etamin == 1.6:
                slope = 7.677
                intercept =16.82
            if etamin == 1.8:
                slope= 7.726
                intercept = 13.36
            if etamin == 2.0:
                slope = 9.621
                intercept =10.62
            if etamin == 2.2:
                slope = 11.23
                intercept = 13.44

        if parity_case ==3 :
            if etamin == 1.6:
                slope = 7.72
                intercept = 16.91
            if etamin == 1.8:
                slope= 8.643
                intercept = 16.45
            if etamin == 2.0:
                slope = 10.02
                intercept =11.83
            if etamin == 2.2:
                slope = 11.83
                intercept = 17.66

    if stx ==2:

        if parity_case == 0: 
            if etamin == 1.2:
                slope = 2.907
                intercept = 5.906
            if etamin == 1.4:
                slope= 2.6
                intercept = 5.191
        if parity_case == 1: 
            if etamin == 1.2:
                slope = 2.409
                intercept = 5.198
            if etamin == 1.4:
                slope= 2.467
                intercept = 4.397

        if parity_case == 2: 
            if etamin == 1.2:
                slope = 2.301
                intercept = 4.929
            if etamin == 1.4:
                slope= 2.23
                intercept = 3.111

        if parity_case == 3: 
            if etamin == 1.2:
                slope = 2.401
                intercept = 4.758
            if etamin == 1.4:
                slope= 2.383
                intercept = 3.782

    return TCut(" (( 1/abs(csc_bending_angle_12) + %f )/%f ) > %f"%(intercept, slope, minptcut))




def pT_from_Positions(etamin = 1.2, stx = 2, parity_case=0, minptcut=10):

    slope = 1.0
    intercept = 1.0
    prop = 0.0

    if stx ==1:
        if parity_case ==0 :
            prop = 0.6484
            if etamin == 1.6:
                slope = 0.05527
                intercept = 0.08944
            if etamin == 1.8:
                slope= 0.08295
                intercept = 0.1279
            if etamin == 2.0:
                slope = 0.166
                intercept = 0.2158
            if etamin == 2.2:
                slope = 0.4952
                intercept = 0.7103

        if parity_case ==1 :
            prop = 0.3542
            if etamin == 1.6:
                slope = 0.1067
                intercept = 0.1957
            if etamin == 1.8:
                slope= 0.1561
                intercept = 0.2654
            if etamin == 2.0:
                slope = 0.3156
                intercept = 0.4514
            if etamin == 2.2:
                slope = 0.8242
                intercept = 1.071

        if parity_case ==2 :
            prop = 0.5636
            if etamin == 1.6:
                slope = 0.05624
                intercept = 0.08417
            if etamin == 1.8:
                slope= 0.08702
                intercept = 0.1426
            if etamin == 2.0:
                slope = 0.1676
                intercept = 0.2198
            if etamin == 2.2:
                slope = 0.4953
                intercept = 0.7272
                
        if parity_case ==3 :
            prop = 0.3217
            if etamin == 1.6:
                slope = 0.1066
                intercept = 0.2026
            if etamin == 1.8:
                slope= 0.1435
                intercept = 0.2118
            if etamin == 2.0:
                slope = 0.2874
                intercept = 0.4055
            if etamin == 2.2:
                slope = 0.7625
                intercept = 1.075
                
    if stx ==2:
        if parity_case ==0 :
            prop = 1.279
            if etamin == 1.2:
                slope = 0.04784
                intercept = 0.1122
            if etamin == 1.4:
                slope= 0.05424
                intercept = 0.09761
      

        if parity_case ==1 :
            prop = 0.6357
            if etamin == 1.2:
                slope = 0.08274
                intercept = 0.2021
            if etamin == 1.4:
                slope= 0.09063
                intercept = 0.1773
          

        if parity_case ==2 :
            prop = 1.001
            if etamin == 1.2:
                slope = 0.038
                intercept = 0.08345
            if etamin == 1.4:
                slope= 0.04157
                intercept = 0.0617
         
                
        if parity_case ==3 :
            prop = 0.5252
            if etamin == 1.2:
                slope = 0.07391
                intercept = 0.1714
            if etamin == 1.4:
                slope= 0.07838
                intercept = 0.1307
          

                


    return TCut(" (( 1/abs((delta_y_gp_23) - %f*(delta_y_gp_12) )  + %f )/%f ) > %f"%(prop, intercept, slope, minptcut))





def Combined_pT(etamin = 1.6, stx = 1, parity_case=0, minptcut=10):


    slope_phi = 0.0
    intercept_phi = 0.0
    sigma_phi = 0.0


    if parity_case ==0 :
        if etamin == 1.6:
            sigma_phi = 0.03583
            slope_phi = 4.53
            intercept_phi = 9.422
        if etamin == 1.8:
            slope_phi = 5.788
            intercept_phi = 9.743
            sigma_phi = 0.04496
        if etamin == 2.0:
            sigma_phi = 0.05328
            slope_phi = 8.367
            intercept_phi = 10.22
        if etamin == 2.2:
            sigma_phi = 0.05832
            slope_phi = 11.02
            intercept_phi = 14.84

    if parity_case ==1 :
        if etamin == 1.6:
            sigma_phi = 0.03948
            slope_phi = 4.779
            intercept_phi = 9.954
        if etamin == 1.8:
            sigma_phi = 0.04381
            slope_phi= 6.273
            intercept_phi = 11.91
        if etamin == 2.0:
            sigma_phi = 0.0557
            slope_phi = 9.315
            intercept_phi =12.21
        if etamin == 2.2:
            sigma_phi = 0.06099
            slope_phi = 10.34
            intercept_phi = 11.02
            
    if parity_case ==2 :
        if etamin == 1.6:
            sigma_phi = 0.05226
            slope_phi = 7.677
            intercept_phi =16.82
        if etamin == 1.8:
            sigma_phi = 0.04987
            slope_phi= 7.726
            intercept_phi = 13.36
        if etamin == 2.0:
            sigma_phi = 0.06021
            slope_phi = 9.621
            intercept_phi =10.62
        if etamin == 2.2:
            sigma_phi = 0.05825
            slope_phi = 11.23
            intercept_phi = 13.44

    if parity_case ==3 :
        if etamin == 1.6:
            sigma_phi = 0.04902
            slope_phi = 7.72
            intercept_phi = 16.91
        if etamin == 1.8:
            sigma_phi = 0.04995
            slope_phi = 8.643
            intercept_phi = 16.45
        if etamin == 2.0:
            sigma_phi = 0.05696
            slope_phi = 10.02
            intercept_phi =11.83
        if etamin == 2.2:
            sigma_phi = 0.05652
            slope_phi = 11.83
            intercept_phi = 17.66


                



    slope = 0.0
    intercept = 0.0
    prop = 0.0
    sigma= 0.0



    if parity_case ==0 :
        prop = 0.6484
        if etamin == 1.6:
            sigma = 0.01319
            slope = 0.05527
            intercept = 0.08944
        if etamin == 1.8:
            sigma = 0.02154
            slope= 0.08295
            intercept = 0.1279
        if etamin == 2.0:
            sigma = 0.03251
            slope = 0.166
            intercept = 0.2158
        if etamin == 2.2:
            sigma = 0.05515
            slope = 0.4952
            intercept = 0.7103

    if parity_case ==1 :
        prop = 0.3542
        if etamin == 1.6:
            sigma = 0.01014
            slope = 0.1067
            intercept = 0.1957
        if etamin == 1.8:
            sigma = 0.01997
            slope= 0.1561
            intercept = 0.2654
        if etamin == 2.0:
            sigma = 0.03619
            slope = 0.3156
            intercept = 0.4514
        if etamin == 2.2:
            sigma = 0.05695
            slope = 0.8242
            intercept = 1.071

    if parity_case ==2 :
        prop = 0.5636
        if etamin == 1.6:
            sigma = 0.008583
            slope = 0.05624
            intercept = 0.08417
        if etamin == 1.8:
            sigma = 0.02352
            slope= 0.08702
            intercept = 0.1426
        if etamin == 2.0:
            sigma = 0.03006
            slope = 0.1676
            intercept = 0.2198
        if etamin == 2.2:
            sigma = 0.05692
            slope = 0.4953
            intercept = 0.7272
            
    if parity_case ==3 :
        prop = 0.3217
        if etamin == 1.6:
            sigma = 0.006731
            slope = 0.1066
            intercept = 0.2026
        if etamin == 1.8:
            sigma = 0.02152
            slope= 0.1435
            intercept = 0.2118
        if etamin == 2.0:
            sigma = 0.03513
            slope = 0.2874
            intercept = 0.4055
        if etamin == 2.2:
            sigma = 0.05173
            slope = 0.7625
            intercept = 1.075



    return TCut("1/(( (%f/(1/abs(csc_bending_angle_12) + %f ))/(%f*%f) +       (%f/( 1/abs((delta_y_gp_23) - %f*(delta_y_gp_12)) + %f))/(%f*%f) )/(1/(%f*%f) + 1/(%f*%f) ) ) > %f "%( slope_phi, intercept_phi, sigma_phi, sigma_phi, slope, prop, intercept, sigma, sigma, sigma, sigma, sigma_phi, sigma_phi, minptcut ))

    

#_____________________________________________________________

def FitHistoFunction68(b1, xBins, xminBin, xmaxBin, yBins, yminBin, ymaxBin, printa): 
    
    r1 = TH1F("r1","1/abs(\Delta\phi) vs p^{Sim}",xBins,xminBin,xmaxBin)
    for x in range(0,xBins):

        if (printa > 0):
            print "*********** For bin x: %d **********************"%x
            

        # Find the total number of frequencies
        totalfreq = b1.Integral(x,x,0,yBins+1)

        # Calculate half of the frequencies
        med = 0
        
        if (totalfreq%2 ==1) :
            med = (totalfreq-1)/2 + 1            # Set the value of the median

        if (totalfreq%2 ==0):
            med = (totalfreq/2) + 0.5               # This might need to be added 0.5

        temporal = 0
        midbin = 0
        
        for m in range (0,yBins+1):
                temporal = b1.Integral(x,x,0,m)

                if (temporal >= med):
                    midbin = m              # Break once I get to the median
                    break

       
        if (printa > 0):

                print "suma: ",totalfreq
                print "Midbin: ",midbin
                print "mediana count: ",temporal


        
        # midbin is the actual value to be stored in (x, midbin) histogram.    
        # Find the error above the median
        
        sumerrup = 0                       # Sum of events up
        binerrup = 0                       # Bin which has the 34% of events
        
        for k in range (midbin, yBins+1): # Looping over the midbin up to 10000 
            sumerrup = b1.Integral(x,x,midbin,k)
  
            if (sumerrup >= 0.33*totalfreq):     # If the summ is bigger or equal to 34% of the total number of entries, break
                binerrup = k
                break
        
        sumerrlow = 0
        binerrlow = 0

        
        for r in range (0, midbin):
            sumerrlow = b1.Integral(x,x,midbin,midbin-r)

            if (sumerrlow >= 0.33*totalfreq):
                binerrlow = midbin-r        # Store the bin which has the 34% value 
                break
        
        # error is the difference averaged on the bin low and bin up
 
        errorbin = abs(binerrup - binerrlow)/2
        if (totalfreq > 0):
            errorbin = errorbin / sqrt (totalfreq)
        
        if (printa > 0):
            print " X position: ",x
            print " Bin y: ",midbin
            print " Error: ",errorbin
            print " Error low: ",binerrlow
            print " Sum err low: ",sumerrlow
            print " Error high: ",binerrup
            print " Sum err high: ",sumerrup

        # Store in a histogram the values of (x, midbin) with an error given by errorbin

        
        if errorbin ==0:
            errorbin == 1

            
        scale = yBins / ymaxBin
        if ymaxBin < yBins:
            scale = ymaxBin / yBins


        r1.SetBinContent(x, midbin*scale)
        r1.SetBinError(x, errorbin*scale)

    return r1                               #Return the histogram 1D 
                             #Return the histogram 1D 


#_______________________________________________________________________________
def get_1D(p, title, h_name, h_bins, to_draw, cut, opt = "", color = kBlue):
    gStyle.SetStatStyle(0)
    gStyle.SetOptStat(11111111)
    #nbins = len(xbins)
    #h = TH1F("h_name", "h_name", nbins, xbins);
    p.tree.Draw(to_draw + ">>" + h_name + h_bins, cut)
    h = TH1F(gDirectory.Get(h_name).Clone(h_name))
    if not h:
        sys.exit('%s does not exist'%(to_draw))
    h.SetTitle(title)
    h.SetLineWidth(2)
    h.SetLineColor(color)    
    h.SetMinimum(0.)
    SetOwnership(h, False)
    return h

#_______________________________________________________________________________
def AND(*arg):
    """AND of any number of TCuts in PyROOT"""
    length = len(arg)
    if length == 0:
        print "ERROR: invalid number of arguments"
        return
    if length == 1:
        return arg[0] 
    if length==2:
        return ANDtwo(arg[0],arg[1])
    if length>2:
        result = arg[0]
        for i in range(1,len(arg)):
            result = ANDtwo(result,arg[i])
        return result

#_______________________________________________________________________________
def OR(*arg):
    """OR of any number of TCuts in PyROOT"""
    length = len(arg)
    if length == 0:
        print "ERROR: invalid number of arguments"
        return
    if length == 1:
        return arg[0] 
    if length==2:
        return ORtwo(arg[0],arg[1])
    if length>2:
        result = arg[0]
        for i in range(1,len(arg)):
            result = ORtwo(result,arg[i])
        return result

#_______________________________________________________________________________
def ANDtwo(cut1,cut2):
    """AND of two TCuts in PyROOT"""
    if cut1.GetTitle() == "":
        return cut2
    if cut2.GetTitle() == "":
        return cut1
    return TCut("(%s) && (%s)"%(cut1.GetTitle(),cut2.GetTitle()))


#_______________________________________________________________________________
def ORtwo(cut1,cut2):
    """OR of two TCuts in PyROOT"""
    if cut1.GetTitle() == "":
        return cut2
    if cut2.GetTitle() == "":
        return cut1
    return TCut("(%s) || (%s)"%(cut1.GetTitle(),cut2.GetTitle()))
