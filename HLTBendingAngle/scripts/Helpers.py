from ROOT import *
import math
import array
from math import log10, floor

# this look-up table is *ONLY* for ME11!!!
# first number is slope; second number if intersept
slope_inter = {
    1.5 : [3.69, 1.01],
    1.6 : [3.127, 4.283], 
    1.7 : [1.34, 7.131],
    1.8 : [5.056, 6.566],
    1.9 : [1.439, 7.424],
    2.0 : [5.368, 3.792], 
    2.1 : [1.156, 5.511],
    2.2 : [4.61, -0.3263], 
    2.3 : [0.6751, -1.342], 
    2.4 : [0.687, 1.013]
}

def has_csc():
    return TCut("nlayerscsc>4")
    
def dxy(dxy_min, dxy_max):
    return TCut("%f < abs(dxy_csc) && abs(dxy_csc)<%f"%(dxy_min, dxy_max))

def eta_sh_st1(eta_min, eta_max):
    return TCut("%f < abs(csc_gp_eta) && abs(csc_gp_eta)< %f"%(eta_min, eta_max))
    
def same_direction_cut():
    return TCut("Lxy_csc <0 && pzvz_csc<0")
    
def has_reco_pt():
    return TCut("calculated_p_csc_10_15==0")
    
def pt_from_povercosh(eta_min, sim_pt, var=12):
    slope = slope_inter[eta_min][0]
    inter = slope_inter[eta_min][1]
    return TCut("(1/abs(csc_bending_angle_%d) + %f) / %f > %f"%(var, inter, slope, sim_pt))

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