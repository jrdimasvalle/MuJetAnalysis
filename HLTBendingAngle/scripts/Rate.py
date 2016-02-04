from ROOT import *
from Helpers import *

varsh = 12
stx = "ME11 - ME2"
test = 15


#____________________________________
def getIntegral(hist):
  nBins = hist.GetXaxis().GetNbins()
  min1 = hist.GetXaxis().GetXmin()
  max1 = hist.GetXaxis().GetXmax()
  integral = TH1F("integral", "integral", nBins, min1, max1)
  integral.SetStats(0)
  hist.Sumw2()
  global1=0
  
  for i in range(min1,nBins):
    total = 0

 
    for k in range(i,nBins):
      
      total = total + hist.GetBinContent(k)

    integral.SetBinContent(i,total);
 
  return integral



#_____________________________________________________________________________________________________________________________________________________________________________________________________________________________
def denominator(dxy_min = 0, dxy_max = 5, eta_min = 1.6, eta_max = 1.7, minpt=0, maxpt = 500, stx=1, parity_case = 0):

        return AND(has_delta_y(), has_csc(1), has_csc(2), has_csc(3), dxy_cut(dxy_min, dxy_max), eta_sh_cut(eta_min, eta_max, 2), has_csc_sh_pairs(12), has_csc_sh_pairs(23), SimTrack_pt_cut(minpt, maxpt), ME1X_only(stx), parity_csc(parity_case))

#_____________________________________________________________________________________________________________________________________________________________________________________________________________________________
def Comparison(file, dir, xaxis, yaxis, x_bins, y_bins, etamin, etamax,  minptcut, maxptcut, ctau, stx, parity_case):
       
    c1 = TCanvas()
    c1.SetGridx()
    c1.SetGridy()
    c1.SetTickx()
    c1.SetTicky()
    gStyle.SetOptFit(0111)
    gStyle.SetOptStat(0)
    #c1.SetLogz()

    
    gStyle.SetStatY(0.25)
    gStyle.SetStatX(0.90)
    gStyle.SetStatW(0.1729)
    gStyle.SetStatH(0.12)

    f = TFile(file)

    tt10 = f.Get(dir)
    tt11 = f.Get(dir)
    tt12 = f.Get(dir)
    tt13 = f.Get(dir)           
    xBins = int(x_bins[1:-1].split(',')[0])
    xminBin = float(x_bins[1:-1].split(',')[1])
    xmaxBin = float(x_bins[1:-1].split(',')[2])
    yBins = int(y_bins[1:-1].split(',')[0])
    yminBin = float(y_bins[1:-1].split(',')[1])
    ymaxBin = float(y_bins[1:-1].split(',')[2])

    todrawb1 = "%s"%yaxis+":"+"%s>>b1"%xaxis
    todrawb2 = "%s"%yaxis+":"+"%s>>b2"%xaxis
    todrawb3 = "%s"%yaxis+":"+"%s>>b3"%xaxis

    e01 = TH1F("e01", "e01", xBins, xminBin, xmaxBin)
    e02 = TH1F("e02", "e02", xBins, xminBin, xmaxBin)

    e1 = TH1F("e1", "e1", xBins, xminBin, xmaxBin)
    e2 = TH1F("e2", "e2", xBins, xminBin, xmaxBin)
 
    
    tt10.Draw(xaxis+">> e01", denominator(0,500,etamin, etamax, 0, 500, 1, parity_case))
    tt11.Draw(yaxis+">> e02", denominator(0,500,etamin, etamax, 0, 500, 1, parity_case))


    
 
    

    #b1.SetMaximum(50)


    c1.SetLogy()
    Nume = 140.0*207552/250000
    e1 = getIntegral(e01)
    e2 = getIntegral(e02)

    print Nume
    e1.Scale(Nume)
    e2.Scale(Nume)

    e1.Sumw2()
    e2.Sumw2()

    stationxx = ""
    if stx == 1:
            stationxx = "ME11" 
    binxmax = xBins


    legend = TLegend(0.6,0.75,0.9,0.9)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetMargin(0.19)
    
    if parity_case == 0:
        legend.SetHeader(" Odd in %s, Even in ME2 and ME3"%stationxx)
    if parity_case == 1:
        legend.SetHeader(" Odd Chambers in %s, ME2 and ME3"%stationxx)

    if parity_case == 2:
        legend.SetHeader(" Even Chambers in %s, ME2 and ME3"%stationxx)
    if parity_case == 3:
        legend.SetHeader(" Even in %s, Odd in ME2 and ME3"%stationxx)



    even_odd = ""
    if parity_case == 0:
        even_odd = "Odd_Even_Even"
    if parity_case == 1:
        even_odd = "Odd_Odd_Odd"
    if parity_case == 2:
        even_odd = "Even_Even_Even"
    if parity_case == 3:
        even_odd = "Even_Odd_Odd"


    
    e1.SetLineColor(kBlack)
    e1.SetLineWidth(6)

    e2.SetLineColor(kRed)
    e2.SetLineWidth(3)
    #b1.Draw()

    e1.GetXaxis().SetTitle("Reco p_{T} using either method ")
    e1.GetYaxis().SetTitle("Count Scaled up ")
    e1.SetTitle(" << MinBias >> Rate on %.2f < | \eta | < %.2f"%(etamin, etamax))    

    e2.GetXaxis().SetTitle("Reco p_{T} using either method")
    e2.GetYaxis().SetTitle("Count Scaled up ")
    e2.SetTitle(" MinBias Rate on %.2f < | \eta | < %.2f"%(etamin, etamax))    

    e1.Draw("h")
    e2.Draw("same h")

    legend.AddEntry(e1, " Direction Based Only (\Delta \phi)", "l")
    legend.AddEntry(e2, " Position Based Only  (\Delta Y)" , "l")
    

    legend.Draw("same")
    etacut = etamin*10


    
    c1.SaveAs("MinBias_eta%d_%s.pdf"%(etacut, even_odd))
    


sta="ME11"
tree = "HLTBendingAngle/trk_eff_csc_"+sta
ctau = 0
  

stx = 1
xaxis2= "pt_SimTrack"


PTdeltaphi = ""
PTdeltay = ""

    
for parity_case in range (0, 4):
        for etamin in (1.6, 1.8, 2.0, 2.2):
                

                sigma_phi = 0.0
                slope_phi = 0.0
                intercept_phi = 0.0
                
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


                x_bins = "(50,0,50)"
                y_bins= "(50,0, 50)"


                PTdeltaphi = "(( 1/abs(csc_bending_angle_12) + %f )/%f )"%(intercept_phi, slope_phi)
                PTdeltay =" (( 1/abs((delta_y_gp_23) - %f*(delta_y_gp_12) )  + %f )/%f )"%(prop, intercept, slope)

                etamax = etamin + 0.2

                Comparison("MinBias.root",tree, PTdeltaphi, PTdeltay, x_bins, y_bins, etamin, etamax, 0, 500, ctau, stx, parity_case)
                


