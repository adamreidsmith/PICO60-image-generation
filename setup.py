import ROOT as rt
import numpy as np
import datetime as dt
from ROOT import TCanvas, TH1D, TH2D

NItoUNIXdelta = dt.date(1970, 1, 1) - dt.date(1904, 1, 1) # 24107 days
NItoUNIXseconds = NItoUNIXdelta.days * 24 * 60 * 60 # 2082844800 seconds


########## SET UP ROOT FILES AND TREES ##########
#mergefile = "/Users/smfallows/Data/recon/30l-16/output/merged_all.root"
'''
mergefile = "/Users/smfallows/Data/recon/30l-16/output/merged_all_all.root"
Aufile = "/Users/smfallows/Data/recon/30l-16/output/AuNess.root"
Qfile = "/Users/smfallows/Data/recon/30l-16/output/GetQ.root"
AP175file = "/Users/smfallows/Data/recon/30l-16/output/AP175.root"
'''
mergefile = "/mnt/dpt/home/pico60_recon/devel/30l-16/merged_all_all.root"
Aufile = "/mnt/dpt/home/pico60_recon/devel/30l-16/AuNess.root"
Qfile = "/mnt/dpt/home/pico60_recon/GetQ.root"
AP175file = "/mnt/dpt/home/pico60_recon/AP175.root"

f = rt.TFile(mergefile, "read")
fAu = rt.TFile(Aufile, "read")
fQ = rt.TFile(Qfile, "read")
fAP175 = rt.TFile(AP175file, "read")
tree = f.Get("T")
Autree = fAu.Get("AuT")
Qtree = fQ.Get("T")
AP175tree = fAP175.Get("T")
tree.AddFriend(Autree)
tree.AddFriend(Qtree)
tree.AddFriend(AP175tree)
########## END SET UP ROOT FILES AND TREES ##########


########## DEFINE BASE SELECTION CUTS ##########
pset17_cut = "(pset==21)"
pset24_cut = "(pset==30)"
pset33_cut = "(pset==30)"
temp17_cut = "(TMath::Abs(ts[2] - 16.05) < 0.5)"
temp24_cut = "(TMath::Abs(ts[2] - 16.05) < 0.5)"
temp33_cut = "(TMath::Abs(ts[2] - 14.05) < 0.5)"
PT17cut = pset17_cut + " && abs(pts[0]-pset)<1 && " + temp17_cut
PT24cut = pset24_cut + " && abs(pts[0]-pset)<1 && " + temp24_cut
PT33cut = pset33_cut + " && abs(pts[0]-pset)<1 && " + temp33_cut

tecut = "(te>25)"

# Data written out correctly, and enforce "good" triggers
# NOT THIS OLD TEXT: (no pressure triggers or time-outs for event search)
write_quality = "(trigger_main==0 && timestamp>0 && pset>0 && pts[0]>0 && ts[2]>0)"
camera_trigger = "(trigger_ctic >= 4)" # NOTE: THIS IS NOT WHAT WE'VE GENERALLY USED, BUT IT'S RIGHT (we think - check DAQ)
# pressure triggers have ctic value of 1
# main DAQ loop (manual, timeout) have ctic value of 2 (cameras are 4, 8, 12, etc.)
# trigger indexing: manual, time-out, auto-relaunch, end run
# trigger_main==0 is a quality cut that means none of the above happened

# Base cut select thermodynamic conditions, good expansion, and good data write-out
base17 = PT17cut + " && " + tecut + " && " + write_quality
base24 = PT24cut + " && " + tecut + " && " + write_quality
base33 = PT33cut + " && " + tecut + " && " + write_quality # "(te>25)"
########## END DEFINE BASE SELECTION CUTS ##########


########## DEFINE FIDUCIAL VOLUME CUTS ##########
# Old fiducial cuts
surf_cut = "(Z<523)"
wall_cut = "( (Dwall>6 && Z>0 && Z<=400) ||" + \
            " (Dwall>6 && Z<=0 && R2<=10000) ||" + \
            " (Dwall>13 && Z<=0 && R2>10000) ||" + \
            " (Dwall>13 && Z>400) )"
Fiducial_cuts = "(" + surf_cut + " && " + wall_cut + ")"


# New fiducial cuts
TazoCuts = "(TazoProcessed!=0 && TazoProcessed!=3 && !(TazoProcessed==2 && cDwall<10 && z_theta<-900 && cZ>0))";

cutFid33 =  "( (cZ<0 && cDwall>5) ||" + \
            " (cZ>0 && cZ<523 && cDwall>10) ||" + \
            " (cZ>0 && cZ<523 && cDwall<10 && cDwall>3 && z_theta < 0.21) )"
cutFid33a = "( (cZ<0 && cDwall>10) ||" + \
            "  (cZ<0 && cDwall<10 && cDwall>5 && z_theta < 0.12 && z_theta > 0) ||" + \
            "  (cZ>=0 && cZ<523 && cDwall>11) ||" + \
            "  (cZ>=0 && cZ<523 && cDwall<11 && cDwall>4 && z_theta < 0.12" + \
                    "&& z_theta_chisq<20 && gFrameDiff<10) )"
cutFiducialNew33 = cutFid33 + " && " + TazoCuts
cutFiducialNew33a = cutFid33a + " && " + TazoCuts

cutFid24 =  "( (cZ<0 && cDwall>5) ||" + \
            " (cZ>0 && cZ<519 && cDwall > 10) ||" + \
            " (cZ>0 && cZ<519 && cDwall<10 && cDwall>3 && z_theta < 0.21) ||" + \
            " (cZ>519 && cZ < 529 && gFrameDiff <10 && cDwall > 6) )"
cutFid24a = "( (cZ<0 && cDwall>5) ||" + \
            "  (cZ>0 && cZ<519 && cDwall>7) ||" + \
            "  (cZ>519 && cZ<529 && gFrameDiff<10 && cDwall>6) ||" + \
            "  (cZ>0 && cZ<519 && cDwall<7 && cDwall>3 && z_theta<0.11" + \
                                    "&& gFrameDiff<10) )"
cutFiducialNew24 = cutFid24 + " && " + TazoCuts
cutFiducialNew24a = cutFid24a + " && " + TazoCuts

cutFid17 =  "( (cZ<0 && cDwall>5) ||" + \
            " (cZ>0 && cZ<531 && cDwall>10) ||" + \
            " (cZ>0 && cZ<531 && cDwall<10 && cDwall>3 && z_theta < 0.21) )"
cutFid17a = "( (cZ<0 && cDwall>5 && gFrameDiff<10) ||" + \
            "  (cZ>0 && cZ<531 && cDwall>10) ||" + \
            "  (cZ>0 && cZ<531 && cDwall<10 && cDwall>3 && z_theta < 0.21) )"
cutFiducialNew17 = cutFid17 + " && " + TazoCuts
cutFiducialNew17a = cutFid17a + " && " + TazoCuts
########## END DEFINE FIDUCIAL VOLUME CUTS ##########


##### Combined base selection cuts with fiducial volume cuts #####
# Original fiducial volume below
bulk17 = base17 + " && " + Fiducial_cuts
bulk24 = base24 + " && " + Fiducial_cuts
bulk33 = base33 + " && " + Fiducial_cuts
# Improved fiducial volume below
bulk17inter = base17 + " && " + cutFiducialNew17
bulk24inter = base24 + " && " + cutFiducialNew24
bulk33inter = base33 + " && " + cutFiducialNew33
bulk17new = base17 + " && " + cutFiducialNew17a
bulk24new = base24 + " && " + cutFiducialNew24a
bulk33new = base33 + " && " + cutFiducialNew33a


########## BEGIN PIEZO QUALITY CUTS ########## (mainly used in singles searches)
piezo3_t0_33 = "(piezo_t0[0] < -0.005 && piezo_t0[0] > -0.045)"
piezo3_t0_24 = "(piezo_t0[0] < -0.005 && piezo_t0[0] > -0.040)"
piezo3_t0_17 = "(piezo_t0[0] < -0.005 && piezo_t0[0] > -0.035)"
# same t0 on both piezo 3 and 7, apparently, so just use piezo 3
#piezo7_t0_33 = "(piezo_t0[2] < -0.05 && piezo_t0[2] > -0.07)"
#piezo7_t0_24 = "(piezo_t0[2] < 0 && piezo_t0[2] > -0.05)"
#piezo7_t0_17 = "(piezo_t0[2] < 0 && piezo_t0[2] > -0.05)"
piezo_t0_33 = piezo3_t0_33# + "&&" + piezo7_t0_33
piezo_t0_24 = piezo3_t0_24# + "&&" + piezo7_t0_24
piezo_t0_17 = piezo3_t0_17# + "&&" + piezo7_t0_17

piezo37_rms_33 = "(piezo_prehit_filtered_rms[0] < 0.013 && piezo_prehit_filtered_rms[2] < 0.011)"
piezo37_rms_24 = "(piezo_prehit_filtered_rms[0] < 0.013 && piezo_prehit_filtered_rms[2] < 0.0105)"
piezo37_rms_17 = "(piezo_prehit_filtered_rms[0] < 0.013 && piezo_prehit_filtered_rms[2] < 0.0105)"

piezo37_band6noise_33 = "(piezo_E[39] < 31500 && piezo_E[41] < 27500)"
piezo37_band6noise_24 = "(piezo_E[39] < 31500 && piezo_E[41] < 23500)"
piezo37_band6noise_17 = "(piezo_E[39] < 31500 && piezo_E[41] < 23500)"

#piezo_quality_orig = "(piezo_t0 < 0 && piezo_t0 > -0.07 && piezo_filtered_rms < 0.013)"
piezo_quality33 = piezo_t0_33 + "&&" + piezo37_rms_33 + "&&" + piezo37_band6noise_33
piezo_quality24 = piezo_t0_24 + "&&" + piezo37_rms_24 + "&&" + piezo37_band6noise_24
piezo_quality17 = piezo_t0_17 + "&&" + piezo37_rms_17 + "&&" + piezo37_band6noise_17
########## END PIEZO QUALITY CUTS ##########


########## DEFINE MULTIPLICITY CUTS ##########
##### Singles cuts: optical, and via dytran #####
singles_base = "(ibub==1 && nbub==1)"
singles_dytranCZT = "(dytranCZT > 0.8 && dytranCZT < 1.2)"
singles33 = singles_base + "&&" + piezo_quality33 + "&&" + singles_dytranCZT #"&& dytranCZ>0.7 && dytranCZ<1.3"
singles24 = singles_base + "&&" + piezo_quality24 + "&&" + singles_dytranCZT #"&& dytranCZ>1.4 && dytranCZ<2.2"
singles17 = singles_base + "&&" + piezo_quality17 + "&&" + singles_dytranCZT #"&& dytranCZ>3. && dytranCZ<4.5"

# Multiplicity cuts from dytran, with permissive (wide) bounds
#  element n is the low edge for a cut selecting multiplicity by dytranCZ, with n+1 as the high edge
#  multiplicity of 1, 2, 3, 4, 5+
#  NOTE: we probably should use only the new dytranCZTmult definition now that it is available
dytranmult33 = [0.5, 1.5, 2.5, 3.5, 4.5]
dytranmult24 = [1, 2.5, 4.75, 6.8, 9]
dytranmult17 = [2.5, 5.5, 9.5, 13, 16.5]
dytranCZTmult = [0, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
########## END DEFINE MULTIPLICITY CUTS ##########


########## DEFINE AP CUTS ##########
APscaling33 = 1.0
APscaling24 = 1.34
APscaling17 = 9.0
AP33 = "(TMath::Log(acoustic_bubnum/%.02f)>-0.69314718 && TMath::Log(acoustic_bubnum/%.02f)<0.4054651081)" % (APscaling33, APscaling33)
AP24 = "(TMath::Log(acoustic_bubnum/%.02f)>-0.69314718 && TMath::Log(acoustic_bubnum/%.02f)<0.4054651081)" % (APscaling24, APscaling24)
AP17 = "(TMath::Log(acoustic_bubnum/%.02f)>-0.69314718 && TMath::Log(acoustic_bubnum/%.02f)<0.4054651081)" % (APscaling17, APscaling17)
#AP33 = "(TMath::Log(acoustic_bubnum)>-0.693 && TMath::Log(acoustic_bubnum)<0.4055)"
#AP24 = "(TMath::Log(acoustic_bubnum/1.34)>-0.693 && TMath::Log(acoustic_bubnum/1.34)<0.4055)"
#AP17 = "(TMath::Log(acoustic_bubnum/9.0)>-0.693 && TMath::Log(acoustic_bubnum/9.0)<0.4055)"
########## END DEFINE AP CUTS ##########


goldencut = "(goldenness == 0)"


# Handy shortcuts from here on out

# Run types
gamma_runs = "(run_type==22 || run_type==23 || run_type==24 || run_type==32 || run_type==41)"
neutron_runs = "(run_type==14 || run_type==15 || run_type==16 || run_type==33)" # or 2, if you want AmBe
bg_runs = "(run_type==0 || run_type==10 || run_type==100)"

# Time- and AP-related aliases
tree.SetAlias("minday", "TMath::Floor(((timestamp-2082844800)%86400)/60)") # minute within the day
tree.SetAlias("UNIXtime", "timestamp-2082844800")
tree.SetAlias("AP33", "TMath::Log(acoustic_bubnum/%.02f)" % APscaling33)
tree.SetAlias("AP24", "TMath::Log(acoustic_bubnum/%.02f)" % APscaling24)
tree.SetAlias("AP17", "TMath::Log(acoustic_bubnum/%.02f)" % APscaling17)
tree.SetAlias("AP18", "TMath::Log(acoustic_bubnum175)")

tree.SetAlias("TE32cal", "ts[2]-0.117")

# Acoustic investigations with
tree.SetAlias("FakeYield","(piezo_E_PosCor[50]+piezo_E_PosCor[53]+piezo_E_PosCor[56]+piezo_E_PosCor[59])/(piezo_E_PosCor[62]+piezo_E_PosCor[65]+piezo_E_PosCor[68]+piezo_E_PosCor[71])")
tree.SetAlias("First4bands","piezo_E_PosCor[50]+piezo_E_PosCor[53]+piezo_E_PosCor[56]+piezo_E_PosCor[59]")
tree.SetAlias("Second4bands", "piezo_E_PosCor[62]+piezo_E_PosCor[65]+piezo_E_PosCor[68]+piezo_E_PosCor[71]")

# Function to set the eastern time zone, with appropriate offset, and configure the input hist's x-axis
def setEST(hist):
        rt.gSystem.Setenv("TZ","EST")
        rt.gStyle.SetTimeOffset(0)
        # Use the below offset if you want to use NI "timestamp" rather than above alias "UNIXtime"
        #rt.gStyle.SetTimeOffset(-1*NItoUNIXseconds)

        hist.GetXaxis().SetTimeDisplay(1)
        hist.GetXaxis().SetTimeFormat("%b-%d %H:%M")

# Create lines that designated the times of PICO-60 chiller events; likely no longer relevant
def getSPlines(hist):
        yaxis = hist.GetYaxis()
        nbinsy = yaxis.GetNbins()
        ymin = yaxis.GetBinLowEdge(1)
        ymax = yaxis.GetBinUpEdge(nbinsy)

        # line148 = rt.TLine(1488846600, ymin, 1488846600, ymax) # ~7:30pm EST Mar 6?
        line148 = rt.TLine(1488844800, ymin, 1488844800, ymax) # ~7pm EST Mar 6?
        line149 = rt.TLine(1488909600, ymin, 1488909600, ymax) # ~1pm EST Mar 7?
        line150 = rt.TLine(1488927240, ymin, 1488927240, ymax) # 5:54pm EST Mar 7

        return [line148, line149, line150]
