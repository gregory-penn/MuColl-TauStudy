from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas, gStyle, TLegend, TLatex
from argparse import ArgumentParser
import os
import math

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('-i', '--inputFile', help='--inputFile input_reco.slcio', 
                  type=str, default='output_reco.slcio')

# Output file
parser.add_argument('-o', '--outputFile', help='--outputFile plot.root', 
                  type=str, default='plot.root')
args = parser.parse_args()

#set global style
gStyle.SetTitleSize(0.042, "X")
gStyle.SetTitleSize(0.042, "Y")
#disabling plot titles - personal preference
gStyle.SetOptTitle(0)
gStyle.SetLineWidth(2)
gStyle.SetOptStat(0) # remove stat box at upper right corner 

#helper functions for calculating theta and phi
def get_theta(px, py, pz):
    pt = math.sqrt(px**2 + py**2)
    return math.atan2(pt, pz)

def get_phi(px, py):
    return math.atan2(py, px)

hists = []

fCls1VsPiEn = TH1F('ratio of cluster 1 and pion energy', 'ClsVsPiEnergy', 14, 0, 1.4)
fCls1VsPiEn.SetXTitle('E^{cluster} / E^{\pi^{-}}')
fCls1VsPiEn.SetMinimum(0.0)
fCls1VsPiEn.SetMaximum(400.0)
pl1Color = 30
fCls1VsPiEn.SetLineColor(pl1Color);
fCls1VsPiEn.SetLineWidth(2)
fCls1VsPiEn.SetMarkerStyle(8)
fCls1VsPiEn.SetMarkerColor(pl1Color)
fCls1VsPiEn.GetYaxis().SetTitle("Entries")

fCls2VsPiEn = TH1F('ratio of cluster 2 and pion energy', 'ClsVsPiEnergy', 14, 0, 1.4)
fCls2VsPiEn.SetXTitle('E^{cluster} / E^{\pi^{-}}')
fCls2VsPiEn.SetMinimum(0.0)
fCls2VsPiEn.SetMaximum(400.0)
pl2Color = 46
fCls2VsPiEn.SetLineColor(pl2Color);
fCls2VsPiEn.SetLineWidth(2)
fCls2VsPiEn.SetMarkerStyle(8)
fCls2VsPiEn.SetMarkerColor(pl2Color)
fCls2VsPiEn.GetYaxis().SetTitle("Entries")

fCls3VsPiEn = TH1F('ratio of cluster 3 and pion energy', 'ClsVsPiEnergy', 14, 0, 1.4)
fCls3VsPiEn.SetXTitle('E^{cluster} / E^{\pi^{-}}')
fCls3VsPiEn.SetMinimum(0.0)
fCls3VsPiEn.SetMaximum(400.0)
pl3Color = 38
fCls3VsPiEn.SetLineColor(pl3Color);
fCls3VsPiEn.SetLineWidth(2)
fCls3VsPiEn.SetMarkerStyle(8)
fCls3VsPiEn.SetMarkerColor(pl3Color)
fCls3VsPiEn.GetYaxis().SetTitle("Entries")

nRecoPion1Cls = TH1F('number of reco pions, events with 1 cluster', 'nPions1Cls', 5, 0, 5)
nRecoPion1Cls.SetXTitle('# of Reco Pions')
nRecoPion1Cls.SetMinimum(0.0)
nRecoPion1Cls.SetMaximum(550.0)
nRecoPion1Cls.SetLineColor(pl1Color);
nRecoPion1Cls.SetLineWidth(2)
nRecoPion1Cls.SetMarkerStyle(8)
nRecoPion1Cls.SetMarkerColor(pl1Color)
nRecoPion1Cls.GetYaxis().SetTitle("Entries")

nRecoPion2Cls = TH1F('number of reco pions, events with 2 clusters', 'nPions2Cls', 5, 0, 5)
nRecoPion2Cls.SetXTitle('# of Reco Pions')
nRecoPion2Cls.SetMinimum(0.0)
nRecoPion2Cls.SetMaximum(550.0)
nRecoPion2Cls.SetLineColor(pl2Color);
nRecoPion2Cls.SetLineWidth(2)
nRecoPion2Cls.SetMarkerStyle(8)
nRecoPion2Cls.SetMarkerColor(pl2Color)
nRecoPion2Cls.GetYaxis().SetTitle("Entries")

nRecoPion3Cls = TH1F('number of reco pions, events with 3 clusters', 'nPions3Cls', 5, 0, 5)
nRecoPion3Cls.SetXTitle('# of Reco Pions')
nRecoPion3Cls.SetMinimum(0.0)
nRecoPion2Cls.SetMaximum(550.0)
nRecoPion3Cls.SetLineColor(pl3Color);
nRecoPion3Cls.SetLineWidth(2)
nRecoPion3Cls.SetMarkerStyle(8)
nRecoPion3Cls.SetMarkerColor(pl3Color)
nRecoPion3Cls.GetYaxis().SetTitle("Entries")

nTrkRecoPion = TH1F('number of tracks per reco pion', 'nTrkPerRecoPion', 5, 0, 5)
nTrkRecoPion.SetXTitle('# of Tracks')
nTrkRecoPion.SetMinimum(0.0)
nTrkRecoPion.SetMaximum(550.0)
nTrkRecoPion.SetLineColor(pl3Color);
nTrkRecoPion.SetLineWidth(2)
nTrkRecoPion.SetMarkerStyle(8)
nTrkRecoPion.SetMarkerColor(pl3Color)
nTrkRecoPion.GetYaxis().SetTitle("Entries")

nPFOs = TH1F('number of PFOs', 'nPFOs', 15, 0, 15)
nPFOs.SetXTitle('# of PFOs')
nPFOs.SetMinimum(0.0)
#nPFOs.SetMaximum(550.0)
nPFOs.SetLineColor(pl1Color);
nPFOs.SetLineWidth(2)
nPFOs.SetMarkerStyle(8)
nPFOs.SetMarkerColor(pl1Color)
nPFOs.GetYaxis().SetTitle("Entries")
hists.append(nPFOs)

nMCPis = TH1F('number of MC Pis', 'nMCPis', 5, 0, 5)
nMCPis.SetXTitle('# of MC Charged Pions')
nMCPis.SetMinimum(0.0)
#nMCPis.SetMaximum(550.0)
nMCPis.SetLineColor(pl1Color);
nMCPis.SetLineWidth(2)
nMCPis.SetMarkerStyle(8)
nMCPis.SetMarkerColor(pl1Color)
nMCPis.GetYaxis().SetTitle("Entries")
hists.append(nMCPis)

nTrks = TH1F('number of tracks', 'nTrks', 5, 0, 5)
nTrks.SetXTitle('# of Tracks')
nTrks.SetMinimum(0.0)
#nTrks.SetMaximum(550.0)
nTrks.SetLineColor(pl1Color);
nTrks.SetLineWidth(2)
nTrks.SetMarkerStyle(8)
nTrks.SetMarkerColor(pl1Color)
nTrks.GetYaxis().SetTitle("Entries")
hists.append(nTrks)

nTrk1PFO = TH1F('number of tracks for 1 PFO', 'nTrk1PFO', 5, 0, 5)
nTrk1PFO.SetXTitle('# of tracks associated to all clusters')
nTrk1PFO.SetMinimum(0.0)
nTrk1PFO.SetMaximum(550.0)
nTrk1PFO.SetLineColor(pl1Color);
nTrk1PFO.SetLineWidth(2)
nTrk1PFO.SetMarkerStyle(8)
nTrk1PFO.SetMarkerColor(pl1Color)
nTrk1PFO.GetYaxis().SetTitle("Entries")
hists.append(nTrk1PFO)

nTrk2PFO = TH1F('number of tracks for 2 PFOs', 'nTrk2PFOs', 5, 0, 5)
nTrk2PFO.SetXTitle('# of tracks associated to all clusters')
nTrk2PFO.SetMinimum(0.0)
nTrk2PFO.SetMaximum(550.0)
nTrk2PFO.SetLineColor(pl2Color);
nTrk2PFO.SetLineWidth(2)
nTrk2PFO.SetMarkerStyle(8)
nTrk2PFO.SetMarkerColor(pl2Color)
nTrk2PFO.GetYaxis().SetTitle("Entries")
hists.append(nTrk2PFO)

nTrk3PFO = TH1F('number of tracks for 3 PFOs', 'nTrk3PFOs', 5, 0, 5)
nTrk3PFO.SetXTitle('# of tracks associated to all clusters')
nTrk3PFO.SetMinimum(0.0)
nTrk3PFO.SetMaximum(550.0)
nTrk3PFO.SetLineColor(pl3Color);
nTrk3PFO.SetLineWidth(2)
nTrk3PFO.SetMarkerStyle(8)
nTrk3PFO.SetMarkerColor(pl3Color)
nTrk3PFO.GetYaxis().SetTitle("Entries")
hists.append(nTrk3PFO)

fSiTrkPtNoFake = TH1F('track p_{T}', 'TrackPtNoFake', 10, 0, 350)
fSiTrkPtNoFake.SetXTitle('track p_{T} (GeV)')
fSiTrkPtNoFake.SetMinimum(0.0)
fSiTrkPtNoFake.GetYaxis().SetTitle("Entries")
hists.append(fSiTrkPtNoFake)

fSiTrkThetaNoFake = TH1F('track theta', 'TrackThetaNoFake', 10, 0, 3.14159)
fSiTrkThetaNoFake.SetXTitle('track #theta (rad)')
fSiTrkThetaNoFake.SetMinimum(0.0)
fSiTrkThetaNoFake.GetYaxis().SetTitle("Entries")
hists.append(fSiTrkThetaNoFake)

fSiTrkPhiNoFake = TH1F('track phi', 'TrackPhiNoFake', 10, -3.14159, 3.14159)
fSiTrkPhiNoFake.SetXTitle('track #phi (rad)')
fSiTrkPhiNoFake.SetMinimum(0.0)
fSiTrkPhiNoFake.GetYaxis().SetTitle("Entries")
hists.append(fSiTrkPhiNoFake)

fMCPiE = TH1F('mc_pi_E', 'MCPionE', 10, 0, 350)
fMCPiE.SetXTitle('Truth Charged Pion E (GeV)')
fMCPiE.SetMinimum(0.0)
fMCPiE.GetYaxis().SetTitle("Entries")

fMCPiPt = TH1F('mc_pi_pt', 'MCPionPt', 10, 0, 350)
fMCPiPt.SetXTitle('Truth Charged Pion p_{T} (GeV)')
fMCPiPt.SetMinimum(0.0)
fMCPiPt.GetYaxis().SetTitle("Entries")
hists.append(fMCPiPt)

fMCPiTheta = TH1F('mc_pi_theta', 'MCPionTheta', 10, 0, 3.14159)
fMCPiTheta.SetXTitle('Truth Charged Pion #theta (rad)')
fMCPiTheta.SetMinimum(0.0)
fMCPiTheta.GetYaxis().SetTitle("Entries")
hists.append(fMCPiTheta)

fMCPiPhi = TH1F('mc_pi_phi', 'MCPionPhi', 10, -3.14159, 3.14159)
fMCPiPhi.SetXTitle('Truth Charged Pion #phi (rad)')
fMCPiPhi.SetMinimum(0.0)
fMCPiPhi.GetYaxis().SetTitle("Entries")
hists.append(fMCPiPhi)

nRecoPFOwTrk = TH1F('PFO_trk', 'PFO_trk', 5, 0, 5)
nRecoPFOwTrk.SetXTitle('Number of PFOs with tracks')
nRecoPFOwTrk.SetMinimum(0.0)
nRecoPFOwTrk.GetYaxis().SetTitle("Entries")
hists.append(nRecoPFOwTrk)

fPFOwTrkPt = TH1F('PFOwTrkPt', 'PFOwTrkPt', 10, 0, 350)
fPFOwTrkPt.SetXTitle('p_{T} of PFOs w/ trks (GeV)')
fPFOwTrkPt.SetMinimum(0.0)
fPFOwTrkPt.GetYaxis().SetTitle("Entries")
hists.append(fPFOwTrkPt)

fPFOwTrkTheta = TH1F('PFOwTrkTheta', 'PFOwTrkTheta', 10, 0, 3.14159)
fPFOwTrkTheta.SetXTitle('#theta of PFOs w/ trks (rad)')
fPFOwTrkTheta.SetMinimum(0.0)
fPFOwTrkTheta.GetYaxis().SetTitle("Entries")
hists.append(fPFOwTrkTheta)

fPFOwTrkPhi = TH1F('PFOwTrkPhi', 'PFOwTrkPhi', 10, -3.14159, 3.14159)
fPFOwTrkPhi.SetXTitle('#phi of PFOs w/ trks (GeV)')
fPFOwTrkPhi.SetMinimum(0.0)
fPFOwTrkPhi.GetYaxis().SetTitle("Entries")
hists.append(fPFOwTrkPhi)

nRecoPerPFOwTrk = TH1F('nPions', 'nPions', 5, 0, 5)
nRecoPerPFOwTrk.SetXTitle('Number of Reco Pions')
nRecoPerPFOwTrk.SetMinimum(0.0)
nRecoPerPFOwTrk.GetYaxis().SetTitle("Entries")
hists.append(nRecoPerPFOwTrk)

fRecoPiPt = TH1F('RecoPiPt', 'RecoPiPt', 10, 0, 350)
fRecoPiPt.SetXTitle('\pi^{-} p_{T} (GeV)')
fRecoPiPt.SetMinimum(0.0)
fRecoPiPt.GetYaxis().SetTitle("Entries")
hists.append(fRecoPiPt)

fRecoPiTheta = TH1F('RecoPiTheta', 'RecoPiTheta', 10, 0, 3.14159)
fRecoPiTheta.SetXTitle('\pi^{-} #theta (rad)')
fRecoPiTheta.SetMinimum(0.0)
fRecoPiTheta.GetYaxis().SetTitle("Entries")
hists.append(fRecoPiTheta)

fRecoPiPhi = TH1F('RecoPiPhi', 'RecoPiPhi', 10, -3.14159, 3.14159)
fRecoPiPhi.SetXTitle('\pi^{-} #phi (rad)')
fRecoPiPhi.SetMinimum(0.0)
fRecoPiPhi.GetYaxis().SetTitle("Entries")
hists.append(fRecoPiPhi)

to_process = []
test1 = 0
totalMCPis = 0
totalTrks = 0

if os.path.isdir(args.inputFile):
  for r, d, f in os.walk(args.inputFile):
    for file in f:
      to_process.append(os.path.join(r, file))
else:
  to_process.append(args.inputFile)

# Open input file(s)
for file in to_process:
  reader = IOIMPL.LCFactory.getInstance().createLCReader()
  reader.open(file)

  # Loop through Events
  for ievt, event in enumerate(reader):
    test1 += 1
    mc_pis = []
    reco_pis = []
    pfo_pi_trks = []
    pfo_w_trks = []
    
    mc_particles = event.getCollection('MCParticle')
    my_cluster = event.getCollection('PandoraClusters')
    pfos = event.getCollection('PandoraPFOs')
    my_tracks = event.getCollection('SiTracks')

    # quick sanity check
    if len(my_cluster) == 0:
      print("There are no clusters in an event!")
      exit()

    nPFOs.Fill(len(pfos))

    # Loop over MC
    for mcp in mc_particles:
        p = mcp.getMomentum()
        px = p[0]
        py = p[1]
        pz = p[2]
        phi = get_phi(px,py)
        theta = get_theta(px,py,pz)
        pt = math.sqrt(px**2 + py**2)

        if abs(mcp.getPDG()) == 211:
          mc_pis.append(mcp)

    # Performing sanity checks early on
    if len(mc_pis) == 0:
       print("Error! There are 0 MC pions in event: ", test1)
       exit()

    totalMCPis += len(mc_pis)

    nMCPis.Fill(len(mc_pis))
    piE = mc_pis[0].getEnergy() 
    p = mc_pis[0].getMomentum()
    px = p[0]
    py = p[1]
    pz = p[2]
    pt = math.sqrt(px**2 + py**2)
    phi = get_phi(px,py)
    theta = get_theta(px,py,pz)

    fMCPiE.Fill(piE)
    fMCPiPt.Fill(pt)
    fMCPiTheta.Fill(theta)
    fMCPiPhi.Fill(phi)

    #requiring from the beginning that there is a track between 17 degrees and 163 for theta (0.3 to 2.85 rad) 
   #  check_centrality = []
   #  for track in my_tracks:
   #    tan_lambda = track.getTanLambda() #"dip angle" in r-z at reference point
   #    theta = (math.pi / 2) - math.atan(tan_lambda)

   #    if theta < 0.35 or theta > 2.75:
   #       check_centrality.append(1)
   #       continue
      
   #  for mc in mc_particles:
   #     px, py, pz = (
   #        mc.getMomentum()[0],
   #        mc.getMomentum()[1],
   #        mc.getMomentum()[2]
   #     )
   #     theta = get_theta(px,py,pz)
   #     if theta < 0.35 or theta > 2.75:
   #        check_centrality.append(1)
   #        continue

   #  test1 += 1 
   #  if len(check_centrality) > 0:
   #     continue
   #  test2 += 1

    if len(pfos) == 3:
       nTrk3PFO.Fill( len(pfos[0].getTracks()) + len(pfos[1].getTracks()) + len(pfos[2].getTracks()) )
    elif len(pfos) == 2:
       nTrk2PFO.Fill( len(pfos[0].getTracks()) + len(pfos[1].getTracks()) )
    elif len(pfos) == 1:
       nTrk1PFO.Fill( len(pfos[0].getTracks()) )                 

    # Loop over PFOs (reconstructed particles)
    for pfo in pfos:
      if len(pfo.getTracks()) == 1:
         pfo_w_trks.append(pfo)
         nRecoPFOwTrk.Fill(1)
         p = pfo.getMomentum()
         px = p[0]
         py = p[1]
         pz = p[2]
         phi = get_phi(px,py)
         theta = get_theta(px,py,pz)
         pt = math.sqrt(px**2 + py**2)
         fPFOwTrkPt.Fill(pt)
         fPFOwTrkPhi.Fill(phi)
         fPFOwTrkTheta.Fill(theta)

         if abs(pfo.getType()) == 211:
            nRecoPerPFOwTrk.Fill(1)
            fRecoPiPt.Fill(pt)
            fRecoPiTheta.Fill(theta)
            fRecoPiPhi.Fill(phi)

    if not pfo_w_trks: #if list is empty
       nRecoPFOwTrk.Fill(0)        
    
    if len(my_cluster) == 1:
        fCls1VsPiEn.Fill(my_cluster[0].getEnergy() / piE)
        nRecoPion1Cls.Fill(len(reco_pis))

    elif len(my_cluster) == 2:
        nRecoPion2Cls.Fill(len(reco_pis))
        for cluster in my_cluster:
            fCls2VsPiEn.Fill(cluster.getEnergy() / piE)

    elif len(my_cluster) == 3:
        nRecoPion3Cls.Fill(len(reco_pis))
        for cluster in my_cluster:
            fCls3VsPiEn.Fill(cluster.getEnergy()/ piE)

    # now I'll remove double counting of tracks when there are fakes. When two tracks are present, I will arbitrarily take the first. 
    # TODO: Maybe there is some chi^2 metric for the track fit which could be compared?  
    # filling histograms for just one track per event. For efficiency plots.

    BFIELD = 3.57 # Taking 5 T for MAIA and 3.57 T for MuColl_v1. I think this is correct?
    FACTOR = 3e-4 #some conversion factor to take T calculation to GeV pT?

    nTrks.Fill(len(my_tracks))

    totalTrks += len(my_tracks)

    if len(my_tracks) != 0:
      omegaNF, tan_lambdaNF, phiNF = (
            my_tracks[0].getOmega(), #signed curvature of track in 1/mm
            my_tracks[0].getTanLambda(), #"dip angle" in r-z at reference point
            my_tracks[0].getPhi(),
            )
      thetaNF = (math.pi / 2) - math.atan(tan_lambdaNF)
      ptNF = BFIELD * FACTOR / abs(omegaNF)
      pzNF = ptNF * tan_lambdaNF
      pNF = math.sqrt(ptNF * ptNF + pzNF * pzNF)

      fSiTrkPtNoFake.Fill(ptNF)
      fSiTrkThetaNoFake.Fill(thetaNF)
      fSiTrkPhiNoFake.Fill(my_tracks[0].getPhi())

  reader.close()

myCanvas = TCanvas("scoobyDoo", "Canvas", 800, 600)
myCanvas.SetFrameLineWidth(2)

fCls1VsPiEn.Draw("E")
fCls2VsPiEn.Draw("E SAME")
fCls3VsPiEn.Draw("E SAME")

legend = TLegend(0.7, 0.7, 0.89, 0.89)
legend.AddEntry(fCls1VsPiEn, "1 cluster", "lep")
legend.AddEntry(fCls2VsPiEn, "2 clusters", "lep")
legend.AddEntry(fCls3VsPiEn, "3 clusters", "lep")
legend.SetBorderSize(0)
legend.Draw()

myCanvas.Update()
myCanvas.SaveAs("e_ratio.png")

myCanvas2 = TCanvas("scoobyDoo2", "Canvas", 800, 600)
myCanvas2.SetFrameLineWidth(2)

nRecoPion1Cls.Draw("E")
nRecoPion2Cls.Draw("E SAME")
nRecoPion3Cls.Draw("E SAME")

legend.Draw()

myCanvas2.Update()
myCanvas2.SaveAs("n_reco_pions.png")

myCanvas3 = TCanvas("scoobyDoo3", "Canvas", 800, 600)
myCanvas3.SetFrameLineWidth(2)

nTrk1PFO.Draw("E")
nTrk2PFO.Draw("E SAME")
nTrk3PFO.Draw("E SAME")

legend.Draw()

myCanvas3.Update()
myCanvas3.SaveAs("n_trk_per_cluster.png")

# Tracking efficiency plots
fTrkPtEff = fSiTrkPtNoFake.Clone('pi_eff')
fTrkPtEff.Divide(fTrkPtEff, fMCPiPt, 1, 1, 'B')
fTrkPtEff.SetLineColor(pl1Color)
fTrkPtEff.SetLineWidth(2)
fTrkPtEff.SetTitle('TrackEfficiencyVsPt')
fTrkPtEff.SetMarkerStyle(8)
fTrkPtEff.SetMarkerColor(pl1Color)
fTrkPtEff.GetXaxis().SetTitle('Truth Charged Pion p_{T} (GeV)')
fTrkPtEff.GetYaxis().SetTitle("Efficiency")
fTrkPtEff.SetMinimum(0.0)
fTrkPtEff.SetMaximum(1.4)

fTrkThetaEff = fSiTrkThetaNoFake.Clone('pi_eff')
fTrkThetaEff.Divide(fTrkThetaEff, fMCPiTheta, 1, 1, 'B')
fTrkThetaEff.SetLineColor(pl1Color)
fTrkThetaEff.SetLineWidth(2)
fTrkThetaEff.SetTitle('TrackEfficiencyVsTheta')
fTrkThetaEff.SetMarkerStyle(8)
fTrkThetaEff.SetMarkerColor(pl1Color)
fTrkThetaEff.GetXaxis().SetTitle('Truth Charged Pion #theta (rad)')
fTrkThetaEff.GetYaxis().SetTitle("Efficiency")
fTrkThetaEff.SetMinimum(0.0)
fTrkThetaEff.SetMaximum(1.4)

fTrkPhiEff = fSiTrkPhiNoFake.Clone('pi_eff')
fTrkPhiEff.Divide(fTrkPhiEff, fMCPiPhi, 1, 1, 'B')
fTrkPhiEff.SetLineColor(pl1Color)
fTrkPhiEff.SetLineWidth(2)
fTrkPhiEff.SetTitle('TrackEfficiencyVsPhi')
fTrkPhiEff.SetMarkerStyle(8)
fTrkPhiEff.SetMarkerColor(pl1Color)
fTrkPhiEff.GetXaxis().SetTitle('Truth Charged Pion #phi (rad)')
fTrkPhiEff.GetYaxis().SetTitle("Efficiency")
fTrkPhiEff.SetMinimum(0.0)
fTrkPhiEff.SetMaximum(1.4)

# Duplicate for individual plotting... want different y-titles. Ugly hack.
fTrkPtEffDupe = fSiTrkPtNoFake.Clone('pi_eff')
fTrkPtEffDupe.Divide(fTrkPtEffDupe, fMCPiPt, 1, 1, 'B')
fTrkPtEffDupe.SetLineColor(pl1Color)
fTrkPtEffDupe.SetLineWidth(2)
fTrkPtEffDupe.SetTitle('TrackEfficiencyVsPt')
fTrkPtEffDupe.SetMarkerStyle(8)
fTrkPtEffDupe.SetMarkerColor(pl1Color)
fTrkPtEffDupe.GetXaxis().SetTitle('Truth Charged Pion p_{T} (GeV)')
fTrkPtEffDupe.GetYaxis().SetTitle("Tracking Efficiency")
fTrkPtEffDupe.SetMinimum(0.0)
fTrkPtEffDupe.SetMaximum(1.4)
hists.append(fTrkPtEffDupe)


# Create pi ID, given PFO with track plot
fPiEffPt = fRecoPiPt.Clone('pi_eff')
fPiEffPt.Divide(fPiEffPt, fPFOwTrkPt, 1, 1, 'B')
fPiEffPt.SetLineColor(pl2Color)
fPiEffPt.SetLineWidth(2)
fPiEffPt.SetTitle('ChargedPionEffVsPt')
fPiEffPt.SetMarkerStyle(8)
fPiEffPt.SetMarkerColor(pl2Color)
fPiEffPt.GetXaxis().SetTitle('PFO (w/ trk) p_{T} (GeV)')
fPiEffPt.GetYaxis().SetTitle("Charged Pion ID Efficiency")
fPiEffPt.SetMinimum(0.0)
fPiEffPt.SetMaximum(1.4)
hists.append(fPiEffPt)

fPiEffTheta = fRecoPiTheta.Clone('pi_eff')
fPiEffTheta.Divide(fPiEffTheta, fPFOwTrkTheta, 1, 1, 'B')
fPiEffTheta.SetLineColor(pl2Color)
fPiEffTheta.SetLineWidth(2)
fPiEffTheta.SetTitle('ChargedPionEffVsTheta')
fPiEffTheta.SetMarkerStyle(8)
fPiEffTheta.SetMarkerColor(pl2Color)
fPiEffTheta.GetXaxis().SetTitle('PFO (w/ trk) #theta (rad)')
fPiEffTheta.GetYaxis().SetTitle("Charged Pion ID Efficiency")
fPiEffTheta.SetMinimum(0.0)
fPiEffTheta.SetMaximum(1.4)
hists.append(fPiEffTheta)

fPiEffPhi = fRecoPiPhi.Clone('pi_eff')
fPiEffPhi.Divide(fPiEffPhi, fPFOwTrkPhi, 1, 1, 'B')
fPiEffPhi.SetLineColor(pl2Color)
fPiEffPhi.SetLineWidth(2)
fPiEffPhi.SetTitle('ChargedPionEffVsPhi')
fPiEffPhi.SetMarkerStyle(8)
fPiEffPhi.SetMarkerColor(pl2Color)
fPiEffPhi.GetXaxis().SetTitle('PFO (w/ trk) #phi (rad)')
fPiEffPhi.GetYaxis().SetTitle("Charged Pion ID Efficiency")
fPiEffPhi.SetMinimum(0.0)
fPiEffPhi.SetMaximum(1.4)
hists.append(fPiEffPhi)

# Efficiency to create a PFO with a track
fPFOTrkEffPt = fPFOwTrkPt.Clone('pi_eff')
fPFOTrkEffPt.Divide(fPFOTrkEffPt, fMCPiPt, 1, 1, 'B')
fPFOTrkEffPt.SetLineColor(pl3Color)
fPFOTrkEffPt.SetLineWidth(2)
fPFOTrkEffPt.SetTitle('PFOwTrkEffVsPt')
fPFOTrkEffPt.SetMarkerStyle(8)
fPFOTrkEffPt.SetMarkerColor(pl3Color)
fPFOTrkEffPt.GetXaxis().SetTitle('MC Charged Pion p_{T} (GeV)')
fPFOTrkEffPt.GetYaxis().SetTitle("Track-Cluster Matching Efficiency")
fPFOTrkEffPt.SetMinimum(0.0)
fPFOTrkEffPt.SetMaximum(1.4)
hists.append(fPFOTrkEffPt)

fPFOTrkEffTheta = fPFOwTrkTheta.Clone('pi_eff')
fPFOTrkEffTheta.Divide(fPFOTrkEffTheta, fMCPiTheta, 1, 1, 'B')
fPFOTrkEffTheta.SetLineColor(pl3Color)
fPFOTrkEffTheta.SetLineWidth(2)
fPFOTrkEffTheta.SetTitle('PFOwTrkEffVsTheta')
fPFOTrkEffTheta.SetMarkerStyle(8)
fPFOTrkEffTheta.SetMarkerColor(pl3Color)
fPFOTrkEffTheta.GetXaxis().SetTitle('MC Charged Pion #theta (rad)')
fPFOTrkEffTheta.GetYaxis().SetTitle("Track-Cluster Matching Efficiency")
fPFOTrkEffTheta.SetMinimum(0.0)
fPFOTrkEffTheta.SetMaximum(1.4)
hists.append(fPFOTrkEffTheta)

fPFOTrkEffPhi = fPFOwTrkPhi.Clone('pi_eff')
fPFOTrkEffPhi.Divide(fPFOTrkEffPhi, fMCPiPhi, 1, 1, 'B')
fPFOTrkEffPhi.SetLineColor(pl3Color)
fPFOTrkEffPhi.SetLineWidth(2)
fPFOTrkEffPhi.SetTitle('PFOwTrkEffVsPhi')
fPFOTrkEffPhi.SetMarkerStyle(8)
fPFOTrkEffPhi.SetMarkerColor(pl3Color)
fPFOTrkEffPt.GetXaxis().SetTitle('MC Charged Pion #phi (rad)')
fPFOTrkEffPhi.GetYaxis().SetTitle("Track-Cluster Matching Efficiency")
fPFOTrkEffPhi.SetMinimum(0.0)
fPFOTrkEffPhi.SetMaximum(1.4)
hists.append(fPFOTrkEffPhi)

geo_label = TLatex()
geo_label.SetNDC()
geo_label.SetTextSize(0.04)

gen_label = TLatex()
gen_label.SetNDC()
gen_label.SetTextSize(0.02)

myCanvas4 = TCanvas("scoobyDoo4", "Canvas", 800, 600)
myCanvas4.SetFrameLineWidth(2)

fTrkPtEff.Draw("E")
fPiEffPt.Draw("E SAME")
fPFOTrkEffPt.Draw("E SAME")

legend2 = TLegend(0.55, 0.7, 0.89, 0.89)
legend2.AddEntry(fTrkPtEff, "Tracking Eff", "lep")
legend2.AddEntry(fPFOTrkEffPt, "Track-Cluster Matching Eff", "lep")
legend2.AddEntry(fPiEffPt, "Charged Pion ID Eff", "lep")
legend2.SetBorderSize(0)
legend2.Draw()

geo_label.DrawLatex(0.135, 0.84, "#it{MuColl_v1}")
gen_label.DrawLatex(0.135, 0.81, "#bf{Charged Pion Gun}")

myCanvas4.Update()
myCanvas4.SaveAs("all_efficiencies_pt.png")

myCanvas5 = TCanvas("scoobyDoo5", "Canvas", 800, 600)
myCanvas5.SetFrameLineWidth(2)

fTrkThetaEff.Draw("E")
fPiEffTheta.Draw("E SAME")
fPFOTrkEffTheta.Draw("E SAME")

legend2.Draw()

geo_label.DrawLatex(0.135, 0.84, "#it{MuColl_v1}")
gen_label.DrawLatex(0.135, 0.81, "#bf{Charged Pion Gun}")

myCanvas5.Update()
myCanvas5.SaveAs("all_efficiencies_theta.png")

myCanvas6 = TCanvas("scoobyDoo6", "Canvas", 800, 600)
myCanvas6.SetFrameLineWidth(2)

fTrkPhiEff.Draw("E")
fPiEffPhi.Draw("E SAME")
fPFOTrkEffPhi.Draw("E SAME")

legend2.Draw()

geo_label.DrawLatex(0.135, 0.84, "#it{MuColl_v1}")
gen_label.DrawLatex(0.135, 0.81, "#bf{Charged Pion Gun}")

myCanvas6.Update()
myCanvas6.SaveAs("all_efficiencies_phi.png")


for hist in hists:
  filename = hist.GetTitle() + '.png'
  canvas = TCanvas()
  hist.Draw("E1")
  canvas.SaveAs(filename)


print("number of total events: ", test1)
print("Number of MC Pis: ", totalMCPis)
print("Number of Tracks: ", totalTrks)

#print("number of total events, after track cleaning: ", test2)
