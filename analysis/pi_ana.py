from pyLCIO import IOIMPL, EVENT, UTIL
from ROOT import TH1F, TFile, TCanvas, gStyle
import math
from argparse import ArgumentParser
from array import array
import os
import fnmatch

# Command line arguments
parser = ArgumentParser()

# Input file
parser.add_argument('-i-', '--inputFile', help='--inputFile output_reco.slcio', 
                  type=str, default='output_reco.slcio')

# Output file
parser.add_argument('-o', '--outputFile', help='--outputFile pi_ana.root', 
                  type=str, default='pi_ana.root')
args = parser.parse_args()

hists = []

#helper functions for calculating theta and phi
def get_theta(px, py, pz):
    pt = math.sqrt(px**2 + py**2)
    return math.atan2(pt, pz)

def get_phi(px, py):
    return math.atan2(py, px)

#set global style
gStyle.SetTitleSize(0.042, "X")
gStyle.SetTitleSize(0.042, "Y")
#disabling plot titles - personal preference
gStyle.SetOptTitle(0)


# Initialize histograms
fPDG = TH1F('pdg', 'MCParticlePDG', 1250, 0, 2500)
fPDG.SetMinimum(0.0)
fPDG.SetXTitle('PDG ID')
fPDG.GetYaxis().SetTitle("Entries")
hists.append(fPDG)
fType = TH1F('type', 'RecoPFOType', 1250, 0, 2500)
fType.SetMinimum(0.0)
fType.SetXTitle('Reco PFO PDG ID')
fType.GetYaxis().SetTitle("Entries")
hists.append(fType)
fNMCPions = TH1F('n_mc_pions', 'NumberOfMCPions', 10, 0, 10)
fNMCPions.SetMinimum(0.0)
fNMCPions.SetXTitle('# of Truth Pions')
fNMCPions.GetYaxis().SetTitle("Entries")
hists.append(fNMCPions)
fNRecoPions = TH1F('n_reco_pions', 'NumberOfRecoPions', 10, 0, 10)
fNRecoPions.SetMinimum(0.0)
fNRecoPions.SetXTitle('# of Reco Pions')
fNRecoPions.GetYaxis().SetTitle("Entries")
hists.append(fNRecoPions)
fMCPiPt = TH1F('mc_pi_pt', 'MCPionPt', 10, 0, 350)
fMCPiPt.SetXTitle('Truth Charged Pion p_{T} (GeV)')
fMCPiPt.SetMinimum(0.0)
fMCPiPt.GetYaxis().SetTitle("Entries")
hists.append(fMCPiPt)
fMCPiPhi = TH1F('mc_pi_phi', 'MCPionPhi', 10, -3.14159, 3.14159)
fMCPiPhi.SetXTitle('Truth Charged Pion #phi (rad)')
fMCPiPhi.SetMinimum(0.0)
fMCPiPhi.GetYaxis().SetTitle("Entries")
hists.append(fMCPiPhi)
fMCPiTheta = TH1F('mc_pi_theta', 'MCPionTheta', 10, 0, 3.14159)
fMCPiTheta.SetXTitle('Truth Charged Pion #theta (rad)')
fMCPiTheta.SetMinimum(0.0)
fMCPiTheta.GetYaxis().SetTitle("Entries")
hists.append(fMCPiTheta)
fMCPi0Pt = TH1F('mc_pi0_pt', 'MCNeutralPionPt', 10, 0, 350)
fMCPi0Pt.SetXTitle('Neutral Pion p_{T} (GeV)')
fMCPi0Pt.SetMinimum(0.0)
fMCPi0Pt.GetYaxis().SetTitle("Entries")
hists.append(fMCPi0Pt)
fRecoPiPt = TH1F('reco_pi_pt', 'RecoPionPt', 10, 0, 350)
fRecoPiPt.SetXTitle('Reco Charged Pion p_{T} (GeV)')
fRecoPiPt.SetMinimum(0.0)
fRecoPiPt.GetYaxis().SetTitle("Entries")
hists.append(fRecoPiPt)
fRecoPiPhi = TH1F('reco_pi_phi', 'RecoPionPhi', 10, -3.14159, 3.14159)
fRecoPiPhi.SetXTitle('Reco Charged Pion #phi (rad)')
fRecoPiPhi.SetMinimum(0.0)
fRecoPiPhi.GetYaxis().SetTitle("Entries")
hists.append(fRecoPiPhi)
fRecoPiTheta = TH1F('reco_pi_theta', 'RecoPionTheta', 10, 0, 3.14159)
fRecoPiTheta.SetXTitle('Reco Charged Pion #theta (rad)')
fRecoPiTheta.SetMinimum(0.0)
fRecoPiTheta.GetYaxis().SetTitle("Entries")
hists.append(fRecoPiTheta)
fRecoPi0Pt = TH1F('reco_pi0_pt', 'RecoNeutralPion Pt', 10, 0, 350)
fRecoPi0Pt.SetXTitle('Reco Neutral Pion p_{T} (GeV)')
fRecoPi0Pt.SetMinimum(0.0)
fRecoPi0Pt.GetYaxis().SetTitle("Entries")
hists.append(fRecoPi0Pt)
fNSiTrks = TH1F('number of tracks', 'NumberOfTracks', 5, 0, 5)
fNSiTrks.SetXTitle('# of Tracks')
fNSiTrks.SetMinimum(0.0)
fNSiTrks.GetYaxis().SetTitle("Entries")
hists.append(fNSiTrks)
fNSiTrksNoFake = TH1F('number of tracks, fakes rejected', 'NumberOfTracksNoFake', 5, 0, 5)
fNSiTrksNoFake.SetXTitle('# of Tracks, no fake')
fNSiTrksNoFake.SetMinimum(0.0)
fNSiTrksNoFake.GetYaxis().SetTitle("Entries")
hists.append(fNSiTrksNoFake)
fSiTrkPt = TH1F('track p_{T}', 'TrackPt', 10, 0, 350)
fSiTrkPt.SetXTitle('track p_{T} (GeV)')
fSiTrkPt.SetMinimum(0.0)
fSiTrkPt.GetYaxis().SetTitle("Entries")
hists.append(fSiTrkPt)
fSiTrkPtNoFake = TH1F('track p_{T}', 'TrackPtNoFake', 10, 0, 350)
fSiTrkPtNoFake.SetXTitle('track p_{T} (GeV)')
fSiTrkPtNoFake.SetMinimum(0.0)
fSiTrkPtNoFake.GetYaxis().SetTitle("Entries")
hists.append(fSiTrkPtNoFake)
fSiTrkPtLeading = TH1F('Leading track p_{T}', 'LeadingTrackPt', 10, 0, 350)
fSiTrkPtLeading.SetXTitle('Leading Track p_{T} (GeV)')
fSiTrkPtLeading.SetMinimum(0.0)
fSiTrkPtLeading.GetYaxis().SetTitle("Entries")
hists.append(fSiTrkPtLeading)
fSiTrkPtSubleading = TH1F('Subleading track p_{T}', 'SubLeadingTrackPt', 10, 0, 350)
fSiTrkPtSubleading.SetXTitle('track p_{T} (GeV)')
fSiTrkPtSubleading.SetMinimum(0.0)
fSiTrkPtSubleading.GetYaxis().SetTitle("Entries")
hists.append(fSiTrkPtSubleading)
fSiTrkTheta = TH1F('track theta', 'TrackTheta', 10, 0, 3.14159)
fSiTrkTheta.SetXTitle('track #theta (rad)')
fSiTrkTheta.SetMinimum(0.0)
fSiTrkTheta.GetYaxis().SetTitle("Entries")
hists.append(fSiTrkTheta)
fSiTrkThetaNoFake = TH1F('track theta', 'TrackThetaNoFake', 10, 0, 3.14159)
fSiTrkThetaNoFake.SetXTitle('track #theta (rad)')
fSiTrkThetaNoFake.SetMinimum(0.0)
fSiTrkThetaNoFake.GetYaxis().SetTitle("Entries")
hists.append(fSiTrkThetaNoFake)
fSiTrkPhi = TH1F('track phi', 'TrackPhi', 10, -3.14159, 3.14159)
fSiTrkPhi.SetXTitle('track #phi (rad)')
fSiTrkPhi.SetMinimum(0.0)
fSiTrkPhi.GetYaxis().SetTitle("Entries")
hists.append(fSiTrkPhi)
fSiTrkPhiNoFake = TH1F('track phi', 'TrackPhiNoFake', 10, -3.14159, 3.14159)
fSiTrkPhiNoFake.SetXTitle('track #phi (rad)')
fSiTrkPhiNoFake.SetMinimum(0.0)
fSiTrkPhiNoFake.GetYaxis().SetTitle("Entries")
hists.append(fSiTrkPhiNoFake)



for hist in hists:
  hist.SetDirectory(0)

pfoTypes = {}
  
to_process = []

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
    # print('Entered event')
    mc_pis = []
    mc_pi0s = []
    reco_pis = []
    reco_pi0s = []
    si_tracks = []
    
    mc_particles = event.getCollection('MCParticle')

    # Loop through MCParticles
    for mc_particle in mc_particles:
      #print('Entered MCParticle')
      pdg = abs(mc_particle.getPDG())
      parents = mc_particle.getParents()
      if (len(parents) != 0):
        motherID = parents[-1].getPDG()
      else:
        motherID = 0

      # Tag charged and neutral pions
#      if (pdg == 211 and motherID == 15):
#        mc_pis.append(mc_particle)
#      elif (pdg == 111 and motherID == 15):
#        mc_pi0s.append(mc_particle)
      #hacking to remove tau dependency, only reading negatively charged pions
      if (pdg == 211 and motherID == 0):
        mc_pis.append(mc_particle)
      #elif (pdg == 111):
      #  mc_pi0s.append(mc_particle)

      # Fill PDG hist
      fPDG.Fill(pdg)

    # Fill NPions hist
    fNMCPions.Fill(len(mc_pis)+len(mc_pi0s))

    # Loop over charged pions
    for mc_pi in mc_pis:
      p = mc_pi.getMomentum()
      px = p[0]
      py = p[1]
      pz = p[2]
      pt = math.sqrt(px**2 + py**2)
      phi = get_phi(px, py)
      theta = get_theta(px,py,pz)

      # Fill MC charged pion pt hist
      fMCPiPt.Fill(pt)
      fMCPiPhi.Fill(phi)
      fMCPiTheta.Fill(theta)

    # Loop over neutral pions
    for mc_pi0 in mc_pi0s:
      p = mc_pi0.getMomentum()
      px = p[0]
      py = p[1]
      pt = math.sqrt(px**2 + py**2)

      # Fill MC neutral pion pt hist
      fMCPi0Pt.Fill(pt)

    pfos = event.getCollection('PandoraPFOs')

    # Loop over PFOs (reconstructed particles)
    for pfo in pfos:
      # print('Entered PFO')
      _type = abs(pfo.getType())

      if str(_type) in pfoTypes:
        pfoTypes[str(_type)] += 1
      else:
        pfoTypes[str(_type)] = 1
      
      # Tag charged and neutral pions
      if _type == 211:
        reco_pis.append(pfo)
      elif _type == 111:
        reco_pi0s.append(pfo)

      # Fill PFO type hist (reconstructed particle type)  
      fType.Fill(_type)

    # Fill NPions hist  
    fNRecoPions.Fill(len(reco_pis) + len(reco_pi0s))

    # Loop over charged pions
    for reco_pi in reco_pis:
      p = reco_pi.getMomentum()
      px = p[0]
      py = p[1]
      pt = math.sqrt(px**2 + py**2)
      phi = get_phi(px,py)
      theta = get_theta(px,py,pz)
      # Fill reco charged pion pt hist
      fRecoPiPt.Fill(pt)
      fRecoPiPhi.Fill(phi)
      fRecoPiTheta.Fill(theta)

    # Loop over neutral pions
    for reco_pi0 in reco_pi0s:
      p = reco_pi0.getMomentum()
      px = p[0]
      py = p[1]
      pt = math.sqrt(px**2 + py**2)

      # Fill reco neutral pion pt hist
      fRecoPi0Pt.Fill(pt)  

    # retrieving tracking information
    my_tracks = event.getCollection('SiTracks')
    fNSiTrks.Fill(len(my_tracks))
    if len(my_tracks) < 1: 
      fNSiTrksNoFake.Fill(len(my_tracks))
    # only counting events with fake tracks as one track
    else: 
      fNSiTrksNoFake.Fill(1)

    BFIELD = 5 # Taking 5 T for MAIA and 3.57 T for MuColl_v1. I think this is correct?
    FACTOR = 3e-4 #some conversion factor to take T calculation to GeV pT?

    # filling histograms for all tracks in the event
    for track in my_tracks:
      omega, tan_lambda, phi = (
        track.getOmega(), #signed curvature of track in 1/mm
        track.getTanLambda(), #"dip angle" in r-z at reference point
        track.getPhi(),
        )
      theta = (math.pi / 2) - math.atan(tan_lambda)
      pt = BFIELD * FACTOR / abs(omega)
      pz = pt * tan_lambda
      p = math.sqrt(pt * pt + pz * pz)

      fSiTrkPt.Fill(pt)
      fSiTrkTheta.Fill(theta)
      fSiTrkPhi.Fill(track.getPhi())

    # test to see kinematics of the tracks when there are fakes present (>1 track for a charged pion gun)
    if len(my_tracks) == 2: 
      trackL = my_tracks[0]
      trackSL = my_tracks[1]
      omegaL, tan_lambdaL, phiL = (
          trackL.getOmega(), #signed curvature of track in 1/mm
          trackL.getTanLambda(), #"dip angle" in r-z at reference point
          trackL.getPhi(),
          )
      omegaSL, tan_lambdaSL, phiSL = (
          trackSL.getOmega(), #signed curvature of track in 1/mm
          trackSL.getTanLambda(), #"dip angle" in r-z at reference point
          trackSL.getPhi(),
          )
      ptL = BFIELD * FACTOR / abs(omegaL)
      ptSL = BFIELD * FACTOR / abs(omegaSL)
      if ptSL > ptL:
        print("The second track is higher pT than the first! Just a heads up.")

      fSiTrkPtLeading.Fill(ptL)
      fSiTrkPtSubleading.Fill(ptSL)

    # now I'll remove double counting of tracks when there are fakes. When two tracks are present, I will arbitrarily take the first. 
    # TODO: Maybe there is some chi^2 metric for the track fit which could be compared?  
    # filling histograms for just one track per event. For efficiency plots.

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

# Create tracking efficiency plots
fTrkPtEff = fSiTrkPtNoFake.Clone('pi_eff')
fTrkPtEff.Divide(fTrkPtEff, fMCPiPt, 1, 1, 'B')
fTrkPtEff.SetLineColor(6)
fTrkPtEff.SetLineWidth(2)
fTrkPtEff.SetTitle('TrackEfficiencyVsPt')
fTrkPtEff.GetXaxis().SetTitle('Truth Charged Pion p_{T} (GeV)')
fTrkPtEff.GetYaxis().SetTitle("Tracking Efficiency")
fTrkPtEff.SetMinimum(0.0)
fTrkPtEff.SetMaximum(1.4)
hists.append(fTrkPtEff)

# Create tracking efficiency plots
fTrkThetaEff = fSiTrkThetaNoFake.Clone('pi_eff')
fTrkThetaEff.Divide(fTrkThetaEff, fMCPiTheta, 1, 1, 'B')
fTrkThetaEff.SetLineColor(6)
fTrkThetaEff.SetLineWidth(2)
fTrkThetaEff.SetTitle('TrackEfficiencyVsTheta')
fTrkThetaEff.GetXaxis().SetTitle('Truth Charged Pion #theta (rad)')
fTrkThetaEff.GetYaxis().SetTitle("Tracking Efficiency")
fTrkThetaEff.SetMinimum(0.0)
fTrkThetaEff.SetMaximum(1.4)
hists.append(fTrkThetaEff)

# Create tracking efficiency plots
fTrkPhiEff = fSiTrkPhiNoFake.Clone('pi_eff')
fTrkPhiEff.Divide(fTrkPhiEff, fMCPiPhi, 1, 1, 'B')
fTrkPhiEff.SetLineColor(6)
fTrkPhiEff.SetLineWidth(2)
fTrkPhiEff.SetTitle('TrackEfficiencyVsPhi')
fTrkPhiEff.GetXaxis().SetTitle('Truth Charged Pion #phi (rad)')
fTrkPhiEff.GetYaxis().SetTitle("Tracking Efficiency")
fTrkPhiEff.SetMinimum(0.0)
fTrkPhiEff.SetMaximum(1.4)
hists.append(fTrkPhiEff)

# Create charged pion pt efficiency hist
fPiPtEff = fRecoPiPt.Clone('pi_eff')
fPiPtEff.Divide(fPiPtEff, fMCPiPt, 1, 1, 'B')
fPiPtEff.SetLineColor(6)
fPiPtEff.SetLineWidth(2)
fPiPtEff.SetTitle('ChargedPionEfficiencyVsPt')
fPiPtEff.GetXaxis().SetTitle('Truth Charged Pion p_{T} (GeV)')
fPiPtEff.GetYaxis().SetTitle("Reco Efficiency")
fPiPtEff.SetMinimum(0.0)
fPiPtEff.SetMaximum(0.5)
hists.append(fPiPtEff)

# vs theta
fPiThetaEff = fRecoPiTheta.Clone('pi_eff')
fPiThetaEff.Divide(fPiThetaEff, fMCPiTheta, 1, 1, 'B')
fPiThetaEff.SetLineColor(6)
fPiThetaEff.SetLineWidth(2)
fPiThetaEff.SetTitle('ChargedPionEfficiencyVsTheta')
fPiThetaEff.GetXaxis().SetTitle('Truth Charged Pion #theta (rad)')
fPiThetaEff.GetYaxis().SetTitle("Reco Efficiency")
fPiThetaEff.SetMinimum(0.0)
fPiThetaEff.SetMaximum(0.5)
hists.append(fPiThetaEff)

# vs phi
fPiPhiEff = fRecoPiPhi.Clone('pi_eff')
fPiPhiEff.Divide(fPiPhiEff, fMCPiPhi, 1, 1, 'B')
fPiPhiEff.SetLineColor(6)
fPiPhiEff.SetLineWidth(2)
fPiPhiEff.SetTitle('ChargedPionEfficiencyVsPhi')
fPiPhiEff.GetXaxis().SetTitle('Truth Charged Pion #phi (rad)')
fPiPhiEff.GetYaxis().SetTitle("Reco Efficiency")
fPiPhiEff.SetMinimum(0.0)
fPiPhiEff.SetMaximum(0.5)
hists.append(fPiPhiEff)


# Create neutral pion pt efficiency hist
fPi0PtEff = fRecoPi0Pt.Clone('pi0_eff')
fPi0PtEff.Divide(fPi0PtEff, fMCPi0Pt, 1, 1, 'B')
fPi0PtEff.SetLineColor(7)
fPi0PtEff.SetLineWidth(2)
fPi0PtEff.SetTitle('NeutralPionEfficiencyVsPt')
fPi0PtEff.GetXaxis().SetTitle('Truth Neutral Pion p_{T} (GeV)')
fPi0PtEff.GetYaxis().SetTitle('Reco Efficiency')
fPi0PtEff.SetMinimum(0.0)
fPi0PtEff.SetMaximum(0.5)
hists.append(fPi0PtEff)

output_file = TFile(args.outputFile, 'RECREATE')
for hist in hists:
  hist.Write()
output_file.Close()


for hist in hists:
  filename = hist.GetTitle() + '.png'
  canvas = TCanvas()
  hist.Draw()
  canvas.SaveAs(filename)

with open('pfo_types.txt', 'w') as file:
  file.write('PDG/Type of Reconstructed PFOs')
  for key, value in pfoTypes.items():
    print(f'PDG: {key}, Total Number: {value}', file=file)
