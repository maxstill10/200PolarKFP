//*-- Author : Yuri Fisyak 02/02/2016
//LETS CHECK IF IT WILL APPEAR!!!!!!!!!!
#include "StKFParticleAnalysisMaker.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TSystem.h"
//--- KF particle classes ---
#include "KFVertex.h"
#include "KFParticle.h"
#include "KFParticleSIMD.h"
#include "KFPTrack.h"
#include "KFParticleTopoReconstructor.h"
#include "StKFParticleInterface.h"
#include "StKFParticlePerformanceInterface.h"
//--- Pico classes ---
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoEpdHit.h"
//--- Mu classes ---
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
//--- TMVA classes ---
#include "TMVA/GeneticAlgorithm.h"
#include "TMVA/GeneticFitter.h"
#include "TMVA/IFitterTarget.h"
#include "TMVA/Factory.h"
//--- StRefMult class ---
#include "StRefMultCorr/StRefMultCorr.h"
#include "StRefMultCorr/CentralityMaker.h"

#include "TMath.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TLorentzVector.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"



#include <cmath>
ClassImp(StKFParticleAnalysisMaker);

const Char_t *gPartName[] = {"Lambda", "AntiLambda", "Omega", "AntiOmega", "OmegaDaughter", "AntiOmegaDaughter" , "Proton_Lam", "AntiProton_AnLam", "AntiProton_Omega", "Proton_Omega"};
const Char_t *gDetName[] = {"BBC", "EPD", "ZDC"};
const float gMass[] = {1.115683, 1.115683, 1.672, 1.672, 1.115683, 1.115683,0.938, 0.938, 0.938, 0.938};
const float gDMass[] = {0.938, 0.938, 1.115683, 1.115683, 0.938, 0.938, 0.938, 0.938, 1.115683, 1.115683, 0.938, 0.938};
const float gD2Mass[] = {0.1349, 0.1349, 0.493, 0.493, 0.1349, 0.1349};
const float gSigma = 0.007;
const int gPID[] = {3122, -3122, 3334, -3334}, gDPID[] = {2212, -2212, -1, -1}, gD2PID[] = {-211, 211, -321, 321}; 
float gLowMinv[] = {1.1,1.1,1.64,1.64, 1.1, 1.1, 0.2, 0.2, 0.2, 0.2}, gTopMinv[] = {1.13,1.13,1.7,1.7, 1.13, 1.13, 1.1, 1.1, 1.1, 1.1};
int H3Fpoint[DetNum], TPpoint[DetNum], TP2Dpoint[DetNum], TP3Dpoint[DetNum], TP3DCentpoint[DetNum];
int H3Fsize, TPsize, TP2Dsize, TP3Dsize, TP3DCentsize;
int PtBins = 10;
//Big Eta
/*
int EtaBins = 18*2;
float EtaCuts[] = {-3., 3.};
*/
//Small Eta
int EtaBins = 18;
float EtaCuts[] = {-1.5, 1.5};
double pt_intervals[6] = {0.15, 0.75, 1.2, 1.55, 2., 3.};
double eta_intervals[10] = {-1.5, -1.2, -1., -0.6, -0.2, 0.2, 0.6, 1., 1.2, 1.5};

double phiCenter[12][31][2];
double deltaPhi = (30.0/180.0)*TMath::Pi();
StRefMultCorr *refmultCorrUtil;



int mEventsDone = 0;
bool mTrackFound = true;
//________________________________________________________________________________
StKFParticleAnalysisMaker::StKFParticleAnalysisMaker(const char *name) : StMaker(name), fNTrackTMVACuts(0), fIsPicoAnalysis(true), fdEdXMode(1), 
  fStoreTmvaNTuples(false), fProcessSignal(false), fCollectTrackHistograms(false), fCollectPIDHistograms(false),fTMVAselection(false), 
  fFlowAnalysis(false), fFlowChain(NULL), fFlowRunId(-1), fFlowEventId(-1), fCentrality(-1), fFlowFiles(), fFlowMap(), 
  fRunCentralityAnalysis(0), fRefmultCorrUtil(0), fCentralityFile(""), fAnalyseDsPhiPi(false)
{
  memset(mBeg,0,mEnd-mBeg+1);
  
  fNTuplePDG[0] = 421;
  fNTuplePDG[1] = 411;
  fNTuplePDG[2] = 431;
  fNTuplePDG[3] = 4122;
  fNTuplePDG[4] = 426;
  fNTuplePDG[5] = 429;
  fNTuplePDG[6] = 521;
  fNTuplePDG[7] = 511;
  
  fNtupleNames[0] = "D0"; 
  fNtupleNames[1] = "DPlus"; 
  fNtupleNames[2] = "Ds"; 
  fNtupleNames[3] = "Lc";
  fNtupleNames[4] = "D0KK";
  fNtupleNames[5] = "D04";
  fNtupleNames[6] = "BPlus";
  fNtupleNames[7] = "B0";
  
  vector<TString> trackCutNames;
  trackCutNames.push_back("chi2Primary_");
  trackCutNames.push_back("dEdXPi_");
  trackCutNames.push_back("dEdXK_");
  trackCutNames.push_back("dEdXP_");
  trackCutNames.push_back("ToFPi_");
  trackCutNames.push_back("ToFK_");
  trackCutNames.push_back("ToFP_");
  fNTrackTMVACuts = trackCutNames.size();
  
  fDaughterNames[0].push_back("K");     fDaughterNames[0].push_back("Pi");                                                                              //D0 -> Kpi
  fDaughterNames[1].push_back("K");     fDaughterNames[1].push_back("Pi1");    fDaughterNames[1].push_back("Pi2");                                      //D+ -> Kpipi
  fDaughterNames[2].push_back("KPlus"); fDaughterNames[2].push_back("KMinus"); fDaughterNames[2].push_back("Pi");                                       //Ds -> KKpi
  fDaughterNames[3].push_back("K");     fDaughterNames[3].push_back("Pi");     fDaughterNames[3].push_back("P");                                        //Lc -> pKpi
  fDaughterNames[4].push_back("KPlus"); fDaughterNames[4].push_back("KMinus");                                                                          //D0 -> KK
  fDaughterNames[5].push_back("K");     fDaughterNames[5].push_back("Pi1");    fDaughterNames[5].push_back("Pi2");  fDaughterNames[5].push_back("Pi3"); //D0 -> Kpipipi
  fDaughterNames[6].push_back("PiD");   fDaughterNames[6].push_back("KD");     fDaughterNames[6].push_back("Pi");                                       //B+ -> D0_bpi
  fDaughterNames[7].push_back("Pi1D");  fDaughterNames[7].push_back("KD");     fDaughterNames[7].push_back("Pi2D"); fDaughterNames[7].push_back("Pi");  //B0 -> D-pi+

  for(int iDecay=0; iDecay<fNNTuples; iDecay++)
  {
    for(unsigned int iDaughter=0; iDaughter<fDaughterNames[iDecay].size(); iDaughter++)
    {
      for(int iTrackTMVACut=0; iTrackTMVACut<fNTrackTMVACuts; iTrackTMVACut++)
      {
        if(iDaughter==0 && iTrackTMVACut==0)
          fNtupleCutNames[iDecay] = trackCutNames[iTrackTMVACut];  
        else
          fNtupleCutNames[iDecay] += trackCutNames[iTrackTMVACut];
        fNtupleCutNames[iDecay] += fDaughterNames[iDecay][iDaughter];
        fNtupleCutNames[iDecay] += ":";
      }
    }
    if(iDecay<6)
      fNtupleCutNames[iDecay] += "Chi2NDF:LdL:Chi2Topo:refMult";
    else if(iDecay>=6 && iDecay<8)
    {
      fNtupleCutNames[iDecay] += "Chi2NDF_D:LdL_D:Chi2Topo_D:Chi2NDF:LdL:Chi2Topo:refMult";
    } 
    
    SetTMVABins(iDecay);
  }
}
//________________________________________________________________________________
StKFParticleAnalysisMaker::~StKFParticleAnalysisMaker() 
{
  SafeDelete(fStKFParticleInterface);
  SafeDelete(fStKFParticlePerformanceInterface);
}

//_____________________________________________________________________________
Int_t StKFParticleAnalysisMaker::Init()
{
  cout << "STARTING INIT " << endl;
  
  TFile *f = GetTFile();
  if(f) 
  {
  //  f->cd();
    BookVertexPlots();
/*    if(fCollectTrackHistograms)
      fStKFParticleInterface->CollectTrackHistograms();
    if(fCollectPIDHistograms)
      fStKFParticleInterface->CollectPIDHistograms();*/
  }
  f->cd(); 
  
    
  int ew = 0;//east

   for(int pp=1; pp<13; pp++){
        double phiPpCenter = TMath::Pi()/2.0 - (pp-0.5)*deltaPhi;
        if(phiPpCenter<0.0) phiPpCenter += 2.0*TMath::Pi();
        phiCenter[pp-1][0][ew] = phiPpCenter;

        for(int tt=2; tt<32; tt+=2){
                phiCenter[pp-1][tt-1][ew] = phiPpCenter - deltaPhi/4.0;
        }

        for(int tt=3; tt<32; tt+=2){
                phiCenter[pp-1][tt-1][ew] = phiPpCenter + deltaPhi/4.0;
        }
   }

   ew = 1;//west 5.89049


   for(int pp=1; pp<13; pp++){
        double phiPpCenter = TMath::Pi()/2.0 + (pp-0.5)*deltaPhi;
        if(phiPpCenter>2.0*TMath::Pi()) phiPpCenter -= 2.0*TMath::Pi();
        phiCenter[pp-1][0][ew] = phiPpCenter;

        for(int tt=2; tt<32; tt+=2){
                phiCenter[pp-1][tt-1][ew] = phiPpCenter + deltaPhi/4.0;
        }
        for(int tt=3; tt<32; tt+=2){
                phiCenter[pp-1][tt-1][ew] = phiPpCenter - deltaPhi/4.0;
        }

  }

  for(int iSub=0; iSub!=nSub; iSub++){
    QWeight_1[nSub] = 0.;
    QWeight_2[nSub] = 0.;
    Qvec_1[2*nSub] = 0.;
    Qvec_2[2*nSub] = 0.;
    Qvec_3[2*nSub] = 0.;
    Psi1[nSub] = 0.;
    Psi2[nSub] = 0.;
    Psi3[nSub] = 0.;
    deltaPsi1[nSub] = 0.;
    deltaPsi2[nSub] = 0.;
    deltaPsi3[nSub] = 0.;
  }
 

  CreateEPDist();
  CreateKFPHists();
  GetWeightCorr();
  GetCentring();
  GetFlattening();

  isFileRead = true;

  refmultCorrUtil = CentralityMaker::instance()->getRefMultCorr() ;
 
  cout << "HERE WE'VE DONE WITH INIT!!! " << endl;
  if(fTMVAselection || fStoreTmvaNTuples)
  {
    for(int iReader=0; iReader<fNNTuples; iReader++)
    {
      TString cutName;
      int firstSymbolOfCutName = 0;
      
      int nCuts = 0;
      while(fNtupleCutNames[iReader].Tokenize(cutName,firstSymbolOfCutName,":"))
        nCuts++;
      fTMVAParticleParameters[iReader].resize(nCuts);
    }
  }
  
  if(fTMVAselection)
  {
    for(int iReader=0; iReader<fNNTuples; iReader++)
    {
      const int nCentralityBins = fTMVACentralityBins[iReader].size() - 1;
      const int nPtBins = fTMVAPtBins[iReader].size() - 1;
      
      for(int iCentralityBin=0; iCentralityBin<nCentralityBins; iCentralityBin++)
      {
        for(int iPtBin=0; iPtBin<nPtBins; iPtBin++)
        {
          fTMVAReader[iReader][iCentralityBin][iPtBin] = new TMVA::Reader("Silent");

          TString cutName;
          int firstSymbolOfCutName = 0;      
          unsigned int iCut = 0;
          while(fNtupleCutNames[iReader].Tokenize(cutName,firstSymbolOfCutName,":"))
          {
            fTMVAReader[iReader][iCentralityBin][iPtBin] -> AddVariable( cutName.Data(), &fTMVAParticleParameters[iReader][iCut] );
            iCut++;
            if(iCut == (fTMVAParticleParameters[iReader].size()-1)) break;
          }
          
          fTMVAReader[iReader][iCentralityBin][iPtBin] -> BookMVA("BDT", fTMVACutFile[iReader][iCentralityBin][iPtBin].Data());
        }
      }
    }
  }
      
  //Create file with NTuples for cut optimization
  if(fStoreTmvaNTuples)
  {  
    TFile* curFile = gFile;
    TDirectory* curDirectory = gDirectory;
    for(int iNtuple=0; iNtuple<fNNTuples; iNtuple++)
    {
      TString SignalPrefix = "_Signal";
      if(!fProcessSignal) SignalPrefix = "_BG";
      TString currentNTupleFileName = fNtupleNames[iNtuple]+SignalPrefix+TString(".root");
      fNTupleFile[iNtuple] = new TFile(currentNTupleFileName.Data(),"RECREATE");
      fCutsNTuple[iNtuple] = new TNtuple(fNtupleNames[iNtuple].Data(), fNtupleNames[iNtuple].Data(), fNtupleCutNames[iNtuple].Data());
    }
    gFile = curFile;
    gDirectory = curDirectory;
  }
  
//  fRefmultCorrUtil = CentralityMaker::instance()->getgRefMultCorr_P16id();
  fRefmultCorrUtil = CentralityMaker::instance()->getRefMultCorr();
  fRefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
  fRefmultCorrUtil->readScaleForWeight("/gpfs01/star/pwg/pfederic/qVectors/StRoot/StRefMultCorr/macros/weight_grefmult_VpdnoVtx_Vpd5_Run16.txt"); //for new StRefMultCorr, Run16, SL16j
  
  //Initialise the chain with files containing centrality and reaction plane
  if(fFlowAnalysis)
  {
    std::cout << "StKFParticleAnalysisMaker: run flow analysis. Flow file list:"<<std::endl;
    
    fFlowChain = new TChain("mTree");
    for(unsigned int iFlowFile=0; iFlowFile<fFlowFiles.size(); iFlowFile++)
    {
      std::cout << "      " << fFlowFiles[iFlowFile] << std::endl;
      fFlowChain->Add(fFlowFiles[iFlowFile].Data());
    }
    
    fFlowChain->SetBranchStatus("*",0);
    fFlowChain->SetBranchAddress("runid",   &fFlowRunId);   fFlowChain->SetBranchStatus("runid", 1);
    fFlowChain->SetBranchAddress("eventid", &fFlowEventId); fFlowChain->SetBranchStatus("eventid", 1);
    fFlowChain->SetBranchAddress("cent", &fCentrality);  fFlowChain->SetBranchStatus("cent", 1);
    
    std::cout << "StKFParticleAnalysisMaker: number of entries in the flow chain" << fFlowChain->GetEntries() << std::endl;
    for(int iEntry=0; iEntry<fFlowChain->GetEntries(); iEntry++)
    {
      fFlowChain->GetEvent(iEntry);
      fFlowMap[GetUniqueEventId(fFlowRunId, fFlowEventId)] = iEntry;
    }
  }
  return kStOK;
}
//________________________________________________________________________________
Int_t StKFParticleAnalysisMaker::InitRun(Int_t runumber) 
{
//   assert(StPicoDstMaker::instance());
//   if (StPicoDstMaker::instance()->IOMode() == StPicoDstMaker::ioRead) {
    //TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO Ask Yuri
//     StPicoDstMaker::instance()->SetStatus("*",0);
//     const Char_t *ActiveBranches[] = {
//       "MuEvent"
//       ,"PrimaryVertices"
//       ,"PrimaryTracks"
//       ,"GlobalTracks"
//       ,"StStMuMcVertex"
//       ,"StStMuMcTrack"
//       ,"CovPrimTrack"
//       ,"CovGlobTrack"
//       ,"StStMuMcVertex"
//       ,"StStMuMcTrack"
//       ,"KFTracks"
//       ,"KFVertices"
//       ,"StBTofHit"
//       ,"StBTofHeader"
//     }; 
//     Int_t Nb = sizeof(ActiveBranches)/sizeof(Char_t *);
//     for (Int_t i = 0; i < Nb; i++) StPicoDstMaker::instance()->SetStatus(ActiveBranches[i],1); // Set Active braches
//   }
  return StMaker::InitRun(runumber);
}
//_____________________________________________________________________________
void StKFParticleAnalysisMaker::PrintMem(const Char_t *opt)
{
  MemInfo_t info;
  gSystem->GetMemInfo(&info);
  cout << opt 
       << "\tMemory : Total = " << info.fMemTotal 
       << "\tUsed = " << info.fMemUsed
       << "\tFree = " << info.fMemFree
       << "\tSwap Total = " << info.fSwapTotal
       << "\tUsed = " << info.fSwapUsed
       << "\tFree = " << info.fSwapFree << endl;
}
//_____________________________________________________________________________
void StKFParticleAnalysisMaker::BookVertexPlots()
{
  TDirectory *dirs[2] = {0};
  dirs[0] = TDirectory::CurrentDirectory(); assert(dirs[0]);
  dirs[0]->cd();
  if (! dirs[0]->GetDirectory("Particles")) {
    dirs[0]->mkdir("Particles");
  }
  dirs[1] = dirs[0]->GetDirectory("Particles"); assert(dirs[1]);
  dirs[1]->cd();
  PrintMem(dirs[1]->GetPath());
  
  fStKFParticleInterface = new StKFParticleInterface;
  bool storeMCHistograms = false;
  if(!fIsPicoAnalysis && fProcessSignal) storeMCHistograms = true;
  fStKFParticlePerformanceInterface = new StKFParticlePerformanceInterface(fStKFParticleInterface->GetTopoReconstructor(), storeMCHistograms);
//  dirs[0]->cd(); 
//  PrintMem(dirs[1]->GetPath());
//  fStKFParticleInterface = new StKFParticleInterface;
}
//_____________________________________________________________________________
Int_t StKFParticleAnalysisMaker::Make()
{  
  if(fIsPicoAnalysis)
  {
    fPicoDst = StPicoDst::instance();
    if(!fPicoDst) return kStOK;
  }
  else
  {  
    fMuDst = StMuDst::instance();
    if(!fMuDst) return kStOK;
    else { if(StMuDst::instance()->numberOfPrimaryVertices() == 0 ) return kStOK; }
  }
  
  //find max global track index
  int maxGBTrackIndex = -1;
  if(fIsPicoAnalysis)
  {
    for(unsigned int iTrack = 0; iTrack < fPicoDst->numberOfTracks(); iTrack++) 
    {
      StPicoTrack *gTrack = fPicoDst->track(iTrack);
      if (! gTrack) continue;
      int index = gTrack->id();
      if(index > maxGBTrackIndex)
        maxGBTrackIndex = index;
    }
  }
  else
  {
    for(unsigned int iTrack = 0; iTrack < fMuDst->numberOfGlobalTracks(); iTrack++) 
    {
      StMuTrack *gTrack = fMuDst->globalTracks(iTrack);
      if (! gTrack) continue;
      int index = gTrack->id();
      if(index > maxGBTrackIndex)
        maxGBTrackIndex = index;
    }
  }
  vector<KFMCTrack> mcTracks(0);
  vector<int> mcIndices(maxGBTrackIndex+1);
  for(unsigned int iIndex=0; iIndex<mcIndices.size(); iIndex++)
    mcIndices[iIndex] = -1;
  
//   fStKFParticleInterface->SetTriggerMode();
//   fStKFParticleInterface->SetSoftKaonPIDMode();
//   fStKFParticleInterface->SetSoftTofPidMode();
//   fStKFParticleInterface->SetChiPrimaryCut(10);
//   
//   fStKFParticleInterface->SetPtCutCharm(0.5);
//   fStKFParticleInterface->SetChiPrimaryCutCharm(8);
//   fStKFParticleInterface->SetLdLCutCharmManybodyDecays(3);
//   fStKFParticleInterface->SetChi2TopoCutCharmManybodyDecays(10);
//   fStKFParticleInterface->SetChi2CutCharmManybodyDecays(3);
//   fStKFParticleInterface->SetLdLCutCharm2D(3);
//   fStKFParticleInterface->SetChi2TopoCutCharm2D(10);
//   fStKFParticleInterface->SetChi2CutCharm2D(3);
  
  vector<int> triggeredTracks;
  bool isGoodEvent = false;
  
  //Process the event
  if(maxGBTrackIndex > 0)
    fStKFParticleInterface->ResizeTrackPidVectors(maxGBTrackIndex+1);
  if(fIsPicoAnalysis)
    isGoodEvent = fStKFParticleInterface->ProcessEvent(fPicoDst, triggeredTracks);
  else
    isGoodEvent = fStKFParticleInterface->ProcessEvent(fMuDst, mcTracks, mcIndices, fProcessSignal);

//   bool openCharmTrigger = false;
//   if(isGoodEvent) openCharmTrigger =  fStKFParticleInterface->OpenCharmTrigger();
//   fStKFParticleInterface->OpenCharmTriggerCompression(triggeredTracks.size(), fPicoDst->numberOfTracks(), openCharmTrigger);
  //collect histograms
 
  if(isGoodEvent)
  {
    int centralityBin = -1;
    float centralityWeight = 0.;
    
    if(fRunCentralityAnalysis)
    {
      fRefmultCorrUtil->init(fPicoDst->event()->runId());
      if(! (fRefmultCorrUtil->isBadRun(fPicoDst->event()->runId())) )
      {
        fRefmultCorrUtil->initEvent(fPicoDst->event()->grefMult(), fPicoDst->event()->primaryVertex().z(), fPicoDst->event()->ZDCx()) ;
        centralityBin = fRefmultCorrUtil->getCentralityBin9();
        centralityWeight = fRefmultCorrUtil->getWeight();
      }
//       refmultCor = fRefmultCorrUtil->getRefMultCorr();
    }
    
    if(fTMVAselection)
    {
      for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
      {
        KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
              
        for(int iReader=0; iReader<fNNTuples; iReader++)
        {
          if( abs(particle.GetPDG()) == fNTuplePDG[iReader] )
          {
            GetParticleParameters(iReader, particle);
            
            const int iTMVACentralityBin = GetTMVACentralityBin(iReader, centralityBin);
            const int iTMVAPtBin = GetTMVAPtBin(iReader, particle.GetPt());
            
            if(iTMVACentralityBin<0 || iTMVAPtBin<0) 
            {
              fStKFParticleInterface->RemoveParticle(iParticle);
              continue;
            }
           
            if(fTMVAReader[iReader][iTMVACentralityBin][iTMVAPtBin]->EvaluateMVA("BDT") < fTMVACut[iReader][iTMVACentralityBin][iTMVAPtBin])
              fStKFParticleInterface->RemoveParticle(iParticle);
            
            if(fAnalyseDsPhiPi && abs(fStKFParticleInterface->GetParticles()[iParticle].GetPDG()) == 431)
            {              
              KFParticle phi;
              if(particle.GetPDG() == 431)
                phi += fStKFParticleInterface->GetParticles()[particle.DaughterIds()[0]];
              else
                phi += fStKFParticleInterface->GetParticles()[particle.DaughterIds()[1]];
              phi += fStKFParticleInterface->GetParticles()[particle.DaughterIds()[2]];
              float mass = 0.f, dmass = 0.f;
              phi.GetMass(mass, dmass);
              if( fabs(mass - 1.01946) > 0.015)
                fStKFParticleInterface->RemoveParticle(iParticle);
            }
          }
        }
      }      
    }
    
    //clean H3L, H4L, Ln, Lnn
    for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
    {
      KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
      if( abs(particle.GetPDG())==3003 || abs(particle.GetPDG())==3103 || abs(particle.GetPDG())==3004 || abs(particle.GetPDG())==3005)
      {
//         if(particle.GetP() < 1.)
//         {
//           fStKFParticleInterface->RemoveParticle(iParticle);
//           continue;
//         }

//         if(particle.GetPhi() > -0.8 && particle.GetPhi() < -0.4)
//         {
//           fStKFParticleInterface->RemoveParticle(iParticle);
//           continue;
//         }
        
        for(int iD=0; iD<particle.NDaughters(); iD++)
        {
          const int daughterId = particle.DaughterIds()[iD];
          const KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterId];
          if(abs(daughter.GetPDG())==211 && daughter.GetP() > 0.5)
            fStKFParticleInterface->RemoveParticle(iParticle);
        }
      }
    }
    
    int eventId = -1;
    int runId = -1;
    
    if(fFlowAnalysis)
    {
      if(fIsPicoAnalysis) 
      {
        runId   = fPicoDst->event()->runId();
        eventId = fPicoDst->event()->eventId();
      }
      else
      {
        runId   = fMuDst->event()->runId();
        eventId = fMuDst->event()->eventId();
      }
    
      long entryId = GetUniqueEventId(runId, eventId);
      std::map<long,int>::iterator flowMapIterator = fFlowMap.find(entryId);
      if (flowMapIterator != fFlowMap.end())
      {
        fFlowChain->GetEvent(fFlowMap[GetUniqueEventId(runId, eventId)]);
        centralityBin = fCentrality;
      }
    }
    
    centralityWeight = 1;
    
    fStKFParticlePerformanceInterface->SetMCTracks(mcTracks);
    fStKFParticlePerformanceInterface->SetMCIndexes(mcIndices);    
    fStKFParticlePerformanceInterface->SetCentralityBin(centralityBin);
    fStKFParticlePerformanceInterface->SetCentralityWeight(centralityWeight);
    Int_t nevent = 100000;
    fStKFParticlePerformanceInterface->SetPrintEffFrequency(nevent);
    fStKFParticlePerformanceInterface->PerformanceAnalysis();
    
    if(fStoreTmvaNTuples)
    {
      for(int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
      {
        KFParticle particle;
        bool isMCParticle = fStKFParticlePerformanceInterface->GetParticle(particle, iParticle);
              
        if( !( (fProcessSignal && isMCParticle) || (!fProcessSignal && !isMCParticle) ) ) continue;
                  
        for(int iNTuple=0; iNTuple<fNNTuples; iNTuple++)
        {
          if( particle.GetPDG() == fNTuplePDG[iNTuple] )
          {
            GetParticleParameters(iNTuple, particle);
            fCutsNTuple[iNTuple]->Fill(fTMVAParticleParameters[iNTuple].data());
          }
        }
      }
    }
  }

  // Here is my Polarization analysis
  double dphi = 0., Phi_ang = 0., phi_Lam = 0.;
  double w = 0.;
    
  

  
  StPicoEvent *myEvent = fPicoDst->event();
  mEventsDone++;
  //cout << "Done event " << mEventsDone << endl;
  
  if (EventCut(myEvent) == false) return kStOk;   
   
  
  //mTrgEff = 1;
  cout << "Done event " << mEventsDone << endl;
  refmultCorrUtil -> init(myEvent->runId());
  refmultCorrUtil -> initEvent(myEvent->refMult(), myEvent->primaryVertex().z(), myEvent->ZDCx());
  Bool_t isBadRun = refmultCorrUtil->isBadRun(myEvent->runId()); //reject bad runs
  Bool_t isPileUpEvt = !refmultCorrUtil->passnTofMatchRefmultCut(1.*myEvent->refMult(), 1.*myEvent->nBTOFMatch()); //reject pileup events

  int cent = refmultCorrUtil->getCentralityBin9() ;
  if (cent < 0 || isBadRun  || isPileUpEvt) return kStOk;

  for(int iSub=0; iSub!=nSub; iSub++){
    Qvec_1[2*iSub] = 0.;
    Qvec_1[2*iSub+1] = 0.;
    Qvec_2[2*iSub] = 0.;
    Qvec_2[2*iSub+1] = 0.;
    Qvec_3[2*iSub] = 0.;
    Qvec_3[2*iSub+1] = 0.;
    QWeight_1[iSub] = 0.;
    QWeight_2[iSub] = 0.;
    Psi1[iSub] = 0.;
    Psi2[iSub] = 0.;
    Psi3[iSub] = 0.;
  }


  Int_t nPicoTracks = fPicoDst->numberOfTracks();
  
  for (int i = 0; i < nPicoTracks; ++i){
    StPicoTrack* trk = (StPicoTrack*)fPicoDst->track(i);
    if ( !trk ) continue;
    if (trk->nHits() < 15 || trk->pPt() < 0.1 || trk->pPt() > 2. || fabs(trk->pMom().Eta()) > 1.5 || fabs(trk->pMom().Eta())<0.1) continue;
    if ((Float_t)trk->nHits()/trk->nHitsPoss() < 0.52 ) continue;

    double eta = trk->pMom().Eta();

    double wEff = trk->pPt();
    Phi_ang = trk->pMom().Phi(); //femtoTrack->isPrimary() ? femtoTrack->pMom().Phi() : 0.;

    if(eta>0) iSide = 0;//West
    else if(eta<0) iSide = 1;//East

    Qvec_2[2*iSide] += wEff*TMath::Cos(2*Phi_ang);
    Qvec_2[2*iSide+1] += wEff*TMath::Sin(2*Phi_ang);
    Qvec_3[2*iSide] += wEff*TMath::Cos(3*Phi_ang);
    Qvec_3[2*iSide+1] += wEff*TMath::Sin(3*Phi_ang);
    QWeight_2[iSide] += wEff;   

  }//for (int i = 0; i < nPicoTracks; ++i)

  //EPD analisys
  Int_t nEpdHits = fPicoDst->numberOfEpdHits();

  //EPD loop
  for(int iEpd=0; iEpd<nEpdHits; iEpd++){
    //Retrieve ith EPDHit
    StPicoEpdHit* femtoEpdHit = (StPicoEpdHit*)fPicoDst->epdHit(iEpd);

    if(!femtoEpdHit->isGood()) continue;

    PP = femtoEpdHit->position();
    TT = femtoEpdHit->tile();
    EW = femtoEpdHit->side();
    w = femtoEpdHit->nMIP();
    row = femtoEpdHit->row();

    double wEff_EPD = w;
    if(wEff_EPD<0.3) continue;
    if(wEff_EPD>3) wEff_EPD = 3;

    if(EW == 1){

      wEff_EPD *= fabs(v1_average[cent][0]->GetBinContent(row));
      Qvec_1[0] += wEff_EPD*TMath::Cos(phiCenter[PP-1][TT-1][EW]);
      Qvec_1[1] += wEff_EPD*TMath::Sin(phiCenter[PP-1][TT-1][EW]);
      QWeight_1[0] += wEff_EPD;
      
    }

    if(EW==-1){
      wEff_EPD *= fabs(v1_average[cent][1]->GetBinContent(row));
      Qvec_1[2] += wEff_EPD*TMath::Cos(phiCenter[PP-1][TT-1][0]);
      Qvec_1[3] += wEff_EPD*TMath::Sin(phiCenter[PP-1][TT-1][0]);
      QWeight_1[1] += wEff_EPD;
      
    }
  }

    

  //.........................................start of RP calculation.....................................
  bool check = true;
  for(int iSub=0; iSub!=nSub-1; iSub++){
    if(QWeight_1[iSub] == 0) check = false;
    if(QWeight_2[iSub] == 0) check = false;
  }
  if(!check) return kStOk;

  //Get Q vectors
  
  for(int iSub=0; iSub!=nSub-1; iSub++){
    Qvec_1[2*iSub] = Qvec_1[2*iSub]/QWeight_1[iSub];
    Qvec_1[2*iSub+1] = Qvec_1[2*iSub+1]/QWeight_1[iSub];
    if(fabs(Qvec_1[2*iSub])>999 || fabs(Qvec_1[2*iSub+1])>999) check = false;

    Qvec_2[2*iSub] = Qvec_2[2*iSub]/QWeight_2[iSub];
    Qvec_2[2*iSub+1] = Qvec_2[2*iSub+1]/QWeight_2[iSub];
    if(fabs(Qvec_2[2*iSub])>999 || fabs(Qvec_2[2*iSub+1])>999) check = false;

    Qvec_3[2*iSub] = Qvec_3[2*iSub]/QWeight_2[iSub];
    Qvec_3[2*iSub+1] = Qvec_3[2*iSub+1]/QWeight_2[iSub];
    if(fabs(Qvec_3[2*iSub])>999 || fabs(Qvec_3[2*iSub+1])>999) check = false;
  }
  if(!check) return kStOk;
  Qvec_1[4] = Qvec_1[0] - Qvec_1[2];
  Qvec_1[5] = Qvec_1[1] - Qvec_1[3];
  Qvec_2[4] = Qvec_2[0] + Qvec_2[2];
  Qvec_2[5] = Qvec_2[1] + Qvec_2[3];
  Qvec_3[4] = Qvec_3[0] - Qvec_3[2];
  Qvec_3[5] = Qvec_3[1] - Qvec_3[3];

   
  //Recentering
  for(int iSub = 0; iSub!=2*nSub; iSub++){
    Qvec_1[iSub] -= Qvec1Prof_TH[iSub]->GetBinContent(cent+1);
    Qvec_2[iSub] -= Qvec2Prof_TH[iSub]->GetBinContent(cent+1);
    Qvec_3[iSub] -= Qvec3Prof_TH[iSub]->GetBinContent(cent+1);
  }

  
  for(int iSub=0; iSub!=nSub; iSub++){
    Psi1[iSub] = GetPsi(1, Qvec_1[2*iSub], Qvec_1[2*iSub+1]);
    if(Psi1[iSub] > 7) check = false;

    Psi2[iSub] = GetPsi(2, Qvec_2[2*iSub], Qvec_2[2*iSub+1]);
    if(Psi2[iSub] > 3.5) check = false;

    Psi3[iSub] = GetPsi(3, Qvec_3[2*iSub], Qvec_3[2*iSub+1]);
    if(Psi3[iSub] > 2.1) check = false;
  }
  if(!check) return kStOk;
  
  
  //Flattaning
  for(int iProf=0; iProf!=10; iProf++){
    for(int iSub=0; iSub!=nSub; iSub++){
      deltaPsi1[iSub] += Coef_A_n_TH_Psi1[iProf][iSub]->GetBinContent(cent+1)*TMath::Cos((iProf+1)*Psi1[iSub]) + \
                         Coef_B_n_TH_Psi1[iProf][iSub]->GetBinContent(cent+1)*TMath::Sin((iProf+1)*Psi1[iSub]);
      
      deltaPsi2[iSub] += Coef_A_n_TH_Psi2[iProf][iSub]->GetBinContent(cent+1)*TMath::Cos((iProf+1)*2*Psi2[iSub]) + \
                         Coef_B_n_TH_Psi2[iProf][iSub]->GetBinContent(cent+1)*TMath::Sin((iProf+1)*2*Psi2[iSub]);

      deltaPsi3[iSub] += Coef_A_n_TH_Psi3[iProf][iSub]->GetBinContent(cent+1)*TMath::Cos((iProf+1)*3*Psi3[iSub]) + \
                         Coef_B_n_TH_Psi3[iProf][iSub]->GetBinContent(cent+1)*TMath::Sin((iProf+1)*3*Psi3[iSub]);
    }
  }

  for(int iSub=0; iSub!=nSub; iSub++){
    Psi1[iSub] += deltaPsi1[iSub];
    while(Psi1[iSub]>2*TMath::Pi()) Psi1[iSub] -= 2*TMath::Pi();
    while(Psi1[iSub]<0.0) Psi1[iSub] += 2*TMath::Pi();
    deltaPsi1[iSub] = 0.0;
    
    Psi2[iSub] += deltaPsi2[iSub]/2.;
    while(Psi2[iSub]>TMath::Pi()) Psi2[iSub] -= TMath::Pi();
    while(Psi2[iSub]<0.0) Psi2[iSub] += TMath::Pi();
    deltaPsi2[iSub] = 0.0;

    Psi3[iSub] += deltaPsi3[iSub]/3.;
    while(Psi3[iSub]>2*TMath::Pi()/3) Psi3[iSub] -= 2*TMath::Pi()/3;
    while(Psi3[iSub]<0.0) Psi3[iSub] += 2*TMath::Pi()/3;
    deltaPsi3[iSub] = 0.0;
  }

    
  //Fill histograms
  for(int iSub = 0; iSub!=nSub; iSub++){
    Qvec1Hist[cent][2*iSub]->Fill(Qvec_1[2*iSub]);
    Qvec1Hist[cent][2*iSub+1]->Fill(Qvec_1[2*iSub+1]);
    Psi1Hist[cent][iSub]->Fill(Psi1[iSub]);

    Qvec2Hist[cent][2*iSub]->Fill(Qvec_2[2*iSub]);
    Qvec2Hist[cent][2*iSub+1]->Fill(Qvec_2[2*iSub+1]);
    Psi2Hist[cent][iSub]->Fill(Psi2[iSub]);

    Qvec3Hist[cent][2*iSub]->Fill(Qvec_3[2*iSub]);
    Qvec3Hist[cent][2*iSub+1]->Fill(Qvec_3[2*iSub+1]);
    Psi3Hist[cent][iSub]->Fill(Psi3[iSub]);
  }

  
  CosOfDiff_1->Fill(cent, TMath::Cos(Psi1[0] - Psi1[1]));
  CosOfDiff_2->Fill(cent, TMath::Cos(2*(Psi2[0] - Psi2[1])));
  CosOfDiff_3->Fill(cent, TMath::Cos(3*(Psi3[1] - Psi3[0])));

  //.........................................end of RP calculation.....................................
  

  //.........................................Pz calculating............................................
  for (int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++){
    KFParticle particle = fStKFParticleInterface->GetParticles()[iParticle];
    TVector3 ParentVec(particle.GetPx(), particle.GetPy(), particle.GetPz()); 
	

    //Lambda research
    if(particle.GetPDG() == 3122){

      if(ParentVec.Eta()<0) iPsi = 0;
      else iPsi = 1;

                          
      //Get daughters of Lambda hyperon
      for (const auto& elem : particle.DaughterIds()) {

        if(elem<0) continue;
        KFParticle DaugParticle = fStKFParticleInterface->GetParticles()[elem];
        TVector3 DaugVec(DaugParticle.GetPx(), DaugParticle.GetPy(), DaugParticle.GetPz());

        if(abs(DaugParticle.GetPDG())!=2212) continue;
        
        TVector3 Lam_mom = ParentVec * (1/particle.GetE());
        TLorentzVector proton_mom(DaugVec, DaugParticle.GetE());
        proton_mom.Boost(-(Lam_mom));

        //Pz studing
        //if(particle.GetMass() <= 1.111 || particle.GetMass() >= 1.121) continue;
        double inv_m = particle.GetMass();
        phi_Lam = ParentVec.Phi();

        sin_diffPsi1Phi = TMath::Sin(Psi1[2] - proton_mom.Phi());
        cos_diffPsi1Phi = TMath::Cos(proton_mom.Phi() - Psi1[2]);
        
        //Psi2 dependences research
        delta_phi = TMath::Pi()/6;
        //research related Psi_comb
        dphi = phi_Lam-Psi2[2];
        while((dphi) < 0.) dphi+=TMath::Pi();//caus of dphi = (-2pi;pi)
        phi_bin = (int)(dphi/delta_phi);
        if(phi_bin == 6) phi_bin = 5;

        prSin_diffPhiPsi1[cent][phi_bin][1]->Fill(inv_m, sin_diffPsi1Phi);
        prCos_diffPhiPsi1[cent][phi_bin][1]->Fill(inv_m, cos_diffPsi1Phi);
        prCos_theta[cent][phi_bin][1]->Fill(inv_m, proton_mom.CosTheta());
        prSin_theta[cent][phi_bin][1]->Fill(inv_m, TMath::Sin(proton_mom.Theta()));
        prSin_diffPhiPsi1Sin_theta[cent][phi_bin][1]->Fill(inv_m, sin_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));
        prCos_diffPhiPsi1Sin_theta[cent][phi_bin][1]->Fill(inv_m, cos_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));

        InvMLamDist[cent][phi_bin][1]->Fill(inv_m);

        //research related Psi_e/w
        dphi = phi_Lam-Psi2[iPsi];
        while((dphi) < 0.) dphi+=TMath::Pi();
        phi_bin = (int)(dphi/delta_phi);
        if(phi_bin == 6) phi_bin = 5;       

        prSin_diffPhiPsi1[cent][phi_bin][0]->Fill(inv_m, sin_diffPsi1Phi);
        prCos_diffPhiPsi1[cent][phi_bin][0]->Fill(inv_m, cos_diffPsi1Phi);
        prCos_theta[cent][phi_bin][0]->Fill(inv_m, proton_mom.CosTheta());
        prSin_theta[cent][phi_bin][0]->Fill(inv_m, TMath::Sin(proton_mom.Theta()));
        prSin_diffPhiPsi1Sin_theta[cent][phi_bin][0]->Fill(inv_m, sin_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));
        prCos_diffPhiPsi1Sin_theta[cent][phi_bin][0]->Fill(inv_m, cos_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));

        InvMLamDist[cent][phi_bin][0]->Fill(inv_m);

        //Psi3 dependences research
        delta_phi = TMath::Pi()/9;
        //research related Psi_comb
        dphi = phi_Lam-Psi3[2];
        while((dphi) < 0.) dphi+=2*TMath::Pi()/3;//caus of dphi = (-5pi/3;pi)
        if((dphi) > 2*TMath::Pi()/3) dphi-=2*TMath::Pi()/3;
        phi_bin = (int)(dphi/delta_phi);
        if(phi_bin == 6) phi_bin = 5;

        prSin_diffPhiPsi1[cent][phi_bin][3]->Fill(inv_m, sin_diffPsi1Phi);
        prCos_diffPhiPsi1[cent][phi_bin][3]->Fill(inv_m, cos_diffPsi1Phi);
        prCos_theta[cent][phi_bin][3]->Fill(inv_m, proton_mom.CosTheta());
        prSin_theta[cent][phi_bin][3]->Fill(inv_m, TMath::Sin(proton_mom.Theta()));
        prSin_diffPhiPsi1Sin_theta[cent][phi_bin][3]->Fill(inv_m, sin_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));
        prCos_diffPhiPsi1Sin_theta[cent][phi_bin][3]->Fill(inv_m, cos_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));

        InvMLamDist[cent][phi_bin][3]->Fill(inv_m);

        //research related Psi_e/w
        dphi = phi_Lam-Psi3[iPsi];
        while((dphi) < 0.) dphi+=2*TMath::Pi()/3;//caus of dphi = (-5pi/3;pi)
        if((dphi) > 2*TMath::Pi()/3) dphi-=2*TMath::Pi()/3;
        phi_bin = (int)(dphi/delta_phi);
        if(phi_bin == 6) phi_bin = 5;

        prSin_diffPhiPsi1[cent][phi_bin][2]->Fill(inv_m, sin_diffPsi1Phi);
        prCos_diffPhiPsi1[cent][phi_bin][2]->Fill(inv_m, cos_diffPsi1Phi);
        prCos_theta[cent][phi_bin][2]->Fill(inv_m, proton_mom.CosTheta());
        prSin_theta[cent][phi_bin][2]->Fill(inv_m, TMath::Sin(proton_mom.Theta()));
        prSin_diffPhiPsi1Sin_theta[cent][phi_bin][2]->Fill(inv_m, sin_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));
        prCos_diffPhiPsi1Sin_theta[cent][phi_bin][2]->Fill(inv_m, cos_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));

        InvMLamDist[cent][phi_bin][2]->Fill(inv_m);
        

      }//for (const auto& elem : particle.DaughterIds())
    }//if(particle.GetPDG() == 3122)

    //AntiLambda research
    if(particle.GetPDG() == -3122){

      if(ParentVec.Eta()<0) iPsi = 0;
      else iPsi = 1;

      //if(particle.GetMass() <= 1.111 || particle.GetMass() >= 1.121) continue;
                          
      //Get daughters of AntiLambda hyperon
      for (const auto& elem : particle.DaughterIds()) {

        if(elem<0) continue;
        KFParticle DaugParticle = fStKFParticleInterface->GetParticles()[elem];
        TVector3 DaugVec(DaugParticle.GetPx(), DaugParticle.GetPy(), DaugParticle.GetPz());

        if(abs(DaugParticle.GetPDG())!=2212) continue;
        
        
        TVector3 Lam_mom = ParentVec * (1/particle.GetE());
        
        TLorentzVector proton_mom(DaugVec, DaugParticle.GetE());
        proton_mom.Boost(-(Lam_mom));

        //if(particle.GetMass() <= 1.111 || particle.GetMass() >= 1.121) continue;

        sin_diffPsi1Phi = TMath::Sin(proton_mom.Phi() - Psi1[2]);
        cos_diffPsi1Phi = TMath::Cos(proton_mom.Phi() - Psi1[2]);
	
        //Pz studing
        phi_Lam = ParentVec.Phi();
        double inv_m = particle.GetMass();
        
        //Psi2 dependences research
        delta_phi = TMath::Pi()/6;
        dphi = phi_Lam-Psi2[2];
        while((dphi) < 0.) dphi+=TMath::Pi();//caus of dphi = (-2pi;pi)
        phi_bin = (int)(dphi/delta_phi);
        if(phi_bin == 6) phi_bin = 5;

        prSin_diffPhiPsi1_LamBar[cent][phi_bin][1]->Fill(inv_m, sin_diffPsi1Phi);
        prCos_diffPhiPsi1_LamBar[cent][phi_bin][1]->Fill(inv_m, cos_diffPsi1Phi);
        prCos_theta_LamBar[cent][phi_bin][1]->Fill(inv_m, proton_mom.CosTheta());
        prSin_theta_LamBar[cent][phi_bin][1]->Fill(inv_m, TMath::Sin(proton_mom.Theta()));
        prSin_diffPhiPsi1Sin_theta_LamBar[cent][phi_bin][1]->Fill(inv_m, sin_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));
        prCos_diffPhiPsi1Sin_theta_LamBar[cent][phi_bin][1]->Fill(inv_m, cos_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));

        InvMLamBarDist[cent][phi_bin][1]->Fill(inv_m);
        

        dphi = phi_Lam-Psi2[iPsi];
        while((dphi) < 0.) dphi+=TMath::Pi();
        phi_bin = (int)(dphi/delta_phi);
        if(phi_bin == 6) phi_bin = 5;

        prSin_diffPhiPsi1_LamBar[cent][phi_bin][0]->Fill(inv_m, sin_diffPsi1Phi);
        prCos_diffPhiPsi1_LamBar[cent][phi_bin][0]->Fill(inv_m, cos_diffPsi1Phi);
        prCos_theta_LamBar[cent][phi_bin][0]->Fill(inv_m, proton_mom.CosTheta());
        prSin_theta_LamBar[cent][phi_bin][0]->Fill(inv_m, TMath::Sin(proton_mom.Theta()));
        prSin_diffPhiPsi1Sin_theta_LamBar[cent][phi_bin][0]->Fill(inv_m, sin_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));
        prCos_diffPhiPsi1Sin_theta_LamBar[cent][phi_bin][0]->Fill(inv_m, cos_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));

        InvMLamBarDist[cent][phi_bin][0]->Fill(inv_m);

        //Psi3 dependences research
        delta_phi = TMath::Pi()/9;
        //research related Psi_comb
        dphi = phi_Lam-Psi3[2];
        while((dphi) < 0.) dphi+=2*TMath::Pi()/3;//caus of dphi = (-5pi/3;pi)
        if((dphi) > 2*TMath::Pi()/3) dphi-=2*TMath::Pi()/3;
        phi_bin = (int)(dphi/delta_phi);
        if(phi_bin == 6) phi_bin = 5;

        prSin_diffPhiPsi1_LamBar[cent][phi_bin][3]->Fill(inv_m, sin_diffPsi1Phi);
        prCos_diffPhiPsi1_LamBar[cent][phi_bin][3]->Fill(inv_m, cos_diffPsi1Phi);
        prCos_theta_LamBar[cent][phi_bin][3]->Fill(inv_m, proton_mom.CosTheta());
        prSin_theta_LamBar[cent][phi_bin][3]->Fill(inv_m, TMath::Sin(proton_mom.Theta()));
        prSin_diffPhiPsi1Sin_theta_LamBar[cent][phi_bin][3]->Fill(inv_m, sin_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));
        prCos_diffPhiPsi1Sin_theta_LamBar[cent][phi_bin][3]->Fill(inv_m, cos_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));

        InvMLamBarDist[cent][phi_bin][3]->Fill(inv_m);

        //research related Psi_e/w
        dphi = phi_Lam-Psi3[iPsi];
        while((dphi) < 0.) dphi+=2*TMath::Pi()/3;//caus of dphi = (-5pi/3;pi)
        if((dphi) > 2*TMath::Pi()/3) dphi-=2*TMath::Pi()/3;
        phi_bin = (int)(dphi/delta_phi);
        if(phi_bin == 6) phi_bin = 5;

        prSin_diffPhiPsi1_LamBar[cent][phi_bin][2]->Fill(inv_m, sin_diffPsi1Phi);
        prCos_diffPhiPsi1_LamBar[cent][phi_bin][2]->Fill(inv_m, cos_diffPsi1Phi);
        prCos_theta_LamBar[cent][phi_bin][2]->Fill(inv_m, proton_mom.CosTheta());
        prSin_theta_LamBar[cent][phi_bin][2]->Fill(inv_m, TMath::Sin(proton_mom.Theta()));
        prSin_diffPhiPsi1Sin_theta_LamBar[cent][phi_bin][2]->Fill(inv_m, sin_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));
        prCos_diffPhiPsi1Sin_theta_LamBar[cent][phi_bin][2]->Fill(inv_m, cos_diffPsi1Phi*TMath::Sin(proton_mom.Theta()));

        InvMLamBarDist[cent][phi_bin][2]->Fill(inv_m);


      }//for (const auto& elem : particle.DaughterIds())
    }//end of AntiLambda research

  }//for (int iParticle=0; iParticle<fStKFParticlePerformanceInterface->GetNReconstructedParticles(); iParticle++)
  
 
  return kStOk; 
 
}

void StKFParticleAnalysisMaker::GetDaughterParameters(const int iReader, int& iDaughterTrack, int& iDaughterParticle, KFParticle& particle)
{
  if(particle.NDaughters() == 1)
  {
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts]   = particle.GetPt();
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+1] = particle.GetDeviationFromVertex(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
    int trackId = particle.DaughterIds()[0];
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+2]   = fStKFParticleInterface->GetdEdXNSigmaPion(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+3]   = fStKFParticleInterface->GetdEdXNSigmaKaon(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+4]   = fStKFParticleInterface->GetdEdXNSigmaProton(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+5]   = fStKFParticleInterface->GetTofNSigmaPion(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+6]   = fStKFParticleInterface->GetTofNSigmaKaon(trackId);
    fTMVAParticleParameters[iReader][iDaughterTrack*fNTrackTMVACuts+7]   = fStKFParticleInterface->GetTofNSigmaProton(trackId);
    
    iDaughterTrack++;
  }
  else if(particle.NDaughters() > 1)
  {
    int order[4] = {0, 1, 2, 3};
    if( particle.GetPDG() == -421 || particle.GetPDG() == -411 || particle.GetPDG() == -431 ||   
        particle.GetPDG() == -429 || particle.GetPDG() == -4122) 
    { 
      order[0] = 1; 
      order[1] = 0; 
    }
    
    for(int iDaughter=0; iDaughter<particle.NDaughters(); iDaughter++)
    {
      const int daughterParticleIndex = particle.DaughterIds()[order[iDaughter]];
      KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterParticleIndex];
      //set pdg for correct order of cuts
      if(particle.GetPDG() == 521 && daughter.GetPDG() == -1) daughter.SetPDG(-421);
      if(particle.GetPDG() ==-521 && daughter.GetPDG() == -1) daughter.SetPDG( 421);
      if(particle.GetPDG() == 511 && daughter.GetPDG() == -1) daughter.SetPDG(-411);
      if(particle.GetPDG() ==-511 && daughter.GetPDG() == -1) daughter.SetPDG( 411);
        
      GetDaughterParameters(iReader, iDaughterTrack, iDaughterParticle, daughter);
    }
    
    fTMVAParticleParameters[iReader][fDaughterNames[iReader].size()*fNTrackTMVACuts + iDaughterParticle*3] = particle.Chi2()/particle.NDF();  
    
    KFParticleSIMD tempSIMDParticle(particle);
    float_v l,dl;
    KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
    tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
    fTMVAParticleParameters[iReader][fDaughterNames[iReader].size()*fNTrackTMVACuts + iDaughterParticle*3 + 1] = l[0]/dl[0];
    
    tempSIMDParticle.SetProductionVertex(pv);
    fTMVAParticleParameters[iReader][fDaughterNames[iReader].size()*fNTrackTMVACuts + iDaughterParticle*3 + 2] = 
      double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);
    
    iDaughterParticle++;
  }
}

void StKFParticleAnalysisMaker::GetParticleParameters(const int iReader, KFParticle& particle)
{
  bool isBMeson = abs(particle.GetPDG()) == 511 || abs(particle.GetPDG()) == 521;
//   if( !isBMeson ) return;
  
  int iDaughterTrack = 0;
  int iDaughterParticle = 0;
  GetDaughterParameters(iReader, iDaughterTrack, iDaughterParticle, particle);

  int nDaughterParticleCut = 0;
  if(isBMeson) nDaughterParticleCut += 3;
  nDaughterParticleCut += fDaughterNames[iReader].size()*fNTrackTMVACuts;
  
  fTMVAParticleParameters[iReader][nDaughterParticleCut]   = particle.Chi2()/particle.NDF();  
  
  KFParticleSIMD tempSIMDParticle(particle);
  float_v l,dl;
  KFParticleSIMD pv(fStKFParticleInterface->GetTopoReconstructor()->GetPrimVertex());
  tempSIMDParticle.GetDistanceToVertexLine(pv, l, dl);
  fTMVAParticleParameters[iReader][nDaughterParticleCut + 1] = l[0]/dl[0];
  
  tempSIMDParticle.SetProductionVertex(pv);
  fTMVAParticleParameters[iReader][nDaughterParticleCut + 2] = double(tempSIMDParticle.Chi2()[0])/double(tempSIMDParticle.NDF()[0]);

  if(fIsPicoAnalysis)
    fTMVAParticleParameters[iReader][nDaughterParticleCut + 3] = fPicoDst->event()->refMult();
  else
    fTMVAParticleParameters[iReader][nDaughterParticleCut + 3] = fMuDst->event()->refMult();
}

Int_t StKFParticleAnalysisMaker::Finish() 
{
  if(fStoreTmvaNTuples)
  {
    TFile* curFile = gFile;
    TDirectory* curDirectory = gDirectory;
    for(int iNtuple=0; iNtuple<fNNTuples; iNtuple++)
    {
      fNTupleFile[iNtuple]->cd();
      fCutsNTuple[iNtuple]->Write();
    }
    gFile = curFile;
    gDirectory = curDirectory;
  }
 
 
  return kStOK;
}

long StKFParticleAnalysisMaker::GetUniqueEventId(const int iRun, const int iEvent) const
{
  long id = 1000000000;
  return id*(iRun%1000) + iEvent;
}

int StKFParticleAnalysisMaker::GetTMVACentralityBin(int iReader, int centrality)
{
  for(unsigned int iBin=0; iBin<fTMVACentralityBins[iReader].size()-1; iBin++)
    if(centrality >= fTMVACentralityBins[iReader][iBin] && centrality < fTMVACentralityBins[iReader][iBin+1])
      return iBin;
  return -1;
}

int StKFParticleAnalysisMaker::GetTMVAPtBin(int iReader, double pt)
{
  for(unsigned int iBin=0; iBin<fTMVAPtBins[iReader].size()-1; iBin++)
    if(pt >= fTMVAPtBins[iReader][iBin] && pt < fTMVAPtBins[iReader][iBin+1])
      return iBin;
  return -1;
}

void StKFParticleAnalysisMaker::SetTMVACentralityBins(int iReader, TString bins)
{
  fTMVACentralityBins[iReader].clear();
  TString value; int firstSymbol = 0;      
  while(bins.Tokenize(value,firstSymbol,":"))
    fTMVACentralityBins[iReader].push_back(value.Atoi());
}

void StKFParticleAnalysisMaker::SetTMVAPtBins(int iReader, TString bins)
{
  fTMVAPtBins[iReader].clear();
  TString value; int firstSymbol = 0;      
  while(bins.Tokenize(value,firstSymbol,":"))
    fTMVAPtBins[iReader].push_back(value.Atof());
}

void StKFParticleAnalysisMaker::SetTMVABins(int iReader, TString centralityBins, TString ptBins)
{
  SetTMVACentralityBins(iReader, centralityBins);
  SetTMVAPtBins(iReader, ptBins);
  
  const int nCentralityBins = fTMVACentralityBins[iReader].size() - 1;
  const int nPtBins = fTMVAPtBins[iReader].size() - 1;
  
  fTMVACutFile[iReader].resize(nCentralityBins);
  fTMVACut[iReader].resize(nCentralityBins);
  fTMVAReader[iReader].resize(nCentralityBins);
  
  for(int iCentralityBin=0; iCentralityBin<nCentralityBins; iCentralityBin++)
  {
    fTMVACutFile[iReader][iCentralityBin].resize(nPtBins);
    fTMVACut[iReader][iCentralityBin].resize(nPtBins);
    fTMVAReader[iReader][iCentralityBin].resize(nPtBins);
  }
}

int StKFParticleAnalysisMaker::CentralityBin(int Multiplicity) {
        int CentId;
        if (Multiplicity >= 290) CentId = 8; //0-5%
        else if (Multiplicity >= 233) CentId = 7; //5-10%
        else if (Multiplicity >= 150) CentId = 6; //10-20%
        else if (Multiplicity >= 94) CentId = 5; //20-30%
        else if (Multiplicity >= 57) CentId = 4; //30-40%
        else if (Multiplicity >= 32) CentId = 3; //40-50%
        else if (Multiplicity >= 17) CentId = 2; //50-60%
        else if (Multiplicity >= 9) CentId = 1; //60-70%
        else if (Multiplicity >= 4) CentId = 0; //70-80%
        else CentId = -1;
        return CentId;
}
/*
int StKFParticleAnalysisMaker::VzBin(float Vz) { //19.6
  if (Vz < -60) return 0;
  else if (Vz < -50) return 1;
  else if (Vz < -40) return 2;
  else if (Vz < -30) return 3;
  else if (Vz < -20) return 4;
  else if (Vz < -10) return 5;
  else if (Vz < 0) return 6;
  else if (Vz < 10) return 7;
  else if (Vz < 20) return 8;
  else if (Vz < 30) return 9;
  else if (Vz < 40) return 10;
  else if (Vz < 50) return 11;
  else if (Vz < 60) return 12;
  else if (Vz < 70) return 13;
  else return 100;
}*/
int StKFParticleAnalysisMaker::VzBin(float Vz) {
  if (Vz < -70) return 0;
  else if (Vz < 0) return 1;
  else if (Vz < 70) return 2;
  else return 3;
}

bool StKFParticleAnalysisMaker::GoodRun(StPicoEvent* event) { //19.6 Preliminary
        int Run = event->runId();
        int badruns[] = {22145028, 22145045, 22145046, 22146001, 22146002, 22146004, 22146007, 22146008,
                         22146009, 22146010, 22146012, 22149004, 22150023, 22155004, 22155007, 22155044};
        bool IsGood = true;
        for (int iR = 0; iR != 16; iR++) {
                if (Run == badruns[iR]) {IsGood = false; cout << Run << endl; return IsGood;}
        }
        return IsGood;
}

bool StKFParticleAnalysisMaker::EventCut(StPicoEvent *event)
{
  bool cut = true;
  if (!event->isTrigger(650000) && !event->isTrigger(650001) && !event->isTrigger(650002) && !event->isTrigger(650003) &&
      !event->isTrigger(650007) && !event->isTrigger(650004) && !event->isTrigger(650005) && !event->isTrigger(650006) &&
      !event->isTrigger(650009) ) cut = false;
  if (event->primaryVertex().Z() < -145. || event->primaryVertex().Z() > 145. || ( pow(event->primaryVertex().X(),2) + pow(event->primaryVertex().Y(),2) > 4)) cut = false;
//  if (!GoodRun(event))  cut = false;
  return cut;



}


void StKFParticleAnalysisMaker::SetMyStuff(Char_t *outFileName, Char_t *runFileName) {
foutFile = outFileName; frunFile = runFileName;	
}

bool StKFParticleAnalysisMaker::PrithwishRun(StPicoEvent *event) {
        int Run = event->runId();
        const int badrunlist[160] = {
19130078, 19130079, 19130085, 19130086, 19131001, 19131003, 19131006, 19131007, 19131009, 19131010, 19131012,
19131014, 19131015, 19131016, 19131019, 19131020, 19131039, 19131045, 19131049, 19131050, 19131052, 19131057,
19131062, 19132016, 19132029, 19132031, 19132035, 19132046, 19132047, 19132063, 19132083, 19133009, 19133012,
19133013, 19133014, 19133018, 19133021, 19133023, 19133032, 19133033, 19133040, 19133041, 19133050, 19133061,
19134005, 19134008, 19134010, 19134013, 19134017, 19134019, 19134025, 19134045, 19134046, 19134047, 19135012,
19135013, 19135014, 19135022, 19135029, 19135039, 19136001, 19136042, 19137001, 19137003, 19137008, 19137009,
19137010, 19137011, 19137013, 19137015, 19137022, 19137027, 19137050, 19137051, 19137052, 19137056, 19138004,
19138008, 19138010, 19138014, 19138025, 19139023, 19139026, 19139027, 19139028, 19139032, 19139033, 19139037,
19139038, 19139041, 19140014, 19140043, 19141004, 19141008, 19141018, 19141019, 19143008, 19143009, 19143010,
19143011, 19143012, 19143013, 19143014, 19143015, 19143016, 19143017, 19144012, 19144013, 19144014, 19144018,
19144019, 19144020, 19144024, 19144025, 19144026, 19144031, 19144032, 19144033, 19144036, 19144042, 19144044,
19144046, 19144047, 19145001, 19145004, 19145005, 19145006, 19145008, 19145009, 19145010, 19145011, 19145013,
19145014, 19145015, 19145017, 19145019, 19145020, 19145031, 19145035, 19145036, 19145040, 19145042, 19145044,
19145050, 19146002, 19146003, 19146007, 19146008, 19146009, 19146012, 19146013, 19146014, 19146016, 19146017,
19146019, 19146020, 19146024, 19146025, 19146026, 19147007};
        bool IsGood = true;
        for (int iR = 0; iR < 19; iR++) {
                if (Run == badrunlist[iR]) {IsGood = false; return IsGood;}
        }
        return IsGood;
}

int StKFParticleAnalysisMaker::JoeyCentrality(int Multiplicity) {
        int CentId;
        if (Multiplicity >= 299) CentId = 8; // 0-5%
        else if (Multiplicity >= 256) CentId = 7; // 5-10%
        else if (Multiplicity >= 194) CentId = 6; // 10-20%
        else if (Multiplicity >= 149) CentId = 5; // 20-30%
        else if (Multiplicity >= 113) CentId = 4; // 30-40%
        else if (Multiplicity >= 85) CentId = 3; // 40-50%
        else if (Multiplicity >= 61) CentId = 2; // 50-60%
        else if (Multiplicity >= 42) CentId = 1; // 60-70%
        else if (Multiplicity >= 27) CentId = 0; // 70-80
        else  CentId = -1;
        return CentId;
}

bool StKFParticleAnalysisMaker::IsKfGoof(KFParticle particle) {
	bool IsGood = true;
/*	for (int iD = 0; iD < particle.NDaughters(); iD++) {
        	const int daughterId = particle.DaughterIds()[iD];
                const KFParticle daughter = fStKFParticleInterface->GetParticles()[daughterId];
                if (TVector3(daughter.Px(),daughter.Py(),daughter.Pz()).Phi() > -1.2 &&
                    TVector3(daughter.Px(),daughter.Py(),daughter.Pz()).Phi() < 0.6) IsGood = false;}
	if (TVector3(particle.GetPx(),particle.GetPy(),particle.GetPz()).Mag() < 0.5) IsGood = false;*/
	if (TVector3(particle.GetPx(),particle.GetPy(),particle.GetPz()).Mag() < 0.2 || TVector3(particle.GetPx(),particle.GetPy(),particle.GetPz()).Perp() > 5.) IsGood = false;
	return IsGood;
/*	StPicoTrack *TmpTrk = fPicoDst->track(FindTrack(fStKFParticleInterface->GetParticles()[particle.DaughterIds()[0]].DaughterIds()[0], fPicoDst));
	if (!TmpTrk) return false;
	else {
        if (TmpTrk->gMom().Phi() > -0.8 && TmpTrk->gMom().Phi() < 0.2) return false;}
        StPicoTrack *TmpTrkTwo = fPicoDst->track(FindTrack(fStKFParticleInterface->GetParticles()[particle.DaughterIds()[1]].DaughterIds()[0], fPicoDst));
	if (!TmpTrkTwo) return false;
	else {
        if (TmpTrkTwo->gMom().Phi() > -0.8 && TmpTrkTwo->gMom().Phi() < 0.2) return false;}*/
//	if (particle.GetPhi() > -0.8 && particle.GetPhi() < 0.2) return false;
//	else return true;
//	return IsGood;
}

void StKFParticleAnalysisMaker::CreateEPDist() {
	for(int iHist=0; iHist!=9; iHist++){
    for(int iSub=0; iSub!=nSub; iSub++){
      Qvec1Hist[iHist][2*iSub] = new TH1F(Form("Qvec1Hist_%i_%i", iHist, 2*iSub), "", 500, -1, 1);
      Qvec1Hist[iHist][2*iSub+1] = new TH1F(Form("Qvec1Hist_%i_%i", iHist, 2*iSub+1), "", 500, -1, 1);
      Qvec2Hist[iHist][2*iSub] = new TH1F(Form("Qvec2Hist_%i_%i", iHist, 2*iSub), "", 500, -1, 1);
      Qvec2Hist[iHist][2*iSub+1] = new TH1F(Form("Qvec2Hist_%i_%i", iHist, 2*iSub+1), "", 500, -1, 1);
      Qvec3Hist[iHist][2*iSub] = new TH1F(Form("Qvec3Hist_%i_%i", iHist, 2*iSub), "", 500, -1, 1);
      Qvec3Hist[iHist][2*iSub+1] = new TH1F(Form("Qvec3Hist_%i_%i", iHist, 2*iSub+1), "", 500, -1, 1);
      
      Psi1Hist[iHist][iSub] = new TH1F(Form("Psi1Hist_%i_%i", iHist, iSub), "", 140, 0, 7);
      Psi2Hist[iHist][iSub] = new TH1F(Form("Psi2Hist_%i_%i", iHist, iSub), "", 70, 0, 3.5);
      Psi3Hist[iHist][iSub] = new TH1F(Form("Psi3Hist_%i_%i", iHist, iSub), "", 56, 0, 2.1);
    }
  }


  CosOfDiff_1 = new TProfile("CosOfDiff_1", "CosOfDiff for Psi_1", 9, 0, 9);
  CosOfDiff_2 = new TProfile("CosOfDiff_2", "CosOfDiff for Psi_2", 9, 0, 9);
  CosOfDiff_3 = new TProfile("CosOfDiff_3", "CosOfDiff for Psi_3", 9, 0, 9);

  for(int iCent=0; iCent!=9; iCent++){
    Qvec1Hist[iCent][0]->SetTitle(Form("Qx West EPD first harm (%i)", iCent));
    Qvec1Hist[iCent][1]->SetTitle(Form("Qy West EPD first harm (%i)", iCent));
    Qvec1Hist[iCent][2]->SetTitle(Form("Qx East EPD first harm (%i)", iCent));
    Qvec1Hist[iCent][3]->SetTitle(Form("Qy East EPD first harm (%i)", iCent));
    Qvec1Hist[iCent][4]->SetTitle(Form("Qx Comb EPD first harm (%i)", iCent));
    Qvec1Hist[iCent][5]->SetTitle(Form("Qy Comb EPD first harm (%i)", iCent));

    Qvec2Hist[iCent][0]->SetTitle(Form("Qx West TPC second harm (%i)", iCent));
    Qvec2Hist[iCent][1]->SetTitle(Form("Qy West TPC second harm (%i)", iCent));
    Qvec2Hist[iCent][2]->SetTitle(Form("Qx East TPC second harm (%i)", iCent));
    Qvec2Hist[iCent][3]->SetTitle(Form("Qy East TPC second harm (%i)", iCent));
    Qvec2Hist[iCent][4]->SetTitle(Form("Qx Comb TPC second harm (%i)", iCent));
    Qvec2Hist[iCent][5]->SetTitle(Form("Qy Comb TPC second harm (%i)", iCent));

    Qvec3Hist[iCent][0]->SetTitle(Form("Qx West TPC third harm (%i)", iCent));
    Qvec3Hist[iCent][1]->SetTitle(Form("Qy West TPC third harm (%i)", iCent));
    Qvec3Hist[iCent][2]->SetTitle(Form("Qx East TPC third harm (%i)", iCent));
    Qvec3Hist[iCent][3]->SetTitle(Form("Qy East TPC third harm (%i)", iCent));
    Qvec3Hist[iCent][4]->SetTitle(Form("Qx Comb TPC third harm (%i)", iCent));
    Qvec3Hist[iCent][5]->SetTitle(Form("Qy Comb TPC third harm (%i)", iCent));

    Psi1Hist[iCent][0]->SetTitle(Form("West Psi_1 EPD (%i)", iCent));
    Psi1Hist[iCent][1]->SetTitle(Form("East Psi_1 EPD (%i)", iCent));
    Psi1Hist[iCent][2]->SetTitle(Form("Comb Psi_1 EPD (%i)", iCent));
    
    Psi2Hist[iCent][0]->SetTitle(Form("West Psi_2 TPC (%i)", iCent));
    Psi2Hist[iCent][1]->SetTitle(Form("East Psi_2 TPC (%i)", iCent));
    Psi2Hist[iCent][2]->SetTitle(Form("Comb Psi_2 TPC (%i)", iCent));

    Psi3Hist[iCent][0]->SetTitle(Form("West Psi_3 TPC (%i)", iCent));
    Psi3Hist[iCent][1]->SetTitle(Form("East Psi_3 TPC (%i)", iCent));
    Psi3Hist[iCent][2]->SetTitle(Form("Comb Psi_3 TPC (%i)", iCent));

  }

}


void StKFParticleAnalysisMaker::CreateKFPHists() {
  for(int iCent=0; iCent!=9; iCent++){
    for(int iphi=0; iphi!=6; iphi++){
      for(int iSub=0; iSub!=4; iSub++){
        prSin_diffPhiPsi1[iCent][iphi][iSub] = new TProfile(Form("prSin_diffPhiPsi1_%i_%i_%i", iCent, iphi, iSub), "", 30, 1.1, 1.13);
        prCos_diffPhiPsi1[iCent][iphi][iSub] = new TProfile(Form("prCos_diffPhiPsi1_%i_%i_%i", iCent, iphi, iSub), "", 30, 1.1, 1.13);
        prCos_theta[iCent][iphi][iSub] = new TProfile(Form("prCos_theta_%i_%i_%i", iCent, iphi, iSub), "", 30, 1.1, 1.13);
        prSin_theta[iCent][iphi][iSub] = new TProfile(Form("prSin_theta_%i_%i_%i", iCent, iphi, iSub), "", 30, 1.1, 1.13);
        prSin_diffPhiPsi1Sin_theta[iCent][iphi][iSub] = new TProfile(Form("prSin_diffPhiPsi1Sin_theta_%i_%i_%i", iCent, iphi, iSub), "", 30, 1.1, 1.13);
        prCos_diffPhiPsi1Sin_theta[iCent][iphi][iSub] = new TProfile(Form("prCos_diffPhiPsi1Sin_theta_%i_%i_%i", iCent, iphi, iSub), "", 30, 1.1, 1.13);

        prSin_diffPhiPsi1_LamBar[iCent][iphi][iSub] = new TProfile(Form("prSin_diffPhiPsi1_LamBar_%i_%i_%i", iCent, iphi, iSub), "", 30, 1.1, 1.13);
        prCos_diffPhiPsi1_LamBar[iCent][iphi][iSub] = new TProfile(Form("prCos_diffPhiPsi1_LamBar_%i_%i_%i", iCent, iphi, iSub), "", 30, 1.1, 1.13);
        prCos_theta_LamBar[iCent][iphi][iSub] = new TProfile(Form("prCos_theta_LamBar_%i_%i_%i", iCent, iphi, iSub), "", 30, 1.1, 1.13);
        prSin_theta_LamBar[iCent][iphi][iSub] = new TProfile(Form("prSin_theta_LamBar_%i_%i_%i", iCent, iphi, iSub), "", 30, 1.1, 1.13);
        prSin_diffPhiPsi1Sin_theta_LamBar[iCent][iphi][iSub] = new TProfile(Form("prSin_diffPhiPsi1Sin_theta_LamBar_%i_%i_%i", iCent, iphi, iSub), "", 30, 1.1, 1.13);
        prCos_diffPhiPsi1Sin_theta_LamBar[iCent][iphi][iSub] = new TProfile(Form("prCos_diffPhiPsi1Sin_theta_LamBar_%i_%i_%i", iCent, iphi, iSub), "", 30, 1.1, 1.13);

        //Lambda invMass distributions
        InvMLamDist[iCent][iphi][iSub] = new TH1F(Form("InvMLamDist_%i_%i_%i", iCent, iphi, iSub), Form("InvMLamDist_%i_%i_%i", iCent, iphi, iSub), 125, 1.1, 1.13);
        InvMLamBarDist[iCent][iphi][iSub] = new TH1F(Form("InvMLamBarDist_%i_%i_%i", iCent, iphi, iSub), Form("InvMLamBarDist_%i_%i_%i", iCent, iphi, iSub), 125, 1.1, 1.13);
      }

      prSin_diffPhiPsi1[iCent][iphi][0]->SetTitle(Form("Sin(Psi1-phi) Vs InvMLam for cent=(%i), phi-Psi2_e/w=(%i)", iCent, iphi));
      prSin_diffPhiPsi1[iCent][iphi][1]->SetTitle(Form("Sin(Psi1-phi) Vs InvMLam for cent=(%i), phi-Psi2_comb=(%i)", iCent, iphi));
      prSin_diffPhiPsi1[iCent][iphi][2]->SetTitle(Form("Sin(Psi1-phi) Vs InvMLam for cent=(%i), phi-Psi3_e/w=(%i)", iCent, iphi));
      prSin_diffPhiPsi1[iCent][iphi][3]->SetTitle(Form("Sin(Psi1-phi) Vs InvMLam for cent=(%i), phi-Psi3_comb=(%i)", iCent, iphi));

      prCos_diffPhiPsi1[iCent][iphi][0]->SetTitle(Form("Cos(Psi1-phi) Vs InvMLam for cent=(%i), phi-Psi2_e/w=(%i)", iCent, iphi));
      prCos_diffPhiPsi1[iCent][iphi][1]->SetTitle(Form("Cos(Psi1-phi) Vs InvMLam for cent=(%i), phi-Psi2_comb=(%i)", iCent, iphi));
      prCos_diffPhiPsi1[iCent][iphi][2]->SetTitle(Form("Cos(Psi1-phi) Vs InvMLam for cent=(%i), phi-Psi3_e/w=(%i)", iCent, iphi));
      prCos_diffPhiPsi1[iCent][iphi][3]->SetTitle(Form("Cos(Psi1-phi) Vs InvMLam for cent=(%i), phi-Psi3_comb=(%i)", iCent, iphi));
      
      prCos_theta[iCent][iphi][0]->SetTitle(Form("Cos(theta) Vs InvMLam for cent=(%i), phi-Psi2_e/w=(%i)", iCent, iphi));
      prCos_theta[iCent][iphi][1]->SetTitle(Form("Cos(theta) Vs InvMLam for cent=(%i), phi-Psi2_comb=(%i)", iCent, iphi));
      prCos_theta[iCent][iphi][2]->SetTitle(Form("Cos(theta) Vs InvMLam for cent=(%i), phi-Psi3_e/w=(%i)", iCent, iphi));
      prCos_theta[iCent][iphi][3]->SetTitle(Form("Cos(theta) Vs InvMLam for cent=(%i), phi-Psi3_comb=(%i)", iCent, iphi));
      
      prSin_theta[iCent][iphi][0]->SetTitle(Form("Sin(theta) Vs InvMLam for cent=(%i), phi-Psi2_e/w=(%i)", iCent, iphi));
      prSin_theta[iCent][iphi][1]->SetTitle(Form("Sin(theta) Vs InvMLam for cent=(%i), phi-Psi2_comb=(%i)", iCent, iphi));
      prSin_theta[iCent][iphi][2]->SetTitle(Form("Sin(theta) Vs InvMLam for cent=(%i), phi-Psi3_e/w=(%i)", iCent, iphi));
      prSin_theta[iCent][iphi][3]->SetTitle(Form("Sin(theta) Vs InvMLam for cent=(%i), phi-Psi3_comb=(%i)", iCent, iphi));
      
      prSin_diffPhiPsi1Sin_theta[iCent][iphi][0]->SetTitle(Form("Sin(Psi1-phi)Sin(Theta) Vs InvMLam for cent=(%i), phi-Psi2_e/w=(%i)", iCent, iphi));
      prSin_diffPhiPsi1Sin_theta[iCent][iphi][1]->SetTitle(Form("Sin(Psi1-phi)Sin(Theta) Vs InvMLam for cent=(%i), phi-Psi2_comb=(%i)", iCent, iphi));
      prSin_diffPhiPsi1Sin_theta[iCent][iphi][2]->SetTitle(Form("Sin(Psi1-phi)Sin(Theta) Vs InvMLam for cent=(%i), phi-Psi3_e/w=(%i)", iCent, iphi));
      prSin_diffPhiPsi1Sin_theta[iCent][iphi][3]->SetTitle(Form("Sin(Psi1-phi)Sin(Theta) Vs InvMLam for cent=(%i), phi-Psi3_comb=(%i)", iCent, iphi));
      
      prCos_diffPhiPsi1Sin_theta[iCent][iphi][0]->SetTitle(Form("Cos(Psi1-phi)Sin(Theta) Vs InvMLam for cent=(%i), phi-Psi2_e/w=(%i)", iCent, iphi));
      prCos_diffPhiPsi1Sin_theta[iCent][iphi][1]->SetTitle(Form("Cos(Psi1-phi)Sin(Theta) Vs InvMLam for cent=(%i), phi-Psi2_comb=(%i)", iCent, iphi));
      prCos_diffPhiPsi1Sin_theta[iCent][iphi][2]->SetTitle(Form("Cos(Psi1-phi)Sin(Theta) Vs InvMLam for cent=(%i), phi-Psi3_e/w=(%i)", iCent, iphi));
      prCos_diffPhiPsi1Sin_theta[iCent][iphi][3]->SetTitle(Form("Cos(Psi1-phi)Sin(Theta) Vs InvMLam for cent=(%i), phi-Psi3_comb=(%i)", iCent, iphi));

      //LamBar
      prSin_diffPhiPsi1_LamBar[iCent][iphi][0]->SetTitle(Form("Sin(Psi1-phi) Vs InvMLamBar for cent=(%i), phi-Psi2_e/w=(%i)", iCent, iphi));
      prSin_diffPhiPsi1_LamBar[iCent][iphi][1]->SetTitle(Form("Sin(Psi1-phi) Vs InvMLamBar for cent=(%i), phi-Psi2_comb=(%i)", iCent, iphi));
      prSin_diffPhiPsi1_LamBar[iCent][iphi][2]->SetTitle(Form("Sin(Psi1-phi) Vs InvMLamBar for cent=(%i), phi-Psi3_e/w=(%i)", iCent, iphi));
      prSin_diffPhiPsi1_LamBar[iCent][iphi][3]->SetTitle(Form("Sin(Psi1-phi) Vs InvMLamBar for cent=(%i), phi-Psi3_comb=(%i)", iCent, iphi));

      prCos_diffPhiPsi1_LamBar[iCent][iphi][0]->SetTitle(Form("Cos(Psi1-phi) Vs InvMLamBar for cent=(%i), phi-Psi2_e/w=(%i)", iCent, iphi));
      prCos_diffPhiPsi1_LamBar[iCent][iphi][1]->SetTitle(Form("Cos(Psi1-phi) Vs InvMLamBar for cent=(%i), phi-Psi2_comb=(%i)", iCent, iphi));
      prCos_diffPhiPsi1_LamBar[iCent][iphi][2]->SetTitle(Form("Cos(Psi1-phi) Vs InvMLamBar for cent=(%i), phi-Psi3_e/w=(%i)", iCent, iphi));
      prCos_diffPhiPsi1_LamBar[iCent][iphi][3]->SetTitle(Form("Cos(Psi1-phi) Vs InvMLamBar for cent=(%i), phi-Psi3_comb=(%i)", iCent, iphi));
      
      prCos_theta_LamBar[iCent][iphi][0]->SetTitle(Form("Cos(theta) Vs InvMLamBar for cent=(%i), phi-Psi2_e/w=(%i)", iCent, iphi));
      prCos_theta_LamBar[iCent][iphi][1]->SetTitle(Form("Cos(theta) Vs InvMLamBar for cent=(%i), phi-Psi2_comb=(%i)", iCent, iphi));
      prCos_theta_LamBar[iCent][iphi][2]->SetTitle(Form("Cos(theta) Vs InvMLamBar for cent=(%i), phi-Psi3_e/w=(%i)", iCent, iphi));
      prCos_theta_LamBar[iCent][iphi][3]->SetTitle(Form("Cos(theta) Vs InvMLamBar for cent=(%i), phi-Psi3_comb=(%i)", iCent, iphi));
      
      prSin_theta_LamBar[iCent][iphi][0]->SetTitle(Form("Sin(theta) Vs InvMLamBar for cent=(%i), phi-Psi2_e/w=(%i)", iCent, iphi));
      prSin_theta_LamBar[iCent][iphi][1]->SetTitle(Form("Sin(theta) Vs InvMLamBar for cent=(%i), phi-Psi2_comb=(%i)", iCent, iphi));
      prSin_theta_LamBar[iCent][iphi][2]->SetTitle(Form("Sin(theta) Vs InvMLamBar for cent=(%i), phi-Psi3_e/w=(%i)", iCent, iphi));
      prSin_theta_LamBar[iCent][iphi][3]->SetTitle(Form("Sin(theta) Vs InvMLamBar for cent=(%i), phi-Psi3_comb=(%i)", iCent, iphi));
      
      prSin_diffPhiPsi1Sin_theta_LamBar[iCent][iphi][0]->SetTitle(Form("Sin(Psi1-phi)Sin(Theta) Vs InvMLamBar for cent=(%i), phi-Psi2_e/w=(%i)", iCent, iphi));
      prSin_diffPhiPsi1Sin_theta_LamBar[iCent][iphi][1]->SetTitle(Form("Sin(Psi1-phi)Sin(Theta) Vs InvMLamBar for cent=(%i), phi-Psi2_comb=(%i)", iCent, iphi));
      prSin_diffPhiPsi1Sin_theta_LamBar[iCent][iphi][2]->SetTitle(Form("Sin(Psi1-phi)Sin(Theta) Vs InvMLamBar for cent=(%i), phi-Psi3_e/w=(%i)", iCent, iphi));
      prSin_diffPhiPsi1Sin_theta_LamBar[iCent][iphi][3]->SetTitle(Form("Sin(Psi1-phi)Sin(Theta) Vs InvMLamBar for cent=(%i), phi-Psi3_comb=(%i)", iCent, iphi));
      
      prCos_diffPhiPsi1Sin_theta_LamBar[iCent][iphi][0]->SetTitle(Form("Cos(Psi1-phi)Sin(Theta) Vs InvMLamBar for cent=(%i), phi-Psi2_e/w=(%i)", iCent, iphi));
      prCos_diffPhiPsi1Sin_theta_LamBar[iCent][iphi][1]->SetTitle(Form("Cos(Psi1-phi)Sin(Theta) Vs InvMLamBar for cent=(%i), phi-Psi2_comb=(%i)", iCent, iphi));
      prCos_diffPhiPsi1Sin_theta_LamBar[iCent][iphi][2]->SetTitle(Form("Cos(Psi1-phi)Sin(Theta) Vs InvMLamBar for cent=(%i), phi-Psi3_e/w=(%i)", iCent, iphi));
      prCos_diffPhiPsi1Sin_theta_LamBar[iCent][iphi][3]->SetTitle(Form("Cos(Psi1-phi)Sin(Theta) Vs InvMLamBar for cent=(%i), phi-Psi3_comb=(%i)", iCent, iphi));
    }
    
  }
}


double StKFParticleAnalysisMaker::GetPsi(int iOrd, double Qx, double Qy){
    double Psi;
    Psi = atan2(Qy, Qx)/iOrd;
    if(Psi<0.) Psi += 2*TMath::Pi() / iOrd;
    return Psi;
}


void StKFParticleAnalysisMaker::GetCentring() {
  TFile *file = new TFile("/star/u/mmorozov/14p5MakeCorrections/correctionsEP_output/output6sem_Psi1_2_3Corr.root", "read");	
  
  for(int iSub=0; iSub!=2*nSub; iSub++){
    Qvec1Prof_TH[iSub] = (TH1F*)file->Get(Form("Qvec1Prof_%i", iSub));
    Qvec2Prof_TH[iSub] = (TH1F*)file->Get(Form("Qvec2Prof_%i", iSub));
    Qvec3Prof_TH[iSub] = (TH1F*)file->Get(Form("Qvec3Prof_%i", iSub));
  }
}


void StKFParticleAnalysisMaker::GetFlattening() {
  TFile *file = new TFile("/star/u/mmorozov/14p5MakeCorrections/correctionsEP_output/outputCentred_Psi1_2_3Corr.root", "read");

  for(int iProf=0; iProf!=10; iProf++){
    for(int iSub=0; iSub!=nSub; iSub++){
      Coef_A_n_TH_Psi1[iProf][iSub] = (TH1F*)file->Get(Form("Coef_A_n_Psi1_%i_%i", iProf, iSub));
      Coef_B_n_TH_Psi1[iProf][iSub] = (TH1F*)file->Get(Form("Coef_B_n_Psi1_%i_%i", iProf, iSub));

      Coef_A_n_TH_Psi2[iProf][iSub] = (TH1F*)file->Get(Form("Coef_A_n_Psi2_%i_%i", iProf, iSub));
      Coef_B_n_TH_Psi2[iProf][iSub] = (TH1F*)file->Get(Form("Coef_B_n_Psi2_%i_%i", iProf, iSub));

      Coef_A_n_TH_Psi3[iProf][iSub] = (TH1F*)file->Get(Form("Coef_A_n_Psi3_%i_%i", iProf, iSub));
      Coef_B_n_TH_Psi3[iProf][iSub] = (TH1F*)file->Get(Form("Coef_B_n_Psi3_%i_%i", iProf, iSub));
    }
  }
}

void StKFParticleAnalysisMaker::GetWeightCorr() {
  TFile *file = new TFile("/star/u/mmorozov/14p5MakeCorrections/correctionsEP_output/outputV1.root", "read");

  for(int iCent=0; iCent!=9; iCent++){
    for(int iSub=0; iSub!=2; iSub++){
      v1_average[iCent][iSub] = (TProfile*)file->Get(Form("v1_average_%i_%i", iCent, iSub));
    }
        
  }
}
