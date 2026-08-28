// $Id: StKFParticleAnalysisMaker.h,v 1.16 2014/08/06 11:43:53 jeromel Exp $
/*!
 * \class  StKFParticleAnalysisMaker
 * \author Maksym Zyzak
 * \date   2017/10/17
 * \brief  class for analysis of PicoDst
 */                                                                      
#ifndef STAR_StKFParticleAnalysisMaker
#define STAR_StKFParticleAnalysisMaker
//#define __DEVT__
#ifndef StMaker_H
#include "StMaker.h"
#endif
#include "TMVA/Reader.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile3D.h"
#include "TLorentzVector.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"




#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEpdHit.h"

const int PolPartNum = 12;
const int DetNum = 3;
const int nSub = 3;//Numbers of Psi
		   //0-West TPC
		   //1-East TPC
		   //2-Comb TPC
		   //3-E+W TPC
		   //4-West EPD
		   //5-East EPD
		   //6-E+W EPD

class StKFParticleInterface;
class StKFParticlePerformanceInterface;
class KFParticle;
class StPicoDst;
class StPicoEvent;
class StMuDst;
class TNtuple;
class TFile;
class TChain;
class StRefMultCorr;

class StKFParticleAnalysisMaker : public StMaker {
 private:
  static const int fNNTuples = 8;
  Char_t                mBeg[1];        //!
  StMuDst                          *fMuDst;
  StPicoDst                        *fPicoDst;                          //!
  StKFParticleInterface            *fStKFParticleInterface;            //!
  StKFParticlePerformanceInterface *fStKFParticlePerformanceInterface; //!
  TNtuple* fCutsNTuple[fNNTuples];
  TFile* fNTupleFile[fNNTuples];
  int fNTuplePDG[fNNTuples];
  TString fNtupleNames[fNNTuples];
  TString fNtupleCutNames[fNNTuples];
  std::vector<TString> fDaughterNames[fNNTuples];
  vector< vector<TString> > fTMVACutFile[fNNTuples];
  vector< vector<double> > fTMVACut[fNNTuples];
  vector< vector<TMVA::Reader*> > fTMVAReader[fNNTuples];
  std::vector<int> fTMVACentralityBins[fNNTuples];
  std::vector<double> fTMVAPtBins[fNNTuples];
  Char_t                mEnd[1];        //!
  std::vector<float> fTMVAParticleParameters[fNNTuples];
  int fNTrackTMVACuts;
  bool fIsPicoAnalysis;
  int fdEdXMode;
  Bool_t fStoreTmvaNTuples;
  Bool_t fProcessSignal;
  Bool_t fCollectTrackHistograms;
  Bool_t fCollectPIDHistograms;
  Bool_t fTMVAselection;
  
  //Centrality and flow
  Bool_t fFlowAnalysis;
  TChain* fFlowChain;
  int fFlowRunId;
  int fFlowEventId;
  int fCentrality;
  std::vector<TString> fFlowFiles;
  std::map<long, int> fFlowMap;
  
  bool fRunCentralityAnalysis;
  StRefMultCorr *fRefmultCorrUtil;
  TString fCentralityFile;
  
  bool fAnalyseDsPhiPi;

  void GetDaughterParameters(const int iReader, int& iDaughterTrack, int& iDaughterParticle, KFParticle& particle);
  void GetParticleParameters(const int iReader, KFParticle& particle);
  long  GetUniqueEventId(const int iRun, const int iEvent) const;
  
  int GetTMVACentralityBin(int iReader, int centrality);
  int GetTMVAPtBin(int iReader, double pt);
  void SetTMVACentralityBins(int iReader, TString bins);
  void SetTMVAPtBins(int iReader, TString bins);
  void SetTMVABins(int iReader, TString centralityBins="-1:1000", TString ptBins="-1.:1000.");
  // Here will be my stuff
  void InitEpFinders();

  //My variables
  const Char_t *runnumber;
  double QWeight_1[nSub], QWeight_2[nSub];
  double Qvec_1[2*nSub], Qvec_2[2*nSub], Qvec_3[2*nSub];
  double Psi1[nSub], Psi2[nSub], Psi3[nSub];
  double deltaPsi1[nSub], deltaPsi2[nSub], deltaPsi3[nSub];
  int PP, TT, EW, iPsi, iSide, row, phi_bin;
  double sin_diffPsi1Phi, cos_diffPsi1Phi, delta_phi;
  
  //My hist
  TProfile *v1_average[9][2];

  TH1F *Coef_A_n_TH_Psi1[10][nSub];
  TH1F *Coef_B_n_TH_Psi1[10][nSub];
  TH1F *Coef_A_n_TH_Psi2[10][nSub];
  TH1F *Coef_B_n_TH_Psi2[10][nSub];
  TH1F *Coef_A_n_TH_Psi3[10][nSub];
  TH1F *Coef_B_n_TH_Psi3[10][nSub];
  
  TH1F *Qvec1Prof_TH[nSub*2];
  TH1F *Qvec2Prof_TH[nSub*2];
  TH1F *Qvec3Prof_TH[nSub*2]; 

  TH1F *Qvec1Hist[9][2*nSub];
  TH1F *Qvec2Hist[9][2*nSub];
  TH1F *Qvec3Hist[9][2*nSub];
  TH1F *Psi1Hist[9][nSub];
  TH1F *Psi2Hist[9][nSub];
  TH1F *Psi3Hist[9][nSub];

  TProfile *CosOfDiff_1;
  TProfile *CosOfDiff_2;
  TProfile *CosOfDiff_3;

  //For Pz
  TProfile *prSin_diffPhiPsi1[9][6][4];
  TProfile *prCos_diffPhiPsi1[9][6][4];
  TProfile *prCos_theta[9][6][4];
  TProfile *prSin_theta[9][6][4];
  TProfile *prSin_diffPhiPsi1Sin_theta[9][6][4];
  TProfile *prCos_diffPhiPsi1Sin_theta[9][6][4];

  TProfile *prSin_diffPhiPsi1_LamBar[9][6][4];
  TProfile *prCos_diffPhiPsi1_LamBar[9][6][4];
  TProfile *prCos_theta_LamBar[9][6][4];
  TProfile *prSin_theta_LamBar[9][6][4];
  TProfile *prSin_diffPhiPsi1Sin_theta_LamBar[9][6][4];
  TProfile *prCos_diffPhiPsi1Sin_theta_LamBar[9][6][4];

  TH1F *InvMLamDist[9][6][4];
  TH1F *InvMLamBarDist[9][6][4];

  StEpdEpFinder *etaFinder;
  StEpdEpFinder *etaVzFinder[14];
  StEpdGeom *mEpdGeom;
  
  bool isFileRead;
  Char_t *frunFile, *foutFile; 
  int CentralityBin(int refMult);
  bool GoodRun(StPicoEvent *event);
  bool EventCut(StPicoEvent *event);
  int VzBin(float Vz);
  bool PrithwishRun(StPicoEvent *event);
  int JoeyCentrality(int Multiplicity);
//  StPicoTrack *FindTrack(int trackId, StPicoDst *dst);
  bool IsKfGoof(KFParticle particle);
  bool isEventPlaneMethod = false, isA2A = true;
  bool isTopo = false;
  float mTrgEff;
  //Its My Functions MAXIM
  void CreateEPDist();
  void CreateKFPHists();
  void GetCentring();
  void GetWeightCorr();
  void GetFlattening();
  double GetPsi(int iOrd, double Qx, double Qy);
  
/*
  void FillTH1F(int PolNPart, StFemtoV0 particle, int CentralityID);
  void FillTH2F(int PolNPart, StFemtoV0 particle, int CentralityID);
  void FillTH3F(int PolNPart, int detNum, float Psi, float phi, StFemtoV0 particle, int iCent);*/
 public: 
  StKFParticleAnalysisMaker(const char *name="KFParticleAnalysis");
  virtual       ~StKFParticleAnalysisMaker();
  virtual Int_t  Init();
  virtual Int_t  InitRun(Int_t runumber);
  void           BookVertexPlots();
  void SetMyStuff(Char_t *outFileName, Char_t *runFileName);
  virtual Int_t  Make();
  virtual Int_t  Finish();
  Bool_t         Check();
  void AnalysePicoDst() { fIsPicoAnalysis = true;  }
  void AnalyseMuDst()   { fIsPicoAnalysis = false; }
  static void    PrintMem(const Char_t *opt = "");
  virtual const char *GetCVS() const {
    static const char cvs[]="Tag $Name:  $ $Id: StKFParticleAnalysisMaker.h,v 1.0 2017/10/07 11:43:53 mzyzak Exp $ built " __DATE__ " " __TIME__ ; 
    return cvs;
  }
  void ProcessSignal() { fProcessSignal = true; }
  void StoreTMVANtuples() { fStoreTmvaNTuples = true; }
  void CollectTrackHistograms() { fCollectTrackHistograms = true; }
  void CollectPIDHistograms() { fCollectPIDHistograms = true; }
  void UseTMVA() { fTMVAselection = true; }
  void SetTMVABinsD0   (TString centralityBins, TString ptBins) { SetTMVABins(0, centralityBins, ptBins); }
  void SetTMVABinsDPlus(TString centralityBins, TString ptBins) { SetTMVABins(1, centralityBins, ptBins); }
  void SetTMVABinsDs   (TString centralityBins, TString ptBins) { SetTMVABins(2, centralityBins, ptBins); }
  void SetTMVABinsLc   (TString centralityBins, TString ptBins) { SetTMVABins(3, centralityBins, ptBins); }
  void SetTMVABinsD0KK (TString centralityBins, TString ptBins) { SetTMVABins(4, centralityBins, ptBins); }
  void SetTMVABinsD04  (TString centralityBins, TString ptBins) { SetTMVABins(5, centralityBins, ptBins); }
  void SetTMVABinsBPlus(TString centralityBins, TString ptBins) { SetTMVABins(6, centralityBins, ptBins); }
  void SetTMVABinsB0   (TString centralityBins, TString ptBins) { SetTMVABins(7, centralityBins, ptBins); }
  void SetTMVAcutsD0   (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[0][iCentralityBin][iPtBin] = file; fTMVACut[0][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsDPlus(TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[1][iCentralityBin][iPtBin] = file; fTMVACut[1][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsDs   (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[2][iCentralityBin][iPtBin] = file; fTMVACut[2][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsLc   (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[3][iCentralityBin][iPtBin] = file; fTMVACut[3][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsD0KK (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[4][iCentralityBin][iPtBin] = file; fTMVACut[4][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsD04  (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[5][iCentralityBin][iPtBin] = file; fTMVACut[5][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsBPlus(TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[6][iCentralityBin][iPtBin] = file; fTMVACut[6][iCentralityBin][iPtBin] = cut; }
  void SetTMVAcutsB0   (TString file, double cut, int iCentralityBin = 0, int iPtBin = 0) { fTMVACutFile[7][iCentralityBin][iPtBin] = file; fTMVACut[7][iCentralityBin][iPtBin] = cut; }
  
  void RunFlowAnalysis()         { fFlowAnalysis = true; }
  void AddFlowFile(TString file) { fFlowFiles.push_back(file); }
  
  void RunCentralityAnalysis() { fRunCentralityAnalysis = true; }
  void SetCentralityFile(TString file) { fCentralityFile = file; }
  
  void AnalyseDsPhiPi() { fAnalyseDsPhiPi = true; }

  //Here will be my stuff
  void SetRunNumber(const Char_t *mrunnumber){ runnumber = mrunnumber; }
  
  ClassDef(StKFParticleAnalysisMaker,0)   //
};
#endif
// $Log: StKFParticleAnalysisMaker.h,v $
