//#include "StRoot/StKFParticleAnalysisMaker/StKFParticlePerformanceInterface.h"
//#include "StRoot/StKFParticleAnalysisMaker/StKFParticleAnalysisMaker.h"

void removeValuesFromFile(const Char_t *filename)
{
    std::ifstream input(filename);
    if (!input.is_open()) {
        std::cerr << "Error: could not open file " << filename << '\n';
        return;
    }

    std::ofstream output(Form("%s.tmp", filename)); // create temporary file for output
    if (!output.is_open()) {
        std::cerr << "Error: could not create temporary file\n";
        return;
    }

    std::string line;
    while (std::getline(input, line)) {
        size_t pos = line.find(' '); // find first space character
        if (pos != std::string::npos) {
            output << line.substr(0, pos) << '\n'; // write filename to temporary file
        }
    }

    input.close();
    output.close();

    if (std::remove(filename)) { // delete original file
        std::cerr << "Error: could not delete file " << filename << '\n';
        return;
    }

    if (std::rename(Form("%s.tmp", filename), filename)) { // rename temporary file to original name
        std::cerr << "Error: could not rename temporary file\n";
        return;
    }

    std::cout << "Values removed from file " << filename << '\n';
}


void analysis(Char_t *inFileName = "../inputfiles/st_physics_20124023_raw_3500008.picoDst.root" ,// "/star/u/alpatov/KFParticle/st_physics_19142027_raw_2500009.picoDst.root",
Char_t *outFileName = "pico.root", Char_t *runFile = "095", bool isPico = true)
{
#if !defined(__CINT__)
  std::cout << "This code cannot be compiled" << std::endl;
#elsev
  //  gSystem->SetFPEMask(kInvalid | kDivByZero | kOverflow );
  gSystem->Load("StEpdUtil");
  
  //removeValuesFromFile(inFileName);

  gROOT->LoadMacro("lMuDst.C");

    
  TString input;
  TString output;
  int year = 2017;//2014
//  Int_t N = 1000;
  Int_t N = 1000000000;
              
  if(isPico)
  {
   // input = "/star/u/alpatov/KFParticle/inputData/*.picoDst.root", 
    input = inFileName;//"/star/u/alpatov/KFParticle/st_physics_19142027_raw_2500009.picoDst.root";
    output = outFileName;//"pico.root";
    lMuDst(-1,input.Data(),"ry2017,RpicoDst,mysql,kfpAna,quiet,nodefault",output);
  }
  else
  {
    input = "/star/u/mzyzak/KFParticle_NewPicoFormat/inputData/*.MuDst.root", 
    output = "mu.root";
    lMuDst(-1,input.Data(),"ry2016,picoEvt,RMuDst,mysql,kfpAna,quiet,nodefault",output.Data());
  }

//  gSystem->Load("StFemtoV0"); //don't forget this stuff!!!
    

  StKFParticleAnalysisMaker* kfpAnalysis = (StKFParticleAnalysisMaker*) StMaker::GetTopChain()->Maker("KFParticleAnalysis");
  kfpAnalysis->SetMyStuff(outFileName, runFile);
//  kfpAnalysis->SetRunFile(runFile);
  if(!isPico) 
  {
    kfpAnalysis->AnalyseMuDst();
    kfpAnalysis->ProcessSignal();
  }
  
  if(year == 2016)
  {
    kfpAnalysis->UseTMVA();
    // D0->Kpi
    kfpAnalysis->SetTMVABinsD0("0:2:3:4:5:6:7:8:9","-1:1000");
    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_0_1_pt0_80_BDT.weights.xml", 0.075, 0);
    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_2_2_pt0_80_BDT.weights.xml", 0.05,  1);
    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_3_3_pt0_80_BDT.weights.xml", 0.05,  2);
    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_4_4_pt0_80_BDT.weights.xml", 0.1,   3);
    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_5_5_pt0_80_BDT.weights.xml", 0.1,   4);
    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_6_6_pt0_80_BDT.weights.xml", 0.125, 5);
    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_7_7_pt0_80_BDT.weights.xml", 0.125, 6);
    kfpAnalysis->SetTMVAcutsD0("/gpfs01/star/pwg/kocmic/TMVA/Centrality/TMVA/D0/weights/TMVAClassification_8_8_pt0_80_BDT.weights.xml", 0.125, 7);
    
    kfpAnalysis->RunCentralityAnalysis();
    kfpAnalysis->SetCentralityFile("/gpfs01/star/pwg/mzyzak/Femto/Template/Centrality/centrality_2016.txt");
  }
  if(year == 2014)
  {
    kfpAnalysis->UseTMVA();
    kfpAnalysis->SetTMVAcutsD0(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/D0.xml",    0.1);
    kfpAnalysis->SetTMVAcutsDPlus("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/DPlus.xml", 0.05);
    kfpAnalysis->SetTMVAcutsDs(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/Ds.xml",    0.075);
    kfpAnalysis->SetTMVAcutsLc(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/Lc.xml",    0.1);
    kfpAnalysis->SetTMVAcutsD0KK( "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/D0KK.xml",  0.125);
    kfpAnalysis->SetTMVAcutsD04(  "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/D04.xml",  -0.05);
    kfpAnalysis->SetTMVAcutsBPlus("/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/BPlus.xml",-0.1);
    kfpAnalysis->SetTMVAcutsB0(   "/gpfs01/star/pwg/mzyzak/Femto/Template/TMVA/2014/B0.xml",   -0.1);
    
//     kfpAnalysis->RunCentralityAnalysis();
//     kfpAnalysis->SetCentralityFile("/gpfs01/star/pwg/mzyzak/Femto/Template/Centrality/centrality_2014.txt");
  }
  
  chain->Init();

//  StKFParticleInterface::instance()->CleanLowPVTrackEvents();
//     StKFParticleInterface::instance()->UseHFTTracksOnly();
//  StKFParticleInterface::instance()->SetOutFileNameMine(outFileName);
//  StKFParticleInterface::instance()->SetRunFileNameMine(runFile);
  StKFParticleInterface::instance()->SetSoftKaonPIDMode();
  StKFParticleInterface::instance()->SetSoftTofPidMode();
//  StKFParticleInterface::instance()->SetChiPrimaryCut(10);
  
//  StKFParticleInterface::instance()->SetPtCutCharm(0.2);
//  StKFParticleInterface::instance()->SetChiPrimaryCutCharm(8);
//  StKFParticleInterface::instance()->SetLdLCut2D(5);
//  StKFParticleInterface::instance()->SetLdLCutCharmManybodyDecays(3);//3
//  StKFParticleInterface::instance()->SetChi2TopoCutCharmManybodyDecays(10);
//  StKFParticleInterface::instance()->SetChi2CutCharmManybodyDecays(3);
//  StKFParticleInterface::instance()->SetLdLCutCharm2D(3);//3
//  StKFParticleInterface::instance()->SetChi2TopoCutCharm2D(10);
//  StKFParticleInterface::instance()->SetChi2CutCharm2D(3);
  StKFParticleInterface::instance()->SetLdLCutXiOmega(5);//10

  StKFParticleInterface::instance()->CleanLowPVTrackEvents();
//  StKFParticleInterface::instance()->SetChiPrimaryCut(10);
  //Add decays to the reconstruction list
//  StKFParticleInterface::instance()->AddDecayToReconstructionList(  310);
  StKFParticleInterface::instance()->AddDecayToReconstructionList( 3122); //Lambda
  StKFParticleInterface::instance()->AddDecayToReconstructionList(-3122); //antiLambda
//  StKFParticleInterface::instance()->AddDecayToReconstructionList( 3312); //Xi
//  StKFParticleInterface::instance()->AddDecayToReconstructionList(-3312);
//  StKFParticleInterface::instance()->AddDecayToReconstructionList(3334);
//  StKFParticleInterface::instance()->AddDecayToReconstructionList(-3334);
  
  Long64_t nevent = N;
  if(isPico)
  {
    StPicoDstMaker* maker = (StPicoDstMaker *) StMaker::GetTopChain()->Maker("PicoDst");
    if (! maker) return;
    maker->SetStatus("*",1);
    TChain *tree = maker->chain();
    Long64_t nentries = tree->GetEntries();
    if (nentries <= 0) return;
    nevent = TMath::Min(nevent,nentries);
    cout << nentries << " events in chain " << nevent << " will be read." << endl;
  }
  
  chain->EventLoop(nevent);
#endif
  
}
