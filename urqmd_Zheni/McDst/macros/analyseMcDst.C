/**
 * \brief Example of how to read a file (list of files) using McDst classes
 *
 * analyseMcDst.C is an example of reading McDst format.
 * One can use either uDst file or a list of mcDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 */

// C++ headers
#include <iostream>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TMath.h"
#include "TF1.h"
#include "TRandom3.h"

// McDst headers
#include "../McDstReader.h"
#include "../McDst.h"
#include "../McEvent.h"
#include "../McParticle.h"
#include "../McRun.h"

// inFile - is a name of name.uDst.root file or a name
//          of a name.lis(t) files that contains a list of
//          name1.uDst.root files
//_________________
void analyseMcDst(const Char_t *inFile,
		  const Char_t *oFileName) {

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;
/*
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  R__LOAD_LIBRARY(../libMcDst)
#else
    gSystem->Load("/star/u/annakraeva/urqmd_hbt_Zheni/StRoot/McDst/libMcDst.so");
#endif
*/
  McDstReader* myReader = new McDstReader(inFile);
  myReader->Init();

  // This is a way if you want to spead up IO
  std::cout << "Explicit read status for some branches" << std::endl;
  myReader->setStatus("*",0);
  myReader->setStatus("Event",1);
  myReader->setStatus("Particle",1);
  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;

  if( !myReader->chain() ) {
    std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInTree = myReader->tree()->GetEntries();
  std::cout << "eventsInTree: "  << eventsInTree << std::endl;
  Long64_t events2read = myReader->chain()->GetEntries();

  std::cout << "Number of events to read: " << events2read << std::endl;

  TFile *oFile = new TFile(oFileName, "RECREATE");

  // Histogramming
  // Event
  TH2D *hImpactParVsNch = new TH2D("hImpactParVsNch","Impact parameter vs. Nch;Nch;Impact parameter (fm)", 300, -0.5, 599.5, 200, 0., 18.);//was  300, -0.5, 599.5, 200, 0., 18.
  TH1D *hNch = new TH1D("hNCh","Number of charged particles;Nch;Entries", 300, -0.5, 599.5); //was 300, -0.5, 599.5
  TH1D *hSqrtSnn = new TH1D("hSqrtSnn","Collision energy;#sqrt{s_{NN}} (GeV);Entries",10, 0., 10. );
  TH1D *hPDG = new TH1D("hPDG","PDG of particles;PDG;Entries",4000,-500.,3500.);
  TH1D *hMult = new TH1D("hMult","hMult",600, -0.5, 599.5);
  TH1D *hMultwoProtonsLess200fm = new TH1D("hMultwoProtonsLess200fm","hMultwoProtonsLess200fm",600, -0.5, 599.5);
  TH1D *hMultEff2D = new TH1D("hMultEff2D","hMultEff2D",600, -0.5, 599.5);
  TH1D *hMultEff2DwoProtonsLess200fm = new TH1D("hMultEff2DwoProtonsLess200fm","hMultEff2DwoProtonsLess200fm",600, -0.5, 599.5);
  TH1D *hEta_LS_afterEtacut = new TH1D("hEta_LS_afterEtacut","#eta;#eta;Entries", 200, -15., 15.);
  TH1D *hEta_LS_afterRapcut = new TH1D("hEta_LS_afterRapcut","#eta;#eta;Entries", 200, -15., 15.);
  TH1D *hEta_LS_afterEtaRapcut = new TH1D("hEta_LS_afterEtaRapcut","#eta;#eta;Entries", 200, -15., 15.);
  TH1D *hRap_LS_afterEtacut = new TH1D("hRap_LS_afterEtacut","Rapidity;y;Entries", 1000,-10.,10.);
  TH1D *hRap_LS_afterRapcut = new TH1D("hRap_LS_afterRapcut","Rapidity;y;Entries", 1000,-10.,10.);
  TH1D *hRap_LS_afterEtaRapcut = new TH1D("hRap_LS_afterEtaRapcut","Rapidity;y;Entries", 1000,-10.,10.);
  TH1D *hP_aftercut = new TH1D("hP_aftercut","p;p (GeV/c); Entries",200,0.,2.5);
  TH1D *hPt_aftercut = new TH1D("hPt_aftercut","p_{T};p_{T} (GeV/c);Entries",200,0.,2.5);
  // Track
  TH1D *hImpact = new TH1D("hImpact","Impact parameter;Impact parameter;Entries",100,0,20);
  TH1D *hCharge = new TH1D("hCharge","Charge of particles;Charge;Entries",5,-2,3);
  TH1D *hX = new TH1D("hX","X position;X (fm);Entries",200,-20.,20.);
  TH1D *hY = new TH1D("hY","Y position;Y (fm);Entries",200,-20.,20.);
  TH1D *hZ = new TH1D("hZ","Z position;Z (fm);Entries",500,-60.,60.); //was -20,20
  TH2D *hXY = new TH2D("hXY","XY position;X (fm);Y (fm)",200,-20.,20.,200,-20.,20.);
  TH2D *hXZ = new TH2D("hXZ","XZ position;X (fm);Z (fm)",200,-20.,20.,1000,-30.,70.);
  TH2D *hYZ = new TH2D("hYZ","YZ position;Y (fm);Z (fm)",200,-20.,20.,1000,-30.,70.);
  TH3F *hXYZ = new TH3F("hXYZ","XYZ position;Z (fm);X (fm);Y (fm)",1000,-10.,70.,200,-20.,20.,200,-20.,20.);
  TH1D *hPx = new TH1D("hPx","p_{x};p_{x} (GeV/c);Entries", 100, -10., 10.);
  TH1D *hPy = new TH1D("hPy","p_{y};p_{y} (GeV/c);Entries", 100, -10., 10.);
  TH1D *hPz = new TH1D("hPz","p_{z};p_{z} (GeV/c);Entries", 100, -10., 10.);
  TH1D *hEta = new TH1D("hEta","#eta;#eta;Entries", 200, -15., 15.);
  TH1D *hP = new TH1D("hP","p;p (GeV/c); Entries",200,0.,2.5);
  TH1D *hPt = new TH1D("hPt","p_{T};p_{T} (GeV/c);Entries",200,0.,2.5);
  TH1D *hPhi = new TH1D("hPhi","#phi;#phi (rad);Entries",200,-5.,5.);
  TH1D *hPDGmass = new TH1D("hPDGmass","Mass according to the PDG code;m^2 (GeV/c^2);Entries",200,-1.,2.);
  TH1D *hMass = new TH1D("hMass","Mass;m^2 (GeV/c^2);Entries",200,-1.,2.);
  TH1D *hEnergy = new TH1D("hEnergy","Energy;E (GeV);Entries",200,0.,10.);
  TH2D *hPtVsEta = new TH2D("hPtVsEta", "p_{T} vs. #eta of primary track;#eta;p_{T} (GeV/c)", 1000, -10, 10, 200, 0., 2.0);
  TH2D *hRapidityP = new TH2D("hRapidityP","Rapidity vs p;y;p (GeV/c)",2000,-2.,2.,1000,0.,2.);
  TH2D *hRapidityPt = new TH2D("hRapidityPt","Rapidity vs p_T;y;p_T (GeV/c)",1000,-10.,10.,200,0.,2.);
  TH1D *hRapidity = new TH1D("hRapidity","Rapidity;y;",1000,-10.,10.);
  TH2D *hBetaP = new TH2D("hBetaP","Beta vs p;#beta;p (GeV/c)",1000, 0., 2.2, 2000, 0., 3.);
  TH2D *hInvBetaP = new TH2D("hInvBetaP","1/#beta vs p;#beta;p (GeV/c)",2000, -2.2, 2.2, 2000, 0., 3.);

//lorentz
  TH2D *hRapidityP_LS = new TH2D("hRapidityP_LS","Rapidity vs p;y;p (GeV/c)",2000,-2.,2.,1000,0.,2.);
  TH2D *hRapidityPt_LS = new TH2D("hRapidityPt_LS","Rapidity vs p_T;y;p_T (GeV/c)",1000,-10.,10.,200,0.,2.);
  TH1D *hRapidity_LS = new TH1D("hRapidity_LS","Rapidity;y;",1000,-10.,10.);
  TH1D *hPz_LS = new TH1D("hPz_LS","p_{z};p_{z} (GeV/c);Entries", 100, -10., 10.);
  TH1D *hPt_LS = new TH1D("hPt_LS","p_{T};p_{T} (GeV/c);Entries",200,0.,2.5); 
  TH1D *hEta_LS = new TH1D("hEta_LS","#eta;#eta;Entries", 200, -15., 15.);
  TH1D *hP_LS = new TH1D("hP_LS","p;p (GeV/c); Entries",200,0.,2.5);

  TH2D *hRapidityP_LS_CMS = new TH2D("hRapidityP_LS_CMS","Rapidity vs p;y;p (GeV/c)",2000,-2.,2.,1000,0.,2.);
  TH2D *hRapidityPt_LS_CMS = new TH2D("hRapidityPt_LS_CMS","Rapidity vs p_T;y;p_T (GeV/c)",1000,-10.,10.,200,0.,2.);
  TH1D *hRapidity_LS_CMS = new TH1D("hRapidity_LS_CMS","Rapidity;y;",1000,-10.,10.);
  TH1D *hPz_LS_CMS = new TH1D("hPz_LS_CMS","p_{z};p_{z} (GeV/c);Entries", 100, -10., 10.);
  TH1D *hPt_LS_CMS = new TH1D("hPt_LS_CMS","p_{T};p_{T} (GeV/c);Entries",200,0.,2.5);
  TH1D *hEta_LS_CMS = new TH1D("hEta_LS_CMS","#eta;#eta;Entries", 200, -15., 15.);
  TH1D *hP_LS_CMS = new TH1D("hP_LS_CMS","p;p (GeV/c); Entries",200,0.,2.5);



  TH2D *hktPairCM = new TH2D("hktPairCM","Rapidity vs k_{T};y;k_{T} (GeV/c)",1000,-10.,10.,200,0.,2.);
  TH2D *hktPairLS = new TH2D("hktPairLS","Rapidity vs k_{T};y;k_{T} (GeV/c)",1000,-10.,10.,200,0.,2.);
 
  //after pion select
  TH1D *hChargePion = new TH1D("hChargePion","Charge of particles;Charge;Entries",5,-2,3);
  TH1D *hPDGPion = new TH1D("hPDGPion","PDG of particles;PDG;Entries",4000,-500.,3500.);
  TH1D *hXPion = new TH1D("hXPion","X position;X (fm);Entries",200,-20.,20.);
  TH1D *hYPion = new TH1D("hYPion","Y position;Y (fm);Entries",200,-20.,20.);
  TH1D *hZPion = new TH1D("hZPion","Z position;Z (fm);Entries",500,-60.,60.);//was -20,20
  TH1D *hPxPion = new TH1D("hPxPion","p_{x};p_{x} (GeV/c);Entries", 100, -10., 10.);
  TH1D *hPyPion = new TH1D("hPyPion","p_{y};p_{y} (GeV/c);Entries", 100, -10., 10.);
  TH1D *hPzPion = new TH1D("hPzPion","p_{z};p_{z} (GeV/c);Entries", 100, -10., 10.);
  TH1D *hEtaPion = new TH1D("hEtaPion","#eta;#eta;Entries", 200, -15., 15.);
  TH1D *hPDGmassPion = new TH1D("hPDGmassPion","Mass according to the PDG code;m^2 (GeV/c^2);Entries",200,-1.,2.);
  TH1D *hMassPion = new TH1D("hMassPion","Mass;m^2 (GeV/c^2);Entries",200,-1.,2.);
  TH1D *hEnergyPion = new TH1D("hEnergyPion","Energy;E (GeV);Entries",200,0.,10.);
  TH1D *hPhiPion = new TH1D("hPhiPion","#phi;#phi (rad);Entries",200,-5.,5.);
  TH1D *hPPion = new TH1D("hPPion","Momentum of #pi;p (GeV/c);Entries",200, 0., 2.5);
  TH1D *hPtPion = new TH1D("hPtPion","p_{T};p_{T} (GeV/c);Entries",200,0.,2.5);


  TH2D *hRapidityPionPlusP = new TH2D("hRapidityPionPlusP","Rapidity vs p pion plus;y;p (GeV/c)",2000,-4.,4.,1000,0.,2.);
  TH2D *hRapidityPionPlusPt = new TH2D("hRapidityPionPlusPt","Rapidity vs p_T pion plus;y;p_T (GeV/c)",2000,-4.,4.,1000,0.,2.);
  TH2D *hRapidityPionMinusP = new TH2D("hRapidityPionMinusP","Rapidity vs p pion minus;y;p (GeV/c)",2000,-4.,4.,1000,0.,2.);
  TH2D *hRapidityPionMinusPt = new TH2D("hRapidityPionMinusPt","Rapidity vs p_T pion minus;y;p_T (GeV/c)",2000,-4.,4.,1000,0.,2.);

  TH2D *hRTauPionPlus = new TH2D("hRTauPionPlus","hRTauPionPlus;r_{T} (fm);#tau (fm/c)",100,0.,30.,100,0.,30.);
  TH2D *hRTauPionMinus = new TH2D("hRTauPionMinus","hRTauPionMinus;r_{T} (fm);#tau (fm/c)",100,0.,30.,100,0.,30.);

  TH2D *hRTPionPlus = new TH2D("hRTPionPlus","hRTPionPlus;r_{T} (fm);t (fm/c)",100,0.,30.,100,0.,30.);
  TH2D *hRTPionMinus = new TH2D("hRTPionMinus","hRTPionMinus;r_{T} (fm);t (fm/c)",100,0.,30.,100,0.,30.);

  TH2D *hXYPionPlus = new TH2D("hXYPionPlus","XY position pion plus;X (fm);Y (fm)",200,-20.,20.,200,-20.,20.);
  TH2D *hXZPionPlus = new TH2D("hXZPionPlus","XZ position pion plus;X (fm);Z (fm)",200,-20.,20.,500,-20.,60.);
  TH2D *hYZPionPlus = new TH2D("hYZPionPlus","YZ position pion plus;Y (fm);Z (fm)",200,-20.,20.,500,-20.,60.);

  TH2D *hXYPionMinus = new TH2D("hXYPionMinus","XY position pion minus;X (fm);Y (fm)",200,-20.,20.,200,-20.,20.);
  TH2D *hXZPionMinus = new TH2D("hXZPionMinus","XZ position pion minus;X (fm);Z (fm)",200,-20.,20.,500,-20.,60.);
  TH2D *hYZPionMinus = new TH2D("hYZPionMinus","YZ position pion minus;Y (fm);Z (fm)",200,-20.,20.,500,-20.,60.);

  TH1D *hXPionPlus = new TH1D("hXPionPlus","X position;X (fm);Entries",200,-20.,20.);
  TH1D *hYPionPlus = new TH1D("hYPionPlus","Y position;Y (fm);Entries",200,-20.,20.);
  TH1D *hZPionPlus = new TH1D("hZPionPlus","Z position;Z (fm);Entries",200,-20.,20.);
  TH1D *hTPionPlus = new TH1D("hTPionPlus","t pion;t (fm/c);N",1500,0.,200.);
  TH1D *hTauPionPlus = new TH1D("hTauPionPlus","#tau pion;#tau (fm/c);N",1000,0.,100.);

  TH1D *hXPionMinus = new TH1D("hXPionMinus","X position;X (fm);Entries",200,-20.,20.);
  TH1D *hYPionMinus = new TH1D("hYPionMinus","Y position;Y (fm);Entries",200,-20.,20.);
  TH1D *hZPionMinus = new TH1D("hZPionMinus","Z position;Z (fm);Entries",200,-20.,20.);
  TH1D *hTPionMinus = new TH1D("hTPionMinus","t pion;t (fm/c);N",1500,0.,200.);
  TH1D *hTauPionMinus = new TH1D("hTauPionMinus","#tau pion;#tau (fm/c);N",1000,0.,100.);


  TH1D *hT = new TH1D("hT","t;t (fm/c);N",1500,0.,200.);
  TH1D *hNProtons200fm = new TH1D("hNProtons200fm",";;N",3,0,3);
  TH1D *hNNeutrons200fm = new TH1D("hNNeutrons200fm",";;N",3,0,3);

  TH1D *hNProtons0_199fm = new TH1D("hNProtons0_199fm",";;N",3,0,3);
  TH1D *hNNeutrons0_199fm = new TH1D("hNNeutrons0_199fm",";;N",3,0,3);

  TH1D *hNProtons = new TH1D("hNProtons",";;N",3,0,3);
  TH1D *hNNeutrons = new TH1D("hNNeutrons",";;N",3,0,3);

  TH1D *hNofParticlesIn1DeffMult = new TH1D("hNofParticlesIn1DeffMult",";;N",3,0,3);
  TH1D *hNofParticlesIn2DeffMult = new TH1D("hNofParticlesIn2DeffMult",";;N",3,0,3);

  TH1D *hNofParticlesIn1DeffMult_woProtons200fm = new TH1D("hNofParticlesIn1DeffMult_woProtons200fm",";;N",3,0,3);
  TH1D *hNofParticlesIn2DeffMult_woProtons200fm = new TH1D("hNofParticlesIn2DeffMult_woProtons200fm",";;N",3,0,3);

  TF1 *eff = new TF1("eff","[0]*exp( - TMath::Power( ([1]/x) , [2] ) )",0.,2.);
  double par0 = 0.766958;
  double par1 = 0.106865;
  double par2 = 2.89999;
  eff->SetParameters(par0, par1, par2);

  // for 2D efficiency
  TF1 *eff2D = new TF1("eff2D","[0]*exp( - TMath::Power( ([1]/x) , [2] ) )",0.,2.);
  vector <double> effPtEta_par0 = {0.209896, 0.572387, 0.67324, 0.722479, 0.781151, 0.826993, 0.854024, 0.86946, 0.887045, 0.888512, 0.887932, 0.895433, 0.89101, 0.890118, 0.890439, 0.895355, 0.892673, 0.89766, 0.905097, 0.901792, 0.906511, 0.905563, 0.907799, 0.904004, 0.914687, 0.907462, 0.906398, 0.906832, 0.913391, 0.906275, 0.90729, 0.909223, 0.899668, 0.900899, 0.900946, 0.896836, 0.905008, 0.900258, 0.895357, 0.891014, 0.885054, 0.892036, 0.874984, 0.879219, 0.849106};
  vector <double> effPtEta_par1 = {0.361997, 0.264772, 0.186217, 0.149996, 0.132803, 0.12689, 0.121994, 0.117229, 0.117883, 0.112186, 0.115001, 0.11486, 0.112434, 0.115041, 0.11263, 0.107745, 0.109371, 0.111086, 0.107953, 0.108031, 0.111217, 0.108695, 0.113001, 0.112749, 0.111384, 0.110075, 0.107814, 0.108848, 0.107642, 0.108907, 0.10976, 0.109058, 0.108591, 0.107734, 0.113506, 0.107831, 0.112075, 0.111591, 0.111581, 0.111599, 0.113503, 0.112723, 0.113835, 0.115637, 0.118168};
  vector <double> effPtEta_par2 = {2.34089, 2.22472, 2.3333, 4.69045, 4.72749, 4.19807, 4.43388, 5.00734, 4.27832, 5.18579, 4.37454, 4.33481, 5.15944, 4.80342, 6.12996, 6.42091, 5.65661, 4.7156, 4.33515, 5.13981, 4.59951, 5.01735, 5.79742, 5.06123, 4.32697, 4.96016, 5.92188, 6.12095, 4.56565, 5.97894, 5.12106, 6.62969, 6.95313, 7.00752, 5.78742, 8.18995, 4.88928, 6.3173, 5.669, 7.61135, 7.26456, 5.97931, 6.38419, 6.26772, 5.62442};
  vector <double> eta_bins = {-1.84, -1.8, -1.76, -1.72, -1.68, -1.64, -1.6, -1.56, -1.52, -1.48, -1.44, -1.4, -1.36, -1.32, -1.28, -1.24, -1.2, -1.16, -1.12, -1.08, -1.04, -1, -0.96, -0.92, -0.88, -0.84, -0.8, -0.76, -0.72, -0.68, -0.64, -0.6, -0.56, -0.52, -0.48, -0.44, -0.4, -0.36, -0.32, -0.28, -0.24, -0.2, -0.16, -0.12, -0.08};

  // check my calc 1D eff
  TH1D *h_1D_eff_func_reg_part = new TH1D("h_1D_eff_func_reg_part",";;N",200,0,2);
  TH1D *h_1D_eff_func_all_part = new TH1D("h_1D_eff_func_all_part",";;N",200,0,2);

  Int_t eventCounter = 0;
  Int_t hundredIter = 0;
  const double MassPion = 0.139570;


  // Run info
  hSqrtSnn->Fill( myReader->run()->nnSqrtS() );
  myReader->run()->print();

  TLorentzVector vecCurCM, vecLastCM, vecCurLS, vecLastLS;

  // Loop over events
  for(Long64_t iEvent=0 ; iEvent<events2read; iEvent++) {

      eventCounter++;
      if( eventCounter >= 1000 ) {
        eventCounter = 0;
        hundredIter++;
        std::cout << "Working on event #[" << (hundredIter * 1000)
      << "/" << events2read << "]" << std::endl;
      }

      Bool_t readEvent = myReader->loadEntry(iEvent);
      if( !readEvent ) {
        std::cout << "Something went wrong, Master! Nothing to analyze..."
      << std::endl;
        break;
      }

      // Retrieve femtoDst
      McDst *dst = myReader->mcDst();

      // Retrieve event information
      McEvent *event = dst->event();
      if( !event ) {
        std::cout << "Something went wrong, Master! Event is hiding from me..."
      << std::endl;
        break;
      }

      // Track analysis
      Int_t nTracks = dst->numberOfParticles();
      Int_t NumOfCharged = 0;
      Int_t mNumberOfTrackMult = 0;
      Int_t mNumberOfTrackMult2D = 0;
      Int_t mNumberOfTrackMultwoProtonsLess200fm = 0;
      Int_t mNumberOfTrackMult2DwoProtonsLess200fm = 0;	

      //if (event->impact()>4.8) continue;
      //if (event->impact()<10) continue;

      // Track loop
      for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

          // Retrieve i-th femto track
          McParticle *particle = dst->particle(iTrk);

          if (!particle) continue;
          //std::cout << "Track #[" << (iTrk+1) << "/" << nTracks << "]"  << std::endl;
          

      
          //variables
          Double_t x = particle->x();
          Double_t y = particle->y();
          Double_t z = particle->z();
          Double_t px = particle->px();
          Double_t py = particle->py();
          Double_t pz = particle->pz();
          Double_t p = particle->ptot();
          Double_t pt = particle->pt();
          Double_t phi = particle->phi();
          Double_t eta = particle->eta();
          Double_t mass = particle->mass();
          //Double_t pdgMass = particle->pdgMass();
          Double_t e = particle->e();
          Double_t t = particle->t();
          Int_t pdg = particle->pdg();


          vecCurCM.SetXYZM(px, py, pz, mass);

          hPDG->Fill( pdg );
          hX->Fill( x );
          hY->Fill( y );
          hZ->Fill( z );
          hXY->Fill( x, y );
          hXZ->Fill( x, z );
          hYZ->Fill( y, z );
          hXYZ->Fill( z, x, y );
          hPx->Fill( px );
          hPy->Fill( py );
          hPz->Fill( pz );
          hEta->Fill( eta );
          //hPDGmass->Fill( pdgMass );
          hMass->Fill( mass );
          hPhi->Fill( phi );
          hP->Fill( p );
          hPt->Fill( pt );
          hEnergy->Fill( e );
          
          TLorentzVector LorVec;
          LorVec.SetXYZM( px, py, pz, mass );

          //if (LorVec.Pt()<0.15 || LorVec.Pt()>1.5) continue;
          //if (LorVec.Rapidity()<0. || LorVec.Rapidity()>1.) continue;

          hRapidityP->Fill(LorVec.Rapidity(),LorVec.P());
          hRapidityPt->Fill(LorVec.Rapidity(),LorVec.Pt());
          hRapidity->Fill(LorVec.Rapidity());
          hBetaP->Fill(LorVec.P(),LorVec.Beta());
          double rap = LorVec.Rapidity();

          //Lorentz
          double MassNuclon = 0.931494; //GeV
          double sNN = 3; 
          double pBeam = TMath::Sqrt(sNN*sNN/4. - MassNuclon*MassNuclon);
          double beta = 2.*pBeam/sNN;
            
          TLorentzVector pcms; 
          pcms.SetXYZM(px, py, pz, mass);
            
          TVector3 boostVect; 
          boostVect.SetXYZ(0., 0., beta);
            
          TLorentzVector plab = pcms;
          plab.Boost(-boostVect);
          vecCurLS = plab;

          double etaLS = plab.PseudoRapidity();
        
          hRapidityP_LS->Fill(plab.Rapidity(),plab.P());
          hRapidityPt_LS->Fill(plab.Rapidity(),plab.Pt()); 
          hPz_LS->Fill(plab.Pz());
          hPt_LS->Fill(plab.Pt());
          hRapidity_LS->Fill(plab.Rapidity());
          hEta_LS->Fill(plab.PseudoRapidity());
          hP_LS->Fill(plab.Mag());
        
          /*plab.Boost(boostVect);

          hRapidityP_LS_CMS->Fill(plab.Rapidity(),plab.P());
          hRapidityPt_LS_CMS->Fill(plab.Rapidity(),plab.Pt());
          hPz_LS_CMS->Fill(plab.Pz());
          hPt_LS_CMS->Fill(plab.Pt());
          hRapidity_LS_CMS->Fill(plab.Rapidity());
          hEta_LS_CMS->Fill(plab.PseudoRapidity());
          hP_LS_CMS->Fill(plab.Mag());*/


          double ktPairCM = 0.5 * (vecCurCM + vecLastCM).Perp();
          double yPairCM = (vecCurCM + vecLastCM).Rapidity();
          double ktPairLS = 0.5 * (vecCurLS + vecLastLS).Perp();   
          double yPairLS = (vecCurLS + vecLastLS).Rapidity();      

          hktPairCM->Fill(yPairCM, ktPairCM);
          hktPairLS->Fill(yPairLS, ktPairLS);
          //Lorentz

          vecLastCM = vecCurCM;
          vecLastLS = vecCurLS;

          if (TMath::Abs(pdg) == 2212 && t>199) {hNProtons200fm->Fill(1);}
          if (TMath::Abs(pdg) == 2112 && t>199) {hNNeutrons200fm->Fill(1);}

          if (TMath::Abs(pdg) == 2212 && t<200) {hNProtons0_199fm->Fill(1);}
          if (TMath::Abs(pdg) == 2112 && t<200) {hNNeutrons0_199fm->Fill(1);}

          if (TMath::Abs(pdg) == 2212) {hNProtons->Fill(1);}
          if (TMath::Abs(pdg) == 2112) {hNNeutrons->Fill(1);}


          // Mult. calc. with efficiency
          double y_eff;
          y_eff = eff->Eval(pt);

          double min = 0;
          double max = 1;
          double frand = (double)rand() / RAND_MAX;
          //double r = min + frand * (max - min);
          double r = gRandom->Uniform(0, 1);
          

          double eta_ls = plab.PseudoRapidity();

          for (int i=0; i<=eta_bins.size(); i++)
          {
              if ( (eta_ls >= eta_bins[i]) && (eta_ls < eta_bins[i+1]) )
              eff2D->SetParameters(effPtEta_par0[i], effPtEta_par1[i], effPtEta_par2[i]);
              if (eta_ls < eta_bins[0]) eff2D->SetParameters(effPtEta_par0[0], effPtEta_par1[0], effPtEta_par2[0]); //std::cout<<"eta < "<<eta_bins[0]<<std::endl;
              if (eta_ls > eta_bins[eta_bins.size()]) eff2D->SetParameters(effPtEta_par0[eta_bins.size()], effPtEta_par1[eta_bins.size()], effPtEta_par2[eta_bins.size()]); //std::cout<<"eta > "<<eta_bins[eta_bins.size()] <<std::endl;
          }

          double y_eff2D;
          y_eff2D = eff2D->Eval(pt);

          if (etaLS>-2. && etaLS<0.) {hRap_LS_afterEtacut->Fill(plab.Rapidity()); hEta_LS_afterEtacut->Fill (etaLS);}
          if (LorVec.Rapidity()>-0.1 && LorVec.Rapidity()<0.) {hRap_LS_afterRapcut->Fill(plab.Rapidity()); hEta_LS_afterRapcut->Fill (etaLS);}
          if (etaLS>-2. && etaLS<0.) {if (LorVec.Rapidity()>-0.1 && LorVec.Rapidity()<0.) {hRap_LS_afterEtaRapcut->Fill(plab.Rapidity()); hEta_LS_afterEtaRapcut->Fill (etaLS);}}

          if (particle->charge() != 0  && etaLS>-2. && etaLS<0. && (r <= y_eff) ) {
              mNumberOfTrackMult++;
              //hEta_LS_aftercut->Fill (etaLS);
              hP_aftercut->Fill ( plab.Mag() );
              hPt_aftercut->Fill ( plab.Pt() );
              hNofParticlesIn1DeffMult->Fill(1);
              h_1D_eff_func_reg_part->Fill(plab.Pt());
              if ((TMath::Abs(pdg) == 2212) && (t<200)) {mNumberOfTrackMultwoProtonsLess200fm++; hNofParticlesIn1DeffMult_woProtons200fm->Fill(1);}
          } //1D eff case
          if (particle->charge() != 0  && etaLS>-2. && etaLS<0.) {
            h_1D_eff_func_all_part->Fill(plab.Pt());
        } //1D eff case


          if (particle->charge() != 0 && etaLS>-2. && etaLS<0. && (r <= y_eff2D) )  //2D eff case   
          {
              mNumberOfTrackMult2D++;  
              hNofParticlesIn2DeffMult->Fill(1);
              
	            if ((TMath::Abs(pdg) == 2212) && (t<200)) {mNumberOfTrackMult2DwoProtonsLess200fm++; hNofParticlesIn2DeffMult_woProtons200fm->Fill(1);}
          } //fxtmult 3gev ls
          if ( particle->charge() ) {
        
              NumOfCharged++;
              if (particle->charge()>0) {hCharge->Fill(1,1); hInvBetaP->Fill(TMath::Abs(LorVec.P()),1./(LorVec.Beta()));}
              if (particle->charge()<0) {hCharge->Fill(-1,1); hInvBetaP->Fill((-1.)*LorVec.P(),1./(LorVec.Beta()));}
              hPtVsEta->Fill( eta, pt );
	        }	
	        if ( pdg == 211 || pdg == -211 && p>0.15 && p<1.5 && pt>0.15 && pt<1.5) {
              hPDGPion->Fill( pdg );
              hXPion->Fill( x );
              hYPion->Fill( y );
              hZPion->Fill( z );
              hPxPion->Fill( px );
              hPyPion->Fill( py );
              hPzPion->Fill( pz );
              hEtaPion->Fill( eta );
              //hPDGmassPion->Fill( pdgMass );
              hMassPion->Fill( mass );
              hPhiPion->Fill( phi );
              hPPion->Fill( p );
              hPtPion->Fill( pt );
              hEnergyPion->Fill( e );}
	        if (pdg == 211 && p>0.15 && p<1.5 && pt>0.15 && pt<1.5) {
		          hChargePion->Fill(1,1);
              hRapidityPionPlusP->Fill( LorVec.Rapidity(),LorVec.P() );
              hRapidityPionPlusPt->Fill( LorVec.Rapidity(),LorVec.Pt() );
              hRTauPionPlus->Fill(TMath::Sqrt( x*x + y*y ), TMath::Sqrt( t*t - z*z ) );
              hRTPionPlus->Fill(TMath::Sqrt( x*x + y*y ), t);
              hXYPionPlus->Fill( x, y );
              hXZPionPlus->Fill( x, z );
              hYZPionPlus->Fill( y, z );
              hTPionPlus->Fill( t );
              hXPionPlus->Fill( x );
              hYPionPlus->Fill( y );
              hZPionPlus->Fill( z );
		          hTauPionPlus->Fill( TMath::Sqrt( t*t - z*z ) );
		      }

        if (pdg == -211 && p>0.15 && p<1.5 && pt>0.15 && pt<1.5) {
              hChargePion->Fill(-1,1);
              hRapidityPionMinusP->Fill( LorVec.Rapidity(),LorVec.P() );
              hRapidityPionMinusPt->Fill( LorVec.Rapidity(),LorVec.Pt() );
              hRTauPionMinus->Fill(TMath::Sqrt( x*x + y*y ), TMath::Sqrt( t*t - z*z ) );
              hRTPionMinus->Fill(TMath::Sqrt( x*x + y*y ), t);
              hXYPionMinus->Fill( x, y );
              hXZPionMinus->Fill( x, z );
              hYZPionMinus->Fill( y, z );
              hTPionMinus->Fill( t );
              hXPionMinus->Fill( x );
              hYPionMinus->Fill( y );
              hZPionMinus->Fill( z );
		          hTauPionMinus->Fill( TMath::Sqrt( t*t - z*z ) );
		        }

	      double imp = event->impact();



    } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    //std::cout<<"iEvent = "<<iEvent<<std::endl;
    //std::cout<<"mNumberOfTrackMult2D = "<<mNumberOfTrackMult2D<<std::endl;

    hImpactParVsNch->Fill( NumOfCharged, event->impact() );
    hNch->Fill( NumOfCharged );
    hImpact->Fill( event->impact() );
    hMult->Fill( mNumberOfTrackMult );
    hMultwoProtonsLess200fm->Fill( mNumberOfTrackMultwoProtonsLess200fm );
    hMultEff2D->Fill( mNumberOfTrackMult2D );
    hMultEff2DwoProtonsLess200fm->Fill( mNumberOfTrackMult2DwoProtonsLess200fm );
  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)


  oFile->Write();
  oFile->Close();

  myReader->Finish();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"<< std::endl;
}
