/**
 * Code is developed to build charged pion, kaons, proton and antiproton
 * spectra from mcDst files. It is designed to work in a standalone mode
 * to compile the code one should use prepared Makefile with the command:
 * ```
 * make 
 * '''
 * Do not forget to adjust path to ROOT and McDst packages, libraries
 * and set environment variables (including McDst).
 * 
 * To run the code simply run next commands from the terminal:
 * ```
 * spectraFromMcDst inputFile outputFile.root
 * '''
 * For the input file one can use either fname.mcDst.root file or
 * a list of mcDst files with the .list or .lis extention.
 */

// C++ headers
#include <iostream>
#include <cmath>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TDatabasePDG.h"

// McDst headers
#include "McDstReader.h"
#include "McDst.h"
#include "McEvent.h"
#include "McParticle.h"
#include "McRun.h"

// Array of momentum spectra bin edges
const Int_t nBinsPtSpectra = 59;
Double_t ptBinEdges[nBinsPtSpectra + 1] =  { 
  0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 
  0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 
  1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 
  1.8, 1.85, 1.9, 1.95, 2.0, 2.1, 2.2, 2.3, 2.4, 2.6, 2.8, 3.0, 3.5, 
  4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0 
};

//________________
Int_t findPtBin(Double_t pt) {
  Int_t bin = 0;
  for (Int_t iBin=0; iBin<nBinsPtSpectra; iBin++) {
    if ( ptBinEdges[iBin] <= pt && pt < ptBinEdges[iBin+1] ) {
      bin = iBin;
      break;
    }
  }
  return bin;
}

//________________
Double_t ptBinWidth(Int_t iBin) {
  return ptBinEdges[iBin+1] - ptBinEdges[iBin];
}

//________________
int main(int argc, char* argv[]) {

// #if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
//   R__LOAD_LIBRARY(libEG)
//   R__LOAD_LIBRARY(libMcDst)
// #else
//   gSystem->Load("libEG");
//   gSystem->Load("libMcDst.so");
// #endif

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  const char* fileName;
  const char* oFileName;

  // Next line is important for the spectra analysis because it
  // will be used to select |rapidity|<rapidityCut
  Double_t rapidityCut = 0.1; 

  switch (argc) {
  case 3:
    fileName = argv[1];
    oFileName = argv[2];
    break;
  default:
    std::cout << "Usage: processMcDstStandalone inputFileName outputFileName.root" << std::endl;
    return -1;
  }
  std::cout << " inputFileName : " << fileName << std::endl;
  std::cout << " outputFileName: " << oFileName << std::endl;

  McDstReader* myReader = new McDstReader(fileName);
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

  TFile *oFile = new TFile(oFileName, "recreate");

  /////////////////////
  //  Histogramming  //
  /////////////////////

  //
  // Event
  //
  int mMult[2] = {-100, 100};
  float mImp[2] = {-100., 100.};

  TH1D *hImpactParameter = new TH1D("hImpactParameter","Impact parameter;b (bm);Entries", 100, 0., 16.);
  TH1D *hRefMult = new TH1D("hRefMult","Reference multiplicity (|#eta|<0.5, p_{T}>0.15 GeV/c); refMult (|#eta|<0.5, p_{T}>0.15 GeV/c);Entries",
                            600, 0, 600.);
  TH2D *hImpactParVsRefMult = new TH2D("hImpactParVsRefMult",
    "Impact parameter vs. refMult (|#eta|<0.5, p_{T}>0.15 GeV/c);refMult (|#eta|<0.5, p_{T}>0.15 GeV/c);Impact parameter (fm)",
    600, 0., 600., 100, 0., 16.);

  // Spectra of charged pion, kaons and protons
  // Charge: 0 - positive, 1 - negative
  TH1D *hPionSpectra_ch[2];
  TH1D *hKaonSpectra_ch[2];
  TH1D *hProtonSpectra_ch[2];

  TH1D *hPionSpectra_ch_mult_[2][9];
  TH1D *hKaonSpectra_ch_mult_[2][9];
  TH1D *hProtonSpectra_ch_mult_[2][9];
  TH1D *hPionSpectra_ch_imp_[2][9];
  TH1D *hKaonSpectra_ch_imp_[2][9];
  TH1D *hProtonSpectra_ch_imp_[2][9];

  TH1D *hPionEta_ch[2];
  TH1D *hPionRapidity_ch[2];

  //norm
  TH1D *hNorm_mult_[9];
  TH1D *hNorm_imp_[9];

  for (Int_t iMult=0; iMult<9; iMult++) {
  	hNorm_mult_[iMult] = new TH1D(Form("hNorm_mult_%d", iMult),
                                     Form("Reference multiplicity (|#eta|<0.5, p_{T}>0.15 GeV/c), mult %d;refMult (|#eta|<0.5, p_{T}>0.15 GeV/c);Entries",iMult),
                                     600, 0., 600.);
	hNorm_imp_[iMult] = new TH1D(Form("hNorm_imp_%d", iMult),
                                     Form("Impact parameter, imp %d;b (bm);Entries",iMult),
                                     100, 0., 16.);
  }

  for (Int_t iCharge=0; iCharge<2; iCharge++) {
    int charge = 0;
    charge = (iCharge == 0) ? 1 : -1;
    //example of wide bins in the end of histogram
    /*hPionSpectra_ch[iCharge] = new TH1D(Form("hPionSpectra_ch%d", iCharge),
                                     Form("p_{T}(#pi) of charge %d;p_{T} GeV/c;#frac{d^{2}N}{2 #pi p_{T} dp_{T} dy}", charge),
                                     nBinsPtSpectra, ptBinEdges);
    hKaonSpectra_ch[iCharge] = new TH1D(Form("hKaonSpectra_ch%d", iCharge),
                                     Form("p_{T}(K) of charge %d;p_{T} GeV/c;#frac{d^{2}N}{2 #pi p_{T} dp_{T} dy}", charge),
                                     nBinsPtSpectra, ptBinEdges);
    hProtonSpectra_ch[iCharge] = new TH1D(Form("hProtonSpectra_ch%d", iCharge),
                                     Form("p_{T}(p) of charge %d;p_{T} GeV/c;#frac{d^{2}N}{2 #pi p_{T} dp_{T} dy}", charge),
                                     nBinsPtSpectra, ptBinEdges);*/

    //bin = 50 MeV
    hPionSpectra_ch[iCharge] = new TH1D(Form("hPionSpectra_ch%d", iCharge),
                                     Form("p_{T}(#pi) of charge %d;p_{T} GeV/c;#frac{d^{2}N}{2 #pi p_{T} dp_{T} dy}", charge),
                                     180, 0., 9.);
    hKaonSpectra_ch[iCharge] = new TH1D(Form("hKaonSpectra_ch%d", iCharge),
                                     Form("p_{T}(K) of charge %d;p_{T} GeV/c;#frac{d^{2}N}{2 #pi p_{T} dp_{T} dy}", charge),
                                     180, 0., 9.);
    hProtonSpectra_ch[iCharge] = new TH1D(Form("hProtonSpectra_ch%d", iCharge),
                                     Form("p_{T}(p) of charge %d;p_{T} GeV/c;#frac{d^{2}N}{2 #pi p_{T} dp_{T} dy}", charge),
                                     180, 0., 9.);

    hPionEta_ch[iCharge] = new TH1D(Form("hPionEta_ch%d", iCharge),
                                     Form("#eta (#pi) of charge %d;#eta;N", charge),
                                     100, -5, 5);
    hPionRapidity_ch[iCharge] = new TH1D(Form("hPionRapidity_ch%d", iCharge),
                                     Form("y (#pi) of charge %d;y;N", charge),
                                     100, -5, 5); 
  
  }

  for (Int_t iMult=0; iMult<9; iMult++) {
   for (Int_t iCharge=0; iCharge<2; iCharge++) {
    int charge = 0;
    charge = (iCharge == 0) ? 1 : -1;
    /*hPionSpectra_ch_mult_[iCharge][iMult] = new TH1D(Form("hPionSpectra_ch_mult_%d%d", iCharge, iMult),
                                     Form("p_{T}(#pi) of charge %d, mult %d;p_{T} GeV/c;#frac{d^{2}N}{2 #pi p_{T} dp_{T} dy}", charge, iMult),
                                     nBinsPtSpectra, ptBinEdges);
    hKaonSpectra_ch_mult_[iCharge][iMult]  = new TH1D(Form("hKaonSpectra_ch_mult_%d%d", iCharge, iMult),
                                     Form("p_{T}(K) of charge %d, mult %d;p_{T} GeV/c;#frac{d^{2}N}{2 #pi p_{T} dp_{T} dy}", charge, iMult),
                                     nBinsPtSpectra, ptBinEdges);
    hProtonSpectra_ch_mult_[iCharge][iMult]  = new TH1D(Form("hProtonSpectra_ch_mult_%d%d", iCharge, iMult),
                                     Form("p_{T}(p) of charge %d, mult %d;p_{T} GeV/c;#frac{d^{2}N}{2 #pi p_{T} dp_{T} dy}", charge, iMult),
                                     nBinsPtSpectra, ptBinEdges);*/
    hPionSpectra_ch_mult_[iCharge][iMult] = new TH1D(Form("hPionSpectra_ch_mult_%d%d", iCharge, iMult),
                                     Form("p_{T}(#pi) of charge %d, mult %d;p_{T} GeV/c;#frac{d^{2}N}{2 #pi p_{T} dp_{T} dy}", charge, iMult),
                                     180, 0., 9.);
    hKaonSpectra_ch_mult_[iCharge][iMult]  = new TH1D(Form("hKaonSpectra_ch_mult_%d%d", iCharge, iMult),
                                     Form("p_{T}(K) of charge %d, mult %d;p_{T} GeV/c;#frac{d^{2}N}{2 #pi p_{T} dp_{T} dy}", charge, iMult),
                                     180, 0., 9.);
    hProtonSpectra_ch_mult_[iCharge][iMult]  = new TH1D(Form("hProtonSpectra_ch_mult_%d%d", iCharge, iMult),
                                     Form("p_{T}(p) of charge %d, mult %d;p_{T} GeV/c;#frac{d^{2}N}{2 #pi p_{T} dp_{T} dy}", charge, iMult),
                                     180, 0., 9.);

    hPionSpectra_ch_imp_[iCharge][iMult] = new TH1D(Form("hPionSpectra_ch_imp_%d%d", iCharge, iMult),
                                     Form("p_{T}(#pi) of charge %d, imp %d;p_{T} GeV/c;#frac{d^{2}N}{2 #pi p_{T} dp_{T} dy}", charge, iMult),
                                     180, 0., 9.);
    hKaonSpectra_ch_imp_[iCharge][iMult]  = new TH1D(Form("hKaonSpectra_ch_imp_%d%d", iCharge, iMult),
                                     Form("p_{T}(K) of charge %d, imp %d;p_{T} GeV/c;#frac{d^{2}N}{2 #pi p_{T} dp_{T} dy}", charge, iMult),
                                     180, 0., 9.);
    hProtonSpectra_ch_imp_[iCharge][iMult]  = new TH1D(Form("hProtonSpectra_ch_imp_%d%d", iCharge, iMult),
                                     Form("p_{T}(p) of charge %d, imp %d;p_{T} GeV/c;#frac{d^{2}N}{2 #pi p_{T} dp_{T} dy}", charge, iMult),
                                     180, 0., 9.);

   }
  }


  /////////////////////
  //     Analysis    //
  /////////////////////

  Int_t eventCounter = 0;
  Int_t clock = 10000;
  Int_t clockIter = 0;

  // Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

    eventCounter++;
    if( eventCounter >= clock ) {
      eventCounter = 0;
      clockIter++;
      std::cout << "Working on event #[" << (clockIter * clock)
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

    Double_t b = event->b();
    hImpactParameter->Fill( b );

    //
    // Particle analysis
    //

    // Retrieve number of particles in the event
    Int_t nParticles = dst->numberOfParticles();
    // Number of charged particles with |eta|<0.5 and pT>0.15 GeV/c
    Int_t refMult{0};
    Int_t refMult_calc{0};

    /////////////////for mult analysis//////////////////
     for(Int_t iTrk=0; iTrk<nParticles; iTrk++) {
      McParticle *particle = (McParticle*)dst->particle(iTrk);
      if (!particle) continue;
      Int_t pdgCode = particle->pdg();
      Double_t pt = particle->pt();
      Double_t eta = particle->eta();
      Int_t charge = particle->charge();
       if ( /*(TMath::Abs(pdgCode) == 211 || TMath::Abs(pdgCode) == 321 || TMath::Abs(pdgCode) == 2212)*/ charge != 0  && pt>0.20 && eta>-2 && eta<0 ) {
          refMult_calc++;
      	}

     }
     ///////////////////////////////////////////////////

    // Track loop
    for(Int_t iTrk=0; iTrk<nParticles; iTrk++) {

      // Retrieve i-th femto track
      McParticle *particle = (McParticle*)dst->particle(iTrk);

      if (!particle) continue;

      // std::cout << "Track #[" << (iTrk+1) << "/" << nParticles << "]"  
      // 		<< std::endl;

      /*
      std::cout << "pdgId: " << particle->pdg() 
             		<< " status: " << particle->status()
       		      << " px/py/pz/E/m/pt/y: " 
       		      << particle->px() << " / " 
       		      << particle->py() << " / "
       		      << particle->pz() << " / " 
       		      << particle->e() << " / "
       		      << particle->mass() << " / "
                << particle->pt() << " / "
                << particle->momentum().Rapidity() << " / "
                << std::endl;
            */

      Int_t pdgCode = particle->pdg();
      Double_t pt = particle->pt();
      Double_t eta = particle->eta();
      Int_t charge = particle->charge();
      Double_t rapidity = particle->momentum().Rapidity();

      

      Int_t iCharge{0};
      if (charge > 0) {
	      iCharge = 0;
	}
      else {
	      iCharge = 1;
	}

	if ( TMath::Abs(pdgCode) == 211) {hPionEta_ch[iCharge]->Fill(eta); hPionRapidity_ch[iCharge]->Fill(rapidity);}	

      // Calculate reference multiplicity
      if ( /*(TMath::Abs(pdgCode) == 211 || TMath::Abs(pdgCode) == 321 || TMath::Abs(pdgCode) == 2212)*/ charge != 0 && pt>0.20 && eta>-2 && eta<0  ) {
          refMult++;
      } // if ( charge != 0 && pt>0.15 && TMath::Abs(eta) < 0.5 )

      if (TMath::Abs(pdgCode) != 211 && TMath::Abs(pdgCode) != 321 && TMath::Abs(pdgCode) != 2212)   continue;

      //if ( TMath::Abs(rapidity)>rapidityCut ) continue;
      if ( pt<0.2 ) continue;
      if ( eta<-2 && eta>0  ) continue; 


      Int_t ptBin = findPtBin( pt );
      //Double_t deltaPt = ptBinWidth( ptBin );
      Double_t deltaPt = 0.05;
      Double_t deltaEta = 2.;
      // Note that 2*rapidity is on purpose, not a mistake
      //Double_t weight = 1. / (2 * TMath::Pi() * pt * (2. * rapidityCut) * deltaPt);
	Double_t weight = 1. / (2 * TMath::Pi() * pt * deltaEta * deltaPt);

      //std::cout << "weight: " << weight << std::endl;
      switch (TMath::Abs(pdgCode)) {
      case 211:
        hPionSpectra_ch[iCharge]->Fill( pt,  weight); 
        break;
      case 321:
        hKaonSpectra_ch[iCharge]->Fill( pt, weight ); 
        break;
      case 2212:
        hProtonSpectra_ch[iCharge]->Fill( pt, weight ); 
        break;
      default:
        break;
      }

      //mult
      for (int iMult=0; iMult < 9; iMult++) {
	
	switch ( iMult ) { //3 gev as in experiment
        case 0: mMult[0] = 139;  mMult[1] = 191; break; //0-5%
        case 1: mMult[0] = 119;  mMult[1] = 139; break; //5-10%
	case 2: mMult[0] = 88;  mMult[1] = 119; break; //10-20%
	case 3: mMult[0] = 65;  mMult[1] = 88; break; //20-30%
	case 4: mMult[0] = 47;  mMult[1] = 65; break;  //30-40%
	case 5: mMult[0] = 33;  mMult[1] = 47; break;   //40-50%
	case 6: mMult[0] = 23;  mMult[1] = 33; break;   //50-60%
	case 7: mMult[0] = 16;  mMult[1] = 23; break;   //60-70%
	case 8: mMult[0] = 11;  mMult[1] = 16; break;    //70-80%
	default: mMult[0] = 0.; mMult[1] = 600; break;
       };



      switch (TMath::Abs(pdgCode)) {
      case 211:
	if (refMult_calc >= mMult[0] && refMult_calc < mMult[1]) {hPionSpectra_ch_mult_[iCharge][iMult]->Fill( pt,  weight);}
        break;
      case 321:
	if (refMult_calc >= mMult[0] && refMult_calc < mMult[1]) {hKaonSpectra_ch_mult_[iCharge][iMult]->Fill( pt, weight );}
        break;
      case 2212:
	if (refMult_calc >= mMult[0] && refMult_calc < mMult[1]) {hProtonSpectra_ch_mult_[iCharge][iMult]->Fill( pt, weight );}
        break;
      default:
        break;
      }

	} //iMult

	// impact
 	for (int imp=0; imp < 9; imp++) {
	
	switch ( imp ) { //5.2 GeV
        case 0: mImp[0] = 0.;  mImp[1] = 3.3; break; //0-5%
	case 1: mImp[0] = 3.3;  mImp[1] = 4.7; break;//5-10%
	case 2: mImp[0] = 4.7;  mImp[1] = 6.6; break; //10-20%
	case 3: mImp[0] = 6.6;  mImp[1] = 8.1; break; //20-30%
	case 4: mImp[0] = 8.1;  mImp[1] = 9.3; break; //30-40%
	case 5: mImp[0] = 9.3;  mImp[1] = 10.4; break; //40-50%
	case 6: mImp[0] = 10.4;  mImp[1] = 11.4; break;//50-60%
	case 7: mImp[0] = 11.4;  mImp[1] = 12.3; break;//60-70%
	case 8: mImp[0] = 12.3;  mImp[1] = 13.2; break;//70-80%
	default: mImp[0] = 0.; mImp[1] = 600; break;
       };

	switch (TMath::Abs(pdgCode)) {
      case 211:
	if (b >= mImp[0] && b < mImp[1]) {hPionSpectra_ch_imp_[iCharge][imp]->Fill( pt,  weight);}
        break;
      case 321:
	if (b >= mImp[0] && b < mImp[1]) {hKaonSpectra_ch_imp_[iCharge][imp]->Fill( pt, weight );}
        break;
      case 2212:
	if (b >= mImp[0] && b < mImp[1]) {hProtonSpectra_ch_imp_[iCharge][imp]->Fill( pt, weight );}
        break;
      default:
        break;
      }

     
	}//impact

    } //for(Int_t iTrk=0; iTrk<nParticles; iTrk++)

    hRefMult->Fill(refMult);
    hImpactParVsRefMult->Fill(refMult, b);

  //////////////////////////////////for norm coefficient///////////////////////////////////////////

  for (Int_t iMult=0; iMult<9; iMult++) {
	switch ( iMult ) { //3 gev as in experiment
        case 0: mMult[0] = 139;  mMult[1] = 191; break; //0-5%
        case 1: mMult[0] = 119;  mMult[1] = 139; break; //5-10%
	case 2: mMult[0] = 88;  mMult[1] = 119; break; //10-20%
	case 3: mMult[0] = 65;  mMult[1] = 88; break; //20-30%
	case 4: mMult[0] = 47;  mMult[1] = 65; break;  //30-40%
	case 5: mMult[0] = 33;  mMult[1] = 47; break;   //40-50%
	case 6: mMult[0] = 23;  mMult[1] = 33; break;   //50-60%
	case 7: mMult[0] = 16;  mMult[1] = 23; break;   //60-70%
	case 8: mMult[0] = 11;  mMult[1] = 16; break;    //70-80%
	default: mMult[0] = 0.; mMult[1] = 600; break;

       };


  if (refMult >= mMult[0] && refMult < mMult[1]) {hNorm_mult_[iMult]->Fill(refMult);}

  }//iMult
    for (int imp=0; imp < 9; imp++) {
	
	switch ( imp ) { //5.2 GeV
	case 0: mImp[0] = 0.;  mImp[1] = 3.3; break; //0-5%
	case 1: mImp[0] = 3.3;  mImp[1] = 4.7; break;//5-10%
	case 2: mImp[0] = 4.7;  mImp[1] = 6.6; break; //10-20%
	case 3: mImp[0] = 6.6;  mImp[1] = 8.1; break; //20-30%
	case 4: mImp[0] = 8.1;  mImp[1] = 9.3; break; //30-40%
	case 5: mImp[0] = 9.3;  mImp[1] = 10.4; break; //40-50%
	case 6: mImp[0] = 10.4;  mImp[1] = 11.4; break;//50-60%
	case 7: mImp[0] = 11.4;  mImp[1] = 12.3; break;//60-70%
	case 8: mImp[0] = 12.3;  mImp[1] = 13.2; break;//70-80%
	default: mImp[0] = 0.; mImp[1] = 600; break;
       };
    if (b >= mImp[0] && b < mImp[1]) {hNorm_imp_[imp]->Fill(b);}

    }//imp
  /////////////////////////////////for norm coefficient////////////////////////////////////////////

  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  oFile->Write();
  oFile->Close();
  myReader->Finish();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!" << std::endl;

  return 0;
}
