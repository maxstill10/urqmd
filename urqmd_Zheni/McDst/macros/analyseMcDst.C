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
#include "TComplex.h"
#include "TProfile.h"

// McDst headers
#include "../McDstReader.h"
#include "../McDst.h"
#include "../McEvent.h"
#include "../McParticle.h"
#include "../McRun.h"

//My functions
int GetCentrality(int refMult);
int Get_refmult(McDst *dst);
double Cacl_2cov_ref(TComplex Qn, int M);
double Cacl_4cov_ref(TComplex Qn, TComplex Q2n, int M);
double Cacl_2cov_dif(TComplex Qn, TComplex pn, int M, int mp, int mq);
double Cacl_4cov_dif(TComplex Qn, TComplex pn, TComplex qn, TComplex Q2n, TComplex q2n, int M, int mp, int mq);

//pt and eta intervals
//double const pt_intervals[] = {0.15, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.4, }


// inFile - is a name of name.uDst.root file or a name
//          of a name.lis(t) files that contains a list of
//          name1.uDst.root files
//_________________
void analyseMcDst(const Char_t *inFile,
		  const Char_t *oFileName) {
  
  gSystem->Load("/star/u/mmorozov/urqmd/urqmd_Zheni/McDst/libMcDst.so");

  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

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

  //Hist initialization
  TH1F *hpT = new TH1F("hpT", "p_{T} of particles", 570, 0.15, 3);

  TProfile *pCor2_ref[3];
  TProfile *pCor4_ref[3];

  for(int n=0; n!=3; n++){
    pCor2_ref[n] = new TProfile(Form("pCor2_ref_%i", n), "", 9, 0, 9);
    pCor4_ref[n] = new TProfile(Form("pCor4_ref_%i", n), "", 9, 0, 9);
  }

  Int_t eventCounter = 0;
  Int_t hundredIter = 0;
  const double MassPion = 0.139570;


  // Run info
  myReader->run()->print();


  // Loop over events
  for(Long64_t iEvent=0 ; iEvent<events2read; iEvent++) {

    eventCounter++;
//std::cout << "Working on event #[" << eventCounter
//    << "/" << events2read << "]" << std::endl;
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

    //Get event centrality
    int refMult = Get_refmult(dst) ;
    int cent = GetCentrality(refMult);
    if(cent<0) continue;
    

    //Variables for particle number counting
    double lQn_calc[15] = {}; // (cos(nphi), sin(nphi), cos(2nphi), sin(2nphi), numPart) * (n=1,2,3)
    double lpn_pos_cent_calc[15] = {};
    double lpn_neg_cent_calc[15] = {};

    // Track analysis
    Int_t nTracks = dst->numberOfParticles();

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

      if(particle->charge() == 0) continue;
      if(particle->pt() < 0.15) continue;

      hpT->Fill(pt);

      //Get pt and eta bins
      //int ipt = 

      //Directed flow measurements
      if(fabs(particle->eta()) > 5.5) continue;

      lQn_calc[0] += TMath::Cos(phi);
      lQn_calc[1] += TMath::Sin(phi);
      lQn_calc[2] += TMath::Cos(2*phi);
      lQn_calc[3] += TMath::Sin(2*phi);
      lQn_calc[4] += 1;

      //Elliptic and triangular flow measurements
      if(fabs(particle->eta()) > 1.) continue;

      for(int n = 1; n!=3; n++){
        lQn_calc[5*n] += TMath::Cos((n+1)*phi);
        lQn_calc[5*n+1] += TMath::Sin((n+1)*phi);
        lQn_calc[5*n+2] += TMath::Cos(2*(n+1)*phi);
        lQn_calc[5*n+3] += TMath::Sin(2*(n+1)*phi);
        lQn_calc[5*n+4] += 1;
      }

    }//for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    
    //Flow calculations by cumulants
    if(cent<3) continue;
    //Reference flow calculation
    for(int n=0; n!=3; n++){
      TComplex Qn(lQn_calc[5*n], lQn_calc[5*n+1]);
      TComplex Q2n(lQn_calc[5*n+2], lQn_calc[5*n+3]);
      int M = int(lQn_calc[5*n+4]);

      double ev_2cor_ref = Cacl_2cov_ref(Qn, M);
      double ev_4cor_ref = Cacl_4cov_ref(Qn, Q2n, M);

      pCor2_ref[n]->Fill(cent, ev_2cor_ref, M*(M-1));
      pCor4_ref[n]->Fill(cent, ev_4cor_ref, M*(M-1)*(M-2)*(M-3));
    }
  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)


  oFile->Write();
  oFile->Close();

  myReader->Finish();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"<< std::endl;
}


int GetCentrality(int refMult){
  if(refMult>48) return -1;
  else if(refMult>=19) return 8; //0-5%
  else if(refMult>=16) return 7; //5-10%
  else if(refMult>=12) return 6; //10-20%
  else if(refMult>=9)  return 5; //20-30%
  else if(refMult>=7)  return 4; //30-40%
  else if(refMult>=5)  return 3; //40-50%
  else if(refMult>=4)  return 2; //50-60%
  else if(refMult>=3)  return 1; //60-70%
  else if(refMult>=2)  return 0; //70-80%
  else                 return -1;
}


int Get_refmult(McDst *dst){
  int refMult = 0;
  
  // Track analysis
  Int_t nTracks = dst->numberOfParticles();

  // Track loop
  for(Int_t iTrk=0; iTrk<nTracks; iTrk++){
    // Retrieve i-th femto track
    McParticle *particle = dst->particle(iTrk);

    if (!particle) continue;

    if(particle->charge() == 0) continue;
    if(particle->pt() < 0.15) continue;
    if(fabs(particle->eta()) > 1.) continue;

    refMult++;
  }

  return refMult;
}





double Cacl_2cov_ref(TComplex Qn, int M){
  return (Qn.Rho2() - M) / M / (M-1);
}


double Cacl_4cov_ref(TComplex Qn, TComplex Q2n, int M){
  TComplex Qn_star = TComplex::Conjugate(Qn);
  double Qn_squared = Qn.Rho2();
  double up_part = Qn_squared*Qn_squared + Q2n.Rho2() - 2*(Q2n*Qn_star*Qn_star).Re() - 4*(M-2)*Qn_squared + 2*M*(M-3);
  return up_part / M / (M-1) / (M-2) / (M-3);
}


double Cacl_2cov_dif(TComplex Qn, TComplex pn, int M, int mp, int mq){
  TComplex Qn_star = TComplex::Conjugate(Qn);
  double up_part = (pn*Qn_star).Re() - mq;
  return up_part / (mp*M-mq);
}


double Cacl_4cov_dif(TComplex Qn, TComplex pn, TComplex qn, TComplex Q2n, TComplex q2n, int M, int mp, int mq){
  TComplex Qn_star = TComplex::Conjugate(Qn);
  TComplex qn_star = TComplex::Conjugate(qn);
  TComplex Q2n_star = TComplex::Conjugate(Q2n);
  double Qn_squared = Qn.Rho2();
  double up_part = (pn*Qn_star*Qn_squared - q2n*Qn_star*Qn_star - pn*Qn_star*Q2n_star - 2*M*pn*Qn_star - 2*mq*Qn_squared +\
                    7*qn*Qn_star - Qn*qn_star + q2n*Q2n_star + 2*pn*Qn_star + 2*mq*M - 6*mq).Re();
  return up_part / (mp*M - 3*mq) / (M-1) / (M-2);
}

