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
#include "/star/u/mmorozov/urqmd/urqmd_Zheni/McDst/McDstReader.h"
#include "/star/u/mmorozov/urqmd/urqmd_Zheni/McDst/McDst.h"
#include "/star/u/mmorozov/urqmd/urqmd_Zheni/McDst/McEvent.h"
#include "/star/u/mmorozov/urqmd/urqmd_Zheni/McDst/McParticle.h"
#include "/star/u/mmorozov/urqmd/urqmd_Zheni/McDst/McRun.h"

//My functions
int GetCentrality(int refMult);
int Get_refmult(McDst *dst);
void Qn_calc(double lQn_calc[], int nOrd, double phi);
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
  TProfile *pCov24_ref[3];
  //diff corr
  TProfile *pCor2_dif_cent[3][2]; //(v1+v2+v3)x(pos, neg)
  TProfile *pCor4_dif_cent[3][2];
  TProfile *pCov_dif_cent[3][2][5];//(v1+v2+v3)x(pos, neg)x(22dif, 24dif, 42dif, 44dif, 2dif4dif)
  TProfile *pCor2_dif_prot_cent_forv1[2];//p+p_bar
  TProfile *pCor4_dif_prot_cent_forv1[2];

  for(int n=0; n!=3; n++){
    pCor2_ref[n] = new TProfile(Form("pCor2_ref_%i", n), "", 9, 0, 9);
    pCor4_ref[n] = new TProfile(Form("pCor4_ref_%i", n), "", 9, 0, 9);
    pCov24_ref[n] = new TProfile(Form("pCov24_ref_%i", n), "", 9, 0, 9);
    //diff corr
    for(int ich=0; ich!=2; ich++){
      pCor2_dif_cent[n][ich] = new TProfile(Form("pCor2_dif_cent_%i_%i", n, ich), "", 9, 0, 9);
      pCor4_dif_cent[n][ich] = new TProfile(Form("pCor4_dif_cent_%i_%i", n, ich), "", 9, 0, 9);

      //covariance
      for(int icov = 0; icov!=5; icov++){
        pCov_dif_cent[n][ich][icov] = new TProfile(Form("pCov_dif_cent_%i_%i_%i", n, ich, icov), "", 9, 0, 9);
      }
    }
    
    
  }

  for(int ich=0; ich!=2; ich++){
    pCor2_dif_prot_cent_forv1[ich] = new TProfile(Form("pCor2_dif_prot_cent_forv1_%i", ich), "", 9, 0, 9);
    pCor4_dif_prot_cent_forv1[ich] = new TProfile(Form("pCor4_dif_prot_cent_forv1_%i", ich), "", 9, 0, 9);
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
    double lpn_prot_cent_calc[5] = {};
    double lpn_protBar_cent_calc[5] = {};

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

      //........................Directed flow measurements........................
      if(fabs(particle->eta()) > 5.5) continue;

      Qn_calc(lQn_calc, 0, phi);

      if(pdg == 211 || pdg == 321){//Positive prticles v1 studing
        Qn_calc(lpn_pos_cent_calc, 0, phi);
      }

      if(pdg == 2212){//Proton v1 studing
        Qn_calc(lpn_prot_cent_calc, 0, phi);
      }

      if(pdg == -211 || pdg == -321){//Negative prticles v1 studing
        Qn_calc(lpn_neg_cent_calc, 0, phi);
      }

      if(pdg == -2212){//Anti Proton v1 studing
        Qn_calc(lpn_protBar_cent_calc, 0, phi);
      }
      //........................end of Directed flow measurements........................


      //..................Elliptic and triangular flow measurements......................
      if(fabs(particle->eta()) > 1.) continue;

      for(int n = 1; n!=3; n++){
        Qn_calc(lQn_calc, n, phi);

        if(pdg == 211 || pdg == 321 || pdg == 2212){//Positive prticles flow studing
          Qn_calc(lpn_pos_cent_calc, n, phi);
        }

        if(pdg == -211 || pdg == -321 || pdg == -2212){//Negative prticles flow studing
          Qn_calc(lpn_neg_cent_calc, n, phi);
        }
      }

      //..................end of Elliptic and triangular flow measurements......................

    }//for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    
    //.......................Flow calculations by cumulants.......................
    if(cent<3) continue;
    
    for(int n=0; n!=3; n++){
      //Reference flow calculation
      TComplex Qn(lQn_calc[5*n], lQn_calc[5*n+1]);
      TComplex Q2n(lQn_calc[5*n+2], lQn_calc[5*n+3]);
      int M = int(lQn_calc[5*n+4]);

      double ev_2cor_ref = Cacl_2cov_ref(Qn, M);
      double ev_4cor_ref = Cacl_4cov_ref(Qn, Q2n, M);

      double weight_ref[2] = {M*(M-1), M*(M-1)*(M-2)*(M-3)};//w2_ref, w4_ref

      pCor2_ref[n]->Fill(cent, ev_2cor_ref, weight_ref[0]);
      pCor4_ref[n]->Fill(cent, ev_4cor_ref, weight_ref[1]);
      pCov24_ref[n]->Fill(cent, ev_2cor_ref*ev_4cor_ref, weight_ref[0]*weight_ref[1]);

      //Diferential flow calculations
      if(lpn_pos_cent_calc[5*n+4]!=0)
      {//Centrality dependence of positive particles flow
        TComplex pn_pos_cent(lpn_pos_cent_calc[5*n], lpn_pos_cent_calc[5*n+1]);
        TComplex p2n_pos_cent(lpn_pos_cent_calc[5*n+2], lpn_pos_cent_calc[5*n+3]);
        int mp_cent = int(lpn_pos_cent_calc[5*n+4]);

        double ev_2cor_dif_cent = Cacl_2cov_dif(Qn, pn_pos_cent, M, mp_cent, mp_cent);
        double ev_4cor_dif_cent = Cacl_4cov_dif(Qn, pn_pos_cent, pn_pos_cent, Q2n, p2n_pos_cent, M, mp_cent, mp_cent);

        double weight_dif[2] = {mp_cent*M - mp_cent, (mp_cent*M - 3*mp_cent)*(M-1)*(M-2)};//wc2_dif, wc4_dif
        
        pCor2_dif_cent[n][0]->Fill(cent, ev_2cor_dif_cent, weight_dif[0]);
        pCor4_dif_cent[n][0]->Fill(cent, ev_4cor_dif_cent, weight_dif[1]);
        //covariance
        pCov_dif_cent[n][0][0]->Fill(cent, ev_2cor_ref*ev_2cor_dif_cent, weight_ref[0]*weight_dif[0]);
        pCov_dif_cent[n][0][1]->Fill(cent, ev_2cor_ref*ev_4cor_dif_cent, weight_ref[0]*weight_dif[1]);
        pCov_dif_cent[n][0][2]->Fill(cent, ev_4cor_ref*ev_2cor_dif_cent, weight_ref[1]*weight_dif[0]);
        pCov_dif_cent[n][0][3]->Fill(cent, ev_4cor_ref*ev_4cor_dif_cent, weight_ref[1]*weight_dif[1]);
        pCov_dif_cent[n][0][4]->Fill(cent, ev_2cor_dif_cent*ev_4cor_dif_cent, weight_dif[0]*weight_dif[1]);
      }//


      if(lpn_neg_cent_calc[5*n+4]!=0)
      {//Centrality dependence of negative particles flow
        TComplex pn_neg_cent(lpn_neg_cent_calc[5*n], lpn_neg_cent_calc[5*n+1]);
        TComplex p2n_neg_cent(lpn_neg_cent_calc[5*n+2], lpn_neg_cent_calc[5*n+3]);
        int mp_cent = int(lpn_neg_cent_calc[5*n+4]);

        double ev_2cor_dif_cent = Cacl_2cov_dif(Qn, pn_neg_cent, M, mp_cent, mp_cent);
        double ev_4cor_dif_cent = Cacl_4cov_dif(Qn, pn_neg_cent, pn_neg_cent, Q2n, p2n_neg_cent, M, mp_cent, mp_cent);

        double weight_dif[2] = {mp_cent*M - mp_cent, (mp_cent*M - 3*mp_cent)*(M-1)*(M-2)};//wc2_dif, wc4_dif

        pCor2_dif_cent[n][1]->Fill(cent, ev_2cor_dif_cent, weight_dif[0]);
        pCor4_dif_cent[n][1]->Fill(cent, ev_4cor_dif_cent, weight_dif[1]);
        //covariance
        pCov_dif_cent[n][1][0]->Fill(cent, ev_2cor_ref*ev_2cor_dif_cent, weight_ref[0]*weight_dif[0]);
        pCov_dif_cent[n][1][1]->Fill(cent, ev_2cor_ref*ev_4cor_dif_cent, weight_ref[0]*weight_dif[1]);
        pCov_dif_cent[n][1][2]->Fill(cent, ev_4cor_ref*ev_2cor_dif_cent, weight_ref[1]*weight_dif[0]);
        pCov_dif_cent[n][1][3]->Fill(cent, ev_4cor_ref*ev_4cor_dif_cent, weight_ref[1]*weight_dif[1]);
        pCov_dif_cent[n][1][4]->Fill(cent, ev_2cor_dif_cent*ev_4cor_dif_cent, weight_dif[0]*weight_dif[1]);

      }//    
      
    }//for(int n=0; n!=3; n++)
  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  //Set TProfile titles
  for(int n=0; n!=3; n++){
    pCor2_ref[n]->SetTitle(Form("Cor2 ref cum for v_%i", n+1));
    pCor4_ref[n]->SetTitle(Form("Cor4 ref cum for v_%i", n+1));
    if(n!=0){//elliptic and triangular
      pCor2_dif_cent[n][0]->SetTitle(Form("Cor2 dif cum for v_%i for (K^{+}, #pi^{+}, p^{+})", n+1));
      pCor4_dif_cent[n][0]->SetTitle(Form("Cor4 dif cum for v_%i for (K^{+}, #pi^{+}, p^{+})", n+1));
      pCor2_dif_cent[n][1]->SetTitle(Form("Cor2 dif cum for v_%i for (K^{-}, #pi^{-}, p^{-})", n+1));
      pCor4_dif_cent[n][1]->SetTitle(Form("Cor4 dif cum for v_%i for (K^{-}, #pi^{-}, p^{-})", n+1));
    }else{//directed
      pCor2_dif_cent[n][0]->SetTitle(Form("Cor2 dif cum for v_%i for (K^{+}, #pi^{+})", n+1));
      pCor4_dif_cent[n][0]->SetTitle(Form("Cor4 dif cum for v_%i for (K^{+}, #pi^{+})", n+1));
      pCor2_dif_cent[n][1]->SetTitle(Form("Cor2 dif cum for v_%i for (K^{-}, #pi^{-})", n+1));
      pCor4_dif_cent[n][1]->SetTitle(Form("Cor4 dif cum for v_%i for (K^{-}, #pi^{-})", n+1));
    }

    //Covariance
    for(int ich=0; ich!=2; ich++){
      pCov_dif_cent[n][ich][0]->SetTitle(Form("Covariance between 22' for v_%i if charge=%i", n+1, 1-2*ich));
      pCov_dif_cent[n][ich][1]->SetTitle(Form("Covariance between 24' for v_%i if charge=%i", n+1, 1-2*ich));
      pCov_dif_cent[n][ich][2]->SetTitle(Form("Covariance between 42' for v_%i if charge=%i", n+1, 1-2*ich));
      pCov_dif_cent[n][ich][3]->SetTitle(Form("Covariance between 44' for v_%i if charge=%i", n+1, 1-2*ich));
      pCov_dif_cent[n][ich][4]->SetTitle(Form("Covariance between 2'4' for v_%i if charge=%i", n+1, 1-2*ich));
    }
  }

  //Hists for proton v1  
  pCor2_dif_prot_cent_forv1[0]->SetTitle("Cor2 dif cum for p^{+} v_{1}");
  pCor2_dif_prot_cent_forv1[1]->SetTitle("Cor2 dif cum for p^{-} v_{1}");
  pCor4_dif_prot_cent_forv1[0]->SetTitle("Cor4 dif cum for p^{+} v_{1}");
  pCor4_dif_prot_cent_forv1[1]->SetTitle("Cor4 dif cum for p^{-} v_{1}");
  


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


void Qn_calc(double lQn_calc[], int nOrd, double phi){
  lQn_calc[5*nOrd] += TMath::Cos((nOrd+1)*phi);
  lQn_calc[5*nOrd+1] += TMath::Sin((nOrd+1)*phi);
  lQn_calc[5*nOrd+2] += TMath::Cos(2*(nOrd+1)*phi);
  lQn_calc[5*nOrd+3] += TMath::Sin(2*(nOrd+1)*phi);
  lQn_calc[5*nOrd +4] += 1;
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
  double up_part = (pn*Qn_star*Qn_squared - q2n*Qn_star*Qn_star - pn*Qn*Q2n_star - 2*M*pn*Qn_star - 2*mq*Qn_squared +\
                    7*qn*Qn_star - Qn*qn_star + q2n*Q2n_star + 2*pn*Qn_star + 2*mq*M - 6*mq).Re();
  return up_part / (mp*M - 3*mq) / (M-1) / (M-2);
}


