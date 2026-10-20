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
#include "TDirectory.h"
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
double GetPsi(int iOrd, double Qx, double Qy);

//pt and eta intervals
const int n_pt_points = 12;
const double pt_intervals[] = {0.15, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2., 2.4, 2.8, 3.};

const int n_eta_points = 28;


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
  oFile->cd();
  TDirectory *ptDep = oFile->mkdir("pt_dependences");
  TDirectory *etaDep = oFile->mkdir("eta_dependences");
  TDirectory *v2_RP_meth = oFile->mkdir("v2_RP_method_dir");
  TDirectory *v1_RP_meth = oFile->mkdir("v1_RP_method_dir");

  //Hist initialization
  TH1F *hpT = new TH1F("hpT", "p_{T} of particles", 570, 0.15, 3);
  TH1F *hEta = new TH1F("hEta", "#eta of particles", 600, -6, 6);
  TH1F *hcent = new TH1F("hcent", "Event Centrality", 9, 0, 9);

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

  
  //pt hists initialization
  ptDep->cd();

  //diff corr
  TProfile *pCor2_dif_pt[3][2][n_pt_points]; //(v1+v2+v3)x(pos, neg)
  TProfile *pCor4_dif_pt[3][2][n_pt_points];
  TProfile *pCov_dif_pt[3][2][n_pt_points][5];//(v1+v2+v3)x(pos, neg)x(22dif, 24dif, 42dif, 44dif, 2dif4dif)

  for(int n=0; n!=3; n++){
    //diff corr
    for(int ich=0; ich!=2; ich++){
      for(int ipt=0; ipt!=n_pt_points; ipt++){
        pCor2_dif_pt[n][ich][ipt] = new TProfile(Form("pCor2_dif_pt_%i_%i_%i", n, ich, ipt), "", 9, 0, 9);
        pCor4_dif_pt[n][ich][ipt] = new TProfile(Form("pCor4_dif_pt_%i_%i_%i", n, ich, ipt), "", 9, 0, 9);

        //covariance
        for(int icov = 0; icov!=5; icov++){
          pCov_dif_pt[n][ich][ipt][icov] = new TProfile(Form("pCov_dif_pt_%i_%i_%i_%i", n, ich, ipt, icov), "", 9, 0, 9);
        }
      }
    }    
  }

  oFile->cd();
  //end of pt hists initialization


  //eta hists initialization
  etaDep->cd();

  TProfile *pCor2_dif_prot_eta[2][n_eta_points]; //proton only (pos, neg)
  TProfile *pCor4_dif_prot_eta[2][n_eta_points];
  TProfile *pCor2_dif_eta[2][n_eta_points]; //pion+kaon (pos, neg)
  TProfile *pCor4_dif_eta[2][n_eta_points];

  TProfile *pCov_dif_eta[2][n_eta_points][5];
  TProfile *pCov_dif_prot_eta[2][n_eta_points][5];
  
  for(int ich=0; ich!=2; ich++){
    for(int ieta=0; ieta!=n_eta_points; ieta++){
      pCor2_dif_prot_eta[ich][ieta] = new TProfile(Form("pCor2_dif_prot_eta_%i_%i", ich, ieta), "", 9, 0, 9);
      pCor4_dif_prot_eta[ich][ieta] = new TProfile(Form("pCor4_dif_prot_eta_%i_%i", ich, ieta), "", 9, 0, 9);

      pCor2_dif_eta[ich][ieta] = new TProfile(Form("pCor2_dif_eta_%i_%i", ich, ieta), "", 9, 0, 9);
      pCor4_dif_eta[ich][ieta] = new TProfile(Form("pCor4_dif_eta_%i_%i", ich, ieta), "", 9, 0, 9);

      //covariance
      for(int icov = 0; icov!=5; icov++){
        pCov_dif_eta[ich][ieta][icov] = new TProfile(Form("pCov_dif_eta_%i_%i_%i", ich, ieta, icov), "", 9, 0, 9);
        pCov_dif_prot_eta[ich][ieta][icov] = new TProfile(Form("pCov_dif_prot_eta_%i_%i_%i", ich, ieta, icov), "", 9, 0, 9);
      }
    }
  }


  oFile->cd();
  //end of eta hists initialization


  //v2 RP dependences initialization
  v2_RP_meth->cd();

  //EP dist
  TH1F *hPsi1[9][3];
  TH1F *hQvec1[9][6];
  TProfile *CosOfDiff_1 = new TProfile("CosOfDiff_1", "CosOfDiff for Psi_1", 9, 0, 9);

  TH1F *hPsi2[9][3];
  TH1F *hQvec2[9][6];
  TProfile *CosOfDiff_2 = new TProfile("CosOfDiff_2", "CosOfDiff for Psi_2", 9, 0, 9);  

  //v2 dist
  TProfile *prCos_2phi_vsCent[2];
  TProfile *prCos_2phi_vsCent_forPt[2][n_pt_points];
  TProfile *prCos_2dphi2_vsCent[2];
  TProfile *prCos_2dphi2_vsCent_forPt[2][n_pt_points];
  TProfile *prCos_2dphi1_vsCent[2];
  TProfile *prCos_2dphi1_vsCent_forPt[2][n_pt_points];

  //Psi, Qvec hists
  for(int iCent=0; iCent!=9; iCent++){
    for(int iSub=0; iSub!=3; iSub++){
      hQvec1[iCent][2*iSub] = new TH1F(Form("hQvec1_%i_%i", iCent, 2*iSub), "", 500, -1, 1);
      hQvec1[iCent][2*iSub+1] = new TH1F(Form("hQvec1_%i_%i", iCent, 2*iSub+1), "", 500, -1, 1);
      hPsi1[iCent][iSub] = new TH1F(Form("hPsi1_%i_%i", iCent, iSub), "", 70, 0, 3.5);

      hQvec2[iCent][2*iSub] = new TH1F(Form("hQvec2_%i_%i", iCent, 2*iSub), "", 500, -1, 1);
      hQvec2[iCent][2*iSub+1] = new TH1F(Form("hQvec2_%i_%i", iCent, 2*iSub+1), "", 500, -1, 1);
      hPsi2[iCent][iSub] = new TH1F(Form("hPsi2_%i_%i", iCent, iSub), "", 70, 0, 3.5);
    }
  }
  
  for(int ich=0; ich!=2; ich++){
    
    //v2
    prCos_2dphi2_vsCent[ich] = new TProfile(Form("prCos_2dphi2_vsCent_%i", ich), "", 9, 0, 9);
    prCos_2dphi1_vsCent[ich] = new TProfile(Form("prCos_2dphi1_vsCent_%i", ich), "", 9, 0, 9);
    prCos_2phi_vsCent[ich] = new TProfile(Form("prCos_2phi_vsCent_%i", ich), "", 9, 0, 9);
    for(int ipt=0; ipt!=n_pt_points; ipt++){
      prCos_2dphi2_vsCent_forPt[ich][ipt] = new TProfile(Form("prCos_2dphi2_vsCent_forPt_%i_%i", ich, ipt), "", 9, 0, 9);
      prCos_2dphi1_vsCent_forPt[ich][ipt] = new TProfile(Form("prCos_2dphi1_vsCent_forPt_%i_%i", ich, ipt), "", 9, 0, 9);
      prCos_2phi_vsCent_forPt[ich][ipt] = new TProfile(Form("prCos_2phi_vsCent_forPt_%i_%i", ich, ipt), "", 9, 0, 9);
    }
  }

  oFile->cd();
  //end of v2 RP dependences initialization


  //v1 RP dependences initialization
  v1_RP_meth->cd();

  TProfile *prCos_phi_vsEta[2][9];
  TProfile *prCos_phi_prot_vsEta[2][9];

  TProfile *prCos_dphi1_vsEta[2][9];
  TProfile *prCos_dphi1_prot_vsEta[2][9];

  for(int ich=0; ich!=2; ich++){
    
    //v1
    for(int iCent=0; iCent!=9; iCent++){
      prCos_phi_vsEta[ich] = new TProfile(Form("prCos_phi_vsEta_%i_%i", ich, iCent), "", n_eta_points, -5.5, 5.5);
      prCos_phi_prot_vsEta[ich] = new TProfile(Form("prCos_phi_prot_vsEta_%i_%i", ich, iCent), "", n_eta_points, -5.5, 5.5);

      prCos_dphi1_vsEta[ich] = new TProfile(Form("prCos_dphi1_vsEta_%i", ich, iCent), "", n_eta_points, -5.5, 5.5);
      prCos_dphi1_prot_vsEta[ich] = new TProfile(Form("prCos_dphi1_prot_vsEta_%i", ich, iCent), "", n_eta_points, -5.5, 5.5);
    }

  }

  oFile->cd();
  //end of v1 RP dependences initialization


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
    hcent->Fill(cent);
    

    //Variables for particle number counting
    double lQn_calc[15] = {}; // (cos(nphi), sin(nphi), cos(2nphi), sin(2nphi), numPart) * (n=1,2,3)
    double lpn_pos_cent_calc[15] = {};
    double lpn_neg_cent_calc[15] = {};
    double lpn_prot_cent_calc[5] = {}; //for proton v1
    double lpn_protBar_cent_calc[5] = {}; //n_pt_points

    //for pt
    double lpn_pos_pt_calc[n_pt_points][15] = {};// (cos(nphi), sin(nphi), cos(2nphi), sin(2nphi), numPart) * (n=1,2,3)
    double lpn_neg_pt_calc[n_pt_points][15] = {};// (cos(nphi), sin(nphi), cos(2nphi), sin(2nphi), numPart) * (n=1,2,3)

    //for eta
    double lpn_pos_eta_calc[n_eta_points][5] = {};// (cos(nphi), sin(nphi), cos(2nphi), sin(2nphi), numPart) for (K+, pi+)
    double lpn_neg_eta_calc[n_eta_points][5] = {};// (cos(nphi), sin(nphi), cos(2nphi), sin(2nphi), numPart) for (K-, pi-)
    double lpn_prot_eta_calc[n_eta_points][5] = {};// (cos(nphi), sin(nphi), cos(2nphi), sin(2nphi), numPart) for (p+)
    double lpn_protBar_eta_calc[n_eta_points][5] = {};// (cos(nphi), sin(nphi), cos(2nphi), sin(2nphi), numPart) for (p-)

    //for EP method
    double lQ1_EP[9] = {}; // (east+west+comb)*(x, y, w)
    double lQ2_EP[9] = {}; // (east+west+comb)*(x, y, w)

    double Psi1_EP[3] = {};    
    double Psi2_EP[3] = {};

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

      hpT->Fill(pt);
      hEta->Fill(eta);

      if(fabs(eta) < 0.05) continue;

      //Get pt and eta bins
      int pt_bin = -1;
      for(int ipt = 0; ipt!=n_pt_points; ipt++){
        if(pt>=pt_intervals[ipt] && pt<pt_intervals[ipt+1]) { pt_bin = ipt; break;}
      }
      if(pt>=pt_intervals[n_pt_points]) pt_bin = n_pt_points;

      //.....................get Psi1_BBC.........................................
      if(fabs(eta) > 2.1 && fabs(eta) < 5.0) {
        
        int iSide = -1;
        if(eta<0) iSide = 0; //East
        else iSide = 1;      //West

        lQ1_EP[3*iSide] += 1.*TMath::Cos(phi);
        lQ1_EP[3*iSide+1] += 1.*TMath::Sin(phi);
        lQ1_EP[3*iSide+2] += 1.;
      }

      //.....................end of get Psi1_BBC.........................................


      //........................Directed flow measurements by cum........................
      if(fabs(particle->eta()) >= 5.5) continue;

      double delta_eta = 11.0 / n_eta_points;
      int eta_bin = (int)((eta + 5.5) / delta_eta);

      Qn_calc(lQn_calc, 0, phi);

      if(pdg == 211 || pdg == 321){//Positive prticles v1 studing
        Qn_calc(lpn_pos_cent_calc, 0, phi);
        Qn_calc(lpn_pos_eta_calc[eta_bin], 0, phi);
      }

      if(pdg == 2212){//Proton v1 studing
        Qn_calc(lpn_prot_cent_calc, 0, phi);
        Qn_calc(lpn_prot_eta_calc[eta_bin], 0, phi);
      }

      if(pdg == -211 || pdg == -321){//Negative prticles v1 studing
        Qn_calc(lpn_neg_cent_calc, 0, phi);
        Qn_calc(lpn_neg_eta_calc[eta_bin], 0, phi);
      }

      if(pdg == -2212){//Anti Proton v1 studing
        Qn_calc(lpn_protBar_cent_calc, 0, phi);
        Qn_calc(lpn_protBar_eta_calc[eta_bin], 0, phi);
      }
      
      //........................end of Directed flow measurements by cum........................

      if(particle->pt() < 0.15) continue;
      if(fabs(particle->eta()) > 1.) continue;
      
      //..................Elliptic and triangular flow measurements by cum......................      

      for(int n = 1; n!=3; n++){
        Qn_calc(lQn_calc, n, phi);

        if(pdg == 211 || pdg == 321 || pdg == 2212){//Positive prticles flow studing
          Qn_calc(lpn_pos_cent_calc, n, phi); //lpn_pos_cent_calc
          Qn_calc(lpn_pos_pt_calc[pt_bin], n, phi);
        }

        if(pdg == -211 || pdg == -321 || pdg == -2212){//Negative prticles flow studing
          Qn_calc(lpn_neg_cent_calc, n, phi);
          Qn_calc(lpn_neg_pt_calc[pt_bin], n, phi);
        }
      }

      //..................end of Elliptic and triangular flow measurements by cum......................


      //.....................get Psi2.........................................

      //if(pt>2.) continue;      

      int iSide = -1;
      if(eta<0) iSide = 0; //East
      else iSide = 1;      //West

      lQ2_EP[3*iSide] += pt*TMath::Cos(2*phi);
      lQ2_EP[3*iSide+1] += pt*TMath::Sin(2*phi);
      lQ2_EP[3*iSide+2] += pt;

    }//for(Int_t iTrk=0; iTrk<nTracks; iTrk++) First track loop

    //=============================================================================================
    //=============================================================================================
    //...............................Flow measurements by EP method................................

    //Get Q_vec for Psi2
    for(int iSub=0; iSub!=2; iSub++){
      if(lQ2_EP[3*iSub+2]!=0){
        lQ2_EP[3*iSub] /= lQ2_EP[3*iSub+2];
        lQ2_EP[3*iSub+1] /= lQ2_EP[3*iSub+2];
      }else{
        lQ2_EP[3*iSub] = -999.;
        lQ2_EP[3*iSub+1] = -999.;
      }
    }

    bool check_q2 = true;
    for(int iSub=0; iSub!=2; iSub++){
      if(lQ2_EP[3*iSub]== -999. || lQ2_EP[3*iSub+1]== -999.) check_q2 = false;
    }

    if(check_q2){
      lQ2_EP[6] = lQ2_EP[0] + lQ2_EP[3];
      lQ2_EP[7] = lQ2_EP[1] + lQ2_EP[4];

      //Psi angle reconstruction
      for(int iSub=0; iSub!=2; iSub++){
        Psi2_EP[iSub] = GetPsi(2, lQ2_EP[3*iSub], lQ2_EP[3*iSub+1]);
      }

      //EP hists filling
      for(int iSub=0; iSub!=3; iSub++){
        hQvec2[cent][2*iSub]->Fill(lQ2_EP[3*iSub]);
        hQvec2[cent][2*iSub+1]->Fill(lQ2_EP[3*iSub+1]);
        hPsi2[cent][iSub]->Fill(Psi2_EP[iSub]);     
      }

      CosOfDiff_2->Fill(cent, TMath::Cos(2*(Psi2_EP[1] - Psi2_EP[0])));
    }
    

    //Get Q_vec for Psi1
    for(int iSub=0; iSub!=2; iSub++){
      if(lQ1_EP[3*iSub+2]!=0){
        lQ1_EP[3*iSub] /= lQ1_EP[3*iSub+2];
        lQ1_EP[3*iSub+1] /= lQ1_EP[3*iSub+2];
      }else{
        lQ1_EP[3*iSub] = -999.;
        lQ1_EP[3*iSub+1] = -999.;
      }
    }
    
    bool check_q1 = true;
    for(int iSub=0; iSub!=2; iSub++){
      if(lQ1_EP[3*iSub]== -999. || lQ1_EP[3*iSub+1]== -999.) check_q1 = false;
    }

    if(check_q1){
      lQ1_EP[6] = lQ1_EP[0] + lQ1_EP[3];
      lQ1_EP[7] = lQ1_EP[1] + lQ1_EP[4];

      //Psi angle reconstruction
      for(int iSub=0; iSub!=3; iSub++){
        Psi1_EP[iSub] = GetPsi(1, lQ1_EP[3*iSub], lQ1_EP[3*iSub+1]);
      }

      //EP hists filling
      for(int iSub=0; iSub!=3; iSub++){
        hQvec1[cent][2*iSub]->Fill(lQ1_EP[3*iSub]);
        hQvec1[cent][2*iSub+1]->Fill(lQ1_EP[3*iSub+1]);
        hPsi1[cent][iSub]->Fill(Psi1_EP[iSub]);     
      }

      CosOfDiff_1->Fill(cent, TMath::Cos(Psi1_EP[1] - Psi1_EP[0]));
    }
    
      

    //flow calculating by EP
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++){
      // Retrieve i-th femto track
      McParticle *particle = dst->particle(iTrk);

      if (!particle) continue;

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
      if(fabs(particle->eta()) >= 5.5) continue;
      
      //Clear v1
      //v1 flow by EP method
      if(pdg == 211 || pdg == 321){
        prCos_phi_vsEta[0][cent]->Fill(eta, TMath::Cos(phi));
      }

      if(pdg == 2212){
        prCos_phi_prot_vsEta[0][cent]->Fill(eta, TMath::Cos(phi));
      }

      if(pdg == -211 || pdg == -321){
        prCos_phi_vsEta[1][cent]->Fill(eta, TMath::Cos(phi));
      }

      if(pdg == -2212){
        prCos_phi_prot_vsEta[1][cent]->Fill(eta, TMath::Cos(phi));
      }


      //v1 by Psi_BBC
      if(check_q1){
        if(pdg == 211 || pdg == 321){
          prCos_dphi1_vsEta[0][cent]->Fill(eta, TMath::Cos(phi - Psi1_EP[2]));
        }

        if(pdg == 2212){
          prCos_dphi1_prot_vsEta[0][cent]->Fill(eta, TMath::Cos(phi - Psi1_EP[2]));
        }

        if(pdg == -211 || pdg == -321){
          prCos_dphi1_vsEta[1][cent]->Fill(eta, TMath::Cos(phi - Psi1_EP[2]));
        }

        if(pdg == -2212){
          prCos_dphi1_prot_vsEta[1][cent]->Fill(eta, TMath::Cos(phi - Psi1_EP[2]));
        }
      }



      if(particle->pt() < 0.15) continue;
      if(fabs(particle->eta()) > 1.) continue;

      //Get pt and eta bins
      int pt_bin = -1;
      for(int ipt = 0; ipt!=n_pt_points; ipt++){
        if(pt>=pt_intervals[ipt] && pt<pt_intervals[ipt+1]) { pt_bin = ipt; break;}
      }
      if(pt>=pt_intervals[n_pt_points]) pt_bin = n_pt_points;      

      
      //Clear v2
      if(pdg == 211 || pdg == 321 || pdg == 2212)
      {//p+K+pi
        double v2 = TMath::Cos(2*(phi));
        prCos_2phi_vsCent[0]->Fill(cent, v2);
        prCos_2phi_vsCent_forPt[0][pt_bin]->Fill(cent, v2);
      }

      if(pdg == -211 || pdg == -321 || pdg == -2212)
      {//p_bar+K_bar+pi_bar
        double v2 = TMath::Cos(2*(phi));
        prCos_2phi_vsCent[1]->Fill(cent, v2);
        prCos_2phi_vsCent_forPt[1][pt_bin]->Fill(cent, v2);
      }

      
      //v2 by Psi1
      if(check_q1){
        
        if(pdg == 211 || pdg == 321 || pdg == 2212)
        {//p+K+pi
          double v2 = TMath::Cos(2*(phi - Psi1_EP[2]));
          prCos_2dphi1_vsCent[0]->Fill(cent, v2);
          prCos_2dphi1_vsCent_forPt[0][pt_bin]->Fill(cent, v2);
        }

        if(pdg == -211 || pdg == -321 || pdg == -2212)
        {//p_bar+K_bar+pi_bar
          double v2 = TMath::Cos(2*(phi - Psi1_EP[2]));
          prCos_2dphi1_vsCent[1]->Fill(cent, v2);
          prCos_2dphi1_vsCent_forPt[1][pt_bin]->Fill(cent, v2);
        }
      }

      
      //v2 by Psi2
      if(check_q2){
        //Get iSide number
        int iSide = -1;
        if(eta<0) iSide = 1; //East
        else iSide = 0;      //West

        if(pdg == 211 || pdg == 321 || pdg == 2212)
        {//p+K+pi
          double v2 = TMath::Cos(2*(phi - Psi2_EP[iSide]));
          prCos_2dphi2_vsCent[0]->Fill(cent, v2);
          prCos_2dphi2_vsCent_forPt[0][pt_bin]->Fill(cent, v2);
        }

        if(pdg == -211 || pdg == -321 || pdg == -2212)
        {//p_bar+K_bar+pi_bar
          double v2 = TMath::Cos(2*(phi - Psi2_EP[iSide]));
          prCos_2dphi2_vsCent[1]->Fill(cent, v2);
          prCos_2dphi2_vsCent_forPt[1][pt_bin]->Fill(cent, v2);
        }
      }


    }//for(Int_t iTrk=0; iTrk<nTracks; iTrk++) Second track loop
    

    //...............................end of flow measurements by EP method.........................
    //=============================================================================================
    //=============================================================================================
    

    
    //=============================================================================================
    //=============================================================================================
    //...............................Flow measurements by cumulants................................
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
      
      
      //............................Pt dependeces...................
      for(int ipt=0; ipt!=n_pt_points; ipt++){
        if(lpn_pos_pt_calc[ipt][5*n+4]!=0)
        {//pT dependence of positive particles flow
          TComplex pn_pos_cent(lpn_pos_pt_calc[ipt][5*n], lpn_pos_pt_calc[ipt][5*n+1]);
          TComplex p2n_pos_cent(lpn_pos_pt_calc[ipt][5*n+2], lpn_pos_pt_calc[ipt][5*n+3]);
          int mp_cent = int(lpn_pos_pt_calc[ipt][5*n+4]);

          double ev_2cor_dif_cent = Cacl_2cov_dif(Qn, pn_pos_cent, M, mp_cent, mp_cent);
          double ev_4cor_dif_cent = Cacl_4cov_dif(Qn, pn_pos_cent, pn_pos_cent, Q2n, p2n_pos_cent, M, mp_cent, mp_cent);

          double weight_dif[2] = {mp_cent*M - mp_cent, (mp_cent*M - 3*mp_cent)*(M-1)*(M-2)};//wc2_dif, wc4_dif
          
          pCor2_dif_pt[n][0][ipt]->Fill(cent, ev_2cor_dif_cent, weight_dif[0]);
          pCor4_dif_pt[n][0][ipt]->Fill(cent, ev_4cor_dif_cent, weight_dif[1]);
          //covariance
          pCov_dif_pt[n][0][ipt][0]->Fill(cent, ev_2cor_ref*ev_2cor_dif_cent, weight_ref[0]*weight_dif[0]);
          pCov_dif_pt[n][0][ipt][1]->Fill(cent, ev_2cor_ref*ev_4cor_dif_cent, weight_ref[0]*weight_dif[1]);
          pCov_dif_pt[n][0][ipt][2]->Fill(cent, ev_4cor_ref*ev_2cor_dif_cent, weight_ref[1]*weight_dif[0]);
          pCov_dif_pt[n][0][ipt][3]->Fill(cent, ev_4cor_ref*ev_4cor_dif_cent, weight_ref[1]*weight_dif[1]);
          pCov_dif_pt[n][0][ipt][4]->Fill(cent, ev_2cor_dif_cent*ev_4cor_dif_cent, weight_dif[0]*weight_dif[1]);
        }//


        if(lpn_neg_pt_calc[ipt][5*n+4]!=0)
        {//pT dependence of negative particles flow
          TComplex pn_pos_cent(lpn_neg_pt_calc[ipt][5*n], lpn_neg_pt_calc[ipt][5*n+1]);
          TComplex p2n_pos_cent(lpn_neg_pt_calc[ipt][5*n+2], lpn_neg_pt_calc[ipt][5*n+3]);
          int mp_cent = int(lpn_neg_pt_calc[ipt][5*n+4]);

          double ev_2cor_dif_cent = Cacl_2cov_dif(Qn, pn_pos_cent, M, mp_cent, mp_cent);
          double ev_4cor_dif_cent = Cacl_4cov_dif(Qn, pn_pos_cent, pn_pos_cent, Q2n, p2n_pos_cent, M, mp_cent, mp_cent);

          double weight_dif[2] = {mp_cent*M - mp_cent, (mp_cent*M - 3*mp_cent)*(M-1)*(M-2)};//wc2_dif, wc4_dif
          
          pCor2_dif_pt[n][1][ipt]->Fill(cent, ev_2cor_dif_cent, weight_dif[0]);
          pCor4_dif_pt[n][1][ipt]->Fill(cent, ev_4cor_dif_cent, weight_dif[1]);
          //covariance
          pCov_dif_pt[n][1][ipt][0]->Fill(cent, ev_2cor_ref*ev_2cor_dif_cent, weight_ref[0]*weight_dif[0]);
          pCov_dif_pt[n][1][ipt][1]->Fill(cent, ev_2cor_ref*ev_4cor_dif_cent, weight_ref[0]*weight_dif[1]);
          pCov_dif_pt[n][1][ipt][2]->Fill(cent, ev_4cor_ref*ev_2cor_dif_cent, weight_ref[1]*weight_dif[0]);
          pCov_dif_pt[n][1][ipt][3]->Fill(cent, ev_4cor_ref*ev_4cor_dif_cent, weight_ref[1]*weight_dif[1]);
          pCov_dif_pt[n][1][ipt][4]->Fill(cent, ev_2cor_dif_cent*ev_4cor_dif_cent, weight_dif[0]*weight_dif[1]);
        }//
      }//for(int ipt=0; ipt!=n_pt_points; ipt++)
      
      //............................end of Pt dependeces...................


      //............................Eta dependeces...................
      if(n == 0){
        for(int ieta=0; ieta!=n_eta_points; ieta++){
          if(lpn_pos_eta_calc[ieta][5*n+4]!=0)
          {//eta dependence of K+, pi+ flow
            TComplex pn_pos_cent(lpn_pos_eta_calc[ieta][5*n], lpn_pos_eta_calc[ieta][5*n+1]);
            TComplex p2n_pos_cent(lpn_pos_eta_calc[ieta][5*n+2], lpn_pos_eta_calc[ieta][5*n+3]);
            int mp_cent = int(lpn_pos_eta_calc[ieta][5*n+4]);

            double ev_2cor_dif_cent = Cacl_2cov_dif(Qn, pn_pos_cent, M, mp_cent, mp_cent);
            double ev_4cor_dif_cent = Cacl_4cov_dif(Qn, pn_pos_cent, pn_pos_cent, Q2n, p2n_pos_cent, M, mp_cent, mp_cent);

            double weight_dif[2] = {mp_cent*M - mp_cent, (mp_cent*M - 3*mp_cent)*(M-1)*(M-2)};//wc2_dif, wc4_dif
            
            pCor2_dif_eta[0][ieta]->Fill(cent, ev_2cor_dif_cent, weight_dif[0]);
            pCor4_dif_eta[0][ieta]->Fill(cent, ev_4cor_dif_cent, weight_dif[1]);
            //covariance
            pCov_dif_eta[0][ieta][0]->Fill(cent, ev_2cor_ref*ev_2cor_dif_cent, weight_ref[0]*weight_dif[0]);
            pCov_dif_eta[0][ieta][1]->Fill(cent, ev_2cor_ref*ev_4cor_dif_cent, weight_ref[0]*weight_dif[1]);
            pCov_dif_eta[0][ieta][2]->Fill(cent, ev_4cor_ref*ev_2cor_dif_cent, weight_ref[1]*weight_dif[0]);
            pCov_dif_eta[0][ieta][3]->Fill(cent, ev_4cor_ref*ev_4cor_dif_cent, weight_ref[1]*weight_dif[1]);
            pCov_dif_eta[0][ieta][4]->Fill(cent, ev_2cor_dif_cent*ev_4cor_dif_cent, weight_dif[0]*weight_dif[1]);
          }//


          if(lpn_neg_eta_calc[ieta][5*n+4]!=0)
          {//eta dependence of K-, pi- flow
            TComplex pn_pos_cent(lpn_neg_eta_calc[ieta][5*n], lpn_neg_eta_calc[ieta][5*n+1]);
            TComplex p2n_pos_cent(lpn_neg_eta_calc[ieta][5*n+2], lpn_neg_eta_calc[ieta][5*n+3]);
            int mp_cent = int(lpn_neg_eta_calc[ieta][5*n+4]);

            double ev_2cor_dif_cent = Cacl_2cov_dif(Qn, pn_pos_cent, M, mp_cent, mp_cent);
            double ev_4cor_dif_cent = Cacl_4cov_dif(Qn, pn_pos_cent, pn_pos_cent, Q2n, p2n_pos_cent, M, mp_cent, mp_cent);

            double weight_dif[2] = {mp_cent*M - mp_cent, (mp_cent*M - 3*mp_cent)*(M-1)*(M-2)};//wc2_dif, wc4_dif
            
            pCor2_dif_eta[1][ieta]->Fill(cent, ev_2cor_dif_cent, weight_dif[0]);
            pCor4_dif_eta[1][ieta]->Fill(cent, ev_4cor_dif_cent, weight_dif[1]);
            //covariance
            pCov_dif_eta[1][ieta][0]->Fill(cent, ev_2cor_ref*ev_2cor_dif_cent, weight_ref[0]*weight_dif[0]);
            pCov_dif_eta[1][ieta][1]->Fill(cent, ev_2cor_ref*ev_4cor_dif_cent, weight_ref[0]*weight_dif[1]);
            pCov_dif_eta[1][ieta][2]->Fill(cent, ev_4cor_ref*ev_2cor_dif_cent, weight_ref[1]*weight_dif[0]);
            pCov_dif_eta[1][ieta][3]->Fill(cent, ev_4cor_ref*ev_4cor_dif_cent, weight_ref[1]*weight_dif[1]);
            pCov_dif_eta[1][ieta][4]->Fill(cent, ev_2cor_dif_cent*ev_4cor_dif_cent, weight_dif[0]*weight_dif[1]);
          }//


          if(lpn_prot_eta_calc[ieta][5*n+4]!=0)
          {//eta dependence of p+ flow
            TComplex pn_pos_cent(lpn_prot_eta_calc[ieta][5*n], lpn_prot_eta_calc[ieta][5*n+1]);
            TComplex p2n_pos_cent(lpn_prot_eta_calc[ieta][5*n+2], lpn_prot_eta_calc[ieta][5*n+3]);
            int mp_cent = int(lpn_prot_eta_calc[ieta][5*n+4]);

            double ev_2cor_dif_cent = Cacl_2cov_dif(Qn, pn_pos_cent, M, mp_cent, mp_cent);
            double ev_4cor_dif_cent = Cacl_4cov_dif(Qn, pn_pos_cent, pn_pos_cent, Q2n, p2n_pos_cent, M, mp_cent, mp_cent);

            double weight_dif[2] = {mp_cent*M - mp_cent, (mp_cent*M - 3*mp_cent)*(M-1)*(M-2)};//wc2_dif, wc4_dif
            
            pCor2_dif_prot_eta[0][ieta]->Fill(cent, ev_2cor_dif_cent, weight_dif[0]);
            pCor4_dif_prot_eta[0][ieta]->Fill(cent, ev_4cor_dif_cent, weight_dif[1]);
            //covariance
            pCov_dif_prot_eta[0][ieta][0]->Fill(cent, ev_2cor_ref*ev_2cor_dif_cent, weight_ref[0]*weight_dif[0]);
            pCov_dif_prot_eta[0][ieta][1]->Fill(cent, ev_2cor_ref*ev_4cor_dif_cent, weight_ref[0]*weight_dif[1]);
            pCov_dif_prot_eta[0][ieta][2]->Fill(cent, ev_4cor_ref*ev_2cor_dif_cent, weight_ref[1]*weight_dif[0]);
            pCov_dif_prot_eta[0][ieta][3]->Fill(cent, ev_4cor_ref*ev_4cor_dif_cent, weight_ref[1]*weight_dif[1]);
            pCov_dif_prot_eta[0][ieta][4]->Fill(cent, ev_2cor_dif_cent*ev_4cor_dif_cent, weight_dif[0]*weight_dif[1]);
          }//


          if(lpn_protBar_eta_calc[ieta][5*n+4]!=0)
          {//eta dependence of p- flow
            TComplex pn_pos_cent(lpn_protBar_eta_calc[ieta][5*n], lpn_protBar_eta_calc[ieta][5*n+1]);
            TComplex p2n_pos_cent(lpn_protBar_eta_calc[ieta][5*n+2], lpn_protBar_eta_calc[ieta][5*n+3]);
            int mp_cent = int(lpn_protBar_eta_calc[ieta][5*n+4]);

            double ev_2cor_dif_cent = Cacl_2cov_dif(Qn, pn_pos_cent, M, mp_cent, mp_cent);
            double ev_4cor_dif_cent = Cacl_4cov_dif(Qn, pn_pos_cent, pn_pos_cent, Q2n, p2n_pos_cent, M, mp_cent, mp_cent);

            double weight_dif[2] = {mp_cent*M - mp_cent, (mp_cent*M - 3*mp_cent)*(M-1)*(M-2)};//wc2_dif, wc4_dif
            
            pCor2_dif_prot_eta[1][ieta]->Fill(cent, ev_2cor_dif_cent, weight_dif[0]);
            pCor4_dif_prot_eta[1][ieta]->Fill(cent, ev_4cor_dif_cent, weight_dif[1]);
            //covariance
            pCov_dif_prot_eta[1][ieta][0]->Fill(cent, ev_2cor_ref*ev_2cor_dif_cent, weight_ref[0]*weight_dif[0]);
            pCov_dif_prot_eta[1][ieta][1]->Fill(cent, ev_2cor_ref*ev_4cor_dif_cent, weight_ref[0]*weight_dif[1]);
            pCov_dif_prot_eta[1][ieta][2]->Fill(cent, ev_4cor_ref*ev_2cor_dif_cent, weight_ref[1]*weight_dif[0]);
            pCov_dif_prot_eta[1][ieta][3]->Fill(cent, ev_4cor_ref*ev_4cor_dif_cent, weight_ref[1]*weight_dif[1]);
            pCov_dif_prot_eta[1][ieta][4]->Fill(cent, ev_2cor_dif_cent*ev_4cor_dif_cent, weight_dif[0]*weight_dif[1]);
          }//
        }//for(int ieta=0; ieta!=n_eta_points; ieta++)
      }//if(n == 0)
      
    }//for(int n=0; n!=3; n++)

    //...............................end of flow measurements by cumulants................................
    //=============================================================================================
    //=============================================================================================
    
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


  //for pt dep
  for(int n=0; n!=3; n++){
    for(int ipt=0; ipt!=n_pt_points; ipt++){
      if(n!=0){//elliptic and triangular
        pCor2_dif_pt[n][0][ipt]->SetTitle(Form("Cor2 dif cum for v_%i for (K^{+}, #pi^{+}, p^{+})", n+1));
        pCor4_dif_pt[n][0][ipt]->SetTitle(Form("Cor4 dif cum for v_%i for (K^{+}, #pi^{+}, p^{+})", n+1));
        pCor2_dif_pt[n][1][ipt]->SetTitle(Form("Cor2 dif cum for v_%i for (K^{-}, #pi^{-}, p^{-})", n+1));
        pCor4_dif_pt[n][1][ipt]->SetTitle(Form("Cor4 dif cum for v_%i for (K^{-}, #pi^{-}, p^{-})", n+1));
      }else{//directed
        pCor2_dif_pt[n][0][ipt]->SetTitle(Form("Cor2 dif cum for v_%i for (K^{+}, #pi^{+})", n+1));
        pCor4_dif_pt[n][0][ipt]->SetTitle(Form("Cor4 dif cum for v_%i for (K^{+}, #pi^{+})", n+1));
        pCor2_dif_pt[n][1][ipt]->SetTitle(Form("Cor2 dif cum for v_%i for (K^{-}, #pi^{-})", n+1));
        pCor4_dif_pt[n][1][ipt]->SetTitle(Form("Cor4 dif cum for v_%i for (K^{-}, #pi^{-})", n+1));
      }

      //Covariance
      for(int ich=0; ich!=2; ich++){
        pCov_dif_pt[n][ich][ipt][0]->SetTitle(Form("Covariance between 22' for v_%i if charge=%i", n+1, 1-2*ich));
        pCov_dif_pt[n][ich][ipt][1]->SetTitle(Form("Covariance between 24' for v_%i if charge=%i", n+1, 1-2*ich));
        pCov_dif_pt[n][ich][ipt][2]->SetTitle(Form("Covariance between 42' for v_%i if charge=%i", n+1, 1-2*ich));
        pCov_dif_pt[n][ich][ipt][3]->SetTitle(Form("Covariance between 44' for v_%i if charge=%i", n+1, 1-2*ich));
        pCov_dif_pt[n][ich][ipt][4]->SetTitle(Form("Covariance between 2'4' for v_%i if charge=%i", n+1, 1-2*ich));
      }
    }
  }


  //for eta dep
  for(int ieta=0; ieta!=n_eta_points; ieta++){
    //directed
    pCor2_dif_eta[0][ieta]->SetTitle("Cor2 dif cum for v_1 for (K^{+}, #pi^{+})");
    pCor4_dif_eta[0][ieta]->SetTitle("Cor4 dif cum for v_1 for (K^{+}, #pi^{+})");
    pCor2_dif_eta[1][ieta]->SetTitle("Cor2 dif cum for v_1 for (K^{-}, #pi^{-})");
    pCor4_dif_eta[1][ieta]->SetTitle("Cor4 dif cum for v_1 for (K^{-}, #pi^{-})");

    pCor2_dif_prot_eta[0][ieta]->SetTitle("Cor2 dif cum for v_1 for (p^{+})");
    pCor4_dif_prot_eta[0][ieta]->SetTitle("Cor4 dif cum for v_1 for (p^{+})");
    pCor2_dif_prot_eta[1][ieta]->SetTitle("Cor2 dif cum for v_1 for (p^{-})");
    pCor4_dif_prot_eta[1][ieta]->SetTitle("Cor4 dif cum for v_1 for (p^{-})");
    

    //Covariance
    for(int ich=0; ich!=2; ich++){
      pCov_dif_eta[ich][ieta][0]->SetTitle(Form("Covariance between 22' for v_1 if charge=%i", 1-2*ich));
      pCov_dif_eta[ich][ieta][1]->SetTitle(Form("Covariance between 24' for v_1 if charge=%i", 1-2*ich));
      pCov_dif_eta[ich][ieta][2]->SetTitle(Form("Covariance between 42' for v_1 if charge=%i", 1-2*ich));
      pCov_dif_eta[ich][ieta][3]->SetTitle(Form("Covariance between 44' for v_1 if charge=%i", 1-2*ich));
      pCov_dif_eta[ich][ieta][4]->SetTitle(Form("Covariance between 2'4' for v_1 if charge=%i", 1-2*ich));

      //for proton
      pCov_dif_prot_eta[ich][ieta][0]->SetTitle(Form("Covariance between 22' for v_1 if charge=%i", 1-2*ich));
      pCov_dif_prot_eta[ich][ieta][1]->SetTitle(Form("Covariance between 24' for v_1 if charge=%i", 1-2*ich));
      pCov_dif_prot_eta[ich][ieta][2]->SetTitle(Form("Covariance between 42' for v_1 if charge=%i", 1-2*ich));
      pCov_dif_prot_eta[ich][ieta][3]->SetTitle(Form("Covariance between 44' for v_1 if charge=%i", 1-2*ich));
      pCov_dif_prot_eta[ich][ieta][4]->SetTitle(Form("Covariance between 2'4' for v_1 if charge=%i", 1-2*ich));
    }
  }
  
  //EP hists
  for(int iCent=0; iCent!=9; iCent++){
    hQvec1[iCent][0]->SetTitle(Form("Qx East for %i centrality", iCent));
    hQvec1[iCent][1]->SetTitle(Form("Qy East for %i centrality", iCent));
    hQvec1[iCent][2]->SetTitle(Form("Qx West for %i centrality", iCent));
    hQvec1[iCent][3]->SetTitle(Form("Qy West for %i centrality", iCent));
    hPsi1[iCent][0]->SetTitle(Form("Psi1 East for %i centrality", iCent));
    hPsi1[iCent][1]->SetTitle(Form("Psi1 West for %i centrality", iCent));
    
    hQvec2[iCent][0]->SetTitle(Form("Qx East for %i centrality", iCent));
    hQvec2[iCent][1]->SetTitle(Form("Qy East for %i centrality", iCent));
    hQvec2[iCent][2]->SetTitle(Form("Qx West for %i centrality", iCent));
    hQvec2[iCent][3]->SetTitle(Form("Qy West for %i centrality", iCent));
    hPsi2[iCent][0]->SetTitle(Form("Psi2 East for %i centrality", iCent));
    hPsi2[iCent][1]->SetTitle(Form("Psi2 West for %i centrality", iCent));
  }
  
  //flow by EP meth hists

  
  for(int iCent=0; iCent!=9; iCent++){
    //Clear v1
    prCos_phi_vsEta[0][iCent]->SetTitle(Form("Cos(#phi) Vs eta for (K^{+}, #pi^{+}) for iCent=(%i)", iCent));
    prCos_phi_vsEta[1][iCent]->SetTitle(Form("Cos(#phi) Vs eta for (K^{-}, #pi^{-}) for iCent=(%i)", iCent));
    prCos_phi_prot_vsEta[0][iCent]->SetTitle(Form("Cos(#phi) Vs eta for (p^{+}) for iCent=(%i)", iCent));
    prCos_phi_prot_vsEta[1][iCent]->SetTitle(Form("Cos(#phi) Vs eta for (p^{-}) for iCent=(%i)", iCent));

    //v1 by Psi_BBC
    prCos_dphi1_vsEta[0][iCent]->SetTitle(Form("Cos(#phi-#Psi_{BBC}) Vs eta for (K^{+}, #pi^{+}) for iCent=(%i)", iCent));
    prCos_dphi1_vsEta[1][iCent]->SetTitle(Form("Cos(#phi-#Psi_{BBC}) Vs eta for (K^{-}, #pi^{-}) for iCent=(%i)", iCent));
    prCos_dphi1_prot_vsEta[0][iCent]->SetTitle(Form("Cos(#phi-#Psi_{BBC}) Vs eta for (p^{+}) for iCent=(%i)", iCent));
    prCos_dphi1_prot_vsEta[1][iCent]->SetTitle(Form("Cos(#phi-#Psi_{BBC}) Vs eta for (p^{-}) for iCent=(%i)", iCent));
  }

  //v2 by Psi_2
  prCos_2dphi2_vsCent[0]->SetTitle("Cos(2#phi-2#Psi_{2}) Vs cent for (K^{+}, #pi^{+}, p^{+})");
  prCos_2dphi2_vsCent[1]->SetTitle("Cos(2#phi-2#Psi_{2}) Vs cent for (K^{+}, #pi^{+}, p^{+})");

  //v2 by Psi_BBC
  prCos_2dphi1_vsCent[0]->SetTitle("Cos(2#phi-2#Psi_{2}) Vs cent for (K^{+}, #pi^{+}, p^{+})");
  prCos_2dphi1_vsCent[1]->SetTitle("Cos(2#phi-2#Psi_{2}) Vs cent for (K^{+}, #pi^{+}, p^{+})");

  //Clear v2
  prCos_2phi_vsCent[0]->SetTitle("Cos(2#phi) Vs cent for (K^{+}, #pi^{+}, p^{+})");
  prCos_2phi_vsCent[1]->SetTitle("Cos(2#phi) Vs cent for (K^{+}, #pi^{+}, p^{+})");

  for(int ipt=0; ipt!=n_pt_points; ipt++){
    //v2 by Psi_2
    prCos_2dphi2_vsCent_forPt[0][ipt]->SetTitle(Form("Cos(2#phi-2#Psi_{2}) Vs cent for (K^{+}, #pi^{+}, p^{+}) for ipt=(%i)", ipt));
    prCos_2dphi2_vsCent_forPt[1][ipt]->SetTitle(Form("Cos(2#phi-2#Psi_{2}) Vs cent for (K^{+}, #pi^{+}, p^{+}) for ipt=(%i)", ipt));

    //v2 by Psi_BBC
    prCos_2dphi1_vsCent_forPt[0][ipt]->SetTitle(Form("Cos(2#phi-2#Psi_{BBC}) Vs cent for (K^{+}, #pi^{+}, p^{+}) for ipt=(%i)", ipt));
    prCos_2dphi1_vsCent_forPt[1][ipt]->SetTitle(Form("Cos(2#phi-2#Psi_{BBC}) Vs cent for (K^{+}, #pi^{+}, p^{+}) for ipt=(%i)", ipt));
    
    //Clear v2
    prCos_2phi_vsCent_forPt[0][ipt]->SetTitle(Form("Cos(2#phi) Vs cent for (K^{+}, #pi^{+}, p^{+}) for ipt=(%i)", ipt));
    prCos_2phi_vsCent_forPt[1][ipt]->SetTitle(Form("Cos(2#phi) Vs cent for (K^{+}, #pi^{+}, p^{+}) for ipt=(%i)", ipt));
  }


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


//Get Psi angle, iOrd starts with 1
double GetPsi(int iOrd, double Qx, double Qy){
    double Psi;
    Psi = atan2(Qy, Qx)/iOrd;
    if(Psi<0.) Psi += 2*TMath::Pi() / iOrd;
    return Psi;

}


