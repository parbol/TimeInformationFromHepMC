//=======================================================================================================================================================================================================================// 
//                                                                                                                                                                                                                       // 
//$$$$$$\ $$$$$$$$\  $$$$$$\   $$$$$$\                      $$\                                    $$\       $$\                            $$\  $$$$$$\                      $$\                     $$                 //
//\_$$  _|$$  _____|$$  __$$\ $$  __$$\                     $$ |                                   $$ |      \__|                           $$ |$$  __$$\                     $$ |                    \__|               //
//  $$ |  $$ |      $$ /  \__|$$ /  $$ |                    $$ |      $$$$$$\  $$$$$$$\   $$$$$$\  $$ |      $$\ $$\    $$\  $$$$$$\   $$$$$$$ |$$ /  $$ |$$$$$$$\   $$$$$$\  $$ |$$\   $$\  $$$$$$$\ $$\  $$$$$$$       //
//  $$ |  $$$$$\    $$ |      $$$$$$$$ |      $$$$$$\       $$ |     $$  __$$\ $$  __$$\ $$  __$$\ $$ |      $$ |\$$\  $$  |$$  __$$\ $$  __$$ |$$$$$$$$ |$$  __$$\  \____$$\ $$ |$$ |  $$ |$$  _____|$$ |$$  _____|     //
//  $$ |  $$  __|   $$ |      $$  __$$ |      \______|      $$ |     $$ /  $$ |$$ |  $$ |$$ /  $$ |$$ |      $$ | \$$\$$  / $$$$$$$$ |$$ /  $$ |$$  __$$ |$$ |  $$ | $$$$$$$ |$$ |$$ |  $$ |\$$$$$$\  $$ |\$$$$$$        // 
//  $$ |  $$ |      $$ |  $$\ $$ |  $$ |                    $$ |     $$ |  $$ |$$ |  $$ |$$ |  $$ |$$ |      $$ |  \$$$  /  $$   ____|$$ |  $$ |$$ |  $$ |$$ |  $$ |$$  __$$ |$$ |$$ |  $$ | \____$$\ $$ | \____$$       //
//$$$$$$\ $$ |      \$$$$$$  |$$ |  $$ |                    $$$$$$$$\\$$$$$$  |$$ |  $$ |\$$$$$$$ |$$$$$$$$\ $$ |   \$  /   \$$$$$$$\ \$$$$$$$ |$$ |  $$ |$$ |  $$ |\$$$$$$$ |$$ |\$$$$$$$ |$$$$$$$  |$$ |$$$$$$$  |     //
//\______|\__|       \______/ \__|  \__|                    \________|\______/ \__|  \__| \____$$ |\________|\__|    \_/     \_______| \_______|\__|  \__|\__|  \__| \_______|\__| \____$$ |\_______/ \__|\_______/      //
//                                                                                       $$\   $$ |                                                                               $$\   $$ |                             // 
//                                                                                       \$$$$$$  |                                                                               \$$$$$$  |                             //
//=======================================================================================================================================================================================================================//
//                                                                                                                                                                                                                       //
// Authors of the code: Celia Fernandez Madrazo                                                                                                                                                                          //
//                      Pablo Martinez Ruiz Del Arbol                                                                                                                                                                    //
//                                                                                                                                                                                                                       //
//=======================================================================================================================================================================================================================//
//                                                                                                                                                                                                                       //
// Description: Main analyzer                                                                                                                                                                                            //
//                                                                                                                                                                                                                       //
//=======================================================================================================================================================================================================================//



#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepPDT/ParticleDataTable.hh"
#include "SimGeneral/HepPDTRecord/interface/PDTRecord.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#include "TLorentzVector.h"
#include "TTree.h"
#include "TFile.h"
#include "TVector3.h"

//=======================================================================================================================================================================================================================//


///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// FUNCTIONS ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////




//=======================================================================================================================================================================================================================//
class LongLivedAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit LongLivedAnalysis(const edm::ParameterSet&);
      ~LongLivedAnalysis();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      void propagate(TLorentzVector *, TLorentzVector *, TLorentzVector *, TLorentzVector *, int);

      edm::ParameterSet parameters;
      edm::EDGetTokenT<edm::HepMCProduct> theHEPMC;   
};
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
LongLivedAnalysis::LongLivedAnalysis(const edm::ParameterSet& iConfig)
{


    parameters = iConfig;

    theHEPMC = consumes<edm::HepMCProduct> (parameters.getParameter<edm::InputTag>("HepMCProduct"));

}
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
LongLivedAnalysis::~LongLivedAnalysis()
{

}
//=======================================================================================================================================================================================================================//



//=======================================================================================================================================================================================================================//
void LongLivedAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   /////////////////////////////////////////////////////////////////////////////////////
   ///////////////////////////////////// MAIN CODE /////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////////



   //////////////////////////////// GET THE COLLECTIONS ////////////////////////////////

   edm::Handle<edm::HepMCProduct> hepProduct;
   iEvent.getByToken(theHEPMC, hepProduct);

   edm::ESHandle<HepPDT::ParticleDataTable> fTable;
   iSetup.get<PDTRecord>().get(fTable);
   const HepPDT::ParticleDataTable* pdt = &(*fTable);


   const HepMC::GenEvent* Evt = hepProduct->GetEvent() ;

   for (HepMC::GenEvent::vertex_const_iterator itVtx=Evt->vertices_begin(); itVtx!=Evt->vertices_end(); ++itVtx) {
       for (HepMC::GenVertex::particles_out_const_iterator itPartOut=(*itVtx)->particles_out_const_begin();
            itPartOut!=(*itVtx)->particles_out_const_end(); ++itPartOut) {
           HepMC::GenVertex *pro_vertex = (*itPartOut)->production_vertex();
           int charge = pdt->particle((*itPartOut)->pdg_id())->charge(); 
           if((*itPartOut)->status() == 1) {
               TLorentzVector inputPosition;
               inputPosition.SetXYZT(pro_vertex->position().x(), pro_vertex->position().y(), pro_vertex->position().z(), pro_vertex->position().t());
               TLorentzVector inputMomentum;
               inputPosition.SetXYZT((*itPartOut)->momentum().x(), (*itPartOut)->momentum().y(), (*itPartOut)->momentum().z(), (*itPartOut)->momentum().t());
               TLorentzVector outputPosition, outputMomentum;
               if(charge != 0) 
               propagate(&inputPosition, &inputMomentum, &outputPosition, &outputMomentum, charge);               

           }            
       }
   }


}
//=======================================================================================================================================================================================================================//


void LongLivedAnalysis::propagate(TLorentzVector *inputPosition, TLorentzVector *inputMomentum, TLorentzVector *outputPosition, TLorentzVector *outputMomentum, int charge) {

 //Constants
    double B = 3.8;
    double PositiveEndcap_length = 3.0;
    double NegativeEndcap_length = -3.0;
    double R1 = 1.25;
    double c_light = 299792458.0;

    double m = inputMomentum->M();
    double E = inputMomentum->E();
    double gamma = inputMomentum->Gamma();
    double pt = inputMomentum->Pt();
    double vt = (pt / (m * gamma)) * c_light;
    double pz = inputMomentum->Pz();
    double vz = (pz / (m * gamma)) * c_light;

    double qoverm = 8.9849e+7 / m; //m in GeV
    double w = charge * qoverm * (B / gamma);
    double delta = atan2(inputMomentum->X(), inputMomentum->Y());
    double x0 = inputPosition->X() + vt/w * cos(delta);
    double y0 = inputPosition->Y() - vt/w * sin(delta);
    double z0 = inputPosition->Z();
    double R2 = fabs(vt / w);


    double tposZ = 1e15;
    double tnegZ = 1e15;

   if(vz != 0) {
        //Cut with positive Z endcap
        tposZ = (PositiveEndcap_length - z0) / vz;
        if(tposZ < 0) tposZ = 1e15;
        //Cut with positive Z endcap
        tnegZ = (NegativeEndcap_length - z0) / vz;
        if(tnegZ < 0) tnegZ = 1e15;
    }
    //Cut with barrel 
    double Delta = R1 * R1 - R2 * R2 + x0 * x0 + y0 * y0;
    double Ac = 4.0 * (x0 * x0 + y0 * y0);
    double Bc = -4.0 * Delta * x0;
    double Cc = Delta * Delta - 4.0 * y0 * y0 * R1 * R1;
    double radicant = Bc*Bc - 4.0* Ac * Cc;
    double tBarrel = 1e15;
    if(radicant >= 0) {
        double x1 = (-Bc + sqrt(radicant))/(2.0*Ac);
        double x2 = (-Bc - sqrt(radicant))/(2.0*Ac);

        double cos1 = -(x1 - x0) * w / vt;
        double cos2 = -(x2 - x0) * w / vt;
        double phase1p = acos(cos1);
        double phase2p = acos(cos2);
        double phase1m = acos(cos1);
        double phase2m = acos(cos2);

        double t1p = (phase1p-delta)/w;
        if(t1p < 0) {
            phase1p = phase1p + 2.0 * TMath::Pi();
            t1p = (phase1p-delta)/w;
        }
        if(t1p < 0) {
            phase1p = phase1p - 4.0 * TMath::Pi();
            t1p = (phase1p-delta)/w;
        }
                                                       
        double t1m = (-phase1m-delta)/w;
        if(t1m < 0) {
            phase1m = phase1m + 2.0 * TMath::Pi();
            t1m = (phase1m-delta)/w;
        }
        if(t1m < 0) {
            phase1m = phase1m - 4.0 * TMath::Pi();
            t1m = (phase1m-delta)/w;
        }

        double t2p = (phase2p-delta)/w;
        if(t2p < 0) {
            phase2p = phase2p + 2.0 * TMath::Pi();
            t2p = (phase2p-delta)/w;
        }
        if(t2p < 0) {
            phase2p = phase2p - 4.0 * TMath::Pi();
            t2p = (phase2p-delta)/w;
        }
        double t2m = (-phase2m-delta)/w;
        if(t2m < 0) {
            phase2m = phase2m + 2.0 * TMath::Pi();
            t2m = (phase2m-delta)/w;
        }
        if(t2m < 0) {
            phase2m = phase2m - 4.0 * TMath::Pi();
            t2m = (phase2m-delta)/w;
        }


        double x1p = - vt / w * cos(w * t1p + delta) + x0;
        double y1p =   vt / w * sin(w * t1p + delta) + y0;
        double diff1p = x1p*x1p + y1p*y1p - R1*R1;

        double x2p = - vt / w * cos(w * t2p + delta) + x0;
        double y2p =   vt / w * sin(w * t2p + delta) + y0;
        double diff2p = x2p*x2p + y2p*y2p - R1*R1;

        double t1, t2;
        if(fabs(diff1p) > 0.00001) {
            t1 = t1m;
        } else {
            t1 = t1p;
        }
        if(fabs(diff2p) > 0.00001) {
            t2 = t2m;
        } else {
            t2 = t2p;
        }

        tBarrel = TMath::Min(t1, t2);
    }
    if(tBarrel == 1e15 && tposZ == 1e15 && tnegZ == 1e15) {
        outputPosition->SetXYZT(0, 0, 0, 0);
        outputMomentum->SetXYZT(0, 0, 0, 0);
        return;
    }
 
    double t = TMath::Min(tposZ, TMath::Min(tnegZ, tBarrel));
    double x = - vt / w * cos(w * t + delta) + x0;
    double y =   vt / w * sin(w * t + delta) + y0;
    double z =   vz * t + z0;
    double px = pt * sin(w * t + delta);
    double py = pt * cos(w * t + delta);

    outputPosition->SetXYZT(x, y, z, t + inputPosition->T());
    outputMomentum->SetXYZT(px, py, pz, E);

}



//=======================================================================================================================================================================================================================//
void LongLivedAnalysis::beginJob()
{

}
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
void LongLivedAnalysis::endJob() 
{

}
//=======================================================================================================================================================================================================================//




//=======================================================================================================================================================================================================================//
void LongLivedAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//=======================================================================================================================================================================================================================//



DEFINE_FWK_MODULE(LongLivedAnalysis);
