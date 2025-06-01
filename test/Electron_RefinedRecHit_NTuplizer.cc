#include "../interface/Electron_RefinedRecHit_NTuplizer.h"

//
// constructors and destructor
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using namespace edm;
using namespace std;
using namespace reco;


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Electron_RefinedRecHit_NTuplizer::Electron_RefinedRecHit_NTuplizer(const edm::ParameterSet& iConfig):
   rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rhoFastJet"))),
   recHitCollectionEBToken_(consumes<EcalRecHitCollection>(edm::InputTag("reducedEcalRecHitsEB"))),
   recHitCollectionEEToken_(consumes<EcalRecHitCollection>(edm::InputTag("reducedEcalRecHitsEE"))),
   recHitCollectionESToken_(consumes<EcalRecHitCollection>(edm::InputTag("reducedEcalRecHitsES"))),
   eleLooseIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"))),
   eleMediumIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"))),
   eleTightIdMapToken_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("eleTightIdMap")))
{
   //now do what ever initialization is needed
   isMC_ = iConfig.getParameter<bool>("isMC");
   miniAODRun_ = iConfig.getParameter<bool>("miniAODRun");
   refinedCluster_ = iConfig.getParameter<bool>("refinedCluster");
   electronsToken_ = mayConsume<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electrons"));
   genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"));
   usesResource("TFileService");

   isRefinedSC = refinedCluster_; // store the type of SC used
}


Electron_RefinedRecHit_NTuplizer::~Electron_RefinedRecHit_NTuplizer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)




}


//
// member functions
//

// ------------ method called for each event  ------------
   void
Electron_RefinedRecHit_NTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   iEvent.getByToken(rhoToken_, rhoHandle);
   iEvent.getByToken(recHitCollectionEBToken_, EBRechitsHandle);
   iEvent.getByToken(recHitCollectionEEToken_, EERechitsHandle);
   iEvent.getByToken(recHitCollectionESToken_, ESRechitsHandle);
   iEvent.getByToken(electronsToken_, electrons);
   iEvent.getByToken(genParticlesToken_, genParticles);
   iEvent.getByToken(eleLooseIdMapToken_, loose_id_decisions);
   iEvent.getByToken(eleMediumIdMapToken_, medium_id_decisions);
   iEvent.getByToken(eleTightIdMapToken_ , tight_id_decisions);


   ESHandle<CaloGeometry> pG;
   iSetup.get<CaloGeometryRecord>().get(pG);
   iSetup.get<EcalPedestalsRcd>().get(_ped);
   const CaloGeometry* geo = pG.product();
   const CaloSubdetectorGeometry* ecalEBGeom = static_cast<const CaloSubdetectorGeometry*>(geo->getSubdetectorGeometry(DetId::Ecal, EcalBarrel));
   const CaloSubdetectorGeometry* ecalEEGeom = static_cast<const CaloSubdetectorGeometry*>(geo->getSubdetectorGeometry(DetId::Ecal, EcalEndcap));

   /////////////////////////////PU related variables////////////////////////////////////////////////////
   rho = *rhoHandle;
   ///////////////////////////////////////////////////////////////////////////////////////////////////////


   clustertools_NoZS = new noZS::EcalClusterLazyTools(iEvent, iSetup, recHitCollectionEBToken_, recHitCollectionEEToken_);

   //Clear all vectors to be written to the tree
   ClearTreeVectors();
   run=0;  event=0;  lumi=0;

   ///////////////////////////Fill Electron/Eleton related stuff/////////////////////////////////////////////////////
   nElectrons_ = 0;
   const EcalRecHitCollection *recHitsEB = clustertools_NoZS->getEcalEBRecHitCollection();
   const EcalRecHitCollection *recHitsEE = clustertools_NoZS->getEcalEERecHitCollection();
   for (size_t i = 0; i < electrons->size(); ++i){
      if (nElectrons_ == 2) break;
      const auto ele = electrons->ptrAt(i);
      if( ele->pt() < 10 ) continue;

      SuperClusterRef sc;

      if (refinedCluster_) sc = ele->superCluster(); // use for refined, comment out for unrefined
      else{
        if (!ele->ecalDriven()) continue;
        if (ele->parentSuperCluster().isNull()) continue;
        sc = ele->parentSuperCluster();
      }

      std::vector< std::pair<DetId, float> > hitsAndFractions = sc->hitsAndFractions();

      isEB = ((*sc->seed()).hitsAndFractions().at(0).first.subdetId() == EcalBarrel);
      isEE = ((*sc->seed()).hitsAndFractions().at(0).first.subdetId() == EcalEndcap);

      EBDetId* DidEB;
      EEDetId* DidEE;

      EcalRecHitCollection::const_iterator oneHit;
      for (const auto&  detitr : hitsAndFractions) {
         if(isEB){
            DidEB = new EBDetId(detitr.first.rawId());
            DetId Did   = detitr.first.rawId();
            shared_ptr<const CaloCellGeometry> geom = ecalEBGeom->getGeometry(Did);
            oneHit = recHitsEB->find( (detitr.first) ) ;
            iEta[nElectrons_].push_back(DidEB->ieta());
            iPhi[nElectrons_].push_back(DidEB->iphi());
            Hit_Eta[nElectrons_].push_back(geom->etaPos());
            Hit_Phi[nElectrons_].push_back(geom->phiPos());
            Hit_X[nElectrons_].push_back(geom->getPosition().x());
            Hit_Y[nElectrons_].push_back(geom->getPosition().y());
            Hit_Z[nElectrons_].push_back(geom->getPosition().z());
         }
         else if(isEE){
            DidEE = new EEDetId(detitr.first.rawId());
            DetId Did   = detitr.first.rawId();
            shared_ptr<const CaloCellGeometry> geom = ecalEEGeom->getGeometry(Did);
            oneHit = recHitsEE->find( (detitr.first) ) ;
            iEta[nElectrons_].push_back(DidEE->ix());
            iPhi[nElectrons_].push_back(DidEE->iy());
            Hit_Eta[nElectrons_].push_back(geom->etaPos());
            Hit_Phi[nElectrons_].push_back(geom->phiPos());
            Hit_X[nElectrons_].push_back(geom->getPosition().x());
            Hit_Y[nElectrons_].push_back(geom->getPosition().y());
            Hit_Z[nElectrons_].push_back(geom->getPosition().z());
         }

         RecHitEn[nElectrons_].push_back(oneHit->energy());
         RecHitFrac[nElectrons_].push_back(detitr.second);
         if(oneHit->checkFlag(EcalRecHit::kGood))	RecHitQuality[nElectrons_].push_back(1);
         else RecHitQuality[nElectrons_].push_back(0);

         for (int iflag=0; iflag<EcalRecHit::kHasSwitchToGain1+1; iflag++){
            RecHitFlag_container[iflag][nElectrons_].push_back(oneHit->checkFlag(iflag));
         }

         if (DEBUG) cout<<endl<<" Reco Flags = "<<oneHit->recoFlag()<<endl;

         if(oneHit->checkFlag(EcalRecHit::kHasSwitchToGain6)) 		RecHitGain[nElectrons_].push_back(6);
         else if(oneHit->checkFlag(EcalRecHit::kHasSwitchToGain1))            RecHitGain[nElectrons_].push_back(1);
         else RecHitGain[nElectrons_].push_back(12);
         HitNoise[nElectrons_].push_back(_ped->find(detitr.first)->rms(1));
      }  

      if(isEE){
         GetESPlaneRecHits(*sc, geo, nElectrons_, 1);     
         GetESPlaneRecHits(*sc, geo, nElectrons_, 2);
      }

      nElectrons_++;
      Ele_pt_.push_back( ele->pt() );
      Ele_eta_.push_back( ele->superCluster()->eta() );
      Ele_phi_.push_back( ele->superCluster()->phi() );
      Ele_energy_.push_back( ele->energy() );
      Ele_energy_error_.push_back( ele->correctedEcalEnergyError() );
      Ele_ecal_mustache_energy_.push_back( ele->correctedEcalEnergy() );
      Ele_R9.push_back(ele->full5x5_r9());
      Ele_SigIEIE.push_back(ele->full5x5_sigmaIetaIeta());
      Ele_SigIPhiIPhi.push_back(ele->full5x5_sigmaIphiIphi());
      Ele_SCEtaW.push_back(ele->superCluster()->etaWidth());
      Ele_SCPhiW.push_back(ele->superCluster()->phiWidth());
      Ele_HadOverEm.push_back(ele->hadronicOverEm());
      const CaloClusterPtr seed_clu = ele->superCluster()->seed();
      //        if (!seed_clu) continue;
      //        Ele_CovIEtaIEta.push_back(clustertools_NoZS->localCovariances(*seed_clu)[0]);
      //        Ele_CovIEtaIPhi.push_back(clustertools_NoZS->localCovariances(*seed_clu)[1]);
      //	Ele_ESSigRR.push_back(clustertools->eseffsirir( *(ele->superCluster()) ) );
      Ele_SCRawE.push_back(ele->superCluster()->rawEnergy());
      Ele_SC_ESEnByRawE.push_back( (ele->superCluster()->preshowerEnergy())/(ele->superCluster()->rawEnergy()) );
      //        Ele_S4.push_back(clustertools_NoZS->e2x2( *seed_clu ) / clustertools_NoZS->e5x5( *seed_clu ) );

      // Fill Isolation variables
      Ele_sumChargedHadronPt.push_back(ele->pfIsolationVariables().sumChargedHadronPt);
      Ele_sumChargedParticlePt.push_back(ele->pfIsolationVariables().sumChargedParticlePt);
      Ele_sumEcalClusterEt.push_back(ele->pfIsolationVariables().sumEcalClusterEt);
      Ele_sumHcalClusterEt.push_back(ele->pfIsolationVariables().sumHcalClusterEt);
      Ele_sumNeutralHadronEt.push_back(ele->pfIsolationVariables().sumNeutralHadronEt);
      Ele_sumPhotonEt.push_back(ele->pfIsolationVariables().sumPhotonEt);
      Ele_sumPUPt.push_back(ele->pfIsolationVariables().sumPUPt);
      Ele_EcalPFClusterIso.push_back(ele->ecalPFClusterIso());
      Ele_HcalPFClusterIso.push_back(ele->hcalPFClusterIso());

      ////// Look up and save the ID decisions
      bool isPassLoose  = (*loose_id_decisions)[ele];
      bool isPassMedium = (*medium_id_decisions)[ele];
      bool isPassTight  = (*tight_id_decisions)[ele];
      passLooseId_.push_back( (int) isPassLoose);
      passMediumId_.push_back( (int) isPassMedium);
      passTightId_.push_back ( (int) isPassTight );

   }

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////


   //////////////////////// Gen Stuff hardcaded for status 1 electrons for now /////////////////////////////////////
   if (isMC_){
      for(edm::View<GenParticle>::const_iterator part = genParticles->begin(); part != genParticles->end(); ++part){
         if( part->status()==1  && abs(part->pdgId())==11 ){
            Ele_Gen_Pt.push_back(part->pt());
            Ele_Gen_Eta.push_back(part->eta());
            Ele_Gen_Phi.push_back(part->phi());
            Ele_Gen_E.push_back(part->energy());
         }
      }
   }

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   /////////////////////////Run, event, lumi//////////////////////////////////
   run=iEvent.id().run();
   event=iEvent.id().event();
   lumi=iEvent.luminosityBlock();
   ///////////////////////////////////////////////////////////////////////////



   T->Fill(); // Write out the events
   delete clustertools_NoZS;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
   void
Electron_RefinedRecHit_NTuplizer::beginJob()
{
   edm::Service<TFileService> fs;
   T=fs->make<TTree>("T","MyTuple");


   T->Branch("iEtaEle1"  ,  &(iEta[0]));
   T->Branch("iPhiEle1"  ,  &(iPhi[0]));
   T->Branch("Hit_ES_Eta_Ele1"  ,  &(Hit_ES_Eta[0]));
   T->Branch("Hit_ES_Phi_Ele1"  ,  &(Hit_ES_Phi[0]));
   T->Branch("Hit_ES_X_Ele1"  ,  &(Hit_ES_X[0]));
   T->Branch("Hit_ES_Y_Ele1"  ,  &(Hit_ES_Y[0]));
   T->Branch("Hit_ES_Z_Ele1"  ,  &(Hit_ES_Z[0]));
   T->Branch("ES_RecHitEnEle1"  ,  &(ES_RecHitEn[0]));

   T->Branch("Hit_Eta_Ele1"  ,  &(Hit_Eta[0]));
   T->Branch("Hit_Phi_Ele1"  ,  &(Hit_Phi[0]));
   T->Branch("Hit_X_Ele1"  ,  &(Hit_X[0]));
   T->Branch("Hit_Y_Ele1"  ,  &(Hit_Y[0]));
   T->Branch("Hit_Z_Ele1"  ,  &(Hit_Z[0]));
   T->Branch("RecHitEnEle1"  ,  &(RecHitEn[0]));
   T->Branch("RecHitFracEle1"  ,  &(RecHitFrac[0]));
   T->Branch("RecHitGain1"  ,  &(RecHitGain[0]));
   T->Branch("RecHitQuality1", &(RecHitQuality[0]));
   T->Branch("HitNoiseEle1", &(HitNoise[0]));

   T->Branch("RecHitFlag_kGood_ele1", &(RecHitFlag_kGood[0]));
   T->Branch("RecHitFlag_kPoorReco_ele1", &(RecHitFlag_kPoorReco[0]));
   T->Branch("RecHitFlag_kOutOfTime_ele1", &(RecHitFlag_kOutOfTime[0]));
   T->Branch("RecHitFlag_kFaultyHardware_ele1", &(RecHitFlag_kFaultyHardware[0]));
   T->Branch("RecHitFlag_kNoisy_ele1", &(RecHitFlag_kNoisy[0]));
   T->Branch("RecHitFlag_kPoorCalib_ele1", &(RecHitFlag_kPoorCalib[0]));
   T->Branch("RecHitFlag_kSaturated_ele1", &(RecHitFlag_kSaturated[0]));
   T->Branch("RecHitFlag_kLeadingEdgeRecovered_ele1", &(RecHitFlag_kLeadingEdgeRecovered[0]));
   T->Branch("RecHitFlag_kNeighboursRecovered_ele1", &(RecHitFlag_kNeighboursRecovered[0]));
   T->Branch("RecHitFlag_kTowerRecovered_ele1", &(RecHitFlag_kTowerRecovered[0]));
   T->Branch("RecHitFlag_kDead_ele1", &(RecHitFlag_kDead[0]));
   T->Branch("RecHitFlag_kKilled_ele1", &(RecHitFlag_kKilled[0]));
   T->Branch("RecHitFlag_kTPSaturated_ele1", &(RecHitFlag_kTPSaturated[0]));
   T->Branch("RecHitFlag_kL1SpikeFlag_ele1", &(RecHitFlag_kL1SpikeFlag[0]));
   T->Branch("RecHitFlag_kWeird_ele1", &(RecHitFlag_kWeird[0]));
   T->Branch("RecHitFlag_kDiWeird_ele1", &(RecHitFlag_kDiWeird[0]));
   T->Branch("RecHitFlag_kHasSwitchToGain6_ele1", &(RecHitFlag_kHasSwitchToGain6[0]));
   T->Branch("RecHitFlag_kHasSwitchToGain1_ele1", &(RecHitFlag_kHasSwitchToGain1[0]));

   T->Branch("RecHitFlag_kESGood_ele1", &(RecHitFlag_kESGood[0]));
   T->Branch("RecHitFlag_kESDead_ele1", &(RecHitFlag_kESDead[0]));
   T->Branch("RecHitFlag_kESHot_ele1", &(RecHitFlag_kESHot[0]));
   T->Branch("RecHitFlag_kESPassBX_ele1", &(RecHitFlag_kESPassBX[0]));
   T->Branch("RecHitFlag_kESTwoGoodRatios_ele1", &(RecHitFlag_kESTwoGoodRatios[0]));
   T->Branch("RecHitFlag_kESBadRatioFor12_ele1", &(RecHitFlag_kESBadRatioFor12[0]));
   T->Branch("RecHitFlag_kESBadRatioFor23Upper_ele1", &(RecHitFlag_kESBadRatioFor23Upper[0]));
   T->Branch("RecHitFlag_kESBadRatioFor23Lower_ele1", &(RecHitFlag_kESBadRatioFor23Lower[0]));
   T->Branch("RecHitFlag_kESTS1Largest_ele1", &(RecHitFlag_kESTS1Largest[0]));
   T->Branch("RecHitFlag_kESTS3Largest_ele1", &(RecHitFlag_kESTS3Largest[0]));
   T->Branch("RecHitFlag_kESTS3Negative_ele1", &(RecHitFlag_kESTS3Negative[0]));
   T->Branch("RecHitFlag_kESSaturated_ele1", &(RecHitFlag_kESSaturated[0]));
   T->Branch("RecHitFlag_kESTS2Saturated_ele1", &(RecHitFlag_kESTS2Saturated[0]));
   T->Branch("RecHitFlag_kESTS3Saturated_ele1", &(RecHitFlag_kESTS3Saturated[0]));
   T->Branch("RecHitFlag_kESTS13Sigmas_ele1", &(RecHitFlag_kESTS13Sigmas[0]));
   T->Branch("RecHitFlag_kESTS15Sigmas_ele1", &(RecHitFlag_kESTS15Sigmas[0]));

   T->Branch("iEtaEle2"  ,  &(iEta[1]));
   T->Branch("iPhiEle2"  ,  &(iPhi[1]));
   T->Branch("Hit_ES_Eta_Ele2"  ,  &(Hit_ES_Eta[1]));
   T->Branch("Hit_ES_Phi_Ele2"  ,  &(Hit_ES_Phi[1]));
   T->Branch("Hit_ES_X_Ele2"  ,  &(Hit_ES_X[1]));
   T->Branch("Hit_ES_Y_Ele2"  ,  &(Hit_ES_Y[1]));
   T->Branch("Hit_ES_Z_Ele2"  ,  &(Hit_ES_Z[1]));
   T->Branch("ES_RecHitEnEle2"  ,  &(ES_RecHitEn[1]));

   T->Branch("Hit_Eta_Ele2"  ,  &(Hit_Eta[1]));
   T->Branch("Hit_Phi_Ele2"  ,  &(Hit_Phi[1]));
   T->Branch("Hit_X_Ele2"  ,  &(Hit_X[1]));
   T->Branch("Hit_Y_Ele2"  ,  &(Hit_Y[1]));
   T->Branch("Hit_Z_Ele2"  ,  &(Hit_Z[1]));
   T->Branch("RecHitEnEle2"  ,  &(RecHitEn[1]));
   T->Branch("RecHitFracEle2"  ,  &(RecHitFrac[1]));
   T->Branch("RecHitGain2"  ,  &(RecHitGain[1]));
   T->Branch("RecHitQuality2", &(RecHitQuality[1]));
   T->Branch("HitNoiseEle2", &(HitNoise[1]));

   T->Branch("RecHitFlag_kGood_ele2", &(RecHitFlag_kGood[1]));
   T->Branch("RecHitFlag_kPoorReco_ele2", &(RecHitFlag_kPoorReco[1]));
   T->Branch("RecHitFlag_kOutOfTime_ele2", &(RecHitFlag_kOutOfTime[1]));
   T->Branch("RecHitFlag_kFaultyHardware_ele2", &(RecHitFlag_kFaultyHardware[1]));
   T->Branch("RecHitFlag_kNoisy_ele2", &(RecHitFlag_kNoisy[1]));
   T->Branch("RecHitFlag_kPoorCalib_ele2", &(RecHitFlag_kPoorCalib[1]));
   T->Branch("RecHitFlag_kSaturated_ele2", &(RecHitFlag_kSaturated[1]));
   T->Branch("RecHitFlag_kLeadingEdgeRecovered_ele2", &(RecHitFlag_kLeadingEdgeRecovered[1]));
   T->Branch("RecHitFlag_kNeighboursRecovered_ele2", &(RecHitFlag_kNeighboursRecovered[1]));
   T->Branch("RecHitFlag_kTowerRecovered_ele2", &(RecHitFlag_kTowerRecovered[1]));
   T->Branch("RecHitFlag_kDead_ele2", &(RecHitFlag_kDead[1]));
   T->Branch("RecHitFlag_kKilled_ele2", &(RecHitFlag_kKilled[1]));
   T->Branch("RecHitFlag_kTPSaturated_ele2", &(RecHitFlag_kTPSaturated[1]));
   T->Branch("RecHitFlag_kL1SpikeFlag_ele2", &(RecHitFlag_kL1SpikeFlag[1]));
   T->Branch("RecHitFlag_kWeird_ele2", &(RecHitFlag_kWeird[1]));
   T->Branch("RecHitFlag_kDiWeird_ele2", &(RecHitFlag_kDiWeird[1]));
   T->Branch("RecHitFlag_kHasSwitchToGain6_ele2", &(RecHitFlag_kHasSwitchToGain6[1]));
   T->Branch("RecHitFlag_kHasSwitchToGain1_ele2", &(RecHitFlag_kHasSwitchToGain1[1]));

   T->Branch("RecHitFlag_kESGood_ele2", &(RecHitFlag_kESGood[1]));
   T->Branch("RecHitFlag_kESDead_ele2", &(RecHitFlag_kESDead[1]));
   T->Branch("RecHitFlag_kESHot_ele2", &(RecHitFlag_kESHot[1]));
   T->Branch("RecHitFlag_kESPassBX_ele2", &(RecHitFlag_kESPassBX[1]));
   T->Branch("RecHitFlag_kESTwoGoodRatios_ele2", &(RecHitFlag_kESTwoGoodRatios[1]));
   T->Branch("RecHitFlag_kESBadRatioFor12_ele2", &(RecHitFlag_kESBadRatioFor12[1]));
   T->Branch("RecHitFlag_kESBadRatioFor23Upper_ele2", &(RecHitFlag_kESBadRatioFor23Upper[1]));
   T->Branch("RecHitFlag_kESBadRatioFor23Lower_ele2", &(RecHitFlag_kESBadRatioFor23Lower[1]));
   T->Branch("RecHitFlag_kESTS1Largest_ele2", &(RecHitFlag_kESTS1Largest[1]));
   T->Branch("RecHitFlag_kESTS3Largest_ele2", &(RecHitFlag_kESTS3Largest[1]));
   T->Branch("RecHitFlag_kESTS3Negative_ele2", &(RecHitFlag_kESTS3Negative[1]));
   T->Branch("RecHitFlag_kESSaturated_ele2", &(RecHitFlag_kESSaturated[1]));
   T->Branch("RecHitFlag_kESTS2Saturated_ele2", &(RecHitFlag_kESTS2Saturated[1]));
   T->Branch("RecHitFlag_kESTS3Saturated_ele2", &(RecHitFlag_kESTS3Saturated[1]));
   T->Branch("RecHitFlag_kESTS13Sigmas_ele2", &(RecHitFlag_kESTS13Sigmas[1]));
   T->Branch("RecHitFlag_kESTS15Sigmas_ele2", &(RecHitFlag_kESTS15Sigmas[1]));

   T->Branch("nElectrons",  &nElectrons_ , "nEle/I");
   T->Branch("pt"  ,  &Ele_pt_);
   T->Branch("eta" ,  &Ele_eta_ );
   T->Branch("phi" ,  &Ele_phi_ );
   T->Branch("energy", &Ele_energy_);
   T->Branch("energy_error", &Ele_energy_error_);
   T->Branch("energy_ecal_mustache", &Ele_ecal_mustache_energy_);

   T->Branch("passLooseId", &passLooseId_);
   T->Branch("passMediumId" ,  &passMediumId_ );
   T->Branch("passTightId"  ,  &passTightId_ );
   //   T->Branch("passMVAMediumId", &passMVAMediumId_);

   T->Branch("Ele_R9"  ,  &Ele_R9);
   T->Branch("Ele_S4"  ,  &Ele_S4);
   T->Branch("Ele_SigIEIE"  ,  &Ele_SigIEIE);
   T->Branch("Ele_SigIPhiIPhi" , &Ele_SigIPhiIPhi);
   T->Branch("Ele_SCEtaW"  ,  &Ele_SCEtaW);
   T->Branch("Ele_SCPhiW"  ,  &Ele_SCPhiW);
   T->Branch("Ele_CovIEtaIEta"  ,  &Ele_CovIEtaIEta);
   T->Branch("Ele_CovIEtaIPhi"  ,  &Ele_CovIEtaIPhi);
   T->Branch("Ele_ESSigRR"  ,  &Ele_ESSigRR);
   T->Branch("Ele_SCRawE"  ,  &Ele_SCRawE);
   T->Branch("Ele_SC_ESEnByRawE"  ,  &Ele_SC_ESEnByRawE);
   T->Branch("Ele_HadOverEm"  ,  &Ele_HadOverEm);

   T->Branch("sumChargedHadronPt", &Ele_sumChargedHadronPt);
   T->Branch("sumChargedParticlePt", &Ele_sumChargedParticlePt);
   T->Branch("sumEcalClusterEt", &Ele_sumEcalClusterEt);
   T->Branch("sumHcalClusterEt", &Ele_sumHcalClusterEt);
   T->Branch("sumNeutralHadronEt", &Ele_sumNeutralHadronEt);
   T->Branch("sumPhotonEt", &Ele_sumPhotonEt);
   T->Branch("sumPUPt", &Ele_sumPUPt);
   T->Branch("Ele_EcalPFClusterIso"	,	&Ele_EcalPFClusterIso);
   T->Branch("Ele_HcalPFClusterIso"	,	 &Ele_HcalPFClusterIso);


   if (isMC_){
      T->Branch("Ele_Gen_Pt" , &Ele_Gen_Pt);
      T->Branch("Ele_Gen_Eta" , &Ele_Gen_Eta);
      T->Branch("Ele_Gen_Phi" , &Ele_Gen_Phi);
      T->Branch("Ele_Gen_E" , &Ele_Gen_E);
   }

   T->Branch("rho", &rho, "rho/F");

   T->Branch("run",&run,"run/I");
   T->Branch("event",&event,"event/I");
   T->Branch("lumi",&lumi,"lumi/I");
   
   T->Branch("isRefinedSC", &isRefinedSC);

}

// ------------ method called once each job just after ending the event loop  ------------
   void
Electron_RefinedRecHit_NTuplizer::endJob()
{
}

/*
//Evaluate if the gen particle dR matched to a reco electron is also a electron
bool Electron_RefinedRecHit_NTuplizer::GetGenMatchType(const reco::GsfElectron& Electron, const reco::GenParticle& GenColl, int pdgId, double dRThresh){


}
*/


// Extract the rechits of the SC from the ES layers
void Electron_RefinedRecHit_NTuplizer::GetESPlaneRecHits(const reco::SuperCluster& sc, const CaloGeometry* &geo, unsigned int elenum, unsigned int planeIndex) {
   double RawenergyPlane = 0.;
   double pfRawenergyPlane = 0.;
   //      if(!_ESRechitsHandle.isValid())
   //              return RawenergyPlane;

   //        if (!sc.preshowerClusters().isAvailable()) //protection for miniAOD
   //                break;

   int NumHits = 0;

   const CaloSubdetectorGeometry* ecalESGeom = static_cast<const CaloSubdetectorGeometry*>(geo->getSubdetectorGeometry(DetId::Ecal, EcalPreshower));


   for(auto iES = sc.preshowerClustersBegin(); iES != sc.preshowerClustersEnd(); ++iES) {//loop over preshower clusters
      const std::vector< std::pair<DetId, float> > hits = (*iES)->hitsAndFractions();
      for(std::vector<std::pair<DetId, float> >::const_iterator rh = hits.begin(); rh != hits.end(); ++rh) { // loop over recHits of the cluster
         //      std::cout << "print = " << (*iES)->printHitAndFraction(iCount);
         //      ++iCount;
         for(ESRecHitCollection::const_iterator esItr = ESRechitsHandle->begin(); esItr != ESRechitsHandle->end(); ++esItr) {//loop over ES rechits to find the one in the cluster
            ESDetId rhid = ESDetId(esItr->id());
            if(rhid == (*rh).first) { // found ESrechit
               //                                        std::cout << " ES energy = " << esItr->energy() << " pf energy = " << (*rh).second << std::endl;
               if((int) rhid.plane() == (int) planeIndex) {
                  std::shared_ptr<const CaloCellGeometry> geom = ecalESGeom->getGeometry(rhid);
                  Hit_ES_Eta[elenum].push_back( geom->etaPos() );
                  Hit_ES_Phi[elenum].push_back( geom->phiPos() );
                  Hit_ES_X[elenum].push_back( geom->getPosition().x() );
                  Hit_ES_Y[elenum].push_back( geom->getPosition().y() );
                  Hit_ES_Z[elenum].push_back( geom->getPosition().z() ) ;
                  ES_RecHitEn[elenum].push_back(esItr->energy());

                  for (int iflag=EcalRecHit::kESGood; iflag<EcalRecHit::kESTS15Sigmas+1; iflag++){
                     bool check_bit = esItr->checkFlag(iflag);
                     RecHitESFlag_container[iflag][elenum].push_back(check_bit);

                     if (DEBUG) cout<< "ES Flag: "<<iflag<<endl;
                  }
                  //						std::cout << "Preshower" <<std::setprecision(4) << " Eta = " <<geom->etaPos() << " : " <<" Phi = "<< geom->phiPos() << " 3D position" << geom->getPosition().z() << std::endl;
                  RawenergyPlane += esItr->energy();
                  pfRawenergyPlane += rh->second;
                  NumHits++;
               }
               break;
            }
         }
      }

      //		std::cout<<std::endl<<" Number of ES hits in plane 1 = "<<NumHits<<std::endl;
   }

   //       return RawenergyPlane;
}


//Clear tree vectors each time analyze method is called
void Electron_RefinedRecHit_NTuplizer::ClearTreeVectors()
{
   nElectrons_ = 0;
   iEta[0].clear();
   iPhi[0].clear();


   Hit_ES_Eta[0].clear();
   Hit_ES_Phi[0].clear();
   Hit_ES_X[0].clear();
   Hit_ES_Y[0].clear();
   Hit_ES_Z[0].clear();
   ES_RecHitEn[0].clear();


   Hit_Eta[0].clear();
   Hit_Phi[0].clear();
   Hit_X[0].clear();
   Hit_Y[0].clear();
   Hit_Z[0].clear();
   RecHitEn[0].clear();
   RecHitFrac[0].clear();
   RecHitGain[0].clear();
   RecHitQuality[0].clear();
   HitNoise[0].clear();
   iEta[1].clear();
   iPhi[1].clear();

   RecHitFlag_kGood[0].clear();
   RecHitFlag_kPoorReco[0].clear();
   RecHitFlag_kOutOfTime[0].clear();
   RecHitFlag_kFaultyHardware[0].clear();
   RecHitFlag_kNoisy[0].clear();
   RecHitFlag_kPoorCalib[0].clear();
   RecHitFlag_kSaturated[0].clear();
   RecHitFlag_kLeadingEdgeRecovered[0].clear();
   RecHitFlag_kNeighboursRecovered[0].clear();
   RecHitFlag_kTowerRecovered[0].clear();
   RecHitFlag_kDead[0].clear();
   RecHitFlag_kKilled[0].clear();
   RecHitFlag_kTPSaturated[0].clear();
   RecHitFlag_kL1SpikeFlag[0].clear();
   RecHitFlag_kWeird[0].clear();
   RecHitFlag_kDiWeird[0].clear();
   RecHitFlag_kHasSwitchToGain6[0].clear();
   RecHitFlag_kHasSwitchToGain1[0].clear();

   RecHitFlag_kESGood[0].clear();
   RecHitFlag_kESDead[0].clear();
   RecHitFlag_kESHot[0].clear();
   RecHitFlag_kESPassBX[0].clear();
   RecHitFlag_kESTwoGoodRatios[0].clear();
   RecHitFlag_kESBadRatioFor12[0].clear();
   RecHitFlag_kESBadRatioFor23Upper[0].clear();
   RecHitFlag_kESBadRatioFor23Lower[0].clear();
   RecHitFlag_kESTS1Largest[0].clear();
   RecHitFlag_kESTS3Largest[0].clear();
   RecHitFlag_kESTS3Negative[0].clear();
   RecHitFlag_kESSaturated[0].clear();
   RecHitFlag_kESTS2Saturated[0].clear();
   RecHitFlag_kESTS3Saturated[0].clear();
   RecHitFlag_kESTS13Sigmas[0].clear();
   RecHitFlag_kESTS15Sigmas[0].clear();

   Hit_ES_Eta[1].clear();
   Hit_ES_Phi[1].clear();
   Hit_ES_X[1].clear();
   Hit_ES_Y[1].clear();
   Hit_ES_Z[1].clear();
   ES_RecHitEn[1].clear();

   Hit_Eta[1].clear();
   Hit_Phi[1].clear();
   Hit_X[1].clear();
   Hit_Y[1].clear();
   Hit_Z[1].clear();
   RecHitEn[1].clear();
   RecHitFrac[1].clear();
   RecHitGain[1].clear();
   RecHitQuality[1].clear();
   HitNoise[1].clear();
   Ele_pt_.clear();
   Ele_eta_.clear();
   Ele_phi_.clear();
   Ele_energy_.clear();
   Ele_energy_error_.clear();
   Ele_ecal_mustache_energy_.clear();
   Ele_R9.clear();
   Ele_S4.clear();
   Ele_SigIEIE.clear();
   Ele_SigIPhiIPhi.clear();
   Ele_SCEtaW.clear();
   Ele_SCPhiW.clear();
   Ele_CovIEtaIEta.clear();
   Ele_CovIEtaIPhi.clear();
   Ele_ESSigRR.clear();
   Ele_SCRawE.clear();
   Ele_SC_ESEnByRawE.clear();
   Ele_HadOverEm.clear();

   Ele_sumChargedHadronPt.clear();
   Ele_sumChargedParticlePt.clear();
   Ele_sumEcalClusterEt.clear();
   Ele_sumHcalClusterEt.clear();
   Ele_sumNeutralHadronEt.clear();
   Ele_sumPhotonEt.clear();
   Ele_sumPUPt.clear();
   Ele_EcalPFClusterIso.clear();
   Ele_HcalPFClusterIso.clear();

   if (isMC_){
      Ele_Gen_Pt.clear();
      Ele_Gen_Eta.clear();
      Ele_Gen_Phi.clear();
      Ele_Gen_E.clear();
   }

   RecHitFlag_kGood[1].clear();
   RecHitFlag_kPoorReco[1].clear();
   RecHitFlag_kOutOfTime[1].clear();
   RecHitFlag_kFaultyHardware[1].clear();
   RecHitFlag_kNoisy[1].clear();
   RecHitFlag_kPoorCalib[1].clear();
   RecHitFlag_kSaturated[1].clear();
   RecHitFlag_kLeadingEdgeRecovered[1].clear();
   RecHitFlag_kNeighboursRecovered[1].clear();
   RecHitFlag_kTowerRecovered[1].clear();
   RecHitFlag_kDead[1].clear();
   RecHitFlag_kKilled[1].clear();
   RecHitFlag_kTPSaturated[1].clear();
   RecHitFlag_kL1SpikeFlag[1].clear();
   RecHitFlag_kWeird[1].clear();
   RecHitFlag_kDiWeird[1].clear();
   RecHitFlag_kHasSwitchToGain6[1].clear();
   RecHitFlag_kHasSwitchToGain1[1].clear();

   RecHitFlag_kESGood[1].clear();
   RecHitFlag_kESDead[1].clear();
   RecHitFlag_kESHot[1].clear();
   RecHitFlag_kESPassBX[1].clear();
   RecHitFlag_kESTwoGoodRatios[1].clear();
   RecHitFlag_kESBadRatioFor12[1].clear();
   RecHitFlag_kESBadRatioFor23Upper[1].clear();
   RecHitFlag_kESBadRatioFor23Lower[1].clear();
   RecHitFlag_kESTS1Largest[1].clear();
   RecHitFlag_kESTS3Largest[1].clear();
   RecHitFlag_kESTS3Negative[1].clear();
   RecHitFlag_kESSaturated[1].clear();
   RecHitFlag_kESTS2Saturated[1].clear();
   RecHitFlag_kESTS3Saturated[1].clear();
   RecHitFlag_kESTS13Sigmas[1].clear();
   RecHitFlag_kESTS15Sigmas[1].clear();

   Ele_pt_.clear();
   Ele_eta_.clear();
   Ele_phi_.clear();
   Ele_energy_.clear();
   Ele_ecal_mustache_energy_.clear();
   Ele_R9.clear();
   Ele_S4.clear();
   Ele_SigIEIE.clear();
   Ele_SigIPhiIPhi.clear();
   Ele_SCEtaW.clear();
   Ele_SCPhiW.clear();
   Ele_CovIEtaIEta.clear();
   Ele_CovIEtaIPhi.clear();
   Ele_ESSigRR.clear();
   Ele_SCRawE.clear();
   Ele_SC_ESEnByRawE.clear();
   Ele_HadOverEm.clear();

   Ele_Gen_Pt.clear();
   Ele_Gen_Eta.clear();
   Ele_Gen_Phi.clear();
   Ele_Gen_E.clear();

   passLooseId_.clear();
   passMediumId_.clear();
   passTightId_ .clear();
   passMVAMediumId_.clear();

   isTrue_.clear();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Electron_RefinedRecHit_NTuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
   //The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);

   //Specify that only 'tracks' is allowed
   //To use, remove the default given above and uncomment below
   //ParameterSetDescription desc;
   //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
   //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Electron_RefinedRecHit_NTuplizer);
