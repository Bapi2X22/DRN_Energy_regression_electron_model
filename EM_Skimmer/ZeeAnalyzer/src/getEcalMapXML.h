#include <TString.h>
#include <TXMLEngine.h>

using namespace reco;

void 
getEBMapXML(const edm::EventSetup& iSetup){

   ESHandle<CaloGeometry> pG;
   iSetup.get<CaloGeometryRecord>().get(pG);
   const CaloGeometry* geo = pG.product();

   TXMLEngine xml;

   const CaloSubdetectorGeometry* ecalEBGeom = static_cast<const CaloSubdetectorGeometry*>(geo->getSubdetectorGeometry(DetId::Ecal, EcalBarrel));
   //----------------------------------------------------
   // Make a lookup of detid and all the detids in R=0.3
   //----------------------------------------------------

   std::map<int, std::vector<int>> EB_neighbour_map;
   const std::vector<DetId>& ids = ecalEBGeom->getValidDetIds(DetId::Ecal, EcalBarrel);
   std::vector<int> neids;

   for (auto id : ids) {

      EBDetId ebid(id.rawId());
      std::shared_ptr<const CaloCellGeometry> xtal_igeom = ecalEBGeom->getGeometry(id);

      float ieta = xtal_igeom->etaPos();
      float iphi = xtal_igeom->phiPos();

      neids.clear();

      for ( auto jd : ids) {

         EBDetId ebjd(jd.rawId());
         std::shared_ptr<const CaloCellGeometry> xtal_jgeom = ecalEBGeom->getGeometry(jd);

         if (ebjd==ebid) continue;

         float jeta = xtal_jgeom->etaPos();
         float jphi = xtal_jgeom->phiPos();

         float dR = reco::deltaR(ieta, iphi, jeta, jphi);

         if (dR<0.3){
            neids.push_back(ebjd.numberByEtaPhi());
         }

      }

      std::sort(neids.begin(), neids.end());
      EB_neighbour_map.insert({ebid.numberByEtaPhi(), neids});
   }


   // Create main node of document treb
   XMLNodePointer_t ebnode = xml.NewChild(0, 0, "EB_xtal_dR0p3_map");

   // Add crystals
   for (std::map<int, std::vector<int>>::iterator it = EB_neighbour_map.begin();
         it != EB_neighbour_map.end(); ++it){

      int tmpxid = it->first; 
      XMLNodePointer_t tmpnode = xml.NewChild(ebnode, 0, "xtal_id", std::to_string(tmpxid).c_str());
      // Add neighbours
      for (unsigned int inn=0; inn<(it->second).size(); inn++){
         int tmpnid = (it->second)[inn];
         xml.NewChild(tmpnode, 0, "neighbour", std::to_string(tmpnid).c_str());
      }
   }

   // now create document and assign main node of document
   XMLDocPointer_t xmldoc = xml.NewDoc();
   xml.DocSetRootElement(xmldoc, ebnode);

   // Save document to file
   xml.SaveDoc(xmldoc, "EB_xtal_dR0p3_map.xml");

   // Release memory before exit
   xml.FreeDoc(xmldoc);
}



void 
getEEMapXML(const edm::EventSetup& iSetup){

   ESHandle<CaloGeometry> pG;
   iSetup.get<CaloGeometryRecord>().get(pG);
   const CaloGeometry* geo = pG.product();

   TXMLEngine xml;

   const CaloSubdetectorGeometry* ecalEEGeom = static_cast<const CaloSubdetectorGeometry*>(geo->getSubdetectorGeometry(DetId::Ecal, EcalEndcap));

   //----------------------------------------------------
   // Make a lookup of detid and all the detids in R=0.3
   //----------------------------------------------------

   std::map<int, std::vector<int>> EE_neighbour_map;
   const std::vector<DetId>& ids = ecalEEGeom->getValidDetIds(DetId::Ecal, EcalBarrel);
   std::vector<int> neids;

   for (auto id : ids) {

      EEDetId eeid(id.rawId());
      std::shared_ptr<const CaloCellGeometry> xtal_igeom = ecalEEGeom->getGeometry(id);
      int xtal_id = eeid.hashedIndex();

      float ieta = xtal_igeom->etaPos();
      float iphi = xtal_igeom->phiPos();

      neids.clear();

      for ( auto jd : ids) {

         EEDetId eejd(jd.rawId());
         std::shared_ptr<const CaloCellGeometry> xtal_jgeom = ecalEEGeom->getGeometry(jd);

         if (eejd==eeid) continue;

         float jeta = xtal_jgeom->etaPos();
         float jphi = xtal_jgeom->phiPos();

         float dR = reco::deltaR(ieta, iphi, jeta, jphi);

         if (dR<0.3){
            int xtal_jd = eejd.hashedIndex();
            neids.push_back(xtal_jd);
         }

      }

      std::sort(neids.begin(), neids.end());
      EE_neighbour_map.insert({xtal_id, neids});
   }



   // Create main node of document tree
   XMLNodePointer_t eenode = xml.NewChild(0, 0, "EE_xtal_dR0p3_map");

   // Add crystals
   for (std::map<int, std::vector<int>>::iterator it = EE_neighbour_map.begin();
         it != EE_neighbour_map.end(); ++it){

      int tmpxid = it->first; 
      XMLNodePointer_t tmpnode = xml.NewChild(eenode, 0, "xtal_id", std::to_string(tmpxid).c_str());
      // Add neighbours
      for (unsigned int inn=0; inn<(it->second).size(); inn++){
         int tmpnid = (it->second)[inn];
         xml.NewChild(tmpnode, 0, "neighbour", std::to_string(tmpnid).c_str());
      }
   }

   // now create document and assign main node of document
   XMLDocPointer_t xmldoc = xml.NewDoc();
   xml.DocSetRootElement(xmldoc, eenode);

   // Save document to file
   xml.SaveDoc(xmldoc, "EE_xtal_dR0p3_map.xml");

   // Release memory before exit
   xml.FreeDoc(xmldoc);
}
