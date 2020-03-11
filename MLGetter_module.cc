////////////////////////////////////////////////////////////////////////
// Class:       MLGetter
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h" 
#include "art/Framework/Principal/Handle.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h" 
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "messagefacility/MessageLogger/MessageLogger.h" 

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
//#include "larcore/Geometry/geo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "larreco/RecoAlg/SpacePointAlg.h"
// #include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

#include "TVector3.h"

//#include <typeinfo>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iterator>

#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventory.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

namespace protoDUNE
{
  class MLGetterAna : public art::EDAnalyzer
  {
    public:
  
      explicit MLGetterAna(fhicl::ParameterSet const & pset);
  
      void analyze(const art::Event & evt) override;
      
      void beginJob();
      //void eventNumberFile(std::string const& eventID);
    
    private:
      std::string HitModuleLabel;
      std::string PFParticleTrackModuleLabel;
//       std::string PFParticleShowerHandle;
      std::string TrackModuleLabel;
      std::string MCModuleLabel;
      std::string ShowerModuleLabel;
      std::string FileName;
      
//       double t0;
      
  };
}


protoDUNE::MLGetterAna::MLGetterAna(fhicl::ParameterSet const & pset):
  EDAnalyzer(pset)
{
  /*this->reconfigure(pset);
  FileName (pset.get< std::string >("MCModuleLabel"))*/
}

/*void protDUNE::reconfigure(fhicl::ParameterSet const& pset)
{
  return;
}*/

void protoDUNE::MLGetterAna::beginJob()
{
}

bool testIfEmpty(std::ifstream& pFile)
{
    return pFile.peek() == std::ifstream::traits_type::eof();
}

void protoDUNE::MLGetterAna::analyze(const art::Event & evt)
{
//  std::ofstream tpcsFile;
//  tpcsFile.open("TPCS_1x2x6.txt",std::ios_base::app);
//  tpcsFile << "This file should contain the coordinates of the tpcs for the 1x2x6 detector geometry" << std::endl;
  
//  auto const* geo = lar::providerFrom<geo::Geometry>();
//  size_t nTPCs{geo->TotalNTPC()};
//  for (size_t i{0}; i < nTPCs ; i++) {
//    tpcsFile << "TPC " << i << ": Xmax: " << geo->TPC(i).MaxX() << ", Xmin: " << geo->TPC(i).MinX() << ", Ymax: " << geo->TPC(i).MaxY() << ", Ymin: " << geo->TPC(i).MinY() << ", Zmax: " << geo->TPC(i).MaxZ() << ", Zmin: " << geo->TPC(i).MinZ() << std::endl;
//  }
//  const calo::CalorimetryAlg fCaloAlg;
  
  std::string runNum    = std::to_string( evt.run()    );
  std::string subRunNum = std::to_string( evt.subRun() );
  std::string evtNum    = std::to_string( evt.event()  );
  
  MCModuleLabel         = "generator";
  PFParticleTrackModuleLabel  = "pandora";
  TrackModuleLabel            = "pandoraTrack";
  ShowerModuleLabel       = "pandoraShower";
//  PFParticleShowerHandle  = "pandora";
  HitModuleLabel        = "linecluster";

  std::ifstream emptyTrackTest;
  emptyTrackTest.open("TrackMLFile.csv");
  std::ofstream TrackMLFile;
  if ( testIfEmpty(emptyTrackTest) ){
    emptyTrackTest.close();
    TrackMLFile.open( "TrackMLFile.csv", std::ios_base::app );
    TrackMLFile << "runNum,subRunNum,evtNum,TrackID,PeakTime,SummedADC,Integral,WireID,matchedTrueTrackPDG,particleOrigin,CCNC,Mode,TrackShower" << std::endl;
  }
  else{
    emptyTrackTest.close();
    TrackMLFile.open( "TrackMLFile.csv", std::ios_base::app );
  }

//   std::ifstream emptyShowerTest;
//   emptyShowerTest.open("ShowerMLFile.csv");
//   std::ofstream ShowerMLFile;
//   if ( testIfEmpty(emptyShowerTest) ){
//     emptyShowerTest.close();
//     ShowerMLFile.open( "ShowerMLFile.csv", std::ios_base::app );
//     ShowerMLFile << "runNum,subRunNum,evtNum,TrackID,PeakTime,SummedADC,Integral,WireID,matchedTrueTrackPDG,particleOrigin" << std::endl;
//   }
//   else{
//     emptyShowerTest.close();
//     ShowerMLFile.open( "ShowerMLFile.csv", std::ios_base::app );
//   }

  /* ====================================================
   * ============== Recover the handles =================
   * ==================================================*/

  // ========= Hits =========
  art::Handle< std::vector<recob::Hit> > hitHandle;
  std::vector< art::Ptr<recob::Hit> > hitList;
  if ( evt.getByLabel( HitModuleLabel, hitHandle) ) {
    art::fill_ptr_vector( hitList, hitHandle );
  }

  // ========= Tracks =========
  art::Handle< std::vector<recob::Track> > trackHandle;
  std::vector< art::Ptr<recob::Track> > trackList;
  if ( evt.getByLabel( TrackModuleLabel, trackHandle ) ) {
    art::fill_ptr_vector( trackList, trackHandle );
  }

  // ========= Showers =========
  art::Handle< std::vector<recob::Shower> > showerHandle;
  std::vector< art::Ptr<recob::Shower> > showerList;
  if ( evt.getByLabel( ShowerModuleLabel, showerHandle ) ) {
    art::fill_ptr_vector( showerList, showerHandle );
  }

  // ========= PFParticlesTrack =========
  art::Handle< std::vector<recob::PFParticle> > PFParticleTrackHandle;
  std::vector< art::Ptr<recob::PFParticle> > PFParticleTrackList;
  if ( evt.getByLabel( PFParticleTrackModuleLabel, PFParticleTrackHandle ) ) {
    art::fill_ptr_vector( PFParticleTrackList, PFParticleTrackHandle );
  }

//   // ========= PFParticlesShower =========
//   art::Handle< std::vector<recob::PFParticle> > PFParticleShowerHandle;
//   std::vector< art::Ptr<recob::PFParticle> > PFParticleShowerList;
//   if ( evt.getByLabel( PFParticleShowerModuleLabel, PFParticleShowerHandle ) ) {
//     art::fill_ptr_vector( PFParticleShowerList, PFParticleShowerHandle );
//   }

  // ========= MC Truth =========
  art::Handle< std::vector<simb::MCTruth> > mcHandle;
  std::vector< art::Ptr<simb::MCTruth> > mcList;
  if ( evt.getByLabel( MCModuleLabel, mcHandle ) ) {
    art::fill_ptr_vector( mcList, mcHandle );
  }
  /* ====================================================
   * ==================================================*/

  /* ====================================================
  * =============== Get Associactions ==================
  * ==================================================*/

  art::FindManyP< recob::Track      > PFParticleToTrackAssns ( PFParticleTrackHandle,  evt,  TrackModuleLabel    );
  art::FindManyP< recob::Shower     > PFParticleToShowerAsns ( PFParticleTrackHandle,  evt,  ShowerModuleLabel   );
  art::FindManyP< recob::Hit        > TrackToHits            ( trackHandle,       evt,  TrackModuleLabel    );
  art::FindManyP< recob::Hit        > ShowerToHits           ( showerHandle,      evt,  ShowerModuleLabel   );
  art::FindManyP< recob::SpacePoint > HitToSpoints           ( hitHandle,         evt,  TrackModuleLabel    );
  art::FindManyP< recob::SpacePoint > HitToSpointsShowers    ( hitHandle,         evt,  ShowerModuleLabel   );

  /* ====================================================
  * ==================================================*/

  // Get the T0 of the event
//   auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
//   t0 = detprop->TriggerOffset();

  // Obtain Track information
  std::cout << "run: " << runNum << " | subrun: " << subRunNum << " | event: " << evtNum << " | trackList.size(): " << trackList.size() << " | showerList.size(): " << showerList.size() << std::endl;
  for ( size_t i{0} ; i < trackList.size() ; i++ ) { // Loop over the list of tracks

    if ( TrackToHits.isValid() ) { // Check that the associaction is valid

      std::vector< art::Ptr<recob::Hit> > vhit = TrackToHits.at(i); // Create a vector of the hits associated to the track

      art::ServiceHandle<cheat::BackTrackerService> bt;
      art::ServiceHandle<cheat::ParticleInventoryService> pi;
      std::map<int,double> trkIDEnergySum;
      bool TrackIDsSize = true;

      for ( unsigned int p = 0; p < vhit.size(); p++) {
        std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackIDEs(vhit[p]);
        if(TrackIDs.size() == 0) {
          TrackIDsSize = false;
        }
        for(size_t k = 0; k < TrackIDs.size(); k++) {
          trkIDEnergySum[TrackIDs[k].trackID] += TrackIDs[k].energy;
        }
      }

      double maxEnergy = -999.9;
      int trkMostEnergyID;

      for(std::map<int,double>::iterator trkIDIt = trkIDEnergySum.begin(); trkIDIt!=trkIDEnergySum.end(); ++trkIDIt) {
        if(trkIDIt->second > maxEnergy) {
          maxEnergy = trkIDIt->second;
          trkMostEnergyID = trkIDIt->first;
        }
      }
      art::Ptr<simb::MCTruth> mc;
      int particleOrigin = 4000;
      int CCNC = 4000;
      int mode = 4000;
      int matchedTrueTrackPDG;
      if(TrackIDsSize) {
        matchedTrueTrackPDG = pi->TrackIdToParticle(trkMostEnergyID).PdgCode();
        mc = pi->TrackIdToMCTruth_P(trkMostEnergyID);
        particleOrigin = mc->Origin();
        CCNC = mc->GetNeutrino().CCNC(); // 0 for CC, 1 for NC
        mode = mc->GetNeutrino().Mode(); // 0 for QE, 1 for Res, 2 for DIS, 3 for Coherent
      }
      else {
        matchedTrueTrackPDG   = 0;
      }

      for ( size_t h = 0 ; h < vhit.size() ; h++ ) {
        TrackMLFile << runNum << "," << subRunNum << "," << evtNum << "," << (i+1) << "," << vhit[h]->PeakTime() << "," << vhit[h]->SummedADC() << "," << vhit[h]->Integral() << "," << vhit[h]->WireID() << "," << matchedTrueTrackPDG << "," << particleOrigin << "," << CCNC << "," << mode << "," << 1 << std::endl;
      }
    }
  }

// Shower space points
  for ( size_t i{0} ; i < showerList.size() ; i++ ) { // Loop over list of showers

    if ( ShowerToHits.isValid() ) { // Check association is valid

      std::vector< art::Ptr<recob::Hit> > vhit = ShowerToHits.at(i); // Create a vector of the hits associated to the shower

      art::ServiceHandle<cheat::BackTrackerService> bt;
      art::ServiceHandle<cheat::ParticleInventoryService> pi;
      std::map<int,double> trkIDEnergySum;
      bool TrackIDsSize = true;

      for ( unsigned int p = 0; p < vhit.size(); p++) {
        std::vector<sim::TrackIDE> TrackIDs = bt->HitToTrackIDEs(vhit[p]);
        if(TrackIDs.size() == 0) {
          TrackIDsSize = false;
        }
        for(size_t k = 0; k < TrackIDs.size(); k++) {
          trkIDEnergySum[TrackIDs[k].trackID] += TrackIDs[k].energy;
        }
      }

      double maxEnergy = -999.9;
      int trkMostEnergyID;

      for(std::map<int,double>::iterator trkIDIt = trkIDEnergySum.begin(); trkIDIt!=trkIDEnergySum.end(); ++trkIDIt) {
        if(trkIDIt->second > maxEnergy) {
          maxEnergy = trkIDIt->second;
          trkMostEnergyID = trkIDIt->first;
        }
      }

      art::Ptr<simb::MCTruth> mc;
      int particleOrigin = 4000;
      int CCNC = 4000;
      int mode = 4000;
      int matchedTrueTrackPDG = 0;
      if(TrackIDsSize) {
        matchedTrueTrackPDG  = pi->TrackIdToParticle(trkMostEnergyID).PdgCode();
        mc = pi->TrackIdToMCTruth_P(trkMostEnergyID);
        particleOrigin = mc->Origin();
        CCNC = mc->GetNeutrino().CCNC(); // 0 for CC, 1 for NC
        mode = mc->GetNeutrino().Mode(); // 0 for QE, 1 for Res, 2 for DIS, 3 for Coherent
      }

      for ( size_t h{0} ; h < vhit.size() ; h++ ) {
        TrackMLFile << runNum << "," << subRunNum << "," << evtNum << "," << (i+1) << "," << vhit[h]->PeakTime() << "," << vhit[h]->SummedADC() << "," << vhit[h]->Integral() << "," << vhit[h]->WireID() << "," << matchedTrueTrackPDG << "," << particleOrigin << "," << CCNC << "," << mode << "," << 0 << std::endl;
//         std::vector< art::Ptr<recob::SpacePoint> > spts = HitToSpointsShowers.at(vhit[h].key());
//         if ( spts.size() > 0 ) {
//           ShowerMLFile << runNum << "," << subRunNum << "," << evtNum << "," << (i+1) << "," << vhit[h]->PeakTime() << "," << vhit[h]->SummedADC() << "," << vhit[h]->Integral() << "," << vhit[h]->WireID() << matchedTrueTrackPDG << "," << matchedTrueTrackPDG << "," << particleOrigin << std::endl;
//         }
      }
    }
  } // ======= End Spacepoints loops =======

  return;

} // end of analyze
//////////////////////////////////////////////////////////////////////////////////////////////////////

DEFINE_ART_MODULE(protoDUNE::MLGetterAna)


