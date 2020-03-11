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
#include "larcore/Geometry/geo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "larreco/RecoAlg/SpacePointAlg.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
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
      std::string PFParticleModuleLabel;
      std::string TrackModuleLabel;
      std::string MCModuleLabel;
      std::string ShowerModuleLabel;
      
      double t0;
      
  };
}


protoDUNE::MLGetterAna::MLGetterAna(fhicl::ParameterSet const & pset):
  EDAnalyzer(pset)
{
}

void protoDUNE::MLGetterAna::beginJob()
{
}

bool testIfEmpty(std::ifstream& pFile)
{
    return pFile.peek() == std::ifstream::traits_type::eof();
}

void protoDUNE::MLGetterAna::analyze(const art::Event & evt)
{                
  const calo::CalorimetryAlg fCaloAlg;
  
  std::string runNum    = std::to_string( evt.run()    );
  std::string subRunNum = std::to_string( evt.subRun() );
  std::string evtNum    = std::to_string( evt.event()  );
  MCModuleLabel         = "generator";  
  PFParticleModuleLabel = "pandora";
  TrackModuleLabel      = "pandora";
  ShowerModuleLabel     = "emshower";
  HitModuleLabel        = "linecluster";
  
  std::ifstream emptyTrackTest;
  emptyTrackTest.open("TrackMLFile.csv");
  std::ofstream TrackMLFile;
  if ( testIfEmpty(emptyTrackTest) ){
    emptyTrackTest.close();
    TrackMLFile.open( "TrackMLFile.csv", std::ios_base::app );
    TrackMLFile << "runNum,subRunNum,evtNum,TrackID,X,Y,Z,StartTick,EndTick,PeakTime,SigmaPeakTime,RMS,PeakAmplitude,SigmaPeakAmplitude,SummedADC,Integral,SigmaIntegral,Multiplicity,LocalIndex,GoodnessOfFit,DegreesOfFreedom,Channel" << std::endl;
  }
  else{
    emptyTrackTest.close();
    TrackMLFile.open( "TrackMLFile.csv", std::ios_base::app );
  }
  
  std::ifstream emptyShowerTest;
  emptyShowerTest.open("ShowerMLFile.csv");
  std::ofstream ShowerMLFile;
  if ( testIfEmpty(emptyShowerTest) ){
    emptyShowerTest.close();
    ShowerMLFile.open( "ShowerMLFile.csv", std::ios_base::app );
    ShowerMLFile << "runNum,subRunNum,evtNum,ShowerID,X,Y,Z,StartTick,EndTick,PeakTime,SigmaPeakTime,RMS,PeakAmplitude,SigmaPeakAmplitude,SummedADC,Integral,SigmaIntegral,Multiplicity,LocalIndex,GoodnessOfFit,DegreesOfFreedom,Channel" << std::endl;
  }
  else{
    emptyShowerTest.close();
    ShowerMLFile.open( "ShowerMLFile.csv", std::ios_base::app );
  }
  
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

  // ========= PFParticles =========
  art::Handle< std::vector<recob::PFParticle> > PFParticleHandle;
  std::vector< art::Ptr<recob::PFParticle> > PFParticleList;
  if ( evt.getByLabel( PFParticleModuleLabel, PFParticleHandle ) ) {
    art::fill_ptr_vector( PFParticleList, PFParticleHandle );
  }

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

  art::FindManyP< recob::Track      > PFParticleToTrackAssns ( PFParticleHandle,  evt,  TrackModuleLabel    );
  art::FindManyP< recob::Shower     > PFParticleToShowerAsns ( PFParticleHandle,  evt,  ShowerModuleLabel   );
  art::FindManyP< recob::Hit        > TrackToHits            ( trackHandle,       evt,  TrackModuleLabel    );
  art::FindManyP< recob::Hit        > ShowerToHits           ( showerHandle,      evt,  ShowerModuleLabel   );
  art::FindManyP< recob::SpacePoint > HitToSpoints           ( hitHandle,         evt,  TrackModuleLabel    );
  art::FindManyP< recob::SpacePoint > HitToSpointsShowers    ( hitHandle,         evt,  ShowerModuleLabel   );

  /* ====================================================
  * ==================================================*/
  
  // Get the T0 of the event
  auto const *detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
  t0 = detprop->TriggerOffset();
  
  // Obtain Track information
  for ( size_t i{0} ; i < trackList.size() ; i++ ) { // Loop over the list of tracks
    
    recob::Track track = trackList[i];
    double tracklength = track.Length()
    
    double length {0}; // Holder for length of track
    
    if ( TrackToHits.isValid() ) { // Check that the associaction is valid
      
      std::vector< art::Ptr<recob::Hit> > vhit = TrackToHits.at(i); // Create a vector of the hits associated to the track
      
      std::vector<double> Xs;
      std::vector<double> Ys;
      std::vector<double> Zs;
      
      for ( size_t h = 0 ; h < vhit.size() ; h++ ) { // Loop over the hits in the track
        
        std::vector< art::Ptr<recob::SpacePoint> > spts = HitToSpoints.at(vhit[h].key()); // Create vector of space points associated to the hits
        
        if ( spts.size() > 0 ) {
          Xs.push_back(spts[0]->XYZ()[0]);
          Ys.push_back(spts[0]->XYZ()[1]);
          Zs.push_back(spts[0]->XYZ()[2]);
        }
      }
      
      length = std::sqrt( ( Xs.back()-Xs.front() )*( Xs.back()-Xs.front() ) + ( Ys.back()-Ys.front() )*( Ys.back()-Ys.front() ) + ( Zs.back()-Zs.front() )*( Zs.back()-Zs.front() ))
      
      std::cout << "track.length = " << tracklength << std::endl;
      std::cout << "homebrew track length = " << length << std::endl << std::endl;

      for ( size_t h = 0 ; h < vhit.size() ; h++ ) {
        std::vector< art::Ptr<recob::SpacePoint> > spts = HitToSpoints.at(vhit[h].key());
        if ( spts.size() > 0 ) {
          TrackMLFile << runNum << "," << subRunNum << "," << evtNum << "," << (i+1) << "," << spts[0]->XYZ()[0] << "," << spts[0]->XYZ()[1] << "," << spts[0]->XYZ()[2] << "," << vhit[h]->StartTick() << "," << vhit[h]->EndTick() << "," << vhit[h]->PeakTime() << "," << vhit[h]->SigmaPeakTime() << "," << vhit[h]->RMS() << "," << vhit[h]->PeakAmplitude() << "," << vhit[h]->SigmaPeakAmplitude() << "," << vhit[h]->SummedADC() << "," << vhit[h]->Integral() << "," << vhit[h]->SigmaIntegral() << "," << vhit[h]->Multiplicity() << "," << vhit[h]->LocalIndex() << "," << vhit[h]->GoodnessOfFit() << "," << vhit[h]->DegreesOfFreedom() << "," << (*vhit[h]).Channel() << std::endl;
        }
      }
    }
  }
  // Shower space points
  for ( size_t i{0} ; i < showerList.size() ; i++ ) {
    if ( ShowerToHits.isValid() ) {
      std::vector< art::Ptr<recob::Hit> > vhit = ShowerToHits.at(i);
      for ( size_t h{0} ; h < vhit.size() ; h++ ) {
        std::vector< art::Ptr<recob::SpacePoint> > spts = HitToSpointsShowers.at(vhit[h].key());
        if ( spts.size() > 0 ) {
          SpacePointsFile << runNum << "," << subRunNum << "," << evtNum << "," << (i+1) << "," << spts[0]->XYZ()[0] << "," << spts[0]->XYZ()[1] << "," << spts[0]->XYZ()[2] << "," << vhit[h]->StartTick() << "," << vhit[h]->EndTick() << "," << vhit[h]->PeakTime() << "," << vhit[h]->SigmaPeakTime() << "," << vhit[h]->RMS() << "," << vhit[h]->PeakAmplitude() << "," << vhit[h]->SigmaPeakAmplitude() << "," << vhit[h]->SummedADC() << "," << vhit[h]->Integral() << "," << vhit[h]->SigmaIntegral() << "," << vhit[h]->Multiplicity() << "," << vhit[h]->LocalIndex() << "," << vhit[h]->GoodnessOfFit() << "," << vhit[h]->DegreesOfFreedom() << "," << (*vhit[h]).Channel() << std::endl;
        }
      }
    }
  } // ======= End Spacepoints loops =======
    
  return;

} // end of analyze
//////////////////////////////////////////////////////////////////////////////////////////////////////

DEFINE_ART_MODULE(protoDUNE::MLGetterAna)
