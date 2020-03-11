////////////////////////////////////////////////////////////////////////
// Class:       DeePID
//
// Author:      James Pillow, University of Warwick, for DUNE
// Email:       j.pillow@warwick.ac.uk
//
// A module to create images of particle tracks to use with a deep learning
// network to enhance the Warwick PID.
//
////////////////////////////////////////////////////////////////////////


// Art includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Larsoft includes
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
//#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/LArPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"

// Larsoft cheating includes
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventory.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

// Root includes
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TVector3.h"

// C++ includes
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iterator>



// ==================================== Class Def ====================================

namespace DeePID {
  class ImageGenAna : public art::EDAnalyzer
  {
    public:
    
    
    private:
    
    std::string TrackModuleLabel;
    std::string ShowerModuleLabel;
    std::string MCModuleLabel;
    std::string GeometryLabel;
    
    
  }; // class ImageGenAna : public art::EDAnalyzer
} // namespace DeePID +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// ==================================== Things to do ====================================

/* Things to do

 - Have function to test if supported geometry is selected
 - Make function for PD corrections
 - Make function for 1x2x6 corrections
 - Make function for FD corrections
 - Make function for producing the images
 - Make function that will create a directory to put the images into
 - Make function that can save the images
 - Make using the MC information optional through the fcl parameters

*/ // Things to do +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// ==================================== DeePID::ImageGenAna::ImageGenAna ====================================

DeePID::ImageGenAna::ImageGenAna(fhicl::ParameterSet const & pset):
  EDAnalyzer(pset),
  TrackModuleLabel  ( pset.get< std::string > ("TrackModuleLabel")  ),
  ShowerModuleLabel ( pset.get< std::string > ("ShowerModuleLabel") ),
  MCModuleLabel     ( pset.get< std::string > ("MCModuleLabel")     )
{
} // DeePID::ImageGenAna::ImageGenAna +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// ==================================== DeePID::ImageGenAna::begingJob ====================================

void DeePID::ImageGenAna::begingJob()
{
  auto const* geo  = lar::providerFrom<geo::Geometry>();
  std::string GeometryLabel = geo->DetectorName();
  std::cout << "\n==================================" << std::endl;
  std::cout << "Detector: " << detector << std::endl;
  std::cout << "==================================\n" << std::endl;
  
} // DeePID::ImageGenAna::begingJob +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// ==================================== DeePID::ImageGenAna::endJob ====================================

void DeePID::ImageGenAna::endJob()
{
} // DeePID::ImageGenAna::endJob +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// ==================================== DeePID::ImageGenAna::analyze ====================================

void DeePID::ImageGenAna::analyze(const art::Event & evt)
{
  
  // Recover our handles and associations
  
  // =============================================
  // ================== Handles ==================
  // =============================================
  
  // ========== Tracks ==========
  art::Handle< std::vector<recob::Track> >  trackHandle;
  std::vector< art::Ptr<recob::Track> >     trackList;
  if ( evt.getByLabel( TrackModuleLabel, trackHandle ) ) {
    art::fill_ptr_vector( trackList, trackHandle );
  }
  
  // ========== Showers ==========
  art::Handle< std::vector<recob::Shower> > showerHandle;
  std::vector< art::Ptr<recob::Shower> >    showerList;
  if ( evt.getByLabel( ShowerModuleLabel, showerHandle ) ) {
    art::fill_ptr_vector( showerList, showerHandle );
  }
  
  // ========== MC Truth ==========
  if ( MCModuleLabel != "" ) {
    art::Handle< std::vector<simb:MCTruth> >  mcHandle;
    std::vector< art::Ptr<simb::MCTruth> >    mcList;
    if ( evt.getByLabel( MCModuleLabel, mcHandle ) ) {
      art::fill_ptr_vector( mcList, mcHandle )
    }
  }
  
  // ========== Associations ==========
  
  art::FindManyP< recob::Hit > TrackToHits  ( trackHandle,  evt, TrackModuleLabel  );
  art::FindManyP< recob::Hit > ShowerToHits ( showerHandle, evt, ShowerModuleLabel );
  
  // ============================================
  // ================== Tracks ==================
  // ============================================
  
  for ( size_t track{0} ; traack < trackList.size() ; track++ ) { // Track Loop
  
    if ( TrackToHits.isValid() ) { // Valid Tracks
      
      // Make a vector of the hits
      std::vector< art::Ptr<recob:Hit> > hitVector = TrackToHits.at(track);
      
      // Find the MC information
      if ( MCModuleLabel != "" ) {
        std::vector<int> mcInfo = MCTruthFinder( hitVector );
      }
      
      //  Make vectors for the track information
      std::vector< int   > tpcVector;
      std::vector< int   > planeVector;
      std::vector< int   > wireVector;
      std::vector< float > adcVector;
      std::vector< float > timeVector;
      
      // Fill vectors of track information
      for ( size_t hit{0} ; hit < hitVector.size() ; hit++ ) { // Hit Loop
        
        std::string wireID = hitVector[hit]->WireID();
        int TPC    = std::stoi ( wireID.substr( wireID.find("T")+2, wireID.find("P") - wireID.find("T") - 3 ) );
        int Plane  = std::stoi ( wireID.substr( wireID.find("P")+2, wireID.find("W") - wireID.find("P") - 3 ) );
        int Wire   = std::stoi ( wireID.substr( wireID.find("W")+2, wireID.size()    - wireID.find("W") ) );
        
        // Flip U and V in certain TPCs
        if (Plane == 0) Plane = (TPC % 4 == 0 || TPC % 4 == 3) ? Plane : Plane + 1;
        if (Plane == 0) Plane = (TPC % 4 == 0 || TPC % 4 == 3) ? Plane : Plane - 1;
        
        if (Plane == 1) Wire = 1148 - wire;
        
        tpcVector.push_back   ( TPC   );
        planeVector.push_back ( Plane );
        wireVector.push_back  ( Wire  );
        adcVector.push_back   ( hitVector[ hit ]->Integral() );
        timeVector.push_back  ( hitVector[ hit ]->PeakTime() );
        
      } // Hit Loop
      
      // Correct the wires
      Corrections( tpcVector, planeVector, wireVector, timeVector, GeometryLabel );
      
      
    } // Valid Tracks
  
  } // Track Loop
  
  // Need to obtain the hit information for the event
  
  // Want this information on a plane separated basis
  
  // We want to loop over each individual track/shower
  
  // We have to apply corrections to wire numbers and to time based on the detector geometry
  
  // We then create an image for each track/shower in each plane
  
  // We then save the images
  
} // DeePID::ImageGenAna::analyze +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// ==================================== DeePID::ImageGenAna::MCTruthFinder ====================================

std::vector< int > DeePID::ImageGenAna::MCTruthFinder( const std::vector< art::Ptr<recob::Hit> > & hitVector )
{
  
  art::ServiceHandle< cheat::BackTrackerService >     backTracker;
  art:ServiceHandle< cheat::ParticleInventoryService> particleInventory;
  
  std::map<int,double> trackIDEnergySum {};
  bool TrackIDsSize {true};

  for (size_t hit{0} ; hit < hitVector.size() ; hit++ ) { // Hit Loop
    
    std::vector< sim::TrackIDE > TrackIDs = backTracker->HitToTrackIDEs( hitVector[hit] );
    
    if ( TrackIDs.size() == 0 ) {
      TrackIDsSize = false;
    }
    
    for ( size_t trackID{0} ; trackID < TrackIDs.size() ; trackID++ ) { // TrackID Loop
      trackIDEnergySum[ TrackIDs[hit].trackID ] += TrackIDs[hit].energy;
    } // TrackID
  } // Hit Loop
  
  double maxEnergy  {-999.9};
  int trackMostEnergyID {};
  
  for ( std::map<int,double>::iterator trackIDIt = trackIDEnergySum.begin() ; trackIDIt != trackIDEnergySum.end() ; trackIDIt++ ) { // TrackIDEnergySum
    if ( trackIDIt->second > maxEnergy ) {
      maxEnergy = trackIDIt->second;
      trackMostEnergyID = trackIDIt->first;
    }
  } // trackIDEnergySum Loop
  
  art::Ptr< simb::MCTruth > mc;
  int particleOrigin  {4000};
  int CCNC            {4000};
  int mode            {4000};
  int PDG             {4000};
  if ( TrackIDsSize ) {
    PDG = particleInventory->TrackIdToParticle( trackMostEnergyID ).PdgCode();
    mc = particleInventory->TrackIdToMCTruth_P( trackMostEnergyID );
    particleOrigin = mc->Origin();
    CCNC = mc->GetNeutrino().CCNC(); // 0 for CC, 1 for NC
    mode = mc->GetNeutrino().Mode(); // 0-QE, 1-Res, 2-DIS, 3-Coherent
  }
  
  std::vector< int > mcInfo = {particleOrigin, CCNC, mode, PDG};
  
  return mcInfo;
  
} // DeePID::ImageGenAna::MCTruthFinder +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// ==================================== DeePID::ImageGenAna::Corrections ====================================

void DeePID::ImageGenAna::Corrections( std::vector<int> & tpcVector, std::vector<int> & planeVector, std::vector<int> & wireVector, std::vector<float> & timeVector, const std::string & GeometryLabel )
{ // This function will corerct the wire positions for events in protoDUNE
  
  // Make corrections to time
  TimeCorrections(tpcVector, timeVector, GeometryLabel);
  
  // Make basic wire corrections based on tpc
  BasicCorrections(tpcVector, planeVector, wireVector);
  
  // Need to get min and max values of wires in each ptpcs
  auto MinMaxWiresVeec = MinMaxWires(tpcVector, planeVector, wireVector, GeometryLabel);
  
  for ( size_t entry{0} ; entry < wireVector.size() ; entry++ ){ // Loop over each wire entry
    const int tpc = tpcVector[entry];
    const int plane = planeVector[entry];
    
    if ( plane == 0 || plane == 1 ) { // Loop for plane 0 and 1
      if ( tpc == 0 ){
        wireVector[entry] += 0
      }
      else {
        wireVector[entry] -= MinMaxWiresVec[plane][tpc-1].second != -999999 ? MinMaxWiresVec[plane][tpc].first - MinMaxWiresVec[plane][tpc-1].second : 0;
      }
    } // Loop for plane 0 and 1
  } // Loop over each wire entry
  
} // DeePID::ImageGenAna::Corrections +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// ==================================== DeePID::ImageGenAna::TimeCorrections ====================================

void DeePID::ImageGenAna::TimeCorrections( const std::vector<int> & tpcVector, std::vector<float> & timeVector, const std::string & GeometryLabel )
{

  if ( GeometryLabel == "protodune" ) { // ProtoDUNE
    
    for ( int entry{0} ; entry < tpcVector.size() ; entry++ ) { // Entry Loop
      
      
    
    } // Entry Loop
    
  } // ProtoDUNE
  
  else if ( GeometryLabel == "dune10kt_v2_1x2x6" ) { // 1x2x6
  
    for ( int entry{0} ; entry < tpcVector.size() ; entry++ ) { // Entry Loop
    
    
    
    } // Entry Loop
  
  } // 1x2x6

} // DeePID::ImageGenAna::TimeCorrections +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// ==================================== DeePID::ImageGenAna::BasicCorrections ====================================

void DeePID::ImageGenAna::BasicCorrections( std::vector<int> & tpcVector, const std::vector<int> & planeVector, std::vector<int> & wireVector )
{ // Function that makes a basic naive correction to wires based only on their tpcs. Also changes to pseudoTPCs
  for (int entry{0} ; entry < wireVector.size() ; entry++){
    const int tpc = tpcVector[entry];
    const int plane = planeVector[entry];
    wireVector[entry] += plane == 2 ? std::floor(tpc/4)*480 : std::floor(tpc/4)*1148;
    tpcVector[entry] = std::floor(tpc/4);
  }
} // DeePID::ImageGenAna::BasicCorrections +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



// ==================================== DeePID::ImageGenAna::MinMaxWires ====================================

std::vector< std::map< int, std::pair<int,int> > > DeePID::ImageGenAna::MinMaxWires( const std::vector<int> & tpcVector, const std::vector<int> & planeVector, const std::vector<int> & wireVector, const std::string & GeometryLabel )
{ // This function will retrieve the min and max wires in each tpc for planes 0 and 1
  
  int NumTPCs = GeometryLabel == "protodune" ? 3 : 6; // How Many ptpcs are there to loop over?
  
  std::map<int, std::pair<int,int> > plane0; // Initialise the maps for the planes
  std::map<int, std::pair<int,int> > plane1;
  
  for ( int tpc{0} ; tpc < NumTPCs ; tpc++ ){ // Fill the maps with the number of tpcs and values for min and max wire
    plane0.emplace(tpc, std::make_pair( 999999, -999999 ) );
    plane1.emplace(tpc, std::make_pair( 999999, -999999 ) );
  }
  
  for ( size_t entry{0} ; entry < wireVector.size() ; entry++ ){ // Loop over each wire entry
    const int tpc  = tpcVector[entry];    // Retrieve the tpc and wire info
    const int wire = wireVector[entry];
    
    if ( planeVector[entry] == 0 ) { // If wire is in plane 0
      plane0[tpc].first  = wire < plane0[tpc].first ? wire : plane0[tpc].first;
      plane0[tpc].second = wire > plane0[tpc].second ? wire : plane0[tpc].second;
    }
    else if ( planeVector[entry] == 1) { // If wire is in plane 1
      plane1[tpc].first  = wire < plane1[tpc].first ? wire : plane1[tpc].first;
      plane1[tpc].second = wire > plane1[tpc].second ? wire : plane1[tpc].second;
    }
  
  } // Loop over each wire entry
  
  std::vector< std::map< int, std::pair<int,int> > > MinMaxs = {plane0,plane1}; // Zip the maps together and return them
  return MinMaxs;
  
} // DeePID::ImageGenAna::MinMaxWires +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// ==================================== DeePID::ImageGenAna:: ====================================

void DeePID::ImageGenAna::SomeKindOfPlottingAlgorithm ( const std::something & fuckKnowsWhat) {

  // Gonna need to use root to plot a 2d histogram or something similar to matplotlib's imshow

} // DeePID::ImageGenAna::SomeKindOfPlottingAlgorithm +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// ==================================== DeePID::ImageGenAna:: ====================================

// ==================================== DeePID::ImageGenAna:: ====================================

// ==================================== DeePID::ImageGenAna:: ====================================

// ==================================== DeePID::ImageGenAna:: ====================================

// ==================================== DeePID::ImageGenAna:: ====================================

// ==================================== DeePID::ImageGenAna:: ====================================

















