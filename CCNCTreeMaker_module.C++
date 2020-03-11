////////////////////////////////////////////////////////////////////////
//
// file CCNCTreeMaker_module.cc
//
// j.pillow@warwick.ac.uk
//
///////////////////////////////////////////////////////////////////////

// Generic C++ includes
#include <iostream>
#include <cstddef> // std::ptrdiff_t
#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional> // std::mem_fun_ref
#include <typeinfo>
#include <memory> // std::unique_ptr<>

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/MCBase/MCStep.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/ExternalTrigger.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RawData/TriggerData.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/Deprecated/BezierTrack.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "larreco/RecoAlg/ShowerEnergyAlg.h" // Shower energy finder
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/MVAPIDResult.h"

#include "TTree.h"
#include "TTimeStamp.h"
#include "TH1D.h"
#include "TH2D.h"

namespace DeePID
{
  
  class CCNCTreeMaker : public art::EDAnalyzer
  {
    public:
    
      explicit CCNCTreeMaker(fhicl::ParameterSet const& pset);
    
      void analyze(const art::Event & evt) override;
    
      void reconfigure(fhicl::ParameterSet const& pset);
    
      std::vector< int > WireCorrector ( int wire, int plane, int tpc );
    
      void PFinfo ( art::FindManyP<recob::Cluster>             const & PFtoClust,
               art::FindManyP<recob::Cluster>             const & HitToClust,
               std::vector< art::Ptr<recob::Hit> >        const & hitList,
               std::vector< art::Ptr<recob::PFParticle> > const & PFPartList,
               bool   const & cheat );
    
      void DaughterInfo ( std::vector< art::Ptr<recob::PFParticle> > const & PFPartList,
                     art::FindManyP<recob::Cluster>             const & PFtoClust,
                     std::map< art::Ptr<recob::Cluster>, std::vector< art::Ptr<recob::Hit> > > const & clustersToHits,
                     std::vector<size_t> const & daughters,
                     std::vector< int >    & wires,
                     std::vector< int >    & planes,
                     std::vector< int >    & tpcs,
                     std::vector< double > & times,
                     std::vector< double > & integrals );
    
      void beginJob();
    
    private:
    
      // Declare root tree
      TTree* CCtree;
      TTree* NCtree;
    
      // Declare variables for root tree
      double initialE;
      int pdgCode;
      int mode;
    
      std::vector< std::vector< int > >    wireList;
      std::vector< std::vector< int > >    planeList;
      std::vector< std::vector< int > >    TPCList;
      std::vector< std::vector< double > > integralList;
      std::vector< std::vector< double > > timeList;
    
      std::vector< std::vector< int > >    wireListCheat;
      std::vector< std::vector< int > >    planeListCheat;
      std::vector< std::vector< int > >    TPCListCheat;
      std::vector< std::vector< double > > integralListCheat;
      std::vector< std::vector< double > > timeListCheat;
    
      std::string MCModuleLabel;
      std::string MCPartModLabel;
      std::string HitLabel;
      std::string pandoraLabel;
      std::string pandoraCheatLabel;
    
    
    
  }; // class CCNCTreeMaker : public art::EDAnalyzer
  
} // namespace DeePID

//---------------------------------------------------------------------------------------------------------

DeePID::CCNCTreeMaker::CCNCTreeMaker(fhicl::ParameterSet const & pset):
  EDAnalyzer(pset)
{
  this->reconfigure(pset);
} // CCNCTreeMaker::CCNCTreeMaker(fhicl::ParameterSet const& pset)

//---------------------------------------------------------------------------------------------------------

void DeePID::CCNCTreeMaker::beginJob()
{
  art::ServiceHandle<art::TFileService> tfs;
  CCtree = tfs->make<TTree>("CCtree","CCtree");
  
  CCtree->Branch( "initEnergy", &initialE );
  CCtree->Branch( "pdgCode", &pdgCode );
  CCtree->Branch( "mode", &mode );
  
  CCtree->Branch( "wireList",     &wireList     );
  CCtree->Branch( "planeList",    &planeList    );
  CCtree->Branch( "TPCList",      &TPCList      );
  CCtree->Branch( "integralList", &integralList );
  CCtree->Branch( "timeList",     &timeList     );

  CCtree->Branch( "wireListCheat",     &wireListCheat     );
  CCtree->Branch( "planeListCheat",    &planeListCheat    );
  CCtree->Branch( "TPCListCheat",      &TPCListCheat      );
  CCtree->Branch( "integralListCheat", &integralListCheat );
  CCtree->Branch( "timeListCheat",     &timeListCheat     );

  NCtree = tfs->make<TTree>("NCtree","NCtree");
  
  NCtree->Branch( "initEnergy", &initialE );
  NCtree->Branch( "pdgCode", &pdgCode );
  NCtree->Branch( "mode", &mode );
  
  NCtree->Branch( "wireList",     &wireList     );
  NCtree->Branch( "planeList",    &planeList    );
  NCtree->Branch( "TPCList",      &TPCList      );
  NCtree->Branch( "integralList", &integralList );
  NCtree->Branch( "timeList",     &timeList     );

  NCtree->Branch( "wireListCheat",     &wireListCheat     );
  NCtree->Branch( "planeListCheat",    &planeListCheat    );
  NCtree->Branch( "TPCListCheat",      &TPCListCheat      );
  NCtree->Branch( "integralListCheat", &integralListCheat );
  NCtree->Branch( "timeListCheat",     &timeListCheat     );

}

void DeePID::CCNCTreeMaker::reconfigure(fhicl::ParameterSet const & p)
{
  // reconfigure without recompiling
  // read in the parameters from the .fcl file
  // allows for interactive changes of the parameter values

  // ==================================================================
  // Data labels
  // ==================================================================
  MCModuleLabel      = p.get< std::string >("MCModuleLabel");     //"generator";
  MCPartModLabel     = p.get< std::string >("MCPartModLabel");    //"largeant";
  HitLabel           = p.get< std::string >("HitLabel");          //"linecluster";
  pandoraLabel       = p.get< std::string >("pandoraLabel");      //"pandora";
  pandoraCheatLabel  = p.get< std::string >("pandoraCheatLabel"); //"cheatPandora";
  
}
//---------------------------------------------------------------------------------------------------------

void DeePID::CCNCTreeMaker::analyze(art::Event const & evt)
{
  // ==================================================================
  // Clear the vectors
  // ==================================================================

  wireList.clear();
  planeList.clear();
  TPCList.clear();
  integralList.clear();
  timeList.clear();

  wireListCheat.clear();
  planeListCheat.clear();
  TPCListCheat.clear();
  integralListCheat.clear();
  timeListCheat.clear();

  // ==================================================================
  // Get the data handles and list of particles
  // ==================================================================
  
  // MCTruth
  art::Handle< std::vector<simb::MCTruth> > mcHandle;
  std::vector< art::Ptr<simb::MCTruth> > mcList;
  if ( evt.getByLabel( MCModuleLabel, mcHandle ) ) art::fill_ptr_vector( mcList, mcHandle );
  
  // MCParticle
  art::Handle< std::vector<simb::MCParticle> > mcPartHandle;
  std::vector< art::Ptr<simb::MCParticle> > mcPartList;
  if ( evt.getByLabel( MCPartModLabel, mcPartHandle ) ) art::fill_ptr_vector( mcPartList, mcPartHandle );
  
  // Hits
  art::Handle< std::vector<recob::Hit> > hitHandle;
  std::vector< art::Ptr<recob::Hit> > hitList;
  if ( evt.getByLabel( HitLabel, hitHandle ) ) art::fill_ptr_vector( hitList, hitHandle );
  
  // Pandora
  art::Handle< std::vector<recob::PFParticle> > PFPartHandle;
  std::vector< art::Ptr<recob::PFParticle> > PFPartList;
  if ( evt.getByLabel( pandoraLabel, PFPartHandle ) ) art::fill_ptr_vector( PFPartList, PFPartHandle );
  
  // Clusters
  art::Handle< std::vector<recob::Cluster> > clusterHandle;
  std::vector< art::Ptr<recob::Cluster> > clusterList;
  if ( evt.getByLabel( pandoraLabel, clusterHandle ) ) art::fill_ptr_vector( clusterList, clusterHandle );
  
  // Pandora Cheated
  art::Handle< std::vector<recob::PFParticle> > PFPartHandleCheat;
  std::vector< art::Ptr<recob::PFParticle> > PFPartListCheat;
  if ( evt.getByLabel( pandoraCheatLabel, PFPartHandleCheat ) ) art::fill_ptr_vector( PFPartListCheat, PFPartHandleCheat );
  
  // Cheated Clusters
  art::Handle< std::vector<recob::Cluster> > clusterHandleCheat;
  std::vector< art::Ptr<recob::Cluster> > clusterListCheat;
  if ( evt.getByLabel( pandoraCheatLabel, clusterHandleCheat ) ) art::fill_ptr_vector( clusterListCheat, clusterHandleCheat );

  // Service handles
  art::ServiceHandle<cheat::BackTrackerService>       bt_serv;
  art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
  
  // ==================================================================
  // Associations
  // ==================================================================
  
  art::FindManyP< recob::Cluster > PFtoClust       ( PFPartHandle,      evt, pandoraLabel      );
  art::FindManyP< recob::Cluster > HitToClust      ( hitHandle,         evt, pandoraLabel      );
  
  art::FindManyP< recob::Cluster > PFtoClustCheat  ( PFPartHandleCheat, evt, pandoraCheatLabel );
  art::FindManyP< recob::Cluster > HitToClustCheat ( hitHandle,         evt, pandoraCheatLabel );
  
  // ==================================================================
  // Get the initial energy of the original MC Particle and T0
  // ==================================================================
  
  // Initial energy
  art::Ptr<simb::MCTruth> mcTruthProto = mcList[0];
  const simb::MCParticle& partp( mcTruthProto->GetParticle(0) );
  initialE = partp.E();
  int CCNC = mcTruthProto->GetNeutrino().CCNC();   // 0 for CC, 1 for NC
  mode = mcTruthProto->GetNeutrino().Mode(); // 0 for QE, 1 for Res, 2 for DIS, 3 for Coherent

  // PFParticles
  if ( PFPartHandle.isValid() || PFPartList.size() > 0 ) {
    PFinfo ( PFtoClust, HitToClust, hitList, PFPartList, 0);
  }
  else {
    std::cerr << "\n\n###################################\n"
              << "\nEvent " << evt.id().event() << "\n"
              << "PFPartHandle is not valid, or empty\n"
              << "PFPartList.size(): " << PFPartList.size() << "\n"
              << "###################################\n\n";
  }
  
  // Cheated PFParticles
  if ( PFPartHandleCheat.isValid() || PFPartListCheat.size() > 0 ) {
    PFinfo ( PFtoClustCheat, HitToClustCheat, hitList, PFPartListCheat, 1);
  }
  else {
    std::cerr << "\n\n########################################\n"
              << "\nEvent " << evt.id().event() << "\n"
              << "PFPartHandleCheat is not valid, or empty\n"
              << "PFPartListCheat.size(): " << PFPartListCheat.size() << "\n"
              << "########################################\n\n";
  }

  if ( CCNC == 0 ) CCtree->Fill();
  else NCtree->Fill();
}

// ------------------------------------------------------------------------------------------------------

std::vector< int > DeePID::CCNCTreeMaker::WireCorrector ( int wire, int plane, int tpc )
{
//  std::cout << "\nThis is a hit on wire " << wire << " on plane " << plane << " in tpc " << tpc << "\n";

  if ( plane == 2 ) {
  
    if ( tpc < 4 ) wire += 0;
    else if ( tpc < 8 ) wire += 480;
    else wire += 960;
    
  }
  
  else {
    
    if ( tpc % 2 == 0 ) {
      plane = ( plane == 1 ) ? 0 : 1;
    }
    
    if ( tpc < 4 ) wire += 0;
    else if ( tpc < 8 ) wire += 1148;
    else wire += 2296;
    
  }
  
//  std::cout << "It has been corrected to a hit on wire " << wire << " on plane " << plane << " in tpc " << tpc << "\n";
  
  std::vector< int > vals = { wire, plane };
  return vals;

}

void DeePID::CCNCTreeMaker::PFinfo ( art::FindManyP<recob::Cluster>             const & PFtoClust,
                                     art::FindManyP<recob::Cluster>             const & HitToClust,
                                     std::vector< art::Ptr<recob::Hit> >        const & hitList,
                                     std::vector< art::Ptr<recob::PFParticle> > const & PFPartList,
                                     bool   const & cheat )
{

  // Associate the hits to clusters
  std::map< art::Ptr<recob::Cluster>, std::vector< art::Ptr<recob::Hit> > > clustersToHits;

  for ( size_t hitIndex{0} ; hitIndex < hitList.size() ; hitIndex++ ) { // hitIndex Loop

    if ( !HitToClust.at(hitIndex).empty() ) {
      const art::Ptr<recob::Cluster> cluster = HitToClust.at(hitIndex)[0];
      const art::Ptr<recob::Hit> hit = hitList[hitIndex];
      clustersToHits[cluster].push_back(hit);
    }

  } // hitIndex Loop
  
  // Primaries
  for ( size_t PFIndex{0} ; PFIndex < PFPartList.size() ; PFIndex++ ) {
    auto particle = PFPartList[PFIndex];
    
    if ( particle->IsPrimary() ) {
    
      // Containers for all the info in the pftree
      std::vector< int > wires;
      std::vector< int > planes;
      std::vector< int > tpcs;
      std::vector< double > times;
      std::vector< double > integrals;
      
      auto daughters = particle->Daughters();
      
      if ( !daughters.empty() ) { // Daughters
        DaughterInfo( PFPartList, PFtoClust, clustersToHits, daughters, wires, planes, tpcs, times, integrals );
      } // Daughters

      std::vector< art::Ptr<recob::Cluster> > allClusters = PFtoClust.at( PFIndex );

      for ( size_t clustIndex{0} ; clustIndex < allClusters.size() ; clustIndex++ ) { // clusters

        const art::Ptr<recob::Cluster> cluster = allClusters.at(clustIndex);
        auto clustIt = clustersToHits.find(cluster);

        std::vector< art::Ptr<recob::Hit> > allHits = clustIt->second;

        for ( size_t h{0} ; h < allHits.size() ; h++ ) { // h loop
          art::Ptr<recob::Hit> hit = allHits[h];
          
          int wire = hit->WireID().Wire;
          int plane = hit->WireID().Plane;
          int tpc = hit->WireID().TPC;
          
          auto wireCorrection = WireCorrector ( wire, plane, tpc );
          
          wire = wireCorrection[0];
          plane = wireCorrection[1];
          
          wires.push_back    ( wire  );
          planes.push_back   ( plane );
          tpcs.push_back     ( tpc   );
          times.push_back    ( hit->PeakTime()     );
          integrals.push_back ( hit->Integral()     );
          
        } // h loop
    
      } // clusters
    
      if ( cheat ) {
        wireListCheat.push_back( wires );
        planeListCheat.push_back( planes );
        TPCListCheat.push_back( tpcs );
        integralListCheat.push_back( integrals );
        timeListCheat.push_back( times );
      }
      else {
        wireList.push_back( wires );
        planeList.push_back( planes );
        TPCList.push_back( tpcs );
        integralList.push_back( integrals );
        timeList.push_back( times );
      }
    
    } // Primary loop
  
  } // PFIndex loop
  
}

void DeePID::CCNCTreeMaker::DaughterInfo ( std::vector< art::Ptr<recob::PFParticle> > const & PFPartList,
                                             art::FindManyP<recob::Cluster>             const & PFtoClust,
                                             std::map< art::Ptr<recob::Cluster>, std::vector< art::Ptr<recob::Hit> > > const & clustersToHits,
                                             std::vector<size_t> const & daughters,
                                             std::vector< int >    & wires,
                                             std::vector< int >    & planes,
                                             std::vector< int >    & tpcs,
                                             std::vector< double > & times,
                                             std::vector< double > & integrals )
{
  
  for ( size_t i{0} ; i < daughters.size() ; i++ ) {
    auto particle = PFPartList[ daughters[i] ];
    auto newDaughters = particle->Daughters();
    
    if ( !newDaughters.empty() ) {
      DaughterInfo( PFPartList, PFtoClust, clustersToHits, newDaughters, wires, planes, tpcs, times, integrals );
    }
    
    std::vector< art::Ptr<recob::Cluster> > allClusters = PFtoClust.at( daughters[i] );

    for ( size_t clustIndex{0} ; clustIndex < allClusters.size() ; clustIndex++ ) {

      const art::Ptr<recob::Cluster> cluster = allClusters.at(clustIndex);
      auto clustIt = clustersToHits.find(cluster);

      std::vector< art::Ptr<recob::Hit> > allHits = clustIt->second;

      for ( size_t h{0} ; h < allHits.size() ; h++ ) {
        art::Ptr<recob::Hit> hit = allHits[h];
        
        int wire = hit->WireID().Wire;
        int plane = hit->WireID().Plane;
        int tpc = hit->WireID().TPC;
        
        auto wireCorrection = WireCorrector ( wire, plane, tpc );
        
        wire = wireCorrection[0];
        plane = wireCorrection[1];
        
        wires.push_back    ( wire  );
        planes.push_back   ( plane );
        tpcs.push_back     ( tpc   );
        times.push_back    ( hit->PeakTime()     );
        integrals.push_back ( hit->Integral()     );
      
      } // h loop
  
    } // This PF Particle Loop
    
  } // iterator loop
  
}
  


DEFINE_ART_MODULE(DeePID::CCNCTreeMaker)

