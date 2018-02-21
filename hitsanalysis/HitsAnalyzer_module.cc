////////////////////////////////////////////////////////////////////////
// Class:       HitsAnalyzer
// Plugin Type: analyzer (art v2_05_01)
// File:        HitsAnalyzer_module.cc
//
// Generated at Wed Feb 14 10:17:20 2018 by Nicolo Foppiani using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////


#include <cstdlib>
#include <iostream>


#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "art/Framework/Services/Optional/TFileService.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"

#include "TTree.h"

class HitsAnalyzer;


class HitsAnalyzer : public art::EDAnalyzer {
public:
  explicit HitsAnalyzer(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HitsAnalyzer(HitsAnalyzer const &) = delete;
  HitsAnalyzer(HitsAnalyzer &&) = delete;
  HitsAnalyzer & operator = (HitsAnalyzer const &) = delete;
  HitsAnalyzer & operator = (HitsAnalyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  std::string _pfp_producer = "pandoraNu";

  TTree * fVerticesTree;
  int fPdgCode;
  double fvx, fvy, fvz, fstartx, fstarty, fstartz;
  int fNhits, fNclusters;
  int fRun, fSubrun, fEvent;

  TTree * fHitsTree;
  int fPlane;
  int fWire;
  double fCharge;

  TTree * fSpacePointsTree;
  double fx, fy, fz, fChargeU, fChargeV, fChargeY;
};


HitsAnalyzer::HitsAnalyzer(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;
  
  fVerticesTree = tfs->make<TTree>("Vertices", "Vertices Tree");
  fVerticesTree->Branch("pdg_code", &fPdgCode, "pdg_code/i");
  fVerticesTree->Branch("vx", &fvx, "vx/d");
  fVerticesTree->Branch("vy", &fvy, "vy/d");
  fVerticesTree->Branch("vz", &fvz, "vz/d");
  fVerticesTree->Branch("start_x", &fstartx, "start_x/d");
  fVerticesTree->Branch("start_y", &fstarty, "start_y/d");
  fVerticesTree->Branch("start_z", &fstartz, "start_z/d");
  fVerticesTree->Branch("n_hits", &fNhits, "n_hits/i");
  fVerticesTree->Branch("n_clusters", &fNclusters, "n_clusters/i");
  fVerticesTree->Branch("event", &fEvent, "event/i");
  fVerticesTree->Branch("run", &fRun, "run/i");
  fVerticesTree->Branch("subrun", &fSubrun, "subrun/i");


  fHitsTree = tfs->make<TTree>("Hits", "Hits Tree");
  fHitsTree->Branch("pdg_code", &fPdgCode, "pdg_code/i");
  fHitsTree->Branch("plane", &fPlane, "plane/i");
  fHitsTree->Branch("wire", &fWire, "wire/i");
  fHitsTree->Branch("charge", &fCharge, "charge/d");
  fHitsTree->Branch("event", &fEvent, "event/i");
  fHitsTree->Branch("run", &fRun, "run/i");
  fHitsTree->Branch("subrun", &fSubrun, "subrun/i");

  fSpacePointsTree = tfs->make<TTree>("SpacePoints", "Space Points Tree");
  fSpacePointsTree->Branch("pdg_code", &fPdgCode, "pdg_code/i");
  fSpacePointsTree->Branch("x", &fx, "x/d");
  fSpacePointsTree->Branch("y", &fy, "y/d");
  fSpacePointsTree->Branch("z", &fz, "z/d");
  fSpacePointsTree->Branch("chargeU", &fChargeU, "chargeU/d");
  fSpacePointsTree->Branch("chargeV", &fChargeV, "chargeV/d");
  fSpacePointsTree->Branch("chargeY", &fChargeY, "chargeY/d");
  fSpacePointsTree->Branch("event", &fEvent, "event/i");
  fSpacePointsTree->Branch("run", &fRun, "run/i");
  fSpacePointsTree->Branch("subrun", &fSubrun, "subrun/i");
}

void HitsAnalyzer::analyze(art::Event const & evt)
{
  fRun = evt.run();
  fSubrun = evt.subRun();
  fEvent = evt.id().event();

  auto const &pfparticle_handle = 
    evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(_pfp_producer);
  auto const &spacepoint_handle = 
    evt.getValidHandle<std::vector<recob::SpacePoint>>(_pfp_producer);

  art::FindOneP<recob::Vertex> vertex_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindOneP<recob::Shower> shower_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_cluster(cluster_handle, evt, _pfp_producer);
  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, _pfp_producer);

  for (size_t i_pfp = 0; i_pfp < pfparticle_handle->size(); i_pfp++)
  {
    recob::PFParticle const &pfparticle = pfparticle_handle->at(i_pfp);
    fPdgCode = pfparticle.PdgCode();

    if (fPdgCode == 11)
    {
      try
      {
        art::Ptr<recob::Shower> pf_obj = shower_per_pfpart.at(i_pfp);
        fstartx = pf_obj->ShowerStart().X();
        fstarty = pf_obj->ShowerStart().Y();
        fstartz = pf_obj->ShowerStart().Z();
      }
      catch (...)
      {
        fstartx = 1000000.;
        fstarty = 1000000.;
        fstartz = 1000000.;
        std::cout << "No start point found for shower " << fPdgCode << std::endl;
      }
    }
    else if (fPdgCode == 13)
    {
      try
      {
        art::Ptr<recob::Track> pf_obj = track_per_pfpart.at(i_pfp);
        fstartx = pf_obj->Start().X();
        fstarty = pf_obj->Start().Y();
        fstartz = pf_obj->Start().Z();
      }
      catch (...)
      {
        fstartx = 1000000.;
        fstarty = 1000000.;
        fstartz = 1000000.;
        std::cout << "No start point found for track " << fPdgCode << std::endl;
      }
    }
    else
    {
      fstartx = 1000000.;
      fstarty = 1000000.;
      fstartz = 1000000.;
    }
      
    fNhits = 0;    
    fNclusters = 0;   
    // Hits
    std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(i_pfp);
    for (art::Ptr<recob::Cluster> &cluster : clusters)
    {
      fNclusters += 1;
      std::vector<art::Ptr<recob::Hit>> hits = hits_per_cluster.at(cluster.key());
      fNhits += hits.size();
      for (art::Ptr<recob::Hit> &hit : hits)
      {
        fPlane = hit->WireID().Plane;
        fWire = hit->WireID().Wire;
        fCharge = hit->Integral();
        fHitsTree->Fill();
      }
    }

    try
    {
      art::Ptr<recob::Vertex> vertex_obj = vertex_per_pfpart.at(i_pfp);
      double reco_neutrino_vertex[3];
      vertex_obj->XYZ(reco_neutrino_vertex);
      fvx = reco_neutrino_vertex[0];
      fvy = reco_neutrino_vertex[1];
      fvz = reco_neutrino_vertex[2];
    }
    catch (...)
    {
      fvx = 1000000.;
      fvy = 1000000.;
      fvz = 1000000.;
      std::cout << "No vertex found for " << fPdgCode << " with " << fNhits << std::endl;
    }

    fVerticesTree->Fill();

    // Space points
    std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(i_pfp);
    for (art::Ptr<recob::SpacePoint> &sps : spcpnts)
    {
      auto xyz = sps->XYZ();
      fx = xyz[0]; 
      fy = xyz[1];
      fz = xyz[2];
      
      std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(sps.key());
      fChargeU = 0;
      fChargeV = 0;
      fChargeY = 0;
      for (art::Ptr<recob::Hit> &hit : hits)
      {
        double hit_plane = hit->WireID().Plane;
        double hit_integral = hit->Integral();
        if (hit_plane == 0) fChargeU += hit_integral;
        else if (hit_plane == 1) fChargeV += hit_integral;
        else if (hit_plane == 2) fChargeY += hit_integral;
        else std::cout << "hit plane != 0, 1, 2, but " << hit_plane << std::endl; 
      }
      fSpacePointsTree->Fill();
    }
  }


}

DEFINE_ART_MODULE(HitsAnalyzer)
