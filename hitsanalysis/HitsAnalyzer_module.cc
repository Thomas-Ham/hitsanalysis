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

class HitsAnalyzer : public art::EDAnalyzer
{
public:
  explicit HitsAnalyzer(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HitsAnalyzer(HitsAnalyzer const &) = delete;
  HitsAnalyzer(HitsAnalyzer &&) = delete;
  HitsAnalyzer &operator=(HitsAnalyzer const &) = delete;
  HitsAnalyzer &operator=(HitsAnalyzer &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;

private:
  std::string _pfp_producer = "pandoraNu";

  TTree *fPFParticlesTree;
  int fPdgCode;
  double fvx, fvy, fvz;
  int fNhits, fNclusters;
  int fNhitsU, fNhitsV, fNhitsY;
  int fRun, fSubrun, fEvent;

  TTree *fClustersTree;
  double fClusterCharge, fClusterWidth, fClusterPosition;
  int fClusterNhits, fClusterPlane;

  TTree *fHitsTree;
  int fPlane;
  int fWire;
  double fCharge;

  TTree *fSpacePointsTree;
  double fx, fy, fz, fChargeU, fChargeV, fChargeY;
};

HitsAnalyzer::HitsAnalyzer(fhicl::ParameterSet const &p)
    : EDAnalyzer(p) // ,
                    // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;

  fPFParticlesTree = tfs->make<TTree>("PFParticles", "PFParticles Tree");
  fPFParticlesTree->Branch("pdg_code", &fPdgCode, "pdg_code/i");
  fPFParticlesTree->Branch("vx", &fvx, "vx/d");
  fPFParticlesTree->Branch("vy", &fvy, "vy/d");
  fPFParticlesTree->Branch("vz", &fvz, "vz/d");
  fPFParticlesTree->Branch("n_hits", &fNhits, "n_hits/i");
  fPFParticlesTree->Branch("n_hitsU", &fNhitsU, "n_hitsU/i");
  fPFParticlesTree->Branch("n_hitsV", &fNhitsV, "n_hitsV/i");
  fPFParticlesTree->Branch("n_hitsY", &fNhitsY, "n_hitsY/i");
  fPFParticlesTree->Branch("n_clusters", &fNclusters, "n_clusters/i");
  fPFParticlesTree->Branch("event", &fEvent, "event/i");
  fPFParticlesTree->Branch("run", &fRun, "run/i");
  fPFParticlesTree->Branch("subrun", &fSubrun, "subrun/i");

  fClustersTree = tfs->make<TTree>("Clusters", "Clusters Tree");
  fClustersTree->Branch("pdg_code", &fPdgCode, "pdg_code/i");
  fClustersTree->Branch("charge", &fClusterCharge, "charge/d");
  fClustersTree->Branch("width", &fClusterWidth, "width/d");
  fClustersTree->Branch("n_hits", &fClusterNhits, "n_hits/i");
  fClustersTree->Branch("plane", &fClusterPlane, "plane/i");
  fClustersTree->Branch("position", &fClusterPosition, "position/d");
  fClustersTree->Branch("event", &fEvent, "event/i");
  fClustersTree->Branch("run", &fRun, "run/i");
  fClustersTree->Branch("subrun", &fSubrun, "subrun/i");

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

void HitsAnalyzer::analyze(art::Event const &evt)
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

    fNhits = 0;
    fNhitsU = 0;
    fNhitsV = 0;
    fNhitsY = 0;
    fNclusters = 0;

    // Clusters and Hits
    std::vector<art::Ptr<recob::Cluster>> clusters = clusters_per_pfpart.at(i_pfp);
    fNclusters = clusters.size();
    for (art::Ptr<recob::Cluster> &cluster : clusters)
    {
      fClusterCharge = cluster->Integral();
      fClusterWidth = cluster->Width();
      fClusterNhits = cluster->NHits();
      fClusterPlane = cluster->View();
      fClusterPosition = (cluster->EndWire() + cluster->StartWire()) / 2.;
      fClustersTree->Fill();

      std::cout << "cluster for " << fPdgCode << " with " << fClusterCharge << "   " << cluster->Integral() << fClusterNhits << "   " << cluster->NHits() << std::endl;

      std::vector<art::Ptr<recob::Hit>> hits = hits_per_cluster.at(cluster.key());
      fNhits += hits.size();
      for (art::Ptr<recob::Hit> &hit : hits)
      {
        fPlane = hit->WireID().Plane;
        fWire = hit->WireID().Wire;
        fCharge = hit->Integral();
        fHitsTree->Fill();

        if (fPlane == 0)
        {
          fNhitsU += 1;
        }
        else if (fPlane == 1)
        {
          fNhitsV += 1;
        }
        else if (fPlane == 2)
        {
          fNhitsY += 1;
        }
      }
    }

    try
    {
      art::Ptr<recob::Vertex> vertex_obj = vertex_per_pfpart.at(i_pfp);
      double neutrino_vertex[3];
      vertex_obj->XYZ(neutrino_vertex);
      fvx = neutrino_vertex[0];
      fvy = neutrino_vertex[1];
      fvz = neutrino_vertex[2];
    }
    catch (...)
    {
      fvx = 1000000.;
      fvy = 1000000.;
      fvz = 1000000.;
      std::cout << "No vertex found for " << fPdgCode << " with " << fNhits << std::endl;
    }

    fPFParticlesTree->Fill();

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
        if (hit_plane == 0)
          fChargeU += hit_integral;
        else if (hit_plane == 1)
          fChargeV += hit_integral;
        else if (hit_plane == 2)
          fChargeY += hit_integral;
        else
          std::cout << "hit plane != 0, 1, 2, but " << hit_plane << std::endl;
      }
      fSpacePointsTree->Fill();
    }
  }
}

DEFINE_ART_MODULE(HitsAnalyzer)
