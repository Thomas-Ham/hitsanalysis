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
  void reconfigure(fhicl::ParameterSet const &p) override;

private:
  std::string m_pfp_producer;
  bool m_is_lite;

  TTree *fPFParticlesTree;
  int fPdgCode;
  double fvx, fvy, fvz;
  double fStartx, fStarty, fStartz;
  int fNhits, fNclusters, fNvertices, fNshowers_or_tracks;
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

  this->reconfigure(p);

  fPFParticlesTree = tfs->make<TTree>("PFParticles", "PFParticles Tree");
  fPFParticlesTree->Branch("pdg_code", &fPdgCode, "pdg_code/i");
  fPFParticlesTree->Branch("vx", &fvx, "vx/d");
  fPFParticlesTree->Branch("vy", &fvy, "vy/d");
  fPFParticlesTree->Branch("vz", &fvz, "vz/d");
  fPFParticlesTree->Branch("start_x", &fStartx, "start_x/d");
  fPFParticlesTree->Branch("start_y", &fStarty, "start_y/d");
  fPFParticlesTree->Branch("start_z", &fStartz, "start_z/d");
  fPFParticlesTree->Branch("n_hits", &fNhits, "n_hits/i");
  fPFParticlesTree->Branch("n_hitsU", &fNhitsU, "n_hitsU/i");
  fPFParticlesTree->Branch("n_hitsV", &fNhitsV, "n_hitsV/i");
  fPFParticlesTree->Branch("n_hitsY", &fNhitsY, "n_hitsY/i");
  fPFParticlesTree->Branch("n_clusters", &fNclusters, "n_clusters/i");
  fPFParticlesTree->Branch("n_vertices", &fNvertices, "n_vertices/i");
  fPFParticlesTree->Branch("n_showers_or_tracks", &fNshowers_or_tracks, "n_showers_or_tracks/i");
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

  if (m_is_lite == false)
  {
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
}

void HitsAnalyzer::reconfigure(fhicl::ParameterSet const &p)
{
  m_pfp_producer = p.get<std::string>("pfp_producer", "pandoraNu");
  m_is_lite = p.get<bool>("is_lite", false);
}

void HitsAnalyzer::analyze(art::Event const &evt)
{
  fRun = evt.run();
  fSubrun = evt.subRun();
  fEvent = evt.id().event();
  // std::cout << "Event: " << fEvent << ", Run: " << fRun << ", subRun: " << fSubrun << std::endl;

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(m_pfp_producer);
  auto const &cluster_handle =
      evt.getValidHandle<std::vector<recob::Cluster>>(m_pfp_producer);

  art::FindManyP<recob::Shower> showers_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
  art::FindManyP<recob::Track> tracks_per_pfpart(pfparticle_handle, evt, m_pfp_producer);

  art::FindManyP<recob::Vertex> vertices_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
  art::FindManyP<recob::Hit> hits_per_cluster(cluster_handle, evt, m_pfp_producer);
  auto const &spacepoint_handle =
        evt.getValidHandle<std::vector<recob::SpacePoint>>(m_pfp_producer);

  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, m_pfp_producer);

  for (size_t i_pfp = 0; i_pfp < pfparticle_handle->size(); i_pfp++)
  {
    recob::PFParticle const &pfparticle = pfparticle_handle->at(i_pfp);
    fPdgCode = pfparticle.PdgCode();

    std::vector<art::Ptr<recob::Vertex>> vertex_obj = vertices_per_pfpart.at(i_pfp);
    fNvertices = vertex_obj.size();
    std::cout << "Found vertices per pf part with lenght " << fNvertices << ", for pdg = " << fPdgCode << std::endl;
    if (fNvertices != 0)
    {
      double vertex_array[3];
      vertex_obj[0]->XYZ(vertex_array);
      fvx = vertex_array[0];
      fvy = vertex_array[1];
      fvz = vertex_array[2];
    }
    else
    {
      std::cout << "Zero vertices, for pdg = " << fPdgCode << std::endl;
      fvx = 1000000.;
      fvy = 1000000.;
      fvz = 1000000.;
    }

    if (fPdgCode == 11)
    {
      std::vector<art::Ptr<recob::Shower>> pf_objs = showers_per_pfpart.at(i_pfp);
      fNshowers_or_tracks = pf_objs.size();
      std::cout << "Found shower or tracks per pf part with lenght " << fNshowers_or_tracks << ", for pdg = " << fPdgCode << std::endl;
      if (fNshowers_or_tracks != 0)
      {
        fStartx = pf_objs[0]->ShowerStart().X();
        fStarty = pf_objs[0]->ShowerStart().Y();
        fStartz = pf_objs[0]->ShowerStart().Z();
      }
      else
      {
        fStartx = 1000000.;
        fStarty = 1000000.;
        fStartz = 1000000.;
        std::cout << "No start point found for shower " << fPdgCode << std::endl;
      }
    }
    else if (fPdgCode == 13)
    {
      std::vector<art::Ptr<recob::Track>> pf_objs = tracks_per_pfpart.at(i_pfp);
      fNshowers_or_tracks = pf_objs.size();
      std::cout << "Found shower or tracks per pf part with lenght " << fNshowers_or_tracks << ", for pdg = " << fPdgCode << std::endl;
      if (fNshowers_or_tracks != 0)
      {
        fStartx = pf_objs[0]->Start().X();
        fStarty = pf_objs[0]->Start().Y();
        fStartz = pf_objs[0]->Start().Z();
      }
      else
      {
        fStartx = 1000000.;
        fStarty = 1000000.;
        fStartz = 1000000.;
        std::cout << "No start point found for track " << fPdgCode << std::endl;
      }
    }
    else
    {
      fNshowers_or_tracks = 0;
      fStartx = 1000000.;
      fStarty = 1000000.;
      fStartz = 1000000.;
    }

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

      fNhits += fClusterNhits;

      std::vector<art::Ptr<recob::Hit>> hits = hits_per_cluster.at(cluster.key());
      for (art::Ptr<recob::Hit> &hit : hits)
      {
        fPlane = hit->WireID().Plane;
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

        if (m_is_lite == false)
        {
          fWire = hit->WireID().Wire;
          fCharge = hit->Integral();
          fHitsTree->Fill();
        }
      }
    }

    fPFParticlesTree->Fill();

    if (m_is_lite == false)
    {
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
}

DEFINE_ART_MODULE(HitsAnalyzer)
