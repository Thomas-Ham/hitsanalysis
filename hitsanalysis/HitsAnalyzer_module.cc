////////////////////////////////////////////////////////////////////////
// Class:       HitsAnalyzer
// Plugin Type: analyzer (art v2_05_01)
// File:        HitsAnalyzer_module.cc
//
// Generated at Wed Feb 14 10:17:20 2018 by Nicolo Foppiani using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

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
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"

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
  double fvx, fvy, fvz;

  TTree * fHitsTree;
  int fPlane;
  int fWire;
  int fCharge;

  TTree * fSpacePointsTree;
  double fx, fy, fz, fChargeU, fChargeV, fChargeY;
  int fEvent;
};


HitsAnalyzer::HitsAnalyzer(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;
  
  fVerticesTree = tfs->make<TTree>("Vertices Tree", "Vertices Tree");
  fVerticesTree->Branch("pdg_code", &fPdgCode, "pdg_code/i");
  fVerticesTree->Branch("vx", &fvx, "vx/d");
  fVerticesTree->Branch("vy", &fvy, "vy/d");
  fVerticesTree->Branch("vz", &fvz, "vz/d");


  fHitsTree = tfs->make<TTree>("HitsTree", "Hits Tree");
  fHitsTree->Branch("pdg_code", &fPdgCode, "pdg_code/i");
  fHitsTree->Branch("plane", &fPlane, "plane/d");
  fHitsTree->Branch("wire", &fWire, "wire/d");
  fHitsTree->Branch("charge", &fCharge, "charge/d");

  fSpacePointsTree = tfs->make<TTree>("SpacePointsTree", "Space Points Tree");
  fSpacePointsTree->Branch("pdg_code", &fPdgCode, "pdg_code/i");
  fSpacePointsTree->Branch("x", &fx, "x/d");
  fSpacePointsTree->Branch("y", &fy, "y/d");
  fSpacePointsTree->Branch("z", &fz, "z/d");
  fSpacePointsTree->Branch("chargeU", &fChargeU, "chargeU/d");
  fSpacePointsTree->Branch("chargeV", &fChargeV, "chargeV/d");
  fSpacePointsTree->Branch("chargeY", &fChargeY, "chargeY/d");
}

void HitsAnalyzer::analyze(art::Event const & e)
{
  // Hits
  // art::Handle< std::vector<recob::Hit> > hitListHandle;
  // std::vector<art::Ptr<recob::Hit> > hitlist;
  // if (evt.getByLabel(fHitsModuleLabel,hitListHandle))
  //   art::fill_ptr_vector(hitlist, hitListHandle);


  auto const &pfparticle_handle = 
    evt.getValidHandle<std::vector<recob::PFParticle>>(_pfp_producer);
  
  art::FindOneP<recob::Vertex> vertex_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::Hit> hits_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  art::FindManyP<recob::SpacePoint> 
    spcpnts_per_pfpart(pfparticle_handle, evt, _pfp_producer);
  auto const &spacepoint_handle = 
    evt.getValidHandle<std::vector<recob::SpacePoint>>(_pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, _pfp_producer);

  for (size_t i_pfp = 0; i_pfp < pfparticle_handle->size(); i_pfp++)
  {
    recob::PFParticle const &pfparticle = pfparticle_handle->at(i_pfp);
    fPdgCode = pfparticle->PdgCode();

    recob::Vertex const &vertex_obj = vertex_per_pfpart.at(ipf_candidate);
    double reco_neutrino_vertex[3];
    vertex_obj->XYZ(reco_neutrino_vertex);
    vx = reco_neutrino_vertex[0];
    vy = reco_neutrino_vertex[1];
    vz = reco_neutrino_vertex[2];
    fVerticesTree->Fill();
  
    // Hits
    std::vector<art::Ptr<recob::Hit>> hits = hits_per_pfpart.at(i_pfp);
    for (art::Ptr<recob::Hit> &hit : hits)
    {
      fPlane = hit->WireID().Plane;
      fWire = hit->WireID().Wire;
      fCharge = hit->Integral();
      fHitsTree->Fill();
    }

    // Space points
    std::vector<art::Ptr<recob::SpacePoint>> spcpnts = spcpnts_per_pfpart.at(i_pfp);
    for (art::Ptr<recob::SpacePoint> &sps : spcpnts)
    {
      double xyz[3] = sps->XYZ();
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
        double hit_wire = hit->WireID().Wire;
        double hit_integral = hit->Integral();
        if (hit_plane == 0) fChargeU += integral;
        else if (hit_plane == 1) fChargeV += integral;
        else if (hit_plane == 2) fChargeY += integral;
        else std::cout << "hit plane != 0, 1, 2, but " << hit_plane << endl; 
      }
      fSpacePointsTree->Fill();
    }
  }


}

DEFINE_ART_MODULE(HitsAnalyzer)
