////////////////////////////////////////////////////////////////////////
// Class:       NuCCFilter
// Plugin Type: filter (art v2_05_01)
// File:        NuCCFilter_module.cc
//
////////////////////////////////////////////////////////////////////////

/**
 * \brief CC-inclusive filter
 *
 * \author Thomas Metter
 *
 * \email thomas.mettler@lhep.unibe.ch
 *
 * \notes This module filters numu CC-inclusive events to fullfill the blinding requirements
 *        based on the mcc9 cc inclusive selection
 *
 */

// Base includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "larcoreobj/SummaryData/POTSummary.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "helpers/PandoraInterfaceHelper.h"
#include "helpers/EnergyHelper.h"
#include "helpers/TrackHelper.h"

#include "TTree.h"

#include <memory>

class NuCCFilter;


class NuCCFilter : public art::EDFilter {
public:
  explicit NuCCFilter(fhicl::ParameterSet const & p);

  NuCCFilter(NuCCFilter const &) = delete;
  NuCCFilter(NuCCFilter &&) = delete;
  NuCCFilter & operator = (NuCCFilter const &) = delete;
  NuCCFilter & operator = (NuCCFilter &&) = delete;

  bool filter(art::Event & e) override;
  void endJob() override;
  
  bool GetReconstructed(art::Event const &evt);
  bool GetMuon(const art::Ptr<recob::PFParticle> &pfp,
                       const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle,
                       const art::FindManyP<anab::ParticleID> &trackPIDAssn);
  void clearEvent();

private:
  
  PandoraInterfaceHelper pandoraInterfaceHelper;
  TrackHelper trackHelper;
  
  // LAr Pandora Helper fields
  lar_pandora::LArPandoraHelper larpandora;
  lar_pandora::PFParticleVector pfparticles;
  lar_pandora::PFParticleVector pfneutrinos;
  lar_pandora::PFParticleVector pfdaughters;
  lar_pandora::TrackVector pftracks;
  lar_pandora::PFParticleMap particleMap;
  lar_pandora::PFParticlesToMetadata particlesToMetadata;
  lar_pandora::PFParticlesToVertices particlesToVertices;
  lar_pandora::PFParticlesToTracks particlesToTracks;

  // producer datalabels
  std::string m_pfp_producer;
  std::string m_hitfinder_producer;

  int verbose_ = 0;
  
  // variables for cut
  int fNu_PDG;
  int fNumPfp;
  double fNu_Score;
  double fVx, fVy, fVz;
  double fNu_Vx,fNu_Vy,fNu_Vz;
  double fTrackScore;
  double fVtxDistance;
  double fTrackPID_chiproton;
  double fTrackPID_chimuon;
  double fTrackLength;
  // cut vaules
  double NuScore_cut;
  int NuPDG_cut;
  int startcontained_cut;
  double track_length_cut;
  double track_score_cut;
  double vrx_dist_cut;
  int longest_track_cut;
  double topologicalScore_cut;
  double track_chi2_muon_per_track_chi2_proton;
  double track_chi2_proton_cut;
  double track_chi2_muon_cut;
  
  double vertex_start_x_cut;
  double vertex_start_y_cut;
  double vertex_start_z_cut;
  
  double vertex_end_x_cut;
  double vertex_end_y_cut;
  double vertex_end_z_cut;
  
  //cout counters
  int counter_nuscore = 0;
  int counter_NuPDG = 0;
  int counter_moun = 0;
  int counter_tracklength = 0;
  int counter_track_score = 0;
  int counter_vrx_dist = 0;
  int counter_chi_muon = 0;
  int counter_chi_proton = 0;
  int counter_chi_ratio = 0;
  int counter_num_nu = 0;
  int counter_fiducial = 0;
  int counter_trackPIDAssn = 0;
  int counter_neutrino_metadata_vec = 0;
  int counter_getPID = 0;
};


NuCCFilter::NuCCFilter(fhicl::ParameterSet const & pset) :
    EDFilter(pset)
{
    //producer datalabels
    m_pfp_producer = pset.get<std::string>("pfp_producer", "pandoraConsolidated");
    m_hitfinder_producer = pset.get<std::string>("hitfinder_producer", "gaushit");
    //cut variables
    NuScore_cut = pset.get<double>("NuScore_cut", 0.0);
    NuPDG_cut = pset.get<int>("NuPDG_cut", 14);
    startcontained_cut = pset.get<int>("startcontained_cut", 0.06);
    track_length_cut = pset.get<double>("track_length_cut", 0.06);
    track_score_cut = pset.get<double>("track_score_cut", 0.06);
    vrx_dist_cut = pset.get<double>("vrx_dist_cut", 0.06);
    longest_track_cut = pset.get<int>("longest_track_cut", 1);
    topologicalScore_cut = pset.get<double>("topologicalScore_cut", 0.06);
    track_chi2_muon_per_track_chi2_proton = pset.get<double>("track_chi2_muon_per_track_chi2_proton", 0.06);
    track_chi2_proton_cut = pset.get<double>("track_chi2_proton_cut", 0.06);
    track_chi2_muon_cut = pset.get<double>("track_chi2_muon_cut", 0.06);
    vertex_start_x_cut = pset.get<double>("vertex_start_x_cut", 0.06);
    vertex_start_y_cut = pset.get<double>("vertex_start_y_cut", 0.06);
    vertex_start_z_cut = pset.get<double>("vertex_start_z_cut", 0.06);
    vertex_end_x_cut = pset.get<double>("vertex_end_x_cut", 0.06);
    vertex_end_y_cut = pset.get<double>("vertex_end_y_cut", 0.06);
    vertex_end_z_cut = pset.get<double>("vertex_end_z_cut", 0.06);

    verbose_ = pset.get<int>("verbose");
}

bool NuCCFilter::filter(art::Event & evt)
{
  clearEvent();
  //get the pfparticles and its metadata
  larpandora.CollectPFParticleMetadata(evt, m_pfp_producer, pfparticles, particlesToMetadata);
  larpandora.BuildPFParticleMap(pfparticles, particleMap);
  if (pfparticles.size() == 0){ //this should never happen
    std::cout << "[NuCCFilter::filter] Event failed: No reconstructed PFParticles in event is " << pfparticles.size() << std::endl;
    return false;
  }
  else{
    //Get the neutrino candidate info
    larpandora.SelectNeutrinoPFParticles(pfparticles, pfneutrinos);
    if (pfneutrinos.size() != 1){ // some events (e.g out of TPC) have no neutrino reconstructed
      std::cout << "[NuCCFilter::filter] Event failed: Number of reconstructed neutrinos in event is " << pfneutrinos.size() << std::endl;
      counter_num_nu++;
      return false;
    }
    else{ //check if neutrino candidate is muon neutrino
      art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();
      fNu_PDG = pfnu->PdgCode(); // has to be 14
      if( fNu_PDG != NuPDG_cut){
        std::cout << "[NuCCFilter::filter] Event failed: NuPDG is " << fNu_PDG << " [ ==" << NuPDG_cut << " ] " << std::endl;
        counter_NuPDG++;
        return false;
      }
      if(!GetReconstructed(evt)){//Try to find its daughter particles for further studies
        std::cout << "[NuCCFilter::filter] Event failed: GetReconstructed not passed" << std::endl;
        return false;
      }
    }
  }
  std::cout << "[NuCCFilter] Event Passes." << std::endl;
  return true;

}
bool NuCCFilter::GetReconstructed(art::Event const &evt)
{
  fNumPfp = pfparticles.size();
  // Get vertex information
  lar_pandora::VertexVector vertexVector_dummy;
  larpandora.CollectVertices(evt, m_pfp_producer, vertexVector_dummy, particlesToVertices);

  // get information of downstream tracks
  larpandora.CollectTracks(evt, m_pfp_producer, pftracks, particlesToTracks);
  art::ValidHandle<std::vector<recob::Track>> trackHandle = evt.getValidHandle<std::vector<recob::Track> >(m_pfp_producer);
  const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle = evt.getValidHandle<std::vector<recob::MCSFitResult>>("pandoraMCSMu");
  const art::FindManyP<anab::ParticleID> trackPIDAssn(trackHandle, evt, "pandoracalipidSCE");
  if (!trackPIDAssn.isValid()){
    std::cout << "[NuCCFilter::FillReconstructed] Event failed: PID is invalid" << std::endl;
    counter_trackPIDAssn++;
    return false;
  } 
  //check the neutrino information
  art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();
  lar_pandora::MetadataVector neutrino_metadata_vec = particlesToMetadata.at(pfnu);
  lar_pandora::VertexVector neutrino_vertex_vec = particlesToVertices.at(pfnu);
  //check if there is a neutrino vertex reconstructed and only one
  if (neutrino_metadata_vec.size() != 1 || neutrino_vertex_vec.size() != 1){
    std::cout << "[NuCCFilter::FillReconstructed] Event failed: Neutrino association failed" << std::endl;
    counter_neutrino_metadata_vec++;
    return false;
  }
  else{
    //check the topological score ov the event
    const larpandoraobj::PFParticleMetadata::PropertiesMap &neutrino_properties = neutrino_metadata_vec.front()->GetPropertiesMap();
    fNu_Score = neutrino_properties.at("NuScore"); // nuscore $$
    if( fNu_Score < NuScore_cut){
      std::cout << "[NuCCFilter::filter] Event failed: NuScore is " << fNu_Score << " [ <" << NuScore_cut << " ] " << std::endl;
      counter_nuscore++;
      return false;
    }
    //check the vertex position, if in fiducial volume
    const recob::Vertex::Point_t &neutrino_vtx = neutrino_vertex_vec.front()->position();
    fNu_Vx = neutrino_vtx.X();
    fNu_Vy = neutrino_vtx.Y();
    fNu_Vz = neutrino_vtx.Z();
    if( fNu_Vx<vertex_start_x_cut || fNu_Vx>vertex_end_x_cut || fNu_Vy<vertex_start_y_cut || fNu_Vy>vertex_end_y_cut || fNu_Vz<vertex_start_z_cut || fNu_Vz>vertex_end_z_cut){
      std::cout << "[NuCCFilter::filter] Event failed: Vertex out of fiducial " << fNu_Vx << " : " << fNu_Vy << " : " << fNu_Vz << std::endl;
      counter_fiducial++;
      return false;
    }
  }
  //now get the muon candidate track for further investigtion (aka the longest track)
  pandoraInterfaceHelper.CollectDownstreamPFParticles(particleMap, pfnu, pfdaughters);
  double max_track_length = 0; //max muon length
  double muon_pfp_key = -1; // which pfp has longest track
  int pfp_counter = 0; // pfp counter
  for (auto const pfp : pfdaughters){  //loop over all daughter pfparticles
    if (!pfp->IsPrimary()){ 
      if (particlesToTracks.find(pfp) != particlesToTracks.end() ){ // get the track like pfp
        const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front(); //get the track
        if( this_track->Length() > max_track_length){ //take the longest track as muon candidate
          max_track_length = this_track->Length();
          muon_pfp_key = pfp_counter;
        }
      }
    }
    pfp_counter++;
  }
  if(muon_pfp_key!=-1){ //check if there was any track like pfp in the downstream 
    if (!GetMuon(pfdaughters.at(muon_pfp_key), MCSMu_handle, trackPIDAssn)){ // check if muon candidate fullfills all the requirements
      std::cout << "[NuCCFilter::FillReconstructed] Event failed: Daughter investigation unsuccessful" << std::endl;
      return false;
    }
    return true; // found muon candidate, all cuts fullfilled, take this event
  }
  else{ // No track like object
    std::cout << "[NuCCFilter::FillReconstructed] Event failed: No Daughter is muon/tracklike" << std::endl;
    counter_moun++;
    return false;
  }
}

bool NuCCFilter::GetMuon(const art::Ptr<recob::PFParticle> &pfp,
                         const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle,
                         const art::FindManyP<anab::ParticleID> &trackPIDAssn){
  //check the track score value
  const larpandoraobj::PFParticleMetadata::PropertiesMap &pfp_properties = particlesToMetadata.at(pfp).front()->GetPropertiesMap();
  fTrackScore = pfp_properties.at("TrackScore");
  if( fTrackScore < track_score_cut){
    std::cout << "[NuCCFilter::FillDaughters] Event failed: track_score is " << fTrackScore << " [ <" << track_score_cut << " ] " << std::endl;
    counter_track_score++;
    return false;
  }
   // get start/vertex position of the track
  const recob::Vertex::Point_t &pfp_vtx = particlesToVertices.at(pfp).front()->position();
  fVx = pfp_vtx.X();//  vtx_distance < 5cm $$
  fVy = pfp_vtx.Y();
  fVz = pfp_vtx.Z();
  //check the distance to the nu vertex
  fVtxDistance = pandoraInterfaceHelper.Distance3D(fVx, fVy, fVz, fNu_Vx, fNu_Vy, fNu_Vz);
  if( fVtxDistance > vrx_dist_cut){
    std::cout << "[NuCCFilter::FillDaughters] Event failed: Track is to far from Vertex: " << fVtxDistance << " [ >" << vrx_dist_cut << " ] " << std::endl;
    counter_vrx_dist++;
    return false;
  }
  // check the track length, chi squares etc
  if (particlesToTracks.find(pfp) != particlesToTracks.end()){ // always true
    // get the recob track
    const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front();
    //check track length
    fTrackLength = this_track->Length();
    if( fTrackLength < track_length_cut){
      std::cout << "[NuCCFilter::FillDaughters] Event failed: Track to short: " << fTrackLength << " [ <" << track_length_cut << " ] " << std::endl;
      counter_tracklength++;
      return false;
    }
    // get the PID infos from association done in GetReconstructed
    std::map<std::string, float> pid_map;
    if(trackHelper.getPID(pid_map, this_track, trackPIDAssn)){
      fTrackPID_chiproton = pid_map.at("chi2_proton");  // chi squqre cuts
      fTrackPID_chimuon = pid_map.at("chi2_muon");
      //check proton chisquare
      if( fTrackPID_chiproton < track_chi2_proton_cut){
        std::cout << "[NuCCFilter::FillDaughters] Event failed: Proton chi square is to small: " << fTrackPID_chiproton << " [ <" << track_chi2_proton_cut << " ] " << std::endl;
        counter_chi_proton++;
        return false;
      }
      //check muon chisquare
      if( fTrackPID_chimuon > track_chi2_muon_cut){
        std::cout << "[NuCCFilter::FillDaughters] Event failed: Muon chi square is to high: " << fTrackPID_chimuon << " [ >" << track_chi2_muon_cut << " ] " << std::endl;
        counter_chi_muon++;
        return false;
      }
      double ratio = fTrackPID_chimuon/fTrackPID_chiproton;
      //check the ratio of them
      if( ratio > track_chi2_muon_per_track_chi2_proton){
        std::cout << "[NuCCFilter::FillDaughters] Event failed: CHi square ration is to high: " << ratio << " [ >" << track_chi2_muon_per_track_chi2_proton << " ] " << std::endl;
        counter_chi_ratio++;
        return false;
      }
    }
    else{
      std::cout << "[NuCCFilter::FillDaughters] Track has no PID attached to it" << std::endl;
      counter_getPID++;
      return false;
    }
  }
  // all requirements passed, to take this event with this muon candidate
  return true;
}

void NuCCFilter::clearEvent(){ // reset all pfp related vectors
    pfparticles.clear();
    pfneutrinos.clear();
    pfdaughters.clear();
    pftracks.clear();
    particleMap.clear();
    particlesToMetadata.clear();
    particlesToVertices.clear();
    particlesToTracks.clear();
}

void NuCCFilter::endJob(){
  // Implementation of optional member function here.
  std::cout << "###### Filter summary  ######################" << std::endl;
  std::cout.width(35); std::cout << "Failed due to has one neutrino: " << std::left << counter_num_nu << std::endl;
  std::cout.width(35); std::cout << "Failed due to Nu PDG: " << std::left << counter_NuPDG << std::endl;
  std::cout.width(35); std::cout << "Failed due to Nuscore/toposcore: " << std::left << counter_nuscore << std::endl;
  std::cout.width(35); std::cout << "Failed due to fiducial: " << std::left << counter_fiducial << std::endl;
  std::cout.width(35); std::cout << "Failed due to no muon candidate: " << std::left << counter_moun << std::endl;
  std::cout.width(35); std::cout << "Failed due to tracklength: " << std::left << counter_tracklength << std::endl;
  std::cout.width(35); std::cout << "Failed due to track_score: " << std::left << counter_track_score << std::endl;
  std::cout.width(35); std::cout << "Failed due to vrx_dist: " << std::left << counter_vrx_dist << std::endl;
  std::cout.width(35); std::cout << "Failed due to chi_muon: " << std::left << counter_chi_muon << std::endl;
  std::cout.width(35); std::cout << "Failed due to chi_proton: " << std::left << counter_chi_proton << std::endl;
  std::cout.width(35); std::cout << "Failed due to chi_ratio: " << std::left << counter_chi_ratio << std::endl;
  std::cout.width(35); std::cout << "Failed due to PID asso: " << std::left << counter_trackPIDAssn << std::endl;
  std::cout.width(35); std::cout << "Failed due to neutrino meta data: " << std::left << counter_neutrino_metadata_vec << std::endl;
  std::cout.width(35); std::cout << "Failed due to getPID: " << std::left << counter_getPID << std::endl;
  std::cout.width(35); std::cout << "###### End summary  #########################" << std::endl;
}

DEFINE_ART_MODULE(NuCCFilter)
