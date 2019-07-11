////////////////////////////////////////////////////////////////////////
// Class:       NumuCCana
// Plugin Type: filter (art v2_05_01)
// File:        NumuCCana_module.cc
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
#include "art/Framework/Core/EDAnalyzer.h"
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

#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/CRT/CRTTrack.hh"
#include "ubobj/CRT/CRTTzero.hh"
#include "ubcrt/CRT/CRTAuxFunctions.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"

#include "lardataobj/RecoBase/OpFlash.h" 
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/AnalysisBase/T0.h"

#include "TTree.h"

#include <memory>

class NumuCCana;


class NumuCCana : public art::EDAnalyzer {
public:
  explicit NumuCCana(fhicl::ParameterSet const & p);

  NumuCCana(NumuCCana const &) = delete;
  NumuCCana(NumuCCana &&) = delete;
  NumuCCana & operator = (NumuCCana const &) = delete;
  NumuCCana & operator = (NumuCCana &&) = delete;

   void analyze(art::Event const & e) override;
  
  void beginJob() override;
  void endJob() override;
  
  void endSubRun(const art::SubRun &subrun);
  
  bool GetVertex(art::Event const &evt);
  bool GetMuon(const art::Ptr<recob::PFParticle> &pfp,
                       const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle,
                       const art::FindManyP<anab::ParticleID> &trackPIDAssn);
  void FillReconTruthMatching(art::Event const &evt);
  void FillTrueNu(art::Event const &evt);
  void FillTrueNuDaughters(art::Event const &evt);
  bool MatchDaughter(art::Event const &evt, const art::Ptr<recob::PFParticle> &pfp);
  
  void clearEvent();
  
  void initialize_tree();
  void initialize_pot();

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
  
  lar_pandora::ShowerVector pfshowers;
  lar_pandora::PFParticlesToShowers particlesToShowers;
  
  // Used for reco truth matching
  lar_pandora::PFParticlesToMCParticles matchedParticles;
  std::set<art::Ptr<simb::MCParticle>> matchedMCParticles;

  // producer datalabels
  std::string m_pfp_producer;
  std::string m_hitfinder_producer;
  //mc
  std::string m_geant_producer;
  std::string m_hit_mcp_producer;
  // MC neutrino daughter
  std::vector<int> fTrueNu_DaughterPDG;
  std::vector<float> fTrueNu_DaughterE;
  std::vector<bool> fTrueNu_DaughterMatched;
  
  std::string                         data_label_DAQHeader_;
  std::string                         data_label_flash_beam_;
  std::string                         data_label_crthit_;
  std::string                         data_label_crtT0asso_;
  std::string                         data_label_crttricorr_;

  int verbose_ = 0;
  int fHardDelay_;
  int fCRTT0off_;
  double beam_start_ = 0;
  double beam_end_ = 0;
  bool is_data_ = false;
  
  double crt_trig_corr_mean = 0;
  double crt_trig_corr_med = 0;
  
  double a_crthit_ts0 = 0;
  double a_crthit_ts1 = 0;
  int a_adc_length = 0;
  double a_crt_adc = 0;
  int a_t0_counter = 0;
  
  art::ServiceHandle<art::TFileService> tfs;
  
  TTree* my_event_;
  int has_neutrino_ = -1;
  double NuScore_ = -1;
  double FlashScore_ = -1;
  double FlashScoreTime_ = -1;
  int NuPDG_ = 0;
  int NumPfp_ = -1;
  double Vx_=-999, Vy_=-999, Vz_=-999;
  double Nu_Vx_=-999,Nu_Vy_=-999,Nu_Vz_=-999;
  double VtxDistance_=-99;
  double TrackScore_=-1;
  double TrackLength_=-99;
  double TrackPID_chiproton_=-99;
  double TrackPID_chimuon_=-99;
  double TrackPID_chipion_=-99;
  double TrackPID_chikaon_=-99;
  
  int NuTracks_ = 0;
  int NuShowers_ = 0;
  
  double TriTim_sec_ = 0;          //event trigger time sec
  double TriTim_nsec_ = 0;          //event trigger time ns
  
  int nr_crthit_ = -1; // # crt hits assigned to a tpc track
  double crthit_ts0_ = -99;
  double crthit_ts1_ = -99;
  int adc_length_ = -1;
  double crt_adc_ = -1;
  
  double TimFla_ = -99;
  double flash_PE_ = -99;
  double flash_y_ = -999;
  double flash_z_ = -999;
  
  double crtt0_time_ = -9999;
  int crtt0_trig_ = -1;
  double crtt0_DCA_ = -1;
  int crtt0_plane_ = -1;
  
  //mc
  // MC neutrino info
  uint NuMCnu; // number of MC neutrinos in event, only one gets saved!
  int MCNu_Interaction;
  int MCNu_CCNC;
  int MCNu_PDG;
  float MCNu_Energy;
  float MCNu_leptonPx, MCNu_leptonPy, MCNu_leptonPz;
  float MCNu_LeptonEnergy;
  float MCNu_Px, MCNu_Py, MCNu_Pz;
  float MCNu_leptonTheta;
  float MCNu_time; // time of the true neutrino interaction
  float MCNu_Vx, MCNu_Vy, MCNu_Vz;
  float MCNu_VxSce, MCNu_VySce, MCNu_VzSce;
  float MCNu_vertexDistance;
  
  // Matched MCParticle info
  bool MCNU_matched;
  bool MCCosmic_matched;
  int MCle_PDG;
  float MCle_Energy;
  float MCle_Vx, MCle_Vy, MCle_Vz;
  float MCle_length;
  float MCle_VxSce, MCle_VySce, MCle_VzSce;
  
  TTree* _sr_tree;
  int _sr_run = -9999;
  int _sr_subrun = -9999;
  double _sr_begintime = -9999;
  double _sr_endtime = -9999;
  double _sr_pot = -9999;
  //int event_counter = 0;
  
  // variables for cut
  
  //cout counters
  
  int counter_moun = 0;
  int counter_num_nu = 0;
  int counter_trackPIDAssn = 0;
  int counter_neutrino_metadata_vec = 0;
  int counter_getPID = 0;
};


NumuCCana::NumuCCana(fhicl::ParameterSet const & pset) :
    EDAnalyzer(pset)
{
    //producer datalabels
    m_pfp_producer = pset.get<std::string>("pfp_producer", "pandoraConsolidated");
    m_hitfinder_producer = pset.get<std::string>("hitfinder_producer", "gaushit");
    m_geant_producer = pset.get<std::string>("geant_producer", "largeant");
    m_hit_mcp_producer = pset.get<std::string>("hit_mcp_producer", "gaushitTruthMatch");
      
    data_label_DAQHeader_ = pset.get<std::string>("data_label_DAQHeader");
    data_label_flash_beam_ = pset.get<std::string>("data_label_flash_beam");
    data_label_crthit_ = pset.get<std::string>("data_label_crthit");
    data_label_crtT0asso_ = pset.get<std::string>("data_label_crtT0asso");
      
    data_label_crttricorr_ = pset.get<std::string>("data_label_crttricorr");
      
    fHardDelay_ = pset.get<int>("fHardDelay",40000);
    fCRTT0off_ = pset.get<int>("fCRTT0off",69000);
    beam_start_ = pset.get<double>("beam_start",3.2);
    beam_end_ = pset.get<double>("beam_end",5);
    is_data_ = pset.get<bool>("is_data",false);
    //cut variables
    //NuScore_cut = pset.get<double>("NuScore_cut", 0.0);

    verbose_ = pset.get<int>("verbose");
}

void NumuCCana::analyze(art::Event const & evt)
{
  clearEvent();
  
  //get the pfparticles and its metadata
  larpandora.CollectPFParticleMetadata(evt, m_pfp_producer, pfparticles, particlesToMetadata);
  larpandora.BuildPFParticleMap(pfparticles, particleMap);
  if (pfparticles.size() == 0){ //this should never happen
    std::cout << "[NumuCCana::filter] Event failed: No reconstructed PFParticles in event is " << pfparticles.size() << std::endl;
    //return false;
    if (!is_data_) // new for >= v08_00_00_16
      {
        FillReconTruthMatching(evt);
        FillTrueNu(evt);
        FillTrueNuDaughters(evt);
      }
  }
  else{
    //Get the neutrino candidate info
    larpandora.SelectNeutrinoPFParticles(pfparticles, pfneutrinos);
    if (pfneutrinos.size() != 1){ // some events (e.g out of TPC) have no neutrino reconstructed
      if(verbose_!=0) std::cout << "[NumuCCana::filter] Event failed: Number of reconstructed neutrinos in event is " << pfneutrinos.size() << std::endl;
      has_neutrino_= 0;
      counter_num_nu++;
    }
    else{ //check if neutrino candidate is muon neutrino
      has_neutrino_ = 1;
      if (!is_data_)
      {
        FillReconTruthMatching(evt);
        FillTrueNu(evt);
        FillTrueNuDaughters(evt);
      }
      if(!GetVertex(evt)){//Try to find its daughter particles for further studies
        if(verbose_!=0) std::cout << "[NumuCCana::filter] Event failed: GetVertex not passed" << std::endl;
        //return false;
      }
    }
  }
  if(verbose_!=0) std::cout << "[NumuCCana] Event Passes." << std::endl;
  // get event time for CRT time calculation
  art::Handle< raw::DAQHeaderTimeUBooNE > rawHandle_DAQHeader;
  evt.getByLabel(data_label_DAQHeader_, rawHandle_DAQHeader);
  raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
  art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();  
  TriTim_sec_ = evtTimeGPS.timeHigh();
  TriTim_nsec_ = evtTimeGPS.timeLow();
  
  art::Handle< std::vector<anab::T0> > rawHandle_CRTTriCorr;
  evt.getByLabel(data_label_crttricorr_, rawHandle_CRTTriCorr); //mergerextra
  std::vector<anab::T0> const& CRTTriCorrCollection(*rawHandle_CRTTriCorr);
  for(std::vector<int>::size_type i = 0; i != CRTTriCorrCollection.size(); i++) {
    crt_trig_corr_med = CRTTriCorrCollection.at(i).fTime;
    crt_trig_corr_mean = CRTTriCorrCollection.at(i).fTriggerConfidence;
  }
  
  // get CRT hits in the beam window
  art::Handle< std::vector<crt::CRTHit> > rawHandle_CRTHit;
  evt.getByLabel(data_label_crthit_, rawHandle_CRTHit); //mergerextra
  std::vector<crt::CRTHit> const& CRTHitCollection(*rawHandle_CRTHit);
  for(std::vector<int>::size_type i = 0; i != CRTHitCollection.size(); i++) {
    //double crthittime = (CRTHitCollection.at(i).ts0_ns + fCRTT0off_ - TriTim_nsec_)/1000.0;
    double crthittime = (CRTHitCollection.at(i).ts0_ns + crt_trig_corr_med - TriTim_nsec_)/1000.0;
    if( abs( crthittime - (beam_end_ + beam_start_)/2 ) < (beam_end_-beam_start_)/2 ){
        if(verbose_!=0){
          std::cout << "################################################" << std::endl;
          std::cout << "CRT Hit in beam window:" << std::endl;
          std::cout << "## Time GPS [us]:\t" << crthittime << std::endl;
          std::cout << "## Time Beam [us]:\t" << crthit_ts1_ << std::endl;
          std::cout << "################################################" << std::endl;
        }
        //fill tree variables
        //crthit_ts0_ = (double)(CRTHitCollection.at(i).ts0_ns + fCRTT0off_ - TriTim_nsec_)/1000.0;
        crthit_ts0_ = (double)(CRTHitCollection.at(i).ts0_ns + crt_trig_corr_med - TriTim_nsec_)/1000.0;
        crthit_ts1_ = ((double)CRTHitCollection.at(i).ts1_ns + fHardDelay_)/1000.0;
        crt_adc_ = CRTHitCollection.at(i).peshit;
        adc_length_ = CRTHitCollection.at(i).pesmap.begin()->second.size();
        if(verbose_!=0) std::cout << "## length of adc vector:\t" << adc_length_ << std::endl;
        nr_crthit_++;
        //has_crthit = 1;
    }
  }
  //get the beam flash
  art::Handle< std::vector<recob::OpFlash> > rawHandle_OpFlashBeam;
  evt.getByLabel(data_label_flash_beam_, rawHandle_OpFlashBeam);
  std::vector<recob::OpFlash> const& OpFlashCollectionBeam(*rawHandle_OpFlashBeam);
  if(verbose_!=0) std::cout << "There are: " << OpFlashCollectionBeam.size() << " in the beamflash collection" << std::endl;
  for(std::vector<int>::size_type i = 0; i != OpFlashCollectionBeam.size(); i++) {
    if( abs( OpFlashCollectionBeam.at(i).Time() - (beam_end_ + beam_start_)/2 ) < (beam_end_-beam_start_)/2 ){
      if(verbose_!=0){
      std::cout << "################################################" << std::endl;
      std::cout << "Flash in beam window:" << std::endl;
      std::cout << "## Time [us]:\t\t" << OpFlashCollectionBeam.at(i).Time() << std::endl;
      std::cout << "## Total PE []:\t" << OpFlashCollectionBeam.at(i).TotalPE() << std::endl;
      std::cout << "## Y [cm]:\t\t" << OpFlashCollectionBeam.at(i).YCenter() << std::endl;
      std::cout << "## Z [cm]:\t\t" << OpFlashCollectionBeam.at(i).ZCenter() << std::endl;
      std::cout << "################################################" << std::endl;
      }
      //fill tree variables
      TimFla_ = OpFlashCollectionBeam.at(i).Time();
      flash_PE_ = OpFlashCollectionBeam.at(i).TotalPE();
      flash_y_ = OpFlashCollectionBeam.at(i).YCenter();
      flash_z_ = OpFlashCollectionBeam.at(i).ZCenter();
    }
  }
  
  
  my_event_->Fill();
  //return true;

}
bool NumuCCana::GetVertex(art::Event const &evt)
{
  
  NumPfp_ = pfparticles.size();
  // Get vertex information
  lar_pandora::VertexVector vertexVector_dummy;
  larpandora.CollectVertices(evt, m_pfp_producer, vertexVector_dummy, particlesToVertices);

  // get information of downstream tracks
  larpandora.CollectTracks(evt, m_pfp_producer, pftracks, particlesToTracks);
  larpandora.CollectShowers(evt, m_pfp_producer, pfshowers, particlesToShowers);
  art::ValidHandle<std::vector<recob::Track>> trackHandle = evt.getValidHandle<std::vector<recob::Track> >(m_pfp_producer);
  const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle = evt.getValidHandle<std::vector<recob::MCSFitResult>>("pandoraMCSMu");
  const art::FindManyP<anab::ParticleID> trackPIDAssn(trackHandle, evt, "pandoracalipidSCE");
  if (!trackPIDAssn.isValid()){
    std::cout << "[NumuCCana::FillReconstructed] Event failed: PID is invalid" << std::endl;
    counter_trackPIDAssn++;
  }
  //check the neutrino information
  art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();
  NuPDG_ = pfnu->PdgCode(); // has to be 14
  
  art::Handle< std::vector<recob::PFParticle> > theParticles;
  evt.getByLabel(m_pfp_producer, theParticles);
  
  art::FindManyP<anab::T0> nuFlashScoreAsso(theParticles, evt, "flashmatch");
  //std::cout << "##################" << pfnu->Self() << std::endl;
  //std::cout << "##################" << nuFlashScoreAsso.size() << std::endl;
  //std::cout << "##################" << pfnu->IsPrimary() << std::endl;
  const std::vector<art::Ptr<anab::T0> > T0_v = nuFlashScoreAsso.at( pfnu.key() );
  //const std::vector<const anab::T0*>& T0_v = nuFlashScoreAsso.at( pfnu.key() );
  //std::cout << "######### Flash Score: " << T0_v.at(0)->Time() << " ########## " << std::endl;
  //std::cout << "######### Flash Score: " << T0_v.at(0)->TriggerType() << " ########## " << std::endl;
  //std::cout << "######### Flash Score: " << T0_v.at(0)->TriggerConfidence() << " ########## " << std::endl;
  //std::cout << "######### Flash Score: " << T0_v.at(0)->TriggerBits() << " ########## " << std::endl;
  //std::cout << "######### Flash Score: " << T0_v.size() << " ########## " << std::endl;
  if(T0_v.size()==1){
  FlashScoreTime_ = T0_v.at(0)->Time();
  FlashScore_ = T0_v.at(0)->TriggerConfidence();
  }
  else std::cout << "[NumuCCana::FillReconstructed] Flash score invalid" << std::endl;
  
  
  lar_pandora::MetadataVector neutrino_metadata_vec = particlesToMetadata.at(pfnu);
  lar_pandora::VertexVector neutrino_vertex_vec = particlesToVertices.at(pfnu);
  //check if there is a neutrino vertex reconstructed and only one
  if (neutrino_metadata_vec.size() != 1 || neutrino_vertex_vec.size() != 1){
    std::cout << "[NumuCCana::FillReconstructed] Event failed: Neutrino association failed" << std::endl;
    counter_neutrino_metadata_vec++;
    return false;
  }
  else{
    //check the topological score ov the event
    const larpandoraobj::PFParticleMetadata::PropertiesMap &neutrino_properties = neutrino_metadata_vec.front()->GetPropertiesMap();
    NuScore_ = neutrino_properties.at("NuScore"); // nuscore $$
    //FlashScore_ = neutrino_properties.at("FlashScore");
    if(verbose_ !=0) for (auto& t : neutrino_properties) std::cout << t.first << " " << t.second << std::endl;
    //for (auto& t : neutrino_properties) std::cout << t.first << " " << t.second << std::endl;
    //check the vertex position, if in fiducial volume
    const recob::Vertex::Point_t &neutrino_vtx = neutrino_vertex_vec.front()->position();
    Nu_Vx_ = neutrino_vtx.X();
    Nu_Vy_ = neutrino_vtx.Y();
    Nu_Vz_ = neutrino_vtx.Z();
  }
  //now get the muon candidate track for further investigtion (aka the longest track)
  pandoraInterfaceHelper.CollectDownstreamPFParticles(particleMap, pfnu, pfdaughters);
  double max_track_length = 0; //max muon length
  double muon_pfp_key = -1; // which pfp has longest track
  int pfp_counter = 0; // pfp counter
  // Implement here smart muon track id selection
  NuShowers_ = 0;
  NuTracks_ = 0;
  art::FindMany<anab::T0> trk_T0_assn_v(trackHandle, evt, data_label_crtT0asso_);
  art::FindMany<crt::CRTHit> trk_crt_assn_v(trackHandle, evt, data_label_crtT0asso_);
  for (auto const pfp : pfdaughters){  //loop over all daughter pfparticles
    if (!pfp->IsPrimary()){ 
      if (particlesToTracks.find(pfp) != particlesToTracks.end() ){ // get the track like pfp
        NuTracks_++;
        const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front(); //get the track
        if( this_track->Length() > max_track_length){ //take the longest track as muon candidate
          max_track_length = this_track->Length();
          muon_pfp_key = pfp_counter;
        }
        
        const std::vector<const anab::T0*>& T0crt = trk_T0_assn_v.at(this_track.key() );
        if(T0crt.size()!=0){
          if(verbose_!=0){
            std::cout << "################################################" << std::endl;
            std::cout << "Found T0 object from crthit - track association:" << std::endl;
            std::cout << "## Time [us]:\t\t" << T0crt.at(0)->Time() << std::endl;
            std::cout << "## Is from CRT []:\t" << T0crt.at(0)->TriggerType() << std::endl;
            std::cout << "## DCA [cm]:\t\t" << T0crt.at(0)->TriggerConfidence() << std::endl;
            std::cout << "## plane:\t\t" << T0crt.at(0)->TriggerBits() << std::endl;
            std::cout << "################################################" << std::endl;
          }
          //fill tree variables
          crtt0_time_ = T0crt.at(0)->Time();
          crtt0_trig_ = T0crt.at(0)->TriggerType();
          crtt0_DCA_ = T0crt.at(0)->TriggerConfidence();
          crtt0_plane_ = T0crt.at(0)->TriggerBits();
        }
        const std::vector<const crt::CRTHit*>& CRTHit_v = trk_crt_assn_v.at(this_track.key());
        for(std::vector<int>::size_type j = 0; j != CRTHit_v.size(); j++) {//CRThitloop
          //nr_crthit++;
          a_crthit_ts0 = CRTHit_v.at(j)->ts0_ns;
          a_crthit_ts1 = CRTHit_v.at(j)->ts1_ns;
          a_adc_length = CRTHit_v.at(j)->pesmap.begin()->second.size();
          a_crt_adc = CRTHit_v.at(j)->peshit;
          a_t0_counter++;
        } //end crthit loop
      }
      if (particlesToShowers.find(pfp) != particlesToShowers.end()){
        NuShowers_++;
      }
    }
    pfp_counter++;
  }
  if(muon_pfp_key!=-1){ //check if there was any track like pfp in the downstream 
    if (!GetMuon(pfdaughters.at(muon_pfp_key), MCSMu_handle, trackPIDAssn)){ // check if muon candidate fullfills all the requirements
      std::cout << "[NumuCCana::FillReconstructed] Event failed: Daughter investigation unsuccessful" << std::endl;
      return false;
    }
    else if (MatchDaughter(evt, pfdaughters.at(muon_pfp_key))){
      std::cout << "[NumuCCana::FillReconstructed] Matched MC daughter successfully" << std::endl;
    }
    return true; // found muon candidate, all cuts fullfilled, take this event
  }
  else{ // No track like object
    std::cout << "[NumuCCana::FillReconstructed] Event failed: No Daughter is muon/tracklike" << std::endl;
    counter_moun++;
    return false;
  }
}

bool NumuCCana::GetMuon(const art::Ptr<recob::PFParticle> &pfp,
                         const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle,
                         const art::FindManyP<anab::ParticleID> &trackPIDAssn){
  //check the track score value
  const larpandoraobj::PFParticleMetadata::PropertiesMap &pfp_properties = particlesToMetadata.at(pfp).front()->GetPropertiesMap();
  TrackScore_ = pfp_properties.at("TrackScore");
   // get start/vertex position of the track
  const recob::Vertex::Point_t &pfp_vtx = particlesToVertices.at(pfp).front()->position();
  Vx_ = pfp_vtx.X();//  vtx_distance < 5cm $$
  Vy_ = pfp_vtx.Y();
  Vz_ = pfp_vtx.Z();
  //check the distance to the nu vertex
  VtxDistance_ = pandoraInterfaceHelper.Distance3D(Vx_, Vy_, Vz_, Nu_Vx_, Nu_Vy_, Nu_Vz_);
  // check the track length, chi squares etc
  if (particlesToTracks.find(pfp) != particlesToTracks.end()){ // always true
    // get the recob track
    const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front();
    //check track length
    TrackLength_ = this_track->Length();
    // get the PID infos from association done in GetVertex
    std::map<std::string, float> pid_map;
    if(trackHelper.getPID(pid_map, this_track, trackPIDAssn)){
      TrackPID_chiproton_ = pid_map.at("chi2_proton");  // chi squqre cuts
      TrackPID_chimuon_ = pid_map.at("chi2_muon");
      TrackPID_chipion_ = pid_map.at("chi2_pion");
      TrackPID_chikaon_ = pid_map.at("chi2_kaon");
    }
    else{
      std::cout << "[NumuCCana::FillDaughters] Track has no PID attached to it" << std::endl;
      counter_getPID++;
      return false;
    }
  }
  // all requirements passed, to take this event with this muon candidate
  return true;
}
bool NumuCCana::MatchDaughter(art::Event const &evt, const art::Ptr<recob::PFParticle> &pfp)
{
  if (is_data_)
    return false;
  art::Ptr<simb::MCParticle> matched_mcp;
  double fGeneration = larpandora.GetGeneration(particleMap, pfp);
  if (fGeneration == 2)
  {
    if (matchedParticles.find(pfp) == matchedParticles.end())
      return false;
    matched_mcp = matchedParticles.at(pfp);
  }
  else if (fGeneration == 3)
  {
    // Generation 3 particle get matched to its parent.
    const auto iter(particleMap.find(pfp->Parent()));
    if (iter == particleMap.end())
      throw cet::exception("NuCC::MatchDaughter") << "Scrambled PFParticle IDs" << std::endl;
    const art::Ptr<recob::PFParticle> &pfp_parent = iter->second;

    if (matchedParticles.find(pfp_parent) == matchedParticles.end())
      return false;
    matched_mcp = matchedParticles.at(pfp_parent);
  }
  else
  {
    std::cout << "[NuCC::MatchDaughter] Generation 4 particle is not matched." << std::endl;
    return false;
  }

  if (!is_data_)
  {
    // Is this MC particle neutrino?
    const art::Ptr<simb::MCTruth> mctruth = pandoraInterfaceHelper.TrackIDToMCTruth(evt, m_geant_producer, matched_mcp->TrackId());
    if (mctruth->Origin() == simb::kBeamNeutrino)
    {
      MCNU_matched = true;
    }
    else
    {
      MCNU_matched = false;
      MCCosmic_matched = true;
    }

    MCle_PDG = matched_mcp->PdgCode();
    MCle_Energy = matched_mcp->E();
    MCle_Vx = matched_mcp->Vx();
    MCle_Vy = matched_mcp->Vy();
    MCle_Vz = matched_mcp->Vz();
    MCle_length = (matched_mcp->Position().Vect() - matched_mcp->EndPosition().Vect()).Mag();

    pandoraInterfaceHelper.SCE(MCle_Vx, MCle_Vy, MCle_Vz, matched_mcp->T(),
                               MCle_VxSce, MCle_VySce, MCle_VzSce);
    std::cout << "[NuCC::MatchDaughter] Daughter matched with PDG: " << MCle_PDG << ", neutrino origin: " << MCNU_matched << std::endl;
  }
  return true;
}
void NumuCCana::FillTrueNuDaughters(art::Event const &evt)
{
  if(verbose_!=0) std::cout << "[NumuCCana::FillTrueNuDaughters] entered function" << std::endl;
  lar_pandora::MCParticleVector mcparticles;
  larpandora.CollectMCParticles(evt, m_geant_producer, mcparticles);
  if(verbose_!=0) std::cout << "[NumuCCana::FillTrueNuDaughters] loop over MC particles of size: " << mcparticles.size() << std::endl;
  for (auto const &mcparticle : mcparticles)
  {
    if (!(mcparticle->Process() == "primary" &&
          mcparticle->T() != 0 &&
          mcparticle->StatusCode() == 1))
      continue;
    if(verbose_!=0) std::cout << "[NumuCCana::FillTrueNuDaughters] in the loop over MC particles" << std::endl;
    if(verbose_!=0) std::cout << "[NumuCCana::FillTrueNuDaughters] "<< m_geant_producer << " - " << mcparticle->TrackId() << std::endl;
    const art::Ptr<simb::MCTruth> mc_truth = pandoraInterfaceHelper.TrackIDToMCTruth(evt, m_geant_producer, mcparticle->TrackId());
    if(verbose_!=0) std::cout << "[NumuCCana::FillTrueNuDaughters] after trackIDToMCTrth function" << std::endl;
    if(verbose_!=0) std::cout << "[NumuCCana::FillTrueNuDaughters] " << mc_truth << std::endl;
    if(0)//if (mc_truth->Origin() == simb::kBeamNeutrino)
    {
      if(verbose_!=0) std::cout << "[NumuCCana::FillTrueNuDaughters] found daughter of MC beam neutrino" << std::endl;
      fTrueNu_DaughterE.push_back(mcparticle->E());
      fTrueNu_DaughterPDG.push_back(mcparticle->PdgCode());

      bool daughter_matched_neutrino_pfp = false;
      if (matchedMCParticles.find(mcparticle) != matchedMCParticles.end())
      {
        // Check if the corresponding pfparticle is also attached to the neutrino:
        for (auto const &[key, val] : matchedParticles)
        {
          if (val->TrackId() == mcparticle->TrackId())
          {
            if (larpandora.IsNeutrino(larpandora.GetParentPFParticle(particleMap, key)))
            {
              daughter_matched_neutrino_pfp = true;
              break;
            }
          }
        }
      }
      fTrueNu_DaughterMatched.push_back(daughter_matched_neutrino_pfp);
      if(verbose_!=0) std::cout << "[NuCC::FillTrueNuDaughters] << PDG: " << mcparticle->PdgCode() << ", E: " << mcparticle->E() << ", was matched? " << fTrueNu_DaughterMatched.back() << std::endl;
    }
    if(verbose_!=0) std::cout << "[NumuCCana::FillTrueNuDaughters] next MC particle" << std::endl;
  }
  if(verbose_!=0) std::cout << "[NumuCCana::FillTrueNuDaughters] end of function" << std::endl;
}
void NumuCCana::FillReconTruthMatching(art::Event const &evt)
{
  if(verbose_!=0) std::cout << "[NumuCCana::FillReconTruthMatching] entered function" << std::endl;
  pandoraInterfaceHelper.Configure(evt, m_pfp_producer, m_pfp_producer, m_hitfinder_producer, m_geant_producer, m_hit_mcp_producer);
   if(verbose_!=0) std::cout << "[NumuCCana::FillReconTruthMatching] Configured" << std::endl;
  pandoraInterfaceHelper.GetRecoToTrueMatches(matchedParticles);
  if(verbose_!=0) std::cout << "[NumuCCana::FillReconTruthMatching] " << std::endl;
  if(verbose_!=0) std::cout << "Number of PFPparticles in event: " << pfparticles.size() << std::endl;
  for (auto it = matchedParticles.begin(); it != matchedParticles.end(); ++it)
  {
    matchedMCParticles.insert(it->second);
  }
  if(verbose_!=0) std::cout << "[NumuCCana::FillReconTruthMatching] ";
  if(verbose_!=0) std::cout << "PFParticlesToMCParticles constructed: Number of PFPparticles matched: " << matchedParticles.size() << std::endl;
}
void NumuCCana::FillTrueNu(art::Event const &evt)
{
  if(verbose_!=0) std::cout << "[NumuCCana::FillTrueNu] entered function" << std::endl;
  if (!is_data_)
  {
    auto const &generator_handle = evt.getValidHandle<std::vector<simb::MCTruth>>("generator");
    auto const &generator(*generator_handle);
    int NuMCnu = generator.size();
    if(verbose_!=0) std::cout << "[NumuCCana::FillTrueNu] True neutrinos found: " << NuMCnu << std::endl;
    if (generator.size() > 0)
    {
      if (generator.front().Origin() != simb::kBeamNeutrino)
      {
        std::cout << "[NumuCCana::FillTrueNu] Origin of generator particle is not kBeamNeutrino." << std::endl;
        return;
      }
      const simb::MCNeutrino &mcnu = generator.front().GetNeutrino();

      MCNu_Interaction = mcnu.InteractionType();
      MCNu_CCNC = mcnu.CCNC();
      MCNu_PDG = mcnu.Nu().PdgCode();
      MCNu_Energy = mcnu.Nu().E();
      MCNu_Px = mcnu.Nu().Px();
      MCNu_Py = mcnu.Nu().Py();
      MCNu_Pz = mcnu.Nu().Pz();
      MCNu_LeptonEnergy = mcnu.Lepton().E();
      MCNu_leptonPx = mcnu.Lepton().Px();
      MCNu_leptonPy = mcnu.Lepton().Py();
      MCNu_leptonPz = mcnu.Lepton().Pz();
      MCNu_leptonTheta = mcnu.Theta();
      MCNu_time = mcnu.Nu().T();
      MCNu_Vx = mcnu.Nu().Vx();
      MCNu_Vy = mcnu.Nu().Vy();
      MCNu_Vz = mcnu.Nu().Vz();
      pandoraInterfaceHelper.SCE(MCNu_Vx, MCNu_Vy, MCNu_Vz, MCNu_time,
                                 MCNu_VxSce, MCNu_VySce, MCNu_VzSce);
      if(verbose_!=0) std::cout << ", CCNC: " << MCNu_CCNC << ", PDG: " << MCNu_PDG << ", E: " << MCNu_Energy << ", z-vertex: " << MCNu_Vz << std::endl;
    }
  }
}
void NumuCCana::clearEvent(){ // reset all pfp related vectors
    pfparticles.clear();
    pfneutrinos.clear();
    pfdaughters.clear();
    pftracks.clear();
    particleMap.clear();
    particlesToMetadata.clear();
    particlesToVertices.clear();
    particlesToTracks.clear();
  //mc
  matchedParticles.clear();
  matchedMCParticles.clear();
  
  has_neutrino_ = -1;
  NuScore_ = -1;
  FlashScore_ = -1;
  FlashScoreTime_ = -1;
  NuPDG_ = 0;
  NumPfp_ = -1;
  Vx_=-999, Vy_=-999, Vz_=-999;
  Nu_Vx_=-999,Nu_Vy_=-999,Nu_Vz_=-999;
  VtxDistance_=-99;
  TrackScore_=-1;
  TrackLength_=-99;
  TrackPID_chiproton_=-99;
  TrackPID_chimuon_=-99;
  TrackPID_chipion_=-99;
  TrackPID_chikaon_=-99;
  NuTracks_ = -1;
  NuShowers_ = -1;
  
  nr_crthit_ = -1; // # crt hits assigned to a tpc track
  crthit_ts0_ = -99;
  crthit_ts1_ = -99;
  adc_length_ = -1;
  crt_adc_ = -1;
  
  TimFla_ = -99;
  flash_PE_ = -99;
  flash_y_ = -999;
  flash_z_ = -999;
  
  crtt0_time_ = -9999;
  crtt0_trig_ = -1;
  crtt0_DCA_ = -1;
  crtt0_plane_ = -1;
  
  crt_trig_corr_mean = 0;
  crt_trig_corr_med = 0;
  
  a_crthit_ts0 = 0;
  a_crthit_ts1 = 0;
  a_adc_length = 0;
  a_crt_adc = 0;
  a_t0_counter = 0;
  
  if(!is_data_){
  NuMCnu = -1; 
  MCNu_Interaction = -1;
  MCNu_CCNC = -1;
  MCNu_PDG = -1;
  MCNu_Energy = -1;
  MCNu_leptonPx = -9, MCNu_leptonPy = -9, MCNu_leptonPz = -9;
  MCNu_LeptonEnergy = -9;
  MCNu_Px = -9, MCNu_Py = -9, MCNu_Pz = -9;
  MCNu_leptonTheta = -9;
  MCNu_time = -99; 
  MCNu_Vx = -999, MCNu_Vy = -999, MCNu_Vz = -999;
  MCNu_VxSce = -999, MCNu_VySce = -999, MCNu_VzSce = -999;
  MCNu_vertexDistance = -9;

  MCNU_matched = false;
  MCCosmic_matched = false;
  MCle_PDG = -1;
  MCle_Energy = -9;
  MCle_Vx = -999, MCle_Vy = -999, MCle_Vz = -999;
  MCle_length = -99;
  MCle_VxSce = -999, MCle_VySce = -999, MCle_VzSce = -999;
  }
}

void NumuCCana::initialize_tree()
{
  // Implementation of required member function here.
  std::cout << "Initialize variables and histograms for event tree" << std::endl;
  //tree stuff for hits: //////////////////////////////////////////////////////////////////////////////////
  my_event_ = tfs->make<TTree>("event","numuCC event tree");
  my_event_->Branch("NuScore",           &NuScore_,               "NuScore/D");
  my_event_->Branch("FlashScore",        &FlashScore_,            "FlashScore/D");
  my_event_->Branch("FlashScoreTime",        &FlashScoreTime_,            "FlashScoreTime/D");
  my_event_->Branch("NuPDG",             &NuPDG_,                 "NuPDG/I");
  my_event_->Branch("NumPfp",           &NumPfp_,                "NumPfp/I");
  
  my_event_->Branch("Nu_Vx",             &Nu_Vx_,                 "Nu_Vx/D");
  my_event_->Branch("Nu_Vy",             &Nu_Vy_,                 "Nu_Vy/D");
  my_event_->Branch("Nu_Vz",             &Nu_Vz_,                 "Nu_Vz/D");
  my_event_->Branch("Vx",                &Vx_,                    "Vx/D");
  my_event_->Branch("Vy",                &Vy_,                    "Vy/D");
  my_event_->Branch("Vz",                &Vz_,                    "Vz/D");
  
  my_event_->Branch("TrackScore",        &TrackScore_,            "TrackScore/D");
  my_event_->Branch("VtxDistance",       &VtxDistance_,           "VtxDistance_/D");
  my_event_->Branch("TrackLength",       &TrackLength_,           "TrackLength_/D");
  my_event_->Branch("TrackPID_chiproton",&TrackPID_chiproton_,    "TrackPID_chiproton_/D");
  my_event_->Branch("TrackPID_chipion",  &TrackPID_chipion_,      "TrackPID_chipion_/D");
  my_event_->Branch("TrackPID_chikaon",  &TrackPID_chikaon_,      "TrackPID_chikaon/D");
  my_event_->Branch("TrackPID_chimuon",  &TrackPID_chimuon_,      "TrackPID_chimuon_/D");
  
  my_event_->Branch("NuTracks",          &NuTracks_,              "NuTracks/I");
  my_event_->Branch("NuShowers",         &NuShowers_,             "NuShowers/I");
  
  my_event_->Branch("nr_crthit",         &nr_crthit_,              "nr_crthit_/I");
  my_event_->Branch("crthit_ts0",        &crthit_ts0_,            "crthit_ts0/D");
  my_event_->Branch("crthit_ts1",        &crthit_ts1_,            "crthit_ts1/D");
  my_event_->Branch("adc_length",        &adc_length_,            "adc_length/I");
  my_event_->Branch("crt_adc",           &crt_adc_,               "crt_adc/D");
  
  my_event_->Branch("TriTim_sec",        &TriTim_sec_,            "TriTim_sec/D");
  my_event_->Branch("TriTim_nsec",       &TriTim_nsec_,           "TriTim_nsec/D");
  
  my_event_->Branch("TimFla",            &TimFla_,                "TimFla/D");
  my_event_->Branch("flash_PE",          &flash_PE_,              "flash_PE/D");
  my_event_->Branch("flash_y",           &flash_y_,               "flash_y/D");
  my_event_->Branch("flash_z",           &flash_z_,               "flash_z/D");
  
  my_event_->Branch("crtt0_time",        &crtt0_time_,            "crtt0_time/D");
  my_event_->Branch("crtt0_trig",        &crtt0_trig_,            "crtt0_trig/I");
  my_event_->Branch("crtt0_DCA",         &crtt0_DCA_,             "crtt0_DCA/D");
  my_event_->Branch("crtt0_plane",       &crtt0_plane_,           "crtt0_plane/I");
  
  my_event_->Branch("crt_trig_corr_mean",&crt_trig_corr_mean,     "crt_trig_corr_mean/D");
  my_event_->Branch("crt_trig_corr_med", &crt_trig_corr_med,      "crt_trig_corr_med/D");
  my_event_->Branch("a_crthit_ts0",      &a_crthit_ts0,           "a_crthit_ts0/D");
  my_event_->Branch("a_crthit_ts1",      &a_crthit_ts1,           "a_crthit_ts1/D");
  my_event_->Branch("a_adc_length",      &a_adc_length,           "a_adc_length/I");
  my_event_->Branch("a_crt_adc",         &a_crt_adc,              "a_crt_adc/D");
  my_event_->Branch("a_t0_counter",      &a_t0_counter,           "a_t0_counter/I");
  
  if(!is_data_){
    my_event_->Branch("MCNu_Interaction",&MCNu_Interaction,"MCNu_Interaction/I");
    my_event_->Branch("MCNu_CCNC",           &MCNu_CCNC,           "MCNu_CCNC/I");
    my_event_->Branch("MCNu_PDG",            &MCNu_PDG,            "MCNu_PDG/I");
    my_event_->Branch("MCNu_Energy",         &MCNu_Energy,         "MCNu_Energy/F");
    my_event_->Branch("MCNu_leptonPx",       &MCNu_leptonPx,       "MCNu_leptonPx/F");
    my_event_->Branch("MCNu_LeptonEnergy",   &MCNu_LeptonEnergy,   "MCNu_LeptonEnergy/F");
    //my_event_->Branch("MCNu_Px",           &MCNu_Px,             "MCNu_Px/D");
    my_event_->Branch("MCNu_leptonTheta",    &MCNu_leptonTheta,    "MCNu_leptonTheta/F");
    my_event_->Branch("MCNu_time",           &MCNu_time,           "MCNu_time/F");
    my_event_->Branch("MCNu_Vx",             &MCNu_Vx,             "MCNu_Vx/F");
    my_event_->Branch("MCNu_Vy",             &MCNu_Vy,             "MCNu_Vy/F");
    my_event_->Branch("MCNu_Vz",             &MCNu_Vz,             "MCNu_Vz/F");
    my_event_->Branch("MCNu_VxSce",          &MCNu_VxSce,          "MCNu_VxSce/F");
    my_event_->Branch("MCNu_VySce",          &MCNu_VySce,          "MCNu_VySce/F");
    my_event_->Branch("MCNu_VzSce",          &MCNu_VzSce,          "MCNu_VzSce/F");
    my_event_->Branch("MCNu_vertexDistance",    &MCNu_vertexDistance,    "MCNu_vertexDistance/F");
    
    my_event_->Branch("MCle_PDG",               &MCle_PDG,               "MCle_PDG/I");
    my_event_->Branch("MCle_Energy",            &MCle_Energy,            "MCle_Energy/F");
    my_event_->Branch("MCle_Vx",                &MCle_Vx,                "MCle_Vx/F");
    my_event_->Branch("MCle_Vy",                &MCle_Vy,                "MCle_Vy/F");
    my_event_->Branch("MCle_Vz",                &MCle_Vz,                "MCle_Vz/F");
    my_event_->Branch("MCle_length",            &MCle_length,            "MCle_length/F");
    my_event_->Branch("MCle_VxSce",             &MCle_VxSce,             "MCle_VxSce/F");
    my_event_->Branch("MCle_VySce",             &MCle_VySce,             "MCle_VySce/F");
    my_event_->Branch("MCle_VzSce",             &MCle_VzSce,             "MCle_VzSce/F");

  }
  
  
  
}
void NumuCCana::initialize_pot()
{
  // Implementation of required member function here.
  std::cout << "Initialize variables and histograms for pot tree" << std::endl;
  //tree stuff for hits: //////////////////////////////////////////////////////////////////////////////////
  _sr_tree = tfs->make<TTree>("pottree","pottree");
  _sr_tree->Branch("run",                &_sr_run,                "run/I");
  _sr_tree->Branch("subrun",             &_sr_subrun,             "subrun/I");
  _sr_tree->Branch("begintime",          &_sr_begintime,          "begintime/D");
  _sr_tree->Branch("endtime",            &_sr_endtime,            "endtime/D");
  _sr_tree->Branch("pot",                &_sr_pot,                "pot/D");
 // _sr_tree->Branch("event_count",        &event_counter,          "event count/I");
  std::cout << "-------Using the following fcl parameters:-------" << std::endl;
  std::cout << "Pandora label:\t\t" << m_pfp_producer << std::endl;
  std::cout << "hitfinder label:\t\t" << m_hitfinder_producer << std::endl;
  std::cout << "geant producer label:\t\t" << m_geant_producer << std::endl;
  std::cout << "mcp producer label:\t\t" << m_hit_mcp_producer << std::endl;
  std::cout << "DAQ label:\t\t" << data_label_DAQHeader_ << std::endl;
  std::cout << "Beam flash label:\t" << data_label_flash_beam_ << std::endl;  std::cout << "mcp producer label:\t\t" << m_hit_mcp_producer << std::endl;
  std::cout << "CRT hit label:\t\t" << data_label_crthit_ << std::endl;
  std::cout << "CRT T0 asso label:\t" << data_label_crtT0asso_ << std::endl;
  
  std::cout << "fHardDelay:\t\t" << fHardDelay_ << std::endl;
  std::cout << "fCRTT0off:\t\t" << fCRTT0off_ << std::endl;
  std::cout << "Beam start:\t\t" << beam_start_ << std::endl;
  std::cout << "Beam end:\t\t" << beam_end_ << std::endl;

  std::cout << "is_data:\t\t\t" << is_data_ << std::endl;
  std::cout << "verbose:\t\t\t" << verbose_ << std::endl;
  std::cout << "------------end fcl parameters-------------------" << std::endl;
  
}
void NumuCCana::endSubRun(art::SubRun const &sr){
  
  if( verbose_ !=0) std::cout << "Write POT infos in to tree..." << std::endl;
  _sr_run       = sr.run();
  _sr_subrun    = sr.subRun();
  _sr_begintime = sr.beginTime().value();
  _sr_endtime   = sr.endTime().value();

  art::Handle<sumdata::POTSummary> potsum_h;
  if ( sr . getByLabel ( "generator" , potsum_h ) ) {
    _sr_pot = potsum_h->totpot;
    //std::cout << "Got POT info" << std::endl;
  }
  _sr_tree->Fill();
  
}

void NumuCCana::beginJob()
{
  // Implementation of optional member function here.
  initialize_tree();
  initialize_pot();
}
void NumuCCana::endJob(){
  // Implementation of optional member function here.
  std::cout.width(35); std::cout << "###### Filter summary  ######################" << std::left << std::endl;
  std::cout.width(35); std::cout << "Failed due to has one neutrino: " << std::left << counter_num_nu << std::endl;
  std::cout.width(35); std::cout << "Failed due to no muon candidate: " << std::left << counter_moun << std::endl;
  std::cout.width(35); std::cout << "Failed due to PID asso: " << std::left << counter_trackPIDAssn << std::endl;
  std::cout.width(35); std::cout << "Failed due to neutrino meta data: " << std::left << counter_neutrino_metadata_vec << std::endl;
  std::cout.width(35); std::cout << "Failed due to getPID: " << std::left << counter_getPID << std::endl;
  std::cout.width(35); std::cout << "###### End summary  #########################" << std::endl;
}

DEFINE_ART_MODULE(NumuCCana)
