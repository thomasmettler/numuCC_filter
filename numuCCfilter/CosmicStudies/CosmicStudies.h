////////////////////////////////////////////////////////////////////////
// Class:       CosmicStudies
// Plugin Type: analyzer (art v2_11_03)
// File:        CosmicStudies_module.cc
//
// Generated at Wed Oct 31 16:11:25 2018 by Wouter Van de pontseele using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#ifndef COSMICMODULE_H
#define COSMICMODULE_H

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larcoreobj/SummaryData/POTSummary.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "helpers/GeometryHelper.h"
#include "helpers/MCParticleHelper.h"
#include "helpers/PandoraInterfaceHelper.h"

#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

namespace constants
{
const float MCP_E_CUT = 0.1;      // Only save MC particles above this energy
const float PFP_LENGTH_CUT = 5.0; // Only save reconstructed tracks above this length
const float MUON_M_MEV = 105.658; // Mass of muons, in MeV
const float CRT_E_CUT = 0.005;     // Minimum particle energy to create a crt hit
} // namespace constants

class CosmicStudies : public art::EDAnalyzer
{
  public:
    explicit CosmicStudies(fhicl::ParameterSet const &p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    CosmicStudies(CosmicStudies const &) = delete;
    CosmicStudies(CosmicStudies &&) = delete;
    CosmicStudies &operator=(CosmicStudies const &) = delete;
    CosmicStudies &operator=(CosmicStudies &&) = delete;

    // Implemented in header
    void clear();
    void clear_MCParticle();
    void clear_PFParticle();
    void clear_Cluster();
    void clear_Flashes();
    // Larpandora way of backtracking
    art::Ptr<simb::MCTruth> TrackIDToMCTruth(art::Event const &e, std::string m_geant_producer, int geant_track_id);

    // Implemented in module.
    void analyze(art::Event const &evt) override;
    void endSubRun(const art::SubRun &sr);
    void reconfigure(fhicl::ParameterSet const &p);

    void fill_flash(art::ValidHandle<std::vector<recob::OpFlash>> const &simple_cosmic_handle, uint number, TTree *tree);
    void fill_MC(art::Event const &evt);
    void fill_TPCreco(art::Event const &evt);

  private:
    // FCL parameters
    std::string m_pfp_producer;
    std::string m_spacepoint_producer;
    std::string m_hitfinder_producer;
    std::string m_geant_producer;
    std::string m_hit_mcp_producer;
    std::string m_cosmic_opflash_producer;
    std::string m_cosmic_simpleflash_producer;
    std::string m_cosmic_ophit_producer;
    std::string m_beam_opflash_producer;
    std::string m_beam_simpleflash_producer;
    std::string m_beam_ophit_producer;

    bool m_is_mcc9;
    bool m_is_data;
    bool m_is_true_nu;
    bool m_verb;
    bool m_slimmed;

    // Other private fields
    PandoraInterfaceHelper pandoraHelper;
    lar_pandora::PFParticlesToMCParticles matchedParticles;
    std::set<art::Ptr<simb::MCParticle>> matchedMCParticles;
    MCParticleHelper mcpHelper;

    // Fields for in the tree!
    TTree *fPOTTree;
    uint fRun_sr, fSubrun_sr;
    float fPot;
    float fDatasetPrescaleFactor; // Prescaling factor for data.
    uint fNevents = 0;

    TTree *fEventTree;
    uint fRun, fSubrun, fEvent;
    uint fNumSimpleBeamFlashes;
    uint fNumOpBeamFlashes;
    uint fNumSimpleCosmicFlashes;
    uint fNumOpCosmicFlashes;
    uint fNumMcp;
    uint fNumMcp_saved; // criteria to keep only relevant particles.
    uint fNumPfp;
    uint fNumPfp_saved;
    // neutrino fields if there is a neutrino:
    uint fNum_nu;
    std::vector<float> fNu_vtx_x;
    std::vector<float> fNu_vtx_y;
    std::vector<float> fNu_vtx_z;
    std::vector<float> fNu_E;
    std::vector<float> fNu_time;
    std::vector<int> fNu_pdg_code;
    std::vector<bool> fNu_ccnc;

    TTree *fMCParticlesTree;
    float fMc_E;
    int fMc_PdgCode;
    float fMc_StartX;
    float fMc_StartY;
    float fMc_StartZ;
    float fMc_EndX;
    float fMc_EndY;
    float fMc_EndZ;
    float fMc_StartMomentumX;
    float fMc_StartMomentumY;
    float fMc_StartMomentumZ;
    int fMc_StatusCode;
    float fMc_Time;
    bool fMc_kBeamNeutrino;
    // using the MCPhelper
    uint fMc_Process; // std::string
    uint fMc_EndProcess;
    bool fMc_Matched; // is this mc particle matched to a PFParticle?
    bool fMc_StartInside;
    bool fMc_EndInside;
    bool fMc_PartInside; // This means that the track is crossing, starts/ends inside, is completely inside.
    bool fCRT_crossed;
    float fMc_StartX_tpc; // if the track starts outside but crosses, this will be stored here.
    float fMc_StartY_tpc;
    float fMc_StartZ_tpc;
    float fMc_EndX_tpc;
    float fMc_EndY_tpc;
    float fMc_EndZ_tpc;
    float fMc_Length;
    float fMc_LengthTPC;

    TTree *fSimpleCosmicFlashesTree;
    TTree *fOpCosmicFlashesTree;
    TTree *fSimpleBeamFlashesTree;
    TTree *fOpBeamFlashesTree;
    float fFlash_Time;
    float fFlash_DiffTime;
    uint fFlash_TotalPE;
    uint fFlash_num10percentPMT; // The number of PMT that is responsible for more than 10% of the total flash.
    float fFlash_Z;
    float fFlash_sigmaZ;
    float fFlash_Y;
    float fFlash_sigmaY;
    float fFlash_Width;
    float fFlash_AbsTime;

    TTree *fPFParticlesTree;
    int fPdgCode;
    uint fNumDaughters;
    bool fIsPrimary;
    float fVx, fVy, fVz;
    uint fNhits, fNclusters;
    uint fNhitsU, fNhitsV, fNhitsY;
    // track
    bool fTrack_Valid;
    bool fTrack_HasMomentum;
    float fTrack_StartX;
    float fTrack_StartY;
    float fTrack_StartZ;
    float fTrack_EndX;
    float fTrack_EndY;
    float fTrack_EndZ;
    float fTrack_Length;
    float fTrack_StartMomentumX;
    float fTrack_StartMomentumY;
    float fTrack_StartMomentumZ;
    float fTrack_EndMomentumX;
    float fTrack_EndMomentumY;
    float fTrack_EndMomentumZ;
    float fTrack_Theta;
    float fTrack_Phi;
    float fTrack_ZenithAngle;
    float fTrack_AzimuthAngle;
    // MC fields
    float fTrack_MCS_mom;
    float fTrack_MCS_err;
    float fTrack_MCS_ll;
    float fTrack_MCS_energy;
    // reco-truth fields
    int fTrack_matched_PdgCode;
    float fTrack_matchedE;
    bool fTrack_matched_kBeamNeutrino;
    uint fTrack_matched_Process;    // std::string
    uint fTrack_matched_EndProcess; // std::string
    float fTrack_matched_Time;
    bool fTrack_matched_StartInside;
    bool fTrack_matched_EndInside;
    bool fTrack_matched_PartInside;  // This means that the track is crossing, starts/ends inside, is completely inside.
    float fTrack_matched_StartX_sce; // spacecharge corrected version of the tpc edge or the inside start end point.
    float fTrack_matched_StartY_sce;
    float fTrack_matched_StartZ_sce;
    float fTrack_matched_EndX_sce;
    float fTrack_matched_EndY_sce;
    float fTrack_matched_EndZ_sce;
    float fTrack_matched_StartX; // spacecharge corrected version of the tpc edge or the inside start end point.
    float fTrack_matched_StartY;
    float fTrack_matched_StartZ;
    float fTrack_matched_EndX;
    float fTrack_matched_EndY;
    float fTrack_matched_EndZ;
    float fTrack_matched_LengthTPC;
    float fTrack_matched_Length_sce;
    float fTrack_matched_StartMomentumX;
    float fTrack_matched_StartMomentumY;
    float fTrack_matched_StartMomentumZ;

    TTree *fClustersTree;
    float fClusterCharge, fClusterWidth, fClusterPosition;
    uint fClusterNhits, fClusterPlane;

    TTree *fCRTcrossTree;
    uint fNumCross;
    float fCrossX;
    float fCrossY;
    float fCrossZ;
    float fCrossE;
    float fCrossT;
};

CosmicStudies::CosmicStudies(fhicl::ParameterSet const &p)
    : EDAnalyzer(p) // ,
                    // More initializers here.
{
    art::ServiceHandle<art::TFileService> tfs;

    this->reconfigure(p);
    mcpHelper.Configure(m_is_mcc9);

    //// Check if things are set up properly:
    std::cout << std::endl;
    std::cout << "[CosmicStudies constructor] Checking set-up" << std::endl;
    std::cout << "[CosmicStudies constructor] verbose_output " << m_verb << std::endl;
    std::cout << "[CosmicStudies constructor] is_mcc9 " << m_is_mcc9 << std::endl;
    std::cout << "[CosmicStudies constructor] is_data " << m_is_data << std::endl;
    std::cout << "[CosmicStudies constructor] is_slimmed " << m_slimmed << std::endl;
    std::cout << "[CosmicStudies constructor] is_true_nu " << m_is_true_nu << std::endl;

    //// Tree for every subrun
    fPOTTree = tfs->make<TTree>("pot", "POT Tree");
    fPOTTree->Branch("run", &fRun_sr, "run/i");
    fPOTTree->Branch("subrun", &fSubrun_sr, "subrun/i");
    fPOTTree->Branch("pot", &fPot, "pot/F");
    fPOTTree->Branch("n_events", &fNevents, "n_events/i");

    //// Tree for every event
    fEventTree = tfs->make<TTree>("Event", "Event Tree");
    fEventTree->Branch("event", &fEvent, "event/i");
    fEventTree->Branch("run", &fRun, "run/i");
    fEventTree->Branch("subrun", &fSubrun, "subrun/i");
    fEventTree->Branch("pot", &fPot, "pot/F");
    fEventTree->Branch("dataset_prescale_factor", &fDatasetPrescaleFactor, "dataset_prescale_factor/F");

    fEventTree->Branch("num_simplebeamflashes", &fNumSimpleBeamFlashes, "num_simplebeamflashes/i");
    fEventTree->Branch("num_opbeamflashes", &fNumOpBeamFlashes, "num_opbeamflashes/i");
    fEventTree->Branch("num_simplecosmicflashes", &fNumSimpleCosmicFlashes, "num_simplecosmicflashes/i");
    fEventTree->Branch("num_opcosmicflashes", &fNumOpCosmicFlashes, "num_opcosmicflashes/i");

    fEventTree->Branch("num_pfp", &fNumPfp, "num_pfp/i");
    fEventTree->Branch("num_pfp_saved", &fNumPfp_saved, "num_pfp_saved/i");
    if (!m_is_data)
    {
        fEventTree->Branch("num_mcp", &fNumMcp, "num_mcp/i");
        fEventTree->Branch("num_mcp_saved", &fNumMcp_saved, "num_mcp_saved/i");
    }
    if (m_is_true_nu)
    {
        fEventTree->Branch("num_nu", &fNum_nu, "num_nu/I");
        fEventTree->Branch("nu_vtx_x", "std::vector< float >", &fNu_vtx_x);
        fEventTree->Branch("nu_vtx_y", "std::vector< float >", &fNu_vtx_y);
        fEventTree->Branch("nu_vtx_z", "std::vector< float >", &fNu_vtx_z);
        fEventTree->Branch("nu_E", "std::vector< float >", &fNu_E);
        fEventTree->Branch("nu_time", "std::vector< float >", &fNu_time);
        fEventTree->Branch("nu_pdg_code", "std::vector< int >", &fNu_pdg_code);
        fEventTree->Branch("nu_ccnc", "std::vector< bool >", &fNu_ccnc);
    }

    //// Tree for MC particles
    if (!m_is_data)
    {
        fMCParticlesTree = tfs->make<TTree>("MCParticles", "MCParticles Tree");
        fMCParticlesTree->Branch("event", &fEvent, "event/i");
        fMCParticlesTree->Branch("run", &fRun, "run/i");
        fMCParticlesTree->Branch("subrun", &fSubrun, "subrun/i");
        fMCParticlesTree->Branch("num_mcp", &fNumMcp, "num_mcp/i");
        fMCParticlesTree->Branch("num_mcp_saved", &fNumMcp_saved, "num_mcp_saved/i");
        fMCParticlesTree->Branch("mc_energy", &fMc_E, "mc_energy/F");
        fMCParticlesTree->Branch("mc_pdg_code", &fMc_PdgCode, "mc_pdg_code/I");
        fMCParticlesTree->Branch("mc_status_code", &fMc_StatusCode, "mc_status_code/I");
        fMCParticlesTree->Branch("mc_process", &fMc_Process, "mc_process/i");
        fMCParticlesTree->Branch("mc_end_process", &fMc_EndProcess, "mc_end_process/i");
        fMCParticlesTree->Branch("mc_is_matched", &fMc_Matched, "mc_is_matched/O");
        fMCParticlesTree->Branch("mc_start_inside", &fMc_StartInside, "mc_start_inside/O");
        fMCParticlesTree->Branch("mc_end_inside", &fMc_EndInside, "mc_end_inside/O");
        fMCParticlesTree->Branch("mc_part_inside", &fMc_PartInside, "mc_part_inside/O");
        fMCParticlesTree->Branch("crt_crossed", &fCRT_crossed, "crt_crossed/O");
        fMCParticlesTree->Branch("mc_neutrino_origin", &fMc_kBeamNeutrino, "mc_neutrino_origin/O");
        fMCParticlesTree->Branch("mc_time", &fMc_Time, "mc_time/F");
        fMCParticlesTree->Branch("mc_startx", &fMc_StartX, "mc_startx/F");
        fMCParticlesTree->Branch("mc_starty", &fMc_StartY, "mc_starty/F");
        fMCParticlesTree->Branch("mc_startz", &fMc_StartZ, "mc_startz/F");
        fMCParticlesTree->Branch("mc_startx_tpc", &fMc_StartX_tpc, "mc_startx_tpc/F");
        fMCParticlesTree->Branch("mc_starty_tpc", &fMc_StartY_tpc, "mc_starty_tpc/F");
        fMCParticlesTree->Branch("mc_startz_tpc", &fMc_StartZ_tpc, "mc_startz_tpc/F");
        fMCParticlesTree->Branch("mc_endx", &fMc_EndX, "mc_endx/F");
        fMCParticlesTree->Branch("mc_endy", &fMc_EndY, "mc_endy/F");
        fMCParticlesTree->Branch("mc_endz", &fMc_EndZ, "mc_endz/F");
        fMCParticlesTree->Branch("mc_endx_tpc", &fMc_EndX_tpc, "mc_endx_tpc/F");
        fMCParticlesTree->Branch("mc_endy_tpc", &fMc_EndY_tpc, "mc_endy_tpc/F");
        fMCParticlesTree->Branch("mc_endz_tpc", &fMc_EndZ_tpc, "mc_endz_tpc/F");
        fMCParticlesTree->Branch("mc_startmomentumx", &fMc_StartMomentumX, "mc_startmomentumx/F");
        fMCParticlesTree->Branch("mc_startmomentumy", &fMc_StartMomentumY, "mc_startmomentumy/F");
        fMCParticlesTree->Branch("mc_startmomentumz", &fMc_StartMomentumZ, "mc_startmomentumz/F");
        fMCParticlesTree->Branch("mc_length", &fMc_Length, "mc_length/F");
        fMCParticlesTree->Branch("mc_length_tpc", &fMc_LengthTPC, "mc_length_tpc/F");
    }

    //// Tree for beam flashes
    fOpBeamFlashesTree = tfs->make<TTree>("OpBeamFlashes", "OpBeamFlashes Tree");
    fOpBeamFlashesTree->Branch("event", &fEvent, "event/i");
    fOpBeamFlashesTree->Branch("run", &fRun, "run/i");
    fOpBeamFlashesTree->Branch("subrun", &fSubrun, "subrun/i");
    fOpBeamFlashesTree->Branch("dataset_prescale_factor", &fDatasetPrescaleFactor, "dataset_prescale_factor/F");
    fOpBeamFlashesTree->Branch("num_flashes", &fNumOpBeamFlashes, "num_flashes/i");
    fOpBeamFlashesTree->Branch("flash_time", &fFlash_Time, "flash_time/F");
    fOpBeamFlashesTree->Branch("flash_difftime", &fFlash_DiffTime, "flash_difftime/F");
    fOpBeamFlashesTree->Branch("flash_totalPE", &fFlash_TotalPE, "flash_total_PE/i");
    fOpBeamFlashesTree->Branch("flash_z", &fFlash_Z, "flash_z/F");
    fOpBeamFlashesTree->Branch("flash_sz", &fFlash_sigmaZ, "flash_sz/F");
    fOpBeamFlashesTree->Branch("flash_y", &fFlash_Y, "flash_y/F");
    fOpBeamFlashesTree->Branch("flash_sy", &fFlash_sigmaY, "flash_sy/F");
    fOpBeamFlashesTree->Branch("flash_width", &fFlash_Width, "flash_width/F");
    fOpBeamFlashesTree->Branch("flash_abstime", &fFlash_AbsTime, "flash_abstime/F");
    fOpBeamFlashesTree->Branch("flash_num_PMT10percent", &fFlash_num10percentPMT, "flash_num_PMT10percent/i");

    fSimpleBeamFlashesTree = tfs->make<TTree>("SimpleBeamFlashes", "SimpleBeamFlashes Tree");
    fSimpleBeamFlashesTree->Branch("event", &fEvent, "event/i");
    fSimpleBeamFlashesTree->Branch("run", &fRun, "run/i");
    fSimpleBeamFlashesTree->Branch("subrun", &fSubrun, "subrun/i");
    fSimpleBeamFlashesTree->Branch("dataset_prescale_factor", &fDatasetPrescaleFactor, "dataset_prescale_factor/F");
    fSimpleBeamFlashesTree->Branch("num_flashes", &fNumSimpleBeamFlashes, "num_flashes/i");
    fSimpleBeamFlashesTree->Branch("flash_time", &fFlash_Time, "flash_time/F");
    fSimpleBeamFlashesTree->Branch("flash_difftime", &fFlash_DiffTime, "flash_difftime/F");
    fSimpleBeamFlashesTree->Branch("flash_totalPE", &fFlash_TotalPE, "flash_total_PE/i");
    fSimpleBeamFlashesTree->Branch("flash_z", &fFlash_Z, "flash_z/F");
    fSimpleBeamFlashesTree->Branch("flash_sz", &fFlash_sigmaZ, "flash_sz/F");
    fSimpleBeamFlashesTree->Branch("flash_y", &fFlash_Y, "flash_y/F");
    fSimpleBeamFlashesTree->Branch("flash_sy", &fFlash_sigmaY, "flash_sy/F");
    fSimpleBeamFlashesTree->Branch("flash_width", &fFlash_Width, "flash_width/F");
    fSimpleBeamFlashesTree->Branch("flash_abstime", &fFlash_AbsTime, "flash_abstime/F");
    fSimpleBeamFlashesTree->Branch("flash_num_PMT10percent", &fFlash_num10percentPMT, "flash_num_PMT10percent/i");

    //// Tree for cosmic flashes
    fOpCosmicFlashesTree = tfs->make<TTree>("OpCosmicFlashes", "OpCosmicFlashes Tree");
    fOpCosmicFlashesTree->Branch("event", &fEvent, "event/i");
    fOpCosmicFlashesTree->Branch("run", &fRun, "run/i");
    fOpCosmicFlashesTree->Branch("subrun", &fSubrun, "subrun/i");
    fOpCosmicFlashesTree->Branch("dataset_prescale_factor", &fDatasetPrescaleFactor, "dataset_prescale_factor/F");
    fOpCosmicFlashesTree->Branch("num_flashes", &fNumOpCosmicFlashes, "num_flashes/i");
    fOpCosmicFlashesTree->Branch("flash_time", &fFlash_Time, "flash_time/F");
    fOpCosmicFlashesTree->Branch("flash_difftime", &fFlash_DiffTime, "flash_difftime/F");
    fOpCosmicFlashesTree->Branch("flash_totalPE", &fFlash_TotalPE, "flash_total_PE/i");
    fOpCosmicFlashesTree->Branch("flash_z", &fFlash_Z, "flash_z/F");
    fOpCosmicFlashesTree->Branch("flash_sz", &fFlash_sigmaZ, "flash_sz/F");
    fOpCosmicFlashesTree->Branch("flash_y", &fFlash_Y, "flash_y/F");
    fOpCosmicFlashesTree->Branch("flash_sy", &fFlash_sigmaY, "flash_sy/F");
    fOpCosmicFlashesTree->Branch("flash_width", &fFlash_Width, "flash_width/F");
    fOpCosmicFlashesTree->Branch("flash_abstime", &fFlash_AbsTime, "flash_abstime/F");
    fOpCosmicFlashesTree->Branch("flash_num_PMT10percent", &fFlash_num10percentPMT, "flash_num_PMT10percent/i");

    fSimpleCosmicFlashesTree = tfs->make<TTree>("SimpleCosmicFlashes", "SimpleCosmicFlashes Tree");
    fSimpleCosmicFlashesTree->Branch("event", &fEvent, "event/i");
    fSimpleCosmicFlashesTree->Branch("run", &fRun, "run/i");
    fSimpleCosmicFlashesTree->Branch("subrun", &fSubrun, "subrun/i");
    fSimpleCosmicFlashesTree->Branch("num_mcp", &fNumMcp, "num_mcp/i"); // This field is needed to distinguish mcc9 events with the same event/subrun/run tag.
    fSimpleCosmicFlashesTree->Branch("dataset_prescale_factor", &fDatasetPrescaleFactor, "dataset_prescale_factor/F");
    fSimpleCosmicFlashesTree->Branch("num_flashes", &fNumSimpleCosmicFlashes, "num_flashes/i");
    fSimpleCosmicFlashesTree->Branch("flash_time", &fFlash_Time, "flash_time/F");
    fSimpleCosmicFlashesTree->Branch("flash_difftime", &fFlash_DiffTime, "flash_difftime/F");
    fSimpleCosmicFlashesTree->Branch("flash_totalPE", &fFlash_TotalPE, "flash_total_PE/i");
    fSimpleCosmicFlashesTree->Branch("flash_z", &fFlash_Z, "flash_z/F");
    fSimpleCosmicFlashesTree->Branch("flash_sz", &fFlash_sigmaZ, "flash_sz/F");
    fSimpleCosmicFlashesTree->Branch("flash_y", &fFlash_Y, "flash_y/F");
    fSimpleCosmicFlashesTree->Branch("flash_sy", &fFlash_sigmaY, "flash_sy/F");
    fSimpleCosmicFlashesTree->Branch("flash_width", &fFlash_Width, "flash_width/F");
    fSimpleCosmicFlashesTree->Branch("flash_abstime", &fFlash_AbsTime, "flash_abstime/F");
    fSimpleCosmicFlashesTree->Branch("flash_num_PMT10percent", &fFlash_num10percentPMT, "flash_num_PMT10percent/i");

    //// Tree for the PF particles
    fPFParticlesTree = tfs->make<TTree>("PFParticles", "PFParticles Tree");
    fPFParticlesTree->Branch("event", &fEvent, "event/i");
    fPFParticlesTree->Branch("run", &fRun, "run/i");
    fPFParticlesTree->Branch("subrun", &fSubrun, "subrun/i");
    fPFParticlesTree->Branch("num_pfp", &fNumPfp, "num_pfp/i");
    fPFParticlesTree->Branch("num_mcp", &fNumMcp, "num_mcp/i");
    fPFParticlesTree->Branch("num_mcp_saved", &fNumMcp_saved, "num_mcp_saved/i");
    fPFParticlesTree->Branch("num_flashes", &fNumSimpleBeamFlashes, "num_flashes/i");
    fPFParticlesTree->Branch("pdg_code", &fPdgCode, "pdg_code/I");
    fPFParticlesTree->Branch("num_daughters", &fNumDaughters, "num_daughters/i");
    fPFParticlesTree->Branch("is_primary", &fIsPrimary, "is_primary/O");
    fPFParticlesTree->Branch("n_hits", &fNhits, "n_hits/i");
    fPFParticlesTree->Branch("n_clusters", &fNclusters, "n_clusters/i");
    fPFParticlesTree->Branch("pfp_vx", &fVx, "pfp_vx/F");
    fPFParticlesTree->Branch("pfp_vy", &fVy, "pfp_vy/F");
    fPFParticlesTree->Branch("pfp_vz", &fVz, "pfp_vz/F");
    //track
    fPFParticlesTree->Branch("track_valid", &fTrack_Valid, "track_valid/O");
    fPFParticlesTree->Branch("track_startx", &fTrack_StartX, "track_startx/F");
    fPFParticlesTree->Branch("track_starty", &fTrack_StartY, "track_starty/F");
    fPFParticlesTree->Branch("track_startz", &fTrack_StartZ, "track_startz/F");
    fPFParticlesTree->Branch("track_endx", &fTrack_EndX, "track_endx/F");
    fPFParticlesTree->Branch("track_endy", &fTrack_EndY, "track_endy/F");
    fPFParticlesTree->Branch("track_endz", &fTrack_EndZ, "track_endz/F");
    fPFParticlesTree->Branch("track_length", &fTrack_Length, "track_length/F");
    fPFParticlesTree->Branch("track_hasmomentum", &fTrack_HasMomentum, "track_hasmomentum/O");
    fPFParticlesTree->Branch("track_startmomentumx", &fTrack_StartMomentumX, "track_startmomentumx/F");
    fPFParticlesTree->Branch("track_startmomentumy", &fTrack_StartMomentumY, "track_startmomentumy/F");
    fPFParticlesTree->Branch("track_startmomentumz", &fTrack_StartMomentumZ, "track_startmomentumz/F");
    fPFParticlesTree->Branch("track_endmomentumx", &fTrack_EndMomentumX, "track_endmomentumx/F");
    fPFParticlesTree->Branch("track_endmomentumy", &fTrack_EndMomentumY, "track_endmomentumy/F");
    fPFParticlesTree->Branch("track_endmomentumz", &fTrack_EndMomentumZ, "track_endmomentumz/F");
    fPFParticlesTree->Branch("track_theta", &fTrack_Theta, "track_theta/F");
    fPFParticlesTree->Branch("track_phi", &fTrack_Phi, "track_phi/F");
    fPFParticlesTree->Branch("track_zenith", &fTrack_ZenithAngle, "track_zeninth/F");
    fPFParticlesTree->Branch("track_azimuth", &fTrack_AzimuthAngle, "track_azimuth/F");
    fPFParticlesTree->Branch("track_mcs_momentum", &fTrack_MCS_mom, "track_mcs_momentum/F");
    fPFParticlesTree->Branch("track_mcs_mom_err", &fTrack_MCS_err, "track_mcs_mom_err/F");
    fPFParticlesTree->Branch("track_mcs_likelihood", &fTrack_MCS_ll, "track_mcs_likelihood/F");
    fPFParticlesTree->Branch("track_mcs_energy", &fTrack_MCS_energy, "track_mcs_energy/F");

    if (!m_is_data)
    {
        fPFParticlesTree->Branch("track_matched_pdgcode", &fTrack_matched_PdgCode, "track_matched_pdgcode/I");
        fPFParticlesTree->Branch("track_matched_energy", &fTrack_matchedE, "track_matched_energy/F");
        fPFParticlesTree->Branch("track_matched_kBeamNeutrino", &fTrack_matched_kBeamNeutrino, "track_matched_kBeamNeutrino/O");
        fPFParticlesTree->Branch("track_matched_time", &fTrack_matched_Time, "track_matched_time/F");
        fPFParticlesTree->Branch("track_matched_process", &fTrack_matched_Process, "track_matched_process/i");
        fPFParticlesTree->Branch("track_matched_end_process", &fTrack_matched_EndProcess, "track_matched_end_process/i");
        fPFParticlesTree->Branch("track_matched_startinside", &fTrack_matched_StartInside, "track_matched_startinside/O");
        fPFParticlesTree->Branch("track_matched_endinside", &fTrack_matched_EndInside, "track_matched_endinside/O");
        fPFParticlesTree->Branch("track_matched_partinside", &fTrack_matched_PartInside, "track_matched_partinside/O");
        fPFParticlesTree->Branch("track_matched_startx", &fTrack_matched_StartX, "track_matched_startx/F");
        fPFParticlesTree->Branch("track_matched_starty", &fTrack_matched_StartY, "track_matched_starty/F");
        fPFParticlesTree->Branch("track_matched_startz", &fTrack_matched_StartZ, "track_matched_startz/F");
        fPFParticlesTree->Branch("track_matched_endx", &fTrack_matched_EndX, "track_matched_endx/F");
        fPFParticlesTree->Branch("track_matched_endy", &fTrack_matched_EndY, "track_matched_endy/F");
        fPFParticlesTree->Branch("track_matched_endz", &fTrack_matched_EndZ, "track_matched_endz/F");
        fPFParticlesTree->Branch("track_matched_startx_sce", &fTrack_matched_StartX_sce, "track_matched_startx_sce/F");
        fPFParticlesTree->Branch("track_matched_starty_sce", &fTrack_matched_StartY_sce, "track_matched_starty_sce/F");
        fPFParticlesTree->Branch("track_matched_startz_sce", &fTrack_matched_StartZ_sce, "track_matched_startz_sce/F");
        fPFParticlesTree->Branch("track_matched_endx_sce", &fTrack_matched_EndX_sce, "track_matched_endx_sce/F");
        fPFParticlesTree->Branch("track_matched_endy_sce", &fTrack_matched_EndY_sce, "track_matched_endy_sce/F");
        fPFParticlesTree->Branch("track_matched_endz_sce", &fTrack_matched_EndZ_sce, "track_matched_endz_sce/F");
        fPFParticlesTree->Branch("track_matched_length_tpc", &fTrack_matched_LengthTPC, "track_matched_length_tpc/F");
        fPFParticlesTree->Branch("track_matched_length_sce", &fTrack_matched_Length_sce, "track_matched_length_sce/F");
        fPFParticlesTree->Branch("track_matched_startmomentumx", &fTrack_matched_StartMomentumX, "track_matched_startmomentumx/F");
        fPFParticlesTree->Branch("track_matched_startmomentumy", &fTrack_matched_StartMomentumY, "track_matched_startmomentumy/F");
        fPFParticlesTree->Branch("track_matched_startmomentumz", &fTrack_matched_StartMomentumZ, "track_matched_startmomentumz/F");
    }

    fClustersTree = tfs->make<TTree>("Clusters", "Clusters Tree");
    fClustersTree->Branch("event", &fEvent, "event/i");
    fClustersTree->Branch("run", &fRun, "run/i");
    fClustersTree->Branch("subrun", &fSubrun, "subrun/i");
    fClustersTree->Branch("pdg_code", &fPdgCode, "pdg_code/I");
    fClustersTree->Branch("charge", &fClusterCharge, "charge/F");
    fClustersTree->Branch("width", &fClusterWidth, "width/F");
    fClustersTree->Branch("n_hits", &fClusterNhits, "n_hits/i");
    fClustersTree->Branch("plane", &fClusterPlane, "plane/i");
    fClustersTree->Branch("position", &fClusterPosition, "position/F");

    if (!m_is_data)
    {
        fCRTcrossTree = tfs->make<TTree>("CRTcross", "CRTcross Tree");
        fCRTcrossTree->Branch("event", &fEvent, "event/i");
        fCRTcrossTree->Branch("run", &fRun, "run/i");
        fCRTcrossTree->Branch("subrun", &fSubrun, "subrun/i");
        fCRTcrossTree->Branch("num_mcp", &fNumMcp, "num_mcp/i"); // This field is needed to distinguish mcc9 events with the same event/subrun/run tag.
        fCRTcrossTree->Branch("cross_x", &fCrossX, "cross_x/F");
        fCRTcrossTree->Branch("cross_y", &fCrossY, "cross_y/F");
        fCRTcrossTree->Branch("cross_z", &fCrossZ, "cross_z/F");
        fCRTcrossTree->Branch("cross_time", &fCrossT, "cross_time/F");
        fCRTcrossTree->Branch("cross_E", &fCrossE, "cross_/F");
        fCRTcrossTree->Branch("mc_time", &fMc_Time, "cross_x/F");
        fCRTcrossTree->Branch("mc_pdg_code", &fMc_PdgCode, "mc_pdg_code/I");
        fCRTcrossTree->Branch("mc_process", &fMc_Process, "mc_process/i");
        fCRTcrossTree->Branch("mc_energy", &fMc_E, "mc_energy/F");
        fCRTcrossTree->Branch("mc_neutrino_origin", &fMc_kBeamNeutrino, "mc_neutrino_origin/O");
    }
}

void CosmicStudies::reconfigure(fhicl::ParameterSet const &p)
{
    m_pfp_producer = p.get<std::string>("pfp_producer", "pandoraCosmic");
    m_spacepoint_producer = p.get<std::string>("spacepoint_producer", "pandoraCosmic");
    m_hitfinder_producer = p.get<std::string>("hitfinder_producer", "gaushit");
    m_geant_producer = p.get<std::string>("geant_producer", "largeant");
    m_hit_mcp_producer = p.get<std::string>("hit_mcp_producer", "gaushitTruthMatch");
    m_cosmic_simpleflash_producer = p.get<std::string>("cosmic_simpleflash_producer", "simpleFlashCosmic");
    m_cosmic_opflash_producer = p.get<std::string>("cosmic_opflash_producer", "opflashCosmic");
    m_cosmic_ophit_producer = p.get<std::string>("cosmic_ophit_producer", "ophitCosmic");
    m_beam_simpleflash_producer = p.get<std::string>("beam_simpleflash_producer", "simpleFlashBeam");
    m_beam_opflash_producer = p.get<std::string>("beam_opflash_producer", "opflashBeam");
    m_beam_ophit_producer = p.get<std::string>("beam_ophit_producer", "ophitBeam");

    m_verb = p.get<bool>("verbose_output", true);
    m_is_mcc9 = p.get<bool>("is_mcc9", false);
    m_is_data = p.get<bool>("is_data", false);
    m_is_true_nu = p.get<bool>("is_true_nu", false);
    m_slimmed = p.get<bool>("is_slimmed", false);
}

// Clear once per event
void CosmicStudies::clear()
{
    fDatasetPrescaleFactor = 1;
    //fNevents = 0; this is cleared every endsubrun
    matchedParticles.clear();
    matchedMCParticles.clear();

    fRun = 0;
    fSubrun = 0;
    fEvent = 0;

    fNumPfp = 0;
    fNumMcp = 0;
    fNumPfp_saved = 0;
    fNumMcp_saved = 0;
    fNumSimpleBeamFlashes = 0;
    fNumOpBeamFlashes = 0;
    fNumSimpleCosmicFlashes = 0;
    fNumOpCosmicFlashes = 0;
    fNum_nu = 0;
    fNumCross = 0;

    fNu_vtx_x.clear();
    fNu_vtx_y.clear();
    fNu_vtx_z.clear();
    fNu_E.clear();
    fNu_time.clear();
    fNu_pdg_code.clear();
    fNu_ccnc.clear();

    // Make sure this is set to false if the event does not contain a neutrino, in that case this field should never be changed
    fMc_kBeamNeutrino = false;
    fTrack_matched_kBeamNeutrino = false;
}

// Clear once per cluster
void CosmicStudies::clear_Cluster()
{
    // fClustersTree;
    fClusterCharge = -9999;
    fClusterWidth = -9999;
    fClusterPosition = -9999;
    fClusterNhits = 0;
    fClusterPlane = 0;
}

// Clear once per PFP
void CosmicStudies::clear_PFParticle()
{
    // fPFParticlesTree;
    fPdgCode = 0;
    fNumDaughters = 0;
    fIsPrimary = false;
    fVx = -9999;
    fVy = -9999;
    fVz = -9999;
    fNhits = 0;
    fNclusters = 0;
    fNhitsU = 0;
    fNhitsV = 0;
    fNhitsY = 0;
    //// track
    fTrack_Valid = false;
    fTrack_HasMomentum = false;
    fTrack_StartX = -9999;
    fTrack_StartY = -9999;
    fTrack_StartZ = -9999;
    fTrack_EndX = -9999;
    fTrack_EndY = -9999;
    fTrack_EndZ = -9999;
    fTrack_Length = -9999;
    fTrack_StartMomentumX = -9999;
    fTrack_StartMomentumY = -9999;
    fTrack_StartMomentumZ = -9999;
    fTrack_EndMomentumX = -9999;
    fTrack_EndMomentumY = -9999;
    fTrack_EndMomentumZ = -9999;
    fTrack_Theta = -9999;
    fTrack_Phi = -9999;
    fTrack_ZenithAngle = -9999;
    fTrack_AzimuthAngle = -9999;
    // Reco-truth matched
    fTrack_matched_PdgCode = 0; // this means that the PFParticel was not matched!
}

// Clear once per beam flash
void CosmicStudies::clear_Flashes()
{
    fFlash_Time = -9999;
    fFlash_DiffTime = -9999;
    fFlash_TotalPE = 0;
    fFlash_Z = -9999;
    fFlash_sigmaZ = -9999;
    fFlash_Y = -9999;
    fFlash_sigmaY = -9999;
    fFlash_num10percentPMT = 0;
    fFlash_Width = -9999;
    fFlash_AbsTime = -9999;
}

art::Ptr<simb::MCTruth> CosmicStudies::TrackIDToMCTruth(art::Event const &e, std::string m_geant_producer, int geant_track_id)
{
    lar_pandora::MCTruthToMCParticles truthToParticles;
    lar_pandora::MCParticlesToMCTruth particlesToTruth;

    lar_pandora::LArPandoraHelper::CollectMCParticles(e, m_geant_producer, truthToParticles, particlesToTruth);

    for (auto iter : particlesToTruth)
    {
        if (iter.first->TrackId() == geant_track_id)
        {
            return iter.second;
        }
    }

    art::Ptr<simb::MCTruth> null_ptr;
    return null_ptr;
}

DEFINE_ART_MODULE(CosmicStudies)

#endif