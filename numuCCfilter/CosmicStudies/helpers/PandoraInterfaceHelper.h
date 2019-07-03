#ifndef PANDORAINTERFACEHELPER_H
#define PANDORAINTERFACEHELPER_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

namespace lar_pandora
{
typedef std::map<art::Ptr<recob::PFParticle>, art::Ptr<simb::MCParticle>> PFParticlesToMCParticles;
}

typedef std::map<art::Ptr<recob::PFParticle>, unsigned int> RecoParticleToNMatchedHits;
typedef std::map<art::Ptr<simb::MCParticle>, RecoParticleToNMatchedHits> ParticleMatchingMap;
typedef std::set<art::Ptr<recob::PFParticle>> PFParticleSet;
typedef std::set<art::Ptr<simb::MCParticle>> MCParticleSet;

class PandoraInterfaceHelper
{
public:
  PandoraInterfaceHelper();
  ~PandoraInterfaceHelper() = default;

  void reconfigure(fhicl::ParameterSet const &pset);

  /**
    * @brief Travers the tree of the daughters of a PFParticle
    *
    * @param pfparticles PFParticles handle
    * @param top_index Index of the parent
    * @param unordered_daugthers Vector of PFParticles daughters
    */
  void traversePFParticleTree(
      art::Handle<std::vector<recob::PFParticle>> pfparticles,
      size_t top_index, std::vector<size_t> &unordered_daugthers,
      std::string m_pfp_producer);

  /**
     * @brief Measures the three-dimensional center of the deposited charge for a
     * PFParticle
     *
     * @param ipf Index of the PFParticle
     * @param pfparticles PFParticles handle
     * @param evt art Event
     * @return vector with: lowest x_sps, center in y, z, and total deposited charge on the collection plane.
     */
  std::vector<double> calculateChargeCenter(
      size_t ipf,
      art::Handle<std::vector<recob::PFParticle>> pfparticles,
      const art::Event &evt,
      std::string m_pfp_producer);

  void get_daughter_tracks(std::vector<size_t> pf_ids, const art::Event &evt,
                           std::vector<art::Ptr<recob::Track>> &tracks,
                           std::string m_pfp_producer);

  void get_daughter_showers(std::vector<size_t> pf_ids, const art::Event &evt,
                            std::vector<art::Ptr<recob::Shower>> &showers,
                            std::string m_pfp_producer);

  /**
    *  @brief Returns matching between true and reconstructed particles
    *
    *  @param matchedParticles the output matches between reconstructed and true particles
    */
  void GetRecoToTrueMatches(lar_pandora::PFParticlesToMCParticles &matchedParticles);

  /**
     *  @brief Configure function parameters (call this function first)
     *
     *  @param e the art::Event
     *  @param m_pfp_producer the PFParticle producer label
     *  @param m_spacepoint_producer the SpacePoint producer label
     *  @param m_hitfinder_producer the Hit producer label
     *  @param m_geant_producer The Geant4 producer label
     */
  void Configure(art::Event const &e,
                 std::string m_pfp_producer,
                 std::string m_spacepoint_producer,
                 std::string m_hitfinder_producer,
                 std::string m_geant_producer,
                 std::string m_hit_mcp_producer);

  void Configure(art::Event const &e,
                 std::string m_pfp_producer,
                 std::string m_spacepoint_producer,
                 std::string m_hitfinder_producer,
                 std::string m_geant_producer);

  art::Ptr<simb::MCTruth> TrackIDToMCTruth(art::Event const &e, std::string _geant_producer, int geant_track_id);
  /**
     *  @brief Collect a vector of MCParticle objects from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the truth information in the event
     *  @param particleVector the output vector of MCParticle objects
     */
  void CollectMCParticles(const art::Event &evt,
                          const std::string &label,
                          lar_pandora::MCParticleVector &particleVector);
  /**
     *  @brief Collect truth information from the ART event record
     *
     *  @param evt the ART event record
     *  @param label the label for the truth information in the event
     *  @param truthToParticles output map from MCTruth to MCParticle objects
     *  @param particlesToTruth output map from MCParticle to MCTruth objects
     */
  void CollectMCParticles(const art::Event &evt,
                          const std::string &label,
                          lar_pandora::MCTruthToMCParticles &truthToParticles,
                          lar_pandora::MCParticlesToMCTruth &particlesToTruth);

protected:
  lar_pandora::HitsToMCParticles m_hit_to_mcps_map; ///< A map from recon hits to MCParticles
  lar_pandora::PFParticlesToHits m_pfp_to_hits_map; ///< A map from PFParticles to recon hits
private:
  bool m_configured = false;
  bool m_isOverlaidSample = false;
  bool m_debug = false;
  bool m_verbose = false;
};

#endif