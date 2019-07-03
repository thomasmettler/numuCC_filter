#ifndef PANDORAINTERFACEHELPER_CXX
#define PANDORAINTERFACEHELPER_CXX

#include "PandoraInterfaceHelper.h"

void PandoraInterfaceHelper::reconfigure(fhicl::ParameterSet const &pset)
{
  m_isOverlaidSample = pset.get<bool>("isOverlaidSample", false);
}

PandoraInterfaceHelper::PandoraInterfaceHelper()
{
  m_configured = false;
}

void PandoraInterfaceHelper::get_daughter_tracks(std::vector<size_t> pf_ids,
                                                 const art::Event &evt,
                                                 std::vector<art::Ptr<recob::Track>> &tracks,
                                                 std::string m_pfp_producer)
{
  try
  {

    auto const &pfparticle_handle =
        evt.getValidHandle<std::vector<recob::PFParticle>>(m_pfp_producer);

    art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt, m_pfp_producer);

    for (auto const &pf_id : pf_ids)
    {
      auto const &track_obj = track_per_pfpart.at(pf_id);
      tracks.push_back(track_obj);
    }
  }
  catch (...)
  {
    std::cout << "[PandoraInterfaceHelper] "
              << "Error getting daughter tracks" << std::endl;
  }
}

// Old version of configure on non-slimmed files
void PandoraInterfaceHelper::Configure(art::Event const &e,
                                       std::string m_pfp_producer,
                                       std::string m_spacepoint_producer,
                                       std::string m_hitfinder_producer,
                                       std::string m_geant_producer)
{

  // Collect hits
  lar_pandora::HitVector hitVector;
  lar_pandora::LArPandoraHelper::CollectHits(e, m_hitfinder_producer, hitVector);

  // Collect PFParticles and match Reco Particles to Hits
  lar_pandora::PFParticleVector recoParticleVector;
  lar_pandora::PFParticleVector recoNeutrinoVector;
  lar_pandora::PFParticlesToHits pfp_to_hits_map;
  lar_pandora::HitsToPFParticles recoHitsToParticles;

  lar_pandora::LArPandoraHelper::CollectPFParticles(e, m_pfp_producer, recoParticleVector);
  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
  lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e,
                                                        m_pfp_producer,
                                                        m_spacepoint_producer,
                                                        pfp_to_hits_map,
                                                        recoHitsToParticles,
                                                        lar_pandora::LArPandoraHelper::kAddDaughters,
                                                        true); // Use clusters to go from pfp to hits

  if (m_verbose)
  {
    std::cout << "[McPfpMatch] RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;
    std::cout << "[McPfpMatch] RecoParticles: " << recoParticleVector.size() << std::endl;
  }

  // Collect MCParticles and match True Particles to Hits
  lar_pandora::MCParticleVector trueParticleVector;
  lar_pandora::MCTruthToMCParticles truthToParticles;
  lar_pandora::MCParticlesToMCTruth particlesToTruth;
  lar_pandora::MCParticlesToHits trueParticlesToHits;
  lar_pandora::HitsToMCParticles hit_to_mcps_map;

  lar_pandora::LArPandoraHelper::CollectMCParticles(e, m_geant_producer, trueParticleVector);
  lar_pandora::LArPandoraHelper::CollectMCParticles(e, m_geant_producer, truthToParticles, particlesToTruth);
  lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(e,
                                                        m_geant_producer,
                                                        hitVector,
                                                        trueParticlesToHits,
                                                        hit_to_mcps_map,
                                                        lar_pandora::LArPandoraHelper::kAddDaughters); // Consider daughters as independent mcps

  if (m_verbose)
  {
    std::cout << "[McPfpMatch] TrueParticles: " << particlesToTruth.size() << std::endl;
    std::cout << "[McPfpMatch] TrueEvents: " << truthToParticles.size() << std::endl;
  }

  // Now set the things we need for the future
  m_hit_to_mcps_map = hit_to_mcps_map;
  m_pfp_to_hits_map = pfp_to_hits_map;

  if (m_debug)
  {
    std::cout << "[McPfpMatch] This is event " << e.id().event() << std::endl;
    //art::ServiceHandle<cheat::BackTracker> bt;
    std::cout << "[McPfpMatch] Number of MCParticles matched to hits: " << trueParticlesToHits.size() << std::endl;
    for (const auto &iter : trueParticlesToHits)
    {
      const art::Ptr<simb::MCTruth> mc_truth = TrackIDToMCTruth(e, m_geant_producer, (iter.first)->TrackId()); //bt->TrackIDToMCTruth((iter.first)->TrackId());
      std::cout
          << "[McPfpMatch] MCParticle with pdg " << (iter.first)->PdgCode()
          << " and origin " << (mc_truth->Origin() == 1 ? "neutrino" : "cosmic")
          << " has " << (iter.second).size() << " hits ass." << std::endl;
      if (mc_truth->Origin() == 1)
      {
        lar_pandora::HitVector hits = (iter.second);
      }
    }
  }

  m_configured = true;
}

// New version of configure on slimmed files
void PandoraInterfaceHelper::Configure(art::Event const &e,
                                       std::string m_pfp_producer,
                                       std::string m_spacepoint_producer,
                                       std::string m_hitfinder_producer,
                                       std::string m_geant_producer,
                                       std::string m_hit_mcp_producer)
{

  lar_pandora::LArPandoraHelper::DaughterMode daughterMode = lar_pandora::LArPandoraHelper::kAddDaughters;

  // Collect hits

  lar_pandora::HitVector hitVector;
  //lar_pandora::LArPandoraHelper::CollectHits(e, m_hitfinder_producer, hitVector);
  art::Handle<std::vector<recob::Hit>> hit_h;

  e.getByLabel(m_hitfinder_producer, hit_h);
  if (!hit_h.isValid())
  {
    std::cout << "[McPfpMatch] Hit Handle is not valid." << std::endl;
    // throw std::exception();
  }

  art::fill_ptr_vector(hitVector, hit_h);

  art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> mcps_from_hit(hit_h, e, m_hit_mcp_producer);

  // Collect PFParticles and match Reco Particles to Hits
  lar_pandora::PFParticleVector recoParticleVector;
  lar_pandora::PFParticleVector recoNeutrinoVector;
  lar_pandora::PFParticlesToHits pfp_to_hits_map;
  lar_pandora::HitsToPFParticles recoHitsToParticles;

  lar_pandora::LArPandoraHelper::CollectPFParticles(e, m_pfp_producer, recoParticleVector);

  lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);

  try
  {
    lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e,
                                                          m_pfp_producer,
                                                          m_spacepoint_producer,
                                                          pfp_to_hits_map,
                                                          recoHitsToParticles,
                                                          daughterMode,
                                                          true); // Use clusters to go from pfp to hits
  }
  catch (...)
  {
    std::cout << "[McPfpMatch] BuildPFParticleHitMaps error" << std::endl;
  }

  m_verbose = false;
  if (m_verbose)
  {
    std::cout << "[McPfpMatch] RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;
    std::cout << "[McPfpMatch] RecoParticles: " << recoParticleVector.size() << std::endl;
  }

  // Collect MCParticles and match True Particles to Hits
  lar_pandora::MCParticleVector trueParticleVector;
  lar_pandora::MCTruthToMCParticles truthToParticles;
  lar_pandora::MCParticlesToMCTruth particlesToTruth;
  lar_pandora::MCParticlesToHits trueParticlesToHits;
  lar_pandora::HitsToMCParticles hit_to_mcps_map;

  if (!e.isRealData() || m_isOverlaidSample)
  {
    CollectMCParticles(e, m_geant_producer, trueParticleVector);
    CollectMCParticles(e, m_geant_producer, truthToParticles, particlesToTruth);

    // Construct a Particle Map (trackID to MCParticle)
    lar_pandora::MCParticleMap particleMap;

    for (lar_pandora::MCTruthToMCParticles::const_iterator iter1 = truthToParticles.begin(), iterEnd1 = truthToParticles.end(); iter1 != iterEnd1; ++iter1)
    {
      const lar_pandora::MCParticleVector &particleVector = iter1->second;
      for (lar_pandora::MCParticleVector::const_iterator iter2 = particleVector.begin(), iterEnd2 = particleVector.end(); iter2 != iterEnd2; ++iter2)
      {
        const art::Ptr<simb::MCParticle> particle = *iter2;
        particleMap[particle->TrackId()] = particle;
      }
    }

    // Loop over the hits, get the ass MCP, and then tru to link

    std::vector<art::Ptr<simb::MCParticle>> mcp_v;
    std::vector<anab::BackTrackerHitMatchingData const *> match_v;
    for (auto hit : hitVector)
    {

      mcp_v.clear();
      match_v.clear();
      mcps_from_hit.get(hit.key(), mcp_v, match_v);

      double max_energy = -1;
      int best_match_id = -1;
      for (size_t m = 0; m < match_v.size(); m++)
      {
        double this_energy = match_v[m]->energy;
        if (this_energy > max_energy)
        {
          best_match_id = m;
          max_energy = this_energy;
        }
      }

      if (best_match_id > -1)
      {
        try
        {
          const art::Ptr<simb::MCParticle> thisParticle = mcp_v.at(best_match_id);
          const art::Ptr<simb::MCParticle> primaryParticle(lar_pandora::LArPandoraHelper::GetFinalStateMCParticle(particleMap, thisParticle));
          const art::Ptr<simb::MCParticle> selectedParticle((lar_pandora::LArPandoraHelper::kAddDaughters == daughterMode) ? primaryParticle : thisParticle);

          if ((lar_pandora::LArPandoraHelper::kIgnoreDaughters == daughterMode) && (selectedParticle != primaryParticle))
            continue;

          if (!(lar_pandora::LArPandoraHelper::IsVisible(selectedParticle)))
            continue;

          hit_to_mcps_map[hit] = selectedParticle;
        }
        catch (...)
        {
          std::cout << "[PandoraInterfaceHelper] "
                    << "Error in the loop of the hits" << std::endl;
          continue;
        }
      }
    }
  }

  // Now set the things we need for the future
  m_hit_to_mcps_map = hit_to_mcps_map;
  m_pfp_to_hits_map = pfp_to_hits_map;

  // std::cout << "hit_to_mcps_map size " << hit_to_mcps_map.size() << std::endl;

  m_configured = true;
}

art::Ptr<simb::MCTruth> PandoraInterfaceHelper::TrackIDToMCTruth(art::Event const &e, std::string m_geant_producer, int geant_track_id)
{

  lar_pandora::MCTruthToMCParticles truthToParticles;
  lar_pandora::MCParticlesToMCTruth particlesToTruth;

  if (!e.isRealData() || m_isOverlaidSample)
  {
    CollectMCParticles(e, m_geant_producer, truthToParticles, particlesToTruth);

    for (auto iter : particlesToTruth)
    {
      if (iter.first->TrackId() == geant_track_id)
      {
        return iter.second;
      }
    }
  }
  art::Ptr<simb::MCTruth> null_ptr;
  return null_ptr;
}

void PandoraInterfaceHelper::CollectMCParticles(const art::Event &evt, const std::string &label, lar_pandora::MCParticleVector &particleVector)
{
  if (evt.isRealData() && !m_isOverlaidSample)
    std::cout << " PandoraCollector::CollectMCParticles --- Trying to access MC truth from real data ";

  art::Handle<std::vector<simb::MCParticle>> theParticles;
  evt.getByLabel(label, theParticles);
  if (!theParticles.isValid())
  {
    std::cout << "[LArPandora]"
              << "  Failed to find MC particles... " << std::endl;
    return;
  }
  else
  {
    if (m_debug)
    {
      std::cout << "[LArPandora]"
                << "  Found: " << theParticles->size() << " MC particles " << std::endl;
    }
  }
  for (unsigned int i = 0; i < theParticles->size(); ++i)
  {
    const art::Ptr<simb::MCParticle> particle(theParticles, i);
    particleVector.push_back(particle);
  }
}

void PandoraInterfaceHelper::get_daughter_showers(
    std::vector<size_t> pf_ids, const art::Event &evt,
    std::vector<art::Ptr<recob::Shower>> &showers, std::string m_pfp_producer)
{

  auto const &pfparticle_handle =
      evt.getValidHandle<std::vector<recob::PFParticle>>(m_pfp_producer);

  art::FindOneP<recob::Shower> shower_per_pfpart(pfparticle_handle, evt, m_pfp_producer);

  for (auto const &pf_id : pf_ids)
  {
    try
    {

      auto const &shower_obj = shower_per_pfpart.at(pf_id);
      showers.push_back(shower_obj);
    }
    catch (...)
    {
      std::cout << "[PandoraLEE] "
                << "Error getting the shower" << std::endl;
    }
  }
}

//___________________________________________________________________________________________________
void PandoraInterfaceHelper::GetRecoToTrueMatches(lar_pandora::PFParticlesToMCParticles &matchedParticles)
{
  bool m_debug = false;

  if (!m_configured)
  {
    std::cout << "Call to " << __PRETTY_FUNCTION__ << " whitout having done configuration. Abort." << std::endl;
    // throw std::exception();
  }

  // Loop over the reco particles
  for (auto iter1 : m_pfp_to_hits_map)
  {

    // The PFParticle
    const art::Ptr<recob::PFParticle> recoParticle = iter1.first;

    if (m_debug)
      std::cout << "[PandoraInterfaceHelper::GetRecoToTrueMatches] Looking at PFP with ID " << recoParticle->Self() << std::endl;

    // The PFParticle's hits
    const lar_pandora::HitVector &hitVector = iter1.second;

    if (m_debug)
      std::cout << "[PandoraInterfaceHelper::GetRecoToTrueMatches] \t This PFP has " << hitVector.size() << " hits." << std::endl;

    lar_pandora::MCParticlesToHits truthContributionMap;
    // Loop over all the hits associated to this reco particle
    for (auto hit : hitVector)
    {

      // Find the MCParticle that share this same hit (if any)
      auto iter3 = m_hit_to_mcps_map.find(hit);
      if (m_hit_to_mcps_map.end() == iter3)
        continue;

      // If exists, get the MCParticle
      const art::Ptr<simb::MCParticle> trueParticle = iter3->second;

      if (m_debug)
        std::cout << "[PandoraInterfaceHelper::GetRecoToTrueMatches] \t Found a hit shared with MCParticle with PDG " << trueParticle->PdgCode() << std::endl;

      // This map will contain all the true particles that match some or all of the hits of the reco particle
      truthContributionMap[trueParticle].push_back(hit);
    }

    // Now we want to find the true particle that has more hits in common with this reco particle than the others
    lar_pandora::MCParticlesToHits::const_iterator mIter = truthContributionMap.end();

    for (lar_pandora::MCParticlesToHits::const_iterator iter4 = truthContributionMap.begin(), iterEnd4 = truthContributionMap.end();
         iter4 != iterEnd4; ++iter4)
    {
      if ((truthContributionMap.end() == mIter) || (iter4->second.size() > mIter->second.size()))
      {
        mIter = iter4;
      }
    }

    if (truthContributionMap.end() != mIter)
    {
      const art::Ptr<simb::MCParticle> trueParticle = mIter->first;

      if (m_debug)
        std::cout << "[PandoraInterfaceHelper::GetRecoToTrueMatches] \t >>> Match found with MCParticle with PDG " << trueParticle->PdgCode() << std::endl;

      // Emplace into the output map
      matchedParticles[recoParticle] = trueParticle;
    }

  } // m_pfp_to_hits_map loop ends
}

void PandoraInterfaceHelper::traversePFParticleTree(
    const art::Handle<std::vector<recob::PFParticle>> pfparticles,
    size_t top_index, std::vector<size_t> &unordered_daugthers,
    std::string m_pfp_producer)
{

  // This is a tree-traversal algorithm.  It returns the index of the top
  // particle, plus the index
  // of all daughter particles.

  // This is a recursive algorithm, so it needs a break clause:
  unordered_daugthers.push_back(top_index);

  if (pfparticles->at(top_index).Daughters().size() == 0)
  {
    return;
  }

  // Else, go through the tree:
  for (size_t i = 0; i < pfparticles->at(top_index).Daughters().size(); i++)
  {
    traversePFParticleTree(pfparticles,
                           pfparticles->at(top_index).Daughters().at(i),
                           unordered_daugthers, m_pfp_producer);
  }

  return;
}

void PandoraInterfaceHelper::CollectMCParticles(const art::Event &evt,
                                                const std::string &label,
                                                lar_pandora::MCTruthToMCParticles &truthToParticles,
                                                lar_pandora::MCParticlesToMCTruth &particlesToTruth)
{
  if (evt.isRealData() && !m_isOverlaidSample)
    std::cout << " PandoraCollector::CollectMCParticles --- Trying to access MC truth from real data ";

  art::Handle<std::vector<simb::MCParticle>> theParticles;
  evt.getByLabel(label, theParticles);

  if (!theParticles.isValid())
  {
    std::cout << "[LArPandora]"
              << "  Failed to find MC particles... " << std::endl;
    return;
  }
  else
  {
    if (m_debug)
    {
      std::cout << "[LArPandora]"
                << "  Found: " << theParticles->size() << " MC particles " << std::endl;
    }
  }

  art::FindOneP<simb::MCTruth> theTruthAssns(theParticles, evt, label);
  for (unsigned int i = 0, iEnd = theParticles->size(); i < iEnd; ++i)
  {
    const art::Ptr<simb::MCParticle> particle(theParticles, i);
    const art::Ptr<simb::MCTruth> truth(theTruthAssns.at(i));
    truthToParticles[truth].push_back(particle);
    particlesToTruth[particle] = truth;
  }
}

// Method to calculate the total the center for a parent particle (index of
// neutrino pfp)
std::vector<double> PandoraInterfaceHelper::calculateChargeCenter(
    size_t top_particle_index,
    const art::Handle<std::vector<recob::PFParticle>> pfparticles,
    const art::Event &evt,
    std::string m_pfp_producer)
{

  // First, get the indexes of pfparticles that are in the hierarchy of this
  // particle:
  std::vector<size_t> daughters;
  daughters.reserve(50);
  traversePFParticleTree(pfparticles, top_particle_index, daughters, m_pfp_producer);

  // Get the associations from pfparticle to spacepoint
  auto const &spacepoint_handle =
      evt.getValidHandle<std::vector<recob::SpacePoint>>(m_pfp_producer);

  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticles, evt, m_pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, m_pfp_producer);

  // Variables for the total weight and center of charge
  double totalweight = 0;
  std::vector<double> chargecenter;
  chargecenter.resize(3);
  double min_x_sps_temp = 9999; //random big value that will be overwritten

  // Loop over the pfparticles, get their space points, and compute the weighted
  // average:

  for (auto &m_i_pfp : daughters)
  {

    // Get the associated spacepoints:
    std::vector<art::Ptr<recob::SpacePoint>> spcpnts =
        spcpnts_per_pfpart.at(m_i_pfp);

    // Loop over the spacepoints and get the associated hits:
    for (auto &m_sps : spcpnts)
    {
      auto xyz = m_sps->XYZ();

      if (xyz[0] > 0 && xyz[0] < min_x_sps_temp)
      {
        min_x_sps_temp = xyz[0];
      }

      std::vector<art::Ptr<recob::Hit>> hits = hits_per_spcpnts.at(m_sps.key());
      // Add the hits to the weighted average, if they are collection hits:
      for (auto &hit : hits)
      {
        if (hit->View() == geo::kZ)
        {
          // Collection hits only
          double weight = hit->Integral();
          chargecenter[0] += (xyz[0]) * weight;
          chargecenter[1] += (xyz[1]) * weight;
          chargecenter[2] += (xyz[2]) * weight;
          totalweight += weight;
          // break; // Exit the loop over hits
        } // if collection

      } // hits

    } // spacepoints

  } // pfparticles

  // Normalize;
  chargecenter[0] /= totalweight;
  chargecenter[1] /= totalweight;
  chargecenter[2] /= totalweight;

  // Store the data:
  std::vector<double> m_center_of_charge(4, 0);
  //_center_of_charge.SetX(chargecenter[0]);
  // IMPORTANT:
  // This function is necessary for optical selection.
  // Flashatching is returning the extimated minimal x position,
  // not the center. Therefore, also here teh minimum is returned and not the center.
  m_center_of_charge[0] = chargecenter[0];
  m_center_of_charge[1] = chargecenter[1];
  m_center_of_charge[2] = chargecenter[2];
  m_center_of_charge[3] = totalweight;

  return m_center_of_charge;
}

#endif