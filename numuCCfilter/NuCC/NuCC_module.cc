#include "NuCC.h"

void NuCC::endSubRun(const art::SubRun &subrun)
{
  if (!m_isData)
  {
    art::Handle<sumdata::POTSummary> potSummaryHandle;
    m_pot = subrun.getByLabel("generator", potSummaryHandle) ? static_cast<float>(potSummaryHandle->totpot) : 0.f;
    std::cout << "[LArPandoraExternalEventBuilding::endSubRun] Storing POT info!" << std::endl;
  }

  m_run = subrun.run();
  m_subrun = subrun.subRun();
  fSubrunTree->Fill();
}

void NuCC::analyze(art::Event const &evt)
{
  clearEvent();
  fRun = evt.run();
  fSubrun = evt.subRun();
  fEvent = evt.id().event();
  art::Timestamp evtTime = evt.time();
  fTimeHigh = evtTime.timeHigh();
  fTimeLow = evtTime.timeLow();
  std::cout << "[NuCC::analyze]: Run " << fRun << ", Subrun " << fSubrun << ", Event " << fEvent << std::endl;
  //std::cout << "[NuCC::analyze]: evt_time_sec " << evtTime.timeHigh() << ", evt_time_nsec " << evtTime.timeLow() << std::endl;

  larpandora.CollectPFParticleMetadata(evt, m_pfp_producer, pfparticles, particlesToMetadata);
  larpandora.BuildPFParticleMap(pfparticles, particleMap);

  if (pfparticles.size() == 0)
    std::cout << "[NuCC::FillReconstructed] No reconstructed PFParticles in event." << std::endl;
  else
  {
    larpandora.SelectNeutrinoPFParticles(pfparticles, pfneutrinos);
    if (pfneutrinos.size() != 1)
      std::cout << "[NuCC::FillReconstructed] Number of reconstructed neutrinos in event is " << pfneutrinos.size() << std::endl;
    else // We have a reconstructed neutrino
    {
      if (!m_isData)
      {
        FillReconTruthMatching(evt);
        FillTrueNu(evt);
        FillTrueNuDaughters(evt);
      }
      FillReconstructed(evt);
      fEventTree->Fill();
    }
  }
  std::cout << "\n\n";
}

void NuCC::FillReconstructed(art::Event const &evt)
{
  fNumPfp = pfparticles.size();
  // Load associations and collections
  lar_pandora::VertexVector vertexVector_dummy;
  lar_pandora::PFParticleVector particleVector_dummy;
  lar_pandora::SpacePointVector spacePointVector_dummy;
  larpandora.CollectVertices(evt, m_pfp_producer, vertexVector_dummy, particlesToVertices);
  larpandora.CollectPFParticles(evt, m_pfp_producer, particleVector_dummy, particlesToClusters);
  larpandora.CollectPFParticles(evt, m_pfp_producer, particleVector_dummy, particlesToSpacePoints);
  larpandora.CollectShowers(evt, m_pfp_producer, pfshowers, particlesToShowers);
  larpandora.CollectTracks(evt, m_pfp_producer, pftracks, particlesToTracks);
  lar_pandora::ClusterVector clusterVector_dummy;
  larpandora.CollectClusters(evt, m_pfp_producer, clusterVector_dummy, clustersToHits);
  //larpandora.CollectSpacePoints(evt, m_pfp_producer, spacePointVector_dummy, spacePointsToHits, hitsToSpacePoints);
  art::ValidHandle<std::vector<recob::Track>> trackHandle = evt.getValidHandle<std::vector<recob::Track> >(m_pfp_producer);
  const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle = evt.getValidHandle<std::vector<recob::MCSFitResult>>("pandoraMCSMu");
  const art::FindManyP<anab::ParticleID> trackPIDAssn(trackHandle, evt, "pandoracalipidSCE");
  if (!trackPIDAssn.isValid()){
    std::cout << "[NuCC::FillReconstructed] trackPIDAssn.isValid() == false" << std::endl;
  } 


  // Start filling information
  art::Ptr<recob::PFParticle> pfnu = pfneutrinos.front();
  fNu_PDG = pfnu->PdgCode();
  fNumPrimaryDaughters = pfnu->NumDaughters();
  lar_pandora::MetadataVector neutrino_metadata_vec = particlesToMetadata.at(pfnu);
  lar_pandora::VertexVector neutrino_vertex_vec = particlesToVertices.at(pfnu);
  if (neutrino_metadata_vec.size() != 1 || neutrino_vertex_vec.size() != 1)
  {
    std::cout << "[NuCC::FillReconstructed] Neutrino association problem." << std::endl;
  }
  else
  {
    const larpandoraobj::PFParticleMetadata::PropertiesMap &neutrino_properties = neutrino_metadata_vec.front()->GetPropertiesMap();
    fNu_Score = neutrino_properties.at("NuScore");
    fNu_SliceIndex = neutrino_properties.at("SliceIndex");
    const recob::Vertex::Point_t &neutrino_vtx = neutrino_vertex_vec.front()->position();
    fNu_Vx = neutrino_vtx.X();
    fNu_Vy = neutrino_vtx.Y();
    fNu_Vz = neutrino_vtx.Z();
    if (m_hasMCNeutrino)
    {
      fTrueNu_VtxDistance = pandoraInterfaceHelper.Distance3D(fTrueNu_VxSce, fTrueNu_VySce, fTrueNu_VzSce, fNu_Vx, fNu_Vy, fNu_Vz);
    }
  }
  pandoraInterfaceHelper.CollectDownstreamPFParticles(particleMap, pfnu, pfdaughters);
  fNumDaughters = pfdaughters.size() - 1; // The neutrino itself is included here.
  std::cout << "[NuCC::FillReconstructed] neutrino PDG: " << fNu_PDG << ", Primary Daughters: " << fNumPrimaryDaughters;
  std::cout << ", Daughters: " << fNumDaughters << ", TopoScore: " << fNu_Score;
  std::cout << ", TrueNu_VtxDistance: " << fTrueNu_VtxDistance << std::endl;

  for (auto const pfp : pfdaughters)
  {
    if (!pfp->IsPrimary())
    {
      if (!FillDaughters(pfp, MCSMu_handle, trackPIDAssn))
      {
        fDaughtersStored = false;
      }
      else
      {
        if (MatchDaughter(evt, pfp))
          fNumMatchedDaughters++;
        fNueDaughtersTree->Fill();
      }
    }
  }
}

bool NuCC::FillDaughters(const art::Ptr<recob::PFParticle> &pfp,
                         const art::ValidHandle<std::vector<recob::MCSFitResult>> &MCSMu_handle,
                         const art::FindManyP<anab::ParticleID> &trackPIDAssn)
{
  clearDaughter();
  const lar_pandora::ClusterVector cluster_vec = particlesToClusters.at(pfp);
  std::vector<uint> nHits;
  std::vector<float> pfenergy;
  energyHelper.energy_from_hits(cluster_vec, nHits, pfenergy);
  fNhitsU = nHits[0];
  fNhitsV = nHits[1];
  fNhitsY = nHits[2];
  fCaloU = pfenergy[0];
  fCaloV = pfenergy[1];
  fCaloY = pfenergy[2];
  fNu_NhitsU += fNhitsU;
  fNu_NhitsV += fNhitsV;
  fNu_NhitsY += fNhitsY;
  fNu_CaloU += fCaloU;
  fNu_CaloV += fCaloV;
  fNu_CaloY += fCaloY;

  if (particlesToSpacePoints.find(pfp) == particlesToSpacePoints.end())
  {
    // If a daughter has no associated spacepoints, count the hits to contribute to the total, but dont save the daughter
    std::cout << "[NuCC::FillDaughters] Daughter had no associated spacepoints." << std::endl;
    return false;
  }
  fNSpacepoints = particlesToSpacePoints.at(pfp).size();
  fNu_NSpacepoints += fNSpacepoints;

  if (particlesToVertices.find(pfp) == particlesToVertices.end())
  {
    // If a daughter has no associated vertex, count the hits to contribute to the total, but dont save the daughter
    std::cout << "[NuCC::FillDaughters] Daughter had no associated vertex." << std::endl;
    return false;
  }

  if (particlesToMetadata.at(pfp).size() != 1 || particlesToVertices.at(pfp).size() != 1)
  {
    std::cout << "[NuCC::FillDaughters] Daughter association problem." << std::endl;
    return false;
  }
  const recob::Vertex::Point_t &pfp_vtx = particlesToVertices.at(pfp).front()->position();
  fVx = pfp_vtx.X();
  fVy = pfp_vtx.Y();
  fVz = pfp_vtx.Z();
  const larpandoraobj::PFParticleMetadata::PropertiesMap &pfp_properties = particlesToMetadata.at(pfp).front()->GetPropertiesMap();
  fTrackScore = pfp_properties.at("TrackScore");
  fVtxDistance = pandoraInterfaceHelper.Distance3D(fVx, fVy, fVz, fNu_Vx, fNu_Vy, fNu_Vz);

  // Track-like fields
  if (particlesToTracks.find(pfp) != particlesToTracks.end())
  {
    fIsTrack = true;
    fNumTracks++;
    const art::Ptr<recob::Track> this_track = particlesToTracks.at(pfp).front();
    const recob::Track::Vector_t &track_dir = this_track->StartDirection();
    fTrackLength = this_track->Length();
    fTrackDirX = track_dir.X();
    fTrackDirY = track_dir.Y();
    fTrackDirZ = track_dir.Z();
    fTrackEndX = this_track->End().X();
    fTrackEndY = this_track->End().Y();
    fTrackEndZ = this_track->End().Z();

    // MCS momentum:
    const recob::MCSFitResult &mcsMu = MCSMu_handle->at(this_track.key());
    fTrackMCS_mom = mcsMu.fwdMomentum();
    fTrackMCS_err = mcsMu.fwdMomUncertainty();
    fTrackMCS_ll = mcsMu.fwdLogLikelihood();

    // PID information:
    std::map<std::string, float> pid_map;
    if(trackHelper.getPID(pid_map, this_track, trackPIDAssn))
    {
      fTrackPID_chiproton = pid_map.at("chi2_proton");
      fTrackPID_chimuon = pid_map.at("chi2_muon");
      std::cout << "[NuCC::FillDaughters] fTrackPID_chiproton: " << fTrackPID_chiproton << ", fTrackPID_chimuon: " << fTrackPID_chimuon << std::endl;
    }
    else
    {
      std::cout << "[NuCC::FillDaughters] Track has no PID attached to it" << std::endl;
    }
  }

  // Shower-like fields
  if (particlesToShowers.find(pfp) != particlesToShowers.end())
  {
    fIsShower = true;
    fNumShowers++;
    const art::Ptr<recob::Shower> this_shower = particlesToShowers.at(pfp).front();
    if (this_shower->has_length() && this_shower->has_open_angle())
    {
      const TVector3 &shower_dir = this_shower->Direction();
      fShowerLength = this_shower->Length();
      fShowerOpenAngle = this_shower->OpenAngle();
      fShowerDirX = shower_dir.X();
      fShowerDirY = shower_dir.Y();
      fShowerDirZ = shower_dir.Z();

      std::vector<float> pitches(3, std::numeric_limits<float>::lowest());
      std::vector<float> dqdx(3, std::numeric_limits<float>::lowest());
      std::vector<std::vector<float>> dqdx_hits(3, std::vector<float>());
      energyHelper.dQdx(shower_dir, cluster_vec, clustersToHits, dqdx, dqdx_hits, pitches);
      std::vector<float> dedx = energyHelper.dEdx_from_dQdx(dqdx);

      fDedxU = dedx[0];
      fDedxV = dedx[1];
      fDedxY = dedx[2];
      fDedxHitsU = dqdx_hits[0].size();
      fDedxHitsV = dqdx_hits[1].size();
      fDedxHitsY = dqdx_hits[2].size();
      fDedxPitchU = pitches[0];
      fDedxPitchV = pitches[1];
      fDedxPitchY = pitches[2];
    }
    else
    {
      std::cout << "[NuCC::FillDaughters] Bad shower, no length or opening angle!" << std::endl;
    }
  }

  // Hierarchy info
  fGeneration = larpandora.GetGeneration(particleMap, pfp);
  if (fNumPrimaryDaughters < fNumDaughters)
  {
    if (particleMap.at(pfp->Parent())->PdgCode() == 13)
    {
      fIsTrackDaughter = true;
    }
    if (pfp->NumDaughters())
    {
      for (const int daughter_id : pfp->Daughters())
      {
        if (particleMap.at(daughter_id)->PdgCode() == 11)
        {
          fHasShowerDaughter = true;
        }
      }
    }
  }
  std::cout << "[NuCC::FillDaughters] Trackscore: " << fTrackScore << ", Generation: " << fGeneration;
  std::cout << ", vtx distance: " << fVtxDistance << std::endl;
  std::cout << "[NuCC::FillDaughters] U Plane: Hits:" << fNhitsU << ", Energy: " << fCaloU << ", dedx hits: " << fDedxHitsU << ", dedx: " << fDedxU << ", pitch: " << fDedxPitchU << std::endl;
  std::cout << "[NuCC::FillDaughters] V Plane: Hits:" << fNhitsV << ", Energy: " << fCaloV << ", dedx hits: " << fDedxHitsV << ", dedx: " << fDedxV << ", pitch: " << fDedxPitchV << std::endl;
  std::cout << "[NuCC::FillDaughters] Y Plane: Hits:" << fNhitsY << ", Energy: " << fCaloY << ", dedx hits: " << fDedxHitsY << ", dedx: " << fDedxY << ", pitch: " << fDedxPitchY << std::endl;
  return true;
}

bool NuCC::MatchDaughter(art::Event const &evt, const art::Ptr<recob::PFParticle> &pfp)
{
  if (m_isData)
    return false;
  art::Ptr<simb::MCParticle> matched_mcp;
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

  if (m_hasMCNeutrino)
  {
    // Is this MC particle neutrino?
    const art::Ptr<simb::MCTruth> mctruth = pandoraInterfaceHelper.TrackIDToMCTruth(evt, m_geant_producer, matched_mcp->TrackId());
    if (mctruth->Origin() == simb::kBeamNeutrino)
    {
      fMatchedNeutrino = true;
    }
    else
    {
      fMatchedNeutrino = false;
      fCosmicMatched = true;
    }

    fTruePDG = matched_mcp->PdgCode();
    fTrueEnergy = matched_mcp->E();
    fTrueVx = matched_mcp->Vx();
    fTrueVy = matched_mcp->Vy();
    fTrueVz = matched_mcp->Vz();
    fTrueLength = (matched_mcp->Position().Vect() - matched_mcp->EndPosition().Vect()).Mag();

    pandoraInterfaceHelper.SCE(fTrueVx, fTrueVy, fTrueVz, matched_mcp->T(),
                               fTrueVxSce, fTrueVySce, fTrueVzSce);
    std::cout << "[NuCC::MatchDaughter] Daughter matched with PDG: " << fTruePDG << ", neutrino origin: " << fMatchedNeutrino << std::endl;
  }
  return true;
}

void NuCC::FillTrueNu(art::Event const &evt)
{
  if (m_hasMCNeutrino)
  {
    auto const &generator_handle = evt.getValidHandle<std::vector<simb::MCTruth>>("generator");
    auto const &generator(*generator_handle);
    fNumNu = generator.size();
    std::cout << "[NuCC::FillTrueNu] True neutrinos found: " << fNumNu;
    if (generator.size() > 0)
    {
      if (generator.front().Origin() != simb::kBeamNeutrino)
      {
        std::cout << "[NuCC::FillTrueNu] Origin of generator particle is not kBeamNeutrino." << std::endl;
        return;
      }
      const simb::MCNeutrino &mcnu = generator.front().GetNeutrino();

      fTrueNu_InteractionType = mcnu.InteractionType();
      fTrueNu_CCNC = mcnu.CCNC();
      fTrueNu_PDG = mcnu.Nu().PdgCode();
      fTrueNu_Energy = mcnu.Nu().E();
      fTrueNu_Px = mcnu.Nu().Px();
      fTrueNu_Py = mcnu.Nu().Py();
      fTrueNu_Pz = mcnu.Nu().Pz();
      fTrueNu_LeptonEnergy = mcnu.Lepton().E();
      fTrueNu_LeptonPx = mcnu.Lepton().Px();
      fTrueNu_LeptonPy = mcnu.Lepton().Py();
      fTrueNu_LeptonPz = mcnu.Lepton().Pz();
      fTrueNu_LeptonTheta = mcnu.Theta();
      fTrueNu_Time = mcnu.Nu().T();
      fTrueNu_Vx = mcnu.Nu().Vx();
      fTrueNu_Vy = mcnu.Nu().Vy();
      fTrueNu_Vz = mcnu.Nu().Vz();
      pandoraInterfaceHelper.SCE(fTrueNu_Vx, fTrueNu_Vy, fTrueNu_Vz, fTrueNu_Time,
                                 fTrueNu_VxSce, fTrueNu_VySce, fTrueNu_VzSce);
      std::cout << ", CCNC: " << fTrueNu_CCNC << ", PDG: " << fTrueNu_PDG << ", E: " << fTrueNu_Energy << ", z-vertex: " << fTrueNu_Vz << std::endl;
    }
  }
}

void NuCC::FillTrueNuDaughters(art::Event const &evt)
{
  lar_pandora::MCParticleVector mcparticles;
  larpandora.CollectMCParticles(evt, m_geant_producer, mcparticles);

  for (auto const &mcparticle : mcparticles)
  {
    if (!(mcparticle->Process() == "primary" &&
          mcparticle->T() != 0 &&
          mcparticle->StatusCode() == 1))
      continue;

    const art::Ptr<simb::MCTruth> mc_truth = pandoraInterfaceHelper.TrackIDToMCTruth(evt, m_geant_producer, mcparticle->TrackId());
    if (mc_truth->Origin() == simb::kBeamNeutrino)
    {
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
      std::cout << "[NuCC::FillTrueNuDaughters] << PDG: " << mcparticle->PdgCode() << ", E: " << mcparticle->E() << ", was matched? " << fTrueNu_DaughterMatched.back() << std::endl;
    }
  }
}

void NuCC::FillReconTruthMatching(art::Event const &evt)
{
  pandoraInterfaceHelper.Configure(evt, m_pfp_producer, m_pfp_producer, m_hitfinder_producer, m_geant_producer, m_hit_mcp_producer);
  pandoraInterfaceHelper.GetRecoToTrueMatches(matchedParticles);
  std::cout << "[NuCC::FillReconTruthMatching] ";
  std::cout << "Number of PFPparticles in event: " << pfparticles.size() << std::endl;
  for (auto it = matchedParticles.begin(); it != matchedParticles.end(); ++it)
  {
    matchedMCParticles.insert(it->second);
  }
  std::cout << "[NuCC::FillReconTruthMatching] ";
  std::cout << "PFParticlesToMCParticles constructed: Number of PFPparticles matched: " << matchedParticles.size() << std::endl;
}
