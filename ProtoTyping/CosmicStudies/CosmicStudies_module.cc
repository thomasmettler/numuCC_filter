#include "CosmicStudies.h"

void CosmicStudies::analyze(art::Event const &evt)
{
  clear();
  fRun = evt.run();
  fSubrun = evt.subRun();
  fEvent = evt.id().event();
  fNevents++;

  //// Filling the MC and TPCreco information -----------------------------------------------------
  if (!m_is_data)
  {
    if (!m_slimmed)
    {
      pandoraHelper.Configure(evt, m_pfp_producer, m_pfp_producer, m_hitfinder_producer, m_geant_producer);
      std::cout << "[CosmicStudies] Reco-Truth matcher configured (full file)." << std::endl;
    }
    else
    {
      pandoraHelper.Configure(evt, m_pfp_producer, m_pfp_producer, m_hitfinder_producer, m_geant_producer, m_hit_mcp_producer);
      std::cout << "[CosmicStudies] Reco-Truth matcher configured (slimmed file)." << std::endl;
    }
    // Get the map PFP->MCP and the set of MCPs
    pandoraHelper.GetRecoToTrueMatches(matchedParticles);
    std::cout << "[CosmicStudies] Reco-Truth matches constructed." << std::endl;
    for (auto it = matchedParticles.begin(); it != matchedParticles.end(); ++it)
    {
      matchedMCParticles.insert(it->second);
    }
    fill_MC(evt);
  }
  fill_TPCreco(evt);
  ////---------------------------------------------------------------------------------

  //// Filling the flashes ------------------------------------------------------------
  art::ValidHandle<std::vector<recob::OpFlash>> const &simple_beam_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(m_beam_simpleflash_producer);
  fNumSimpleBeamFlashes = simple_beam_handle->size();
  std::cout << "[CosmicStudies] fNumSimpleBeamFlashes: " << fNumSimpleBeamFlashes << std::endl;
  fill_flash(simple_beam_handle, fNumSimpleBeamFlashes, fSimpleBeamFlashesTree);

  art::ValidHandle<std::vector<recob::OpFlash>> const &op_beam_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(m_beam_opflash_producer);
  fNumOpBeamFlashes = op_beam_handle->size();
  std::cout << "[CosmicStudies] fNumOpBeamFlashes: " << fNumOpBeamFlashes << std::endl;
  fill_flash(op_beam_handle, fNumOpBeamFlashes, fOpBeamFlashesTree);

  art::ValidHandle<std::vector<recob::OpFlash>> const &simple_cosmic_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(m_cosmic_simpleflash_producer);
  fNumSimpleCosmicFlashes = simple_cosmic_handle->size();
  std::cout << "[CosmicStudies] fNumSimpleCosmicFlashes: " << fNumSimpleCosmicFlashes << std::endl;
  fill_flash(simple_cosmic_handle, fNumSimpleCosmicFlashes, fSimpleCosmicFlashesTree);

  art::ValidHandle<std::vector<recob::OpFlash>> const &op_cosmic_handle = evt.getValidHandle<std::vector<recob::OpFlash>>(m_cosmic_opflash_producer);
  fNumOpCosmicFlashes = op_cosmic_handle->size();
  std::cout << "[CosmicStudies] fNumOpCosmicFlashes: " << fNumOpCosmicFlashes << std::endl;
  fill_flash(op_cosmic_handle, fNumOpCosmicFlashes, fOpCosmicFlashesTree);
  ////---------------------------------------------------------------------------------

  fEventTree->Fill();
}

void CosmicStudies::endSubRun(const art::SubRun &sr)
{
  fRun_sr = sr.run();
  fSubrun_sr = sr.subRun();

  art::Handle<sumdata::POTSummary> potListHandle;
  if (!m_is_data)
  {
    if (sr.getByLabel("generator", potListHandle))
    {
      fPot = potListHandle->totpot;
      std::cout << "[CosmicStudies::endSubRun] POT for SubRun: " << fPot << std::endl;
    }
    else
      fPot = 0.;
  }
  else
  {
    if (sr.getByLabel("beamdata", "bnbETOR860", potListHandle))
      fPot = potListHandle->totpot;
    else
      fPot = 0.;
  }

  fPOTTree->Fill();

  if (m_verb)
  {
    std::set<std::string> string_process = mcpHelper.getProcesses();
    std::cout << "[CosmicStudies::endSubRun] string_process has members: " << string_process.size() << std::endl;
    for (auto elem : string_process)
    {
      std::cout << elem << ", ";
    }
  }

  // Reset events counter at the end of the subrun
  fNevents = 0;
}

void CosmicStudies::fill_flash(art::ValidHandle<std::vector<recob::OpFlash>> const &flash_handle, uint number, TTree *tree)
{
  float prevTime = 0; // For the first flash, the flashtime will be identical to the flashtimediff
  for (uint ifl = 0; ifl < number; ++ifl)
  {
    clear_Flashes();

    recob::OpFlash const &flash = flash_handle->at(ifl);
    fFlash_TotalPE = flash.TotalPE();
    fFlash_Time = flash.Time();
    fFlash_DiffTime = fFlash_Time - prevTime; //This will not work if the flashes are sorted on PE instead of time!
    prevTime = fFlash_Time;
    fFlash_Y = flash.YCenter();
    fFlash_Z = flash.ZCenter();
    fFlash_sigmaY = flash.YWidth();
    fFlash_sigmaZ = flash.ZWidth();
    fFlash_AbsTime = flash.AbsTime();
    fFlash_Width = flash.TimeWidth();

    for (uint i_pmt = 0; i_pmt < 32; i_pmt++)
    {
      if (flash.PE(i_pmt) > (fFlash_TotalPE / 10.0))
      {
        fFlash_num10percentPMT++;
      }
    }
    tree->Fill();
  }
}

void CosmicStudies::fill_MC(art::Event const &evt)
{
  if (m_is_true_nu)
  {
    auto const &generator_handle = evt.getValidHandle<std::vector<simb::MCTruth>>("generator");
    auto const &generator(*generator_handle);
    if( generator.size() == 1){
      fNum_nu = generator[0].NeutrinoSet(); 
    }
    else{
      fNum_nu = generator.size();
    }

    if (m_verb)
    {
      std::cout << "[CosmicStudies] True neutrinos found: " << fNum_nu << std::endl;
    }

    for (auto &gen : generator)
    {

      if (gen.Origin() == simb::kBeamNeutrino)
      {
        fNu_pdg_code.push_back(gen.GetNeutrino().Nu().PdgCode());
        fNu_time.push_back(gen.GetNeutrino().Nu().T());
        fNu_E.push_back(gen.GetNeutrino().Nu().E());
        fNu_ccnc.push_back(gen.GetNeutrino().CCNC());
        fNu_vtx_x.push_back(gen.GetNeutrino().Nu().Vx());
        fNu_vtx_y.push_back(gen.GetNeutrino().Nu().Vy());
        fNu_vtx_z.push_back(gen.GetNeutrino().Nu().Vz());
      }
    }
  }

  auto const &mcparticles_handle = evt.getValidHandle<std::vector<simb::MCParticle>>(m_geant_producer);
  std::vector<art::Ptr<simb::MCParticle>> mcp_v;
  art::fill_ptr_vector(mcp_v, mcparticles_handle);

  fNumMcp = mcparticles_handle->size();

  for (size_t i_mcp = 0; i_mcp < fNumMcp; i_mcp++)
  {
    simb::MCParticle const &mcparticle = mcparticles_handle->at(i_mcp);
    // Important, only save MC particles with energy over 100MeV, THIS WILL SAVE ALL MUONS.
    uint pdg = abs(mcparticle.PdgCode());
    bool pdg_ok = (pdg == 11 or pdg == 13 or pdg == 211 or pdg == 111 or pdg == 22 or pdg == 2112 or pdg == 2212);
    if (mcparticle.E() > constants::MCP_E_CUT && pdg_ok)
    {
      fMc_Time = mcparticle.T();
      fMc_StatusCode = mcparticle.StatusCode();
      fMc_E = mcparticle.E();
      fMc_PdgCode = mcparticle.PdgCode();
      fMc_StartMomentumX = mcparticle.Px();
      fMc_StartMomentumY = mcparticle.Py();
      fMc_StartMomentumZ = mcparticle.Pz();
      fMc_StartX = mcparticle.Vx();
      fMc_StartY = mcparticle.Vy();
      fMc_StartZ = mcparticle.Vz();
      fMc_EndX = mcparticle.EndX();
      fMc_EndY = mcparticle.EndY();
      fMc_EndZ = mcparticle.EndZ();

      if (m_is_true_nu)
      {
        // Is this MC particle neutrino?
        const art::Ptr<simb::MCTruth> mctruth = pandoraHelper.TrackIDToMCTruth(evt, m_geant_producer, mcparticle.TrackId());
        if (mctruth->Origin() == simb::kBeamNeutrino)
        {
          fMc_kBeamNeutrino = true;
        }
        else
        {
          fMc_kBeamNeutrino = false;
        }
      }

      MCParticleInfo this_mcp = mcpHelper.fillMCP(mcparticle);
      fMc_Process = this_mcp.process;
      fMc_EndProcess = this_mcp.end_process;
      fMc_StartInside = this_mcp.startInside;
      fMc_EndInside = this_mcp.endInside;
      fMc_PartInside = this_mcp.partInside;
      fMc_StartX_tpc = this_mcp.startX_tpc;
      fMc_StartY_tpc = this_mcp.startY_tpc;
      fMc_StartZ_tpc = this_mcp.startZ_tpc;
      fMc_EndX_tpc = this_mcp.endX_tpc;
      fMc_EndY_tpc = this_mcp.endY_tpc;
      fMc_EndZ_tpc = this_mcp.endZ_tpc;
      fMc_Length = this_mcp.length;
      fMc_LengthTPC = this_mcp.lengthTPC;

      // is this mcp matched to a pfp? Only check if primary!
      if (fMc_Process == 23)
      {
        fMc_Matched = matchedMCParticles.find(mcp_v.at(i_mcp)) != matchedMCParticles.end();
      }

      // Check if we have a CRT crossing:
      // Require a charged particle, a minimum energy and an end point below the CRT in y:
      bool pdg_crt_ok = pdg == 11 or pdg == 13 or pdg == 211 or pdg == 2212 or pdg == 2112;
      if (pdg_crt_ok && fMc_E > constants::CRT_E_CUT && fMc_EndY < constants::BY)
      {
        CRTcrossing this_crossing = mcpHelper.isCrossing(mcparticle);
        fCRT_crossed = this_crossing.CRT_cross;
        if (fCRT_crossed)
        {
          fNumCross++;
          fCrossE = this_crossing.crossE;
          fCrossT = this_crossing.crossT;
          fCrossX = this_crossing.crossX;
          fCrossY = this_crossing.crossY;
          fCrossZ = this_crossing.crossZ;
          fCRTcrossTree->Fill();
        }
      }
      fMCParticlesTree->Fill();
      fNumMcp_saved++; //Counter
    }
  }
  std::cout << "[CosmicStudies::fill_MC] ";
  std::cout << "Number of MC CRT crossing particles in event: " << fNumCross << std::endl;
}

void CosmicStudies::fill_TPCreco(art::Event const &evt)
{
  auto const &pfparticle_handle = evt.getValidHandle<std::vector<recob::PFParticle>>(m_pfp_producer);
  std::vector<art::Ptr<recob::PFParticle>> pfp_v;
  art::fill_ptr_vector(pfp_v, pfparticle_handle);

  auto const &cluster_handle = evt.getValidHandle<std::vector<recob::Cluster>>(m_pfp_producer);
  auto const &spacepoint_handle = evt.getValidHandle<std::vector<recob::SpacePoint>>(m_pfp_producer);
  auto const &MCSMu_handle = evt.getValidHandle<std::vector<recob::MCSFitResult>>("pandoraCosmicMCSMu");

  art::FindOneP<recob::Vertex> vertex_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
  art::FindManyP<recob::Cluster> clusters_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
  art::FindManyP<recob::Hit> hits_per_cluster(cluster_handle, evt, m_pfp_producer);
  art::FindManyP<recob::SpacePoint> spcpnts_per_pfpart(pfparticle_handle, evt, m_pfp_producer);
  art::FindManyP<recob::Hit> hits_per_spcpnts(spacepoint_handle, evt, m_pfp_producer);
  art::FindOneP<recob::Track> track_per_pfpart(pfparticle_handle, evt, m_pfp_producer);

  if (m_verb && !m_is_data)
  {
    std::cout << "[CosmicStudies::fill_TPCreco] ";
    std::cout << "PFParticlesToMCParticles constructed: Number of PFPparticles matched: " << matchedParticles.size() << std::endl;
    std::cout << "[CosmicStudies::fill_TPCreco] ";
    std::cout << "Number of PFPparticles in event: " << pfp_v.size() << std::endl;
  }

  /// PFParticle loop
  fNumPfp = pfp_v.size();
  for (uint i_pfp = 0; i_pfp < fNumPfp; i_pfp++)
  {
    clear_PFParticle();

    recob::PFParticle const &pfparticle = pfparticle_handle->at(i_pfp);
    fPdgCode = pfparticle.PdgCode();
    fNumDaughters = pfparticle.NumDaughters();
    fIsPrimary = pfparticle.IsPrimary();

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
      clear_Cluster();

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
        uint plane = hit->WireID().Plane;
        if (plane == 0)
        {
          fNhitsU += 1;
        }
        else if (plane == 1)
        {
          fNhitsV += 1;
        }
        else if (plane == 2)
        {
          fNhitsY += 1;
        }
      }
    }

    // The Pfparticle is a TRACK
    if (fPdgCode == 13)
    {
      art::Ptr<recob::Track> const &track_obj = track_per_pfpart.at(i_pfp);
      if (track_obj.isNull())
      {
        std::cout << "[CosmicStudies::fill_TPCreco] track is Null" << std::endl;
        fTrack_Valid = false;
      }
      else
      {
        fTrack_Valid = true;
        fTrack_HasMomentum = track_obj->HasMomentum();

        fTrack_StartX = track_obj->Start().X();
        fTrack_StartY = track_obj->Start().Y();
        fTrack_StartZ = track_obj->Start().Z();
        fTrack_EndX = track_obj->Trajectory().End().X();
        fTrack_EndY = track_obj->Trajectory().End().Y();
        fTrack_EndZ = track_obj->Trajectory().End().Z();
        fTrack_Length = track_obj->Length();
        fTrack_StartMomentumX = track_obj->StartMomentumVector().X();
        fTrack_StartMomentumY = track_obj->StartMomentumVector().Y();
        fTrack_StartMomentumZ = track_obj->StartMomentumVector().Z();
        fTrack_EndMomentumX = track_obj->EndMomentumVector().X();
        fTrack_EndMomentumY = track_obj->EndMomentumVector().Y();
        fTrack_EndMomentumZ = track_obj->EndMomentumVector().Z();
        fTrack_Theta = track_obj->Theta();
        fTrack_Phi = track_obj->Phi();
        fTrack_ZenithAngle = track_obj->ZenithAngle();
        fTrack_AzimuthAngle = track_obj->AzimuthAngle();

        // MCS momentum:
        const recob::MCSFitResult &mcsMu = MCSMu_handle->at(track_obj.key());
        fTrack_MCS_mom = mcsMu.fwdMomentum();
        fTrack_MCS_err = mcsMu.fwdMomUncertainty();
        fTrack_MCS_ll = mcsMu.fwdLogLikelihood();
        // Muon energy hypothesis
        fTrack_MCS_energy = (std::sqrt(std::pow(fTrack_MCS_mom * 1000, 2) + std::pow(constants::MUON_M_MEV, 2)) - constants::MUON_M_MEV) / 1000.;
      }
    }

    try
    {
      art::Ptr<recob::Vertex> vertex_obj = vertex_per_pfpart.at(i_pfp);
      double vertex[3];
      vertex_obj->XYZ(vertex);
      fVx = vertex[0];
      fVy = vertex[1];
      fVz = vertex[2];
    }
    catch (...)
    {
      std::cout << "[CosmicStudies::fill_TPCreco] No vertex found for " << fPdgCode << " with " << fNhits << std::endl;
    }

    // Fill the reco truth matched fields
    if (!m_is_data)
    {
      if (matchedParticles.find(pfp_v.at(i_pfp)) != matchedParticles.end())
      {
        art::Ptr<simb::MCParticle> matched_mcp = matchedParticles[pfp_v.at(i_pfp)];

        if (m_is_true_nu)
        {
          // Is this MC particle neutrino?
          const art::Ptr<simb::MCTruth> mctruth = pandoraHelper.TrackIDToMCTruth(evt, m_geant_producer, matched_mcp->TrackId());
          if (mctruth->Origin() == simb::kBeamNeutrino)
          {
            fTrack_matched_kBeamNeutrino = true;
          }
          else
          {
            fTrack_matched_kBeamNeutrino = false;
          }
        }

        fTrack_matched_PdgCode = matched_mcp->PdgCode();
        fTrack_matchedE = matched_mcp->E();
        fTrack_matched_Time = matched_mcp->T();
        fTrack_matched_StartMomentumX = matched_mcp->Px();
        fTrack_matched_StartMomentumY = matched_mcp->Py();
        fTrack_matched_StartMomentumZ = matched_mcp->Pz();

        MCParticleInfo this_mcp = mcpHelper.fillMCP(*matched_mcp);
        fTrack_matched_Process = this_mcp.process;
        fTrack_matched_EndProcess = this_mcp.end_process;
        fTrack_matched_StartInside = this_mcp.startInside;
        fTrack_matched_EndInside = this_mcp.endInside;
        fTrack_matched_PartInside = this_mcp.partInside;

        fTrack_matched_StartX = this_mcp.startX_tpc;
        fTrack_matched_StartY = this_mcp.startY_tpc;
        fTrack_matched_StartZ = this_mcp.startZ_tpc;
        fTrack_matched_EndX = this_mcp.endX_tpc;
        fTrack_matched_EndY = this_mcp.endY_tpc;
        fTrack_matched_EndZ = this_mcp.endZ_tpc;

        fTrack_matched_StartX_sce = this_mcp.startX_sce;
        fTrack_matched_StartY_sce = this_mcp.startY_sce;
        fTrack_matched_StartZ_sce = this_mcp.startZ_sce;
        fTrack_matched_EndX_sce = this_mcp.endX_sce;
        fTrack_matched_EndY_sce = this_mcp.endY_sce;
        fTrack_matched_EndZ_sce = this_mcp.endZ_sce;
        fTrack_matched_LengthTPC = this_mcp.lengthTPC;
        fTrack_matched_Length_sce = this_mcp.length_sce;
      }
    }
    // Save only PF tracks with a length over 5cm
    if (fTrack_Valid && fTrack_Length > constants::PFP_LENGTH_CUT)
    {
      fPFParticlesTree->Fill();
      fNumPfp_saved++;
    }
  }
}
