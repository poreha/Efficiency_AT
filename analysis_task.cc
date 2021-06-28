//
// Created by mikhail on 6/16/20.
// Edited by alexey on 5/21/21
//

#include "analysis_task.h"
#include <TCanvas.h>
#include <math.h>

namespace AnalysisTree
{
  void AnalysisTask::Init(std::map<std::string, void *> &branch_map)
  {
    // linking pointers with branch fields
    event_header_ = static_cast<EventHeader *>(branch_map.at("event_header"));
    mdc_vtx_tracks_ = static_cast<Particles *>(branch_map.at("mdc_vtx_tracks"));
    meta_hits_ = static_cast<HitDetector *>(branch_map.at("meta_hits"));
    mdc_meta_matching_ = static_cast<Matching *>(branch_map.at("mdc_vtx_tracks2meta_hits"));

    // for simulated
    sim_header_ = static_cast<EventHeader *>(branch_map.at("sim_header"));
    sim_tracks_ = static_cast<Particles *>(branch_map.at("sim_tracks"));
    track_matching_ = static_cast<Matching *>(branch_map.at("mdc_vtx_tracks2sim_tracks"));

    // getting branch configurations, which store information about fields in branches
    auto event_header_config = config_->GetBranchConfig("event_header");
    auto mdc_vtx_tracks_config = config_->GetBranchConfig("mdc_vtx_tracks");
    auto meta_hits_config = config_->GetBranchConfig("meta_hits");

    auto sim_header_config = config_->GetBranchConfig("sim_header");
    auto sim_tracks_config = config_->GetBranchConfig("sim_tracks");

    // linking necessary for analysis fields with enumerator for fast access to them
    fields_id_.insert(std::make_pair(FIELDS::HITS_TOF, event_header_config.GetFieldId("selected_tof_hits")));
    fields_id_.insert(std::make_pair(FIELDS::BETA, meta_hits_config.GetFieldId("beta")));
    fields_id_.insert(std::make_pair(FIELDS::M2, meta_hits_config.GetFieldId("mass2")));
    fields_id_.insert(std::make_pair(FIELDS::CHARGE, meta_hits_config.GetFieldId("charge")));
    fields_id_.insert(std::make_pair(FIELDS::PT2, event_header_config.GetFieldId("physical_trigger_2")));
    fields_id_.insert(std::make_pair(FIELDS::TOF, meta_hits_config.GetFieldId("time_of_flight")));
    fields_id_.insert(std::make_pair(FIELDS::PATH_LEN, meta_hits_config.GetFieldId("path_length")));

    mc_particle_id_.insert(std::make_pair(uPARTICLES::PROTON, 2212));
    mc_particle_id_.insert(std::make_pair(uPARTICLES::KAONp, 321));
    mc_particle_id_.insert(std::make_pair(uPARTICLES::ELECTRON, 11));

    // initializing histograms
    tof_multiplicity_distribution_ = new TH1F("tof_multiplicity", ";TOF hits;counts", 100, 0, 100);

    // use legend
    PID_RECO_ = new TH2F("PID_reco", "All reconstructed protons;Y;p_{T}, [GeV/c]", 100, 0., 2.5, 100, 0., 2.);
    PID_PRIM_ = new TH2F("PID_prim", "Reconstructed primary protons matched with MC-primary particles;Y;p_{T}, [GeV/c]", 100, 0., 2.5, 100, 0., 2.);
    PID_SEC_ = new TH2F("PID_sec", "Reconstructed secondary protons matched with MC-secondary particles;Y;p_{T}, [GeV/c]", 100, 0., 2.5, 100, 0., 2.);

    PDG_PRIM_RECO_ = new TH2F("PDG_prim_RECO", "Particles matched with MC-primary proton;Y;p_{T}, [GeV/c]", 100, 0., 2.5, 100, 0., 2.5);
    PDG_SEC_RECO_ = new TH2F("PDG_sec_RECO", "Particles matched with MC-secondary proton;Y;p_{T}, [GeV/c]", 250, 0., 2.5, 250, 0., 2.);

    MISMATCH_ = new TH2F("Mismatch", "Protons matched with all MC-particles except protons;Y;p_{T}, GeV/c", 100, 0., 2.5, 100, 0., 2.);
    LOSS_ = new TH2F("Loss", "All MC-particles, except protons, matched with MC-protons;Y;p_{T}, GeV/c", 100, 0, 2.5, 100, 0, 2.);

    PDG_ = new TH2F("PDG", "All MC-protons;Y;p_{T}, GeV/c", 100, 0., 2.5, 100, 0., 2.);
    PDG_PRIM_ = new TH2F("PDG_prim", "All primary MC-protons;Y;p_{T}, GeV/c", 100, 0., 2.5, 100, 0., 2.);
    PDG_SEC_ = new TH2F("PDG_sec", "All secondary MC-protons;Y;p_{T}, GeV/c", 250, 0, 2.5, 250, 0, 2.);

    PHIvPT_ = new TH2F("Phi_vs_pT", ";Phi;p_{T}, GeV", 256, -3.5, 3.5, 250, 0, 5);
    MOMENTUMvBETA_ = new TH2F("Momentum_vs_beta", ";p/q, GeV/c; beta", 256, -3, 6.2, 256, 0, 1.3);
    M2vMOMENTUM_ = new TH2F("Momentum_vs_Mass2", ";p/q, GeV/c;mass^2, GeV", 256, -6, 6, 512, -10, 60); // check it
    RAPIDITYvPHI_ = new TH2F("Rapidity_vs_Phi", ";Phi ;Rapidity", 256, -3.5, 3.5, 256, 0, 2);

    Efficiency_Acceptance_ = new TH2F("Eff", "Efficiency (x) Acceptance; Y; p_{T}", 100, 0., 2.5, 100, 0, 2.);
    PID_SIM_ = new TH2F("PID_SIM", "PID for simulated momentum; Y, p_{T}", 100, 0., 2.5, 100, 0., 2.4);

    DCA = new TH2F("dca", "DCA xy vs z; dca_xy, mm; dca_z, mm", 100, -100, 100, 100, -100, 150);
    VTX_XZ = new TH2F("vtx_xz", "Vertex x vs z; vtx_x, mm; vtx_z, mm", 100, -30, 30, 100, -100, 100);
    
    DCA_XY = new TH1F("dca_xy", "DCA xy; dca_xy", 1024, -20, 20);
    DCA_Z = new TH1F("dca_z", "DCA z; dca_z", 1024, -50, 50);
    VTX_X = new TH1F("vtx_x", "Vertex x; vtx, mm", 1024, -20, 20);
    VTX_Y = new TH1F("vtx_y", "Vertex y; vtx, mm", 1024, -20, 20);
    VTX_Z = new TH1F("vtx_z", "Vertex z; vtx, mm", 1024, -50, 50);

    RESOLUTION_X = new TH1F("resolution_x", "Resolution; mm", 1024, -20, 20);
  }

  void AnalysisTask::Exec()
  {
    PdgCode_t proton = 2212; //Particle Data Group code == int // not needed
    //auto hits_tof = event_header_->GetField<int>(fields_id_.at(FIELDS::HITS_TOF)); // getting multiplicity from event header
    //tof_multiplicity_distribution_->Fill(hits_tof); // filling histogram

//auto im_param = sim_header_->GetField<float>(1); // impact parameter // посмотреть по событиям, а не по трекам

    int n_tracks = mdc_vtx_tracks_->GetNumberOfChannels(); // number of tracks in current event
    for (int i = 0; i < n_tracks; ++i)
    {                                              // loop over all tracks if current event
      auto track = mdc_vtx_tracks_->GetChannel(i); // getting track from track detector
      PdgCode_t pid = track.GetPid();

      Int_t matched_track = track_matching_->GetMatchDirect(i); // matching simulated track's ID to current reconstructed track

      //  для no_match(id = -999) ввести обработку
      if (matched_track == -999)
      {
        continue;
      }
      auto sim_track = sim_tracks_->GetChannel(matched_track);    // getting matched simulated track // changed i->matched_track
      PdgCode_t pdg = sim_track.GetPid();                           // getting PID of matched simulated track
      Integer_t match_meta_hit = mdc_meta_matching_->GetMatchDirect(i); // getting index of matched with track TOF-system hit

      auto hit = meta_hits_->GetChannel(match_meta_hit); // getting matched with track hit in TOF-system

      auto pT = track.GetPt();         // getting transverse momentum
      auto sim_pT = sim_track.GetPt(); // getting transverse momentum of matched simulated track

      auto eta = track.GetEta();                   // getting pseudorapidity
      auto rapidity = track.GetRapidity();         // getting rapidity
      auto sim_rapidity = sim_track.GetRapidity(); // getting rapidity of matched simulated track

      auto sim_p = sim_track.GetP();
      auto sim_m2 = sim_track.GetMass();

      auto charge = hit.GetField<int>(fields_id_.at(FIELDS::CHARGE)); // getting charge from meta_hits
      auto beta = hit.GetField<float>(fields_id_.at(FIELDS::BETA));   // getting beta from meta_hits
      auto m2 = hit.GetField<float>(fields_id_.at(FIELDS::M2));       // getting mass squarred from meta_hits

      // Filters and Cuts
      bool pt2 = event_header_->GetField<bool>(5);             // physical_trigger_2 (id=5)
      bool kGoodVertexCand = event_header_->GetField<bool>(1); // good_vertex_candidate (id=1)
      bool kGoodTrigger = event_header_->GetField<bool>(13);   // good_trigger (id=13)
      bool kNoVeto = event_header_->GetField<bool>(15);        // no_veto (id=15)
      bool kNoPileUpMeta = event_header_->GetField<bool>(7);

      bool kGoodStart = event_header_->GetField<bool>(2);
      bool kGoodStartMeta = event_header_->GetField<bool>(17);

      float dca_xy = track.GetField<float>(3);
      float dca_z = track.GetField<float>(4);
      float chi2 = track.GetField<float>(0);
      
      // Vertexes
      float reco_vtx_x = sim_track.GetField<float>(0);
      float reco_vtx_y = sim_track.GetField<float>(1);
      float reco_vtx_z = sim_track.GetField<float>(2);
      // filling distributions
      /*
      filter: pie+ 211
              pie0 111
              proton 2212
              K0 311
              K+ 321
              e- 11
    */
      PdgCode_t cert_particle[5] = {211, 111, 2212, 321, 11}; 
      bool is_prim = sim_track.GetField<bool>(0); // sim_particle is a primary

      PdgCode_t particle = proton; // proton
      PdgCode_t mc_pdg = mc_particle_id_.at(uPARTICLES::PROTON);
      // Filling histograms

      // PID
      if (pid == mc_pdg) 
      { 
        PID_RECO_->Fill(rapidity, pT); // all reconstructed certain particles
        if (dca_xy < 15. && dca_xy > -15. && dca_z < 15. && chi2 < 100.)
        {
          if (is_prim)
            {                                
              PID_PRIM_->Fill(rapidity, pT); // any matched certain particle is primary

            }
        }
        else
        {
          if (!is_prim)
          {
            PID_SEC_->Fill(rapidity, pT); // any matched certain particle is secondary

          }
        }
      }

      // MISMATCH if pid == particle , pdg != particle
      else
      {
        MISMATCH_->Fill(sim_rapidity, sim_pT);
        
      }

      // LOSS
      if (pid != particle && pdg == particle)
      {
        LOSS_->Fill(rapidity, pT);
       
      }
    
      // PDG_RECO
      if (pdg == mc_pdg)
      {
        
        if (is_prim)
        {
          PDG_PRIM_RECO_->Fill(sim_rapidity, sim_pT); // matched particle is primary // changed to sim_
        }
        else
        {
          PDG_SEC_RECO_->Fill(sim_rapidity, sim_pT); // matched particle is secondary
          
        }
      }

      RESOLUTION_X->Fill(reco_vtx_x - dca_xy);
      DCA->Fill(dca_xy, dca_z);
      DCA_XY->Fill(dca_xy);
      DCA_Z->Fill(dca_z);
      //PHIvPT_->Fill(phi, pT);
      //MOMENTUMvBETA_->Fill(p, beta); // charge to see negatively charged
      //M2vMOMENTUM_->Fill(p, m2);
      //RAPIDITYvPHI_->Fill(phi, rapidity);
      //mass2_distribution_branch_->Fill(m2);
    } //end of reco-cycle


    // GEN -> PDG
    int sim_ntracks = sim_tracks_->GetNumberOfChannels();
    for (int sim_i = 0; sim_i < sim_ntracks; sim_i++)
    {
      auto sim_track = sim_tracks_->GetChannel(sim_i); // getting a simulated track 
      int pdg = sim_track.GetPid();                // getting PID of the simulated track   

      auto vtx_x = sim_track.GetField<float>(0);
      auto vtx_y = sim_track.GetField<float>(1);
      auto vtx_z = sim_track.GetField<float>(2);
      

      auto sim_pT = sim_track.GetPt();             // getting transverse momentum of the simulated track
      auto sim_rapidity = sim_track.GetRapidity(); // getting rapidity of the simulated track

      auto sim_p = sim_track.GetP();
      auto sim_m2 = sim_track.GetMass();

      bool is_prim = sim_track.GetField<bool>(0); // sim_particle is a primary

      if (pdg == 2212)
      {

        PDG_->Fill(sim_rapidity, sim_pT);
        VTX_X->Fill(vtx_x);
        VTX_Y->Fill(vtx_y);
        VTX_Z->Fill(vtx_z);
        VTX_XZ->Fill(vtx_x, vtx_z);
        if (is_prim)
        {
          PDG_PRIM_->Fill(sim_rapidity, sim_pT);
          
        }
        else
        {
          PDG_SEC_->Fill(sim_rapidity, sim_pT);
          
        }
      }
    }//end of mc-cycle

        Efficiency_Acceptance_->Divide(PID_PRIM_,PDG_PRIM_);
  }

  void AnalysisTask::Finish()
  {
    // Writing histograms to file
    //selected_tof_rpc_centrality_->Write();
    DCA->SetMinimum(10);
    VTX_XZ->SetMinimum(10);
    DCA->Write();
    VTX_XZ->Write();
    DCA_XY->Write();
    DCA_Z->Write();
    VTX_X->Write();
    VTX_Y->Write();
    VTX_Z->Write();
    RESOLUTION_X->Write();
    Efficiency_Acceptance_->Write();

    PID_RECO_->Write();
    PID_PRIM_->Write();
    PID_SEC_->Write();
    PDG_PRIM_RECO_->Write();
    PDG_SEC_RECO_->Write();
    MISMATCH_->Write();
    //LOSS_->Write();
    //PDG_->Write();
    PDG_PRIM_->Write();
    //PDG_SEC_->Write();

    //PHIvPT_->Write();
    //MOMENTUMvBETA_->Write();
    M2vMOMENTUM_->Write();
    //RAPIDITYvPHI_->Write();
    
  }
} // namespace AnalysisTree

/*
    TO DO
    make a Factory to create multiple histos
*/