/**
 *  @file   larpandora/LArPandoraEventBuilding/CRTTriggerProducer_module.cc
 *
 *  @brief  module for lar pandora external event building
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/make_tool.h"

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "larpandora/LArPandoraEventBuilding/Slice.h"
#include "larpandora/LArPandoraEventBuilding/LArPandoraEvent.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

//#include "larpandora/LArPandoraObjects/PFParticleMetadata.h"

#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/CRT/CRTTrack.hh"
#include "ubobj/CRT/CRTTzero.hh"
#include "ubcrt/CRT/CRTAuxFunctions.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"

//#include "ubcrt/CRTXSEC/CRTAnaFun.hh"

#include "lardataobj/RecoBase/Track.h"                                                                
#include "lardataobj/RecoBase/Hit.h"     
#include "lardataobj/AnalysisBase/T0.h" 
#include "lardataobj/AnalysisBase/ParticleID.h" 
#include "lardataobj/AnalysisBase/CosmicTag.h"                                                        
#include "lardataobj/AnalysisBase/Calorimetry.h"                                                                 
#include "lardataobj/RecoBase/OpFlash.h"                                                              
#include "lardata/Utilities/AssociationUtil.h"
#include "lardataobj/RecoBase/PFParticle.h"

#include "Pandora/PdgTable.h"

#include "TTree.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TClonesArray.h"
#include "TObject.h"


namespace lar_pandora
{

class CRTTriggerProducer : public art::EDProducer
{
public:
    explicit CRTTriggerProducer(fhicl::ParameterSet const & pset);
    
    CRTTriggerProducer(CRTTriggerProducer const &) = delete;
    CRTTriggerProducer(CRTTriggerProducer &&) = delete;
    CRTTriggerProducer & operator = (CRTTriggerProducer const &) = delete;
    CRTTriggerProducer & operator = (CRTTriggerProducer &&) = delete;

    void produce(art::Event &evt) override;
    
    // Selected optional functions.
    void beginJob() override;
    void endJob() override;
    
private:
    // data types for pandora functions
    art::ServiceHandle<art::TFileService> tfs;
    TTree * my_event_;
    
    int event_counter = 0;
    uint32_t fEvtNum;                //Number of current event                       
    uint32_t frunNum;                //Run Number taken from event  
    uint32_t fsubRunNum;             //Subrun Number taken from event 
    double fTriTim_sec = 0;          //event trigger time sec
    double fTriTim_nsec = 0;          //event trigger time ns
    UInt_t fTimeHigh = 0;
    UInt_t fTimeLow = 0;
    
    double fAbsTimFla = 0;
    double flash_PE = 0;
    double flash_y = 0;
    double flash_z = 0;
    int flash_inTime = -1;
    
    int nr_crthit = -1; // # crt hits assigned to a tpc track
    double crthit_ts0 = 0;
    double crthit_ts1 = 0;
    int adc_length = 0;
    double crt_adc = 0;
    
    double crtt0_time = -1;
    int crtt0_trig = -1;
    double crtt0_DCA = -1;
    int crtt0_plane = -1;
    
    TTree * my_event_trigger;
    
    int flash_tot = 0;
    int flash_counter = 0;
    int t0_counter = 0;
    int t0_outlayer = 0;
    double crthit_correction = 0;
    double crthit_corr_med = 0;
    
    std::vector<double> flash_crthit_diff;
    
    void initialize_tmyevent();
    void reset_tree();
    void reset_tree_event();
    
    int verbose_;
    int saveTTree_ = 0;
    int run_MC = 0;
    int store_t0_ = 1;
    //double flash_crt_window = 1.0;

    std::string                         m_inputProducerLabel;  ///< Label for the Pandora instance that produced the collections we want to consolidated
    std::string                         data_label_DAQHeader_;
    std::string                         data_label_flash_beam_;
    std::string                         data_label_flash_cosmic_;
    std::string                         data_label_crthit_;
    std::string                         data_label_crtT0asso_;
    std::string                         trackProducerLabel;
    // pandora functions
};

DEFINE_ART_MODULE(CRTTriggerProducer)

} // namespace lar_pandora

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

#include "Pandora/PdgTable.h"

namespace lar_pandora
{

CRTTriggerProducer::CRTTriggerProducer(fhicl::ParameterSet const &pset) :
  EDProducer(pset)
{
    //control variables
    verbose_ = pset.get<int>("verbose");
    saveTTree_ = pset.get<int>("saveTTree");
    //datalabels
    data_label_DAQHeader_ = pset.get<std::string>("data_label_DAQHeader");
    data_label_flash_beam_ = pset.get<std::string>("data_label_flash_beam");
    data_label_flash_cosmic_ = pset.get<std::string>("data_label_flash_cosmic");
    data_label_crthit_ = pset.get<std::string>("data_label_crthit");
    data_label_crtT0asso_ = pset.get<std::string>("data_label_crtT0asso");
    
    trackProducerLabel = pset.get<std::string>("TrackProducerLabel");
    // crt variables
    store_t0_ = pset.get<int>("store_t0",0);

    if(store_t0_ == 1) produces< std::vector<anab::T0>   >();
    
    
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CRTTriggerProducer::produce(art::Event &evt)
{
  reset_tree_event();
  std::cout << "Prozessing event nr: " << event_counter << std::endl;
  if(verbose_!=0) std::cout << "Run " << evt.run() << ", subrun " << evt.subRun() << std::endl;
  frunNum    = evt.run();
  fsubRunNum = evt.subRun();
  fEvtNum = evt.event();
  event_counter++;
  
  art::Handle< raw::DAQHeaderTimeUBooNE > rawHandle_DAQHeader;
  evt.getByLabel(data_label_DAQHeader_, rawHandle_DAQHeader);
  raw::DAQHeaderTimeUBooNE const& my_DAQHeader(*rawHandle_DAQHeader);
  art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();  
  fTriTim_sec = evtTimeGPS.timeHigh();
  fTriTim_nsec = evtTimeGPS.timeLow();
  
  art::Timestamp evtTime = evt.time();
  fTimeHigh = evtTime.timeHigh();
  fTimeLow = evtTime.timeLow();
  
  art::Handle< std::vector<recob::OpFlash> > rawHandle_OpFlashBeam;
  evt.getByLabel(data_label_flash_beam_, rawHandle_OpFlashBeam);
  std::vector<recob::OpFlash> const& OpFlashCollectionBeam(*rawHandle_OpFlashBeam);
  if(verbose_!=0) std::cout << "There are: " << OpFlashCollectionBeam.size() << " in the beamflash collection" << std::endl;
  
  art::Handle< std::vector<recob::OpFlash> > rawHandle_OpFlashCosmic;
  evt.getByLabel(data_label_flash_cosmic_, rawHandle_OpFlashCosmic);
  std::vector<recob::OpFlash> const& OpFlashCollectionCosmic(*rawHandle_OpFlashCosmic);
  if(verbose_!=0) std::cout << "There are: " << OpFlashCollectionCosmic.size() << " in the cosmic flash collection" << std::endl;

  flash_tot = OpFlashCollectionCosmic.size() + OpFlashCollectionBeam.size();
  
  art::Handle<std::vector<recob::Track> > rawHandle_TPCtrack;
  evt.getByLabel(trackProducerLabel, rawHandle_TPCtrack);
  
  std::vector<recob::Track> const& TrackCollection(*rawHandle_TPCtrack); 
  
  art::FindMany<crt::CRTHit> trk_crt_assn_v(rawHandle_TPCtrack, evt, data_label_crtT0asso_);
  art::FindMany<anab::T0> trk_T0_assn_v(rawHandle_TPCtrack, evt, data_label_crtT0asso_);
  
  
  for(std::vector<int>::size_type i = 0; i != TrackCollection.size(); i++) {    //start loop over tpc tracks
    reset_tree();
    const std::vector<const crt::CRTHit*>& CRTHit_v = trk_crt_assn_v.at(i);
    for(std::vector<int>::size_type j = 0; j != CRTHit_v.size(); j++) {//CRThitloop
      nr_crthit++;
      crthit_ts0 = CRTHit_v.at(j)->ts0_ns;
      crthit_ts1 = CRTHit_v.at(j)->ts1_ns;
      adc_length = CRTHit_v.at(j)->pesmap.begin()->second.size();
      crt_adc = CRTHit_v.at(j)->peshit;
      t0_counter++;
    } //end crthit loop
    const std::vector<const anab::T0*>& T0_v = trk_T0_assn_v.at(i);
    for(std::vector<int>::size_type j = 0; j != T0_v.size(); j++) {//CRTt0 loop
      crtt0_time = T0_v.at(j)->Time();
      crtt0_trig = T0_v.at(j)->TriggerType();
      crtt0_DCA = T0_v.at(j)->TriggerConfidence();
      crtt0_plane = T0_v.at(j)->TriggerBits();
      
      for(std::vector<int>::size_type k = 0; k != OpFlashCollectionBeam.size(); k++) {//B
        recob::OpFlash my_flash = OpFlashCollectionBeam.at(k);
        if( abs( my_flash.Time() - crtt0_time)<5.0){
          fAbsTimFla = my_flash.Time();
          flash_PE = my_flash.TotalPE();
          flash_y = my_flash.YCenter();
          flash_z = my_flash.ZCenter();
          flash_inTime = 1;
          if(adc_length == 32){
            double time_diff = fAbsTimFla*1000 - (crthit_ts0 - fTriTim_nsec);
            flash_crthit_diff.push_back(time_diff);
            flash_counter++;
          }
        }
      } //B
      for(std::vector<int>::size_type k = 0; k != OpFlashCollectionCosmic.size(); k++) {//B
        recob::OpFlash my_flash = OpFlashCollectionCosmic.at(k);
        if( abs( my_flash.Time() - crtt0_time)<5.0){
          fAbsTimFla = my_flash.Time();
          flash_PE = my_flash.TotalPE();
          flash_y = my_flash.YCenter();
          flash_z = my_flash.ZCenter();
          flash_inTime = 0;
          if(adc_length == 32){
            double time_diff = fAbsTimFla*1000 - (crthit_ts0 - fTriTim_nsec);
            flash_crthit_diff.push_back(time_diff);
            flash_counter++;
          }
          
        }
      } //B
    //if(saveTTree_ == 1) my_event_->Fill();
    } //end crtt0 loop
  }//end loop over tpc tracks
  // calculation using mean, cut everything outside 1 std
  if(flash_crthit_diff.size()>1){
    double mean_corr = 0;
    double std_corr = 0;
    for(std::vector<int>::size_type i = 0; i != flash_crthit_diff.size(); i++) {//mean difference calculation
      mean_corr = mean_corr+flash_crthit_diff.at(i);
    }
    mean_corr = mean_corr/flash_crthit_diff.size();
    for(std::vector<int>::size_type i = 0; i != flash_crthit_diff.size(); i++) {//mean difference calculation
      std_corr = std_corr + (flash_crthit_diff.at(i)-mean_corr)*(flash_crthit_diff.at(i)-mean_corr);
    }
    std_corr = sqrt( std_corr/(flash_crthit_diff.size()-1) );
    int counter = 0;
    crthit_correction = 0;
    t0_outlayer = 0;
    for(std::vector<int>::size_type i = 0; i != flash_crthit_diff.size(); i++) {//mean difference calculation
      if( abs(flash_crthit_diff.at(i) - mean_corr)< std_corr){
        counter++;
        crthit_correction+=flash_crthit_diff.at(i);
      }
      else t0_outlayer++;
    }
    crthit_correction = crthit_correction/counter;
    if(counter==0) crthit_correction = mean_corr;
    if(flash_crthit_diff.size()==1) crthit_correction = flash_crthit_diff.at(0);
  }
  ////////////////////////////////////////////////////
  //calculation using median//////////////////////////
  if(flash_crthit_diff.size()>1){
    std::sort(flash_crthit_diff.begin(), flash_crthit_diff.end()); 
    for(std::vector<int>::size_type i = 0; i != flash_crthit_diff.size(); i++) {//mean difference calculation
      if(verbose_>0) std::cout << "i = " << i << " value: " << flash_crthit_diff.at(i) << std::endl;
    }
    int middle = (flash_crthit_diff.size()/2.0)+0.5;
    crthit_corr_med = flash_crthit_diff.at(middle-1);
    if(verbose_>0) std::cout << "Middle = " << middle << " median: " << crthit_corr_med << std::endl;
  }
  if(flash_crthit_diff.size()==1){
    crthit_corr_med = flash_crthit_diff.at(0);
  }
  /////////////////////////////////////////////////////
  my_event_trigger->Fill();
  
  // rerun over all t0 etc... now with the correction
  for(std::vector<int>::size_type i = 0; i != TrackCollection.size(); i++) {    //start loop over tpc tracks
    reset_tree();
    const std::vector<const crt::CRTHit*>& CRTHit_v = trk_crt_assn_v.at(i);
    for(std::vector<int>::size_type j = 0; j != CRTHit_v.size(); j++) {//CRThitloop
      nr_crthit++;
      crthit_ts0 = CRTHit_v.at(j)->ts0_ns;
      crthit_ts1 = CRTHit_v.at(j)->ts1_ns;
      adc_length = CRTHit_v.at(j)->pesmap.begin()->second.size();
      crt_adc = CRTHit_v.at(j)->peshit;
      t0_counter++;
    } //end crthit loop
    const std::vector<const anab::T0*>& T0_v = trk_T0_assn_v.at(i);
    for(std::vector<int>::size_type j = 0; j != T0_v.size(); j++) {//CRTt0 loop
      crtt0_time = T0_v.at(j)->Time();
      crtt0_trig = T0_v.at(j)->TriggerType();
      crtt0_DCA = T0_v.at(j)->TriggerConfidence();
      crtt0_plane = T0_v.at(j)->TriggerBits();
      
      for(std::vector<int>::size_type k = 0; k != OpFlashCollectionBeam.size(); k++) {//B
        recob::OpFlash my_flash = OpFlashCollectionBeam.at(k);
        if( abs( my_flash.Time() - crtt0_time)<5.0){
          fAbsTimFla = my_flash.Time();
          flash_PE = my_flash.TotalPE();
          flash_y = my_flash.YCenter();
          flash_z = my_flash.ZCenter();
          flash_inTime = 1;
        }
      } //B
      for(std::vector<int>::size_type k = 0; k != OpFlashCollectionCosmic.size(); k++) {//B
        recob::OpFlash my_flash = OpFlashCollectionCosmic.at(k);
        if( abs( my_flash.Time() - crtt0_time)<5.0){
          fAbsTimFla = my_flash.Time();
          flash_PE = my_flash.TotalPE();
          flash_y = my_flash.YCenter();
          flash_z = my_flash.ZCenter();
          flash_inTime = 0;
          double time_diff = fAbsTimFla*1000 - (crthit_ts0 - fTriTim_nsec);
          flash_crthit_diff.push_back(time_diff);
          flash_counter++;
          
        }
      } //B
    if(saveTTree_ == 1) my_event_->Fill();
    } //end crtt0 loop
  }//end loop over tpc tracks
  
  
  
  // produce anab::t0 object to tag event
  std::unique_ptr<std::vector<anab::T0> > T0_collection(new std::vector<anab::T0>);
  
  if(store_t0_ == 1){
    anab::T0 my_t0;
    my_t0.fTime = crthit_corr_med;  // double, corrected time used median
    my_t0.fTriggerType = t0_counter;    //unsigned int # of associated CRT hits to tracks
    my_t0.fTriggerBits = flash_counter;    //int # of matched flashes used...
    my_t0.fID = flash_tot;             //int # of flashes 
    my_t0.fTriggerConfidence = crthit_correction; //double, corrected time used mean with outlayer cut
    T0_collection->emplace_back(my_t0);
    evt.put(std::move(T0_collection));
  }
  
}


//------------------------------------------------------------------------------------------------------
void CRTTriggerProducer::initialize_tmyevent()
{
  // Implementation of optional member function here.
  std::cout << "Initialize variables and histograms for root tree output" << std::endl;
  // tree stuff for tracks: ////////////////////////////////////////////////////////////////////////////////// 
  my_event_ = tfs->make<TTree>("my_event","my_event");
  // event infos
  my_event_->Branch("run_MC", &run_MC, "run_MC/I");
  my_event_->Branch("event_counter", &event_counter, "event_counter/I");
  my_event_->Branch("frunNum", &frunNum, "Run Number/i");
  my_event_->Branch("fsubRunNum", &fsubRunNum, "SubRun Number/i");
  my_event_->Branch("fEvtNum", &fEvtNum, "Event Number/i");
  // DAQ time info
  my_event_->Branch("fTriTim_sec", &fTriTim_sec, "fTriTim_sec s/D");
  my_event_->Branch("fTriTim_nsec", &fTriTim_nsec, "fTriTim_nsec ns/D");
  
  my_event_->Branch("evt_time_sec", &fTimeHigh, "evt_time_sec/i");
  my_event_->Branch("evt_time_nsec", &fTimeLow, "evt_time_nsec/i");
  // Beam flash info
  my_event_->Branch("flash_time", &fAbsTimFla, "flash time us/D");
  my_event_->Branch("flash_PE", &flash_PE, "flash_PE ns/D");
  my_event_->Branch("flash_y", &flash_y, "flash_y cm/D");
  my_event_->Branch("flash_z", &flash_z, "flash_z cm/D");
  my_event_->Branch("flash_inTime", &flash_inTime, "is in BNB/I");
  // crt hit info
  my_event_->Branch("nr_crthit", &nr_crthit, "nr_crthit in BNB/I");
  my_event_->Branch("crthit_ts0", &crthit_ts0, "crthit_ts0 ns/D");
  my_event_->Branch("crthit_ts1", &crthit_ts1, "crthit_ts1 ns/D");
  my_event_->Branch("adc_length", &adc_length, "adc_length/I");
  my_event_->Branch("peshit", &crt_adc, "peshit/D");
  // asso tcrt t0 
  my_event_->Branch("crtt0_time", &crtt0_time, "crtt0_time us/D");
  my_event_->Branch("crtt0_trig", &crtt0_trig, "crtt0_trig = 2/I");
  my_event_->Branch("crtt0_DCA", &crtt0_DCA, "crtt0_DCA ncm/D");
  my_event_->Branch("crtt0_plane", &crtt0_plane, "crtt0_plane 3=top/I");
  
  //my_event_->Branch("fTriTim_sec", &fTriTim_sec, "fTriTim_sec s/D");
  //my_event_->Branch("fTriTim_nsec", &fTriTim_nsec, "fTriTim_nsec ns/D");
  //my_event_->Branch("event_counter", &event_counter, "event_counter/I");
  my_event_->Branch("flash_tot",      &flash_tot,          "flash_tot/I");
  my_event_->Branch("flash_counter",  &flash_counter,      "flash_counter/I");
  my_event_->Branch("t0_counter",     &t0_counter,         "t0_counter/I");
  my_event_->Branch("t0_outlayer",    &t0_outlayer,        "t0_outlayer/I");
  my_event_->Branch("crthit_corr",    &crthit_correction,  "crthit_corr/D");
  my_event_->Branch("crthit_corr_med",    &crthit_corr_med,  "crthit_corr_med/D");
  //general event info
  
  my_event_trigger = tfs->make<TTree>("my_trigger","my_trigger");
  my_event_trigger->Branch("fTriTim_sec", &fTriTim_sec, "fTriTim_sec s/D");
  my_event_trigger->Branch("fTriTim_nsec", &fTriTim_nsec, "fTriTim_nsec ns/D");
  my_event_trigger->Branch("event_counter", &event_counter, "event_counter/I");
  my_event_trigger->Branch("flash_tot",      &flash_tot,          "flash_tot/I");
  my_event_trigger->Branch("flash_counter",  &flash_counter,      "flash_counter/I");
  my_event_trigger->Branch("t0_counter",     &t0_counter,         "t0_counter/I");
  my_event_trigger->Branch("t0_outlayer",    &t0_outlayer,        "t0_outlayer/I");
  my_event_trigger->Branch("crthit_corr",    &crthit_correction,  "crthit_corr/D");
  my_event_trigger->Branch("crthit_corr_med",    &crthit_corr_med,  "crthit_corr_med/D");

}


void CRTTriggerProducer::reset_tree()
{
  fAbsTimFla = -1;
  flash_PE = -1;
  flash_y = -1;
  flash_z = -1;
  flash_inTime = 0;
  nr_crthit = 0;
  crthit_ts0 = -1;
  crthit_ts1 = -1;
  adc_length = -1;
  crt_adc = -1;
  crtt0_time = -1;
  crtt0_trig = -1;
  crtt0_DCA = -1;
  crtt0_plane = -1;
}
void CRTTriggerProducer::reset_tree_event()
{  
  flash_crthit_diff.clear();
  flash_tot = -1;
  flash_counter = 0;
  t0_counter = 0;
  t0_outlayer = 0;
  crthit_correction = -1;
  crthit_corr_med = -1;
}  
  
void CRTTriggerProducer::beginJob()
{
  // Implementation of optional member function here.
  std::cout << "Starting CRT trigger producer module" << std::endl;
  initialize_tmyevent();
  std::cout << "-------Using the following fcl parameters:-------" << std::endl;
  std::cout << "Pandora label:\t\t" << m_inputProducerLabel << std::endl;
  std::cout << "CRT T0 asso label:\t" << data_label_crtT0asso_ << std::endl;
  std::cout << "Beam flash label:\t" << data_label_flash_beam_ << std::endl;
  std::cout << "CRT hit label:\t\t" << data_label_crthit_ << std::endl;
  std::cout << "saveTTree:\t\t" << saveTTree_ << std::endl;
  std::cout << "run_MC:\t\t\t" << run_MC << std::endl;
  std::cout << "verbose:\t\t\t" << verbose_ << std::endl;
  std::cout << "store_t0:\t\t\t" << store_t0_ << std::endl;
  std::cout << "------------end fcl parameters-------------------" << std::endl;

}
void CRTTriggerProducer::endJob()
{
  std::cout << "Ending CRT trigger producer module" << std::endl;
  std::cout << "-------Using the following fcl parameters:-------" << std::endl;
  std::cout << "Pandora label:\t\t" << m_inputProducerLabel << std::endl;
  std::cout << "CRT T0 asso label:\t" << data_label_crtT0asso_ << std::endl;
  std::cout << "Beam flash label:\t" << data_label_flash_beam_ << std::endl;
  std::cout << "CRT hit label:\t\t" << data_label_crthit_ << std::endl;
  std::cout << "saveTTree:\t\t" << saveTTree_ << std::endl;
  std::cout << "run_MC:\t\t\t" << run_MC << std::endl;
  std::cout << "verbose:\t\t\t" << verbose_ << std::endl;
  std::cout << "store_t0:\t\t\t" << store_t0_ << std::endl;
  std::cout << "------------end fcl parameters-------------------" << std::endl;
}
} // namespace lar_pandora

