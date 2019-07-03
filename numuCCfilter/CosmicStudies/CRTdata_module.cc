////////////////////////////////////////////////////////////////////////
// Class:       CRTdata
// Plugin Type: analyzer (art v2_11_03)
// File:        CRTdata_module.cc
// Mon Nov 12 2018 by Wouter Van De Pontseele
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

#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/CRT/CRTSimData.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"
#include <TTree.h>

class CRTdata;

class CRTdata : public art::EDAnalyzer
{
public:
  explicit CRTdata(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTdata(CRTdata const &) = delete;
  CRTdata(CRTdata &&) = delete;
  CRTdata &operator=(CRTdata const &) = delete;
  CRTdata &operator=(CRTdata &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;

private:
  bool m_is_data;
  std::string m_CRTHit_producer;
  const std::string m_DAQHeaderProducer = "daq";
  const float m_DTOffset = 68600.; // in microseconds

  // TTree Declaration.
  TTree *fCRTTree;

  // Event info.
  uint fRun, fSubrun, fEvent;
  uint fNhits_in_event;

  // CRT hit info.
  uint fPlane;
  float fTime;
  float fPE;
  float fX;
  float fX_err;
  float fY;
  float fY_err;
  float fZ;
  float fZ_err;
  double fT0_ns;
  double fT1_ns;
};

CRTdata::CRTdata(fhicl::ParameterSet const &p)
    : EDAnalyzer(p) // ,
                    // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;

  m_is_data = p.get<bool>("is_data", false);
  m_CRTHit_producer = p.get<std::string>("CRTHit_producer", "merger");

  std::cout << std::endl;
  std::cout << "[CRTdata constructor] CRTHit_producer: " << m_CRTHit_producer << std::endl;
  std::cout << "[CRTdata constructor] is_data: " << m_is_data << std::endl;

  //// Tree for every CRT hit
  fCRTTree = tfs->make<TTree>("CRT", "CRT Tree");

  fCRTTree->Branch("event", &fEvent, "event/i");
  fCRTTree->Branch("run", &fRun, "run/i");
  fCRTTree->Branch("subrun", &fSubrun, "subrun/i");
  fCRTTree->Branch("nhits", &fNhits_in_event, "nhits/i");

  fCRTTree->Branch("plane", &fPlane, "plane/i");
  fCRTTree->Branch("time", &fTime, "time/F");
  fCRTTree->Branch("pe", &fPE, "pe/F");
  fCRTTree->Branch("x", &fX, "x/F");
  fCRTTree->Branch("x_err", &fX_err, "x_err/F");
  fCRTTree->Branch("y", &fY, "y/F");
  fCRTTree->Branch("y_err", &fY_err, "y_err/F");
  fCRTTree->Branch("z", &fZ, "z/F");
  fCRTTree->Branch("z_err", &fZ_err, "z_err/F");

  fCRTTree->Branch("t0_ns", &fT0_ns, "time/D");
  fCRTTree->Branch("t1_ns", &fT1_ns, "time/D");
}

void CRTdata::analyze(art::Event const &e)
{
  // Fill the event info.
  fEvent = e.event();
  fSubrun = e.subRun();
  fRun = e.run();

  float evt_timeGPS_nsec = 0;

  if (m_is_data)
  {
    // Declare an object for the GPS timestamp of the event so that you can offset the cosmic t0 times.
    art::Handle<raw::DAQHeaderTimeUBooNE> rawHandle_DAQHeader;
    e.getByLabel(m_DAQHeaderProducer, rawHandle_DAQHeader);

    raw::DAQHeaderTimeUBooNE const &my_DAQHeader(*rawHandle_DAQHeader);
    art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
    evt_timeGPS_nsec = evtTimeGPS.timeLow();
  }

  // Load CRT hits.
  art::Handle<std::vector<crt::CRTHit>> crthit_h;
  e.getByLabel(m_CRTHit_producer, crthit_h);

  art::Handle<std::vector<crt::CRTSimData>> crtsimdata_h;
  e.getByLabel("crtdetsim", crtsimdata_h);

  if (!crtsimdata_h.isValid())
  {
    std::cout << "[CRTdata::analyze] crtsimdata_h is not valid!" << std::endl;
  }

  if (!crthit_h.isValid())
  {
    std::cout << "[CRTdata::analyze] ... could not locate CRT Hits, Event skipped!" << std::endl;
    return;
  }

  // Set the variable for the number of CRT hits in the event.
  fNhits_in_event = crthit_h->size();
  std::cout << "[CRTdata::analyze] Number of CRT hits in event:" << fNhits_in_event << std::endl;
  std::cout << std::endl;

  for (size_t j = 0; j < fNhits_in_event; j++)
  {
    fPlane = crthit_h->at(j).plane;
    fPE = crthit_h->at(j).peshit;
    fX = crthit_h->at(j).x_pos;
    fY = crthit_h->at(j).y_pos;
    fZ = crthit_h->at(j).z_pos;
    fX_err = crthit_h->at(j).x_err;
    fY_err = crthit_h->at(j).y_err;
    fZ_err = crthit_h->at(j).z_err;

    fT0_ns = crthit_h->at(j).ts0_ns;
    fT1_ns = crthit_h->at(j).ts1_ns;
    // Time of the CRT Hit wrt the event timestamp
    fTime = ((fT0_ns - evt_timeGPS_nsec + m_DTOffset) / 1000.);

    std::cout  <<  "fTime " << fTime;
    std::cout  <<  "\tfT1_ns " << fT1_ns;
    std::cout  <<  "\tevt_timeGPS_nsec " << evt_timeGPS_nsec << std::endl;

    fCRTTree->Fill();
  }
}

DEFINE_ART_MODULE(CRTdata)
