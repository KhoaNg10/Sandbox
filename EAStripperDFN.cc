//TODO: Look at new electronics in all Sd Stations
//Use as starting point

#include "EAStripperDFN.h"
//#include "sde_trigger_defs.h"
//#include "trigger_check.h"
//#include "TriggerParam.h"
#include <TStyle.h>
#include <TH1.h>
#include <TFile.h>
#include <TProfile.h>
#include "IoSdData.h"
//#include <THist.h>
//#include <TH1I> //placeholder
#include "IoSdData.h"


#include <stdio.h>
#include <sstream>

//#include "sde_trigger_defs.h"
//#include "trigger_check.h"
//#include "TriggerParam.h"

using namespace det;
using namespace evt;
using namespace sevt;
using namespace std;
using namespace utl;
using namespace fwk;

bool Enable[3];
int adcs[5][2048];
int FiltTrace[3][768];
double Base[3];
int Baseline[5];
double VEMpk[3];
int TargetStation;
int FirstMonth;
int Year;
int LastMonth;
//int LastYear;
int FirstDay;
int FirstHour;
int FirstMinute;
int LastDay;
int LastHour;
int LastMinute;
//int dt;


// public:
//     bool* getEnable() { return Enable; }
//     int* getAdcs(int index) { return adcs[index]; }

sevt::Station *CurrentStation;
sevt::SEvent::StationIterator stationIt;

//EventHistogramModule eventHistogramModule;
EAStripperDFN::EAStripperDFN() {}

EAStripperDFN::~EAStripperDFN() {}



TFile *outputFile; //= nullptr;
TH1I* eventHist; // Histogram for event numbers
std::map<int, TH1I*> eventHists;

//



/*EAStripperDFN & EAStripperDFN::operator =(const EAStripperDFN & st) {
    if (this == &st)
    return *this;
  Trigger = st.Trigger;
}

EAStripperDFN::EAStripperDFN(const EAStripperDFN & st) {
  Trigger = st.Trigger;
}

IoSdT2Trigger EAStripperDFN::trigger() const {
  return Trigger;
}

int EAStripperDFN::TriggerType() {
  if (Trigger.IsTOTd())
    return 4;
}*/




VModule::ResultFlag EAStripperDFN::Init()
{
  // Initialize your module here. This method 
  // is called once at the beginning of the run.
  // The eSuccess flag indicates the method ended
  // successfully.  For other possible return types, 
  // see the VModule documentation.
  //outputTree->Branch("sevent", &event.GetSEvent(), 32000, 99);
  //eventHistogramModule.Init();
  //  char name [20];
 // char name2 [80];
  
  INFO("EAStripperDFN::Init()");
  gStyle->SetOptStat(0);

  // Get parameters from the XML configuration file.

  CentralConfig *theConfig = CentralConfig::GetInstance();
  Branch topBranch = theConfig->GetTopBranch("EAStripperDFN"),
    dataBranch;
 
  //topBranch.GetChild("EnableDumpHistos").GetData(DumpHistos);
  topBranch.GetChild("TargetStation").GetData(TargetStation);

  topBranch.GetChild("FirstMonth").GetData(FirstMonth);
  topBranch.GetChild("Year").GetData(Year);
  topBranch.GetChild("LastMonth").GetData(LastMonth);
  //topBranch.GetChild("LastYear").GetData(FirstYear);
  topBranch.GetChild("FirstDay").GetData(FirstDay);
  topBranch.GetChild("FirstHour").GetData(FirstHour);
  topBranch.GetChild("FirstMinute").GetData(FirstMinute);
  topBranch.GetChild("LastDay").GetData(LastDay);
  topBranch.GetChild("LastHour").GetData(LastHour);
  topBranch.GetChild("LastMinute").GetData(LastMinute);

  cout_selected.open("SelectedEvents.txt", std::fstream::out);
  cout_rejected.open("RejectedEvents.txt", std::fstream::out);
  outputFile = new TFile("output_hist.root", "RECREATE");
  //eventHist = new TH1I("eventHist", "Histogram of Event Numbers", 100, 0, 20000); // adjust binning as needed
  return eSuccess;
}

VModule::ResultFlag EAStripperDFN::Run(Event& event)
{ 
  ostringstream info;
  bool reject;
  double theta;
  int nPMTs;
  int i, p;
  int traceLength;
  // Run EventHistogramModule
  //eventHistogramModule.Run(event);

  
  if (!event.HasSEvent())
    {
      INFO("Event doesn't have surface detector event");
      return eFailure;
    }
if (!event.HasSEvent ()) return eContinueLoop;
sevt::SEvent& theSEvent = event.GetSEvent();
fEventNumber = theSEvent.GetHeader().GetId();
const utl::TimeStamp& eventTime = theSEvent.GetHeader().GetTime();
const utl::TimeStamp& firstTime = UTCDateTime(Year,FirstMonth,FirstDay,FirstHour,FirstMinute).GetTimeStamp();
const utl::TimeStamp& lasttime = UTCDateTime(Year,LastMonth,LastDay,LastHour,LastMinute).GetTimeStamp();

unsigned long t = eventTime.GetGPSSecond();
unsigned long t0 = firstTime.GetGPSSecond();
unsigned long t1 = lasttime.GetGPSSecond();

//double dt = (t - t0)/60./60.;

if ((t < t0) || (t >t1))
  {
    info.str("");
    info << "Event out of time window";
    INFO(info);
    return eContinueLoop;
  }

  // Pick out some pieces of the data.

  theta = 0.;
  sevt::SEvent& sevent = event.GetSEvent();

  if (event.HasRecShower()) 
    {
  // Pick out some pieces of the data.

  const ShowerSRecData& showersrec = event.GetRecShower().GetSRecShower();
  const Vector& axis = showersrec.GetAxis();
  const CoordinateSystemPtr localCS =
  LocalCoordinateSystem::Create(showersrec.GetCorePosition()); 
  theta = axis.GetTheta(localCS)/degree;

      info.str("");
      info << "New reconstruced event: Zenith=" << theta;
      INFO(info);
    }

  //Trigger = ToTD;
  isUUB = false;
  //IsToTD = true;
  int eventId = sevent.GetHeader().GetId();

  // Create a new histogram for this event number if one doesn't exist
  if(eventHists.find(eventId) == eventHists.end()) {
    eventHists[eventId] = new TH1I(("eventHist_" + std::to_string(eventId)).c_str(), 
                                    ("Histogram of Event " + std::to_string(eventId)).c_str(), 
                                    100, 0, 2048); // adjust binning as needed
  }

  // Now use the histogram
  eventHist = eventHists[eventId];
  
  for (stationIt = sevent.StationsBegin();
       stationIt != sevent.StationsEnd(); stationIt++)
    {
  fCurrentStation = &(*stationIt);
  fStationID = fCurrentStation->GetId();
  //if (fStationID != 56 && fStationID != 1733 && fStationID != 1736 && fStationID != 1737 && fStationID != 1738 && fStationID != 1742 && fStationID != 1744 && fStationID != 734) continue;
  if (!Trigger.IsTOTd()) continue;         
  const int firstPMT = sdet::Station::GetFirstPMTId();
  nPMTs = 0;

  for (p = 0; p < 3; p++) 
  {
    const int pmtId = p + firstPMT;
    const PMT& pmt = stationIt->GetPMT(pmtId);
    Enable[p] = false;

    if (pmt.HasFADCTrace())
    {
      nPMTs++;
      const TraceI& trace = pmt.GetFADCTrace(sdet::PMTConstants::eHighGain);
      traceLength = trace.GetSize();

      // create a new histogram for this station and PMT if one doesn't exist
      int histId = eventId*10000 + fStationID*100 + pmtId;  // Create a unique ID for each histogram
      if(eventHists.find(histId) == eventHists.end()) {
        eventHists[histId] = new TH1I(("eventHist_" + std::to_string(histId)).c_str(), 
                                      ("Histogram of Event " + std::to_string(eventId) + ", Station " + std::to_string(fStationID) + ", PMT " + std::to_string(pmtId)).c_str(), 
                                      2048, 0, 2048); // adjust binning as needed
      }

      // Now use the histogram
      TH1I* eventHist = eventHists[histId];

      for (i=0; i<traceLength; i++)
      {
        //eventHist->Reset();
        adcs[p][i] = (int) trace[i];
        eventHist->Fill(i, adcs[p][i]);  
        //eventHist->Write();
        std::cout << "Filling histogram for event " << eventId << " with value " << adcs[p][i] <<  " corresponding" << i  << std::endl;
      }
	      eventHist->Write();
      Base[p] = pmt.GetCalibData().GetBaseline();
	    Baseline[p] = Base[p]*16 + .5;
	    VEMpk[p] = pmt.GetCalibData().GetVEMPeak();
	    info << ":" << VEMpk[p];
	    if (pmt.GetCalibData().IsTubeOk()) Enable[p] = true;
	    if (!pmt.GetCalibData().IsTubeOk()) info << "bad:";
    }
      std::cout << "Filling histogram for event " << eventId << " with value " << adcs[p][i] <<  " corresponding" << i  << std::endl;

      }
      isUUB = isUUB || fCurrentStation->IsUUB();
      //IsToTD = IsToTD || Trigger.IsTOTd();
      //Trigger.IsTOTd;
    }
  //not necessary if the selected stations already has UUB
  if (isUUB)
    reject = false;
    //std::cout << "has UUB and ToTD" <<
  else
    reject = true;

  if (reject)
    {
      cout_rejected << "Event " << sevent.GetHeader().GetId()
                    << " Theta=" << theta << endl;
      //std::cout << "station" << std::to_string(fStationID) + "has UUB" << std::endl;
      return eContinueLoop;
      
      //return eSuccess;
    }
  else 
    {
      cout_selected << "Event " << sevent.GetHeader().GetId()
                    << " Theta=" << theta << endl;
  return eSuccess;
    }
  //add a new module that will selects stations that have ToTD
  
}


VModule::ResultFlag EAStripperDFN::Finish()
{
outputFile->Write();
//std::cout << "Histogram is in file: " << (outputFile->Get("eventHist_" + std::to_string(eventId)) != nullptr) << std::endl;
outputFile->Close();
//TFile *checkFile = new TFile("output_hist.root", "READ");
//std::cout << "Number of keys in file: " << checkFile->GetNkeys() << std::endl;
//checkFile->Close();
//delete checkFile;
return eSuccess;
}

// Configure (x)emacs for this file ...
// Local Variables:
// mode: c++
// End:
