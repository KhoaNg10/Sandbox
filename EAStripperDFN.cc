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
//#include "IoSdData.h"


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
//int TargetStation;
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
int trig_id;

sevt::Station *CurrentStation;
sevt::SEvent::StationIterator stationIt;

//EventHistogramModule eventHistogramModule;
EAStripperDFN::EAStripperDFN() {}

EAStripperDFN::~EAStripperDFN() {}



TFile *outputFile = nullptr;
TH1I* eventHist; // Histogram for event numbers
std::map<int, TH1I*> eventHists;

VModule::ResultFlag EAStripperDFN::Init()
{
  // Initialize your module here. This method 
  // is called once at the beginning of the run.
  // The eSuccess flag indicates the method ended
  // successfully.  For other possible return types, 
  // see the VModule documentation.
  INFO("EAStripperDFN::Init()");
  gStyle->SetOptStat(0);

  // Get parameters from the XML configuration file.

  CentralConfig *theConfig = CentralConfig::GetInstance();
  Branch topBranch = theConfig->GetTopBranch("EAStripperDFN"),
    dataBranch;
 
  //topBranch.GetChild("EnableDumpHistos").GetData(DumpHistos);
  //topBranch.GetChild("CurrentStation").GetData(CurrentStation);

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

  return eSuccess;
}

VModule::ResultFlag EAStripperDFN::Run(Event& event)
{ 
  ostringstream info;
  //bool reject;
  double theta;
  int nPMTs;
  int i, p;
  int traceLength;

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

  // Now use the histogram
  eventHist = eventHists[eventId];
  fStationID = 0;
  for (stationIt = sevent.StationsBegin();
       stationIt != sevent.StationsEnd(); stationIt++)
    {
  fCurrentStation = &(*stationIt);
  fStationID = fCurrentStation->GetId();
  bool isCurrentStationUUB = fCurrentStation->IsUUB();
  const int firstPMT = sdet::Station::GetFirstPMTId();
  nPMTs = 0;
  
  // Accessing the trigger data
  const sevt::StationTriggerData& trig = fCurrentStation->GetTriggerData();
  trig_id = trig.GetPLDTrigger();
  
  // Define the ToTD triggers condition (adjust as per your need)
  int ToTDTriggers = StationTriggerData::ePLDLatchTOTB | StationTriggerData::ePLDTOTB;

  if (!isCurrentStationUUB || (trig_id & ~ToTDTriggers) != 0) continue;
  
  //std::cout << "Event " << eventId << " is ToTD" << std::endl;
  //if (isUUB) 
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
                                      2048, 0, 2048); // Create a new histogram for this station and PMT
      }

      // Now use the histogram
      TH1I* eventHist = eventHists[histId];
      //selecting specific RMS range
      for (i=0; i<traceLength; i++)
      {
        //eventHist->Reset();
        adcs[p][i] = (int) trace[i];
        eventHist->Fill(i, adcs[p][i]);  // Fill the histogram with the ADC value at this time bin
        //eventHist->Write(); 
        //std::cout << "Filling histogram for event " << eventId << " with value " << adcs[p][i] <<  " corresponding" << i  << std::endl;
      }

      // Create a histogram for the second half of the bins
      int startIndex = traceLength/2;
      TH1I* secondHalfHist = new TH1I("secondHalf", "secondHalf", traceLength - startIndex, startIndex, traceLength);
      for (i=startIndex; i<traceLength; i++)
      {
        secondHalfHist->Fill(i, adcs[p][i]);
      }
      double rms = secondHalfHist->GetRMS(); // Compute RMS of the trace amplitude

          if (rms >= 141 && rms <= 213)
    {
      // ... (existing code for processing the PMT) ...
      Base[p] = pmt.GetCalibData().GetBaseline();
      Baseline[p] = Base[p] * 16 + .5; // convert to 12-bit ADC counts
      VEMpk[p] = pmt.GetCalibData().GetVEMPeak();
      info << ":" << VEMpk[p]; 
      if (pmt.GetCalibData().IsTubeOk()) Enable[p] = true;
      if (!pmt.GetCalibData().IsTubeOk()) info << "bad:";
      eventHist->Write();
      delete secondHalfHist;  // clean up the second half histogram
      return eSuccess;
    }
    else
    {
      // Histogram doesn't meet the criteria, so delete it
      delete eventHists[histId]; // Delete the histogram object
      eventHists.erase(histId); // Remove the entry from the map
      std::cout << "Deleting histogram for event " << eventId << ", Station " << fStationID << ", PMT " << pmtId << " with RMS out of range: " << rms << std::endl;
      continue; // Skip the remaining processing for this PMT
      delete secondHalfHist; // Delete second half of the histogram object
      return eContinueLoop;       
    }
      // old part
      /*Base[p] = pmt.GetCalibData().GetBaseline();  
	    Baseline[p] = Base[p]*16 + .5; // convert to 12-bit ADC counts
	    VEMpk[p] = pmt.GetCalibData().GetVEMPeak(); 
	    info << ":" << VEMpk[p]; 
	    if (pmt.GetCalibData().IsTubeOk()) Enable[p] = true;
	    if (!pmt.GetCalibData().IsTubeOk()) info << "bad:";
      //const PMTCalibData& pmtCalib = pmt.GetCalibData();
      eventHist->Write();*/
      std::cout << "Filling histogram for event " << eventId << " with value " << adcs[p][i] <<  " corresponding" << i  << std::endl;
    }
    }
    //std::cout << "Histogram is in file: " << (outputFile->Get("eventHist_" + std::to_string(eventId)) != nullptr) << std::endl;
  }
  return eContinueLoop;
}
VModule::ResultFlag EAStripperDFN::Finish()
{
outputFile->Write();
//std::cout << "Histogram is in file: " << (outputFile->Get("eventHist_" + std::to_string(eventId)) != nullptr) << std::endl;
outputFile->Close();
return eSuccess;
}
