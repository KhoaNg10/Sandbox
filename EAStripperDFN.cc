//TODO: Look at new electronics in all Sd Stations
//Use as starting point

#include "EAStripperDFN.h"
#include "sde_trigger_defs.h"
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

    int histId = eventId*10000 + fStationID*100 + pmtId;  // Create a unique ID for each histogram

    if(eventHists.find(histId) == eventHists.end()) 
    {
        eventHists[histId] = new TH1I(("eventHist_" + std::to_string(histId)).c_str(), 
                                      ("Histogram of Event " + std::to_string(eventId) + ", Station " + std::to_string(fStationID) + ", PMT " + std::to_string(pmtId)).c_str(), 
                                      2048, 0, 2048); // Create a new histogram for this station and PMT
    }

    TH1I* eventHist = eventHists[histId];
    
    int startIndex = traceLength/2;
    std::string secondHalfHistName = "secondHalf_event_" + std::to_string(eventId) + "_PMT_" + std::to_string(pmtId);
    TH1I* secondHalfHist = new TH1I(secondHalfHistName.c_str(), "secondHalf", traceLength - startIndex, startIndex, traceLength);

    for (i=0; i<traceLength; i++)
    {
        eventHist->Fill(i, trace[i]);
        
        if(i >= startIndex) 
        {
            secondHalfHist->Fill(i, trace[i]);
        }
    }
    
double rms = secondHalfHist->GetRMS(); // Compute RMS of the trace amplitude

// For debugging
std::cout << "RMS for event " << eventId << ", PMT " << pmtId << " is: " << rms << std::endl;

if (rms >= 1 && rms <= 5)
{
    std::cout << "Writing histogram for event " << eventId << ", PMT " << pmtId << std::endl;
    eventHist->Write();
}
else
{
    std::cout << "Skipping histogram for event " << eventId << ", PMT " << pmtId << std::endl;
}

// Deleting the second half histogram
delete secondHalfHist;

// Applying FIR filter to get the filtered trace
int fir[21] = {5,0,12,22,0,-61,-96,0,256,551,681,551,256,0,-96,-61,0,22,12,0,5};
int filt[2048];
for (i=21; i<2048; i++)
{
    filt[i] = 0;
    for (int j=0; j<21; j++)
    {
        filt[i] += trace[i-21+j] * fir[j];
    }
}

TFile filteredOutputFile("filteredtr.root", "UPDATE");
TH1I filteredHist(("filteredHist_" + std::to_string(histId)).c_str(), 
                 ("Filtered Histogram of Event " + std::to_string(eventId) + ", Station " + std::to_string(fStationID) + ", PMT " + std::to_string(pmtId)).c_str(), 
                 2048, 0, 2048);
for (i=0; i<2048; i++)
{
    filteredHist.Fill(i, filt[i]);
}
filteredHist.Write();
filteredOutputFile.Close();

// Downsampling the filtered trace
TFile downsampledOutputFile("downfilt.root", "UPDATE");
TH1I downsampledHist(("downsampledHist_" + std::to_string(histId)).c_str(), 
                    ("Downsampled Histogram of Event " + std::to_string(eventId) + ", Station " + std::to_string(fStationID) + ", PMT " + std::to_string(pmtId)).c_str(), 
                    2048/3, 0, 2048/3);
for (i=0; i<2048; i+=3)
{
    downsampledHist.Fill(i/3, filt[i]);
}
downsampledHist.Write();
downsampledOutputFile.Close();

    }
    }
    //std::cout << "Histogram is in file: " << (outputFile->Get("eventHist_" + std::to_string(eventId)) != nullptr) << std::endl;
  }
  return eSuccess;
}
VModule::ResultFlag EAStripperDFN::Finish()
{
outputFile->Write();
//std::cout << "Histogram is in file: " << (outputFile->Get("eventHist_" + std::to_string(eventId)) != nullptr) << std::endl;
outputFile->Close();
return eSuccess;
}
