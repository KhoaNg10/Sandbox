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
#include <cmath>
#include <vector>
#include <TH1F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <csignal>   // Added header for signal handling
#include <THStack.h>  // Added header for THStack

#include <stdio.h>
#include <sstream>

#include "sde_trigger_defs.h"
//#include "trigger_check.h"
//#include "TriggerParam.h"

/*#define SEQUENCE_MASK 0xFF           // Mask for the lowest 8 bits
#define ADC_SAMPLE_MASK 0xFFF        // Mask for 12 bits
#define FILTERED_TRACE_MASK 0x3      // Mask for the 2 least significant bits */


using namespace det;
using namespace evt;
using namespace sevt;
using namespace std;
using namespace utl;
using namespace fwk;

bool Enable[3];
int adcs[5][2048];
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
int filt[2048];


sevt::Station *CurrentStation;
sevt::SEvent::StationIterator stationIt;

//EventHistogramModule eventHistogramModule;
EAStripperDFN::EAStripperDFN() {}

EAStripperDFN::~EAStripperDFN() {}

TH1F* rmsHist = nullptr;
TFile *outputFile = nullptr;
TH1I* eventHist; // Histogram for event numbers
// Histogram for RMS values

std::map<int, TH1I*> eventHists;
std::vector<double> rmsValues; // Store RMS values
std::vector<double> frequencies; // Store frequencies
std::vector<std::pair<double, int>> rmsValuesWithIDs; // Store RMS values with histogram IDs

void SaveRMSGraph()
{
  if (!outputFile) return;

  // Debug: Print the number of RMS values and some sample values
  std::cout << "Number of RMS values: " << rmsValues.size() << std::endl;
  if (!rmsValues.empty()) {
    std::cout << "Sample RMS values: " << rmsValues[0] << ", " << rmsValues[1] << ", " << rmsValues[2] << std::endl;
  }

  // Find the range of RMS values
  double rmsMin = *std::min_element(rmsValues.begin(), rmsValues.end());
  double rmsMax = *std::max_element(rmsValues.begin(), rmsValues.end());

  std::cout << "RMS value range: [" << rmsMin << ", " << rmsMax << "]" << std::endl;

  // Create histogram for RMS values
  int nBins = 50; // Adjust number of bins as needed
  TH1F* rmsHist = new TH1F("rmsHist", "RMS vs Frequency", nBins, rmsMin, rmsMax);

  for (double rms : rmsValues) {
    rmsHist->Fill(rms);
  }

  rmsHist->SetFillColor(kBlue); // Set bar color to blue
  rmsHist->SetBarWidth(0.8); // Set bar width to 80% of bin width

  TCanvas *c1 = new TCanvas("c1", "RMS vs Frequency", 800, 600);
  rmsHist->GetXaxis()->SetTitle("RMS Value");
  rmsHist->GetYaxis()->SetTitle("Frequency");
  rmsHist->Draw("bar");

  c1->SaveAs("RMS_vs_Frequency.png");

  outputFile->cd();
  rmsHist->Write();
  outputFile->Write();
  outputFile->Close();
  delete c1;
}

void SignalHandler(int signal)
{
  if (signal == SIGINT || signal == SIGTSTP) {
    std::cout << "Interrupt signal received. Saving RMS vs Frequency graph." << std::endl;
    SaveRMSGraph();
    exit(signal);
  }
}

VModule::ResultFlag EAStripperDFN::Init()
{
  // Initialize your module here. This method 
  // is called once at the beginning of the run.
  // The eSuccess flag indicates the method ended
  // successfully.  For other possible return types, 
  // see the VModule documentation.
  INFO("EAStripperDFN::Init()");
  gStyle->SetOptStat(1111111);

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

  // Register signal handler
  std::signal(SIGINT, SignalHandler);
  std::signal(SIGTSTP, SignalHandler);

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
outputFile->cd();
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

    //TH1I* totalHistogram = new TH1I("totalHist", "Total Histogram", /* binning parameters */);
    //int startIndex = traceLength;

    for (i=0; i<traceLength; i++)
    {
    adcs[p][i] = (int) trace[i];
    std::cout << "Filling histogram for event " << eventId << " with value " << adcs[p][i] << std::endl;
    eventHist->Fill(i, adcs[p][i]);  // Fill the histogram with the ADC value at this time bin
    std::cout << "Number of entries for event " << eventId << ": " << eventHist->GetEntries() << std::endl;
    }

    eventHist->Write();


        // Calculate RMS for the first 300 bins
        double sum = 0;
        double sumSq = 0;
        int binCount = 300;
        for (int bin = 1; bin <= binCount; ++bin) 
        {
          double value = eventHist->GetBinContent(bin);
          sum += value;
          sumSq += value * value;
        }
        double mean = sum / binCount;
        double rms = std::sqrt(sumSq / binCount - mean * mean);

        std::cout << "Calculated RMS: " << rms << std::endl; // Debug: Print each calculated RMS


        // Store RMS value and corresponding frequency
        rmsValues.push_back(rms);
        frequencies.push_back(fStationID * 100 + pmtId);
    
    /* for  (j=0; j<nevents; ++j)
      {
        // Print time of event.
        fwrite(&linux_time[j],sizeof(int),1,output_file);
        fwrite(&secondsu[j],sizeof(int),1,output_file);
        fwrite(&delta_tics[j],sizeof(int),1,output_file);
        fwrite(&trig_bits[j],sizeof(int),1,output_file);
        for (i=0; i<3; i++)
          fwrite(&uncleaned_VEMs[j][i],sizeof(double),1,output_file);
        for (i=0; i<3; i++)
          fwrite(&cleaned_VEMs[j][i],sizeof(double),1,output_file);
        for (i=0; i<3; i++)
          fwrite(&baseline[j][i],sizeof(short),1,output_file);
        for (i=0; i<3; i++)
          fwrite(&totd_thres[j][i],sizeof(short),1,output_file);

        for (i=0; i<SHWR_MEM_WORDS; i++)
          {
        for (k=0; k<4; k++) // ssd is in k=3
          fwrite(&adcs[j][k][i],sizeof(u16),1,output_file);
          }
      }
      */

  }
  
    }
    
    } 
     return eSuccess;
}

VModule::ResultFlag EAStripperDFN::Finish()
{
  SaveRMSGraph();

  outputFile->Write();
  outputFile->Close();


  return eSuccess;
}