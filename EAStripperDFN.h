#ifndef _EAStripperDFN_h
#define _EAStripperDFN_h

#include <det/Detector.h>

#include <evt/Event.h>
#include <evt/ShowerRecData.h>
#include <evt/ShowerSRecData.h>
#include <evt/ShowerSimData.h>

#include <fwk/VModule.h>
#include <fwk/RunController.h>
#include <fwk/CentralConfig.h>
#include <fwk/LocalCoordinateSystem.h>

#include <sevt/SEvent.h>
#include <sevt/Header.h>
#include <sdet/SDetector.h>
#include <sevt/Station.h>
#include <sevt/PMT.h>
#include <sevt/PMTRecData.h>
#include <sevt/PMTCalibData.h>
#include <sevt/StationRecData.h>
#include <sevt/PMTSimData.h>
#include <sevt/StationSimData.h>
#include <sevt/StationTriggerData.h>

#include <utl/AugerUnits.h>
#include <utl/MathConstants.h>
#include <utl/PhysicalConstants.h>
#include <utl/TimeStamp.h>
#include <utl/TimeDistribution.h>
#include <utl/Particle.h>
#include <utl/TabulatedFunction.h>
#include <utl/ErrorLogger.h>
#include <utl/UTCDateTime.h>
#include <utl/TraceAlgorithm.h>

#include <fstream>
#include <string>
#include <vector>


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <unistd.h>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>

#include "TH1.h"
#include "TObject.h"
#include "IoSdData.h"

//#include "sde_trigger_defs.h"
//#include "trigger_check.h"
//#include "TriggerParam.h"


namespace evt {
  class Event;
}

namespace sevt {
  class Station;
  class SEvent;
}

/*class IoSdT2Trigger {
public:
  IoSdT2Trigger() {
    Type = Window = Offset = 0;
    Name = "";
  }
  virtual ~ IoSdT2Trigger() {}
  UInt_t type() const;              ///< Returns trigger type (see PLD documentation)
  UInt_t window() const;            ///< Returns window around the trigger time
  Int_t offset() const;             ///< Returns window offset with respect of T3 time

  UInt_t Type;                ///< Trigger type, 16 bits of 2 8 bits words
  bool IsTOTd() const {
    return 4;
  }
  UInt_t Window;              ///< Window around the trigger where stations should look for a trace to send back. Window=0 means that the station participated to the T3 trigger, an therefore has a T2
  Int_t Offset;               ///< Offset for the window. Contains the offset of the station with respect to the earliest one for Window=0 stations. Contains the average offset of all Window=0 stations for the others
  string Name;           ///< Trigger name.

  ClassDef(IoSdT2Trigger, 5)
} */



class EAStripperDFN : public fwk::VModule {

public:
  EAStripperDFN();
  // Beware to modify these functions if you add new fields to the IoSdStation class
  EAStripperDFN(const EAStripperDFN & st);
  EAStripperDFN & operator =(const EAStripperDFN & st);
  IoSdT2Trigger Trigger;

  virtual ~EAStripperDFN();

  fwk::VModule::ResultFlag Init();
  fwk::VModule::ResultFlag Run(evt::Event& event);
  fwk::VModule::ResultFlag Finish();


private:
  std::fstream cout_selected;
  std::fstream cout_rejected;
  long fEventNumber;
  long fStationID;
  bool isUUB;

  sevt::Station *fCurrentStation;
  sevt::SEvent::StationIterator stationIt;

  REGISTER_MODULE("EAStripperDFN", EAStripperDFN);

};



using namespace det;
using namespace evt;
using namespace fwk;
using namespace sevt;
using namespace utl;
using namespace std;

#endif

