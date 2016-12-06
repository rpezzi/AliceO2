//
//  DigitizerTask.h
//  ALICEO2
//
//  Created by Markus Fasel on 16.07.15.
//
//

#ifndef __ALICEO2__DigitizerTask__
#define __ALICEO2__DigitizerTask__

#include <stdio.h>
#include "FairTask.h"  // for FairTask, InitStatus
#include "Rtypes.h"    // for DigitizerTask::Class, ClassDef, etc

#include "ITSSimulation/Digitizer.h"

class TClonesArray;
namespace AliceO2 { namespace ITS { class Digitizer; }}  // lines 19-19

namespace AliceO2 {
namespace ITS {

class Digitizer;

class DigitizerTask : public FairTask
{
  public:
    DigitizerTask(Bool_t useAlpide=kFALSE);

    virtual ~DigitizerTask();

    virtual InitStatus Init();

    virtual void Exec(Option_t *option);

    Digitizer &GetDigitizer() { return fDigitizer; }

  private:
    Bool_t fUseAlpideSim;
    Digitizer fDigitizer;

    TClonesArray *fPointsArray;
    TClonesArray *fDigitsArray;

  ClassDef(DigitizerTask, 1)
};
}
}

#endif /* defined(__ALICEO2__DigitizerTask__) */
