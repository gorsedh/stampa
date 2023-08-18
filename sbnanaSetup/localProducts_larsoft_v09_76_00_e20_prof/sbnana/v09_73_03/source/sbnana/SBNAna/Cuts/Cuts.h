#pragma once

// Definition of the generic Cut object
#include "sbnana/CAFAna/Core/Cut.h"

namespace ana{

  // Select beam mode
  // extern const SpillCut kIsRHC;
  extern const SpillCut kFirstEvents;
  extern const SpillCut kFlashTrigger;

  extern const SpillCut kCRTHitVetoND;
  extern const SpillCut kCRTHitVetoFD;

  extern const Cut kActiveVolumeND;
  extern const Cut kFiducialVolumeND;

  extern const Cut kActiveVolumeFDCryo1;
  extern const Cut kActiveVolumeFDCryo2;
  extern const Cut kFiducialVolumeFDCryo1;
  extern const Cut kFiducialVolumeFDCryo2;

  extern const Cut kSlcIsRecoNu;
  extern const Cut kSlcNuScoreCut;

  extern const Cut kSlcHasFlashMatch;
  extern const Cut kSlcFlashMatchCut;

  const Cut kContainedFD  = kFiducialVolumeFDCryo1 || kFiducialVolumeFDCryo2;

} // namespace
