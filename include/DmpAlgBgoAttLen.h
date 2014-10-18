#ifndef DmpAlgBgoAttLen_H
#define DmpAlgBgoAttLen_H

#include "DmpVAlg.h"
#include "DmpEvtHeader.h"
#include "DmpEvtBgoRaw.h"
#include "DmpEvtBgoDyCoe.h"
#include "DmpEvtBgoMips.h"
#include "DmpEvtBgoAttLen.h"
#include "TProfile.h"
#include <map>

class DmpEvtHeader;
class DmpEvtBgoRaw;
class DmpEvtBgoMips;
class DmpEvtBgoAttLen;

class DmpAlgBgoAttLen : public DmpVAlg{
/*
 *  DmpAlgBgoAttLen
 *
 */
public:
  DmpAlgBgoAttLen();
  ~DmpAlgBgoAttLen();

  //void Set(const std::string &type,const std::string &value);
  // if you need to set some options for your algorithm at run time. Overload Set()
  bool GetMipsPar();
  bool Initialize();
  bool ProcessThisEvent();    // only for algorithm
  bool Finalize();
  bool Reset();

private:
  DmpEvtHeader          *fEvtHeader;
  DmpEvtBgoRaw          *fBgoRaw;
  DmpEvtBgoMips         *fBgoMips;
  DmpEvtBgoAttLen       *fBgoAttLen;
  std::map<short,TProfile*>  fAttLenHist;

  double MipsPar[14][22][3][3];//layer,bar,side :2 mean combined value,0:MPV 1:Gsigma 2:Lwidth
  double adc_Buf[14][22][2];
  double max_adc[14];
  short barID_max_adc[14];
};

#endif
