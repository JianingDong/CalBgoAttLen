/*   $Id: DmpEvtBgoAttLen.cc, 2014-09-29 20:28:57+08:00 DAMPE $
 *--------------------------------------------------------
 *  Author(s):
 *    Jianing Dong (jndong@mail.ustc.edu.cn) 29/09/2014
 *--------------------------------------------------------
*/

#include "DmpEvtBgoAttLen.h"

ClassImp(DmpEvtBgoAttLen)

DmpEvtBgoAttLen::DmpEvtBgoAttLen(){
  Reset();
}

DmpEvtBgoAttLen::DmpEvtBgoAttLen(const DmpEvtBgoAttLen &r){
  Reset();
  UsedFileName = r.UsedFileName;
  StartTime = r.StartTime;
  StopTime = r.StopTime;
  short n = GlobalBarID.size();
  for(size_t i = 0;i<n;++i){
    GlobalBarID.push_back(r.GlobalBarID[i]);
    Slope.push_back(r.Slope[i]);
    Intercept.push_back(r.Intercept[i]);
    Slp_Err.push_back(r.Slp_Err[i]);
    Inc_Err.push_back(r.Inc_Err[i]);
    ChiS.push_back(r.ChiS[i]);
  }
}

DmpEvtBgoAttLen::DmpEvtBgoAttLen(const DmpEvtBgoAttLen *&r){
  Reset();
  UsedFileName = r->UsedFileName;
  StartTime = r->StartTime;
  StopTime = r->StopTime;
  short n = GlobalBarID.size();
  for(size_t i = 0;i<n;++i){
    GlobalBarID.push_back(r->GlobalBarID[i]);
    Slope.push_back(r->Slope[i]);
    Intercept.push_back(r->Intercept[i]);
    Slp_Err.push_back(r->Slp_Err[i]);
    Inc_Err.push_back(r->Inc_Err[i]);
    ChiS.push_back(r->ChiS[i]);
  }
}

DmpEvtBgoAttLen::~DmpEvtBgoAttLen()
{
}

void DmpEvtBgoAttLen::Reset()
{
  UsedFileName = "NO";
  StartTime = 0;
  StopTime = 0xafffffff;
  GlobalBarID.clear();
  Slope.clear();
  Intercept.clear();
  Slp_Err.clear();
  Inc_Err.clear();
  ChiS.clear();
}
