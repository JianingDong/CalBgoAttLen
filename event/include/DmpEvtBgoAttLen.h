/*   $Id: DmpEvtBgoAttLen.h, 2014-09-29 20:28:57+08:00 DAMPE $
 *--------------------------------------------------------
 *  Author(s):
 *    Jianing Dong (jndong@mail.ustc.edu.cn) 29/09/2014
 *--------------------------------------------------------
*/

#ifndef DmpEvtBgoAttLen_H
#define DmpEvtBgoAttLen_H

#include "TObject.h"

class DmpEvtBgoAttLen : public TObject{ 
/* 
 *  DmpEvtBgoAttLen
 *
 */
public:
  DmpEvtBgoAttLen();
  DmpEvtBgoAttLen(const DmpEvtBgoAttLen &r);
  DmpEvtBgoAttLen(const DmpEvtBgoAttLen *&r);
  ~DmpEvtBgoAttLen();

  void Reset();

public:
  std::string UsedFileName;
  int StartTime;
  int StopTime;
  std::vector<short> GlobalBarID;
  std::vector<double> Slope;
  std::vector<double> Intercept;

  std::vector<double> Slp_Err;
  std::vector<double> Inc_Err;

  std::vector<double> ChiS;

  ClassDef(DmpEvtBgoAttLen,1)

};

#endif
