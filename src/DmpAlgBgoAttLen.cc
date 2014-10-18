#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "DmpAlgBgoAttLen.h"
#include "DmpDataBuffer.h"
#include "DmpBgoBase.h"
#include "DmpCore.h"
#include <iostream>
#include <fstream>

//-------------------------------------------------------------------
DmpAlgBgoAttLen::DmpAlgBgoAttLen()
 :DmpVAlg("Cal/Bgo/AttLen"),
  fEvtHeader(0),
  fBgoRaw(0),
  fBgoMips(0),
  fBgoAttLen(0)
{
//  Reset();	
}

//-------------------------------------------------------------------
bool DmpAlgBgoAttLen::GetMipsPar(){
  TFile *fMipPar=new TFile("./MIPs/MIPsPar.root");
  if(fMipPar==0){
    std::cout<<"Can not open MIPs Par root file!"<<std::endl;
    exit(0);
  }
  TTree *Miptree=(TTree*)fMipPar->Get("Calibration/Bgo");
  TBranch *b_fBgoMip;
  Miptree->SetBranchAddress("Mips",&fBgoMips,&b_fBgoMip);
  Miptree->GetEntry(0);
  //prepare parameters
  short gid=0,l=0,b=0,s=0;
  //std::cout<<"MipsPar---------------------"<<std::endl;
  short nBars=(short)fBgoMips->GlobalBarID.size();
  for(short i=0;i<nBars;++i){
    gid=fBgoMips->GlobalBarID[i];
    l=DmpBgoBase::GetLayerID(gid);
    b=DmpBgoBase::GetBarID(gid);
    s=fBgoMips->BgoSide[i];//s=0,1,2
    MipsPar[l][b][s][0]=fBgoMips->MPV[i]+0.222783*fBgoMips->Lwidth[i];//corrected MPV
    MipsPar[l][b][s][1]=fBgoMips->Gsigma[i];
    //std::cout<<"MipsPar="<<MipsPar[l][b][s][0]<<std::endl;
    MipsPar[l][b][s][2]=fBgoMips->Lwidth[i];
  }
  delete Miptree;
  delete fMipPar;
  //usage: QdcCoe[fGidOrder[gid]];//Quadratic Coefficients
  //       Slope[...],Cst[...] are same.
  return true;
}

//-------------------------------------------------------------------
DmpAlgBgoAttLen::~DmpAlgBgoAttLen(){
}
//-------------------------------------------------------------------
bool DmpAlgBgoAttLen::Initialize(){
  //read input data
  fEvtHeader = dynamic_cast<DmpEvtHeader*>(gDataBuffer->ReadObject("Event/Cutped/EventHeader"));
  fBgoRaw = dynamic_cast<DmpEvtBgoRaw*>(gDataBuffer->ReadObject("Event/Cutped/Bgo"));
  //create output data holder
  fBgoAttLen = new DmpEvtBgoAttLen();
  gDataBuffer->RegisterObject("Calibration/Bgo/AttLen",fBgoAttLen,"DmpEvtBgoAttLen");
  
  fBgoAttLen->UsedFileName = gRootIOSvc->GetInputFileName();
  gRootIOSvc->PrepareEvent(gCore->GetCurrentEventID());
  fBgoAttLen->StartTime = fEvtHeader->GetSecond();

  bool prepareMipPar = GetMipsPar();
  if(!prepareMipPar){
    std::cout<<"Error:Can not read Mips Par!"<<std::endl;
    return false;
  }
  //create Hist map
  for(short l=0;l<14;++l){
    for(short b=0;b<22;++b){
      char name[50];
      short gid_bar = DmpBgoBase::ConstructGlobalBarID(l,b);
      snprintf(name,50,"BgoAttLen_%05d-L%02d_B%02d",gid_bar,l,b);
      fAttLenHist.insert(std::make_pair(gid_bar,new TProfile(name,name,22,0,60)));
    }
  }
  return true;
}

//-------------------------------------------------------------------
bool DmpAlgBgoAttLen::Reset(){
  for(short layer=0;layer<14;++layer){
    for(short bar=0;bar<22;++bar){
      for(short side=0;side<2;++side){
        adc_Buf[layer][bar][side] = 0.;
      }
    }
  }
  for(short layer=0;layer<14;++layer){
    max_adc[layer] = 0.;
    //DmpLogInfo<<"layer="<<layer<<" max_adc[layer]="<<max_adc[layer-1];
    barID_max_adc[layer] = 22;
  }
  return true;
}

//-------------------------------------------------------------------
bool DmpAlgBgoAttLen::ProcessThisEvent(){
  fBgoAttLen->Reset();
  Reset();
  short gid = 0,l = -1,b = -1,s = -1,d = -1;
  double adc = 0.;
  //DmpLogInfo<<"layer="<<layerNo<<" max_adc[layer]="<<max_adc[layerNo-1];
  short nSignal = fBgoRaw->fGlobalDynodeID.size();
  for(short i=0;i<nSignal;++i){
    gid = fBgoRaw->fGlobalDynodeID[i];
    adc = fBgoRaw->fADC[i];
    //l = DmpBgoBase::GetLayerID(gid);
    //b = DmpBgoBase::GetBarID(gid);
    //s = DmpBgoBase::GetSideID(gid);
    DmpBgoBase::LoadLBSDID(gid,l,b,s,d);
    //DmpLogInfo<<"layer="<<layerNo<<" max_adc[layer]="<<max_adc[layerNo-1];
    if(b > 21){
      continue;
    }
    if(d == 8){
      adc_Buf[l][b][s] = adc;
      if(s == 0){
        if(adc > max_adc[l]){
	  barID_max_adc[l] = b;
	  max_adc[l] = adc;
	}
      }
    }
  }
  if((max_adc[0]>200||max_adc[1]>200)&&(max_adc[12]>200||max_adc[13]>200)){
    TF1 *linear = new TF1("linear","[0]*x+[1]",-1,24);
    linear->SetParNames("Par0","Par1");
    TH2F *fTrackHist[2];
    fTrackHist[0] = new TH2F("fTrack_X","fTrack_X",160,-1,15,240,-1,23);
    fTrackHist[1] = new TH2F("fTrack_Y","fTrack_Y",160,-1,15,240,-1,23);
    for(short l=0;l<14;++l){
      if(l%2 == 0){
        fTrackHist[1]->Fill(l+1,barID_max_adc[l]);
      }
      else{
	fTrackHist[0]->Fill(l+1,barID_max_adc[l]);
      }
    } 
    double Par0 = 0.,Par1 = 0.;
    double HitPos[14];
    memset(HitPos,0,sizeof(HitPos));
    for(short idim=0;idim<2;++idim){
      fTrackHist[idim]->Fit(linear,"RQ0");
      Par0=linear->GetParameter(0);
      Par1=linear->GetParameter(1);
      linear->SetParameter(0,Par0);
      linear->SetParameter(1,Par1);
      fTrackHist[idim]->Fit(linear,"RQ0");
      Par0=linear->GetParameter(0);
      Par1=linear->GetParameter(1);
      short ChiS=linear->GetChisquare();
      if(ChiS<10){
        for(short l=0;l<14;++l){
	  if(l%2 == idim){
	    HitPos[l]=((Par0*(l+1)+Par1)-0.5)*2.75;
	  }
        }
      }
      else{
        memset(HitPos,0,sizeof(HitPos));
      }
    }
    delete linear;
    delete fTrackHist[0];
    delete fTrackHist[1];

    double AttPar0[14];
    double AttPar1[14];
    memset(AttPar0,0,sizeof(AttPar0));
    memset(AttPar1,0,sizeof(AttPar1));
    for(short l=0;l<14;++l){
      short ib = barID_max_adc[l];
      if(ib == 22)continue;
      short gid_bar = DmpBgoBase::ConstructGlobalBarID(l,ib);
      AttPar0[l] = adc_Buf[l][ib][0]/MipsPar[l][ib][0][0];
      AttPar1[l] = adc_Buf[l][ib][1]/MipsPar[l][ib][1][0];
      //double AttPar = AttPar0/AttPar1;
      if(AttPar1[l] != 0 && (AttPar0[l]/AttPar1[l]) > 0 && HitPos[l] != 0){
	fAttLenHist[gid_bar]->Fill(HitPos[l],TMath::Log(AttPar0[l]/AttPar1[l]));
      }
    }
  }
  Reset();
  return true;
}

//-------------------------------------------------------------------
bool DmpAlgBgoAttLen::Finalize(){
  TF1 *LinearFit = new TF1("LinearFit","[0]*x+[1]",0,60);
  LinearFit->SetParNames("Intercept","Slope");
  LinearFit->SetRange(15,45);
  std::string histFileName = "./AttLen/Histograms/"+gRootIOSvc->GetOutputStem()+"_Hist.root";
  TFile *histFile = new TFile(histFileName.c_str(),"RECREATE");
  fBgoAttLen->StopTime = fEvtHeader->GetSecond();
  for(std::map<short,TProfile*>::iterator aHist=fAttLenHist.begin();aHist!=fAttLenHist.end();++aHist){
    fBgoAttLen->GlobalBarID.push_back(aHist->first);
    aHist->second->Fit(LinearFit,"R");
    double Slope = LinearFit->GetParameter(0),Intercept = LinearFit->GetParameter(1);
    double Slp_Err = LinearFit->GetParError(0),Inc_Err = LinearFit->GetParError(1);
    double ChiS = LinearFit->GetChisquare();
    fBgoAttLen->Slope.push_back(Slope);
    fBgoAttLen->Slp_Err.push_back(Slp_Err);
    fBgoAttLen->Intercept.push_back(Intercept);
    fBgoAttLen->Inc_Err.push_back(Inc_Err);
    fBgoAttLen->ChiS.push_back(ChiS);

    aHist->second->Write();
    delete aHist->second;
  }
  delete histFile;
  return true;
}

