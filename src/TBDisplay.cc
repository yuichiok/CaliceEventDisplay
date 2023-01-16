#define TBDisplay_cxx

#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveVSD.h>
#include <TEveVSDStructs.h>

#include <TEveTrack.h>
#include <TEveTrackPropagator.h>
#include <TEveGeoShape.h>
#include <TEventList.h>

#include <TGTab.h>
#include <TGButton.h>

#include <TFile.h>
#include <TKey.h>
#include <TSystem.h>
#include <TPRegexp.h>

#include <algorithm> // for std::find
#include <iterator> // for std::begin, std::end
#include <string>
#include <map>

#include "../include/TBDisplay.hh"
#include "../include/Fit3D.hh"

using std::cout;
using std::endl;

const bool debug = false;
const int nslabs = 15;
const int nscas = 15;
const float beamX = 20.0, beamY = 15.0;
const float MARKER_SIZE = 3.5;

// define the parametric line equation
void line(double t, const double *p, double &x, double &y, double &z) {
   // a parametric line is define from 6 parameters but 4 are independent
   // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
   // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
  //  x = p[0] + p[1]*t;
  //  y = p[2] + p[3]*t;
   x = p[0];
   z = p[1];
   y = t;
}

void TBDisplay::Next()
{
  GotoEvent(fCurEv + 1);
}

void TBDisplay::Prev()
{
  GotoEvent(fCurEv - 1);
}

void TBDisplay::GoTo() {

  int goto_eventNum;
  cout << "Go To: ";
  cin >> goto_eventNum;
  cout << endl;

  GotoEvent(goto_eventNum);

}

void TBDisplay::DropEvent()
{   
  gEve->GetViewers()->DeleteAnnotations();
  gEve->GetCurrentEvent()->DestroyElements();
}

void TBDisplay::ClearObjects()
{   
  TSeqCollection* canvases = gROOT->GetListOfCanvases();
  TIter next(gROOT->GetListOfCanvases());
  while(TCanvas *c = (TCanvas*)next())
  {
    delete c;
  }

  for (int islab=0; islab<nslabs; islab++){
    TString h2name = TString::Format("hitmap layer%i",islab);
    delete gROOT->FindObject(h2name);
  }
}

Bool_t TBDisplay::GotoEvent(Int_t ev)
{
  Long64_t nb = 0;
  if (fChain == 0) return kFALSE;

  Int_t nentries = fMaxEv;

  if (ev < 0 || ev >= nentries)
  {
    Warning("GotoEvent", "Invalid event id %d.", ev);
    return kFALSE;
  }

  DropEvent();
  ClearObjects();

  cout << endl;
  cout << "Going to " << ev << "..." << endl;
  cout << endl;

  fCurEv = ev;

  Long64_t ientry = LoadTree( evlist->GetEntry(ev) );
  if (ientry == 0) {
    Warning("GotoEvent", "Entry is empty");
    return kFALSE;
  }

  nb = fChain->GetEntry(evlist->GetEntry(ev));

  // Profile
  Float_t X0s[nslabs] = {1.198630137, 2.397260274, 3.595890411, 4.794520548, 5.993150685, 7.191780822, 8.390410959, 9.589041096, 10.78767123, 12.38584475, 13.98401826, 15.58219178, 17.1803653, 18.77853881, 20.37671233};
  Float_t E_layers[nslabs] = {0};
  vector<Float_t> layer_hit_x[nslabs];
  vector<Float_t> layer_hit_y[nslabs];

  // 2D Hit Maps
  TH2F *hitmap[nslabs];
  for (int islab=0; islab<nslabs; islab++){
    TString h2name = TString::Format("hitmap layer%i",islab);
    hitmap[islab] = new TH2F(h2name,h2name,32,-90,90,32,-90,90);
  }

  // Load event data into visualization structures.
  cout << "number of hits: " << nhit_len   << endl;
  cout << "sum energy    : " << sum_energy << endl;
  for (int ihit=0; ihit<nhit_len; ihit++){
    
    E_layers[hit_slab[ihit]] += hit_energy[ihit];
    layer_hit_x[hit_slab[ihit]].push_back(hit_x[ihit]);
    layer_hit_y[hit_slab[ihit]].push_back(hit_y[ihit]);

    hitmap[hit_slab[ihit]]->Fill(hit_x[ihit], hit_y[ihit], hit_energy[ihit]);

    LoadHits_Box(fHits_Box,ihit);

  } // hit loop

  // Add overlayed color bar
    ColorBar();


  // Hit Map
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(4,4);
  for (int islab=0; islab<nslabs; islab++){
    c1->cd(islab+1);
    hitmap[islab]->Draw("col");
  }
  c1->Draw();

  // Shower Profile
  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  TGraph *gr_profile = new TGraph(nslabs,X0s,E_layers);
  gr_profile->Draw("AC*");
  c2->Draw();

  vector<vector<Float_t>> Mean_SD_x;
  vector<vector<Float_t>> Mean_SD_y;
  for (int islab=0; islab<nslabs; islab++){
    vector<Float_t> iMean_SD_x = Mean_SD(islab, layer_hit_x[islab]);
    vector<Float_t> iMean_SD_y = Mean_SD(islab, layer_hit_y[islab]);

    Mean_SD_x.push_back(iMean_SD_x);
    Mean_SD_y.push_back(iMean_SD_y);
  }


  Int_t n_valid_slab = Mean_SD_x.size();
  vector<Int_t> valid_slabs;
  if(n_valid_slab){
    for (int islab=0; islab<n_valid_slab; islab++){

      Bool_t is_nhit   =  10 < Mean_SD_x.at(islab).at(1);
      Bool_t is_sigw_x = ( 0 < Mean_SD_x.at(islab).at(3) ) && ( Mean_SD_x.at(islab).at(3) < 20 );
      Bool_t is_sigw_y = ( 0 < Mean_SD_y.at(islab).at(3) ) && ( Mean_SD_y.at(islab).at(3) < 20 );

      if( is_nhit && is_sigw_x && is_sigw_y ){
        valid_slabs.push_back(Mean_SD_x.at(islab).at(0));
        // cout << "layer " << Mean_SD_x.at(islab).at(0) << " | " << Mean_SD_x.at(islab).at(2) << " " << Mean_SD_x.at(islab).at(3) << " | " << Mean_SD_y.at(islab).at(2) << " " << Mean_SD_y.at(islab).at(3) << " :: " << Mean_SD_x.at(islab).at(1)<< endl;
      }

    }
  }

  // Fit
  TCanvas  *c0 = new TCanvas("c0","c0",800,800);
  TGraph2D *gr = new TGraph2D();
  FindCoG(valid_slabs,gr);

  cout << gr->GetN() << endl;

  ROOT::Fit::Fitter  fitter;
  // make the functor objet
  Fit3D sdist(gr);
  ROOT::Math::Functor fcn(sdist,2);
  // set the function and the initial parameter values
  double pStart[2] = {-50,-100};
  fitter.SetFCN(fcn,pStart);
  // set step sizes different than default ones (0.3 times parameter values)
  for (int i = 0; i < 2; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);

  bool fit_ok = fitter.FitFCN();
  if (!fit_ok) {
    Error("line3Dfit","Line3D Fit failed");
    return false;
  }
  const ROOT::Fit::FitResult & result = fitter.Result();

  std::cout << "Total final distance square " << result.MinFcnValue() << std::endl;
  result.Print(std::cout);

  gr->SetTitle(";x;z;y");
  gr->GetXaxis()->SetRangeUser(-95,95);
  gr->Draw("p0");

  // get fit parameters
  const double * parFit = result.GetParams();

  // draw the fitted line
  int n = 1000;
  double t0 = 0;
  double dt = 210;
  TPolyLine3D *l = new TPolyLine3D(n);
  for (int i = 0; i <n;++i) {
    double t = t0+ dt*i/n;
    double x,y,z;
    line(t,parFit,x,y,z);
    l->SetPoint(i,x,y,z);
  }
  l->SetLineColor(kRed);
  l->Draw("same");

  c0->Draw();

  if(fit_ok && 2 < gr->GetN()){
    BeamAxisLine(parFit);
  }


  // MultiView
  // ProjectView();

  gEve->Redraw3D(kFALSE, kTRUE);

  return kTRUE;
}

void TBDisplay::FindCoG(vector<Int_t> arr, TGraph2D *gr)
{
  Float_t Xcm[nslabs] = {0};
  Float_t Ycm[nslabs] = {0};
  Float_t Zcm[nslabs] = {0};
  Float_t Wsum[nslabs] = {0};

  // Load event data into visualization structures.
  for (int ihit=0; ihit<nhit_len; ihit++){

    if (std::binary_search(arr.begin(), arr.end(), hit_slab[ihit]))
    {
    // CoG based on hit energy
      Xcm[hit_slab[ihit]]  += hit_x[ihit] * hit_energy[ihit];
      Ycm[hit_slab[ihit]]  += hit_y[ihit] * hit_energy[ihit];
      Wsum[hit_slab[ihit]] += hit_energy[ihit];

      Zcm[hit_slab[ihit]]  =  hit_z[ihit];
    }

  } // hit loop

  Int_t N = 0;
  for (int islab = 0; islab < nslabs; islab++)
  {
    if(Xcm[islab]==0 || Ycm[islab]==0) continue;

    Xcm[islab] = Xcm[islab] / Wsum[islab];
    Ycm[islab] = Ycm[islab] / Wsum[islab];

    cout << Xcm[islab] << " | " << Ycm[islab] << endl;

    gr->SetPoint(N,Xcm[islab],Zcm[islab],Ycm[islab]);
    N++;

  } // match slab
}

vector<Float_t> TBDisplay::Mean_SD(int slab, vector<Float_t> arr)
{
  static vector<Float_t> Mean_SD_vec;
  Float_t sum   = 0;
  Int_t nhits = arr.size();
  Mean_SD_vec = {static_cast<float>(slab), static_cast<float>(nhits), -1000.0, -1000.0};

  if (nhits == 0) {
    return Mean_SD_vec;
  }

  for (int i=0; i<nhits; i++){
    sum += arr.at(i);
  }

  Float_t mean  = sum / nhits;
  Float_t sigma = 0;
  for (int i=0; i<nhits; i++){
    sigma += pow(arr.at(i) - mean, 2);
  }

  sigma = sqrt(sigma / nhits);

  Mean_SD_vec[0] = slab;
  Mean_SD_vec[1] = nhits;
  Mean_SD_vec[2] = mean;
  Mean_SD_vec[3] = sigma;

  return Mean_SD_vec;

}

void TBDisplay::BeamAxisLine(const double *p)
{
  double x0,y0,z0;
  double x1,y1,z1;
  line(0,p,x0,z0,y0);
  line(210,p,x1,z1,y1);

  auto BeamLine = new TEveStraightLineSet();
  BeamLine->AddLine( x0,y0,z0,
                     x1,y1,z1);
  BeamLine->SetLineColor(kRed);

  gEve->AddElement(BeamLine);
  gEve->Redraw3D();

}

//______________________________________________________________________________
void TBDisplay::MakeViewerScene(TEveWindowSlot* slot, TEveViewer*& v, TEveScene*& s)
{
  // Create a scene and a viewer in the given slot.

  v = new TEveViewer("Viewer");
  v->SpawnGLViewer(gEve->GetEditor());
  slot->ReplaceWindow(v);
  gEve->GetViewers()->AddElement(v);
  s = gEve->SpawnNewScene("Scene");
  v->AddScene(s);
}

void TBDisplay::LoadHits(TEvePointSet*& ps, int i)
{
  ps = new TEvePointSet("Hit");
  ps->SetOwnIds(kTRUE);
  ps->SetMarkerSize(MARKER_SIZE);
  ps->SetMarkerStyle(54);
  ps->IncDenyDestroy();

  ps->SetNextPoint(hit_x[i],hit_y[i],hit_z[i]);
  // ps->SetMainColor(TColor::GetColorPalette
  //                  (hit_adc_high[i]));
  ps->SetMainColor(TColor::GetColorPalette
                  (hit_energy[i]));
  if(hit_adc_high[i] < hit_energy[i]) ps->SetMainColor(kRed);
  // if(hit_energy[i] > 700) ps->SetMainColor(kRed);
  // ps->SetMainColor(TColor::GetColorPalette
  //                  (hit_energy[i]));
  ps->SetPointId(new TNamed(Form("Point %d", i), ""));
  ps->SetTitle(TString::Format("hit_adc_high=%i\n hit_energy=%f\n hit_isHit=%i\n (%i,%i,%i,%i)",
                                hit_adc_high[i],
                                hit_energy[i],
                                hit_isHit[i],
                                hit_slab[i], hit_chip[i], hit_chan[i], hit_sca[i]));

  gEve->AddElement(ps);
}

void TBDisplay::LoadHits_Box(TEveBoxSet*& bs, int i)
{
  bs = new TEveBoxSet("BoxSet");
  TEveRGBAPalette *pal = new TEveRGBAPalette(0, 10);
  pal->SetupColorArray();
  bs->SetPalette(pal);

  bs->Reset(TEveBoxSet::kBT_AABox, kFALSE, 64);

  bs->AddBox(hit_x[i], hit_y[i], hit_z[i],
              5, 5, 0.5);
  bs->DigitValue(hit_energy[i]);

  bs->SetName(TString::Format("hit_adc_high=%i\n hit_energy=%f\n hit_isHit=%i\n (%i,%i,%i,%i)",
                                hit_adc_high[i],
                                hit_energy[i],
                                hit_isHit[i],
                                hit_slab[i], hit_chip[i], hit_chan[i], hit_sca[i]));


  bs->RefitPlex();

  TEveTrans& t = bs->RefMainTrans();
  t.SetPos(0,0,0);

  // Uncomment these two lines to get internal highlight / selection.
  bs->SetPickable(1);
  bs->SetAlwaysSecSelect(1);

  gEve->AddElement(bs);

}

void TBDisplay::ColorBar()
{
  TEveRGBAPalette *pal = new TEveRGBAPalette(0, 10);
  pal->SetupColorArray();
  auto po = new TEveRGBAPaletteOverlay(pal, 0.55, 0.1, 0.4, 0.05);
  auto v  = gEve->GetDefaultGLViewer();
  v->AddOverlayElement(po);
}

void TBDisplay::MultiDisplay(TEveElement* gentle_geom)
{
  gMultiView = new MultiView;
  gMultiView->f3DView->GetGLViewer()->SetStyle(TGLRnrCtx::kOutline);

  gMultiView->SetDepth(-10);
  gMultiView->ImportGeomRPhi(gentle_geom);
  gMultiView->ImportGeomRhoZ(gentle_geom);
  gMultiView->SetDepth(0);

}

void TBDisplay::ProjectView()
{
  // Fill projected views.
  auto top = gEve->GetCurrentEvent();

  gMultiView->DestroyEventRPhi();
  gMultiView->ImportEventRPhi(top);

  gMultiView->DestroyEventRhoZ();
  gMultiView->ImportEventRhoZ(top);
}

void TBDisplay::Debug(bool debug=false, Long64_t entry=0)
{

  if(!debug) return;
  else if( !(entry % 10000) ) return;

  for (int ihit = 0; ihit < nhit_len; ++ihit)
  {

    cout << hit_sca[ihit] << " " ;

  }

  cout << endl;
  
}
