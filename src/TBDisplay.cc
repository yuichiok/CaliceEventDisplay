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
  // TGraph2D *gr_profile = new TGraph2D();
  Float_t X0s[nslabs] = {1.198630137, 2.397260274, 3.595890411, 4.794520548, 5.993150685, 7.191780822, 8.390410959, 9.589041096, 10.78767123, 12.38584475, 13.98401826, 15.58219178, 17.1803653, 18.77853881, 20.37671233};
  Float_t E_layers[nslabs] = {0};

  // 2D Hit Maps
  TH2F *hitmap[nslabs];
  for (int islab=0; islab<nslabs; islab++){
    TString h2name = TString::Format("hitmap layer%i",islab);
    hitmap[islab] = new TH2F(h2name,h2name,32,-90,90,32,-90,90);
  }

  // Load event data into visualization structures.
  for (int ihit=0; ihit<nhit_len; ihit++){
    
    E_layers[hit_slab[ihit]] += hit_energy[ihit];

    hitmap[hit_slab[ihit]]->Fill(hit_x[ihit], hit_y[ihit], hit_energy[ihit]);

    LoadHits_Box(fHits_Box,ihit);

  } // hit loop

  // Add overlayed color bar
    ColorBar();

  // Fit
  TCanvas  *c0 = new TCanvas("c0","c0",800,800);
  TGraph2D *gr = new TGraph2D();
  FindCoG(gr);
  gr->SetTitle(";x;z;y");
  gr->GetXaxis()->SetRangeUser(-95,95);
  gr->Draw("p0");
  c0->Draw();

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





  // MultiView
  // ProjectView();

  gEve->Redraw3D(kFALSE, kTRUE);

  return kTRUE;
}

void TBDisplay::FindCoG(TGraph2D *gr)
{
  Float_t Xcm[nslabs] = {0};
  Float_t Ycm[nslabs] = {0};
  Float_t Zcm[nslabs] = {0};
  Float_t Wsum[nslabs] = {0};

  // Load event data into visualization structures.
  for (int ihit=0; ihit<nhit_len; ihit++){
    
    for (int islab = 0; islab < nslabs; islab++)
    {
        if (hit_slab[ihit]==islab)
        {
        // CoG based on hit energy
          Xcm[islab]  += hit_x[ihit] * hit_energy[ihit];
          Ycm[islab]  += hit_y[ihit] * hit_energy[ihit];
          Wsum[islab] += hit_energy[ihit];
        // CoG based on adc high
          // Xcm[islab]  += hit_x[ihit] * hit_adc_high[ihit];
          // Ycm[islab]  += hit_y[ihit] * hit_adc_high[ihit];
          // Wsum[islab] += hit_adc_high[ihit];

          Zcm[islab]  =  hit_z[ihit];
        }

    } // match slab

  } // hit loop

  Int_t N = 0;
  for (int islab = 0; islab < nslabs; islab++)
  {
    if(Xcm[islab]==0 || Ycm[islab]==0) continue;

    Xcm[islab] = Xcm[islab] / Wsum[islab];
    Ycm[islab] = Ycm[islab] / Wsum[islab];

    gr->SetPoint(islab,Xcm[islab],Zcm[islab],Ycm[islab]);

  } // match slab
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
