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
#include "../include/MultiView.hh"
MultiView* gMultiView = 0;

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

   cout << endl;
   cout << "Going to " << ev << "..." << endl;
   cout << endl;

   TGraph2D *gr = new TGraph2D();
   fCurEv = ev;

   Long64_t ientry = LoadTree( evlist->GetEntry(ev) );
   if (ientry == 0) {
      Warning("GotoEvent", "Entry is empty");
      return kFALSE;
   }

   nb = fChain->GetEntry(evlist->GetEntry(ev));

   float Xcm[nslabs] = {0};
   float Ycm[nslabs] = {0};
   float Zcm[nslabs] = {0};
   float Wsum[nslabs] = {0};

   // Load event data into visualization structures.
   for (int ihit=0; ihit<nhit_len; ihit++){
      
      // LoadHits(fHits,ihit);
      LoadHits_Box(fHits_Box,ihit);

      for (int islab = 0; islab < nslabs; islab++)
      {
         if (hit_slab[ihit]==islab)
         {
            // Xcm[islab]  += hit_x[ihit] * hit_energy[ihit];
            // Ycm[islab]  += hit_y[ihit] * hit_energy[ihit];
            Xcm[islab]  += hit_x[ihit] * hit_adc_high[ihit];
            Ycm[islab]  += hit_y[ihit] * hit_adc_high[ihit];
            Zcm[islab]  =  hit_z[ihit];
            Wsum[islab] += hit_adc_high[ihit];
         }

      } // match slab

   } // hit loop

   int N = 0;
   for (int islab = 0; islab < nslabs; islab++)
   {
      if(Xcm[islab]==0 || Ycm[islab]==0) continue;

      Xcm[islab] = Xcm[islab] / Wsum[islab];
      Ycm[islab] = Ycm[islab] / Wsum[islab];

      gr->SetPoint(N,Xcm[islab],Ycm[islab],Zcm[islab]);
      N++;

   } // match slab

   // Add overlayed color bar
   ColorBar();

   gEve->Redraw3D(kFALSE, kTRUE);

   return kTRUE;
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

void TBDisplay::Display()
{
   Long64_t nbytes = 0, nb = 0;
   int ievent = 0;
   int cnt_event = 0;

   cout << "Done.\n";

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
