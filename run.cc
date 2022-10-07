#include "TROOT.h"
#include "TFile.h"
#include "TApplication.h"

#include "src/TBDisplay.cc"
#include "include/MultiView.hh"
#include "src/FileSelector.cc"

TBDisplay *gDisplay = 0;

void make_gui();

void run(TString setting = "reco e- 10"){

   FileSelector fs(setting);
   TString filein = fs.GetRunName_with_path();

	cout << "Input: " << filein << endl; 

	gDisplay = new TBDisplay(filein);

	TFile::SetCacheFileDir(".");

	TEveManager::Create();

	const int nslabs = 15;
	TEveGeoShape *layer[nslabs];
	for(int i=0; i<nslabs; i++){
		layer[i] = new TEveGeoShape;
		layer[i]->SetShape(new TGeoBBox(0.95e+02, 0.95e+02, 0.5)); // dx, dy, dz
		layer[i]->SetNSegments(100); // number of vertices
		layer[i]->SetMainColor(kGreen);
		layer[i]->SetMainAlpha(0.2);
		layer[i]->RefMainTrans().SetPos(0, 0, 0.5+i*15.0); // set position
		gEve->AddGlobalElement(layer[i]);
	}

	gStyle->SetOptStat(0);

	// gDisplay->MultiDisplay(layer[0]);

   // Set GLViewer Style
   //=====================

	TGLViewer *tglv = gEve->GetDefaultGLViewer();
	tglv->SetGuideState(TGLUtil::kAxesEdge, kTRUE, kFALSE, 0);
	tglv->SetStyle(TGLRnrCtx::kOutline);

	gEve->GetBrowser()->GetTabRight()->SetTab(1);

	make_gui();

	gEve->AddEvent(new TEveEventManager("Event", "SiWECAL VSD Event"));

	gDisplay->GotoEvent(0);

	gEve->Redraw3D(kTRUE); // Reset camera after the first event has been shown.


	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 200 );

}

//______________________________________________________________________________
void make_gui()
{
   // Create minimal GUI for event navigation.

   auto browser = gEve->GetBrowser();
   browser->StartEmbedding(TRootBrowser::kLeft);

   auto frmMain = new TGMainFrame(gClient->GetRoot(), 1000, 600);
   frmMain->SetWindowName("XX GUI");
   frmMain->SetCleanup(kDeepCleanup);

   auto hf = new TGHorizontalFrame(frmMain);
   {
      TString icondir("./icons/");
      TGPictureButton* b = 0;

      b = new TGPictureButton(hf, gClient->GetPicture(icondir+"back.png"));
      hf->AddFrame(b);
      b->Connect("Clicked()", "TBDisplay", gDisplay, "Prev()");

      b = new TGPictureButton(hf, gClient->GetPicture(icondir+"search.png"));
      hf->AddFrame(b);
      b->Connect("Clicked()", "TBDisplay", gDisplay, "GoTo()");

      b = new TGPictureButton(hf, gClient->GetPicture(icondir+"next.png"));
      hf->AddFrame(b);
      b->Connect("Clicked()", "TBDisplay", gDisplay, "Next()");
   }
   frmMain->AddFrame(hf);

   frmMain->MapSubwindows();
   frmMain->Resize();
   frmMain->MapWindow();

   browser->StopEmbedding();
   browser->SetTabTitle("Event Control", 0);
}
