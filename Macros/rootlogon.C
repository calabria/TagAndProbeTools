{

 TString roofitbase = getenv("ROOFITSYS");

 cout<<"Loading RooFit libraries..."<<endl;

 cout<<roofitbase + "lib/libRooFit.so"<<endl;
 cout<<roofitbase + "lib/libRooFitCore.so"<<endl;
 cout<<"-I" + roofitbase + "include"<<endl;

  if (roofitbase.Length() > 0) {
	gSystem->Load(roofitbase + "lib/libRooFit.so");
	gSystem->Load(roofitbase + "lib/libRooFitCore.so");
	using namespace RooFit;
	gSystem->AddIncludePath("-I" + roofitbase + "include");
  }

 cout<<"RooFit libraries loaded."<<endl;

}
