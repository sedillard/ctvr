#include "MainWindow.hpp"

#include <string>
#include <cstdlib>
#include <fstream>

#include <QApplication>

using namespace std;


int main( int argc, char* argv[] )
{
  QApplication app( argc, argv );

  string filename  = argv[1];
  uint32_t ncols   = atoi(argv[2]); 
  uint32_t nrows   = atoi(argv[3]); 
  uint32_t nstacks = atoi(argv[4]); 
  uint32_t nvoxels = ncols*nrows*nstacks;

  uint8_t *image = new uint8_t[nvoxels];
  {
    ifstream file( filename.c_str(), ios::in | ios::binary );
    file.read( reinterpret_cast<char*>(image), nvoxels );
    size_t nread = file.gcount();
    if ( nread != nvoxels ) {
      cerr << "file was shorter than expected" << endl; 
      app.exit();
    }
  }

  cout << "building contour tree" << endl;
  ContourTree *ct = new ContourTree(image,ncols,nrows,nstacks); 
  ct->build();

  cout << "setting up volume renderer" << endl;
  ContourTreeVolumeRenderer *ctvr = 
    new ContourTreeVolumeRenderer(*ct,image,ncols,nrows,nstacks);

  MainWindow *main_window = new MainWindow(ct,ctvr);

  main_window->show();
  app.exec();
}


MainWindow::MainWindow
( 
  ContourTree *ct_, 
  ContourTreeVolumeRenderer *ctvr_ 
) :
  ct(ct_),
  ctvr(ctvr_)
{
  vr = new VolumeRendererWidget(ctvr,this);
  setCentralWidget(vr);
}


