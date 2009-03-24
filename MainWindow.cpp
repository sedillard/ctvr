#include "MainWindow.hpp"

#include <string>
#include <cstdlib>
#include <fstream>

#include <QApplication>
#include <QFile>

using namespace std;


int main( int argc, char* argv[] )
{
  QApplication app( argc, argv );

  if ( argc < 5 ) {
    cerr << "usage: ctvr volume_file width height depth" << endl;
    app.exit();
  }

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

  QFile file(":/RayCaster.glsl");
  file.open(QIODevice::ReadOnly | QIODevice::Text );
  ctvr->fshader_src = string(file.readAll().data());

  MainWindow *main_window = new MainWindow(ct,ctvr);

  main_window->show();
  main_window->resize(700,600);
  app.exec();
}


MainWindow::MainWindow
( 
  ContourTree *ct_, 
  ContourTreeVolumeRenderer *ctvr_ ,
  QWidget *parent
) :
  QMainWindow(parent),
  ct(ct_),
  ctvr(ctvr_)
{
  vr = new VolumeRendererWidget(ctvr,this);
  setCentralWidget(vr);
}


