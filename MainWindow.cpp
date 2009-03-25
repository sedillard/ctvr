#include "MainWindow.hpp"

#include <string>
#include <cstdlib>
#include <fstream>

#include <QApplication>
#include <QFile>

#include "ColorMapEditor.hpp"

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

  if ( argc > 5 ) {
    ctvr->read_tracker_file(argv[5]);
    ctvr->cluster_tf();
  }

  MainWindow *main_window = new MainWindow(ct,ctvr);

  main_window->show();
  main_window->resize(700,600);

  ColorMapEditor *color_map_ed = new ColorMapEditor;
  color_map_ed->show();
  app.connect( color_map_ed, SIGNAL(edited(RGBA8*)), main_window, SLOT(tf_edited(RGBA8*)) );

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

void MainWindow::tf_edited( RGBA8 *colors )
{
  vr->makeCurrent();
  ctvr->set_global_tf( colors );
  ctvr->update_global_tf_tex();
  vr->update(); 
}
