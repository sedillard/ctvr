#include "MainWindow.hpp"

#include <string>
#include <cstdlib>
#include <fstream>

#include <QApplication>
#include <QFile>
#include <QSplitter>
#include <QToolBar>
#include <QSpinBox>
#include <QSlider>
#include <QVBoxLayout>

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
    ctvr->read_events_file(argv[5]);
  }

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

  QSplitter *split = new QSplitter(this);
  split->setOrientation(Qt::Vertical);
  setCentralWidget(split);
  
  vr = new VolumeRendererWidget(ctvr,split);

  QWidget *controls = new QWidget(split);
  QVBoxLayout *controls_layout = new QVBoxLayout;
  controls->setLayout(controls_layout);
  
  cme = new ColorMapEditor(controls);
  isoval = new QSlider(Qt::Horizontal);
  isoval->setRange(-1,256);

  controls_layout->addWidget(cme);
  controls_layout->addWidget(isoval);

  //tw = new TreeWidget(ctvr);
  //tw->show();

  split->setStretchFactor(0,2);
  split->setStretchFactor(1,1);

  QToolBar *toolbar = addToolBar(tr("Stuff"));

  smooth_amt = new QSpinBox;
  toolbar->addWidget(smooth_amt);
  smooth_amt->setRange(0,4000000);
  smooth_amt->setValue(0);
  
  connect( smooth_amt, SIGNAL(editingFinished()), this, SLOT(smooth_amt_changed()) );
  connect( cme, SIGNAL(edited(RGBA8*)), this, SLOT(tf_edited(RGBA8*)) );
  connect( isoval, SIGNAL(sliderReleased()), this, SLOT(isoval_changed()));
  
}

void MainWindow::tf_edited( RGBA8 *colors )
{
  vr->makeCurrent();
  ctvr->set_global_tf( colors );
  ctvr->update_global_tf_tex();
  vr->update(); 
}

void MainWindow::smooth_amt_changed()
{
  int t = smooth_amt->value();
  ctvr->do_merge_events(t); 
  ctvr->cluster_tf();
  vr->makeCurrent();
  ctvr->update_tf_tex();
  vr->update();
  //tw->update();
}

void MainWindow::isoval_changed()
{
  cout << "isovalue = " << this->isoval->value() << endl;
  ctvr->isoval = this->isoval->value() / 255.0;
  vr->update();
}
