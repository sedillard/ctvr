#ifndef MAINWINDOW_HPP_INCLUDED
#define MAINWINDOW_HPP_INCLUDED

#include <QMainWindow>
#include <QSlider>
#include <QSpinBox>

#include "ContourTree.hpp"
#include "ContourTreeVolumeRenderer.hpp"
#include "VolumeRendererWidget.hpp"
#include "Color.hpp"
#include "ColorMapEditor.hpp"
#include "TreeWidget.hpp"

class MainWindow : public QMainWindow
{
  Q_OBJECT

  ContourTree *ct;
  ContourTreeVolumeRenderer *ctvr; 
  VolumeRendererWidget *vr;
  ColorMapEditor *cme;
  TreeWidget *tw;
  QSlider *isoval;
  QSpinBox *smooth_amt;

  public:
    MainWindow( ContourTree *ct, ContourTreeVolumeRenderer *ctvr,
                QWidget *parent=0);

  public Q_SLOTS:
    void tf_edited( RGBA8* ); 
    void smooth_amt_changed();
    void isoval_changed();
};





#endif
