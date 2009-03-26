#ifndef MAINWINDOW_HPP_INCLUDED
#define MAINWINDOW_HPP_INCLUDED

#include <QMainWindow>

#include "ContourTree.hpp"
#include "ContourTreeVolumeRenderer.hpp"
#include "VolumeRendererWidget.hpp"
#include "Color.hpp"
#include "ColorMapEditor.hpp"

class MainWindow : public QMainWindow
{
  Q_OBJECT

  ContourTree *ct;
  ContourTreeVolumeRenderer *ctvr; 
  VolumeRendererWidget *vr;
  ColorMapEditor *cme;

  public:
    MainWindow( ContourTree *ct, ContourTreeVolumeRenderer *ctvr,
                QWidget *parent=0);

  public slots:
    void tf_edited( RGBA8* ); 
    void lod_changed(int);
};





#endif
