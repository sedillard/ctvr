#ifndef MAINWINDOW_HPP_INCLUDED
#define MAINWINDOW_HPP_INCLUDED

#include <QMainWindow>

#include "ContourTree.hpp"
#include "ContourTreeVolumeRenderer.hpp"
#include "VolumeRendererWidget.hpp"

class MainWindow : public QMainWindow
{
  Q_OBJECT

  ContourTree *ct;
  ContourTreeVolumeRenderer *ctvr; 
  VolumeRendererWidget *vr;

  public:
    MainWindow( ContourTree *ct, ContourTreeVolumeRenderer *ctvr,
                QWidget *parent=0);

};





#endif
