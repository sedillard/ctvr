#ifndef MAINWINDOW_HPP_INCLUDED
#define MAINWINDOW_HPP_INCLUDED

#include <QMainWindow>
#include <QSlider>
#include <QMenu>
#include <QThread>
#include <QAction>

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
    MainWindow( ContourTree *ct, ContourTreeVolumeRenderer *ctvr);

};





#endif
