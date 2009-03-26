
TEMPLATE = app

QT += opengl
CONFIG += warn_on build_all debug_and_release
CONFIG(debug, debug|release) {
  TARGET = ctvr_dbg
  MOC_DIR = debug
  OBJECTS_DIR = debug
} else {
   TARGET = ctvr
   MOC_DIR = release
   OBJECTS_DIR = release
}

LIBS += -fopenmp
QMAKE_CXXFLAGS += -fopenmp

RESOURCES = ctvr.qrc

HEADERS	= \
  HexMeshRayTracer.hpp \
  ColorMapEditor.hpp \
  ContourTree.hpp \
  ContourTreeVolumeRenderer.hpp \
  MainWindow.hpp \
  VolumeRendererWidget.hpp \

SOURCES	= \
  ColorMapEditor.cpp \
  ContourTree.cpp \
  ContourTreeVolumeRenderer.cpp \
  MainWindow.cpp \
  VolumeRendererWidget.cpp \




