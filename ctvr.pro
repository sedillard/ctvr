
TEMPLATE = app

QT += opengl
CONFIG += warn_on build_all debug_and_release  no_keywords
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
  ColorMapEditor.hpp \
  ContourTree.hpp \
  ContourTreeVolumeRenderer.hpp \
  HexMeshRayTracer.hpp \
  MainWindow.hpp \
  TreeWidget.hpp \
  VolumeRendererWidget.hpp \

SOURCES	= \
  ColorMapEditor.cpp \
  ContourTree.cpp \
  ContourTreeVolumeRenderer.cpp \
  MainWindow.cpp \
  TreeWidget.cpp \
  VolumeRendererWidget.cpp \




