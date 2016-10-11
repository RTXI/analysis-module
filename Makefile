PLUGIN_NAME = analysis_tools

RTXI_INCLUDES=

HEADERS = analysis-tools.h

SOURCES = analysis-tools.cpp\
          moc_analysis-tools.cpp

LIBS = -lqwt -lgsl -lhdf5 -lrtplot -lrtdsp -lrtmath

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile
