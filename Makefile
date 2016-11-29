PLUGIN_NAME = analysis_tools

HEADERS = analysis-tools.h

SOURCES = analysis-tools.cpp\
          moc_analysis-tools.cpp

LIBS = -lgsl -lhdf5 -lrtplot -lrtdsp

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile
