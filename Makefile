PLUGIN_NAME = analysis_tools

RTXI_INCLUDES=/usr/local/lib/rtxi_includes

HEADERS = analysis-tools.h\
          $(RTXI_INCLUDES)/scatterplot.h\
          $(RTXI_INCLUDES)/incrementalplot.h\
          $(RTXI_INCLUDES)/basicplot.h\
          $(RTXI_INCLUDES)/scrollzoomer.h\
		  $(RTXI_INCLUDES)/runningstat.h\
		  $(RTXI_INCLUDES)/scrollbar.h\
		  $(RTXI_INCLUDES)/DSP/gen_win.h\
		  $(RTXI_INCLUDES)/DSP/rectnglr.h\
		  $(RTXI_INCLUDES)/DSP/trianglr.h\
		  $(RTXI_INCLUDES)/DSP/hamming.h\
		  $(RTXI_INCLUDES)/DSP/hann.h\
		  $(RTXI_INCLUDES)/DSP/dolph.h\
		  $(RTXI_INCLUDES)/DSP/kaiser.h\

SOURCES = analysis-tools.cpp\
		  moc_analysis-tools.cpp\
		  $(RTXI_INCLUDES)/scatterplot.cpp\
		  $(RTXI_INCLUDES)/incrementalplot.cpp\
		  $(RTXI_INCLUDES)/basicplot.cpp\
		  $(RTXI_INCLUDES)/scrollzoomer.cpp\
		  $(RTXI_INCLUDES)/runningstat.cpp\
		  $(RTXI_INCLUDES)/scrollbar.cpp\
		  $(RTXI_INCLUDES)/moc_scatterplot.cpp\
		  $(RTXI_INCLUDES)/moc_incrementalplot.cpp\
		  $(RTXI_INCLUDES)/moc_basicplot.cpp\
		  $(RTXI_INCLUDES)/moc_scrollzoomer.cpp\
		  $(RTXI_INCLUDES)/moc_runningstat.cpp\
		  $(RTXI_INCLUDES)/moc_scrollbar.cpp\
		  
			 
LIBS = -lqwt -lgsl -lhdf5

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile
