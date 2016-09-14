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
          $(RTXI_INCLUDES)/DSP/fft.h\
          $(RTXI_INCLUDES)/DSP/dit_sino.h\
          $(RTXI_INCLUDES)/DSP/complex.h\
          $(RTXI_INCLUDES)/DSP/cbitrev.h\
          $(RTXI_INCLUDES)/DSP/log2.h\
          $(RTXI_INCLUDES)/DSP/misdefs.h\
          $(RTXI_INCLUDES)/DSP/acosh.h\

SOURCES = analysis-tools.cpp\
          moc_analysis-tools.cpp\
          $(RTXI_INCLUDES)/scatterplot.cpp\
          $(RTXI_INCLUDES)/incrementalplot.cpp\
          $(RTXI_INCLUDES)/basicplot.cpp\
          $(RTXI_INCLUDES)/scrollzoomer.cpp\
          $(RTXI_INCLUDES)/runningstat.cpp\
          $(RTXI_INCLUDES)/scrollbar.cpp\
          $(RTXI_INCLUDES)/DSP/gen_win.cpp\
          $(RTXI_INCLUDES)/DSP/rectnglr.cpp\
          $(RTXI_INCLUDES)/DSP/trianglr.cpp\
          $(RTXI_INCLUDES)/DSP/hamming.cpp\
          $(RTXI_INCLUDES)/DSP/hann.cpp\
          $(RTXI_INCLUDES)/DSP/dolph.cpp\
          $(RTXI_INCLUDES)/DSP/kaiser.cpp\
          $(RTXI_INCLUDES)/DSP/fft.cpp\
          $(RTXI_INCLUDES)/DSP/dit_sino.cpp\
          $(RTXI_INCLUDES)/DSP/complex.cpp\
          $(RTXI_INCLUDES)/DSP/cbitrev.cpp\
          $(RTXI_INCLUDES)/DSP/log2.cpp\
          $(RTXI_INCLUDES)/DSP/acosh.cpp\
          $(RTXI_INCLUDES)/moc_scatterplot.cpp\
          $(RTXI_INCLUDES)/moc_incrementalplot.cpp\
          $(RTXI_INCLUDES)/moc_basicplot.cpp\
          $(RTXI_INCLUDES)/moc_scrollzoomer.cpp\
          $(RTXI_INCLUDES)/moc_runningstat.cpp\
          $(RTXI_INCLUDES)/moc_scrollbar.cpp\
          $(RTXI_INCLUDES)/DSP/moc_gen_win.cpp\
          $(RTXI_INCLUDES)/DSP/moc_rectnglr.cpp\
          $(RTXI_INCLUDES)/DSP/moc_trianglr.cpp\
          $(RTXI_INCLUDES)/DSP/moc_hamming.cpp\
          $(RTXI_INCLUDES)/DSP/moc_hann.cpp\
          $(RTXI_INCLUDES)/DSP/moc_dolph.cpp\
          $(RTXI_INCLUDES)/DSP/moc_kaiser.cpp\
          $(RTXI_INCLUDES)/DSP/moc_fft.cpp\
          $(RTXI_INCLUDES)/DSP/moc_dit_sino.cpp\
          $(RTXI_INCLUDES)/DSP/moc_complex.cpp\
          $(RTXI_INCLUDES)/DSP/moc_cbitrev.cpp\
          $(RTXI_INCLUDES)/DSP/moc_log2.cpp\
          $(RTXI_INCLUDES)/DSP/moc_acosh.cpp\

LIBS = -lqwt -lgsl -lhdf5

### Do not edit below this line ###

include $(shell rtxi_plugin_config --pkgdata-dir)/Makefile.plugin_compile
