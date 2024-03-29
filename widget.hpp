/*
Copyright (C) 2011 Georgia Institute of Technology

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <QCheckBox>
#include <QComboBox>
#include <QSpinBox>
#include <QTreeWidget>

#include <hdf5.h>
#include <hdf5_hl.h>
#include <qwt_plot_curve.h>
#include <rtxi/dsp/dolph.h>
#include <rtxi/dsp/fft.h>
#include <rtxi/dsp/gen_win.h>
#include <rtxi/dsp/hamming.h>
#include <rtxi/dsp/hann.h>
#include <rtxi/dsp/kaiser.h>
#include <rtxi/dsp/rectnglr.h>
#include <rtxi/dsp/trianglr.h>
#include <rtxi/plot/basicplot.h>
#include <rtxi/widgets.hpp>

// This is an generated header file. You may change the namespace, but
// make sure to do the same in implementation (.cpp) file
namespace analysis_module
{

constexpr std::string_view MODULE_NAME = "analysis-module";

// This format is taken from data recorder header file.
typedef struct data_sample
{
  int64_t time;
  double value;
} data_sample;

inline std::vector<Widgets::Variable::Info> get_default_vars()
{
  return {};
}

inline std::vector<IO::channel_t> get_default_channels()
{
  return {};
}

herr_t op_func(hid_t loc_id,
               const char* name,
               const H5O_info_t* info,
               void* operator_data);

class Panel : public Widgets::Panel
{
  Q_OBJECT
public:
  Panel(const Panel&) = delete;
  Panel(Panel&&) = delete;
  Panel& operator=(const Panel&) = delete;
  Panel& operator=(Panel&&) = delete;
  Panel(QMainWindow* main_window, Event::Manager* ev_manager);
  ~Panel() override;

  void customizeGUI();

  enum window_t
  {
    RECT = 0,
    TRI,
    HAMM,
    HANN,
    CHEBY,
    KAISER
  };
  enum plot_t
  {
    TIMESERIES,
    SCATTER,
    FFT
  };

private:
  // inputs, states, flags, calculated values
  bool fwrChecked=false;
  GenericWindow* disc_window=nullptr;
  plot_t plot_mode;
  window_t window_shape;
  double Kalpha;  // Kaiser window sidelobe attenuation parameter
  double Calpha;  // Chebyshev window sidelobe attenuation parameter

  // File I/O
  hid_t file_id;
  hid_t dataset_id;

  // GUI components
  BasicPlot* omniplot = nullptr;
  QwtPlotCurve* tscurve = nullptr;
  QwtPlotCurve* fftcurve = nullptr;
  double xmin, xmax, ymin, ymax;
  // QwtText xAxisTitle, yAxisTitle;

  QGroupBox* plotControls = nullptr;
  QGroupBox* plotOptions = nullptr;
  QLineEdit* fileNameEdit = nullptr;
  QPushButton* resetPlotButton = nullptr;
  QPushButton* savePlotButton = nullptr;
  QPushButton* exportSeriesButton = nullptr;
  QComboBox* windowShape = nullptr;
  QComboBox* plotType = nullptr;
  QButtonGroup* plotOptionsButtons = nullptr;
  QCheckBox* FWRCheckBox = nullptr;
  QDoubleSpinBox* kalphaEdit = nullptr;
  QDoubleSpinBox* calphaEdit = nullptr;

  // custom functions
  void initParameters();
  int openFile(QString& filename);
  void closeFile();
  void makeWindow(int);

  // Tree traversal variables for hdf5 datasets
  QTreeWidget* treeViewer = nullptr;
  QTreeWidgetItem* treeParent = nullptr;
  std::vector<data_sample> data_buffer;
  std::vector<double> channel_data;
  std::vector<double> time_buffer;
  std::vector<double> fft_output_y;
  std::vector<double> fft_output_x;
  std::vector<complex> fft_input;
  std::vector<complex> fft_buffer;

private slots:
  void getTrialData();
  void changeChannel(QModelIndex);
  void updatePlotMode(int);
  void updatePlot();
  void screenshot();
  void toggleFWR(bool);
  void resetAxes();
  void changeDataFile();
  void updateWindow(int);
  void updateKalpha(double);
  void updateCalpha(double);
  void exportData();

signals:
  void setPlotRange(double, double, double, double);

  // Any functions and data related to the GUI are to be placed here
};

class Plugin : public Widgets::Plugin
{
public:
  explicit Plugin(Event::Manager* ev_manager);
};

}  // namespace analysis_module
