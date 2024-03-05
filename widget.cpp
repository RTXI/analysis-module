/*
Copyright (C) 2015 Georgia Institute of Technology

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

#include <QButtonGroup>
#include <QFileDialog>
#include <QMessageBox>
#include <QSettings>
#include <QTimer>
#include <cmath>

#include "widget.hpp"

#include <qnamespace.h>

analysis_module::Plugin::Plugin(Event::Manager* ev_manager)
    : Widgets::Plugin(ev_manager, std::string(analysis_module::MODULE_NAME))
{
}

analysis_module::Panel::Panel(QMainWindow* main_window,
                              Event::Manager* ev_manager)
    : Widgets::Panel(
        std::string(analysis_module::MODULE_NAME), main_window, ev_manager)
{
  setWhatsThis(
      "<p><b>Analysis Tools</b></p><p>View HDF data recorded by the data "
      "recorded. Simply open the file in the File Control Menu and click "
      "specific channels in the HDF5 Viewer list. Under 'Plotting Options,' "
      "you can change from view from time series to an FFT. To zoom in to the "
      "plot, use the mouse to click-and-drag an area. Plots can be saved to "
      "PDF (via the Plot button), and data used in the plot (both x and y "
      "axis) can be saved to a text file (via the Export button)</p>");
  initParameters();
  customizeGUI();
  QTimer::singleShot(0, this, SLOT(resizeMe()));
}

analysis_module::Panel::~Panel()
{
  if (file_id != 0) {
    H5Fclose(file_id);
    file_id = 0;
  }
}

void analysis_module::Panel::initParameters()
{
  fwrChecked = false;
  file_id = 0;
  dataset_id = 0;
  plot_mode = TIMESERIES;
  window_shape = RECT;
  Kalpha = 1.5;
  Calpha = 70;
  data_buffer.clear();
  channel_data.clear();
  time_buffer.clear();
  fft_output_y.clear();
  fft_output_x.clear();
  fft_input.clear();
  fft_buffer.clear();
}

void analysis_module::Panel::customizeGUI()
{
  auto* customlayout = new QHBoxLayout(this);

  // File control
  auto* fileColumnLayout = new QVBoxLayout;
  auto* fileBox = new QGroupBox(tr("File Control"));
  auto* fileLayout = new QHBoxLayout;
  fileBox->setLayout(fileLayout);
  fileNameEdit = new QLineEdit;
  fileNameEdit->setReadOnly(true);
  fileLayout->addWidget(fileNameEdit);
  auto* fileChangeButton = new QPushButton("Open");
  fileLayout->addWidget(fileChangeButton);
  QObject::connect(
      fileChangeButton, SIGNAL(released()), this, SLOT(changeDataFile()));
  fileColumnLayout->addWidget(fileBox);

  // Plot controls
  auto* plotColumnLayout = new QVBoxLayout;
  plotControls = new QGroupBox("Plot Controls");
  auto* plotControlsLayout = new QHBoxLayout;
  plotControls->setLayout(plotControlsLayout);
  resetPlotButton = new QPushButton("Reset");
  QObject::connect(resetPlotButton, SIGNAL(clicked()), this, SLOT(resetAxes()));
  savePlotButton = new QPushButton("Save");
  QObject::connect(savePlotButton, SIGNAL(clicked()), this, SLOT(screenshot()));
  savePlotButton->setToolTip("Save screenshot of the plot");
  exportSeriesButton = new QPushButton("Export");
  QObject::connect(
      exportSeriesButton, SIGNAL(clicked()), this, SLOT(exportData()));
  plotControlsLayout->addWidget(resetPlotButton);
  plotControlsLayout->addSpacerItem(
      new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum));
  plotControlsLayout->addWidget(savePlotButton);
  plotControlsLayout->addWidget(exportSeriesButton);
  plotColumnLayout->addWidget(plotControls);
  plotControls->setEnabled(false);

  // Put plot under plot controls
  omniplot = new BasicPlot(this);
  /*
   * xAxisTitle.setText("X Axis");
   * yAxisTitle.setText("Y Axis");
   * omniplot->setAxisTitle(QwtPlot::xBottom, xAxisTitle);
   * omniplot->setAxisTitle(QwtPlot::yLeft, yAxisTitle);
   */
  plotColumnLayout->addWidget(omniplot);
  tscurve = new QwtPlotCurve;
  tscurve->setPen(Qt::white);
  fftcurve = new QwtPlotCurve;
  fftcurve->setPen(Qt::white);

  // Plot options / file controls
  plotOptions = new QGroupBox(tr("Plotting Options"));
  auto* plotOptionsLayout = new QGridLayout;
  plotOptions->setLayout(plotOptionsLayout);

  auto* plotTypeLabel = new QLabel("Plot Type");
  plotType = new QComboBox;
  plotType->insertItem(1, "Time Series");
  plotType->insertItem(2, "Scatter");
  plotType->insertItem(3, "FFT");
  QObject::connect(plotType,
                   SIGNAL(currentIndexChanged(int)),
                   this,
                   SLOT(updatePlotMode(int)));
  plotOptionsLayout->addWidget(plotTypeLabel, 1, 0);
  plotOptionsLayout->addWidget(plotType, 1, 1);

  auto* plotOptionsVerticalLayout = new QVBoxLayout;
  plotOptionsButtons = new QButtonGroup;
  plotOptionsButtons->setExclusive(false);
  FWRCheckBox = new QCheckBox("TS: Full wave rectify");
  plotOptionsVerticalLayout->addWidget(FWRCheckBox);
  plotOptionsButtons->addButton(FWRCheckBox);
  FWRCheckBox->setChecked(false);
  QObject::connect(
      FWRCheckBox, SIGNAL(toggled(bool)), this, SLOT(toggleFWR(bool)));
  FWRCheckBox->setToolTip("Enable full wave rectification of time series plot");
  plotOptionsLayout->addLayout(plotOptionsVerticalLayout, 2, 0, 1, 2);

  auto* windowLabel = new QLabel("FFT window shape:");
  windowShape = new QComboBox;
  windowShape->insertItem(1, "Rectangular");
  windowShape->insertItem(2, "Triangular (Bartlett)");
  windowShape->insertItem(3, "Hamming");
  windowShape->insertItem(4, "Hann");
  windowShape->insertItem(5, "Chebyshev");
  windowShape->insertItem(6, "Kaiser");
  QObject::connect(
      windowShape, SIGNAL(activated(int)), this, SLOT(updateWindow(int)));
  windowShape->setToolTip(
      "Choose a window to apply for the FFT plot. For no window, choose "
      "Rectangular.");
  plotOptionsLayout->addWidget(windowLabel, 3, 0);
  plotOptionsLayout->addWidget(windowShape, 3, 1);

  auto* kalphaLabel = new QLabel("Kaiser Alpha");
  plotOptionsLayout->addWidget(kalphaLabel, 4, 0);
  auto* kalphaEdit = new QDoubleSpinBox(plotOptions);
  kalphaEdit->setValue(Kalpha);
  QObject::connect(kalphaEdit,
                   SIGNAL(valueChanged(double)),
                   this,
                   SLOT(updateKalpha(double)));
  kalphaEdit->setToolTip("Attenuation parameter for Kaiser window");
  plotOptionsLayout->addWidget(kalphaEdit, 4, 1);

  auto* calphaLabel = new QLabel("Chebyshev (dB)");
  plotOptionsLayout->addWidget(calphaLabel, 5, 0);
  calphaEdit = new QDoubleSpinBox(plotOptions);
  calphaEdit->setValue(Calpha);
  QObject::connect(calphaEdit,
                   SIGNAL(valueChanged(double)),
                   this,
                   SLOT(updateCalpha(double)));
  calphaEdit->setToolTip("Attenuation parameter for Chebyshev window");
  plotOptionsLayout->addWidget(calphaEdit, 5, 1);
  fileColumnLayout->addWidget(plotOptions);
  plotOptions->setEnabled(false);

  // HDF5 viewer
  treeViewer = new QTreeWidget;
  treeViewer->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
  treeViewer->setHeaderLabels(QStringList("HDF Tree"));
  QObject::connect(treeViewer,
                   SIGNAL(clicked(QModelIndex)),
                   this,
                   SLOT(changeChannel(QModelIndex)));
  fileColumnLayout->addWidget(treeViewer);

  customlayout->addLayout(fileColumnLayout);
  customlayout->addLayout(plotColumnLayout);

  setLayout(customlayout);
  QObject::connect(this,
                   SIGNAL(setPlotRange(double, double, double, double)),
                   omniplot,
                   SLOT(setAxes(double, double, double, double)));
}

void analysis_module::Panel::exportData()
{
  auto* fd = new QFileDialog(this);
  fd->setFileMode(QFileDialog::AnyFile);
  fd->setViewMode(QFileDialog::Detail);
  fd->setDefaultSuffix("txt");

  QString fileName =
      QFileDialog::getSaveFileName(this,
                                   "Save the plotted data",
                                   "~/",
                                   "Text Files (*.txt);;All Files (*.*)");

  // If filename does not include .csp extension, add extension
  if (!(fileName.endsWith(".txt"))) {
    fileName.append(".txt");
  }

  // If filename exists, warn user
  if (QFileInfo(fileName).exists()
      && QMessageBox::warning(this,
                              "File Exists",
                              "Do you wish to overwrite " + fileName + "?",
                              QMessageBox::Yes | QMessageBox::No)
          != QMessageBox::Yes)
  {
    return;  // Return if answer is no
  }

  // Save protocol to file
  QFile file(fileName);  // Open file
  if (!file.open(QIODevice::WriteOnly))
  {  // Open file, return error if unable to do so
    QMessageBox::warning(
        this, "Error", "Unable to save file: Please check folder permissions.");
    return;
  }

  QTextStream stream(&file);
  size_t n = 0;
  switch (plot_mode) {
    case TIMESERIES:
      n = tscurve->dataSize();
      for (int i = 0; i < n; i++) {
        stream << time_buffer[i] << " " << channel_data[i] << "\n";
      }
      break;
    case SCATTER:
      break;
    case FFT:
      n = fftcurve->dataSize();
      for (int i = 0; i < n; i++) {
        stream << fft_output_x[i] << " " << fft_output_y[i] << "\n";
      }
      break;
    default:
      break;
  }
  file.close();
}

void analysis_module::Panel::changeChannel(QModelIndex /*id*/)
{
  getTrialData();
}

void analysis_module::Panel::resetAxes()
{
  emit setPlotRange(xmin, xmax, ymin, ymax);
}

void analysis_module::Panel::updateWindow(int index)
{
  window_shape = static_cast<window_t>(index);
  if (plot_mode == FFT) {
    getTrialData();
  }
}

void analysis_module::Panel::updateKalpha(double KalphaInput)
{
  Kalpha = KalphaInput;
  if (plot_mode == FFT) {
    getTrialData();
  }
}

void analysis_module::Panel::updateCalpha(double CalphaInput)
{
  Calpha = CalphaInput;
  if (plot_mode == FFT) {
    getTrialData();
  }
}

void analysis_module::Panel::makeWindow(int num_points)
{
  switch (window_shape) {
    case RECT:  // rectangular
      disc_window = new RectangularWindow(num_points);
      break;
    case TRI:  // triangular
      disc_window = new TriangularWindow(num_points, 1);
      break;
    case HAMM:  // Hamming
      disc_window = new HammingWindow(num_points);
      break;
    case HANN:  // Hann
      disc_window = new HannWindow(num_points, 1);
      break;
    case CHEBY:  // Dolph-Chebyshev
      disc_window = new DolphChebyWindow(num_points, Calpha);
      break;
    case KAISER:
      disc_window = new KaiserWindow(num_points, Kalpha);
      break;
  }  // end of switch on window_shape
}

void analysis_module::Panel::screenshot()
{
  QwtPlotRenderer renderer;
  renderer.exportTo(omniplot, "screenshot.pdf");
}

void analysis_module::Panel::toggleFWR(bool fwrStatus)
{
  fwrChecked = fwrStatus;
  if (plot_mode == FFT) {
    getTrialData();
  }
}

// TODO: may need to restore toggle functions to allow plots to be cleared when
// deselected
void analysis_module::Panel::changeDataFile()
{
  QFileDialog fileDialog(this);
  fileDialog.setFileMode(QFileDialog::AnyFile);
  fileDialog.setWindowTitle("Select Data File");

  QSettings userprefs;
  QSettings::setPath(QSettings::NativeFormat,
                     QSettings::SystemScope,
                     "/usr/local/share/rtxi/");
  fileDialog.setDirectory(
      userprefs.value("/dirs/data", getenv("HOME")).toString());

  QStringList filterList;
  filterList << "*.h5";
  fileDialog.setNameFilters(filterList);

  QStringList files;
  QString filename;
  if (fileDialog.exec() != 0) {
    closeFile();  // close previous file
    treeViewer->clear();
    files = fileDialog.selectedFiles();
    filename = files[0];
    fileNameEdit->setText(filename);
    openFile(filename);
  }
}

// TODO: populate HDF5, attribute, and parameter viewer contents
//        enable any scatter/FFT specific options
int analysis_module::Panel::openFile(QString& filename)
{
  herr_t status = -1;
  currentTrialFlag = 0;
  currentTrial = "";
  firstChannelSelected = 0;

  if (QFile::exists(filename)) {
    file_id =
        H5Fopen(filename.toLatin1().constData(), H5F_ACC_RDONLY, H5P_DEFAULT);
    // Iterate through file
    status = H5Ovisit(
        file_id, H5_INDEX_NAME, H5_ITER_NATIVE, op_func, this->treeViewer);
    if (status == 0) {
      plotControls->setEnabled(true);
      plotOptions->setEnabled(true);
      getTrialData();
      updatePlot();
    } else {
      closeFile();
    }
  }
  return status;
}

// TODO: erase HDF5, attribute, and parameter viewer contents
//        disable plot button and any scatter/FFT specific options
void analysis_module::Panel::closeFile()
{
  if (file_id != 0) {
    H5Fclose(file_id);
    file_id = 0;
  }
}

void analysis_module::Panel::getTrialData()
{
  QTreeWidgetItem* current_tree_item = this->treeViewer->currentItem();
  if (current_tree_item == nullptr) {
    return;
  }
  QVariant data_id_variant = current_tree_item->data(0, Qt::UserRole);
  if (!data_id_variant.isValid()) {
    return;
  }
  herr_t status = 0;
  hsize_t nrecords = 0;
  hid_t packettable_id = 0;
  double channelDataSum = 0;
  double channelDataMean = NAN;

  // Open packet table
  packettable_id =
      H5PTopen(file_id, data_id_variant.value<QString>().toLatin1());
  if (packettable_id < 0) {
    ERROR_MSG("analysis_module::Panel::getTrialData : H5PTopen error {}",
              packettable_id);
    return;
  }

  // Get packet count
  status = H5PTget_num_packets(packettable_id, &nrecords);
  if (status < 0) {
    ERROR_MSG(
        "analysis_module::Panel::getTrialData : H5PTget_num_packets error {}",
        status);
    return;
  }
  if (nrecords <= 1) {
    return;
  }

  // Initialize data buffer -- module will crash for large trials...
  data_buffer.resize(nrecords);
  channel_data.resize(nrecords);
  time_buffer.resize(nrecords);

  // Read data
  status = H5PTread_packets(
      packettable_id, 0, static_cast<int>(nrecords), data_buffer.data());
  if (status < 0) {
    printf("Throw error - H5PTread_packets error %d\n", status);
  }

  for (int i = 0; i < nrecords; i++) {
    channel_data[i] = data_buffer[i].value;
    channelDataSum = channelDataSum + channel_data[i];
  }
  // find DC offset -- for long signals, it's approx. the mean
  channelDataMean = channelDataSum / (static_cast<double>(nrecords));

  const auto start_time = static_cast<double>(data_buffer[0].time);
  // Build time vector
  for (int i = 1; i < nrecords; i++) {
    time_buffer[i] = static_cast<double>(data_buffer[i].time) - start_time;
  }

  // FFT plot
  // TODO: check if FFT plot is selected
  double windowedChannelVal = NAN;
  int fft_length = 2;
  while (fft_length < nrecords) {
    fft_length = fft_length * 2;
  }
  fft_input.resize(fft_length);
  fft_buffer.resize(fft_length);
  fft_output_y.resize(fft_length);
  fft_output_x.resize(fft_length);
  makeWindow(static_cast<int>(nrecords));
  for (int i = 0; i < fft_length; i++) {
    if (i < nrecords) {
      windowedChannelVal = disc_window->GetDataWinCoeff(i) * channel_data[i];
      fft_input[i] = complex(windowedChannelVal, 0.0);
    } else {
      fft_input[i] = 0;  // pad the rest of the array with 0's
    }
  }
  fft(fft_input.data(), fft_buffer.data(), fft_length);
  // We have to use a consistent period so we'll just use the average. This will give
  // really wrong answers if too many data points are skipped or lost
  const double avg_period = (time_buffer.back() - time_buffer.front())
      / static_cast<double>(time_buffer.size());
  for (size_t i = 0; i < fft_length; i++) {
    fft_output_y[i] = fabs(real(fft_buffer[i]));
    fft_output_x[i] = static_cast<double>(i) / ((avg_period/1e9) * fft_length);
  }

  // Full wave rectification
  if (fwrChecked) {
    // subtract DC offset from each value and take absolute value
    for (int i = 0; i < static_cast<int>(nrecords); i++) {
      channel_data[i] = std::abs(channel_data[i] - channelDataMean);
    }
  }

  // Plot
  tscurve->setRawSamples(
      time_buffer.data(), channel_data.data(), static_cast<int>(nrecords));
  fftcurve->setRawSamples(fft_output_x.data(), fft_output_y.data(), fft_length);
  updatePlot();

  // Close identifiers
  H5PTclose(packettable_id);
}

// Temporary function for validating data access
void analysis_module::Panel::updatePlot()
{
  tscurve->detach();
  fftcurve->detach();

  switch (plot_mode) {
    case TIMESERIES:
      tscurve->attach(omniplot);
      break;
    case SCATTER:
      break;
    case FFT:
      fftcurve->attach(omniplot);
      break;
    default:
      break;
  }

  omniplot->setAxisAutoScale(BasicPlot::yLeft, true);
  omniplot->setAxisAutoScale(BasicPlot::xBottom, true);
  omniplot->replot();
  xmin = omniplot->axisScaleDiv(QwtPlot::xBottom).lowerBound();
  xmax = omniplot->axisScaleDiv(QwtPlot::xBottom).upperBound();
  ymin = omniplot->axisScaleDiv(QwtPlot::yLeft).lowerBound();
  ymax = omniplot->axisScaleDiv(QwtPlot::yLeft).upperBound();
  emit setPlotRange(xmin, xmax, ymin, ymax);
}

void analysis_module::Panel::updatePlotMode(int mode)
{
  plot_mode = static_cast<plot_t>(mode);
  updatePlot();
}

// Temporary function for validating data access
void analysis_module::Panel::dump_vals(double* data, hsize_t* /*ndims*/)
{
  // Only printing first value out or else the printf will block
  for (size_t i = 0; i < 10; i++) {
    printf("value is %f\n", data[i]);
  }
}

herr_t analysis_module::op_func(hid_t loc_id,
                                const char* name,
                                const H5O_info_t* info,
                                void* operator_data)
{
  auto* tree = reinterpret_cast<QTreeWidget*>(operator_data);
  QTreeWidgetItem* tree_item = nullptr;
  QTreeWidgetItem* parent_item = nullptr;
  QList<QTreeWidgetItem*> matching_groups;
  const QString qName = QString(name);
  QStringList split_path = qName.split("/", Qt::SkipEmptyParts);
  switch (info->type) {
    case H5O_TYPE_GROUP:
      tree_item = new QTreeWidgetItem;
      tree_item->setText(0, split_path.back());
      if (split_path.size() == 1) {
        tree->addTopLevelItem(tree_item);
      } else {
        split_path.pop_back();
        matching_groups = tree->findItems(split_path.back(), Qt::MatchExactly);
        matching_groups[0]->addChild(tree_item);
      }
      break;
    case H5O_TYPE_DATASET:
      tree_item = new QTreeWidgetItem;
      tree_item->setText(0, split_path.back());
      tree_item->setData(0, Qt::UserRole, QVariant::fromValue<QString>(qName));
      split_path.pop_back();
      matching_groups = tree->findItems(split_path.back(),
                                        Qt::MatchExactly | Qt::MatchRecursive);
      matching_groups.back()->addChild(tree_item);
      break;
    default:
      break;
  }
  return 0;
}

///////// DO NOT MODIFY BELOW //////////
// The exception is if your plugin is not going to need real-time functionality.
// For this case just replace the craeteRTXIComponent return type to nullptr.
// RTXI will automatically handle that case and won't attach a component to the
// real time thread for your plugin.

std::unique_ptr<Widgets::Plugin> createRTXIPlugin(Event::Manager* ev_manager)
{
  return std::make_unique<analysis_module::Plugin>(ev_manager);
}

Widgets::Panel* createRTXIPanel(QMainWindow* main_window,
                                Event::Manager* ev_manager)
{
  return new analysis_module::Panel(main_window, ev_manager);
}

std::unique_ptr<Widgets::Component> createRTXIComponent(
    Widgets::Plugin* host_plugin)
{
  return nullptr;
}

Widgets::FactoryMethods fact;

extern "C"
{
Widgets::FactoryMethods* getFactories()
{
  fact.createPanel = &createRTXIPanel;
  fact.createComponent = &createRTXIComponent;
  fact.createPlugin = &createRTXIPlugin;
  return &fact;
}
};

//////////// END //////////////////////
