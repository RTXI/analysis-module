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
  period_buffer.clear();
  fft_output_y.clear();
  fft_output_x.clear();
  fft_input.clear();
  fft_buffer.clear();
}

void analysis_module::Panel::customizeGUI()
{
  auto* customlayout = new QVBoxLayout(this);

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
  customlayout->addLayout(plotColumnLayout);
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
    status = H5Ovisit(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, op_func, nullptr);
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

// TODO: think through error cases here (e.g. when one of the top-level groups
// are selected, etc.)
void analysis_module::Panel::getTrialData()
{
  // TODO: check that current item is a dataset (and not a group), only
  // open/plot if a dataset is selected (maybe display an warning otherwise?)
  herr_t status = 0;
  hsize_t nrecords = 0;
  hsize_t ntrials = 0;
  hsize_t nchannels = 0;
  hid_t packettable_id = 0;
  hid_t trial_id = 0;
  hid_t period_id = 0;
  double channelDataSum = 0;
  double channelDataMean = NAN;

  // Check if selected channel is actually a channel (should have 0 children in
  // the treeViewer)
  if (treeViewer->currentItem()->childCount() != 0
      || treeViewer->currentItem()->text(0) == "Asynchronous Data"
      || treeViewer->currentItem()->text(0) == "Synchronous Data")
  {
    printf("getTrialData error: selected channel is not a dataset\n");
    return;
  }

  // Reset data buffers if something is already plotted

  // Get elements from GUI
  QString selectedTrial = treeViewer->currentItem()->text(0);
  QString channelNum = selectedTrial.split(" ").at(0);
  int channelNumInt = channelNum.toInt();
  QString trialToRead = treeViewer->currentItem()->parent()->parent()->text(0)
      + "/" + treeViewer->currentItem()->parent()->text(0);
  QString channelToRead = treeViewer->currentItem()->parent()->parent()->text(0)
      + "/" + treeViewer->currentItem()->parent()->text(0) + "/Channel Data";
  QString periodToRead =
      treeViewer->currentItem()->parent()->parent()->text(0) + "/Period (ns)";
  QString trialLengthToRead =
      treeViewer->currentItem()->parent()->parent()->text(0)
      + "/Trial Length (ns)";

  // Open packet table
  packettable_id = H5PTopen(file_id, channelToRead.toLatin1().constData());
  if (packettable_id < 0) {
    printf("Throw error - H5PTopen error %ld\n", packettable_id);
  }

  // Get packet count
  status = H5PTget_num_packets(packettable_id, &nrecords);
  if (status < 0) {
    printf("Throw error - H5PTget_num_packets error %d\n", status);
  }

  // Get number of trials
  status = H5Gget_num_objs(file_id, &ntrials);
  if (status < 0) {
    printf("Throw error - H5Gget_num_objs %d\n", status);
  }

  // Get identifier for trial and group
  printf("%s\n", trialToRead.toStdString().c_str());
  trial_id = H5Gopen1(file_id, trialToRead.toLatin1().constData());
  if (trial_id < 0) {
    printf("Throw error - H5Gopen1 %ld\n", trial_id);
  }

  // Get number of channels from trial
  // Returned number of objects includes "Channel Data" struct, so we subtract 1
  status = H5Gget_num_objs(trial_id, &nchannels);
  if (status < 0) {
    printf("Throw error - H5Gget_num_objs %d\n", status);
  }
  nchannels--;

  // Initialize data buffer -- module will crash for large trials...
  data_buffer =
      static_cast<double*>(malloc(sizeof(double) * static_cast<int>(nrecords)
                                  * static_cast<int>(nchannels)));
  channel_data =
      static_cast<double*>(malloc(sizeof(double) * static_cast<int>(nrecords)));
  time_buffer =
      static_cast<double*>(malloc(sizeof(double) * static_cast<int>(nrecords)));
  period_buffer = static_cast<double*>(malloc(sizeof(double)));

  // Read data
  status = H5PTread_packets(
      packettable_id, 0, static_cast<int>(nrecords), data_buffer);
  if (status < 0) {
    printf("Throw error - H5PTread_packets error %d\n", status);
  }

  // Build channel_data array -- pretty messy but seems unavoidable with the
  // current Data Recorder structure. To eliminate a loop later on, find
  // channelDataSum here if full wave rectification option is checked
  if (fwrChecked) {
    int j = 0;
    for (int i = channelNumInt - 1;
         i < static_cast<int>(nrecords) * static_cast<int>(nchannels);
         i = i + static_cast<int>(nchannels))
    {
      channel_data[j] = data_buffer[i];
      channelDataSum = channelDataSum + channel_data[j];
      j++;
    }
    // find DC offset -- for long signals, it's approx. the mean
    channelDataMean = channelDataSum / ((int)nrecords);
  } else {
    int j = 0;
    for (int i = channelNumInt - 1; i < (int)nrecords * (int)nchannels;
         i = i + nchannels)
    {
      channel_data[j] = data_buffer[i];
      j++;
    }
  }

  // Build time vector
  period_id =
      H5Dopen2(file_id, periodToRead.toLatin1().constData(), H5P_DEFAULT);
  if (period_id < 0) {
    printf("Throw error - H5Dopen2 error %ld\n", period_id);
  }
  status = H5Dread(period_id,
                   H5T_NATIVE_DOUBLE,
                   H5S_ALL,
                   H5S_ALL,
                   H5P_DEFAULT,
                   period_buffer);
  if (status < 0) {
    printf("Throw error - H5Dread error %d\n", status);
  }
  time_buffer[0] = 0;
  for (int i = 1; i < static_cast<int>(nrecords); i++) {
    time_buffer[i] = time_buffer[i - 1] + (*period_buffer / 1e9);
  }

  // FFT plot
  // TODO: check if FFT plot is selected
  double windowedChannelVal = NAN;
  int fft_length = 2;
  while (fft_length < (int)nrecords) {
    fft_length = fft_length * 2;
  }
  fft_input = static_cast<complex*>(malloc(sizeof(complex) * fft_length));
  fft_buffer = static_cast<complex*>(malloc(sizeof(complex) * fft_length));
  fft_output_y = static_cast<double*>(malloc(sizeof(double) * fft_length));
  fft_output_x = static_cast<double*>(malloc(sizeof(double) * fft_length));
  makeWindow(static_cast<int>(nrecords));
  for (int i = 0; i < fft_length; i++) {
    if (i < static_cast<int>(nrecords)) {
      windowedChannelVal = disc_window->GetDataWinCoeff(i) * channel_data[i];
      fft_input[i] = complex(windowedChannelVal, 0.0);
    } else {
      fft_input[i] = 0;  // pad the rest of the array with 0's
    }
  }
  fft(fft_input, fft_buffer, fft_length);
  for (int i = 0; i < fft_length; i++) {
    fft_output_y[i] = fabs(real(fft_buffer[i]));
    fft_output_x[i] = i / ((*period_buffer / 1e9) * fft_length);
  }

  // Full wave rectification
  if (fwrChecked) {
    // subtract DC offset from each value and take absolute value
    for (int i = 0; i < static_cast<int>(nrecords); i++) {
      channel_data[i] = std::abs(channel_data[i] - channelDataMean);
    }
  }

  // Plot
  tscurve->setRawSamples(time_buffer, channel_data, static_cast<int>(nrecords));
  fftcurve->setRawSamples(fft_output_x, fft_output_y, fft_length);
  updatePlot();

  // Close identifiers
  H5PTclose(packettable_id);
  H5Gclose(trial_id);
  H5Dclose(period_id);
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

// TODO: clean up treeViewer -- no need to list full path for each parent/child
herr_t analysis_module::op_func(hid_t loc_id,
                                const char* name,
                                const H5O_info_t* info,
                                void* operator_data)
{
  QString qName = QString(name);
  if (name[0] == '.') {
  } else {
    switch (info->type) {
      case H5O_TYPE_GROUP:
        // printf ("%s  (Group)\n", name);
        if (!qName.startsWith(currentTrial) || currentTrial == "") {
          currentTrial = qName;
          treeParent = new QTreeWidgetItem;
          treeParent->setText(0, qName);
          treeViewer->addTopLevelItem(treeParent);
        } else if (qName.startsWith(currentTrial)
                   && qName.endsWith("Asynchronous Data"))
        {
          treeChild1 = new QTreeWidgetItem;
          currentGroup = "Asynchronous Data";
          treeChild1->setText(0, "Asynchronous Data");
          treeParent->addChild(treeChild1);
        } else if (qName.startsWith(currentTrial)
                   && qName.endsWith("Synchronous Data"))
        {
          treeChild1 = new QTreeWidgetItem;
          currentGroup = "Synchronous Data";
          treeChild1->setText(0, "Synchronous Data");
          treeParent->addChild(treeChild1);
        }
        break;
      case H5O_TYPE_DATASET:
        // printf ("%s  (Dataset)\n", name);
        if (qName.startsWith(currentTrial + "/Asynchronous Data")
            && !qName.endsWith("Channel Data"))
        {
          treeChild2 = new QTreeWidgetItem;
          treeChild2->setText(0,
                              qName.right(qName.length()
                                          - (currentTrial.length()
                                             + currentGroup.length() + 2)));
          treeChild1->addChild(treeChild2);
          treeChild2->setToolTip(0, qName);
          if (firstChannelSelected == 0) {
            firstChannelSelected = 1;
            treeViewer->setCurrentItem(treeChild2);
          }
        } else if (qName.startsWith(currentTrial + "/Synchronous Data")
                   && !qName.endsWith("Channel Data"))
        {
          treeChild2 = new QTreeWidgetItem;
          treeChild2->setText(0,
                              qName.right(qName.length()
                                          - (currentTrial.length()
                                             + currentGroup.length() + 2)));
          treeChild1->addChild(treeChild2);
          treeChild2->setToolTip(0, qName);
          if (firstChannelSelected == 0) {
            firstChannelSelected = 1;
            treeViewer->setCurrentItem(treeChild2);
          }
        }
        // else statement -- add all remaining dataset contents to parameters
        break;
      case H5O_TYPE_NAMED_DATATYPE:
      default:
        break;
    }
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
