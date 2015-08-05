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

/*
 * HDF viewer and analysis tools
 */

#include <QtGui>
#include <QtGlobal>
#include <algorithm>

#include <main_window.h>
#include <qwt_plot_renderer.h>
#include "analysis-tools.h"

#include <time.h>
#include <gsl/gsl_math.h>

#include <iostream>

extern "C" Plugin::Object *createRTXIPlugin(void) {
	return new AnalysisTools();
}

// Globals
QTreeWidget *treeViewer;
QTreeWidgetItem *treeParent;
QTreeWidgetItem *treeChild1;
QTreeWidgetItem *treeChild2;
QString currentTrial, currentGroup;
int currentTrialFlag;
int firstChannelSelected;
double *data_buffer, *channel_data, *time_buffer, *period_buffer, *fft_output_y, *fft_output_x;
complex *fft_input, *fft_buffer;

static DefaultGUIModel::variable_t vars[] = {
	{ "Input", "Input", DefaultGUIModel::INPUT, },
	{ "Output", "Output", DefaultGUIModel::OUTPUT, }, 
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

AnalysisTools::AnalysisTools(void) :  DefaultGUIModel("Analysis Tools", ::vars, ::num_vars) {
	setWhatsThis(
			"<p><b>Analysis Tools</b></p><p>Analysis tools</p>"); // TO-DO: add detail here
	initParameters();
	DefaultGUIModel::createGUI(vars, num_vars); // this is required to create the GUI
	customizeGUI();
	update(INIT);
	refresh(); // this is required to update the GUI with parameter and state values
	QTimer::singleShot(0, this, SLOT(resizeMe()));
}

AnalysisTools::~AnalysisTools(void)
{
	if(file_id) {
		H5Fclose(file_id);
		file_id = NULL;
	}
	freePlotBuffers();
}

void AnalysisTools::execute(void) {
	return;
}

void AnalysisTools::update(DefaultGUIModel::update_flags_t flag) {
	switch (flag) {
		case INIT:
			break;
		case MODIFY:
			break;
		case UNPAUSE:
			break;
		case PAUSE:
			break;
		case PERIOD:
			break;
		default:
			break;
	}
};

void AnalysisTools::initParameters() {
	fwrChecked = 0;
	file_id = NULL;
	data_buffer = NULL;
	channel_data = NULL;
	time_buffer = NULL;
	period_buffer = NULL;
	dataset_id = NULL;
	
	window_shape = RECT;
	Kalpha = 1.5;
	Calpha = 70;
}

void AnalysisTools::customizeGUI(void) {
	QGridLayout *customlayout = DefaultGUIModel::getLayout(); 

	// Screenshot buttons
	QGroupBox *plotBox = new QGroupBox("Save Screenshot");
	QHBoxLayout *plotBoxLayout = new QHBoxLayout;
	plotBox->setLayout(plotBoxLayout);
	QPushButton *saveTSPlotButton = new QPushButton("Time Series Plot");
	saveTSPlotButton->setToolTip("Save screenshot of the time series plot");
	plotBoxLayout->addWidget(saveTSPlotButton);
	QPushButton *saveScatterPlotButton = new QPushButton("Scatter Plot");
	saveScatterPlotButton->setToolTip("Save screenshot of the scatter plot");
	plotBoxLayout->addWidget(saveScatterPlotButton);
	QPushButton *saveFFTPlotButton = new QPushButton("FFT Plot");
	saveFFTPlotButton->setToolTip("Save screenshot of the FFT plot");
	plotBoxLayout->addWidget(saveFFTPlotButton);
	customlayout->addWidget(plotBox, 0, 2, 1, 4);
	plotBox->setMinimumSize(910, 0);

	// Initialize plots
	QGroupBox *tsplotBox = new QGroupBox("Time Series Plot");
	QHBoxLayout *tsplotBoxLayout = new QHBoxLayout;
	tsplotBox->setLayout(tsplotBoxLayout);
	tsplot = new BasicPlot(this);
	tsplot->setFixedSize(450, 270);
	tsplotBoxLayout->addWidget(tsplot);
	customlayout->addWidget(tsplotBox, 1, 2, 3, 2);
	tscurve = new QwtPlotCurve;
	tscurve->attach(tsplot);

	QGroupBox *scatterplotBox = new QGroupBox("Scatter Plot");
	QHBoxLayout *scatterplotBoxLayout = new QHBoxLayout;
	scatterplotBox->setLayout(scatterplotBoxLayout);
	splot = new ScatterPlot(this);
	splot->setFixedSize(450, 270);
	scatterplotBoxLayout->addWidget(splot);
	customlayout->addWidget(scatterplotBox, 1, 4, 3, 2);

	QGroupBox *fftplotBox = new QGroupBox("FFT Plot");
	QHBoxLayout *fftplotBoxLayout = new QHBoxLayout;
	fftplotBox->setLayout(fftplotBoxLayout);
	fftplot = new BasicPlot(this);
	fftplot->setFixedSize(675, 270);
	fftplotBoxLayout->addWidget(fftplot);
	customlayout->addWidget(fftplotBox, 4, 3, 3, 3);
	fftcurve = new QwtPlotCurve;
	fftcurve->attach(fftplot);

	// Connect screenshot buttons to functions
	QObject::connect(saveTSPlotButton, SIGNAL(clicked()), this, SLOT(screenshotTS()));
	QObject::connect(saveScatterPlotButton, SIGNAL(clicked()), this, SLOT(screenshotScatter()));
	QObject::connect(saveFFTPlotButton, SIGNAL(clicked()), this, SLOT(screenshotFFT()));

	// Plot options
	QGroupBox *plotOptionsBox = new QGroupBox(tr("Plot Options"));
	QGridLayout *plotOptionsBoxLayout = new QGridLayout;
	plotOptionsBox->setLayout(plotOptionsBoxLayout);
	QVBoxLayout *plotOptionsVerticalLayout = new QVBoxLayout;
	QButtonGroup *plotOptionsButtons = new QButtonGroup;
	plotOptionsButtons->setExclusive(false);
	QCheckBox *FWRCheckBox = new QCheckBox("TS: Full wave rectify");
	plotOptionsVerticalLayout->addWidget(FWRCheckBox);
	plotOptionsButtons->addButton(FWRCheckBox);
	FWRCheckBox->setChecked(false);
	QObject::connect(FWRCheckBox, SIGNAL(toggled(bool)),this,SLOT(toggleFWR(bool)));
	FWRCheckBox->setToolTip("Enable full wave rectification of time series plot");
	plotOptionsBoxLayout->addLayout(plotOptionsVerticalLayout, 0, 0);
	QLabel *windowLabel = new QLabel("FFT window shape:");
	windowShape = new QComboBox;
	windowShape->insertItem(1, "Rectangular");
	windowShape->insertItem(2, "Triangular (Bartlett)");
	windowShape->insertItem(3, "Hamming");
	windowShape->insertItem(4, "Hann");
	windowShape->insertItem(5, "Chebyshev");
	windowShape->insertItem(6, "Kaiser");
	QObject::connect(windowShape, SIGNAL(activated(int)), this, SLOT(updateWindow(int)));
	windowShape->setToolTip("Choose a window to apply for the FFT plot. For no window, choose Rectangular.");
	plotOptionsBoxLayout->addWidget(windowLabel, 1, 0);
	plotOptionsBoxLayout->addWidget(windowShape, 1, 1);
	QLabel *kalphaLabel = new QLabel("Kaiser Alpha");
	plotOptionsBoxLayout->addWidget(kalphaLabel, 2, 0);
	QDoubleSpinBox *kalphaEdit = new QDoubleSpinBox(plotOptionsBox);
	kalphaEdit->setValue(Kalpha);
	QObject::connect(kalphaEdit, SIGNAL(valueChanged(double)), this, SLOT(updateKalpha(double)));
	kalphaEdit->setToolTip("Attenuation parameter for Kaiser window");
	plotOptionsBoxLayout->addWidget(kalphaEdit, 2, 1);
	QLabel *calphaLabel = new QLabel("Chebyshev (dB)");
	plotOptionsBoxLayout->addWidget(calphaLabel, 3, 0);
	QDoubleSpinBox *calphaEdit = new QDoubleSpinBox(plotOptionsBox);
	calphaEdit->setValue(Calpha);
	QObject::connect(calphaEdit, SIGNAL(valueChanged(double)), this, SLOT(updateCalpha(double)));
	calphaEdit->setToolTip("Attenuation parameter for Chebyshev window");
	plotOptionsBoxLayout->addWidget(calphaEdit, 3, 1);
	customlayout->addWidget(plotOptionsBox, 2, 0, 2, 1);

	// Global plot options
	QGroupBox *optionBox = new QGroupBox;
	QGridLayout *optionBoxLayout = new QGridLayout;
	optionBox->setLayout(optionBoxLayout);
	QHBoxLayout *plotSelectionLayout = new QHBoxLayout;
	QButtonGroup *optionButtons = new QButtonGroup;
	optionButtons->setExclusive(false);
	QCheckBox *plotTSCheckBox = new QCheckBox("Time Series");
	plotSelectionLayout->addWidget(plotTSCheckBox);
	optionButtons->addButton(plotTSCheckBox);
	plotTSCheckBox->setChecked(true);
	QCheckBox *plotScatterCheckBox = new QCheckBox("Scatter");
	plotSelectionLayout->addWidget(plotScatterCheckBox);
	optionButtons->addButton(plotScatterCheckBox);
	plotScatterCheckBox->setChecked(true);
	QCheckBox *plotFFTCheckBox = new QCheckBox("FFT");
	plotSelectionLayout->addWidget(plotFFTCheckBox);
	optionButtons->addButton(plotFFTCheckBox);
	plotFFTCheckBox->setChecked(true);
	// TO-DO: add toggles for all related plot options
	QObject::connect(plotTSCheckBox,SIGNAL(toggled(bool)),tsplot,SLOT(setEnabled(bool)));
	QObject::connect(plotTSCheckBox,SIGNAL(toggled(bool)),saveTSPlotButton,SLOT(setEnabled(bool)));
	QObject::connect(plotTSCheckBox,SIGNAL(toggled(bool)),FWRCheckBox,SLOT(setEnabled(bool)));
	//QObject::connect(plotTSCheckBox,SIGNAL(toggled(bool)),this,SLOT(toggleTSplot(bool)));
	QObject::connect(plotScatterCheckBox,SIGNAL(toggled(bool)),splot,SLOT(setEnabled(bool)));
	QObject::connect(plotScatterCheckBox,SIGNAL(toggled(bool)),saveScatterPlotButton,SLOT(setEnabled(bool)));
	//QObject::connect(plotScatterCheckBox,SIGNAL(toggled(bool)),this,SLOT(toggleScatterplot(bool)));
	QObject::connect(plotFFTCheckBox,SIGNAL(toggled(bool)),fftplot,SLOT(setEnabled(bool)));
	QObject::connect(plotFFTCheckBox,SIGNAL(toggled(bool)),saveFFTPlotButton,SLOT(setEnabled(bool)));
	//QObject::connect(plotFFTCheckBox,SIGNAL(toggled(bool)),this,SLOT(toggleFFTplot(bool)));
	plotTSCheckBox->setToolTip("Enable time series plot");
	plotScatterCheckBox->setToolTip("Enable scatter plot");
	plotFFTCheckBox->setToolTip("Enable FFT plot");
	optionBoxLayout->addLayout(plotSelectionLayout, 0, 0);

	QVBoxLayout *plotButtonLayout = new QVBoxLayout;
	plotButton = new QPushButton("Plot");
	plotButtonLayout->addWidget(plotButton);
	plotButton->setEnabled(false);
	QObject::connect(plotButton,SIGNAL(released(void)), this, SLOT(plotTrial(void)));
	plotButton->setToolTip("Plot data for selected trial and channel");
	optionBoxLayout->addLayout(plotButtonLayout, 1, 0);
	customlayout->addWidget(optionBox, 1, 0, 1, 1);

	// File control
	QGroupBox *fileBox = new QGroupBox(tr("File Control"));
	QHBoxLayout *fileLayout = new QHBoxLayout;
	fileBox->setLayout(fileLayout);
	fileLayout->addWidget(new QLabel(tr("File Name")));
	fileNameEdit = new QLineEdit;
	fileNameEdit->setReadOnly(true);
	fileLayout->addWidget(fileNameEdit);
	QPushButton *fileChangeButton = new QPushButton("Choose File");
	fileLayout->addWidget(fileChangeButton);
	QObject::connect(fileChangeButton,SIGNAL(released(void)),this,SLOT(changeDataFile(void)));
	customlayout->addWidget(fileBox, 0, 0, 1, 1);

	// HDF5 viewer
	treeViewer = new QTreeWidget;
	treeViewer->setHeaderLabels(QStringList("HDF5 Viewer"));
	customlayout->addWidget(treeViewer, 4, 0, 3, 1);

	// Parameters
	QGroupBox *paramsBox = new QGroupBox(tr("Parameters"));
	QHBoxLayout *paramsLayout = new QHBoxLayout;
	paramsBox->setLayout(paramsLayout);
	parameterView = new QTextEdit;
	parameterView->setFixedWidth(220);
	paramsLayout->addWidget(parameterView);
	paramsBox->setMaximumSize(250, 500);
	customlayout->addWidget(paramsBox, 4, 2, 3, 1);

	// Standard module buttons
	DefaultGUIModel::pauseButton->setEnabled(false);
	DefaultGUIModel::modifyButton->setEnabled(false);
	DefaultGUIModel::unloadButton->setToolTip("Close module");

	setLayout(customlayout);
}

void AnalysisTools::updateWindow(int index) {
	if (index == 0) {
		window_shape = RECT;
	} else if (index == 1) {
		window_shape = TRI;
	} else if (index == 2) {
		window_shape = HAMM;
	} else if (index == 3) {
		window_shape = HANN;
	} else if (index == 4) {
		window_shape = CHEBY;
	} else if (index == 5) {
		window_shape = KAISER;
	}
}

void AnalysisTools::updateKalpha(double KalphaInput) {
	Kalpha = KalphaInput;
}

void AnalysisTools::updateCalpha(double CalphaInput) {
	Calpha = CalphaInput;
}

void AnalysisTools::makeWindow(int num_points) {
	switch (window_shape) {
		case RECT: // rectangular
			disc_window = new RectangularWindow(num_points);
			break;
		case TRI: // triangular
			disc_window = new TriangularWindow(num_points, 1);
			break;
		case HAMM: // Hamming
			disc_window = new HammingWindow(num_points);
			break;
		case HANN: // Hann
			disc_window = new HannWindow(num_points, 1);
			break;
		case CHEBY: // Dolph-Chebyshev
			disc_window = new DolphChebyWindow(num_points, Calpha);
			break;
		case KAISER:
			disc_window = new KaiserWindow(num_points, Kalpha);
			break;
	} // end of switch on window_shape
}

void AnalysisTools::screenshotTS() {
	QwtPlotRenderer renderer;
	renderer.exportTo(tsplot,"screenshot.pdf");
}

void AnalysisTools::screenshotScatter() {
	QwtPlotRenderer renderer;
	renderer.exportTo(splot,"screenshot.pdf");
}

void AnalysisTools::screenshotFFT() {
	QwtPlotRenderer renderer;
	renderer.exportTo(fftplot,"screenshot.pdf");
}

void AnalysisTools::clearData() {
}

void AnalysisTools::toggleFWR(bool fwrStatus) {
	fwrChecked = fwrStatus;
}

// TO-DO: may need to restore toggle functions to allow plots to be cleared when deselected

void AnalysisTools::changeDataFile(void) {
	QFileDialog fileDialog(this);
	fileDialog.setFileMode(QFileDialog::AnyFile);
	fileDialog.setWindowTitle("Select Data File");

	QSettings userprefs;
	userprefs.setPath(QSettings::NativeFormat, QSettings::SystemScope, "/usr/local/share/rtxi/");
	fileDialog.setDirectory(userprefs.value("/dirs/data", getenv("HOME")).toString());

	QStringList filterList;
	filterList << "*.h5";
	fileDialog.setNameFilters(filterList);

	QStringList files;
	QString filename;
	if(fileDialog.exec()) {
		closeFile(); // close previous file
		treeViewer->clear();
		files = fileDialog.selectedFiles();
		filename = files[0];
		fileNameEdit->setText(filename);
		openFile(filename);
	}
}

// TO-DO: populate HDF5, attribute, and parameter viewer contents
//        enable any scatter/FFT specific options
int AnalysisTools::openFile(QString &filename) {
	herr_t status = -1;
	currentTrialFlag = 0;
	currentTrial = "";
	firstChannelSelected = 0;

	if (QFile::exists(filename)) {
		file_id = H5Fopen(filename.toLatin1().constData(), H5F_ACC_RDONLY, H5P_DEFAULT);
		// Iterate through file
		status = H5Ovisit(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, op_func, NULL);
		if (!status)
			plotButton->setEnabled(true);
		else
			closeFile();
	}
	return status;
}

// TO-DO: erase HDF5, attribute, and parameter viewer contents
//        disable plot button and any scatter/FFT specific options
void AnalysisTools::closeFile()
{
	if(file_id) {
		H5Fclose(file_id);
		file_id = NULL;
	}
	freePlotBuffers();
}

// TO-DO: clean up treeViewer -- no need to list full path for each parent/child
herr_t op_func(hid_t loc_id, const char *name, const H5O_info_t *info, void *operator_data) {
	QString qName = QString(name);
	//printf ("/");
	if (name[0] == '.') {
		//printf ("  (Group)\n");
	} else {
		switch (info->type) {
			case H5O_TYPE_GROUP:
				//printf ("%s  (Group)\n", name);
				if (!currentTrialFlag) {
					currentTrialFlag = 1;
					currentTrial = qName;
					treeParent = new QTreeWidgetItem;
					treeParent->setText(0, qName);
					treeViewer->addTopLevelItem(treeParent);
				} else if (qName.startsWith(currentTrial) && qName.endsWith("Asynchronous Data")) {
					treeChild1 = new QTreeWidgetItem;
					currentGroup = "Asynchronous Data";
					treeChild1->setText(0, "Asynchronous Data");
					treeParent->addChild(treeChild1);
				} else if (qName.startsWith(currentTrial) && qName.endsWith("Synchronous Data")) {
					currentTrialFlag = 0; // reset flag to start parsing next trial
					treeChild1 = new QTreeWidgetItem;
					currentGroup = "Synchronous Data";
					treeChild1->setText(0, "Synchronous Data");
					treeParent->addChild(treeChild1);
				}
				break;
			case H5O_TYPE_DATASET:
				//printf ("%s  (Dataset)\n", name);
				if (qName.startsWith(currentTrial + "/Asynchronous Data") && !qName.endsWith("Channel Data")) {
					treeChild2 = new QTreeWidgetItem;
					treeChild2->setText(0, qName.right(qName.length() - (currentTrial.length()+currentGroup.length()+2)));
					treeChild1->addChild(treeChild2);
					treeChild2->setToolTip(0, qName);
					if (!firstChannelSelected) {
						firstChannelSelected = 1;
						treeViewer->setCurrentItem(treeChild2);
					}
				} else if (qName.startsWith(currentTrial + "/Synchronous Data") && !qName.endsWith("Channel Data")) {
					treeChild2 = new QTreeWidgetItem;
					treeChild2->setText(0, qName.right(qName.length() - (currentTrial.length()+currentGroup.length()+2)));
					treeChild1->addChild(treeChild2);
					treeChild2->setToolTip(0, qName);
					if (!firstChannelSelected) {
						firstChannelSelected = 1;
						treeViewer->setCurrentItem(treeChild2);
					}
				}
				// else statement -- add all remaining dataset contents to parameters
				break;
			case H5O_TYPE_NAMED_DATATYPE:
				//printf ("%s  (Datatype)\n", name);
				break;
			default:
				//printf ("%s  (Unknown)\n", name);
				break;
		}
	}

	return 0;
}

// TO-DO: think through error cases here (e.g. when one of the top-level groups are selected, etc.)
void AnalysisTools::plotTrial() {
	// TO-DO: check that current item is a dataset (and not a group), only open/plot if a dataset is selected (maybe display an warning otherwise?)
	herr_t status;
	hsize_t nrecords, ntrials, nchannels;
	hid_t packettable_id, trial_id, period_id;
	double channelDataSum = 0;
	double channelDataMean;
	
	// Check if selected channel is actually a channel (should have 0 children in the treeViewer)
	if (treeViewer->currentItem()->childCount() != 0) {
		printf("plotTrial error: selected channel is not a dataset\n");
		return;
	}
	
	// Reset data buffers if something is already plotted
	freePlotBuffers();

	// Get elements from GUI
	QString selectedTrial = treeViewer->currentItem()->text(0);
	QString channelNum = selectedTrial.at(0);
	int channelNumInt = channelNum.toInt();
	QString trialToRead = treeViewer->currentItem()->parent()->parent()->text(0) + "/" +
			treeViewer->currentItem()->parent()->text(0);
	QString channelToRead = treeViewer->currentItem()->parent()->parent()->text(0) + "/" +
			treeViewer->currentItem()->parent()->text(0) + "/Channel Data";
	QString periodToRead = treeViewer->currentItem()->parent()->parent()->text(0) + "/Period (ns)";
	QString trialLengthToRead = treeViewer->currentItem()->parent()->parent()->text(0) + "/Trial Length (ns)";

	// Open packet table
	packettable_id = H5PTopen(file_id, channelToRead.toLatin1().constData());
	if(packettable_id < 0)
		printf("Throw error - H5PTopen error %d\n", packettable_id);

	// Get packet count
	status = H5PTget_num_packets(packettable_id, &nrecords);
	if(status < 0)
		printf("Throw error - H5PTget_num_packets error %d\n", status);

	// Get number of trials
	status = H5Gget_num_objs(file_id, &ntrials);
	if(status < 0)
		printf("Throw error - H5Gget_num_objs %d\n", status);

	// Get identifier for trial and group
	printf("%s\n", trialToRead.toStdString().c_str());
	trial_id = H5Gopen1(file_id, trialToRead.toLatin1().constData());
	if(trial_id < 0)
		printf("Throw error - H5Gopen1 %d\n", trial_id);

	// Get number of channels from trial
	// Returned number of objects includes "Channel Data" struct, so we subtract 1
	status = H5Gget_num_objs(trial_id, &nchannels);
	if(status < 0)
		printf("Throw error - H5Gget_num_objs %d\n", status);
	nchannels--;

	// Initialize data buffer -- module will crash for large trials...
	data_buffer = (double *)malloc(sizeof(double)*(int)(nrecords)*(int)nchannels);
	channel_data = (double *)malloc(sizeof(double)*(int)(nrecords));
	time_buffer = (double *)malloc(sizeof(double)*(int)(nrecords));
	period_buffer = (double *)malloc(sizeof(double));

	// Read data
	status = H5PTread_packets(packettable_id, 0, (int)nrecords, data_buffer);
	if(status < 0)
		printf("Throw error - H5PTread_packets error %d\n", status);

	// Build channel_data array -- pretty messy but seems unavoidable with the current Data Recorder structure
	// To eliminate a loop later on, find channelDataSum here if full wave rectification option is checked
	if (fwrChecked) {
		int j = 0;
		for (int i = channelNumInt-1; i < (int)nrecords*(int)nchannels; i = i + (int)nchannels) {
			channel_data[j] = data_buffer[i];
			channelDataSum = channelDataSum + channel_data[j];
			j++;
		}
		// find DC offset -- for long signals, it's approx. the mean
		channelDataMean = channelDataSum / ((int)nrecords);
	} else {
		int j = 0;
		for (int i = channelNumInt-1; i < (int)nrecords*(int)nchannels; i = i + nchannels) {
			channel_data[j] = data_buffer[i];
			j++;
		}
	}

	// Build time vector
	period_id = H5Dopen2(file_id, periodToRead.toLatin1().constData(), H5P_DEFAULT);
	if(period_id < 0)
		printf("Throw error - H5Dopen2 error %d\n", period_id);
	status = H5Dread(period_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, period_buffer);
	if(status < 0)
		printf("Throw error - H5Dread error %d\n", status);
	time_buffer[0] = 0;
	for (int i = 1; i < (int)nrecords; i++) {
		time_buffer[i] = time_buffer[i-1] + (*period_buffer / 1e9);
	}
	
	// FFT plot
	// TO-DO: check if FFT plot is selected
	double windowedChannelVal;
	int fft_length = 2;
	while (fft_length < (int)nrecords) {
		fft_length = fft_length * 2;
	}
	fft_input = (complex *)malloc(sizeof(complex)*fft_length);
	fft_buffer = (complex *)malloc(sizeof(complex)*fft_length);
	fft_output_y = (double *)malloc(sizeof(double)*fft_length);
	fft_output_x = (double *)malloc(sizeof(double)*fft_length);
	makeWindow((int)nrecords);
	for (int i = 0; i < fft_length; i++) {
		if (i < (int)nrecords) {
			windowedChannelVal = disc_window->GetDataWinCoeff(i) * channel_data[i];
			fft_input[i] = complex(windowedChannelVal, 0.0);
		} else {
			fft_input[i] = 0; // pad the rest of the array with 0's
		}
	}
	fft(fft_input, fft_buffer, fft_length);
	for (int i = 0; i < fft_length; i++) {
		fft_output_y[i] = abs(real(fft_buffer[i]));
		fft_output_x[i] = i / ((*period_buffer / 1e9) * fft_length);
	}
	
	// Full wave rectification
	if (fwrChecked) {
		// subtract DC offset from each value and take absolute value
		for (int i = 0; i < (int)nrecords; i++) {
			channel_data[i] = std::abs(channel_data[i] - channelDataMean);
		}
	}
	
	// Plot
	tscurve->setRawSamples(time_buffer, channel_data, (int)nrecords);
	fftcurve->setRawSamples(fft_output_x, fft_output_y, fft_length);
	// Refresh enabled plots (turn autoscaling on in case the zoom feature was used)
	tsplot->setAxisAutoScale(tsplot->yLeft, true);
	tsplot->setAxisAutoScale(tsplot->xBottom, true);
	tsplot->replot();
	fftplot->setAxisAutoScale(fftplot->yLeft, true);
	fftplot->setAxisAutoScale(fftplot->xBottom, true);
	fftplot->replot();

	// Close identifiers
	H5PTclose(packettable_id);
	H5Gclose(trial_id);
	H5Dclose(period_id);
}

// Temporary function for validating data access
void AnalysisTools::dump_vals(double *data, hsize_t *ndims)
{
	// Only printing first value out or else the printf will block
	for(size_t i=0; i<10; i++)
		printf("value is %f\n", data[i]);
}

void AnalysisTools::freePlotBuffers()
{
	if(data_buffer) {
		free(data_buffer);
		data_buffer = NULL;
	}
	if(channel_data) {
		free(channel_data);
		channel_data = NULL;
	}
	if(time_buffer) {
		free(time_buffer);
		time_buffer = NULL;
	}
	if(period_buffer) {
		free(period_buffer);
		period_buffer = NULL;
	}
	if(fft_buffer) {
		free(fft_buffer);
		fft_buffer = NULL;
	}
	if(fft_input) {
		free(fft_input);
		fft_input = NULL;
	}
	if(fft_output_x) {
		free(fft_output_x);
		fft_output_x = NULL;
	}
	if(fft_output_y) {
		free(fft_output_y);
		fft_output_y = NULL;
	}
}
