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

// global variables (necessary for op_func)
QTreeWidget *treeViewer;
QTreeWidgetItem *treeParent;
QTreeWidgetItem *treeChild1;
QTreeWidgetItem *treeChild2;
QString currentTrial;
int currentTrialFlag;
int firstChannelSelected;

extern "C" Plugin::Object *createRTXIPlugin(void) {
	return new AnalysisTools();
}

static DefaultGUIModel::variable_t vars[] = {
	{ "Input", "Input", DefaultGUIModel::INPUT, },
	{ "Output", "Output", DefaultGUIModel::OUTPUT, }, 
};

static size_t num_vars = sizeof(vars) / sizeof(DefaultGUIModel::variable_t);

AnalysisTools::AnalysisTools(void) :  DefaultGUIModel("Analysis Tools", ::vars, ::num_vars) {
	setWhatsThis(
			"<p><b>Analysis Tools</b></p><p>Analysis tools</p>"); // TO-DO: add detail here
	DefaultGUIModel::createGUI(vars, num_vars); // this is required to create the GUI
	customizeGUI();
	update(INIT);
	refresh(); // this is required to update the GUI with parameter and state values
	QTimer::singleShot(0, this, SLOT(resizeMe()));
}

AnalysisTools::~AnalysisTools(void) {}

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
	fftplot->setFixedSize(450, 270);
	fftplotBoxLayout->addWidget(fftplot);
	customlayout->addWidget(fftplotBox, 4, 4, 3, 2);
	// Connect screenshot buttons to functions
	QObject::connect(saveTSPlotButton, SIGNAL(clicked()), this, SLOT(screenshotTS()));
	QObject::connect(saveScatterPlotButton, SIGNAL(clicked()), this, SLOT(screenshotScatter()));
	QObject::connect(saveFFTPlotButton, SIGNAL(clicked()), this, SLOT(screenshotFFT()));

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
	QObject::connect(plotTSCheckBox,SIGNAL(toggled(bool)),tsplot,SLOT(setEnabled(bool)));
	QObject::connect(plotTSCheckBox,SIGNAL(toggled(bool)),saveTSPlotButton,SLOT(setEnabled(bool)));
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

	// Scatter/FFT plot options
	QGroupBox *plotOptionsBox = new QGroupBox(tr("Scatter/FFT Plot Options"));
	// TO-DO: add detail here (later)
	customlayout->addWidget(plotOptionsBox, 2, 0, 1, 1);

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
	customlayout->addWidget(treeViewer, 3, 0, 4, 1);

	// Attributes
	QGroupBox *attributesBox = new QGroupBox(tr("Attributes"));
	QHBoxLayout *attributesLayout = new QHBoxLayout;
	attributesBox->setLayout(attributesLayout);
	customlayout->addWidget(attributesBox, 4, 2, 3, 1);

	// Parameters
	QGroupBox *paramsBox = new QGroupBox(tr("Parameters"));
	QHBoxLayout *paramsLayout = new QHBoxLayout;
	paramsBox->setLayout(paramsLayout);
	customlayout->addWidget(paramsBox, 4, 3, 3, 1);

	// Standard module buttons
	QObject::connect(DefaultGUIModel::pauseButton,SIGNAL(toggled(bool)),DefaultGUIModel::modifyButton,SLOT(setEnabled(bool)));
	DefaultGUIModel::pauseButton->setToolTip("Start/Stop analysis tools");
	DefaultGUIModel::modifyButton->setToolTip("Commit changes to parameter values");
	DefaultGUIModel::unloadButton->setToolTip("Close module");

	setLayout(customlayout);
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

// TO-DO: may need to restore toggle functions to allow plots to be cleared when deselected

void AnalysisTools::changeDataFile(void) {
	QFileDialog fileDialog(this);
	fileDialog.setFileMode(QFileDialog::AnyFile);
	fileDialog.setWindowTitle("Select Data File");

	QSettings userprefs;
	userprefs.setPath(QSettings::NativeFormat, QSettings::SystemScope, "/usr/local/share/rtxi/");
	fileDialog.setDirectory(userprefs.value("/dirs/data", getenv("HOME")).toString());

	QStringList filterList;
	filterList.push_back("HDF5 files (*.h5)");
	fileDialog.setFilters(filterList);
	fileDialog.selectNameFilter("HDF5 files (*.h5)");

	QStringList files;
	QString filename;
	if(fileDialog.exec()) {
		files = fileDialog.selectedFiles();
		filename = files[0];
		fileNameEdit->setText(filename);
		openFile(filename);
	}
}

// TO-DO: populate HDF5, attribute, and parameter viewer contents
//        enable any scatter/FFT specific options
//        call closeFile if fileLoaded flag is set to 1
int AnalysisTools::openFile(QString &filename) {
	herr_t status = -1;
	currentTrialFlag = 0;
	currentTrial = "";
	firstChannelSelected = 0;
	if (QFile::exists(filename)) {
		file_id = H5Fopen(filename.toLatin1().constData(), H5F_ACC_RDONLY, H5P_DEFAULT);

		// Iterate through file
		//printf ("Objects in the file:\n");
		status = H5Ovisit(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, op_func, NULL);

		if (!status) plotButton->setEnabled(true);
	}

	return status;
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
					treeChild1->setText(0, qName);
					treeParent->addChild(treeChild1);
				} else if (qName.startsWith(currentTrial) && qName.endsWith("Synchronous Data")) {
					currentTrialFlag = 0; // reset flag to start parsing next trial
					treeChild1 = new QTreeWidgetItem;
					treeChild1->setText(0, qName);
					treeParent->addChild(treeChild1);
				}
				break;
			case H5O_TYPE_DATASET:
				//printf ("%s  (Dataset)\n", name);
				if (qName.startsWith(currentTrial + "/Asynchronous Data") && !qName.endsWith("Channel Data")) {
					treeChild2 = new QTreeWidgetItem;
					treeChild2->setText(0, qName);
					treeChild1->addChild(treeChild2);
					treeChild2->setToolTip(0, qName);
					if (!firstChannelSelected) {
						firstChannelSelected = 1;
						treeViewer->setCurrentItem(treeChild2);
					}
				} else if (qName.startsWith(currentTrial + "/Synchronous Data") && !qName.endsWith("Channel Data")) {
					treeChild2 = new QTreeWidgetItem;
					treeChild2->setText(0, qName);
					treeChild1->addChild(treeChild2);
					treeChild2->setToolTip(0, qName);
					if (!firstChannelSelected) {
						firstChannelSelected = 1;
						treeViewer->setCurrentItem(treeChild2);
					}
				}
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

// TO-DO: erase HDF5, attribute, and parameter viewer contents
//        disable plot button and any scatter/FFT specific options
void AnalysisTools::closeFile(bool shutdown)
{
	//#ifdef DEBUG
	//if(!pthread_equal(pthread_self(),thread))
	//{
	//ERROR_MSG("DataRecorder::Panel::closeFile : called by invalid thread\n");
	//PRINT_BACKTRACE();
	//}
	//#endif

	//H5Fclose(file.id);
	//if (!shutdown) {
	//CustomEvent *event = new CustomEvent(static_cast<QEvent::Type>QSetFileNameEditEvent);
	//SetFileNameEditEventData data;
	//data.filename = "";
	//event->setData(static_cast<void*>(&data));

	//QApplication::postEvent(this, event);
	//data.done.wait(&mutex);
	//}
}

// TO-DO: think through error cases here (e.g. when one of the top-level groups are selected, etc.)
void AnalysisTools::plotTrial() {
	// TO-DO: check that current item is a dataset (and not a group), only open/plot if a dataset is selected (maybe display an warning otherwise?)
	//        need to open appropriate column in Channel Data, not the header dataset
	//        plot using QWT and setData/setSamples -- look into their assocated warnings
	//        only plot if check boxes are selected
	herr_t status;
	hsize_t dims[2];
	hid_t dataset_id, datatype_id, dataspace_id, memspace_id;
	H5T_class_t class_id;
	int rank, dim_status;
	double *data_buffer;

	// Get elements from GUI
	QString selectedTrial = treeViewer->currentItem()->text(0);
	QString channelNum = selectedTrial.at(selectedTrial.size()-1);
	QString channelToRead = treeViewer->currentItem()->parent()->text(0) + "/Channel Data";

	// Open dataset
	dataset_id = H5Dopen2(file_id, channelToRead.toLatin1().constData(), H5P_DEFAULT);
	if(dataset_id < 0)
		printf("Throw error - H5Dopen2 error %d\n", dataset_id);
	
	// Get dataset info
	datatype_id = H5Dget_type(dataset_id);
	if(datatype_id < 0)
		printf("Throw error - H5Dget_type %d\n", datatype_id);

	// Get dataset class
	class_id = H5Tget_class(datatype_id);
	if(class_id < 0)
		printf("Throw error - H5Tget_class %d\n", class_id);

	// Get dataset space in bytes
	dataspace_id = H5Dget_space(dataset_id);
	if(dataspace_id < 0)
		printf("Throw error - H5Dget_space error %d\n", dataspace_id);

	// Get rank of dataspace
	rank = H5Sget_simple_extent_ndims(dataspace_id);
	if(rank < 0)
		printf("Throw error - H5Sget_simple_extent_ndims error %d\n", rank);

	// Get dimensions of data
	// If the second dimension is zero, then the dataset is a scalar
	dim_status = H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
	if(dim_status < 0)
		printf("Throw error - H5Sget_simple_extent_dims error %d\n", dim_status);

	// Define memory dataspace
	memspace_id = H5Screate_simple(rank, dims, NULL);
	if(memspace_id < 0)
		printf("Throw error - H5Screate_simple error %d\n", memspace_id);

	// Print for debug
	printf("dataset_id: %d\n" "datatype_id: %d\n" "dataspace_id: %d\n" "memspace_id: %d\n" "rank: %d\n" "dimensions: %lu x %lu \n\n",
			dataset_id, datatype_id, dataspace_id, memspace_id, rank, (unsigned long)(dims[0]), (unsigned long)(dims[1]));

	// Initialize data buffer using rank 
	if(rank == 1)
		data_buffer = (double *)malloc(sizeof(double)*dims[0]);
	else if(rank == 2)
		data_buffer = (double *)malloc(sizeof(double)*dims[0]*dims[1]);

	// Read data
	status = H5Dread(dataset_id, datatype_id, memspace_id, dataspace_id, H5P_DEFAULT, data_buffer);
	if(status < 0)
		printf("Throw error - H5Dread error %d\n", status);

	// Temporarily use print to check if data is really there
	dump_vals(data_buffer, rank, dims);

	// The rest is easier
	// TO-DO: read time duration and period, build time vector for plotting

	// TO-DO: plot selected trial and channel
	//tsplot->setData(

	// Close everything
	H5Dclose(dataset_id);
	H5Tclose(datatype_id);
	H5Sclose(memspace_id);
	H5Sclose(dataspace_id);
}

// Temporary function for validating data access
void AnalysisTools::dump_vals(double *data, int rank, hsize_t *ndims)
{
	printf("Rank: %d\n"
			"Dimensions: %lu x %lu \n\n", 
			rank, (unsigned long)(ndims[0]),
			(unsigned long)(ndims[1]));

	for(size_t i=0; i<ndims[0]; i++)
		if(ndims[1] == 0)
			printf("value is %f\n", data[i]);
		else
			for(size_t j=0; j<ndims[1]; j++)
				printf("value is %f\n", data[i*ndims[1]+j]);
}
