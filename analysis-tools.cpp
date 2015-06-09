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
QTreeWidgetItem *treeChild;
QString currentTrial;
int currentTrialFlag;
int firstChannelSelected;

extern "C" Plugin::Object *createRTXIPlugin(void) {
	return new AnalysisTools();
}

static DefaultGUIModel::variable_t vars[] = {
	{ "Input", "Input", DefaultGUIModel::INPUT, },
	{ "Output", "Output", DefaultGUIModel::OUTPUT, },
	//{ "Time (s)", "Time (s)", DefaultGUIModel::STATE, }, 
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

AnalysisTools::~AnalysisTools(void) {}

void AnalysisTools::execute(void) {
	systime = count * dt; // current time, s

	count++; // increment count to measure time
	return;
}

void AnalysisTools::update(DefaultGUIModel::update_flags_t flag) {
	switch (flag) {
	case INIT:
		setState("Time (s)", systime);
		break;
	case MODIFY:
		bookkeep();
		break;
	case UNPAUSE:
		bookkeep();
		break;
	case PAUSE:
		output(0) = 0; // stop command in case pause occurs in the middle of command
		break;
	case PERIOD:
		dt = RT::System::getInstance()->getPeriod() * 1e-9;
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
	// TO-DO: gray out plot button before file is loaded
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

void AnalysisTools::initParameters() {
	dt = RT::System::getInstance()->getPeriod() * 1e-9; // s
	bookkeep();
}

void AnalysisTools::bookkeep() {
	count = 0;
	systime = 0;
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
//        update and enable channelSelection (channelSelection->insertItem(1, "Ch0");)
//        update and enable trialSelection (trialSelection->insertItem(1, "trial0");)
//        enable any scatter/FFT specific options
int AnalysisTools::openFile(QString &filename) {
	currentTrialFlag = 0;
	currentTrial = "";
	firstChannelSelected = 0;
	if (QFile::exists(filename)) {
		hid_t file_id;
		herr_t status;
		file_id = H5Fopen(filename.toLatin1().constData(), H5F_ACC_RDONLY, H5P_DEFAULT);
		
		// Iterate through file
		//printf ("Objects in the file:\n");
		status = H5Ovisit(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, op_func, NULL);
		
		// TO-DO: if file opened successfully, enable plot button and select the first channel of first trial (to block that potential error case)
		if (!status) {
			plotButton->setEnabled(true);
		}
			
		
		// TO-DO: remove commented code
		//treeViewer->addTopLevelItems(treeList);
		
		//size_t trial_num;
		//QString trial_name;
		//for (trial_num = 1;; ++trial_num) {
			//trial_name = "/Trial" + QString::number(trial_num);
			//file.trial = H5Gopen(file.id, trial_name.toLatin1().constData(), H5P_DEFAULT);
			//if (file.trial < 0) {
				//H5Eclear(H5E_DEFAULT);
				//break;
			//} else
				//H5Gclose(file.trial);
		//}
		//trialNum->setNum(int(trial_num)-1);
	}
	//if (file.id < 0) {
		//H5E_type_t error_type;
		//size_t error_size;
		//error_size = H5Eget_msg(file.id, &error_type, NULL, 0);
		//char error_msg[error_size + 1];
		//H5Eget_msg(file.id, &error_type, error_msg, error_size);
		//error_msg[error_size] = 0;
		//H5Eclear(file.id);

		//ERROR_MSG("DataRecorder::Panel::processData : failed to open \"%s\" for writing with error : %s\n", filename.toStdString().c_str(),error_msg);
		//return -1;
	//}

	//CustomEvent *event = new CustomEvent(static_cast<QEvent::Type>QSetFileNameEditEvent);
	//SetFileNameEditEventData data;
	//data.filename = filename;
	//event->setData(static_cast<void*>(&data));

	//QApplication::postEvent(this, event);
	//data.done.wait(&mutex);

	return 0;
}

herr_t op_func(hid_t loc_id, const char *name, const H5O_info_t *info, void *operator_data) {
	QString qName = QString(name);
    //printf ("/"); // Print root group in object path
    if (name[0] == '.') { // Root group, do not print '.'
        //printf ("  (Group)\n");
    } else {
        switch (info->type) {
            case H5O_TYPE_GROUP:
                //printf ("%s  (Group)\n", name);
                if (!currentTrialFlag) {
					currentTrialFlag = 1;
					currentTrial = qName;
				} else if (qName.startsWith(currentTrial) && qName.endsWith("Synchronous Data")) {
					currentTrialFlag = 0;
					treeParent = new QTreeWidgetItem;
					treeParent->setText(0, qName);
					treeViewer->addTopLevelItem(treeParent);
				}
                break;
            case H5O_TYPE_DATASET: // TO-DO: reset synchDataFlag
                //printf ("%s  (Dataset)\n", name);
                if (qName.startsWith(currentTrial + "/Synchronous Data") && !qName.endsWith("Channel Data")) {
					//printf ("Found synchronous dataset\n");
					treeChild = new QTreeWidgetItem;
					treeChild->setText(0, qName);
					treeParent->addChild(treeChild);
					if (!firstChannelSelected) {
						firstChannelSelected = 1;
						treeViewer->setCurrentItem(treeChild);
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
//        erase contents of channelSelection and trialSelection
//        disable channelSelection, trialSelection and any scatter/FFT specific options
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

void AnalysisTools::plotTrial() {
}
