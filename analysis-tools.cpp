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

static DefaultGUIModel::variable_t vars[] = {
	{ "Input", "Input", DefaultGUIModel::INPUT, },
	{ "Output", "Output", DefaultGUIModel::OUTPUT, },
	{ "Time (s)", "Time (s)", DefaultGUIModel::STATE, }, 
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
	// TO-DO: put QGroupBox around each plot to make it look neater, will eliminate resizing issues too
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

	// Plot options
	// TO-DO: disable scatter and FFT plot / gray scatter and FFT screenshot buttons by default
	QGroupBox *optionBox = new QGroupBox("Plot Selection");
	QHBoxLayout *optionBoxLayout = new QHBoxLayout;
	optionBox->setLayout(optionBoxLayout);
	QButtonGroup *optionButtons = new QButtonGroup;
	optionButtons->setExclusive(false);
	QCheckBox *plotTSCheckBox = new QCheckBox("Time Series");
	optionBoxLayout->addWidget(plotTSCheckBox);
	optionButtons->addButton(plotTSCheckBox);
	plotTSCheckBox->setChecked(true);
	plotTSCheckBox->setEnabled(true);
	QCheckBox *plotScatterCheckBox = new QCheckBox("Scatter");
	optionBoxLayout->addWidget(plotScatterCheckBox);
	optionButtons->addButton(plotScatterCheckBox);
	plotScatterCheckBox->setChecked(true);
	plotScatterCheckBox->setEnabled(true);
	QCheckBox *plotFFTCheckBox = new QCheckBox("FFT");
	optionBoxLayout->addWidget(plotFFTCheckBox);
	optionButtons->addButton(plotFFTCheckBox);
	plotFFTCheckBox->setChecked(true);
	plotFFTCheckBox->setEnabled(true);
	QObject::connect(plotTSCheckBox,SIGNAL(toggled(bool)),tsplot,SLOT(setEnabled(bool)));
	QObject::connect(plotTSCheckBox,SIGNAL(toggled(bool)),saveTSPlotButton,SLOT(setEnabled(bool)));
	QObject::connect(plotTSCheckBox,SIGNAL(toggled(bool)),this,SLOT(toggleTSplot(bool)));
	QObject::connect(plotScatterCheckBox,SIGNAL(toggled(bool)),splot,SLOT(setEnabled(bool)));
	QObject::connect(plotScatterCheckBox,SIGNAL(toggled(bool)),saveScatterPlotButton,SLOT(setEnabled(bool)));
	QObject::connect(plotScatterCheckBox,SIGNAL(toggled(bool)),this,SLOT(toggleScatterplot(bool)));
	QObject::connect(plotFFTCheckBox,SIGNAL(toggled(bool)),fftplot,SLOT(setEnabled(bool)));
	QObject::connect(plotFFTCheckBox,SIGNAL(toggled(bool)),saveFFTPlotButton,SLOT(setEnabled(bool)));
	QObject::connect(plotFFTCheckBox,SIGNAL(toggled(bool)),this,SLOT(toggleFFTplot(bool)));
	plotTSCheckBox->setToolTip("Enable time series plot");
	plotScatterCheckBox->setToolTip("Enable scatter plot");
	plotFFTCheckBox->setToolTip("Enable FFT plot");
	customlayout->addWidget(optionBox, 0, 0, 1, 1);
	
	// File control
	QGroupBox *fileBox = new QGroupBox(tr("File Control"));
	QHBoxLayout *fileLayout = new QHBoxLayout;
	fileBox->setLayout(fileLayout);
	fileLayout->addWidget(new QLabel(tr("File Name")));
	QLineEdit *fileNameEdit = new QLineEdit;
	fileNameEdit->setReadOnly(true);
	fileLayout->addWidget(fileNameEdit);
	QPushButton *fileChangeButton = new QPushButton("Choose File");
	fileLayout->addWidget(fileChangeButton);
	QObject::connect(fileChangeButton,SIGNAL(released(void)),this,SLOT(changeDataFile(void)));
	customlayout->addWidget(fileBox, 2, 0, 1, 1);

	// HDF5 viewer
	QGroupBox *viewBox = new QGroupBox(tr("HDF5 Viewer"));
	QHBoxLayout *viewLayout = new QHBoxLayout;
	viewBox->setLayout(viewLayout);
	customlayout->addWidget(viewBox, 3, 0, 4, 1);
	
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

void AnalysisTools::toggleTSplot(bool on) {
}

void AnalysisTools::toggleScatterplot(bool on) {
}

void AnalysisTools::toggleFFTplot(bool on) {
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
	filterList.push_back("All files (*.*)");
	fileDialog.setFilters(filterList);
	fileDialog.selectNameFilter("HDF5 files (*.h5)");

	QStringList files;
	if(fileDialog.exec())
		files = fileDialog.selectedFiles();

	QString filename;
	if(files.isEmpty() || files[0] == NULL || files[0] == "/" )
		return;
	else
		filename = files[0];

	if (!filename.toLower().endsWith(QString(".h5")))
		filename += ".h5";

	// Write this directory to the user prefs as most recently used
	//userprefs.setValue("/dirs/data", fileDialog.directory().path());

	// Post to event queue
	// TO-DO: is this necessary? Taken from data recorder
	//OpenFileEvent RTevent(filename, fifo);
	//RT::System::getInstance()->postEvent(&RTevent);
}

int AnalysisTools::openFile(QString &filename) {
//#ifdef DEBUG
	//if(!pthread_equal(pthread_self(),thread))
	//{
		//ERROR_MSG("DataRecorder::Panel::openFile : called by invalid thread\n");
		//PRINT_BACKTRACE();
	//}
//#endif

	//if (QFile::exists(filename)) {
		//CustomEvent *event = new CustomEvent(static_cast<QEvent::Type>QFileExistsEvent);
		//FileExistsEventData data;

		//event->setData(static_cast<void *>(&data));
		//data.filename = filename;

		//QApplication::postEvent(this, event);
		//data.done.wait(&mutex);
		//if (data.response == 0) { // append
			//file.id = H5Fopen(filename.toLatin1().constData(), H5F_ACC_RDWR, H5P_DEFAULT);
			//size_t trial_num;
			//QString trial_name;
			//H5Eset_auto(H5E_DEFAULT, NULL, NULL);
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
		//} else if (data.response == 1) { //overwrite
			//file.id = H5Fcreate(filename.toLatin1().constData(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
			//trialNum->setText("0");
		//} else {
			//return -1;
		//}
	//} else {
		//file.id = H5Fcreate(filename.toLatin1().constData(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		//trialNum->setText("0");
	//}
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

//bool AnalysisTools::OpenFile(QString FName) {
	//return true;
//}
