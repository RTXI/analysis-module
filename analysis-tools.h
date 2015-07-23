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

#include <QtGui>
#include <QMessageBox>
#include <QFileInfo>

#include <qwt_plot_curve.h>
#include <default_gui_model.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <scatterplot.h>

class AnalysisTools : public DefaultGUIModel {

	Q_OBJECT

	public:
		AnalysisTools(void);
		virtual ~AnalysisTools(void);
		void execute(void);
		void customizeGUI(void);

	public slots:

	signals: // custom signals

	protected:
		virtual void update(DefaultGUIModel::update_flags_t);

	private:
		// inputs, states, flags, calculated values
		bool fwrChecked;
		
		// File I/O
		hid_t file_id;
		hid_t dataset_id;
		
		// GUI components
		BasicPlot *tsplot;
		QwtPlotCurve *tscurve;
		ScatterPlot *splot;
		BasicPlot *fftplot;
		QLineEdit *fileNameEdit;
		QPushButton *plotButton;

		// custom functions
		void initParameters(void);
		int openFile(QString &filename);
		void closeFile(void);
		void dump_vals(double *, hsize_t*);

	private slots:
		void plotTrial(void);
		void screenshotTS(void);
		void screenshotScatter(void);
		void screenshotFFT(void);
		void clearData(void);
		void toggleFWR(bool);
		void changeDataFile(void);
};

herr_t op_func(hid_t, const char*, const H5O_info_t*, void*);
