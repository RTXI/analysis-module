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

#include <QMessageBox>
#include <QFileInfo>

#include <DSP/gen_win.h>
#include <DSP/rectnglr.h>
#include <DSP/trianglr.h>
#include <DSP/hamming.h>
#include <DSP/hann.h>
#include <DSP/dolph.h>
#include <DSP/kaiser.h>
#include <DSP/fft.h>

#include <qwt_plot_curve.h>
#include <default_gui_model.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <scatterplot.h>

namespace HdfViewerUtils {

class HdfViewer : public DefaultGUIModel {

	Q_OBJECT

	public:
		HdfViewer(void);
		virtual ~HdfViewer(void);
		void execute(void);
		void customizeGUI(void);
		
		enum window_t	{
			RECT=0, TRI, HAMM, HANN, CHEBY, KAISER
		};
		enum plot_t	{
			TIMESERIES, SCATTER, FFT
		};

	public slots:

	signals: // custom signals

	protected:
		virtual void update(DefaultGUIModel::update_flags_t);

	private:
		// inputs, states, flags, calculated values
		bool fwrChecked;
		GenericWindow *disc_window;
		plot_t plot_mode;
		window_t window_shape;
		double Kalpha; // Kaiser window sidelobe attenuation parameter
		double Calpha; // Chebyshev window sidelobe attenuation parameter
		
		// File I/O
		hid_t file_id;
		hid_t dataset_id;
		
		// GUI components
		BasicPlot *omniplot;
		QwtPlotCurve *tscurve, *fftcurve;
		// QwtText xAxisTitle, yAxisTitle;

		QGroupBox *plotControls, *plotOptions;
		QLineEdit *fileNameEdit;
		QPushButton *resetPlotButton, *savePlotButton, 
		            *exportSeriesButton;
		QComboBox *windowShape, *plotType;
		QButtonGroup *plotOptionsButtons;
		QCheckBox *FWRCheckBox;
		QDoubleSpinBox *kalphaEdit, *calphaEdit;

		// custom functions
		void initParameters(void);
		int openFile(QString &filename);
		void closeFile(void);
		void dump_vals(double *, hsize_t*);
		void makeWindow(int);
		void freePlotBuffers(void);

	private slots:
		void getTrialData(void);
		void changeChannel(QModelIndex);
		void updatePlotMode(int);
		void updatePlot(void);
		void screenshot(void);
		void clearData(void);
		void toggleFWR(bool);
		void resetAxes(void);
		void changeDataFile(void);
		void updateWindow(int);
		void updateKalpha(double);
		void updateCalpha(double);
};
};
