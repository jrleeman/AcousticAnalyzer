"""


Many concepts modeled after: http://eli.thegreenplace.net/files/prog_code/qt_mpl_bars.py.txt
"""

import sys, os, random
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from datetime import datetime
import numpy as np
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

class AppForm(QMainWindow):

    def __init__(self,parent=None):
        QMainWindow.__init__(self,parent)
        self.setWindowTitle('Acoustic Analyzer')
        
        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()
        
        self.on_draw()
        
        self.first_plot = True
        
    def save_plot(self):
        file_choices = "PNG (*.png)|*.png"
        
        path = unicode(QFileDialog.getSaveFileName(self,'Save file', '', file_choices))
        
        if path:
            self.canvas.print_figure(path,dpi=self.dpi)
            self.statusBar().showMessage('Saved to %s' %path, 2000)
            
    def draw_waveform(self,listname):
        self.stored_xlimits = self.ax1.get_xlim()
        self.stored_ylimits = self.ax1.get_ylim()
        self.ax1.clear()
        self.ax1.set_xlabel(r'Time [$\mu s$]')
        self.ax1.set_ylabel(r'Amplitude [V]')
        self.ax1.grid(self.grid_cb.isChecked())
        seismogram = self.experiment.seismograms[str(listname)]
        self.ax1.plot(self.experiment.time,seismogram.amplitude,color='k',linewidth=1)
        
        if self.first_plot == True:
            self.first_plot = False
        else:
            self.ax1.set_xlim(self.stored_xlimits)
            self.ax1.set_ylim(self.stored_ylimits)
            
        self.trigger_show()
        self.points_show()
        
        self.canvas.draw()
    
    def points_show(self):
        seismogram = self.experiment.seismograms[str(self.loaded_list.currentItem().text())]
        if self.show_points_cb.isChecked():
            self.scatter_points = self.ax1.scatter(self.experiment.time,seismogram.amplitude,color='k',s=10)
        else:
            try:
                self.scatter_points.remove()
            except:
                pass
        self.canvas.draw()
        
    def trigger_show(self):
        if self.trigger_cb.isChecked():
            self.trigger_line = self.ax1.axvline(x=self.experiment.trigger_time,color='r')
        else:
            try:
                self.trigger_line.remove()
            except:
                pass
        self.canvas.draw()
    
    def on_trigger_zero(self):
        self.stored_xlimits = self.ax1.get_xlim()
        if self.zero_trigger_cb.isChecked():
            self.experiment.time = self.experiment.time - self.experiment.trigger_time
            self.stored_xlimits = [item - self.experiment.trigger_time for item in self.stored_xlimits]
            self.experiment.trigger_time = 0.0
            self.ax1.set_xlim(self.stored_xlimits)
            
            
        else:
            self.experiment.time = np.arange(self.experiment.total_recs)/float(self.experiment.trigger_rate)
            self.experiment.trigger_time = (self.experiment.pre_trig_recs) * (1/float(self.experiment.trigger_rate))
            self.stored_xlimits = [item + self.experiment.trigger_time for item in self.stored_xlimits]
            self.ax1.set_xlim(self.stored_xlimits)
        self.draw_waveform(self.loaded_list.currentItem().text())
            
    def on_about(self):
        msg = """ ENTER AN ABOUT MESSAGE HERE
        IT CAN EVEN BE MULTI-LINE
        """
        QMessageBox.about(self, "About Acoustic Analyzer", msg.strip())
        
        
    def on_pick(self):
        box_points = event.artist.get_bbox().get_points()
        msg = "Coordinates: %s" %box_points
        QMessageBox.information(self, "Marked Cursor", msg)
        
    def on_draw(self):
        self.ax1.grid(self.grid_cb.isChecked())
        self.canvas.draw()
        
    def on_clear(self):
        self.ax1.clear()
        self.ax1.grid(self.grid_cb.isChecked())
        self.canvas.draw()
        
    def create_meta_box(self):
        ### META BOX
        # Make the boxes we'll use and do some initial setup
        meta_left_vbox = QVBoxLayout()
        meta_right_vbox = QVBoxLayout()
        meta_hbox = QHBoxLayout()
        meta_hbox.addLayout(meta_left_vbox)
        meta_hbox.addSpacing(300)
        meta_hbox.addLayout(meta_right_vbox)
        meta_hbox.addSpacing(300)
        
        # Make the experiment meta label objects
        self.exp_date_label = QLabel(self)
        self.exp_recrate_label = QLabel(self)
        self.exp_pulserate_label = QLabel(self)
        self.exp_pretrigrecs_label = QLabel(self)
        self.exp_posttrigrecs_label = QLabel(self)
        self.exp_totalrecs_label = QLabel(self)
        # Add them to the left vbox
        meta_left_vbox.addWidget(self.exp_date_label)
        meta_left_vbox.addWidget(self.exp_recrate_label)
        meta_left_vbox.addWidget(self.exp_pulserate_label)
        meta_left_vbox.addWidget(self.exp_pretrigrecs_label)
        meta_left_vbox.addWidget(self.exp_posttrigrecs_label)
        meta_left_vbox.addWidget(self.exp_totalrecs_label)
        # Set default labels
        self.exp_date_label.setText('Date:')
        self.exp_recrate_label.setText('Recording Rate [MHz]:')
        self.exp_pulserate_label.setText('Pulse Rate [pps]:')
        self.exp_posttrigrecs_label.setText('Post-Trigger Samples:')
        self.exp_pretrigrecs_label.setText('Pre-Trigger Samples:')
        self.exp_totalrecs_label.setText('Total Samples:')
        
        # Make the seismogram meta label objects
        self.seismogram_recnum_label = QLabel(self)
        self.seismogram_syncvoltage_label = QLabel(self)
        self.seismogram_maxamplitude_label = QLabel(self)
        self.seismogram_timestamp_label = QLabel(self)
        # Add them to the right vbox
        meta_right_vbox.addWidget(self.seismogram_recnum_label)
        meta_right_vbox.addWidget(self.seismogram_syncvoltage_label)
        meta_right_vbox.addWidget(self.seismogram_maxamplitude_label)
        meta_right_vbox.addWidget(self.seismogram_timestamp_label)
        # Set default labels
        self.seismogram_recnum_label.setText('Record Number: ')
        self.seismogram_syncvoltage_label.setText('Sync Voltage: ')
        self.seismogram_maxamplitude_label.setText('Maximum Amplitude: ')
        self.seismogram_timestamp_label.setText('Time Stamp: ')
        
                

        
        return meta_hbox
        
    def create_main_frame(self):
        self.main_frame = QWidget()
        
        #Create the plotting area
        self.dpi = 100
        self.fig = Figure((10.0,3.0), dpi=self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)
        self.ax1 = self.fig.add_subplot(111)
        self.fig.subplots_adjust(bottom=0.2)
        
        # Bind the picker event
        self.canvas.mpl_connect('pick_event', self.on_pick)
        
        # Create and bind the toolbar for the plot
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        
        # Other GUI Cotrols will go here
        self.loaded_list = QListWidget()
        
        self.grid_cb = QCheckBox("Show &Grid")
        self.grid_cb.setChecked(False)
        self.connect(self.grid_cb, SIGNAL('stateChanged(int)'), self.on_draw)
        
        self.trigger_cb = QCheckBox("Show &Trigger")
        self.trigger_cb.setChecked(False)
        self.connect(self.trigger_cb, SIGNAL('stateChanged(int)'), self.trigger_show)
        
        self.zero_trigger_cb = QCheckBox("Zero &Trigger")
        self.zero_trigger_cb.setChecked(False)
        self.connect(self.zero_trigger_cb, SIGNAL('stateChanged(int)'), self.on_trigger_zero)
        
        self.show_points_cb = QCheckBox("Show &Points")
        self.show_points_cb.setChecked(False)
        self.connect(self.show_points_cb, SIGNAL('stateChanged(int)'), self.points_show)
        
        self.plot_multiple_button = QPushButton("&Plot Multiple")
        #self.connect(self.force_draw_button, SIGNAL('clicked()'), self.on_force_draw)
        
        self.force_clear_button = QPushButton("&Force Clear")
        self.connect(self.force_clear_button, SIGNAL('clicked()'), self.on_clear)

        self.loaded_list.currentTextChanged.connect(self.draw_waveform)
        
        # Do layout with box sizers
        
        # Make the left Vbox that contains the loaded list
        left_hbox = QVBoxLayout()
        left_hbox.addWidget(self.loaded_list)
        
        # Make the plot Hbox that contains the plotting and plot toolbar
        plot_vbox = QVBoxLayout()
        plot_vbox.addWidget(self.canvas)
        plot_vbox.addWidget(self.mpl_toolbar)
        
        central_vbox = QVBoxLayout()
        options_hbox = QHBoxLayout()
        
        
        meta_hbox = self.create_meta_box()
        
        # Plot options vbox
        plot_options_vbox = QVBoxLayout()
        plot_options_vbox.addWidget(self.plot_multiple_button)
        plot_options_vbox.addWidget(self.force_clear_button)
        plot_options_vbox.addWidget(self.grid_cb)
        plot_options_vbox.addWidget(self.trigger_cb)
        plot_options_vbox.addWidget(self.zero_trigger_cb)
        plot_options_vbox.addWidget(self.show_points_cb)
        
        options_hbox.addLayout(plot_options_vbox)
        options_hbox.addStretch(1)
        options_hbox.addLayout(meta_hbox)
        
        central_vbox.addLayout(plot_vbox)
        central_vbox.addLayout(options_hbox)
        
        # Make the main box
        main_box = QHBoxLayout()
        main_box.addLayout(left_hbox)
        main_box.addLayout(central_vbox)
        self.main_frame.setLayout(main_box)
        self.setCentralWidget(self.main_frame)

    def create_status_bar(self):
        self.status_text = QLabel("Ready")
        self.statusBar().addWidget(self.status_text, 1)
        
    def create_menu(self): 
    
        #
        # File menu
        #
           
        self.file_menu = self.menuBar().addMenu("&File")
        
        load_file_action = self.create_action("&Save plot",
            shortcut="Ctrl+S", slot=self.save_plot, 
            tip="Save the plot")
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, 
            (load_file_action, None, quit_action))
        
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About", 
            shortcut='F1', slot=self.on_about, 
            tip='About the demo')
        
        self.add_actions(self.help_menu, (about_action,))
        
        #
        # Analyze Menu
        #
        self.analyze_menu = self.menuBar().addMenu("&Analyze")
        
        remove_mean_action = self.create_action("&Remove Mean (all)", shortcut="", slot=self.save_plot,tip="Remove the Mean of all the seismogram")
        
        detrend_action = self.create_action("&Detrend", shortcut="", slot=self.save_plot,tip="Linearly detrend the seismogram")
        
        calculate_velocities_action = self.create_action("&Calculate Velocities", shortcut="", slot=self.save_plot,tip="Calculate velocities")
        
        manual_picker_action = self.create_action("&Manual Picker", shortcut="", slot=self.save_plot,tip="Use manual picking")
        
        self.add_actions(self.analyze_menu,
            (remove_mean_action, detrend_action,calculate_velocities_action,manual_picker_action,))
        
        #
        # Cross-Correlation Sub-Menu
        #
        self.analyze_menu = self.analyze_menu.addMenu("&Cross-Correlation")
        pick_pattern_action = self.create_action("&Pick Pattern", shortcut="", slot=self.save_plot,tip="Pick the pattern used for cross-correlation")
        
        run_cross_correlation_action = self.create_action("&Run Cross-Correlation", shortcut="", slot=self.save_plot,tip="Run the Cross-Correlation Algorithm")
        
        self.add_actions(self.analyze_menu,
            (pick_pattern_action,run_cross_correlation_action))
            
        #
        # Data Menu
        #
        self.data_menu = self.menuBar().addMenu("&Data")
        load_experiment_action = self.create_action("&Load Experiment", shortcut="", slot=self.load_experiment,tip="Load an experiment from a log file.")
        
        load_waveform_action = self.create_action("&Load Waveform", shortcut="", slot=self.save_plot,tip="Load an individual waveform")
        
        set_calibrations_action = self.create_action("&Set Calibrations", shortcut="", slot=self.save_plot,tip="Set correction and calibration factors.")
        
        export_seismogram_data_action = self.create_action("&Export Seismogram Data", shortcut="", slot=self.save_plot,tip="Export the meta-data for all seismograms.")
        
        self.add_actions(self.data_menu,
            (load_experiment_action,load_waveform_action, set_calibrations_action,export_seismogram_data_action))
        
        
    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action
        
    def load_experiment(self):
        experiment = SeismicExperiment()
        log_file_name = QFileDialog.getOpenFileName(self,'Open Experiment - Select Log File')
        path = str(log_file_name)
        experiment.seismic_path = os.path.split(path)[0] + '/'
        experiment.seismic_prefix = os.path.split(path)[1].strip('log.txt')
        experiment.ReadLogFile()
        experiment.ReadSeismograms()
        
        # Put seismograms in loaded_list 
        rows_seismograms = []
        for key in experiment.seismograms:
            fname = experiment.seismograms[key].filename
            fname = os.path.split(fname)[1]
            row = experiment.seismograms[key].rec_number
            rows_seismograms.append([row,fname])

        rows_seismograms.sort()
        for item in rows_seismograms:
            QListWidgetItem(item[1], self.loaded_list)
            
        self.experiment = experiment
        
        self.update_experiment_meta()
            
        return self.experiment
        
    def update_experiment_meta(self):
        string_date = self.experiment.date
        self.exp_date_label.setText('Date: %s' %string_date)
        self.exp_recrate_label.setText('Recording Rate [MHz]: %.2f' %(self.experiment.record_rate/1e6))
        self.exp_pulserate_label.setText('Pulse Rate [pps]: %d' %self.experiment.trigger_rate)
        self.exp_posttrigrecs_label.setText('Post-Trigger Samples: %d' %self.experiment.post_trig_recs)
        self.exp_pretrigrecs_label.setText('Pre-Trigger Samples: %d' %self.experiment.pre_trig_recs)
        self.exp_totalrecs_label.setText('Total Samples: %d' %self.experiment.total_recs) 

class Seismogram:
    """
    Class to hold a seismogram and all the data pertaining 
    to that seismogram.
    """
    
    def __init__(self):
        self.rec_number    = 0     # Seismic record number [int]
        self.sync_voltage  = 0.0   # Sync voltage recorded at seismic computer [float, volts]
        self.amplitude     = None  # Stores the amplitude array [float array, volts]
        self.timestamp     = None  # Timestamp from the log file as a datetime object 
        self.filename      = ''    # Filename of the record
        self.maxamp        = 0.0   # Maximum amplitude (biggest value - smallest value)
        self.cross_correlation_picks = None  # Row of first arrival [int, row number]
        self.manual_picks  = None  # Array of rows of picks [array, row number[
        self.format = None
        
    def set_file(self,fname):
        """
        Used to read in seismogram.
        """
        self.filename    = fname # Set filename
        self.rec_number = fname.split('rec')[1]
        self.rec_number = int(self.rec_number.split('.')[0])
        
        GoodRead = False
        
        # Fist we will try to read the binary format.  If that errors out then
        # we try to read as ASCII.
        
        try:
            self.amplitude,length = self.ReadBinSeismo(self.filename) # Read waveform
            GoodRead = True
        except:
            GoodRead = False
            
        if GoodRead == False:
            try:
                self.amplitude,length = self.ReadAsciiSeismo(self.filename) # Read waveform
                GoodRead = True
            except:
                GoodRead = False
      
        self.maxamp = self.get_maxdeltaamplitude() # Calculate and store max amplitude
        
    def get_maxdeltaamplitude(self):
        """
        Return the maximum amplitude difference over the entire waveform.
        """ 
        return max(self.amplitude) - min(self.amplitude)
    
    def RemoveMean(self):
        """
        Remove the mean of the seismogram from the amplitude array.
        """
        self.amplitude = self.amplitude - np.mean(self.amplitude[0:250])

    def ReadBinSeismo(self,fname):
        """
        Read the binary seismogram file written from the vi.  Kudos to 
        Ryan May from Univeristy of Oklahoma for decoding the format.
        """
        f = open(fname,'rb')
        seismo = np.fromfile(f,dtype='<i4',count=-1,sep='')
        
        # For 2D array, size information is encoded in 2 4-byte integers
        size = seismo[0:2]
        
        # After grabbing the size, take the remaining data and re-interpret as
        # doubles
        seismo = seismo[2:]
        seismo.dtype = '<d'

        return seismo.reshape(size[1]),size[1]
        
    def ReadAsciiSeismo(self,fname):
        """
        Read the ASCII seismogram file written from the vi.
        """
        f = open(fname,'r')
        seismo = np.loadtxt(f)
        return seismo, np.size(seismo)

            
class SeismicExperiment:
    """
    Class to store the experiment.  Holds names and data from both the seismic 
    setup and the biax.  Performs data alignment operations and amplitude mapping.
    """
    def __init__(self):
        self.seismic_path  = '' # Path to directory containing the seismic records
        self.number_seismograms = 0 # Number of seismograms in the folder
        self.seismic_prefix = '' # Prefix to the seismic files
        self.rec_data = None # Stores data from the seismic log file about each record
        self.footer = None # Stores data from the footer of the seismic file
        self.seismograms = {} # Dictionary of seismogram file names avaliable
        self.time = None #Stores time in microseconds for the seismogram x axis
        self.trigger_time = 0
    
    def ReadLogFile(self):
        """
        Read the logfile for the seismic records and parse data from it.
        """
        # The logfile has the siesmic prefix followed by log.txt
        fname = self.seismic_path + self.seismic_prefix + 'log.txt'
        f = open(fname,'r')
        
        # We store this data in an array 
        self.rec_data = np.array([],dtype='S')
        self.rec_data = self.rec_data.reshape([0,3])
        
        for line in f:
        
            if 'c:\\' in line:
                # This is a file record
                line    = line.split()
                fname   = line[0].split('\\')[-1]
                time    = line[1]
                ampm    = line[2]
                self.rec_data = np.append(self.rec_data,[[fname,time,ampm]],axis=0 )
                
            elif '/' in line:
                # This is probably the footer line
                print "\n\nFooter reads:"
                line = line.split()
                self.date = line[0]
                self.trigger_rate = int(line[1])
                self.record_rate  = int(line[2])
                self.post_trig_recs = int(line[3])
                self.total_recs = int(line[4])
                self.pre_trig_recs = self.total_recs - self.post_trig_recs
                #self.footer = [date, trigger_rate,record_rate,post_trig_recs,total_recs]
                    
                # Show footer data
                print "-----------------------------------------"
                print "|%25s|%13s|" %('Date', self.date)
                print "|%25s|%13s|" %('Trigger Rate', self.trigger_rate)
                print "|%25s|%13s|" %('Recording Rate', self.record_rate)
                print "|%25s|%13s|" %('Post Trigger Data Points', self.post_trig_recs)
                print "|%25s|%13s|" %('Total Data Points', self.total_recs)
                print "-----------------------------------------"
                
            else:
                print "ERROR: Unknown line encountered in log file."

        self.time = np.arange(self.total_recs)/float(self.trigger_rate)
        self.trigger_time = (self.pre_trig_recs) * (1/float(self.trigger_rate))

    def ReadSeismograms(self):
        """
        By using file names from the log file we read the seismogram files
        into a dictionary of objects.
        """
        print "Reading seismograms to dictionary"
        FirstSeismogram = True
        for rec in self.rec_data:
            f = rec[0]
            self.seismograms[f] = Seismogram()
            self.seismograms[f].set_file(self.seismic_path+f)
            timestr = self.date + ' ' + rec[1] + ' ' + rec[2]    
            self.seismograms[f].timestamp = datetime.strptime(timestr,'%m/%d/%Y %I:%M:%S %p')

class SeismicProcessing:
            
    def CrossCorrelate(self,pattern):
        for seismogram in self.seismograms:
            seismogram = self.seismograms[seismogram]
            xcorr = np.correlate(seismogram.amplitude,pattern,"valid")
            max_corr_row = xcorr.argmax() 
            seismogram.first_arrival = max_corr_row
            seismogram.first_arrival_time = (max_corr_row - self.pre_trig_recs) * (1./self.record_rate) * 1e6   
        
def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()