import numpy as np
import matplotlib.pyplot as plt
import pyqtgraph
import pyqtgraph.opengl as gl
from PyQt5 import QtCore, QtGui, uic
import os,sys,time

home_dir = os.path.dirname(os.path.realpath(__file__))

sys.path.insert(0,home_dir)
from plot_galaxy import elem_to_state


star_file = os.path.abspath(os.path.join(home_dir,'../data/stars.txt'))
ui_file = os.path.join(home_dir,'template.ui')

star_data = np.loadtxt(star_file, delimiter=',', skiprows=1)
ids = star_data[:, 0].astype(int)
star_elem = star_data[:, 1:]

class GTOC(QtGui.QMainWindow):
    def __init__(self):
        super(GTOC,self).__init__()
        uic.loadUi(ui_file,self)
        self.time = 0
        self.color_set = plt.get_cmap('jet')
        self.init_plots()
        self.update_time_plot()

    def init_plots(self):
        self.time_plot = gl.GLScatterPlotItem(pos=np.array([0,0,0]).reshape(-1,3))
        self.time_graphics_view.addItem(self.time_plot)
        self.speeds = None

    def update_time_from_click_box(self,time):
        self.time_slider.setValue(int(time))
        self.time = self.time_click_box.value()
        self.update_time_plot()

    def update_time_from_slider(self,time):
        self.time_click_box.setValue(time)
        self.time = self.time_click_box.value()
        self.update_time_plot()

    def update_time_plot(self):
        start = time.time()
        states = elem_to_state(star_elem,self.time)
        if self.speeds is None:
            self.speeds = np.sqrt(np.sum(states[:,3:]**2,axis=1))
            self.speeds -= np.mean(self.speeds)
            self.speeds /= np.std(self.speeds)*3
            self.speeds += 1
            self.speeds /=2
            self.speeds = np.clip(self.speeds,0,1)
        color = self.color_set(self.speeds)
        self.time_plot.setData(pos=states[:,:3],color=color,size=1)
        
                
if __name__=='__main__':
    app = QtGui.QApplication(sys.argv)
    window = GTOC()
    window.show()
    sys.exit(app.exec_())
