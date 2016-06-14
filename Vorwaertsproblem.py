import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import vorwaerts as vw

class Vorwaertsproblem():
    
    def __init__(self, l_e=10, l_r=3, h_m=30, h_b=90, phi_0=135, delta_phi=45, psi=90, show_legend=True):
        self.l_e = l_e
        self.l_r = l_r
        self.h_m = h_m
        self.h_b = h_b        
        self.phi = phi_0
        self.delta_phi = delta_phi 
        self.psi = psi 
        self.show_legend = show_legend
        
    def initUI(self):
        self.fig, self.ax = plt.subplots()
        self.fig.canvas.set_window_title('Numerik DGL 2 - Schattenprojektion')
        self.fig.suptitle('Vorwaertsproblem')
        plt.subplots_adjust(bottom=0.3)
        plt.axis([0, 200, -100, 100])
        plt.axis('equal')
        axColor = 'lightgoldenrodyellow'
        
        # Slider - Phi
        axPhi = plt.axes([0.18, 0.2, 0.65, 0.03], axisbg=axColor)                
        self.sliderPhi = Slider(axPhi, 'Phi', 0, 360, valinit=self.phi, valfmt='%1d')
        self.sliderPhi.on_changed(self.onUpdate)
        # Slider - Delta Phi
        axDeltaPhi = plt.axes([0.18, 0.15, 0.65, 0.03], axisbg=axColor)        
        self.sliderDeltaPhi = Slider(axDeltaPhi, 'Delta Phi', 5, 360, valinit=self.delta_phi, valfmt='%1d')        
        self.sliderDeltaPhi.on_changed(self.onUpdate)
        # Slider - Psi
        axPsi = plt.axes([0.18, 0.1, 0.65, 0.03], axisbg=axColor)        
        self.sliderPsi = Slider(axPsi, 'Psi', 0, 180, valinit=self.psi, valfmt='%1d')        
        self.sliderPsi.on_changed(self.onUpdate)
        # Button - Previous
        axPrev = plt.axes([0.18, 0.03, 0.1, 0.05])
        self.buttonPrev = Button(axPrev, 'Previous')
        self.buttonPrev.on_clicked(self.onPrevious)
        # Button - Next
        axNext = plt.axes([0.29, 0.03, 0.1, 0.05])
        self.buttonNext = Button(axNext, 'Next')
        self.buttonNext.on_clicked(self.onNext)        
        # Button - Reset        
        axReset = plt.axes([0.73, 0.03, 0.1, 0.05])
        self.buttonReset = Button(axReset, 'Reset')
        self.buttonReset.on_clicked(self.onReset)
        
        self.onDraw()
    
    def onPrevious(self, event):        
        self.phi -= self.delta_phi    
        self.phi = self.phi % 360
        self.onDraw()
    
    def onNext(self, event):    
        self.phi += self.delta_phi
        self.phi = self.phi % 360
        self.onDraw()
    
    def onReset(self, event):
        self.sliderPhi.reset()
        self.sliderDeltaPhi.reset()    
        self.sliderPsi.reset()    
        self.onDraw()
        
    def onUpdate(self, val):
        self.phi = int(self.sliderPhi.val)
        self.delta_phi = int(self.sliderDeltaPhi.val)
        self.psi = int(self.sliderPsi.val)
        self.onDraw()        
        self.fig.canvas.draw_idle()
        
    def onDraw(self):
        self.ax.cla()  # clear the axes 
        psi = np.radians(self.psi)
        phi = np.radians(self.phi)
        # Werte berechnen
        h_r = vw.get_h_r(self.l_e, self.h_m, phi)
        alpha = vw.get_alpha(self.l_e, phi, h_r)
        beta = vw.get_beta(alpha, psi)
        h = vw.get_h(psi, self.h_b, beta)    
        gamma = vw.get_gamma(self.l_r, h_r)
        b_bottom = vw.get_b_bottom(h, beta, gamma)
        b_top = vw.get_b_top(h, beta, gamma)
        b = vw.get_b(h, beta, gamma)
        # Abstand bis Drehteller vom Ursprung
        point_phi_x, point_phi_y = self.plot_line(0, 0, np.radians(0), self.h_m, name="hm", color='k-')
        self.ax.plot(self.h_m, 0, 'kx')
        self.ax.text(self.h_m, -0.3, 'h_m')
        # Untere Bodenlinie
        point_psi_x, point_psi_y = self.plot_line(0, 0, np.radians(0), self.h_b, name="hb", color='k-')
        self.ax.plot(self.h_b, 0, 'kx')
        self.ax.text(self.h_b, -0.3, 'h_b')
        # Drehteller
        circlele = plt.Circle((point_phi_x, point_phi_y), self.l_e, color='k', linestyle='dashed', fill=False)
        self.ax.add_artist(circlele)    
        # Winkelhalbierende
        point_b_x, point_b_y = self.plot_line(0,0, alpha, h, name="h", color='g-')
        # Winkelhalbierende bis Kreismittelpunkt
        point_r_x, point_r_y = self.plot_line(0, 0, alpha, h_r, name="hr", color = 'b-')
        # Abstand Mitte Drehteller, Mitte Zylinder
        self.ax.plot((self.h_m, point_r_x), (0, point_r_y), 'y-', label='l_e')
        self.ax.plot(point_r_x, point_r_y, 'kx')
        self.ax.text(point_r_x, point_r_y-0.3, 'M')
        # Zylinder
        circler = plt.Circle((point_r_x,point_r_y), self.l_r, color='k', fill=False)
        self.ax.add_artist(circler)
        # x,y Koordinate der Punkt b_top und b_bottom bestimmen
        b_top_x, b_top_y = self.plot_line(point_b_x, point_b_y, np.pi-psi, b_top, plot=False)
        b_bottom_x, b_bottom_y = self.plot_line(point_b_x, point_b_y, 2*np.pi-psi, b_bottom, plot=False)
        # Wand
        if b_bottom_y > point_psi_y and b_top_y > point_psi_y:  # Schatten oberhalb der Mittellinie
            kathete_top_a = np.abs(b_bottom_x-point_psi_x)
            kathete_top_b = np.abs(b_bottom_y-point_psi_y)
            len_top = 2*np.sqrt(kathete_top_a**2+kathete_top_b**2)+b
            len_bottom = 1
        elif b_bottom_y < point_psi_y and b_top_y < point_psi_y:  # Schatten unterhalb der Mittellinie            
            kathete_bottom_a = np.abs(b_top_x-point_psi_x)
            kathete_bottom_b = np.abs(b_top_y-point_psi_y)
            len_bottom = 2*np.sqrt(kathete_bottom_a**2+kathete_bottom_b**2)+b
            len_top = 1
        else:  # Schatten ober- und unterhalb der Mittellinie
            kathete_top_a = np.abs(b_top_x-point_psi_x)
            kathete_top_b = np.abs(b_top_y-point_psi_y)
            len_top = np.sqrt(kathete_top_a**2+kathete_top_b**2) + 1
            kathete_bottom_a = np.abs(b_bottom_x-point_psi_x)
            kathete_bottom_b = np.abs(b_bottom_y-point_psi_y)
            len_bottom = np.sqrt(kathete_bottom_a**2+kathete_bottom_b**2) + 1 
        self.plot_line(point_psi_x, point_psi_y, np.pi-psi, len_top, color="k-")  # obere Haelfte
        self.plot_line(point_psi_x, point_psi_y, 2*np.pi-psi, len_bottom, color="k-") # untere Haelfte
        # b an Wand
        self.plot_line(b_top_x, b_top_y, 2*np.pi-psi, b, name="b", color='r-')
        # Unterer Scheitel
        self.ax.plot((0, b_bottom_x), (0, b_bottom_y), 'c-.')
        self.ax.plot(b_bottom_x, b_bottom_y, 'kx')
        # Oberer Scheitel
        self.ax.plot((0, b_top_x), (0, b_top_y), 'c-.')
        self.ax.plot(b_top_x, b_top_y, 'kx')
        self.ax.axis('equal')
        if self.show_legend:
            self.ax.legend()
        
    def plot_line(self, pointx, pointy, angle, distance, plot=True, name=None, color='-'):
        x = np.cos(angle)*distance+pointx
        y = np.sin(angle)*distance+pointy
        if plot:
            self.ax.plot((pointx,x),(pointy,y),color, label=name)
        return x,y
        
    def plot(self):
        self.initUI()
        plt.show()


if __name__ == '__main__':
     vorwaertsproblem = Vorwaertsproblem(show_legend=False)
     vorwaertsproblem.plot()
