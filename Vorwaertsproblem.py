import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import vorwaerts as vw

class Vorwaertsproblem():
    
    def __init__(self, l_e=10, l_r=3, h_m=30, h_b=90, phi_0=135, delta_phi=45, psi=90):
        self.l_e = l_e
        self.l_r = l_r
        self.h_m = h_m
        self.h_b = h_b        
        self.phi = phi_0
        self.delta_phi = delta_phi 
        self.psi = psi 
        
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
        pointphix, pointphiy = self.plot_line(0, 0, np.radians(0), self.h_m, name="hm", color='k-')
        self.ax.plot(self.h_m, 0, 'kx')
        self.ax.text(self.h_m, -0.3, 'h_m')
        # Untere Bodenlinie
        pointpsix, pointpsiy = self.plot_line(0, 0, np.radians(0), self.h_b, name="hb", color='k-')
        self.ax.plot(self.h_b, 0, 'kx')
        self.ax.text(self.h_b, -0.3, 'h_b')
        # Wand
        topx, topy = self.plot_line(pointpsix, pointpsiy, np.pi - psi, 1.3*b, plot=False)
        self.plot_line(topx, topy, 2*np.pi-psi, 2.6*b, color='k-')
        # Drehteller
        circlele = plt.Circle((pointphix, pointphiy), self.l_e, color='k', linestyle='dashed', fill=False)
        self.ax.add_artist(circlele)    
        # Winkelhalbierende
        pointbx, pointby = self.plot_line(0,0, alpha, h, name="h", color='g-')
        # Winkelhalbierende bis Kreismittelpunkt
        pointrx, pointry = self.plot_line(0, 0, alpha, h_r, name="hr", color = 'b-')
        # Abstand Mitte Drehteller, Mitte Zylinder
        self.ax.plot((self.h_m, pointrx), (0, pointry), 'y-', label='l_e')
        self.ax.plot(pointrx, pointry, 'kx')
        self.ax.text(pointrx, pointry-0.3, 'M')
        # Zylinder
        circler = plt.Circle((pointrx,pointry), self.l_r, color='k', fill=False)
        self.ax.add_artist(circler)
        # b an Wand
        btopx, btopy = self.plot_line(pointbx, pointby, np.pi-psi, b_top, plot=False)
        self.plot_line(btopx, btopy, 2*np.pi-psi, b, name="b", color='r-')
        bbottomx, bbottomy = self.plot_line(pointbx, pointby, 2*np.pi-psi, b_bottom, plot=False)
        # Unterer Scheitel
        self.ax.plot((0, bbottomx), (0, bbottomy), 'c-.')
        self.ax.plot(bbottomx, bbottomy, 'kx')
        # Oberer Scheitel
        self.ax.plot((0, btopx), (0, btopy), 'c-.')
        self.ax.plot(btopx, btopy, 'kx')
        self.ax.axis('equal')
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
     vorwaertsproblem = Vorwaertsproblem()
     vorwaertsproblem.plot()
