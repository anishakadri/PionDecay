#help: z position is being weird, if i shift it to all positive.
#separate simulation and data collection.
#class variables-->local variable 'pion_list' referenced before assignment.

import Particles as pa
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from mpl_toolkits.mplot3d import Axes3D

class Chamber:
    
    """
    A class defining a cylindrical decay chamber, defined by the following 
    arguments:
        Emin; the minimum energy of a pion source,
        Emax; the maximum energy of a pion source,
        d; the diameter of the cylinder,
        l; the length of the cyclinder.
    Other data attributes include a list of electrons, muons and pions created in a 
    simulation (electron_list, pion_list, muon_list).
    
    """
    def __init__(self, Emin= 500.0, Emax= 10000.0, d= 5.0, l= 100.0):
        self.e_min = Emin
        self.e_max = Emax
        self.d = d
        self.l = l
        self.electrons = np.array([])
        self.mu_electrons = np.array([])
        self.pi_electrons = np.array([])
        self.pions = np.array([])
        self.muons = np.array([])
    
    def __repr__(self):
        return "%s( Emin=%r, Emax=%r, diameter=%r, length=%r)" % ("Chamber", self.e_min, self.e_max, self.d, self.l)
        
    def electron_num(self):
        return len(self.electrons)
    
    def MUe_num(self):
        return len(self.mu_electrons)
    
    def PIe_num(self):
        return len(self.pi_electrons)
        
    def pion_num(self):
        print len(self.pions)
        
    def muon_num(self):
        return len(self.muons)
    
    def clear_particles(self):
        self.electrons = np.array([])
        self.mu_electrons = np.array([])
        self.pi_electrons = np.array([])
        self.pions = np.array([])
        self.muons = np.array([])
        
    def create_pions(self, E=500.0, n=1000):
        """A function that creates pions and adds it to 
        a list of muons.for a chamber 
        running at constant energy.
        It takes 2 arguments:
        E ; the total energy (in MeV) of the pions produced by the source.
        n ; the number of pions produced by the source.
        """
        for y in range(0, n):
            pi = pa.Pion(E)
            self.pions = np.append(pi, self.pions)

    def create_muons(self):
        """A function that creates muons from the list of pre-existing pions
        and adds it to a list of muons.
        (for a chamber running at constant energy.)"""
        if len(self.pions) == 0:
            raise Exception("Create pions first.")
        else:
            for i in self.pions:
                mu = i.decay_muon()
                if mu == None:
                    mu = 0.0
                else:
                    self.muons = np.append(mu, self.muons)
    
    def create_PIe(self):
        """A function that creates electrons via pion decay and adds it to 
        a list of pi- electrons and electrons.(for a chamber 
        running at constant energy.)
        It takes 2 arguments:
        E ; the total energy (in MeV) of the pions produced by the source.
        n ; the number of pions produced by the source."""
        if len(self.pions) == 0:
            raise Exception("Create pions first.")
        else:
            for i in self.pions:
                mu = i.decay_electron()
                if mu == None:
                    mu = 0.0
                else:
                    self.electrons = np.append(mu, self.electrons)
                    self.pi_electrons = np.append(mu, self.pi_electrons)
            
    def create_MUe(self):
        """A function that creates electrons via muon decay and adds it to 
        a list of mu- electrons and electrons.for a chamber 
        running at constant energy.
        It takes 2 arguments:
        E ; the total energy (in MeV) of the pions produced by the source.
        n ; the number of pions produced by the source."""
        if len(self.muons) == 0:
            raise Exception("Create muons first.")
        else:
            for i in self.muons:
                el = i.decay_electron()
                if el == None:
                    el = 0.0
                else:
                    self.electrons = np.append(el, self.electrons)
                    self.mu_electrons = np.append(el, self.mu_electrons)
                    
    def momentum(self, particle):
        x = np.array([])
        y = np.array([])
        z = np.array([])
        for z in electron_list:
            x = np.append(z.momentum()[0],x)
            y = np.append(z.momentum()[1],y)
            z = np.append(z.momentum()[2],z)
        
        return x,y,z
        
    def pielectron_scatter(self, E=500.0, n=1000):
        """A function that returns a scatter plot of the final positions of 
        electrons decaying from pions, for a chamber running at constant energy.
        It takes 2 arguments:
        E ; the total energy (in MeV) of the pions produced by the source.
        n ; the number of pions produced by the source."""
    
        self.clear_particles()
        self.create_pions(E,n)
        self.create_PIe()
        
        electron_finalx = np.array([])
        electron_finaly = np.array([])
        electron_finalz = np.array([])
        
        for z in self.electrons:
            electron_finalx = np.append(z.final_position()[0],electron_finalx)
            electron_finaly = np.append(z.final_position()[1],electron_finaly)
            electron_finalz = np.append(z.final_position()[2],electron_finalz)
                
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    
        ax.scatter(electron_finalx, electron_finaly, electron_finalz, c= 'k')

        grid_x = np.linspace(-1, 1, 100)
        grid_z = np.linspace(0, 4, 100)
        x, z = np.meshgrid(grid_x, grid_z)
        y = np.sqrt(1-x**2)
        x = x*2.5
        y = y*2.5
        z = z*25    
        rstride = 20
        cstride = 10
        ax.plot_surface(x, y, z, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(x, -y, z, alpha=0.2, rstride=rstride, cstride=cstride)
    
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        plt.show() 
    
    def muelectron_scatter(self, E=500.0, n=1000):
        """A function that returns a scatter plot of the final positions of 
        electrons decaying from muons, for a chamber running at constant energy.
        It takes 2 arguments:
        E ; the total energy (in MeV) of the pions produced by the source.
        n ; the number of pions produced by the source."""
    
        self.clear_particles()
        self.create_pions(E,n)
        self.create_muons()
        self.create_MUe()
        
        electron_finalx = np.array([])
        electron_finaly = np.array([])
        electron_finalz = np.array([])
        
        for z in self.electrons:
            electron_finalx = np.append(z.final_position()[0],electron_finalx)
            electron_finaly = np.append(z.final_position()[1],electron_finaly)
            electron_finalz = np.append(z.final_position()[2],electron_finalz)
                
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    
        ax.scatter(electron_finalx, electron_finaly, electron_finalz, c= 'w')

        grid_x = np.linspace(-1, 1, 100)
        grid_z = np.linspace(0, 4, 100)
        x, z = np.meshgrid(grid_x, grid_z)
        y = np.sqrt(1-x**2)
        x = x*2.5
        y = y*2.5
        z = z*25    
        rstride = 20
        cstride = 10
        ax.plot_surface(x, y, z, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(x, -y, z, alpha=0.2, rstride=rstride, cstride=cstride)
    
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        plt.show() 
        
    def muon_scatter(self, E=500.0, n=1):
        """A function that returns a scatter plot of the final positions of muons
        for a chamber running at constant energy.
        It takes 2 arguments:
        E ; the total energy (in MeV) of the pions produced by the source.
        n ; the number of pions produced by the source."""
    
        muon_finalx = np.array([])
        muon_finaly = np.array([])
        muon_finalz = np.array([])
        self.clear_particles()
        self.create_pions(E,n)
        self.create_muons()
        
        for n in self.muons: 
            muon_finalx = np.append(n.final_position()[0],muon_finalx)
            muon_finaly = np.append(n.final_position()[1],muon_finaly)
            muon_finalz = np.append(n.final_position()[2],muon_finalz)
                
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    
        ax.scatter(muon_finalx, muon_finaly, muon_finalz, c='r')

        grid_x = np.linspace(-1, 1, 100)
        grid_z = np.linspace(0, 4, 100)
        x, z = np.meshgrid(grid_x, grid_z)
        y = np.sqrt(1-x**2)
        x = x*2.5
        y = y*2.5
        z = z*25    

        rstride = 20
        cstride = 10
        ax.plot_surface(x, y, z, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(x, -y, z, alpha=0.2, rstride=rstride, cstride=cstride)
    
        ax.set_xlabel("x distance/(m)")
        ax.set_ylabel("y distance/ (m)")
        ax.set_zlabel("z distance/ (m)")
        plt.show()     
         
    def pion_scatter(self, E=500.0, n=1):
        """A function that returns a scatter plot of the final positions of pions
        for a chamber running at constant energy.
        It takes 2 arguments:
        E ; the total energy (in MeV) of the pions produced by the source.
        n ; the number of pions produced by the source."""
        
        self.clear_particles()
        self.create_pions(E,n)
        pion_list = np.array([])
        pion_finalx = np.array([])
        pion_finaly = np.array([])
        pion_finalz = np.array([])
        
        for i in self.pions: 
            pion_finalx = np.append(i.final_position()[0],pion_finalx)
            pion_finaly = np.append(i.final_position()[1],pion_finaly)
            pion_finalz = np.append(i.final_position()[2],pion_finalz)
                
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
    
        ax.scatter(pion_finalx, pion_finaly, pion_finalz, c='g')

        grid_x=np.linspace(-1, 1, 100)
        grid_z=np.linspace(0, 4, 100)
        x, z=np.meshgrid(grid_x, grid_z)
        y = np.sqrt(1-x**2)
        x = x*2.5
        y = y*2.5
        z = z*25    

        rstride = 20
        cstride = 10
        ax.plot_surface(x, y, z, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(x, -y, z, alpha=0.2, rstride=rstride, cstride=cstride)
    
        ax.set_xlabel("x distance/(m)")
        ax.set_ylabel("y distance/ (m)")
        ax.set_zlabel("z distance/ (m)")
        plt.show()     
    
    def escape_ratio(self, E= 550):
        """ This function finds  the proportion of muons 
        leaving the chamber at a certain energy E."""
        self.clear_particles()
        self.create_pions(E, 5000)
        self.create_muons()
        self.create_MUe()                   
        decay_ratio = np.float64(len(self.electrons))/np.float64(len(self.muons))
        esc_ratio = 1.0 - decay_ratio
        return esc_ratio
        
    def particle_tracks(self, electron_initial, electron_final, muon_initial, muon_final, E = 500.0):
        mpl.rcParams['legend.fontsize'] = 10
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        z1 = np.linspace(0, 50, 100)
        x1 = np.linspace(0, 0, 100)
        y1 = np.linspace(0, 0, 100)
        ax.plot(x1, y1, z1, 'b', label='pion')


        z2 = np.linspace(50, 100, 100)
        x2 = np.linspace(0, 2.5, 100)
        y2 = np.linspace(0, 2.5, 100)
        ax.plot(x2, y2, z2, 'r', label='muon')

        ax.legend()
        grid_x=np.linspace(-1, 1, 100)
        grid_z=np.linspace(0, 4, 100)
        x, z=np.meshgrid(grid_x, grid_z)
        y = np.sqrt(1-x**2)
        x = x*2.5
        y = y*2.5
        z = z*25    
        rstride = 20
        cstride = 10
        ax.plot_surface(x, y, z, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(x, -y, z, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.set_xlabel("x distance/(m)")
        ax.set_ylabel("y distance/ (m)")
        ax.set_zlabel("z distance/ (m)")
        plt.show() 
      
        