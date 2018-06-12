#HELP:
# 1) case of -ve velocity for muons and electrons
# 2) compare decay time and exit time for muons
# try a check decay pos, raise exception, use try to call loop
#3) what to do when lines exceed the width??-- google

import properties as pp
import fourmomentum as rel
import numpy as np 
import scipy.optimize as spo

class Particle:
    """
    A class definining an instance of a particle;
    defined by its energy-momentum four vector, and the initial and 
    final vertices of its path. All quantities are required to be 
    in natural units.
    """
    def __init__(self, E=0.00, p= [0.00,0.00,0.00], i=[0.00,0.00,0.00],f=[0.00,0.00,0.00]):
        self.EP= rel.Fourvector(E, p)
        self.initial = np.array(i)
        self.final = np.array(f)
        if len(p) > 3:
            raise Exception("Error: Momentum Vector parameter size")
        if len(i) > 3:
            raise Exception("Error: Initial Vector parameter size")
        if len(f) > 3:
            raise Exception("Error: Final Vector parameter size")
            
    def __repr__(self):
        return "%s(Energy-Momentum %r, intial=%r, final=%r)" % ("Particle", self.EP, self.initial, self.final)
    
    def __str__(self):
        return "[%g, %r]" % (self.__E, self.__p)    
            
    def momentum(self):
        """Returns the momentum of a particle in the lab frame."""
        return self.EP.momentum()
        
    def mass(self):     
        """Returns the rest mass of a particle, 
        which is equivalent to the rest energy in natural units."""
        return np.sqrt((self.energy()**2)-self.EP.magp_sq())
    
    def energy(self):
        """Returns the total energy of the particle."""
        return self.EP.energy()
        
    def kin_energy(self):
        """Returns the kinetic energy of a particle."""
        return self.energy() - self.rest_energy()
        
    def beta(self):
        """Returns the beta value of a particle.
        For a particle travelling at a speed v, Beta = v/c ."""
        if self.energy() == 0.0:
            return "Total Energy = 0, particle does not exist."
        elif self.EP.magp()/self.energy() >= 1:
            return "Error: Particle travelling faster than light."
        else:
            return self.EP.magp()/self.energy()
            
    def gamma(self):
        """Returns the gamma factor for a particle. 
        In particular, the gamma factor of the particles rest frame 
        with respect to the lab frame."""
        return (self.energy()/self.mass())
        
class Pion(Particle):
    """
    Creates a Pion, with a given total energy (E), a vector momentum (p),
    and initial (i) and final (f) vertices of its path in cartesian coordinates. 
    Although a momentum can be entered upon initialisation, the momentum for the 
    pion is automatically generated to agree with the Energy-Momentum relation.
    """
    def __init__(self, E=0.00, p= [0.00,0.00,0.00], i=[0.00,0.00,0.00],f=[0.00,0.00,0.00]):
        
        Particle.__init__(self, E, p, i, f)
        p_mag = np.sqrt(E**2-(pp.pion_mass)**2)
        p = p_mag*np.array([0.0, 0.0, 1.0])
        self.EP = rel.Fourvector(E, p)
        self.initial = i
        self.final = f
        self.__decay_time = 0
        self.__exit_time = 0
        if len(p) > 3:
            raise Exception("Momentum Vector parameter size")
        if len(i) > 3:
            raise Exception("Initial Vector parameter size")
        if len(f) > 3:
            raise Exception("Final Vector parameter size")
        if E < pp.pion_mass:
            raise Exception("Total energy is less than Pion rest energy.") 
          
    def random_decay(self):
        """Decays a pion into either an electron or a muon, using the predefined
        branching ratio."""
        x= np.random.random_integers(0, pp.electron_ratio, 1)                    #binomial not as accurate at sampling.
        for i in x:
            if i == 1:
                return self.decay_electron()
            else:
                return self.decay_muon()
            
    def final_position(self):
        """Calculates the final position of a particle, whether it decays within
        the chamber or exits one of the chamber sides."""
        exit_pos = self.exit_position()                    #Lab note: have to define as set variables
        decay_pos = self.decay_position()
        if self.__exit_time > self.__decay_time:
            self.final = decay_pos
            return self.final
        else:
            self.final = exit_pos
            return self.final
            
    def exit_position(self):
        """Calculates the final position of a particle, 
        if it exits one of the chamber sides."""
        velocity = (self.momentum()*(1/pp.pion_mass))*pp.c
        exit_time = (100.00-self.initial[2])/(velocity[2])
        self.__exit_time = exit_time
        return self.initial + (exit_time* velocity)
        
    def decay_position(self):
        """Calculates the final position of a particle, if it decays."""
        velocity = (self.momentum()*(1/pp.pion_mass))*pp.c
        pion_time =  np.random.exponential(pp.pion_tau)    # Lab note: exponential(poisson) distribution of mean lifetime
        lab_time = self.gamma()*pion_time
        self.__decay_time = lab_time
        return self.initial + (lab_time*(velocity))
        
    def decay_muon(self):
        """Decays a pion into a muon."""
        mu_start = self.final_position()
        if self.__decay_time > self.__exit_time:
            return None
        else:
            muon_energy = ((pp.pion_mass**2) + (pp.muon_mass**2))/(2*pp.pion_mass)
            p_mag = np.sqrt(muon_energy**2 - pp.muon_mass**2)
            phi = np.random.uniform(0.0, 2*np.pi)    # Lab note: particles uniformly distributed in phi
            theta = np.arccos(np.random.uniform(-1.0, 1.0))    # Lab note: particles uniformly distributed in cos theta(since azimuthal angle)
            p_dir = np.array([np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta), np.cos(phi)])
            EP_pionframe = rel.Fourvector(muon_energy,p_mag*p_dir)
            EP_labframe = EP_pionframe.boost(-self.beta())
            return Muon(EP_labframe.energy(), EP_labframe.momentum(), i=mu_start, f=[])
               
    def decay_electron(self):
        """Decays a pion into an electron."""
        el_start = self.final_position()
        if self.__decay_time > self.__exit_time:
            return None
        else:
            electron_energy = ((pp.pion_mass**2) + (pp.electron_mass**2))/(2*pp.pion_mass)  
            p_mag = np.sqrt(electron_energy**2-pp.electron_mass**2)
            phi = np.random.uniform(0.0, 2*np.pi)    
            theta = np.arccos(np.random.uniform(-1.0,1.0))   
            p_dir = np.array([np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta), np.cos(phi)])
            EP_pionframe = rel.Fourvector(electron_energy,p_mag*p_dir)
            EP_labframe = EP_pionframe.boost(-self.beta())
            return Electron(EP_labframe.energy(), p=EP_labframe.momentum(), i=el_start, f=[])
            
    def energy_deposit(self, x):
        """Calculates the energy deposited in an NaI scintillator by a pion, 
        for every x cm travelled through it. Takes one argument, x."""
        return 4.8*x

class Muon(Particle):
    """
    Creates a Muon, with a given total energy (E), a vector momentum (p),
    and initial (i) and final (f) vertices of its path in cartesian coordinates.
    """
    def __init__(self, E=0.00, p= [0.00,0.00,0.00], i=[0.00,0.00,0.00],f=[0.00,0.00,0.00]):
        Particle.__init__(self, E, p, i, f)
        self.EP = rel.Fourvector(E, p)
        self.initial = i
        self.final = f
        self.__decay_time = 0
        self.__exit_time = 0
        if len(p) > 3:
            raise Exception("Momentum Vector parameter size")
        if len(i) > 3:
            raise Exception("Initial Vector parameter size")
        if len(f) > 3:
            raise Exception("Final Vector parameter size")
        if E < pp.muon_mass:
            raise Exception("Total energy is less than Muon rest energy.")
                     
    def exit_position(self):
        """Calculates the final position of a particle, 
        if it exits one of the chamber sides."""
        v = (self.momentum()/pp.muon_mass)*pp.c
        Vx, Vy, Vz = v[0], v[1], v[2] 
        x0, y0, z0 = self.initial[0], self.initial[1], self.initial[2]
        if Vz > 0:                                                                # for particle travelling in +z direction
            if (Vx)**2 + (Vy)**2 == 0.0:                                         #if no radial velocity, then exits parallel to z axis
                t = (100.00-z0)/(Vz)
                self.__exit_time = t
                return self.initial + (v*t)
            else:                                                                #if it has radial velocity, it exits when radius of position is > 2.5
                def func(t):
                    return (x0 + Vx*np.fabs(t))**2 + (y0 + Vy*np.fabs(t))**2 -6.25
                t = spo.fsolve(func, 1)
                exit_point = self.initial + (v*np.fabs(t))
                if exit_point[2] > 100.0:                                        #if  leaving circular face, then exit is at z=100
                    t = (100.00-z0)/(Vz)
                    self.__exit_time = t
                    return self.initial + (v*t)
                else:
                    self.__exit_time = t
                    return exit_point
        elif Vz < 0:                                                            #for particle travelling in -z direction
            if (Vx)**2 + (Vy)**2 == 0.0:                                        #if no radial velocity, exits parallel to z axis.
                t = z0/(Vz)
                self.__exit_time = t
                return self.initial + (v*t)
            else: 
                def func(t):
                    return (x0 + Vx*np.fabs(t))**2 + (y0 + Vy*np.fabs(t))**2 -6.25
                t = spo.fsolve(func, 1) 
                exit_point = self.initial + (v*np.fabs(t))
                if exit_point[2] > 100.0:                                                 #if  leaving circular face, then exit is at z=0
                    t = z0/(Vz)
                    self.__exit_time = t
                    return self.initial + (v*t)
                else:
                    self.__exit_time = t
                    return exit_point
        else:                                                                     #if no velocity in z-axis
            if (Vx)**2 + (Vy)**2 == 0.0:                                         #if no velocity at all, it doesn't exit
                return self.initial
            else:
                def func(t):                                                    #if only radial velocity, it exits parallel to x-y plane
                    return (x0 + Vx*np.fabs(t))**2 + (y0 + Vy*np.fabs(t))**2 -6.25
                t = spo.fsolve(func, 3)  
                self.__exit_time = t
                return self.initial + (v*np.fabs(t))
            
    def decay_position(self):
        """Calculates the final position of a particle, if it decays."""
        velocity = (self.momentum()*(1/pp.muon_mass))*pp.c
        muon_time = np.random.exponential(pp.muon_tau)
        lab_time = self.gamma()*muon_time
        self.__decay_time = lab_time
        return self.initial + (lab_time*(velocity))
        
    def final_position(self):
        """Calculates the final position of a particle, whether it decays within
        the chamber or exits one of the chamber sides."""
        exit_pos = self.exit_position()                                            #Lab note: have to define as set variables, since decay position is randomly generated
        decay_pos = self.decay_position()
        if self.__exit_time > self.__decay_time:
            self.final = decay_pos
            return self.final
        else:
            self.final = exit_pos
            return self.final
        
    def michel(self):
        x, y = np.random.random(2)
        if x > y:
            return x*53.
        if y> x:
            return y*53.
    
    def decay_electron(self):
        """Decays a muon into a lectron if it is still in the chamber."""
        el_start = self.final_position()
        if self.__decay_time > self.__exit_time:
            return None
        else:
            electron_energy = self.michel()  
            p_mag = np.sqrt(electron_energy**2 - pp.electron_mass**2)
            phi = np.random.uniform(0.0, 2*np.pi)    
            theta = np.arccos(np.random.uniform(-1.0,1.0))    
            p_dir = np.array([np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta), np.cos(phi)])
            EP_pionframe = rel.Fourvector(electron_energy,p_mag*p_dir)
            EP_labframe = EP_pionframe.boost(-self.beta())
            return Electron(E= EP_labframe.energy(), p=EP_labframe.momentum(), i=el_start, f=[])
            
    def energy_deposit(self, x):
        """Calculates the energy deposited in an NaI scintillator by a muon, 
        for every x cm travelled through it. Takes one argument, x."""
        return 4.8*x
 
class Electron(Particle):
    """
    Creates an electron, with a given total energy (E), a vector momentum (p),
    and initial (i) and final (f) vertices of its path in cartesian coordinates.
    """
    def __init__(self, E=0.00, p=[0.00,0.00,0.00], i=[0.00,0.00,0.00],f=[0.00,0.00,0.00]):
        Particle.__init__(self, E, p, i, f)
        self.EP = rel.Fourvector(E, p)
        self.initial = i
        self.final = f
        self.__exit_time = 0
        if len(p) > 3:
            raise Exception("Error: Momentum Vector parameter size")
        if len(i) > 3:
            raise Exception("Error: Initial Vector parameter size")
        if len(f) > 3:
            raise Exception("Error: Final Vector parameter size")
        if E < pp.electron_mass:
            raise Exception("Total energy is less than electron rest energy.")

    def exit_position(self):
        """Calculates the final position of a particle, 
        if it exits one of the chamber sides."""
        v = (self.momentum()/pp.electron_mass)*pp.c
        Vx, Vy, Vz = v[0], v[1], v[2] 
        x0, y0, z0 = self.initial[0], self.initial[1], self.initial[2]
        if Vz > 0:                                                                # for particle travelling in +z direction
            if (Vx)**2 + (Vy)**2 == 0.0:                                         #if no radial velocity, then exits parallel to z axis
                t = (100.00-z0)/(Vz)
                self.__exit_time = t
                return self.initial + (v*t)
            else:                                                                #if it has radial velocity, it exits when radius of position is > 2.5
                def func(t):
                    return (x0 + Vx*np.fabs(t))**2 + (y0 + Vy*np.fabs(t))**2 -6.25
                t = spo.fsolve(func, 1)
                exit_point = self.initial + (v*np.fabs(t))
                if exit_point[2] > 100.0:                                        #if  leaving circular face, then exit is at z=100
                    t = (100.00-z0)/(Vz)
                    self.__exit_time = t
                    return self.initial + (v*t)
                else:
                    self.__exit_time = t
                    return exit_point
        elif Vz < 0:                                                            #for particle travelling in -z direction
            if (Vx)**2 + (Vy)**2 == 0.0:                                        #if no radial velocity, exits parallel to z axis.
                t = z0/(Vz)
                self.__exit_time = t
                return self.initial + (v*t)
            else: 
                def func(t):
                    return (x0 + Vx*np.fabs(t))**2 + (y0 + Vy*np.fabs(t))**2 -6.25
                t = spo.fsolve(func, 1) 
                exit_point = self.initial + (v*np.fabs(t))
                if exit_point[2] > 100.0:                                                 #if  leaving circular face, then exit is at z=0
                    t = z0/(Vz)
                    self.__exit_time = t
                    return self.initial + (v*t)
                else:
                    self.__exit_time = t
                    return exit_point
        else:                                                                     #if no velocity in z-axis
            if (Vx)**2 + (Vy)**2 == 0.0:                                         #if no velocity at all, it doesn't exit
                return self.initial
            else:
                def func(t):                                                    #if only radial velocity, it exits parallel to x-y plane
                    return (x0 + Vx*np.fabs(t))**2 + (y0 + Vy*np.fabs(t))**2 -6.25
                t = spo.fsolve(func, 3)  
                self.__exit_time = t
                return self.initial + (v*np.fabs(t))
            
    def final_position(self):
        """where the electron leaves the chamber"""
        self.final = self.exit_position()
        return self.final
            
    def energy_deposit(self, x):
        """Calculates the energy deposited in an NaI scintillator by 
        an electron, for every x cm travelled through it. 
        Takes one argument, x."""
        return self.energy()*np.exp(-x/2.6)                                        #N.B: DISTANCE X IN UNITS OF CM.
        

