import numpy as np

class Fourvector:
    
    """
    A class defining a template for Energy-Momentum four vectors.
    Two arguments define the fourvector; 
    E; the energy component of the fourvector in MeV, and 
    float
    default value = 0.0
    p=[p(x), p(y), p(z)] ; the momentum vector in MeV/c. 
    array like
    default value = [0.00,0.00,0.00]
    
    """
    
    def __init__(self, E=0.00, p=[0.00,0.00,0.00]):
        self.__E = E
        self.__p = np.array(p) #convert to natural units by dividing by c??
        self._fourv = np.append(self.__E,self.__p)
        if len(p) > 3:
            raise Exception("Error: Four Vector parameter size")
        
    def __repr__(self):
        return "%s([E, P] =%r)" % ("Four Vector", self._fourv)
 
    def __str__(self):
        return "[%g, %r]" % (self.__E, self.__p)
    
    def fourvec(self):
        """Returns the full Four Vector."""
        return self._fourv
    
    def energy(self):
          """Returns the Energy of a Four Vector, E."""
          return self.__E
 
    def momentum(self):
         """ Returns the momentum attribute (p) of a Four Vector, 
         which is a vector."""
         return self.__p
       
    def copy(self):
        """The function returns a copy of a Four Vector instance as another 
        independent instance."""
        return Fourvector(self.__E, self.__p)
    
    def __add__(self,other):
        return Fourvector(self.__E+other.__E, self.__p+other.__p)
        
    def __iadd__(self,other):
        self.__E += other.__E
        self.__p += other.__p
        return self
    
    def __sub__(self,other):
        return Fourvector(self.__E-other.__E, self.__p-other.__p)
        
    def __isub__(self,other):
        self.__E -= other.__E
        self.__p -= other.__p
        return self
        
    def inner(self,other):
        """The function takes two four vector arguments, and returns the inner 
        product of the two four vectors. """ 
        return (self.__E*other.__E) - np.dot(self.__p, other.__p)
    
    def magp_sq(self):
        """ The function takes one four vector as an argument and returns the 
        magnitude of its momentum squared"""
        return np.dot(self.__p, self.__p)
        
    def magp(self):
        """ The function takes one four vector as an argument and returns the 
        magnitude of its momentum"""
        return np.sqrt(np.dot(self.__p, self.__p))

    def boost(self, beta=0.0):
        """The function takes two arguments, a four vector and a beta value, 
        (where beta= v/c). The function returns a four vector with a boost in 
        the z-direction."""
        gamma = 1/np.sqrt(1-(beta)**2)
        return Fourvector(gamma*(self._fourv[0]-beta*self._fourv[3]), [self._fourv[1], self._fourv[2], gamma*(self._fourv[3]-beta*self._fourv[0])])
        #lorentz transform, does it still apply for p if you use E instead of E/c????
        
        
    