import numpy as np

class Detector:
    def __init__(self, d= [0.0, 0.0, 0.0], l= 0.3):
        """
        A class of cubic detectors. A detector is defined by two arguments:
            d; the co-ordinates of its centre from the origin. 
            l; the length of its vertices. 
        """
        self.c = np.array(d + l/2)
        self.x_min = np.array(d - l/2)
        self.x_max = np.array(d + l/2)
        self.y_min = np.array(d - l/2)
        self.y_max = np.array(d + l/2)
        self.z_min = np.array(d - 1/2)
        self.z_max = np.array(d + 1/2)
        if len(d) > 3:
            raise Exception("Error: Four Vector parameter size")
        
    def detect(self)

        

#smearing response
res = (0.02*E)/np.sqrt(E)
n = len(res)
#smear using a random number generator
noise = np.random.uniform(0.0,res,n)
E_detected = E + noise
plt.plot(cosO,E,'r+')