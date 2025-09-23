import numpy as np
const = {
    "G": 6.6743e-11, # m^3 / kg*yr^2
    "AU": 149597870700.0, #m
}

class Bodies:
    # Class variable
    species = "name"

    # Constructor
    def __init__(self, name, mass, position, period):  #As much data as you want!
        # Instance variables
        self.name = name
        self.mass = mass
        self.position = position
        
        if period == 0:
            self.velocity = np.array([0, 0, 0])
        else:
            v = (6.179/period)*self.position #could just be position divided by period automatically
            self.velocity = v
        
        
# Creating instances of the Planet class from May 13
Sun = Bodies("Sun",1.99e+30, np.array([0,0,0]), 0)

Earth = Bodies("Earth", 5.97217e+24, np.array([-9.346304131972046e+07, -1.201287277696097e+08, 3.241015186724067e+04]), 1)
Venus = Bodies("Venus",4.86731e+24, np.array([-1.230819115060806e+07, -1.088684474163899e+08, -7.913773451080099e+05]), 0.62)
Mercury = Bodies("Mercury",0.330103e+24, np.array([5.220311855303314e+07, -2.208024748280561e+07, -6.568944458960462e+06]), 0.24)
Mars = Bodies("Mars", 0.641691e+24, np.array([-2.435119918362688e+08, 5.275791756072435e+07, 7.100638156769440e+06]), 1.88)
Jupiter = Bodies("Jupiter", 1898.125e+24, np.array([7.524924043986979e+06, 7.662396025180806e+08, -3.346132848762810e+06]), 11.86)
Saturn = Bodies("Saturn", 568.317e+24, np.array([1.424526634558748e+09, -1.559981411335355e+08, -5.400533237864482e+07]), 29.46)
Uranus = Bodies("Uranus", 86.8099e+24, np.array([1.595150229232406e+09, 2.446606849745856e+09, -1.157881742462659e+07]), 84.01)
Neptune = Bodies("Neptune", 102.409e+24, np.array([4.469600768360281e+09, -3.325898805858069e+07, -1.023216026265953e+08]), 164.79)
#Asteroids
Ceres = Bodies("Ceres", 938.416e+18, 2.77, 4.60)
Vesta = Bodies("Vesta", 259.076e+18, 2.36, 3.63)
Pallas = Bodies("Pallas", 204e+18, 2.77, 4.61)
Hygiea = Bodies("Hygiea",87e+18, 3.14, 5.57)
Interamnia = Bodies("Interamnia",35e+18, 3.06, 5.34)
Eunomia = Bodies("Eunnomia", 30e+18, 2.64, 4.30)
Juno = Bodies("Juno", 27e+18, 2.67, 4.36)
Davida = Bodies("Davida", 27e+18, 3.17, 5.61)
Europa = Bodies("Europa", 24e+18, 3.10, 3.55)
Psyche = Bodies("Psyche",23e+18, 2.92, 5.01)

