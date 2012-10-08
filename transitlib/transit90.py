# -*- coding: utf-8 -*-

import scipy
import scipy.integrate
import scipy.constants

class Transit90(object):
    
    orbit_radius        = None
    star_radius         = None
    planet_radius       = None
    star_temperature    = None
    planet_temperature  = None
    star_darkening      = None
    
    phase_start     = float(0)
    phase_end       = float(1)
    phase_step      = float(0.0001)
    _phase          = float(0)
    completed       = 0

    phase1 = float(0)
    phase2 = float(0)
    phase6 = float(0)
    stopped = False
    
    def __init__(self):
        super(Transit90, self).__init__()
        return
        
    def set_orbit_radius(self, radius):
        self.orbit_radius = float(radius)
        
    def set_star_radius(self, radius):
        self.star_radius = float(radius)
        
    def set_planet_radius(self, radius):
        self.planet_radius = float(radius)
        
    def set_star_temperature(self, temperature):
        self.star_temperature = float(temperature)
        
    def set_planet_temperature(self, temperature):
        self.planet_temperature = float(temperature)
        
    def set_star_darkening(self, darkening):
        self.star_darkening = float(darkening)

    def set_phase_start(self, start):
        self.phase_start = float(start)
        
    def set_phase_end(self, end):
        self.phase_end = float(end)
        
    def set_phase_step(self,step):
        self.phase_step = float(step)

    def stop(self):
        self.stopped = True
    
    def star_intensity(self):
        return scipy.constants.sigma * self.star_temperature**4
        
    def planet_intensity(self):
        return scipy.constants.sigma * self.planet_temperature**4
        
    def star_luminosity(self):
        return scipy.pi * self.star_radius**2 * self.star_intensity() * ( 1 - (self.star_darkening/3))
        
    def planet_luminosity(self):
        return scipy.pi * self.planet_radius**2 * self.planet_intensity()
        
    def x0(self):
        return self.orbit_radius * scipy.sin(2 * scipy.pi * self._phase)
        
    def x1(self, x):
        x0 = self.x0()

        if x0 == 0 :
            return 0
        
        return (x0**2 + x**2 - self.planet_radius**2)/(2*x0)
        
    def gamma(self, x):
        return scipy.arccos(self.x1(x)/x)
        
    def eq1(self, x):
        return self.star_intensity() * 2 * x * self.gamma(x) * ( 1 - self.star_darkening + self.star_darkening * scipy.sqrt((1-(x/self.star_radius)**2)))
        
    def eq2(self, x):
        return self.star_intensity() * 2 * scipy.pi * x * ( 1 - self.star_darkening + self.star_darkening * scipy.sqrt((1-(x/self.star_radius)**2)))
        
    def onStop(self):
        pass
    
    def onProgress(self, progress):
        pass
    
    def onComplete(self, phases, values):
        pass
        
    def run(self):

        if self.orbit_radius is None : raise Exception("Invalid orbit radius!")
        if self.star_radius is None : raise Exception("Invalid star radius!")
        if self.planet_radius is None : raise Exception("Invalid planet radius!")
        if self.star_temperature is None : raise Exception("Invalid start temperature!")
        if self.planet_temperature is None : raise Exception("Invalid planet temperature!")
        if self.star_darkening is None : raise Exception("Invalid star darkening coeficient!")
        
        self.phase1 = (1/(2*scipy.pi)) * scipy.arcsin((self.star_radius+self.planet_radius)/self.orbit_radius)
        self.phase2 = (1/(2*scipy.pi)) * scipy.arcsin((self.star_radius-self.planet_radius)/self.orbit_radius)
        self.phase6 = (1/(2*scipy.pi)) * scipy.arcsin(self.planet_radius/self.orbit_radius)

        maxPhase = max(self.phase1, self.phase2, self.phase6)
        self._phase = self.phase_start

        result = []
        phases = []
        
        while( self._phase <= self.phase_end and self._phase <= maxPhase ):
            
            if self.stopped:
                self.onStop()
                return
                
            current_result = None
            if self._phase >= min(self.phase1, self.phase2) and self._phase <= max(self.phase1, self.phase2) :
                a = self.x0() - self.planet_radius
                b = self.star_radius
                current_result = scipy.integrate.quad(self.eq1, a, b)
                current_result = current_result[0]
            elif self._phase >= min(self.phase2, self.phase6) and self._phase <= max(self.phase2, self.phase6):
                a = self.x0() - self.planet_radius
                b = self.x0() + self.planet_radius
                current_result = scipy.integrate.quad(self.eq1, a, b)
                current_result = current_result[0]
            elif self._phase <= self.phase6:
                a1 = 0
                b1 = self.planet_radius - self.x0() 
                a2 = self.planet_radius - self.x0()
                b2 = self.x0() + self.planet_radius
                
                result1 = scipy.integrate.quad(self.eq2, a1, b1)
                result2 = scipy.integrate.quad(self.eq1, a2, b2)
                current_result =  result1[0] + result2[0]
            else:
                current_result = None
                
                            
            if current_result is not None :
                phases.append(self._phase)
                result.append( 1 - ( -2.5*scipy.log10(((self.star_luminosity()+self.planet_luminosity())-current_result)/(self.star_luminosity()+self.planet_luminosity()))) )
            

            self._phase += self.phase_step

            newcompleted = min(int((self._phase/maxPhase).real*100),100)
            
            if newcompleted != self.completed :
                self.completed = newcompleted
                self.onProgress(self.completed)
        
        self.onComplete(phases, result)