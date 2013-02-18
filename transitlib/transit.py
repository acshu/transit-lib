# -*- coding: utf-8 -*-
#

from numpy import radians, degrees, sin, cos, arcsin, arctan, finfo
from numpy import arccos as acos
from scipy import sqrt, log
from scipy.integrate import quad
from scipy.constants import pi, sigma
import time, sys, math


class Transit(object):

    def __init__(self):
        super(Transit, self).__init__()
        self.semi_major_axis = None
        self.star_radius = None
        self.planet_radius = None
        self.star_temperature = None
        self.planet_temperature = None
        self.star_darkening_1 = None
        self.star_darkening_2 = None
        self.inclination = None
        self.phases_injection = []
        self.star_darkening_type = 'linear'
        self.star_darkening_law = None

        self.phase_start = float(0)
        self.phase_end = float(1)
        self.phase_step = float(0.0001)
        self.precision = 1e-10
        self.completed = 0

        self.stopped = False
        pass

    def set_semi_major_axis(self, semi_major_axis):
        self.semi_major_axis = float(semi_major_axis)

    def set_star_radius(self, radius):
        self.star_radius = float(radius)

    def set_planet_radius(self, radius):
        self.planet_radius = float(radius)

    def set_star_temperature(self, temperature):
        self.star_temperature = float(temperature)

    def set_planet_temperature(self, temperature):
        self.planet_temperature = float(temperature)

    def set_inclination(self, inclination):
        self.inclination = float(inclination)
        self.inclination_rad = radians(self.inclination)
        
    def set_star_darkening_type(self, type):
        self.star_darkening_type = type

    def set_star_darkening_1(self, darkening):
        self.star_darkening_1 = float(darkening)
        
    def set_star_darkening_2(self, darkening):
        self.star_darkening_2 = float(darkening)

    def set_phase_start(self, start):
        self.phase_start = float(start)

    def set_phase_end(self, end):
        self.phase_end = float(end)

    def set_phase_step(self, step):
        self.phase_step = float(step)

    def set_precision(self, precision):
        self.precision = precision

    def set_phases_injection(self, phases):
        self.phases_injection = phases

    def stop(self):
        self.stopped = True

    def onStop(self):
        pass

    def onProgress(self, progress):
        pass

    def onComplete(self, phases, values):
        pass

    def quad(self, integral, a, b, args=None):
        if a == b or abs(a - b) < finfo(float).eps:
            return 0, 0
        return quad(integral, a, b, args=args, epsrel=self.precision, epsabs=0)

    def run(self):
        start = time.time()
        rs = self.star_radius / self.semi_major_axis
        rp = self.planet_radius / self.semi_major_axis
        k = (self.planet_temperature / self.star_temperature) ** 4

        i1 = acos(rs + rp)
        i2 = acos(rs - rp)
        i5 = acos(rp)

        #Phase limits
        f1 = (1 / (2 * pi)) * arcsin(sqrt((rs + rp) ** 2 - cos(self.inclination_rad) ** 2) / sin(self.inclination_rad))
        f2 = (1 / (2 * pi)) * arcsin(sqrt((rs - rp) ** 2 - cos(self.inclination_rad) ** 2) / sin(self.inclination_rad))
        f5 = (1 / (2 * pi)) * arctan(cos(self.inclination_rad))
        f6 = (1 / (2 * pi)) * arcsin(sqrt(rp ** 2 - cos(self.inclination_rad) ** 2) / sin(self.inclination_rad))
        
        
        phases = []
        phi = self.phase_start
        while phi < self.phase_end :
            phases.append(phi)
            phi += self.phase_step

        if len(self.phases_injection):
            for phase_injection in self.phases_injection:
                phases.append(math.fabs(phase_injection))


        # final point may missing due delta
        phases.append(self.phase_end)
        # make unique and sort phases
        phases = sorted(list(set(phases)))

        # darkening laws
        star_luminosity_darkening_law = None
           
        if self.star_darkening_type == 'linear' :
            self.star_darkening_law = lambda q: 1 - self.star_darkening_1 + self.star_darkening_1 * sqrt(1 - (q / rs) ** 2)
            star_luminosity_darkening_law = 1 - self.star_darkening_1 / 3
            
        elif self.star_darkening_type == 'quadratic' :
            self.star_darkening_law = lambda q: 1 - self.star_darkening_1 * (1 - sqrt(1 - (q / rs) ** 2)) - self.star_darkening_2 * (1 - sqrt(1 - (q / rs) ** 2)) ** 2
            star_luminosity_darkening_law = 1 - self.star_darkening_1 / 3 - self.star_darkening_2 / 6
            
        elif self.star_darkening_type == 'squareroot' :
            self.star_darkening_law = lambda q: 1 - self.star_darkening_1 * (1 - sqrt(1 - (q / rs) ** 2)) - self.star_darkening_2 * (1 - (1 - (q / rs) ** 2) ** (1.0 / 4.0))
            star_luminosity_darkening_law = 1 - self.star_darkening_1 / 3 - self.star_darkening_2 / 5
            
        elif self.star_darkening_type == 'logarithmic' :
            self.star_darkening_law = lambda q: 1 - self.star_darkening_1 * (1 - sqrt(1 - (q / rs) ** 2)) - self.star_darkening_2 * sqrt(1 - (q / rs) ** 2) * log(sqrt(1 - (q / rs) ** 2))
            star_luminosity_darkening_law = 1 - self.star_darkening_1 / 3 - 2 * self.star_darkening_2 / 9
            
        
        star_intensity = sigma * self.star_temperature ** 4
        planet_intensity = sigma * self.planet_temperature ** 4
        
        star_luminosity = pi * self.star_radius ** 2 * star_intensity * star_luminosity_darkening_law
        planet_luminosity = pi * self.planet_radius ** 2 * planet_intensity
        

        int1 = lambda q, rs, b, x0, y0, rp : self.star_darkening_law(q) * q * (acos((x0 * (b ** 2 + q ** 2 - rp ** 2) - y0 * sqrt(4 * b ** 2 * q ** 2 - (b ** 2 + q ** 2 - rp ** 2) ** 2)) / (2 * b ** 2 * q)) - acos((x0 * (b ** 2 + q ** 2 - rp ** 2) + y0 * sqrt(4 * b ** 2 * q ** 2 - (b ** 2 + q ** 2 - rp ** 2) ** 2)) / (2 * b ** 2 * q)))
        int2 = lambda q, rs, b, x0, y0, rp : self.star_darkening_law(q) * q * (acos((x0 * (b ** 2 + q ** 2 - rp ** 2) - y0 * sqrt(4 * b ** 2 * q ** 2 - (b ** 2 + q ** 2 - rp ** 2) ** 2)) / (2 * b ** 2 * q)) + acos((x0 * (b ** 2 + q ** 2 - rp ** 2) + y0 * sqrt(4 * b ** 2 * q ** 2 - (b ** 2 + q ** 2 - rp ** 2) ** 2)) / (2 * b ** 2 * q)))
        int3 = lambda q, rs, b, x0, y0, rp : self.star_darkening_law(q) * q * ((2 * pi) - acos(abs((x0 * (b ** 2 + q ** 2 - rp ** 2) + y0 * sqrt(4 * b ** 2 * q ** 2 - (b ** 2 + q ** 2 - rp ** 2) ** 2)) / (2 * b ** 2)) / q) + acos(abs((x0 * (b ** 2 + q ** 2 - rp ** 2) - y0 * sqrt(4 * b ** 2 * q ** 2 - (b ** 2 + q ** 2 - rp ** 2) ** 2)) / (2 * b ** 2)) / q))
        int4 = lambda q, rs, b, x0, y0, rp : self.star_darkening_law(q) * q * (pi + acos(abs((x0 * (b ** 2 + q ** 2 - rp ** 2) + y0 * sqrt(4 * b ** 2 * q ** 2 - (b ** 2 + q ** 2 - rp ** 2) ** 2)) / (2 * b ** 2)) / q) + acos(abs((x0 * (b ** 2 + q ** 2 - rp ** 2) - y0 * sqrt(4 * b ** 2 * q ** 2 - (b ** 2 + q ** 2 - rp ** 2) ** 2)) / (2 * b ** 2)) / q))

        mag = []
        iteration = 0
        for phase in phases:

            if self.stopped:
                self.onStop()
                return

            result = 1
            intResult = None

            if phase >= 0 and phase < f1:
                if self.inclination_rad > i1 and self.inclination_rad <= i2:
                    x0 = sin(2 * pi * phase)
                    y0 = cos(self.inclination_rad) * cos(2 * pi * phase)
                    b = sqrt(x0 ** 2 + y0 ** 2)
                    qmin = b - rp
                    qmax = rs

                    if y0 >= rp:
                        integr = int1
                    else:
                        integr = int2

                    intResult, err = self.quad(integr, qmin, qmax, args=(rs, b, x0, y0, rp))
                    result = intResult

            if self.inclination_rad > i2 and self.inclination_rad <= i5:
                if phase >= 0 and phase < f1:
                    x0 = sin(2 * pi * phase)
                    y0 = cos(self.inclination_rad) * cos(2 * pi * phase)
                    b = sqrt(x0 ** 2 + y0 ** 2)
                    if phase <= f2:
                        qmin = b - rp
                        qmax = b + rp
                        if y0 >= rp:
                            integr = int1
                            intResult, err = self.quad(integr, qmin, qmax, args=(rs, b, x0, y0, rp))
                        else:
                            q1 = x0 + sqrt(rp ** 2 - y0 ** 2)
                            q2 = x0 - sqrt(rp ** 2 - y0 ** 2)
                            j1, err = self.quad(int1, qmin, q2, args=(rs, b, x0, y0, rp))
                            j2, err = self.quad(int2, q2, q1, args=(rs, b, x0, y0, rp))
                            j3, err = self.quad(int1, q1, qmax, args=(rs, b, x0, y0, rp))
                            intResult = (j1 + j2 + j3)
                        
                    else:
                        qmin = b - rp
                        qmax = rs
                        integr = int1
                        intResult, err = self.quad(integr, qmin, qmax, args=(rs, b, x0, y0, rp))

                    result = intResult

            if self.inclination_rad > i5 and self.inclination_rad <= radians(90):
                if f5 < f6:
                    limit1 = f5
                    limit2 = f6
                else:
                    limit1 = f6
                    limit2 = 0

                if phase >= 0 and phase < f1:
                    x0 = sin(2 * pi * phase)
                    y0 = cos(self.inclination_rad) * cos(2 * pi * phase)
                    b = sqrt(x0 ** 2 + y0 ** 2)

                    if phase < limit1:
                        q1 = sqrt(rp ** 2 - y0 ** 2) + x0
                        q4 = sqrt(rp ** 2 - y0 ** 2) - x0
                        q2 = sqrt(rp ** 2 - x0 ** 2) + y0
                        q3 = sqrt(rp ** 2 - x0 ** 2) - y0

                        j1, err = quad(lambda q: 2 * pi * q * self.star_darkening_law(q), 0, (rp - b))
                        j2, err = self.quad(int3, (rp - b), q3, args=(rs, b, x0, y0, rp))
                        j3, err = self.quad(int4, q3, q4, args=(rs, b, x0, y0, rp))
                        j4, err = self.quad(int2, q4, q1, args=(rs, b, x0, y0, rp))
                        j5, err = self.quad(int1, q1, (rp + b), args=(rs, b, x0, y0, rp))
                        intResult = (j1 + j2 + j3 + j4 + j5)

                    elif phase < limit2:
                        q1 = sqrt(rp ** 2 - y0 ** 2) + x0
                        q3 = sqrt(rp ** 2 - x0 ** 2) - y0
                        q4 = sqrt(rp ** 2 - y0 ** 2) - x0

                        j1, err = quad(lambda q: 2 * pi * q * self.star_darkening_law(q), 0, (rp - b))
                        j2, err = self.quad(int3, (rp - b), q4, args=(rs, b, x0, y0, rp))
                        j4, err = self.quad(int2, q4, q1, args=(rs, b, x0, y0, rp))
                        j5, err = self.quad(int1, q1, (rp + b), args=(rs, b, x0, y0, rp))
                        intResult = (j1 + j2 + j4 + j5)
                        # nan

                    else:
                        # f6<= ph <f1
                        qmin = b - rp
                        if phase <= f2:
                            # f6<= ph <=f2
                            qmin = b - rp
                            qmax = b + rp
                            if y0 >= rp:
                                integr = int1
                                intResult, err = self.quad(integr, qmin, qmax, args=(rs, b, x0, y0, rp))

                            else:
                                q1 = x0 + sqrt(rp ** 2 - y0 ** 2)
                                q2 = x0 - sqrt(rp ** 2 - y0 ** 2)

                                j1, err = self.quad(int1, qmin, q2, args=(rs, b, x0, y0, rp))
                                j2, err = self.quad(int2, q2, q1, args=(rs, b, x0, y0, rp))
                                j3, err = self.quad(int1, q1, qmax, args=(rs, b, x0, y0, rp))
                                intResult = (j1 + j2 + j3)
                                #nan

                        else:
                            #f2 < ph < f1
                            qmin = b - rp
                            qmax = rs
                            if y0 >= rp:
                                integr = int1
                                intResult, err = self.quad(integr, qmin, qmax, args=(rs, b, x0, y0, rp))
                            else:
                                q2 = x0 - sqrt(rp ** 2 - y0 ** 2)
                                if q2 >= rs:
                                    integr = int1
                                    intResult, err = self.quad(integr, qmin, qmax, args=(rs, b, x0, y0, rp))
                                else:
                                    q1 = x0 + sqrt(rp ** 2 - y0 ** 2)
                                    if q1 < qmax:
                                        q1 = x0 + sqrt(rp ** 2 - y0 ** 2)
                                        q2 = x0 - sqrt(rp ** 2 - y0 ** 2)

                                        j1, err = self.quad(int1, qmin, q2, args=(rs, b, x0, y0, rp))
                                        j2, err = self.quad(int2, q2, q1, args=(rs, b, x0, y0, rp))
                                        j3, err = self.quad(int1, q1, qmax, args=(rs, b, x0, y0, rp))
                                        intResult = (j1 + j2 + j3)

                                    else:

                                        j1, err = self.quad(int1, qmin, q2, args=(rs, b, x0, y0, rp))
                                        j2, err = self.quad(int2, q2, qmax, args=(rs, b, x0, y0, rp))
                                        intResult = (j1 + j2)

            # final Flux
            if intResult is not None :
                result = 1 - intResult / ( pi * k * rp ** 2 + (star_luminosity / (star_intensity * self.semi_major_axis ** 2 )))
            else:
                result = 1
                        
            mag.append(result)
            iteration += 1
            self.completed = min(int((float(iteration) / len(phases)).real * 100), 100)
            self.onProgress(self.completed)

        #print time.time() - start, 's'
        self.onComplete(phases, mag)
