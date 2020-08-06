#!/usr/bin/python3
'''
Simon Kueppers
7/29/2020

Mike Dawson-Haggerty
1/30/2014

Inspired by:
http://www.zincland.com/hypocycloid

Implemented from:
'Gear geometry of cycloid drives', Chen BingKui et al.
http://download.springer.com/static/pdf/96/art%253A10.1007%252Fs11431-008-0055-3.pdf
Equations 10, 11
'''

import numpy as np
import matplotlib.pyplot as plt
from math import *
from mathutils import *
import ezdxf

class CycloidalDrive:
    def __init__(self, count_pin, count_cam, eccentricity, radius_pin, radius_pattern, resolution, reduction_tolerance):
        self.count_pin = count_pin
        self.count_cam = count_cam
        self.eccentricity = eccentricity
        self.radius_pin = radius_pin
        self.radius_pattern = radius_pattern
        self.resolution = resolution
        self.reduction_tolerance = radians(reduction_tolerance)

    def reduce_polyline(self, polyline):
        vertices = [[],[]]
        last_vertex = [polyline[0][0], polyline[1][0]]

        # Look through all vertices except start and end vertex
        # Calculate by how much the lines before and after the vertex
        # deviate from a straight path.
        # If the deviation angle exceeds the specification, store it
        for vertex_idx in range(1, len(polyline[0])-1):
            next_slope = np.arctan2(    polyline[1][vertex_idx+1] - polyline[1][vertex_idx+0],
                                        polyline[0][vertex_idx+1] - polyline[0][vertex_idx+0]   )
            prev_slope = np.arctan2(    polyline[1][vertex_idx-0] - last_vertex[1],
                                        polyline[0][vertex_idx-0] - last_vertex[0]   )

            deviation_angle = abs(prev_slope - next_slope)

            if (deviation_angle > self.reduction_tolerance):
                vertices[0] += [polyline[0][vertex_idx]]
                vertices[1] += [polyline[1][vertex_idx]]
                last_vertex = [polyline[0][vertex_idx], polyline[1][vertex_idx]]

        # Return vertices along with first and last point of the original polyline
        return np.array([
            np.concatenate([ [polyline[0][0]], vertices[0], [polyline[0][-1]] ]),
            np.concatenate([ [polyline[1][0]], vertices[1], [polyline[1][-1]] ])
        ])

    def generate_cam(self):
        Rz = self.radius_pattern  # radius of pin pattern
        rz = self.radius_pin      # radius of pin
        e =  self.eccentricity    # eccentricity
        Zb = self.count_pin       # number of pins
        Zg = self.count_cam       # tooth count on gear

        Ze = Zb / (Zb - Zg)
        Zd = Zg / (Zb - Zg)
        K1 = (e * Zb) / (Rz * (Zb - Zg))

        # in the paper they say you should calculate this numerically...
        psi = np.linspace(0, np.pi * 2, int(360 * self.resolution), endpoint=False)

        denom_B = np.sqrt(1 + K1**2 - 2 * K1 * np.cos(Zd * psi))
        cos_B = np.sign(Zb - Zg) * ((K1 * np.sin(Ze * psi)) -
                                    np.sin(psi)) / denom_B
        sin_B = np.sign(Zb - Zg) * ((-K1 * np.cos(Ze * psi)) +
                                    np.cos(psi)) / denom_B

        x = Rz * np.sin(psi) - e * np.sin(Ze * psi) + rz * cos_B
        y = Rz * np.cos(psi) - e * np.cos(Ze * psi) - rz * sin_B

        return self.reduce_polyline(np.array((x, y + self.eccentricity)))

    def generate_pins(self):
        psi = np.linspace(0, np.pi * 2, int(360 * self.resolution), endpoint=False)
        x = self.radius_pin * np.sin(psi)
        y = self.radius_pin * np.cos(psi) + self.radius_pattern
        pin = self.reduce_polyline(np.array((x , y)))
        pins = [ np.dot(rotation_matrix(np.pi * 2 / self.count_pin * n), pin) for n in range(0, self.count_pin) ]
        return np.array(pins)

    def get_dxf(self):
        doc = ezdxf.new('R2010')
        doc.header['$MEASUREMENT'] = 1
        doc.header['$INSUNITS'] = 4
        msp = doc.modelspace()

        doc.layers.new(name='cam')
        msp.add_lwpolyline(self.generate_cam().transpose(), dxfattribs = {'closed': True, 'layer': 'cam'})
        msp.add_point( (0, self.eccentricity), dxfattribs = {'layer': 'cam'})

        doc.layers.new(name='pins')
        for n in range(0, self.count_pin):
            msp.add_circle(
                (0, self.radius_pattern),
                self.radius_pin, dxfattribs = {'layer': 'pins'}
            ).transform(ezdxf.math.Matrix44.z_rotate(np.pi * 2 / self.count_pin * n))

        return doc

if __name__ == '__main__':
    drive = CycloidalDrive(   count_pin = 34,
                              count_cam = 33,
                              eccentricity = 1,
                              radius_pin = 2,
                              radius_pattern = 37,
                              resolution = 100,
                              reduction_tolerance = 2.0)

    cam = drive.generate_cam()
    pins = drive.generate_pins()

    plt.axis('equal')
    plt.grid(True)
    plt.plot(cam[0,:], cam[1,:])
    for pin in pins:
        plt.plot(pin[0,:], pin[1,:], 'r')
    plt.show()

    drive.get_dxf().saveas('cycloidal.dxf')
