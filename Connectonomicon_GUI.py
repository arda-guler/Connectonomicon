## Connectonomicon
## Copyright (C) 2024  H. A. Guler
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see
## <https://www.gnu.org/licenses/>.

import numpy as np
import spiceypy as spice
from datetime import datetime, timedelta
import re
import sys
import matplotlib.pyplot as plt
import subprocess
import os
import mplcursors

import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk

title = "== == == Connectonomicon == == =="
version = "v20241215.00"
authors = "Authors: H. A. Guler\n\nConnectonomicon uses find_orb for orbit determination, developed by Bill Gray. (https://www.projectpluto.com/)"
license_text = "Connectonomicon is licensed under GNU General Public License v2.\nfind_orb is licensed under GNU General Public License v2, see https://www.projectpluto.com/find_orb.htm#License.\n" +\
               "SPICE kernels are not shared alongside this tool, but their use may be credited as described in https://naif.jpl.nasa.gov/naif/credit.html"

help_text =\
"""
Connectonomicon is a tool that helps with the identification of separate astrometric measurements belonging to the same object. Below is a short version of README.txt, read that file for more info.

Connectonomicon makes use of find_orb for orbit determination. It then uses its gravitationally perturbed orbit propagator to predict the future trajectory of the object. """ +\
"""There is also a built-in (native) orbit determination algorithm based on the Laplace and Herget methods, but it is not as robust as find_orb and is quite a lot slower. """ +\
"""For short arcs, it also attempts to figure out the possible area where other sightings at a given date may be, accounting for possible measurement errors. """ +\
"""Lastly, it classifies potential related astrometric observations from a given list of observations.

Connectonomicon requires the console version of find_orb (fo64.exe) to be in the same directory as itself. You can get fo64.exe from the following URL on Project Pluto website:

https://www.projectpluto.com/fo_usage.htm

You need to define observations in Obs80 (MPC 80-column format) in two files - a list of known observations which will be used for the orbit fit and a list of undetermined observations to search within.

You also need SPICE kernels for the orbit propagator, which you can get through NAIF (https://naif.jpl.nasa.gov/naif/). Place them in a 'data' folder:

data/de440.bsp
data/naif0012.tls
"""

# Constants
km_per_AU = 149597870.7
AU = 149597870.7
km_per_mAU = 149597870.7 * 1e-3
s_per_d = 86400
day = 86400
arcsecond = 0.00027778 # deg

## CLASS DEFINITIONS
class MainBody: # for the Sun and planets
    def __init__(self, name, pos, GM):
        self.name = name
        self.pos = pos
        self.GM = GM # a.k.a. mu, standard gravitational parameter

class MP: # minor planets, comets
    def __init__(self, des, pos, vel):
        self.des = des
        self.pos = pos
        self.vel = vel

class Obs: # MPC 80-column observations
    def __init__(self, obs_str, debug=False):
        if type(obs_str) == str:
            self.obs_str = obs_str
            
            packed_perm = obs_str[0:5]
            packed_prov = obs_str[5:12]
            date_str = obs_str[15:31]
            RA_str = obs_str[32:43]
            DEC_str = obs_str[44:55]
            mag_str = obs_str[65:69]
            obscode_str = obs_str[77:80]

            # date handling
            obs_date_parts = date_str.split()
            year = int(obs_date_parts[0])
            month = int(obs_date_parts[1])
            day = int(float(obs_date_parts[2]))

            decimal_day = float(obs_date_parts[2]) - day
            day_seconds = 86400
            decimal_secs = day_seconds * decimal_day

            obs_datetime = datetime(year, month, day)
            obs_datetime = obs_datetime + timedelta(seconds=decimal_secs)

            # RA - DEC handling
            hour2deg = 360/24
            minute2deg = 360/(24*60)
            second2deg = 360/(24*60*60)

            RA_parts = RA_str.split()
            RA_deg = float(RA_parts[0]) * hour2deg + float(RA_parts[1]) * minute2deg + float(RA_parts[2]) * second2deg

            DEC_parts = DEC_str.split()
            DEC_deg = (abs(float(DEC_parts[0])) + float(DEC_parts[1]) / 60 + float(DEC_parts[2]) / 3600) * sign(float(DEC_parts[0]))

            # mag handling
            try:
                mag = float(mag_str)
            except ValueError:
                mag = None

            self.perm = packed_perm
            self.prov = packed_prov
            self.date = obs_datetime
            self.RA = RA_deg
            self.DEC = DEC_deg
            self.mag = mag
            self.obs_code = obscode_str
            
        elif type(obs_str) == list:
            self.RA = obs_str[0]
            self.DEC = obs_str[1]
            self.date = obs_str[2]

    def __str__(self):
        return self.obs_str

    def __repr__(self):
        return str(self)

class ObsPair: # A pair is any 2 astrometric observations of the same object
    def __init__(self, o1, o2, pix_error=1, deg_per_pixel=9.793873680970562e-05):
        self.o1 = o1
        self.o2 = o2

        self.delta_t = (o2.date - o1.date).total_seconds()
        self.delta_RA = (o2.RA - o1.RA)
        self.delta_DEC = (o2.DEC - o1.DEC)
        self.RA_vel = self.delta_RA / self.delta_t
        self.DEC_vel = self.delta_DEC / self.delta_t

        self.main_dir = [self.delta_RA, self.delta_DEC]

        # some "numerical" error finding
        self.o1_p = [[self.o1.RA + deg_per_pixel * pix_error, self.o1.DEC],
                      [self.o1.RA - deg_per_pixel * pix_error, self.o1.DEC],
                      [self.o1.RA, self.o1.DEC + deg_per_pixel * pix_error],
                      [self.o1.RA, self.o1.DEC - deg_per_pixel * pix_error],
                      [self.o1.RA + deg_per_pixel * pix_error, self.o1.DEC + deg_per_pixel * pix_error],
                      [self.o1.RA + deg_per_pixel * pix_error, self.o1.DEC - deg_per_pixel * pix_error],
                      [self.o1.RA - deg_per_pixel * pix_error, self.o1.DEC + deg_per_pixel * pix_error],
                      [self.o1.RA - deg_per_pixel * pix_error, self.o1.DEC - deg_per_pixel * pix_error]]

        self.o2_p = [[self.o2.RA + deg_per_pixel * pix_error, self.o2.DEC],
                      [self.o2.RA - deg_per_pixel * pix_error, self.o2.DEC],
                      [self.o2.RA, self.o2.DEC + deg_per_pixel * pix_error],
                      [self.o2.RA, self.o2.DEC - deg_per_pixel * pix_error],
                      [self.o2.RA + deg_per_pixel * pix_error, self.o2.DEC + deg_per_pixel * pix_error],
                      [self.o2.RA + deg_per_pixel * pix_error, self.o2.DEC - deg_per_pixel * pix_error],
                      [self.o2.RA - deg_per_pixel * pix_error, self.o2.DEC + deg_per_pixel * pix_error],
                      [self.o2.RA - deg_per_pixel * pix_error, self.o2.DEC - deg_per_pixel * pix_error]]

        all_lines = []
        all_dists = []
        for i in range(8):
            for j in range(8):
                c_line = [self.o2_p[j][0] - self.o1_p[i][0], self.o2_p[j][1] - self.o1_p[i][1]]
                all_lines.append(c_line)
                all_dists.append([abs(self.o2_p[j][0] - self.o1_p[i][0]), abs(self.o2_p[j][1] - self.o1_p[i][1])])

        # find max difference in slope
        max_diff = 0

        if self.delta_RA != 0:
            main_slope = self.delta_DEC / self.delta_RA
        else:
            main_slope = 1e10 # very large finite number
            
        for i in range(len(all_lines)):
            if all_lines[i][0] != 0:
                slope = all_lines[i][1] / all_lines[i][0]
            else:
                slope = 1e10 # very large finite number
                
            slope_diff = abs(slope - main_slope)
            if slope_diff > max_diff:
                max_diff = slope_diff

        min_RA_vel = float("Inf")
        max_RA_vel = -float("Inf")
        min_DEC_vel = float("Inf")
        max_DEC_vel = -float("Inf")
        for i in range(len(all_dists)):
            if all_dists[i][0] / self.delta_t < min_RA_vel:
                min_RA_vel = all_dists[i][0] / self.delta_t
            if all_dists[i][0] / self.delta_t > max_RA_vel:
                max_RA_vel = all_dists[i][0] / self.delta_t

            if all_dists[i][1] / self.delta_t < min_DEC_vel:
                min_DEC_vel = all_dists[i][1] / self.delta_t
            if all_dists[i][1] / self.delta_t > max_DEC_vel:
                max_DEC_vel = all_dists[i][1] / self.delta_t

        self.max_vel = (max_RA_vel**2 + max_DEC_vel**2)**0.5
        self.min_vel = (min_RA_vel**2 + min_DEC_vel**2)**0.5
        self.slope_diff = max_diff
        self.main_slope = main_slope

    def __str__(self):
        output = str(o1) + "\n" + str(o2) + "\n\n"
        output += "Max. slope diff: " + str(self.slope_diff) + "\n"
        output += "Max. vel: " + str(self.max_vel) + " deg s-1,\tMin. vel: " + str(self.min_vel) + " deg s-1\n"
        return output

    def __repr__(self):
        return str(self)

    def guessPositions(self, new_dates, plot=True):
        # sanitize dates
        new_delta_ts = []
        for new_date in new_dates:
            obs_date_parts = new_date.split()
            year = int(obs_date_parts[0])
            month = int(obs_date_parts[1])
            day = int(float(obs_date_parts[2]))

            decimal_day = float(obs_date_parts[2]) - day
            day_seconds = 86400
            decimal_secs = day_seconds * decimal_day

            obs_datetime = datetime(year, month, day)
            obs_datetime = obs_datetime + timedelta(seconds=decimal_secs)

            new_delta_t = (obs_datetime - self.o1.date).total_seconds()
            new_delta_ts.append(new_delta_t)

        position_sets = []
        for i in range(len(new_delta_ts)):
            dt = new_delta_ts[i]
            p_center = [self.o1.RA + self.RA_vel * dt, self.o1.DEC + self.DEC_vel * dt]

            max_slope = self.main_slope + self.slope_diff
            min_slope = self.main_slope - self.slope_diff
            max_normalized_RA_rate = 1 / (1 + max_slope**2)**0.5
            max_normalized_DEC_rate = max_slope / (1 + max_slope**2)**0.5
            min_normalized_RA_rate = 1 / (1 + min_slope**2)**0.5
            min_normalized_DEC_rate = min_slope / (1 + min_slope**2)**0.5

            if max_normalized_RA_rate * self.RA_vel < 0:
                max_normalized_RA_rate = -max_normalized_RA_rate
                min_normalized_RA_rate = -min_normalized_RA_rate
                max_normalized_DEC_rate = -max_normalized_DEC_rate
                min_normalized_DEC_rate = -min_normalized_DEC_rate
            
            p1 = [self.o1.RA + max_normalized_RA_rate * self.max_vel * dt, self.o1.DEC + max_normalized_DEC_rate * self.max_vel * dt]
            p2 = [self.o1.RA + max_normalized_RA_rate * self.min_vel * dt, self.o1.DEC + max_normalized_DEC_rate * self.min_vel * dt]
            p3 = [self.o1.RA + min_normalized_RA_rate * self.max_vel * dt, self.o1.DEC + min_normalized_DEC_rate * self.max_vel * dt]
            p4 = [self.o1.RA + min_normalized_RA_rate * self.min_vel * dt, self.o1.DEC + min_normalized_DEC_rate * self.min_vel * dt]

            position_sets.append([p_center, p1, p2, p3, p4])

        plt.scatter([o1.RA, o2.RA], [o1.DEC, o2.DEC], color="b", label="Previous Observations")
        # plt.plot([o1.RA, o2.RA], [o1.DEC, o2.DEC], "--")

        for i in range(len(new_delta_ts)):
            p_center = position_sets[i][0]
            p1 = position_sets[i][1]
            p2 = position_sets[i][2]
            p3 = position_sets[i][3]
            p4 = position_sets[i][4]
            
            plt.scatter([p_center[0]], [[p_center[1]]], color="r", label="Perfect Observation")
            plt.plot([p1[0], p2[0], p4[0], p3[0], p1[0]], [p1[1], p2[1], p4[1], p3[1], p1[1]], color="yellow", label="Prediction Range")
            plt.plot([o1.RA, p_center[0]], [o1.DEC, p_center[1]], "--")

        plt.grid()
        plt.show()

        return position_sets

    def checkDateObs(self, new_date, plot=True, checked_obs=None):
        dt = (new_date - self.o1.date).total_seconds()

        p_center = [self.o1.RA + self.RA_vel * dt, self.o1.DEC + self.DEC_vel * dt]

        max_slope = self.main_slope + self.slope_diff
        min_slope = self.main_slope - self.slope_diff
        max_normalized_RA_rate = 1 / (1 + max_slope**2)**0.5
        max_normalized_DEC_rate = max_slope / (1 + max_slope**2)**0.5
        min_normalized_RA_rate = 1 / (1 + min_slope**2)**0.5
        min_normalized_DEC_rate = min_slope / (1 + min_slope**2)**0.5

        if max_normalized_RA_rate * self.RA_vel < 0:
            max_normalized_RA_rate = -max_normalized_RA_rate
            min_normalized_RA_rate = -min_normalized_RA_rate
            max_normalized_DEC_rate = -max_normalized_DEC_rate
            min_normalized_DEC_rate = -min_normalized_DEC_rate
        
        p1 = [self.o1.RA + max_normalized_RA_rate * self.max_vel * dt, self.o1.DEC + max_normalized_DEC_rate * self.max_vel * dt]
        p2 = [self.o1.RA + max_normalized_RA_rate * self.min_vel * dt, self.o1.DEC + max_normalized_DEC_rate * self.min_vel * dt]
        p3 = [self.o1.RA + min_normalized_RA_rate * self.max_vel * dt, self.o1.DEC + min_normalized_DEC_rate * self.max_vel * dt]
        p4 = [self.o1.RA + min_normalized_RA_rate * self.min_vel * dt, self.o1.DEC + min_normalized_DEC_rate * self.min_vel * dt]

        if not plot:
            return p_center, p1, p2, p3, p4

        plt.scatter([self.o1.RA, self.o2.RA], [self.o1.DEC, self.o2.DEC], color="b", label="Previous Observations")
            
        plt.scatter([p_center[0]], [[p_center[1]]], color="r", label="Perfect Observation")
        plt.plot([p1[0], p2[0], p4[0], p3[0], p1[0]], [p1[1], p2[1], p4[1], p3[1], p1[1]], color="yellow", label="Prediction Range")
        plt.plot([self.o1.RA, p_center[0]], [self.o1.DEC, p_center[1]], "--")

        if checked_obs:
            plt.scatter([new_obs.RA], [new_obs.DEC], color="g", label="Checked Observation")

        plt.title("Observation Prediction")
        plt.xlabel("RA (deg)")
        plt.ylabel("DEC (deg)")
        plt.legend(loc='upper left')
        plt.grid()
        plt.show()

def performPairAnalysis(knwon_pair, oc, pix, res):
    known_obses = text1.split("\n")
    o1 = Obs(known_pair[0])
    o2 = Obs(known_pair[1])

    pair = ObsPair(o1, o2, pix, res)
    check = pair.check_obs(oc)

def spherical2cartezian(d, RA, DEC):
    x = d * np.cos(DEC) * np.cos(RA)
    y = d * np.cos(DEC) * np.sin(RA)
    z = d * np.sin(DEC)

    return np.array([x, y, z])

def cartezian2spherical(vec):
    d = (vec[0]**2 + vec[1]**2 + vec[2]**2)**0.5
    RA = np.arctan2(vec[1], vec[0])
    DEC = np.arcsin(vec[2] / d)

    RA_deg = np.rad2deg(RA)
    if RA_deg < 0:
        RA_deg += 360

    return d, RA_deg, np.rad2deg(DEC)
    # return d, np.rad2deg(RA), np.rad2deg(DEC)

def sign(x):
    if x >= 0:
        return 1
    return -1

def gravAccel(mp, bodies):
    accel = np.array([0, 0, 0])
    
    for body in bodies:
        dist = np.linalg.norm(body.pos - mp.pos)
        grav_dir = (body.pos - mp.pos) / dist
        grav_mag = body.GM / dist**2
        
        accel = accel + grav_mag * grav_dir

    return accel

def stepYoshida8(mp, bodies, dt):
    # - - - CONSTANTS - - -
    w1 = 0.311790812418427e0
    w2 = -0.155946803821447e1
    w3 = -0.167896928259640e1
    w4 = 0.166335809963315e1
    w5 = -0.106458714789183e1
    w6 = 0.136934946416871e1
    w7 = 0.629030650210433e0
    w0 = 1.65899088454396

    ds = [w7, w6, w5, w4, w3, w2, w1, w0, w1, w2, w3, w4, w5, w6, w7]

    cs = [0.3145153251052165, 0.9991900571895715, 0.15238115813844, 0.29938547587066, -0.007805591481624963,
          -1.619218660405435, -0.6238386128980216, 0.9853908484811935, 0.9853908484811935, -0.6238386128980216,
          -1.619218660405435, -0.007805591481624963, 0.29938547587066, 0.15238115813844, 0.9991900571895715,
          0.3145153251052165]
    # - - -   - - -   - - -

    for i in range(15):
        mp.pos = mp.pos + mp.vel * cs[i] * dt
        accel = gravAccel(mp, bodies)
        mp.vel = mp.vel + accel * ds[i] * dt

    mp.pos = mp.pos + mp.vel * cs[15] * dt

def computeKepler(r, v, mu=1.3271244004193938e11):
    global km_per_AU
    
    r_mag = np.linalg.norm(r)
    v_mag = np.linalg.norm(v)

    h = np.cross(r, v)
    h_mag = np.linalg.norm(h)

    inclination = np.degrees(np.arccos(h[2] / h_mag))

    k = np.array([0, 0, 1])
    n = np.cross(k, h)
    n_mag = np.linalg.norm(n)

    if n_mag != 0:
        omega = np.degrees(np.arccos(n[0] / n_mag))
        if n[1] < 0:
            omega = 360 - omega
    else:
        omega = 0

    e_vec = (1 / mu) * (np.cross(v, h) - mu * r / r_mag)
    eccentricity = np.linalg.norm(e_vec)

    if n_mag != 0:
        if eccentricity != 0:
            arg_periapsis = np.degrees(np.arccos(np.dot(n, e_vec) / (n_mag * eccentricity)))
            if e_vec[2] < 0:
                arg_periapsis = 360 - arg_periapsis
        else:
            arg_periapsis = 0
    else:
        arg_periapsis = 0

    if eccentricity != 0:
        true_anomaly = np.degrees(np.arccos(np.dot(e_vec, r) / (eccentricity * r_mag)))
        if np.dot(r, v) < 0:
            true_anomaly = 360 - true_anomaly
    else:
        true_anomaly = np.degrees(np.arccos(np.dot(r / r_mag, v / v_mag)))

    specific_energy = v_mag**2 / 2 - mu / r_mag
    if abs(eccentricity - 1) > 1e-8:
        semi_major_axis = -mu / (2 * specific_energy)
    else:
        semi_major_axis = np.inf

    if eccentricity < 1:
        E = 2 * np.arctan(np.tan(np.radians(true_anomaly) / 2) * np.sqrt((1 - eccentricity) / (1 + eccentricity)))
        if E < 0:
            E += 2 * np.pi
        mean_anomaly = np.degrees(E - eccentricity * np.sin(E))
        
    elif eccentricity > 1:
        F = 2 * np.arctanh(np.tan(np.radians(true_anomaly) / 2) * np.sqrt((eccentricity - 1) / (eccentricity + 1)))
        mean_anomaly = np.degrees(eccentricity * np.sinh(F) - F)
        
    else:
        mean_anomaly = None

    return {
        "a": semi_major_axis / AU,
        "e": eccentricity,
        "i": inclination,
        "lon_asc": omega,
        "arg_peri": arg_periapsis,
        "true_anomaly": true_anomaly,
        "mean anomaly": mean_anomaly
    }

def propagate(p0, v0, date_init, date_final, mark_date=None, dt=None):
    global AU, day

    # generate bodies
    mp = MP("", p0, v0)
    
    bodies = []
    body_names = ["MERCURY BARYCENTER",
                  "VENUS BARYCENTER",
                  "EARTH BARYCENTER",
                  "MARS BARYCENTER",
                  "JUPITER BARYCENTER",
                  "SATURN BARYCENTER",
                  "URANUS BARYCENTER",
                  "NEPTUNE BARYCENTER",
                  "SUN"]
    
    body_GMs = [2.2031780000000021E+04,
                3.2485859200000006E+05,
                4.0350323550225981E+05,
                4.2828375214000022E+04,
                1.2671276480000021E+08,
                3.7940585200000003E+07,
                5.7945486000000080E+06,
                6.8365271005800236E+06,
                1.3271244004193938e11]

    for i in range(9):
        new_body = MainBody(body_names[i], np.array([0, 0, 0]), body_GMs[i])
        bodies.append(new_body)

    # set time parameters
    time_interval = (date_final - date_init).total_seconds()

    if mark_date:
        md_time_interval = (mark_date - date_init).total_seconds()

    if not dt:
        if mark_date:
            dt = min([md_time_interval / 200, 10 * day])
        else:
            dt = min([time_interval / 200, 10 * day])
    
    N_cycles = int(time_interval // dt) + 1
    date_final_actual = date_init + timedelta(seconds=N_cycles * dt)

    rhos = []
    RAs = []
    DECs = []
    mark_RA = None
    mark_DEC = None

    # numerically propagate orbit
    for cycle in range(N_cycles):
        cycle_date = date_init + timedelta(seconds=cycle * dt)
        cycle_date_str = cycle_date.strftime('%Y-%m-%dT%H:%M:%S')
        t = spice.str2et(cycle_date_str)
        
        for ib, body in enumerate(bodies):
            state, _ = spice.spkezr(body.name, t, 'ECLIPJ2000', 'NONE', 'SOLAR SYSTEM BARYCENTER')
            body.pos = state[:3]

        stepYoshida8(mp, bodies, dt)

        if N_cycles > 5000 and cycle % 1000 == 0:
            percent_done = round(cycle / N_cycles * 100, 2)
            print(f"Propagating: {percent_done}%")

        rho, RA, DEC = getRADEC(cycle_date, mp.pos)
        rhos.append(rho)
        RAs.append(RA)
        DECs.append(DEC)

        if mark_date and abs((mark_date - cycle_date).total_seconds()) < 2 * dt:
            mark_RA = RA
            mark_DEC = DEC

    return mp.pos, mp.vel, date_final_actual, rhos, RAs, DECs, mark_RA, mark_DEC

def getEarthPos(obs_date, coord='ecliptic'):
    t = spice.str2et(obs_date.strftime('%Y-%m-%dT%H:%M:%S'))
    if coord == 'ecliptic':
        earth_state, _ = spice.spkezr("EARTH", t, 'ECLIPJ2000', 'NONE', 'SOLAR SYSTEM BARYCENTER')
    else:
        earth_state, _ = spice.spkezr("EARTH", t, 'J2000', 'NONE', 'SOLAR SYSTEM BARYCENTER')
    earth_pos = earth_state[:3]
    return earth_pos

def getRADEC(obs_date, pos):
    earth_pos = getEarthPos(obs_date, 'equatorial')
    return cartezian2spherical(ecliptic2equatorial(pos) - earth_pos)

def equatorial2ecliptic(eq_pos):
    epsilon = np.deg2rad(23.439281)
    
    rotation_matrix = np.array([
        [1, 0, 0],
        [0, np.cos(epsilon), np.sin(epsilon)],
        [0, -np.sin(epsilon), np.cos(epsilon)]
    ])
    
    ecliptic_coords = np.dot(rotation_matrix, eq_pos)
    
    return ecliptic_coords

def ecliptic2equatorial(ec_pos):
    epsilon = np.deg2rad(23.439281)
    
    rotation_matrix = np.array([
        [1, 0, 0],
        [0, np.cos(epsilon), -np.sin(epsilon)],
        [0, np.sin(epsilon), np.cos(epsilon)]
    ])
    
    equatorial_coords = np.dot(rotation_matrix, ec_pos)
    
    return equatorial_coords

# === === === ORBIT DETERMINATION FUNCTIONS === === ===
def placeRelToEarth(obs_date, R, RA, DEC, coord='ecliptic'):
    t = spice.str2et(obs_date.strftime('%Y-%m-%dT%H:%M:%S'))
    if coord == 'ecliptic':
        earth_state, _ = spice.spkezr("EARTH", t, 'ECLIPJ2000', 'NONE', 'SOLAR SYSTEM BARYCENTER')
    else:
        earth_state, _ = spice.spkezr("EARTH", t, 'J2000', 'NONE', 'SOLAR SYSTEM BARYCENTER')
    earth_pos = earth_state[:3]
    p = earth_pos + spherical2cartezian(R, np.deg2rad(RA), np.deg2rad(DEC))
    return p

def constructUnitVector(RA, DEC):
    x = np.cos(np.deg2rad(DEC)) * np.cos(np.deg2rad(RA))
    y = np.cos(np.deg2rad(DEC)) * np.sin(np.deg2rad(RA))
    z = np.sin(np.deg2rad(DEC))

    return np.array([x, y, z]) / np.linalg.norm(np.array([x, y, z]))

def get_rho2s(o1, o2, o3):
    mu = 1.3271244004193938e11

    rhat_1 = constructUnitVector(o1.RA, o1.DEC)
    rhat_2 = constructUnitVector(o2.RA, o2.DEC)
    rhat_3 = constructUnitVector(o3.RA, o3.DEC)

    R_1 = getEarthPos(o1.date, 'equatorial')
    R_2 = getEarthPos(o2.date, 'equatorial')
    R_3 = getEarthPos(o3.date, 'equatorial')

    t3mt1 = (o3.date - o1.date).total_seconds()
    rhatprime_2 = (rhat_3 - rhat_1) / t3mt1
    rhatprimeprime_2 = (rhat_3 - 2 * rhat_2 + rhat_1) / (t3mt1**2)

    tau_1 = (o1.date - o2.date).total_seconds()
    tau_3 = (o3.date - o2.date).total_seconds()
    tau = (o3.date - o1.date).total_seconds()

    p_1 = np.cross(rhat_2, rhat_3)
    p_2 = np.cross(rhat_1, rhat_3)
    p_3 = np.cross(rhat_1, rhat_2)

    D_0 = np.dot(rhat_1, p_1)
    D_11 = np.dot(R_1, p_1)
    D_12 = np.dot(R_1, p_2)
    D_13 = np.dot(R_1, p_3)
    D_21 = np.dot(R_2, p_1)
    D_22 = np.dot(R_2, p_2)
    D_23 = np.dot(R_2, p_3)
    D_31 = np.dot(R_3, p_1)
    D_32 = np.dot(R_3, p_2)
    D_33 = np.dot(R_3, p_3)

    A = 1/D_0 * (-D_12 * tau_3 / tau + D_22 + D_32 * tau_1 / tau)
    B = 1 / (6 * D_0) * (D_12 * (tau_3**2 - tau**2) * tau_3 / tau + D_32 * (tau**2 - tau_1**2) * tau_1 / tau)
    E = np.dot(R_2, rhat_2)

    R_2sq = np.dot(R_2, R_2)
    
    a = -(A**2 + 2*A*E + R_2sq)
    b = -2*mu*B*(A+E)
    c = -mu**2 * B**2

    rho_2s = np.roots([1, 0, a, 0, 0, b, 0, 0, c]).real

    return rho_2s

def perfectV0(R1, deltaR, o1, o2):
    R2 = R1 + deltaR

    p0 = equatorial2ecliptic(placeRelToEarth(o1.date, R1, o1.RA, o1.DEC, 'equatorial'))
    pf = equatorial2ecliptic(placeRelToEarth(o2.date, R2, o2.RA, o2.DEC, 'equatorial'))

    date_final = o2.date
    date_init = o1.date

    delta_time = (date_final - date_init).total_seconds()
    v0 = (pf - p0) / delta_time # km s-1

    error = float('Inf')
    tol = 1e-4
    while error > tol:
        p_final, v_final, date_final_actual, _, _, _, _, _ = propagate(p0, v0, date_init, date_final)
        R_final, RA_final, DEC_final = getRADEC(o2.date, pf)

        p_err = pf - p_final
        v0 = v0 + p_err / delta_time

        error = np.linalg.norm(p_err)

    # print(f"Perfected v0 with {error} km of error.")
    return v0

def determineOrbit(obs_all):
    mu = 1.3271244004193938e11
    
    o1 = obs_all[0]
    o3 = obs_all[len(obs_all) - 1]
    o2 = obs_all[int(len(obs_all) / 2)]

    date_init = o1.date
    date_final = o2.date
    date_check = o3.date

    rho_2s = abs(get_rho2s(o1, o2, o3))

    # --- initial observation guess ---
    # initially assuming perfect observation with no errors
    R1 = max(rho_2s) # max. is usually the closest one
    deltaR = 0 * AU
    R2 = R1 + deltaR

    p0 = equatorial2ecliptic(placeRelToEarth(o1.date, R1, o1.RA, o1.DEC, 'equatorial'))
    pf = equatorial2ecliptic(placeRelToEarth(o2.date, R2, o2.RA, o2.DEC, 'equatorial'))

    delta_time = (date_final - date_init).total_seconds()
    v0 = perfectV0(R1, deltaR, o1, o2)
    v0 = (mu / np.linalg.norm(p0))**0.5 * np.array([-p0[1] / np.linalg.norm(p0), p0[0] / np.linalg.norm(p0), 0])

    orbital_elems = computeKepler(p0, v0)
    print("Initial guess:")
    print(orbital_elems)
    # --- --- --- --- --- --- --- --- ---

    adjust_factor = 1
    adjust_pfactor = 1

    good_fit = False
    retry_count = 0
    max_retry = 100
    while (not good_fit) and retry_count <= max_retry:
        err_val = 0
        for idx_o, o in enumerate(obs_all):
            if o.date != date_init:
                p_check, v_check, date_check_actual, rhos, RAs, DECs, _, _ = propagate(p0, v0, date_init, o.date)
                d_prop, RA_prop, DEC_prop = getRADEC(o.date, p_check)

                RA_err = o.RA - RA_prop
                DEC_err = o.DEC - DEC_prop

                err_val += RA_err**2 + DEC_err**2

        print(f"Iter: {retry_count}, errRA: {RA_err}, errDEC: {DEC_err}")
        orbital_elems = computeKepler(p0, v0)
        # print(orbital_elems)

        if abs(RA_err) < 5 * arcsecond and abs(DEC_err) < 5 * arcsecond:
            good_fit = True
        else:
            adjust_vals = [0, 0, 0, 0, 0, 0]
            adjust_vecs = [np.array([0.01, 0, 0]),
                           np.array([-0.01, 0, 0]),
                           np.array([0, 0.01, 0]),
                           np.array([0, -0.01, 0]),
                           np.array([0, 0, 0.01]),
                           np.array([0, 0, -0.01])]

            for i in range(6):
                adjust_vecs[i] *= adjust_factor

            # adjust vel
            for i in range(6):
                v0_1 = v0 + adjust_vecs[i]

                adjust_vals[i] = 0
                for idx_o, o in enumerate(obs_all):
                    if o.date != date_init:
                        p_check_1, v_check_1, _, _, _, _, _, _ = propagate(p0, v0_1, date_init, o.date)
                        _, RA_prop_1, DEC_prop_1 = getRADEC(o.date, p_check_1)

                        RA_err_1 = o.RA - RA_prop_1
                        DEC_err_1 = o.DEC - DEC_prop_1

                        adjust_vals[i] += RA_err_1**2 + DEC_err_1**2

            idx_min = np.argmin(adjust_vals)
            if adjust_vals[idx_min] < err_val:
                v0 = v0 + adjust_vecs[idx_min]
                adjust_factor *= 2
            else:
                adjust_factor *= 0.1

            # adjust pos
            adjust_pvals = [0, 0, 0, 0, 0, 0]
            adjust_pvecs = [np.array([0.01, 0, 0]),
                           np.array([-0.01, 0, 0]),
                           np.array([0, 0.01, 0]),
                           np.array([0, -0.01, 0]),
                           np.array([0, 0, 0.01]),
                           np.array([0, 0, -0.01])]

            for i in range(6):
                adjust_pvecs[i] *= adjust_pfactor
                
            for i in range(6):
                p0_1 = p0 + adjust_pvecs[i]
                adjust_pvals[i] = 0
                
                for idx_o, o in enumerate(obs_all):
                    if o.date != date_init:
                        p_check_1, v_check_1, _, _, _, _, _, _ = propagate(p0_1, v0, date_init, o.date)
                        _, RA_prop_1, DEC_prop_1 = getRADEC(o.date, p_check_1)

                        RA_err_1 = o.RA - RA_prop_1
                        DEC_err_1 = o.DEC - DEC_prop_1

                        adjust_pvals[i] += RA_err_1**2 + DEC_err_1**2

            idx_min = np.argmin(adjust_pvals)
            if adjust_pvals[idx_min] < err_val:
                p0 = p0 + adjust_pvecs[idx_min]
                adjust_pfactor *= 2
            else:
                adjust_pfactor *= 0.1

        retry_count += 1

    pf, vf, date_check_actual, rhos, RAs, DECs, _, _ = propagate(p0, v0, date_init, obs_all[-1].date)
    orbital_elems = computeKepler(p0, v0)
    print("Final fit:")
    print(p0, v0)
    print(orbital_elems)

    return p0, v0, date_init

# === === === === === === === === === === === === === === === ===

def extractStateVector(filename="elements.txt"):
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    for i, line in enumerate(lines):
        if line.startswith("# State vector"):
            pos_line = lines[i+1][3:53]
            vel_line = lines[i+2][3:53]
            position = np.array([float(x) for x in pos_line.split()])
            velocity = np.array([float(x) for x in vel_line.split()])
        
        if line.startswith("Epoch"):
            epoch_str = line[6:20].strip()
            epoch_date = datetime.strptime(epoch_str, '%Y %b %d.%f')
            
    return position, velocity, epoch_date

def classifyObsNearCurve(X, Y, obs_list, tolerance, known_obs, final_date, downsample_factor=None):
    curve = np.column_stack((X, Y))
    
    close_obs = []
    far_obs = []

    known_mag = known_obs.mag
    
    if downsample_factor is not None:
        curve = curve[::downsample_factor]
    
    for obs in obs_list:
        point = np.array([obs.RA, obs.DEC])
        
        # Find the closest point on the curve
        distances = np.linalg.norm(curve - point, axis=1)
        min_distance = np.min(distances)
        
        if (min_distance < tolerance and abs(obs.mag - known_mag) < 2 # distance and magnitude
            and known_obs.date < obs.date < final_date + timedelta(seconds=(final_date - known_obs.date).total_seconds() * 0.2)): # obs dates
            close_obs.append(obs)
        else:
            far_obs.append(obs)
    
    return close_obs, far_obs

def readObsFile(filename='primary.obs'):
    obses = []
    with open(filename, "r") as f:
        lines = f.readlines()
        for line in lines:
            new_o = Obs(line)
            obses.append(new_o)

    return obses

class AstrometryApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Connectonomicon " + str(version))

        current_row = 0
        
        icon_path = "connectonomicon.png"
        if os.path.exists(icon_path):
            self.ico_img = tk.PhotoImage(file=icon_path).subsample(6, 6)
            self.root.iconphoto(False, self.ico_img)
            header_label = tk.Label(root, image=self.ico_img, bg="#212331")
            header_label.grid(row=current_row, column=0, columnspan=2, sticky="n", pady=10)
            current_row += 1
            
        self.root.configure(bg="#212331")

        # Input fields
        tk.Label(root, text="Known Observations File:", bg="#212331", fg="#eeeeee").grid(row=current_row, column=0, padx=10, pady=5, sticky="w")
        self.primary_file_entry = tk.Entry(root, width=40, bg="#212350", fg="#eeeeee")
        self.primary_file_entry.insert(0, "primary.obs")
        self.primary_file_entry.grid(row=current_row, column=1, padx=10, pady=5)
        current_row += 1

        tk.Label(root, text="Potential Observations File:", bg="#212331", fg="#eeeeee").grid(row=current_row, column=0, padx=10, pady=5, sticky="w")
        self.secondary_file_entry = tk.Entry(root, width=40, bg="#212350", fg="#eeeeee")
        self.secondary_file_entry.insert(0, "secondary.obs")
        self.secondary_file_entry.grid(row=current_row, column=1, padx=10, pady=5)
        current_row += 1

        tk.Label(root, text="Propagation Time (days):", bg="#212331", fg="#eeeeee").grid(row=current_row, column=0, padx=10, pady=5, sticky="w")
        self.prop_time_entry = tk.Entry(root, width=40, bg="#212350", fg="#eeeeee")
        self.prop_time_entry.insert(0, '5')
        self.prop_time_entry.grid(row=current_row, column=1, padx=10, pady=5)
        current_row += 1

        tk.Label(root, text="Max. Measurement Error (pixel):", bg="#212331", fg="#eeeeee").grid(row=current_row, column=0, padx=10, pady=5, sticky="w")
        self.pix_error_entry = tk.Entry(root, width=40, bg="#212350", fg="#eeeeee")
        self.pix_error_entry.insert(0, '1')
        self.pix_error_entry.grid(row=current_row, column=1, padx=10, pady=5)
        current_row += 1

        tk.Label(root, text="Pixel Resolution (deg/pixel):", bg="#212331", fg="#eeeeee").grid(row=current_row, column=0, padx=10, pady=5, sticky="w")
        self.resolution_entry = tk.Entry(root, width=40, bg="#212350", fg="#eeeeee")
        self.resolution_entry.insert(0, '9.793873680970562e-05')
        self.resolution_entry.grid(row=current_row, column=1, padx=10, pady=5)
        current_row += 1

        s = ttk.Style() # silly ttk needs styles to style widgets
        s.configure('Connect.TRadiobutton',
                    background='#212331',
                    foreground='#eeeeee')
        
        tk.Label(root, text="Orbit Determination:", bg="#212331", fg="#eeeeee").grid(row=current_row, column=0, padx=10, pady=5, sticky="w")
        self.OD_var = tk.StringVar()
        self.OD_var.set('find_orb')
        self.OD_native = ttk.Radiobutton(
            root,
            text="Native (Experimental)",
            value='native',
            variable=self.OD_var,
            style='Connect.TRadiobutton')
        self.OD_native.grid(row=current_row, column=1, padx=10, pady=5, sticky="w")
        current_row += 1
        self.OD_find_orb = ttk.Radiobutton(
            root,
            text="find_orb (Recommended)",
            value='find_orb',
            variable=self.OD_var,
            style='Connect.TRadiobutton')
        self.OD_find_orb.grid(row=current_row, column=1, padx=10, pady=5, sticky="w")
        current_row += 1

        self.unlikely_var = tk.BooleanVar()
        self.unlikely_checkbox = tk.Checkbutton(
            root, 
            text="Plot Unlikely Observations", 
            variable=self.unlikely_var, 
            bg="#212331", 
            fg="#ffffff",
            selectcolor="#212331",  
            indicatoron=True,     
            highlightbackground="#212331",
            activebackground="#2123aa",
            activeforeground="#eeeeee"
        )
        self.unlikely_checkbox.grid(row=current_row, column=0, columnspan=2, pady=10)
        current_row += 1

        self.run_button = tk.Button(root, text="RUN", command=self.runAnalysis, bg="#0078d7", fg="#fff", width=20)
        self.run_button.grid(row=current_row, column=0, columnspan=2, pady=20)
        current_row += 1

        self.about_button = tk.Button(root, text="About / Help", command=self.showAbout, bg="#aa1100", fg="#eee", width=20)
        self.about_button.grid(row=current_row, column=0, columnspan=2, pady=10)
        current_row += 1

        print(title)
        print("Connectonomicon", version, "initialized.")

    def showAbout(self):
        about_window = tk.Toplevel(self.root)
        about_window.title("About")
        about_window.geometry("900x600")
        about_window.resizable(False, False)

        tk.Label(about_window, text="About Connectonomicon", font=("Arial", 14, "bold")).pack(pady=10)
        tk.Label(
            about_window,
            text=(help_text + "\n" + authors + "\n\n" + license_text),
            wraplength=800,
            justify="left",
        ).pack(pady=10)

        tk.Button(about_window, text="Close", command=about_window.destroy, bg="#0078d7", fg="#fff", width=10).pack(pady=10)

    def runAnalysis(self):
        print("Connectonomicon is running...")
        
        # Get inputs
        primary_file = self.primary_file_entry.get()
        secondary_file = self.secondary_file_entry.get()
        prop_time = float(self.prop_time_entry.get()) * day
        pix_error = float(self.pix_error_entry.get())
        pix_resolution = float(self.resolution_entry.get())

        plot_unlikely_obs = self.unlikely_var.get()
        OD = self.OD_var.get()

        #try:
        # load SPICE kernels
        print("Loading SPICE kernels...")
        spice.furnsh('data/naif0012.tls')
        spice.furnsh('data/de440.bsp')
        
        # orbit determination
        if OD == 'find_orb':
            print("Orbit determination via find_orb...")
            subprocess.run(["fo64", primary_file], timeout=10)

            # extract state vector
            print("Getting state vectors...")
            FO_pos, FO_vel, FO_epoch = extractStateVector()
            p0 = equatorial2ecliptic(FO_pos * km_per_AU)
            v0 = equatorial2ecliptic(FO_vel * km_per_mAU / s_per_d)

            # read observation files
            print("Reading input observations from file...")
            obs_all = readObsFile(primary_file)
            obs_sec = readObsFile(secondary_file)
        else:
            # read observation files
            print("Reading input observations from file...")
            obs_all = readObsFile(primary_file)
            obs_sec = readObsFile(secondary_file)
            
            print("Orbit determination via built-in algorithm...")
            p0, v0, FO_epoch = determineOrbit(obs_all)

        # orbit propagation
        print("Propagating estimated orbit...")
        epoch_compensation = abs(obs_all[0].date - FO_epoch).total_seconds() * 2
        pf, vf, date_final_actual, rhos, RAs, DECs, mark_RA, mark_DEC = propagate(p0, v0, FO_epoch, FO_epoch + timedelta(seconds=prop_time + epoch_compensation), obs_all[-1].date + timedelta(seconds=prop_time))

        # linear prediction
        print("Short arc linear prediction...")
        final_obs_pair = ObsPair(obs_all[-2], obs_all[-1], pix_error, pix_resolution)
        p_center, p1, p2, p3, p4 = final_obs_pair.checkDateObs(obs_all[-1].date + timedelta(seconds=prop_time), False)

        # classify observations
        print("Classifying observations...")
        final_date = obs_all[-1].date + timedelta(seconds=prop_time)
        obs_close, obs_far = classifyObsNearCurve(RAs, DECs, obs_sec, 0.1, obs_all[0], final_date)

        # plot
        print("Plotting...")
        self.plot_results(obs_all, obs_close, obs_far, RAs, DECs, p1, p2, p3, p4, final_obs_pair, p_center, mark_RA, mark_DEC, plot_unlikely_obs)

        print("Done!")

        #except Exception as e:
            #messagebox.showerror("Error", f"An error occurred: {e}")

    def plot_results(self, obs_all, obs_close, obs_far, RAs, DECs, p1, p2, p3, p4, final_obs_pair, p_center, mark_RA, mark_DEC, plot_unlikely_obs):
        plt.plot(RAs, DECs, label="Estimated Orbit", linestyle="--")
        sc1 = plt.scatter([o.RA for o in obs_all], [o.DEC for o in obs_all], label="Known Obs.", color="blue")
        sc = plt.scatter([o.RA for o in obs_close], [o.DEC for o in obs_close], label="Potential Obs.", color="green", marker='+')

        if plot_unlikely_obs:
            sc2 = plt.scatter([o.RA for o in obs_far], [o.DEC for o in obs_far], label="Unlikely Obs.", color="red", marker='*')
            mplcursors.cursor(sc2, hover=False).connect("add", lambda sel: sel.annotation.set_text(str(obs_far[sel.index])))
        
        plt.scatter([mark_RA], [mark_DEC], label="Est. Pos. on Orbit", color="yellow", marker='^')
        plt.plot([p1[0], p2[0], p4[0], p3[0], p1[0]], [p1[1], p2[1], p4[1], p3[1], p1[1]], color="yellow", label="Linear Error Range")
        plt.plot([final_obs_pair.o1.RA, p_center[0]], [final_obs_pair.o1.DEC, p_center[1]], "--", label="Ideal Linear Path")
        plt.grid()
        plt.title("Astrometry Chart")
        plt.legend()
        plt.xlabel("RA (deg)")
        plt.ylabel("DEC (deg)")
        plt.axis('equal')
        mplcursors.cursor(sc, hover=False).connect("add", lambda sel: sel.annotation.set_text(str(obs_close[sel.index])))
        mplcursors.cursor(sc1, hover=False).connect("add", lambda sel: sel.annotation.set_text(str(obs_all[sel.index])))
        plt.show()

if __name__ == "__main__":
    root = tk.Tk()
    app = AstrometryApp(root)
    root.mainloop()
