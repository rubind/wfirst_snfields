from numpy import * # Yup
import numpy as np
from matplotlib import use
use("PDF")
import matplotlib.pyplot as plt
import ephem
import sfdmap
from astropy.io import ascii
import sys


def add_point(RA, Dec, flipRA = 1, color = 'b', txtlabel = None, symb = '.', label = "", txtang = 0, deg2 = 0):
    assert abs(RA) <= 2.*pi
    assert abs(Dec) <= pi

    if RA <= pi:
        plotRA = RA
    else:
        plotRA = RA - 2.*pi

    if deg2 == 0:
        plt.plot(plotRA*(1. - 2*flipRA), Dec, symb, color = color, label = label)
        plotr = 0
    else:
        assert abs(Dec) < 0.9*pi
        
        plottheta = linspace(0., 2.*pi, 40)
        plotr = sqrt(deg2/pi)*(pi/180.)

        plt.plot(plotRA*(1. - 2*flipRA) + plotr*cos(plottheta)/cos(Dec), Dec + plotr*sin(plottheta), color = color)

    if txtlabel != None:
        plt.text(plotRA*(1. - 2*flipRA) + plotr, Dec, " " + txtlabel, color = color, va = 'top'*int(txtang <= 0)  + 'bottom'*int(txtang > 0), rotation = txtang)

def rad_to_hours(val, flipRA = 1):
    print(val)
    degval = val*(1. - 2*flipRA)*(180./pi)

    return "%02i:00" % around(degval/15. % 24)

def add_bright_stars(threshold = -1):
    star_data = loadtxt("jcat_prime2p_brightest.txt")
    mags = star_data[:,2]
    inds = where(mags < threshold)[0]
    
    for i in range(len(inds)):
        add_point(star_data[inds[i],0]/degrad, star_data[inds[i],1]/degrad, color = 'c', symb = '*', label = "Bright Star, $J < " + str(threshold) + "$" if i == 0 else "")
        print("Bright star!", star_data[inds[i],0], star_data[inds[i],1])

def add_sfd():
    m = sfdmap.SFDMap("sfddata-master") # Must be installed separately

    x1d = linspace(-pi, pi, 400)
    y1d = linspace(-pi/2., pi/2, 200)
    xvals, yvals = meshgrid(x1d, y1d)

    #xvals += random.normal(size = xvals.shape)*0.0001
    #yvals += random.normal(size = xvals.shape)*0.0001

    zvals = xvals*0.
    for i in range(len(x1d)):
        for j in range(len(y1d)):
            zvals[j,i] = m.ebv(x1d[i]*degrad, y1d[j]*degrad)

    #CS = plt.contourf(-xvals, yvals, zvals, levels = [0.01, 0.02, 0.03, 100], colors = [(1, 0, 0), (1, 0.5, 0.5), (1, 0.75, 0.75), (1., 0.9, 0.9)][::-1], alpha = 0.3)
    CS = plt.contourf(-xvals, yvals, zvals, levels = [0, 0.01, 0.02, 0.03], colors = [(0, 0, 1), (0.5, 0.5, 1), (0.75, 0.75, 1), (1, 1, 1)], alpha = 0.3)
    #plt.clabel(CS)

def add_CVZ():
    x1d = linspace(-pi, pi, 100)
    y1d = linspace(-pi/2., pi/2, 50)
    xvals, yvals = meshgrid(x1d, y1d)
    zvals = xvals*0.

    
    
    for i in range(len(x1d)):
        for j in range(len(y1d)):
            eq = ephem.Equatorial(x1d[i], y1d[j], epoch = 2000)
            ec = ephem.Ecliptic(eq, epoch = 2000)
            
            zvals[j,i] = float(ec.lat)*degrad
            print(x1d[i]*degrad, y1d[j]*degrad, zvals[j,i]*degrad)


    CS = plt.contour(-xvals, yvals, zvals, levels = [-ecliptic_limit, 0, ecliptic_limit], colors = 'b', linestyles = "solid")
    plt.clabel(CS)

def add_HLS():
    RAs, Decs = loadtxt("hls_vertices.txt", unpack = True)
    plt.plot(RAs*pi/180., Decs*pi/180., color = 'gray')
    


def add_EP():
    ec = ephem.Ecliptic(0, 90./degrad, epoch = 2000)
    eq = ephem.Equatorial(ec, epoch = 2000)
    add_point(float(eq.ra), float(eq.dec), flipRA = 1, color = 'b', label = "NEP/SEP", symb = 's')

    ec = ephem.Ecliptic(0, -90./degrad, epoch = 2000)
    eq = ephem.Equatorial(ec, epoch = 2000)
    add_point(float(eq.ra), float(eq.dec), flipRA = 1, color = 'b', symb = 's')

def add_fields():
    fields_to_plot = ascii.read("fields_to_plot.txt")

    for i in range(len(fields_to_plot["RA"])):
        eq = ephem.Equatorial(fields_to_plot["RA"][i]/degrad, fields_to_plot["Dec"][i]/degrad, epoch = 2000)
        ec = ephem.Ecliptic(eq, epoch = 2000)
        
        if abs(float(ec.lat)*degrad) >= ecliptic_limit:
            pltcolor = 'g'
        else:
            pltcolor = 'r'
        
        add_point(fields_to_plot["RA"][i]/degrad, fields_to_plot["Dec"][i]/degrad, color = pltcolor, txtlabel = fields_to_plot["Field"][i],
                  txtang = fields_to_plot["PlotAngle"][i], deg2 = fields_to_plot["Deg2"][i])


def add_plane():
    for longi in linspace(0, 360., 100):
        ga = ephem.Galactic(longi/degrad, 0./degrad, epoch = 2000)
        eq = ephem.Equatorial(ga, epoch = 2000)
        add_point(float(eq.ra), float(eq.dec), flipRA = 1, color = 'orange', label = "Galactic Plane" if (longi==0) else "")            



if len(sys.argv) == 1:
    ecliptic_limit = 54
else:
    ecliptic_limit = float(sys.argv[1])
        
fig = plt.figure(figsize=(12, 6), dpi = 288)
#fig = plt.figure(figsize=(18, 9), dpi = 288)
ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)

degrad = 180./pi


add_CVZ()
add_sfd()
add_EP()
add_plane()
#add_HLS()
#add_bright_stars()
add_fields()

plt.xticks(plt.xticks()[0], [rad_to_hours(val) for val in plt.xticks()[0]])

plt.legend(loc = "lower right")
plt.colorbar(label = "MW $E(B-V)$")
plt.savefig("field_considerations_%i.png" % int(np.around(ecliptic_limit)), bbox_inches = "tight", dpi = 288)
plt.close()
