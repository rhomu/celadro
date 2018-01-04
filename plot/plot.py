# This file is part of CELADRO, Copyright (C) 2016-17, Romain Mueller
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#
# Plotting routines
#

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, pi
from matplotlib.colors import LinearSegmentedColormap
from scipy import ndimage

def get_velocity_field(phases, vel):
    """Compute the collective velocity field from a collection of phase-fields
     and their velocities"""
    v = []
    for k in range(len(phases[0])):
        v = v + [ [ sum([ vel[n][0]*phases[n][k] for n in range(len(phases)) ]),
                    sum([ vel[n][1]*phases[n][k] for n in range(len(phases)) ]) ] ]
    return np.array(v)

def phasefields(frame, engine=plt):
    """Plot all phase fields"""
    for p in frame.phi:
        engine.contour(np.arange(0, frame.parameters['Size'][0]),
                       np.arange(0, frame.parameters['Size'][1]),
                       p.T,
                       #levels = [1e-10, 1e-5, .5])
                       levels = [.5],
                       #color='mediumblue'
                       colors='k')

def interfaces(frame, engine=plt):
    """Plot the interfaces density"""
    totphi = np.zeros(frame.parameters['Size'])
    for i in range(len(frame.phi)):
        totphi += frame.phi[i]*frame.parameters['walls']
        for j in range(i+1, len(frame.phi)):
            totphi += frame.phi[i]*frame.phi[j]

    cmap = LinearSegmentedColormap.from_list('mycmap', ['grey', 'white'])
    engine.imshow(totphi.T, interpolation='lanczos', cmap=cmap, origin='lower')

def solidarea(frame, engine=plt):
    """Plot all phase fields with solid colours corresponding to individual areas"""
    for i in range(len(frame.phi)):
        p = frame.phi[i]
        engine.contourf(np.arange(0, frame.parameters['Size'][0]),
                        np.arange(0, frame.parameters['Size'][1]),
                        frame.phi[i],
                        #levels = [1e-10, 1e-5, .5])
                        levels = [.5, 10.],
                        #color='mediumblue'
                        colors=str(frame.area[i]/(np.pi*frame.parameters['R']**2)))

def com(frame, engine=plt):
    """Plot the center-of-mass of each cell"""
    for c in frame.com:
        engine.plot(c[0], c[1], 'ro')

def ellipses(frame, engine=plt):
    """Plot the shape-ellipses of each cell"""
    for n in range(frame.nphases):
        radius = np.sqrt(frame.area[n]/np.pi/(1-frame.S_order[n]**2))
        print frame.S_order[n], radius
        omega  = frame.S_angle[n]
        p = frame.phi[n].reshape(frame.parameters['Size'])
        c = frame.com[i]
        an = np.linspace(-omega, 2*np.pi-omega, 100)
        engine.plot(c[0] + radius*(1+10*frame.S_order[n])*np.cos(an),
                    c[1] + radius*(1-10*frame.S_order[n])*np.sin(an))

def velc(frame, engine=plt):
    """Print contractile part of the velocity"""
    for i in range(frame.nphases):
        p = frame.phi[i].reshape(frame.parameters['Size'])
        c = frame.com[i]
        v = frame.velc[i]
        # correction factor
        a = frame.parameters['ninfo']*frame.parameters['nsubsteps']
        engine.arrow(c[0], c[1], a*v[0], a*v[1], color='g')


def velp(frame, engine=plt):
    """Print inactive part of the velocity"""
    for i in range(frame.nphases):
        p = frame.phi[i].reshape(frame.parameters['Size'])
        c = frame.com[i]
        v = frame.velp[i]
        # correction factor
        a = frame.parameters['ninfo']*frame.parameters['nsubsteps']
        engine.arrow(c[0], c[1], a*v[0], a*v[1], color='b')

def pol(frame, engine=plt):
    """Print polarisation of each cell"""
    for i in range(frame.nphases):
        p = frame.phi[i].reshape(frame.parameters['Size'])
        c = frame.com[i]
        v = frame.pol[i]
        a = 4#frame.parameters['ninfo']*frame.parameters['nsubsteps']
        engine.arrow(c[0], c[1],  a*v[0],  a*v[1], color='k')
        #engine.arrow(c[0], c[1], -a*v[0], -a*v[1], color='k')

def velf(frame, engine=plt):
    """Print active part of the velocity"""
    for i in range(frame.nphases):
        p = frame.phi[i].reshape(frame.parameters['Size'])
        c = frame.com[i]
        v = frame.velf[i]
        # correction factor
        a = frame.parameters['ninfo']*frame.parameters['nsubsteps']
        engine.arrow(c[0], c[1], a*v[0], a*v[1], color='brown')

def vela(frame, engine=plt):
    """Print active part of the velocity"""
    for i in range(frame.nphases):
        p = frame.phi[i].reshape(frame.parameters['Size'])
        c = frame.com[i]
        v = frame.pol[i]
        a = frame.parameters['alpha']/frame.parameters['xi']*frame.parameters['ninfo']*frame.parameters['nsubsteps']
        engine.arrow(c[0], c[1], a*v[0], a*v[1], color='r')

def vel(frame, engine=plt):
    """Print active part of the velocity"""
    for i in range(frame.nphases):
        p = frame.phi[i].reshape(frame.parameters['Size'])
        c = frame.com[i]
        v = frame.vela[i] + frame.veli[i] + frame.velf[i] + frame.velc[i]
        a = frame.parameters['ninfo']*frame.parameters['nsubsteps']
        engine.arrow(c[0], c[1], a*v[0], a*v[1], color='k', head_width=1, zorder=10)


def phase(frame, n, engine=plt):
    """Plot single phase as a density plot"""
    cax = engine.imshow(frame.phi[n].T, interpolation='lanczos', cmap='Greys', origin='lower'
                        #, clim=(0., 1.)
                        )
    cbar = plt.colorbar(cax)

def velocityfield(frame, size=15, engine=plt):
    """Plot the total veloctity field assiciated with the cells"""
    v = get_velocity_field(frame.phi, frame.parameters['alpha']*frame.pol+frame.velp)
    vx = np.array([ i[0] for i in v ])
    vy = np.array([ i[1] for i in v ])
    vx = ndimage.filters.uniform_filter(vx, size=size, mode='constant')
    vy = ndimage.filters.uniform_filter(vy, size=size, mode='constant')
    cax = engine.quiver(np.arange(0, frame.parameters['Size'][0]),
                        np.arange(0, frame.parameters['Size'][1]),
                        vx.T, vy.T,
                        pivot='tail', units='dots', scale_units='dots')

def walls(frame, engine=plt):
    """Plot the wall phase-field"""
    cax = engine.imshow(frame.parameters['walls'], cmap='Greys', origin='lower', clim=(0., 1.))

def patch(frame, n, engine=plt):
    """Plot the restricted patch of a single cell"""
    plot = lambda m, M: engine.fill([ m[0], M[0], M[0], m[0], m[0], None ],
                                    [ m[1], m[1], M[1], M[1], m[1], None ],
                                    color = 'b', alpha=0.04)
    LX, LY = frame.parameters['Size']
    m = frame.patch_min[n]
    M = frame.patch_max[n]

    if(m[0]==M[0]):
        m[0] += 1e-1
        M[0] -= 1e-1
    if(m[1]==M[1]):
        m[1] += 1e-1
        M[1] -= 1e-1

    if(m[0]>M[0] and m[1]>M[1]):
        plot(m, [ LX, LY ])
        plot([ 0, 0 ], M)
        plot([ m[0], 0 ], [ LX, M[1] ])
        plot([0, m[1] ], [ M[0], LY ])
    elif(m[0]>M[0]):
        plot(m, [ LX, M[1] ])
        plot([ 0, m[1] ], M)
    elif(m[1]>M[1]):
        plot(m, [ M[0], LY ])
        plot([ m[0], 0 ], M)
    else:
        plot(m, M)

def patches(frame, engine=plt):
    """Plot the restricted patches of each cell"""
    for n in range(frame.nphases):
        patch(frame, n, engine)

def masks(frame, engine=plt):
    """Plot division/death masks"""
    m1 = np.array([ 1 if i else 0 for i in frame.division_mask ])
    m2 = np.array([ 1 if i else 0 for i in frame.death_mask ])

    engine.contour(np.arange(0, frame.parameters['LX']),
                   np.arange(0, frame.parameters['LY']),
                   m1.reshape(frame.parameters['Size']).T,
                   levels = [.5], colors = ['b'])

    engine.contour(np.arange(0, frame.parameters['LX']),
                   np.arange(0, frame.parameters['LY']),
                   m2.reshape(frame.parameters['Size']).T,
                   levels = [.5], colors = [ 'r' ])
