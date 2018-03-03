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
from math import sqrt, pi, atan2, floor, atan, cos, sin
from matplotlib.colors import LinearSegmentedColormap
from scipy import ndimage, signal
from itertools import product
from Queue import Queue

def get_velocity_field(phases, vel, size=1):
    """Compute the coarse grained collective velocity field from a collection
    of phase-fields and their velocities"""
    # glue all velocities together
    vx = np.zeros(phases[0].shape)
    vy = np.zeros(phases[0].shape)
    for n in range(len(phases)):
        vx += vel[n][0]*phases[n]
        vy += vel[n][1]*phases[n]

    # coarse grain
    vx = ndimage.filters.uniform_filter(vx, size=size, mode='wrap')
    vy = ndimage.filters.uniform_filter(vy, size=size, mode='wrap')

    return vx, vy

def get_Qtensor(phases, Qxx, Qxy, size=1):
    """Compute the coarse grained tissue nematic field from individual cells"""
    # glue all Q-tensors together
    QQxx = np.zeros(phases[0].shape)
    QQxy = np.zeros(phases[0].shape)
    for n in range(len(phases)):
        QQxx += Qxx[n]*phases[n]
        QQxy += Qxy[n]*phases[n]

    # coarse grain
    QQxx = ndimage.filters.uniform_filter(QQxx, size=size, mode='wrap')
    QQxy = ndimage.filters.uniform_filter(QQxy, size=size, mode='wrap')

    return QQxx, QQxy

def get_vorticity_field(ux, uy):
    """Compute vorticity from velocity field"""
    dyux = np.gradient(ux, axis=1)
    dxuy = np.gradient(uy, axis=0)
    return dxuy - dyux

def get_corr(u):
    """Compute the correlation of a real two dimensional scalar field"""
    # get 2d correlations
    c = np.fft.rfft2(u)
    c = np.fft.irfft2(np.multiply(c, np.conj(c)))
    # go to polar coords
    s = int(sqrt(c.size)/2)
    r = np.zeros(s)
    n = np.zeros(s)
    k = 0
    for (i, j), v in np.ndenumerate(c):
        k = int(sqrt(i**2 + j**2))
        if k>=s:
            continue
        r[k] += v
        n[k] += 1
    r = np.divide(r, n)
    r /= r[0]
    return r

def get_corr2(ux, uy):
    """Compute the correlation of a real two dimensional vector field"""
    # get 2d correlations
    cx = np.fft.rfft2(ux)
    cx = np.fft.irfft2(np.multiply(cx, np.conj(cx)))
    cy = np.fft.rfft2(uy)
    cy = np.fft.irfft2(np.multiply(cy, np.conj(cy)))
    c  = cx + cy
    # go to polar coords
    s = int(sqrt(c.size)/2)
    r = np.zeros(s)
    n = np.zeros(s)
    k = 0
    for (i, j), v in np.ndenumerate(c):
        k = int(sqrt(i**2 + j**2))
        if k>=s:
            continue
        r[k] += v
        n[k] += 1
    r = np.divide(r, n)
    r /= r[0]
    return r

def charge_array(Q00, Q01):
    """Compute the charge array associated with a Q-tensor field. The defects
    show up as small regions of non-zero charge (typically 2x2)."""

    # compute angle
    def wang(a, b):
        """Infamous chinese function"""
        ang = atan2(abs(a[0]*b[1]-a[1]*b[0]), a[0]*b[0]+a[1]*b[1])
        if ang>pi/2.:
            b = [ -i for i in b ]
        m = a[0]*b[1]-a[1]*b[0]
        return -np.sign(m)*atan2(abs(m), a[0]*b[0]+a[1]*b[1])

    # get shape and init charge array
    (LX, LY) = Q00.shape
    w = np.zeros((LX, LY))

    # we use the director field instead of Q
    S  = np.vectorize(sqrt)(Q00**2 + Q01**2)
    nx = np.vectorize(sqrt)((1 + Q00/S)/2)
    ny = np.sign(Q01)*np.vectorize(sqrt)((1 - Q00/S)/2)

    # This mysterious part was stolen from Amin's code.
    for i in range(LX):
        for j in range(LY):
            ax1 = [ nx[(i+1)%LX,j],              ny[(i+1)%LX,j]              ]
            ax2 = [ nx[(i-1+LX)%LX,j],           ny[(i-1+LX)%LX,j]           ]
            ax3 = [ nx[i,(j-1+LY)%LY],           ny[i,(j-1+LY)%LY]           ]
            ax4 = [ nx[i,(j+1)%LY],              ny[i,(j+1)%LY]              ]
            ax5 = [ nx[(i+1)%LX,(j-1+LY)%LY],    ny[(i+1)%LX,(j-1+LY)%LY]    ]
            ax6 = [ nx[(i-1+LX)%LX,(j-1+LY)%LY], ny[(i-1+LX)%LX,(j-1+LY)%LY] ]
            ax7 = [ nx[(i+1)%LX,(j+1)%LY],       ny[(i+1)%LX,(j+1)%LY]       ]
            ax8 = [ nx[(i-1+LX)%LX,(j+1)%LY],    ny[(i-1+LX)%LX,(j+1)%LY]    ]

            w[i,j]  = wang(ax1, ax5)
            w[i,j] += wang(ax5, ax3);
            w[i,j] += wang(ax3, ax6);
            w[i,j] += wang(ax6, ax2);
            w[i,j] += wang(ax2, ax8);
            w[i,j] += wang(ax8, ax4);
            w[i,j] += wang(ax4, ax7);
            w[i,j] += wang(ax7, ax1);
            w[i,j] /= 2.*pi

    return w

def get_defects(w, Qxx, Qxy):
    """Returns list of defects from charge array.

    Input:
        w   charge array
    Ouptut:
        list of the form [ [ (x, y), charge] ]
    """

    # defects show up as 2x2 regions in the charge array w and must be
    # collapsed to a single point by taking the average position of neighbouring
    # points with the same charge (by breath first search).

    # bfs recursive function
    def collapse(i, j, s, x=0, y=0, n=0):
        if s*w[i,j]>.4:
            x += i + 1.5
            y += j + 1.5
            n += 1
            w[i,j] = 0
            collapse((i+1)%LX, j, s, x, y, n)
            collapse((i-1+LX)%LX, j, s, x, y, n)
            collapse(i, (j+1)%LY, s, x, y, n)
            collapse(i, (j-1+LY)%LY, s, x, y, n)
        return x/n, y/n

    (LX, LY) = w.shape
    d = []

    for i in range(LX):
        for j in range(LY):
            if abs(w[i,j])>0.4:
                # charge sign
                s = np.sign(w[i,j])
                # bfs
                x, y = collapse(i, j, s)
                # compute angle, see doi:10.1039/c6sm01146b
                num = 0
                den = 0
                for (dx, dy) in [ (0,0), (0,1), (1,1), (1,0) ]:
                    # coordinates of nodes around the defect
                    k = (int(x) + LX + dx) % LX
                    l = (int(y) + LY + dy) % LY
                    # derivative at these points
                    dxQxx = .5*(Qxx[(k+1)%LX, l] - Qxx[(k-1+LX)%LX, l])
                    dxQxy = .5*(Qxy[(k+1)%LX, l] - Qxy[(k-1+LX)%LX, l])
                    dyQxx = .5*(Qxx[k, (l+1)%LY] - Qxx[k, (l-1+LY)%LY])
                    dyQxy = .5*(Qxy[k, (l+1)%LY] - Qxy[k, (l-1+LY)%LY])
                    # accumulate numerator and denominator
                    num += s*dxQxy - dyQxx
                    den += dxQxx + s*dyQxy
                psi = s/(2.-s)*atan2(num, den)
                # add defect to list
                d.append({ "pos": np.array([x, y]),
                           "charge": .5*s,
                           "angle": psi
                        })

    return d

class defect_tracker:
    """Defect tracker

    Honestly this is not good code... if you are reading this I feel for you!
    """

    def __init__(self, max_dst = 5):
        # all defects in id order
        self.defects = []
        # list of id of active defects by charge
        self.active  = {}
        # cutoff distance for matching
        self.max_dst = max_dst

    # wrapping distance on pbc domain (ugly)
    def norm(self, v1, v2, LX, LY):
        return min(
                [ np.linalg.norm(v1-v2+[x*LX, y*LY])
                  for x, y in product([-1, 0, 1], [-1, 0, 1]) ]
                )

    def add_frame(self, Q00, Q01, t):
        """Add a new frame and detect defects"""
        # get dimension of pbc domain
        (LX, LY) = Q00.shape
        # get defects from new frame
        w       = charge_array(Q00, Q01)
        defects = get_defects(w, Q00, Q01)
        # new list of active defects
        active  = {}
        # compare with old ones
        for new in defects:
            # distance from all active defects with same charge
            dist = [ self.norm(new["pos"], self.defects[i]["pos"][-1], LX, LY)
                     for i in self.active.get(new['charge'], []) ]
            if len(dist)!=0:
                # find minimum and index
                v, i = min((v, i) for (i, v) in enumerate(dist))
                # if min dist is smaller than cutoff then we have found our defect
                if v < self.max_dst:
                    # add position and angle to stored defect
                    self.defects[self.active[new['charge']][i]]['pos'].append(new['pos'])
                    # register it as active for next round
                    active[new['charge']] = active.get(new['charge'], []) + [
                            self.active[new['charge']][i] ]
                    # delete from current list
                    del self.active[new['charge']][i]
                    continue

            # we have found a new defect
            self.defects.append({
                'charge': new['charge'],
                'birth': t,
                'pos': [ new['pos'] ]
                })
            # register it as active for next round
            active[new['charge']] = active.get(new['charge'], []) + [ len(self.defects)-1 ]

        # save list of acive defects for next run (all the remaining defects are inactive)
        self.active = active

    def plot(self, engine=plt):
        # get all active defects with a certain charge
        for chg, dft in self.active.iteritems():
            # get all indices of these defects
            for ind in dft:
                p = self.defects[ind]['pos'][-1]
                engine.plot(p[0], p[1],
                            'b' if chg==0.5 else 'g',
                            marker=r"$ {} $".format(ind),
                            markersize=15)

def defects(Q00, Q01, engine=plt):
    """Plot single defects of the nematic field Q"""
    w = charge_array(Q00, Q01)
    defects = get_defects(w, Q00, Q01)
    for d in defects:
        if d['charge']==0.5:
            engine.plot(d["pos"][0], d["pos"][1], 'go')
            a = 5 # arrow length
            engine.arrow(d['pos'][0], d['pos'][1],
                         a*cos(d['angle']),  a*sin(d['angle']),
                         color='k', head_width=2, head_length=3)
        else:
            engine.plot(d["pos"][0], d["pos"][1], 'b^')


def cell(frame, i, engine=plt):
    """Plot phase field defining one cells contour"""
    p = frame.phi[i]
    engine.contour(np.arange(0, frame.parameters['Size'][0]),
                   np.arange(0, frame.parameters['Size'][1]),
                   p.T,
                   #levels = [1e-10, 1e-5, .5])
                   levels = [.5],
                   #color='mediumblue'
                   colors='k')


def cells(frame, engine=plt):
    """Plot all phase fields defining the cells contours"""
    for i in range(len(frame.phi)):
        cell(frame, i, engine)

def solid_area(frame, engine=plt):
    """Plot all phase fields with solid colours corresponding to individual areas"""
    for i in range(len(frame.phi)):
        engine.contourf(np.arange(0, frame.parameters['Size'][0]),
                        np.arange(0, frame.parameters['Size'][1]),
                        frame.phi[i].T,
                        #levels = [1e-10, 1e-5, .5])
                        levels = [.5, 10.],
                        colors='mediumblue')
                        #colors=str(frame.area[i]/(np.pi*frame.parameters['R']**2)))

def com(frame, engine=plt):
    """Plot the center-of-mass of each cell as a red dot"""
    for c in frame.com:
        engine.plot(c[0], c[1], 'ro')

def shape(frame, engine=plt):
    """Print shape tensor of each cell as a nematic vector"""

    for i in range(frame.nphases):
        Q00 = frame.S00[i]
        Q01 = frame.S01[i]
        S = sqrt(Q00**2 + Q01**2)
        nx = sqrt((1 + Q00/S)/2)
        ny = np.sign(Q01)*sqrt((1 - Q00/S)/2)
        c = frame.com[i]
        a = S
        engine.arrow(c[0], c[1],  a*nx,  a*ny, color='k')
        engine.arrow(c[0], c[1], -a*nx, -a*ny, color='k')

def director(Qxx, Qxy, avg=1, scale=False, engine=plt):
    """Plot director field associated with nematic tensor with components Q00, Q01"""

    # obtain S, nx, and ny
    S  = np.vectorize(sqrt)(Qxy**2 + Qxx**2)
    nx = np.vectorize(sqrt)((1 + Qxx/S)/2)
    ny = np.sign(Qxy)*np.vectorize(sqrt)((1 - Qxx/S)/2)

    # coarse grain
    S  = ndimage.generic_filter(S , np.mean, size=avg)
    nx = ndimage.generic_filter(nx, np.mean, size=avg)
    ny = ndimage.generic_filter(ny, np.mean, size=avg)

    (LX, LY) = S.shape

    # construct nematic lines
    x = []
    y = []
    for i, j in product(np.arange(LX, step=avg),
                        np.arange(LY, step=avg)):
        f = avg*(S[i,j] if scale else 1.)
        x.append(i + .5 - f*nx[i,j]/2.)
        x.append(i + .5 + f*nx[i,j]/2.)
        x.append(None)
        y.append(j + .5 - f*ny[i,j]/2.)
        y.append(j + .5 + f*ny[i,j]/2.)
        y.append(None)

    engine.plot(x, y, color='k', linestyle='-', linewidth=1)

def nematic_field(frame, size=1, avg=1, show_def=False, engine=plt):
    """Plot nematic field associated with the internal degree of freedom"""

    # get field
    (Qxx, Qxy) = get_Qtensor(frame.phi, frame.Q00, frame.Q01, size=size)
    # plot
    director(Qxx, Qxy, avg=avg, engine=engine)
    # defects
    if show_def: defects(Qxx, Qxy, engine)

def shape_field(frame, size=1, avg=1, show_def=False, engine=plt):
    """Plot nematic field associated with the shape of each cell"""

    # get field
    (Qxx, Qxy) = get_Qtensor(frame.phi, frame.S00, frame.S01, size=size)
    # plot
    director(Qxx, Qxy, avg=avg, engine=engine)
    # defects
    if show_def: defects(Qxx, Qxy, engine)

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

def velf(frame, engine=plt):
    """Print active part of the velocity"""
    for i in range(frame.nphases):
        p = frame.phi[i].reshape(frame.parameters['Size'])
        c = frame.com[i]
        v = frame.velf[i]
        # correction factor
        a = frame.parameters['ninfo']*frame.parameters['nsubsteps']
        engine.arrow(c[0], c[1], a*v[0], a*v[1], color='brown')

def vel(frame, engine=plt):
    """Print active part of the velocity"""
    for i in range(frame.nphases):
        p = frame.phi[i].reshape(frame.parameters['Size'])
        c = frame.com[i]
        v = frame.vel
        a = frame.parameters['ninfo']*frame.parameters['nsubsteps']
        engine.arrow(c[0], c[1], a*v[0], a*v[1], color='k', head_width=1, zorder=10)

def nematic(frame, engine=plt):
    """Print director of each cell"""
    for i in range(frame.nphases):
        Q00 = frame.Q00[i]
        Q01 = frame.Q01[i]
        S = sqrt(Q00**2 + Q01**2)
        nx = sqrt((1 + Q00/S)/2)
        ny = np.sign(Q01)*sqrt((1 - Q00/S)/2)
        c = frame.com[i]
        a = frame.parameters['R'][i]/2.5*S
        #print S
        engine.arrow(c[0], c[1],  a*nx,  a*ny, color='k')
        engine.arrow(c[0], c[1], -a*nx, -a*ny, color='k')

def phase(frame, n, engine=plt):
    """Plot single phase as a density plot"""
    cax = engine.imshow(frame.phi[n].T, interpolation='lanczos', cmap='Greys', origin='lower'
                        #, clim=(0., 1.)
                        )
    cbar = plt.colorbar(cax)

def velocity_field(frame, size=15, engine=plt, magn=True, avg=1):
    """Plot the total veloctity field assiciated with the cells"""
    vx, vy = get_velocity_field(frame.phi, frame.velp + frame.velc + frame.velf, size)

    if magn:
        m = np.sqrt(vx**2 + vy**2)
        cax = engine.imshow(m.T, interpolation='lanczos', cmap='plasma', origin='lower')
        plt.colorbar(cax)

    vx = vx.reshape((vx.shape[0]//avg, avg, vx.shape[1]//avg, avg))
    vx = np.mean(vx, axis=(1,3))
    vy = vy.reshape((vy.shape[0]//avg, avg, vy.shape[1]//avg, avg))
    vy = np.mean(vy, axis=(1,3))

    cax = engine.quiver(np.arange(0, frame.parameters['Size'][0], step=avg),
                        np.arange(0, frame.parameters['Size'][1], step=avg),
                        vx.T, vy.T,
                        pivot='tail', units='dots', scale_units='dots')

def vorticity_field(frame, size=15, engine=plt, cbar=True):
    """Plot the total veloctity field assiciated with the cells"""
    vx, vy = get_velocity_field(frame.phi, frame.velp + frame.velc + frame.velf, size)
    w = get_vorticity_field(vx, vy)
    cax = engine.imshow(w.T, interpolation='lanczos', cmap='viridis', origin='lower')
    if cbar: plt.colorbar(cax)

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
