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

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, pi, atan2, cos, sin
from matplotlib.colors import LinearSegmentedColormap, to_rgba
from scipy import ndimage
from itertools import product
from scipy.spatial import Voronoi,voronoi_plot_2d


def _get_field(phases, vals, size=1, mode='wrap'):
    """
    Compute the coarse grained field from a collection of phase-fields and
    associated values: ret[i] = sum_j phases[j]*values[i, j].

    Args:
        phases: List of phase-fields.
        vals: List of lists of size (None, len(phases)) of values to be
            associated with each phase-field.
        size: Coarse-graining size.
        mode: How to treat the boundaries, see
            scipy.ndimage.filters.uniform_filter.

    Returns:
        A list of fields, each corresponding to the individual values.
    """
    ret = []
    for vlist in vals:
        assert len(vlist) == len(phases)
        field = np.zeros(phases[0].shape)
        for n in range(len(phases)):
            field += vlist[n]*phases[n]
        field = ndimage.filters.uniform_filter(field, size=size, mode=mode)
        #field = ndimage.filters.gaussian_filter(field, sigma=size, mode=mode)
        ret.append(field)
    return ret

def get_velocity_field(phases, vel, size=1, mode='wrap'):
    """
    Compute coarse-grained nematic field.

    Args:
        phases: List of phase fields.
        vel: List of 2d velocities associated with each phase field.
        size: Coarse-graining size.
        mode: How to treat the boundaries, see
            scipy.ndimage.filters.uniform_filter.
    """
    v0 = [v[0] for v in vel]
    v1 = [v[1] for v in vel]
    return _get_field(phases, [v0, v1], size, mode)

def get_polarity_field(phases, pol, size=1, mode='wrap'):
    """
    Compute coarse-grained nematic field.

    Args:
        phases: List of phase fields.
        vel: List of 2d velocities associated with each phase field.
        size: Coarse-graining size.
        mode: How to treat the boundaries, see
            scipy.ndimage.filters.uniform_filter.
    """
    p0 = [p[0] for p in pol]
    p1 = [p[1] for p in pol]
    return _get_field(phases, [p0, p1], size, mode)


def get_nematic_field(phases, qxx, qxy, size=1, mode='wrap'):
    """
    Compute coarse-grained nematic field.

    Args:
        phases: List of phase fields.
        qxx, qxy: Components of the nematic field of the individual phase
            fields.
        size: Coarse-graining size.
        mode: How to treat the boundaries, see
            scipy.ndimage.filters.uniform_filter.
    """
    return _get_field(phases, [qxx, qxy], size, mode)


def get_vorticity_field(ux, uy, pbc=True):
    """
    Compute vorticity field from velocity field

    Args:
        ux, uy: the individual components of the velocity field.
        pbc: How to treat boundaries, set to true if using pbc.

    Returns:
        Vorticity field.
    """
    if pbc:
        dxuy = np.gradient(np.pad(uy, 2, mode='wrap'), axis=0)[2:-2, 2:-2]
        dyux = np.gradient(np.pad(ux, 2, mode='wrap'), axis=1)[2:-2, 2:-2]
    else:
        dxuy = np.gradient(uy, axis=0)
        dyux = np.gradient(ux, axis=1)
    return dxuy - dyux


def get_gradient_field(ux, uy, pbc=True):
    """
    Compute gradient field from velocity field

    Args:
        ux, uy: the individual components of the velocity field.
        pbc: How to treat boundaries, set to true if using pbc.

    Returns:
        Gradient field.
    """
    if pbc:
        dxux = np.gradient(np.pad(ux, 2, mode='wrap'), axis=0)[2:-2, 2:-2]
        dyuy = np.gradient(np.pad(uy, 2, mode='wrap'), axis=1)[2:-2, 2:-2]
    else:
        dxux = np.gradient(ux, axis=0)
        dyuy = np.gradient(uy, axis=1)
    return dxux + dyuy


def get_corr(u):
    """
    Compute the cross-correlation (as a function of distance) of a real two-
    dimensional scalar field.

    Arguments:
        u: The scalar field.

    Returns:
        The cross-correlation of u as an array.
    """
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
        if k >= s:
            continue
        r[k] += v
        n[k] += 1
    r = np.divide(r, n)
    r /= r[0]
    return r


def get_corr2(ux, uy):
    """
    Compute the correlation (as a function of distance) of two real two-
    dimensional scalar fields.

    Arguments:
        ux, uy: The scalar fields.

    Returns:
        The correlation of ux and uy as an array.
    """
    # get 2d correlations
    cx = np.fft.rfft2(ux)
    cx = np.fft.irfft2(np.multiply(cx, np.conj(cx)))
    cy = np.fft.rfft2(uy)
    cy = np.fft.irfft2(np.multiply(cy, np.conj(cy)))
    c = cx + cy
    # go to polar coords
    s = int(sqrt(c.size)/2)
    r = np.zeros(s)
    n = np.zeros(s)
    k = 0
    for (i, j), v in np.ndenumerate(c):
        k = int(sqrt(i**2 + j**2))
        if k >= s:
            continue
        r[k] += v
        n[k] += 1
    r = np.divide(r, n)
    r /= r[0]
    return r

def nematic_corr(frame,engine = plt,size = 2,show = True):
    # get field
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    (Qxx, Qxy) = get_nematic_field(frame.phi, frame.Q00, frame.Q01,
                                   size=size, mode=mode)
    Qxx *= (1.-frame.parameters['walls'])
    Qxy *= (1.-frame.parameters['walls'])
    (LX,LY) = Qxx.shape
    corr = get_corr2(Qxx, Qxy) 
    radius = int(sqrt(LX*LY)/2)
    if show:
        engine.plot(range(radius), corr,'r')
    return corr

def shape_corr(frame,engine = plt,size = 2,show = True):
    # get field
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    (Qxx, Qxy) = get_nematic_field(frame.phi, frame.S00, frame.S01,
                                   size=size, mode=mode)
    Qxx *= (1.-frame.parameters['walls'])
    Qxy *= (1.-frame.parameters['walls'])
    (LX,LY) = Qxx.shape
    corr = get_corr2(Qxx, Qxy) 
    radius = int(sqrt(LX*LY)/2)
    if show:
        engine.plot(range(radius), corr,'g')
    return corr
 
def velocity_corr(frame,engine = plt,size =1,show = True):
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    vx, vy = get_velocity_field(frame.phi, frame.com_velocity, size, mode=mode)
    vx *= (1.-frame.parameters['walls'])
    vy *= (1.-frame.parameters['walls'])
    corr = get_corr2(vx, vy)
    (LX,LY) = vx.shape
    radius = int(sqrt(LX*LY)/2)
    if show:
        engine.plot(range(radius), corr,'k')
    return corr
    

def vorticity_corr(frame,engine = plt,size = 1,show = True):
    vx, vy = get_velocity_field(frame.phi, frame.com_velocity, size, mode=mode)
    vx *= (1.-frame.parameters['walls'])
    vy *= (1.-frame.parameters['walls'])
    w = get_vorticity_field(vx, vy)
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    w = ndimage.filters.uniform_filter(w, size=size, mode=mode)
    corr = get_corr(w)
    (LX,LY) = vx.shape
    radius = int(sqrt(LX*LY)/2)
    if show:
        engine.plot(range(radius), corr,'r')
    return corr

def _charge_array(vx, vy,Type = 'nematic'):
    """
    Compute the charge array of vector field (vx,vy)

    Args:
        vx: x component of vector field
        vy: y component of vector field
        Type: nematic -- (vx,vy) is head-tail symmetrical
              polar-- (vx,vy) is not head-tail symmetical
    Returns:
        Field of the charge distribution with the same shape as vx and vy
    """
    # compute angle
    def wang(a, b):
        """Infamous chinese function"""
        '''Type = 'nematic' or 'polar' '''
        ang = atan2(abs(a[0]*b[1]-a[1]*b[0]), a[0]*b[0]+a[1]*b[1])	
        if (Type == 'nematic') and (ang > pi/2.):
            b = [-i for i in b]
        m = a[0]*b[1]-a[1]*b[0]
        return -np.sign(m)*atan2(abs(m), a[0]*b[0]+a[1]*b[1])
    (LX, LY) = vx.shape
    w = np.zeros((LX,LY))
    # This mysterious part was stolen from Amin's code.
    # (calculate the winding number)
    for i in range(LX):
        for j in range(LY):
            ax1 = [vx[(i+1) % LX, j],
                   vy[(i+1) % LX, j]]
            ax2 = [vx[(i-1+LX) % LX, j],
                   vy[(i-1+LX) % LX, j]]
            ax3 = [vx[i, (j-1+LY) % LY],
                   vy[i, (j-1+LY) % LY]]
            ax4 = [vx[i, (j+1) % LY],
                   vy[i, (j+1) % LY]]
            ax5 = [vx[(i+1) % LX, (j-1+LY) % LY],
                   vy[(i+1) % LX, (j-1+LY) % LY]]
            ax6 = [vx[(i-1+LX) % LX, (j-1+LY) % LY],
                   vy[(i-1+LX) % LX, (j-1+LY) % LY]]
            ax7 = [vx[(i+1) % LX, (j+1) % LY],
                   vy[(i+1) % LX, (j+1) % LY]]
            ax8 = [vx[(i-1+LX) % LX, (j+1) % LY],
                   vy[(i-1+LX) % LX, (j+1) % LY]]

            w[i, j] = wang(ax1, ax5)
            w[i, j] += wang(ax5, ax3)
            w[i, j] += wang(ax3, ax6)
            w[i, j] += wang(ax6, ax2)
            w[i, j] += wang(ax2, ax8)
            w[i, j] += wang(ax8, ax4)
            w[i, j] += wang(ax4, ax7)
            w[i, j] += wang(ax7, ax1)
            w[i, j] /= 2.*pi

    return w

def Q_charge_array(Q00,Q01):
    '''
    Compute the charge array associated with a Q-tensor field. The defects
    then show up as small regions of non-zero charge (typically 3x3).
    Args:
        Q00: xx component of  Q-tensor field
        Q01: xx component of  Q-tensor field
    Returns:
        Field of the charge distribution with the same shape as Qxx and Qxy
    '''
    (LX, LY) = Q00.shape
    w = np.zeros((LX, LY))
    # get shape and init charge array
    # we use the director field instead of Q
    S = np.vectorize(sqrt)(Q00**2 + Q01**2)
    # here vx and vy represents nematic director
    nx = np.sqrt(2.0*S)*np.vectorize(sqrt)((1 + Q00/S)/2)
    ny = np.sqrt(2.0*S)*np.sign(Q01)*np.vectorize(sqrt)((1 - Q00/S)/2)
    w = _charge_array(nx,ny,Type = 'nematic')
    return w

def polar_charge_array(px,py):
    '''
    Compute the charge array associated with a polarity field. The defects
    then show up as small regions of non-zero charge (typically 3x3).
    Args:
        px: x component of polarisation field
        py: y component of polarisation field
    Returns:
        Field of the charge distribution with the same shape as px and py
    '''
    return _charge_array(px,py,Type = 'polar')
    


def get_defects(w, Qxx = 0, Qxy = 0,cal_angle = False):
    """
    Returns list of defects from charge array.

    Args:
        w: Charge array.
        if Type == 'nematic' then Qxx,Qxy is useful
        if Type == 'polar' then Qxx,Qxy can be omitted

    Returns:
        List of the form [ [ (x, y), charge] ].
    """
    # defects show up as 2x2 regions in the charge array w and must be
    # collapsed to a single point by taking the average position of
    # neighbouring points with the same charge (by breath first search).

    # bfs recursive function
    def collapse(i, j, s, x=0, y=0, n=0,rng = [0.4,0.6]):
        (LX,LY) = w.shape
        if (s*w[i, j] > rng[0]) and (s*w[i, j] < rng[1]):
            w[i, j] = 0
            x1,y1,n1 = collapse((i+1) % LX, j, s, x, y, n,rng)
            x2,y2,n2 = collapse((i-1+LX) % LX, j, s, x, y, n, rng)
            x3,y3,n3 = collapse(i, (j+1) % LY, s, x, y, n, rng)
            x4,y4,n4 = collapse(i, (j-1+LY) % LY, s, x, y, n, rng)
            x = i + x1 + x2 +x3 +x4
            y = j + y1 + y2 +y3 +y4
            n = 1 + n1 + n2 +n3 +n4
            return x,y,n
        else:
            return 0,0,0 

    (LX, LY) = w.shape
    d = []

    for i in range(LX):
        for j in range(LY):
            # detect 
            if  (abs(w[i, j]) > 0.4) and (abs(w[i,j])<0.6):
                # charge sign
                s = np.sign(w[i, j])
                # bfs
                sum_x, sum_y, n = collapse(i, j, s,rng = [0.4,0.6])
                x,y = sum_x/n,sum_y/n
                # compute angle, see doi:10.1039/c6sm01146b
                if cal_angle:
                    num = 0
                    den = 0
                    for (dx, dy) in [(0, 0), (0, 1), (1, 1), (1, 0)]:
                        # coordinates of nodes around the defect
                        kk = (int(x) + LX + dx) % LX
                        ll = (int(y) + LY + dy) % LY
                        # derivative at these points
                        dxQxx = .5*(Qxx[(kk+1) % LX, ll] - Qxx[(kk-1+LX) % LX, ll])
                        dxQxy = .5*(Qxy[(kk+1) % LX, ll] - Qxy[(kk-1+LX) % LX, ll])
                        dyQxx = .5*(Qxx[kk, (ll+1) % LY] - Qxx[kk, (ll-1+LY) % LY])
                        dyQxy = .5*(Qxy[kk, (ll+1) % LY] - Qxy[kk, (ll-1+LY) % LY])
                        # accumulate numerator and denominator
                        num += s*dxQxy - dyQxx
                        den += dxQxx + s*dyQxy
                    psi = s/(2.-s)*atan2(num, den)
                else:
                    psi = 0 

                # add defect to list
                d.append({"pos": np.array([x, y]),
                          "charge": .5*s,
                          "angle": psi})
            elif (abs(w[i, j]) > 0.9) and (abs(w[i,j])<1.1):
		# charge sign
                s = np.sign(w[i, j])
                # bfs
                sum_x, sum_y, n = collapse(i, j, s,rng = [0.9,1.1])
                x,y = sum_x/n,sum_y/n
                # add defect to list
                d.append({"pos": np.array([x, y]),
                          "charge": 1*s,
                          "angle": 0})
		
    return d


def Q_defects(Q00, Q01, engine=plt, arrow_len=0):
    """
    Plot defects of the nematic field Q.

    Args:
        Q00, Q01: Components of the nematic field.
        engine: Plotting engine or axis.
        arrow_len: If non-zero plot speed of defects as well.
    """
    w = Q_charge_array(Q00, Q01)
    defects = get_defects(w, Q00, Q01)
    for d in defects:
        if d['charge'] == 0.5:
            engine.plot(d["pos"][0], d["pos"][1], 'go')
            # plot direction of pos defects
            if not arrow_len == 0:
                engine.arrow(d['pos'][0], d['pos'][1],
                             -arrow_len*cos(d['angle']),
                             -arrow_len*sin(d['angle']),
                             color='r', head_width=3, head_length=3)
        elif d['charge'] == -0.5:
            engine.plot(d["pos"][0], d["pos"][1], 'b^')
        elif d['charge'] == 1:
            engine.plot(d["pos"][0], d["pos"][1], 'r*')
        elif d['charge'] == -1:
            engine.plot(d["pos"][0], d["pos"][1], 'kX')

def polar_defects(px, py, engine=plt, arrow_len=0):
    """
    Plot defects of the polar field (px,py)

    Args:
        px, py: Components of the nematic field.
        engine: Plotting engine or axis.
        arrow_len: If non-zero plot speed of defects as well.
    """
    w = polar_charge_array(px, py)
    defects = get_defects(w)
    for d in defects:
        if d['charge'] == 0.5:
            engine.plot(d["pos"][0], d["pos"][1], 'go')
            # plot direction of pos defects
            if not arrow_len == 0:
                engine.arrow(d['pos'][0], d['pos'][1],
                             -arrow_len*cos(d['angle']),
                             -arrow_len*sin(d['angle']),
                             color='r', head_width=3, head_length=3)
        elif d['charge'] == -0.5:
            engine.plot(d["pos"][0], d["pos"][1], 'b^')
        elif d['charge'] == 1:
            engine.plot(d["pos"][0], d["pos"][1], 'r*')
        elif d['charge'] == -1:
            engine.plot(d["pos"][0], d["pos"][1], 'kX')

  


def cell(frame, i, engine=plt, color='k'):
    """
    Plot a single phase field as a contour.

    Args:
        frame: Frame to plot, from archive module.
        i: Index of the cell to plot.
        engine: Plotting engine or axis.
        color: Color to use for the contour.
    """
    p = frame.phi[i]
    engine.contour(np.arange(0, frame.parameters['Size'][0]),
                   np.arange(0, frame.parameters['Size'][1]),
                   p.T,
                   levels=[.5],
                   colors=color)


def cells(frame, engine=plt, colors='k'):
    """
    Plot all cells as contours.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        colors: Colors to use for the contour. Can also be a list of colors,
            one for each cell.
    """
    if not isinstance(colors, list):
        colors = len(frame.phi)*[colors]

    for i in range(len(frame.phi)):
        cell(frame, i, engine, color=colors[i])


def interfaces(frame, engine=plt):
    """
    Plot the overlap between cells as heatmap using beautiful shades of gray
    for an absolutely photorealistic effect that will impress all your friends.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    totphi = np.zeros(frame.parameters['Size'])
    for i in range(len(frame.phi)):
        totphi += frame.phi[i]*frame.parameters['walls']
        for j in range(i+1, len(frame.phi)):
            totphi += frame.phi[i]*frame.phi[j]

    cmap = LinearSegmentedColormap.from_list('mycmap', ['grey', 'white'])
    engine.imshow(totphi.T, interpolation='lanczos', cmap=cmap, origin='lower')


def interfaces2(frame, engine=plt):
    """
    Plot the overlap between cells as heatmap in a different but also beatiful
    way.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    totphi = [np.zeros(frame.parameters['Size']),
              np.zeros(frame.parameters['Size'])]
    for i in range(len(frame.phi)):
        k = 0 if i < 64 else 1
        totphi[k] += frame.phi[i]*frame.parameters['walls']
        for j in range(i+1, len(frame.phi)):
            totphi[k] += frame.phi[i]*frame.phi[j]

    cmap0 = LinearSegmentedColormap.from_list('mycmap', ['grey', 'white'])
    cmap1 = LinearSegmentedColormap.from_list('mycmap', ['blue', 'white'])
    engine.imshow(totphi[0].T, interpolation='lanczos',
                  cmap=cmap0, origin='lower')
    engine.imshow(totphi[1].T, interpolation='lanczos',
                  cmap=cmap1, origin='lower')


def solidarea(frame, engine=plt, colors='gray'):
    """
    Plot all phase fields with solid colours corresponding to individual areas.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        colors: Colors to use for the contour. Can also be a list of colors,
            one for each cell.
    """
    if not isinstance(colors, list):
        colors = len(frame.phi)*[colors]

    for i in range(len(frame.phi)):
        alpha = 1-min(1, frame.area[i]/(pi*frame.parameters['R'][i]**2))
        engine.contourf(np.arange(0, frame.parameters['Size'][0]),
                        np.arange(0, frame.parameters['Size'][1]),
                        frame.phi[i].T,
                        levels=[.5, 10.],
                        colors=colors[i],
                        alpha=alpha)

def show_neighbour(frame,engine = plt,cell_index = 0):
    for n in frame.nbr_cells[cell_index]:
        engine.text(5+2*n,25,str(n),fontsize = 15)

def com(frame, engine=plt,plotIndex = False):
    """
    Plot the center-of-mass of each cell as a red dot. Not really
    photorealistic.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    for c in frame.com:
        engine.plot(c[0], c[1], 'ro')
    if plotIndex:
        for i in range(len(frame.com)):
            engine.text(frame.com[i][0],frame.com[i][1],str(i),fontsize = 15)
    


def shape(frame, engine=plt, color='k'):
    """
    Print shape tensor of each cell as the director of a nematic tensor.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    for i in range(frame.nphases):
        Q00 = frame.S00[i]
        Q01 = frame.S01[i]
        S = sqrt(Q00**2 + Q01**2)
        w = atan2(Q01, Q00)/2
        nx = cos(w)
        ny = sin(w)
        c = frame.com[i]
        engine.arrow(c[0], c[1],  S*nx,  S*ny, color=color)
        engine.arrow(c[0], c[1], -S*nx, -S*ny, color=color)


def director(Qxx, Qxy, avg=1, scale=False, engine=plt):
    """
    Plot director field associated with a given nematic field.

    Args:
        Qxx, Qxy: Components of the nematic field.
        avg: Coarse-graining size.
        scale: Scale factor that controls the size of the director.
        engine: Plotting engine or axis.
    """

    # obtain S, nx, and ny
    S = np.vectorize(sqrt)(Qxy**2 + Qxx**2)
    nx = np.vectorize(sqrt)((1 + Qxx/S)/2)
    ny = np.sign(Qxy)*np.vectorize(sqrt)((1 - Qxx/S)/2)

    (LX, LY) = S.shape

    # construct nematic lines
    x = []
    y = []
    Smax = np.amax(S)
    for i, j in product(np.arange(LX, step=avg),
                        np.arange(LY, step=avg)):
        f = avg*(S[i, j]/Smax if scale else 1.)
        x.append(i + .5 - f*nx[i, j]/2.)
        x.append(i + .5 + f*nx[i, j]/2.)
        x.append(None)
        y.append(j + .5 - f*ny[i, j]/2.)
        y.append(j + .5 + f*ny[i, j]/2.)
        y.append(None)

    engine.plot(x, y, color='k', linestyle='-', linewidth=1)

def shape_angle(frame, size=1,show_def=False, engine=plt):
    """
    Plot angle of director field associated with a given nematic field.

    Args:
        Qxx, Qxy: Components of the nematic field.
        size: size of uniform filter
        engine: Plotting engine or axis.
        show_def: show defect or not
    """
    # get field
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    (Qxx, Qxy) = get_nematic_field(frame.phi, frame.S00, frame.S01,
                                   size=size, mode=mode)
    Qxx *= (1.-frame.parameters['walls'])
    Qxy *= (1.-frame.parameters['walls'])
    # defects
    if show_def:
        Q_defects(Qxx, Qxy, engine=engine)
    # obtain S, nx, and ny
    S = np.vectorize(sqrt)(Qxy**2 + Qxx**2)
    nx = np.vectorize(sqrt)((1 + Qxx/S)/2)
    ny = np.sign(Qxy)*np.vectorize(sqrt)((1 - Qxx/S)/2)
    # get angle
    angle = np.arctan2(ny,nx) # -pi < angle < pi
    angle[np.where(angle<1e-3)] += np.pi   # make sure angle is between 0,pi
    (LX, LY) = S.shape
    x = np.arange(0,LX)
    y = np.arange(0,LY)
    xv,yv = np.meshgrid(x,y)
    cb = engine.pcolor(angle.T)
    fig = engine.get_figure()
    fig.colorbar(cb, ax=engine)


def nematic_field(frame, size=1, avg=1, show_def=False, arrow_len=0,
                  engine=plt):
    """
    Plot nematic field associated with the internal degree of freedom

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        avg: Average size (reduces the number of points plotted)
        show_def: If true, show defects.
        arrow_len: If non-zero, prints defect speed.
        engine: Plotting engine or axis.
    """
    # get field
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    (Qxx, Qxy) = get_nematic_field(frame.phi, frame.Q00, frame.Q01,
                                   size=size, mode=mode)
    Qxx *= (1.-frame.parameters['walls'])
    Qxy *= (1.-frame.parameters['walls'])
    # defects
    if show_def:
        Q_defects(Qxx, Qxy, engine=engine, arrow_len=arrow_len)
    # plot
    director(Qxx, Qxy, avg=avg, engine=engine,scale = True)


def shape_field(frame, size=1, avg=1, show_def=False, arrow_len=0,
                engine=plt):
    """
    Plot nematic field associated with the shape tensor of each cell.

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        avg: Average size (reduces the number of points plotted)
        show_def: If true, show defects.
        arrow_len: If non-zero, prints defect speed.
        engine: Plotting engine or axis.
    """
    # get field
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    (Qxx, Qxy) = get_nematic_field(frame.phi, frame.S00, frame.S01,
                                   size=size, mode=mode)
    Qxx *= (1.-frame.parameters['walls'])
    Qxy *= (1.-frame.parameters['walls'])
    # defects
    if show_def:
        Q_defects(Qxx, Qxy, engine=engine, arrow_len=arrow_len)
    # plot
    director(Qxx, Qxy, avg=avg, engine=engine,scale = True)

def polarity_field(frame, size=2, show_def = False,engine=plt,magn=True, cbar=True, avg=2,width = 0.2):
    """
    Plot polarity field associated with the internal degree of freedom

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        avg: Average size (reduces the number of points plotted)
        show_def: If true, show defects.
        arrow_len: If non-zero, prints defect speed.
        engine: Plotting engine or axis.
    """
    # get field
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    (px, py) = get_polarity_field(frame.phi, frame.polarization, size=size, mode=mode)
    px *= (1.-frame.parameters['walls'])
    py *= (1.-frame.parameters['walls'])
    # defects
    if show_def:
        polar_defects(px, py, engine=engine)
    # plot
    if magn:
        m = np.sqrt(px**2 + py**2)
        im = engine.imshow(m.T, interpolation='lanczos', cmap='plasma',
                            origin='lower')
        if cbar:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(engine)
            ax_cb = divider.new_horizontal(size="3%", pad=0.03)
            fig = engine.get_figure()
            fig.add_axes(ax_cb)
            plt.colorbar(im, cax=ax_cb)

    px = px.reshape((px.shape[0]//avg, avg, px.shape[1]//avg, avg))
    px = size*np.mean(px, axis=(1, 3))
    py = py.reshape((py.shape[0]//avg, avg, py.shape[1]//avg, avg))
    py = size*np.mean(py, axis=(1, 3))
    
    cax = engine.quiver(np.arange(0, frame.parameters['Size'][0], step=avg),
                        np.arange(0, frame.parameters['Size'][1], step=avg),
                        px.T, py.T,
                        pivot='tail', units='xy',scale_units='xy',width = width)
def velocity_field(frame, size=2, engine=plt,magn=False, cbar=True, avg=2,width = 0.2):
    """
    Plot velocity field associated with the shape tensor of each cell.

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        engine: Plotting engine or axis.
        magn: Plot velocity magnitude as a heatmap?
        cbar: Show color bar?
        avg: Size of the averaging (drops points)
    """

    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    vx, vy = get_velocity_field(frame.phi, frame.com_velocity, size, mode=mode)
    vx *= (1.-frame.parameters['walls'])
    vy *= (1.-frame.parameters['walls'])

    if magn:
        m = np.sqrt(vx**2 + vy**2)
        im = engine.imshow(m.T, interpolation='lanczos', cmap='plasma',
                            origin='lower')
        if cbar:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(engine)
            ax_cb = divider.new_horizontal(size="3%", pad=0.03)
            fig = engine.get_figure()
            fig.add_axes(ax_cb)
            plt.colorbar(im, cax=ax_cb)

    vx = vx.reshape((vx.shape[0]//avg, avg, vx.shape[1]//avg, avg))
    vx = np.mean(vx, axis=(1, 3))
    vy = vy.reshape((vy.shape[0]//avg, avg, vy.shape[1]//avg, avg))
    vy = np.mean(vy, axis=(1, 3))

    cax = engine.quiver(np.arange(0, frame.parameters['Size'][0], step=avg),
                        np.arange(0, frame.parameters['Size'][1], step=avg),
                        vx.T, vy.T,width = width,
                        pivot='tail', units='dots', scale_units='dots')


def __modu(a,b):
    while a >= b:
        a -= b
    while a < 0:
        a += b 
    return a

      
def vorticity_field(frame, size=15, engine=plt, cbar=True):
    """
    Plot nematic field associated with the shape tensor of each cell.

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        engine: Plotting engine or axis.
        cbar: Show color bar?
    """
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    vx, vy = get_velocity_field(frame.phi, frame.com_velocity, size, mode=mode)
    vx *= (1.-frame.parameters['walls'])
    vy *= (1.-frame.parameters['walls'])
    w = get_vorticity_field(vx, vy)
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    w = ndimage.filters.uniform_filter(w, size=size, mode=mode)
    cax = engine.imshow(w.T, interpolation='lanczos', cmap='viridis',
                        origin='lower')

    if cbar:
        plt.colorbar(cax)
    return w
    
def force_density(frame, engine=plt,force_type = 'all',magn=True, cbar=True, avg=1,width = 1,scale = None):
    """
    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        type: all, polar,dipole,passive
        magn: Plot forcedistribution  magnitude as a heatmap?
        cbar: Show color bar?
        avg: Size of the averaging (drops points)
    """
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    # get passive and active force density    
    if force_type == 'all': 
        fx = frame.fpol_field_x + frame.fdipole_field_x + frame.fp_field_x
        fy = frame.fpol_field_y + frame.fdipole_field_y + frame.fp_field_y
    elif force_type == 'polar':
        fx,fy = frame.fpol_field_x,frame.fpol_field_y
    elif force_type == 'dipole': 
        fx,fy = frame.fdipole_field_x,frame.fdipole_field_y
    elif force_type == 'passive':
        fx,fy = frame.fp_field_x,frame.fp_field_y

    fx = fx.reshape((fx.shape[0]//avg, avg, fx.shape[1]//avg, avg))
    fx = np.mean(fx, axis=(1, 3))
    fy = fy.reshape((fy.shape[0]//avg, avg, fy.shape[1]//avg, avg))
    fy = np.mean(fy, axis=(1, 3))

    if magn:
        m = np.sqrt(fx**2 + fy**2)
        im = engine.imshow(m.T, interpolation='lanczos', cmap='plasma',
                            origin='lower')

    if cbar:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(engine)
        ax_cb = divider.new_horizontal(size="3%", pad=0.03)
        fig = engine.get_figure()
        fig.add_axes(ax_cb)
        plt.colorbar(im, cax=ax_cb)

    cax = engine.quiver(np.arange(0, frame.parameters['Size'][0], step=avg),
                        np.arange(0, frame.parameters['Size'][1], step=avg),
                        fx.T, fy.T,width = width,scale = scale,
                        pivot='tail', units='dots', scale_units='dots')

def _force(frame, i, v, engine=plt, **kwargs):
    """
    Helper function to plot forces.
    """
    c = frame.com[i]
    l = sqrt(v[0]*v[0]+v[1]*v[1])
    engine.arrow(c[0], c[1], v[0], v[1],head_width=l/6, head_length=l/3,**kwargs)


def velocity(frame, engine=plt, color='r'):
    """
    Plot total velocity of each cell.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        color: Color of the arrow.
    """
    scale = 0.5*frame.parameters['ninfo']*frame.parameters['nsubsteps']
    for i in range(frame.nphases):
        _force(frame, i,
               1*scale*frame.com_velocity[i],
               engine=engine,
               color=color)

def polar_force(frame, engine=plt, color='b'):
    """
    Plot polar force of each cell.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        color: Color of the arrow.
    """
    scale = 0.2*frame.parameters['ninfo']*frame.parameters['nsubsteps']
    for i in range(frame.nphases):
        _force(frame, i,
               scale*frame.Fpol[i],
               engine=engine,
               color=color)
def shape_force(frame, engine=plt, color='b'):
    """
    Plot shape force of each cell stemming from deformation tensor

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        color: Color of the arrow.
    """
    scale = frame.parameters['ninfo']*frame.parameters['nsubsteps']
    for i in range(frame.nphases):
        _force(frame, i,
               scale*frame.Fshape[i],
               engine=engine,
               color=color)

def passive_force(frame, engine=plt, color='b'):
    """
    Plot passive force of each cell.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        color: Color of the arrow.
    """
    scale = frame.parameters['ninfo']*frame.parameters['nsubsteps']
    for i in range(frame.nphases):
        _force(frame, i,
               scale*frame.Fpassive[i],
               engine=engine,
               color=color)

def interaction_force(frame, engine=plt, color='b'):
    """
    Plot interaction force of each cell.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        color: Color of the arrow.
    """
    scale = frame.parameters['ninfo']*frame.parameters['nsubsteps']*0.01
    for i in range(frame.nphases):
        _force(frame, i,
               scale*frame.Fint[i],
               engine=engine,
               color=color)

def nematic(frame, engine=plt):
    """
    Print director of each cell as a line at their center.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    for i in range(frame.nphases):
        Q00 = frame.Q00[i]
        Q01 = frame.Q01[i]
        S = sqrt(Q00**2 + Q01**2)
        nx = sqrt((1 + Q00/S)/2)
        ny = np.sign(Q01)*sqrt((1 - Q00/S)/2)
        c = frame.com[i]
        a = frame.parameters['R'][i]/2.5*S
        engine.arrow(c[0], c[1],  a*nx,  a*ny, color='k')
        engine.arrow(c[0], c[1], -a*nx, -a*ny, color='k')

def polarization(frame, engine=plt):
    """
    Print poalrization of each cell as a line at their center.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    for i in range(frame.nphases):
        px = frame.polarization[i][0]
        py = frame.polarization[i][1]
        S = sqrt(px**2 + py**2)
        c = frame.com[i]
        engine.arrow(c[0], c[1],  px,  py, head_width=2/3, head_length=0.8,color='g')

def phase(frame, n, engine=plt, cbar=False):
    """
    Plot single phase as a density plot.

    Args:
        frame: Frame to plot, from archive module.
        n: Index of the cell to plot.
        engine: Plotting engine or axis.
        cbar: Display cbar?
    """
    cax = engine.imshow(frame.phi[n].T, interpolation='lanczos', cmap='Greys',
                        origin='lower')
    if cbar:
        plt.colorbar(cax)


def walls(frame, engine=plt, cbar=False):
    """
    Plot walls.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        cbar: Display cbar?
    """
    cax = engine.imshow(frame.parameters['walls'].T, cmap='Greys',
                        origin='lower', clim=(0., 1.))
    if cbar:
        plt.colorbar(cax)


def patch(frame, n, engine=plt):
    """Plot the restricted patch of a single cell

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    def plot(m, M): engine.fill([m[0], M[0], M[0], m[0], m[0], None],
                                [m[1], m[1], M[1], M[1], m[1], None],
                                color='b', alpha=0.04)
    LX, LY = frame.parameters['Size']
    m = frame.patch_min[n]
    M = frame.patch_max[n]

    if(m[0] == M[0]):
        m[0] += 1e-1
        M[0] -= 1e-1
    if(m[1] == M[1]):
        m[1] += 1e-1
        M[1] -= 1e-1

    if(m[0] > M[0] and m[1] > M[1]):
        plot(m, [LX, LY])
        plot([0, 0], M)
        plot([m[0], 0], [LX, M[1]])
        plot([0, m[1]], [M[0], LY])
    elif(m[0] > M[0]):
        plot(m, [LX, M[1]])
        plot([0, m[1]], M)
    elif(m[1] > M[1]):
        plot(m, [M[0], LY])
        plot([m[0], 0], M)
    else:
        plot(m, M)


def patches(frame, engine=plt):
    """
    Plot the subdomain patches of each cell.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    for n in range(frame.nphases):
        patch(frame, n, engine)


def masks(frame, engine=plt):
    """
    Plot division/death masks.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    m1 = np.array([1 if i else 0 for i in frame.division_mask])
    m2 = np.array([1 if i else 0 for i in frame.death_mask])

    engine.contour(np.arange(0, frame.parameters['LX']),
                   np.arange(0, frame.parameters['LY']),
                   m1.reshape(frame.parameters['Size']).T,
                   levels=[.5], colors=['b'])

    engine.contour(np.arange(0, frame.parameters['LX']),
                   np.arange(0, frame.parameters['LY']),
                   m2.reshape(frame.parameters['Size']).T,
                   levels=[.5], colors=['r'])
def trajectories(frame,engine = plt,color = 'm'):
    import animation
    for pos in animation.position: # pos is the postition for each time step until now
        x = pos[:,0]
        y = pos[:,1]
        engine.scatter(x,y,s = 0.3,c = color)
    
def voronoi_lattice(frame,engine = plt):
    """
    plot voronoi lattice polygon
    """
    vor = Voronoi(frame.com)
    voronoi_plot_2d(vor,ax = engine)

