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
from matplotlib.colors import LinearSegmentedColormap
from scipy import ndimage
from itertools import product


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
        ret.append(field)
    return ret


def get_velocity_field(phases, vel, qxy, size=1, mode='wrap'):
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


def charge_array(Q00, Q01):
    """
    Compute the charge array associated with a Q-tensor field. The defects
    then show up as small regions of non-zero charge (typically 2x2).

    Args:
        Q00, Q01: The components of the nematic field.

    Returns:
        Field of the same shape as Q00 and Q01.
    """
    # compute angle
    def wang(a, b):
        """Infamous chinese function"""
        ang = atan2(abs(a[0]*b[1]-a[1]*b[0]), a[0]*b[0]+a[1]*b[1])
        if ang > pi/2.:
            b = [-i for i in b]
        m = a[0]*b[1]-a[1]*b[0]
        return -np.sign(m)*atan2(abs(m), a[0]*b[0]+a[1]*b[1])

    # get shape and init charge array
    (LX, LY) = Q00.shape
    w = np.zeros((LX, LY))

    # we use the director field instead of Q
    S = np.vectorize(sqrt)(Q00**2 + Q01**2)
    nx = np.vectorize(sqrt)((1 + Q00/S)/2)
    ny = np.sign(Q01)*np.vectorize(sqrt)((1 - Q00/S)/2)

    # This mysterious part was stolen from Amin's code.
    for i in range(LX):
        for j in range(LY):
            ax1 = [nx[(i+1) % LX, j],
                   ny[(i+1) % LX, j]]
            ax2 = [nx[(i-1+LX) % LX, j],
                   ny[(i-1+LX) % LX, j]]
            ax3 = [nx[i, (j-1+LY) % LY],
                   ny[i, (j-1+LY) % LY]]
            ax4 = [nx[i, (j+1) % LY],
                   ny[i, (j+1) % LY]]
            ax5 = [nx[(i+1) % LX, (j-1+LY) % LY],
                   ny[(i+1) % LX, (j-1+LY) % LY]]
            ax6 = [nx[(i-1+LX) % LX, (j-1+LY) % LY],
                   ny[(i-1+LX) % LX, (j-1+LY) % LY]]
            ax7 = [nx[(i+1) % LX, (j+1) % LY],
                   ny[(i+1) % LX, (j+1) % LY]]
            ax8 = [nx[(i-1+LX) % LX, (j+1) % LY],
                   ny[(i-1+LX) % LX, (j+1) % LY]]

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


def get_defects(w, Qxx, Qxy):
    """
    Returns list of defects from charge array.

    Args:
        w: Charge array.

    Returns:
        List of the form [ [ (x, y), charge] ].
    """
    # defects show up as 2x2 regions in the charge array w and must be
    # collapsed to a single point by taking the average position of
    # neighbouring points with the same charge (by breath first search).

    # bfs recursive function
    def collapse(i, j, s, x=0, y=0, n=0):
        if s*w[i, j] > .4:
            x += i + 1.5
            y += j + 1.5
            n += 1
            w[i, j] = 0
            collapse((i+1) % LX, j, s, x, y, n)
            collapse((i-1+LX) % LX, j, s, x, y, n)
            collapse(i, (j+1) % LY, s, x, y, n)
            collapse(i, (j-1+LY) % LY, s, x, y, n)
        return x/n, y/n

    (LX, LY) = w.shape
    d = []

    for i in range(LX):
        for j in range(LY):
            if abs(w[i, j]) > 0.4:
                # charge sign
                s = np.sign(w[i, j])
                # bfs
                x, y = collapse(i, j, s)
                # compute angle, see doi:10.1039/c6sm01146b
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
                # add defect to list
                d.append({"pos": np.array([x, y]),
                          "charge": .5*s,
                          "angle": psi})
    return d


def defects(Q00, Q01, engine=plt, arrow_len=0):
    """
    Plot defects of the nematic field Q.

    Args:
        Q00, Q01: Components of the nematic field.
        engine: Plotting engine or axis.
        arrow_len: If non-zero plot speed of defects as well.
    """
    w = charge_array(Q00, Q01)
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
        else:
            engine.plot(d["pos"][0], d["pos"][1], 'b^')


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


def solidarea(frame, engine=plt):
    """
    Plot all phase fields with solid colours corresponding to individual areas.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    for i in range(len(frame.phi)):
        color = str(min(1, frame.area[i]/(pi*frame.parameters['R']**2)))
        engine.contourf(np.arange(0, frame.parameters['Size'][0]),
                        np.arange(0, frame.parameters['Size'][1]),
                        frame.phi[i].T,
                        levels=[.5, 10.],
                        colors=color)


def com(frame, engine=plt):
    """
    Plot the center-of-mass of each cell as a red dot. Not really
    photorealistic.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
    """
    for c in frame.com:
        engine.plot(c[0], c[1], 'ro')


def shape(frame, engine=plt):
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
        a = S
        engine.arrow(c[0], c[1],  a*nx,  a*ny, color='k')
        engine.arrow(c[0], c[1], -a*nx, -a*ny, color='k')


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

    # coarse grain
    S = ndimage.generic_filter(S, np.mean, size=avg)
    nx = ndimage.generic_filter(nx, np.mean, size=avg)
    ny = ndimage.generic_filter(ny, np.mean, size=avg)

    (LX, LY) = S.shape

    # construct nematic lines
    x = []
    y = []
    for i, j in product(np.arange(LX, step=avg),
                        np.arange(LY, step=avg)):
        f = avg*(S[i, j] if scale else 1.)
        x.append(i + .5 - f*nx[i, j]/2.)
        x.append(i + .5 + f*nx[i, j]/2.)
        x.append(None)
        y.append(j + .5 - f*ny[i, j]/2.)
        y.append(j + .5 + f*ny[i, j]/2.)
        y.append(None)

    engine.plot(x, y, color='k', linestyle='-', linewidth=1)


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
        defects(Qxx, Qxy, engine=engine, arrow_len=arrow_len)
    # plot
    director(Qxx, Qxy, avg=avg, engine=engine)


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
        defects(Qxx, Qxy, engine=engine, arrow_len=arrow_len)
    # plot
    director(Qxx, Qxy, avg=avg, engine=engine)


def velocity_field(frame, size=15, engine=plt, magn=True, cbar=True, avg=1):
    """
    Plot nematic field associated with the shape tensor of each cell.

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        engine: Plotting engine or axis.
        magn: Plot velocity magnitude as a heatmap?
        cbar: Show color bar?
        avg: Size of the averaging (drops points)
    """
    mode = 'wrap' if frame.parameters['BC'] == 0 else 'constant'
    vx, vy = get_velocity_field(frame.phi, frame.velocity, size, mode=mode)
    vx *= (1.-frame.parameters['walls'])
    vy *= (1.-frame.parameters['walls'])

    if magn:
        m = np.sqrt(vx**2 + vy**2)
        cax = engine.imshow(m.T, interpolation='lanczos', cmap='plasma',
                            origin='lower')
        if cbar:
            plt.colorbar(cax)

    vx = vx.reshape((vx.shape[0]//avg, avg, vx.shape[1]//avg, avg))
    vx = np.mean(vx, axis=(1, 3))
    vy = vy.reshape((vy.shape[0]//avg, avg, vy.shape[1]//avg, avg))
    vy = np.mean(vy, axis=(1, 3))

    cax = engine.quiver(np.arange(0, frame.parameters['Size'][0], step=avg),
                        np.arange(0, frame.parameters['Size'][1], step=avg),
                        vx.T, vy.T,
                        pivot='tail', units='dots', scale_units='dots')


def vorticity_field(frame, size=15, engine=plt, cbar=True):
    """
    Plot nematic field associated with the shape tensor of each cell.

    Args:
        frame: Frame to plot, from archive module.
        size: Coarse-graining size.
        engine: Plotting engine or axis.
        cbar: Show color bar?
    """
    vx, vy = get_velocity_field(frame.phi, frame.velocity, size)
    w = get_vorticity_field(vx, vy)
    cax = engine.imshow(w.T, interpolation='lanczos', cmap='viridis',
                        origin='lower')
    if cbar:
        plt.colorbar(cax)


def _force(frame, i, v, engine=plt, **kwargs):
    """
    Helper function to plot forces.
    """
    c = frame.com[i]
    engine.arrow(c[0], c[1], v[0], v[1], **kwargs)


def velocity(frame, engine=plt, color='r'):
    """
    Plot total velocity of each cell.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        color: Color of the arrow.
    """
    scale = frame.parameters['ninfo']*frame.parameters['nsubsteps']
    for i in range(frame.nphases):
        _force(frame, i,
               scale*frame.velocity[i],
               engine=engine,
               color=color)


def pressure_force(frame, engine=plt, color='b'):
    """
    Plot pressure force of each cell.

    Args:
        frame: Frame to plot, from archive module.
        engine: Plotting engine or axis.
        color: Color of the arrow.
    """
    scale = frame.parameters['ninfo']*frame.parameters['nsubsteps']
    for i in range(frame.nphases):
        _force(frame, i,
               scale*frame.Fpressure[i],
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
        Q00 = frame.S00[i]
        Q01 = frame.S01[i]
        S = sqrt(Q00**2 + Q01**2)
        nx = sqrt((1 + Q00/S)/2)
        ny = np.sign(Q01)*sqrt((1 - Q00/S)/2)
        c = frame.com[i]
        a = frame.parameters['R']/2.5*S
        engine.arrow(c[0], c[1],  a*nx,  a*ny, color='k')
        engine.arrow(c[0], c[1], -a*nx, -a*ny, color='k')


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
