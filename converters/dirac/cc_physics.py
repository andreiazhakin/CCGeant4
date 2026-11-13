"""One positive spectrum interpolation shared by rates and all energy moments.

E is the native numeric MeV/u table coordinate. T is electron kinetic energy
in MeV. No isotope/finite-mass transport mapping is performed here.
"""
import numpy as np

T0 = 0.001
ME_MEV = 0.51099895
MU_MEV = 931.494


def tmax(energy):
    return 2*ME_MEV*((1+np.asarray(energy)/MU_MEV)**2-1)


def bracket(grid, query):
    query = np.atleast_1d(np.asarray(query, dtype=float))
    if query.ndim != 1 or np.any(~np.isfinite(query)) or np.any(query < grid[0]) or np.any(query > grid[-1]):
        raise ValueError('Energy outside table coverage')
    hi = np.clip(np.searchsorted(grid, query, side='right'), 1, len(grid)-1)
    lo = hi-1
    w = np.log(query/grid[lo])/np.log(grid[hi]/grid[lo])
    return query, lo, hi, w


def powerlaw_cells(t, density, order):
    """Integral T**order D(T)dT in each log-log linear cell."""
    log_t = np.log(t[1:]/t[:-1])
    v = np.log(density[..., 1:]/density[..., :-1])+(order+1)*log_t
    factor = np.ones_like(v)
    np.divide(np.expm1(v), v, out=factor, where=np.abs(v)>1e-14)
    return density[..., :-1]*t[:-1]**(order+1)*log_t*factor


def integrate_spectrum(t, density):
    result = []
    for order in range(3):
        cells = powerlaw_cells(t, density, order)
        tail = np.zeros_like(density)
        tail[..., :-1] = np.cumsum(cells[..., ::-1], axis=-1)[..., ::-1]
        result.append(tail)
    return result


def delta_moments(data, energy, zp):
    """Return [query,N=0..27,moment=0..2] in cm2*MeV**moment."""
    e, lo, hi, w = bracket(data['energy_MeV_u'], energy)
    idx = data['offsets'][:-1]
    native = np.stack([data[k][:, idx].T for k in
                       ['tail0_cm2', 'tail1_MeV_cm2', 'tail2_MeV2_cm2']], axis=-1)
    lengths = tmax(data['energy_MeV_u'])-T0
    means_u = np.empty_like(native)
    means_u[..., 0] = native[..., 0]
    means_u[..., 1] = (native[..., 1]-T0*native[..., 0])/lengths[:, None]
    means_u[..., 2] = (native[..., 2]-2*T0*native[..., 1]+T0*T0*native[..., 0])/lengths[:, None]**2
    mixed = (1-w[:, None, None])*means_u[lo]+w[:, None, None]*means_u[hi]
    L = (tmax(e)-T0)[:, None]
    out = np.empty((len(e), 28, 3))
    ncalc = native.shape[1]
    out[:, :ncalc, 0] = mixed[..., 0]
    out[:, :ncalc, 1] = T0*mixed[..., 0]+L*mixed[..., 1]
    out[:, :ncalc, 2] = T0*T0*mixed[..., 0]+2*T0*L*mixed[..., 1]+L*L*mixed[..., 2]
    for n in range(ncalc, 28):
        out[:, n] = out[:, ncalc-1]*((zp-n)/(zp-(ncalc-1)))**2
    return out


def node_tail(data, ei, n, cut, order):
    sl = slice(int(data['offsets'][ei]), int(data['offsets'][ei+1]))
    t = data['T_MeV'][sl]
    density = data['D_cm2_MeV'][n, sl]
    tail = data[['tail0_cm2','tail1_MeV_cm2','tail2_MeV2_cm2'][order]][n, sl]
    if cut <= t[0]:
        return float(tail[0])
    if cut >= t[-1]:
        return 0.
    i = int(np.searchsorted(t, cut, side='right')-1)
    power = np.log(density[i+1]/density[i])/np.log(t[i+1]/t[i])
    dc = density[i]*np.exp(power*np.log(cut/t[i]))
    cell = powerlaw_cells(np.array([cut,t[i+1]]), np.array([dc,density[i+1]]), order)[0]
    return float(tail[i+1]+cell)


def delta_tail(data, energy, n, cut, zp):
    """[M0,M1,M2] above cut, with the exact same affine mixture law."""
    if (not isinstance(n, (int, np.integer)) or not 0 <= n <= 27
            or not np.isscalar(cut) or not np.isfinite(cut) or cut < T0):
        raise ValueError('Require N0..27 and electron cut >=1keV')
    e, lo, hi, w = bracket(data['energy_MeV_u'], energy)
    if len(e) != 1:
        raise ValueError('Scalar energy required')
    L = float(tmax(e[0])-T0)
    if cut >= T0+L:
        return np.zeros(3)
    ncalc = data['D_cm2_MeV'].shape[0]
    nn = min(n,ncalc-1)
    scale = ((zp-n)/(zp-nn))**2
    out = np.zeros(3)
    for j, weight in [(int(lo[0]),1-float(w[0])),(int(hi[0]),float(w[0]))]:
        Lj = float(tmax(data['energy_MeV_u'][j])-T0)
        a = L/Lj
        c = T0+(cut-T0)/a
        m0,m1,m2 = [node_tail(data,j,nn,c,k) for k in range(3)]
        b = T0*(1-a)
        out += weight*np.array([m0,a*m1+b*m0,a*a*m2+2*a*b*m1+b*b*m0])
    return scale*out


def interpolate_B(stop, energy):
    result = []
    for ti in range(len(stop['target_Z'])):
        sl = slice(int(stop['offsets'][ti]),int(stop['offsets'][ti+1]))
        grid = stop['energy_MeV_u'][sl]
        query, lo, hi, w = bracket(grid,energy)
        right = stop['B_MeV_cm2'][sl]
        left = stop['B_left_MeV_cm2'][sl]
        values = (1-w)*right[lo]+w*left[hi]
        values[query==grid[-1]] = right[-1]
        result.append(values)
    return np.array(result)


def stopping_at_energy(delta, stop, energy, zp):
    hard = delta_moments(delta,energy,zp)[...,1]
    q2 = (zp-np.arange(28,dtype=float))**2
    return (interpolate_B(stop,energy)[...,None]*q2
            +stop['target_Z'][:,None,None]*hard[None,:,:])
