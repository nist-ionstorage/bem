import numpy as np
import matplotlib.pyplot as plt
from scipy import special


def spherical_to_cartesian(r, theta, phi):
    """transform spherical coordinates (r, theta, phi) to cartesian
    coordinates (x, y, z)"""
    cis_theta = np.exp(1j*theta)
    cis_phi = np.exp(1j*phi)
    x = r*cis_theta.imag*cis_phi.real
    y = r*cis_theta.imag*cis_phi.imag
    z = r*cis_theta.real
    return np.broadcast_arrays(x, y, z)


def cartesian_to_spherical(x, y, z):
    """transform cartesian coordinates (x, y, z) to spherical
    coordinates (r, theta, phi)"""
    xy = np.square(x) + np.square(y)
    r = np.sqrt(xy + np.square(z))
    # theta from z-axis down
    # (flip args for xy-plane up)
    theta = np.arctan2(np.sqrt(xy), z)
    phi = np.arctan2(y, x)
    return np.broadcast_arrays(r, theta, phi)


def derivatives_spherical_cartesian(d, r, theta, phi,
        direction="to_spherical"):
    """transform gradient d (3, p, q) at theta (p, q) and phi (p, q)
    between spherical and cartesian coordinates. The components of d are
    (dx, dy, dz) or (dr, dtheta, dphi) respectively.
    """
    cis_theta = np.exp(1j*theta)
    ct, st = cis_theta.real, cis_theta.imag
    cis_phi = np.exp(1j*phi)
    cp, sp = cis_phi.real, cis_phi.imag
    if direction == "to_spherical":
        m = [[cp*st, sp*st, ct],
             [-sp/(r*st), cp/(r*st), np.zeros_like(r)],
             [cp*ct/r, sp*ct/r, -st/r]]
    elif direction == "to_cartesian":
        m = [[cp*st, -r*sp*st, r*cp*ct],
             [sp*st, r*cp*st, r*sp*ct],
             [ct, np.zeros_like(r), -r*st]]
    return [sum(m[i][j]*k for j, k in enumerate(d))
            for i in range(3)]



class SphericalHarmonics(object):
    """Spherical harmonic analysis and synthesis.

    spherical coordinates are a (nmax, 2*nmax) array:

    theta: angle from the z axis
        samples are nmax roots of the Legendre
        polynomials in cos(theta) weighted accordingly
        (Gauss-Legendre quadrature, exact)
    phi: angle from the x axis to the xy projection
        samples are 2*nmax evenly spaced in [0, 2*pi)

    spherical harmonics coefficients are a (nmax, nmax) array.
    the (m, n) component (with 0 <= n < nmax and |m| <= n) is
    located at (       m,            n) for m >=0
               (nmax + m, nmax - 1 - n) for m < 0, like:

    [[( 0, 0), ( 0, 1), ( 0, 2)]
     [(-2, 2), ( 1, 1), ( 1, 2)]
     [(-1, 2), (-1, 1), ( 2, 2)]]

    For m >= 0 (the upper triangle) this is the same order
    as in scipy.special.lpmn().

    The sampled Ymn are normalized such that

    (Ymn * Ym'n' * weight * measure).sum() == delta(m, m') * delta(n, n')

    With Ymn = ptmn * exp(1j * m * phi) where ptmn is the (nmax, nmax,
    nmax) (theta, m, n) array of theta-sampled scaled, associated
    Legendre polynomials.

    The choice of phase and normalization of the Ymn is that of Jackson,
    Magnus, Condon, Shortley etc.
    """

    def __init__(self, nmax):
        self.nmax = nmax
        # precompute legendre roots, weights, theta, phi, ptmn, n, m
        x, w = special.orthogonal.p_roots(nmax)
        self.theta = np.arccos(x)[:, None]
        # self.measure = np.sin(self.theta)
        self.weight = w[:, None]
        self.phi = np.arange(0, 2 * np.pi, np.pi / nmax)[None, :]
        self.n, self.m = np.empty((2, nmax, nmax), np.int)
        for n, m, p, q in self.indices():
            self.n[p, q], self.m[p, q] = n, m
        ptmn = np.array([special.lpmn(nmax - 1, nmax - 1, xi)[0] for xi in x])
        # scaling as in scipy.special.sph_harm
        b = (2 * self.n + 1) / (4 * np.pi)
        b *= special.gamma(self.n - np.absolute(self.m) + 1)
        b /= special.gamma(self.n + np.absolute(self.m) + 1)
        c = (-1) ** np.where(self.m < 0, self.m, 0)
        self.ptmn = np.sqrt(b) * c * ptmn[:, np.absolute(self.m), self.n]
        # prepare weighted and measured ptmn for analysis()
        self.ptmnw = self.ptmn * (self.weight[:, :, None] * np.pi / nmax)

    def spherical_grid(self):
        return np.broadcast_arrays(self.theta, self.phi)

    def cartesian_grid(self, r=1.):
        return spherical_to_cartesian(r, self.theta, self.phi)

    def indices(self):
        """yields n, m spherical harmonic indices
        and corresponding coefficient array indices p, q"""
        for n in range(self.nmax):
            for m in range(-n, 0):
                yield n, m, self.nmax + m, self.nmax - n - 1
            for m in range(n + 1):
                yield n, m, m, n

    def check_orthogonality(self):
        for n1, m1, p1, q1 in self.indices():
            y1 = self.ptmnw[:, p1, q1][:, None] * np.exp(1j * m1 * self.phi) 
            for n2, m2, p2, q2 in self.indices():
                if n2 > n1 or (n2 == n1 and m2 > m1):
                    continue
                y2 = self.ptmn[:, p2, q2][:, None] * np.exp(1j * m2 * self.phi)
                y12 = np.absolute((y1 * y2.conj()).sum())
                if n1 == n2 and m1 == m2:
                    assert abs(y12-1) < 1e-6, ((n1, m1), (n2, m2), y12)
                else:
                    assert y12 < 1e-6, ((n1, m1), (n2, m2), y12)

    def check_identity(self):
        for n, m, p, q in self.indices():
            y = self.ptmn[:, p, q][:, None] * np.exp(1j * m * self.phi)
            ys = np.absolute(y).sum()
            a = self.analysis(y)
            assert np.allclose(np.absolute(a).sum(), 1), (n, m, a)
            assert np.allclose(a[p, q], 1.), (n, m, a)
            y1 = self.synthesis(a)
            y1s = np.absolute(y1).sum()
            assert np.allclose(ys, y1s), (n, m, ys, y1s)
            assert np.allclose(y1, y), (n, m, y, y1)

    def analysis(self, f): # theta, phi
        """Spherical Harmonic Analysis.

        Expand the function f on [self.theta, self.phi]
        in spherical harmonic coefficients a on [self.m, self.n].
        f may be complex or real, a is complex
        This is O(nmax^3) in space and time."""
        f = np.atleast_2d(f)
        assert f.shape == (self.nmax, 2 * self.nmax)
        ft = np.fft.fft(f, axis=1) # theta, m
        # negative, positive frequency lookup works nicely
        # with negative indexing
        a = np.einsum("ijk,ijk->jk", self.ptmnw, ft[:, self.m])
        #a = np.tensordot(self.ptmnw, ft[:, self.m], [[0], [0]])
        #a = (self.ptmnw * ft[:, self.m]).sum(0)
        return a # m, n

    def synthesis(self, a): # m, n
        """Spherical Harmonic Synthesis.

        Expand the spherical harmonic coefficients a on
        [self.m, self.n] in spherical coordinates
        [self.theta, self.phi].
        a may be complex or real, f is complex
        This is O(nmax^3) in time."""
        a = np.atleast_2d(a)
        assert a.shape == (self.nmax, self.nmax)
        # TODO: maybe use np.histogramdd or np.bincount here
        #at = a[None, :, :]*self.ptmn # theta, m, n
        #triui = np.triu_indices(self.nmax, 0)
        #trili = np.tril_indices(self.nmax, -1)
        #pi = np.r_[triui[0], trili[0]]
        #qi = np.r_[triui[1], trili[1]]
        #ft = at[:, pi, qi].sum(2) # theta, m
        ft = np.zeros((self.nmax, 2 * self.nmax), np.complex128)
        for i in range(self.nmax):
            at = self.ptmn[i] * a
            al = np.triu(at, 0).sum(1) # sum over n
            ar = np.tril(at, -1).sum(1) # sum over n
            ft[i] = np.concatenate((al, ar)) # m
        f = np.fft.ifft(ft, axis=1) * self.nmax * 2
        return f # theta, phi


def test_one_mn():
    s = SphericalHarmonics(4)
    s.check_orthogonality()
    s.check_identity()
    n, m = 3, 2
    f = special.sph_harm(m, n, s.phi, s.theta)
    #f = s.ptmn[:, m, n][:, None] * np.exp(1j * m * s.phi)
    a = s.analysis(f)
    f1 = s.synthesis(a)
    a1 = s.analysis(f1)
    fig, ax = plt.subplots(1, 4, figsize=(18, 7),
            subplot_kw=dict(aspect="equal"))
    for axi, yi, xm in zip(ax, (f, a, f1, a1), (2, 1, 2, 1)):
        axi.quiver(yi.real, yi.imag, angles="xy", scale_units="xy", scale=1)
        axi.set_xlim(-1, xm * s.nmax)
        axi.set_ylim(-1, s.nmax)


def test_random():
    s = SphericalHarmonics(16)
    a = np.random.randn(2, s.nmax, s.nmax)
    a = a[0] + 1j * a[1]
    y = s.synthesis(a)
    a1 = s.analysis(y)
    y1 = s.synthesis(a1)
    assert np.allclose(a1, a)
    assert np.allclose(y, y1)
    fig, ax = plt.subplots(1, 4, figsize=(18, 7),
            subplot_kw=dict(aspect="equal"))
    for axi, yi in zip(ax, (y, a, y1, a1)):
        axi.quiver(yi.real, yi.imag, angles="xy", scale_units="xy", scale=2)


def test_gaussian():
    s = SphericalHarmonics(8)
    #y = np.zeros((s.nmax, 2*s.nmax))
    #y[3, 4] = 1.
    y = np.exp(-((s.theta - np.pi / 3) / .5) ** 2
               -((s.phi - np.pi / 3) / .5) ** 2)
    a = s.analysis(y)
    y1 = s.synthesis(a)
    print(np.absolute(y).sum())
    print(np.absolute(y1).sum())
    fig, ax = plt.subplots(1, 4, figsize=(18, 7),
            subplot_kw=dict(aspect="equal"))
    for axi, yi in zip(ax, (y, 3 * a, y1, 20 * (y - y1))):
        axi.quiver(yi.real, yi.imag, angles="xy", scale_units="xy", scale=1)


if __name__ == "__main__":
    test_one_mn()
    test_random()
    test_gaussian()
    plt.show()
