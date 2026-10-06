"""The redshift as the separation of two neighbouring null geodesics.

Light is a bundle of null geodesics; the period is carried by the separation
between two neighbouring rays (two successive crests), which is Lie
transported along the bundle.  This script measures it directly, with no use
of k.U, of Killing charges, or of any expansion:

  1. ray 1 leaves the emitter E at proper time 0 and reaches the observer O;
  2. ray 2 leaves the emitter's worldline a proper time tau_i later and is
     shot (by root finding on its initial direction and affine length) so that
     it lands on the observer's worldline;
  3. tau_f is the proper time along the observer's worldline between the two
     arrivals, i.e. the length of the separation vector between the two rays
     at O, chosen along the observer's worldline.

tau_f / tau_i - 1 is then compared with the closed formula of Proposition 5
(static: alpha_O/alpha_E; ZAMO: (N_O/N_E)(1 - w_E b)/(1 - w_O b)) for all
three metrics, both patches, both observer families.  A central difference in
tau_i removes the O(tau_i) term.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from expansion_check import build, observer, frame


def killing_dir(obs, x, g, xi, iph):
    eph = np.zeros(len(x)); eph[iph] = 1.
    if obs == 'static':
        K = xi.copy()
    else:
        w = -(xi @ g @ eph)/(eph @ g @ eph)
        K = xi + w*eph
    return K, np.sqrt(-(K @ g @ K))


def geodesic(fG, n, x0, k0, L):
    def rhs(l, y):
        x, k = y[:n], y[n:]
        return np.concatenate([k, -np.einsum('abc,b,c->a', fG(x), k, k)])
    s = solve_ivp(rhs, (0, L), np.concatenate([x0, k0]), method='DOP853', rtol=1e-13, atol=1e-15)
    return s.y[:n, -1], s.y[n:, -1]


def run(kind, pars, patch, obs, rO=60., step=1.5, nhat=(0.04, 0.6, 0.8), theta_off=0.01, h=2e-4):
    n, idx, xi, fg, fG = build(kind, pars)
    th0 = np.pi/2 + theta_off if patch == 'equator' else 3*theta_off
    xE = np.zeros(n); xE[idx['r']] = rO; xE[idx['th']] = th0
    gE = fg(xE)
    uE = observer(obs, xE, gE, xi, idx['ph'])
    KE, NE = killing_dir(obs, xE, gE, xi, idx['ph'])
    er, et, ep = frame(uE, gE, idx, n)
    n0 = np.dot(np.array(nhat)/np.linalg.norm(nhat), [er, et, ep])
    # directions to vary the ray: two orthogonal to n0 in the (r,th,ph) span, plus,
    # in geodesic monism, the remaining spatial direction of the (t,u) plane
    basis = [er, et, ep]
    extra = []
    for v in basis:
        w_ = v - (n0 @ gE @ v)*n0
        for f in extra:
            w_ = w_ - (f @ gE @ w_)*f
        if np.sqrt(abs(w_ @ gE @ w_)) > 1e-6:
            extra.append(w_/np.sqrt(w_ @ gE @ w_))
        if len(extra) == 2:
            break
    if kind in ('gm', 'gmrot'):
        e4 = np.zeros(n); e4[1] = 1.          # d_u, made orthogonal to u, e_r, e_th, e_ph
        e4 = e4 + (uE @ gE @ e4)*uE
        for f in (er, et, ep):
            e4 = e4 - (f @ gE @ e4)*f
        extra.append(e4/np.sqrt(e4 @ gE @ e4))
    kE = uE + n0
    xO, kO = geodesic(fG, n, xE, kE, step)
    gO = fg(xO)
    KO, NO = killing_dir(obs, xO, gO, xi, idx['ph'])
    orbit = [idx['ph']] + ([1] if kind in ('gm', 'gmrot') else [])   # coordinates that move along the worldline

    def landing(p, sigma):
        x2 = xE + sigma*KE
        nn = n0 + sum(c*v for c, v in zip(p[:-1], extra))
        g2 = fg(x2)
        u2 = observer(obs, x2, g2, xi, idx['ph'])
        nn = nn + (u2 @ g2 @ nn)*u2          # spatial in the emitter's frame
        nn = nn/np.sqrt(nn @ g2 @ nn)
        x, _ = geodesic(fG, n, x2, u2 + nn, p[-1])
        s = (x[0] - xO[0])/KO[0]
        res = [x[idx['r']] - xO[idx['r']], x[idx['th']] - xO[idx['th']]]
        res += [x[i] - xO[i] - s*KO[i] for i in orbit]
        return np.array(res), s

    ratios = []
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    for sigma in (h, -h):
        p0 = np.zeros(len(extra) + 1); p0[-1] = step
        p = fsolve(lambda q: landing(q, sigma)[0], p0, xtol=1e-14)
        res, s = landing(p, sigma)
        assert np.max(np.abs(res)) < 1e-11, res
        tau_i = sigma*NE
        tau_f = s*NO
        ratios.append(tau_f/tau_i)
    z_rays = 0.5*(ratios[0] + ratios[1]) - 1
    # closed formula (Proposition 5)
    eph = np.zeros(n); eph[idx['ph']] = 1.
    E = -(kE @ gE @ xi); Lphi = kE @ gE @ eph
    b = Lphi/E
    wE = -(xi @ gE @ eph)/(eph @ gE @ eph) if obs == 'ZAMO' else 0.
    wO = -(xi @ gO @ eph)/(eph @ gO @ eph) if obs == 'ZAMO' else 0.
    z_formula = (NO/NE)*(1 - wE*b)/(1 - wO*b) - 1
    return z_rays, z_formula


if __name__ == '__main__':
    cases = [('kerr', dict(M=1., a=0.9)), ('mp5', dict(mu=4., a=0.9)),
             ('gm', dict(c0=-0.5, c1=-0.05, c2=1e-4, c3=0., J=0.8)),
             ('gmrot', dict(c0=-0.5, c1=-0.05, c2=1e-4, J=0.8, K=0.002, B=2e-5, Q=1e-9))]
    worst = 0.
    for kind, pars in cases:
        for patch in ('equator', 'pole'):
            for obs in ('static', 'ZAMO'):
                zr, zf = run(kind, pars, patch, obs)
                rel = abs(zr - zf)/abs(zf)
                worst = max(worst, rel)
                print(f'{kind:4s} {patch:7s} {obs:6s} two rays: z = {zr:+.9e}   formula: z = {zf:+.9e}   rel. diff {rel:.1e}')
    print('max relative difference', worst)
    assert worst < 1e-4
    print('PASS Lie transport: the separation of two neighbouring null geodesics gives the redshift of Proposition 5')
