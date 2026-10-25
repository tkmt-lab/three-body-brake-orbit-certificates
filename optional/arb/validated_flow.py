"""Arb Taylor enclosures with a moving (verified-invertible) QR basis.

State: q1x,q1y,q2x,q2y,p1x,p1y,p2x,p2y; q3=-q1-q2.
All decisions use Arb comparisons. Floating point is used ONLY to choose a basis;
the chosen dyadic matrix and its Arb inverse are used in the enclosure formulas.
"""
from flint import arb, arb_mat, arb_series, ctx
from fractions import Fraction
from collections import Counter
import math


class Reject(Exception):
    pass


def sym(x):
    if not x.is_finite():
        raise Reject('nonfinite symmetric enclosure')
    return arb(0, x.abs_upper().man_exp())


def stack(a, b):
    return arb_mat(a.tolist() + b.tolist())


def eye(n):
    return arb_mat([[int(i == j) for j in range(n)] for i in range(n)])


def norm(a):
    return max(sum(a[i, j].abs_upper() for j in range(a.ncols())).upper() for i in range(a.nrows()))


def norm_data(A, weighted=False):
    """Point diagonal D, and a >= ||D^-1 A D||_inf for A=[[0,I],[H,0]]."""
    if not finite(A):
        raise Reject('nonfinite Jacobian enclosure')
    raw = norm(A)
    if not raw.is_finite():
        raise Reject('nonfinite Jacobian norm bound')
    if not weighted:
        return [arb(1)] * 8, raw
    if weighted == 'hessian':
        raw = max(sum(A[i,j].abs_upper() for j in range(4)).upper() for i in range(4,8))
    if raw == 0:
        return [arb(1)] * 8, arb(1)
    m, e = raw.abs_upper().man_exp()
    # A point power of two at least sqrt(raw); for the legacy mode raw>=1.
    exponent = (int(m).bit_length() + int(e) + 1) // 2
    w = arb(2) ** exponent
    return [arb(1)] * 4 + [w] * 4, max(w, (raw / w).upper())


def weighted_norm(W, d):
    return max(sum(W[i, j].abs_upper() for j in range(W.ncols())).upper() / d[i]
               for i in range(W.nrows())).upper()


def finite(a):
    return all(x.is_finite() for x in a.entries())


def rhs(y, jac=False):
    x, u, z, v = [y[i, 0] for i in range(4)]
    ds = [(z-x, v-u), (-2*x-z, -2*u-v), (-x-2*z, -u-2*v)]
    forces, bs, distances = [], [], []
    for dx, dy in ds:
        r2 = dx*dx + dy*dy
        if not r2 > 0:
            raise Reject('distance box contains collision')
        distances.append(r2.lower())
        g = 1 / (r2 * r2.sqrt())
        forces.append((dx*g, dy*g))
        if jac:
            w = g/r2
            bs.append(arb_mat([[g-3*dx*dx*w, -3*dx*dy*w],
                               [-3*dx*dy*w, g-3*dy*dy*w]]))
    a, b, c = forces
    f = arb_mat([[y[i, 0]] for i in range(4,8)] + [[a[0]+b[0]], [a[1]+b[1]],
                                                              [-a[0]+c[0]], [-a[1]+c[1]]])
    if not jac:
        return f
    blocks = [[-bs[0]-2*bs[1], bs[0]-bs[1]], [bs[0]-bs[2], -bs[0]-2*bs[2]]]
    H = arb_mat([[blocks[i//2][j//2][i%2,j%2] for j in range(4)] for i in range(4)])
    A = arb_mat([[0]*8 for _ in range(8)])
    for i in range(4):
        A[i,4+i] = 1
        for j in range(4):
            A[4+i,j] = H[i,j]
    return f, A, min(distances)


def series_rhs(q, jac):
    x,u,z,v = q
    ds = [(z-x,v-u),(-2*x-z,-2*u-v),(-x-2*z,-u-2*v)]
    fs, bs = [], []
    for dx,dy in ds:
        r2 = dx*dx+dy*dy
        g = (r2.sqrt())**(-3)
        fs.append((dx*g,dy*g))
        if jac:
            w = g/r2
            bs.append([[g-3*dx*dx*w,-3*dx*dy*w],[-3*dx*dy*w,g-3*dy*dy*w]])
    a,b,c = fs
    force = [a[0]+b[0],a[1]+b[1],-a[0]+c[0],-a[1]+c[1]]
    if not jac:
        return force, None
    H = [[None]*4 for _ in range(4)]
    for i in range(2):
        for j in range(2):
            H[i][j] = -bs[0][i][j]-2*bs[1][i][j]
            H[i][2+j] = bs[0][i][j]-bs[1][i][j]
            H[2+i][j] = bs[0][i][j]-bs[2][i][j]
            H[2+i][2+j] = -bs[0][i][j]-2*bs[2][i][j]
    return force,H


def coefficients(y, N, M=None):
    """Normalized derivatives through N for y'=f(y), M'=Df(y) M."""
    ys = [[y[i,0]] for i in range(8)]
    if M is not None:
        mq = [arb_mat(M.tolist()[:4])]
        mp = [arb_mat(M.tolist()[4:])]
        hs = []
    for n in range(N):
        q = [arb_series(ys[i],prec=n+1) for i in range(4)]
        f,H = series_rhs(q,M is not None)
        for i in range(4):
            ys[i].append(ys[4+i][n]/(n+1))
            ys[4+i].append(f[i][n]/(n+1))
        if M is not None:
            hs.append(arb_mat([[H[i][j][n] for j in range(4)] for i in range(4)]))
            s = arb_mat(4,M.ncols())
            for k in range(n+1):
                s += hs[k]*mq[n-k]
            mq.append(mp[n]/(n+1))
            mp.append(s/(n+1))
    yc = [arb_mat([[ys[i][n]] for i in range(8)]) for n in range(N+1)]
    if M is None:
        return yc,None
    return yc,[stack(mq[n],mp[n]) for n in range(N+1)]


def evaluate(cs,h):
    a = cs[-1]
    for c in reversed(cs[:-1]):
        a = a*h+c
    return a


def tube(Y,h,both=False,weighted=False):
    """Picard self-map; h is nonnegative. Enclose the WHOLE time slab."""
    span = arb(0,h) if both else arb(h/2,h/2)
    eps = arb(2)**(-ctx.prec+20)
    f = rhs(Y)
    B = Y + span*f
    B = arb_mat([[arb(v.mid(),(v.rad()*arb('1.25')+eps).upper().man_exp())] for v in B.entries()])
    for _ in range(16):
        f,A,d = rhs(B,True)
        diagonal, a = norm_data(A, weighted)
        if not h*a < arb('0.25'):
            raise Reject('Picard contraction bound too large')
        image = Y+span*f
        if B.contains(image):
            E = h*a/(1-h*a)
            MB = eye(8)
            for i in range(8):
                for j in range(8):
                    MB[i,j] += sym(E * diagonal[i] / diagonal[j])
            return B,MB,d
        # Rebuild a symmetric inflation around Y. This contains both candidates.
        B = arb_mat([[Y[i,0]+sym((image[i,0]-Y[i,0]).abs_upper()*2+eps)] for i in range(8)])
    raise Reject('Picard self-map not obtained')


def qr_basis(A):
    """Numerical proposal only. No orthogonality assumption enters the proof."""
    cols=[]
    for j in range(8):
        v=[float(A[i,j].mid()) for i in range(8)]
        for _ in range(2):
            for q in cols:
                d=sum(x*y for x,y in zip(v,q))
                v=[x-d*y for x,y in zip(v,q)]
        n=math.sqrt(sum(x*x for x in v))
        if not math.isfinite(n) or n<1e-100:
            raise Reject('numerical basis proposal failed')
        cols.append([x/n for x in v])
    Q=arb_mat([[arb(cols[j][i]) for j in range(8)] for i in range(8)])
    Qi=Q.inv()
    if not finite(Qi):
        raise Reject('basis inverse not validated')
    return Q,Qi


def encode(x):
    return {'mid':[str(x.mid().man_exp()[0]),int(x.mid().man_exp()[1])],
            'rad':[str(x.rad().man_exp()[0]),int(x.rad().man_exp()[1])]}


def encode_mat(A):
    return [[encode(A[i,j]) for j in range(A.ncols())] for i in range(A.nrows())]


def fracarb(q):
    return arb(q.numerator)/q.denominator


def fraction(x):
    m,e=x.man_exp()
    return Fraction(int(m))*Fraction(2)**int(e)


def initial(vx,vy,r):
    Y=arb_mat([[-1],[0],[1],[0],[vx],[vy],[vx],[vy]])
    c=Y.mid();Q=eye(8)
    rp=arb_mat([[sym(Y[i,0]-c[i,0])] for i in range(8)])
    rb=arb_mat([[sym(rp[i,0]+(sym(r) if i>=4 else 0))] for i in range(8)])
    D=arb_mat([[0,0],[0,0],[0,0],[0,0],[1,0],[0,1],[1,0],[0,1]])
    return c,Q,rp,rb,D


def integrate(vx,vy,T,r,order=24,tolbits=170,maxstep=Fraction(1,16),
              progress=None,stop_time=None,checkpoint=None,weighted=False,minstep_exp=-28):
    c,Q,rp,rb,R=initial(vx,vy,r)
    end=fraction((T/4).mid())
    if stop_time is not None:
        end=min(end,Fraction(stop_time))
    t=Fraction(0);h0=maxstep;tol=arb(2)**(-tolbits)
    mind=arb(1);steps=[];rejected=0;reasons=Counter()
    while t<end:
        hf=min(h0,end-t);h=fracarb(hf)
        Y=c+Q*rb
        try:
            B,MB,d=tube(Y,h,weighted=weighted)
            yr,mr=coefficients(B,order,MB)
            hp=h**order
            error=max(norm(yr[-1]*hp),norm(mr[-1]*hp))
            if not error < tol:
                raise Reject('Taylor remainder exceeds local tolerance')
            yp,_=coefficients(c,order-1)
            _,ms=coefficients(Y,order-1,eye(8))
            C=evaluate(yp,h)+yr[-1]*hp
            M=evaluate(ms,h)+mr[-1]*hp
            Qn,Qi=qr_basis(M*Q)
            cn=C.mid();delta=C-cn
            Z=Qi*M*Q
            pn=Qi*delta+Z*rp
            bn=Qi*delta+Z*rb
            rn=Z*R
            candidate=[cn,Qn,arb_mat([[sym(x)] for x in pn.entries()]),arb_mat([[sym(x)] for x in bn.entries()]),rn]
            if not all(finite(z) for z in candidate):
                raise Reject('nonfinite endpoint enclosure')
            c,Q,rp,rb,R=candidate
            mind=min(mind,d)
            t+=hf;steps.append([str(hf.numerator),str(hf.denominator)])
            h0=min(maxstep,hf*2) if error < tol/(arb(2)**(order+2)) else hf
            if progress:
                progress(dict(t=float(t),end=float(end),steps=len(steps),rejected=rejected,
                              point_radius=float(norm(rp)),box_radius=float(norm(rb)),
                              min_distance_squared_lower=float(mind),rejection_reasons=dict(reasons)))
            if checkpoint and len(steps)%100==0:
                checkpoint(dict(c=encode_mat(c),Q=encode_mat(Q),rp=encode_mat(rp),rb=encode_mat(rb),
                                R=encode_mat(R),t=[str(t.numerator),str(t.denominator)],steps=steps))
        except Reject as e:
            reasons[str(e)]+=1
            rejected+=1;h0=hf/2
            if h0<Fraction(2)**minstep_exp:
                raise RuntimeError(f'step size underflow at guard 2^{minstep_exp}; last rejection: {e}')
    return dict(c=c,Q=Q,rp=rp,rb=rb,R=R,time=end,steps=steps,
                min_distance_squared_lower=mind,rejected=rejected,rejection_reasons=dict(reasons))


def krawczyk(F,J,r):
    C=J.mid().inv().mid()
    if not finite(C) or C.det().contains(0):
        raise Reject('preconditioner not verified nonsingular')
    Z=eye(3)-C*J
    K=-C*F+Z*arb_mat([[sym(r)] for _ in range(3)])
    ratios=[(K[i,0].abs_upper()/r).upper() for i in range(3)]
    q=norm(Z)
    return dict(passed=all(v<1 for v in ratios) and q<1,
                C=C,K=K,Z=Z,ratios=ratios,contraction=q)


def finish(flow,T,r,weighted=False):
    c,Q,rp,rb,R=[flow[k] for k in ['c','Q','rp','rb','R']]
    dt=T/4-fracarb(flow['time'])
    delta=(dt.abs_upper()+r/4).upper()
    B,MB,dist=tube(c+Q*rb,delta,both=True,weighted=weighted)
    f,A,_=rhs(B,True)
    W=Q*R
    diagonal,a=norm_data(A,weighted);E=delta*a/(1-delta*a)
    wn=weighted_norm(W,diagonal)
    WB=W+arb_mat([[sym(diagonal[i]*E*wn) for j in range(2)] for i in range(8)])
    WE=W+arb(0,delta)*A*WB
    yp=c+Q*rp+dt*f
    ids=[4,5,6]
    F=arb_mat([[yp[i,0]] for i in ids])
    J=arb_mat([[WE[i,0],WE[i,1],f[i,0]/4] for i in ids])
    k=krawczyk(F,J,r)
    g=B[0,0]+2*B[2,0]
    k['nondegenerate']=not g.contains(0)
    k['passed']=k['passed'] and k['nondegenerate']
    k.update(F=F,J=J,endpoint=B,g=g,
             distance_squared_lower=min(flow['min_distance_squared_lower'],dist))
    return k
