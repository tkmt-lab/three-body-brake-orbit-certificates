"""Exact dyadic helpers reused from own prior R015 audit; no producer imports."""
from fractions import Fraction as Q

class D:
    __slots__=('a','b','e')
    def __init__(self,a,b=None,e=0):
        self.a=int(a);self.b=int(a if b is None else b);self.e=e
        assert self.a<=self.b
    @classmethod
    def read(cls,x):
        x=[{'n':v['n'],'e':v['e'] if int(v['n']) else 0} for v in x]
        e=min(x[0]['e'],x[1]['e'])
        return cls(int(x[0]['n'])<<(x[0]['e']-e),int(x[1]['n'])<<(x[1]['e']-e),e)
    @classmethod
    def pairs(cls,x):
        a,b=(Q(int(v[0]),int(v[1])) for v in x)
        assert a.denominator&(a.denominator-1)==0 and b.denominator&(b.denominator-1)==0
        e=-max(a.denominator.bit_length()-1,b.denominator.bit_length()-1)
        return cls(int(a*(1<<(-e))),int(b*(1<<(-e))),e)
    def aligned(self,e):return self.a<<(self.e-e),self.b<<(self.e-e)
    def __add__(self,o):
        e=min(self.e,o.e);a,b=self.aligned(e);c,d=o.aligned(e)
        return D(a+c,b+d,e)
    def __neg__(self):return D(-self.b,-self.a,self.e)
    def __sub__(self,o):return self+-o
    def __mul__(self,o):
        if self.a==self.b==0 or o.a==o.b==0:return D(0)
        p=(self.a*o.a,self.a*o.b,self.b*o.a,self.b*o.b)
        return D(min(p),max(p),self.e+o.e)
    def square(self):
        return D(0 if self.a<=0<=self.b else min(self.a*self.a,self.b*self.b),max(self.a*self.a,self.b*self.b),2*self.e)
    def rounded(self,bits):
        shift=self.e+bits
        if shift>=0:return D(self.a<<shift,self.b<<shift,-bits)
        return D(self.a>>(-shift),-((-self.b)>>(-shift)),-bits)
    def meet(self,o):
        e=min(self.e,o.e);a,b=self.aligned(e);c,d=o.aligned(e)
        return D(max(a,c),min(b,d),e)
    def contains(self,o):
        e=min(self.e,o.e);a,b=self.aligned(e);c,d=o.aligned(e)
        return a<=c<=d<=b
    def __eq__(self,o):return self.contains(o) and o.contains(self)
    def q(self,n):return Q(n<<self.e) if self.e>=0 else Q(n,1<<(-self.e))
    def endpoints(self):return self.q(self.a),self.q(self.b)
    def pairs_out(self):return [[str(q.numerator),str(q.denominator)] for q in self.endpoints()]

def matrix(raw):return [[D.read(x) for x in row] for row in raw]

def product(a,b,bits=None):
    result=[]
    for row in a:
        line=[]
        for j in range(len(b[0])):
            value=D(0)
            for k,v in enumerate(row):value=value+v*b[k][j]
            line.append(value if bits is None else value.rounded(bits))
        result.append(line)
    return result
