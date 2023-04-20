from addlib import *
from pyqist import pq
defaultQIST = {
     "t0"       : 0.00000000000000000000000000000000000,
     "tf"       : 1.51111111111111111111282993860588617,
     "filepath" : "/home/david/wrk/nstgro/pycurve/denseqists/",
     "trajfile" : "denseQISTGW_16_0.0.strm",
     }

class Itraj(object):
    def __init__(self, **kwargs):
        """
        kwargs dict must include the following keys with types:
        t0              : real64 (fortran real*8)
        tf              : real64 (fortran real*8)
        filepath        : character
        trajfile        : character
        """
        for k,v in kwargs.items():
            setattr(self,k,v)
        pq.pw_init_i(self.t0,self.tf, self.filepath, self.trajfile)
    def state(self,t):
        return pq.pw_state(t)[:6]
    def stm(self,t):
        return pq.pw_stm(t)[:6,:6]
    def stt(self,t):
        return pq.pw_stt(t)[:6,:6,:6]
    def stmInv(self,t):
        return pq.pw_stm_i(t)[:6,:6]
    def sttInv(self,t):
        return pq.pw_stt_i(t)[:6,:6,:6]
    def prop(self,ta,tb,xa,order=None):
        ordr= 2
        func = {1: pq.pw_prop_once, 2: pq.pw_prop_many}
        r = xa.ndim
        if order is not None : ordr=order
        if tb>self.tf:
            return norbit_push(self,ta,tb,xa, ordr)
        else:
            return func[r](ta,tb,xa, ordr)[:6]
    def sttsAToB(self,ta,tb):
        return pq.pw_stts_ab(ta,tb)
    def zMap(self,t,order):
        return pq.pw_zmap(t,order)
