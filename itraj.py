from addlib import *
import numpy as np
from pyqist import pq
defaultQIST = {"path" : "/home/david/wrk/nstgro/qist/data/europa/",
               "nml"  : "europa_config_namelist_2024Aug142100002024Aug20120000.qist"
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
        pq.pw_init_n(self.path+self.nml)
    def state(self,t):
        return pq.pw_state(t)
    def stm(self,t):
        return pq.pw_stm(t)
    def stt(self,t):
        return pq.pw_stt(t)
    def stmInv(self,t):
        return pq.pw_stm_i(t)
    def sttInv(self,t):
        return pq.pw_stt_i(t)
    def prop(self,ta,tb,xa,order=None):
        ordr= 2
        func = {1: pq.pw_prop_once, 2: pq.pw_prop_many}
        r = xa.ndim
        if order is not None : ordr=order
        return func[r](ta,tb,xa, ordr)
    def sttsAToB(self,ta,tb):
        scratchm, scratcht = pq.pw_stts_ab(ta,tb)
        stm = scratchm.T
        stt = np.einsum('kji -> ijk', scratcht)
        return stm,stt
    def sttUpdate(self,ta,tb,xa):
        stm, stt = pq.pw_stt_update(ta,tb,xa)
        return stm,stt
    def tensorChangeBasis(self,RNewOldf, ROldNew0, old_stm, old_stt):
        return pq.pw_tensor_change_of_basis(RNewOldf, ROldNew0, old_stm, old_stt)
    def zMap(self,t,order):
        return pq.pw_zmap(t,order)
