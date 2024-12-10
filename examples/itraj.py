import numpy as np
from pyqist import pq

class Itraj(object):
    def __init__(self, namelist):
        pq.pw_init_n(namelist)
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
    def prop_back(self,tb,ta,xb,order=None):
        ordr= 2
        if order is not None : ordr=order
        return pq.pw_prop_back(tb,ta,xb,ordr)
    def sttsAToB(self,ta,tb):
        scratchm, scratcht = pq.pw_stts_ab(ta,tb)
        stm = scratchm
        stt = scratcht #np.einsum('kji -> ijk', scratcht)
        return stm,stt
    def sttUpdate(self,ta,tb,xa):
        stm, stt = pq.pw_stt_update(ta,tb,xa)
        return stm,stt
    def tensorChangeBasis(self,RNewOldf, ROldNew0, old_stm, old_stt):
        return pq.pw_tensor_change_of_basis(RNewOldf, ROldNew0, old_stm, old_stt)
    def zMap(self,t,order):
        return pq.pw_zmap(t,order)
