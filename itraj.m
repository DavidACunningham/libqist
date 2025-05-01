from addlib import *
import numpy as np
from pyqist import pq

class Itraj(object):
    """
    A lightweight wrapper class for running QIST models in Python.
    

    Parameters
    ----------
    None

    Examples
    --------
    >>> from itraj import Itraj
    >>> it = Itraj(path="./examples/",nml="curve_deimos_config_namelist_2026Nov26120000002026Nov2617450000.nml")
    """
    def __init__(self, namelist):
        """
        namelist             : character
        """
        pq.pw_init_n(namelist)
    def state(self,t):
        """
        Returns the reference state at a time t.
        For this version of QIST, the state is just the time in seconds
        past J2000, since the reference state is never computed by QIST
        """
        return pq.pw_state(t)
    def stm(self,t):
        """
        Returns the STM at a time t in seconds past J2000
        """
        return pq.pw_stm(t)
    def stt(self,t):
        """
        Returns the STT at a time t in seconds past J2000
        """
        return pq.pw_stt(t)
    def stmInv(self,t):
        """
        Returns the inverse STM at a time t in seconds past J2000
        """
        return pq.pw_stm_i(t)
    def sttInv(self,t):
        """
        Returns the inverse STT at a time t in seconds past J2000
        """
        return pq.pw_stt_i(t)
    def prop(self,ta,tb,xa,order=None):
        """
        Propagates relative states.
        
        Inputs
        ------
        ta: float, time in seconds past J2000 TDB to start propagation
        tb: float, time in seconds past J2000 TDB to end propagation
        xa: float, dimension (8,:), initial relative state in km and km/s with two padding zeros

        Returns the state at a time tb in seconds past J2000 with two padding zeros
        """
        ordr= 2
        func = {1: pq.pw_prop_once, 2: pq.pw_prop_many}
        r = xa.ndim
        if order is not None : ordr=order
        return func[r](ta,tb,xa, ordr)
    def sttsAToB(self,ta,tb):
        """
        Gets stm and stt from a time ta to a time tb
        
        Inputs
        ------
        ta: float, time in seconds past J2000 TDB to start propagation
        tb: float, time in seconds past J2000 TDB to end propagation

        Returns
        ------
        stm: float (8,8), first order STM of the reference trajectory
        stt: float (8,8,8), second order STT of the reference trajectory
        """
        stm, stt = pq.pw_stts_ab(ta,tb)
        return stm,stt
    def sttUpdate(self,ta,tb,xa):
        """
        Returns approximate STM and STT of the chaser with state xa from time ta to time tb
        Inputs
        ------
        ta: float, time in seconds past J2000 TDB to start propagation
        tb: float, time in seconds past J2000 TDB to end propagation
        xa: float, dimension (8,:), initial relative state in km and km/s with two padding zeros

        Returns
        ------
        stm: float (8,8), approximate first order STM of the CHASER trajectory
        stt: float (8,8,8), approximate second order STT of the CHASER trajectory
        """
        stm, stt = pq.pw_stt_update(ta,tb,xa)
        return stm,stt
    def tensorChangeBasis(self,RNewOldf, ROldNew0, old_stm, old_stt):
        """
        Transform STM and STT from an old basis to a new basis. Rotation matrix epochs
        should correspond to the times at which the original STM and STT were evaluated.
        Inputs
        ------
        RNewOldf: float, (8,8) Matrix with rotation matrix in the upper left 6x6 block
                               to transform from the new frame to the old frame at time tf
        ROldNew0: float, (8,8) Matrix with rotation matrix in the upper left 6x6 block
                               to transform from the old frame to the new frame at time t0
        old_stm: float, dimension (8,8), STM in old coordinate frame from time t0 to tf
        old_stt: float, dimension (8,8,8), STT in old coordinate frame from time t0 to tf

        Returns
        ------
        new_stm: float (8,8),  STM in new coordinate frame from time t0 to tf
        new_stt: float (8,8,8),  STT in new coordinate frame from time t0 to tf
        """
        return pq.pw_tensor_change_of_basis(RNewOldf, ROldNew0, old_stm, old_stt)
    def zMap(self,t,order):
        """
        Propagates a test state to a time after the model epoch and then back to the
        origin. The result should be zero.

        Inputs
        ------
        t: float, in seconds past J2000 TDB
        order: 1 or 2

        Returns
        ------
        some numbers very close to zero
        """
        return pq.pw_zmap(t,order)
