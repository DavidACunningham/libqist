# Itraj.py
# Thin wrapper interface for QIST Python
# Author: David Cunningham
import numpy as np
from pyqist import pq

class Itraj(object):
    """
    class Itraj
    Description:
        Object-oriented wrapper for the QIST propagation runtime.
        Basic usage: initialize from a namelist file as specified in the users
        guide and then use the functions below to propagate relative motion
        from the model specified by the namelist.

        NOTE: Currently Python can only load one QIST model at a time per process, 
        so something like 
        it1 = Itraj(nml1)
        it2 = Itraj(nml2)
        will result in the statement it1==it2 evaluating to True.
    """
    def __init__(self, namelist):
        """
        Itraj initializer
        inputs: 
               namelist (string) -- a file on disk containing the ITRAJ_CONFIG namelist
        Description: Loads the specified namelist into memory for use by QIST.

        """
        pq.pw_init_n(namelist)
    def state(self,t):
        """
        Itraj.state(t)
        inputs: 
               t (float) -- time in seconds past J2000
        returns: 
               state (float)
        Description: returns the ``state'' stored in the QIST model. In the current version
            of QIST, this is just the value of an independent variable.

        """
        return pq.pw_state(t)
    def stm(self,t):
        """
        Itraj.stm(t)
        inputs: 
               t (float) -- time in seconds past J2000
        returns: 
               stm (8x8 ndarray)
        Description: returns the STM around the reference trajectory from time t0 to time t
            where t0 is the initial time specified in the namelist.
            The 8 rows/columns of the STM correspond to:
            [x y z xdot ydot zdot t0 TOF]

            If you don't need time partials then just do:
            this_stm = it.stm(t)[:6,:6]

        """
        return pq.pw_stm(t)
    def stt(self,t):
        """
        Itraj.stt(t)
        inputs: 
               t (float) -- time in seconds past J2000
        returns: 
               stt (8x8x8 ndarray)
        Description: returns the STT around the reference trajectory from time t0 to time t
            where t0 is the initial time specified in the namelist.
            The 8 pages/rows/columns of the STM correspond to:
            [x y z xdot ydot zdot t0 TOF]

            If you don't need time partials then just do:
            this_stt = it.stt(t)[:6,:6,:6]

        """
        return pq.pw_stt(t)
    def stmInv(self,t):
        """
        Itraj.stmInv(t)
        inputs: 
               t (float) -- time in seconds past J2000
        returns: 
               invstm (8x8 ndarray)
        Description: returns the inverse STM around the reference trajectory from time t to time t0
            where t0 is the initial time specified in the namelist.
        """
        return pq.pw_stm_i(t)
    def sttInv(self,t):
        """
        Itraj.sttInv(t)
        inputs: 
               t (float) -- time in seconds past J2000
        returns: 
               invstt (8x8x8 ndarray)
        Description: returns the inverse STT around the reference trajectory from time t to time t0
            where t0 is the initial time specified in the namelist.
        """
        return pq.pw_stt_i(t)
    def prop(self,ta,tb,xa,order=None):
        """
        Itraj.prop(ta,tb,xa,order=2)
        inputs: 
                ta (float) -- initial time in seconds past J2000
                tb (float) -- final time in seconds past J2000
                xa (8 x n ndarray) -- Initial relative states to propagate between times ta and tb
                                      one column should contain a vector 
                                      [dx, dy, dz, dxdot, dydot, dzdot, dt0, dTOF]
                                      to be propagated to time tb.
                                      The reference frame must be in SPICE J2000 (ICRF), 
                                      km, km/s, and s, relative to the reference trajectory.
                order (integer) -- Order of propagation, should be 1 or 2. 
                                   Other values will throw an error. 
                                   1 will cause QIST to propagate linearly, 
                                   i.e. with the STM only. 
                                   2 will cause QIST to propagate a quadratic
                                   solution, i.e. using the STM and STT.
        returns: 
                xb (8 x n ndarray) -- the final states at time tb.
        Description: The core function of QIST: propagates relative states between times.
        """
        ordr= 2
        func = {1: pq.pw_prop_once, 2: pq.pw_prop_many}
        r = xa.ndim
        if order is not None : ordr=order
        return func[r](ta,tb,xa, ordr)
    def sttsAToB(self,ta,tb):
        """
        Itraj.sttsAToB(ta,tb)
        inputs: 
                ta (float) -- initial time in seconds past J2000
                tb (float) -- final time in seconds past J2000
        returns: 
                (STMab, STTab) :: ((8x8 ndarray), (8x8x8 ndarray)) -- STM and STT from 
                                   time ta to time tb in SPICE J2000 frame in km, km/s, and s.
        Description: Returns the "chained" STM and STT between any two times in the simulation for targeting,
                     optimization, uncertainty propagation, or any other application requiring partials of the flow.
        """
        scratchm, scratcht = pq.pw_stts_ab(ta,tb)
        stm = scratchm
        stt = scratcht #np.einsum('kji -> ijk', scratcht)
        return stm,stt
    def sttUpdate(self,ta,tb,xa):
        """
        Itraj.sttUpdate(ta,tb,xa)
        inputs: 
                ta (float) -- initial time in seconds past J2000
                tb (float) -- final time in seconds past J2000
                xa (8 x n ndarray) -- Initial relative states to propagate between times ta and tb
                                      one column should contain a vector 
                                      [dx, dy, dz, dxdot, dydot, dzdot, dt0, dTOF]
                                      to be propagated to time tb.
                                      The reference frame must be in SPICE J2000 (ICRF), 
                                      km, km/s, and s, relative to the reference trajectory.
        returns: (STMab_new, STTab_new) :: ((8x8 ndarray), (8x8x8 ndarray)) -- STM and STT from 
                                   time ta to time tb in SPICE J2000 frame in km, km/s, and s.
        Description: Returns an approximation of the STM and STT of the _relative_ trajectory using the initial
                     relative state and the _reference_ trajectory STM and STT.
        """
        stm, stt = pq.pw_stt_update(ta,tb,xa)
        return stm,stt
    def tensorChangeBasis(self,RNewOldf, ROldNew0, old_stm, old_stt):
        """
        Itraj.tensorChangeBasis(RNewOldf, ROldNew0, old_stm, old_stt)
        inputs: 
                RNewOldf (8 x 8 ndarray) -- Transformation matrix from the old basis to the new basis at time tf
                ROldNew0 (8 x 8 ndarray) -- Transformation matrix from the new basis to the old basis at time t0
                old_stm  (8 x 8 ndarray) -- STM from time ta to tb in old basis
                old_stt  (8 x 8 x 8 ndarray) -- STT from time ta to tb in old basis
        returns: (STMab_new, STTab_new) :: ((8x8 ndarray), (8x8x8 ndarray)) -- STM and STT from 
                                   time ta to time tb in km, km/s, and s in new basis.
        Description: This function applies a rotation matrix or other linear transformation to the STM and STT
                     to change their basis. Useful for when a user wants to propagate or otherwise operate 
                     in a frame other than J2000.
                    
        """
        return pq.pw_tensor_change_of_basis(RNewOldf, ROldNew0, old_stm, old_stt)
