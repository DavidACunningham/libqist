from addlib import *
from abc import ABC, abstractmethod
from dac3B import LU
import attrs
import numpy as np
from coordrots import rswfromijk
import numpy.linalg as l
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itraj import Itraj, defaultQIST
from rel_state_methods import rel_state_methods as rs
from copy import deepcopy

def rfunc(b): return np.array([*b, 0., 0., 0.])
def vfunc(b): return np.array([0., 0., 0.,*b])
def randvfunc(b,c): return np.array([*b,*c])
km = 50
eps = km/(LU*np.sqrt(3))
defaults = {"order" : 2,
            "x_rel_a" : eps*np.array([1, 1, 1, 0, 0, 0]),
            "t_a" : 1e-8}

@attrs.define(kw_only=True)
class RelState(ABC):
    """
    relstate: abstract class for propagating relative states

    relstate: abstract class for propagating relative states
    with one of several propagation models. This class contains
    some support functions, but note that the method that does
    most of the heavy lifting for each class (propagate) is
    an abstract method and needs definition in the subclasses.

    """
    order: int
    x_rel_a: np.ndarray
    t_a: float
    label: str = attrs.field(init=False, default="")
    it = Itraj(**defaultQIST)
    def __attrs_post_init__(self):
        # Load itraj here!
        self.x_targ_a = self.it.state(self.t_a)
    @abstractmethod
    def propagate(self,t_b):
        """
        propagate: propagate the relative state given a t
        """
    def compare_truth(self,t_b):
        """
        compare_truth: compares current relative state to truth model relative state, if it exists.
        """
    def read_truth(self,filename):
        """
        read_truth: read in truth model from file
        """
    def GWState(self, t):
        """
        GWState: Return the Gateway state at a given value of t (time)

        """
        return self.it.state(t)

    def generate_tag(self):
        """
        generate_tag: generate unique id tag for initial conditions for later saving and retrieval

        MUST BE STATIC after array is instantiated
        """
        tag_str = ''.join([*self.x_rel_a[:6],self.t_a])
        return hex(hash(tag_str))

@attrs.define(kw_only=True,slots=False)
class RelTruth(RelState):
    """
    RelTruth: class for truth model relative states

    RelTruth: class implementing truth model relative states.
    Will look up Gateway state, but must integrate relative states not previously
    used. If relative state is previously used (as determined by UID), will load
    previously used state.
    """
    rtol = 1e-15
    atol = 1e-22

    def save_truth(self):
        """ save_truth: saves data to a truth trajectory file """
    def propagate(self,t_b):
        """ propagate: Propagate the truth trajectory """
        # Check if previously integrated
        # If previously integrated, get chaser state from file
        # If not previously integrated, numerically integrate and save state to file
        dx_b_one = rs.prop_truth_t
        dx_b_many = rs.prop_ts
        dx_b = {1: dx_b_one, 2: dx_b_many}
        args = [t_b,self.t_a,self.x_targ_a[:-1], self.x_rel_a]
        # return chaser state - GW state
        return dx_b[self.x_rel_a.ndim](*args).T
    def check_save_tb(self):
        """ check_save_tb: see if the current trajectory is saved through tb """

@attrs.define(kw_only=True,slots=False)
class RelQist(RelState):
    # RelQist: Propagate QIST trajectories, but with a uniform interface
    def __init__(self,order,x_rel_a,t_a):
        # Load itraj here!
        super().__init__(order=order,x_rel_a=x_rel_a, t_a=t_a)
    def propagate(self,t_b,order=2):
        """ Propagate to t_b.
        """
        return self.it.prop(self.t_a, t_b, self.x_rel_a, order)[:6].T
    def truth(self):
        return RelTruth(order=2,x_rel_a=self.x_rel_a, t_a=self.t_a)
@attrs.define(kw_only=True,slots=False)
class RelTimeSTT(RelState):
    # RelTimeSTT: Propagate time STM and STT for comparison
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
    def propagate(self,t_b, *args, **kwargs):
        # propagate: Propagate the truth trajectory
        # Check if previously integrated
        # If previously integrated, get chaser state from file
        # If not previously integrated, numerically integrate and save state to file
        try:
            order = kwargs["order"]
        except:
            order = 2
        pfun = rs.prop_fo_so_t
        foout, soout = pfun(t_b,self.t_a,self.x_targ_a,self.x_rel_a)
        if order==1:
            return foout.T
        else:
            return soout.T
class TrajBall:
    def __init__(self, initsep, npoints, t_a, states=None, initvel=None, order=2,rv="r", mode="qist",_fromtruediv=False):
        self.knd = {"qist" : RelQist, 
               "truth": RelTruth,
               "time" : RelTimeSTT}
        self.rv = rv
        self.rvs = {"r": rfunc,
                    "v": vfunc,
                    'randv': randvfunc}
        if states is not None:
            self.ics = np.squeeze(np.array(states))
            self.rels = [self.knd[mode](order=order,
                                   x_rel_a=self.ics.T,
                                   t_a=t_a) 
                    for s in states]
        else:
            if self.rv=='randv' and isinstance(initsep,tuple):
                if len(initsep) != 2:
                    print("error: need two values for randv mode")
                self.ics = np.array([self.rvs[rv](i,j) 
                           for i,j in 
                           zip(self.makeball(initsep[0],npoints),
                               self.makeball(initsep[1],npoints))])
            else:
                self.ics = np.array([self.rvs[rv](i) 
                           for i in 
                           self.makeball(initsep,npoints)])
            # ics has shape (npoints,6)
            self.rels = [self.knd[mode](order=order,
                                   x_rel_a=np.array(self.ics).T,
                                   t_a=t_a) ]
        self.initsep = initsep
        self._fromtruediv = _fromtruediv
        self.t_a = t_a
        self.mode=mode
        self.order = order
        self.npoints=npoints
        if mode == "time" or mode == "truth" : 
            self.label = mode
        else:
            self.label=''.join([mode,", order ",str(order)])

    def __truediv__(me,you):
        mestates = me.ics
        youstates = you.ics
        ourstates = []
        for m,y in zip(mestates,youstates):
            rnorm = l.norm(y[:3])
            vnorm = l.norm(y[3:6])
            ourstates.append(np.r_[m[:3]/rnorm, m[3:6]/vnorm])
        return TrajBall(me.initsep,me.npoints,me.t_a, 
                        states=ourstates,
                        order=me.order, 
                        rv=me.rv,
                        mode=me.mode,
                        _fromtruediv=True)
    def __sub__(me,you):
        mestates = me.ics
        youstates = you.ics
        ourstates = [m[:6] - y[:6] for m,y in zip(mestates,youstates)]
        return TrajBall(me.initsep,me.npoints,me.t_a, 
                        states=ourstates,
                        order=me.order, 
                        rv = me.rv,
                        mode=me.mode)
    def makeball(self, r, npoints):
        # Fibbonaci spiral method adjusted
        # final ball shape: (npoints,3)
        # See Martin Roberts blog
        # http://extremelearning.com.au/evenly-distributing-points-on-a-sphere/
        i = np.arange(0,npoints,dtype=float) + 0.5
        phi = np.arccos(1-2*(i+0.5)/npoints)
        gr = (1+np.sqrt(5))/2
        theta = 2*np.pi*i/gr
        ball = []
        ball.append(np.cos(theta)*np.sin(phi))
        ball.append(np.sin(theta)*np.sin(phi))
        ball.append(np.cos(phi))
        ball = r*np.array(ball).T
        return ball
    
    def get_ivs(self,t_b, **kwargs):
        # push: Propagate a ball of states forward in time to tb
        # using the propagate method in the underlying Rel**** instances.
        # Uses multiprocessing.
        # kwargs dict catches the optional arguments for the Rel****
        # method.
        # Returns a new TrajBall
        states_at_tb = [r.propagate(t_b,**kwargs) for r in self.rels]
        states_at_tb = np.squeeze(np.asarray(states_at_tb))
        times = [s[-1] for s in states_at_tb]
        return np.array(times)

    def push(self, t_b, **kwargs):
        # push: Propagate a ball of states forward in time to tb
        # using the propagate method in the underlying Rel**** instances.
        # Uses multiprocessing.
        # kwargs dict catches the optional arguments for the Rel****
        # method.
        # Returns a new TrajBall
        states_at_tb = [r.propagate(t_b,**kwargs) for r in self.rels]
        return TrajBall(self.initsep,self.npoints,t_b, 
                        states=states_at_tb, 
                        rv = self.rv,
                        order=self.order, 
                        mode=self.mode)
                        
    def plot(self,mode="3d",rotmat=None,show=False,ax=None,title=False,timeintitle=False,alpha=1,color=None,marker=None,markersize=None,units=False,label=False):
        
        # plot: plot the current ball of states on a 3D axis
        if rotmat is not None:
            rot = rotmat
        else:
            rot = np.eye(3)

        rotstates = rot@self.ics.T[:3]
        un = ""
        if units and not self._fromtruediv:
            rotstates = rotstates*LU
            un = " (km)"
        elif units:
            un = " (normalized)"
        if ax is None:
            fig = plt.figure()
            if mode=="3d":
                ax = fig.add_subplot(111, projection='3d')
            else:
                ax = fig.add_subplot(111)
        if mode =="3d":
            ax.scatter(*rotstates,alpha=alpha,color=color,marker=marker,s=markersize)
            ax.elev=10
            ax.azim = 0
            if label and rotmat is not None:
                ax.set_xlabel("R"+ un)
                ax.set_ylabel("S" + un)
                ax.set_zlabel("W" + un)
            elif label:
                ax.set_xlabel("X" + un)
                ax.set_ylabel("Y" + un)
                ax.set_zlabel("Z" + un)
        elif mode=="xz":
            ax.scatter(*rotstates[::2],alpha=alpha,color=color,marker=marker,s=markersize)
            if label and rotmat is not None:
                ax.set_xlabel("R" + un)
                ax.set_ylabel("W" + un)
            elif label:
                ax.set_xlabel("X" + un)
                ax.set_ylabel("Z" + un)
        elif mode=="xy":
            ax.scatter(*rotstates[:2],alpha=alpha,color=color,marker=marker,s=markersize)
            if label and rotmat is not None:
                ax.set_xlabel("R" + un)
                ax.set_ylabel("S" + un)
            elif label:
                ax.set_xlabel("X" + un)
                ax.set_ylabel("Y" + un)
        elif mode=="yz":
            ax.scatter(*rotstates[1:],alpha=alpha,color=color,marker=marker,s=markersize)
            if label and rotmat is not None:
                ax.set_xlabel("S" + un)
                ax.set_ylabel("W" + un)
            elif label:
                ax.set_xlabel("Y" + un)
                ax.set_ylabel("Z" + un)
        if timeintitle:
            ax.set_title(''.join([self.label," $t = ", str(self.t_a), "$"]))
        elif title:
            ax.set_title(self.label)
        if ax is None:
            if show:
                plt.show()
            else:
                return (fig,ax)


    def clone(self,newmode,order=2):
        """
        clone: Makes a ball of identical IC's with a new rel mode

        """ 
        return TrajBall(self.r,self.npoints,self.t_a, 
                        states=self.ics,
                        order=order,
                        mode=newmode)

    def stats(self,mode='pos'):
        if mode == 'pos':
            norm = l.norm(self.ics[:,:3],axis=1)
        elif mode=='vel':
            norm = l.norm(self.ics[:,3:],axis=1)
        else:
            norm = l.norm(self.ics,axis=1)
        # return mean magnitude
        avg = np.mean(norm)
        std = np.std(norm)
        maximum = norm.max()
        minimum = norm.min()
        # return max-min magnitude
        # return stdev magnitude?
        return (avg, 3*std, maximum, minimum)

if __name__=="__main__":
    from time import time
    from dac3B import TU

    sep= 500
    proj = "yz"
    num = 500
    inter_times = [10,  65, 77.838, 78.5, 79.665 ]
    itimes = [t*3600/TU for t in inter_times]
    times =  itimes + [1.5]
    # STM, STM+STT, TRUTH, TRUTH-(STM), TRUTH-(STT+STM) 
    times_h = [t*TU/3600 for t in times]
    b_true = TrajBall(sep/LU,num,0.01, mode='truth')
    b_qist = TrajBall(sep/LU,num,0.01, mode='qist')
    
    allballs = {"truth": b_true, 
                "qist1_2": b_qist, 
                "qist2_2": b_qist}
    nmodels = len(list(allballs.keys()))
    pushes = {}
    for k,ball in allballs.items():
        print("push", k)
        tstart = time()
        w = False
        if k=="qist1_2":
            kwargs = {"order":1}
            w = True
        if w:
            pushes[k]   = [ball.push(t, **kwargs) for t in times]
        else:
            pushes[k]   = [ball.push(t) for t in times]
        tstop = time()
        print(tstop-tstart,"s elapsed")
        for s in pushes[k]:
            s.label = k
    tobediffed = {k:v for k,v in pushes.items() if k!="truth"}
    diffs = {}
    for k,approx in tobediffed.items():
        diffs[k] = [(a-truth)/truth for a,truth in zip(approx,pushes["truth"])]
        for s in diffs[k]:
            s.label = ''.join([k,'-truth'])
    print("plotting")
    fig,axs = plt.subplots(len(times),1,figsize=(4*(2*nmodels-1),4*len(times)),constrained_layout=True,dpi=120)
    for ax in axs:
        ax.remove()
    gridspec = axs[0].get_subplotspec().get_gridspec()
    subfigs = [fig.add_subfigure(gs) for gs in gridspec]
    axs = []
    for t, subfig in zip(times_h,subfigs):
        subfig.suptitle(f't = {t} hrs')
        axs.append(subfig.subplots(1,2*nmodels-1))
    gw = RelTruth(**defaults)

    for i,q in enumerate(zip(*list(pushes.values()),*list(diffs.values()))):
        for a,x in zip(axs[i],q):
            T = rswfromijk(gw.GWState(times[i]))
            x.plot(rotmat=T,show=False,mode=proj,ax=a)
    savestring = f"spheres_{sep}km_" + proj + ".png"
    plt.savefig(savestring,format="png", pad_inches=0.0001)
