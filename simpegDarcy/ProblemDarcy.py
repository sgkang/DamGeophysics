from SimPEG import Problem, Utils, Maps, Mesh
from SimPEG.Utils import sdiag
import numpy as np
from SimPEG.Utils import Zero
from SimPEG.EM.Static.DC import getxBCyBC_CC
from SimPEG.EM.Base import BaseEMSurvey
from SimPEG import Props, Survey, Solver as SimpegSolver

class BaseDarcyProblem(Problem.BaseProblem):

    K, KMap, KDeriv = Props.Invertible(
        "Hydraulic conductivity (m/s)"
    )

    Ki, KiMap, KDeriv = Props.Invertible(
        "Hydraulic resistivity (s/m)"
    )

    Props.Reciprocal(K, Ki)

    surveyPair = Survey.BaseSurvey  #: The survey to pair with.
    dataPair = Survey.Data  #: The data to pair with.

    mapPair = Maps.IdentityMap

    Solver = SimpegSolver
    solverOpts = {}

    verbose = False

    Ainv = None
    sigma = None
    f = None
    Ainv = None

    def fields(self, m=None):
        if m is not None:
            self.model = m

        if self.Ainv is not None:
            self.Ainv.clean()

        A = self.getA()
        self.Ainv = self.Solver(A, **self.solverOpts)
        RHS = self.getRHS()
        u = self.Ainv * RHS
        return u

    @property
    def deleteTheseOnModelUpdate(self):
        toDelete = []
        return toDelete

    @property
    def MfK(self):
        """
            Face inner product matrix for \\(\\Ki\\). Used in the H-J formulation
        """
        if getattr(self, '_MfK', None) is None:
            self._MfK = self.mesh.getFaceInnerProduct(self.K)
        return self._MfK

    @property
    def MfKi(self):
        """
            Face inner product matrix for \\(\\Ki\\). Used in the H-J formulation
        """
        if getattr(self, '_MfKi', None) is None:
            self._MfKi = self.mesh.getFaceInnerProduct(self.Ki)
        return self._MfKi

    # TODO: This should take a vector
    def MfKiDeriv(self, u):
        """
        Derivative of :code:`MfKi` with respect to the model.
        """
        return self.mesh.getFaceInnerProductDeriv(self.Ki)(u) * self.KiDeriv

    @property
    def MfKiI(self):
        """
        Inverse of :code:`MfKi`
        """
        if getattr(self, '_MfKiI', None) is None:
            self._MfKiI = self.mesh.getFaceInnerProduct(self.Ki, invMat=True)
        return self._MfKiI

    def MfKiIDeriv(self, u):
        """
            Derivative of :code:`MfKiI` with respect to the model.
        """
        dMfKiI_dI = -self.MfKiI**2
        dMf_dKi = self.mesh.getFaceInnerProductDeriv(self.Ki)(u)
        return dMfKiI_dI * (dMf_dKi * self.KiDeriv)

    def gradh (self, u):
        return -self.mesh.cellGrad*u

    def vel (self, u):
        MfI = self.mesh.getFaceInnerProduct(invMat=True)
        return self.MfKiI*self.Grad*u

    def p (self, u):
        return u - self.mesh.gridCC[:, 2]

    def divgradh(self, u):
        return -self.mesh.faceDiv*self.mesh.cellGrad*u


class Problem_CC(BaseDarcyProblem):

    _solutionType = 'hSolution'
    _formulation = 'HJ'  # CC potentials means J is on faces
    modelType = None

    def __init__(self, mesh, **kwargs):
        BaseDarcyProblem.__init__(self, mesh, **kwargs)
        self.mesh.setCellGradBC("neumann")
        # if self.Ki is None:
        #     raise Exception("Resistivity:Ki needs to set when \
        #                      initializing SPproblem")
        # self.setBC()

    def getA(self):
        """

        Make the A matrix for the cell centered DC resistivity problem

        A = D MfKiI G

        """

        D = self.Div
        G = self.Grad
        MfKiI = self.MfKiI
        A = D * MfKiI * G
        return A

    def getADeriv(self, u, v, adjoint= False):
        # We assume conductivity is known
        return Zero()

    def getRHS(self):
        """
        RHS for the Darcy problem

        """
        RHS = self.Div*self.MfKiI*self.P_BC*self.x_BC
        return RHS

    def getRHSDeriv(self, src, v, adjoint=False):
        """
        Derivative of the right hand side with respect to the model
        """
        return Zero()

    def mapBC(self, bc):
        alpha = []
        beta = []
        for bc_temp in bc:
            if bc_temp == "neumann":
                alpha.append(0.)
                beta.append(1.)
            elif bc_temp == "dirichlet":
                alpha.append(1.)
                beta.append(0.)
            else:
                raise Exception("Only n or d")
        return alpha, beta

    #This needs to generalized
    def setBC(self, bc, h):
        """
        Set up boundary conditions of 2D and 3D problem

        """
        if self.mesh.dim == 3:
            if len(bc) and len(h) != 3:
                raise Exception("Lenth of bc and h should be 3 for 3D")

            # Set gamma value
            hx = h[0]
            hy = h[1]
            hz = h[2]

            # Set alpha and beta value
            alphax, betax = self.mapBC(bc[0])
            alphay, betay = self.mapBC(bc[1])
            alphaz, betaz = self.mapBC(bc[2])

            fxm, fxp, fym, fyp, fzm, fzp = self.mesh.faceBoundaryInd
            gBFxm = self.mesh.gridFx[fxm,:]
            gBFxp = self.mesh.gridFx[fxp,:]
            gBFym = self.mesh.gridFy[fym,:]
            gBFyp = self.mesh.gridFy[fyp,:]
            gBFzm = self.mesh.gridFz[fzm,:]
            gBFzp = self.mesh.gridFz[fzp,:]

            # Setup Mixed B.C (alpha, beta, gamma)
            temp_xm, temp_xp = np.ones_like(gBFxm[:,0]), np.ones_like(gBFxp[:,0])
            temp_ym, temp_yp = np.ones_like(gBFym[:,1]), np.ones_like(gBFyp[:,1])
            temp_zm, temp_zp = np.ones_like(gBFzm[:,2]), np.ones_like(gBFzp[:,2])

            alpha_xm, alpha_xp = temp_xm*alphax[0], temp_xp*alphax[1]
            alpha_ym, alpha_yp = temp_ym*alphay[0], temp_yp*alphay[1]
            alpha_zm, alpha_zp = temp_zm*alphaz[0], temp_zp*alphaz[1]

            beta_xm, beta_xp = temp_xm*betax[0], temp_xp*betax[1]
            beta_ym, beta_yp = temp_ym*betay[0], temp_yp*betay[1]
            beta_zm, beta_zp = temp_zm*betaz[0], temp_zp*betaz[1]

            gamma_xm, gamma_xp = temp_xm*hx[0], temp_xp*hx[1]
            gamma_ym, gamma_yp = temp_ym*hy[0], temp_yp*hy[1]
            gamma_zm, gamma_zp = temp_zm*hz[0], temp_zp*hz[1]

            alpha = [alpha_xm, alpha_xp, alpha_ym, alpha_yp, alpha_zm, alpha_zp]
            beta =  [beta_xm, beta_xp, beta_ym, beta_yp, beta_zm, beta_zp]
            gamma = [gamma_xm, gamma_xp, gamma_ym, gamma_yp, gamma_zm, gamma_zp]

        elif self.mesh.dim==2:
            if len(bc) and len(h) != 2:
                raise Exception("Lenth of bc and h should be 2 for 2D")

            # Set gamma value
            hx = h[0]
            hy = h[1]

            # Set alpha and beta value
            alphax, betax = self.mapBC(bc[0])
            alphay, betay = self.mapBC(bc[1])

            fxm, fxp, fym, fyp = self.mesh.faceBoundaryInd
            gBFxm = self.mesh.gridFx[fxm,:]
            gBFxp = self.mesh.gridFx[fxp,:]
            gBFym = self.mesh.gridFy[fym,:]
            gBFyp = self.mesh.gridFy[fyp,:]

            # Setup Mixed B.C (alpha, beta, gamma)
            temp_xm, temp_xp = np.ones_like(gBFxm[:,0]), np.ones_like(gBFxp[:,0])
            temp_ym, temp_yp = np.ones_like(gBFym[:,1]), np.ones_like(gBFyp[:,1])

            alpha_xm, alpha_xp = temp_xm*alphax[0], temp_xp*alphax[1]
            alpha_ym, alpha_yp = temp_ym*alphay[0], temp_yp*alphay[1]

            beta_xm, beta_xp = temp_xm*betax[0], temp_xp*betax[1]
            beta_ym, beta_yp = temp_ym*betay[0], temp_yp*betay[1]

            gamma_xm, gamma_xp = temp_xm*hx[0], temp_xp*hx[1]
            gamma_ym, gamma_yp = temp_ym*hy[0], temp_yp*hy[1]

            alpha = [alpha_xm, alpha_xp, alpha_ym, alpha_yp]
            beta =  [beta_xm, beta_xp, beta_ym, beta_yp]
            gamma = [gamma_xm, gamma_xp, gamma_ym, gamma_yp]

        x_BC, y_BC = getxBCyBC_CC(self.mesh, alpha, beta, gamma)
        V = Utils.sdiag(self.mesh.vol)
        self.Div = V * self.mesh.faceDiv
        P_BC, B = self.mesh.getBCProjWF_simple()
        M = B*self.mesh.aveCC2F
        self.Grad = self.Div.T - P_BC*Utils.sdiag(y_BC)*M
        self.x_BC = x_BC
        self.P_BC = P_BC

class DarcySurvey(Survey.BaseSurvey):
    """docstring for RichardsSurvey"""

    srcList = None

    def __init__(self, rxList, **kwargs):
        self.rxList = rxList
        self.srcList = [Survey.BaseSrc(rxList)]
        Survey.BaseSurvey.__init__(self, **kwargs)

    @property
    def nD(self):
        return rxList[0].nD


class DarcyRx(Survey.BaseRx):
    """
    Base DC receiver
    """
    locs = None
    rxType = None

    knownRxTypes = {
        'h': ['h', None],
        'gradhx': ['gradh', 'x'],
        'gradhy': ['gradh', 'y'],
        'gradhz': ['gradh', 'z'],
        'velx': ['vel', 'x'],
        'vely': ['vel', 'y'],
        'velz': ['vel', 'z'],
    }

    def __init__(self, locs, rxType='h', **kwargs):
        Survey.BaseRx.__init__(self, locs, rxType, **kwargs)
        locs = np.atleast_2d(locs)

    @property
    def projField(self):
        """Field Type projection (e.g. e b ...)"""
        return self.knownRxTypes[self.rxType][0]

    def projGLoc(self, f):
        """Grid Location projection (e.g. Ex Fy ...)"""
        comp = self.knownRxTypes[self.rxType][1]
        if comp is not None:
            return f._GLoc(self.rxType) + comp
        return f._GLoc(self.rxType)

    def eval(self, src, mesh, f):
        P = self.getP(mesh, self.projGLoc(f))
        return P*f[src, self.projField]

    def evalDeriv(self, src, mesh, f, v, adjoint=False):
        P = self.getP(mesh, self.projGLoc(f))
        if not adjoint:
            return P*v
        elif adjoint:
            return P.T*v

    @property
    def nD(self):
        """Number of data in the receiver."""
        return self.locs.shape[0]

    def getP(self, mesh, Gloc):
        if mesh in self._Ps:
            return self._Ps[mesh]

        P = mesh.getInterpolationMat(self.locs, Gloc)

        if self.storeProjections:
            self._Ps[mesh] = P
        return P
