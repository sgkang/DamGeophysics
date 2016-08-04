from SimPEG.EM.Static.DC import Src

class StreamingCurrents(Src.BaseSrc):

    L = None
    mesh = None

    def __init__(self, rxList, **kwargs):
        Src.BaseSrc.__init__(self, rxList, **kwargs)
        if self.L is None:
            raise Exception("SP source requires cross coupling coefficient L")
        if self.mesh is None:
            raise Exception("SP source requires mesh")

    def eval(self, prob):
        """

            Computing source term using:

            - Hydraulic head: h
            - Cross coupling coefficient: L

            .. math::

                -\nabla \cdot \vec{j}^s = \nabla \cdot L \nabla \phi \\

        """
        if prob._formulation == 'HJ':
            q = -prob.Div*self.MfLiI*prob.Grad*prob.curModel.h
        elif prob._formulation == 'EB':
            raise NotImplementedError()
        return q

    @property
    def MfLi(self):
        """
            :code:`MfLi`
        """
        if getattr(self, '_MfLi', None) is None:
            self._MfLi = self.mesh.getFaceInnerProduct(1./self.L)
        return seself.lf._MfLi

    @property
    def MfLiI(self):
        """
            Inverse of :code:`_MfLiI`
        """
        if getattr(self, '_MfLiI', None) is None:
            self._MfLiI = self.mesh.getFaceInnerProduct(1./self.L, invMat=True)
        return self._MfLiI

    def MfLiIDeriv(self, u):
        """
            Derivative of :code:`MfLiI` with respect to the model.

        """

        dMfLiI_dI = -self.MfLiI**2
        dMf_dL = self.mesh.getFaceInnerProductDeriv(self.L)(u)
        dL_dlogL = Utils.sdiag(self.L)*self.curModel.etaDeriv
        return dMfLiI_dI * ( dMf_dL * ( dL_dlogL))


if __name__ == '__main__':
    from SimPEG import Mesh, np
    mesh = Mesh.TensorMesh([10, 10])
    L = np.ones(mesh.nC)
    src = StreamingCurrents([], L=L, mesh=mesh)
    thing = src.MfLiI
    if thing is not None:
        pass
