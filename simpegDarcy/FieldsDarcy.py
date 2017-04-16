from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import SimPEG
from SimPEG.Utils import Identity, Zero
import numpy as np


class FieldsDarcy(SimPEG.Problem.Fields):
    knownFields = {}
    dtype = float

    def _hDeriv(self, src, du_dm_v, v, adjoint=False):
        if (
            getattr(self, '_hDeriv_u', None) is None or
            getattr(self, '_hDeriv_m', None) is None
        ):
            raise NotImplementedError(
                'Getting hDerivs from {0!s} is not '
                'implemented'.format(self.knownFields.keys()[0])
            )

        if adjoint:
            return (self._hDeriv_u(src, v, adjoint=adjoint),
                    self._hDeriv_m(src, v, adjoint=adjoint))

        return (np.array(self._hDeriv_u(src, du_dm_v, adjoint) +
                         self._hDeriv_m(src, v, adjoint), dtype=float))


class Fields_CC(FieldsDarcy):
    knownFields = {'hSolution': 'CC'}
    aliasFields = {
        'h': ['hSolution', 'CC', '_h'],
        'vel': ['hSolution', 'F', '_vel'],
        'gradh': ['hSolution', 'F', '_gradh']
    }
    # primary - secondary
    # CC variables

    def __init__(self, mesh, survey, **kwargs):
        FieldsDarcy.__init__(self, mesh, survey, **kwargs)
        mesh.setCellGradBC("neumann")
        cellGrad = mesh.cellGrad

    def startup(self):
        self.prob = self.survey.prob

    def _GLoc(self, fieldType):
        if fieldType == 'h':
            return 'CC'
        elif fieldType == 'e' or fieldType == 'j':
            return 'F'
        else:
            raise Exception('Field type must be h, e, j')

    def _h(self, hSolution, srcList):
        return hSolution

    def _hDeriv_u(self, src, v, adjoint=False):
        return Identity()*v

    def _hDeriv_m(self, src, v, adjoint=False):
        return Zero()

    def _vel(self, hSolution, srcList):
        """
            .. math::
                \mathbf{j} = \mathbf{M}^{f \ -1}_{\Ki} \mathbf{G} \h
        """
        return self.prob.MfKiI*self.prob.Grad*hSolution

    def _gradh(self, hSolution, srcList):
        """
            In HJ formulation e is not well-defined!!
            .. math::
                \vec{e} = -\nabla \h
        """
        return -self.mesh.cellGrad*hSolution
