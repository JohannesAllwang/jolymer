"""
"""


import sys
import os
import numpy as np

from .original_regals.regals import *
from .original_regals.efa import *

class bioComponent(component):
    def __init__(self, concentration, profile,
                 uv_scale=0.01):
        self.concentration = concentration
        self.profile = profile
        self.uv_scale = uv_scale

    def to_dict(self):
        return {
            "class": "bioComponent",
            "concentration": self.concentration.to_dict(),
            "profile": self.profile.to_dict(),
            "uv_scale": self.uv_scale,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            concentration=ConcentrationClass.from_dict(d["concentration"]),
            profile=ProfileClass.from_dict(d["profile"]),
            uv_scale=d.get("uv_scale", 0.01),
        )


class bioMIXTURE(mixture):

    def __init__(self, components, lambda_concentration=np.array([]), lambda_profile = np.array([]),
                 u_concentration = [], u_profile = [],
                 uv_meas=None, uv_err=None, uv_weight=1.0):
        super().__init__(components,
                         lambda_concentration=np.array([]),
                         lambda_profile=np.array([]),
                         u_concentration=[],
                         u_profile=[])
        self.uv_meas = uv_meas           # array length Nx (frames)
        self.uv_err = uv_err             # error or weight
        self.uv_weight = uv_weight

    def concentration_problem(self, I, err, calc_Ab=True):
        sout = super().concentration_problem(I, err, calc_Ab=calc_Ab)
        if self.uv_meas is None:
            return sout
        A = [comp.concentration.A for comp in self.components]
        if calc_Ab:
            AA, Ab = sout
        else:
            AA = sout
        Abs = self.uv_meas                    # shape (Nt,)
        uv_scale = [comp.uv_scale for comp in self.components]
        # --- UV normal-equation contribution ---
        w = 1 / np.mean(err,1)
        y = self.profiles
        y = w[:,np.newaxis] * y
        AA = [
            [
                (y[:,k1] @ y[:,k2]) * (A[k1].T @ A[k2])
                + self.uv_weight * self.components[k1].uv_scale * \
                        self.components[k2].uv_scale * (A[k1].T @ A[k2])
                for k2 in range(self.Nc)
            ]
            for k1 in range(self.Nc)
        ]
        D = w[:,np.newaxis] * I
        if calc_Ab:
            Ab = [
                A[k].T @ (D.T @ y[:,k])
                + self.uv_weight * uv_scale[k] * (A[k].T @ Abs)
                for k in range(self.Nc)
            ]

        # Assemble blocks
        AA = sp.vstack(tuple(sp.hstack(tuple(row)) for row in AA))
        if calc_Ab:
            Ab = np.hstack(tuple(Ab))
            return AA, Ab
        else:
            return AA

    def to_dict(self):
        return {
            "class": "bioMIXTURE",
            "components": [c.to_dict() for c in self.components],
            # REGALS / mixture hyperparameters
            "lambda_concentration": list(self.lambda_concentration),
            "lambda_profile": list(self.lambda_profile),
            "u_concentration": list(self.u_concentration),
            "u_profile": list(self.u_profile),
            # UV-specific additions
            "uv_meas": None if self.uv_meas is None else list(self.uv_meas),
            "uv_err": None if self.uv_err is None else list(self.uv_err),
            "uv_weight": self.uv_weight,
        }

    @classmethod
    def from_dict(cls, d):
        components = [
            bioComponent.from_dict(cd)
            for cd in d["components"]
        ]
        return cls(
            components=components,
            lambda_concentration=np.array(d.get("lambda_concentration", [])),
            lambda_profile=np.array(d.get("lambda_profile", [])),
            u_concentration=d.get("u_concentration", []),
            u_profile=d.get("u_profile", []),
            uv_meas=None if d["uv_meas"] is None else np.array(d["uv_meas"]),
            uv_err=None if d["uv_err"] is None else np.array(d["uv_err"]),
            uv_weight=d.get("uv_weight", 1.0),
        )



class bioREGALS(regals):
    def __init__(self, I=None,
                 sigma=None,
                 I_waxs=None,
                 err_waxs=None,
                 onlineUV=None):
        self.I = I
        self.sigma = sigma
        self.err = sigma
        self.I_waxs = I_waxs
        self.err_waxs = err_waxs
        self.onlineUV = onlineUV

