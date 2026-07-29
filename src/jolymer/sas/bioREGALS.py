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

class bioEFA(efa):
    def __init__(self, I, sigma=np.array([]), uv=None):
        super().__init__(I, sigma)
        if uv is not None:
            self.wl = uv.iloc[:,0].values
            self.uv = uv.iloc[:,1:].values
            if self.uv.shape[1] != self.Nc:
                raise ValueError('uv frame count does not match I column count')
        else:
            self.wl = np.array([])
            self.uv = np.array([])

    @property
    def Nuv(self):
        return self.uv.shape[0]

    @property
    def Nframes(self):
        return self.uv.shape[1]

    def select_uv_range(self, wl_min=220, wl_max=300):
        mask = np.logical_and(self.wl >= wl_min, self.wl <= wl_max)
        self.wl = self.wl[mask]
        self.uv = self.uv[mask,:]

    def uv_svd(self, k=None, cols=None):
        if cols is None:
            cols = np.arange(self.Nframes)
        B = self.uv[:,cols]
        if k is None:
            k = min(B.shape)
        U, s, Vt = np.linalg.svd(B, full_matrices=False)
        k = min(k, len(s))
        return U[:,:k], s[:k], Vt[:k,:]

    def constrained_factors(self, k, cols=None):
        if cols is None:
            cols = np.arange(self.Nc)
        A = self.I[:,cols]
        if len(self.sigma) != 0:
            if self.sigma.shape[0] != self.Nr:
                raise ValueError('sigma does not match number of rows in I')
            if len(self.sigma.shape) == 1 or self.sigma.shape[1] == 1:
                A = np.diag(1/self.sigma) @ A
            else:
                A = np.diag(1/np.mean(self.sigma[:,cols],1)) @ A
        Us, ss, Vst = np.linalg.svd(A, full_matrices=False)
        assert ss.shape[0] == Vst.shape[0]
        r = min(k, ss.shape[0])
        Us, ss, Vst = Us[:,:r], ss[:r], Vst[:r,:]
        C_saxs = np.diag(ss) @ Vst
        Uuv, suv, Vuv = self.uv_svd(k=r, cols=cols)
        C_uv = Vuv
        scale = np.zeros(r)
        for i in range(r):
            x = C_uv[i]
            y = C_saxs[i]
            scale[i] = np.dot(x,y) / np.dot(x,x)
        C = np.diag(scale) @ C_uv
        F = Us
        return {
            "F": F,
            "C": C,
            "scale": scale,
            "uv_spectra": Uuv,
            "uv_singular": suv,
            "cols": cols,
            "k": r
        }

    def evolving_uv_factors(self, k, direction="forward", skip=1):
        ncols = self.Nc
        factors = np.full((k,ncols), np.nan)
        for j in range(0,ncols,skip):
            if direction.lower() == "forward":
                cols = np.arange(j+1)
            elif direction.lower() == "reverse":
                cols = np.arange(j,ncols)
            else:
                raise ValueError("direction must be forward/reverse")
            if len(cols) < k:
                continue
            result = self.constrained_factors(k, cols=cols)
            factors[:,j] = np.linalg.norm(result["C"], axis=1)
        return factors
