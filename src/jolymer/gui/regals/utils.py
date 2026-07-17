import numpy as np
from jolymer.sas.bioREGALS import bioMIXTURE, bioComponent
from copy import deepcopy


def clone_mixture_for_waxs(
    saxs_mix: bioMIXTURE,
    q_waxs,
    *,
    scale_Nw: float | None = None,
    add_water: bool = True,
):
    """
    Clone a SAXS mixture for WAXS:
      - keep concentrations and x-grid
      - replace q in all profiles
      - optionally rescale Nw
      - always start unsolved
      - optionally add a smooth water component
    """
    d = saxs_mix.to_dict()
    new_components = []
    # -------------------------------------------------
    # clone existing components
    # -------------------------------------------------
    for comp in d["components"]:
        comp = deepcopy(comp)
        prof = comp["profile"]
        if "params" in prof:
            prof["params"]["q"] = np.asarray(q_waxs).tolist()
            if scale_Nw is not None and "Nw" in prof["params"]:
                prof["params"]["Nw"] = max(
                    1, int(scale_Nw * prof["params"]["Nw"])
                )
        new_components.append(comp)
    # -------------------------------------------------
    # add explicit water component
    # -------------------------------------------------
    if add_water:
        x = np.asarray(d["x"])
        xmin, xmax = x.min(), x.max()
        water = {
            "concentration": {
                "kind": "smooth",
                "params": {
                    "x": x.tolist(),
                    "xmin": float(xmin),
                    "xmax": float(xmax),
                    "Nw": 10,
                    "is_zero_at_xmin": False,
                    "is_zero_at_xmax": False,
                },
            },
            "profile": {
                "kind": "smooth",
                "params": {
                    "q": np.asarray(q_waxs).tolist(),
                    "Nw": 10,
                },
            },
            "uv_scale": 0.0,
        }
        new_components.append(water)
    # -------------------------------------------------
    # rebuild mixture dict
    # -------------------------------------------------
    d["components"] = new_components
    # explicitly remove solution state if present
    d.pop("u_profile", None)
    d.pop("u_concentration", None)
    return bioMIXTURE.from_dict(d)

