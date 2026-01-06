from __future__ import annotations

import json
import warnings
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Type, Dict

from dataclasses import dataclass, field, asdict
from typing import Optional
from jolymer.samples.bioMOLECULE import bioMOLECULE, DNA
from jolymer.sas.ms import Ms
from jolymer.sas.CoupledMeasurement import CoupledMeasurement
from jolymer.sas.bioREGALS import *
from jolymer.uv.onlineUV import onlineUV
from jolymer.sas.SAXS_Measurement import SAXS_Measurement

RegalsResult = dataclass



class StateEncoder(json.JSONEncoder):
    def default(self, obj: Any):
        if isinstance(obj, Path):
            return {
                "_type": "Path",
                "value": str(obj),
            }
        return super().default(obj)


# ---------- serialization ----------
def dataclass_to_dict(obj: Any) -> Dict[str, Any] | None:
    if obj is None:
        return None
    if not is_dataclass(obj):
        warnings.warn(
            f"Object of type {type(obj).__name__} is not a dataclass; skipping",
            RuntimeWarning,
        )
        return None
    d = asdict(obj)
    d["_type"] = obj.__class__.__name__
    return d


def dict_to_dataclass(
    d: Dict[str, Any],
    registry: Dict[str, Type],
):
    if d is None:
        return None
    if "_type" not in d:
        warnings.warn(
                "Missing _type field in JSON object",
                RuntimeWarning)
        return None
    if d["_type"] == "Ms":
        print("processing Ms")
        mlist = []
        for md in d['ms']:
            path = Path(md.pop('path'))
            m = SAXS_Measurement(path=path, **md)
            mlist.append(m)
        out =  Ms(mlist)
        out.name = d['name']
        return out
    type_name = d.pop("_type")
    print("type ", type_name)
    if type_name not in registry:
        warnings.warn(f"Unknown dataclass type '{type_name}'",
                      RuntimeWarning)
        return None
    try:
        cls = registry[type_name]
        return cls(**d)
    except Exception as e:
        warnings.warn(f"Loading dataclass type '{type_name}' failed",
                      RuntimeWarning)
        return None


# ---------- Path restoration ----------

def restore_special_types(obj: Any):
    """
    Recursively restore special JSON objects like Path.
    """
    if isinstance(obj, dict):
        if obj.get("_type") == "Path":
            return Path(obj["value"])
        return {k: restore_special_types(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [restore_special_types(v) for v in obj]
    else:
        return obj


@dataclass
class BioREGALSState:

    DATACLASS_REGISTRY = {
        # "Sample": Sample, ## this needs to be fixed
        "bioMOLECULE": bioMOLECULE,
        "DNA": DNA,
        "SAXS_Measurement": SAXS_Measurement,
        "Ms": Ms,
        "CoupledMeasurement": CoupledMeasurement,
        "onlineUV": onlineUV,
        "bioMIXTURE": bioMIXTURE,
        "bioREGALS": bioREGALS,
    }
    sample: Optional[bioMOLECULE] = None
    saxs: Optional[Ms] = None
    cm_saxs: Optional[CoupledMeasurement] = None
    uv: Optional[onlineUV] = None
    mixture: bioMIXTURE | None = None

    to_regals: dict = field(default_factory=dict)
    bioREGALS: Optional[bioREGALS] = None
    regals_result: dict = field(default_factory=dict)

    figures: dict = field(default_factory=dict)
    axes: dict = field(default_factory=dict)

    def to_json(self, path: Path):
        outdict = {
            "sample": dataclass_to_dict(self.sample),
            "saxs": dataclass_to_dict(self.saxs),
            "uv": dataclass_to_dict(self.uv),
            "mixture": dataclass_to_dict(self.mixture),
            "bioREGALS": dataclass_to_dict(self.bioREGALS),
            "to_regals": self.to_regals,
        }

        with open(path, "w") as f:
            json.dump(outdict, f, indent=2, cls=StateEncoder)

    def load_from_json( self, path: Path):
        registry = BioREGALSState.DATACLASS_REGISTRY
        with open(path) as f:
            raw = json.load(f)
        raw = restore_special_types(raw)
        print('test 1')
        sample=dict_to_dataclass(raw["sample"], registry)
        if not sample is None:
            self.sample=sample
        saxs=dict_to_dataclass(raw["saxs"], registry)
        if not saxs is None:
            self.saxs=saxs
        uv=dict_to_dataclass(raw["uv"], registry)
        if not uv is None:
            self.uv=uv
        mixture=dict_to_dataclass(raw["mixture"], registry)
        if not mixture is None:
            self.mixture=mixture
        bioREGALS=dict_to_dataclass(raw["bioREGALS"], registry)
        if not bioREGALS is None:
            self.bioREGALS=bioREGALS
        if "to_regals" in raw:
            to_regals = raw['to_regals']
            if not to_regals is None:
                self.to_regals=to_regals
        print('test 2')
        new_cm = CoupledMeasurement(saxs_list=self.saxs, uv=self.uv, sample=self.sample)
        if self.cm_saxs is None:
            self.cm_saxs = new_cm
        else:
            self.cm_saxs.__dict__.update(new_cm.__dict__)
        return self
