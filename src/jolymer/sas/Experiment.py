"""
TODO: all
"""

from dataclasses import dataclass

from . import desy
from .. import sasview_utility


@dataclass
class SASExperiment:

    model: sasview_utility.SasModel
    dataloader: desy.Desy

    def save(self):
        pass

    def to_markdown():
        pass
