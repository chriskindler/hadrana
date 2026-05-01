# cora/ensembles/specs.py
from __future__ import annotations
from numpy import s_
import numpy as np

from .ensemble_specs import (
    CONFIGURATION_NUMBERS,
    EXCEPTIONAL_CONFIGURATION_NUMBERS,
    LATTICE_DIMENSIONS,
    REPLICA_ORDER,
    ENSEMBLES_SYM,
    ENSEMBLES_TRM,
    ENSEMBLES_MSC,
    MEASUREMENTS_C2PT,
    MEASUREMENTS_C3PT,
)

class EnsembleHelpers:
    def __init__(self, ensemble: str):
        self.ensemble = ensemble

    def get_spatial_lattice_dimension(self) -> int:
        return LATTICE_DIMENSIONS[self.ensemble].n_space

    def get_temporal_lattice_dimension(self) -> int:
        return LATTICE_DIMENSIONS[self.ensemble].n_time

    def get_trajectory_type(self) -> str:
        if self.ensemble in ENSEMBLES_SYM:
            return "SYM"
        if self.ensemble in ENSEMBLES_TRM:
            return "TRM"
        if self.ensemble in ENSEMBLES_MSC:
            return "MSC"

    def get_config_list(self, replica: int) -> list[int]:
        r"""
        @replica: int,
            integer id of replica as found in r000, r001, ...
        """
        info = CONFIGURATION_NUMBERS[self.ensemble][f"r{replica:0>3}"]
        start: int = info.get("start_cfg", 1)
        step: int = info.get("step_cfg", 1)
        end: int = info["end_cfg"]
        cfg_list = list(range(start, end + 1, step))
        for excl in info["broken_cfg"]:
            cfg_list.remove(excl)
        return cfg_list

    def get_config_number_per_replica(self) -> dict[str, int]:
        n_cfg_per_replica: dict[str, int] = {}
        for replica in CONFIGURATION_NUMBERS[self.ensemble].keys():
            n_cfg_per_replica[replica] = len(
                self.get_config_list(int(replica[1:]))
            )
        return n_cfg_per_replica

    def get_replicas(self) -> list[tuple[str, int]]:
        return REPLICA_ORDER[self.ensemble]

    def get_replica_slices(self) -> dict:
        r"""
        Returns the np.s_[start:stop:step] slice for each replica of an ensemble.
        For replicas with order == -1 the slice is reversed.
        """
        n_cfg_per_replica = self.get_config_number_per_replica()
        replicas = self.get_replicas()
        slices: dict = {}
        n_cfg_acc: int = 0
        for replica, order in replicas:
            start = n_cfg_acc
            stop = n_cfg_acc + n_cfg_per_replica[replica]
            n_cfg_acc = stop
            if order == -1:
                # Reverse: walk from stop-1 down to start
                slices[replica] = s_[stop - 1 : start - 1 if start > 0 else None : -1]
            else:
                slices[replica] = s_[start:stop:order]
        return slices

    def get_total_config_list(self) -> list[str]:
        r"""
        Returns a list of all configs in analysis order as
            f"{replica}:{n_cfg}"
        e.g. "r000:12" — replica r000, configuration number 12 (1-indexed, not array index).
        Look up the analysis index of a given config via:
            cfg_list.index("r000:12")
        """
        return [
            f"{replica}:{n_cfg}"
            for replica, order in self.get_replicas()
            for n_cfg in self.get_config_list(int(replica[1:]))[::order]
        ]

    def get_replica_order_indices(self) -> list[int]:
        return [order_index for _ , order_index in REPLICA_ORDER[self.ensemble]]

    def get_configuration_numbers_per_replica(self, replica: str) -> int:
        return CONFIGURATION_NUMBERS[self.ensemble][replica]['end_cfg']
        
    def get_total_configuration_number(self) -> int:
        return sum(self.get_config_number_per_replica().values())

    def get_c2pt_source_sets(self) -> list[str]:
        return list(MEASUREMENTS_C2PT[self.ensemble]['c2pt'].keys())

    def get_c3pt_source_sets(self) -> list[str]:
        return list(MEASUREMENTS_C3PT[self.ensemble]['c3pt'].keys())

    def get_c2pt_measurement_ids(self, source_set: str) -> list[int]:
        return list(MEASUREMENTS_C2PT[self.ensemble]['c2pt'][source_set]['src_ids'])

    def get_c3pt_measurement_ids(self, source_set: str, tsep: int) -> list[int]:
        return list(MEASUREMENTS_C3PT[self.ensemble]['c3pt'][source_set][tsep]['src_ids'])

    def get_source_sink_separations(self, source_set: str) -> list[int]:
        return list(MEASUREMENTS_C3PT[self.ensemble]['c3pt'][source_set].keys())

    def get_exceptionals(self) -> np.ndarray:
        entry = EXCEPTIONAL_CONFIGURATION_NUMBERS.get(self.ensemble)
        if entry is None:
            return np.array([], dtype=int)
        combined = sorted(set(entry["rwfs"]) | set(entry["cfgs"]))
        return np.array(combined, dtype=int)   
