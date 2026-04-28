# cora/ensembles/specs.py
from __future__ import annotations

from .ensemble_specs import (
    CONFIGURATION_NUMBERS,
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
    
    def get_replicas(self) -> list[str]:
        return [rep for rep , _ in REPLICA_ORDER[self.ensemble]]

    def get_replica_order(self) -> list[tuple[str, int]]:
        return list(REPLICA_ORDER[self.ensemble])

    def get_replica_order_indices(self) -> list[int]:
        return [order_index for _ , order_index in REPLICA_ORDER[self.ensemble]]

    def get_configuration_numbers_per_replica(self, replica: str) -> int:
        return CONFIGURATION_NUMBERS[self.ensemble][replica]['end_cfg']
        
    def get_total_configuration_number(self) -> int:
        return sum(self.get_configuration_numbers_per_replica(rep) for rep in self.get_replicas())

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

    
