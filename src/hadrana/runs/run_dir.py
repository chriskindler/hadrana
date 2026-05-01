# hadrana/runs/run_dir.py
from __future__ import annotations

import datetime as dt
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

def create_run_directory(
    analysis_path: Path,
    ensemble: str,
    observable: Literal["c2pt", "c2pt_ratio", "ratio", "formfactor"],
    fit_identifier: str
) -> Path:
    """
        INPUT:
            analysis_root: Path
                Root path for fit results /home/ck/phd/nucleon_analysis

            ensemble: str
                Ensemble identifier

            observable: str
                We consider four observables: (p) = finite momentum, (0) = zero momentum
                Two-point functions       c2pt(p)
                Two-point function ratios c2pt(p) / c2pt(0)
                Ratios                    c3pt(p) / c2pt(0)
                Form factors              GA axial, GP pseudoscalar, GtildeP induced pseudoscalar

            label: str
                The UUID of the fit run that specifies the fit type with crucial parameters.
                It will be combined with UTC timestamp and ensemble identifier.

        RETURNS:
            The path to the fit results
            Path(analysis_path / run_id) 

            UTC TIME: https://www.utctime.net/
            run_id is derived from UTC timestamp + ensemble + label, e.g.
            observable = "c2pt"
            model      = "one-state"
            2026-04-30-T22-17-UTC-D251-c2pt-one-state-nsquare01-tmin04-StN05
            
    """

    time_format: str = "%Y-%m-%d-T%H-%M-%Z"
    time_stamp:  str = datetime.now(timezone.utc).strftime(time_format)

    run_id:      str = f"{time_stamp}-{ensemble}-{observable}-{fit_identifier}"

    # intentional failure of directory creation if it already exists to prevent accidental overwrite
    run_dir: Path = Path(analysis_path / "runs" / run_id)
    run_dir.mkdir(parents=True, exist_ok=False)

    return run_dir


if __name__ == "__main__":
    ensemble   = "D251"
    observable = "c2pt"

    analysis_path = Path("/home/ck/phd/nucleon_analysis")

    nsquare = 1
    model = "one-state"
    label = f"c2pt-{model}-nsquare0{nsquare}-tmin04-StN05"
    test = create_run_directory(analysis_path, ensemble, observable, label)


