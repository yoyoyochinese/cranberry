"""M6 review: do the *_nominal presets get the values they write themselves?

D74-7 records "a file's own keys are overwritten by its own `overlays:`". Every F_*_nominal.yaml applies
points/F_common_nominal.yaml as its own overlay, so any offline.params key written in BOTH files resolves to the
F_common_nominal value. This lists every such key and the resolved value (asserts that F_bor6_nominal's p_pr_avg_w
— the D50 "same-frequency probe <= 1-3 nW for the OPO mode" choice — is lost).
Usage: uv run python reviews/scripts/m6_nominal_override.py
"""
from pathlib import Path

from cranberry import config
from cranberry.config.yaml_io import load_yaml

ROOT = Path(__file__).resolve().parents[2]
common = load_yaml(ROOT / "presets/future/points/F_common_nominal.yaml")["offline"]["params"]
lost = {}
for b in ("F_wg100", "F_wg6", "F_bor6", "F_liu20"):
    own = load_yaml(ROOT / f"presets/future/{b}_nominal.yaml")["offline"]["params"]
    spec = config.load(ROOT / f"presets/future/{b}_nominal.yaml")
    res = spec.offline["params"]
    for k, v in own.items():
        if k.endswith("_db"):   # D56-1: <x>_db -> <x>_lin at load
            k, v = k[:-3] + "_lin", 10 ** (float(v) / 10)
        if abs(float(res[k]) - float(v)) > 1e-12 * abs(float(v)):
            lost[(b, k)] = (v, res[k])
            print(f"{b}_nominal: {k} written {v!r} -> resolved {res[k]!r} (F_common_nominal: {common.get(k)!r})")
    print(f"{b}_nominal: {len(own)} own keys, overlap with F_common_nominal: {sorted(set(own) & set(common))}")
assert lost.get(("F_bor6", "p_pr_avg_w")) == (1.0e-9, 1.0e-6), lost
print("confirmed: F_bor6_nominal runs with P̄ = 1 µW, not the 1 nW it writes")
