"""Phase 5 M1 review: D61-5 follow-up — do coil_birefringence / per_budget accept YAML-loaded dB keys?

``config.load`` converts ``<x>_db`` → ``<x>_lin`` (D56-1). D61-5 says per_budget and coil_birefringence
need the same ``_lin`` support as polarizer_projection / dc_fringe. This builds a minimal spec through
``config.load`` with both blocks offline and reports what happens.

Run: uv run python reviews/scripts/m1_refD_lin_names.py
"""

from cranberry import config

COIL = {"n_turns": 38, "turn_length_m": 0.22, "bend_radius_m": 0.015, "cladding_radius_m": 62.5e-6,
        "bend_coeff": 0.136, "beat_length_m": 4.2e-3, "wavelength_m": 1.5488e-6, "K_T_rad_per_m_K": -2.05,
        "bend_loss_db": 0.5, "orientation": "random", "segment": "c", "applies_to": ["probe"], "source": "x",
        "toggles": {"pol_bend_coil": True, "bend_loss": True}}
PERB = {"elements": [{"id": "cp", "type": "coupler", "per_db": 25.0, "source": "x"}], "loss": None,
        "toggles": {"pol_crosstalk": True}}
SEG = {"id": "c", "kind": "pm", "L_m": 8.4, "n_g": 1.468, "beta2_s2_per_m": 0.0, "beta3_s3_per_m": 0.0,
       "zone": "lab", "n_sub": 1, "spatial": "lumped", "pol_model": "jones", "pol_params": {"J_slot": "J.c"},
       "drivers": {}}


def main():
    for kind, params in (("coil_birefringence", COIL), ("per_budget", PERB)):
        tree = {"run": {"duration_s": 1e-4, "f_sim_hz": 1e6, "chunk_samples": 100, "seed": 1, "seed_mode": "common"},
                "topology": {"segments": [SEG], "beams": [], "observers": [], "locks": [], "ports": [], "coef": []},
                "offline": {"blocks": [{"block": kind, "name": "b", "params": params}]}, "pipeline": []}
        try:
            spec = config.load(tree)
            p = spec.offline["blocks"][0]["params"] if isinstance(spec.offline["blocks"][0], dict) else None
            print(kind, "loaded; params after load:", {k: v for k, v in (p or {}).items() if "db" in k or "lin" in k}
                  or (p or {}).get("elements"))
            from cranberry.runner.run import run
            run(spec)
            print(kind, "ran OK")
        except Exception as e:  # noqa: BLE001
            print(kind, "ERROR:", type(e).__name__, str(e)[:300])


if __name__ == "__main__":
    main()
