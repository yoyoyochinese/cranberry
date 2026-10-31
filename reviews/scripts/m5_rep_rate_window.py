"""M5 review: rep_rate_tail_check with tail_model = source when 1/f_rep lies beyond the source window.

Run: uv run python reviews/scripts/m5_rep_rate_window.py
"""
from cranberry.blocks.chk.rep_rate_tail import make_rep_rate_tail_check
from cranberry.blocks.source.pulsed_testing import BOREALIS

if __name__ == "__main__":
    for win in (12.0, 3.0):
        src = {"block": "pulsed_opo_source", "params": {**BOREALIS, "x_ss": 2.4, "window_kappa": win, "rep_rate_hz": 60e6}}
        for model in ("source", "analytic"):
            p = {"rep_rate_hz": 60e6, "threshold": 1e-3, "tail_model": model, "toggles": {"enabled": True}}
            p.update({"source": src} if model == "source" else {"fwhm_hz": 66e6})
            c = make_rep_rate_tail_check(p)
            print(f"window_kappa={win:4.1f} tail_model={model:8s}: tail={c.result['tail']:.3e} ok={c.result['ok']}")
            if win == 3.0 and model == "source":  # finding M5-08: warning silently missed
                assert c.result["tail"] == 0.0 and c.result["ok"]
