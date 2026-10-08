"""
位相安定化シミュレーションの最小プロトタイプ（2026-10-07）

目的:
  Moku + EOM/AOM + 高 QE PD を使った「ファイバ伝送スクイーズド光の位相ロック」で、
  ファイバ長 L に対して残留位相ジッタ θ_rms と測定スクイージング量がどう変わるかを、
  周波数領域の雑音バジェットと時間領域の簡易シミュレーションで見積もる。

注意:
  - ファイバ位相雑音 PSD は「仮置き」のモデル（S_REF_PER_M など）。実測値に置き換えること。
  - Moku の遅延は Liquid Instruments の KB の値（LLB fast, 1 MHz LPF で約 1 µs）を丸めて使う。
  - AOM の音響遅延・EOM の Vπ・ループ構成は仮定。ノートの「仮定一覧」を参照。

依存: numpy, matplotlib のみ（遅延 e^{-sτ} は周波数領域で厳密に扱う）
"""
from __future__ import annotations

import dataclasses as dc
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------- 定数
C0 = 299_792_458.0          # 光速 [m/s]
H = 6.62607015e-34          # プランク定数 [J s]
LAM = 1550e-9               # 波長 [m]
NU = C0 / LAM               # 光周波数 [Hz]
N_G = 1.468                 # SMF-28 の群屈折率（概数）

# ---------------------------------------------------------------- 仮定パラメータ
@dc.dataclass
class Params:
    # ファイバ位相雑音（仮置き）:
    #   S_d(f) = (L/1 m) * S_REF * (F_REF/f)^2 / (1 + (f/F_KNEE)^2)  [rad^2/Hz]
    #   音響帯（〜kHz）までは f^-2、それより上は f^-4 で落ちると仮定
    S_REF_PER_M: float = 1e-4   # 1 m あたり, f = F_REF での値
    F_REF: float = 100.0
    F_KNEE: float = 5e3
    # 検出・誤差信号
    P_PROBE: float = 1e-6       # BHD に届く probe パワー [W]
    ETA_PD: float = 0.99        # PD 量子効率
    CLEARANCE_DB: float = 15.0  # BHD の shot noise / 電子雑音 クリアランス [dB]
    # 遅延
    TAU_MOKU: float = 1.0e-6    # 復調 + PID（LLB fast, LPF 1 MHz 相当）[s]
    TAU_AOM: float = 0.5e-6     # AOM 音響遅延（仮定）[s]
    TAU_EOM: float = 20e-9      # EOM ドライバ等 [s]
    # 安定性の要求
    PM_MIN_DEG: float = 45.0
    GM_MIN_DB: float = 6.0
    # 積分帯域
    F_MIN: float = 1.0
    F_MAX: float = 2e7
    # スクイージング源と損失
    SQZ_GEN_DB: float = 15.0    # 生成時（純粋）スクイージング [dB]
    ETA_SOURCE: float = 0.95    # 共振器/導波路の脱出効率など
    ETA_MM: float = 0.98        # LO とのモード整合（パワー）
    FIBER_DB_PER_KM: float = 0.2
    CONNECTOR_LOSS: float = 0.97  # コネクタ・スプライス合計（仮定）


P = Params()
f = np.logspace(np.log10(P.F_MIN), np.log10(P.F_MAX), 6000)
s = 2j * np.pi * f


# ---------------------------------------------------------------- 物理モデル
def tau_fiber(L):
    """片道の伝搬遅延 [s]"""
    return N_G * L / C0


def S_fiber(L, ff=None):
    """ファイバ位相雑音 PSD（片側, rad^2/Hz）。長さに比例（独立な区間の和）と仮定。"""
    ff = f if ff is None else ff
    return (L / 1.0) * P.S_REF_PER_M * (P.F_REF / ff) ** 2 / (1 + (ff / P.F_KNEE) ** 2)


def S_sensor():
    """ヘテロダイン（AOM でずらした probe と LO のビート）の位相換算 shot noise + 電子雑音"""
    shot = H * NU / (P.ETA_PD * P.P_PROBE)          # heterodyne 復調: hν/(ηP)
    return shot * (1 + 10 ** (-P.CLEARANCE_DB / 10)) * np.ones_like(f)


def loop_shape(kind: str, tau: float, f_c: float):
    """
    kind = 'flat'  : EOM（位相アクチュエータ, 平坦）＋ I・PI（= PI²: 積分を制御器側で用意）
    kind = 'integ' : AOM（周波数アクチュエータ → 位相に対して 1/s）＋ PI
    どちらも UGF 付近で L ∝ 1/s、f_c/4 以下で 1/s² のブーストになる。
    （平坦プラントに PI だけだと高域で |L| が下がらず、遅延で発振する）
    f_c を UGF の目標として |L(f_c)| = 1 に規格化した L(f) を返す。
    """
    wI = 2 * np.pi * f_c / 4.0                  # ブーストの零点を f_c/4 に置く
    C = 1 + wI / s
    if kind == "flat":
        C = C / s
    A = np.ones_like(s) if kind == "flat" else 1 / s
    L0 = C * A * np.exp(-s * tau)
    k = 1 / np.abs(np.interp(f_c, f, np.abs(L0)))
    return k * L0


def margins(L):
    mag = np.abs(L)
    ph = np.unwrap(np.angle(L))
    i_c = np.where(np.diff(np.sign(mag - 1)) != 0)[0]
    if len(i_c) == 0:
        return -np.inf, -np.inf
    pm = 180 + np.degrees(ph[i_c[-1]])
    # UGF より上で位相が -180° を下向きに横切る点でのゲイン余裕
    #（AOM + PI は低域で位相がちょうど -180° になるので、UGF より上だけを見る）
    above = np.arange(len(f)) > i_c[-1]
    i_p = np.where(above[:-1] & (ph[:-1] > -np.pi) & (ph[1:] <= -np.pi))[0]
    gm = np.inf if len(i_p) == 0 else -20 * np.log10(mag[i_p[0]])
    return pm, gm


def residual(L, Sd, Sn, extra=None):
    """閉ループの残留位相 PSD: |S|^2 S_d + |T|^2 S_n (+ 抑圧できない項)"""
    S = 1 / (1 + L)
    T = L / (1 + L)
    out = np.abs(S) ** 2 * Sd + np.abs(T) ** 2 * Sn
    if extra is not None:
        out = out + np.abs(T) ** 2 * extra
    return out


def rms(psd):
    return np.sqrt(np.trapezoid(psd, f))


def best_loop(kind, tau, Sd, Sn, extra=None):
    """UGF を走査し、PM/GM 条件を満たす中で θ_rms 最小の設定を返す（BO の代わりの全探索）"""
    best = None
    for f_c in np.logspace(2, 6.5, 400):
        L = loop_shape(kind, tau, f_c)
        pm, gm = margins(L)
        if pm < P.PM_MIN_DEG or gm < P.GM_MIN_DB:
            continue
        th = rms(residual(L, Sd, Sn, extra))
        if best is None or th < best[0]:
            best = (th, f_c, pm, gm, L)
    return best


# ---------------------------------------------------------------- スクイージング
def eta_total(L):
    eta_fiber = 10 ** (-P.FIBER_DB_PER_KM * L / 1000 / 10) * P.CONNECTOR_LOSS
    return P.ETA_SOURCE * P.ETA_MM * P.ETA_PD * eta_fiber


def measured_sqz_db(L, theta_rms):
    """V± = η e^{±2r} + 1 - η（shot noise 単位）, ガウス位相揺らぎで平均"""
    r = P.SQZ_GEN_DB / (20 * np.log10(np.e))       # 10log10(e^{2r}) = dB
    eta = eta_total(L)
    Vm = eta * np.exp(-2 * r) + 1 - eta
    Vp = eta * np.exp(+2 * r) + 1 - eta
    V = Vm + (Vp - Vm) * (1 - np.exp(-2 * theta_rms ** 2)) / 2
    return 10 * np.log10(V), 10 * np.log10(Vm)


# ---------------------------------------------------------------- 構成の比較
CONFIGS = {
    # 受信側ロック: probe はスクイーズド光と同じファイバを伝搬し、受信側の LO を EOM で動かす
    "A-EOM: 受信側ロック (EOM)": dict(kind="flat", tau=lambda L: P.TAU_MOKU + P.TAU_EOM, fnc=False),
    "A-AOM: 受信側ロック (AOM)": dict(kind="integ", tau=lambda L: P.TAU_MOKU + P.TAU_AOM, fnc=False),
    # 送信側 FNC: 往復光で検出し、入口の AOM で補正（Ma 1994 型）。ループに往復遅延が入る
    "B-FNC: 送信側 FNC (AOM, 往復)": dict(kind="integ", tau=lambda L: P.TAU_MOKU + P.TAU_AOM + 2 * tau_fiber(L), fnc=True),
}


def evaluate(L):
    Sd = S_fiber(L)
    Sn = S_sensor()
    rows = {}
    for name, c in CONFIGS.items():
        extra = None
        if c["fnc"]:
            # Williams, Swann & Newbury (2008): 遅延で抑圧できない片道残差 ≈ (1/3)(2π f τ)^2 S_fiber
            extra = (1 / 3) * (2 * np.pi * f * tau_fiber(L)) ** 2 * Sd
        b = best_loop(c["kind"], c["tau"](L), Sd, Sn, extra)
        rows[name] = b
    return rows


def main(outdir="."):
    lengths = np.logspace(0, 4, 17)    # 1 m – 10 km
    res = {k: [] for k in CONFIGS}
    for L in lengths:
        rows = evaluate(L)
        for k, b in rows.items():
            res[k].append((b[0], b[1]))

    # ---- 表示用テキスト
    print(f"hν/(ηP) sensor floor = {S_sensor()[0]:.2e} rad^2/Hz  (P_probe={P.P_PROBE*1e6:.1f} µW)")
    print(f"{'L [m]':>8} | " + " | ".join(f"{k.split(':')[0]:>18}" for k in CONFIGS))
    for i, L in enumerate(lengths):
        cells = []
        for k in CONFIGS:
            th, fc = res[k][i]
            sq, sq0 = measured_sqz_db(L, th)
            cells.append(f"{th*1e3:6.1f} mrad {sq:6.2f} dB")
        print(f"{L:8.0f} | " + " | ".join(cells) + f" | (θ=0: {measured_sqz_db(L,0)[0]:.2f} dB)")

    # ---- 図 1: θ_rms と測定スクイージング vs L
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for k in CONFIGS:
        th = np.array([r[0] for r in res[k]])
        ax[0].loglog(lengths, th * 1e3, "o-", ms=3, label=k)
        # θ_rms > 0.3 rad は小角近似もロックの前提も崩れるので表示しない
        sq = np.array([measured_sqz_db(L, t)[0] for L, t in zip(lengths, th)])
        sq[th > 0.3] = np.nan
        ax[1].semilogx(lengths, sq, "o-", ms=3, label=k)
    ax[1].semilogx(lengths, [measured_sqz_db(L, 0)[0] for L in lengths], "k--", lw=1, label="位相揺らぎなし（損失のみ）")
    ax[0].set_xlabel("ファイバ長 L [m]"); ax[0].set_ylabel("残留位相 θ_rms [mrad]")
    ax[1].set_xlabel("ファイバ長 L [m]"); ax[1].set_ylabel("測定スクイージング [dB]（θ_rms ≤ 0.3 rad のみ）")
    for a in ax:
        a.grid(True, which="both", alpha=0.3)
    ax[0].legend(fontsize=7); ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(f"{outdir}/2026-10-07_budget_vs_length.png", dpi=150)

    # ---- 図 2: L = 1 km での PSD 内訳（A-EOM と B-FNC）
    L = 1000.0
    rows = evaluate(L)
    Sd = S_fiber(L)
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.loglog(f, Sd, "k", lw=1, label="自由走行 S_d (1 km, 仮置き)")
    ax.loglog(f, S_sensor(), "gray", ls=":", label="センサ雑音 hν/(ηP)")
    for k in ["A-EOM: 受信側ロック (EOM)", "B-FNC: 送信側 FNC (AOM, 往復)"]:
        th, fc, pm, gm, Lf = rows[k]
        extra = (1 / 3) * (2 * np.pi * f * tau_fiber(L)) ** 2 * Sd if CONFIGS[k]["fnc"] else None
        ax.loglog(f, residual(Lf, Sd, S_sensor(), extra), label=f"{k}  UGF={fc/1e3:.0f} kHz, θ={th*1e3:.1f} mrad")
    ax.set_ylim(1e-16, 1e2); ax.set_xlabel("f [Hz]"); ax.set_ylabel("位相 PSD [rad²/Hz]")
    ax.grid(True, which="both", alpha=0.3); ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(f"{outdir}/2026-10-07_psd_1km.png", dpi=150)
    return res


# ---------------------------------------------------------------- 時間領域: EOM のレンジ飽和とオフロード
def run_td(phi_d, dt, d_samp, kp, ki, k_off, use_aom, eom_range=np.pi):
    """離散時間の受信側ロック。誤差 = sin(Δφ)、d_samp サンプルの遅延つき。"""
    n = len(phi_d)
    buf = np.zeros(d_samp + 1)
    integ = 0.0
    phi_eom = 0.0
    phi_aom = 0.0
    res = np.empty(n)
    eom_hist = np.empty(n)
    for i in range(n):
        dphi = phi_d[i] - phi_eom - phi_aom
        res[i] = dphi
        buf[i % (d_samp + 1)] = np.sin(dphi)
        e = buf[(i - d_samp) % (d_samp + 1)] if i >= d_samp else 0.0
        integ = np.clip(integ + ki * e * dt, -eom_range, eom_range)   # anti-windup
        phi_eom = np.clip(kp * e + integ, -eom_range, eom_range)
        if use_aom:
            phi_aom += k_off * phi_eom * dt      # EOM のずれを AOM 周波数（位相の積分）で肩代わり
        eom_hist[i] = phi_eom
    return res, eom_hist


def td_prediction(kp, ki, d_samp, dt, psd):
    """同じ離散ループの周波数領域予測（ZOH なし, 小角近似）: θ_rms^2 = ∫ |S(e^{jωdt})|^2 S_d df"""
    ff = np.logspace(0, np.log10(0.5 / dt), 4000)
    z = np.exp(1j * 2 * np.pi * ff * dt)
    C = kp + ki * dt / (1 - 1 / z)
    Lz = C * z ** (-d_samp - 1)          # 計算順序による 1 サンプルの遅れを含める
    return np.sqrt(np.trapezoid(np.abs(1 / (1 + Lz)) ** 2 * psd(ff), ff))


def time_domain_demo(outdir=".", T=0.25, dt=1e-6, seed=1, L_m=100.0):
    rng = np.random.default_rng(seed)
    n = int(T / dt)
    t = np.arange(n) * dt
    psd = lambda ff: S_fiber(L_m, ff)
    # FFT 合成（周期境界）。分散を PSD 積分で検算する
    fk = np.fft.rfftfreq(n, dt)
    S = np.zeros_like(fk); S[1:] = psd(fk[1:])
    X = np.sqrt(S / 4 / dt * n) * (rng.standard_normal(len(fk)) + 1j * rng.standard_normal(len(fk)))
    X[0] = 0
    acoustic = np.fft.irfft(X, n)
    var_target = np.trapezoid(S[1:], fk[1:])
    print(f"合成雑音の rms = {acoustic.std():.3f} rad（PSD 積分からの期待値 {np.sqrt(var_target):.3f} rad）")
    drift = 30.0 * t                     # 熱ドリフト（仮置き, 30 rad/s）
    d_samp = int(round(P.TAU_MOKU / dt))
    kp, ki, k_off = 0.3, 2 * np.pi * 30e3, 2 * np.pi * 300.0

    # (1) 検算: ドリフトなし・レンジ制限なし → 周波数領域の予測と一致するか
    res0, _ = run_td(acoustic, dt, d_samp, kp, ki, 0.0, False, eom_range=1e9)
    pred = td_prediction(kp, ki, d_samp, dt, psd)
    tail = slice(n // 5, None)
    print(f"[検算] 時間領域 θ_rms = {res0[tail].std()*1e3:.2f} mrad / 周波数領域予測 = {pred*1e3:.2f} mrad")

    # (2) ドリフトあり: EOM のみ vs EOM + AOM オフロード
    out = {}
    for case, use_aom in [("EOM のみ", False), ("EOM + AOM オフロード", True)]:
        res, eom = run_td(drift + acoustic, dt, d_samp, kp, ki, k_off, use_aom)
        res_w = np.angle(np.exp(1j * res))
        slips = int(np.round(np.abs(res[-1] - res_w[-1]) / (2 * np.pi)))
        out[case] = (res_w, eom)
        print(f"{case:>20}: θ_rms（±π に畳んだ値, 後半）= {res_w[n//2:].std()*1e3:8.1f} mrad, 2π スリップ ≈ {slips} 回")

    fig, ax = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
    for case, (res_w, eom) in out.items():
        ax[0].plot(t * 1e3, res_w, lw=0.6, label=case)
        ax[1].plot(t * 1e3, eom, lw=0.6, label=case)
    ax[0].set_ylabel("残留位相 [rad]（±π に畳む）")
    ax[1].set_ylabel("EOM の位相 [rad]"); ax[1].set_xlabel("t [ms]")
    ax[0].set_title(f"受信側ロック, L = {L_m:.0f} m 相当の音響雑音 + 30 rad/s ドリフト（仮置き）", fontsize=9)
    for a in ax:
        a.grid(alpha=0.3); a.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(f"{outdir}/2026-10-07_range_offload.png", dpi=150)


if __name__ == "__main__":
    import sys
    od = sys.argv[1] if len(sys.argv) > 1 else "."
    plt.rcParams["font.family"] = ["Hiragino Sans", "Noto Sans CJK JP", "IPAexGothic", "DejaVu Sans"]
    main(od)
    time_domain_demo(od)
