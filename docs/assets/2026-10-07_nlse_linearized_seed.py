"""
線形化 NLSE ＋ 共分散行列によるファイバ中の量子雑音計算 — 最小の種コード（2026-10-07）

計画ノート [[2026-10-07 phase-stab-sim_implementation-plan]] の M2–M4 の「動く最小例」。
  1. RK4IP で NLSE（分散 + Kerr）を伝搬
  2. jax.jacfwd で古典場まわりの線形化行列 S（= 揺らぎの伝達行列）を厳密に得る
  3. 真空入力の共分散 V_in = I（真空単位）→ V_out = S S^T
  4. 検証:
     (a) 分散なし・平坦パルスで、Kerr スクイージングの解析式 1+2Φ²−2Φ√(1+Φ²) と一致
     (b) 損失なしでは S がシンプレクティック（S^T Ω S = Ω）
     (c) 1 次ソリトンの形が保たれる
     (d) Bloch–Messiah（S の特異値）から最良のスクイージングを得る

単位: 時間 ps, 長さ m, パワー W, β2 [ps²/m], γ [1/(W m)]
NLSE: ∂A/∂z = -i(β2/2)∂²A/∂T² + iγ|A|²A （Agrawal の符号）
FFT: numpy/jax の fft（順方向 e^{-iωt}）。∂/∂T → -iω なので ∂²/∂T² → -ω²
"""
import jax
import jax.numpy as jnp
import numpy as np

jax.config.update("jax_enable_x64", True)


# ------------------------------------------------------------------ 伝搬
def make_propagator(N, dt, beta2, gamma, L, nsteps):
    w = 2 * np.pi * np.fft.fftfreq(N, dt)                 # 角周波数 [rad/ps]
    Dw = 1j * (beta2 / 2) * w**2                          # 線形演算子（周波数領域）
    h = L / nsteps
    half = jnp.exp(Dw * h / 2)

    def N_op(A):                                          # 非線形項（時間領域）
        return 1j * gamma * jnp.abs(A) ** 2 * A

    def rk4ip_step(A, _):
        # Hult (2007) の RK4IP: 相互作用描像で RK4
        AI = jnp.fft.ifft(half * jnp.fft.fft(A))
        k1 = jnp.fft.ifft(half * jnp.fft.fft(h * N_op(A)))
        k2 = h * N_op(AI + k1 / 2)
        k3 = h * N_op(AI + k2 / 2)
        k4 = h * N_op(jnp.fft.ifft(half * jnp.fft.fft(AI + k3)))
        A_new = jnp.fft.ifft(half * jnp.fft.fft(AI + k1 / 6 + k2 / 3 + k3 / 3)) + k4 / 6
        return A_new, None

    @jax.jit
    def propagate(A0):
        A, _ = jax.lax.scan(rk4ip_step, A0, None, length=nsteps)
        return A

    return propagate


def real_map(propagate):
    """複素場の写像を実ベクトル (Re, Im) の写像にする（jacfwd 用）"""
    def f(x):
        n = x.shape[0] // 2
        A = propagate(x[:n] + 1j * x[n:])
        return jnp.concatenate([A.real, A.imag])
    return f


# ------------------------------------------------------------------ 測定
def homodyne_vector(u, theta, dt):
    """時間モード u（Σ|u|²dt = 1）, LO 位相 θ の quadrature を、真空単位の (Re,Im) 座標で表すベクトル w
    揺らぎ δa_j = (y_re + i y_im)/sqrt(4 dt) とおくと、真空の y は共分散 I。X = w·y, |w| = 1"""
    c = jnp.exp(-1j * theta) * jnp.conj(u) * 2 * dt / jnp.sqrt(4 * dt)
    # X = Re(Σ c_j (y_re + i y_im)) = Σ Re(c) y_re - Im(c) y_im
    return jnp.concatenate([c.real, -c.imag])


def scan_theta(V, u, dt, n=2001):
    th = np.linspace(-np.pi / 2, np.pi / 2, n)
    vals = np.array([float(homodyne_vector(u, t, dt) @ V @ homodyne_vector(u, t, dt)) for t in th])
    return vals.min(), vals.max()


def db(x):
    return 10 * np.log10(x)


# ------------------------------------------------------------------ 検証
def test_kerr_cw_limit():
    N, T = 512, 400.0                      # 時間窓 400 ps
    dt = T / N
    t = (np.arange(N) - N / 2) * dt
    P0, gamma, L = 1.0, 1.3e-3, 1000.0     # 1 W, γ = 1.3 /W/km, 1 km → Φ = 1.3 rad
    A0 = np.sqrt(P0) * np.exp(-(t / 120) ** 20)            # 平坦な超ガウス
    prop = make_propagator(N, dt, beta2=0.0, gamma=gamma, L=L, nsteps=200)
    x0 = jnp.concatenate([jnp.array(A0.real), jnp.array(A0.imag)])
    S = jax.jacfwd(real_map(prop))(x0)
    V = S @ S.T
    A_out = prop(jnp.array(A0, dtype=complex))
    top = np.abs(t) < 40
    u = np.where(top, np.exp(1j * np.angle(np.array(A_out))), 0)
    u = u / np.sqrt(np.sum(np.abs(u) ** 2) * dt)
    vmin, vmax = scan_theta(V, jnp.array(u), dt)
    Phi = gamma * P0 * L
    ana_min = 1 + 2 * Phi**2 - 2 * Phi * np.sqrt(1 + Phi**2)
    ana_max = 1 + 2 * Phi**2 + 2 * Phi * np.sqrt(1 + Phi**2)
    print(f"(a) Kerr Φ={Phi:.2f}: 数値 {db(vmin):.3f}/{db(vmax):.3f} dB, 解析 {db(ana_min):.3f}/{db(ana_max):.3f} dB")
    return abs(db(vmin) - db(ana_min)) < 0.05


def test_symplectic_and_soliton():
    N, T = 256, 40.0
    dt = T / N
    t = (np.arange(N) - N / 2) * dt
    beta2, gamma = -22e-3, 1.3e-3          # SMF: -22 ps²/km, 1.3 /W/km（m 単位に換算）
    T0 = 1.0
    P0 = abs(beta2) / (gamma * T0**2)      # 1 次ソリトン
    LD = T0**2 / abs(beta2)
    L = 2 * LD
    A0 = np.sqrt(P0) / np.cosh(t / T0)
    prop = make_propagator(N, dt, beta2, gamma, L, nsteps=2000)
    A_out = np.array(prop(jnp.array(A0, dtype=complex)))
    shape_err = np.max(np.abs(np.abs(A_out) - np.abs(A0))) / np.sqrt(P0)
    print(f"(c) 1 次ソリトン: L = {L:.1f} m 伝搬後の |A| の最大相対誤差 = {shape_err:.2e}")

    x0 = jnp.concatenate([jnp.array(A0), jnp.zeros(N)])
    S = np.array(jax.jacfwd(real_map(prop))(x0))
    Om = np.block([[np.zeros((N, N)), np.eye(N)], [-np.eye(N), np.zeros((N, N))]])
    sym_err = np.max(np.abs(S.T @ Om @ S - Om))
    print(f"(b) シンプレクティック誤差 max|SᵀΩS − Ω| = {sym_err:.2e}")

    # (d) Bloch–Messiah: S の特異値 = 各スーパーモードのスクイージング因子
    sv = np.linalg.svd(S, compute_uv=False)
    print(f"(d) 最良スーパーモード: {db(sv.min()**2):.2f} dB（反スクイージング {db(sv.max()**2):.2f} dB）")
    # 参考: LO = 出力パルス形状のときの最良値
    u = A_out / np.sqrt(np.sum(np.abs(A_out) ** 2) * dt)
    vmin, vmax = scan_theta(S @ S.T, jnp.array(u), dt, n=721)
    print(f"    LO = 出力パルス形状: {db(vmin):.2f} dB / {db(vmax):.2f} dB")
    return sym_err < 1e-8 and shape_err < 1e-3


if __name__ == "__main__":
    ok_a = test_kerr_cw_limit()
    ok_b = test_symplectic_and_soliton()
    print("PASS" if (ok_a and ok_b) else "FAIL")
