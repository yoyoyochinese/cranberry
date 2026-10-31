# T17 パルス励起 OPO 光源（pulsed-opo-source）

> **訂正（2026-11-01、DECISIONS D60-5・D70-2）:** x_ss = 2.5 で r₀ = 1.038（1.029 ではない）、Gaussian の LO の重なり 0.955、矩形 0.866、R± は Borealis −7.72 / +8.82 dB、Liu −8.95 / +9.22 dB。


- Phase: 1（追加の deep dive、D43）／ 日付: 2026-10-31
- 担当: Borealis 型（Madsen et al. 2022）と Liu 2026 型の「パルスで励起するしきい値以下の PPKTP ヘミリシック OPO」。
  時間モード、FSR モードの除去、OPO の長さロック、ポンプの安定化と位相ロック。T02・T03・T06・T07・T10・T11・T12 の結論がどう変わるか。
- 入力: `docs/squeezing_sim_prompt.md`、`AGENTS.md`、`DECISIONS.md`（D01–D49）、`notes/conventions.md`（2026-10-31 改訂）、`notes/survey.md`、
  T02・T03・T06・T07・T10・T11・T12 のノート、`refs/xanadu_2022.pdf`（**Supplementary Information だけ、37 ページ。本文は入っていない**）、
  `refs/Liu et al. - 2026 - … w Supplementary Information.pdf`（SI、64 ページ）と同じ論文の本文 PDF。
- 値の印: **[M-SI p.n]** = Madsen 2022 SI の該当ページ（PDF のページ下の番号）、**[M-SI Fig. S7]** = 図の中の表（PDF を画像にして読んだ。pdftotext では表が抜ける）、
  **[L-SI §n]**・**[L-Methods]**・**[L-本文]** = Liu 2026、**[AE]** = この担当の概算、**[sim]** = この担当の数値計算（scratchpad の使い捨てスクリプト、コミットしない。式と入力は §2・§7 にすべて書いた）、
  **[UNVERIFIED]** = 原典で確認していない値・書誌。

## 0. 結論（先に）

1. **時間モードは κ_c と T_p で決まり、Borealis と Liu はどちらも T_p κ_c ≈ 0.62–0.63 に設計されている。**
   Borealis: 線幅 66 MHz FWHM（κ_c = 2.07×10⁸ s⁻¹、1/κ_c = 4.8 ns）、ポンプ 3 ns 矩形 [M-SI p.15, Fig. S7]。
   Liu: FSR 7.4 GHz・フィネス 59 → 125 MHz FWHM [AE]（κ_c = 3.94×10⁸ s⁻¹、1/κ_c = 2.5 ns）、ポンプ 1.6 ns 矩形 [L-SI §1.1]。
   出てくる時間モードは「立ち上がり ~T_p、減衰 e^{−κ_c t}」の形で、実効モード時間 τ_eff = (∫u dt)² ≈ 3.0/κ_c（**Borealis 14.5 ns、Liu 7.5 ns**）[sim]。
   D07 の 1 ps（τ_eff = 2.1 ps）より 3.5×10³–6.8×10³ 倍長い。
2. **パルス励起 OPO の量子 Langevin モデル（§2）は、文献の数値を合わせ込みなしで再現する。** Borealis の値（T_p κ_c = 0.622、κ_p/κ_c = 1.924）で r = 0.533 のとき Schmidt 数 K = 1.10（Borealis の測定 K ≈ 1.1、g⁽²⁾ = 1.91）、
   Liu の値で r = 1.8 のとき K = 1.026（純度 0.975、Liu の「r = 1.8 で識別不能性 97%」）[sim]。K は r とともに下がる（Liu SI §1.2 の「ポンプを強くすると純度が上がる」と同じ向き）。
   T_p → ∞ では CW の T02 の結果（r → 2 artanh x、離調の回転 −(Δ/κ_c)/(1 + x²)）に戻る [sim]。
3. **ns の時間モードでは、T03・T06・T11 のタイミング・分散の要求はほぼ消える。** 重なり 99% の時間ずれの許容は 0.090/κ_c = **436 ps（Borealis）／227 ps（Liu）**（ps パルスでは 8.5 fs–0.85 ps）。
   分散長 L_D ~ 10⁸ m、熱によるタイミングドリフト（63 m で 2.5 ps/K）は 99% の許容に対して約 170 K 分の余裕。Borealis の「250 ps の精度で合わせた」[M-SI p.18] と整合。
   代わりに **新しい制約が 2 つ出る**: (a) モードの裾 ≈ 8.8 e^{−2κ_c T}（T はポンプの立ち上がりから）が次のスロットに入るので、**f_rep ≲ 45 MHz（Borealis 型）／≲ 87 MHz（Liu 型）**（裾 10⁻³）。
   Borealis 型の OPO を 100 MHz で使うと 13.8% が次のスロットに入る（悪い構成として残す）。(b) LO の時間波形: ガウス形の LO では重なりが最大 0.951、矩形のゲートでは 0.864（§6.3）。
4. **散乱雑音の「モードあたりの光子数」は、狭帯域の雑音（Rayleigh、GAWBS のピーク）についてだけ τ_eff に比例して増え、広帯域の雑音（Raman）は変わらない。**
   同じ周波数・同じ偏光の逆向き probe 1 mW・63 m の Rayleigh は、1 ps モードの 0.03–0.06 から **435（Borealis）／226（Liu）光子/モード** になる [AE]。
   −10 dB を 0.05 dB 以内に保つには P̄_pr ε_pol ≲ 1.3 nW（63 m）。T10 の積の式では、ロックの位相誤差の下限が σ_lock ≥ 1.1 mrad（63 m、10 kHz）／3.5 mrad（630 m）になる。
   一方、**モードの帯域が 100 MHz 級に狭いので、周波数をずらした probe（Borealis は +7 GHz）と Brillouin（10.8 GHz）はスペクトルで分けられる**（ps パルスではできなかった）。
5. **ロックの役割**（D42）:
   - OPO の長さロック（PDH）は squeezed 光の位相を直接は決めない。決めるのはロック光のモードの共振で、squeezed 光の共振とのずれ（離調 Δ_s）は**分散・複屈折の温度変化**で残る。
     Borealis はロック光を「+7 GHz の側帯を SHG した 775 nm 付近、ポンプと直交偏光、CW」にしている [M-SI p.15–16]。Liu はロック光を「ポンプと同じ偏光・同じ周波数の 1 µs パルス、時間をずらす」にしている [L-SI §1.1]。
   - スクイーズ角は**キャビティ内のポンプ位相**で決まる: θ_sq = φ_p,ic/2 + π/2 + δθ_sig、φ_p,ic = φ_p,in − arctan(Δ_p/κ_p)（定常）。
   - 離調の効き方（Borealis の値、パルス励起）[sim]: δθ_sq ≈ −0.59 – −0.66 × Δ_s/κ_c（信号の離調だけ）、長さの変化（Δ_p = 2Δ_s）では ≈ −1.05 – −1.09 × Δ_s/κ_c（x_ss = 0.5–2.4）。
     **Liu 型の「透過したポンプを参照光とホモダインしてポンプ光路の PZT を動かす」ロック [L-Methods] はポンプ側の項だけを消し、−0.53 – −0.57 × Δ_s/κ_c が残る。**
     これは T02 の「同軸の probe は離調による回転を完全には追わない」と同じ種類の残差で、Δ_s/κ_c = 2F δl/λ0 から **約 40 mrad/nm（往復の光路長）**（Borealis・Liu とも）。
   - 単一の光源を時間で分ける構成では、光源のスクイーズ角の揺らぎは指標 (ii) には 4 sin²(πf(τ_k − τ_0)) の重みでしか入らない（遅いドリフトは共通）。指標 (i)（LO との相対位相）にはそのまま入る。
6. **Borealis の SI にはポンプ位相のロックの記述がない。** ポンプパワーは OPO を透過したポンプ（PD2）を見て AOM で安定化している [M-SI p.15]。透過パワーはキャビティ内のパワーなので、このループは Δ_p のドリフトによる x の変化も補正する [AE]。

## 1. 記号（symbols, English）

conventions.md の記号（t, Ω, ω0, λ0, f, κ_c, Δν_c, x, θ_sq, φ_p, τ_eff, N）をそのまま使う。追加分:

- `κ_c`: signal-mode field decay rate (HWHM) [rad/s] (conventions §4). `Δν_c = κ_c/π`: FWHM [Hz]. `κ_T = ρ κ_c`: decay through the output coupler. `κ_L = κ_c − κ_T`.
- `ρ = κ_T/κ_c`: escape efficiency. `T_oc`: output-coupler power transmittance at λ_s. `L_ic`: other round-trip loss (HR leak, AR, absorption).
- `FSR = c/l` [Hz], `l`: round-trip optical path length [m]. `F = FSR/Δν_c`: finesse. Then `κ_c = π c/(l F)`.
- `κ_p`: pump-mode field decay rate (HWHM) [rad/s] of the doubly resonant cavity; `None` for a single-resonant (pump non-resonant) cavity. `Δ_p`: pump-mode detuning [rad/s].
- `Δ_s` (= `Δ` of T02): signal-mode detuning `ω_cav − ω0` [rad/s]; `d = Δ_s/κ_c`.
- `s(t)`: normalized input pump envelope (rectangular, 1 for 0 ≤ t < T_p). `T_p`: pump pulse duration [s]. `f_rep = 1/T_rep`: pulse repetition rate.
- `b(t)`: normalized intracavity pump amplitude, `db/dt = −(κ_p + iΔ_p) b + κ_p s(t)` (steady state `b = 1` at Δ_p = 0; `b = s` if non-resonant).
- `ε(t) = κ_c x_ss b(t) e^{iφ_p,in}`: parametric coupling rate [rad/s] (= `κ|β_p,ic|` in the notation of conventions §4 H_int).
- `x_ss = √(P_pk/P_th,CW)`: normalized pump amplitude that the pulse's peak power would give in CW steady state (`x` of conventions §4 for a CW pump). `P_pk`: peak input pump power, `P_th,CW`: CW threshold.
- `q̂(t), p̂(t)`: intracavity quadratures; `X = (q, p)ᵀ`; `X_in`, `X_L`: vacuum inputs through the coupler and the internal loss.
- `M(t)`: 2×2 drift matrix (§2.1). `G(t, t′)`: its time-ordered propagator (`G = 0` for t < t′).
- `C(t, t′)`: shot-noise-normalized 2×2 covariance kernel of the output field (vacuum = `δ(t − t′) I`); its discretization on bins of width δt is the 2N×2N matrix `C` (vacuum = identity).
- `u(t)`: real temporal mode, `∫u² dt = 1`; `w = (w_q, w_p)`: general measured mode (complex LO waveform). `R_−`, `R_+`: min / max of `wᵀ C w` (squeezed / anti-squeezed variance normalized to shot noise).
- `r_n`: supermode squeeze parameters (from eigenvalues `e^{∓2r_n}` of `C` at ρ = 1). `r_0 = max r_n`. `K = (Σ sinh² r_n)²/Σ sinh⁴ r_n`: Schmidt number (= `1/(g⁽²⁾ − 1)` of the thermal state used by Borealis; `1/K` = purity `P` of Liu SI Eq. S10).
- `τ_eff = (∫u dt)²` (conventions §2) of the optimal mode. `ũ(f) = ∫u(t) e^{+i2πft} dt` (conventions §8 sign), `∫|ũ|² df = 1`, `|ũ(0)|² = τ_eff`.
- `ε_tail(T)`: energy of the optimal mode beyond time T after pump turn-on.
- `r_m`: squeeze parameter of the non-degenerate pair at `ω0 ± m·FSR` (m ≥ 1); `T_f(ν)`: power transmission of the FSR filter.
- `δl`: change of round-trip optical path length [m] (`δl > 0` → `Δ_s < 0`, T02). `T_c`: crystal temperature. `L_c`: crystal length.
- `Δ(dn/dT)`: difference of thermo-optic coefficients between the lock-beam mode and the signal mode (§2.6) [1/K].
- `Φ_st` [1/s]: stray photon flux into the squeezed spatial/polarization mode (T10). `n_mode`: stray photons in the squeezed temporal mode. `N_stray`: per detection gate.

**文献の記号との変換**:
- Borealis・Liu は線幅を FWHM [Hz] で書く → κ_c = π Δν_c。Liu はフィネスと FSR だけを書く → Δν_c = FSR/F。
- Liu SI Eq. (S1)–(S6) の r, t は鏡の振幅反射・透過、γ は非線形係数、β_k(ω) はポンプのスペクトル。このノートのモデルは、その「キャビティ内を伝搬する場＋鏡の境界条件」を、フィネス ≫ 1 の単一縦モード（平均場）に縮めたもの。
  対応は ε(t) ↔ γ L_c × （キャビティ内ポンプ振幅）× c/l の程度で、数値の係数は照合していない（x_ss を直接パラメータにするので不要）。
- Borealis の K = 1/(g⁽²⁾_min − 1) [M-SI Fig. S12]、Liu の純度 P = (n̄² + m̄²)/(n̄ + m̄)² [L-SI Eq. S10]（2 モードのモデル）。どちらも多モードの熱的状態の g⁽²⁾ − 1 = Σn_i²/(Σn_i)² = 1/K と同じ量として比べる。
- φ_OPA,0（conventions §4）: このノートの Langevin 方程式では 0（Δ = 0、φ_p,in = 0 で p̂ が squeezed、θ_sq = π/2）。共振するポンプでは −arctan(Δ_p/κ_p)/2 が定常のオフセットとして加わる（§2.5）。

## 2. 支配方程式

### 2.1 パルス励起の量子 Langevin 方程式

T02 §2.1 の式で、ε を時間の関数にし、ポンプも（二重共振なら）キャビティの方程式に従わせる:

  db/dt = −(κ_p + iΔ_p) b + κ_p s(t)
  dâ/dt = −(κ_c + iΔ_s) â + ε(t) â† + √(2κ_T) â_in + √(2κ_L) b̂_in,  ε(t) = κ_c x_ss b(t) e^{iφ_p,in}
  â_out = √(2κ_T) â − â_in

ポンプの枯渇なし（squeezed vacuum のパルスあたりの光子数 ≲ 3、ポンプは nJ）、単一縦モード（ω0 の近くの (0,0) の帯域）、good-cavity（T_oc = 0.1）。
直交位相 X = (q, p)ᵀ で書くと、ε = ε_r + iε_i として（conventions §2 の q, p。T02 の M(Ω) と同じ符号）

  dX/dt = M(t) X + √(2κ_T) X_in + √(2κ_L) X_L,  M(t) = [[−κ_c + ε_r, Δ_s + ε_i], [−Δ_s + ε_i, −κ_c − ε_r]]

### 2.2 出力の共分散（時間領域）

時間順序の伝搬子 G(t, t′)（dG/dt = M G、G(t′, t′) = I）を使うと、出力の共分散の核は（ショット雑音で規格化、真空 = δ I）

  C(t, t′) = δ(t − t′) I − 2κ_T [G(t, t′) + G(t′, t)ᵀ] + 4κ_T κ_c ∫ G(t, s) G(t′, s)ᵀ ds + 2κ_T G(t, t_0) G(t′, t_0)ᵀ

最後の項は、時刻 t_0（計算の窓の始まり）にキャビティの中にあった真空が出てくる分。**これを落とすと、ρ = 1 で固有値 0（ρ < 1 で 1 − ρ）の偽の「完全な squeezing」が出る**（試作で実際に起きた。§7 の V-P3）。
数値では幅 δt のビンに分け、G を区間ごとの行列指数の積で作る（2N×2N の行列 C）。δt を半分にすると誤差は 1/4（§7 V-P1）。

- 任意の LO 波形（振幅と位相の時間変化）で測る直交位相は w = (w_q, w_p)、Σw² = 1 で、分散は wᵀCw。**最良の squeezing R_− は C の最小固有値**、その固有ベクトルが最適な LO のモード。
- ρ = 1 では状態は純粋で、C の固有値は e^{±2r_n} の対になる（Bloch–Messiah、T03 §2.2）。r_n から K と平均光子数 Σ sinh² r_n = (tr C − 2N)/4 が出る。
- Δ_s = 0、φ_p,in = 0、ε 実数では q と p が分かれ、G はスカラー h_∓(t, t′) = exp[−∫_{t′}^{t} (κ_c ∓ ε(s)) ds]。

### 2.3 CW への縮退（T02 との一致）

ε = κ_c x（一定）、Δ_s = 0 では h_∓ = e^{−λ_∓|t − t′|}、λ_∓ = κ_c(1 ∓ x) で

  C_q,p(τ) = δ(τ) + 2κ_T (κ_c/λ_∓ − 1) e^{−λ_∓|τ|}

これを Fourier 変換（conventions §8）すると S(Ω) = 1 ± ρ·4x/((1 ∓ x)² + (Ω/κ_c)²)、すなわち T02 §2.2 の Ref B Eq. (2)（2π を入れた形）そのもの。
長さ T の平らな窓のモードでは R = 1 + 2κ_T(κ_c/λ − 1)(2/λ)[1 − (1 − e^{−λT})/(λT)]（V-P1 の解析値）。
パルス励起の計算は T_p → ∞ で連続に CW に近づく: x = 0.5 で T_p κ_c = 10, 30 のとき最良のモードの r_0 = 0.996, 1.078（CW の Ω = 0 で 2 artanh 0.5 = 1.099）、離調の回転 −0.775 d, −0.794 d（CW −0.8 d）[sim]。
このとき K は T_p に比例して増える（3.5、9.9）: 長いポンプは多モード（CW の側帯波の描像）になる。

### 2.4 単一時間モードになる条件

ε(t) の時間幅とキャビティの寿命 1/κ_c の比で決まる [sim、ρ = 1、r_0 を固定して x_ss を合わせた値]:

| T_p κ_c | 二重共振（κ_p = 1.924 κ_c）: r_0 = 0.5 / 1.0 の K | 必要な x_ss | ポンプ非共振: K | 必要な x_ss |
|---|---|---|---|---|
| 0.1 | 1.086 / 1.060 | 7.3 / 14.4 | 1.001 / 1.001 | 5.2 / 10.3 |
| 0.3 | 1.089 / 1.062 | 2.5 / 4.8 | 1.008 / 1.005 | 1.8 / 3.7 |
| **0.62** | **1.103 / 1.072** | **1.22 / 2.40** | 1.030 / 1.020 | 0.97 / 1.92 |
| 1.0 | 1.134 / 1.093 | 0.80 / 1.58 | 1.072 / 1.049 | 0.67 / 1.33 |
| 2.0 | 1.289 / 1.202 | 0.47 / 0.93 | 1.256 / 1.179 | 0.43 / 0.84 |
| 5.0 | 2.14 / 1.85 | 0.31 / 0.59 | 2.17 / 1.88 | 0.30 / 0.57 |

- 短いポンプ（T_p κ_c ≪ 1）では、キャビティの中の場が一度に squeeze され、e^{−κ_c t} で出ていく。単一モードに近づくが、同じ r に要るピークのポンプ（x_ss² 倍の CW しきい値）が T_p に反比例して増える。
- 二重共振ではキャビティ内のポンプが 1/κ_p で立ち上がり・減衰するので、入力のパルスを短くしても ε(t) は ~1/κ_p より短くならず、K ≈ 1.06–1.09 で止まる。
- 両文献の T_p κ_c ≈ 0.62 は、K ≈ 1.1 と x_ss ≈ 1–2.5（CW しきい値の 1.5–6 倍のピーク）の折り合いの点にあたる [AE]。Borealis は「単一時間モードになるように幅を選んだ」[M-SI p.15、文献 [11] Vaidya 2020]、Liu も同様 [L-Methods]。
- ピークでは x_ss > 1（CW なら発振する）だが、パルスが短いので真空からの増幅は有限で、しきい値以下として扱える（このモデルは線形のまま）。

### 2.5 スクイーズ角を決めるもの

Δ = 0 の定常では conventions §4 のとおり θ_sq = φ_p,ic/2 + π/2。共振するポンプでは

  φ_p,ic = φ_p,in − arctan(Δ_p/κ_p)（定常）、パルスでは b(t) の位相が立ち上がりの間に変わる（チャープ）

信号の離調 Δ_s はさらに楕円を回す。CW（T02）では δθ_sig = −(Δ_s/κ_c)/(1 + x²)。パルス励起での値（Borealis の値、ρ = 0.95、Δ = 0 の最適モード u_0 を LO に固定して角度だけ最適化、d = Δ_s/κ_c = 0.01 と 0.05 で線形）[sim]:

| 条件（Δ_p/κ_p = 1.04 d は長さの変化に対応） | x_ss = 0.5 | 1.3 | 2.4 | CW の参考 |
|---|---|---|---|---|
| 信号の離調だけ（Δ_p = 0）: δθ/d | −0.660 | −0.638 | −0.591 | −1/(1 + x²) = −0.80 / −0.37 / −0.15（x < 1 のみ意味あり） |
| ポンプの離調だけ（Δ_p = 2Δ_s）: δθ/d | −0.394 | −0.425 | −0.499 | −½ arctan(Δ_p/κ_p)/d = −0.52 |
| 長さの変化（両方）、ポンプ位相のロックなし | −1.054 | −1.063 | −1.090 | — |
| 同、透過ポンプのロックあり（arg ∫b dt を固定） | −0.534 | −0.543 | −0.571 | −1/(1 + x²) |

- 符号: d > 0（キャビティが短くなる）で θ_sq が減る（T02 V-B7 と同じ）。squeezing の量の劣化は d の 2 次（d = 0.05、x_ss = 2.4 で −7.431 → −7.405 dB）。
- パルス励起では信号の離調の回転が x にあまり依存しない（−0.6 d 前後）。CW の −1/(1 + x²) は使えない。

**長さと温度への換算**:
- Δ_s/κ_c = −2F δl/λ0（δl: 往復の光路長）。Borealis F = 57: 0.0735/nm、Liu F = 59: 0.0761/nm。
  ポンプのロックありの残差は |δθ| ≈ 0.55 × 2F δl/λ0 ≈ **40 mrad/nm**（Borealis）、**42 mrad/nm**（Liu）。10 mrad 以下には δl ≲ 0.25 nm（往復）[AE]。
- ロック光が squeezed 光と違うモード（波長・偏光）を共振させるとき、PDH は (n_L L_c + D) を一定に保つだけなので、信号の共振は結晶の温度で動く:
  Δ_s/κ_c ≈ −4F L_c Δ(dn/dT) δT_c/λ0、Δ(dn/dT) = (dn/dT)_signal − (dn/dT)_lock（熱膨張の項は省略）[AE]。
  - Borealis: ロック光は 775 nm 付近の直交偏光（n_y）、信号は 1550 nm の n_z。Liu: ロック光は 775.7 nm の同じ偏光（n_z）、信号は 1551.4 nm の n_z。
  - KTP の dn/dT の値はこの担当は確認していない **[UNVERIFIED]**。Δ(dn/dT) = 10⁻⁶–10⁻⁵ /K を仮に置くと Δ_s/κ_c = 1.5–15 /K（L_c = 10 mm、F = 57）、δθ ≈ 0.8–8 rad/K、10 mrad には 1–12 mK の温度安定度 [AE]。
  - ポンプの共振（Δ_p）も同じ形で温度に依存する。Liu ではロック光がポンプと同じモードなので Δ_p は帯域内で 0。Borealis ではロック光と直交偏光なので Δ_p も温度で動き、二重共振は結晶の温度で保つ [M-SI p.15]。

### 2.6 ポンプパワーの揺らぎ

r_0 と x_ss の関係（Borealis の値）: r_0 = 0.533 ↔ x_ss = 1.297、0.886 ↔ 2.136、1.0 ↔ 2.401、1.5 ↔ 3.531、1.8 ↔ 4.180（Liu の値ではそれぞれ 1.279、2.107、2.369、3.486、4.129）[sim]。
したがって d ln r_0/d ln P ≈ 0.51（r = 0.5–1）、平均光子数は d ln(sinh² r)/d ln P ≈ r coth r（r = 1 で 1.31）[AE]。
Liu は大きい光子数の実験で分布がポンプのドリフトで広がったと書いている [L-SI §3]。Borealis は透過ポンプのパワー（PD2）を AOM で安定化する [M-SI p.15]。

### 2.7 FSR の櫛（非縮退のモードの組）

キャビティの共振の櫛により、ω0 ± m·FSR の組（m = 1, 2, …）がそれぞれ 2 モード squeezing（各モード単独では熱的、平均 sinh² r_m）を作る。ポンプの帯域（~κ_p）≪ FSR なので、組は対角（(m, −m)）だけ [L-SI §1.2]。
位相整合の帯域は FSR よりずっと広いので r_m ≈ r_0 が多数の m について続く [L-SI §1.2]（二重共振の FSR の分散で m の大きい組が外れる効果は文献で確認していない、パラメータ m_max）。

- BHD（LO は ω0 の近くだけ）には見えない。光子計数（GBS）には、フィルタの漏れ T_f(m) の分だけ区別できる熱的な光子 n_FSR = Σ_m 2 T_f(m) sinh² r_m が加わる。
- Borealis: 2 枚の FBG（平らな反射帯 約 5 GHz、互いに 約 1 GHz ずらす）で ±1, ±2, … を 30 dB 以上抑圧、10 nm のバンドパス（室内光・ポンプを約 60 dB）、フィルタとファイバ結合で透過 約 80% [M-SI p.19]。
- Liu: 不平衡 MZI 3 段（光路差 L_c, 2L_c, 4L_c。次数 ±8, ±24, ±40, … 以外を除く）、各 99.8% 透過・40 dB 以上、その後 0.4 nm のバンドパス 3 層（平均 98.5% 超）、モード整合 99.95% 超、合計 95% 超 [L-Methods、L-本文]。
  MZI の式 T(ν) = ½ + ½ cos(2πL(ν − ν′)/c) [L-Methods Eq. 2]。

## 3. 仮定と適用範囲

- しきい値以下・ポンプの枯渇なし・ガウス状態・線形化。ピークで x_ss > 1 でも、パルスの間の増幅が有限なら線形のまま（§2.4）。パルスあたり sinh² r_0 ≲ 3（r ≤ 1.8）の範囲。
- フィネス ≫ 1（57–59）の単一縦モード（平均場）。結晶内の伝搬・GVM（10 mm で ps 級）・ガウスビームの発散は ns のモードには効かない（Liu SI §1.2 も発散を無視）。
- ポンプのキャビティは線形のフィルタ（b の方程式）。ポンプの入力は理想的な矩形（実際の立ち上がり時間はパラメータ。Borealis は 4 GHz、Liu は 20 GHz の LN 変調器 [M-SI p.15、L-SI §1.1]）。
- 離調は 1 パルスの間は一定（離調の雑音は DC–kHz、パルスは 10–20 ns）。パルスごとに Δ(t_j) をサンプルする（準静的）。
- 「r」の値: Borealis Table S2 の r（0.533、0.886）は SI で参照面（光源の直後か）が書かれていない。216 × sinh²(0.533) = 68 光子と測定の ⟨N⟩ = 21.1 の比 0.31 は、共通の透過 0.482 × デマルチ 0.93 × ループの透過（0.85–0.92 を数回）と矛盾しないので、光源の面の値と読んだ [AE]。
- P_th,CW は両文献とも書いていない。x_ss を直接のパラメータにする。Borealis の平均 3.7 mW・6 MHz・3 ns から P_pk ≈ 206 mW [AE]。これが r = 0.533 の実験なら P_th,CW ≈ 123 mW、r = 0.886 なら 45 mW（どちらの実験の値かは SI に書かれていない）[AE]。
- 最適な LO モードは理想的に整形できると仮定した値。実際の LO の波形（§6.3）は別のパラメータ。

## 4. 周波数帯

| 効果 | 帯域 | 備考 |
|---|---|---|
| squeezed のモード（(0,0) の帯域） | 光の周波数で ±(数 κ_c)。|ũ(f)|²/τ_eff = 0.46（30.5 MHz）、0.10（66 MHz）、9.8×10⁻³（125 MHz）、2.5×10⁻⁶（322 MHz）（Borealis の最適モード）[sim] | Liu ではそれぞれ 0.78、0.39、0.11、2.1×10⁻³ |
| FSR の櫛 | ±m × 3.7 GHz（Borealis）、±m × 7.4 GHz（Liu） | FBG・MZI で除く |
| ポンプのキャビティ | 127 MHz（Borealis）、255 MHz（Liu [AE]）FWHM | ε(t) の立ち上がり |
| 離調・長さ雑音・温度 | DC–kHz（時間ループ段） | パルスごとの準静的な角度の回転 |
| ポンプパワーのドリフト | 分–時間 [M-SI p.15] | AOM で補正 |
| Borealis のロック光 | +7 GHz（1550 nm、ループ用）、2 × (ν0 + 7 GHz)（OPO 用、775 nm 付近）、PDH の変調 40 MHz | ループ用は FBG でさらに 50 dB 落とす [M-SI p.22] |
| Liu の OPO ロック光 | ポンプと同じ周波数、100 MHz の PDH 変調、1 µs のパルス（計算に使う時間の外） | [L-SI §1.1、Fig. S4] |
| GAWBS の R0,m | 30.5、82.1、…、323 MHz（125 µm、T07） | 低次の 1–2 本だけが OPO のモードの帯域に入る |
| Brillouin | ±10.8 GHz | OPO のモードの帯域の外（§5.4） |

## 5. スケーリングと、我々の構成での桁の見積もり

### 5.1 パラメータ（Borealis と Liu のプリセット）

| 量 | Borealis（Madsen 2022 SI） | Liu 2026 |
|---|---|---|
| 波長 λ_s / λ_p | 1550 / 775 nm [M-SI p.15] | 1551.4 / 775.7 nm [L-SI §1.1] |
| 結晶 | PPKTP 10 mm、type-0、ヘミリシック、二重共振 [M-SI p.15] | 同（10 mm）[L-SI Fig. S2] |
| 鏡 | RoC1 12 mm（HR_s 99.9%、HR_p 90%）、AR < 0.2%、RoC2 25 mm（S2_s 90%、S2_p 90%）、D 22 mm（可変）[M-SI Fig. S7] | RoC1 12 mm（HR_ω 99.978%、HR_2ω 90%）、AR < 0.05%、RoC2 7 mm（S2_ω 90%、S2_2ω 90%）、D 2.1 mm [L-SI Fig. S2] |
| ウエスト | 40 µm（1550）、28 µm（775）[M-SI Fig. S7] | 記載なし |
| 往復の光路長 l | 81 mm [AE: c/FSR、n_KTP ≈ 1.816 で 2(1.816·10 + 22) = 80.3 mm と整合] | 40.4 mm [L-SI §1.1] |
| FSR | 3.7 GHz（両波長）[M-SI Fig. S7] | 7.4 GHz [L-SI §1.1] |
| 線幅 FWHM（信号）| **66 MHz**、フィネス 57 [M-SI Fig. S7, p.18] | 125 MHz [AE: FSR/F]、フィネス 59 [L-SI §1.1] |
| 線幅 FWHM（ポンプ）| 127 MHz、フィネス 29 [M-SI Fig. S7] | 255 MHz [AE]、フィネス 29 [L-SI §1.1] |
| κ_c、1/κ_c | 2.07×10⁸ s⁻¹、4.8 ns [AE] | 3.94×10⁸ s⁻¹、2.5 ns [AE] |
| 取り出し効率 ρ | 鏡の値から 0.955–0.991（AR を上限 0.2% と 0 で）、線幅と FSR から 0.93–0.95 [AE]。参照値 0.95 | 0.988–0.998（AR の上限と 0）[AE]。参照値 0.99 |
| ポンプ | 3 ns 矩形、6 MHz、平均 3.7 mW @775 nm（0.62 nJ/パルス、ピーク 約 206 mW [AE]）。4 GHz LN 強度変調器＋AWG → EDFA → ファイバ結合 MgO:LN リッジ導波路で SHG [M-SI p.15] | 1.6 ns 矩形、間隔 50 ns（20 MHz）。20 GHz LN 振幅変調器 → EDFA → AOM（パルス数の切り替え）→ バルク PPLN で SHG [L-SI §1.1] |
| T_p κ_c | 0.622 [AE] | 0.630 [AE] |
| ポンプパワーの安定化 | OPO の前の AOM、OPO を透過したポンプを PD2 で監視 [M-SI p.15] | 記載なし（ドリフトの補正を解析で入れた [L-SI §3]） |
| OPO の長さロック | PDH。ロック光 = 主レーザーを約 7 GHz で位相変調 → エタロン（線幅 2 GHz、FSR 60 GHz）で +1 次だけ → 増幅 → 40 MHz の位相変調（RedPitaya）→ 温度制御したリッジ導波路で SHG（約 775 nm）。ポンプと直交偏光、CW、ポンプと独立。透過を PD1 で検出。アクチュエータは PZT の鏡 [M-SI p.15–16, Fig. S6] | PDH。100 MHz の位相変調 → EDFA → 1 µs の矩形 → PPLN 導波路で SHG。ポンプと同じ偏光・同じ中心周波数（AOM を同じ周波数・次数で駆動）、時間をずらす（回路の最大遅延の外）、可変の BS でポンプと合わせる [L-SI §1.1, Fig. S4] |
| ポンプ位相のロック | 記載なし | ポンプ光路の PZT 鏡。OPO を透過したポンプを二色鏡で分け、参照光と干渉（平衡ホモダイン）。精度 1 nm [L-Methods] |
| 出力 | 二色鏡でポンプと分け、SMF に結合 [M-SI p.15] | — |
| FSR フィルタ | FBG 2 枚（5 GHz の平らな帯、ずれ 1 GHz）30 dB 以上、10 nm バンドパス 60 dB。フィルタ＋結合で 約 80% [M-SI p.19] | 不平衡 MZI 3 段（99.8%、> 40 dB）＋ 0.4 nm フィルタ 3 層。合計 > 95% [L-Methods] |
| 光源の効率 | 記載なし（共通の透過 0.482(9) に含まれる [M-SI Table S3]） | 92% [L-本文] |
| r | 0.533（216 モード、⟨N⟩ = 21.12）、0.886（72 モード、⟨N⟩ = 22.39）[M-SI Table S2] | 0.9–1.8 [L-本文] |
| 単一モード性 | g⁽²⁾_min = 1.91 → K ≈ 1.1 [M-SI Fig. S12] | 識別不能性 97%（r = 1.8）[L-本文]。遅延 50 ns、800 ns、12.8 µs で変わらない [L-SI §1.5] |
| squeezed パルスの幅 | 「10 ns」、ゲートの評価窓 25 ns、ゲート 167 ns [M-SI p.16, p.20] | — |
| 時間ビンの合わせ | 250 ps の精度 [M-SI p.18] | 波束の重なり > 99.7% [L-Methods] |
| 遅延ループ | τ ≈ 167 ns（34 m）、6τ、36τ。ループの損失 5/7/12%（結合＋ファイバ）、結合 約 97%。透過（往復＋VBS 2 回目）0.918(9)、0.889(5)、0.850(5)。分散 < 18 ps/(nm·km)、PMD < 0.04 ps/√km [M-SI p.18, Table S3] | τ = 50 ns、16τ、最大遅延 12,800 ns。干渉計 99%/モード [L-Methods] |
| ループのロック | 1550 nm、+7 GHz、逆向き、99.5% の BS で入れ、自由空間サーキュレータで取り出す（squeezed 光に 5% の損失）、1 µW。約 10⁻³ が PNR 側に戻り、FBG でさらに 50 dB。AOM でゲート: ロック光 65 µs、squeezed 光 35 µs（10 kHz）。VBS で 1.3/3.5/1.8 MHz の振幅変調（深さ 約 V_π/2/10）、5 MHz の PD、lock-in、PyRPL の PI、PZT 鏡（ループ 1）とファイバストレッチャ（ループ 2、3）、コイルは温度安定化した発泡スチロールの箱 [M-SI p.21–23] | 逆向き、同じ周波数、直交偏光、AOM でパルス化して時間をずらす。PBS で入れ出し、参照光と干渉。PZT 鏡＋TEC。精度 3 nm、10 nm、全体 λ/200 [L-Methods, L-SI §4.2.3] |
| 振り分け・検出 | 1→16、共振型 EOM 15 個（3/1.5/0.75/0.375 MHz）、消光比 > 200:1、透過 > 90%（出口の結合 95% を含む）、相対 0.930(3) [M-SI p.23, Table S3]。PNR（TES）、迷光 0.003 光子/パルス [M-SI p.26] | SNSPD、検出効率 93%、回復 43 ns。全体の効率 51% [L-本文] |

注: T11 は Madsen 2022 の本文（PMC）から「振り分けの損失 約 15%」「ループの位相雑音 0.02/0.03/0.15 rad」を引いている。SI の「振り分けの透過 > 90%」とは定義（結合・検出を含むか）が違う可能性がある。
また T11 §2 の「測定 65 µs と安定化 35 µs を交互」は SI と逆（**ロック光 65 µs、squeezed 光 35 µs**、[M-SI p.22]）。

### 5.2 時間モードの数値（最適な LO モード、ρ = プリセットの値、x_ss = 2.4、r_0 ≈ 1）[sim]

| 量 | Borealis | Liu | D07 の ps パルス（比較） |
|---|---|---|---|
| τ_eff = (∫u dt)² | 3.00/κ_c = **14.5 ns** | 2.97/κ_c = **7.5 ns** | 2.13 τ_FWHM（1 ps で 2.1 ps） |
| モードの中心（ポンプの立ち上がりから） | 1.43/κ_c = 6.9 ns | 1.43/κ_c = 3.6 ns | — |
| 強度の rms 幅 | 0.62/κ_c = 3.0 ns | 0.62/κ_c = 1.6 ns | — |
| 重なり 0.99 の時間ずれ | **436 ps** | **227 ps** | 8.5 fs（100 fs）– 0.85 ps（10 ps） |
| 重なり 0.999 | 137 ps | 71 ps | — |
| 裾のエネルギー ε_tail(10 ns) | 0.138 | 3.4×10⁻³ | — |
| ε_tail(20 ns) | 2.3×10⁻³ | 1.2×10⁻⁶ | — |
| 近似 ε_tail(T) | ≈ 8.8 e^{−2κ_c T}（両方、κ_c T ≥ 2 で 5% 以内） | | |
| R_−/R_+（ρ = 0.95 / 0.99、最適 LO） | x_ss = 2.5 で −7.66 / +8.75 dB | x_ss = 2.5 で −8.99 / +9.26 dB | — |

最適モードの squeezing は、単一モード r_0 に純損失 ρ をかけた値とほぼ同じ（Borealis x_ss = 2.5: r_0 = 1.029、0.95 e^{−2r_0} + 0.05 → −7.67 dB）。**K ≈ 1.1 の多モード性は最適モードの squeezing をほとんど下げず、GBS の純度（識別不能性）に効く。**

### 5.3 将来の 32 ポート構成（D08・D39）への当てはめ

- **繰り返しの上限**（裾 ε_tail ≤ 10⁻³ を次のスロットまでに）: T_rep ≥ ln(8.8/ε)/(2κ_c) + （ガード）= 21.9 ns（Borealis 型）、11.5 ns（Liu 型）→ f_rep ≲ 45 MHz、87 MHz。ε = 10⁻⁴ なら 27.4 ns、14.4 ns。
  100 MHz（T_rep = 10 ns）は Borealis 型で 13.8%、Liu 型で 0.34% が次のスロットに入る（T11 §3.3 の例外 (a)「隣のパルスと重なる」に当たり、ε は損失ではなく**コヒーレントな混合**になる）。
- **振り分けのエッジ条件**（T11 §3.4）: t_r + 6σ_j + （モードの長さ ~ε_tail で決まる 22 ns）+ t_ring ≤ N_b T_rep。PC（t_r 8–10 ns）でパルスごと（N_b = 1）なら T_rep ≳ 35 ns だが、PC の 3 MHz の切り替えのほうが先に効く（f_rep ≤ 6 MHz）。
  Liu の 20 MHz を PC の木で振り分けるなら N_b ≥ 20/(2·3) → 4。L_31 = 31·4·c·50 ns/n_g = 1.27 km [AE]。
- **遅延ファイバのタイミング**: 6 MHz、N_b = 1 で L_31 = 1.06 km（T11）、熱による 41 ps/K。重なり 99% の 436 ps に対して約 10 K の余裕。T_rep のドリフト δT/T: δt_31 = 31 T_rep δT/T、436 ps には δT/T ≤ 8×10⁻⁵（6 MHz）[AE]。電子回路で作るパルス（AWG）なら問題ない（T11）。
- **分散**: Borealis の見積もり（7.3 km で < 66 fs、PMD < 110 fs）[M-SI p.19] のとおり無視できる。L_D = T0²/|β2| ≈ (3 ns)²/(21.7 ps²/km) ≈ 4×10⁸ m [AE]。
- **散乱雑音**（逆向きの probe、D30。T10 の Φ_st = 3.0×10¹⁰ /s @ 1 mW・63 m・同じ偏光・同じ周波数の Rayleigh を使う）[AE]:

| | Borealis 型（τ_eff = 14.5 ns） | Liu 型（7.5 ns） | 1 ps パルス（τ_eff = 2.1 ps、T10/T07） |
|---|---|---|---|
| n_mode（Rayleigh、1 mW、63 m） | 435 | 226 | 0.064 |
| n_mode ≤ 5.8×10⁻⁴（−10 dB → −9.95 dB）の P̄_pr ε_pol（63 m / 630 m） | 1.3 nW / 0.13 nW | 2.6 nW / 0.26 nW | 9 µW / 0.9 µW |
| PM・直交偏光（ε_pol = 2×10⁻³）なら | 0.67 µW / 67 nW | 1.3 µW / 0.13 µW | — |
| T10 式 (11) の n_mode σ_lock²（10 kHz、同じ偏光） | 7.1×10⁻¹⁰ rad²（63 m）→ σ_lock ≥ 1.1 mrad、630 m で ≥ 3.5 mrad | 3.7×10⁻¹⁰ rad² → ≥ 0.8 mrad、2.5 mrad | 1 ns ゲートで 4.9×10⁻¹¹ |
| 周波数をずらした probe（+7 GHz、FBG なし）の重み |ũ(7 GHz)|²/τ_eff | ≤ 2.2×10⁻⁵（片側指数のモードの Lorentz 形、保守的）。最適モードの計算では 10⁻⁹ 以下（数値の床）[sim] | ≤ 8.0×10⁻⁵ | 約 1（スペクトルの内側） |
| Brillouin（10.8 GHz）の重み | ≤ 9.3×10⁻⁶ → N_B ≈ 3×10⁻⁷（1 mW、63 m）| ≤ 3.4×10⁻⁵ | 0.998（N_B = 4.5×10⁻³） |
| Raman（広帯域） | 1 ps と同じ N_R ≈ 9×10⁻⁶（1 mW、63 m） | 同じ | 9×10⁻⁶ |

- 一般式: 狭帯域の雑音（周波数 f_0 に集中した光子流 Φ）は N = Φ |ũ(f_0)|²、f_0 = 0 で Φ τ_eff。広帯域（モードの帯域の中で平ら、密度 n_0）は N = n_0（τ によらない）。conventions §2 の N = ∫n(f)|ũ(f)|² df の 2 つの極限。
- 同じ向きの明るい光の GAWBS: N_G = Σ_m (P σ²_G,m/(hν0)) |ũ(f_m)|²（T07 の N_G ≈ τ_eff P σ_G²/(hν0) を f_m τ ~ 1 に拡げた形）[AE]。最強の R0,7（322 MHz）は Borealis の最適モードで |ũ|²/τ_eff = 2.5×10⁻⁶ に落ちる。
- 検出ゲートあたりの迷光（GBS のクリック、D46）はゲート長に比例: ゲートはモードの長さ（~25 ns、Borealis の評価窓）以上が要るので、T10 の 1 ns の値の約 25 倍 [AE]。
  Borealis はロック光を時間で分ける（ロック 65 µs と測定 35 µs を交互）ことと、+7 GHz の周波数のずれ＋FBG（50 dB）で避けている。

### 5.4 スケーリングのまとめ

| 量 | κ_c | T_p | x_ss / r | ρ | f_rep | ファイバ長 L |
|---|---|---|---|---|---|---|
| τ_eff、時間ずれの許容 | ∝ 1/κ_c（3.0/κ_c、0.090/κ_c） | 弱い（T_p κ_c ≲ 1） | 弱い | — | — | — |
| K | T_p κ_c の関数（§2.4） | 大きい T_p で ∝ T_p | r が大きいほど 1 に近づく | — | — | — |
| 必要なピークのポンプ | — | 短いほど増える | x_ss ≈ 2.4 r（T_p κ_c = 0.62） | — | 平均 ∝ f_rep | — |
| 裾による f_rep の上限 | ∝ κ_c | — | — | — | — | — |
| 離調の感度 δθ/δl | ∝ F/λ0（≈ 0.55 × 2F/λ0） | — | 弱い | — | — | — |
| 狭帯域の散乱雑音 n_mode | ∝ 1/κ_c | — | — | — | — | ∝ L |
| R_−（最適 LO） | — | — | r | ρ e^{−2r} + 1 − ρ | — | 損失 |

### 5.5 Ref D

Ref D（D02）には squeezer がない。このトピックは Ref D には関係しない。Ref B・C（CW の OPO）は T02 のまま（§2.3 のとおりこのモデルの CW の極限）。

## 6. Txx の結論の変わり方（100 fs–10 ps → OPO の ns モード）

### 6.1 一覧

| トピック | 100 fs–10 ps（D07）での結論 | OPO の ns モード（Borealis・Liu 型）での結論 |
|---|---|---|
| T02 | CW の OPO（Eq. 2）と D07 のパルスは両立しない。離調の回転 −(Δ/κ)/(1 + x²)。同軸 probe の残差 | パルス励起の Langevin（§2）で扱える。CW は T_p → ∞ の極限。離調の回転は −0.6 Δ_s/κ_c 前後で x にあまり依存しない。ポンプも共振するので −½ arctan(Δ_p/κ_p) が加わる。透過ポンプのロックでも −0.55 Δ_s/κ_c が残る |
| T03 | 強い多モード（K ≳ 200）。時間ずれの許容 8.5 fs–0.85 ps。熱ドリフト 3.5 mK–0.35 K。分散でポート間の重なりが壊れる（100 MHz で τ < 3.7 ps は補償が必要） | K ≈ 1.03–1.1（T_p κ_c ≈ 0.62）。許容 0.23–0.44 ns。熱ドリフトは約 10 K まで無関係。分散は無関係。**新しく LO の時間波形が問題**（ガウス 0.951、矩形 0.864、片側指数 0.84） |
| T06 | 100 fs–1 ps では LO の SPM・ソリトン・SSFS が効く。分散補償が必要 | GNLSE は不要（時間ごとの SPM だけの準 CW、φ_NL(T) = γ P(T) L_eff）。例: 1 mW 平均・6 MHz・15 ns の LO でピーク 約 11 mW、63 m で 0.9 mrad [AE]。ただし **SBS はパルスが音響の寿命（~10 ns）と同程度なので抑えられない**（T07 の CW の閾値をピークパワーで使う） |
| T07 | Brillouin は ps のスペクトルの内側でフィルタできない。Rayleigh の約 1/12。モードあたりの雑音は τ_eff に比例 | Brillouin（10.8 GHz）はモードの帯域の外（重み ≤ 10⁻⁵）。Rayleigh（狭帯域）は τ_eff に比例して 3.5×10³–6.8×10³ 倍。Raman（広帯域）は変わらない。GAWBS は低次のピークだけが帯域に入る |
| T10 | 1 mW の probe は 6 桁大きい。nW–µW で足りる。直交偏光の PM は温度ドリフト | 同じ周波数の probe は ≲ 1–3 nW（同じ偏光、63 m）に下げる必要。ロックの位相誤差に 1–3.5 mrad の下限（10 kHz）。**周波数をずらした probe（Borealis の +7 GHz）＋FBG が使えるようになる**（帯域 100 MHz ≪ 7 GHz）。ずらした分の温度ドリフト 2πΔν·(dτ_g/dT)L ≈ 0.11 rad/K（7 GHz、63 m）[AE] |
| T11 | τ_p ≪ T_rep は常に成り立つ。漏れは損失。タイミングの許容は τ_p に比例して厳しい。PC で f_rep ≤ 6 MHz | **τ_mode ≪ T_rep が成り立たない場合がある**（裾 8.8 e^{−2κ_c T}）→ f_rep ≲ 45/87 MHz。ゲートの長さ ≥ 約 25 ns。タイミングの許容は緩い。パルスごとに読む BHD は 200–500 MHz 級の帯域で足りる [AE]（100 MHz の ps パルスでは GHz 級） |
| T12 | 上流は導波路 2%＋光学部品 1%＋結合 2% の仮の値（η_up = 0.951） | 上流を ρ（0.93–0.99）× FSR フィルタ（0.80 Borealis / > 0.95 Liu）× 結合に置き換える。Borealis 型では η_up ≈ 0.76 で、10 dB の目標（全損失 ≤ 8.4%）は不可能、Liu 型（0.92–0.94）でも上流だけで 6–8% を使う。LO の時間波形の整合を別の損失として持つ |

### 6.2 なぜ変わるか（要点）

- モードの時間が 3–4 桁長くなり、帯域が 3–4 桁狭くなる。時間の許容は長さに比例して緩み、狭帯域の雑音の光子数は長さに比例して増え、帯域の外の雑音は |ũ(f)|² で落ちる。
- 多モード性は「結晶の GVM と位相整合」ではなく「ポンプの時間幅とキャビティの寿命の比」で決まる。
- パルスの時間幅（D07 の τ_p）は OPO では設定値ではない。**OPO のプリセットの掃引のパラメータは (Δν_c, T_p, κ_p, x_ss, ρ) にするべき**（§9 の提案）。

### 6.3 LO の時間モード（T12・T13・D45 への入力）

最適モードとの重なり（エネルギー）[sim、両プリセットでほぼ同じ]: 片側指数 e^{−κ_c t}（遅延を最適化）0.841–0.845、ガウス（中心と幅を最適化）0.951、矩形の窓（幅 2/κ_c が最良）0.864。
- LO を EOM で最適モードの形に切り出す（GHz の変調器で ns の形は作れる）か、CW の LO と速い BHD でデジタルに重みをかける（BHD の帯域 ≳ 数 κ_c/π）。
- 重なりの不足は純損失（T12 §1.3）。ガウスの LO だけで 4.9% の損失。

## 7. 検証テスト（実装が通すべきもの）

数値はこの担当の計算 [sim]（δt = 1/(40–160 κ_c)）。許容差は離散化の誤差（δt² に比例）を見て決めた。

| # | テスト | 合格条件 |
|---|---|---|
| V-P1 | CW の極限（ε 一定、窓 [0, 40/κ_c]）: 中央の幅 20/κ_c の平らなモードで R を §2.3 の解析式と比べる（x = 0.5、ρ = 0.95）: R_p = 0.183704、R_q = 7.84003 | (R − 1) の相対誤差 < 10⁻⁴（δt = 1/(80κ_c) で 2×10⁻⁵）。δt を半分にして誤差が 1/4（20→40→80: 3.4×10⁻⁴、8.5×10⁻⁵、2.1×10⁻⁵）。核の要素 C(t, t′) も解析式と一致 |
| V-P2 | 真空: x_ss = 0 で C の固有値がすべて 1（prompt §7-4 の OPO 版） | \|λ − 1\| ≤ 2 (κ_c δt)² |
| V-P3 | 初期のキャビティ内の真空の項（§2.2 の最後の項）: これを除くと ρ = 1 で固有値 ≈ 0 が出る（回帰テスト） | 項ありで最小固有値が x_ss = 0 で 1 |
| V-P4 | 純粋さ: ρ = 1（離調・複素のポンプを含む: x_ss = 2.4、d = 0.05、Δ_p = 0.1κ_c）で固有値が λ_i λ_{2N+1−i} = 1 の対 | δt = 1/(40κ_c) で \|λλ′ − 1\| ≤ 1.3×10⁻³、δt² で減る |
| V-P5 | **Borealis の K**: T_p κ_c = 0.622、κ_p/κ_c = 1.924、ρ = 1、r_0 = 0.533 → x_ss = 1.297、K = 1.101 | K = 1.10 ± 0.01（測定 g⁽²⁾ = 1.91 → K ≈ 1.1 [M-SI Fig. S12] と ±0.05 で一致） |
| V-P6 | **Liu の純度**: T_p κ_c = 0.630、κ_p/κ_c = 2.034、r_0 = 1.8 → x_ss = 4.129、K = 1.026 | 1/K = 0.975 ± 0.005（Liu の 97% と ±0.01） |
| V-P7 | 長いポンプ → CW（ポンプ非共振、x = 0.5）: T_p κ_c = 10、30 で r_0 = 0.996、1.078 → 2 artanh x = 1.099 に単調に近づく。離調の回転 −0.775 d、−0.794 d → −0.8 d | 単調に近づくこと、T_p κ_c = 30 で r_0 の差 < 3%、回転の差 < 1% |
| V-P8 | 位相の規約: Δ = 0、φ_p,in = 0 で squeezed 直交位相が p（θ_sq = π/2）。φ_p,in = 0.4 で θ_sq = π/2 + 0.2（mod π）（conventions §4） | 10⁻⁴ rad |
| V-P9 | 離調（Borealis の値、ρ = 0.95、x_ss = 1.3、d = 0.01）: 信号だけ −0.638 d、ポンプだけ（Δ_p = 2Δ_s）−0.425 d、両方 −1.063 d、透過ポンプのロック −0.543 d。符号: d > 0 で θ が減る | ±0.005 d。d = 0.05 でも係数が ±0.002 d で同じ（線形） |
| V-P10 | 時間ずれ: 最適モードの重なり 0.99 / 0.999 が \|Δt\| = 0.090/κ_c / 0.028/κ_c（両プリセット） | ±3% |
| V-P11 | 裾: ε_tail(T) ≈ 8.8 e^{−2κ_c T}（κ_c T ≥ 2）。Borealis で ε_tail(10 ns) = 0.138、ε_tail(20 ns) = 2.3×10⁻³ | ±10% |
| V-P12 | 短いポンプ（非共振、T_p κ_c = 0.1、r_0 = 0.5）で K → 1.001。二重共振では K ≥ 1.06 で止まる（§2.4 の表） | ±0.005 |
| V-P13 | 純損失: ρ を変えたとき最適モードの R_− が ρ e^{−2r_0} + 1 − ρ にほぼ一致（Borealis x_ss = 2.5: −7.66 dB 対 −7.67 dB） | ±0.05 dB（多モード性の分だけずれてよい） |
| V-P14 | 狭帯域・広帯域の雑音: n(f) = Φ δ(f − f_0) で N = Φ\|ũ(f_0)\|²、平らな n_0 で N = n_0（τ_eff によらない）を、最適モードと 1 ps のガウスで確認 | 10⁻⁶ |
| V-P15 | FSR の櫛: T_f(m) = 0 で n_FSR = 0、T_f = 1 で各組が平均 sinh² r_m の熱的な光子を 2 つのモードに出す | 厳密 |

## 8. シミュレータに必要なもの（すべて切り替えられる設定、D34）

### 8.1 ブロックと表

| 名前（案） | 段 | 入力 | 出力 | 切り替え |
|---|---|---|---|---|
| `pulsed_opo_source` | O | Δν_c（または FSR と F、または T_oc・L_ic・l）、ρ、κ_p（None = 非共振）、ポンプの波形（矩形 T_p、立ち上がり、任意の配列）、x_ss（または P_pk と P_th,CW）、Δ_s、Δ_p、φ_p,in、時間格子 | C（2N×2N）、r_n、K、最適 LO モード u_opt(t)、τ_eff、\|ũ(f)\|²、ε_tail(T)、R_±（任意の LO 波形で）、表 (V_−, V_+, θ_off)(Δ_s, Δ_p, x_ss) | `opo.pump_resonant`、`opo.pulsed`（off → T02 の CW ブロック）、`opo.include_initial_vacuum`（検証用。常に on） |
| `fsr_comb_filter` | O | FSR、r_m（callable、m_max）、フィルタの透過 T_f(ν)（FBG の組、MZI の段、バンドパス）、通過帯の損失 | (0,0) の損失、n_FSR | `opo.fsr_modes`（off = 理想フィルタ）、プリセット `fbg_pair`・`mzi_cascade`・`none` |
| `opo_length_lock` | T | ロック光の種類（Borealis: 2(ν0 + 7 GHz)、直交偏光、CW。Liu: 2ν0、同じ偏光、1 µs パルス）、PDH の残差 S_δl(f)、T_c(t)、Δ(dn/dT)（信号・ポンプそれぞれ）、L_c、F | Δ_s(t)、Δ_p(t)（パルスごとにサンプル） | `opo.length_noise`、`opo.thermal_detuning` |
| `opo_angle_map` | T | Δ_s、Δ_p、x_ss（表）、ロックの状態 | δθ_sq(t)（パルスごと）、R_− の劣化 | — |
| `pump_power_stabilization` | T | ポンプの RIN・ドリフト、AOM のループ（透過ポンプを監視）、Δ_p | x_ss(t) | `pump.power_lock`（Borealis: on、Liu: off） |
| `pump_phase_lock_transmitted` | T | 透過ポンプの位相 arg∫b dt、参照光（2ω0）の位相、PZT（ポンプ光路） | φ_p,in の補正 | `pump.phase_lock`（Liu: on、Borealis: off）。D42 の比較の候補 |
| `lo_temporal_shape` | O | LO の波形（最適、ガウス、片側指数、矩形、デジタルの重み）、BHD の帯域 | 重なり η_t、θ_off | D45 と共通 |
| 雑音のモード重み | O | n(f)（T07・T10 の各源）、\|ũ(f)\|² | n_mode、N_stray（ゲート長 T_gate） | T07・T10 と共通 |
| 繰り返しの整合性チェック | O | f_rep、N_b、ε_tail、t_r | 警告（裾の混合 > 10⁻³） | 悪い構成でも計算は続ける（D34） |

### 8.2 プリセット（同格、D34）

- `source.borealis_opo`（§5.1 の左列）、`source.liu2026_opo`（右列）、`source.refA_waveguide_opa`（T01/T03、既存）、`source.refB_cw_opo`・`source.refC_cw_opo`（T02）。
- 比較用（悪い・極端な構成）: Borealis 型 OPO を 100 MHz で使う、FSR フィルタなし、ガウス・矩形の LO、長いポンプ（T_p κ_c = 5、K ≈ 2）、ポンプ非共振、直交偏光のロック光で結晶の温度制御なし、同じ周波数の 1 mW の逆向き probe。
- ロック系のプリセット: Borealis（ポンプ位相のロックなし、パワーの AOM、+7 GHz の逆向きロック光を時間で分ける）、Liu（透過ポンプのホモダインロック、同じ周波数・直交偏光の逆向き probe、TEC）。

### 8.3 State・光路のモデルへの要求

- OPO を「区間」として持つとき、係数は「区間 × ビーム」ごと（conventions §9 の既存の規則）: squeezed 光には §2.5 の係数、透過ポンプには −arctan(Δ_p/κ_p)、OPO を通らない逆向き probe には 0。
- ポンプのキャビティをポンプ経路の区間とし、その位相 −arctan(Δ_p/κ_p) が θ_sq に ½ で入る（T01 の「結晶をポンプ経路の区間として扱う」と同じ扱い）。
- パルスごとに Δ_s(t_j)、Δ_p(t_j)、x_ss(t_j) をサンプルする（時刻 t_j はポンプの出射時刻、conventions §7）。

## 9. DECISIONS・conventions との食い違い、提案（調整役が判断）

1. **conventions §4 への追加の提案**: κ_p（ポンプのキャビティの HWHM）、Δ_p、Δ_s（= T02 の Δ）、ε(t)、x_ss = √(P_pk/P_th,CW)、T_p。「二重共振の OPO のポンプのキャビティはポンプ経路の区間で、定常の位相 −arctan(Δ_p/κ_p)」。
   T02 の x = √(P/P_th) は CW の場合、パルスでは x_ss を使う。
2. **conventions §6 の σ² の積分の上限**: 「光源のスペクトルが平らな範囲」は OPO の ns モードでは ~κ_c（数十 MHz）。それより速い位相雑音は 1 パルスの中で角度を変え、遅いジッタとは同じに効かない。OPO のときはモードの帯域を上限として明記する。
3. **D07・D43**: OPO のプリセットでは τ_p を掃引しない。代わりに (Δν_c, T_p, κ_p, x_ss, ρ)。f_rep の上限は裾で決まる（§5.3）。
4. **D42**: 比較の構成に「(d) Liu 型: OPO を透過したポンプを 2ω0 の参照光とホモダインし、ポンプ光路の PZT で固定」を加える。残差は信号の離調の回転（≈ −0.55 Δ_s/κ_c）。
   また「単一光源を時間で分ける構成では、光源のスクイーズ角の遅いドリフトは指標 (ii) には 4 sin²(πf(τ_k − τ_0)) でしか入らず、指標 (i) にはそのまま入る」ことを明記する。
5. **D46**: ゲートの長さはモードの長さ（~25 ns）以上。ゲートあたりの迷光は T10 の 1 ns の値の約 25 倍。
6. **T11 の訂正**: Borealis の時間配分は「ロック光 65 µs、squeezed 光 35 µs」（T11 §2 は逆）。振り分けの透過は SI では > 90%（本文の「約 15% の損失」との違いは定義の差の可能性）。
   T11 §3.3 の「τ_p ≪ T_rep は常に成り立つ」は OPO では成り立たない（§5.3）。
7. **T02 §6.4 の更新**: 「CW 励起の OPO を時間ビンで切り出すなら 1/(2κ) ≪ T_bin」は、パルス励起では「ε_tail(T_rep) ≈ 8.8 e^{−2κ_c T_rep} ≪ 1」。
8. **T12 §4.2**: OPO のプリセットでは η_up を ρ × η_filter × η_coupling に置き換える（Borealis 型 約 0.76、Liu 型 約 0.92–0.94）。
9. **survey の構成 F** の光源の行に「パルス励起 OPO（D43）」を加える（Phase 2 の統合担当へ）。

## 10. ユーザの決定が必要な点（推測しない）

1. D43 の確認（Open のまま）: 「Xanadu 2022 の Supplementary と同等の OPA」が §5.1 左列の OPO を指すか。
2. 将来の OPO の線幅（κ_c）、T_p、f_rep、N_b。Borealis 型（66 MHz）なら f_rep ≲ 45 MHz、32 ポートをパルスごとに PC で分けるなら 6 MHz（L_31 ≈ 1.06 km）。
3. OPO のロック光の方式: Borealis 型（直交偏光・+7 GHz・CW、結晶の温度で二重共振）か Liu 型（同じ偏光・同じ周波数・時間をずらす）か。
4. ポンプ位相のロックを入れるか（Liu 型の透過ポンプのロック、T14 の案 A/B、または入れない Borealis 型）。指標 (i) を使うか (ii) だけかで要否が変わる。
5. FSR フィルタ（FBG の組か MZI の段か）、その透過。
6. LO の作り方: 最適モードに整形したパルス LO か、CW の LO＋速い BHD でデジタルの重みか（D45）。
7. 逆向き probe の周波数: 同じ周波数（Liu）か、ずらす（Borealis の +7 GHz、FBG で除く）か（D40 の選択肢の中で）。
8. OPO を複数台使うか（Liu は 4 台。複数台なら OPO ごとのポンプ位相の絶対的なロックが要る）。

## 11. 文献

読んだもの:

- L. S. Madsen, F. Laudenbach, M. F. Askarani, F. Rortais, T. Vincent, J. F. F. Bulmer, F. M. Miatto, L. Neuhaus, L. G. Helt, M. J. Collins, A. E. Lita, T. Gerrits, S. W. Nam, V. D. Vaidya, M. Menotti, I. Dhand, Z. Vernon, N. Quesada, J. Lavoie,
  "Quantum computational advantage with a programmable photonic processor," Nature (2022)。DOI: 10.1038/s41586-022-04725-x（`refs/xanadu_2022.pdf` の Supplementary Information の表紙から）。
  **ローカルの PDF は SI だけ**（§2.1–2.6 を全部、Fig. S6・S7 は画像で読んだ）。巻・ページ（606, 75–81）は T11 が PMC から引いたもので、この担当は確認していない **[UNVERIFIED]**。本文の値（ループの位相雑音など）は T11 経由。
- H.-L. Liu, H. Su, Y.-H. Deng, S.-Q. Gong, et al., "Gaussian boson sampling with 1,024 squeezed states in 8,176 modes," Nature **653** (2026)。DOI: 10.1038/s41586-026-10523-6（SI の表紙）。
  本文（Methods の「Squeezed light sources and spectral filtering」「Phase stabilization of the system」、Extended Data Fig. 1）と SI §1.1–1.5、§2、§4.2 を読んだ。ページ（本文の柱に 692）と著者の全リストは確認していない **[UNVERIFIED]**。
- T02 で確認済みの入出力理論（Gardiner & Collett 1985、Phys. Rev. A 31, 3761、DOI 10.1103/PhysRevA.31.3761）。§2 の時間領域の式はこの担当がそこから導き、CW の極限で T02 の式と数値で照合した（V-P1）。

書誌だけ（SI の文献リストから引き写し、本文は読んでいない **[UNVERIFIED]**）:

- M. Mehmet et al., "Squeezed light at 1550 nm with a quantum noise reduction of 12.3 dB," Opt. Express 19, 25763–25772 (2011)（Borealis のヘミリシック OPO の文献 [10]）。
- V. D. Vaidya et al., "Broadband quadrature-squeezed vacuum and nonclassical photon number correlations from a nanophotonic device," Sci. Adv. 6, eaba9186 (2020)（単一時間モードの条件の文献 [11]）。
- S. Takeda, K. Takase, A. Furusawa, "On-demand photonic entanglement synthesizer," Sci. Adv. 5, eaaw4530 (2019)（Pockels セルの文献 [12]）。
- L. Neuhaus et al., "PyRPL (Python Red Pitaya Lockbox) …," European Quantum Electronics Conference, EA_P_8 (2017)（ループのロックの文献 [13]）。
- Liu SI の文献 [1]（時間順序の効果）と [2]（Heisenberg 方程式による数値モデル）: 書誌をこの担当は確認していない。

値として **UNVERIFIED** のもの: KTP の dn/dT（Δ(dn/dT) は仮の範囲）、P_th,CW、Borealis の r の参照面（§3 の推定）、二重共振の FSR の分散による r_m の打ち切り。
