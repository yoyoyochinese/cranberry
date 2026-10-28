# T07 ファイバの量子雑音層と散乱雑音（fiber-quantum-noise-scattering）

- Phase: 1（deep-dive、T07）
- 日付: 2026-10-28
- 担当する効果: E32（自発ラマン雑音・熱フォノン）、E33（GAWBS）、E34（SBS 閾値）、E35（自発ブリルアン後方散乱）、
  prompt §5.1 の量子雑音層の方式（線形化 vs truncated Wigner）、prompt §7-4・§7-6
- 入力: `docs/squeezing_sim_prompt.md`、`DECISIONS.md`（D01–D33）、`notes/conventions.md`、`notes/survey.md`（E32–E35、T07 行）、
  Ref C §1・§2、Liu et al. 2026（Methods「Phase stabilization of the system」）、下記の一次文献
- 値の印は survey §0.3 に従う: **[AE]** = この担当の概算（未検証）、**[UNVERIFIED]** = 原典で確認していない値・記述、**[Ref X]**・**[文献名]** = 原典で読んだ値。
- 概算の計算は使い捨てスクリプトで行った（コミットしない）。式と入力はすべて本文に書いたので再現できる。

---

## 0. 結論（先に）

1. **量子雑音層は「線形化（古典解のまわりの Jacobian ＋ 共分散）」を推奨する。** truncated Wigner（TW）は、明るい LO パルスの φ_NL が大きい場合の
   交差検証用の任意ブロックとして残す（§2.5）。理由: ファイバ中の squeezed 光は squeezed vacuum で平均場がほぼ 0（φ_NL ~ 10⁻⁷ rad、§2.5）なので、
   squeezed 光の経路の写像は線形であり、ガウス状態の共分散計算は近似なしに正しい。このトピックの雑音（損失、ラマン、GAWBS、ブリルアン）は
   どれも「古典的な明るい光 × 揺らぎ場」なので、線形化の範囲で入れられる。さらに、決定論的に計算でき（Phase 6 の掃引・BoTorch が Monte Carlo 誤差と戦わない）、
   §7-4 と symplectic 性を機械精度で検証できる。
2. **probe がない構成（squeezed 光だけが入るファイバ、LO は別ファイバ、D09）では、T07 の効果は無視できる。**
   GAWBS による squeezed 光の位相ジッタは 63 m で 約 0.1–0.2 mrad、630 m で 約 0.3–0.65 mrad rms [AE]。−10 dB/+20 dB の状態で squeezing の劣化は 10⁻³ dB 未満。
   squeezed vacuum 自身によるラマン・ブリルアン散乱は、平均光子数がスロットあたり ~2 なので無視できる。
3. **同じファイバを逆向きに伝搬する CW probe（Liu 型、D30）は、自発ブリルアン後方散乱が squeezed 光と同じ向き・GBS 側に光子を入れる。**
   1 mW、1 ps のモードあたり N_B ≈ 4.5×10⁻³（63 m）、4.5×10⁻²（630 m）[AE]。それぞれ −10 dB → −9.6 dB、−7.2 dB。
   10.8 GHz ずれているが、ps パルスの squeezed 光のスペクトル幅（~0.4 THz）の内側なので、スペクトルのフィルタでは取れない。
   後方散乱は往復時間 2n_gL/c にわたって広がるので、**probe をパルスにしてスロットとずらしても除けない（デューティ比の分だけ減る）**。
   同じ周波数の Rayleigh 後方散乱（T10、E44）はこの約 12 倍 [AE]。
   −10 dB での劣化を 0.05 dB 以下に抑えるには、ブリルアンだけで P_pr·D ≲ 130 µW（63 m、1 ps）/ 13 µW（630 m、1 ps）、Rayleigh を含めると約 1/12 になる（§6.3）。
4. **自発ラマン雑音は、probe が squeezed 光と同じ周波数なら小さい。** 1 mW・1 ps で N_R ≈ 9×10⁻⁶（63 m）、9×10⁻⁵（630 m）、295 K [AE]。
   Ω → 0 では g_R(Ω) ∝ Ω と n_th ∝ 1/Ω が打ち消し合い、熱による雑音密度は有限の値 g_R′·k_BT/ħ になる。T → 0 では熱の部分が消え、
   残るのは Stokes 側の自発散乱（∝ g_R(Ω)、Ω → 0 で 0）だけになる（§3、§8 V6）。
5. **GAWBS が問題になるのは、明るい光が squeezed 光と同じファイバを同じ向きに通るときだけである**（Ref C の文献 [27] の構成。Ref C はこれを別ファイバで避けた）。
   CW の LO 1 mW を 63 m 通した場合、R0,7 ピーク（321.75 MHz）の側帯波でショット雑音の約 10 倍 [AE]。
   ps パルスのモードあたりに直すと 約 1.4×10⁻⁴（63 m）と小さいが、10 ns のゲートあたりでは 約 0.65 光子になる（GBS の計数では効く）。
   ピーク周波数は R0m の式から 30.5、82.1、130.7、…、323.3 MHz（クラッド 125 µm）と出て、SMF の測定値 321.75 MHz と 0.5% で合う。80 µm クラッドではすべて ×1.5625 になる（§4）。
6. **SBS 閾値**（Smith 1972 の式、g_B = 5×10⁻¹¹ m/W、A_eff = 80 µm²、いずれも UNVERIFIED）: 63 m で 0.53 W、630 m で 54 mW、1 km（Ref C）で 34 mW。
   CW の LO 15.4 mW（Ref A の値）を 630 m に通すと G = g_B P L_eff/A_eff ≈ 6 で、閾値（G ≈ 21）より下だが、自発散乱の増幅は無視できない。
   パルスの LO・probe では SBS は抑えられる（§5.1）。

---

## 1. 記号（symbols, English）

conventions.md §1 の記号（t, T, z, f, Ω, ω0, λ0）をそのまま使う。追加分:

- `A(z,T)`: classical complex envelope [√W], `|A|²` = power. Field ∝ A e^{−iω0 t} (conventions §4); propagation A → A e^{+iψ}.
- `Â(z,T)`: field operator, `Â = A + √(hν0) â`. `â(z,T)`: fluctuation operator in photon-flux units, `[â(z,T), â†(z,T′)] = δ(T − T′)`.
- `hν0 = ħω0`: photon energy (1.285×10⁻¹⁹ J at 1545.32 nm).
- `q̂, p̂`: quadratures, `â = (q̂ + i p̂)/√2`, vacuum variance 1/2 (conventions §2).
- `x`: real 2N vector of discretized quadratures `(q_1…q_N, p_1…p_N)`, `â_n = √Δt â(T_n)`.
- `V`: symmetrized covariance matrix of x, `V_ij = ⟨{x_i, x_j}⟩/2 − ⟨x_i⟩⟨x_j⟩`. Vacuum: `V = I/2`.
- `K(z)`: real 2N×2N generator of the linearized equation. `M(z2,z1)`: transfer (Jacobian) matrix, `dM/dz = K M`.
- `D(z)`: diffusion (noise) matrix per unit length. `Ω_s = [[0, I_N], [−I_N, 0]]`: symplectic form (`[x_i, x_j] = i Ω_s,ij`).
- `u(T)`: LO temporal mode, `∫|u|² dT = 1`. `ũ(f)`: its spectrum, `∫|ũ|² df = 1`. `τ_eff = |ũ(0)|² = (∫u dT)²`; Gaussian: `τ_eff ≈ 2.13 τ_FWHM` (intensity FWHM).
- `N`: mean number of added (phase-insensitive) noise photons in a mode. Added variance per quadrature = N, i.e. `ΔV/V_SN = 2N`.
- `n(f)`: added photon-flux spectral density [photons s⁻¹ Hz⁻¹] = photons per mode at frequency offset f. `N = ∫ n(f) |ũ(f)|² df`.
- `n_th(Ω, T) = 1/(exp(ħΩ/k_B T) − 1)`: thermal phonon occupation. `k_B`: Boltzmann constant. `T_K`: fiber temperature [K]（式中で T と区別が必要なときだけ）。
- `γ`: Kerr coefficient [1/(W m)]. `f_R`: Raman fraction. `h_R(T)`: Raman response, `∫h_R dT = 1`. `h̃_R(Ω) = ∫ h_R(T) e^{+iΩT} dT` (conventions §8).
- `R(T) = (1 − f_R) δ(T) + f_R h_R(T)`: total nonlinear response.
- `g_R(Ω) = 2 γ f_R |Im h̃_R(Ω)|`: Raman power-gain coefficient [1/(W m)] (`Ω` = pump − scattered angular frequency; Stokes for Ω > 0).
- `g_R′ = dg_R/dΩ|_{Ω→0} = 2 γ f_R ⟨T⟩_h`, `⟨T⟩_h = ∫ T h_R(T) dT`.
- `φ_R(z,T)`, `φ_G(z,T)`: Raman / GAWBS phase-noise fields [rad/m] (phase per unit length).
- `S_φ,j(f)`: one-sided PSD of accumulated phase from mechanism j [rad²/Hz] (conventions §7).
- `g_B`: peak Brillouin gain [m/W]. `ν_B`: Brillouin shift [Hz]. `Δν_B`: Brillouin FWHM [Hz]. `A_eff`: effective area [m²].
- `L_eff = (1 − e^{−αL})/α`: effective length. `α`: power loss coefficient [1/m].
- `G_B = g_B P L_eff / A_eff`: Brillouin single-pass gain exponent.
- `R0,m`, `TR2,m`: radial / torsional-radial guided acoustic modes of the cladding. `f_m`: cut-off (peak) frequency. `Γ_m`: acoustic FWHM [rad/s].
- `γ_OM^{(m)}`: peak opto-mechanical (forward-SBS) nonlinear coefficient of mode m [1/(W m)] (London et al. 2018).
- `V_L`, `V_S`: longitudinal (dilatational) and shear sound velocities. `d`: cladding diameter, `a = d/2`.
- `P_pr`: probe power, `D_pr`: probe duty cycle. `S_Ray`: Rayleigh capture fraction, `α_s`: Rayleigh scattering loss (T10).
- `N_traj`: number of truncated-Wigner trajectories.

---

## 2. 量子雑音層の方式（prompt §5.1）: 比較と推奨

### 2.1 古典の方程式（T06 と共通、規約での形）

conventions §8（Agrawal の Fourier 規約、∂_T → −iΩ）と §4（A → A e^{+iψ}）での GNLSE:

```
∂_z A = −(α/2) A + i Σ_{k≥2} (β_k / k!) (i ∂_T)^k A
        + i γ (1 + i τ_sh ∂_T) [ A(z,T) ∫ R(T′) |A(z, T − T′)|² dT′ ]
```

周波数領域では分散の演算子は `D̂(Ω) = i Σ β_k Ω^k / k! − α/2` で、A → A e^{+iβ(ω)z} の符号と一致する（τ_sh = 1/ω0）。T06 が実装する。

### 2.2 方式 1: 線形化（Jacobian ＋ 共分散）

**揺らぎの方程式。** Â = A + √(hν0) â を入れて â について 1 次まで残す（Haus & Lai 1990 の線形化の考え方。Raman と雑音の項は Drummond & Corney 2001 を線形化したもの）:

```
∂_z â = −(α/2) â + i Σ (β_k/k!) (i∂_T)^k â
        + i γ (1 − f_R) [ 2|A|² â + A² â† ]
        + i γ f_R [ (h_R ∗ |A|²) â + A · (h_R ∗ (A* â + A â†)) ]
        + √α ŝ(z,T)                                  (loss: vacuum input, [ŝ, ŝ†] = δ(z−z′)δ(T−T′))
        + i A(z,T) φ_R(z,T) / √(hν0)                 (Raman phonon noise)
        + i A(z,T) φ_G(z,T) / √(hν0)                 (GAWBS thermal phase noise)
        + ŝ_bs(z,T)                                  (back-scattered noise from a counter-propagating bright field, §5.2)
```

（`∗` は T についての畳み込み。shock 項 τ_sh は同じ形で入るので省略。）

雑音場の相関（対称順序、片側にしない「2 側」スペクトル、単位長さあたり）:

```
⟨φ_R(z,T) φ_R(z′,T′)⟩ → spectrum  S_R^(2)(Ω) = hν0 · g_R(|Ω|) · [ n_th(|Ω|) + 1/2 ] · δ(z − z′)     [rad² s / m]
⟨φ_G(z,T) φ_G(z′,T′)⟩ → spectrum  S_G^(2)(Ω) = hν0 · g_G(|Ω|) · [ n_th(|Ω|) + 1/2 ] · δ(z − z′)
```

- S_R は Drummond & Corney 2001 の Eq. (77)（Wigner の Raman 雑音 `⟨Γ^R Γ^R*⟩ = α^R(|Ω|)[n_th + 1/2]/n̄`、`α^R = 2|Im h̃|`）を物理単位に直したもの。
  検算: 明るい場 A に掛かると、側帯波の光子束の密度は |A|² S_R^(2)/(hν0) = P g_R (n_th + 1/2)。真空の 1/2 と合わせて正規順序で Stokes 側 n_th + 1、anti-Stokes 側 n_th になり、自発ラマン散乱の標準形（§3）に一致する。
- GAWBS は同じ形で、g_G(Ω) は前方 SBS（F-SBS）の利得スペクトル Σ_m 2 γ_OM^{(m)} · Lorentz_m(Ω)。「γ_OM → 利得 2γ_OM」の因子 2 は London et al. の定義を原文で確認していない [UNVERIFIED]。
  GAWBS の周波数（< 1 GHz）では n_th ≈ k_BT/(hf) ≈ 2×10⁴ ≫ 1 で、古典の熱雑音として扱ってよい。
- z について δ 相関とするのは、音響モードの軸方向の相関長が L より十分短いという仮定（GAWBS の雑音パワーが L に比例するという Shelby et al. 1985 の観測と整合 [UNVERIFIED: 原文未読、Elser 2006 の要約による]）。

**離散化と共分散。** 時間格子 T_n（N 点、Δt）で â_n = √Δt â(T_n)、x = (q, p) とすると

```
dx/dz = K(z) x + ξ(z),     ⟨ξ(z) ξ(z′)ᵀ⟩_sym = D(z) δ(z − z′)
dV/dz = K V + V Kᵀ + D,    V(L) = M V(0) Mᵀ + ∫_0^L M(L,z) D(z) M(L,z)ᵀ dz
```

- 損失: K ⊃ −(α/2) I_{2N}、D ⊃ (α/2) I_{2N}。このとき V = I/2 は dV/dz = 0 を満たす（vacuum stays vacuum）。
- â† が入るので複素線形ではない。実 2N ベクトル（または (â, â†) の組）で扱う。分散の相互作用表示は â に e^{D̂h}、â† にその複素共役をかける。
- Raman・GAWBS の D は、ξ = (√2 Re, √2 Im) of [i A φ √Δt / √(hν0)] の相関から作る。A ≈ 0 の squeezed vacuum の経路では 0 になる（自分自身の散乱は無視できる、§2.5）。

**測る量。** LO モード u、角度 θ の直交位相 x̂_θ[u] = ∫ [u*(T) â(T) e^{−iθ} + h.c.]/√2 dT は、実ベクトル w で x̂_θ = wᵀx と書ける。

```
Var(x̂_θ[u]) = w_inᵀ V(0) w_in + ∫_0^L w(z)ᵀ D(z) w(z) dz,     w(z) = M(L,z)ᵀ w,  w_in = M(L,0)ᵀ w
```

w(z) は**随伴（adjoint）方程式** dw/dz = −K(z)ᵀ w を z = L から 0 へ解けば得られる。**1 つの LO モードと 1 つの θ に対し随伴の伝搬 1 回**で済み、2N×2N の M 全体は要らない。
θ を回した曲線（R±、反スクイーズ角）には θ = 0, π/2 と交差項の 3 つの 2 次形式があれば足りる。

**Jacobian の作り方**（優先順）:
- (a) 接線（tangent-linear）RK4IP: 古典と同じ積分器で線形化した右辺を積む。古典解 A(z,T) は z ごとに保存するか、同時に積分する。
- (b) 有限差分 Jacobian（中心差分、摂動 ε ~ 10⁻⁶|A|）: (a) の検算用（§8 V2）。
- (c) 自動微分（JAX など）は新しい依存になるので使わない。

**CW のとき（LASER_MODE = CW、§5.1）。** GNLSE は不要で、側帯波 ±Ω の組ごとの 4×4（または 2×2）の解析的な伝達になる。
損失は η = e^{−αL}、GAWBS・ラマン・ブリルアンは §3–§5 の n(f) を加える雑音として入れる。

### 2.3 方式 2: truncated Wigner の確率 GNLSE

Drummond & Corney 2001, Eq. (64)（φ* 形）を物理単位に直す:

```
∂_z Ψ = [GNLSE right-hand side of §2.1 evaluated with Ψ]  + Γ(z,T)  + i Ψ(z,T) φ_R^W(z,T)
Ψ(0,T) = A_in(T) + √(hν0) w(T),   ⟨w(T) w*(T′)⟩ = δ(T − T′)/2   (vacuum; squeezed input: sample from V(0))
⟨Γ Γ*⟩ = (α/2) hν0 δ(z−z′) δ(T−T′),   φ_R^W : spectrum hν0 g_R(|Ω|)[n_th + 1/2]  (D&C Eq. 75, 77)
```

対称順序の期待値がそのまま出るので、Var(x̂_θ) は規約の V（真空 1/2）と直接比べられる。光子数は ⟨|Ψ|²⟩/hν0 からモード数/2 を引く（D&C）。

**文献の規約からの変換**: D&C は ζ = z/x0、τ = (t − z/v)/t0、φ = Ψ √(v t0 / n̄)、分散の項 ±(i/2)∂²φ/∂τ²（+ が異常分散）。
これは Agrawal 形の −(iβ2/2)∂_T² A（β2 < 0 で +i|β2|/2）と符号が一致する。D&C の Fourier 変換は場に (1/√2π)∫dτ e^{+iΩτ}、応答関数に ∫dτ e^{+iΩτ}（係数なし）で、
符号は conventions §8 と同じ、場の規格化だけ √(2π) が違う。

### 2.4 比較

| 項目 | 線形化（Jacobian ＋ 共分散） | truncated Wigner |
|---|---|---|
| 計算量 | 古典 1 回 ＋ 随伴 2–3 回 / LO モード（M 全体なら接線 2N 回）。決定論的 | 古典と同じ計算を N_traj 回。統計誤差 √(2/N_traj)（分散の相対誤差）。0.1 dB（2.3%）の分解能に N_traj ≈ 4×10³、余裕をみて 10⁴ |
| 適用範囲 | 明るい部分の揺らぎ ≪ 平均場。squeezed vacuum のような線形の伝搬では厳密。Kerr squeezing の「曲がり」（位相拡散で反スクイーズが squeezing 側に回り込む領域）や、光子数が小さいソリトンでは破れる | 光子数（モードあたりの占有）が大きいこと（3 次・4 次微分の切り捨て、D&C）。1 次を超える揺らぎの非線形性を含む。占有が小さいモードで長く伝搬すると破れる（Sinatra et al. 2002）。格子の高周波の真空雑音が非線形項に入る（帯域の打ち切りが必要） |
| ラマン雑音 | D への加法的なガウス項（S_R） | 乗法的な雑音 iΨφ_R^W |
| GAWBS | 同上（古典熱雑音、S_G） | 同上 |
| ブリルアン・Rayleigh 後方散乱 | 加法的な位相によらない雑音 n(f)（§5.2） | 同じ（逆向きの場はどちらの方式でも同じ GNLSE の枠に入らない） |
| §7-4 vacuum stays vacuum | 厳密（V = I/2 が丸め誤差まで保たれる） | 統計的にのみ（許容幅 ∝ 1/√N_traj） |
| symplectic 性 | M Ω_s Mᵀ = Ω_s を直接検査できる | 直接は検査できない |
| CW | 側帯波ごとの解析的な行列 | 不要 |
| 出力 | 共分散・伝達行列・雑音表がそのまま出る | 標本から推定 |

### 2.5 推奨: 線形化を既定、TW は任意の検証ブロック

理由:

1. **squeezed 光の経路は線形。** 10 dB の squeezed vacuum でも平均光子数は sinh²r ≈ 2（r = 1.15）/ パルス。1 ps なら尖頭パワー 2 hν0/τ_p ≈ 2.6×10⁻⁷ W、
   φ_NL = γ P L ≈ 1.3×10⁻³ × 2.6×10⁻⁷ × 630 ≈ 2×10⁻⁷ rad [AE]。1 次を超える項はこの程度で、線形化は事実上厳密。
2. **雑音はすべて「古典の明るい光 × 揺らぎ場」。** 損失、ラマン、GAWBS、ブリルアン・Rayleigh 後方散乱は、線形化の加法的なガウス雑音 D で表せる。
3. **時間ループ段が要るのは共分散・雑音の表。** prompt §5.1 の出力（squeezing/反 squeezing の伝達、過剰雑音スペクトル、モード重なり）はそのまま 2 次形式で出る。
   Phase 6 の掃引・BoTorch で Monte Carlo の揺らぎが目的関数に乗らない。
4. **検証が厳密。** §7-4 と symplectic 性を機械精度で確認でき、§7-3（純損失）も同じ枠で厳密に出る。
5. **TW が要るのは明るい LO・probe 自身の量子雑音で φ_NL が大きいとき**（survey E29: 100 fs、1 mW、63 m で φ_NL ≈ 8 rad [AE]）。
   その場合でも LO の振幅雑音 → SPM → LO 位相雑音は線形化で十分なことが多い。TW は「線形化 vs TW を小さい問題で比べる」検証ブロック（§8 V9）として持ち、既定ではオフにする。

OPA（T01/T03）の中の高利得パルス squeezing は別のトピック。ここではファイバの区間だけを扱い、入力共分散 V(0) を T01/T03 から受け取る。

---

## 3. 自発ラマン散乱と熱フォノン（E32）

### 3.1 式

明るい光（パワー P、周波数 ν_p）から、周波数 ν_p − Ω/2π の squeezed 光のモードへ散乱される光子の密度（同じ向き、長さ L、利得が小さい g_R P L ≪ 1）:

```
n_R(Ω) = g_R(|Ω|) · P · L · { n_th(|Ω|, T) + 1   (Stokes,  Ω > 0: scattered frequency below the bright field)
                               n_th(|Ω|, T)       (anti-Stokes, Ω < 0) }
N_R = ∫ n_R(Ω) |ũ(f)|² df
```

- g_R(Ω) = 2γ f_R |Im h̃_R(Ω)|（Agrawal 形: g_R = (2ω0/c) n2 f_R Im h̃_R を A_eff で割ったもの [UNVERIFIED: 教科書の式、原文未確認]）。
- **Ω → 0 の極限**: Im h̃_R(Ω) ≈ Ω ⟨T⟩_h、n_th ≈ k_BT/(ħΩ) なので、n_R(0) = g_R′ · k_BT/ħ · P L（Stokes・anti-Stokes で同じ、有限）。
  probe が squeezed 光と同じ周波数（Liu 型）のとき、squeezed 光の帯域に入るのはこの値。
- **T → 0**: 熱の部分 ∝ n_th → 0。残るのは Stokes 側の g_R(Ω) P L だけで、これも Ω → 0 で 0。帯域 B ≪ k_BT/h のモードでは N_R(T→0)/N_R(295 K) ~ hB/(k_BT)。
- **詳細釣り合い**: anti-Stokes/Stokes = n_th/(n_th + 1) = exp(−ħΩ/k_BT)。
- 自発散乱による光子は位相によらない雑音として Var に N を足す（ΔV/V_SN = 2N）。
- 逆向きの probe の後方ラマン散乱も同じ式で、捕捉率が前方とほぼ同じと仮定する [AE、UNVERIFIED。QKD との共存の実験（Eraerds et al. 2010）は前方と後方の両方を測っている、内容は未読]。

**ラマン応答のモデル**: Blow & Wood 1989 の 1 振動子モデル h_R(T) = (τ1² + τ2²)/(τ1 τ2²) e^{−T/τ2} sin(T/τ1)、τ1 = 12.2 fs、τ2 = 32 fs、f_R = 0.18 [UNVERIFIED: よく引用される値、原文未読]。
これで γ = 1.3 /(W km) のとき g_R のピーク = 0.68 /(W km)（13.1 THz）、⟨T⟩_h = 8.7 fs、g_R′ = 4.1×10⁻¹⁸ /(W m)/(rad/s) [AE、計算]。
**1 THz 以下では 1 振動子モデルは不確か**（boson peak、準弾性散乱: Perlmutter et al. 1990 の主題 [内容 UNVERIFIED]）。
Hollenbeck & Cantrell 2002 の多振動子モデル、Lin & Agrawal 2006 のモデル、または D&C の 10 Lorentz の当てはめを Phase 4 で比べ、Ω → 0 の値の不確かさ（係数 2–3 程度と予想 [AE]）を表に残す。

### 3.2 周波数帯

0 – 約 40 THz（ピーク 13 THz）。squeezed 光のモードの帯域（1 ps で ~0.44 THz、100 fs で ~4.4 THz）と probe の周波数の差の範囲に入る部分だけが効く。
CW の squeezing（Ref A の 3 MHz、Ref C の GHz）では Ω ≈ 0 の値 n_R(0) が白色の雑音として全帯域に入る。

### 3.3 スケーリング

∝ P（明るい光）× L × [n_th(Ω,T) + 1 or n_th]。Ω ≈ 0 では ∝ T。パルス幅: 雑音は THz にわたって平坦なので、モードあたりの N_R は τ_p にほとんどよらない（100 fs で少し変わる）。

### 3.4 概算（295 K、P = 1 mW、λ0 = 1545.32 nm）

n_R(Ω)/(P L) [1/(W m)]: 0.1 THz で 1.5×10⁻⁴（S）/1.5×10⁻⁴（aS）、1 THz で 1.6×10⁻⁴/1.4×10⁻⁴、5 THz で 2.6×10⁻⁴/1.2×10⁻⁴、13.2 THz で 7.7×10⁻⁴/9.0×10⁻⁵ [AE]。

| τ_p | L | N_R（295 K） | N_R（1 mK） |
|---|---|---|---|
| 100 fs | 63 m | 9.6×10⁻⁶ | 1.2×10⁻⁶ |
| 100 fs | 630 m | 9.6×10⁻⁵ | 1.2×10⁻⁵ |
| 1 ps | 63 m | 9.3×10⁻⁶ | 1.1×10⁻⁷ |
| 1 ps | 630 m | 9.3×10⁻⁵ | 1.1×10⁻⁶ |
| 10 ps | 63 m / 630 m | 9.3×10⁻⁶ / 9.3×10⁻⁵ | 1.1×10⁻⁸ / 1.1×10⁻⁷ |

−10 dB に対し 630 m でも劣化 0.02 dB 未満。**probe と squeezed 光の周波数が同じなら自発ラマンは小さい。**
probe の波長を数 THz 以上ずらす案（survey E46、WDM で分離）では、squeezed 光の帯域が Stokes/anti-Stokes の強い部分に入るので、上の表の値の数倍になる（13 THz で 5 倍）。

---

## 4. GAWBS（E33）

### 4.1 機構と式

クラッドの横方向の音響モード（軸方向の波数がほぼ 0 の、カットオフ付近の案内音響波）が熱で揺れ、屈折率を変える。前方にだけ位相整合する（後方は音響の波数 2β が必要で、それは通常の後方 SBS（§5）になる）。

- R0,m（放射状）: 屈折率の等方的な変化 → **位相雑音**（polarized GAWBS）。
- TR2,m（ねじれ・放射の混合）: 位相雑音と**複屈折の揺らぎ → 偏光雑音・直交偏光への散乱**（depolarized GAWBS）（Elser et al. 2006 の Fig. 1 の説明）。

**R0,m のピーク周波数**（自由な円柱、表面に拘束なし）:

```
(1 − α_v²) J0(y_m) − α_v² J2(y_m) = 0,     α_v = V_S / V_L,     f_m = y_m V_L / (2π a)
```

この式は Shelby et al. 1985 の式とされるが原文で確認していない [UNVERIFIED]。下の数値が SMF の測定値と 0.5% で合うことで間接的に確かめた。
V_L = 5996 m/s、V_S = 3740 m/s（溶融石英の典型値 [UNVERIFIED]）で:

| クラッド径 d | R0,m のピーク [MHz]（m = 1, 2, …） |
|---|---|
| 125 µm | 30.5, 82.1, 130.7, 179.0, 227.2, 275.3, **323.3**, 371.3, 419.4, 467.4, 515.4, 563.4 |
| 80 µm | 47.7, 128.2, 204.3, 279.7, 355.0, 430.1, 505.2, 580.2, 655.2, 730.3, 805.3, 880.2 |

- 測定値: 標準 SMF（125 µm、アクリレート 2 層被覆）で R0,7 = 321.75 MHz、線幅 4.3 MHz、γ_OM^{(0,7)} ≈ 2 /(W km)、F-SBS が最も強いのは R0,7（London et al. 2018、Table 1 と §2.1）。
- **添字の流儀が文献で違う**: London et al. は 30.5 MHz を m = 1 と数えて 321.75 MHz を R0,7 と呼ぶ。Elser et al. は最低次を (0,0) と書く。実装ではどちらかに決めて記録する。
- 間隔は大きい m で V_L/(2a)（125 µm で 48 MHz）。f_m ∝ 1/d。
- TR2,m の周波数は 2×2 の行列式の根（Shelby 1985）で、ここでは計算していない。Elser et al. は 80 µm の PM ファイバで TR2,m（m = 0…12）を 400 MHz 以下に観測。Phase 4 で式を原文から取って実装する。

**PM と SMF の違い**:
- PANDA などの応力付与部はクラッドの円対称を崩し、モードの分裂・追加のピーク・R0m と TR2m の混合を起こす。
  Elser et al. 2006 は 80 µm の PM ファイバ（Fibercore HB800G）で「理論にないピークを PM の幾何学的な非対称によると推定」と書いている。
  Nishizawa et al.（Opt. Lett. 1994、Opt. Rev. 1996）が PM ファイバの GAWBS を扱う [内容 UNVERIFIED、書誌のみ確認]。
- PM ファイバでも depolarized GAWBS による 2 軸間の散乱は起こる（Nishizawa 1994 が PM ファイバでループミラーを使って観測、[内容 UNVERIFIED]）。
  2 軸の伝搬定数の差 Δβ ≈ 2π B/λ に対応する音響の波数は小さく、カットオフ付近の音響モードの周波数はほとんど変わらないので、位相整合は崩れないと見積もる [AE]。
- 細径（80 µm）クラッドの PM ファイバではピークが 1.56 倍の周波数に上がり、雑音の結合も変わる（コア ↔ クラッド表面の距離が近い）。

**位相雑音の大きさ（揺動散逸からの概算）** [AE]: §2.2 の S_G から、モード m の片側 PSD のピークと分散は

```
S_φ,G^(1)(f_m) ≈ 4 k_B T ω0 γ_OM^{(m)} L / Ω_m,          σ_G,m² ≈ S_φ,G^(1)(f_m) · (π/2) · Γ_m/(2π)
```

R0,7（321.75 MHz、4.3 MHz、2 /(W km)）、295 K で S_φ ≈ 1.96×10⁻¹⁷ rad²/Hz/m、σ² ≈ 1.33×10⁻¹⁰ rad²/m。
全モード（R0m ＋ TR2m）の合計は R0,7 だけの 1–5 倍と見る [AE、根拠は弱い]。

**整合の確認** [AE]:
- Elser et al. 2006: 8 m の PM ファイバ（80 µm、MFD 4.5 µm、1064 nm）、0.9 mW、可視度 70%、検出効率 86% で、位相雑音のピークがショット雑音の上 1 dB 前後（図の縦軸範囲 0–1.25 dB から読み取った目安 [UNVERIFIED]）。
  上の式で 8 m・0.9 mW の P·S_φ/(hν) ≈ 1.1 → 効率 0.42 を掛けて 0.46（1.6 dB）。ファイバと波長が違うが、桁は合う。
- Takefushi et al. 2020: SSMF で 1024-QAM の誤りなし伝送距離が GAWBS で 340 km に制限される。1024-QAM の許容位相誤差を ~0.01 rad rms と仮定すると σ² ≈ 3×10⁻¹⁰ rad²/m [AE、仮定込み]。上の R0,7 だけの値の約 2 倍で整合。

**明るい光が同じ向きに通るとき、squeezed 光のモードへ散乱される光子**:

```
n_G(f) = P · S_φ,G^(1)(f) / (2 hν0)       (per sideband; phase-modulation sidebands, correlated between ±f)
N_G(pulsed mode) ≈ τ_eff · P σ_G² / (hν0)     (narrow-band noise, τ_p ≪ 1/f_m)
```

### 4.2 周波数帯

R0m: 約 30 MHz – 1 GHz の離散ピーク（125 µm）。強いのは 100–400 MHz 付近。TR2m: 数十 MHz – GHz。ピークの幅は被覆による減衰で数 MHz（London et al. の R0,7 で 4.3 MHz）。
- 時間ループ段（kHz 帯）のロックには入らない。ロック帯域より上の位相雑音として「θ_rms に足す定数」で十分。
- パルスの squeezed 光（τ_p ≪ 1 ns）にはパルスごとのランダムな位相（相関時間 ~1/Γ ~ 70 ns）。f_rep = 10–100 MHz では隣のパルスと部分的に相関する。

### 4.3 スケーリング

σ_G² ∝ L · T（古典の熱、hf_m/k_B ≈ 15 mK ≪ 295 K）。f_m ∝ 1/d。モードあたりの散乱光子 ∝ P_bright · L · T · τ_eff（パルス）。ファイバの A_eff が小さいほど広帯域に広がる（Takefushi et al. 2020 の要旨）。

### 4.4 概算

| 状況 | 値 [AE] | squeezing への影響（−10 dB / +20 dB） |
|---|---|---|
| squeezed 光自身の位相ジッタ、63 m | σ_G ≈ 0.09（R0,7）– 0.2 mrad | 10⁻⁴ dB 未満 |
| 同、630 m | σ_G ≈ 0.29 – 0.65 mrad | 2×10⁻³ dB 未満 |
| Ref C、1 km、別ファイバ | σ_G ≈ 0.36 – 0.8 mrad | Ref C の全体 0.1 rad に対し無視できる |
| 同じ向きの CW LO 1 mW、63 m（R0,7 ピークの側帯波で） | P S_φ/(hν) ≈ 10（ショット雑音の約 10 倍） | CW での測定は不可能（Ref C の文献 [27] が LO をチョップした理由と整合） |
| 同、1 km（Ref C の文献 [27] の型） | ≈ 150 | 同上 |
| 同じ向きの probe 1 mW、63 m、1 ps モードあたり | N_G ≈ 1.4×10⁻⁴（×5 で 7×10⁻⁴） | −10 dB → −9.99 / −9.94 dB |
| 同、630 m、10 ps モードあたり | N_G ≈ 1.4×10⁻²（×5 で 7×10⁻²） | −10 dB → −8.9 / −6.2 dB |
| 同、10 ns の計数ゲートあたり（GBS） | 0.65 光子（63 m）– 33 光子（630 m、×5） | GBS の許容量は Open（D28） |

注意:
- 同じ向き・同じ偏光の probe の R0m 散乱は squeezed 光の偏光に入る。**同じ向き・直交偏光の probe では TR2m（depolarized）だけが squeezed 光の偏光に入る**。大きさは同程度と見る（Elser の Fig. 4・5 で位相雑音と偏光雑音が同程度 [目安]）。
- **同じ向きの probe はパルスにして時間をずらせば GAWBS の散乱光も一緒にずれる**（GAWBS は遅い変調なので、散乱光は probe パルスの時間窓に留まる）。逆向きの probe（D30）では GAWBS は squeezed 光の向きに散乱しない。
- Ref D（CW 3.33 mW、PM の MZI）: GAWBS は 30 MHz 以上の位相雑音で、DSC1 の帯域（DC–100 kHz）の外。σ_G ~ 0.1 mrad 程度で、ロックへの影響はない（既定でオフ）。

---

## 5. ブリルアン散乱（E34、E35）

### 5.1 SBS の閾値（E34）

Smith 1972 の臨界パワー（CW、狭線幅、後方 Stokes が入力の程度に育つ条件）:

```
P_th ≈ 21 · A_eff / (g_B · L_eff)        (G_B = g_B P L_eff / A_eff ≈ 21 at threshold)
```

- 線幅 Δν_L の光では g_B が Δν_B/(Δν_B + Δν_L) に下がる（Lorentz の重ね合わせ、[UNVERIFIED: Agrawal・Kobyakov 2010 にある形、原文未確認]）。
- パラメータ: g_B = 5×10⁻¹¹ m/W、A_eff = 80 µm²、損失 0.2 dB/km（Ref C の SMF）。g_B と A_eff は典型値 [UNVERIFIED]、PM ファイバの値は未確認（D23 の代用）。

| L | L_eff | P_th | CW LO 15.4 mW（Ref A の値）の G_B | probe 1 mW の G_B |
|---|---|---|---|---|
| 63 m | 62.9 m | 0.53 W | 0.61 | 0.039 |
| 630 m | 621 m | 54 mW | 6.0 | 0.39 |
| 1 km（Ref C） | 978 m | 34 mW | 9.6 | 0.63 |

- 630 m・1 km に CW の LO を十数 mW 通すと、閾値より下でも自発 Stokes が e^{G}/G 程度増幅される。後方へ出るので squeezed 光の向きには来ないが、
  LO の透過パワーの揺らぎ（閾値の近くで SBS が作る RIN）と後方反射によるレーザーへの戻りに注意（T05/T13 に渡す。G_B > 5 で警告を出す）。
- パルス: パルス幅が音響の寿命 1/(πΔν_B) ≈ 10 ns より十分短いと、SBS の利得は減る。繰り返し f_rep > Δν_B の櫛では、櫛の 1 本ごとのパワー
  （平均パワーより十分小さい）で考える [AE]。f_rep = 10 MHz < Δν_B（~30 MHz）では数本が同じ Brillouin 線幅に入る。ps パルスの LO・probe では実質的に問題にならない [AE]。
- Ref D（3.33 mW、PM ファイバ数 m ＋ 915B の巻線）: 閾値に遠く及ばない。

### 5.2 自発ブリルアン後方散乱（E35）: 逆向きの probe から squeezed 光の向きへ

逆向きに伝搬する probe（D30、Liu 型）は、後方の自発ブリルアン散乱で squeezed 光と**同じ向き**の光を作る。散乱光は GBS 側（probe の入口）から出る。

```
Φ_B = (g_B / A_eff) · P_pr · L_eff · (π/2) Δν_B · [ (n_th(ν_B) + 1) + n_th(ν_B) ]     [photons/s]   (Stokes at ν − ν_B, anti-Stokes at ν + ν_B)
ν_B = 2 n_eff V_A / λ0,   n_th(ν_B) ≈ k_B T / (h ν_B) ≈ 570  (295 K, 10.8 GHz)
N_B(mode) = Φ_B · τ_eff · exp[−4 ln2 (ν_B / Δν_u)²]       (Gaussian mode with spectral FWHM Δν_u = 0.441/τ_FWHM)
```

- 「自発放出 = 1 モードあたり n_th + 1（Stokes）の雑音が利得で増幅される」という扱い（Boyd et al. 1990 の雑音起動の考え方）で、G_B ≪ 1 の場合。
  G_B が 1 に近いときは Φ_B に (e^{G_B} − 1)/G_B を掛ける [AE]。
- ν_B ≈ 10.8 GHz（SMF、1550 nm）、Δν_B ≈ 30 MHz [どちらも典型値 UNVERIFIED、Kobyakov et al. 2010 で要確認]。
- 単位長さ・単位パワーあたりの散乱割合 ≈ 4.3×10⁻⁹ /m [AE]。Rayleigh の捕捉割合 S_Ray α_s ≈ 5.3×10⁻⁸ /m（survey の値、UNVERIFIED）の約 1/12。
  Rayleigh と Brillouin の比（Landau–Placzek 比）の典型値 ~20 前後と同じ桁 [AE、比較の数値は UNVERIFIED]。
- **時間ゲートで除けない**: 位置 z で散乱した光は、probe が入った時刻から 2 n_g (L − z)/c 遅れて GBS 側に出る。散乱光は 0 – 2n_gL/c に広がる（OTDR と同じ）ので、
  probe をパルスにしてスロットからずらしても、squeezed パルスに重なる量はデューティ比 D_pr の分しか減らない。
- **周波数で除けない（パルスの場合）**: 10.8 GHz は ps パルスの帯域（1 ps で 441 GHz）の内側。CW の狭帯域 squeezing（OPO の線幅 < 10 GHz）なら一部はフィルタできる。

### 5.3 周波数帯とスケーリング

±ν_B（約 10.8 GHz）、幅 Δν_B（数十 MHz）。∝ P_pr · D_pr · L_eff · T（n_th ∝ T）· τ_eff。SBS の閾値は ∝ A_eff/(g_B L_eff)。

---

## 6. 概算のまとめ（自分たちの構成）

共通の値: λ0 = 1545.32 nm（D07 による Ref A の値の代用）、T = 295 K（D27）、γ = 1.3 /(W km)（London et al. 2018 の SMF の値。PM ファイバの値は未確認で同じ値を代用）、
損失 0.2 dB/km（Ref C の SMF の値で代用）、A_eff = 80 µm²・g_B = 5×10⁻¹¹ m/W・ν_B = 10.8 GHz・Δν_B = 30 MHz（典型値、UNVERIFIED）。
probe の値は **D28 で Open**。ここでは P_pr = 1 mW の CW を仮の値とし、すべて P_pr·D_pr に比例する。τ_p は D07 の 100 fs – 10 ps。
L_k ≤ 63 m（f_rep = 100 MHz）と ≤ 630 m（10 MHz）（survey §3）。依存トピック（T06 の γ・A_eff、T09 の PM の値、T10 の Rayleigh）の値が出たら差し替える。

### 6.1 Ref D（D02）

CW 3.33 mW、PM ファイバの MZI。T07 の効果は全部無視できる（GAWBS は 30 MHz 以上で σ ~ 0.1 mrad、SBS の閾値に遠い、量子雑音の評価は不要）。**既定でオフ。**

### 6.2 将来の構成 F: probe なし（squeezed 光だけが通るファイバ、LO は別ファイバ）

| 効果 | 63 m | 630 m |
|---|---|---|
| GAWBS 位相ジッタ（squeezed 光、ロック帯域の外） | 0.09–0.2 mrad | 0.29–0.65 mrad |
| 同、LO 側のファイバ（同じ長さ） | 同じ値（独立） | 同じ値 |
| ラマン・ブリルアン（squeezed vacuum 自身） | 無視できる（φ_NL ~ 10⁻⁸） | 無視できる（~10⁻⁷） |
| SBS（LO を CW 15.4 mW で送る場合） | G_B = 0.6 | G_B = 6（警告） |

GAWBS は LO 側と squeezed 側で独立なので、相対位相のジッタは √2 倍（63 m で ≲ 0.3 mrad、630 m で ≲ 0.9 mrad）。Ref A の θ = 9 mrad に比べて小さい。

### 6.3 将来の構成 F: 逆向きの CW probe（Liu 型、D30）、P_pr = 1 mW（仮）

| τ_p | L | N_B（ブリルアン） | −10 dB → | N_R（ラマン） | 参考: N_Ray（Rayleigh、T10） | −10 dB →（B ＋ Ray） |
|---|---|---|---|---|---|---|
| 100 fs | 63 m | 4.5×10⁻⁴ | −9.96 dB | 9.6×10⁻⁶ | 5.5×10⁻³ | −9.51 dB |
| 100 fs | 630 m | 4.5×10⁻³ | −9.63 dB | 9.6×10⁻⁵ | 5.5×10⁻² | −6.60 dB |
| 1 ps | 63 m | 4.5×10⁻³ | −9.63 dB | 9.3×10⁻⁶ | 5.5×10⁻² | −6.60 dB |
| 1 ps | 630 m | 4.5×10⁻² | −7.22 dB | 9.3×10⁻⁵ | 0.55 | +1.1 dB |
| 10 ps | 63 m | 3.8×10⁻² | −7.54 dB | 9.3×10⁻⁶ | 0.55 | +1.0 dB |
| 10 ps | 630 m | 0.38 | −0.65 dB | 9.3×10⁻⁵ | 5.5 | +10.7 dB |

（Rayleigh は同じ周波数のコヒーレントな散乱で、ゆっくり揺らぐ変位として入る。ここでは同じ N の雑音として並べた。正しい扱いは T10。）

**−10 dB での劣化を 0.05 dB 以下（N ≤ 5.8×10⁻⁴）にする probe の上限 P_pr·D_pr** [AE]:

| τ_p | L | ブリルアンだけ | Rayleigh だけ |
|---|---|---|---|
| 100 fs | 63 m / 630 m | 1.3 mW / 130 µW | 106 µW / 11 µW |
| 1 ps | 63 m / 630 m | 130 µW / 13 µW | 11 µW / 1.1 µW |
| 10 ps | 63 m / 630 m | 15 µW / 1.5 µW | 1.1 µW / 0.11 µW |

参考: 1 µW の probe のショット雑音で決まる位相の PSD は hν/P ≈ 1.3×10⁻¹³ rad²/Hz、10 kHz の帯域で 約 36 µrad rms [AE、係数は T14 で確認]。
**probe を µW 以下にしてもショット雑音の上ではロックに足りる**。検出器の雑音で足りるかは T13/T14 の問題。

**GBS（光子計数）の 10 ns ゲートあたり**（1 mW CW）: ブリルアン 21 光子（63 m）/ 210 光子（630 m）、Rayleigh 260 / 2600（survey と同じ）。許容量は Open（D28）。

### 6.4 将来の構成 F: 同じ向きの probe（比較用）

- 同じ偏光: GAWBS（R0m）・前方ラマン・XPM が squeezed 光に入る。パルスで時間をずらせば、散乱光は probe と一緒にずれるので大部分を除ける（§4.4）。
- 直交偏光: TR2m（depolarized GAWBS）が squeezed 光の偏光に入る（§4.4 の表）。PER による漏れは T09/T10。

### 6.5 Ref C（1 km SMF × 2 本）

- 別ファイバの LO で GAWBS を避けた（Ref C §1）ことは、§4.4 の「同じ向きの CW LO 1 mW、1 km でショット雑音の ~150 倍」と整合。
- 別ファイバでの GAWBS 位相ジッタ 0.4–0.8 mrad（各ファイバ）は、Ref C の全体 0.1 rad に対して無視できる。
- SBS: P_th ≈ 34 mW。Ref C の LO パワーは本文に書かれていない（grep で確認）ので、G_B は出せない。
- 注: Ref C 本文は「Ref. [27] は LO を同じファイバで送り、GAWBS を抑えるため LO をチョップ、1.6 km 後に 2.5 dB」と書くが、文献リストの [27] は Qin et al.「20 km fiber」（Appl. Sci. 2019）で、距離が合わない。
  番号のずれの可能性がある [UNVERIFIED]。

---

## 7. 効果の一覧（周波数帯とスケーリング）

| 効果 | 周波数帯 | L | P | T | τ_p | 段 | 既定 |
|---|---|---|---|---|---|---|---|
| ファイバ損失の真空雑音 | 全帯域 | e^{−αL} | — | — | — | Q | オン |
| 自発ラマン（明るい光から） | 0–40 THz（Ω ≈ 0 で平坦） | ∝ L | ∝ P_bright | ∝ n_th + 1 / n_th（Ω ≈ 0 で ∝ T） | モードあたりほぼ一定 | Q | オン（probe があるとき） |
| ラマンの熱位相雑音（明るい光自身） | THz | ∝ L | — | 同上 | — | Q | オン（GNLSE の LO・probe） |
| GAWBS 位相ジッタ | 30 MHz – 1 GHz（離散） | ∝ L（分散） | — | ∝ T | パルスごと | Q（表）＋ T（定数の位相分散） | オン |
| GAWBS 散乱光（同じ向きの明るい光から） | 同上 | ∝ L | ∝ P | ∝ T | ∝ τ_eff | Q / O | オン（同じ向きの光があるとき） |
| SBS 閾値・増幅 | — | 1/L_eff | 閾値 | — | パルスで抑制 | O（チェック） | 警告のみ |
| 自発ブリルアン後方散乱 | ±10.8 GHz | ∝ L_eff | ∝ P_pr·D_pr | ∝ T | ∝ τ_eff | Q / O | オン（逆向きの probe があるとき） |

---

## 8. 検証テスト（実装が通すべきもの）

| # | テスト | 合格の基準 |
|---|---|---|
| V1 | §7-4 vacuum stays vacuum: 非線形・損失・雑音なし（分散だけ）で V(0) = I/2 → V(L) | max\|V(L) − I/2\| < 10⁻¹² |
| V1b | 損失だけ（雑音 D_loss 込み） | 同上（損失の真空雑音がちょうど補う） |
| V2 | symplectic 性: 損失・雑音なし、Kerr ＋ 分散 ＋ 明るいパルス（φ_NL = 1、ソリトン次数 1） | ‖M Ω_s Mᵀ − Ω_s‖ / ‖Ω_s‖ < 10⁻⁹。接線 RK4IP の M と中心差分の M の差 < 10⁻⁶。刻みを半分にして誤差が約 1/16（§7-5 と同じ） |
| V3 | 物理性: 損失・雑音あり | V + iΩ_s/2 の最小固有値 ≥ −10⁻¹² |
| V4 | Kerr squeezing の解析解（CW または十分長いパルス、分散なし、線形化）: V = (1/2)[[1, 2φ_NL], [2φ_NL, 1 + 4φ_NL²]]、V_min = (1/2)(√(1+φ_NL²) − φ_NL)² | φ_NL = 0.5、1、2 で −4.18、−7.66、−12.54 dB（反 squeezing は符号反転）に 0.01 dB 以内。Kitagawa & Yamamoto 1986 の線形の極限 |
| V5 | §7-3 純損失: squeezed 入力（−10 dB）、η = e^{−αL} | V_out = ηV_in + (1 − η)/2 を 10⁻¹² で満たす |
| V6 | §7-6 ラマン、T → 0: Ω/2π = 100 GHz、P L = 1 W·m | 熱の部分 g_R·n_th が T = 1 mK で 295 K の値の 10⁻⁶ 未満。T = 1 mK で anti-Stokes = 0、Stokes = g_R(Ω)·P L（295 K の値 1.48×10⁻⁴ に対し 2.39×10⁻⁶ /(W m)、Blow–Wood モデル）。1 ps モード、63 m、1 mW で N_R(1 mK)/N_R(295 K) ≈ 0.012（§3.4） |
| V6b | 詳細釣り合い | aS/S = exp(−ħΩ/k_BT) を相対 10⁻⁶ で（13.2 THz、295 K で 0.117） |
| V6c | ラマン利得のピーク | Blow–Wood、γ = 1.3 /(W km)、f_R = 0.18 で g_R の最大 0.68 /(W km)、13.1 THz（±2%） |
| V7 | §7-6 GAWBS のピーク | 125 µm で R0,m = 30.5, 82.1, 130.7, …, 323.3 MHz（±0.1 MHz、V_L = 5996、V_S = 3740 m/s）。7 番目が測定値 321.75 MHz（London et al.）と 1% 以内。80 µm で全ピークがちょうど 125/80 倍 |
| V7b | GAWBS の分散のスケーリング | L を 2 倍で σ_G² が 2 倍、T を 2 倍で 2 倍（相対 10⁻¹²）、T → 0 で熱の部分（n_th）が 0（対称順序の 1/2 は残るが、295 K の値に対し hf/k_BT ≈ 5×10⁻⁵ の割合） |
| V8 | SBS の閾値 | §5.1 のパラメータで 63 m → 0.534 W、630 m → 54.1 mW、1 km → 34.4 mW（±1%）。G_B > 5 で警告が出る |
| V9 | 線形化 vs TW（任意） | φ_NL = 0.5 の Kerr（V4 と同じ問題）を TW の N_traj = 10⁴ で解き、線形化の V_min と統計誤差の 3σ 以内。光子数を下げると差が出始めることも記録する |
| V10 | 自発ブリルアン後方散乱 | §5.2 のパラメータで散乱割合 4.3×10⁻⁹ /m、1 mW・63 m で Φ_B = 2.1×10⁹ /s（±1%）。N_B は P_pr・L・D_pr に線形、τ_eff に比例 |
| V11 | 位相のない雑音の加わり方 | N を加えたとき、すべての θ で Var(x̂_θ) が N だけ増える（ΔV/V_SN = 2N） |

---

## 9. シミュレータへの要求

### 9.1 ブロック・表

| 名前（案） | 段 | 入力 | 出力 | トグル |
|---|---|---|---|---|
| `fiber_quantum_linearized` | Q | 古典解 A(z,T)（T06）、入力共分散 V(0)（T01/T03）、LO モード u(T)、区間の値（L、α、β_k、γ、f_R、h_R のモデル、T、d、ファイバの種類）、雑音源の表 | V(L) または LO モードの 2 次形式（Var(x̂_θ)、R±、反スクイーズ角のずれ）、伝達行列 M（必要なら）、加わった雑音 N | `loss_noise`、`raman_noise`、`gawbs_noise`、`backscatter_noise`（ブリルアン・ラマン・Rayleigh）、`kerr_fluct`（揺らぎへの Kerr） |
| `fiber_quantum_tw` | Q | 同上 ＋ N_traj、seed | 同上（統計誤差つき） | 既定オフ（検証用） |
| `raman_noise_table` | O | h_R のモデル、γ、f_R、T、P、向き、Δν（probe と squeezed 光の周波数差）、τ_p または u | n_R(f)、N_R | `raman_model` ∈ {blow_wood, hollenbeck_cantrell, lin_agrawal, dc_10lorentz} |
| `gawbs_table` | O | d、V_L、V_S、ファイバの種類（SMF/PM）、T、L、γ_OM と Γ のモード表 | f_m の表、S_φ,G(f)、σ_G²、同じ向きの明るい光からの n_G(f)、N_G | `gawbs_modes` ∈ {R0m, R0m+TR2m} |
| `brillouin_table` | O | g_B、ν_B、Δν_B、A_eff、α、L、P、線幅、D_pr、τ_eff | P_th、G_B、Φ_B、N_B（モードあたり・ゲートあたり）、警告 | `sbs_check`、`brillouin_backscatter` |
| `gawbs_phase_jitter` | T | σ_G²（区間ごと、軸ごと） | 区間の位相への白色（パルスごと）の雑音、または θ_rms に足す定数 | `gawbs_noise` |

- 区間ごとに持つ（D11）。逆向きの probe の後方散乱は「probe が通る区間 → squeezed 光が通る同じ区間」に、時刻をずらさない定常な雑音として入れる（2n_gL/c に広がるので）。
- CW のときは `fiber_quantum_linearized` が側帯波ごとの解析的な行列に縮退する（GNLSE なし、§5.1）。
- 出力の雑音 N は T04 の写像（R′± と θ_rms）に「位相によらない雑音」として渡す。GAWBS のジッタは T04 の θ_rms の高周波側に足す。

### 9.2 他のトピックとの関係

- T06: 古典の GNLSE、γ・A_eff・β_k・損失の値。XPM（逆向き・同じ向きの probe → squeezed 光）は T06。
- T01/T03: 入力共分散 V(0)、LO の時間モード u(T)。
- T04: 加わった雑音 N と位相ジッタ → 測定値。
- T09: PM ファイバの値（A_eff、クラッド径、2 軸）、depolarized GAWBS の 2 軸間の扱い。
- T10: Rayleigh 後方散乱（E44）、PER による漏れ（E43）、ここのブリルアン・ラマンの数値を合わせた「逆向き probe の混入光子」の合計。
- T13/T14: µW の probe でロックの SNR が足りるか。

---

## 10. 規約・DECISIONS との食い違い、ユーザの決定が必要な点

**食い違い・補足（調整役への提案）**

1. **prompt §7-6「Raman noise → 0 as temperature → 0」は文字どおりには成り立たない。** 熱の部分（∝ n_th）は 0 になるが、Stokes 側の自発散乱（∝ g_R(Ω)、T = 0 でも残る。D&C の (n_th + 1/2)、Voss & Kumar 2004 の T = 0 の雑音指数の限界）は残る。
   0 になるのは「熱の部分」と「Ω → 0（probe と squeezed 光が同じ周波数）の極限」。テストは V6 の形（熱の部分 → 0、T = 0 の残りが解析値と一致）を提案する。§7 の文言の修正は調整役の判断。
2. survey E33 の「ピーク周波数 約 20 MHz – 1 GHz」: R0m の最低次は 125 µm で 30.5 MHz。TR2m の最低次は計算していないので、20 MHz の下限は未確認。
3. survey E35 の「約 11 GHz [UNVERIFIED]」: 2n_eff V_A/λ で 10.8–11.2 GHz（V_A による）。値は Kobyakov 2010 で確認が必要（未確認のまま）。
4. survey E34 の「1 km で約 30 mW、63 m で約 0.5 W」: 同じ式で 34 mW、0.53 W（g_B・A_eff は同じく UNVERIFIED）。整合。
5. conventions.md への提案（追加）: 「位相によらない雑音 N 光子 → Var に N（ΔV/V_SN = 2N）」と「モードあたりの光子数 N = ∫ n(f)|ũ(f)|² df、τ_eff = (∫u dT)²」を §2 か §7 に置くと、T07・T10・T11 で数え方が揃う。
6. Ref C の文献 [27] の番号のずれ（§6.5）。

**ユーザの決定が必要な点（推測しない）**

1. **将来の probe のパワー・デューティ比・パルス幅**（D28、Open）: §6.3 のとおり、mW の逆向き CW probe は 630 m・ps パルスで squeezing を大きく下げる（ブリルアンだけで −10 → −7.2 dB、Rayleigh を含めると消える）。
   µW 以下にできるか（検出側の都合）、あるいは probe の波長・方式を変えるか（E46）。
2. **GBS の許容混入光子数**（GBS_REQUIREMENTS、Open）: 10 ns ゲートあたりの値（§6.3）を判定する基準が要る。
3. **将来のパルス幅 τ_p**（D07、Open）: 狭帯域の雑音（ブリルアン・Rayleigh・GAWBS の散乱光）はモードあたり τ_eff に比例するので、10 ps と 100 fs で 100 倍違う。
4. **LO を CW で送るかパルスで送るか、パワー**（D09、Open）: 630 m 以上で CW の LO を ≳10 mW にすると SBS の G_B が 5 を超える。
5. **PM ファイバの具体的な型番（クラッド径 125/80 µm、応力付与部の形）**（D32、Phase 2）: GAWBS のピーク位置と強さが変わる。

---

## 11. 文献

一次文献（書誌は Crossref または PDF で確認。本文を読んだかどうかを併記）:

1. R. M. Shelby, M. D. Levenson, P. W. Bayer, "Guided acoustic-wave Brillouin scattering," Phys. Rev. B 31, 5244–5252 (1985). doi:10.1103/PhysRevB.31.5244. **本文未読**（R0m の周波数の式・TR2m の式・L 比例は原文で確認が必要）。
2. R. M. Shelby, M. D. Levenson, P. W. Bayer, "Resolved forward Brillouin scattering in optical fibers," Phys. Rev. Lett. 54, 939–942 (1985). doi:10.1103/PhysRevLett.54.939. 本文未読。
3. R. M. Shelby, M. D. Levenson, S. H. Perlmutter, R. G. DeVoe, D. F. Walls, "Broad-band parametric deamplification of quantum noise in an optical fiber," Phys. Rev. Lett. 57, 691–694 (1986). doi:10.1103/PhysRevLett.57.691. 本文未読。
4. S. H. Perlmutter, M. D. Levenson, R. M. Shelby, M. B. Weissman, "Polarization properties of quasielastic light scattering in fused-silica optical fiber," Phys. Rev. B 42, 5294–5305 (1990). doi:10.1103/PhysRevB.42.5294. 本文未読（低周波の準弾性散乱、ラマンの Ω → 0 に関係）。
5. A. J. Poustie, "Bandwidth and mode intensities of guided acoustic-wave Brillouin scattering in optical fibers," J. Opt. Soc. Am. B 10, 691 (1993). doi:10.1364/JOSAB.10.000691. 本文未読。
6. A. J. Poustie, "Guided acoustic-wave Brillouin scattering with optical pulses," Opt. Lett. 17, 574 (1992). doi:10.1364/OL.17.000574. 本文未読。
7. P. D. Townsend, A. J. Poustie, "Measured reduction of guided-acoustic-wave Brillouin scattering in a fiber interferometer by time-delayed pulses," Opt. Lett. 20, 37 (1995). doi:10.1364/OL.20.000037. 本文未読。
8. D. Elser, U. L. Andersen, A. Korn, O. Glöckl, S. Lorenz, Ch. Marquardt, G. Leuchs, "Reduction of guided acoustic wave Brillouin scattering in photonic crystal fibers," Phys. Rev. Lett. 97, 133901 (2006). doi:10.1103/PhysRevLett.97.133901（arXiv:quant-ph/0512044 を読んだ）。
9. Y. London, H. H. Diamandi, G. Bashan, A. Zadok, "Invited Article: Distributed analysis of nonlinear wave mixing in fiber due to forward Brillouin scattering and Kerr effects," APL Photonics 3 (2018). doi:10.1063/1.5042635（arXiv:1811.06282 を読んだ。R0,7 = 321.75 MHz、4.3 MHz、γ_OM ≈ 2 /(W km)、γ = 1.3 /(W km)）。論文番号は未確認。
10. N. Takefushi, M. Yoshida, K. Kasai, T. Hirooka, M. Nakazawa, "Theoretical and experimental analyses of GAWBS phase noise in various optical fibers for digital coherent transmission," Opt. Express 28, 2873 (2020). doi:10.1364/OE.384505（要旨のみ読んだ）。
11. Nishizawa, Miyauchi, Kume, Mori, Goto（名の頭文字は未確認）, "Symmetric and asymmetric fiber loop mirrors for observing guided-acoustic-wave Brillouin scattering in polarization-maintaining fibers," Opt. Lett. 19, 1424 (1994). doi:10.1364/OL.19.001424. 本文未読。
12. Nishizawa, Kume, Mori, Goto, Miyauchi（名の頭文字は未確認）, "Characteristics of guided acoustic wave Brillouin scattering in polarization maintaining fibers," Opt. Rev. 3, 29–33 (1996). doi:10.1007/s10043-996-0029-2. 本文未読。
13. P. D. Drummond, J. F. Corney, "Quantum noise in optical fibers. I. Stochastic equations," J. Opt. Soc. Am. B 18, 139 (2001). doi:10.1364/JOSAB.18.000139（arXiv:quant-ph/9912095 の式 (64)、(74)–(77) を確認）。
14. J. F. Corney, P. D. Drummond, "Quantum noise in optical fibers. II. Raman jitter in soliton communications," J. Opt. Soc. Am. B 18, 153 (2001). doi:10.1364/JOSAB.18.000153. 本文未読。
15. S. J. Carter, P. D. Drummond, "Squeezed quantum solitons and Raman noise," Phys. Rev. Lett. 67, 3757–3760 (1991). doi:10.1103/PhysRevLett.67.3757. 本文未読。
16. H. A. Haus, Y. Lai, "Quantum theory of soliton squeezing: a linearized approach," J. Opt. Soc. Am. B 7, 386 (1990). doi:10.1364/JOSAB.7.000386. 本文未読。
17. F. X. Kärtner, D. J. Dougherty, H. A. Haus, E. P. Ippen, "Raman noise and soliton squeezing," J. Opt. Soc. Am. B 11, 1267 (1994). doi:10.1364/JOSAB.11.001267. 本文未読。
18. M. Kitagawa, Y. Yamamoto, "Number-phase minimum-uncertainty state with reduced number uncertainty in a Kerr nonlinear interferometer," Phys. Rev. A 34, 3974–3988 (1986). doi:10.1103/PhysRevA.34.3974. 本文未読（V4 の式はこの担当が線形化から導いたもの）。
19. A. Sinatra, C. Lobo, Y. Castin, "The truncated Wigner method for Bose-condensed gases: limits of validity and applications," J. Phys. B 35, 3599–3631 (2002). doi:10.1088/0953-4075/35/17/301. 本文未読（TW の適用範囲の一般論）。
20. P. L. Voss, P. Kumar, "Raman-noise-induced noise-figure limit for χ(3) parametric amplifiers," Opt. Lett. 29, 445 (2004). doi:10.1364/OL.29.000445（要旨相当のみ）。
21. P. L. Voss, K. G. Köprülü, P. Kumar, "Raman-noise-induced quantum limits for χ(3) nondegenerate phase-sensitive amplification and quadrature squeezing," J. Opt. Soc. Am. B 23, 598 (2006). doi:10.1364/JOSAB.23.000598. 本文未読。
22. K. J. Blow, D. Wood, "Theoretical description of transient stimulated Raman scattering in optical fibers," IEEE J. Quantum Electron. 25, 2665–2673 (1989). doi:10.1109/3.40655. 本文未読（τ1、τ2、f_R は UNVERIFIED）。
23. R. H. Stolen, J. P. Gordon, W. J. Tomlinson, H. A. Haus, "Raman response function of silica-core fibers," J. Opt. Soc. Am. B 6, 1159 (1989). doi:10.1364/JOSAB.6.001159. 本文未読。
24. D. Hollenbeck, C. D. Cantrell, "Multiple-vibrational-mode model for fiber-optic Raman gain spectrum and response function," J. Opt. Soc. Am. B 19, 2886 (2002). doi:10.1364/JOSAB.19.002886. 本文未読。
25. Q. Lin, G. P. Agrawal, "Raman response function for silica fibers," Opt. Lett. 31, 3086 (2006). doi:10.1364/OL.31.003086. 本文未読。
26. R. G. Smith, "Optical power handling capacity of low loss optical fibers as determined by stimulated Raman and Brillouin scattering," Appl. Opt. 11, 2489 (1972). doi:10.1364/AO.11.002489. 本文未読（係数 21 は広く引用される値、UNVERIFIED）。
27. R. W. Boyd, K. Rzążewski, P. Narum, "Noise initiation of stimulated Brillouin scattering," Phys. Rev. A 42, 5514–5521 (1990). doi:10.1103/PhysRevA.42.5514. 本文未読。
28. A. Kobyakov, M. Sauer, D. Chowdhury, "Stimulated Brillouin scattering in optical fibers," Adv. Opt. Photon. 2, 1 (2010). doi:10.1364/AOP.2.000001. 本文未読（g_B、ν_B、Δν_B の確認先）。
29. P. Eraerds, N. Walenta, M. Legré, N. Gisin, H. Zbinden, "Quantum key distribution and 1 Gbps data encryption over a single fibre," New J. Phys. 12, 063027 (2010). doi:10.1088/1367-2630/12/6/063027. 本文未読（前方・後方ラマン雑音の測定の確認先）。
30. K. Bergman, H. A. Haus, "Squeezing in fibers with optical pulses," Opt. Lett. 16, 663 (1991). doi:10.1364/OL.16.000663. 本文未読。

リファレンス（`refs/`）:

- Ref C: S. Verclas, C. Darsow-Fromm, B. Tohermes, R. Schnabel, "Fiber distribution of phase-stabilized GHz-bandwidth squeezed vacuum states of light between two buildings," Optica Quantum 4, 1 (2026). doi:10.1364/OPTICAQ.567418（§1 の GAWBS の記述、§2 を読んだ）。
- H.-L. Liu et al., "Gaussian boson sampling with 1,024 squeezed states in 8,176 modes," Nature (2026). doi:10.1038/s41586-026-10523-6（Methods「Phase stabilization of the system」を読んだ: 逆向き・同じ周波数・直交偏光の probe、AOM でパルス化して squeezed パルス列から時間をずらす、TEC のフィードバック。probe のパワーは本文に書かれていない）。

教科書（UNVERIFIED、このセッションでは読んでいない）: G. P. Agrawal, *Nonlinear Fiber Optics*（GNLSE の形、g_R の式、SBS の線幅の式）。
