# T01 導波路 OPA による squeezing（opa-waveguide-squeezing）

- Phase: 1（deep dive）、トピック T01（`notes/survey.md` の E01–E07）
- 日付: 2026-10-28
- 入力: `docs/squeezing_sim_prompt.md`、`DECISIONS.md`（D01–D33）、`notes/conventions.md`、`notes/survey.md`、Ref A 本文と図（`refs/`）、Ref A の文献 [19]、[23]、[28] の本文（arXiv 版）、Liu et al. 2026（`refs/`）
- 依存: なし。T03（パルス）、T04（位相雑音 → squeezing）、T12（損失・モード整合）、T14（誤差信号）と数値をやりとりする。依存先の値はリファレンス値で代用し、そう書く。
- 値の印は survey と同じ: **[Ref A §x]**、**[AE]**（この担当による概算・導出、未検証）、**[UNVERIFIED]**（原典で確認していない）。数値の検算スクリプトはスクラッチ（コミットしない）で、結果だけをここに書く。

## 0. 結論の要約

1. prompt §7-1 の単一モードモデルは L = 8 %、α = 9.06 /W、P = 0.64 W、θ = 9 mrad で **R′− = −10.148 dB、R′+ = +20.556 dB**。conventions §6 の厳密なガウス平均でも 10⁻⁴ dB 以内で同じ。
2. **640 mW は squeezer OPA の「出口」のポンプパワー（推定値）** [Ref A §3.2「The output pump power immediately after the squeezer OPA is estimated to be 640 mW」、§4.1「pump power after the squeezer OPA」]。α = 906 %/W はこの P の定義と集中損失モデルに結びついた**当てはめの値**で、導波路の SHG 効率の測定値（[19] で 820 %/W）とは別物として扱う。
3. 分布損失（[28] の付録）を入れると、0.1 dB/cm・45 mm・640 mW で導波路の実効損失は **2.03 %**（Ref A の「2 %」を再現）。ただし同じ α で比べると反 squeezing が集中損失モデルより **0.27 dB 低く**なる（20.29 dB。測定値 20.2 ± 0.2 dB に近い）。完全なシミュレーションの許容幅はこの差を含めて決める（§7.3）。
4. 位相の符号は conventions §4 と整合する（θ_amp = φ_p/2、θ_sq = φ_p/2 + π/2 を数値で確認）。結合波方程式の形（dA/dz = iκ A_p A* e^{iΔkz}）からは π/2 の定数オフセットが出るので、OPA ごとの定数 φ_OPA,0 として吸収する（§2.2）。
5. Ref A の lock 2（BHD の 1 MHz ビート）は、squeezer を通った probe を使うので、**lock 1 の誤差 ε は LO の誤差に e^{−2r} ε（Ref A で 0.008 ε）しか伝わらない** [AE]。lock 1 の役割は、主に probe を増幅側の直交位相に置いて lock 2 の信号を大きく保ち、復調位相を固定することにある。逆に D30 の統合案のように squeezer に probe を通さない構成では、**squeezer の結晶そのものの温度による位相（~6 rad/K、[UNVERIFIED] の dn/dT による）と位相不整合による回転がどのループにも入らない**（§8）。
6. 将来のパルス構成（D07）では、45 mm の PPLN の群速度不整合による走り抜け（walk-off）は **約 13 ps** [AE、[19] の SHG の波長許容幅から]。掃引範囲 100 fs–10 ps はすべて準 CW 近似の外にあり、45 mm の導波路では T03 の多モード（時間モード）の扱いが要る。

## 1. 記号（symbols, English）

- `a_s` [1/m]: power propagation-loss coefficient of the signal (1545 nm) mode in the waveguide. `a_s ℓ = (dB/cm)·(ℓ/cm)·ln10/10`.
- `a_p` [1/m]: power propagation-loss coefficient of the pump (772.7 nm) mode.
- `ℓ` [m]: waveguide (crystal) length. Ref A: 0.045 m.
- `η_SHG` [1/W]: normalized SHG efficiency of the waveguide, `P_SH = η_SHG P_FH²` (undepleted, lossless, phase matched). [19]: 820 %/W = 8.20 /W measured.
- `α` [1/W]: fitted "SHG efficiency" in the Ref A model (906 %/W = 9.06 /W), used with `P` = pump power at the waveguide output.
- `P` [W]: pump power used in the Ref A model (= estimated output pump power). `P_in`, `P_out`: pump power at the input / output face.
- `g` [1/m]: amplitude parametric gain per length, `g = r/ℓ` for uniform pump.
- `r`: squeeze parameter (conventions §3). Lossless, phase matched: `r = √(η_SHG P)`.
- `R± = V±/V_vac`: anti-squeezed (+) / squeezed (−) quadrature variance normalized to shot noise (= vacuum, `V_vac = 1/2`).
- `R′±`: same, including phase jitter.
- `L`: total effective loss (lumped). `η = 1 − L`: total efficiency.
- `L_wg`: effective loss of the waveguide alone (defined by `R− = L_wg + (1 − L_wg) e^{−2r}`).
- `θ`, `σ`: RMS phase error between LO and the squeezed quadrature, `Δ = φ_LO − θ_sq` (conventions §5, §6). Ref A writes it `θ̃`.
- `φ_p`: pump phase (conventions §4). `θ_sq = φ_p/2 + π/2`, `θ_amp = φ_p/2`.
- `φ_OPA,0`: constant offset mapping the physical pump phase in a given crystal to conventions' `φ_p` (§2.2).
- `Δk` [1/m]: phase mismatch `k_p − 2k_s − 2π/Λ` (Λ: poling period). `δ = Δk/2`.
- `ΔT_FWHM` [K]: FWHM of the SHG temperature tuning curve (sinc²) of the crystal. `Δλ_FWHM` [m]: FWHM of the SHG wavelength tuning curve (fundamental wavelength).
- `n_g,s`, `n_g,p`: group indices at signal / pump. `Δn_g = n_g,p − n_g,s`. `GVM = Δn_g/c` [s/m].
- `μ`, `ν`: Bogoliubov coefficients of the OPA, `â_out = μ â_in + ν â_in†` (convention: phases of the pump absorbed as in conventions §4).
- `Ω_m` [rad/s]: AOM net frequency shift of the probe (Ref A: 2π × 1 MHz).
- `r_pd`: squeeze (gain) parameter of the phase-detection OPA (PD-OPA).
- `ε`: lock-1 error (deviation of `φ_pr − Ω_m t − φ_p/2` from its set point).
- `v`: `φ_LO − φ_p/2` (LO relative to half pump phase); squeezing measured at `v = π/2`.
- `τ_p`, `f_rep`, `P_peak`, `P_avg`: pulse duration, repetition rate, peak / average pump power (D07).
- `w_j`, `r_j`, `Δθ_j`: weight of LO overlap, squeeze parameter, and squeezing-angle offset of spatial mode j (§2.9).

## 2. 支配方程式

### 2.1 単一モード・集中損失のモデル（Ref A §4.1、prompt §7-1）

Ref A（[23]、[19] も同じ形）:

- R± = L + (1 − L)·exp(±2√(αP))
- R′±(θ) ≈ R± cos²θ + R∓ sin²θ

規約との対応: R± = V±/V_vac（conventions §2 の S_dB = 10 log10 R）。損失は conventions §2 の純損失 V_out = ηV_in + (1−η)V_vac を R で書いたもの。
位相ジッタは conventions §6 の厳密な式 R′ = (R− + R+)/2 ± (R− − R+)/2 · e^{−2σ²} があり、Ref A の式はその σ ≪ 1 の近似。θ = 9 mrad では差は 10⁻⁴ dB 未満（§7.1）。
[23] は η = 1 − L を「実効透過率」と書いて同じ式を使う（[23] Eq. (2)）。

最適ポンプ（[AE]、近似式から導出）: R′− = L + (1 − L)(e^{−2r} cos²θ + e^{2r} sin²θ) を r で最小化すると
- e^{4r_opt} = cot²θ、すなわち **r_opt = ½ ln(cot θ)**、P_opt = r_opt²/α
- **R′−,min = L + (1 − L) sin 2θ**
θ = 9 mrad、α = 9.06 /W で r_opt = 2.355、P_opt = 612 mW、R′−,min = −10.15 dB。Ref A は最適点（640 mW）の近くで測っている。
この式は T04 の「目標 × ファイバ長 → 必要な θ_rms」表（D16）の解析的な基準になる（§6.3）。

### 2.2 r = √(ηP) の出どころと位相の符号

パワーで規格化した包絡線 A（|A|² = パワー [W]）の縮退な結合波方程式（損失なし。ポンプ枯渇を含む一般形）:

- dA_s/dz = iκ A_p A_s* e^{iΔkz}
- dA_p/dz = iκ A_s² e^{−iΔkz}

この係数の組でエネルギーが保存する（d(|A_s|² + |A_p|²)/dz = 0 を確認）[AE]。
- SHG（A_p(0) = 0、ポンプが減らない、Δk = 0）: A_p(ℓ) = iκℓ A_s², P_SH = κ²ℓ² P_FH² → **η_SHG = κ²ℓ²**。
- OPA（A_p 一定、Δk = 0）: 直交位相の振幅が e^{±κ|A_p|ℓ} → **r = κ|A_p|ℓ = √(η_SHG P)**、分散は e^{±2r}。
  これが Ref A の exp(±2√(αP))。

符号の変換（文献の規約 → conventions §4）:
- conventions の H_int = i(κ/2)(β_p â†² − β_p* â²) から Heisenberg 方程式は dâ/dt = κβ_p â†（i が付かない）。増幅されるのは 2φ_s − φ_p = 0、θ_amp = φ_p/2、θ_sq = φ_p/2 + π/2。
  数値で確認した（φ_p = 0.7 rad の 2×2 シンプレクティック行列で分散最小の角度 = 1.92080 rad = φ_p/2 + π/2、V_min = ½e^{−2κt} と一致）。
- 結合波方程式は iκA_p を持つので、増幅されるのは 2φ_s − φ_p = π/2。つまり **φ_p(conventions) = φ_p(結合波) + π/2**。
- さらに分極反転の符号や SHG 結晶の違いで定数がずれる。そこで **OPA（結晶）ごとに定数 φ_OPA,0 を持たせ、φ_p = （その結晶に入るポンプの複素振幅の偏角）+ φ_OPA,0 とする**ことを提案する。ロックの設定値で吸収される定数なので、物理量には効かないが、符号のテストで必要になる（§9 V-T01-5）。

### 2.3 分布損失（E01、[28] の付録）

導波路の中で生成と減衰が同時に起きる。[28]（Yamashima et al. 2025、Appendix）の結果を、この規約に直したもの:

- dR±/dz = ±2g R± − a_s (R± − 1)
- 入力が真空（R±(0) = 1）、g が一定のとき
  **R± = [±2g e^{(±2g − a_s)ℓ} − a_s] / (±2g − a_s)**

規約の変換: [28] は振幅の消光係数 α_Y（パワーで 2α_Y）と振幅の利得 g を使い、x̂_out = e^{−(g+α_Y)L} x̂_in + √(α_Y(1 − e^{−2(g+α_Y)L})/(g+α_Y)) x̂_vac と書く。a_s = 2α_Y と置けば上の式と一致する（確認済み）。
[28] の x̂ は squeezed 側、p̂ は増幅側。
ODE を数値で解いた結果と閉じた式は 10⁻⁸ の相対差で一致した。

実効損失: L_wg = (R− − e^{−2r})/(1 − e^{−2r})。R+ のほうは同じ L_wg の集中損失より小さい（利得そのものが e^{−a_s ℓ/2} 程度下がるため）。これは「分布損失は集中損失 L_wg と等価ではない」ことを意味する。

| a_s（dB/cm） | P [W] | L_wg | R− [dB] | R+ [dB] | 同じ L_wg の集中損失での R+ |
|---|---|---|---|---|---|
| 0.1 | 0.1 | 3.53 % | −7.47 | 8.02 | 8.14 |
| 0.1 | 0.3 | 2.68 % | −12.02 | 14.00 | 14.21 |
| 0.1 | **0.64** | **2.03 %** | −15.50 | 20.56 | 20.83 |
| 0.1 | 1.0 | 1.67 % | −17.20 | 25.77 | 26.07 |
| 0.1 | 1.2 | 1.53 % | −17.78 | 28.26 | 28.57 |
| 0.07 | 0.64 | 1.43 % | −16.53 | 20.67 | 20.85 |

（r = √(9.06 P)、ℓ = 45 mm。0.1 dB/cm は Ref A の値で、[19] の測定は「45 mm で 7 %、すなわち多くとも 0.1 dB/cm」[19]。7 %/45 mm は 0.07 dB/cm に当たる。）
**L_wg はポンプパワーに依存する**（0.1 W で 3.5 %、1.2 W で 1.5 %）。Ref A の Fig. 5 は L = 8 % を一定として当てはめているので、ここにモデルの差が出る（§7.2）。

### 2.4 ポンプの損失と「P」の定義

ポンプも導波路内で減衰する（a_p、772.7 nm での値は Ref A・[19] に記載なし）。P(z) = P_in e^{−a_p z} とすると
- r = κ∫₀^ℓ |A_p(z)| dz = (√η_SHG/ℓ)∫₀^ℓ √P(z) dz = √(η_SHG P_in) · (2/(a_p ℓ))(1 − e^{−a_p ℓ/2})
  = √(η_SHG P_out) · (2/(a_p ℓ))(e^{a_p ℓ/2} − 1) ≈ √(η_SHG P_out)(1 + a_p ℓ/4)。

Ref A は P = P_out（出口のポンプの推定値）で当てはめたので、**α ≈ η_SHG(1 + a_p ℓ/2)** になる [AE]。
906/820 = 1.105 から逆算すると a_p ℓ ≈ 0.2（772.7 nm で約 0.2 dB/cm）になるが、Ref A の導波路は [19] と「同じ方法で作った」別の素子である可能性があり、パワーの校正の差でも説明できるので、**これは整合の確認にすぎない**（[AE]、根拠にしない）。
シミュレーションでは α と P の組を「当てはめの値の組」として扱い、P の基準（input/output）をパラメータで明示する。

### 2.5 位相不整合（E05）

ポンプが一定、Δk ≠ 0 のとき、B = A e^{−iΔkz/2}（出口のポンプの位相を基準にした枠）で
- d/dz (q, p)ᵀ = M (q, p)ᵀ、M = [[g, δ], [−δ, −g]]、δ = Δk/2（q が増幅側になる位相の選び方）
- 固有値 ±√(g² − δ²)。δ > g では増幅しない。
- 転送行列 T = exp(Mℓ) = cosh(sℓ) I + sinh(sℓ) M/s、s = √(g² − δ²)。det T = 1。
- 実効的な squeeze parameter: cosh(2r_eff) = cosh²(sℓ) + sinh²(sℓ)(g² + δ²)/(g² − δ²)。
- squeezed 軸の回転（出口のポンプ/2 に対して）: **δθ_mm ≈ −δℓ/(2r) = −Δkℓ/(4r)**（r ≫ 1、|δ| ≪ g）。

数値（r = 2.408、Ref A の値）[AE]:

| Δkℓ/2 | r_eff | δθ_mm |
|---|---|---|
| 0.01 | 2.4080 | −2.04 mrad |
| 0.1 | 2.4068 | −20.4 mrad |
| 0.3 | 2.3969 | −61.4 mrad |
| 1.392（SHG の半値） | 2.1617 | −299 mrad |

符号は M の取り方で決まるので、実装では V-T01-7 のテストで固定する。

温度との対応: SHG の sinc²(Δkℓ/2) の半値は Δkℓ/2 = ±1.392 なので、**Δkℓ/2 = 2.784 (T − T_pm)/ΔT_FWHM**。
45 mm 素子の ΔT_FWHM は Ref A・[19] に記載がない。[19] Fig. 1(d) の SHG の波長許容幅は図から読んで Δλ_FWHM ≈ 0.27 nm（[AE]、±30 %）。dλ_pm/dT ≈ 0.1 nm/K（[UNVERIFIED]、PPLN の 1550 nm SHG の典型値）とすれば ΔT_FWHM ≈ 2.7 K [AE]。
- r = 2.4 で、温度ずれ 10 mK → Δkℓ/2 ≈ 0.010 → δθ_mm ≈ 2 mrad、r_eff の低下は無視できる。

**結晶そのものの位相**（出口のポンプの位相）: 出口の squeezing 角は
**θ_sq,out = [φ_p,in + ψ_p,xtal]/2 + π/2 + δθ_mm**、ψ_p,xtal = k_p(T)ℓ（分極反転格子の定数を含む）。
すなわち、**squeezer の結晶は「長さ ℓ、屈折率 n_e(2ω0) のポンプの光路の区間」として conventions §9 の区間モデルに入れ、θ_sq にはその 1/2 が入る**。
温度係数は dψ_p,xtal/dT ≈ (2ω0/c)(dn_e/dT + n_e α_th)ℓ。dn_e/dT ≈ (3–4)×10⁻⁵ /K（[UNVERIFIED]、ZnO 添加の値は不明）、熱膨張を無視すると θ_sq には ½ × 2 × (2π/1.545 µm) × 3.5×10⁻⁵ × 0.045 ≈ **6 rad/K（6 mrad/mK）** [AE]。
Ref A では probe も同じ結晶を通るので lock 2 で打ち消される（§8.1）。

### 2.6 ポンプ枯渇（E02）

結合波方程式（§2.2）をそのまま解けば入る。squeezed vacuum では生成されるパワーは
P_sq ≈ ħω0 · sinh²r · B（B: 位相整合帯域）≈ 1.3×10⁻¹⁹ J × 30 × 6 THz ≈ **2×10⁻⁵ W** [AE]（B = 6 THz は [19] の値）。
640 mW のポンプに対して 4×10⁻⁵ なので、squeezer では無視できる。
ただし PD-OPA（Ref A）は mW の probe を入れ、利得があるので、PD-OPA のポンプパワー（Ref A に記載なし）によっては枯渇が効く。これは誤差信号の傾きの非線形性に効くだけで、squeezing には効かない。

### 2.7 ポンプの強度雑音と位相雑音（E03、E04）

- 強度: δr/r = δP/(2P)。反 squeezing の傾き dR′+_dB/dP ≈ 8.686 r/(2P) = **16.3 dB/W**（Ref A の点）。
  P が +1 % 変わると R′+ は +0.104 dB、R′− は +0.002 dB（最適点の近くなので平ら）。+5 % で +0.52 dB / +0.015 dB。
  MHz の RIN は r(t) を変調するが、真空の入力なので平均の分散への効果は 2 次（⟨e^{−2r(t)}⟩ ≈ e^{−2r}(1 + 2⟨δr²⟩)）で、RIN の積分値 10⁻³ なら δr ≈ 1×10⁻³ と小さい [AE]。
  遅い変動は R′+ の時系列の揺れとロックの誤差信号の傾き（Ref A の lock 2 の傾き ∝ e^r）に効く。
- 位相: θ_sq = φ_p/2 + π/2 なので、**ポンプの光路の位相 ψ_p の揺れは θ_sq に ½ψ_p で入る**。同じ長さの変化 δL なら ψ_p = 2ω0 n_p δL/c なので、θ_sq の変化は基本波の光路 ω0 n_p δL/c と同じ大きさになる。
  Ref A ではポンプは光源の SHG（φ_p = 2φ_L + 定数）で、LO・probe と同じレーザーから作るので、レーザー位相雑音は遅延差の形で入る（T05）。経路ごとの位相（conventions §9）で扱う。

### 2.8 ポンプ誘起の吸収・光屈折（E06）

[19] の導波路は ZnO 添加の PPLN（直接接合、機械研削）で、「数百 mW のポンプでも当てはめの SHG 係数（823 %/W）が測定値（820 %/W）と一致し、ポンプありの損失がポンプなしの挿入損失と一致した」ので、ポンプ誘起の劣化は小さいと述べている [19]。
Ref A の Fig. 5 は 40–985 mW で L 一定の当てはめに図から読んだ値で ±0.3 dB 以内で合う（§7.2）。L が 1 W で 1 % 増えると R′− は +0.04 dB しか動かないので、**この当てはめから言える上限は「1 W までで ΔL ≲ 数 %」程度で、きつい制限ではない** [AE]。
モデルとしては L_wg を P の関数（呼び出し可能なパラメータ、Ref B の L(x) と同じ扱い）にし、既定は一定とする。
ポンプの吸収による発熱で位相整合温度がずれる効果（熱的なずれ）は、Δk(P, T) として §2.5 に入れる（値はなし、既定オフ）。

### 2.9 空間多モード性（E07）

[19] の導波路は「準単一モード」で、1545 nm で 2 つのモードしか持たず、2 番目のモードは非対称なので Gauss 型のポンプとほとんど結合しないと述べている [19]。
Ref A の残りのモード不整合 4 % のうち、probe–LO の可視度 99 %（効率 0.99² = 0.98、損失 2 %）は測定値で、残り 2 % は原因不明（「squeezed 光と probe の空間モードの僅かな違い」と推定）[Ref A §4.2]。

モデル [AE]: 出口の squeezed 光を空間モード u_j（squeeze parameter r_j、squeezing 角のずれ Δθ_j。導波路モードの伝搬定数が違えば Δθ_j ≠ 0）の和とし、
R = Σ_j w_j R_j(Δ + Δθ_j) + (1 − Σ_j w_j)、w_j = |⟨LO|u_j⟩|²。
既定は単一モード（w_0 = η_mm）。Ref A の値: η_mm = 0.96、うち probe–LO の重なりで 0.98、squeezed–probe の差で 0.98。値の担当は T12。

## 3. 仮定と適用範囲

- 単一の空間モード、単一の時間モード（CW、または準 CW: τ_p ≫ GVM·ℓ）。
- ポンプは古典的で、雑音は遅い（パラメータの揺れとして扱う）。squeezer ではポンプ枯渇なし（§2.6）。
- 位相整合帯域（THz）は検出帯域（35 MHz）より十分に広いので、ESA の周波数では R± は周波数に依らない（Ref A Fig. 4(b) で 100 MHz まで 8 dB 以上、低下は PD の帯域による [Ref A §3.2]）。周波数特性は T13 の検出器の応答で入れる。
- 集中損失モデル（§2.1）は 1 点の P では分布損失と等価に合わせられるが、P を変えると L_wg(P) が変わる（§2.3）。
- 小角近似 R′ ≈ R cos²θ + R sin²θ は θ ≲ 0.1 rad で厳密なガウス平均と 0.01 dB 以内（θ = 9 mrad で 10⁻⁴ dB）。大きな θ（D16 の 3 dB 目標で 0.2 rad 級）では conventions §6 の厳密な式を使う。
- 位相不整合の回転 δθ_mm ≈ −Δkℓ/(4r) は |δ| ≪ g、r ≫ 1 で有効。一般には §2.5 の行列を数値で使う。

## 4. 効果が働く周波数帯

| 効果 | 帯域 | 備考 |
|---|---|---|
| E01 OPA squeezing（分布損失） | DC–THz（位相整合帯域 2.5–6 THz [18][19]） | ESA の 3 MHz を含む全帯域で一定 |
| E02 ポンプ枯渇 | — | 静的 |
| E03 ポンプ強度雑音 | DC–検出帯域 | 遅い成分は R′+ の揺れ、速い成分は 2 次の効果 |
| E04 ポンプ位相雑音 | DC–MHz | ½ψ_p として θ_sq へ。ループ帯域内はロックで補正 |
| E05 位相整合の温度 | < 1 Hz（TEC の時定数） | 回転 δθ_mm と、結晶の光路の位相 ½ψ_p,xtal |
| E06 ポンプ誘起の吸収 | s–h | L(P, t) |
| E07 空間多モード | 静的 | η_mm |

## 5. スケーリング

- **ポンプパワー**: r ∝ √P、R+ ∝ e^{2√(αP)}。R′− は P_opt = (½ ln cot θ)²/α で最小。P_opt は θ だけで決まり L に依らない。L_wg は P とともに下がる（§2.3 の表）。
- **導波路長**: η_SHG ∝ ℓ²（[19] で 40 %/W·cm² × 4.5² cm² ≈ 810 %/W）、r ∝ ℓ√P。伝搬損失 a_sℓ ∝ ℓ。走り抜け GVM·ℓ ∝ ℓ。位相整合の温度許容幅 ΔT_FWHM ∝ 1/ℓ。結晶の熱的な位相 ∝ ℓ。
- **温度**: Δkℓ/2 = 2.784 δT/ΔT_FWHM。回転 ≈ Δkℓ/(4r)。結晶の位相 ≈ 6 rad/K（45 mm、[UNVERIFIED] の dn/dT）。
- **ファイバ長**: T01 の効果はファイバの手前で起きるので、ファイバ長には依らない。ただし R+（反 squeezing の大きさ）がファイバによる位相雑音の罰則を決め、P_opt はファイバを含めた θ_tot で決まる（§6.3）。
- **パルス幅**: 準 CW なら r(t) = √(αP(t)) がポンプの包絡線に従う。τ_p ≲ GVM·ℓ（45 mm で約 13 ps）では準 CW が破れ、時間モードの多モード性（Schmidt 数）を T03 で扱う。ピークパワーで r が決まるので、P_avg = P_peak f_rep τ_p（矩形の場合）。

## 6. 概算

### 6.1 Ref A（リファレンスの値）

- r = √(9.06 × 0.64) = 2.408。R− = −10.58 dB（θ = 0）、R′− = −10.15 dB、R′+ = 20.56 dB。
- L の内訳 [Ref A §4.2]: 導波路 2 %（§2.3 で 2.03 % を再現）、光学部品 1 %、光電変換 1 %、回路雑音 < 0.2 %、モード不整合 4 %。掛け算で 1 − 0.98 × 0.99 × 0.99 × 0.998 × 0.96 = 7.98 %。
- 感度: L +0.1 % → R′− +0.044 dB。θ +1 mrad → R′− +0.096 dB。P または α +1 % → R′+ +0.104 dB。
- 生成される squeezed 光のパワー ≈ 2×10⁻⁵ W（§2.6）。

### 6.2 Ref D（D02）

Ref D には OPA がない（coherent 光の MZI）。**T01 のブロックは Ref D の構成ではすべてオフ**。Ref D の検証には関係しない。

### 6.3 将来の構成 F（D07–D13、D30–D32）

光源の仕様は Open（D07）なので、OPA は Ref A の値（PPLN 45 mm、α = 9.06 /W、0.1 dB/cm、1545.32 nm）で代用する。

- **目標 squeezing と許される θ**（§2.1 の R′−,min = L + (1 − L) sin 2θ、α = 9.06 /W）[AE]:

| 目標（GBS 入力） | L_tot | θ_max | P_opt | そのときの R′+ |
|---|---|---|---|---|
| 3 dB | 8 % | 238 mrad | 55 mW | 5.9 dB |
| 3 dB | 20 % | 193 mrad | 74 mW | 6.3 dB |
| 3 dB | 40 % | 85 mrad | 168 mW | 8.7 dB |
| 6 dB | 8 % | 94 mrad | 154 mW | 9.9 dB |
| 6 dB | 20 % | 32 mrad | 327 mW | 14.0 dB |
| 6 dB | ≥ 25 % | 不可能 | — | — |
| 10 dB | 8 % | 10.9 mrad | 564 mW | 19.3 dB |
| 10 dB | > 10 % | 不可能 | — | — |

  L_tot は OPA から GBS 入力までの全損失（ファイバ、振り分け素子、コネクタを含む。T12・T11）。θ が 0.1 rad を超える行は小角近似なので、T04 が厳密な式で作り直すこと。
  θ_max の中でファイバの位相雑音に使える分は、ここから OPA 側の残差を引いたもの。**32 ポートすべてで同じ θ_max を満たす必要がある**（ポートごとに L_k が違えば L_tot も違う）。
- **パルス**: [19] の SHG の波長許容幅 Δλ_FWHM ≈ 0.27 nm（図から読んだ値）から Δn_g = 5.567 λ²/(4π ℓ Δλ_FWHM) ≈ 0.087、GVM ≈ 0.29 ps/mm、**ℓ = 45 mm で走り抜け 13 ps** [AE]。D07 の掃引範囲 100 fs–10 ps はすべてこれより短いか同程度で、45 mm の導波路では準 CW のモデルは使えない（T03）。
  例: 準 CW を保つために ℓ = 1 cm（走り抜け 2.9 ps）にすると η_SHG ≈ 40 %/W（[19] の 40 %/W·cm²、[AE]）、r = 2.4 に P_peak ≈ 14.5 W、f_rep = 100 MHz・τ_p = 10 ps（デューティ 10⁻³）で P_avg ≈ 15 mW [AE]。
  45 mm のまま準 CW の式を使えば、r = 2.4 に要る P_peak は 0.64 W で、100 MHz・1 ps なら P_avg = 64 µW（矩形）だが、これは準 CW の範囲外の数字で根拠にしない。
- ポンプ枯渇は F でも無視できる（1 パルスあたりポンプ光子 10⁶ 以上に対し、squeezed 光の光子は sinh²r × モード数 ~ 10³ 以下）[AE]。

## 7. §7-1・§7-2 の検算と許容幅

### 7.1 §7-1（単一モードモデル）

| 量 | 小角近似 | 厳密なガウス平均（conventions §6） | Ref A の測定値 |
|---|---|---|---|
| R′− | −10.148 dB | −10.148 dB | −10.1 ± 0.2 dB |
| R′+ | +20.556 dB | +20.556 dB | +20.2 ± 0.2 dB |

prompt の「約 −10.1 / +20.5 dB」と一致する。

### 7.2 §7-2（Fig. 5 の 4 本の曲線）

R′−（dB）/ R′+（dB）。α = 9.06 /W、集中損失モデル。

| P [mW] | 40 | 100 | 200 | 400 | 640 | 800 | 1000 | 1200 |
|---|---|---|---|---|---|---|---|---|
| this work（θ = 9 mrad、L = 8 %） | −4.48 / 4.98 | −6.62 / 7.96 | −8.43 / 11.36 | −9.84 / 16.18 | −10.15 / 20.56 | −9.98 / 23.02 | −9.47 / 25.78 | −8.67 / 28.28 |
| 前の研究（14 mrad、12 %） | −4.15 / 4.85 | −5.98 / 7.80 | −7.40 / 11.18 | −8.32 / 15.99 | −8.29 / 20.36 | −7.91 / 22.83 | −7.14 / 25.59 | −6.07 / 28.08 |
| 青（0 mrad、8 %） | −4.49 / 4.98 | −6.63 / 7.96 | −8.47 / 11.36 | −9.98 / 16.18 | −10.58 / 20.56 | −10.75 / 23.02 | −10.85 / 25.78 | −10.90 / 28.28 |
| 緑（14 mrad、8 %） | −4.48 / 4.98 | −6.61 / 7.96 | −8.39 / 11.35 | −9.64 / 16.18 | −9.60 / 20.56 | −9.08 / 23.02 | −8.06 / 25.78 | −6.71 / 28.28 |

Fig. 5（400 dpi で描画して読んだもの、±0.3 dB・±10 mW 程度）の 1200 mW の端の値（黒 −8.6、破線 −4.8、緑 −6.7、青 −10.75）と合う。
**「前の研究」の破線は 1200 mW で約 −4.8 dB に見え、計算の −6.07 dB と合わない**（他の 3 本は合う）。破線の図は α = 820 %/W（[23] の値）で描かれた可能性があるが、確認できない。§7-2 の「形の再現」の対象は 3 本＋破線の最小値の位置（P_opt = 503 mW、−8.40 dB）とし、破線の右端は使わないことを提案する（§11）。

Fig. 5 のデータ点（読み取り、11 点）と this work の曲線の差（モデル − データ）は反 squeezing で −0.17〜+0.30 dB、squeezing で −0.20〜+0.22 dB。

### 7.3 完全なシミュレーションの許容幅（提案）

完全なシミュレーション = 分布損失の導波路（a_s = 0.1 dB/cm）＋下流の集中損失（1 %、1 %、0.2 %、4 %）＋ロックの時間ループが作る位相ジッタ（θ_rms = 9 mrad に調整）。

| モデル | R′− | R′+ |
|---|---|---|
| 単一モード（§7.1） | −10.15 | 20.56 |
| 分布損失、α = 9.06 /W のまま | −10.17 | 20.29 |
| 分布損失、R′+ を合わせて α を当てはめ直す（9.29 /W、L_wg = 2.01 %） | −10.18 | 20.56 |

提案する許容幅:
- **R′−: −10.15 dB ± 0.10 dB**。根拠: 損失の内訳の取り方による差（上表）0.03 dB、ジッタがガウスでない（ロックの残差は狭帯域の成分を持ちうる）ことによる差は 9 mrad では ⟨sin²δ⟩ が同じなら 1 次で消える、θ_rms の調整精度 ±0.5 mrad が ±0.05 dB に相当。測定の不確かさ ±0.2 dB（Ref A、3σ）より狭くしておけば、実装の誤り（係数 2 など）を見つけられる。
- **R′+: 20.4 dB ± 0.3 dB**（20.1–20.7 dB）。根拠: 集中損失と分布損失の差 0.27 dB はモデルの選び方の不定性で、どちらも許す。測定値 20.2 ± 0.2 dB もこの範囲に入る。
- 単一モードのブロックそのもの（§9 V-T01-1）は解析値に 0.01 dB 以内。
- ESA の zero-span の時系列（VBW 100 Hz）で、ロックしたときの squeezing の揺れの標準偏差は Ref A で 0.07 dB [Ref A §3.2]。これは VBW より遅い位相の揺れと ESA の radiometer 雑音の和で、T13・T04 の照合点に使える（T01 の許容幅には入れない）。

## 8. Ref A の phase-detection OPA に要るもの（D30 の統合案に関係）

### 8.1 Ref A の 2 つのロックの仕組み（[AE] の導出、T14 で確認すること）

構成 [Ref A §3.1、Fig. 2]: 光源（Laser + EDFA + SHG、Fig. 2 に EDFA あり。本文は NKT Koheras Harmonik）→ 基本波は BPF → カプラで LO と probe に分ける。probe は AOM 2 台（正味 +1 MHz）＋ EOM（control 1）→ ダイクロイックミラーでポンプと合わせ → **tapping mirror で squeezer OPA と PD-OPA に分ける**（pump と probe を同時に）。PD-OPA の出口はダイクロイックでポンプと probe を分け、probe をファイバ結合の PD で受けて phase-locking signal 1。LO は EOM（control 2）＋フィルタキャビティ。BHD の出力を ESA と phase-locking signal 2 に分ける。

- **lock 1（PD-OPA）**: probe（ω0 + Ω_m）が PD-OPA で idler（ω0 − Ω_m）を作り、PD の光電流に **2Ω_m（2 MHz）のビート**が出る。
  I(t) ∝ P_pr [cosh 2r_pd + sinh 2r_pd cos(2Ω_m t − (2φ_pr − φ_p) − arg(μν*))]（conventions の符号で、2φ_pr − φ_p = 0 で増幅）。
  2 MHz で復調することは [23] に明記（「probe light induces a differential-frequency-generation process in the OPA, the beat frequency becomes twice the original shift frequency」[23]）。Ref A は周波数を書いていないが、同じ AOM（1 MHz）を使うので同じと考える [AE]。
  ビートの振幅 ∝ P_pr sinh 2r_pd（r_pd ≪ 1 で ≈ 2 r_pd P_pr）。Ref A は PD-OPA の probe を「数 mW」まで上げられることを利点としている [Ref A §2]。
- **lock 2（BHD の 1 MHz）**: squeezer を通った probe（信号 + idler）と LO のビート。復調した値は（lock 1 で決まる位相 ε を使って）
  e(v) = Re[e^{−iε}(e^{r} cos v + i e^{−r} sin v)]、v = φ_LO − φ_p/2。
  ε = 0 なら e = e^{r} cos v で、**零点 v = π/2 が squeezing の測定点、傾きは e^{r}**（増幅された probe なので大きい）。
  ε ≠ 0 のとき零点は v − π/2 ≈ e^{−2r} ε にずれる。数値: r = 2.408、ε = 0.01 rad → 8.1×10⁻⁵ rad、ε = 0.1 rad → 8.1×10⁻⁴ rad。
  すなわち **Ref A では Δ = φ_LO − θ_sq は lock 2 がほぼ決め、lock 1 の誤差は e^{−2r} ≈ 0.008 倍に抑えられる** [AE]。
  lock 1 は probe を増幅側に置いて lock 2 の信号を保つ役割（ε = π/2 なら信号は e^{−r} に落ち、零点は反 squeezing 側に移る）。
- このため、squeezer の結晶の温度による位相（§2.5 の ½ψ_p,xtal、~6 rad/K）と回転 δθ_mm は、probe と squeezed 光に共通なので lock 2 で消える。PD-OPA の結晶の Δk による lock 1 のずれ（arg(μν*)）も e^{−2r} で抑えられる。
- 一方、lock 2 の精度は BHD に入る probe の量（GBS では混入光子になる）と SNR で決まる。Ref A の従来方式で probe が µW に制限されるのは、BHD に入る probe を減らすため [Ref A §2]。

### 8.2 PD-OPA に要るもの

1. squeezer と PD-OPA への tapping mirror から先の、ポンプ・probe の**非共通の光路**が短く安定していること。lock 1 の位相は 2φ_pr − φ_p なので、非共通区間の誤差は 2δψ_pr − δψ_p（θ には半分）。Ref A は「ほぼ同じ光路を通るので保たれる」としている [Ref A §2]。
2. PD-OPA の位相整合（温度）の安定。PD-OPA の出力のビートの位相は Δk_pd により arg(μν*) だけずれる（Δk だけに依り、絶対的な屈折率には依らない）[AE]:

| r_pd | Δkℓ/2 = 0.01 | 0.1 | 0.3 |
|---|---|---|---|
| 0.05 | −0.0100 | −0.0999 | −0.2998 |
| 0.5 | −0.0092 | −0.0924 | −0.2776 |
| 1.0 | −0.0076 | −0.0762 | −0.2291 |
| 2.408 | −0.0041 | −0.0409 | −0.1228 |

  （rad。低利得で ≈ −Δkℓ/2、高利得で ≈ −Δkℓ/(2r)。θ への換算は半分。）
3. PD-OPA の利得 r_pd とポンプパワー、probe パワー（SNR）。高利得・mW の probe ではポンプ枯渇が誤差信号の傾きを変えうる（§2.6）。
4. パルス光源（F）では、PD-OPA でポンプパルスと probe（CW かパルスか）の時間の重なり。CW の probe ならビートはポンプのデューティで平均される。

### 8.3 Liu 型（D30）との統合への示唆（ユーザの検討用。決定はしない）

- Liu 型の逆向き probe は、**振り分け点（参照光と干渉させる点）から GBS 側の分岐点まで**のポートの光路を固定する（D30、D31）。squeezer の結晶とその前後の光路は Liu 型のループの外にある。
- Ref A の lock 2 は「squeezer を通った probe」を BHD で LO と比べるので、squeezer の中の位相がすべて共通になる。**統合案で squeezer に順方向の probe を通さないなら、PD-OPA だけでは squeezer の結晶の位相（~6 mrad/mK）と非共通のポンプ光路は補正されない**。結晶の温度安定に ≲ 1 mK（θ の予算 9 mrad の場合）が要ることになる [AE、dn/dT が UNVERIFIED]。
- 参考: Liu et al. の最初の位相ロック（OPO の光源）は、**OPO を通ったポンプ（775 nm）をダイクロイックで分け、参照光と BHD で干渉させて PZT を制御する**（Liu et al., Methods「Because the pump beam co-propagates with the squeezed light through the cavities, it serves as an ideal probe」）。
  単一パスの導波路 OPA でも、出口のポンプの位相 φ_p,in + ψ_p,xtal は §2.5 の式で θ_sq を（δθ_mm を除いて）決めるので、**squeezer を出たポンプを 2ω0 の参照光（参照光の SHG）と比べる方式**は、1550 nm の probe を GBS 側に入れずに squeezer の結晶の位相を拾える候補になる。残るのは δθ_mm ≈ −Δkℓ/(4r)（10 mK で約 2 mrad）だけ [AE]。
  ただし 2ω0 の参照光を作る SHG とその光路の位相が新しい非共通区間になる。これは D14・D30 の比較候補に加えるかどうか、ユーザの判断（§12）。

## 9. 検証テスト（実装が通すべきもの）

- **V-T01-1（§7-1）**: 単一モードのブロックで L = 0.08、α = 9.06 /W、P = 0.64 W、θ = 0.009 rad → R′− = −10.148 ± 0.01 dB、R′+ = +20.556 ± 0.01 dB。conventions §6 の厳密な式でも同じ（差 < 10⁻³ dB）。
- **V-T01-2（§7-2）**: §7.2 の表の 4 本 × 8 点を ±0.01 dB で再現（解析式の実装の確認）。this work の曲線と Fig. 5 のデータ点の差は ±0.35 dB 以内。最小点: this work P_opt = 612 mW・−10.15 dB、前の研究 P_opt = 503 mW・−8.40 dB。
- **V-T01-3（分布損失）**: ODE（dR±/dz = ±2gR± − a_s(R± − 1)）を数値で解いた値と閉じた式が相対 10⁻⁶ で一致。0.1 dB/cm・45 mm・r = 2.408 で L_wg = 2.03 ± 0.01 %、R− = 0.02821、R+ = 113.74。極限: a_s → 0 で e^{±2r}、g → 0 で R = 1（真空は真空のまま、§7-4）、入力が R_in ≠ 1 で g = 0 なら R_out = e^{−a_sℓ}R_in + 1 − e^{−a_sℓ}（§7-3 の純損失）。
- **V-T01-4（不確定性）**: 損失なしで R+R− = 1、損失ありで R+R− > 1、すべての P と a_s で。
- **V-T01-5（符号）**: conventions §4 の H_int で φ_p = 0.7 rad のとき、分散が最小の角度 = φ_p/2 + π/2（mod π）を 10⁻⁴ rad で、θ = φ_p/2 のシードが e^{r} 倍に増幅。LO を φ_LO = θ_sq に置いた BHD の模擬で R− が出る。φ_OPA,0 を変えたとき θ_sq が φ_OPA,0/2 だけ動く。
- **V-T01-6（最適ポンプ）**: R′− を P で数値最小化した結果が r_opt = ½ ln cot θ、R′−,min = L + (1 − L) sin 2θ と 10⁻⁶ で一致（θ = 9 mrad で P_opt = 612.3 mW）。
- **V-T01-7（位相不整合）**: δ = 0 で §2.1 に戻る。|δ| > g で増幅しない（R+ が有界で振動）。r_eff が §2.5 の cosh(2r_eff) の式と一致。r = 2.408、Δkℓ/2 = 0.01 で |δθ_mm| = 2.04 ± 0.05 mrad（符号はこのテストで固定）。SHG の効率が sinc²(Δkℓ/2) に従い、半値が Δkℓ/2 = 1.392。
- **V-T01-8（ポンプ枯渇）**: 結合波方程式（§2.2）を数値で解き、(a) シードが弱いとき増幅率 e^{2r}、r = √(η_SHG P) を 10⁻⁴ で再現、(b) |A_s|² + |A_p|² の保存を 10⁻⁹ で、(c) SHG（A_p(0) = 0）で P_SH = η_SHG P_FH² を低変換で再現。
- **V-T01-9（ポンプパワーの感度）**: P を +1 % → R′+ +0.104 ± 0.005 dB、R′− +0.002 dB（数値微分と 16.3 dB/W の解析値の一致）。
- **V-T01-10（Ref A の 2 つのロック、T14 と共同）**: PD-OPA の光電流のビートが 2Ω_m に出て、その位相が 2φ_pr − φ_p に従う。lock 2 の零点の ε 依存が e^{−2r} ε（r = 2.408、ε = 0.1 rad で 8.1×10⁻⁴ rad）。squeezer の結晶の温度を変えても（probe が squeezer を通る構成では）Δ が変わらない。PD-OPA だけで squeezer に probe を通さない構成では Δ が ½ψ_p,xtal で動く。

## 10. シミュレータが必要とするもの

| ブロック / 表 | 段 | 入力 | 出力 | 切り替え |
|---|---|---|---|---|
| `opa_waveguide`（単一パス OPA の伝達） | O | η_SHG または α、P と基準（input/output）、ℓ、a_s、a_p、Δk（または T と T_pm、ΔT_FWHM）、φ_OPA,0 | Gaussian チャネル（2×2 のシンプレクティック行列 S と付加雑音行列 N）、r_eff、L_wg、δθ_mm、ψ_p,xtal | `distributed_loss`（オフで集中損失 L_wg）、`pump_loss`、`phase_mismatch`、`pump_depletion`（結合波方程式で解く）、`multimode_spatial` |
| `single_mode_reference`（§7-1 の式） | O | L、α、P、θ | R′± | 検証用。常に残す（保護リスト） |
| OPA の表（P、T の格子） | O | 上と同じ | R±(P, T)、r_eff、δθ_mm | — |
| ポンプの強度雑音源 | T | RIN の PSD、遅いドリフト | P(t) → r(t) | `pump_rin` |
| 結晶の温度 | T | ドリフトの速さ [K/s]（D27）、TEC の制御 | T(t) → Δk(t)、ψ_p,xtal(t) | `crystal_thermal` |
| squeezer の結晶 = ポンプの光路の区間 | T | ℓ、n_e(2ω0)、dn/dT | ψ_p,xtal、θ_sq,out = (φ_p,in + ψ_p,xtal)/2 + π/2 + δθ_mm | （区間モデルの一部） |
| 損失の P 依存 | O | L(P) の関数（Ref B と同じ仕組み） | L_wg(P) | `pump_induced_loss` |
| PD-OPA の誤差信号（Ref A の lock 1） | T | probe のパワー・位相（AOM のずれを含む）、PD-OPA のポンプのパワー・位相、r_pd、Δk_pd | 2Ω_m のビートを含む光電流 | T14 の誤差信号 #1 の一部 |
| 空間モード | O | w_j、r_j、Δθ_j | 実効の R | 既定は単一モード（T12） |

出力の必須項目: R±（θ = 0）、R′±、L_wg、r_eff、θ_sq の内訳（ポンプの光路 / 結晶 / 回転）。D12 の指標 (i)・(ii) を作るのに θ_sq の内訳が要る。

## 11. DECISIONS・conventions との食い違いと提案（調整役の判断用）

1. **PUMP_POWER の定義**: prompt §2 の「pump power (ref: 640 mW)」は出口のポンプの推定値 [Ref A §3.2]。PROCESS_LOG・DECISIONS に「P は OPA 出口のパワー、α = 906 %/W はこの P と集中損失モデルで当てはめた値」と書くことを提案。survey §6-9 の指摘を確認した。
2. **survey §6-9 の EDFA**: Ref A の Fig. 2 には「Laser + EDFA + SHG」と BPF がある。prompt §1 の「CW seed + EDFA + SHG」と整合する。
3. **conventions §4 への追加の提案**: (a) OPA ごとの定数 φ_OPA,0（§2.2）。(b) OPA の結晶は θ_sq に ½ψ_p,xtal で入るポンプの光路の区間として扱う（§2.5）。どちらも符号の約束は変えない。
4. **§7-2 の破線**: 「前の研究」（14 mrad、12 %）の曲線は α = 9.06 /W では 1200 mW で −6.07 dB で、図の破線（約 −4.8 dB）と合わない。破線の右端は §7-2 の検証対象から外すか、α を別に当てはめる。
5. **§7-1 の許容幅**: R′− ±0.10 dB、R′+ 20.4 ± 0.3 dB を提案（§7.3）。

## 12. ユーザの決定が必要な点（推測しない）

1. 完全なシミュレーションで導波路を**分布損失（物理的）と集中損失（Ref A の当てはめと同じ）のどちらを既定にするか**。分布損失を既定にするなら α を 9.29 /W に当てはめ直すか、9.06 /W のまま R′+ が 0.27 dB 低くなるのを受け入れるか。
2. D30 の統合案で、**squeezer の結晶の位相（~6 mrad/mK）をどう扱うか**: (a) 結晶の温度安定（≲ 1 mK）で済ませる、(b) squeezer を出たポンプを 2ω0 の参照光と比べるロック（Liu の OPO のロックに近い）を候補に加える、(c) 順方向の probe を squeezer に通す（GBS への混入光子と引き換え、D28）。D14 の比較候補に (b) を加えるか。
3. 将来の OPA の導波路長（D07 のパルス幅の掃引は 45 mm では準 CW の範囲外）。T03 で多モードとして扱うか、短い導波路を候補にするか。
4. 772.7 nm でのポンプの伝搬損失 a_p と、45 mm 素子の ΔT_FWHM の値（Ref A・[19] に記載なし）。メーカー（NTT）の資料があれば `refs/` に置いてほしい。なければパラメータとして掃引する。

## 13. 文献

読んだもの（本文を確認）:
- **Ref A**: K. Hirota, T. Kashiwazaki, G. Ha, T. Yamashima, P. Jaturaphagorn, T. Suzuki, K. Takahashi, A. Kawasaki, A. Inoue, W. Asavanant, M. Endo, T. Umeki, A. Furusawa, "Generation of 10-dB squeezed light from a broadband waveguide optical parametric amplifier with improved phase locking method," arXiv:2511.15082（2025-11-19 投稿）。**掲載誌・DOI は UNVERIFIED**（`refs/` の PDF にも記載なし）。
- **[19]** T. Kashiwazaki, T. Yamashima, N. Takanashi, A. Inoue, T. Umeki, A. Furusawa, "Fabrication of low-loss quasi-single-mode PPLN waveguide and its application to a modularized broadband high-level squeezer," Appl. Phys. Lett. 119, 251104 (2021). doi:10.1063/5.0063118（arXiv:2201.01457 の PDF に記載）。
- **[23]** T. Kashiwazaki, T. Yamashima, K. Enbutsu, T. Kazama, A. Inoue, K. Fukui, M. Endo, T. Umeki, A. Furusawa, "Over-8-dB squeezed light generation by a broadband waveguide optical parametric amplifier toward fault-tolerant ultra-fast quantum computers," Appl. Phys. Lett. 122, 234003 (2023). doi:10.1063/5.0144385（arXiv:2301.12658 の書誌情報。論文番号 234003 は Ref A の要旨から）。
- **[28]** T. Yamashima, T. Kashiwazaki, T. Suzuki, R. Nehra, T. Nakamura, A. Inoue, T. Umeki, K. Takase, W. Asavanant, M. Endo, A. Furusawa, "All-optical measurement-device-free feedforward enabling ultra-fast quantum information processing," Opt. Express 33(3), 5769–5780 (2025). doi:10.1364/OE.550832（arXiv:2410.20693 の書誌情報。付録の式は arXiv 版で確認）。
- **Liu et al.**: H.-L. Liu et al., "Gaussian boson sampling with 1,024 squeezed states in 8,176 modes," Nature (2026), published online 13 May 2026. doi:10.1038/s41586-026-10523-6（`refs/` の PDF に記載）。Methods「Phase stabilization」を参照。

書誌情報だけを確認したもの（本文は未読）:
- D. H. Jundt, "Temperature-dependent Sellmeier equation for the index of refraction, n_e, in congruent lithium niobate," Opt. Lett. 22, 1553–1555 (1997). DOI は **UNVERIFIED**。dn_e/dT の値は読んでいない（本文の (3–4)×10⁻⁵ /K は二次情報、**UNVERIFIED**）。

Ref A の文献リストからの引用で、本文を読んでいないもの（**UNVERIFIED**）:
- [18] T. Kashiwazaki et al., "Continuous-wave 6-dB-squeezed light with 2.5-THz-bandwidth from single-mode PPLN waveguide," APL Photonics 5, 036104 (2020).
- [16] D. K. Serkland, M. M. Fejer, R. L. Byer, Y. Yamamoto, "Squeezing in a quasi-phase-matched LiNbO3 waveguide," Opt. Lett. 20, 1649 (1995).
- [25]–[27]（R′± の式の出典）: S. Suzuki et al., Appl. Phys. Lett. 89, 061116 (2006); T. C. Zhang et al., Phys. Rev. A 67, 033802 (2003); T. Aoki, G. Takahashi, A. Furusawa, Opt. Express 14, 6930 (2006).
