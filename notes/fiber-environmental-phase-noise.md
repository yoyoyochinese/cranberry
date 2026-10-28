# T08 ファイバの環境による位相雑音（fiber-environmental-phase-noise）

- Phase: 1（deep dive、トピック T08）
- 日付: 2026-10-28
- 担当する効果: `notes/survey.md` の E22（温度による位相ドリフト）、E23（音響・振動）、E24（熱力学的な基底雑音）、E25（ファイバ間の雑音相関）、E26（群遅延の温度ドリフトの係数）、E27 の熱的な部分
- 入力: `docs/squeezing_sim_prompt.md`、`DECISIONS.md`（D02、D04、D07–D13、D27、D29–D32）、`notes/conventions.md`（§1、§4、§7、§9、§10）、Ref C、Liu et al. 2026、`refs/manuals/914-915B-Specifications.pdf`、下記 §12 の一次文献
- 値の印: **[文献名]** = その文献の本文から読んだ値。**[AE]** = この担当による概算（agent estimate、未検証）。**UNVERIFIED** = 原典の本文で確認できていないもの。

## 0. 要約

1. **熱位相係数**: k_T = (2π/λ)(dn/dT + n·α_L) は、Dong et al. 2016 の Table I の物性値（dn/dT = 9.2×10⁻⁶ /K、n = 1.468、α_L = 5.5×10⁻⁷ /K）で 1548.8 nm のとき **40.6 rad/(m·K)**。サーベイの dn/dT = 1.0×10⁻⁵ /K では 43.8 rad/(m·K)。Slavík et al. 2015 は SMF-28 の群遅延の熱係数を 39 ps/(km·K)（文献値）、37.4 ps/(km·K)（実測、1535 nm）とし、位相に直すと「約 48 rad/(m·K)」と書いている。**サーベイの約 44 rad/(m·K) は文献の幅 40–48 rad/(m·K) に入っていて妥当**。既定値は 44 rad/(m·K)、掃引は 40–48 rad/(m·K)。被覆やケーブルの構造によってはさらに大きくなりうる（Hilweg et al. 2022: 「数百 rad/(m·K) もありうる」）。
2. **群遅延の熱係数**: dτ_g/dT ≈ 37.4–39 ps/(km·K) [Slavík 2015]。L = 63 m で 2.4–2.5 ps/K。100 fs のパルスでは 0.04 K で 1 パルス幅ずれる（サーベイ E26 の値を確認）。
3. **機械的な雑音（音響・振動）の実測の目安**: Williams et al. 2008 では、ロックしていないファイバの位相雑音は f < f_c ≈ 1 kHz で S_ψ(f) ≈ h/f²。地中の実リンク（Boulder の BRAN）で h/L ≈ 4 rad²·Hz/km、**巻いたファイバ（spool）は同じ長さで 20–25 dB 低い** → h/L ≈ 0.013–0.04 rad²·Hz/km。これを 63 m の spool に当てると、1 Hz より上の rms は 28–50 mrad（ロックなし）。f_u = 100 Hz の 1 次のロックで残差は約 6 mrad [AE]。
4. **熱力学的な基底雑音**（Wanser の熱伝導雑音 ＋ Duan の熱機械雑音）は、Dong et al. の SMF-28 の実測と合う。1 m あたり S ≈ 1.9×10⁻¹³/f rad²/Hz（熱機械、約 150 Hz より下で支配的）＋ 2.2×10⁻¹⁵ rad²/Hz（熱伝導、約 190 Hz まで平坦）。63 m・0.1 Hz–10 kHz の積分で **27 µrad rms**。どの構成でも無視できる（既定はオフ）。
5. **長さによるロック持続時間**: 温度が r_T [K/s] で変わると位相は k_T·L·r_T [rad/s] で流れる。915B の全レンジ（16 µm の光路長 = 64.9 rad）は、63 m では 0.023 K で使い切り、r_T = 10⁻⁴ K/s なら約 4 分。630 m（f_rep = 10 MHz の L_31）では 23 s。**長い L_k では PZT ストレッチャだけではレンジが足りず、Liu et al. のような TEC の帰還（遅いアクチュエータ）が必要になる**。Ref C でもロックの持続時間は数分で、遅いアクチュエータのレンジで決まっていた。
6. **ファイバ間の相関**: 同じケーブルの 2 本のファイバ（各 700 km、地中）で、補正信号の相関係数は Δt = 0 で 0.98 を超えた。コヒーレンス C(f) ≈ 1 なのは f < 0.002 Hz と 0.04–0.2 Hz で、それ以外の帯域では下がる [Koke et al. 2019]。**相関は周波数に依存する ρ(f) としてモデルに持たせる**。差の PSD は 2S(1 − ρ)（2 本の S が等しいとき）。Ref C は相関を定性的にしか書いていない。
7. **逆向き伝搬との関係（T10 に渡す）**: 区間の中で雑音が位置について無相関なら、同じ時刻に両端から出る順方向と逆方向の位相の差は S_D = 2S[1 − sinc(2πfτ)] ≈ (1/3)(2πfτ)² S（τ = n_g L/c）。これは Williams の遅延で抑圧しきれない雑音の係数 a = 1/3 と一致する。区間の全体で一様にコヒーレントな擾乱（spool 全体の温度など）なら S_D = 0。63 m・1 kHz では抑圧の係数が 1.3×10⁻⁶ なので、**逆向き伝搬による時刻のずれは、我々の長さでは環境雑音についてはほぼ効かない** [AE]。

## 1. 記号（symbols, English）

- `λ`: vacuum wavelength [m]. Ref D: 1548.8 nm (D02). Future F: 1545.32 nm until decided (D07).
- `ω0 = 2πc/λ`: optical carrier angular frequency [rad/s]. `c`: speed of light in vacuum [m/s].
- `n_a`, `n_g,a`: effective phase / group index of polarization axis `a` [-]. Default n = n_g = 1.468 for both axes (UNVERIFIED for PM15-U25D; only used to set scale).
- `β_a = n_a ω0 / c`: propagation constant of axis `a` [rad/m].
- `L_s`: length of path segment `s` [m] (conventions §9). `d_s ∈ {+1, −1}`: propagation direction. `a ∈ {slow, fast}` (PM) or `{x, y}` (SMF).
- `ψ_s,a(t)`: phase accumulated by light on axis `a` traversing segment `s` [rad] (conventions §4: increases when the segment gets longer or hotter).
- `T0`: nominal temperature, default 295 K (D27). `T_amb,z(t)`: ambient temperature of thermal zone `z` [K]. `ΔT_s(t)`: temperature change of segment `s` relative to T0 [K].
- `r_T`: ambient temperature drift rate [K/s] (D27 parameter).
- `τ_th,s`: thermal time constant between ambient and fiber of segment `s` [s].
- `dn_a/dT`: thermo-optic coefficient [1/K]. `α_L = (1/L) dL/dT`: linear thermal expansion coefficient of the fiber [1/K].
- `k_T,a = (2π/λ)(dn_a/dT + n_a α_L)`: thermal phase coefficient per unit length [rad/(m·K)].
- `κ_τ = dτ_g/dT / L = (dn_g/dT + n_g α_L)/c`: thermal coefficient of group delay per unit length [s/(m·K)].
- `ε_s(t)`: axial strain of segment `s` (dimensionless), the driver for mechanical (acoustic/vibration) noise.
- `ξ_eo`: photo-elastic correction factor, δψ = ξ_eo β δL (conventions §4, value owned by T16). `k_ε,a = ξ_eo,a β_a`: strain phase coefficient [rad/m per unit strain].
- `S_x(f)`: one-sided PSD of x (conventions §7). `f`: Fourier frequency [Hz]. `ω = 2πf` only inside formulas quoted from sources that use ω.
- `h`: coefficient of the f⁻² mechanical phase-noise model, S_ψ = h/f² [rad²·Hz]. `h_L = h/L` [rad²·Hz/m] (per-length form, γ = 1).
- `L_ref`: reference length used to normalize `h` when the length exponent `γ ≠ 1` [m]. `γ ∈ [1, 2]`: length exponent of the PSD (1 = spatially uncorrelated noise along the fiber, 2 = spatially coherent).
- `f_c`: upper corner frequency of the f⁻² region [Hz]. `f_lo`: lower limit of the mechanical model; below it the thermal model takes over [Hz].
- `ρ_ij(f)`: correlation coefficient (real, normalized cross-spectrum) between drivers of segments `i` and `j`; `C_ij(f) = |S_ij|²/(S_ii S_jj) = ρ_ij²`: magnitude-squared coherence.
- `τ_s = n_g L_s / c`: one-way transit time of segment `s` [s].
- Thermodynamic noise: `k_B` Boltzmann constant; `κ_th` thermal conductivity [W/(m·K)]; `D_th` thermal diffusivity [m²/s]; `k_max = 2/r_mf`, `k_min ≈ 2.405/a_clad` [1/m] (`r_mf` mode-field radius, `a_clad` fiber outer radius); `E0` bulk modulus [Pa]; `A` fiber cross-section [m²]; `φ0` mechanical loss angle [-].
- `R_act`: phase range of an actuator [rad]. `t_sat`: time until a drift exhausts `R_act` [s].
- `f_u`: unity-gain frequency of the phase lock [Hz] (owned by T16; used here only for estimates).

## 2. 支配方程式

### 2.1 区間の位相を環境のドライバで表す（conventions §9、D11）

区間 s・偏光の軸 a の位相を、区間で共通の物理的なドライバ（温度 ΔT_s、軸方向の歪み ε_s）と、軸ごとの係数の積で書く。

```
ψ_s,a(t) = ψ̄_s,a + k_T,a L_s ΔT_s(t) + k_ε,a L_s ε_s(t) + δψ^thd_s,a(t)        (2.1)
```

- ψ̄_s,a = β_a L_s: 名目値（T0 のとき）。
- 第 2 項: 温度（§2.2、§2.4）。符号は conventions §4 のとおり、温度が上がると ψ は増える（k_T > 0）。
- 第 3 項: 機械的な雑音（音響・振動・地面の揺れ）を歪み ε_s で表す（§2.5）。ストレッチャの作用（T16）と同じ係数 k_ε,a を使うので、補正と雑音の比が自動的に正しくなる。
- 第 4 項: 熱力学的な基底雑音（§2.6）。軸ごとに独立に生成してよい（オプション）。
- 「温度は区間で共通、応答は軸で異なる」（D11）は、ΔT_s と ε_s を区間で 1 本だけ生成し、軸ごとに k_T,a・k_ε,a を変えることで実現する。

文献の位相雑音 PSD（rad²/Hz）から歪みの PSD への変換は S_ε = S_ψ / (k_ε L)²。文献の値は「位相」で与えられているので、実装では基準の軸（slow、SMF では x）について位相で PSD を与え、ε に直してから軸ごとの係数をかけてもよい。

### 2.2 熱位相係数（E22）

```
k_T,a = (1/L) dψ_a/dT = (2π/λ) (dn_a/dT + n_a α_L)        (2.2)
```

- 仮定: ファイバの温度が一様、被覆による付加的な伸び（コーティングの熱膨張がガラスを引っ張る効果）は α_L に含める（有効な α_L）。
- 数値（λ = 1548.8 nm、n = 1.468）:

| 出典 | dn/dT [1/K] | α_L [1/K] | k_T [rad/(m·K)] |
|---|---|---|---|
| Dong et al. 2016 の Table I（理論計算に使った値） | 9.2×10⁻⁶ | 5.5×10⁻⁷ | **40.6** |
| サーベイ（UNVERIFIED の値） | 1.0×10⁻⁵ | 5.5×10⁻⁷ | 43.8 |
| Slavík et al. 2015: 39 ps/(km·K) ≈ 48 rad/(m·K)（本文の換算） | — | — | 約 48（2π c κ_τ/λ = 47.4） |
| 同、SMF-28 の実測 37.4 ps/(km·K)（1535 nm）を同じ換算で | — | — | 45.5 |

- Slavík の換算は、位相の係数と群遅延の係数を同じとみなしている（dn_g/dT ≈ dn/dT）。両者の差は dn/dT の波長分散で決まり、ここでは確認していない（UNVERIFIED、差は数 % と見込む [AE]）。
- **結論**: 既定値 k_T = 44 rad/(m·K)、掃引の範囲 40–48 rad/(m·K)。被覆やケーブルによっては「数百 rad/(m·K)」もありうる [Hilweg 2022]。被覆つき PM ファイバ（PM15-U25D）の値は見つけられなかった（UNVERIFIED）ので、同じ範囲を使う。
- サーベイの約 44 rad/(m·K) は**妥当**。L = 63 m では 2.8×10³ rad/K。

軸ごとの差（PM の複屈折の温度依存）:

```
k_T,slow − k_T,fast = (2π/λ) (dB/dT + B α_L),   B = n_slow − n_fast        (2.3)
```

値は T09 の担当。D11 の概算 2.4 rad/(m·K)（UNVERIFIED）を T09 の結果が出るまでの仮の値とする。k_T,a は k_T,slow = 44 rad/(m·K)、k_T,fast = k_T,slow − 2.4 rad/(m·K)（符号は T09 で確認する）と置く。

### 2.3 群遅延の熱係数（E26）

```
dτ_g,s/dT = κ_τ L_s,   κ_τ = (dn_g/dT + n_g α_L)/c ≈ 37.4–39 ps/(km·K)        (2.4)
```

- [Slavík 2015]: SMF-28 で 39 ps/(km·K)（文献値として引用）、37.4 ps/(km·K)（同論文の実測、1535 nm）。既定値は 39 ps/(km·K)。
- 区間の群遅延 τ_g,s(t) = n_g L_s/c + κ_τ L_s ΔT_s(t) を、時間ループ段のタイミングの状態として出す。帰結（モードの重なり・パルス到着時刻）は T03 の担当。
- **位相のロックと群遅延の関係** [AE]: 熱による擾乱では Δτ_g/Δψ = κ_τ/k_T ≈ 1/ω0。PZT ストレッチャ（歪み）による補正では Δτ_g/Δψ ≈ n_g/(n ξ_eo ω0) ≈ 1.28/ω0（ξ_eo = 0.78 のとき）。したがって熱ドリフトをストレッチャで位相だけ打ち消すと、群遅延は (1 − 1/ξ_eo) Δψ/ω0 ≈ −0.28 Δψ/ω0 だけ残る。63 m・1 K（Δψ ≈ 2.8×10³ rad）で約 0.6 ps。TEC による補正は擾乱と同じ熱の機構なので、位相と群遅延を同時に打ち消す。100 fs 級のパルスではこの差が効きうる（T03・T16 に渡す）。

### 2.4 温度の時間モデル（D27）

ゾーン z（同じ空気・同じ箱にある区間の集まり）の周囲の温度:

```
T_amb,z(t) = T0 + r_T,z t + w_z(t),   S_w(f) = b_T / f²   (f ≥ f_min)        (2.5)
```

区間 s の温度は 1 次の遅れで周囲に追従する:

```
τ_th,s dΔT_s/dt = (T_amb,z(s)(t) − T0) − ΔT_s(t)        (2.6)
```

- r_T [K/s] は D27 の掃引パラメータ。既定 10⁻⁴ K/s、範囲 10⁻⁶–10⁻³ K/s [AE、一次資料なし]。w_z はランダムウォーク（b_T [K²·Hz]、既定 0 = オフ）。周期的な空調の成分（振幅・周期）を足せるようにしておくとよい [AE]。
- 定常のドリフトでは ΔT_s は r_T (t − τ_th,s) に追従し、位相は k_T,a L_s r_T [rad/s] で流れる。
- τ_th,s は巻き方・ケーブル・箱で大きく変わる（秒〜時間）[AE]。既定 300 s、範囲 1–3000 s。
- 2 本の区間が同じゾーンにあり τ_th も同じなら、ΔT は低周波で完全に相関する（§2.8）。τ_th が違うと、共通の温度変化から差が生まれる（|H_i − H_j|² の分だけ）。

### 2.5 機械的な雑音（音響・振動、E23）

一次資料には実験室のファイバに使える普遍的なモデルはない（Hilweg 2022: 10–100 Hz は人の活動が支配的で、時間と場所で大きく変わり、モデル化しにくい）。そこで、実測で確かめられている形を基本にした**パラメトリックなモデル**を提案する。

```
S_ψ,mech(f; L) = h_ref (L / L_ref)^γ · f⁻² · 1 / [1 + (f/f_c)²] + Σ_i P_i(f)    (f_lo ≤ f)        (2.7)
```

- f⁻² の領域と f_c: Williams et al. 2008 は、ロックしていない 251 km のリンクで S_fiber(f) ≈ h/f²（f < f_c ≈ 1 kHz）、f > f_c ではより急に下がる、と報告している。ここでは f_c より上を f⁻⁴ にした（形は [AE]）。
- 振幅の目安（Williams 2008）:
  - 地中の実リンク（BRAN、建物のダクトを約 1 km 通る）: h/L ≈ 4 rad²·Hz/km（251 km で h ≈ 10³ rad²·Hz）。
  - **巻いたファイバ: 同じ長さで 20–25 dB 低い** → h/L ≈ 0.013–0.04 rad²·Hz/km。
  - この値は 50 km 級の spool（多数の spool をつないだもの）での値。長さのスケーリング γ = 1 は「位置について無相関な雑音」を仮定している（Williams の Appendix A、式 (A4)）。
- 既定値: γ = 1、L_ref = 1 km、h_ref = 0.04 rad²·Hz（spool の上限）、f_c = 1 kHz、f_lo = 0.1 Hz（これより下は §2.4 の熱モデルが担う）。
- 掃引: h_ref = 10⁻³–4 rad²·Hz（静かな spool から実リンクまで）、γ ∈ {1, 2}、f_c = 100 Hz–10 kHz。
  - **γ = 2 を残す理由** [AE]: 実験室の 1 つの spool は音の波長（1 kHz で 34 cm）や振動の空間スケールに比べて小さいので、spool 全体が同じ擾乱を受ける（空間的にコヒーレント）。このとき位相の振幅は L に比例し、PSD は L² に比例する。Williams の値は L_ref = 1 km で γ = 1 と γ = 2 が一致するように正規化して使う。
- P_i(f): 機械的な共振や電源周波数のピーク（Lorentz 型、中心 f_i・Q_i・面積 A_i）。既定はなし。Dong 2016 は約 1 Hz（地面）と 40 Hz（音）にピークを見ている。
- 補助の物理モデル（下限・上限の見積もりに使う）:
  - 圧力（音圧）の直接の感度: δψ/(L δp) ≈ 3×10⁻⁵ rad/(Pa·m)（λ = 1.55 µm）[Hilweg 2022]。63 m・20 mPa rms（約 60 dB SPL）で 38 µrad rms。実際の spool の音響雑音は、巻き枠・台の振動を通る経路の方がずっと大きい [AE]。
  - 縦置きの spool の重力加速度による歪み: δL/L = (ν_S ρ_S h_S / 4E_S) δa、S_ψ = 2C_cmrr (k n L ν_S ρ_S h_S/(4E_S))² S_a(f) [Hilweg 2022 の式 (26)]。巻き枠の値がないので既定では使わない。

### 2.6 熱力学的な基底雑音（E24）

熱伝導雑音（Wanser の式）[Hilweg 2022 の式 (5)・(7)、Dong 2016 の式 (1)]:

```
S_ψ^TC(f) = (2π k_B T² L / (κ_th λ²)) (dn/dT + n α_L)² ln[ (k_max⁴ + (ω/D_th)²) / (k_min⁴ + (ω/D_th)²) ],  ω = 2πf        (2.8)
```

熱機械雑音（Duan、低周波の近似）[Dong 2016 の式 (2)、Hilweg 2022 の式 (9) を位相・f に直したもの]:

```
S_ψ^TM(f) = (2π n/λ)² · 2 k_B T L φ0 / (3π E0 A f)        (2.9)
```

- 物性値は Dong 2016 の Table I: T = 300 K、dn/dT = 9.2×10⁻⁶ /K、n = 1.468、α_L = 5.5×10⁻⁷ /K、κ_th = 1.37 W/(m·K)、D_th = 8.2×10⁻⁷ m²/s、k_max = 3.846×10⁵ /m、k_min = 3.848×10⁴ /m、E0 = 1.9×10¹⁰ Pa、A = 4.91×10⁻⁸ m²、φ0 = 0.01。
- 規約の確認: 上の 2 式は片側・Hz あたりの PSD として扱う。根拠: (2.9) に Table I の値を入れると 1 m あたり 3.34×10⁻¹³/f rad²/Hz。Dong は理論から 2.5 dB を引くと実測に合うと書いており、引いた値は 1.88×10⁻¹³/f。一方、実測（ポリイミド被覆 103 m の Michelson、往復と 2 本の腕で実効 412 m、2.5×10⁻¹¹/f）を SMF-28 に戻す（+5 dB）と 1 m あたり 1.92×10⁻¹³/f で、0.1 dB の範囲で一致する。（2.8）の両側/片側は本文から確認できていない（UNVERIFIED）が、どちらでも結論（無視できる）は変わらない。
- 1 m あたりの値（300 K、1550 nm）: S^TC = 2.2×10⁻¹⁵ rad²/Hz（f ≲ k_min² D_th/2π ≈ 190 Hz で平坦、1 kHz で 1.4×10⁻¹⁵、10 kHz で 3.7×10⁻¹⁶）、S^TM = 3.3×10⁻¹³/f rad²/Hz（理論）、1.9×10⁻¹³/f（実測に合わせた値）。両者が等しくなるのは約 150 Hz。
- 長さのスケーリング: L に比例（位置について無相関）。Dong は 327/134/68/34 m で 0、−3.9、−6.8、−9.8 dB を確認している。

### 2.7 区間の中の空間構造と逆向き伝搬（T10 への入力）

区間 s の中の位置 x ∈ [0, L] での単位長さあたりの位相の擾乱を δψ'(x, t) とする。順方向（d = +1）の光が時刻 t に x = L から出るときと、逆方向（d = −1）の光が同じ時刻 t に x = 0 から出るときの位相は（conventions §9）

```
ψ_fwd(t) = ∫_0^L δψ'(x, t − τ + x n_g/c) dx,    ψ_bwd(t) = ∫_0^L δψ'(x, t − x n_g/c) dx        (2.10)
```

- 位置について無相関で一様（S_δψ'(f) = S(f)/L）なら、相互スペクトルは S(f)·sinc(2πfτ)（sinc x = sin x / x）。差の PSD は

```
S_D(f) = 2 S(f) [1 − sinc(2πfτ)] ≈ (1/3) (2πfτ)² S(f)   (2πfτ ≪ 1)        (2.11)
```

  低周波の係数 1/3 は、Williams 2008 の「遅延で抑圧されない雑音」S_D ≈ a(2πfτ)² S_fiber、a = 1/3（雑音が一様に分布するとき）と一致する。（2.10）からの導出はこの担当によるもの [AE] で、Williams の式 (A4)・(A10)–(A11) と整合する。
- 区間全体で一様にコヒーレントな擾乱（δψ'(x, t) = g(t)）なら、ψ_fwd(t) と ψ_bwd(t) はどちらも g の同じ時間窓 [t − τ, t] の平均なので、S_D = 0。
- 1 点 x0 だけの擾乱なら S_D = 2S[1 − cos(2πf(τ − 2x0 n_g/c))]。
- したがって、逆向きの probe による相殺の崩れは「全体の PSD」だけでなく「空間の構造」に依存する。実装では、区間を N_sub 個の小区間に分け、小区間ごとにドライバを持たせ（空間の相関は「lumped（全部同じ）」か「distributed（独立）」か、その混合）、逆向きの光は小区間ごとに異なる時刻で評価する（遅延バッファ）。
- 数値の目安: τ = 0.31 µs（63 m）、3.1 µs（630 m）、4.9 µs（1 km）。f = 1 kHz での抑圧 (1/3)(2πfτ)² は 1.3×10⁻⁶（63 m）、1.3×10⁻⁴（630 m）、3.2×10⁻⁴（1 km）。

### 2.8 ファイバ間の相関（E25、D09、D13）

区間 i と j（同じ束・同じ電線管・同じ spool にあるファイバ）のドライバを、共通の成分と独立な成分に分けて生成する:

```
X_i(f) = √S_i(f) [ √ρ_ij(f) G_c(f) + √(1 − ρ_ij(f)) G_i(f) ]        (2.12)
```

G は独立な単位の白色ガウス雑音（周波数領域）。このとき S_ij = ρ_ij √(S_i S_j)、コヒーレンス C_ij = ρ_ij²。2 本の差（例: squeezed と LO_k）の PSD は

```
S_Δ(f) = S_i + S_j − 2 ρ_ij(f) √(S_i S_j)   →   2S (1 − ρ)   (S_i = S_j = S)        (2.13)
```

- ρ(f) のモデル（提案）: ρ(f) = ρ_0 / [1 + (f/f_ρ)²]^{1/2}（低周波で ρ_0、f_ρ より上で下がる）、ρ_0 ∈ [0, 0.999]、f_ρ ∈ [0.01 Hz, 1 kHz]。形は [AE]。
- 実測の目安: Koke et al. 2019（同じ地中ケーブルの 2 本、各約 700 km）: Δt = 0 での相関係数 > 0.98、C(f) ≈ 1 は f < 0.002 Hz と 0.04–0.2 Hz、0.002–0.04 Hz と 0.2 Hz より上で下がる（後者の原因は不明と書かれている）。距離・環境が我々とまったく違うので、**形（周波数とともに相関が崩れる）だけを参考にする**。
- Ref C: 2 本のファイバは別々の被覆で、空気に囲まれて硬い電線管の中に並べ、断熱チューブで一部を覆った。地震の揺れ・温度変化・音で両方に余分な位相雑音が乗り、「相関した位相雑音は平衡検出で打ち消される」が残りの差を量子雑音ロックで補償した。通信用のファイバ束ならファイバが機械的にずっと強く結びついているので、差の雑音は大幅に小さいはず、とも書いている。**数値はない**。
- 熱のドライバの相関は §2.4 のゾーンと τ_th で表し、機械のドライバの相関は（2.12）で表す。長さが違う 2 本（2 ポートの組、D13）は、並んで走る部分だけを相関する区間にし、残りを独立な区間にする。

### 2.9 Ref C・Liu 型のロックとの関係（規約の変換）

- Ref C の BHD が測る量は Δ = φ_LO − θ_sq（conventions §5、§10 (i)）。squeezed 光と LO はそれぞれ 1 km のファイバを順方向に通るので、環境の寄与は δΔ(t) = ψ_LO-fiber(t) − ψ_sq-fiber(t)（θ_sq はファイバの位相 ψ をそのまま受け継ぐ。縮退の squeezed 光の直交位相の角度は場の位相と同じだけ回る。符号は T04 で確認すること）。PSD は（2.13）。
- Liu 型（D30）: probe は逆向き（d = −1）・直交偏光で同じ区間を通る。順方向の squeezed 光（軸 a_sq）と逆向きの probe（軸 a_pr）の差は、（2.1）で軸の係数の差 (k_T,a_sq − k_T,a_pr) L ΔT（T09・T10）と、（2.11）の時刻のずれの両方から生まれる。

## 3. 仮定と適用範囲

- 線形: 環境の擾乱は光の位相に線形に足される。パワーに依存しない（吸収による加熱は無視 [AE]）。
- ガウス性: 機械・熱力学の雑音はガウス過程として生成する。実際の実験室の雑音は非定常（人の出入り、空調の周期）で、ピークがある。Hilweg 2022 は「都市環境のファイバの短い時間スケールの雑音はガウス分布」と書いている。
- k_T は 295 K 付近で一定とする（±数 K の範囲）。T0 を大きく変える設定では dn/dT(T) が必要（UNVERIFIED、範囲外）。
- Williams の数値は通信用 SMF・km 級の spool と地中リンクのもの。我々の 1–600 m の実験室のファイバ（PM、FC/APC のコネクタ、光学定盤の上）にそのまま当てはまる保証はない。**h_ref は 3 桁以上掃引する前提のパラメータ**で、Williams の値は「目安の 1 点」。
- 熱力学の式（2.8）・（2.9）は SMF-28 の物性値。PM ファイバ（応力付与部をもつ）の E0、A、φ0 は違うが、桁は同じと見込む [AE]。
- 偏波のドリフト（SMF の偏光状態の揺れ、可視度の低下）は T09 の担当。ここでは各軸の位相だけを扱う。

## 4. 効く周波数帯

| 成分 | 帯域 | 形 | 出典 |
|---|---|---|---|
| 温度のドリフト（E22） | DC – 1/(2π τ_th)（τ_th = 300 s なら 0.5 mHz） | ランプ ＋ 1/f² の温度に 1 次 LPF | §2.4 [AE] |
| 温度の速い揺れ（空気の流れ、裸のファイバ） | 〜0.01–1 Hz | 同上で τ_th が小さい場合 | [AE] |
| 機械（音響・振動、E23） | 0.1 Hz – f_c（≈1 kHz）は f⁻²、それより上は急に下がる | （2.7） | Williams 2008 |
| 地面の揺れ（microseism など） | 0.02–1 Hz、人の活動は 10–100 Hz | ピーク | Hilweg 2022 |
| 熱力学（E24） | 熱機械は ≲150 Hz で 1/f、熱伝導は ≲190 Hz で平坦、それより上で下がる | （2.8）（2.9） | Dong 2016 |
| 群遅延のドリフト（E26） | 温度と同じ（< 1 Hz） | — | Slavík 2015 |
| ファイバ間の相関の崩れ（E25） | 相関は低周波で高く、f_ρ より上で下がる | （2.13） | Koke 2019（形のみ） |

## 5. スケーリング

| 量 | ファイバ長 L | パワー | 温度 T0 | パルス幅 τ_p |
|---|---|---|---|---|
| 熱ドリフトの位相 k_T L ΔT | ∝ L | 依存しない | dn/dT(T) を通じて弱く | 依存しない |
| 群遅延のドリフト κ_τ L ΔT | ∝ L | 依存しない | 弱く | 許容される ΔT ∝ τ_p / (κ_τ L) |
| 機械的な雑音 PSD | ∝ L^γ（γ = 1: 位置について無相関、γ = 2: spool 全体がコヒーレント） | 依存しない | 依存しない | 依存しない（kHz 以下の雑音はパルス列に共通） |
| 熱伝導雑音（Wanser） | ∝ L | 依存しない | ∝ T²·(dn/dT + nα_L)² | 依存しない |
| 熱機械雑音（Duan） | ∝ L | 依存しない | ∝ T | 依存しない |
| 逆向き伝搬の抑圧されない雑音 | (1/3)(2πfτ)² S ∝ L²·S（γ = 1 なら ∝ L³） | 依存しない | — | — |
| アクチュエータのレンジを使い切る時間 t_sat | ∝ 1/L | — | — | — |

## 6. 我々のパラメータでの概算

既定値: k_T = 44 rad/(m·K)、κ_τ = 39 ps/(km·K)、spool の機械雑音 h/L = 1.3×10⁻⁵–4×10⁻⁵ rad²·Hz/m（γ = 1）、熱力学は 295 K、915B のレンジ R = 2π × 16 µm / 1548.8 nm = **64.9 rad**（serial 8669-1。8669-2 は 17 µm で 69.0 rad）[914-915B-Specifications.pdf]。

ロックの残差の目安には、1 次のループ（|1/(1 + G)|² = f²/(f² + f_u²)）で S = h/f² を抑えた σ² = π h/(2 f_u) を使う [AE]。

### 6.1 Ref D（D02: PM ファイバの MZI、片腕に 915B）

- 915B の巻線: 38 巻、PM15-U25D、FC/APC（slow 軸に keyed）[manual]。巻線部のファイバ長は資料に書かれていない。Photonics.com の製品掲載の検索結果の抜粋に「9 m のファイバで最大 16 µm の伸び」とある（**UNVERIFIED**、ページは 403 で開けなかった）。ここでは L_915B ≈ 9 m を仮の値とし、腕の長さの差 ΔL は D02 どおり 0–10 m で掃引する。
- 熱: 2 本の腕が同じ温度変化を受けると、MZI の位相は k_T ΔL ΔT。ΔL = 9 m で 395 rad/K。915B のレンジ 64.9 rad は一様な 0.16 K で使い切る。r_T = 10⁻⁴ K/s で t_sat ≈ 27 分、10⁻³ K/s で約 2.7 分。2 本の腕の温度が違えば k_T (L_A ΔT_A − L_B ΔT_B) で、ΔL = 0 でも残る。
- 機械（ロックなし）: 9 m の spool なら h = 1.1×10⁻⁴–3.6×10⁻⁴ rad²·Hz、1 Hz より上で 11–19 mrad rms、0.1 Hz より上で 34–60 mrad rms。2 本の腕が独立なら両方の長さを足す。定盤の上の裸の PM ファイバは spool より雑音が大きい可能性がある（掃引）。
- 熱力学: 9 m で 10 µrad rms（0.1 Hz–10 kHz）、無視できる。
- **D04（ブラインド検証）**: ここの値は文献の値で決めた予測で、Ref D の実測には合わせない。

### 6.2 Ref C（1 km × 2 本、うち約 180 m は建物間に敷設、残りは研究室で巻いてある）

- 熱: 1 km で 4.4×10⁴ rad/K（1 本あたり）。2 本の差は ΔT の差と τ_th の差で決まる（値は論文にない）。
- ロックの結果: 全体の残留位相雑音は約 0.1 rad（フィットの値）、うち dither が ±4.5° → 0.056 rad rms。二乗で引くと、**dither 以外（環境の残差＋その他）は 0.083 rad 以下**。Ref C のプリセットはこの範囲に入る必要がある（§9 V8）。
- ロックの持続時間: 数分。遅いアクチュエータ（ファイバの位相シフタ）のレンジで制限された [Ref C §4]。位相シフタのレンジは論文に書かれていないので、差のドリフトの速さ ≈ R_fs / t_lock（t_lock = 数分）という関係だけが言える。R_fs は T16 で調べる。
- 機械（ロックなし）の目安: 敷設部 180 m を Williams の実リンク（4 rad²·Hz/km）、巻いた 820 m を spool（0.04 rad²·Hz/km）とし、2 本が無相関（ρ = 0）なら差の h ≈ 2 × (0.72 + 0.033) = 1.5 rad²·Hz。f_u = 1 kHz で 49 mrad、10 kHz で 15 mrad [AE]。これは 0.083 rad の上限と矛盾しない。ρ を上げれば小さくなる。
- 熱力学: 1 km で 0.1 mrad rms。

### 6.3 将来の構成 F（D07–D13、D30–D32）

L_k = k c T_rep / n_g（D08）。L_31 = 6.3 m（f_rep = 1 GHz）、63 m（100 MHz）、633 m（10 MHz）。f_rep は Open なので 63 m と 633 m を代表にする。

| 量 | L = 63 m | L = 633 m |
|---|---|---|
| 熱位相 k_T L | 2.8×10³ rad/K | 2.8×10⁴ rad/K |
| r_T = 10⁻⁴ K/s での位相の流れ | 0.28 rad/s | 2.8 rad/s |
| 915B 級のレンジ（64.9 rad）を使い切る ΔT | 0.023 K | 0.0023 K |
| 同、t_sat（r_T = 10⁻⁴ / 10⁻⁵ K/s） | 235 s / 2350 s | 23 s / 235 s |
| 群遅延のドリフト κ_τ L | 2.5 ps/K | 25 ps/K |
| 100 fs のパルスを 1 パルス幅ずらす ΔT | 0.04 K | 0.004 K |
| 機械（spool、ロックなし、f > 1 Hz） | 28–50 mrad | 89–160 mrad |
| 同、f_u = 100 Hz / 1 kHz のロック後 | 6.3 / 2.0 mrad | 20 / 6.3 mrad |
| 熱力学（0.1 Hz–10 kHz） | 27 µrad | 85 µrad |
| 逆向き伝搬の抑圧の係数 (1/3)(2πfτ)² @ 1 kHz | 1.3×10⁻⁶ | 1.3×10⁻⁴ |

- 比較: Ref A の θ = 9 mrad（prompt §7-1）。機械的な雑音は f_u ≳ 100 Hz のロックで数 mrad に下がる見込み。**支配的なのは熱ドリフトがアクチュエータのレンジを超えること（ロック外れ）と、直交偏光の probe のときの軸間の熱ドリフト（T09・T10）**。
- squeezed のファイバと LO_k のファイバ（D09）: 長さが揃っていて同じ環境なら、差は（2.13）で 2S(1 − ρ)。熱の差は ΔT と τ_th の差で決まる。ρ はパラメータ（D09）。
- 指標 (ii)（D31）: ポートごとにロックするなら、ポート間の差は各ポートのロックの残差の差（ポート間の相関は残差の相関を通じてだけ効く）。共通点の揺れ（probe の分岐、参照光の配送路、D31 の 2）も区間として同じモデルで雑音を持たせる。

### 6.4 Liu et al. 2026 の扱い（D29）

- 構成: 時間間隔 τ = 50 ns、1 段目の遅延ループは t, t + τ, …, t + 15τ、2 段目は t, t + 16τ, …, t + 15 × 16τ [Liu Methods]。全部をファイバ（n_g = 1.468）とすると 1 段あたり約 10.2 m、1 段目の最長 15τ で約 153 m、2 段目の最長 240τ で約 2.45 km（**この換算は担当によるもの [AE]**。論文はファイバ長を書いていない）。
- ロックの方法: 空間モードごとに PZT ミラー。probe は周波数が同じで偏光が直交、逆向きに伝搬し、AOM でパルス化して squeezed のパルス列と時間をずらす。1 段目の干渉計の出口で PBS で取り出し、参照光と干渉させて PZT を駆動する [Liu Methods]。
- **長いファイバの温度ドリフト**: 「温度の揺らぎによる長いファイバの大きな位相ドリフトは、しばしば 1 つの PZT の補償レンジを超える。そこで遅い熱ドリフトを抑えるために熱電素子（TEC）の帰還を併用した」[Liu Methods]。TEC の帯域・制御方法・TEC で温めるファイバの長さは書かれていない（→ パラメータ）。
- 達成値: 位相ロックの精度は OPO 光源の部分で 1 nm、1 段目の遅延ループで 3 nm、2 段目で 10 nm（1550 nm で 4.1、12、41 mrad）、数十時間維持。系全体で約 λ/200（31 mrad）[Liu 本文、Methods]。
- Ref C も同じ方向の改良（「2 本のうち一方の数 m を温度で動かす」）を提案している [Ref C §4]。
- シミュレータへの帰結: 遅いアクチュエータとして**熱アクチュエータ**（TEC: 入力 → 1 次の遅れ τ_TEC → ΔT_TEC → k_T L_TEC ΔT_TEC の位相、レンジは ΔT_TEC の上下限）を用意する。係数は（2.2）（2.4）、ブロックの所有は T16。

## 7. シミュレータに必要なもの

### 7.1 ブロックと表

| 名前（案） | 段 | 役割 | 入力 | 出力 |
|---|---|---|---|---|
| `fiber_env_coeffs` | O | ファイバの種類・軸ごとの係数の表 | fiber type（PM15-U25D / SMF-28 / custom）、λ | k_T,a、k_ε,a、κ_τ、n_a、n_g,a、出典の印 |
| `ambient_temperature` | T（雑音源） | ゾーンごとの周囲温度（2.5） | T0 = 295 K、r_T、b_T、（周期成分）、seed | T_amb,z(t) |
| `segment_thermal` | T（雑音源） | 区間の温度（2.6） | ゾーンの割り当て、τ_th,s | ΔT_s(t) |
| `segment_mech_noise` | T（雑音源） | 機械的な雑音（2.7）。区間間の相関（2.12）、区間内の空間構造（N_sub、lumped/distributed） | h_ref、L_ref、γ、f_c、f_lo、ピークの表、ρ_ij(f)（ρ_0、f_ρ）、seed | ε_s(t)（または小区間ごとの ε_s,m(t)） |
| `segment_thermodynamic` | T（雑音源） | 熱力学の基底（2.8）（2.9） | 物性値（Dong の Table I が既定）、T0 | δψ^thd_s,a(t) |
| `segment_phase` | T | （2.1）で ψ_s,a(t) と τ_g,s(t) を組み立てる | 上の出力、係数の表 | ψ_s,a(t)（区間 × 軸）、τ_g,s(t)（タイミングの状態、T03 へ） |
| `thermal_actuator`（TEC） | T（アクチュエータ） | 遅い熱の補正（Liu 型、Ref C の提案） | 制御入力、τ_TEC、L_TEC、ΔT_TEC の範囲 | ΔT_TEC → 区間の ψ と τ_g（所有は T16） |

### 7.2 トグル（すべて独立にオン・オフ、保護リスト）

- `env.thermal.drift`（ランプ r_T）、`env.thermal.random`（b_T）、`env.mech`、`env.mech.peaks`、`env.thermodynamic`（既定オフ）、`env.correlation`（オフなら区間間は独立）、`env.spatial`（lumped / distributed / mixed）、`env.group_delay`（τ_g の更新）。

### 7.3 数値の注意（実装担当へ）

- 時間スケールが 6 桁以上違う（熱は 10²–10³ s、機械は〜kHz）。熱の成分は低いレート（例: 10–100 Hz）で積分して補間し、機械の成分は FFT で色付けしたガウス雑音をチャンクごとに生成する（多重レート）。乱数の seed は記録する（prompt §6）。
- f⁻² の雑音は f_lo で必ず切る（積分が発散する）。σ² は conventions §6 に従い「測定時間の逆数から帯域の上限まで」で評価する。
- ロックの持続時間の試験（§9 V7）は数分〜数十分の模擬時間が必要。

### 7.4 他のトピックへ渡すもの

- T03: κ_τ、τ_g,s(t)、ストレッチャ補正後に残る群遅延（§2.3）。
- T04: ψ の PSD と σ²（指標 (i)(ii) への写像は T04）。
- T09: k_T,slow − k_T,fast（2.3）の値を返してもらう。
- T10: （2.10）（2.11）、空間構造のパラメータ、軸ごとの係数。
- T16: R_act の値（915B: 64.9 rad）、熱アクチュエータの係数、ドリフトの速さ k_T L r_T、ξ_eo の扱い（§10 の 2）。

## 8. 入力・出力の一覧

- 入力: λ、ファイバの種類、区間の表（L_s、d_s、軸、ゾーン、τ_th,s、N_sub）、T0、r_T、b_T、h_ref、L_ref、γ、f_c、f_lo、ピークの表、ρ_ij(f) のパラメータ（ρ_0、f_ρ）、熱力学の物性値、seed。
- 出力: ψ_s,a(t)、τ_g,s(t)、各成分の PSD（解析式と推定値）、区間間の相関（推定値）。

## 9. 検証テスト（実装が通すべきもの）

- **V1 熱係数と符号**: L = 1 m、一様に ΔT = +1 K、λ = 1548.8 nm、n = 1.468、α_L = 5.5×10⁻⁷ /K。dn/dT = 9.2×10⁻⁶ /K で Δψ = +40.60 rad、1.0×10⁻⁵ /K で +43.84 rad（相対誤差 10⁻⁶ 以下）。符号は正（conventions §4）。
- **V2 群遅延**: L = 63 m、ΔT = 1 K、κ_τ = 39 ps/(km·K) → Δτ_g = 2.457 ps。
- **V3 機械雑音の PSD**: h = 10⁻³ rad²·Hz、f_lo = 0.1 Hz、f_c = 1 kHz で生成した系列の Welch 推定が、1–300 Hz の各 1/3 オクターブ帯で（2.7）と ±10% 以内（十分な平均回数で）。分散が ∫ S df と ±5% 以内。
- **V4 長さのスケーリング**: 長さ 10 倍で、γ = 1（小区間が独立）なら PSD は 10 倍（±10%）、γ = 2（lumped）なら 100 倍。
- **V5 相関**: S_i = S_j、ρ = 0.9（周波数によらない）で、差の PSD = 0.2 S（±10%）、コヒーレンスの推定 ≈ 0.81。ρ = 1 で差が 0（数値の床まで）。`env.correlation` をオフにすると差は 2S。
- **V6 逆向き伝搬**: N_sub ≥ 100 の独立な小区間（distributed、白色の擾乱）で、同じ時刻に両端から出る順方向・逆方向の位相の差の PSD が 2S[1 − sinc(2πfτ)] に一致（2πfτ = 0.3、1、π、2π で ±10%。2πfτ = π では S_D = 2S）。lumped では S_D = 0（数値の床まで）。
- **V7 レンジの使い切り**: lumped の区間 L = 63 m、τ_th = 10 s、r_T = 10⁻⁴ K/s、k_T = 44 rad/(m·K)（このテストでは 43.84 を使う）、レンジ 64.9 rad の積分器のアクチュエータで、定常のドリフトになってから飽和までの時間 = 64.9/(43.84 × 63 × 10⁻⁴) = 235 s（±1%、τ_th の遅れを補正したうえで）。
- **V8 Ref C の目安**: Ref C のプリセットで、dither 以外の残留位相雑音（環境の残差）が 0.083 rad 以下。dither を足した全体が約 0.1 rad（prompt §7-10）。これは合わせ込みではなく、パラメータの範囲が報告と矛盾しないかの確認。
- **V9 熱力学**: Dong の Table I と L = 1308 m で、S^TM(1 Hz) = 4.36×10⁻¹⁰ rad²/Hz、S^TC(1 Hz) = 2.86×10⁻¹² rad²/Hz（±1%）。L を 327 → 134、68、34 m（Michelson の実効長はその 4 倍）にすると −3.9、−6.8、−9.8 dB（Dong 2016 の記述）。
- **V10 1 次の遅れ**: 周囲温度のステップ ΔT = 1 K に対し、ΔT_s(t) = 1 − e^{−t/τ_th}（±10⁻³）。
- **V11 トグル**: すべての `env.*` をオフにすると ψ_s,a(t) = ψ̄_s,a（厳密）、τ_g,s 一定。
- **V12 軸の係数**: 同じ ΔT_s で、slow と fast の位相の差が (k_T,slow − k_T,fast) L ΔT（T09 の値を入れて厳密に一致）。ε_s による差は k_ε,a の差だけ。

## 10. 規約・DECISIONS・サーベイとの関係（調整役への提案）

1. **サーベイ E22 の 44 rad/(m·K)** は文献の範囲（40–48）で確認できた。E26 の 39 ps/(km·K) も Slavík 2015 で確認できた。サーベイの「UNVERIFIED」の印を外し、出典を Dong 2016・Slavík 2015 にしてよい（サーベイの編集は調整役）。
2. **conventions §4 の ξ_eo と 915B**: 915B の資料の「16 µm」は「optical path length measurement」と書かれている（光路長として測った値）。この場合は光弾性の補正をすでに含んでいて、δψ = 2π δOPL/λ（ξ_eo をもう一度かけない）。T16 で確認し、conventions §4 に注記を足すことを提案する。
3. **conventions §7 の片側 PSD**: Williams 2008 は PSD が片側か両側かを明記していない（読んだ範囲で）。ここでは片側として使った（UNVERIFIED）。h_ref は 3 桁掃引するので結論には効かない。
4. **D11 の 2.4 rad/(m·K)**: T08 では検証していない（T09 の担当）。（2.3）の形と仮の値だけを使った。
5. **サーベイ E24 の「Wanser 型」**: 低周波では Wanser の熱伝導雑音より Duan の熱機械雑音（1/f）が支配的（約 150 Hz より下）。E24 の名前を「熱力学的な基底（熱伝導＋熱機械）」にすることを提案する。どちらも無視できる大きさ。
6. **サーベイ E25 の「C は定性的な記述のみ」**: 確認した。Ref C に数値はない。
7. DECISIONS との食い違いは見つからなかった。

## 11. ユーザの決定が必要な点（推測で埋めていない）

1. **温度の環境**: r_T の既定値（ここでは 10⁻⁴ K/s [AE]）と範囲、空調の周期成分を入れるか。ポートのファイバを箱に入れる・温度制御するか（τ_th と b_T が大きく変わる）。D27 によりユーザの実験室の測定値は Phase 5 の終わりまで使わない。
2. **ポートのファイバの置き方**（D13 の相関に直結）: 32 本を 1 つの spool・束にまとめるか、別々の spool か。squeezed と LO_k を同じ束に入れるか（ρ_0 の既定値が 0.9 か 0 かで結果が変わる）。
3. **遅いアクチュエータ**: TEC（Liu 型）を候補に入れるか（D14・D30 の比較に加える）。入れるなら L_TEC と ΔT_TEC の範囲。
4. **Ref D の 915B のファイバ長**: 資料に書かれていない。メーカーへの確認か、ユーザが知っている値があれば使う（D02 は「メーカー資料の値だけ」）。

## 12. 文献

本文を読んで値を取ったもの:

- Ref C: S. Verclas, C. Darsow-Fromm, B. Tohermes, R. Schnabel, "Fiber distribution of phase-stabilized GHz-bandwidth squeezed vacuum states of light between two buildings," Optica Quantum 4, 1 (2026). DOI 10.1364/OPTICAQ.567418.（`refs/`）
- H.-L. Liu et al., "Gaussian boson sampling with 1,024 squeezed states in 8,176 modes," Nature 653, 687 (2026). DOI 10.1038/s41586-026-10523-6.（`refs/`。開始ページは PDF の柱 "Nature | Vol 653 | 21 May 2026 | 687" から読んだ）
- Evanescent Optics, Model 914-2 Piezo Controller / 915B Piezo Assembly の仕様書（serial 1527、8669-1、8669-2）。`refs/manuals/914-915B-Specifications.pdf`（スキャン画像を目で読んだ）。
- P. A. Williams, W. C. Swann, N. R. Newbury, "High-stability transfer of an optical frequency over long fiber-optic links," J. Opt. Soc. Am. B 25, 1284–1293 (2008). DOI 10.1364/JOSAB.25.001284（**DOI は書誌から組み立てたもので UNVERIFIED**、巻・号・開始ページは PDF の柱で確認）。PDF: https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=32873（2026-10-28 閲覧）。
- R. Slavík, G. Marra, E. Numkam Fokoua, N. Baddela, N. V. Wheeler, M. Petrovich, F. Poletti, D. J. Richardson, "Ultralow thermal sensitivity of phase and propagation delay in hollow core optical fibres," Sci. Rep. 5, 15447 (2015). DOI 10.1038/srep15447.
- J. Dong, J. Huang, T. Li, L. Liu, "Observation of fundamental thermal noise in optical fibers down to infrasonic frequencies," arXiv:1512.02799 (2015)。arXiv のページにある関連 DOI は 10.1063/1.4939918（Appl. Phys. Lett. と思われるが**掲載誌・巻・ページは UNVERIFIED**）。値は arXiv 版から取った。
- C. Hilweg, D. Shadmany, P. Walther, N. Mavalvala, V. Sudhir, "Limits and prospects for long-baseline optical fiber interferometry," Optica 9, 1238–1252 (2022). DOI 10.1364/OPTICA.470430. 値は arXiv:2208.09247 から取った。
- S. Koke, A. Kuhl, T. Waterholter, S. M. F. Raupach, O. Lopez, E. Cantin, N. Quintin, A. Amy-Klein, P.-E. Pottie, G. Grosche, "Combining fiber Brillouin amplification with a repeater laser station for fiber-based optical frequency dissemination over 1400 km," arXiv:1911.01215 (2019)。**掲載誌は UNVERIFIED**（プレプリントの体裁は IOP）。

書誌だけを確認し、本文は読んでいないもの:

- L.-S. Ma, P. Jungner, J. Ye, J. L. Hall, "Delivering the same optical frequency at two places: accurate cancellation of phase noise introduced by an optical fiber or other time-varying path," Opt. Lett. 19, 1777–1779 (1994). DOI 10.1364/OL.19.001777. ファイバ雑音補償の原典。数値は使っていない（本文を入手できなかった）。
- K. H. Wanser, "Fundamental phase noise limit in optical fibres due to temperature fluctuations," Electron. Lett. 28, 53–54 (1992). DOI 10.1049/el:19920033. 式は Hilweg 2022・Dong 2016 の引用から取った。
- L. Z. Duan, "General treatment of the thermal noises in optical fibers," Phys. Rev. A 86, 023817 (2012)（タイトルは UNVERIFIED、書誌は Dong 2016 の文献リストから）。L. Z. Duan, Electron. Lett. 46, 1515 (2010)（同）。
- S. Foster, A. Tikhomirov, M. Milnes, IEEE J. Quantum Electron. 43, 378 (2007)（Dong 2016 の文献リストから、UNVERIFIED）。
- N. Lagakos, J. A. Bucaro, J. Jarzynski, "Temperature-induced optical phase shifts in fibers," Appl. Opt. 20, 2305–2308 (1981). DOI 10.1364/AO.20.002305（著者の Jarzynski は UNVERIFIED）。被覆による熱位相係数の増加の文献として挙げるだけで、値は使っていない。

UNVERIFIED のまま使った値:

- 915B の巻線部のファイバ長 約 9 m（Photonics.com の製品掲載の検索結果の抜粋。ページは開けなかった）。
- 温度の時間モデルの既定値（r_T、τ_th、b_T）、ρ(f) の形、f_c より上の f⁻⁴、γ = 2 の正規化: すべて [AE]。
