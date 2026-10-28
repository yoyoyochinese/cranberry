# T16 アクチュエータ・制御器・デジタル部（actuators-control-digital）

- Phase: 1（deep dive）、トピック T16（`notes/survey.md` §4）
- 日付: 2026-10-28
- 担当する効果: E27 の機械的な部分、E77–E84、機器プリセット（Thorlabs DSC1、Moku:Delta、Moku:Pro の要否）、prompt §5.3 の抽象制御器インターフェース、§5.4 の時間ループ実行方式、§7-7・§7-10 の検証。
- 入力: `docs/squeezing_sim_prompt.md`、`DECISIONS.md`（D01–D33）、`notes/conventions.md`、`notes/survey.md`、
  `refs/manuals/914-915B-Specifications.pdf`、`refs/manuals/thorlabs-dsc1.pdf`、`refs/manuals/MokuDelta-Specifications.pdf`、Ref A・B・C、Liu et al. 2026。
- 値の印: **[915B doc p.N]** などは `refs/manuals/` の資料のページ。**[AE]** は本ノート担当の概算（未検証）。**UNVERIFIED** は原典で確認できなかった値・文献。
- `blind/` は読んでいない（D04）。Ref D の値はメーカー資料と物理だけから作った。

---

## 0. 要点（先に結論）

1. **915B + 914 の位相係数**: 資料の「full range の駆動電圧で 16 µm（S/N 8669-1）/ 17 µm（8669-2）、optical path length measurement」[915B doc p.1] と「入力 ±5 V で full swing」[914 doc p.1] から、
   入力電圧あたり **K_s = 2π ΔOPL_full /(λ0 ΔV_in,full) = 6.49 rad/V（8669-1）/ 6.90 rad/V（8669-2）**（λ0 = 1548.8 nm、ΔOPL を光路長の変化と読んだ場合）。
   16 µm を物理的な伸びと読むと ξ_eo n_eff ≈ 1.16 倍の 7.5 / 8.0 rad/V。**どちらの読みかは資料から決められない**（§11 の質問 1）。
2. **Ref D のレンジと量子化は DSC1 で決まる**: DSC1 の出力は ±4 V [DSC1 p.2] なので 914 の入力 ±5 V の 80% しか使えず、位相のレンジは **±26.0 rad（全幅 51.9 rad、約 8.3 縞）**。
   出力 12 bit の 1 LSB = 1.95 mV は **12.7 mrad**。白色量子化雑音として 914 の帯域までで約 1.5 mrad rms、limit cycle なら ±6 mrad 程度 [AE]。Ref A の θ = 9 mrad と同じ桁なので、量子化は既定でオンにする。
   ADC 側（16 bit、122 µV）は K_D ≈ 3.3 V/rad で 37 µrad と小さい。
3. **ループ遅延**: DSC1 は「入力→出力の −180° が typ. 60 kHz（>58 kHz）」[DSC1 p.1–2]。入力帯域 120 kHz・出力帯域 100 kHz を 1 次の低域通過とみなすと、残りの純遅延は **τ ≈ 5.7 µs**（全部を遅延とみなした上限は 8.3 µs）。
   Moku:Delta の PID は「400 kHz で位相遅れ 30°」[Moku p.65] から **τ ≲ 0.21 µs**。ポート光路の光の遅延（63 m で 0.31 µs、630 m で 3.1 µs）も同じ桁になりうるので、遅延バッファに含める。
4. **将来の構成 F では、squeezed 光路に置けるのは PZT ストレッチャ（＋熱アクチュエータ）だけ**: LiNbO3 ファイバ EOM の挿入損失は 2.5–3.5 dB（販売店の表示、UNVERIFIED）で、squeezed 光路に入れると squeezing の上限が 約 2.6–3.6 dB になり、D16 の目標 3–10 dB と両立しない。
   一方、ポート長 L_31 = 63 m の熱ドリフトは 約 2.8×10³ rad/K [T08 の値待ち、survey の AE] で、915B 1 台のレンジ（65 rad）では 約 0.02 K しか追えない。TEC（Liu et al.）、スプール全体の温度制御、2π の整数倍のリセット（fringe hop）のどれかが要る。
5. **Moku:Pro のプリセットは不要**（D25 への回答）: Ref C の §7-10 の検証項目（dither 0.0555 rad rms、全体 約 0.1 rad、遅いアクチュエータのレンジでロック時間が決まること）は、論文の値と汎用の FPGA 制御器の値（遅延 0.1–2 µs、16 bit）で再現でき、kHz 帯のループでは遅延の差が数度以下しか効かない。
   足りないのは PID の利得・各アクチュエータのレンジ・雑音 PSD で、これらは Moku:Pro の仕様書にも載っていない（利用者の設定値）。
6. **§5.4 の推奨**: 閉ループは **(c) コンパイルした 1 サンプルごとのループ**（numba など）、f_sim ≈ 0.5–10 MHz。それより速い RF 搬送波（ESA の 3 MHz、Ref C の 8–300 MHz の雑音電力検出、必要なら Ref A の 1 MHz ヘテロダイン）は **(a) baseband/統計モデル**。
   (b) のチャンク方式は「チャンク長 ≤ ループ遅延」が条件で、Moku 系の 0.2 µs では 1 サンプルになり、数分のロック時間（Ref C、Ref D）を Python のオーバーヘッドで回せない。開ループの部分（雑音生成、ESA の後処理）だけに使う。

---

## 1. 記号（symbols, English）

共通:
- `λ0`: vacuum carrier wavelength [m]. `k0 = 2π/λ0`. `β = n_eff k0`: propagation constant [rad/m].
- `ψ_seg,a`: phase of path segment `seg` on polarization axis `a` [rad] (conventions §9). `δψ`: actuator-induced change of `ψ_seg,a`.
- `u(t)`: actuator drive voltage at the actuator-chain input (controller DAC output after the load divider) [V].
- `s`: Laplace variable [rad/s]; `f`: Fourier frequency [Hz]; `ω = 2πf`.
- `σ_act ∈ {+1, −1}`: polarity (wiring) sign of an actuator; calibration knob.
- `K_D`: error-signal slope at the lock point [V/rad] (from T14/T15).
- `L(s)`: open-loop gain; `S(s) = 1/(1 + L)`: sensitivity; `T(s) = L/(1 + L)`: complementary sensitivity.
- `f_u`: unity-gain frequency of `L` [Hz]. `PM`: phase margin [deg].

PZT ファイバストレッチャ:
- `δL`: physical elongation of the fiber [m]. `ΔOPL`: change of optical path length [m], `ΔOPL = ξ_eo n_eff δL`.
- `ξ_eo = 1 − (n²/2)[p12 − ν(p11 + p12)]`: elasto-optic (strain-optic) correction factor [–].
- `p11`, `p12`: Pockels strain-optic coefficients of silica [–]. `ν`: Poisson ratio of the fiber [–]. `n`: core index [–].
- `ΔOPL_full`: ΔOPL for the full driver swing [m]. `ΔV_in,full`: full input swing of the driver (914: 10 V = ±5 V) [V].
- `K_s = 2π ΔOPL_full / (λ0 ΔV_in,full)`: DC stretcher gain referred to the driver input [rad/V].
- `H_drv(s)`: small-signal driver (914) transfer function, normalized to 1 at DC. `f_d1`, `f_d2`: its poles [Hz]. `SR`: driver slew-rate limit [V/s].
- `H_mech(s) = Σ_i r_i ω_i² / (s² + (ω_i/Q_i) s + ω_i²)`: modal mechanical response, `Σ_i r_i = 1`; `ω_i = 2π f_i` resonance, `Q_i` quality factor, `r_i` modal residue (may be negative → anti-resonance).
- `N_w`: number of fiber wraps (38). `ℓ_w`: fiber length per wrap [m]. `L_w = N_w ℓ_w`: wound fiber length [m].
- `h_hys`: hysteresis width as a fraction of the actuator range [–]. `γ_cr`: creep coefficient per decade [–].
- `C_p`: PZT capacitance [F]. `R_o`: driver output resistance [Ω]. `R_in,914`: 914 input resistance [Ω].

EOM:
- `V_π`: half-wave voltage [V]; `V_π,LF` at low (servo) frequencies, `V_π,RF` at RF. `δψ = π u / V_π`.
- `IL_EOM`: insertion loss [dB]; `η_EOM = 10^(−IL_EOM/10)` [–].

PZT ミラー:
- `δx`: mirror displacement along its normal [m]. `θ_i`: angle of incidence [rad]. `d_x`: displacement per volt [m/V].

TEC / 熱アクチュエータ:
- `ΔT_h`: temperature change of the heated/cooled fiber section [K]. `L_h`: length of that section [m].
- `c_T = (1/L) dψ/dT = k0 (dn/dT + n α_L)`: thermal phase coefficient per length [rad/(m·K)] (value from T08).
- `K_T`: DC thermal gain [K/V] (or [K/A]). `τ_th`: thermal time constant [s]. `τ_tr`: thermal transport (diffusion) delay [s].

AOM:
- `Ω_RF = 2π f_RF`: acoustic (RF) angular frequency. `φ_RF`: RF drive phase [rad]. `m ∈ {+1, −1}`: diffraction order. `τ_a`: acoustic transit time across the beam [s].

デジタル部:
- `N_ADC`, `N_DAC`: bits. `V_FS`: full-scale span (e.g. 8 V for ±4 V) [V]. `q = V_FS / 2^N`: LSB [V].
- `f_s = 1/h`: controller update rate [Hz] (`h`: update period, DSC1 notation). `f_sim = 1/Δt`: simulation rate.
- `τ_L`: total pure loop delay [s] (`τ_ADC + τ_proc + τ_DAC + τ_opt + τ_cable`). `τ_opt = n_g L_opt / c`: optical propagation delay from actuator to detector.
- `S_q = q² / (6 f_s)`: one-sided PSD of white quantization noise [V²/Hz].
- `K_p` [–], `K_i` [1/s], `K_d` [s]: canonical PID coefficients. `f_I`, `f_D`: integrator / differentiator crossover frequencies (Moku convention).
- `e(t)`: error signal. In this project `e = (linear combination of beam phases) − setpoint` (conventions §9), mapped to volts by `K_D`.
- `R_k`: one-sided phase range of actuator k [rad]. `r_drift`: drift rate [rad/s]. `D_ψ`: phase diffusion coefficient, `Var = 2 D_ψ t` [rad²/s].

---

## 2. 規約との対応（conventions §4, §9）

- **アクチュエータの符号**: 「駆動電圧 u を上げると区間の位相 ψ が増える」を正（σ_act = +1）と定義する。ストレッチャで伸ばす・温度を上げる・ミラーを遠ざける（光路が長くなる）は conventions §4 により ψ を増やす。
  実機の配線の極性は資料にないので、σ_act を較正用の値として残す。
- **誤差の符号**: conventions §9 は e = (線形結合) − 設定値。DSC1 は **e_DSC1 = S − v**（S: 設定電圧、v: 入力電圧）[DSC1 p.18]。
  入力電圧が v = V_0 + K_D·(線形結合 − φ_set) なら e_DSC1 = −K_D·e_conv。プリセットはこの符号反転を明示的に持つ（P・I・D の符号は DSC1 では連動して反転する [DSC1 p.19]）。
- **ストレッチャの位相**: conventions §4 の δψ = ξ_eo β δL をそのまま使う。資料が ΔOPL を与える場合は δψ = k0 ΔOPL（ξ_eo と n_eff は ΔOPL に含まれている）。§3.1.2 で変換を示す。
- **PSD は片側**（conventions §7）。量子化雑音の S_q = q²/(6 f_s) も片側。

---

## 3. アクチュエータ

### 3.1 PZT ファイバストレッチャ（Evanescent Optics 915B ＋ Model 914-2）

#### 3.1.1 資料にある値（`refs/manuals/914-915B-Specifications.pdf`、スキャン 3 ページ、販売代理店ハナムラオプティクス経由）

| 項目 | 値 | 出典 |
|---|---|---|
| ドライバ | Model 914-2 Piezo Controller、S/N 1527、2 チャンネル | [914 doc p.1] |
| 入力 | BNC、**±5 V で full piezo voltage swing**、シェルは接地 | [914 doc p.1] |
| 出力 | 915B へ telco ケーブル、DB-9 のモニタ: ピン 2・3（Ch1）、7・8（Ch2）が push-pull 出力 | [914 doc p.1, p.3] |
| 周波数応答 | Fig. 1: ±5 V 入力でのピエゾ電圧。DC–約 7 kHz で 52 V、以下 §3.1.4 の読み取り | [914 doc p.2] |
| 915B S/N 8669-1 | 38 巻、PM15-U25D、**16 µm**（full range、optical path length measurement）、FC/APC（slow 軸） | [914 doc p.1] |
| 915B S/N 8669-2 | 38 巻、PM15-U25D、**17 µm**、FC/APC（slow 軸） | [914 doc p.1] |
| 機械共振 | **記載なし** | — |
| PZT 容量、ドライバの出力抵抗・入力抵抗・雑音・小信号帯域 | **記載なし** | — |
| ヒステリシス・クリープ | **記載なし** | — |
| 巻き 1 周のファイバ長 | 220 mm、「9 m のファイバで最大 16 µm」、「DC–10 kHz で平坦、100 kHz で 0」、楕円形の型に PZT 板を弾性的に付け、直線部だけを伸ばす（PM の曲げによるクロストークを避ける） | photonics.com の製品掲載（検索結果の抜粋のみ、ページ本体は HTTP 403 で取得できず、2026-10-28）。**UNVERIFIED** |

注: D02 は「915B の rad/V・共振・レンジはメーカー資料の値だけを使う」。rad/V とレンジは上の資料から出せるが、**共振は資料にない**（§11 の質問 2）。

#### 3.1.2 支配方程式

伸び δL に対する位相（conventions §4）:

  δψ = ξ_eo β δL = k0 ξ_eo n_eff δL = k0 ΔOPL.

光弾性の補正（Butter & Hocker 1978 の式の形）:

  ξ_eo = 1 − (n²/2) [p12 − ν (p11 + p12)].

**ξ_eo の確認（conventions §4 の「≈ 0.78、T16 で確認」）**:
- バルクの石英でよく使われる p11 = 0.121、p12 = 0.270、ν = 0.17、n = 1.458（この組は原典を確認していない、UNVERIFIED）→ **ξ_eo = 0.784**。
- Bertholds & Dändliker 1988 がファイバで測った p11 = 0.113、p12 = 0.252（検索で書誌と値を確認）と ν = 0.16–0.17、n = 1.444–1.468 → **ξ_eo = 0.795–0.798**。
- よって **ξ_eo = 0.78–0.80**。conventions の 0.78 は 2% 以内で妥当。既定値 0.78、掃引範囲 0.78–0.80 を提案。
- PM ファイバ（PM15-U25D、応力付与部あり）では軸ごとに ξ_eo がわずかに違い、伸ばすと 2 軸の位相差も変わる。軸ごとの係数は T09 の担当（D11 の「係数は軸ごと」に入れる）。

915B では資料が ΔOPL_full を与えるので、ドライバ入力あたりの直流の係数は

  K_s = 2π ΔOPL_full / (λ0 ΔV_in,full).

動特性まで含めたブロック:

  δψ(s) = σ_act K_s H_drv(s) H_mech(s) u(s),  |u| ≤ 5 V（914 の入力）、さらに制御器の出力レンジでクリップ。

- H_drv: 914 の小信号応答（§3.1.4）。大振幅ではスルーレート制限 |dV_p/dt| ≤ SR を非線形のクリップとして別に持つ（オン・オフ可）。
- H_mech: モード和の形にする（反共振（零点）を表せるので、極の積より実機に近い）。既定は 1 モード（r_1 = 1）。
- ヒステリシス（オプション、既定オフ）: 最も簡単な play operator（バックラッシュ）
  y[n] = min(u[n] + w/2, max(u[n] − w/2, y[n−1])), w = h_hys × (u のレンジ)。
  h_hys は資料にない。開ループの積層 PZT で一般に 10–15% と言われる値（UNVERIFIED、メーカー一般の記述）を掃引の目安とする。
  ループの中では積分器が補正するので、効くのは「局所的な傾き K_s の変化（ループ利得の変化）」と、ランプでの再引き込みの位置ずれ。
- クリープ（オプション、既定オフ）: δψ(t) = δψ_0 [1 + γ_cr log10(t/t_0)]、γ_cr ≈ 0.01–0.02（UNVERIFIED）。Hz 以下の遅い変化で、積分器で消える。

#### 3.1.3 静的な値（Ref D）

λ0 = 1548.8 nm（D02）で計算 [AE]:

| 量 | 8669-1（16 µm） | 8669-2（17 µm） | 備考 |
|---|---|---|---|
| 全レンジ（914 ±5 V） | 64.9 rad（10.3 縞） | 69.0 rad | ΔOPL と読んだ場合 |
| K_s | **6.49 rad/V** | 6.90 rad/V | 914 の入力電圧あたり |
| DSC1 ±4 V でのレンジ | **±26.0 rad**（全幅 51.9 rad） | ±27.6 rad | 914 入力の 80% |
| 2π に相当する入力電圧 u_2π = 2π/K_s | 0.968 V | 0.911 V | fringe hop の幅 |
| 物理的な伸びと読んだ場合の K_s | 7.53 rad/V | 8.00 rad/V | × ξ_eo n_eff ≈ 0.79 × 1.468 |

- 巻いたファイバの長さ: L_w = 38 × 0.22 m ≈ **8.4 m**（photonics.com の「9 m」と整合、UNVERIFIED）。D02 の「ΔL の範囲は 915B の巻線ファイバ長から見積もる」に使える: **ΔL ≳ 8 m（＋ピッグテールの差）**。T05・T15 に渡す。
- ひずみ: 16 µm / (0.784 × 1.47 × 8.4 m) ≈ 1.7×10⁻⁶。弾性の範囲で十分小さい。
- DSC1 の出力抵抗 220 Ω と 914 の入力抵抗 R_in,914（資料にない）で分圧: u = u_DAC R_in,914/(R_in,914 + 220 Ω)。R_in,914 ≥ 10 kΩ なら 2% 以下（DSC1 は >100 kΩ の負荷を想定 [DSC1 p.2 脚注 b]）。較正用の値として残す。

#### 3.1.4 914 の周波数応答（Fig. 1 の読み取り）

[914 doc p.2] の Fig. 1 を 300 dpi で読み取った値（±5 V 入力時のピエゾ電圧、縦軸の 52 V が振幅か p-p かは書かれていない）:

| f [kHz] | 7.06 | 10.0 | 19.7 | 29.9 | 50.0 | 100 |
|---|---|---|---|---|---|---|
| V_p [V] | 51.7 | 50.1 | 42.3 | 36.2 | 23.7 | 11.0 |
| V_p/52.1 | 0.992 | 0.962 | 0.812 | 0.695 | 0.455 | 0.211 |

- −3 dB は 約 28 kHz。
- **線形モデル（最小位相を仮定）**: 2 つの実の極 f_d1 = 29.5 kHz、f_d2 = 112 kHz で全点を 2% 以内で再現 [AE、対数誤差の最小二乗]。1 極（24.4 kHz）では 10–30 kHz で 5–10% ずれる。
  この場合の位相遅れ: 1 kHz で 2.5°、5 kHz で 12°、10 kHz で 24°。
- **スルーレートとしての読み**: 30–100 kHz の裾は V_p ∝ 1/f で、2π f V_p = 6.8–7.4 V/µs（52 V を振幅と読んだ場合。p-p なら半分）。
  位相に直すと、どちらの読みでも最大の位相変化速度は 約 **4.3 rad/µs**（100 kHz 点: 64.9/2 × 0.211 × 2π × 10⁵）。ロックで要る速度より何桁も大きく、効くのは大振幅のランプや fringe hop の過渡だけ。
- Fig. 1 は大振幅（full swing）の測定なので、小信号の帯域はこれ以上ある。**既定は保守側の線形 2 極モデル**、スルーレート制限はオプション。
- PZT 容量 C_p とドライバの出力抵抗 R_o が分かれば f_RC = 1/(2π R_o C_p) で置き換える。資料にない（§11 の質問 2）。
- 出力雑音: 資料にない。入力換算の雑音電圧密度 e_n,914 [V/√Hz] をパラメータにし、位相雑音 S_φ = (K_s |H_drv H_mech|)² e_n,914²。

#### 3.1.5 機械共振（資料にない）

- 資料（914/915B）には共振の記載がない。販売店の掲載の「DC–10 kHz 平坦、100 kHz で 0」は 914 の Fig. 1 と同じ形で、ドライバの応答の可能性が高い（UNVERIFIED）。
- 物理の目安 [AE]: 円筒 PZT のフープ（呼吸）モードは f ≈ c_PZT/(2π r)、c_PZT = √(Y/ρ) ≈ √(6×10¹⁰/7.75×10³) ≈ 2.8 km/s。r = 10–20 mm で 22–45 kHz。
  915B は楕円の型に PZT 板を貼る構造（UNVERIFIED）で、型・ファイバ・取り付けの複合系の共振はこれより低いことが多い（Ref B の PZT ミラーでも、単体より複合系の固有振動数で決まると述べている [Ref B §3]）。
- **シミュレータの扱い**: f_m と Q を掃引するパラメータとし、既定値は置かない（D02 に従う）。掃引範囲の提案: f_m = 5–50 kHz、Q = 3–100。
  ループの帯域の目安は、ノッチなしの PI で f_u ≲ f_m/5–f_m/10（Q が大きいほど低く）。

### 3.2 ファイバ EOM（位相変調器）

- 支配方程式: δψ = σ_act π u / V_π。帯域は GHz（LiNbO3 進行波型）なので、時間ループでは静的な利得として扱い、帯域はドライバ・制御器側で決まる。
- レンジ: ±π V_max/V_π。V_max はドライバの出力（DSC1 ±4 V、Moku:Delta PID ±0.5 V into 50 Ω [Moku p.65]）か EOM の許容入力の小さい方。
  V_π = 4 V なら DSC1 で ±π、Moku:Delta PID で ±0.39 rad（50 Ω 負荷）。**ドリフトを追うレンジはない**ので、Ref C のように速い段だけを受け持つ。
- 低周波（servo 帯）の V_π,LF と RF の V_π,RF は一般に異なる（LiNbO3 の圧電・焦電の寄与）。パラメータを分ける。値は資料がないので UNVERIFIED。
- 値:
  - Ref C: exail MPZ-LN-10、**U_π = 4 V**、200 mVpp → |±4.5°| [Ref C §3]。検算: π × 0.1 V / 4 V = 0.0785 rad = 4.50°、RMS 0.0555 rad（prompt §7-10 の 約 0.056 rad と一致）。
  - Ref A: iXblue MPX-LN-0.1 を LO と probe の両方に使用 [Ref A §3]。V_π は本文に記載なし（UNVERIFIED、必要なら製品資料で確認）。
  - 挿入損失: MPZ-LN-10 で 2.5–3.5 dB（販売店 gophotonics の掲載、検索結果の抜粋、2026-10-28、UNVERIFIED）。exail の MPX/MPZ シリーズ資料の 2.5 dB という表示も検索で見たが、型番との対応は確認できなかった。
- **配置の規則（prompt §5.5「squeezed 光路の変調器は挿入損失を squeezed 状態にかける」）**: η_EOM = 10^(−0.25…−0.35) = 0.45–0.56。
  無限の squeezing でも V_out ≥ (1 − η) V_vac なので、上限は −10 log10(1 − η) = **2.6–3.6 dB**。D16 の目標 3–10 dB と両立しないので、F では EOM は pump・LO・probe・参照光の光路だけに置く（D10 の「スロット同期 EOM でポンプ位相を補正」は pump 側なので可）。
- 残留振幅変調（RAM）は誤差信号のオフセットになる。T14 の担当（このブロックは RAM 係数の受け口だけ持つ）。

### 3.3 PZT ミラー（Ref B、Ref C）

- 支配方程式: 反射で光路が往復するので δψ = σ_act 2 k0 cos θ_i δx、δx = d_x H_drv(s) H_mech(s) u。
- Ref B: 単層 PZT から積層 PZT（Thorlabs AE0203D04）に替え、**PZT1 の共振 14 kHz、PZT3 の共振 22 kHz**（単層では 約 2 kHz）。共振は PZT・ミラー・マウントの複合系の固有振動数で決まり、2 つの差は個体差と著者は書いている [Ref B §3]。Q は記載なし。
  PZT2 は 295 kHz の変調に使い、そこが PZT2 の共振でもある [Ref B §3]。d_x（AE0203D04 の変位/電圧）は本文にない（UNVERIFIED、Thorlabs の製品ページで確認できるが、Ref B の再現には共振周波数とループ利得があれば足りる）。
- Ref C: 「piezo-mounted mirror」＋ HV アンプ、ランプ発生器を加算 [Ref C §2、Fig. 3]。値の記載なし → Ref B と同じ形で f_m = 10–20 kHz を掃引の目安（UNVERIFIED）。
- Liu et al.: 各空間モードの光路に PZT ミラーを置き、逆向き probe の誤差信号で駆動。精度 1 nm（光源）、3 nm（1 段目の遅延線）、10 nm（2 段目）を数十時間維持 [Liu et al. Methods]。
  2 段目の 10 nm は 2π × 10/1550 = 0.041 rad（光路長の誤差として。位相 = k0 × 光路長）。

### 3.4 TEC・熱アクチュエータ（Liu et al.、Ref C の提案）

- Liu et al.: 長いファイバの温度ドリフトが 1 つの PZT のレンジを超えることが多いので、TEC のフィードバックで遅い熱ドリフトを抑える [Liu et al. Methods]。数値（利得、帯域、長さ）の記載はない。
- Ref C: 「2 本のファイバの一方の数 m を温度で動かせば、ロック時間を任意に延ばせる」と提案 [Ref C §4]。
- モデル（集中定数の 1 次＋輸送遅延）:
  ΔT_h(s) = K_T e^{−s τ_tr} / (1 + s τ_th) · u(s)、 δψ = σ_act c_T L_h ΔT_h、 |ΔT_h| ≤ ΔT_max。
  TEC のペルチェ・ジュール項による電流と熱量の非線形性は、動作点のまわりで線形化して K_T に入れる（必要なら 2 次の項を追加）。
- c_T は T08 の値を使う。T08 が出るまで survey の [AE] 値 **c_T ≈ 44 rad/(m·K)** を仮に使う。
- 概算 [AE]: L_h = 3 m、ΔT_max = ±10 K → レンジ ±1.3×10³ rad。τ_th はファイバの被覆の熱拡散だけなら (125 µm)²/(10⁻⁷ m²/s) ≈ 0.16 s だが、実際は巻き枠・板の熱容量で決まり 10–100 s 程度（設計次第、UNVERIFIED）。帯域は mHz 帯。
- 役割: 遅い段（Ref C の 3 段目、または F の PZT の下の段）。PZT との受け渡し周波数は mHz。

### 3.5 AOM（周波数・位相のアクチュエータ）

- conventions §4 では E ∝ A e^{−iω0 t}。周波数を Δω 上げる（ω0 → ω0 + Δω）ことは A ∝ e^{−iΔω t}、すなわち **φ(t) = φ(0) − ∫Δω dt**。
- RF 駆動 ∝ cos(Ω_RF t + φ_RF) のとき、m 次の回折光の複素振幅は A_m = √η_AO A e^{−i m (Ω_RF t + φ_RF)}（この規約での導出。物理的な符号は機器の向きで変わるので σ_act を残す）。
  よって**位相アクチュエータとして δψ = −m φ_RF**、**周波数アクチュエータとして dψ/dt = −m 2π δf_RF**。
- 特長: RF 位相は 2π で巻き戻せるので**レンジの制限がない**。周波数で駆動すれば積分器（1/s）の伝達関数になる。帯域は音響波がビームを横切る時間 τ_a ≈ w/v_a で決まり、MHz 程度（例: w = 0.5 mm、TeO2 の縦波 約 4.2 km/s で 0.12 µs、UNVERIFIED）。
- 欠点: 回折効率（ファイバ結合で一般に 50–80%、UNVERIFIED）。**squeezed 光路には置けない**。参照光・probe の光路では使える（Ref A は 2 台の AOM で probe を 1 MHz ずらしているが、アクチュエータとしては EOM を使う [Ref A §3]。Liu et al. は AOM で probe をパルス化 [Liu et al. Methods]）。
- D10 の注意: 参照光（LO）側の AOM で誤差を 0 にしても、GBS に入る squeezed 光の位相は直らない。指標 (ii) には効かない。

### 3.6 レンジ飽和、ロック外れの検出、再引き込み（E83）

- **飽和**: 各段の出力を [u_min, u_max] でクリップ（制御器の出力レンジ、ドライバの入力レンジ、ユーザが設定した CV の制限 [DSC1 p.9]）。
  積分器の anti-windup は機器ごとに不明（DSC1・Moku:Delta とも資料に記載なし）→ モード `none | clamp | back_calculation` を持たせ、既定は `none`（最も素朴な実装）、感度を Phase 6 で見る。
- **ロック時間**（線形ドリフト r_drift、中央から開始）: t_lock = R/|r_drift|（R は片側のレンジ）。
  ランダムウォーク（Var = 2 D_ψ t）なら、区間 (−R, R) を中央から出るまでの平均時間 **⟨t⟩ = R²/(2 D_ψ)**。
- **ロック外れの検出**（シミュレータの真値と、実機で見える量の両方）:
  1. どれかの段が t_sat 以上続けて飽和。
  2. |e| > e_th が t_e 以上続く。DC fringe lock（Ref D）では誤差は sin 型なので、真の位相誤差が |δ| > π/2 を越えると隣の安定点（2π 先）へ滑る。真値として unwrap した位相で「縞の滑り」を数える。
  3. Ref C の量子雑音ロックでは、測った雑音分散の急変（T14 の誤差信号から）。
- **再引き込み**（モード、既定 `none`）:
  1. `reset_center`: 積分器を中央に戻し、自然に引き込む。
  2. `fringe_hop`: u ← u − n u_2π（n は中央に最も近づく整数、915B では u_2π = 0.968 V）。Ref D や F のように位相の関係が 2π を法として同じなら、ロックを保ったまま範囲を戻せる。
     ストレッチャの過渡（共振、スルーレート）の間 t_hop はデータを捨てる。指標 (ii)（ポート間）は 2π を法とするので可。
  3. `ramp_search`: ループを止め、u を [u_min, u_max] でランプ（周期 T_ramp）し、|e| < e_cap かつ傾きの符号が正しいところで積分器を入れる。
     DSC1 の RAMP モード（10–400 ms の鋸歯状波、振幅・オフセット設定）[DSC1 p.7, p.20] は手動の掃引だけで、自動の引き込みは資料にない。Moku:Delta の Laser Lock Box には「Lock Assist」（誤差の任意の零交差へ引き込む）[Moku p.42] がある。
- 出力: ロック区間の長さの分布、デューティ比、滑り・リセットの回数と時刻。

---

## 4. 制御器とデジタル部

### 4.1 ADC/DAC の量子化とレンジ（E81）

- q = V_FS/2^N。白色雑音とみなせる条件（信号が多くの LSB にわたって動く）では S_q = q²/(6 f_s)（片側、0–f_s/2）。
- ループの中では DAC の量子化は**アクチュエータ入力への外乱**: δψ_q(f) = K_s H(f) q(f)/(1 + L(f))。残差の分散
  σ_q² = ∫_0^{f_s/2} |K_s H(f)/(1 + L(f))|² S_q df。
  外乱が小さく信号が静かなときは白色近似が破れ、±1 LSB の limit cycle（振幅 〜K_s q/2–K_s q）になる。両方のモデルを検証テストで確かめる（§8 V11）。
- ADC の量子化は**センサ雑音**: δψ = q_ADC/K_D 相当。ループの帯域内ではそのまま残差になる。
- **レンジと分解能の両立**: 1 つの DAC でレンジ R_tot と分解能 δψ_min を同時に満たすには N_DAC ≥ log2(R_tot/δψ_min)。
  915B（65 rad）で 1 mrad なら 16 bit。12 bit の DSC1 では足りないので、粗い段と細かい段に分ける（TEC＋PZT、EOM＋PZT、2 チャンネルの DAC を抵抗で加算する、など）。

### 4.2 サンプリングとループ遅延（E82）

- 連続時間のモデル: L(s) = K_D C(s) Σ_k A_k(s) e^{−s τ_L}、ZOH の平均遅延 h/2 を τ_L に含めるか、離散化で厳密に扱う。
- 位相余裕の上限: 平らなプラント＋P 制御で f_u,max = (180° − PM)/(360° τ_L)。PI の積分の −90° が効く帯域ではさらに下がる。
- **DSC1**: 入力 16 bit・1 MHz、2 点平均で実効 500 kHz、PID 更新 500 kHz、入力帯域 typ. 120 kHz、出力帯域 typ. 100 kHz（低域フィルタで量子化の段差を均す）、−180° typ. 60 kHz（>58 kHz）[DSC1 p.1–2]。
  入出力の帯域を 1 次とみなすと 60 kHz での遅れは 26.6° + 31.0° = 57.5°、残りの 122.5° を純遅延とすると **τ_L = 5.7 µs**。すべてを遅延とみなす上限は 8.3 µs。
  プリセットは「1 次 LPF 120 kHz（入力）＋ 1 次 LPF 100 kHz（出力）＋ τ_L = 5.7 µs」を既定とし、τ_L = 5.7–8.3 µs を不確かさとして掃引する。
  1 kHz での遅れ 2.1–3.0°、10 kHz で 20–30°。Ref D のループ帯域を決めるのは、DSC1 の遅延より 915B の機械共振（不明）と 914 の極（§3.1.4）である可能性が高い。
- **Moku:Delta**: PID Controller の「最大帯域 400 kHz で位相遅れ 30°」[Moku p.65] → τ_L ≲ 0.21 µs。Laser Lock Box の速い PID は 78 MHz、遅い PID は 1.22 MHz でサンプル [Moku p.45–46]。
- **光の伝搬遅延**: アクチュエータから検出器までの光路 τ_opt = n_g L_opt/c（約 4.9 ns/m）。F の逆向き probe では最大 n_g L_k/c: L_31 = 63 m で 0.31 µs、630 m で 3.1 µs。Moku 系では支配的、DSC1 でも無視できない。
  Ref C では 1 km のファイバは誤差の検出点とアクチュエータ（どちらも建物 B）の上流にあるので、ループ遅延には入らない [Ref C Fig. 3]。

### 4.3 PID/IIR の構造（E84）

正準形（シミュレータ内部、連続時間、各項を独立に書く）:

  C(s) = K_p + K_i / (s + ω_IS) + K_d s / (1 + s/ω_DS)  （+ 必要なら二重積分 K_i2/s²）

- ω_IS: 積分の飽和（低周波の利得の上限、Moku の IS）。0 なら純積分。ω_DS: 微分の飽和（Moku の DS）。
- 離散化は機器ごとに切り替える:
  - **DSC1**: u[n] = K_p e[n] + h Σ_{k≤n} K_i[k] e[k] + (K_d/h)(e[n] − e[n−1])、h は ADC サンプル間隔 [DSC1 p.19]。微分のフィルタはない。
    **bumpless**: K_i を積分の中に入れるので、K_i = 0 の間は積分値が 0 のまま、K_i を入れたときに段差が出ない。I を無効にすると K_i = 0 として扱う [DSC1 p.18]。
    P・I・D の符号は連動 [DSC1 p.19]。K_p、K_i、K_d は ±10⁻⁵–10⁵ [DSC1 p.6]、GUI では K_i [1/s]、K_d [s] [DSC1 p.11]。
    **食い違い**: Ziegler–Nichols の表の値を DSC1 に入れるときは「積分項に 2×10⁻⁶ を掛け、微分項を 2×10⁻⁶ で割る」[DSC1 p.20] とあり、入力値が 1 サンプルあたり（K_i h）であることを示唆する。p.11・p.19 の「K_i [1/s] で h を掛ける」と矛盾する。
    また p.19 は「h は ADC サンプル間の時間」とするが、2×10⁻⁶ は PID 更新（500 kHz）の周期。→ プリセットは h = 2 µs、K_i の単位を `per_second | per_sample` の較正用の切り替えで持つ（§11 の質問 3）。
  - **Moku:Delta**: P は ±60 dB、積分・微分は「交差周波数」で指定。交差周波数 f_I は積分の利得が比例の利得と等しくなる周波数（P がオフなら 1）（Liquid Instruments のナレッジベース、検索結果の抜粋、2026-10-28、伝達関数の式は画像で未確認、UNVERIFIED）。
    よって K_i = K_p 2π f_I、K_d = K_p/(2π f_D)。範囲: f_I = 3.125 Hz–312.5 kHz、二重積分の交差 3.125 Hz–f_I、積分の飽和は K_p と +60 dB の間、f_D = 31.25 Hz–3.125 MHz [Moku p.65]。
  - **汎用（Ref A の自作 PID、Ref B の自作 servo、Ref C）**: Tustin（双一次）変換。
- IIR は 2 次の section（biquad）の縦続（SOS）で表す。Moku:Delta の Digital Filter Box は 305.18 kHz、4.8828 MHz、39.063 MHz で IIR、Butterworth・Chebyshev・Elliptic など、次数 2–16 [Moku p.34–36]。Laser Lock Box の低域フィルタは 78.125 MHz の IIR、2・4 次、2.601 kHz–35.16 MHz [Moku p.44]。

### 4.4 並列の帯域分割と、積み重ね（stacked cascade）

アクチュエータ k の伝達関数を A_k(s)（[rad/V]、§3）、誤差信号の傾きを K_D とする。

- **並列の帯域分割**: u_k = F_k(s) e、F_k = (HPF または LPF) · PID_k。
  L(s) = K_D Σ_k A_k F_k e^{−s τ_L}。
- **積み重ね（Ref C）**: u_1 = C_1 e、u_k = C_k u_{k−1}（k ≥ 2）、各段は自分のアクチュエータを駆動しつつ次の段の入力になる。
  L(s) = K_D [A_1 C_1 + A_2 C_2 C_1 + A_3 C_3 C_2 C_1] e^{−s τ_L}。
  C_k（k ≥ 2）が積分を持つと、定常で u_{k−1} → 0: 速い段のアクチュエータが中央に戻される（offload）。
- 線形の範囲では、積み重ねは F_1 = C_1、F_2 = C_2 C_1、F_3 = C_3 C_2 C_1 の並列と**等価**。違いが出るのは飽和・量子化・再引き込みなど非線形の場面（例: 遅い段が飽和すると、その上の段が順に飽和する）。
- **受け渡し（crossover）の条件**: 隣り合う段の受け渡し周波数 f_{k,k+1} を |A_{k+1} F_{k+1}| = |A_k F_k| で定義する。
  - 低周波で A_{k+1} F_{k+1} と A_k F_k が同じ向き（同位相）になるよう C_{k+1} の符号を選ぶ（積み重ねでは sign(C_{k+1}) = sign(A_k/A_{k+1})）。逆だと offload が正帰還になる。
  - f_{k,k+1} での 2 つの位相差が ±180° に近いと、合成したアクチュエータの応答に深い谷（零点）ができ、ループが不安定になる。位相差を 135° 以下に保つ（目安、Bechhoefer 2005 の一般論からの導出 [AE]）。積分 1 つで受け渡すと位相差は 約 90°。
- **Ref C の対応**: C_1 → ファイバ EOM（kHz 帯）、C_2 → ピエゾミラー（その下、HV アンプ、ランプを加算）、C_3 → ファイバストレッチャ（Hz–sub-Hz、レンジでロック時間が決まる）[Ref C §2]。
  1 MHz の dither は外部の発振器から EOM の駆動に加算 [Ref C Fig. 3]。図の配線は線の重なりで読みにくいが、本文の帯域の割り当てに従う。
- 構造は prompt §5.2 の入れ子のパイプラインで両方表せる: 並列 = `[hpf >> pid1 >> eom] ‖ [lpf >> pid2 >> stretcher]`、積み重ね = `pid1 >> (eom ‖ (pid2 >> (mirror ‖ (pid3 >> stretcher))))`。

### 4.5 機器プリセット

#### Thorlabs DSC1（`refs/manuals/thorlabs-dsc1.pdf`、Rev. B, 2025-11-03, DOC-103674）

| 項目 | 値 | 出典 |
|---|---|---|
| 入力 | 16 bit、±4 V、1 MΩ、1 MHz サンプル、2 点平均で実効 500 kHz、帯域 typ. 120 kHz | p.1–2 |
| 入力雑音 | 0.3 mV rms（100 Hz–250 kHz）、平均の雑音床 −120 dB V²/Hz（1 µV/√Hz）、ピーク −105 dB V²/Hz | p.2 |
| 出力 | 12 bit、±4 V、出力抵抗 220 Ω（>100 kΩ 負荷を想定、R_dev/(R_dev + 220 Ω) で縮む）、帯域 typ. 100 kHz（LPF） | p.1–2 |
| LSB | 入力 122 µV、出力 1.95 mV [計算] | — |
| 系の帯域 | DC–100 kHz、入力→出力 −180° は >58 kHz（typ. 60 kHz） | p.2 |
| 推定の純遅延 | 5.7 µs（1 次 LPF 2 つを差し引く）、上限 8.3 µs [AE] | §4.2 |
| PID | 更新 500 kHz、u[n] = K_p e + h Σ K_i e + (K_d/h) Δe、bumpless、e = S − v、符号連動、係数 ±10⁻⁵–10⁵ | p.2, 6, 18–19 |
| 設定値・オフセット | ±4 V、CV の制限を設定可（範囲外はクリップ） | p.6, 9 |
| RAMP | 鋸歯状波、振幅・オフセット ±4 V、周期 10–400 ms | p.7 |
| PEAK（dither 極値ロック） | 100 Hz–100 kHz（100 Hz 刻み）、HPF → 復調 → LPF → 積分 K [V/s] | p.2, 7–8, 20–22 |
| 自動の再引き込み | 記載なし | — |
| 外部制御 | USB、GUI と SDK あり（SDK の資料は refs にない） | p.1 |

#### Moku:Delta（`refs/manuals/MokuDelta-Specifications.pdf`、v26-0817）

| 項目 | 値 | 出典 |
|---|---|---|
| アナログ入力 | 8 ch、5 GSa/s、14 bit（高帯域）/ 20 bit（低帯域）の blended ADC、ENOB 10.8 | p.6 |
| 入力レンジ | 100 mVpp（+20 dB）、1 Vpp、10 Vpp（50 Ω）、40 Vpp（1 MΩ） | p.6 |
| 入力帯域 | 2 GHz（50 Ω）、**1 MHz（1 MΩ 入力）**、AC 結合の角 16 kHz（50 Ω）/ 1.6 Hz（1 MΩ） | p.6 |
| 入力雑音 | 10 nV/√Hz @ 100 Hz、<10 nV/√Hz > 2 MHz | p.6 |
| アナログ出力 | 8 ch、5 GSa/s、14 bit、1 Vpp（DC–200 MHz）/ 10 Vpp（DC–100 MHz）into 50 Ω、出力抵抗 50 Ω、DC 結合 | p.7 |
| 出力の LSB | 61 µV（1 Vpp）、610 µV（10 Vpp）[計算]。915B では 0.40 / 4.0 mrad | — |
| PID Controller | 4 ch、入力の制御行列 −20–+20、開ループ帯域 35.38 MHz、制御の周波数範囲 DC–40 MHz、最大帯域 400 kHz（位相遅れ 30°）、P ±60 dB、I/I+/IS/D/DS、**出力の制限 ±500 mV into 50 Ω** | p.64–65 |
| Laser Lock Box の PID | 速い: 78 MHz サンプル、f_I 3.125 Hz–312.5 kHz。遅い: 1.22 MHz、f_I 48.83 mHz–4.883 kHz、f_D 488.3 mHz–48.83 kHz。設定値 ±500 mV（100 µV 刻み）、走査（ランプ・三角波）1 mHz–10 MHz、±500 mV | p.45–46 |
| 復調 | 内部・外部参照、PLL（帯域 1 Hz–1 MHz）、LPF は 2・4 次 IIR | p.43–44, 48 |
| Multi-Instrument | 3/5/8 スロット、入出力 5 / 1.25 / 0.625 GSa/s、スロット間 16 bit、出力レンジ 1 Vpp / 10 Vpp | p.78–80 |
| 遅延 | 明記なし。400 kHz で 30° から τ ≲ 0.21 µs [AE] | p.65 |

注:
- PID Controller の出力制限 ±0.5 V（50 Ω）は 914 の入力 ±5 V の 10%（高インピーダンス負荷で電圧が 2 倍になるとしても 20%）。そのままでは 915B のレンジは ±3.2–6.5 rad。PID 計器で 10 Vpp の出力レンジが使えるかは資料から読めない（§11 の質問 5）。
- 32 ポートには出力 8 ch × 4 台、または別の多チャンネル DAC が要る（DSC1 は 1 台 1 チャンネル）。

#### Moku:Pro（Ref C、D25）

- Ref C は Moku:Pro で復調・3 つの PID・加算器を実装 [Ref C §2]。論文には PID の利得・交差周波数・遅延・サンプル周波数・各アクチュエータのレンジの記載がない。
- §7-10 の検証で必要なのは、(i) dither の 0.0555 rad rms が残差に出ること、(ii) 全体で 約 0.1 rad、(iii) ロック時間が遅い段のレンジで決まること。どれも制御器の機種には依存しない。
  kHz 帯の EOM ループで、遅延 0.1–2 µs の差は 10 kHz でも 0.4–7° しか変えない。
- **結論: Moku:Pro のプリセットは作らない**。Ref C は「汎用 FPGA」プリセット（τ_L = 0.2–1 µs、16 bit 相当、レンジは Moku:Delta の値で代用）で設定し、検証テストに「τ_L を 0.1–2 µs で変えても §7-10 の結果が 10% 以内で変わらない」を加える（§8 V19）。

### 4.6 抽象制御器インターフェース（prompt §5.3）

シミュレータ（`SimController`）と実機のラッパ（`MokuController`、`DSC1Controller`）が同じインターフェースを実装し、BoTorch のループはどちらも同じ呼び方で動かす。

1. 構成（`configure(config)`、機種に依存しない値だけ）
   - チャンネルごとの入力: f_s、N_ADC、レンジ、結合（AC/DC）、入力インピーダンス、雑音密度、アンチエイリアス LPF。
   - チャンネルごとの出力: N_DAC、レンジ、出力インピーダンスと負荷、出力 LPF、CV の上下限、オフセット。
   - 処理: 更新周期、純遅延 τ_L、経路（入力 → フィルタ → PID → 出力の行列、積み重ねの接続、加算器）。
   - 各 PID の正準係数（K_p、K_i [1/s]、K_d [s]、ω_IS、ω_DS、二重積分）、設定値、符号、anti-windup のモード。
   - 補助: dither（周波数、振幅、加える出力）、復調（参照位相、LPF）、ランプ（周期、振幅、オフセット）。
   - ロック監視: 飽和・誤差のしきい値、待ち時間、再引き込みのモードとそのパラメータ（§3.6）。
2. 機種への写像（プリセットの YAML ＋ 変換関数）: 正準係数 → DSC1（K_i の単位の扱い、e = S − v による符号反転、h = 2 µs）、Moku（dB と交差周波数）。写像できない値（範囲外）は警告して最も近い値に丸め、丸めた値を返す。
3. 実行時の操作: `set_params(dict)` / `get_params()`、`engage()` / `disengage()`、`reset_integrators()`、`ramp(...)`。
4. 測定: `acquire(duration, signals)` → 誤差、各段の出力、飽和フラグ、ロック状態の時系列。`measure_transfer_function(injection_point, freqs)`（掃引正弦の注入、§7-7 の検証にも使う）。
5. 状態: `status()` → locked/unlocked、飽和中の段、滑り・リセットの回数。
6. 記録: 構成、機種、ファームウェア（実機）、シード、git commit（prompt §6）。

---

## 5. 時間ループの実行方式（prompt §5.4）の推奨

**推奨: (c) コンパイルした 1 サンプルごとの閉ループ ＋ (a) 速い RF 搬送波の baseband/統計モデル ＋ 遅い段の多重レート化。(b) は開ループ部分にだけ使う。**

理由:
1. **ループ遅延の正確さ**: (b) は「チャンク長 ≤ 最小のループ遅延」でないとチャンク内の帰還を落とす。DSC1（τ_L ≈ 5.7 µs）なら f_sim = 1 MHz で 5 サンプル、Moku 系（≲ 0.21 µs）では 1 サンプル未満になり、(b) は (c) と同じになる。
   (c) なら遅延は「整数サンプルのリングバッファ ＋ 端数の補間（線形補間、または 1 次の Thiran 全域通過）」で任意の τ を表せる。
2. **必要な時間の長さ**: Ref C のロック時間は数分、Ref D も §7.1 の概算で数分–数十分。f_sim = 1 MHz で 60 s は 6×10⁷ ステップ。
   Python で 1 ステップ・1 ブロックあたり 5–20 µs かかると 10 ブロックで 1 時間以上、numba などで 1 ステップ 50–200 ns なら 3–12 s [AE、Phase 4 で計測して確かめる]。
3. **MHz 以上の搬送波は (a) でしか扱えない**: Ref C の雑音電力検出は 8–300 MHz、ESA の中心 3 MHz・RBW 1 MHz。これを直接サンプルすると f_sim ≥ 600 MHz で現実的でない。
   復調後（または包絡線検出・ログアンプの後）の量を、瞬時の位相の関数（T14、T04）＋有限帯域・有限時間の推定の揺らぎ（radiometer 雑音、T14）として baseband で作る。ESA のゼロスパンも同じ（VBW で平均した雑音電力）。
   Ref A の 1 MHz ヘテロダインと Ref C の 1 MHz dither は、f_sim = 8–10 MHz なら直接表せる（numba で 60 s が 30–120 s [AE]）。直接と baseband の両方を持てば、互いの検証に使える。
4. **f_sim の選び方**: f_sim ≥ max(プリセットの PID 更新周期の逆数（DSC1: 500 kHz、ADC は 1 MHz で 2 点平均）、最高のモデル化した機械共振の 20 倍、f_u の 10 倍)。
   DSC1 プリセットでは f_sim = 1 MHz（ADC）で、PID は 2 サンプルごとに更新するとビット単位で同じ離散則になる。アクチュエータの連続時間モデルは ZOH の厳密な離散化（`scipy.signal.cont2discrete`）。
5. **多重レート**: TEC・遅い PID（mHz–Hz）は間引いたレートで回す（Moku:Delta の Laser Lock Box も遅い PID は 1.22 MHz で別レート [Moku p.46]）。
6. **保護リスト（AGENTS.md）との関係**: ブロックの factory と `step(state, params) -> state` の形、入れ子のパイプライン、遅延バッファは残す。コンパイルするのは各ブロックの内部と、パイプラインを組んだ後の内側のループ。
   同じ構成を純 Python の参照実装でも短時間動かし、一致を確かめる（§8 V20）。この設計は後から変えると高くつくので、Phase 3 の計画の確認（grill）の対象にすることを提案する（§11 の質問 6）。

---

## 6. 周波数帯とスケーリング

| 効果 | 周波数帯 | ファイバ長 L | パワー P | 温度 T | パルス幅 τ_p |
|---|---|---|---|---|---|
| ストレッチャの利得 K_s | DC–f_m | 巻線長 L_w に比例（光源ポートの L_k とは独立） | なし | c_T 経由で巻線部の熱ドリフト ∝ L_w | なし |
| 機械共振 | 5–50 kHz（不明、掃引） | 巻き枠の大きさ | なし | 弱い | なし |
| ドライバの帯域・スルー | 30–110 kHz（914）、大振幅で 1/f | なし | なし | なし | なし |
| EOM | DC–GHz（servo 帯は制御器が上限） | なし | 挿入損失は squeezed 光で致命的 | 焦電・ドリフト（UNVERIFIED） | パルスでも同じ（帯域内） |
| PZT ミラー | DC–14/22 kHz（Ref B） | なし | なし | 弱い | なし |
| TEC | mHz | L_h に比例したレンジ | なし | c_T（T08） | なし |
| AOM | DC–MHz | なし | 回折効率 | なし | probe のパルス化に使う（Liu） |
| DAC 量子化 | 白色 0–f_s/2、ループで整形 | 必要なレンジ ∝ c_T L_k ΔT → 必要ビット数が log2 で増える | 誤差の傾き K_D ∝ P（ADC 側） | — | — |
| ループ遅延 | f_u 付近 | 光の遅延 n_g L/c が加わる | — | — | — |
| レンジ飽和 | < 1 Hz | ドリフト ∝ L、ロック時間 ∝ R/L | — | ドリフトの速さ ∝ dT/dt | — |

---

## 7. 概算（自分たちのパラメータ）

### 7.1 Ref D（Cobrite DX1 → PM ファイバ MZI、915B ＋ 914、PDB450C、DSC1）

- K_s = 6.49 rad/V（8669-1 と仮定。どちらの個体かは不明）。DSC1 でのレンジ ±26.0 rad、914 本来のレンジ ±32.5 rad。
- DAC の LSB は 12.7 mrad。白色量子化雑音: S_q = (1.95 mV)²/(6 × 500 kHz) = 1.27×10⁻¹² V²/Hz → 位相で 5.4×10⁻¹¹ rad²/Hz。914 の帯域（等価雑音帯域 ≈ 28 kHz × π/2）までで **約 1.5 mrad rms**、量子化の段の幅そのままなら最大 ±6.3 mrad [AE]。
- ADC の LSB は 37 µrad、入力雑音 0.3 mV rms は 0.09 mrad（K_D = 3.3 V/rad、survey の上限）。
- DSC1 の出力雑音を −120 dB V²/Hz（1 µV/√Hz、入力側か出力側かは資料では不明）とすると、位相で 4×10⁻¹¹ rad²/Hz、量子化と同じ桁。
- 遅延 5.7–8.3 µs、914 の 2 極（29.5 kHz、112 kHz）: 合わせた位相遅れは 1 kHz で 4.6–5.5°、5 kHz で 22–27°。915B の共振が 20 kHz、Q = 10 なら PI の f_u は 2–4 kHz 程度が上限の目安 [AE]。
- **熱ドリフトとロック時間** [AE、c_T は T08 の値待ち]: 腕の長さの差 ΔL ≈ 8.4 m（巻線長）として、2 つの腕が同じだけ温まるときの位相差の係数は c_T ΔL ≈ 44 × 8.4 ≈ **3.7×10² rad/K**。
  レンジ ±26 rad は ±0.070 K に当たる。室温のドリフト 0.1 K/h なら 約 42 分、1 K/h なら 約 4 分でレンジを使い切る（D27 のドリフト速度を掃引）。
  915B の巻き枠そのものの熱膨張によるひずみ（ファイバの c_T に加わる）は資料にないので、T08 で扱う。

### 7.2 Ref C

- dither: 振幅 0.0785 rad（4.50°）、RMS 0.0555 rad。全体 約 0.1 rad から、dither 以外の寄与は √(0.1² − 0.0555²) ≈ 0.083 rad [AE]。
- EOM のレンジ: V_π = 4 V で 0.785 rad/V。出力 ±1 V なら ±0.79 rad。遅いドリフトは持てないので、積み重ねで下の段に渡すのが必須。
- 遅い段（ファイバストレッチャ）のレンジ R_3 は記載がない。「数分」のロック時間から r_drift ≈ R_3/(数分)。仮に 915B と同じ 65 rad（±32 rad）なら r_drift ≈ 0.1–0.2 rad/s（中央から始め、2–5 分）。T08 の 1 km（敷設 180 m）の差動ドリフトの値と照合する。
- 遅延 0.2–1 µs: 10 kHz の EOM ループで 0.7–3.6°。

### 7.3 将来の構成 F（32 ポート、D07–D13、D30–D32）

値がないもの（f_rep、パルス幅など）は survey と同じ仮の値（f_rep = 100 MHz → L_31 = 63 m）を使う。
- 位相の熱係数: c_T L_31 ≈ 44 × 63 ≈ **2.8×10³ rad/K**（T08 待ち）。915B 1 台（±32 rad）で追えるのは ±0.012 K。
  - TEC（L_h = 3 m、±10 K）: ±1.3×10³ rad。L_31 では ±0.47 K に相当。スプール全体を温度制御すれば L_k 全体で ±c_T L_k ΔT_max。
  - fringe hop: ドリフト 0.1 K/h なら 2.8×10² rad/h、±32 rad のレンジで 約 7 分ごとにリセット。指標 (ii) は 2π を法とするので可、ただしリセット中（t_hop）のデータを捨てる。
- 分解能: 指標 (ii) の目標が θ_rms ≈ 数 mrad（D16 の表は T04）なら、DAC の 1 LSB は ≲ 1 mrad。915B のレンジ全体で 16 bit。12 bit の DSC1 は不可、Moku:Delta の 14 bit（1 Vpp で 0.40 mrad/LSB、ただしレンジは ±3.2 rad）は細かい段向き。
- 挿入損失: squeezed 光路に EOM は置けない（§3.2）。915B は融着で入れればファイバの曲げ損失だけ（値は資料にない、UNVERIFIED）。FC/APC のコネクタ 2 つ（D08）を挟むと各 0.1–0.3 dB（一般値、UNVERIFIED）で、T12 の損失収支に入れる。
- 光の遅延: 0.31 µs（63 m）、3.1 µs（f_rep = 10 MHz で 630 m）。
- 逆向き probe の誤差信号（Liu 型、D30）は、検出点がアクチュエータの上流にあるので、ループ遅延 = τ_ctrl + n_g L_k/c（アクチュエータから先のポート光路の長さぶん）。

---

## 8. 検証テスト

数値は §3–§7 の値。許容誤差は理由を付けて Phase 3 で確定する。

- **V1 ストレッチャの静的利得**: λ0 = 1548.8 nm、ΔOPL_full = 16 µm、ΔV_in,full = 10 V → u = +1 V で δψ = +6.491 rad（相対 10⁻⁹）。914 で ±5 V、DSC1 で ±4 V にクリップ → レンジ 51.93 rad。
- **V2 ξ_eo**: 物理的な伸びを入力する形で、p11 = 0.121、p12 = 0.270、ν = 0.17、n = 1.458 → ξ_eo = 0.7837、δψ/δL = ξ_eo β。
- **V3 機械共振の離散化**: 1 モード、f_m = 20 kHz、Q = 10、f_sim = 1 MHz: |H_mech(f_m)| = Q ± 1%、位相 −90° ± 0.5°。
- **V4 914 の応答**: ±5 V の正弦入力で、§3.1.4 の 6 点（51.7、50.1、42.3、36.2、23.7、11.0 V）を ±5% で再現（線形 2 極モデル）。スルーのモデルでは 30–100 kHz の 3 点を ±10%。
- **V5 遅延バッファ**: τ = 5.7 µs、f_sim = 1 MHz（端数 0.7 サンプル）: 伝達関数の位相が −360° f τ に対して f ≤ f_sim/10 で ±1°、群遅延 5.7 ± 0.05 µs。
- **V6 DSC1 プリセット**: 入力 → 出力の位相が −180° になる周波数が 60 ± 3 kHz（仕様 >58 kHz、typ. 60 kHz）。
- **V7 Moku:Delta プリセット**: 400 kHz での位相遅れ ≤ 30°。
- **V8 離散 PID**: DSC1 の離散則と連続の PID の周波数応答が f ≤ f_s/20 で 1% 以内。
- **V9 bumpless**: K_i を 0 から有限値に切り替えた直後の |Δu| < 1 LSB。
- **V10 閉ループ（prompt §7-7）**: 量子化・飽和なしの線形ループ（P = K_s H_drv H_mech e^{−sτ}、PI）で、白色雑音を注入して求めた S(f) = 1/(1 + L) が解析解と f ≤ f_sim/10 で 0.5 dB・3° 以内。フィードバックを切ると残差 = free-running の雑音（分散の比 1、統計誤差 ±2/√N_avg 以内）。
- **V11 量子化**: (a) 雑音なし・一定の外乱で、定常の位相誤差 |δ| ≤ K_s q = 12.7 mrad。(b) 外乱が q より十分大きいとき、量子化を入れた分の分散の増加が ∫|K_s H/(1 + L)|² q²/(6 f_s) df と 10% 以内で一致。
- **V12 レンジ飽和とロック時間（§7-10 の後半）**: 線形ドリフト r = 0.1 rad/s、R = 25.96 rad、中央から → ロック外れが t = 259.6 s ± 1 サンプル。ランダムウォーク（D_ψ）で、200 回以上の試行の平均が R²/(2 D_ψ) と 5% 以内。
- **V13 積み重ね = 並列（線形）**: C_1、C_2、C_3 の積み重ねと F_k = C_k ⋯ C_1 の並列で、出力が 10⁻⁹ 以内で一致。一定の外乱で、C_2 が積分を持つと u_1 → 0（offload）。
- **V14 符号の検査**: C_2 の符号を反転すると offload が発散する（または構成の検査が警告する）。
- **V15 Ref C の dither（§7-10 の前半）**: V_π = 4 V、200 mVpp、1 MHz → 位相変調の振幅 0.0785 rad（4.50°）、RMS 0.0555 rad。他の雑音を切った閉ループで、dither が残差に 0.0555 rad rms として現れる（ループは 1 MHz では効かない）。
- **V16 再引き込み**: fringe hop で位相が n·2π だけ変わる（線形モデルで誤差 < 10⁻⁶ rad）、t_hop + 5/f_u 以内にロックが戻る。ramp_search で、正しい傾きの零交差でだけ引き込む。
- **V17 AOM**: δf_RF = 1 Hz（m = +1）で位相が −2π t rad。φ_RF の段差で光の位相が −m φ_RF だけ跳ぶ。
- **V18 TEC**: ステップ応答の時定数 τ_th、直流の利得 K_T c_T L_h。PZT ＋ TEC の積み重ねで、ランプ状のドリフトの下で PZT の出力が中央に戻る。
- **V19 Ref C の制御器の機種に依存しないこと**: τ_L を 0.1–2 µs、DAC を 14–16 bit で変えても §7-10 の残差が 10% 以内で同じ（Moku:Pro プリセット不要の根拠）。
- **V20 実行方式**: コンパイルした実装と純 Python の参照実装で、同じシード・同じ構成の 10 ms の結果が 10⁻¹² 以内で一致。f_sim を 2 倍にして残差の PSD が f ≤ f_sim/10 で 5% 以内で同じ（数値の収束）。

---

## 9. シミュレータに必要なもの

| ブロック | 段 | 入力 | 出力 | 主なパラメータ | オン・オフ |
|---|---|---|---|---|---|
| `pzt_fiber_stretcher` | T（挿入損失は O） | u（遅延バッファ経由） | 区間の δψ（軸ごと）、飽和フラグ | K_s（または ΔOPL_full、ΔV_in,full、λ0）、σ_act、モード（f_i、Q_i、r_i）、ドライバ（f_d1、f_d2、SR、入力レンジ、R_in）、h_hys、γ_cr、区間 id、軸ごとの係数（T09） | resonance、driver、slew、hysteresis、creep、saturation |
| `eom_phase` | T（損失は O） | u | δψ | V_π,LF、V_π,RF、入力レンジ、IL、RAM 係数（T14） | — |
| `pzt_mirror` | T | u | δψ | d_x、θ_i、モード、ドライバ | resonance |
| `thermal_actuator`（TEC） | T | u | δψ、ΔT_h | K_T、τ_th、τ_tr、ΔT_max、L_h、c_T（T08） | — |
| `aom_shift` | T | δf_RF または φ_RF | δψ | m、f_RF、τ_a、η_AO | — |
| `adc` / `dac` | T | 電圧 | 量子化した電圧 | N、レンジ、f_s、平均の点数、雑音密度、LPF、出力抵抗と負荷 | quantization、noise |
| `delay_buffer` | T | 任意 | 遅らせた値 | τ（整数＋端数）、補間の方式 | — |
| `pid` | T | e | u | 正準係数、離散化の方式（dsc1/moku/tustin）、設定値、符号、制限、anti-windup | P、I、D 個別 |
| `iir` / `hpf` / `lpf` | T | 任意 | 任意 | SOS 係数 | — |
| `adder` / `ramp` / `dither` | T | — | — | 周波数、振幅、周期 | — |
| `lock_supervisor` | T | e、各段の u・飽和、真の位相 | ロック状態、イベント | しきい値、待ち時間、再引き込みのモード | relock |
| 機器プリセット | 設定ファイル | — | — | §4.5 の表（`dsc1.yaml`、`moku_delta.yaml`、`fpga_generic.yaml`） | — |
| オフラインの表 | O | — | — | 914 の Fig. 1 の読み取り値（V4 用）、ξ_eo の係数の組 | — |

出力（per run）: 各アクチュエータの u(t)・δψ(t)・飽和の割合、各ロックの残差の時系列と PSD、注入で測った L(f)・S(f)、ロック区間とイベント。

---

## 10. 他のトピックへの依存と仮の値

- **T08**（環境の位相雑音）: c_T と温度ドリフトのスペクトル。仮に c_T = 44 rad/(m·K)（survey の AE）、ドリフト 0.1–1 K/h（D27 の掃引）を使った。
- **T09**（偏光・PM）: ストレッチャの軸ごとの係数、伸ばしたときの 2 軸の位相差。
- **T12**（損失）: ストレッチャ・コネクタ・EOM の挿入損失の収支。
- **T14**（誤差信号）: K_D、RAM、Ref C の radiometer 雑音、dither と復調。§5 の baseband 化はこのモデルに依存。
- **T15**（Ref D）: 本ノートの 915B/914/DSC1 の値を使う。ΔL ≳ 8 m の見積もり（§3.1.3）を渡す。
- **T04**: 目標 squeezing → 必要な θ_rms（§7.3 の分解能の要求の基準）。

---

## 11. 調整役への提案と、ユーザの決定が必要な点

### 提案（coordinator-only のファイルへの変更案）

- `notes/conventions.md` §4: 「ξ_eo ≈ 0.78（値は T16 で確認）」→「ξ_eo = 0.78–0.80（バルクの係数で 0.784、Bertholds & Dändliker 1988 の係数で 0.795–0.798）。資料が光路長の変化 ΔOPL を与える場合は δψ = k0 ΔOPL を使い、ξ_eo は掛けない」。
- `notes/conventions.md` に追記: 「アクチュエータの正の向き = 駆動電圧を上げると ψ が増える（σ_act で較正）」、「DSC1 の誤差は e = S − v で、§9 の e と符号が逆」。
- `notes/survey.md` E78・E80: 「資料待ち」→ 本ノートの値。E82: −180° @ 60 kHz からの純遅延は 5.7 µs（LPF を差し引く）–8.3 µs。
- `DECISIONS.md` D25: 「T16 の結果、Moku:Pro プリセットは不要。Ref C は汎用 FPGA プリセットで設定」を追記。

### ユーザの決定が必要な点（推測で埋めていない）

1. **915B の「16 µm / 17 µm」の意味**: 光路長の変化（K_s = 6.49 rad/V）か物理的な伸び（7.53 rad/V）か。資料の「optical path length measurement」からは前者と読むのが自然だが確定できない。また Ref D に入っているのは S/N 8669-1（16 µm）か 8669-2（17 µm）か。販売代理店（ハナムラオプティクス）または Evanescent Optics への問い合わせが確実。
2. **資料にない 915B/914 の値**（機械共振、PZT 容量、914 の小信号帯域・出力抵抗・入力抵抗・出力雑音、ヒステリシス）: (a) メーカーに問い合わせる、(b) photonics.com の掲載（取得できず UNVERIFIED）を使う、(c) 掃引パラメータのまま進める、(d) 実測する（ただし D04 のブラインド検証と衝突しないか）のどれにするか。当面は (c) で進める。
3. **DSC1 の K_i・K_d の単位**: p.11・p.19（K_i [1/s]、h を掛ける）と p.20（Ziegler–Nichols の値に 2×10⁻⁶ を掛ける）が食い違う。SDK の資料かメーカーへの問い合わせで決める。当面は切り替え可能な較正値にする。
4. 914 の入力インピーダンス（DSC1 の 220 Ω との分圧）。影響は小さいと見込むが資料にない。
5. **Moku:Delta の PID 計器の出力 ±0.5 V（50 Ω）**: 10 Vpp の出力レンジが PID 計器で使えるか。使えないと、Moku:Delta で 915B を駆動したときのレンジが 1/10–1/5 になる。将来の機材選び（D15）に効く。
6. **§5.4 の方式（コンパイルしたループ、numba などの依存の追加、ブロック内部をコンパイル可能な形に制限すること）**は、後から変えると高くつく設計判断なので、Phase 3 の計画の確認で決めてほしい。
7. F の大レンジの手段（TEC・スプールの温度制御・fringe hop）は D14・D30 の比較（Phase 6）の候補に加えることを提案。

---

## 12. 文献

一次文献・論文:
- C. D. Butter and G. B. Hocker, "Fiber optics strain gauge," Appl. Opt. 17(18), 2867–2869 (1978). doi:10.1364/AO.17.002867（書誌は検索で確認、本文は未読）。
- A. Bertholds and R. Dändliker, "Determination of the individual strain-optic coefficients in single-mode optical fibers," J. Lightwave Technol. 6(1), 17–20 (1988). doi:10.1109/50.3956（書誌と p11 = 0.113、p12 = 0.252 は EPFL Infoscience の記録で確認、本文は未読）。
- J. Bechhoefer, "Feedback for physicists: A tutorial essay on control," Rev. Mod. Phys. 77, 783–836 (2005). doi:10.1103/RevModPhys.77.783（書誌は確認、§4.4 の受け渡しの条件は一般論からの本ノートの導出）。
- Ref B: Y. Takeno, M. Yukawa, H. Yonezawa, A. Furusawa, Opt. Express 15, 4321–4327 (2007)（`refs/`、§3 を使用）。
- Ref C: S. Verclas, C. Darsow-Fromm, B. Tohermes, R. Schnabel, Optica Quantum 4, 1 (2026). doi:10.1364/OPTICAQ.567418（`refs/`、§2–4、Fig. 3 を使用）。
- Ref A: K. Hirota, T. Kashiwazaki et al.（`refs/`、§2–3 を使用。掲載誌・DOI は survey と同じく UNVERIFIED）。
- H.-L. Liu et al., "Gaussian boson sampling with 1,024 squeezed states in 8,176 modes," Nature (2026). doi:10.1038/s41586-026-10523-6（`refs/`、Methods「Phase stabilization of the system」を使用）。
- T. I. Laakso, V. Välimäki, M. Karjalainen, U. K. Laine, "Splitting the unit delay," IEEE Signal Process. Mag. 13(1), 30–60 (1996)（端数遅延・Thiran フィルタ。**UNVERIFIED**、DOI 未確認）。

メーカー資料（D24）:
- Evanescent Optics / ハナムラオプティクス, "Model 914-2 Piezo Controller" と 915B Piezo Assembly S/N 8669-1, 8669-2 の仕様（スキャン 3 ページ）、`refs/manuals/914-915B-Specifications.pdf`。
- Thorlabs, "DSC1 Compact Digital Servo Controller User Guide," Rev. B, 2025-11-03, DOC-103674、`refs/manuals/thorlabs-dsc1.pdf`。
- Liquid Instruments, "Moku:Delta Specifications," v26-0817、`refs/manuals/MokuDelta-Specifications.pdf`。
- photonics.com, "Model 915B Piezo-electric Fiber Stretcher"（Evanescent Optics の掲載）, https://www.photonics.com/Buyers-Guide/Products/Model-915B-Piezo-electric-Fiber-Stretcher/psp9424 、2026-10-28 閲覧を試みたが HTTP 403。値は検索結果の抜粋のみ（**UNVERIFIED**）。
- gophotonics, "MPZ-LN-10-00-P-P-FA-FA"（iXblue/exail の販売店掲載）, https://www.gophotonics.com/products/electro-optic-phase-modulators/ixblue-photonics/174-854-mpz-ln-10-00-p-p-fa-fa 、2026-10-28、検索結果の抜粋のみ（**UNVERIFIED**）。
- Liquid Instruments Knowledge Base, "How do I convert Ki, Kp, Kd to crossover frequencies or vice versa," https://knowledge.liquidinstruments.com/pid-controller/how-do-i-convert-ki-kp-kd-to-crossover-frequencies-or-vice-versa 、2026-10-28、検索結果の抜粋のみ（**UNVERIFIED**）。

記憶からの一般値で原典を確認していないもの（すべて **UNVERIFIED**）: バルク石英の p11・p12・ν の組、PZT のヒステリシス 10–15%・クリープ係数、PZT の音速と円筒のフープモードの式の数値、TeO2 の音速、AOM の回折効率、コネクタの損失、被覆の熱拡散率。
