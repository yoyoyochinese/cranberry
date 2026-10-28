# T05 レーザー雑音と遅延差のある干渉（laser-noise-delay-interferometry）

- Phase: 1（deep dive）、トピック T05（`notes/survey.md` の E12–E18、DECISIONS D18-1、D31-3）
- 日付: 2026-10-28
- 入力: `docs/squeezing_sim_prompt.md`、`DECISIONS.md`（D01–D33）、`notes/conventions.md`、`notes/survey.md`、
  Ref A（Hirota 2026）、Ref B（Takeno 2007）、Liu 2026、`refs/manuals/914-915B-Specifications.pdf`、
  ID Photonics の CoBrite データシートと DX1 マニュアル（URL は §14）
- 値の印: **[manual]** = メーカー資料、**[Ref X]** = リファレンス本文、**[AE]** = このノートの担当による概算（未検証）、
  **[UNVERIFIED]** = 原典で確認していない値・文献情報

## 0. 要点

1. **Cobrite DX1 の値（メーカー資料）**: 線幅 FWHM < 100 kHz（N・GC タイプ）、80 kHz typ.（SC タイプ）。
   「instantaneous（3.5 µs）」線幅 25 kHz typ.（N・GC）、< 100 kHz（SC）。RIN −145 dB/Hz（N・GC、10 MHz–22 GHz）、
   −140 dB/Hz（SC、100 kHz–20 MHz）。**周波数雑音 PSD そのものは公表されていない。100 kHz より下の RIN も公表されていない。**
   1 時間の周波数確度 0.3 GHz。周波数を固定するための dither tone があり（振幅・周波数は非公開）、
   SBS 抑圧機能（振幅 0–1 GHz）がある。ユーザの DX1 がどのタイプかは分かっていない（§13）。
2. **Ref D の ΔL**: 915B の PDF には巻数（38）と伸び（16 µm または 17 µm、光路長）だけがあり、巻いたファイバの長さは書かれていない。
   メーカーの製品掲載文には「9 m のファイバで最大 16 µm」とある（直接は取得できず、検索結果の抜粋でのみ確認）。9 m / 38 巻 = 0.237 m/巻。
   **ΔL の名目値は 10 m、掃引範囲は 5–15 m** を提案する。ΔL = 10 m なら τ = 49 ns で、白色 FM（25 kHz）による位相差雑音の床は
   2.7×10⁻⁵ rad/√Hz、100 kHz までの積分で 8.7 mrad rms（100 kHz 線幅なら 17 mrad）。
   さらに、仕様と矛盾しない flicker FM（h₋₁ ≈ 3–4×10⁸ Hz²、[AE]）があると、10 Hz–100 kHz で約 20 mrad になる。
   1 時間の周波数確度 0.3 GHz が実際のドリフトとして現れると、MZI の位相は 92 rad 動き、915B の全レンジ（16 µm → 65 rad）を超える。
3. **ポート間の相関の崩れ（D31-3）**: 白色 FM で σ² = 2πΔν·k·T_rep。32 ポートで最悪のポート（k = 31）を σ_t 以下にするには、
   **Δν ≤ σ_t²·f_rep / (2π·31)** が必要。f_rep = 100 MHz、σ_t = 10 mrad で Δν ≤ 51 Hz。DX1（25 kHz）では σ = 0.22 rad で、
   Ref A の損失のもとでは squeezing が完全に消える。
4. **D31-3 の修正の提案**: Liu 型ロック（D30）で probe が同じレーザーから同じ周波数で作られていれば、ループ帯域より下のレーザー周波数雑音
   （ドリフト、flicker の大部分）はロックで**打ち消される**（§3.4 で導出）。打ち消されないのはループ帯域より上の成分で、白色 FM では分散の
   ほぼすべて（帯域より下にあるのは 2τf_BW ≈ 0.06% だけ）がそこにある。したがって線幅への要求は白色 FM の床（h0）に対する要求になる。
   一方、その同じロックは、遅延を合わせた LO_k に対する指標 (i) に 2πδν(τ_k + τ_p − τ_r) のずれを書き込む。
   **指標 (i) と (ii) は、レーザー周波数のドリフトに同時には強くできない。**
5. **ASE**: LO の ASE は 2 つの経路で効く。(a) ASE が広帯域の squeezed 光とビートする経路は、等価損失 ε に (V₊ + V₋)/2 を掛けた量が混ざる効果になる。
   Ref A の R± では、ε = 10⁻³ で −10.6 dB が −8.4 dB に落ちる。(b) LO と ASE のビートは CMRR を通って漏れる。
   フィルタキャビティは (a) を π/(2F) 倍に抑える（F はフィネス）。probe に乗った ASE と、probe のレーザー位相雑音は、
   BHD で squeezed 光と同じモードに入る。この効果はサーベイの表にないので、追加を提案する（§12）。
6. **フィルタキャビティのロック残差**: LO の位相は arg T = arctan(2δ/Δν_FC) ずれる。Δν_FC = 1 MHz なら、1 kHz rms の離調で 2 mrad になる。
   キャビティは LO の位相雑音に対しては低域通過フィルタ（遅延 τ_eff = 1/(πΔν_FC) の干渉計と同じ形）として働く。

## 1. 記号（symbols）

- `t`: lab time [s]. `f`: Fourier frequency of slow noise [Hz]. `ν0 = c/λ0`: optical carrier frequency [Hz].
- `φ_L(t)`: phase of the laser complex amplitude at the laser output, `A_L = |A_L| e^{iφ_L}` (conventions §4).
- `δν(t)`: instantaneous optical frequency deviation from ν0 [Hz]. With `E ∝ A e^{−iω0 t}`, **`φ_L(t) = φ_L(0) − 2π ∫_0^t δν(t') dt'`**.
- `S_ν(f)`: one-sided PSD of δν [Hz²/Hz]. `S_φ(f) = S_ν(f)/f²` [rad²/Hz] (conventions §7).
- `h0` [Hz²/Hz]: white-FM level. `h₋₁` [Hz²]: flicker-FM coefficient (`S_ν ⊃ h₋₁/f`). `Δν = π h0`: Lorentzian FWHM of the white-FM part.
- `ℒ(f)` [dBc/Hz]: single-sideband phase noise, `ℒ = 10 log10(S_φ/2)` (one-sided S_φ, small-angle).
- `S_RIN(f)` [1/Hz]: one-sided relative intensity noise PSD, `δP/P`.
- `τ_b`: group delay from the laser (or the common split point) to the detection point along beam path b [s]. `τ = τ_a − τ_b`.
- `n_g`: group index (placeholder 1.468, SMF-28-like, **UNVERIFIED for PM15-U25D**; T09 should supply). `ΔL`: Ref D arm length difference. `δL_k`: residual mismatch between squeezed fiber k and LO fiber k (D09).
- `S_Δφ(f)`: one-sided PSD of the phase difference of two beams derived from the same laser [rad²/Hz].
- `f_rep = 1/T_rep`: repetition rate. `N = 32`: number of ports. `τ_k ≈ k T_rep`: emission-time offset of port k relative to port 0 (D08).
- `Θ_k = θ_sq,k − θ_sq,0`: metric (ii) (conventions §10). `Δ_k = φ_LO,k − θ_sq,k`: metric (i).
- `τ_p`: group delay laser → GBS-side probe split point (Liu-type, D30). `τ_r`: group delay laser → reference (LO) interference point.
- `G(f)`: open-loop gain of a lock. `S_lock(f) = 1/(1 + G(f))`: sensitivity function. `f_BW`: lock bandwidth (unity-gain frequency).
- `σ²`: variance of a phase difference over a stated band [rad²]. `σ_t`: target RMS phase error allocated to a mechanism.
- `ν_n = n f_rep + f_ceo`: comb line n. `f_ceo`: carrier-envelope offset frequency. `S_ν,n(f)`: frequency-noise PSD of comb line n.
- `Δν_FC`: filter-cavity FWHM linewidth [Hz]. `FSR`: free spectral range. `F = FSR/Δν_FC`: finesse. `δ = ν_L − ν_res`: cavity detuning [Hz].
- `T_FC(δ)`: complex amplitude transmission of the filter cavity. `H_FC(f)`: transfer function of laser phase fluctuations through the cavity.
- `S_ASE` [W/Hz]: one-sided ASE PSD **in the polarization and spatial mode of the LO**, near the carrier. `n_ASE = S_ASE/(hν0)`: ASE photons per mode.
- `P_LO`, `P_pr`: LO and residual-probe powers at the BHD. `ℜ = η e/(hν0)`: responsivity. `η`: detector QE. `e`: elementary charge. `h`: Planck constant.
- `ε`: normalized common-mode leakage of the balanced detector (`i_diff = ε·i_sum` for common-mode signals; CMRR = −20 log10 ε).
- `ε_ASE`: equivalent loss due to ASE beating with broadband squeezed light (§5).
- `V₊`, `V₋`: anti-squeezed / squeezed variance normalized to shot noise (`R±` of Ref A). `B_sq`: one-sided optical bandwidth of the squeezed light [Hz].
- `Δν_d`, `f_d`: peak frequency deviation and frequency of a sinusoidal FM dither tone.

## 2. レーザー雑音のモデル（E12、E15、E16）

### 2.1 周波数雑音と位相

周波数雑音の PSD を次の和で表す。各項は個別にオン・オフできるようにする。

```
S_ν(f) = h0                                   (white FM)
       + h₋₁ / f                              (flicker FM)
       + Σ_i A_i · L(f; f_i, γ_i)              (technical peaks; L = unit-area Lorentzian or Gaussian)
       + S_drift(f)                           (random-walk FM / slow drift, e.g. h₋₂/f²)
```

さらに決定的な成分として、正弦波の FM dither tone δν(t) = Δν_d sin(2π f_d t) と、線形のドリフト δν(t) = ν̇ t を加える。

- 位相の符号: conventions §4 の E ∝ A e^{−iω0 t} では、周波数が δν だけ上がると A ∝ e^{−i2πδν t} になる。
  したがって **φ_L(t) = −2π∫δν dt**。δν を使うブロックはこの負号を必ず入れる。
- 白色 FM の Lorentz 線幅: Δν = π h0（conventions §7）。文献が両側 PSD S^(2)_ν を使う場合は h0 = 2 S^(2)_ν、Δν = 2π S^(2)_ν。
  dBc/Hz の ℒ(f) から変換するときは S_φ = 2·10^{ℒ/10}。
- 白色 FM 以外を含む線幅は観測時間 T0 によって変わる。β-separation line の近似（Di Domenico et al. 2010）:
  S_ν(f) > 8 ln2·f/π² となる周波数の範囲（下限は 1/T0）で S_ν を積分して A とすると、FWHM ≈ √(8 ln2 · A)。
- 一部のメーカー（NKT/Hamamatsu の Koheras）は、位相雑音を「アーム長差 1 m の干渉計で測った値」[dB(rad/√Hz)/m] で示している。
  低周波（f ≪ c/(n_g·1 m)）では √S_ν = 10^{X/20} / (2π n_g · 1 m / c) [Hz/√Hz] で変換できる。

### 2.2 Cobrite DX1（Ref D の光源、D02）[manual]

出典: ID Photonics「Tunable Laser Product Line Data Sheet」（CoBrite シリーズ、2026-07-22 版、p. 5 の表 "Tunable Laser Specification"）、
同社の製品ページ（線幅の比較図）、「CoBriteDX1 User's Guide」V1.24（§5.5.4.5、§5.5.4.7、SCPI の SBS コマンド）。URL と閲覧日は §14。

| 項目 | N タイプ（NC/NL/NX） | SC タイプ | GC タイプ | 単位 |
|---|---|---|---|---|
| 線幅 FWHM | < 100 | 80 typ. | < 100 | kHz |
| instantaneous 線幅（3.5 µs） | 25 typ. | < 100（Pout < 16 dBm）、< 150（それ以外と読める。表の配置からの読み取りで UNVERIFIED） | 25 typ. | kHz |
| 平均 RIN | −145（10 MHz–22 GHz、11 dBm） | −140（100 kHz–20 MHz）、−150（20 MHz–1 GHz） | −145（10 MHz–22 GHz、11 dBm） | dB/Hz |
| 周波数確度（寿命 / 1 時間） | ±2.5 / 0.3 | ±1.5 / 0.3 | ±2.5 / 0.3 | GHz |
| 微調整の分解能 / 範囲 | 1 MHz / ±6 GHz | 10 MHz / ±10 GHz | 1 MHz / ±6 GHz | — |
| 出力パワー（C バンド） | 10.0–16.0 | 8.8–17.8 | 9.5–15.5 | dBm |
| 出力パワーの安定度 | ±0.03 dB / 1 h、±0.05 dB / 24 h（typ.、温度が安定なとき） | 同左 | 同左 | — |
| 出力ファイバ | PM（PANDA）、PER > 18 dB、25 dB typ. | 同左 | 同左 | — |

製品ページの比較図では、典型的な線幅は LX・NL・NX・NC・GC が 0–25 kHz、SC が 0–100 kHz。

マニュアルの記述:
- **Dither（§5.5.4.5）**: 「レーザー周波数を設定値に保つため、キャリア周波数に dither tone をかけている」。オフにできるのは特別版だけで、
  オフにすると安定化ループが止まり、周波数やパワーがドリフトし、マルチモードになることもある。
  **dither の周波数と振幅はマニュアルにもデータシートにも書かれていない。**
- **SBS 抑圧（§5.5.4.7、`SBS:FREQ`）**: 振幅を 0–1 GHz（100 MHz 刻み）で設定できる。変調の周波数は書かれていない。

このノートでの解釈（[AE]）:
- 「instantaneous 3.5 µs」の線幅は、観測時間 3.5 µs、すなわち約 1/(3.5 µs) ≈ 290 kHz より上の雑音で決まる線幅とみなす。
  これを白色 FM の床とし、**h0 = 25 kHz/π ≈ 8.0×10³ Hz²/Hz**（√h0 ≈ 89 Hz/√Hz、N・GC）とする。
  SC タイプは上限の 100 kHz で h0 ≈ 3.2×10⁴ Hz²/Hz。
- FWHM < 100 kHz と 25 kHz の差は、低周波の雑音（flicker FM、dither、ドリフト）による。β-separation の近似で FWHM = 100 kHz になる flicker の係数は、
  h₋₁ ≈ 4.3×10⁸ Hz²（T0 = 1 ms）または 2.7×10⁸ Hz²（T0 = 10 ms）。データシートには T0 が書かれていないので、**h₋₁ は 0 から 5×10⁸ Hz² まで掃引する**。
  標準値は 3×10⁸ Hz² とする（S_ν(1 kHz) ≈ 3×10⁵ Hz²/Hz）。
- 100 kHz より下の RIN はメーカーの値がない。D23 により、ここでは −140 dB/Hz の平坦な値を仮に使い、PROCESS_LOG に記録するよう調整役に依頼する。
  Ref D への影響は §4.2 のとおり小さい。
- 1 時間の周波数確度 0.3 GHz は確度の上限であり、典型的なドリフトではない。ドリフトの速さはパラメータ ν̇ とし、0.3 GHz/h を上限として掃引する。

### 2.3 Ref A の光源（Koheras Harmonik ＋ EDFA ＋ SHG）

- 本文（§3.1）には「NKT Photonics, Koheras harmonik」とだけ書かれている。図 2 を画像で確認したところ、ブロックは「Laser + EDFA + SHG」で、
  その後に BPF、カプラ（LO と probe に分岐）、LO 経路に EOM とフィルタキャビティがある（サーベイ §6-9 の疑問は、EDFA があることで解消）。
- 型番が書かれていないので、雑音スペクトルは分からない。参考として、同じシリーズで 770–785 nm を出す Koheras HARMONIK HP H78 の製品ページには、
  線幅 < 0.4 kHz、RIN のピーク約 0.7 MHz（ピークで < −97 dBc/Hz）、10 MHz で < −132 dBc/Hz とある（Hamamatsu の製品ページ、URL は §14）。
  **これが Ref A の装置と同じ型である保証はない（UNVERIFIED）。** Ref A を再現するときの標準値は「白色 FM の線幅 0.2 kHz」とし、
  RIN のピーク（0.7 MHz）は技術的ピークの項で表す。いずれも代用として PROCESS_LOG に記録する。
- Ref A §4.1 の結びには「より精密な位相ロックのためには、高フィネスの光キャビティを入れて、シードレーザーから出る光の雑音を減らすべき」とある。
  Ref A ではレーザー雑音が θ = 9 mrad の残りの一因と考えられていることになる。

### 2.4 Ref B の光源（E16）

- Ti:Sapphire（Coherent MBR110、860 nm）。本文 §3:「変調周波数は、89 kHz と 177 kHz 付近にある元々のレーザー雑音の影響を受けないよう十分高くした。
  295 kHz ではレーザー雑音がショット雑音レベルとみなせるので、そこに設定した。」
- **雑音が周波数雑音か強度雑音かは書かれていない。高さと幅も書かれていない。** 177 ≈ 2 × 89 kHz なので、基本波と 2 次高調波の可能性があるが、これは推測である。
- モデル: 技術的ピーク 2 本（f_1 = 89 kHz、f_2 = 177 kHz）を、S_ν と S_RIN のどちらに入れるかを選べるパラメータとして持つ。
  高さと幅はパラメータとし、標準値は「295 kHz でショット雑音以下になる」という条件だけを満たすように選ぶ。

## 3. 遅延差のある干渉計による位相差雑音（E13、E14）

### 3.1 一般式

conventions §9 に従い、検出点（または GBS 入力）でのビーム b の位相は φ_b(t) = φ_L(t − τ_b) + Σ_seg ψ_seg である。
同じレーザーから出た 2 本のビームの位相差 Δφ_ab = φ_a − φ_b のうち、レーザーから来る部分は

```
Δφ_ab^(L)(t) = φ_L(t − τ_a) − φ_L(t − τ_b)
S_Δφ(f)      = |e^{−i2πfτ_a} − e^{−i2πfτ_b}|² S_φ(f) = 4 sin²(π f τ) S_ν(f) / f²,   τ = τ_a − τ_b
```

（conventions §7 と同じ式。）

- 低周波（f ≪ 1/(πτ)）: S_Δφ ≈ (2πτ)² S_ν(f)。位相差は周波数ゆらぎに比例する: Δφ^(L) ≈ +2π τ δν(t)。
  符号: φ_L = −2π∫δν なので、Δφ = φ_L(t−τ_a) − φ_L(t−τ_b) = +2πδν(τ_a − τ_b)。長い経路 a ほど、周波数が上がると位相が進む。
  これは ψ = n ω L/c が ω とともに増えることと整合する（伝搬位相は名目値に固定し、周波数による変化は放出時刻の項に入れる。二重に数えないこと）。
- 高周波（f ≥ 1/τ）: 4 sin² は 0 から 4 の間で振動し、f = m/τ では 0 になる。
- 白色 FM だけのとき、全帯域の分散は **σ² = 2π² h0 τ = 2π Δν τ**（τ ≪ 1/Δν のとき。厳密には ⟨e^{iΔφ}⟩ = exp(−πΔν|τ|)）。
- 帯域 [f_1, f_2] に限ったとき（f_2 ≪ 1/τ）: σ² ≈ (2πτ)² [h0 (f_2 − f_1) + h₋₁ ln(f_2/f_1)]。
- ロックがある場合: ロックの誤差信号に入るのがこの Δφ そのものなら、残差は |S_lock(f)|² S_Δφ(f)。
  squeezed 光や LO に対する影響は、どのビームが誤差信号に入るかで決まるので、経路ごとの位相モデルで計算する（D11 の「仮定せず計算する」）。
- 仮定: レーザーの位相雑音は定常でガウス分布に従う。ファイバの分散による周波数ごとの遅延の違いは無視する（τ には群遅延を使う）。
  遅延 τ はゆっくりとしか変わらない（熱による変化は T08 の領域）。

### 3.2 (a) Ref D の腕の長さの差 ΔL

**ΔL の見積もり。**
- `914-915B-Specifications.pdf` [manual]: 915B は 2 台（Serial 8669-1、8669-2）。いずれも 38 巻、ファイバは PM15-U25D、
  ドライバの全レンジでの伸び（光路長として測定）は 16 µm（8669-1）と 17 µm（8669-2）、FC/APC で slow 軸にキー合わせ。
  **巻いたファイバの長さと巻き枠の寸法は書かれていない。** Ref D でどちらの個体を使っているかも分からない（§13）。
- Evanescent Optics の製品掲載（photonics.com Buyers' Guide「Model 915B Piezo-electric Fiber Stretcher」）: 「9 m のファイバで最大 16 µm の伸び、
  DC–10 kHz、100 kHz で 0 に落ちる」「楕円形の枠の側面にピエゾ板を付け、直線部分だけでファイバを伸ばす（PM ファイバの曲げによる結合を避けるため）」。
  **ページを直接取得すると HTTP 403 になり、検索エンジンの抜粋でしか確認できていない（部分的に UNVERIFIED）。**
- 整合性の確認 [AE]: 9 m / 38 巻 = 0.237 m/巻。掲載にある「Length 4.0 in.」（約 10 cm）を楕円の直線部分の長さとみなすと、
  周長は 2 × 10 cm + π × d で、d ≈ 1.2 cm になる。ありうる寸法なので、9 m は矛盾しない。
- 915B のない腕は「かなり短い」（D02）。カプラ PNH1550R5A2 のピグテールが両腕にあり、相殺する。915B のピグテールの長さは分からない。
- **提案**: ΔL の名目値 10 m（巻いた 9 m ＋ 915B のピグテール分）、掃引範囲 5–15 m。τ = n_g ΔL / c = 49.0 ns（ΔL = 10 m）。
  MZI の FSR は c/(n_g ΔL) = 20.4 MHz。

**位相差雑音（ロックなし、Ref D の DSC1 帯域 DC–100 kHz を意識した帯域）** [AE、h0 は §2.2 の解釈]:

| ΔL [m] | τ [ns] | √S_Δφ 低周波（25 kHz）[rad/√Hz] | σ 10 mHz–10 kHz（25 kHz）[mrad] | σ 10 mHz–100 kHz（25 kHz / 100 kHz）[mrad] | σ 全帯域（25 kHz / 100 kHz）[mrad] | 0.3 GHz のドリフト [rad] |
|---|---|---|---|---|---|---|
| 5 | 24.5 | 1.4×10⁻⁵ | 1.4 | 4.3 / 8.7 | 62 / 124 | 46 |
| 9 | 44.1 | 2.5×10⁻⁵ | 2.5 | 7.8 / 15.6 | 83 / 166 | 83 |
| 10 | 49.0 | 2.7×10⁻⁵ | 2.7 | 8.7 / 17.4 | 88 / 175 | 92 |
| 12 | 58.8 | 3.3×10⁻⁵ | 3.3 | 10.4 / 20.8 | 96 / 192 | 111 |
| 15 | 73.5 | 4.1×10⁻⁵ | 4.1 | 13.0 / 26.0 | 107 / 215 | 138 |

- flicker FM を加える（h₋₁ = 5×10⁸ Hz²、ΔL = 10 m）と、10 Hz–100 kHz で σ ≈ 21 mrad、1 kHz–100 kHz で 15 mrad。
  白色の 8.7 mrad より大きい。低周波側はロックで抑えられるので、残差はループの形しだい（T15・T16）。
- **周波数の変化に対する感度**: dΔφ/dν = 2πτ = 3.1×10⁻⁷ rad/Hz（ΔL = 10 m）。915B の全レンジ 16 µm（光路長）は
  2π × 16 µm / 1.5488 µm = **65 rad**（17 µm の個体では 69 rad）。これはレーザー周波数の約 ±100 MHz（全幅 211 MHz）に相当する。
  DX1 の「1 時間の確度 0.3 GHz」が実際のドリフトとして出ると、レンジを超えてロックが外れる。
  実際には温度ドリフト（T08、ΔL = 10 m の腕で数百 rad/K、サーベイの概算）のほうが大きい可能性が高い。
- **dither tone**: 振幅 Δν_d、周波数 f_d の FM は、MZI の位相を振幅 2Δν_d|sin(πf_dτ)|/f_d ≈ 2πτΔν_d で変調する。
  例えば Δν_d = 1 MHz なら 0.31 rad、10 MHz なら 3.1 rad（f_d ≪ 1/τ）。**DX1 の dither の値は分からないので、このモデルで最大の不確かさになりうる。**
  ループ帯域より上にあれば誤差信号に正弦波が残り、帯域内ならストレッチャが追従する。
- **SBS 抑圧**: オンにすると振幅 0.1–1 GHz の周波数変調がかかり、MZI の位相を 31–310 rad 変調する。**Ref D ではオフであることが前提**（設定のフラグとして持つ）。
- RIN の寄与は §4.2（無視できる）。

### 3.3 (b) 将来の構成: squeezed 光と LO_k の遅延の合わせ残り δL_k（指標 (i)）

- τ_δ = n_g δL_k / c。白色 FM だけなら σ_i² = 2πΔν τ_δ。指標 (i) に σ_t を割り当てると、許される合わせ残りは **δL_k ≤ σ_t² c / (2π Δν n_g)**。

| Δν | σ_t = 3 mrad | 10 mrad | 30 mrad |
|---|---|---|---|
| 25 kHz（DX1 N/GC の instantaneous） | 1.2 cm | 13 cm | 1.2 m |
| 100 kHz | 0.29 cm | 3.3 cm | 29 cm |

- パルス光源（D07、100 fs–10 ps）では、時間モードの重なり（T03）のために δL_k はパルス長より十分短く（サブ mm）なければならない。
  したがってパルスの場合は、この要求は自動的に満たされる。CW の LO 配送（Ref C 型）では、この表がそのまま要求になる。
- δL_k の熱による変化（T08 の群遅延の熱係数）は τ_δ をゆっくり変えるだけなので、ここでは静的なパラメータとして扱う。
- Ref C（squeezed 光と LO が別々の 1 km ファイバ）: ファイバ長の差は論文に書かれていない。δL をパラメータとし、名目値は 0 とする（代用として記録）。
  Ref C の光源の線幅も書かれていない（調べた範囲では見つからなかった。T02 の担当範囲）。

### 3.4 (c) 時間多重によるポート間の相関の崩れ（指標 (ii)、D31-3）

**定式化。** squeezer のポンプはレーザーの SHG なので、φ_p = 2φ_L(t_e) + φ_SHG（t_e はポンプが squeezer にいる時刻、φ_SHG は定数）。
conventions §4 から θ_sq = φ_p/2 + π/2 = φ_L(t_e) + 定数。ポート k のパルスは、ポート 0 より τ_k = k T_rep だけ早く生まれ、ファイバ L_k で τ_k だけ遅らされて、
同じ時刻 t に GBS に入る。したがって

```
Θ_k(t) = [φ_L(t − τ_k) − φ_L(t)] + [ψ_k(t) − ψ_0(t)] + const
         \_______ D31-3 ______/   \___ path phases (D31-1, -2) ___/
S_Θ^(3)(f) = 4 sin²(π f k T_rep) S_ν(f)/f²,  σ² = 2πΔν k T_rep   (white FM, no lock)
```

**ロックによる補償（D31-3 への修正の提案）。** Liu 型ロック（D30）では、probe を GBS 側で分けてポート k を逆向きに通し、上流で参照光と干渉させる。
静的な周波数のずれ δν（φ_L = −2πδν t）で確かめる。probe は時刻 t' − τ_k − τ_p にレーザーを出て、時刻 t' に上流の参照点に着く。
参照光は時刻 t' − τ_r にレーザーを出ている。ロックの条件 φ_L(t' − τ_k − τ_p) + ψ_k = φ_L(t' − τ_r) + ψ_r + c_k から

```
ψ_k = −2πδν (τ_k + τ_p − τ_r) + ψ_r + c_k   ⇒   ψ_k − ψ_0 = −2πδν τ_k + const
Θ_k = [−2πδν(t − τ_k) + 2πδν t] + [−2πδν τ_k] + const = const
```

ロックが書き込む補正は、放出時刻の差による項とちょうど打ち消し合う。動的には、この打ち消しは周波数成分のうち f ≪ f_BW かつ f ≪ 1/τ_k の部分で成り立つ。
ループの外の成分は打ち消されない。おおよそ

```
S_Θ^(3),locked(f) ≈ |S_lock(f)|² · 4 sin²(π f τ_k) S_ν(f)/f²   (f ≪ 1/τ_k for the cancellation; above f_BW S_lock → 1)
```

正確な値は、経路ごとの位相モデルで probe の往路と復路の時刻を区別して計算する（解析式は検証テスト V8 の目安）。
**成り立つ条件**: probe が同じレーザーから、同じ（またはずれが既知で一定の）光周波数で作られていること。
Liu 2026 のフィルタ用 MZI のように probe が別波長（1064 nm）なら、補償されない。AOM で f_AOM だけずらした probe では、
2π f_AOM τ_k の一定のずれと、その熱による変化（f_AOM = 80 MHz、群遅延の熱係数 39 ps/(km·K)（サーベイの UNVERIFIED 値）、63 m で約 1.2 mrad/K、[AE]）が残るだけ。

帰結:
- 白色 FM の分散は、ほとんどが f ~ 1/(2τ_k) 付近にある。f_BW より下にあるのは割合で 2τ_k f_BW（τ = 310 ns、f_BW = 1 kHz で 0.06%）だけ。
  したがって **線幅（h0）への要求はロックで緩まない**。D31-3 の「レーザー線幅への要求として出す」は白色 FM については正しい。
- 一方、ドリフトと低周波の flicker FM はロックで打ち消される。もし打ち消されなかったら、0.3 GHz のドリフトでポート 31（τ = 310 ns）の位相は 584 rad、
  1 MHz の周波数のずれでも 1.95 rad 動く。**probe を同じレーザーから作ることは要求として明記すべき。**
- **指標 (i) は逆にずれる**: LO_k の遅延を squeezed 光の経路に合わせた場合（D09）、ロックが入れる補正 ψ_k は LO_k にはかからない。
  そのため Δ_k は +2πδν(τ_k + τ_p − τ_r) だけずれる（conventions の符号で導出、§3.1 と同じ手順）。
  BHD でポート k を評価するときは、レーザーのドリフトがそのまま squeezing 角のドリフトとして見える
  （τ_k + τ_p − τ_r = 310 ns なら、1 MHz あたり 1.95 rad）。指標 (i) と (ii) は、レーザー周波数のドリフトに同時には強くできない。
  LO_k の作り方（D09 Open）の判断材料になる。

**線幅への要求（白色 FM、ロックなしと同じ）。** 最悪のポート（k = N − 1 = 31）で σ ≤ σ_t にするには

```
Δν_max = σ_t² f_rep / (2π (N − 1)),    h0,max = Δν_max/π
```

ポートで平均した RMS（k = 0…31）で決めるなら、平均の k は (N−1)/2 なので Δν_max は 2 倍になる（σ は 1/√2 倍）。

| f_rep | τ_31 | L_31（n_g = 1.468） | Δν_max（σ_t = 3 / 10 / 30 mrad） | σ_31（Δν = 25 kHz / 1 kHz / 1 Hz） |
|---|---|---|---|---|
| 10 MHz | 3.1 µs | 633 m | 0.46 / 5.1 / 46 Hz | 0.70 / 0.14 / 0.0044 rad |
| 20 MHz（Liu 2026 と同じ） | 1.55 µs | 317 m | 0.92 / 10 / 92 Hz | 0.49 / 0.099 / 0.0031 rad |
| 100 MHz | 310 ns | 63 m | 4.6 / 51 / 462 Hz | 0.22 / 0.044 / 0.0014 rad |
| 1 GHz | 31 ns | 6.3 m | 46 / 513 / 4620 Hz | 0.070 / 0.014 / 0.00044 rad |

- 参考: Ref A の R±（L = 8%、640 mW）に厳密式（conventions §6）を当てると、σ = 10 mrad で −10.05 dB、30 mrad で −7.2 dB、
  44 mrad で −5.1 dB、88 mrad で −0.2 dB、0.22 rad で +7.3 dB。この写像は T04 の担当で、ここでは目安。
- flicker FM（h₋₁ = 5×10⁸ Hz²、τ = 310 ns）のうち、ループ帯域より上（10 kHz–50 MHz）の部分だけで 0.096 rad。
  線幅の数字だけでなく、S_ν(f) の形で要求を書く必要がある: ∫_{f_BW}^{f_rep/2} 4 sin²(πfτ_31) S_ν(f)/f² df ≤ σ_t²。
- 積分の上限: GBS はショットごとの測定なので、f_rep/2 までのすべての成分が効く。白色 FM では sin² の因子のために f ~ 1/τ より上の寄与が小さく、
  上限を f_rep/2 で切った誤差は小さい。

**モード同期レーザーの櫛の線幅と CW の線幅。**
- パルスを CW から切り出す場合（Liu 2026: シードに EOM、ポンプに AOM でゲート）、φ_L は CW レーザーの位相そのもの。上の表がそのまま CW の線幅への要求になる。
- モード同期レーザーの場合、k T_rep 離れた 2 つのパルスのキャリアの位相差は、名目値 2πk f_ceo/f_rep（mod 2π、決定的）と、
  キャリア付近の櫛の線 ν_n（n ≈ ν0/f_rep。100 MHz なら約 2×10⁶）の周波数雑音による部分からなる。
  後者は上と同じ式で、S_ν を S_ν,n に置き換えればよい。δν_n = δf_ceo + n δf_rep なので、f_ceo のゆらぎと、n 倍された f_rep のゆらぎの両方が入る。
- 「elastic tape」（固定点）モデル（Newbury & Swann 2007）では、櫛の線の雑音は固定点からの距離に比例して増える。
  繰り返しを RF 基準にロックすると、RF の位相雑音が n² 倍されてキャリア付近の線幅は広くなる。櫛の線を狭線幅の CW レーザー（キャビティ安定化）にロックすると、
  キャリア付近の線幅はその CW レーザーの線幅になる。**自由発振の櫛の具体的な線幅は機種によって大きく違うので、ここでは数値を挙げない（UNVERIFIED）。**
- 要求（f_rep = 100 MHz、σ_t = 10 mrad で約 50 Hz）は、自由発振の半導体 ECL（DX1 25 kHz）や、おそらく自由発振のファイバ櫛でも満たせない。
  キャビティ安定化した CW レーザー（Hz レベル）か、それに櫛の線をロックした構成が必要になる。
- 他の対策の候補: f_rep を上げる（τ_31 が短くなる）。レーザーの位相を遅延自己ヘテロダインで測り、スロットに同期した速い EOM でフィードフォワードする
  （D10 の比較用構成「スロット同期の EOM」と組み合わせる）。いずれも Phase 6 の比較対象。
- 時間モードの重なり（タイミングジッタ、δL_k による包絡線のずれ）は T03 の担当。

### 3.5 シミュレーションでの扱い（遅延がサンプル間隔より短い場合）

時間ループのサンプル間隔 dt（DSC1 で 1 µs、他の構成でも 0.1–1 µs 程度）は、τ（Ref D で 49 ns、ポート間で 31 ns–3 µs）と同じか、より長い。
遅延バッファで φ_L(t − τ) を整数サンプルで取ると、τ が丸められてしまう。

1. レーザーの雑音源ブロックは δν(t) をループのサンプル周波数で作る（f_Nyq = 1/(2dt) までに帯域制限）。φ_L(t) = −2π∫δν dt も保持する。
2. 経路ごとの位相は φ_L を **端数の遅延で線形補間して** 取る。τ ≪ dt では、これは Δφ ≈ +2πτ δν(t) と一致する（帯域内では正しい）。
3. f_Nyq より上の寄与 σ²_oob = ∫_{f_Nyq}^{∞} 4 sin²(πfτ) S_ν/f² df はオフライン（O）の表で計算し、T04 の squeezing への写像に分散として渡す。
   白色 FM で f_Nyq ≪ 1/τ なら、σ²_oob ≈ 2π²h0τ − 4π²τ²h0 f_Nyq。ロックは f_Nyq より上の成分を見ないので、この扱いで矛盾しない。
4. 同じ φ_L(t) をすべての経路が共有する（雑音源は 1 つ）。経路の間の相関は、各経路の τ_b から自動的に出てくる。

## 4. 強度雑音（RIN）の伝わり方（E15）

### 4.1 BHD（Ref A・F）: CMRR を通る漏れ

LO の RIN は同相信号として両方の PD に入る。差電流への漏れを ε とすると、ショット雑音に対する過剰雑音の比は

```
N_RIN / N_shot = ε² ℜ P_LO S_RIN(f) / (2e)
```

- Ref A の LO 15.4 mW、η = 0.99（ℜ = 1.24 A/W）、RIN −145 dB/Hz（DX1 N の値で代用）: ε = 0.01（CMRR 40 dB）で 1.9×10⁻²（+0.08 dB）。
  ε = 0.056（PDB450C の保証値 25 dB）なら 0.59（+2.0 dB）。RIN が −140 dB/Hz なら、それぞれ 5.9×10⁻² と 1.9。
- LO と ASE のビートも、RIN S_RIN,ASE = 4 S_ASE/P_LO として同じ式で漏れる（§5.2）。
- CMRR の値と周波数依存性は T13 の担当。このノートは S_RIN(f) を渡す。

### 4.2 Ref D の DC side-fringe lock

差信号 ∝ P(ε + V cos φ)（サーベイ E75、T15）。φ = π/2 の零点では、RIN は ε を通してだけ入る: δφ_eq = ε (δP/P) / V。
- RIN −140 dB/Hz（100 kHz より下は仮の値）、ε = 0.056、V = 1 → 5.6×10⁻⁹ rad/√Hz。レーザーの周波数雑音による床（2.7×10⁻⁵ rad/√Hz）より 4 桁小さい。
- パワーの安定度 ±0.03 dB/h（0.7%）× ε = 0.056 → 0.39 mrad のオフセット。小さい。T15 で確認すること。

### 4.3 probe とコヒーレント光の古典雑音（サーベイの表にない、追加を提案）

BHD の signal 側に残るコヒーレント光（Ref A の残留 probe、Ref B の同軸 probe、F の Liu 型 probe の漏れ）は、
測定周波数 f_s でのレーザー雑音を squeezed 光と同じモードに持ち込む。パワー P_pr の光が測定の直交位相に加える過剰雑音（ショット雑音に対する比）は

```
amplitude quadrature:  P_pr S_RIN(f_s) / (2 hν0)
phase quadrature:      2 P_pr S_Δφ(f_s) / (hν0),   S_Δφ: probe−LO phase-difference noise (§3.1, incl. filter cavity §6)
```

（どちらも η を掛ける。ショット雑音レベルの RIN は 2hν0/P、位相雑音は hν0/(2P)。）
- 例 [AE]: P_pr = 1 µW、RIN −145 dB/Hz → 0.012。位相雑音（DX1 の h0、f_s = 3 MHz、LO のレーザー位相雑音と打ち消し合わない場合）→ 1.4×10⁴。
  probe と LO の遅延差が 1 m なら、打ち消しの因子 4 sin²(πf_sτ) = 8.5×10⁻³ が掛かって 117、10 cm なら 1.2。
  LO だけがフィルタキャビティを通る Ref A の構成では、Δν_FC/2 より上で打ち消しが成り立たない（§6）。
- 帰結: 位相雑音の床が高いレーザー（ECL）では、BHD に届く残留 probe は pW–nW 程度に抑える必要がある。
  Ref A の probe の残留パワーは本文に書かれていない（1 MHz のピークとして見えているだけ）。
- 担当の提案: この効果を新しい行（E-new「残留 probe の古典雑音」、段 O ＋ T、担当 T05 ＋ T12/T14）として survey に加える（§12）。

### 4.4 SPM による RIN → 位相

パルスの LO や probe がファイバで受ける SPM の位相 φ_NL ∝ P は、RIN を位相雑音に変える（δφ = φ_NL δP/P）。
係数は T06 の担当。このノートは S_RIN(f) を渡すところまでを受け持つ。

## 5. ASE と LO のフィルタキャビティ（E17）

### 5.1 ASE と広帯域 squeezed 光のビート（Ref A §3.1 の機構）

LO 側の入力に、強い LO と弱い広帯域 ASE（位相がランダムな熱的な場）があるとする。差電流には「ASE × signal 側の場」の項が入る。
RF 周波数 f の雑音は、光周波数で f 離れたすべての組から集まる。ASE は位相がランダムなので、signal 側の各モードの直交位相の分散を位相で平均した値
(V₊ + V₋)/2（真空なら 1）を拾う（2 モード squeezed 状態の片側サイドバンドの光子数 sinh²r から、2 sinh²r + 1 = cosh 2r = (V₊ + V₋)/2）。
ショット雑音の基準（signal 側を遮って測る）は P_LO + P_ASE,eff に比例する。したがって [AE、このノートでの導出]:

```
V_meas/V_SN = (1 − ε_ASE) V_θ + ε_ASE (V₊ + V₋)/2,     ε_ASE = P_ASE,eff / (P_LO + P_ASE,eff)
P_ASE,eff  = ∫ S_ASE(ν) |T_FC(ν)|² |T_BPF(ν)|² dν   over |ν − ν0| ≤ B_sq (LO mode and polarization only)
```

- V₊ と V₋ は、ASE がある光周波数での squeezed 光の値（帯域 B_sq の中でほぼ一定とみなす。周波数による違いは T01 の表から取る）。
- Ref A の R±（−10.58 dB / +20.56 dB、θ = 0）で: ε_ASE = 10⁻⁴ → −10.31 dB、10⁻³ → −8.41 dB、10⁻² → −1.84 dB。
  **反 squeezing が大きいほど ASE に敏感になる**ので、10 dB 級では ε_ASE ≲ 10⁻⁴ が必要。
- 例 [AE、OSNR は仮定]: EDFA の出力で OSNR 40 dB（0.1 nm = 12.5 GHz、両偏光）→ 1 偏光で S_ASE/P = 4×10⁻¹⁵ /Hz。
  B_sq = 2.5 THz（Kashiwazaki 2020 の帯域を代用）でフィルタがないと ε_ASE = 0.02 になり、squeezing は消える。
- フィルタキャビティ（FSR ごとに透過のピーク、Lorentz 形）: 帯域 B_sq ≫ FSR では、各ピークが (π/2)Δν_FC 分の ASE を通す。
  そのため P_ASE,eff はフィルタなしの **π/(2F)** 倍になる。F = 100 / 1000 / 10⁴ で ε_ASE = 3.1×10⁻⁴ / 3.1×10⁻⁵ / 3.1×10⁻⁶。
  Ref A の BPF（図 2、帯域は不明）があれば、積分範囲はさらに狭くなる。
- **Ref A のフィルタキャビティの線幅（FILTER_CAVITY_LW）と ASE レベル（ASE_LEVEL）は論文に書かれていない。リファレンス値がないので、ユーザに確認が必要**（D23 の規則、§13）。
  それまではパラメータとして掃引し、標準値は「ASE なし」とする。

### 5.2 LO と ASE のビート（CMRR を通る）と、signal 側の ASE

- LO と ASE のビート（Olsson 1989 の signal–spontaneous ビート）は同相なので、ε を通してだけ漏れる。§4.1 の式に S_RIN,ASE = 4S_ASE/P_LO を入れると、
  N/N_shot = 2ε² η n_ASE,LO。n_ASE,LO = S_ASE/(hν0) は、LO ビーム中の 1 モードあたりの ASE 光子数。
  例: P_LO = 15.4 mW、S_ASE/P = 4×10⁻¹⁵ /Hz → n = 480、ε = 0.01 で 0.096。この漏れも、フィルタキャビティが Δν_FC/2 より上の周波数で取り除く。
- probe に乗った ASE は signal 側のモードに入り、熱的な雑音として 2n_ASE,pr を加える。n_ASE,pr = (S_ASE/P)·P_pr/(hν0)。
  例: P_pr = 1 µW（BHD に届く量）→ n = 0.031、2n = 0.062（ショット雑音比）。−10.6 dB の squeezing が −8.3 dB 程度に落ちる。
  Ref A では、図 2 の BPF のあとに probe の経路にフィルタキャビティがない。BHD に届く probe のパワーが小さいことが前提になっていると考えられる（推測）。
- F（Ref D）では EDFA がないので、この節の効果は無視できる。将来の構成で EDFA を使う場合に有効にする。

## 6. フィルタキャビティ / MCC のロック残差（E18）

conventions（E ∝ A e^{−iω0 t}、伝搬で A → A e^{+iβL}）では、共振付近の Fabry–Perot キャビティの振幅透過は

```
T_FC(δ) = T0 / (1 − i 2δ/Δν_FC),    arg T_FC = +arctan(2δ/Δν_FC),    |T_FC|² = T0² / (1 + (2δ/Δν_FC)²)
```

（往復位相 φ_rt = 2π(ν/FSR) を 1/(1 − r e^{iφ_rt}) に入れ、r ≈ 1 で展開。定数の伝搬位相は除く。符号は実装で conventions と照合すること。）

- **離調による LO の位相**: δφ_LO ≈ 2δ/Δν_FC。キャビティのロック残差 δ_rms から σ_LO = 2δ_rms/Δν_FC。
  Δν_FC = 0.1 / 1 / 10 MHz、δ_rms = 1 kHz で 20 / 2 / 0.2 mrad。Ref A の θ = 9 mrad と比べられる大きさなので、Δν_FC が小さいほどロックを厳しくする必要がある。
- **離調による強度**: |T|² の変化は δ の 2 次で小さい。ただし δ の平均値がずれていると 1 次になり、LO の RIN になる。
- **レーザー位相雑音の低域通過**: 透過光の位相ゆらぎは H_FC(f) = 1/(1 + i 2f/Δν_FC) で伝わる。キャビティを通らない経路（probe、squeezer のポンプ）との位相差の雑音は
  S_Δφ = |1 − H_FC(f)|² S_φ(f)。低周波では (2/Δν_FC)² S_ν で、§3.1 の τ を τ_eff = 1/(πΔν_FC) に置き換えたものと同じ（Δν_FC = 1 MHz で τ_eff = 318 ns）。
  **フィルタキャビティは LO にとって、長さ c·τ_eff/n_g の遅延線と同じ働きをする。** Ref A のロック 2（BHD の probe–LO ビート）は、
  帯域内ではこれを LO の EOM で補償する。帯域外は残る。
- キャビティ長の雑音: δ = −ν0 δL_cav/L_cav（長さが伸びると共振周波数が下がる）。キャビティ長の雑音と、キャビティのロック（PDH など、Black 2001）の残差は T16 の担当と重なる。
  このノートは、δ(t) を入力として LO の位相と振幅を出すブロックを定める。
- Ref B の MCC（LO 経路）も同じブロックで表す。線幅は論文に書かれていない（パラメータ）。

## 7. 仮定と適用範囲

- レーザーの雑音は定常でガウス分布に従う（dither・SBS・ドリフトは決定的な項として別に持つ）。モードホップとマルチモード化は扱わない（DX1 のマニュアルは dither をオフにするとマルチモードになりうると書いている）。
- 位相差の式は線形（Δφ を位相として足し合わせる）。縞の非線形性（cos）は誤差信号のブロック（T14・T15）が扱う。
- 遅延は群遅延で表し、分散（GVD）による周波数ごとの遅延の違いは無視する。パルスでは T03・T06 が扱う。
- §3.4 のロックによる補償は、probe と squeezed 光のポンプが同じレーザーから同じ光周波数で作られ、ロックの帯域が 1/τ_k より十分低い場合に成り立つ。
- §5.1 の ASE のモデルは、ASE が熱的（位相がランダム）で、LO と同じ空間モード・偏光の成分だけがビートし、ASE と squeezed 光が統計的に独立という仮定に基づく。導出は [AE] で、文献で確認していない。
- DX1 の h0 と h₋₁ は、データシートの線幅からの解釈（§2.2）であり、測定値ではない。

## 8. 周波数帯とスケーリング

| 効果 | 周波数帯 | ファイバ長・遅延 | パワー | 温度 | パルス幅 |
|---|---|---|---|---|---|
| E12 周波数雑音 | DC–数十 MHz（白色の床は 1/τ まで効く） | — | 依存しない（線幅は出力によって少し変わる。SC タイプは 16 dBm 以上で広がる） | レーザー本体の温度でドリフト | パルスでは櫛の線の雑音 |
| E13 遅延差の干渉計 | DC – 1/τ（それより上は sin² で振動） | 低周波では S_Δφ ∝ τ²、白色 FM の全分散は ∝ τ | — | τ の熱変化はゆっくり（T08） | δL_k の要求はパルスではモード重なりで決まる |
| E14 ポート間 | f_BW – f_rep/2 が線幅への要求になる（それより下はロックで補償） | σ² ∝ k T_rep = k/f_rep | — | — | 櫛の構成による |
| E15 RIN | DC–GHz（DX1 の値は 100 kHz 以上だけ） | — | CMRR の漏れは ∝ P_LO（ショット雑音比）、probe の古典雑音は ∝ P_pr | — | SPM を通る分は ∝ P_peak ∝ 1/τ_p |
| E16 技術的ピーク | 89・177 kHz（Ref B）、0.7 MHz（Koheras の RIN ピーク） | 周波数雑音のピークなら ∝ τ² | — | — | — |
| E17 ASE | 光の帯域（THz）。RF では検出帯域全体 | — | ε_ASE ∝ S_ASE/P_LO（ASE とビートする経路） | — | — |
| E18 フィルタキャビティ | DC – Δν_FC | τ_eff = 1/(πΔν_FC) | — | キャビティ長の熱ドリフト | — |

## 9. 概算のまとめ（自分たちのパラメータ）

| 対象 | 条件 | 値 | 印 |
|---|---|---|---|
| Ref D 位相差の床 | ΔL = 10 m、h0 = 25 kHz/π | 2.7×10⁻⁵ rad/√Hz、100 kHz まで 8.7 mrad | AE（h0 は §2.2 の解釈） |
| Ref D、flicker あり | h₋₁ = 5×10⁸ Hz²、10 Hz–100 kHz | 約 21 mrad | AE |
| Ref D、周波数ドリフト | 0.3 GHz、ΔL = 10 m | 92 rad（915B のレンジ 65 rad を超える） | manual ＋ AE |
| Ref D、dither | Δν_d = 1 MHz（仮定） | 0.31 rad の正弦波 | AE（dither の値は非公開） |
| (b) δL_k | Δν = 25 kHz、σ_t = 10 mrad | ≤ 13 cm | AE |
| (c) ポート 31 | f_rep = 100 MHz、DX1 25 kHz | 0.22 rad | AE |
| (c) 要求 | f_rep = 100 MHz、σ_t = 10 mrad | Δν ≤ 51 Hz | AE |
| ASE（ビート） | ε_ASE = 10⁻³、Ref A の R± | −10.58 → −8.41 dB | AE |
| RIN（CMRR） | 15.4 mW、−145 dB/Hz、ε = 0.01 | ショット雑音の 1.9% | AE |
| フィルタキャビティ | Δν_FC = 1 MHz、δ_rms = 1 kHz | 2 mrad | AE（Δν_FC は仮定） |

## 10. 検証テスト（実装が通すこと）

数値は `est.py` 相当の計算（このノートの担当が scratchpad で実行。コードはコミットしていない）による。許容誤差は統計誤差の目安。

- **V1 遅延の伝達関数**: 白色 FM（Δν = 25 kHz）、dt = 1 ns、τ = 49 dt = 49 ns で φ_L を作り、Δφ = φ_L(t) − φ_L(t − τ) の Welch PSD が
  4 sin²(πfτ)·h0/f² に ±10% で一致する（零点 f = m/τ の近くを除く、f < 0.4/dt）。分散は 2πΔντ = 7.70×10⁻³ rad²（σ = 87.7 mrad）に ±5% で一致する。
- **V2 コヒーレンス**: 同じ条件で τ = 10 µs のとき ⟨cos Δφ⟩ = exp(−πΔντ) = 0.456 ± 0.01。
- **V3 共通の経路**: τ_a = τ_b なら、レーザーから来る Δφ は厳密に 0（|Δφ| < 10⁻¹²）。prompt §7-11 の co-propagation テストと同じ考え方。
- **V4 サンプル間隔より短い遅延（§3.5）**: dt = 1 µs、τ = 49 ns。帯域内（f < f_Nyq/2）の PSD が (2πτ)²S_ν に ±10% で一致し、
  帯域内の分散とオフラインの σ²_oob の和が解析値 2πΔντ に ±2% で一致する。
- **V5 Ref D の周波数ステップ**: ΔL = 10 m（τ = 49.0 ns）で δν = +1 MHz のステップ。915B の腕を a とする。
  ロックなし: Δφ がちょうど +0.3077 rad 変わる。ロックあり: 定常状態で、ストレッチャの位相が −0.3077 rad 変わる。
- **V6 Ref D のロックなしの床**: 音響・熱の雑音源をオフにすると、誤差信号の位相換算 PSD は (2πτ)²h0 = 7.5×10⁻¹⁰ rad²/Hz（ΔL = 10 m、Δν = 25 kHz）に ±10% で一致する（prompt §7-7 の後半）。
- **V7 dither tone**: Δν_d = 1 MHz、f_d = 1 kHz、τ = 49 ns で、f_d に振幅 0.308 rad の線が出る。SBS 抑圧のフラグを 1 GHz で入れると振幅は 308 rad になる。設定の確認に使う。
- **V8 ポート間（ロックなし）**: Δν = 1 kHz、T_rep = 10 ns、k = 31 で Var(Θ_31^(3)) = 1.95×10⁻³ rad²（σ = 44.1 mrad）、±5%。
- **V9 Liu 型ロックの静的な補償（2 ポートの設定、D13）**: f_rep = 100 MHz、τ_p = τ_r。δν = +1 MHz のステップのあと定常状態で
  Θ_1 の変化は 0（< 10⁻⁶ rad）、Δ_1 の変化は +2π·10⁶·10 ns = +0.0628 rad（LO_1 の遅延を合わせた場合）。
  probe の周波数を Δf = 80 MHz ずらす設定では、Θ_1 に一定のずれ 2πΔf τ_1 が出て、δν のステップでは変化しない。
  probe を別のレーザーにする設定では補償されず、Θ_1 は放出時刻の項 φ_L(t − τ_1) − φ_L(t) の変化 +2πδντ_1 = +0.0628 rad だけ変わる。
  （符号は conventions の φ_L = −2π∫δν から導出。実装者は自分で導出し直して照合すること。）
- **V10 ロックありのポート間の残差**: 線形のループで、白色 FM による Θ_k の PSD が |S_lock(f)|²·4 sin²(πfτ_k)S_ν/f² に、f < 1/(10τ_k) の範囲で ±10% で一致する。
- **V11 フィルタキャビティ**: 静的な離調 δ = 0.1Δν_FC で arg T = +0.1974 rad、|T|²/T0² = 0.9615。
  位相雑音の伝達 |H_FC(f)|² = 1/(1 + (2f/Δν_FC)²) が f = Δν_FC/2 で 0.5 になる。
- **V12 ASE**: S_ASE = 0 なら結果は厳密に変わらない。ε_ASE = 10⁻³、Ref A の R± で −8.41 dB（±0.01 dB）。
  フィネス F のキャビティを通すと、帯域の広い ASE の P_ASE,eff が π/(2F) 倍になる（B_sq ≫ FSR のとき、±2%）。
- **V13 RIN（CMRR）**: P_LO = 15.4 mW、η = 0.99、RIN −145 dB/Hz、ε = 0.01 で、過剰雑音がショット雑音の 1.88×10⁻² になる。
- **V14 Ref B のピーク**（傾向のテスト）: 89・177 kHz にピークを入れ、dither を 295 kHz にしたときと 89 kHz にしたときで、
  ロック 2 の残留位相雑音が後者で大きくなる。前者はピークがないときと（統計誤差の範囲で）同じになる。

## 11. シミュレータがこのトピックから必要とするもの

| ブロック / 表 | 段 | 入力 | 出力 | 切り替えのフラグ |
|---|---|---|---|---|
| `laser_noise_source`（光源ごとに 1 つ、全経路で共有） | T | h0、h₋₁、ピークの一覧 {f_i, γ_i, A_i, 周波数雑音か RIN か}、ドリフト ν̇、dither {Δν_d, f_d}、SBS {振幅, 周波数}、RIN のモデル、seed | δν(t)、φ_L(t)（−2π∫δν）、δP/P(t) を state に書く | `laser_fm_white`、`laser_fm_flicker`、`laser_peaks`、`laser_drift`、`laser_dither`、`laser_sbs`、`laser_rin` |
| 経路ごとの放出時刻（遅延バッファ、端数の遅延の補間） | T | 各ビーム経路の群遅延 τ_b（区間の L・n_g の和、D11）、ループの dt | φ_L(t − τ_b)（ビームの位相の光源の項） | — （経路の定義に含まれる） |
| `laser_oob_variance` | O | S_ν(f) のモデル、経路の組と τ、f_Nyq、ロックの帯域 | ループで表せない帯域の σ²_oob（T04 へ） | `laser_fm_*` と同じ |
| `filter_cavity` | O ＋ T | Δν_FC、FSR、T0、δ(t)（キャビティのロック残差、T16 の出力でもよい） | LO の位相 arg T(δ)、透過 |T|²、レーザー位相の低域通過 H_FC | `filter_cavity`、`filter_cavity_lock_residual` |
| `ase_equivalent_loss` | O | S_ASE(ν)（OSNR から）、BPF・フィルタキャビティの透過、B_sq、V±（T01 の表）、P_LO、P_pr | ε_ASE、LO の RIN_ASE、probe の n_ASE | `ase_lo`、`ase_probe` |
| `coherent_beam_classical_noise`（§4.3） | O ＋ M | P_pr、S_RIN(f_s)、S_Δφ(f_s)、η | 測定周波数でのショット雑音比の過剰雑音 | `probe_classical_noise` |

- 出力の指標: (i) と (ii) の両方（D12）。(ii) は D31 の 3 つの成分に分けて出す。成分 (3) は、ロックで補償された分と補償されない分に分けて出すことを提案する（§3.4）。
- 掃引するパラメータ: ΔL（5–15 m）、h0（DX1 の N/GC/SC、Koheras 相当）、h₋₁（0–5×10⁸ Hz²）、f_rep、ν̇、Δν_d・f_d（不明なので広く）、Δν_FC、S_ASE。
- 機器プリセット: `cobrite_dx1_N`、`cobrite_dx1_SC`、`cobrite_dx1_GC`（§2.2 の表）、`koheras_harmonik_refA`（UNVERIFIED）、`tisa_mbr110_refB`（ピーク 2 本）。

## 12. DECISIONS・conventions・survey との食い違いと提案（調整役向け）

1. **D31-3 の文言**: 「ロックでは補正できない」は、Liu 型で probe を同じレーザーから同じ周波数で作る場合は、ループ帯域より下では成り立たない（§3.4）。
   提案: 「ループ帯域より上の成分は補正できないので線幅への要求として出す。帯域内の成分（ドリフト、flicker）は、probe を同じレーザーから作れば補償される。
   これは計算で確かめる」と書き換える。D30 に「probe は squeezer のポンプと同じレーザーから、同じ（またはずれが既知の）光周波数で作る」という要求を加える。
2. **指標 (i) と (ii) のトレードオフ**（§3.4）: D09（LO_k の遅延を合わせる）と D30 の組み合わせでは、レーザー周波数のドリフトが指標 (i) に出る。D09 の Open の判断材料として記録を提案する。
3. **survey への追加**: E-new「残留 probe（コヒーレント光）の古典雑音（RIN・レーザー位相雑音・ASE）が測定周波数で BHD に入る」（§4.3、§5.2）。段 O ＋ T、担当 T05 ＋ T12/T14。
4. **survey の書誌の訂正**: Ref A は Optics Express 34, 7958 (2026)、DOI 10.1364/OE.585323（Crossref で確認。PDF には書かれていない）。
   Ref C は Crossref の issued が 2025（survey は 2026。オンライン公開日と号の年の違いかもしれない。UNVERIFIED）。
5. **survey §6-9**: Ref A 図 2 で「Laser + EDFA + SHG」→ BPF を確認した。EDFA はある。
6. **survey E13 の概算の訂正**: 「ΔL = 10 m・Δν = 100 kHz で全帯域 0.18 rad、100 kHz 以下 17 mrad」は、このノートの計算（175 mrad、17.4 mrad）と一致する。
   ただし DX1 の N/GC タイプなら、白色の床は 25 kHz 相当で、その半分になる。
7. conventions に追加を提案: §7 に「φ_L(t) = −2π∫δν dt（E ∝ A e^{−iω0t} から）」「ℒ(f) = S_φ/2」「Δν = π h0 は片側 PSD」の 3 行。
   §9 に「ビームの位相の光源の項は、群遅延で決まる放出時刻で評価し、周波数による伝搬位相の変化と二重に数えない」の 1 行。

## 13. ユーザの決定・情報が必要な点（推測しない）

1. **Ref D の DX1 のタイプ**（型番 CBDX1-1-**XX**-FA の XX: NC / GC / SC など）。線幅の標準値が 25 kHz か 100 kHz かを決める。
2. **DX1 の dither** が標準版（常時オン）か特別版（オフにできる）か。オンの場合、dither の振幅と周波数をメーカーに問い合わせるかどうか。Ref D で最大の不確かさになりうる（§3.2）。
3. **SBS 抑圧** が Ref D の運用でオフになっているか（`SBS:STAT?` の値）。
4. **915B の個体**: Ref D で使っているのは 8669-1（16 µm）か 8669-2（17 µm）か。915B のピグテールの長さ。
5. **ΔL の測り方**: DX1 の微調整（1 MHz 刻み、±6 GHz）で MZI の縞を数えれば ΔL を 1% 未満で決められる（FSR は 10 m で 20.4 MHz）。これは D04 の「実測値を使わない」に当たるか（装置の寸法の測定とみなせるか）。
6. **Ref A のフィルタキャビティの線幅と ASE レベル**: 論文に値がない（D23 では「リファレンス値がなければ止めて聞く」）。著者への問い合わせ、文献 [23] の確認、または「パラメータとして掃引する」でよいか。
7. **100 kHz より下の DX1 の RIN**: メーカー値がない。仮に −140 dB/Hz の平坦な値で続けてよいか（Ref D への影響は小さい、§4.2）。
8. **将来の光源の位相雑音の要求の割り当て**: σ_t（D31-3 に割り当てる位相誤差）をいくつにするか。T04 の「目標 × L → 必要な θ_rms」表と合わせて決める必要がある。f_rep = 100 MHz で 10 mrad を割り当てると、約 50 Hz の線幅（白色 FM）が必要になる。

## 14. 文献

一次文献（Crossref で書誌と DOI を確認済み、2026-10-28）:
- K. Hirota, T. Kashiwazaki, G. Ha, T. Yamashima, P. Jaturaphagorn, T. Suzuki, et al., "Generation of 10-dB squeezed light from a broadband waveguide optical parametric amplifier with improved phase locking method," Opt. Express 34, 7958 (2026). doi:10.1364/OE.585323（本文は `refs/` の PDF で読んだ）
- Y. Takeno, M. Yukawa, H. Yonezawa, A. Furusawa, "Observation of −9 dB quadrature squeezing with improvement of phase stability in homodyne measurement," Opt. Express 15, 4321 (2007). doi:10.1364/OE.15.004321（本文は `refs/`）
- H.-L. Liu, H. Su, Y.-H. Deng, S.-Q. Gong, Y.-C. Gu, et al., "Gaussian boson sampling with 1,024 squeezed states in 8,176 modes," Nature 653, 687–692 (2026). doi:10.1038/s41586-026-10523-6（本文は `refs/`。Methods「Phase stabilization of the system」）
- S. Verclas, C. Darsow-Fromm, B. Tohermes, R. Schnabel, "Fiber distribution of phase-stabilized GHz-bandwidth squeezed vacuum states of light between two buildings," Optica Quantum 4, 1. doi:10.1364/OPTICAQ.567418（このノートでは本文を参照していない）
- G. Di Domenico, S. Schilt, P. Thomann, "Simple approach to the relation between laser frequency noise and laser line shape," Appl. Opt. 49, 4801 (2010). doi:10.1364/AO.49.004801（β-separation line。式の形は記憶による。本文は未読）
- L. E. Richter, H. I. Mandelberg, M. S. Kruger, P. A. McGrath, "Linewidth determination from self-heterodyne measurements with subcoherence delay times," IEEE J. Quantum Electron. 22, 2070–2074 (1986). doi:10.1109/JQE.1986.1072909（遅延差のある干渉計の位相差雑音。本文は未読）
- T. Okoshi, K. Kikuchi, A. Nakayama, "Novel method for high resolution measurement of laser output spectrum," Electron. Lett. 16, 630–631 (1980). doi:10.1049/el:19800437（遅延自己ヘテロダイン法。本文は未読）
- L. B. Mercer, "1/f frequency noise effects on self-heterodyne linewidth measurements," J. Lightwave Technol. 9, 485–493 (1991). doi:10.1109/50.76663（flicker FM の扱い。本文は未読）
- D. S. Elliott, R. Roy, S. J. Smith, "Extracavity laser band-shape and bandwidth modification," Phys. Rev. A 26, 12–18 (1982). doi:10.1103/PhysRevA.26.12（周波数雑音と線形の関係。本文は未読）
- N. R. Newbury, W. C. Swann, "Low-noise fiber-laser frequency combs (Invited)," J. Opt. Soc. Am. B 24, 1756 (2007). doi:10.1364/JOSAB.24.001756（elastic tape / 固定点モデル。本文は未読）
- E. D. Black, "An introduction to Pound–Drever–Hall laser frequency stabilization," Am. J. Phys. 69, 79–87 (2001). doi:10.1119/1.1286663（キャビティの透過と PDH。本文は未読）
- N. A. Olsson, "Lightwave systems with optical amplifiers," J. Lightwave Technol. 7, 1071–1082 (1989). doi:10.1109/50.29634（signal–spontaneous ビート雑音。本文は未読）
- L.-S. Ma, P. Jungner, J. Ye, J. L. Hall, "Delivering the same optical frequency at two places: accurate cancellation of phase noise introduced by an optical fiber or other time-varying path," Opt. Lett. 19, 1777 (1994). doi:10.1364/OL.19.001777（LO 配送のロック、D12。T08 と共有）
- T. Kashiwazaki, N. Takanashi, T. Yamashima, T. Kazama, K. Enbutsu, R. Kasahara, et al., "Continuous-wave 6-dB-squeezed light with 2.5-THz-bandwidth from single-mode PPLN waveguide," APL Photonics 5, 036104 (2020). doi:10.1063/1.5142437（B_sq の代用値として題名の帯域だけを使った）

メーカー資料（D24）:
- Evanescent Optics, "Model 914-2 Piezo Controller / 915B Piezo Assembly"（ハナムラオプティクスの書類のスキャン）、`refs/manuals/914-915B-Specifications.pdf`。
- ID Photonics GmbH, "Tunable Laser Product Line Data Sheet"（CoBrite シリーズ、2026-07-22 版）、
  https://id-photonics.com/download/cobrite-series-datasheet/ （データシート一覧 https://id-photonics.com/resources/datasheets/ から取得、閲覧 2026-10-28）。
- ID Photonics GmbH, CoBrite 製品ページ（線幅の比較図）、https://id-photonics.com/products-solutions/cobrite-tunable-laser/ （閲覧 2026-10-28）。
- ID Photonics GmbH, "CoBriteDX1 User's Guide" V1.24。販売代理店（Laser Lab Source）のミラー:
  https://shop.laserdiodesource.com/wp-content/uploads/2022/08/Tunable_Laser_Manual_CoBrite_DX1-1659720345.pdf （閲覧 2026-10-28。メーカーのサイトに直接の置き場所は見つからなかった）。
- Evanescent Optics, "Model 915B Piezo-electric Fiber Stretcher"（photonics.com Buyers' Guide の掲載）、
  https://www.photonics.com/Buyers-Guide/Products/Model-915B-Piezo-electric-Fiber-Stretcher/psp9424 （閲覧 2026-10-28。直接の取得は HTTP 403 で、内容は検索結果の抜粋でのみ確認。**部分的に UNVERIFIED**）。
- Hamamatsu（旧 NKT Photonics）, "Koheras HARMONIK HP H78" 製品ページ、https://lf.hamamatsu.com/products/single-frequency-fiber-lasers/koheras-harmonik-h78/ （閲覧 2026-10-28。Ref A の装置と同じ型かは UNVERIFIED）。
