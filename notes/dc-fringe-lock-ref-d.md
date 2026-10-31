# T15 `dc-fringe-lock-ref-d`: Ref D の DC side-fringe lock の解析モデル

> **訂正（2026-10-31、DECISIONS D65-3）:** §2.11・§6.1・§7 test 10 の「DAC のリミットサイクルは ±1/2 LSB」は誤り。±LSB/2 で抑えられるのは rms だけで、遅延と動特性があるとリミットサイクルは f₁₈₀ の振動になる（ptp ≈ (4/π)kq|P(i2πf₁₈₀)|）。K_i・K_d の単位はマニュアル p.20（D53）、915B の共振の範囲は 5–60 kHz・Q 3–100（D51-11）。


- Phase: 1（deep dive）／ 日付: 2026-10-28
- 担当範囲（survey §4）: E39 の D 部分、E75、E76、DSC1 プリセット（D05・D06）、Ref D 全体の解析モデル（prompt §7-7、D04）
- 依存: T05（Cobrite DX1 の周波数雑音・RIN）、T16（915B の機械共振、914 の小信号応答、DSC1 の離散化の詳細）、T08（環境雑音）、T09（PER）。**本ノートではこれらをプレースホルダで置き、そう明記する。**
- 入力資料: `docs/squeezing_sim_prompt.md`、`DECISIONS.md`（D01–D33）、`notes/conventions.md`、`notes/survey.md`、
  `refs/manuals/thorlabs-dsc1.pdf`、`refs/manuals/thorlabs-pdb.pdf`、`refs/manuals/914-915B-Specifications.pdf`。
- D04 に従い、Ref D の観測値（`blind/`）は一切読んでいない。値はメーカー資料と物理だけから取った。

### 頁の引き方

- DSC1: 「DSC1 p.N」は**マニュアルの印刷頁**（Rev. B, 2025-11-03, DOC-103674）。PDF の頁は印刷頁 + 5（例: 仕様表 p.2 = PDF p.7）。
- PDB: 「PDB p.N」は**印刷頁**（© 2019）。PDF の頁は印刷頁 + 2（例: PDB450x 個別仕様 p.19 = PDF p.21）。
- 914/915B: 「EO p.N」は PDF の頁（全 3 頁、スキャン画像。テキスト層なし。画像から読んだ）。

---

## 1. 系と記号

光路（D02）: Cobrite DX1（CW、λ0 = 1548.8 nm）→ MN9605C アッテネータ → Optoquest PC →（以降 PM、FC/APC）→ Opneti 偏光子 → PNH1550R5A2（2×2、カプラ 1）→
腕 1: 915B（38 巻）／腕 2: 短い PM ファイバ → PNH1550R5A2（カプラ 2）→ PDB450C（INPUT+ / INPUT−）→ RF OUTPUT → DSC1 入力 → DSC1 出力 → Model 914 入力 → 915B。

### Symbols (English)

| Symbol | Definition | Unit |
|---|---|---|
| `P_in` | optical power into coupler 1 (after attenuator + polarizer), 3.33 mW (D02) | W |
| `T_opt` | power transmission from coupler 1 input to the two PDs (couplers' excess loss, 915B insertion loss, connectors). **Unknown → parameter, default 1 (upper bound)** | – |
| `P_det = T_opt·P_in` | total power reaching the two photodiodes, `P_1 + P_2` at coupler 2 input | W |
| `P_1, P_2` | powers in arm 1 (915B) and arm 2 at coupler 2 input | W |
| `κ_1, κ_2` | power coupling ratios of coupler 1 and coupler 2 (nominal 0.5) | – |
| `ℜ` | PD responsivity at λ0; `ℜ_± = ℜ(1 ± ε_R)` for INPUT± | A/W |
| `ε_R` | responsivity imbalance of the two PDs | – |
| `G` | PDB450C transimpedance gain, 10³ V/A into high-Z load (PDB p.19 note 1) | V/A |
| `V_0 = G·ℜ·P_det` | full-scale difference signal | V |
| `ε` | normalized DC imbalance of the difference signal (Eq. 3) | – |
| `𝒱` | effective fringe visibility of the difference signal (Eq. 3) | – |
| `μ = |μ|e^{iφ_μ}` | complex overlap of the two arm fields at coupler 2 (polarization × temporal coherence) | – |
| `x_j` | fast-axis power fraction in arm j (`x = 10^{−PER/10}` for a single cross-coupling point) | – |
| `Δχ` | differential birefringent phase (slow − fast) accumulated between the cross-coupling point and coupler 2, arm 2 minus arm 1 | rad |
| `ψ_1, ψ_2` | segment phases of arm 1 / arm 2 (conventions §4, `ψ = n_eff ω0 L/c`) | rad |
| `φ_D = φ_2 − φ_1` | **Ref D lock variable**: phase of the arm-2 field minus phase of the arm-1 (915B) field at coupler 2 (includes laser phase at the two emission times) | rad |
| `v(t)` | PDB450C RF OUTPUT voltage = DSC1 input (process variable PV) | V |
| `V_off` | electrical DC offset at the DSC1 input (PDB offset + DSC1 input offset) | V |
| `S` | DSC1 setpoint voltage (D02: initial S = 0, i.e. difference-signal zero) | V |
| `e_DSC1 = S − v` | DSC1 error (DSC1 p.18, §5.4.1). **Sign opposite to conventions §9**, see §10 | V |
| `φ_lock` | value of `φ_D` at which `v = S` (lock point) | rad |
| `K_d = ∂v/∂φ_D |_{φ_lock}` | discriminator slope | V/rad |
| `δp(t)` | relative power fluctuation, `P_det → P_det(1+δp)`; `S_RIN(f)` its one-sided PSD | 1/Hz |
| `u(t)` | DSC1 output voltage (control variable CV), ±4 V | V |
| `K_p, K_i, K_d^{PID}` | DSC1 gains: dimensionless, 1/s, s (DSC1 p.12, p.17). `K_d^{PID}` written with superscript to avoid clash with `K_d` | – |
| `h` | DSC1 controller sample interval, 2 µs (PID update 500 kHz, DSC1 p.2) | s |
| `τ_DSC1` | DSC1 effective pure delay (fitted, §2.7) | s |
| `k_act` | quasi-static actuator gain, optical phase per volt at the 914 input | rad/V |
| `A_914(s)` | 914 normalized transfer function (=1 at DC) | – |
| `H_915B(s)` | 915B normalized mechanical response (=1 at DC) | – |
| `L(s)` | open-loop transfer function (closed loop = 1/(1+L)) | – |
| `f_c` | unity-gain (crossover) frequency of `L` | Hz |
| `ΔL = L_1 − L_2` | arm length difference (fiber), swept (D02) | m |
| `τ = n_g ΔL / c` | differential delay of the MZI | s |
| `n_g` | group index of PM15-U25D, placeholder 1.468 (UNVERIFIED, T09) | – |
| `S_ν(f)`, `Δν` | laser frequency-noise PSD, Lorentzian linewidth (`S_ν = Δν/π` for white FM; conventions §7) | Hz²/Hz, Hz |
| `S_env(f)` | free-running environmental PSD of `φ_D` (T08 placeholder) | rad²/Hz |

---

## 2. 支配方程式

### 2.1 干渉信号

カプラ 2 の入力での各腕の場を `E_j = √P_j · ê_j · e^{iφ_j}`（`ê_j` は単位偏光ベクトル）とし、カプラ 2 を損失のない 2×2 とする。
conventions §5 のビームスプリッタ（`b_1 = (a_s + a_LO)/√2`, `b_2 = (a_s − a_LO)/√2`）を一般の結合比 `κ_2` に拡げ、`b_1` を INPUT+、`b_2` を INPUT− に入れると

```
P_± = κ_2^{(±)} P_1 + (1 − κ_2^{(±)}) P_2 ± 2√(κ_2(1−κ_2) P_1 P_2) · |μ| cos(φ_D + φ_μ)          (1)
      with κ_2^{(+)} = κ_2, κ_2^{(−)} = 1 − κ_2,   P_+ + P_− = P_1 + P_2 = P_det
μ = ⟨ê_1, ê_2⟩ · γ_coh                                                                         (2)
```

PDB450C の RF OUTPUT は `v = G(ℜ_+ P_+ − ℜ_− P_−)`（PDB p.8、§3.3）。これに電気オフセットを足すと

```
v = V_0 [ ε + 𝒱 cos(φ_D + φ_μ) ] + V_off                                                       (3)
V_0 = G ℜ P_det
ε   = (2κ_2 − 1)(P_1 − P_2)/P_det + ε_R
𝒱   = 2√(κ_2(1−κ_2)) · [2√(P_1 P_2)/(P_1 + P_2)] · |μ|
```

- 理想（`κ_1 = κ_2 = 1/2`、損失なし、`ε_R = 0`、`μ = 1`）では `ε = 0`、`𝒱 = 1`、`v = V_0 cos φ_D`。conventions §5 の `⟨n̂_1 − n̂_2⟩ = 2|α_LO||α_pr| cos(φ_LO − φ_pr)` と同じ形。
- **カプラ 2 の結合比だけのずれは `ε` に入るが、腕のパワー比 `P_1/P_2` のずれと掛け算になるので 2 次**。1 次で効くのは PD の感度差 `ε_R` と、`κ_2 ≠ 1/2` かつ `P_1 ≠ P_2`（915B の挿入損失、`κ_1` のずれ）の組み合わせ。
- `𝒱` は `κ_2` と `P_1/P_2` のずれに対して 2 次で鈍感（例: `P_1/P_2 = 0.8` で `2√(P_1P_2)/(P_1+P_2) = 0.994`）。`𝒱` を大きく下げるのは `|μ|`（偏波と可干渉性）。

### 2.2 conventions との対応（提案）

- 腕 2（915B なし）を conventions §5 の「LO」、腕 1（915B、アクチュエータ側）を「probe（pr）」に対応させる: `φ_D = φ_LO − φ_pr`。
  将来系でアクチュエータを probe・squeezed 側に置く（D10）のと向きを揃えるための選択。conventions §9 の言葉では「ロック = 線形結合 `φ_2 − φ_1` を設定値に保つ」。
- ストレッチャで腕 1 を伸ばすと `ψ_1` が増え（conventions §4）、`φ_D` は**減る**: `∂φ_D/∂u = −k_act`（`k_act > 0` と定義。実機の PZT 極性が逆なら符号パラメータで反転）。
- 設定値の初期値 `S = 0` に対応する `φ_lock` は、理想系で `+π/2`（D02、conventions §5）。`−π/2` も零点だが、傾きの符号が逆。DSC1 は P・I・D の符号を一括で反転させる（DSC1 p.19）ので、どちらにもロックできる。
  シミュレータでは `loop_sign ∈ {+1, −1}`（PDB のどちらの出力を INPUT+ に入れたか × PZT 極性 × DSC1 のゲインの符号）を明示的なパラメータにし、`φ_lock = +π/2` で負帰還になる値を既定とする。

### 2.3 ロック点と傾き

ロック条件 `v = S` から（厳密）

```
cos(φ_lock + φ_μ) = c_lock ≡ (S − V_off)/(V_0 𝒱) − ε/𝒱                                         (4)
φ_lock = arccos(c_lock) − φ_μ          (the +π/2 branch)
K_d = ∂v/∂φ_D |_{φ_lock} = −V_0 𝒱 sin(φ_lock + φ_μ) = −V_0 𝒱 √(1 − c_lock²)                      (5)
```

理想（`ε = V_off = S = 0`、`𝒱 = 1`）では `φ_lock = π/2`、`K_d = −V_0`。`|c_lock| < 1` でなければロック点がない（S をフリンジの外に置いた場合）。

### 2.4 ロック点と傾きの摂動（D06、E75、E76）

`c_lock` を 1 次で変分すると、ループが `v = S` を保つことで真の `φ_D` が次だけずれる（`φ_lock` の近傍、`sin ≈ 1`）:

```
δφ_lock = [ δV_off/(V_0 𝒱) + δε/𝒱 + c_lock·δ𝒱/𝒱 + (c_lock + ε/𝒱)·δp ] / sin(φ_lock + φ_μ) − δφ_μ      (6)
        with (c_lock + ε/𝒱) = (S − V_off)/(V_0 𝒱)
```

各項の意味:

1. **DC オフセットのドリフト** `δV_off`: `δφ = δV_off/|K_d|`。PDB の RF OUTPUT オフセット < ±15 mV（PDB p.19）、工場で 23 °C に零調整、暖機 約 5 分で小さくドリフト（PDB p.8）。DSC1 の入力オフセットは資料に記載なし（プレースホルダ）。
2. **光学的な不平衡のドリフト** `δε`（PD 感度の温度変化、915B の挿入損失の変化）: `δφ = δε/𝒱`。
3. **可視度の変化** `δ𝒱`: `c_lock = 0`（真の直交点）では 1 次で効かない。`ε` や `V_off` で直交点から外れているときだけ `c_lock·δ𝒱/𝒱` で効く。
4. **パワー変化・RIN** `δp`: 係数は `(S − V_off)/(V_0𝒱) = (S − V_off)/|K_d|`。
   **すなわち強度雑音がロック点を動かすのは、設定値が電気的なゼロ（`V_off`）からずれている分だけ**。`ε` は式 (4) でロック点に吸収されるので、光学的な不平衡そのものは RIN の結合係数に入らない（`v` のうち `δp` に比例する部分は `V_0(ε + 𝒱 cos) δp = (S − V_off) δp`）。
   これは survey E75 の「φ = π/2 では強度雑音は ε を通してだけ入る」を、ロック点が `ε` で動くことまで含めて正確にしたもの。
5. **偏波の重なりの位相** `δφ_μ`: ロック点を直接 1 対 1 で動かす（§2.5）。

`v` の雑音として見ると、RIN は `δv = (S − V_off)·δp` の**センサ雑音**（ループで補正できない）。閉ループで真の位相に移る量は §2.8 の `T(s)` を掛けたもの。

傾きの変化（ループ利得の変化、E76）: `δK_d/K_d = δp + δ𝒱/𝒱 + O(c_lock δc)`。パワーが 10% 下がれば `f_c` も（積分器型のループでは）約 10% 下がる。

### 2.5 偏波と PER による `μ`（E39 の D 部分、T09 と分担）

偏光子の後、カプラ 1 の入力で fast 軸に `x` のパワーが入る（偏光子の消光比、コネクタのキーずれ角 θ_k で `x ≈ sin²θ_k`）と、両腕に同じ `x` が乗る。2 軸の成分はそれぞれ自分どうしで干渉するので

```
μ = (1 − x) + x · e^{iΔχ},       Δχ = (2π/λ0) · B · (L_2 − L_1)  (+ terms from connectors/coupler leads)   (7)
|μ| ∈ [1 − 2x, 1],               |φ_μ| ≤ arcsin( x/(1 − x) )
```

腕ごとに別の点で `x_1`、`x_2` が生じる（915B の巻きの応力、腕の中のコネクタ）場合は `μ = √((1−x_1)(1−x_2)) + √(x_1 x_2) e^{iΔχ}`、`|φ_μ| ≲ √(x_1 x_2)`。

- `Δχ` は 2 軸の位相差の温度ドリフト（D11 の概算 2.4 rad/(m·K)、未検証、T09 が確認）で回る。ΔL = 9 m なら 約 22 rad/K。すなわち **0.1 K 程度の温度変化でロック点が ±x の範囲を 1 周する遅い揺らぎ**になる。
- ストレッチャでの伸びも `B` に比例して `Δχ` をわずかに変える（`δΔχ/δψ ~ B/n ~ 3×10⁻⁴`、[AE]）。無視してよい。
- 偏光子の手前（DX1 → アッテネータ → PC の SMF 部分）の偏波ドリフトは、偏光子で**パワーの変化** `δp` に変わる（§2.4 の 4）。量は T09 の SMF の値を使う。

### 2.6 腕の長さの差によるレーザー位相雑音（E13、T05 と分担）

`φ_D` はレーザーの位相を 2 つの出射時刻で含む: `φ_D ⊃ φ_L(t − τ_2) − φ_L(t − τ_1)`。conventions §7 より

```
S_{φ_D, laser}(f) = 4 sin²(π f τ) · S_ν(f)/f²,      τ = n_g ΔL / c                               (8)
f ≪ 1/τ:   S ≈ (2πτ)² S_ν(f)
white FM:  σ²_total = ∫_0^∞ S df = 2π Δν τ
```

- これは**真の位相の揺らぎ**（外乱）で、ループで抑えられる（センサ雑音ではない）。ただし DSC1 の帯域（DC–100 kHz）より上は抑えられず、検出器の帯域（150 MHz）には現れる。
- `1/τ` より上まで含めた時間平均は、時間ループ段で明示的に扱わないとき `𝒱_eff = 𝒱 · exp(−σ²_HF/2)`（`σ²_HF` はモデル化した帯域の外の分）として可視度に入れる。**明示的にモデル化した帯域を可視度に二重に入れないこと**（`γ_coh` は時間ループ段では 1 とし、帯域外の分だけを入れる）。

### 2.7 ループの各要素（値はマニュアルのみ）

**(a) 検出器 PDB450C**（PDB p.14、p.15、p.19、p.21）
- 利得 `G = 10³ V/A`（高インピーダンス負荷のとき。50 Ω 負荷では 1/2。DSC1 の入力終端は 1 MΩ（DSC1 p.2）なので 10³ V/A を使う）。
- 帯域 DC–150 MHz（10³ V/A）、ループに対しては定数とみなせる。
- 感度: 表は「Max. Responsivity typ. 1.0 A/W」（p.19）。p.15 の典型曲線では 1550 nm 付近で 約 1.03 A/W（図から読んだ値）。既定 `ℜ = 1.0 A/W`。
- CW 飽和 4.5 mW @ 1550 nm、平衡時は 2 入力の**差**がこれ未満であること（p.7、p.8）。RF OUTPUT の最大振幅 ±4.6 V（高 Z、p.14）。
- RF OUTPUT の DC オフセット < ±15 mV（p.19）。CMRR > 25 dB（typ. > 30 dB）（p.19）。
- 最小 NEP 68.6 pW/√Hz（DC–150 MHz、10³ V/A、PDB450C、p.19）。

**(b) DSC1**（DSC1 p.1、p.2、p.5–6、p.12、p.17–20）
- 入力: 16 bit、±4 V、1 MHz サンプル、2 点平均で実効 500 kHz、入力帯域 typ. 120 kHz、入力終端 1 MΩ（p.1、p.2）。
- 出力: 12 bit、±4 V、出力インピーダンス 220 Ω（高 Z 負荷向け。負荷 `R_dev` ではレンジが `R_dev/(R_dev + 220 Ω)` 倍）、出力の LPF で帯域 100 kHz（p.2 注 b・d）。
- システム帯域 DC–100 kHz、入力→出力の −180° 周波数 > 58 kHz（typ. 60 kHz）（p.2）。
- 入力 RMS 雑音 0.3 mV（100 Hz–250 kHz）、平均雑音フロア −120 dB V²/Hz、ピーク −105 dB V²/Hz（p.2）。
- PID 更新 500 kHz（p.2）。誤差 `e = S − v`（p.18）。制御則（p.19）:

```
u[n] = K_p e[n] + h Σ_{m=1}^{n} K_i e[m] + (K_d^{PID}/h)(e[n] − e[n−1])  (+ output offset O, clipped to CV limits)   (9)
C(z) = K_p + K_i h · z/(z − 1) + (K_d^{PID}/h)(1 − z^{−1})
C(s) ≈ K_p + K_i/s + K_d^{PID} s     (continuous, p.17)
```

  - 積分は `K_i` を積分の中に入れる bumpless 形（p.18）。I を無効にすると `K_i = 0`。P・I・D の符号は連動（p.19）。ゲインの設定範囲は絶対値 10⁻⁵–10⁵（p.6）。
  - **`h` の値**: p.19 は「ADC サンプルの間隔」と書くが、ADC は 1 MHz、PID 更新は 500 kHz。p.20 の「Ziegler–Nichols の積分項に 2·10⁻⁶ を掛けて正規化」から `h = 2 µs` とする。
  - **`K_i` の単位の食い違い**: p.12 は「Ki は 1/s」、p.20 は「表の Ki（1/s）に 2·10⁻⁶ を掛けて入れる」。後者なら画面に入れる値は 1 サンプルあたりの無次元量。**どちらかで積分のゲインが 5×10⁵ 倍変わる**。プリセットに `ki_entry: {per_second, per_sample}` を持たせ、既定は p.12・p.17 に従い `per_second`。§11 で確認を依頼。
- **ループ遅延の当てはめ** [AE]: 入力帯域・出力帯域をそれぞれ 1 次の極（120 kHz、100 kHz。次数は資料にない）とし、残りを純遅延とすると、60 kHz で −180° になる遅延は

```
τ_DSC1 = [180° − atan(60/120) − atan(60/100)] / (360° · 60 kHz) = 5.67 µs     (58 kHz: 5.94 µs)
```

  極を入れず純遅延だけなら 8.33 µs（58 kHz: 8.62 µs、survey E82 の値）。既定は「極 2 つ ＋ 5.67 µs」、比較用に「純遅延 8.33 µs」。
  ここで −180° は P のみの設定で、`e = S − v` による符号反転を除いた位相と解釈した（p.2 注 a は定義を書いていない。UNVERIFIED）。
  2 点平均（1 MHz で 0.5 µs）と 500 kHz の ZOH（1 µs）はこの 5.67 µs に含まれる。

**(c) Model 914-2 Piezo Controller**（EO p.1、p.2）
- 入力 ±5 V で 915B のフルスイング（EO p.1）。DSC1 の出力は ±4 V なので、**使えるのはフルレンジの 80%**（914 の入力インピーダンスは資料になく、高 Z と仮定）。
- 出力はプッシュプル（DB-9 の 2・3 ピン、EO p.1、p.3）。
- 周波数応答は EO p.2 Fig. 1 の振幅だけ（位相なし）。**±5 V のフルスイング入力**での圧電電圧。画像から読んだ値（読み取り誤差 約 ±1 V）:

| f [kHz] | 1 | 5 | 7 | 10 | 20 | 27 | 40 | 50 | 100 |
|---|---|---|---|---|---|---|---|---|---|
| 圧電電圧 [V] | 52 | 52 | 52 | 50 | 42 | 38 | 29–30 | 23–24 | 11 |

  - −3 dB（36.8 V）は 約 28–29 kHz。
  - 40–100 kHz では `f × V ≈ 一定`（2π f V = 7.4、7.4、6.9 V/µs、52 V を片振幅とした場合）で、**スルーレート制限の形**。フルスイングの大信号の測定なので、ロック中の小さな補正に対する小信号帯域はこれより広い可能性がある。
  - 既定モデル（プレースホルダ、T16 で確定）: `A_914(s) = 1/(1 + s/(2π·28 kHz))`（小信号）、加えて大信号の制限として Fig. 1 の振幅表（オフライン表）またはスルーレート 約 7 V/µs（52 V が片振幅の場合。両振幅なら半分）。
  - 914 の出力雑音・オフセット・入力インピーダンスは資料にない（T16 のプレースホルダ）。

**(d) 915B ファイバストレッチャ**（EO p.1）
- 2 台分の記載: S/N 8669-1、8669-2。いずれも 38 巻、ファイバ PM15-U25D、FC/APC（slow 軸にキー）。
- 「Fiber stretch (optical path length measurement) for full range of driver voltage」: **8669-1 = 16 µm、8669-2 = 17 µm**。
- 「光路長として測った」とあるので、位相への換算に光弾性係数 ξ_eo は掛けない:

```
k_act = 2π · ΔOPL_full / (λ0 · ΔV_in,full) = 2π · 16 µm / (1548.8 nm · 10 V) = 6.49 rad/V   (8669-1)
                                                                       = 6.90 rad/V   (8669-2)
```

  もし 16 µm が機械的な伸びなら `× n ξ_eo ≈ 1.45 × 0.78 = 1.13` で 7.3 rad/V（conventions §4、ξ_eo は T16）。既定は光路長の解釈。
- DSC1 の ±4 V で使えるレンジ: `8 V × 6.49 rad/V = 51.9 rad`（8.3 フリンジ）、中心から ±26 rad。
- **機械共振・容量・巻線ファイバ長は 914/915B の資料にない。** Photonics.com の製品掲載（検索結果の抜粋のみ。ページ本体は 403 で開けず UNVERIFIED）に「9 m のファイバで最大 16 µm」「1 巻 220 mm」「DC–10 kHz で使え、100 kHz で 0 に落ちる」とある。38 × 0.22 m = 8.4 m。
  共振は 2 次系 `H_915B(s) = ω_r²/(s² + s ω_r/Q + ω_r²)` のプレースホルダとし、`f_r ∈ [10, 60] kHz`、`Q ∈ [5, 30]` を掃引（T16 が値を決めるまで）。

**(e) 光の伝搬遅延**: 915B から PD まで 約 10 m で 約 50 ns。`τ_DSC1` に比べ無視できるが、遅延バッファとしては区間の長さから自動で入る。

### 2.8 開ループ伝達関数と閉ループ

外乱 `φ_n`（環境、レーザー、914 の出力雑音、DAC の量子化）を `φ_D` に、センサ雑音 `v_n`（DSC1 入力雑音、PDB の NEP、ショット雑音、RIN 由来の `(S − V_off) δp`、オフセットドリフト）を `v` に足す。線形化すると

```
L(s) = −loop_sign · K_d · C(s) · H_in(s) H_out(s) · e^{−sτ_DSC1} · (−k_act) · A_914(s) · H_915B(s)
     ≡ K_0 · C(s) · H_in H_out e^{−sτ_DSC1} A_914 H_915B ,     K_0 = |K_d| · k_act                         (10)
δφ_D = 1/(1+L) · φ_n  −  L/(1+L) · v_n/K_d                                                                       (11)
e_DSC1 (in-loop, measurable) = −K_d/(1+L) · φ_n  −  1/(1+L) · v_n
```

`H_in = 1/(1 + s/(2π·120 kHz))`、`H_out = 1/(1 + s/(2π·100 kHz))`。`K_d < 0`（+π/2）、`k_act > 0` で、正の DSC1 ゲインが負帰還になるよう `loop_sign` を取ると `K_0 > 0`。

- `K_0 = 3.33 V/rad × 6.49 rad/V = 21.6`（`T_opt = 1`、`𝒱 = 1`、8669-1）。**DSC1 のゲインは 10⁻² 程度（P）、10³ s⁻¹ 程度（I）で十分**ということになる。
- Ref D は検出器が 1 つなので、実機で測れるのは in-loop の `e_DSC1` だけ。シミュレータは真の `δφ_D`（out-of-loop）と `e_DSC1` の両方を出す。センサ雑音は in-loop 誤差では押し込まれて見えるが、真の位相には `T = L/(1+L)` で乗る。

### 2.9 安定余裕と PID ゲイン（数値、[AE]、§2.7 の既定モデル）

`fs` を 1 Hz–500 kHz にとり、(10) から数値で求めた（スクリプトはノートに含めない。式は上の通り）。915B は剛体（`H_915B = 1`）と共振プレースホルダの両方。

| 構成 | 量 | 値 |
|---|---|---|
| P のみ、915B 剛体、914 = 28 kHz 極 | −180° の周波数 / 発振する `K_p,u` | 40.9 kHz / 0.093 |
| P のみ、914 の極なし（DSC1 だけ） | 同上 | 60.0 kHz / 0.060 |
| P のみ、純遅延 8.33 µs のみ | 同上 | 45.1 kHz / 0.054 |
| P のみ、915B 共振 20 kHz・Q = 10 | 同上 | 19.9 kHz / 0.0059 |
| I のみ、剛体、位相余裕 60° | `K_i` → `f_c`、利得余裕 | 1.75×10³ s⁻¹ → 5.9 kHz、11.5 dB |
| I のみ、剛体、位相余裕 45° | 同上 | 2.73×10³ s⁻¹ → 8.9 kHz、7.7 dB |
| I のみ、利得余裕 6 dB を満たす最大 | 剛体 / 共振 40 kHz Q10 / 20 kHz Q10 / 20 kHz Q30 / 10 kHz Q10 | `K_i` = 3.3×10³ / 2.5×10³ / 1.0×10³ / 0.81×10³ / 0.21×10³ s⁻¹、`f_c` = 10.5 / 8.7 / 3.7 / 2.8 / 0.7 kHz |
| PI（`K_i/K_p = 2π·1 kHz`）、剛体 | `K_p` = 0.002 / 0.01 / 0.03 | `f_c` = 43 / 221 / 851 Hz、利得余裕 33 / 19 / 10 dB |
| 同、共振 20 kHz Q10 | `K_p` = 0.002 / 0.005 / 0.01 | 利得余裕 9.4 / 1.4 / −4.6 dB（0.01 で不安定） |

要点:
- **積分器だけのループの `f_c` の上限は剛体で 約 10 kHz、915B の共振で決まる**（`f_c,max ≈` 共振周波数 / (数 × Q)）。共振の値（T16）が Ref D の予測を最も大きく左右する未知量。
- P ゲインは共振を直接持ち上げるので、PI の `K_p` は共振ありで 10⁻³ 台に制限される。
- 積分器だけのループの近似: `L ≈ f_c/(i f)`（`f ≪ f_c` 付近まで）、`f_c = K_0 K_i/(2π)`。

### 2.10 閉ループの残差（パラメトリックな環境雑音、T08 のプレースホルダ）

```
σ²_φ = ∫_{f_min}^{f_max} [ |1/(1+L)|² (S_env + S_laser + S_act) + |L/(1+L)|² S_v/K_d² ] df                    (12)
S_env(f) = S_0 (f / 1 Hz)^{−α}      (T08 placeholder: S_0 ∝ ΔL for uncorrelated segment noise, α ∈ {0, 1, 2})
```

積分器だけのループ（`|1/(1+L)|² = f²/(f² + f_c²)`）での閉じた形:
- `α = 2`: `σ² = π S_0/(2 f_c)`（f_min → 0、f_max → ∞）。`f_c = 5.9 kHz` で `S_0 = 10⁻⁴ / 10⁻² / 1 rad²·Hz` → `σ = 0.16 / 1.6 / 16 mrad`。
- `α = 1`: `σ² = (S_0/2) ln(1 + f_max²/f_c²)`。
- 温度のランプ（D27 の速さ `r_T` [K/s]）: `r_φ = (dφ/dT)·ΔL·r_T`。積分器 1 段のループでは定常誤差 `r_φ/(2π f_c)`。
  `dφ/dT ≈ (2π/λ0)(dn/dT + n α_L) ≈ 44 rad/(m·K)`（survey E22 [AE]、`dn/dT = 1.0×10⁻⁵ /K`、`α_L = 0.55×10⁻⁶ /K` を仮定、T08 で確認）。

### 2.11 レンジと量子化

- **レンジ**: ±26 rad（中心から）。一様な温度変化に対しては ΔL の分だけ効くので、`ΔT_span = 51.9 rad / (44 rad/(m·K) · ΔL)` = 0.13 K（ΔL = 9 m）。ロック時間 ≈ 26 rad / r_φ = 0.066 K / r_T（ΔL = 9 m、中心から開始）。例: 1 mK/s で 約 66 s。
  CV リミット（DSC1 p.9）に達した後、積分器が止まるか（anti-windup）は資料にない（T16、UNVERIFIED）。
- **ADC**: LSB = 8 V/2¹⁶ = 122 µV → 位相 37 µrad。入力雑音 0.3 mV より小さい。
- **DAC**: LSB = 8 V/2¹² = 1.95 mV → **1 LSB = 12.7 mrad の光位相**（8669-1）。ループの外乱（レーザー雑音など）で DAC のコードがランダムに切り替わるなら量子化誤差は一様分布で `k_act·LSB/√12 = 3.7 mrad rms`（500 kHz の更新で 250 kHz まで白色）。DSC1 の出力 LPF（100 kHz）と 914（28 kHz）で絞られ、`f_c` 以下はループで抑えられるので、真の位相への寄与は概算 1–2 mrad rms [AE]。
  外乱が小さいとコードが固定されたまま積分器が溜まって跳ぶ**リミットサイクル**（振幅 ±1/2 LSB = ±6.3 mrad）になりうる。どちらになるかは時間ループのシミュレーションで判定する（T16 と共有）。
  DSC1 と 914 の間に分圧器を入れると LSB とレンジが同じ比で縮む（トレードオフ、Phase 6 の探索対象）。

---

## 3. 仮定と適用範囲

- CW、単一周波数。カプラは損失なし 2×2（損失は `T_opt` と腕ごとの透過率にまとめる）。結合比は波長に依存しない（1548.8 nm 近傍）。
- PDB は線形（2 入力の差 < 4.5 mW、出力 < ±4.6 V）。DSC1 入力 ±4 V に収まる条件: `V_0(|ε| + 𝒱) + |V_off| < 4 V` → `P_det ≲ 4 mW`。3.33 mW は収まるが余裕は 20%。
- 線形化（式 5、6、10）は `|δφ_D| ≪ 1 rad` で有効。ロック外れ・引き込みは時間ループで非線形の式 (3) をそのまま使う。
- 914 の小信号応答は 1 次の極（28 kHz）と仮定（Fig. 1 は大信号）。915B の機械応答は不明（プレースホルダ）。
- `k_act` は静的に線形。ヒステリシス・クリープ（PZT）は T16 の範囲で、ここでは入れない。
- `μ` の式 (7) は、交差結合点が 1 か所（または腕ごとに 1 か所）で、2 軸の成分が互いに独立に干渉するとした近似。
- DSC1 の内部構造（フィルタの次数、遅延の配分、`h`、`K_i` の単位）は資料が不十分で、当てはめと解釈を含む（§2.7）。

## 4. 効く周波数帯

| 効果 | 帯域 | 備考 |
|---|---|---|
| 温度による ΔL 位相ドリフト、レンジ飽和 | < 1 Hz | T08。ロック時間を決める |
| PDB オフセットのドリフト、偏波 Δχ のドリフト、パワードリフト | < 1 Hz | ロック点の遅い移動（センサ側、補正不能） |
| 音響・振動（ファイバ、915B の巻き） | Hz – kHz | T08。ループで抑える |
| レーザー周波数雑音 × ΔL | DC – 1/τ（ΔL = 9 m で 23 MHz） | 低周波では白色（白色 FM のとき）。ループで抑えられるのは `f_c` 以下だけ |
| RIN（`(S − V_off)` 経由） | DC – 100 kHz（DSC1） | センサ雑音 |
| DSC1 入力雑音・ADC 量子化 | 100 Hz – 250 kHz（仕様の積分帯域） | センサ雑音 |
| DAC 量子化、914 の雑音 | DC – 100 kHz | アクチュエータ側の外乱 |
| ループの帯域 | `f_c` ≈ 0.7–10 kHz | 915B の共振しだい |

## 5. スケーリング

- **ΔL（ファイバ長）**: レーザー雑音 `S ∝ τ² ∝ ΔL²`（低周波）、全帯域 `σ² = 2πΔντ ∝ ΔL`。温度ドリフトの速さ ∝ ΔL、レンジが持つ温度幅 ∝ 1/ΔL。`Δχ` のドリフト ∝ ΔL（腕の長さの差）。
- **パワー `P_det`**: `K_d ∝ P`、`K_0 ∝ P` → 同じゲインで `f_c ∝ P`、`K_p,u ∝ 1/P`。センサ雑音（DSC1 入力、PDB）の位相換算 ∝ 1/P、ショット雑音 ∝ 1/√P。オフセット経由の RIN 結合 ∝ 1/P。`ε` によるロック点のずれは P によらない。上限は DSC1 入力 ±4 V と PDB 飽和 4.5 mW で、`P_det ≲ 4 mW`。
- **温度**: 環境の温度は ΔL と Δχ を通して効く。PDB のオフセットは 23 °C で零調整（PDB p.8）、温度係数は資料にない。
- **パルス幅**（将来系で同じ DC 干渉ロックを使う場合）: 2 つのパルスの重なりで可視度が決まる。ガウス型の場 `A ∝ exp(−t²/(2T_0²))` で遅延差 δτ のとき `|γ| = exp(−δτ²/(4T_0²))`。δτ には ΔL の群遅延の温度ドリフト（survey E26、T03・T08）が入る。分散の差でも重なりが下がる（T03・T06）。

## 6. 概算

### 6.1 Ref D（D02）。`T_opt = 1`、`𝒱 = 1`、`ℜ = 1.0 A/W`、8669-1 を仮定

| 量 | 値 | 根拠 |
|---|---|---|
| `V_0 = Gℜ P_det` | 3.33 V（`T_opt` を掛ける。0.8 なら 2.66 V） | D02、PDB p.19 |
| `|K_d|` | 3.33 V/rad | 式 (5) |
| `k_act` | 6.49 rad/V（8669-2: 6.90） | EO p.1 |
| `K_0` | 21.6 | 式 (10) |
| DSC1 の CV レンジ → 位相レンジ | 51.9 rad（±26 rad） | DSC1 p.2、EO p.1 |
| 1 K あたりの ΔL 位相（ΔL = 9 m） | 約 390 rad/K → レンジは 0.13 K 分 | [AE]、T08 |
| `ε = 0.05` のときのロック点のずれ | 50.0 mrad | 式 (4) |
| `V_off = 15 mV`（PDB の最大）のときのずれ | 4.5 mrad | 式 (4)、PDB p.19 |
| RIN の結合（`S − V_off = 15 mV`） | `δφ = 4.5×10⁻³ · δp`。RIN = −140 dB/Hz（プレースホルダ、T05）で 4.5×10⁻¹⁰ rad/√Hz → 無視できる | 式 (6) |
| パワーが 1% 変わったときのロック点（`V_off = 15 mV`） | 45 µrad | 式 (6) |
| 入力 PER 20 dB（x = 0.01）/ 30 dB（x = 0.001）でのロック点の揺れ | ±10.1 mrad / ±1.0 mrad（Δχ の 1 周期で） | 式 (7)、T09 |
| センサ雑音（v 換算） | DSC1 0.6–1 µV/√Hz（0.3 mV rms / 250 kHz、−120 dB V²/Hz）≫ PDB NEP 69 nV/√Hz > ショット 33 nV/√Hz | DSC1 p.2、PDB p.19 |
| センサ雑音（位相換算） | 3×10⁻⁷ rad/√Hz、帯域全体で 0.09 mrad rms | — |
| ADC / DAC の 1 LSB（位相） | 37 µrad / 12.7 mrad | §2.11 |
| DAC 量子化の寄与 | 1–2 mrad rms（ランダム）、リミットサイクルなら ±6 mrad | §2.11 [AE] |
| レーザー × ΔL（式 8）、低周波の ASD | Δν = 25 / 100 / 300 kHz、ΔL = 9 m（τ = 44 ns）で 2.5 / 4.9 / 8.6 ×10⁻⁵ rad/√Hz | T05 のプレースホルダ |
| 同、DC–100 kHz の rms | 7.8 / 15.6 / 27 mrad（ΔL = 5–15 m で ×0.56–1.67） | — |
| 同、全帯域 rms（`√(2πΔντ)`） | 0.083 / 0.17 / 0.29 rad | — |
| 帯域外の位相雑音による可視度の低下 | 約 1.4%（Δν = 100 kHz、ΔL = 9 m） | §2.6 |
| ループ: 積分器のみ、PM 60° | `K_i` = 1.75×10³ s⁻¹、`f_c` = 5.9 kHz（剛体） | §2.9 |
| 環境雑音 `S_0/f²` の残差（同ループ） | `σ = √(πS_0/(2f_c))`。`S_0 = 10⁻²` rad²·Hz で 1.6 mrad | §2.10、T08 |

Δν と RIN は Cobrite DX1 の資料が `refs/manuals/` にないため T05 待ち。参考までに、代理店に置かれたデータシートの写し（メーカーの原本ではない、UNVERIFIED）には線幅 typ. 25 kHz（狭線幅版）/ 80 kHz（標準版）、RIN 約 −140 dB/Hz（100 kHz–20 MHz）とあった。**100 kHz より下の RIN は見当たらなかった。** 表の 25 / 100 / 300 kHz は掃引の点で、DX1 の値ではない。

**Ref D の予測の主な結論（暫定）**: (i) ΔL が数 m 以上なら、DSC1 の帯域内の残差位相はレーザー周波数雑音 × ΔL で決まり、`f_c` より上の分は抑えられない（Δν = 25–300 kHz で 8–27 mrad rms、DC–100 kHz）。
(ii) ロック点の遅い移動は PDB のオフセット（≤ 4.5 mrad）と PER（入力 PER 20 dB で ±10 mrad）が支配的で、RIN・パワードリフトは無視できる。(iii) ループの帯域は 915B の共振で決まる。(iv) DAC の分解能（12.7 mrad/LSB）が mrad 級の残差の床を作りうる。(v) ロック時間は ΔL × 温度ドリフト速度で決まる（ΔL = 9 m、1 mK/s で 約 1 分）。

### 6.2 将来の 32 ポート構成（D07–D13、D30–D32）への含意

- Liu 型（D30）で probe と参照光の干渉を DC で読む場合、同じ式 (3)–(6) が使える。違いは `P_det`（probe は µW 級、D28 は Open）と遅延差 τ（ポート k の往復光路と参照光の差）。
  センサ雑音の位相換算は `v_n/(GℜP𝒱)` ∝ 1/P。例（プレースホルダ）: P = 10 µW、G = 10⁵ V/A（PDB450C の 4 MHz 設定）なら `V_0 = 1 V`、DSC1 の雑音で 1 µrad/√Hz 程度 [AE]。
- 式 (8) の τ を `k·T_rep` やポート間の遅延差に置き換えると、D31-3（ロックで補正できない部分、レーザー線幅への要求）の評価そのものになる（T05 の担当）。
- パルス光で DC 干渉を使う場合は §5 の重なり `|γ|` が `𝒱` に掛かる。

---

## 7. 検証テスト（実装が満たすべきもの）

既定のパラメータ: `P_in = 3.33 mW`、`T_opt = 1`、`ℜ = 1.0`、`G = 10³`、`k_act = 6.49 rad/V`、§2.7 の既定モデル。テストは決定的な部分は数値の許容差、確率的な部分は推定誤差で判定する。

1. **静的な縞（RAMP 相当）**: DSC1 の出力を ±4 V で掃引し、雑音を切ると `v(φ_D)` が式 (3) と相対誤差 10⁻¹² で一致。`ε = 0`、`𝒱 = 1` で振幅 3.33 V、フリンジ数 = 51.9/(2π) = 8.26。
2. **ロック点と傾き**: `ε = 0.05` で `φ_lock − π/2 = 50.02 mrad`、`V_off = 15 mV` で `4.505 mrad`、`K_d` が式 (5) と 10⁻⁶ で一致（数値微分）。
3. **RIN の結合の消失**: `S = V_off`（`ε` は任意）で、白色の RIN を入れても真の位相の PSD に RIN 由来の成分が出ない（`S − V_off = 15 mV` のときの値の 10⁻⁶ 未満）。`S − V_off = 15 mV` では `f ≪ f_c` で `S_φ = (15 mV/|K_d|)² S_RIN` に Welch 推定の誤差内で一致。
4. **閉ループ伝達関数（prompt §7-7）**: φ_D に正弦波の外乱を f = 10 Hz、100 Hz、1 kHz、`f_c/2`、`f_c`、`2f_c`、20 kHz で入れ、`δφ_D/φ_n` の振幅・位相が
   (a) 実装と同じ離散化の z 領域の解析式（式 9 の `C(z)`、2 点平均、ZOH、遅延サンプル数、`A_914`・`H_915B` を同じ方法で離散化）と相対 10⁻³・0.1° で一致、
   (b) 連続の `1/(1+L(i2πf))`（式 10）と `f ≤ f_c` で 5%・3° 以内で一致。
   P のみで `K_p` を上げたときに発振が始まる `K_p,u` と発振周波数が解析値（剛体で 0.093、40.9 kHz）と 2% 以内で一致（DSC1 p.19–20 の Ziegler–Nichols の手順の再現）。
5. **フィードバックなし → free-running（prompt §7-7）**: `K_p = K_i = K_d^{PID} = 0`（または STOP）で、真の φ_D の PSD が入力した外乱の和 `S_env + S_laser + S_act` と、50 区間以上の Welch 平均で各ビン ±3/√N_avg（相対）以内で一致。in-loop の `e_DSC1` は非線形の式 (3) の時系列と一致。
6. **レーザー × ΔL**: フィードバックなし、白色 FM（Δν = 100 kHz）、ΔL = 9 m で、φ_D の PSD が式 (8) に一致（零点 `f = m/τ` を含む、検出器帯域まで出力した場合）。全帯域の分散が `2πΔντ = 0.0277 rad²` に 5% 以内。
7. **積分器ループの残差の閉じた形**: `L = f_c/(if)`（理想積分器、遅延なし、テスト用の設定）と `S_env = S_0/f²` で `σ² = πS_0/(2f_c)` に 5% 以内（有限の f_min・f_max の補正は解析的に入れる）。
8. **レンジ飽和**: 一定の位相ランプ `r_φ` で、ロックが外れる時刻が `(26 rad − φ_0)/r_φ`（CV リミット ±4 V）と 1% 以内で一致し、外れた後は縞を滑る。
9. **PER**: `x = 0.01`、Δχ を 0–2π で掃引したとき、`|μ|` が [0.98, 1]、ロック点のずれの最大が `arcsin(x/(1−x)) = 10.10 mrad` に 10⁻⁶ で一致。
10. **DAC 量子化**: 外乱なし、ゆっくりした位相ランプで、`e_DSC1` の量子化による真の位相の揺れが ±k_act·LSB/2 = ±6.3 mrad の範囲に収まる（リミットサイクルの振幅の上限）。外乱を入れた場合の rms は §2.11 の概算と比べて記録する（合否なし、T16 と照合）。
11. **符号**: `loop_sign` を反転すると `φ_lock = +π/2` で正帰還になり発散し、`−π/2` に引き込まれる。

## 8. シミュレータが必要とするもの

### 8.1 ブロック（段 T）と入出力

| ブロック | 入力 | 出力 | パラメータ（既定） | 出典 |
|---|---|---|---|---|
| `cw_laser`（Ref D: DX1） | — | 出射時刻ごとの φ_L(t)、δp(t) | `S_ν(f)`、`S_RIN(f)`、パワードリフト | T05 |
| 光路の区間（腕 1: 915B ＋ピグテール、腕 2、共通部） | 環境雑音、温度、ストレッチャの伸び | ψ_seg,slow/fast(t)、遅延 | 長さ、n_g、B、雑音係数 | T08、T09、conventions §9 |
| `mzi_dc_fringe_detector`（PDB450C） | 2 腕の場（位相、パワー、偏光成分） | `v(t)` | `κ_1, κ_2, T_opt, ℜ, ε_R, G, V_off(t)`、NEP、ショット雑音、±4.6 V の飽和 | 式 (1)–(3)、PDB p.14, 19 |
| `adc`（DSC1 入力） | `v` | 500 kHz の PV | 120 kHz の入力帯域、1 MHz サンプル＋2 点平均、16 bit、±4 V、入力雑音 0.3 mV rms | DSC1 p.1–2 |
| `pid_dsc1` | PV | CV | 式 (9)、`S`、`O`、CV リミット、`h = 2 µs`、`ki_entry`、符号連動 | DSC1 p.5–6, 9, 17–20 |
| `dac`（DSC1 出力） | CV | `u(t)` | 12 bit、±4 V、100 kHz の LPF、220 Ω と負荷の分圧 | DSC1 p.2 |
| 遅延バッファ | — | — | `τ_DSC1` = 5.67 µs（既定）/ 8.33 µs（比較） | §2.7 |
| `hv_amp_914` | `u` | 圧電の駆動量（正規化） | 入力 ±5 V でフルスイング、小信号の極 28 kHz（プレースホルダ）、大信号の制限（Fig. 1 の表）、出力雑音（プレースホルダ） | EO p.1–2、T16 |
| `pzt_stretcher_915B` | 駆動量 | 腕 1 の ψ の変化 | `k_act` = 6.49 rad/V（8669-1）/ 6.90（8669-2）、`f_r`、`Q`（プレースホルダ） | EO p.1、T16 |
| 出力（段 M） | — | 真の φ_D(t) と PSD、`e_DSC1(t)` と PSD、CV、ロック状態 | — | prompt §6 |

### 8.2 オフラインの表・静的値（段 O）

- 914 の Fig. 1（大信号の振幅）の表（§2.7(c)）。
- `μ` の静的な値（入力 PER、腕ごとの PER）。

### 8.3 機器プリセット `presets/dsc1.yaml`（D05）

```yaml
# Thorlabs DSC1, values from User Guide Rev. B (2025-11-03), DOC-103674 (printed page in comments)
adc: {fs_hz: 1.0e6, avg_samples: 2, bits: 16, range_v: [-4, 4], input_bw_hz: 120e3,   # p.1, p.2
      input_noise_vrms: 0.3e-3, noise_band_hz: [100, 250e3], termination_ohm: 1.0e6}  # p.2
pid: {update_hz: 500e3, error: "S - v", integral: bumpless, gain_sign_locked: true,  # p.2, p.18, p.19
      gain_abs_range: [1.0e-5, 1.0e5], ki_entry: per_second}                         # p.6, p.12 (vs p.20: see note)
dac: {bits: 12, range_v: [-4, 4], output_bw_hz: 100e3, r_out_ohm: 220}               # p.2
loop: {phase_m180_hz_typ: 60e3, phase_m180_hz_min: 58e3,                             # p.2
       tau_s_fitted: 5.67e-6}                                                          # AE, this note §2.7
```

### 8.4 オン・オフのフラグ

`laser_phase_noise`、`rin`、`power_drift`、`offset_drift`、`optical_imbalance`、`per_crosscoupling`、`visibility_hf_averaging`、`adc_noise`、`adc_quantization`、`dac_quantization`、`amp914_large_signal`、`stretcher_resonance`、`env_noise`、`feedback`（§7-5 用）。

### 8.5 時間刻み

DSC1 の更新に合わせ `Δt = 2 µs`（またはその約数）。prompt §5.4(b) のチャンク処理なら、チャンク長 ≤ 実際のループ遅延（約 6 µs）なので 2 µs × 3 サンプル以下。レーザー雑音を DSC1 の帯域外まで出すなら、その成分は §2.6 の可視度の補正に回し、時間ループでは 500 kHz までに留めるのが現実的。

## 9. 他トピックへの依存とプレースホルダ

| 量 | 既定（プレースホルダ） | 決めるトピック |
|---|---|---|
| DX1 の `S_ν(f)`・Δν | 白色 FM、Δν を 25–300 kHz で掃引 | T05 |
| DX1 の RIN（特に 100 kHz 以下）、パワードリフト | −140 dB/Hz 白色、±0.01 dB/h（代理店の写し、UNVERIFIED） | T05 |
| 915B の共振 `f_r, Q`、容量 | 10–60 kHz、Q 5–30 を掃引 | T16 |
| 914 の小信号応答・雑音・入力インピーダンス | 1 次の極 28 kHz、雑音 0、高 Z | T16 |
| `ΔL`、各腕の長さ | 9 m（38 × 0.22 m = 8.4 m ＋ピグテール）、5–15 m を掃引 | D02（掃引）。巻き長は Photonics.com の抜粋（UNVERIFIED） |
| `S_env` | `S_0 f^{−α}`、`S_0 ∝ ΔL` | T08 |
| `dφ/dT`、`dB/dT`、`n_g` | 44 rad/(m·K)、2.4 rad/(m·K)、1.468 | T08、T09 |
| 入力 PER、腕の PER、PNH1550R5A2 の結合比の公差・過剰損失 | x = 0.01（20 dB）、κ = 0.5、`T_opt = 1` | T09（PNH1550R5A2 の資料なし） |

## 10. DECISIONS・conventions との食い違い、提案

1. **誤差の符号**: DSC1 は `e = S − v`（DSC1 p.18）、conventions §9 は `e = （線形結合）− 設定値`。物理量（位相）での誤差は conventions に従い、DSC1 ブロックの内部だけマニュアルの定義を使う、と conventions §9 に一文を足すことを提案。
2. **Ref D の腕と conventions §5 の対応**: 「腕 2（915B なし）↔ LO、腕 1（915B）↔ pr、`b_1` ↔ PDB INPUT+」を conventions §5 に追記することを提案（§2.2）。これで「設定値 +π/2」の意味が一意になる。ただし実機の極性は不明なので `loop_sign` をパラメータに残す。
3. **D02「915B の資料にある巻線ファイバ長からレンジを見積もる」**: 914/915B の資料（EO p.1）には巻線長がなく、代わりにフルレンジの光路長変化（16 / 17 µm）が直接書かれている。レンジはこれから求めた。ΔL の範囲の根拠には Photonics.com の掲載（UNVERIFIED）を使った。
4. **survey E82**: 「純遅延なら 約 8 µs」→ 入力・出力の帯域を極として入れると 5.7 µs（§2.7）。どちらも当てはめ。
5. **survey E75**: 「φ = π/2 では強度雑音は ε を通してだけ入る」→ 正確には `(S − V_off)/|K_d|` を通して入る（ε はロック点に吸収される、§2.4）。
6. **D11 の「PM の 2 軸の位相差 2.4 rad/(m·K)」**は、Ref D では PER と組んでロック点の揺れ（§2.5）として効く。T09 の確認待ち。
7. **survey §1 Ref D**: 「PDB450C 利得 10³ V/A（高インピーダンス負荷時）」は DSC1 の 1 MΩ 入力で満たされる（確認済み）。

## 11. ユーザの決定・情報が必要な点（推測しない）

1. **Ref D で使っている 915B はどちらか**（S/N 8669-1: 16 µm、8669-2: 17 µm）。`k_act` が 6% 変わる。
2. **DSC1 の設定値**（`K_p`、`K_i`、`K_d^{PID}`、`S`、`O`、CV リミット、どの項を有効にしているか）。これは測定結果ではなく構成の入力なので D04 と衝突しないと考えるが、ブラインド予測の入力にしてよいか判断してほしい。なければ §2.9 のゲインの掃引で予測を出す。
3. **DSC1 の `K_i` の単位**（1/s か、1 サンプルあたりか、§2.7(b)）。系を測るのではなく、Thorlabs の SDK 資料かサポートへの問い合わせで決めるのが D04 上安全。
4. **DSC1 出力と 914 入力の間**に分圧器・ケーブル以外のものがあるか。914 の入力インピーダンス。
5. **資料の追加**: PNH1550R5A2 の仕様書（結合比の公差、過剰損失、PER）、Opneti 偏光子の消光比、Cobrite DX1 のメーカーのデータシート、915B の機械共振と容量（Evanescent Optics / ハナムラオプティクスに問い合わせ可能か）。Thorlabs の製品ページは JavaScript で描画され、ここからは取得できなかった。
6. **DX1 と偏光子の間のファイバが SMF か PM か**（偏光子の手前の偏波ドリフトがパワードリフトになる量が変わる、§2.5）。

## 12. 文献

一次資料（このノートで読んだもの）:
- Thorlabs, "DSC1 Compact Digital Servo Controller User Guide," Rev. B, November 3, 2025, DOC-103674. `refs/manuals/thorlabs-dsc1.pdf`。
- Thorlabs, "PDB440x and PDB450x Balanced Amplified Photodetectors, Operation Manual," © 2019. `refs/manuals/thorlabs-pdb.pdf`。
- Evanescent Optics（販売: 株式会社ハナムラオプティクス）, "Model 914-2 Piezo Controller, Serial no. 1527" および "915B Piezo Assembly Serial 8669-1 / 8669-2"（3 頁、スキャン、日付なし）。`refs/manuals/914-915B-Specifications.pdf`。

書誌情報を doi.org で確認したもの（本文は未読。式 (8) と制御の標準的な扱いの出典として挙げる）:
- T. Okoshi, K. Kikuchi, A. Nakayama, "Novel method for high resolution measurement of laser output spectrum," Electron. Lett. 16(16), 630–631 (1980). doi:10.1049/el:19800437
- L. E. Richter, H. I. Mandelberg, M. S. Kruger, P. A. McGrath, "Linewidth determination from self-heterodyne measurements with subcoherence delay times," IEEE J. Quantum Electron. 22(11), 2070–2074 (1986). doi:10.1109/JQE.1986.1072909
- J. Bechhoefer, "Feedback for physicists: A tutorial essay on control," Rev. Mod. Phys. 77, 783–836 (2005). doi:10.1103/RevModPhys.77.783
- J. G. Ziegler, N. B. Nichols, "Optimum settings for automatic controllers," Trans. ASME 64, 759–765 (1942). doi:10.1115/1.4019264（DSC1 p.19–20 の調整則の元）

ウェブ（D24）:
- Photonics.com Buyers' Guide, "Model 915B Piezo-electric Fiber Stretcher"（Evanescent Optics Inc.）、https://www.photonics.com/Buyers-Guide/Products/Model-915B-Piezo-electric-Fiber-Stretcher/psp9424 、2026-10-28 に検索結果の抜粋のみ取得（ページ本体は HTTP 403）。**UNVERIFIED**（「9 m で 16 µm」「220 mm/巻」「DC–10 kHz、100 kHz で 0」）。メーカーの公式ページではない。
- Cobrite DX1 の代理店によるデータシートの写し（例: https://shop.laserdiodesource.com/ の CBDX1 データシート）、2026-10-28 に検索結果の抜粋のみ。**UNVERIFIED、メーカー原本ではないので D02 の値としては使わない**（T05 の参考）。

[AE] と書いた値はこの担当の概算で、一次資料による確認が必要。
