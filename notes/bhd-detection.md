# T13 `bhd-detection`: 平衡ホモダイン検出・平衡検出のモデル

- Phase: 1（deep-dive）／ 日付: 2026-10-28 ／ 担当効果: E57–E62（`notes/survey.md`）
- 規約: `notes/conventions.md`（§2 雑音レベル、§5 BHD と DC fringe、§7 片側 PSD）。DECISIONS D01–D33 を前提とする。
- 値の印: **[Ref A]** 等＝文献本文、**[manual]**＝`refs/manuals/thorlabs-pdb.pdf`（PDB440x/PDB450x Operation Manual, v2.4, 08-Apr-2019, M0009-510-1043）、
  **[manual-fig]**＝同マニュアルの図から目で読んだ値（±数 %、UNVERIFIED 扱い）、**[AE]**＝本ノート担当の概算・導出（未検証）、**[UNVERIFIED]**＝原典で確認できていない。

## 0. 記号（symbols, English）

- `ℜ`: photodiode responsivity [A/W]. `η_q = ℜ h ν / e`: photodiode quantum efficiency. `e`: elementary charge. `h ν`: photon energy (`ν = c/λ0`).
- `G`: transimpedance gain [V/A] (PDB450C: specified for high-impedance load; halved into 50 Ω [manual]).
- `P_1, P_2`: optical powers on the two photodiodes (PDB450C: `P_1` → INPUT+, `P_2` → INPUT−). `P_LO`: LO power before the BHD beam splitter.
- `I_tot = ℜ (P_1 + P_2)`: total DC photocurrent. `i_−(t)`: difference current. `V_RF(t)`: RF OUTPUT voltage.
- `S_SN = 2 e I_tot`: one-sided shot-noise PSD of the difference current [A²/Hz]. `S_el(f)`: one-sided input-referred electronic (circuit/dark) noise PSD [A²/Hz].
- `C(f) = |H(f)|² S_SN / S_el(f)`: shot-noise clearance (linear; `C_dB = 10 log10 C`).
- `S_App = 1 + C`: Appel et al.'s "signal-to-noise ratio" (vacuum noise incl. electronic noise over electronic noise).
- `η_e`: equivalent efficiency of electronic noise. `η_det = η_q · η_e · η_bal`: total detector efficiency (multiplies optical efficiency, conventions §2).
- `H(f)`: complex frequency response of the difference channel (normalized to 1 at low f). `H_1(f), H_2(f)`: responses of the two photodiode channels.
- `δ(f) = (H_1 − H_2)/(H_1 + H_2) + δ_BS`: complex imbalance (detector + beam-splitter split error `δ_BS = (T − R)/(T + R)`). `CMRR = |H_1| / |H_1 − H_2|` (Thorlabs definition, [manual] §3.6).
- `S_RIN(f)`: one-sided LO relative intensity noise [1/Hz]. `N_ex(f)`: LO amplitude excess noise in shot-noise units.
- `R = V/V_SN`: true normalized quadrature variance of the signal mode (conventions §2). `R_meas`: what the instrument reports normalized to its own measured shot noise.
- Ref D: `P_in` = power into the MZI (3.33 mW, D02), `η_opt` = MZI power transmission to the detector, `V` = fringe visibility, `ε` = DC imbalance (fraction of `η_opt P_in`), `Δφ = φ_LO − φ_pr` (conventions §5), `s = ±1` = sign set by which coupler output goes to INPUT+, `V_off` = electronic DC offset, `K_D` = error slope [V/rad].
- Pulsed: `T_rep = 1/f_rep`, `h(t)` = impulse response of the difference channel, `w(t)` = per-pulse integration window, `c_m` = inter-pulse crosstalk coefficient, `N_LO` = LO photons per pulse, `σ_el` = electronic noise per pulse [electrons rms].

## 1. 支配方程式

### 1.1 理想 BHD（規約 §5 の確認）

規約 §5 のとおり、強い LO で `n̂_1 − n̂_2 ≈ √2 |α_LO| x̂_θ`、`θ = φ_LO`。差電流の片側 PSD は
`S_i−(f) = S_SN · R(f)`、`S_SN = 2 e I_tot`（`I_tot ≈ ℜ P_LO`）。
`R(f)` は信号モードの直交位相の分散を `V_SN = 1/2` で割ったもの（サイドバンド周波数 `f_s = f` で評価）。
真空入力で `R = 1`（ショット雑音）。2 つの PD のショット雑音は独立に足されるので、平衡でも `S_SN = 2e(I_1 + I_2)`。

### 1.2 量子効率（E57）

PD の量子効率は純損失として扱う（規約 §2）: `V → η_q V + (1 − η_q) V_vac`。
`η_q = ℜ h ν / e`。1548.8 nm で `e/(hν) = 1.249 A/W` [AE] なので、
PDB450C の `ℜ = 1.0 A/W`（表の typ. 値 [manual]）→ `η_q ≈ 0.80`、図の読み `ℜ ≈ 1.03 A/W @ 1549 nm` [manual-fig] → `η_q ≈ 0.82`。
**PDB450C は squeezing の測定には使えない（損失 約 20%）**。Ref D（古典の DC 干渉）では効率は信号の大きさにしか効かない。
窓の反射、受光面からのはみ出し（PDB の受光径 0.3 mm）も「光電変換損失」としてここに含める（Ref A は 約 1% [Ref A §4.2]）。

### 1.3 回路雑音（暗雑音）と等価損失（E59）

電子雑音は差電流に足される: `S_meas(f) = G²|H(f)|² S_SN R(f) + G² S_el(f)`（`S_el` は入力換算、`H` の後で加わると仮定）。
装置は「LO だけ（信号ポートは真空）」で測った値をショット雑音の基準とするので、
`R_meas = (C R + 1)/(C + 1) = η_e R + (1 − η_e)`、`η_e = C/(1 + C)`。
すなわち**電子雑音は効率 η_e の純損失と厳密に等価**（基準化の仕方による）。Appel et al. (2007) Eq. (14) は
`η_eq = (S_App − 1)/S_App`、`S_App = (α² + T)/T`（LO だけの雑音（電子雑音込み）/ 電子雑音）で、`S_App = 1 + C` とおけば同じ式になる。

- 「クリアランス X dB」の定義が文献で `C` か `S_App` かは書かれていないことが多い。X ≥ 20 dB では差は無視できるが、Ref C の 5 dB では等価損失が 24%（`C` 定義）と 32%（`S_App` 定義）で大きく違う [AE]。
- `C ∝ P_LO`（CW）、パルスでは `C ∝ N_LO`（1 パルスあたりの LO 光子数）。等価損失 `1 − η_e ≈ 1/C ∝ 1/P_LO`。
- 電子雑音を差し引く（"dark noise corrected"、Ref A Fig. 4(b)、Ref C の全値）と、期待値では `R_corr = (S_meas − G²S_el)/(G²|H|²S_SN) = R`。残るのは統計誤差だけ（T04 の ESA モデルで扱う）。
  Ref A の 10.1 dB は**補正なし**、Ref C の 4–4.7 dB は**補正あり** [Ref A §3.2、Ref C Fig. 4]。

### 1.4 周波数応答と干渉ディップ（E58）

`H(f)` は差のチャネル全体の応答。モデル（パラメトリック、[AE]）:
`H(f) = H_TIA(f) · H_dip(f)`、
- `H_TIA(f)`: 1 次または 2 次の低域通過（`f_3dB`、Q）。PDB450C の G = 10³ では DC–150 MHz（−3 dB）[manual]、図では 100 MHz まで平坦、130 MHz 付近に約 +1 dB のピーク、200 MHz で約 −15 dB [manual-fig]。
- `H_dip(f) = (1 + r e^{−i 2π f τ_d})/(1 + r)`: 遅延 `τ_d` の 2 経路の干渉による凹み。最小は `f_dip = 1/(2τ_d)`、深さ `(1 − r)/(1 + r)`。
  Ref C は自作 BHD で「信号出力の destructive interference」を 0.96 GHz（BHD_A）、0.83 GHz（BHD_B）に報告 [Ref C Fig. 5 caption] → `τ_d = 0.52 ns`、`0.60 ns` [AE]。
  機構（RF の反射、PD 間の配線遅延など）と深さ `r` は Ref C に書かれていないので、`r` は掃引するパラメータ（UNVERIFIED）。

`H` は信号とショット雑音に同じようにかかるので、**正規化した squeezing には `C(f)` を通してだけ効く**（`H` が小さい周波数ではクリアランスが下がり、等価損失が増える）。
ただし `H_1 ≠ H_2`（2 つの PD チャネルの応答差）は 1.5 の不平衡として別に効く。
Ref A は 100 MHz での squeezing の低下を「帯域 35 MHz の PD の検出効率の劣化」によるとしている [Ref A §3.2]。線形フィルタだけなら正規化で打ち消すはずなので、
これは周波数に依存する効率（キャリア走行など）か不平衡を意味する。本ノートでは `η_q(f)` を任意のパラメータとして持たせ、既定は一定とする [AE]。

### 1.5 不平衡・CMRR と LO の古典雑音の漏れ（E60）

2 つの PD チャネルの重み `H_1, H_2` とビームスプリッタの分岐誤差を入れると [AE]
`i_− ∝ (H_1 + H_2)/2 · [ x̂_θ(signal) + δ(f) · X̂_LO,amp ]`、
`X̂_LO,amp` は LO の振幅直交位相の揺らぎ（コヒーレントなら真空、分散 1/2）。古典雑音があると
`Var(X̂_LO,amp)/V_vac = 1 + N_ex(f)`、`N_ex(f) = η_q P_LO S_RIN(f)/(2 h ν)`（片側 RIN、規約 §7）。
LO だけの基準にも同じ漏れが入るので
`R_meas = (R + x)/(1 + x)`、`x(f) = |δ(f)|² (1 + N_ex(f))`、
すなわち**等価効率 `η_bal = 1/(1 + x)` の純損失**。`x` の真空部分 `|δ|²` だけでも損失になる。

- CMRR との関係: `|H_1 − H_2|/|H_1| = 10^{−CMRR/20}` → `|δ_det| ≈ 10^{−CMRR/20}/2`、`|δ_det|² ≈ 10^{−CMRR/10}/4`。
  PDB450C の CMRR >25 dB（typ. >30 dB）[manual] → `|δ|² ≤ 7.9×10⁻⁴`（typ. 2.5×10⁻⁴）。
- 光路長差による CMRR の上限 [AE]: 2 入力の遅延差 `Δt = n_g Δl / c` で `|H_1 − H_2|/|H| = 2|sin(π f Δt)|`。
  ファイバで `Δl = 1 cm` なら 100 MHz で 30 dB、1 GHz で 10 dB。マニュアルは「光路長差の許容値の図」を参照しているが、その図は PDF に**入っていない**（§3.7 が参照する "next page" は保守の章）。
- EDFA の ASE（T05/E17）: ASE と LO の搬送波のビートは同相成分で、上の `N_ex` に入る（CMRR で抑えられる）。
  一方、LO に乗った広帯域 ASE が広帯域の squeezed 光（遠く離調した反 squeezing 成分）とビートする項は差の信号そのものなので、**CMRR では抑えられない**。Ref A がフィルタキャビティを入れた理由 [Ref A §3.1]。この項の大きさは T05 が出す [AE]。

### 1.6 飽和・線形性（E61）

- PDB450C の RF 出力の最大振幅 ±4.6 V（高インピーダンス負荷、G = 10³–10⁶）、50 Ω 負荷では半分 [manual]。CW 飽和パワー 4.5 mW @ 1550 nm（**2 入力のパワーの差**に対する値）[manual]。
  `4.6 V/(ℜ G) = 4.47 mW`（ℜ = 1.03）で、飽和は TIA の出力振幅で決まると読める [AE]。
- MONITOR 出力: 10 V/mW、最大 10 V（高インピーダンス）→ **1 PD あたり 約 1 mW で飽和** [manual §3.3]。
- 損傷しきい値 20 mW（1 PD）[manual]。線形性のために 4 W/cm² 以下を推奨（受光面全体を照らすこと）[manual §3.2]。
  0.3 mm 径の面全体で 4 W/cm² は 2.8 mW。ファイバ出射の実際のスポット径は書かれていないので、FC 入力で線形性の条件を満たしているかは不明（UNVERIFIED）。
- モデル: 出力のクリップ `V_RF → clip(V_RF, ±V_sat)` と、任意の弱い圧縮 `ℜ_k(P) = ℜ_k (1 − P/P_nl,k)`（`P_nl,k` は未知、既定 ∞）。
  2 つの PD で `P_nl` が違うと、強い LO で DC の差と CMRR が変わる。

### 1.7 Ref D: 平衡 PD による DC 干渉信号（T15 への入力）

2 つの 2×2 カプラの間の MZI。規約 §5 の DC 干渉の式を、腕 A（915B のある腕）を `pr`、腕 B を `LO` の役として書く:
`P_1 − P_2 = s · η_opt P_in [ε + V cos Δφ]`、`P_1 + P_2 = η_opt P_in`、
`V_RF(t) = ℜ G (P_1 − P_2)(1 + δ_P(t)) + V_off(t) + n_V(t)`、クリップ ±4.6 V。

- `s = ±1`: 実際の対称カプラ（交差で i）を 2 回通ると、規約 §5 の実数のビームスプリッタとの差は符号だけになる（`P_out2 − P_out1 = +η P cos(ψ_A − ψ_B)` [AE]）。どちらの出力を INPUT+ に入れたかで `s` が決まる。
  `s` はロックが安定になる零点（+π/2 か −π/2）を選ぶので、DSC1 の符号と組で T15 が使う。
- `ε`: 2 段目のカプラの分岐比の誤差、出力側の損失差、2 PD の応答の差 `(ℜ_1 − ℜ_2)/(ℜ_1 + ℜ_2)` の和。差信号の零点は `cos Δφ = −ε/V`、すなわち `Δφ = π/2 + ε/V`（小さい ε）にずれる。
- **RIN の結合** [AE]: `ε` が光学的な静的不平衡なら、差信号全体が `P_in` に比例するので、零点でのロックは 1 次で RIN に感じない。
  RIN が効くのは (i) `V_off` のような `P_in` に比例しない項があるとき（ロック点 `V_RF = 0` が `Δφ` で `V_off/K_D` ずれ、その量が `δP/P` で揺らぐ）、
  (ii) 設定値が零点以外のとき、(iii) 周波数に依存する不平衡 `δ(f)`（kHz 帯では小さい）。詳細は T15。
- `K_D = ℜ G η_opt P_in V`（`Δφ = π/2`）[V/rad]。

### 1.8 パルス BHD（E62、将来の構成 F）

1 パルスの直交位相の推定値 `X_j = ∫ w(t − j T_rep) (h ∗ i_−)(t) dt`。パルス m の電荷 `q_m` が推定値 j に入る割合
`c_{j−m} = ∫ w(t) h(t + (j − m) T_rep) dt / ∫ w(t) h(t) dt`（`c_0 = 1`）。
`X_j^meas = Σ_m c_m q_{j−m} + n_el,j`、`Var(n_el,j) = σ_el² e²`。

- ショット雑音（1 パルス）: `Var(q_j)|_SN = e² η_q N_LO`。クリアランス `C_pulse = η_q N_LO / σ_el²`（`σ_el` は電子数 rms）。等価損失は §1.3 と同じ形。
- 1 次の低域通過（時定数 `τ = 1/(2π B)`）で `w` がパルスより十分短いとき `c_1 = exp(−2π B / f_rep)` [AE]。`c_1 < 10⁻²` には `B/f_rep > 0.73`。
  リンギング（2 次で Q > 0.5）があると `c_m` は振動し、符号が変わる。
- 隣のスロットが全部同じ状態で独立なら、`c_m` は正規化した分散を変えない（ショット雑音の基準も同じように混ざる）。効くのは
  (a) 隣のスロットの中身が違うとき（真空のスロット、明るい probe パルス、振り分け素子の有限の消光比で漏れたパルス）: 平均値の漏れ `c_m ⟨q_{j−m}⟩` と雑音の混入、
  (b) パルス間の相関を使う測定（時間多重のモード間の相関、GBS の入力の評価）: 見かけの相関 `ρ_adj = Σ_m c_m c_{m+1} / Σ_m c_m²`。
- `B ≪ f_rep` のとき BHD はパルス列を平均し、`f < f_rep/2` のサイドバンドで CW と同じように squeezing スペクトルを測る（測るモードは LO のパルス形で決まる、T03）。
  将来の BHD がパルス分解型か周波数領域型かは未定（D15、§8 の問い）。

## 2. 仮定と適用範囲

- 強い LO（`|α_LO|² ≫` 信号の光子数）、線形化した BHD。LO と信号の空間・時間モードの重なりは T12・T03 が扱い、ここでは効率として受け取る。
- 電子雑音は信号と独立なガウス雑音で、`H` の後で加わる（TIA の入力換算雑音 `S_el` を `|H|²` で割ってもよい。どちらかを設定で選ぶ）。
- 不平衡は小さい（`|δ| ≪ 1`）。`δ` の 2 次まで。
- PDB450C の値はマニュアルの typ./min. 値で、個体差がある（ℜ 曲線は「典型値、個体でずれる」と明記 [manual §3.3]）。D04 により実測で合わせ込まない。
- PDB450C のスペクトル雑音の図の縦軸は「Amplitude [dB]」とだけ書かれている。ESA（RBW 100 kHz、VBW 10 kHz [manual §5]）の 50 Ω 入力での dBm と仮定して換算した（UNVERIFIED）。
  この仮定で得た平均（§4.1）は NEP の値と矛盾しない。
- パルス BHD の式（§1.8、§4.3）は、検出器が各パルスの電荷を線形に積分し、応答 `h(t)` が時不変であると仮定する。

## 3. 効果が効く周波数帯

| ID | 効果 | 周波数帯 |
|---|---|---|
| E57 | 量子効率 | 全帯域（静的。Ref A の解釈では PD の帯域付近で低下） |
| E58 | 周波数応答・ディップ | `f_3dB` 付近（Ref A: 35 MHz、PDB450C: 100–200 MHz、Ref C: 0.83・0.96 GHz） |
| E59 | 回路雑音 | 測定周波数（3 MHz、30 MHz–1.2 GHz）。Ref D の制御では DC–120 kHz（DSC1 の入力帯域） |
| E60 | 不平衡・CMRR | LO の RIN・ASE があるところ（kHz–MHz の緩和振動など、T05）。CMRR は光路長差で高周波ほど悪化 |
| E61 | 飽和・非線形 | DC と出力振幅の最大値（Ref D の縞の山） |
| E62 | パルス BHD | `f_rep` とその高調波。`B` と `f_rep` の比 |

## 4. 概算

### 4.1 Ref D: PDB450C、G = 10³ V/A、P_in = 3.33 mW（D02）

前提: `ℜ = 1.0 A/W`（表）／ 1.03 A/W（図の読み）、`η_opt = V = 1`、`ε = 0`（**上限**。カプラ PNH1550R5A2 の分岐比・過剰損失は資料がなく未確認 → T12/T15 の値に置き換える）。
DSC1 の入力終端は 1 MΩ [`thorlabs-dsc1.pdf` 仕様表] なので、高インピーダンス負荷の `G = 10³ V/A` を使う。**途中に 50 Ω 終端があれば G = 500 V/A で、以下の電圧は半分**。

| 量 | 値 | 出典 |
|---|---|---|
| 帯域（−3 dB） | DC–150 MHz（100 MHz まで平坦、130 MHz に +1 dB のピーク） | [manual] / [manual-fig] |
| 変換利得 | 1.0×10³ V/W（ピーク応答度） | [manual] |
| NEP（DC–150 MHz） | 68.6 pW/√Hz → 入力換算 `i_n ≈ 69 pA/√Hz`（帯域平均） | [manual] |
| 雑音スペクトル（暗、G = 10³） | 低周波 約 −88.5 dB、130 MHz で 約 −81 dB（基準の床 −90 dB）→ 入力換算 約 29 pA/√Hz（低周波）、約 118 pA/√Hz（130 MHz） | [manual-fig]、dBm/100 kHz・50 Ω と仮定 [AE] |
| CMRR | >25 dB（typ. >30 dB） | [manual] |
| DC オフセット | < ±15 mV（23 °C で 0 に調整、5 分程度の暖機ドリフト） | [manual] |
| 飽和 | 出力 ±4.6 V（高 Z）、差のパワー 4.5 mW @ 1550 nm | [manual] |
| DC fringe の振幅 | `ℜ G η_opt P_in V` = **±3.33 V**（ℜ = 1.0）／ ±3.43 V（ℜ = 1.03） | [AE] |
| 傾き `K_D`（Δφ = π/2） | **3.3–3.4 V/rad** | [AE] |
| 余裕 | PDB の ±4.6 V まで 26–28%、**DSC1 の入力 ±4 V まで 14–17%**（3.9 mW の差で DSC1 が先にクリップ） | [AE] |
| 1 PD あたりのパワー | 縞の山で最大 3.33 mW、ロック点（π/2）で 1.67 mW。損傷 20 mW、飽和 4.5 mW より下 | [AE] |
| MONITOR 出力 | ロック点で 1.67 mW → **飽和**（1 mW で 10 V）。規格化・強度モニタには使えない | [manual] / [AE] |
| ショット雑音 | `I_tot = 3.33–3.43 mA` → `√(2eI) = 33 pA/√Hz` | [AE] |
| ショット雑音 / 電子雑音 | 低周波（DSC1 の帯域）**+1.2 dB**、NEP の帯域平均 −6.3 dB、130 MHz −11 dB | [AE] |
| 位相換算の雑音（PDB のみ） | `√(S_SN + S_el)·G / K_D ≈ 1.3×10⁻⁸ rad/√Hz` → 120 kHz で 4.4 µrad rms | [AE] |
| 比較: DSC1 の入力雑音 | 0.3 mV rms（100 Hz–250 kHz）→ 87 µrad。DSC1 の雑音床 −120 dB V²/Hz → 0.29 µrad/√Hz | [`thorlabs-dsc1.pdf`] / [AE] |
| 比較: PDB の DC オフセット | 15 mV → ロック点のずれ 最大 4.4 mrad（静的。ドリフトの温度係数は資料にない） | [AE] |

結論 [AE]: Ref D では、検出器の量子雑音・電子雑音はロックのセンサ雑音として DSC1 の入力雑音より約 20 倍（振幅）小さく、無視できる。
効くのは (i) 縞の振幅＝傾き（`η_opt`、`V`、`ℜ`、負荷インピーダンス）、(ii) DC オフセットとそのドリフト、(iii) クリップ（DSC1 の ±4 V が先）。
ただし PDB の出力の広帯域雑音は 150 MHz 全体で `69 nV/√Hz × √(150 MHz) ≈ 0.84 mV rms` あり、DSC1 の 0.3 mV より大きい。
**DSC1 の入力に 120 kHz のアナログのアンチエイリアスフィルタがあるかは DSC1 の資料に書かれていない**（「入力の典型帯域 120 kHz」とだけある）。なければ折り返しで雑音が 1 桁近く増えうる → T16 で確認する点。

### 4.2 CW の高 QE BHD（Ref A/B/C と、将来の構成 F の候補、D15 で Open）

| 値 | Ref A | Ref B | Ref C |
|---|---|---|---|
| PD の QE | 約 99%（Laser Components IGHQEX0500-1550、2 個）| 0.998（Si、Hamamatsu S3590-06、860 nm 付近） | 97% |
| 帯域 | PD の帯域 35 MHz | — | GHz（自作、Tohermes et al. 2024 の型） |
| クリアランス | 約 28 dB @ 3 MHz、LO 15.4 mW | 回路雑音 −21.7 dB @ LO 5.5 mW | 11 / 7.4 / 5 / 5.7 dB @ 30 MHz / 0.25 / 0.69 / 1.21 GHz |
| 等価損失 `1/(1+C)` | **0.16%**（本文「<0.2%」と一致） | 0.67% | 7.4% / 15% / 24% / 21%（`S_App` 定義なら 7.9 / 18 / 32 / 27%） |
| 特徴 | 電子雑音の補正なしで 10.1 dB | — | 0.96 GHz（BHD_A）、0.83 GHz（BHD_B）の干渉ディップ。値は暗雑音補正あり |

出典: [Ref A §2–4.2]、[Ref B §4]、[Ref C §2–3、Fig. 4–5]。

- 将来の構成 F の既定値（Open、D15）: QE 0.99、クリアランス 28 dB @ 15.4 mW を `C ∝ P_LO` で換算（Ref A の値の代用）。GHz 帯の評価には Ref C の値を使う。
  10 dB 目標で全体の損失収支が 8% 程度のとき、検出器（QE＋電子雑音＋不平衡）に許せるのは 1–2% のオーダー [AE]（T04/T12 の損失収支で確定する）。
- LO の RIN の漏れ（§1.5）[AE]: LO 15.4 mW、`η_q = 1` で `N_ex = P_LO S_RIN/(2hν) = 60`（RIN −150 dB/Hz）、6（−160）、0.6（−170）。
  CMRR 30 dB（`|δ|² = 2.5×10⁻⁴`）なら `x = 0.015 / 0.0015 / 0.0004`。RIN −150 dB/Hz では **x = 1.5% の等価損失**で、10 dB の squeezing が 9.45 dB に下がる。
  CMRR 25 dB なら 4.7%。必要な CMRR は測定周波数での LO の RIN（T05）で決まる。
- PDB450C を量子測定に流用した場合（参考）[AE]: QE 約 0.8、G = 10³ の NEP 69 pA/√Hz に対して LO 15 mW のショット雑音は 約 70 pA/√Hz で、クリアランス 約 0 dB。使えない。

### 4.3 パルス BHD（F、パラメータ）

- 繰り返し `f_rep` は Open（D07）。パルス分解には `B ≳ f_rep` が要る（1 次の応答で `c_1 = exp(−2πB/f_rep)`、`B = f_rep` で `c_1 = 1.9×10⁻³` [AE]）。
- 文献の例: Hansen et al. 2001 は繰り返し 最大 1 MHz、CMRR >85 dB、電子雑音 730 電子/パルス、個々のパルスの量子雑音を 14 dB の SNR で測定（abstract から）。
  この `σ_el = 730` で 14 dB のクリアランスには `η_q N_LO ≈ 1.3×10⁷` 光子/パルス [AE]。
  Kouadou et al. 2025（arXiv、査読前）は 150 MHz までのパルス分解 BHD（1550 nm 帯）を報告（abstract のみ確認、性能の数値は未読）。
  76 MHz のパルス分解測定（Opt. Lett. 33, 1458 (2008)、隣のパルスとの相関でパルスの独立性を確認）は書誌の一部だけ確認（著者・題名 UNVERIFIED）。
- ファイバ長・パルス幅への依存: 検出器は電荷を積分するので、`τ_p ≪ 1/B` ならパルス幅に依存しない。パルス幅が効くのは LO とのモード整合（T03）と、短パルスの尖頭電流による PD の飽和（定量値なし、UNVERIFIED）。

## 5. スケーリング

| 変数 | 依存 |
|---|---|
| ファイバ長 | 直接はない。LO の配送ファイバ（D09）の損失で `P_LO` が下がると `C ∝ P_LO` が下がる（0.2 dB/km 程度なら 100 m 以下で無視できる）。光路長差は CMRR を下げる（§1.5） |
| パワー | `C ∝ P_LO`（CW）、`∝ N_LO`（パルス）。LO の RIN の漏れ `x ∝ |δ|² P_LO S_RIN`（**P_LO を上げると増える**）。飽和は差のパワーと出力振幅で決まる。Ref D の `K_D ∝ P_in` |
| 温度 | PDB の DC オフセット（23 °C で調整、暖機ドリフト。係数は資料なし）。応答度の温度依存は資料なし |
| パルス幅 | §4.3 |
| 繰り返し | パルス間クロストーク `c_1 = exp(−2πB/f_rep)`（1 次の場合） |

## 6. 検証テスト（実装が通すべきもの）

1. **電子雑音の等価損失**: クリアランス C で `R_meas = (C R + 1)/(C + 1)` が厳密に成り立つ。C = 28 dB で等価損失 0.158%（Ref A の「<0.2%」）、21.7 dB で 0.67%。
   `subtract_dark = on` で `R_corr = R`（期待値）、モンテカルロでは相対誤差 `√(2/N)` 以内。
2. **純損失としての QE**: `η_q = 0.99` で `V_out = η_q V_in + (1 − η_q)/2`（prompt §7-3 と同じ実装を通ること）。
3. **不平衡の漏れ**: `δ = 0` なら LO の RIN をどれだけ大きくしても `R_meas = R`。`|δ|² = 2.5×10⁻⁴`、`N_ex = 60` で R = 0.1（−10 dB）が `R_meas = 0.1135`（**−9.45 dB**）。真空入力では常に `R_meas = 1`。
4. **周波数応答**: `H_dip` を `f_dip = 0.96 GHz`、`r = 0.9` で入れると、|H| の最小が 0.96 GHz で深さ `(1−r)/(1+r)`。電子雑音なしでは正規化した R(f) は変わらず、電子雑音ありではディップの位置でだけ等価損失が増える。
5. **Ref D の DC fringe**: `P_in = 3.33 mW`、`ℜ = 1.0 A/W`、`G = 10³`、`η_opt = V = 1`、`ε = 0`、`V_off = 0` で、Δφ を 0 → 2π と回すと `V_RF = ±3.33 V`、Δφ = π/2 での傾き −3.33 V/rad（s = +1）。
   `P_in = 5 mW` にすると縞の山が PDB で ±4.6 V、DSC1 の入力で ±4.0 V に平らに切られる。50 Ω 負荷の設定で振幅と傾きがちょうど半分。
6. **零点のずれと RIN**: `ε = 0.03`、`V = 1` で零点が `Δφ = π/2 + 0.030 rad`（arcsin の厳密値と 10⁻⁴ rad 以内）。`V_off = 0` なら `P_in` に 1% のステップを入れても零点は動かない。`V_off = 15 mV` なら零点は 4.5 mrad ずれ、1% のステップで 45 µrad 動く。
7. **パルス間クロストーク**: 1 次の `h(t)` と短い `w(t)` で、独立な真空パルス列の隣接相関が解析値 `Σ c_m c_{m+1} / Σ c_m²`（1 次では `c_1`）と一致し、`B → ∞` で 0。各パルスの正規化分散は `c_m` によらず R（全スロットが同じ状態のとき）。
8. **ショット雑音のスケーリング**: 真空入力で、出力の PSD が `P_LO` に比例（電子雑音を引いた後）、電子雑音を足すと C が `P_LO` に比例。

## 7. シミュレータがこのトピックから必要とするもの

| ブロック / 表 | 段 | 入力 | 出力 | トグル |
|---|---|---|---|---|
| `detector_efficiency`（静的） | O | `η_q`（または ℜ, λ0）、`C(f)` または `S_el(f)`, `P_LO`、`δ(f)`、`S_RIN(f)`（T05） | `η_det(f) = η_q η_e η_bal`、等価損失の内訳 | `det.qe`, `det.circuit_noise`, `det.imbalance_leak` |
| `bhd_frontend` | M | 信号の R(f) または直交位相の時系列（T04/T01–T03）、`H(f)`（TIA＋dip）、`S_el(f)`、`G` | 差電流・電圧の時系列 / PSD（ESA 模擬 T04 へ）、ショット雑音の基準 | `det.freq_response`, `det.dip`, `det.circuit_noise`, `det.subtract_dark`, `det.saturation` |
| `balanced_pd_dc`（Ref D） | T | `P_in(t)`（RIN を含む、T05）、`Δφ(t)`（区間の位相、規約 §9）、`η_opt`, `V`, `ε`, `s`, `ℜ`, `G`, 負荷（high-Z/50 Ω）、`V_off(t)`、`S_el` | `V_RF(t)` → 誤差信号（T15）、クリップのフラグ | `det.offset_drift`, `det.saturation`, `det.circuit_noise`（既定 on、影響は小） |
| `pulsed_bhd_kernel` | O | `h(t)` または `B, Q`、`w(t)`、`f_rep`、`σ_el`、`N_LO` | `c_m`、`C_pulse`、隣接相関 | `det.pulsed_crosstalk` |
| 機器プリセット | 設定 | — | `pdb450c_g1e3`（表 4.1）、`bhd_refA`、`bhd_refB`、`bhd_refC_A`（dip 0.96 GHz）、`bhd_refC_B`（0.83 GHz）、`bhd_generic_highqe`（Open） | — |

PDB450C の 150 MHz の応答は、DSC1 の 500 kHz の更新では「定数＋数 ns の遅延」とみなせる。時間ループ段では既定で `H = 1`、必要なら 2 次の低域通過をオプションで持つ。

## 8. ほかのトピックとの関係

- T04: `R_meas` の ESA 模擬（RBW/VBW、ログ平均）と、電子雑音の補正の有無。T13 は `H(f)`, `S_el(f)`, `η_det(f)` を渡す。
- T05: LO の RIN `S_RIN(f)`・ASE（`N_ex` と、CMRR で抑えられない ASE×squeezed のビート）。Ref D の `δ_P(t)`。
- T12: QE 以外の損失収支、カプラ（PNH1550R5A2）の分岐比と損失 → `η_opt`, `ε`。T09: 偏光のずれ → `V`。
- T14: Ref A/B/C の誤差信号は BHD の出力から作る（帯域、Ref C の 8–300 MHz の包絡線検出の帯域で `H(f)` と `C(f)` を使う）。
- T15: §1.7 と表 4.1 をそのまま使う。T16: DSC1 の入力帯域・アンチエイリアス。
- T03/T11: パルスの時間モード、振り分け素子の消光比で隣のスロットに漏れるパルス（§1.8 の (a)）。

## 9. 調整役への提案と未決の点

提案（調整役のファイルは編集していない）:
1. `notes/survey.md` E61 の「PD 1 個あたり 約 1.7–3.3 mW、飽和 4.5 mW に近い」は、**飽和の規定が 2 入力の差に対するもの**であることを補う。Ref D で先にクリップするのは DSC1 の ±4 V（差 3.9 mW）で、余裕は 14–17%。MONITOR 出力は 1 mW で飽和するので Ref D では使えない。
2. `notes/conventions.md` §5 に「物理的な 2×2 カプラは対称（交差で i）。MZI では規約の実数の BS との差は符号 `s = ±1` だけで、`s` はどちらの出力を正の入力に入れたかで決まる」を加えるとよい。
3. クリアランスの dB 値は `C = S_SN/S_el` として扱う（Ref A の値と本文の <0.2% が一致）。Ref C の値は暗雑音補正済みなので再現には影響しないが、補正なしの再現をする場合は定義を明記する。

ユーザの決定が必要な点（推測しない）:
1. **PDB450C の RF 出力から DSC1 までの終端**: DSC1 に直結（1 MΩ、G = 10³ V/A）か、途中に 50 Ω 終端・T 分岐・オシロスコープの 50 Ω 入力があるか（G = 500 V/A、縞の振幅と傾きが半分）。
2. 使っている PDB450C が標準の DC 結合版か（D02 の DC lock からは DC 結合と推定）。
3. 応答度の既定値を表の 1.0 A/W とするか、曲線の読み 約 1.03 A/W とするか（D02 は「メーカー資料の値」。本ノートは 1.0 を既定、1.03 を感度として扱った）。
4. 将来の BHD をパルス分解型（`B ≳ f_rep`）にするか、周波数領域型（`B ≪ f_rep`）にするか（D15 の Open。`f_rep` も Open）。

## 10. 文献

一次資料（本文・該当部分を読んだもの）:
- [Ref A] K. Hirota, T. Kashiwazaki, et al., "Generation of 10-dB squeezed light from a broadband waveguide optical parametric amplifier with improved phase locking method." §2–4.2。**掲載誌・年・DOI は PDF から読めない（UNVERIFIED、survey と同じ）。**
- [Ref B] Y. Takeno, M. Yukawa, H. Yonezawa, A. Furusawa, "Observation of −9 dB quadrature squeezing with improvement of phase stability in homodyne measurement," Opt. Express 15, 4321–4327 (2007). DOI は PDF に記載なし（UNVERIFIED）。§2、§4。
- [Ref C] S. Verclas, C. Darsow-Fromm, B. Tohermes, R. Schnabel, "Fiber distribution of phase-stabilized GHz-bandwidth squeezed vacuum states of light between two buildings," Optica Quantum 4, 1 (2026). DOI 10.1364/OPTICAQ.567418. §2–3、Fig. 4–5。
- J. Appel, D. Hoffman, E. Figueroa, A. I. Lvovsky, "Electronic noise in optical homodyne tomography," Phys. Rev. A 75, 035802 (2007). DOI 10.1103/PhysRevA.75.035802（arXiv:quant-ph/0610116 の本文で Eq. (12)–(14) を確認）。
- Thorlabs, "PDB440A/C, PDB450A/C Balanced Amplified Photodetectors Operation Manual," v2.4, 08-Apr-2019, M0009-510-1043（`refs/manuals/thorlabs-pdb.pdf`）。§3.2–3.7、§5.1、§5.1.2 と図。
- Thorlabs, "DSC1 Compact Digital Servo Controller" manual, Rev. B, November 3, 2025（`refs/manuals/thorlabs-dsc1.pdf`）。仕様表（入力終端 1 MΩ、入力 ±4 V、入力 RMS 雑音 0.3 mV）。

abstract・書誌だけ確認したもの（本文未読）:
- H. Hansen, T. Aichele, C. Hettich, P. Lodahl, A. I. Lvovsky, J. Mlynek, S. Schiller, "Ultrasensitive pulsed, balanced homodyne detector: application to time-domain quantum measurements," Opt. Lett. 26, 1714–1716 (2001). DOI 10.1364/OL.26.001714（arXiv:quant-ph/0104084 の abstract から数値）。
- R. Kumar, E. Barrios, A. MacRae, E. Cairns, E. H. Huntington, A. I. Lvovsky, "Versatile wideband balanced detector for quantum optical homodyne tomography," Opt. Commun. 285, 5259–5267 (2012). DOI 10.1016/j.optcom.2012.07.103（著者の並びは検索結果と記憶による。UNVERIFIED）。
- T. Kouadou et al., "Homodyne detection for pulse-by-pulse squeezing measurements," arXiv:2511.04578 (2025)、査読前。
- 76 MHz のパルス分解ホモダイン測定、Opt. Lett. 33, 1458 (2008), DOI 10.1364/OL.33.001458。**著者（Okubo, Hirano ら？）と正確な題名は UNVERIFIED**。
- B. Tohermes, S. Verclas, R. Schnabel, "Directly measured squeeze factors over GHz bandwidth from monolithic PPKTP resonators," New J. Phys. 27, 094506 (2024)（Ref C の文献 [20]、Ref C の BHD の元。書誌は Ref C の文献リストのまま、未読）。

確認できなかったもの:
- Thorlabs PNH1550R5A2 の分岐比・過剰損失（製品ページは 2026-10-28 に取得できなかった）。`R5` が 50:50 を意味するかも UNVERIFIED。T12 で確認する。
