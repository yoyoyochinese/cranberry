# T04 位相雑音 → 測定される squeezing（phase-noise-to-squeezing）

> **訂正（2026-10-31、DECISIONS D65-2）:** ESA のログ平均のトレースのばらつき 0.093–0.099 dB は独立な標本でだけ成り立つ。RBW 1 MHz・VBW 100 Hz では矩形 0.0834–0.0886 dB、Gaussian 0.0706–0.0750 dB（パワー平均は 0.0725–0.0770 / 0.0609–0.0647 dB）。C24 は論文の側の食い違い（D65-6）。


- Phase: 1（deep-dive）／トピック T04（survey の E63–E65、D12・D31 の指標、D16 の等高線、prompt §7-1・§7-2）
- 日付: 2026-10-28
- 規約: `notes/conventions.md`（以下「conv §n」）に従う。文献の規約との対応は §3.5。
- 数値の印: **[Ref A]** などはリファレンス本文の値、**[calc]** はこのノートで計算した値（スクリプトは §12 のテストとして再現できる形で書いた）、**[placeholder]** は Open の値の代用、**UNVERIFIED** は原典で確認していないもの。

## 1. 結論（要約）

1. **準静的な位相揺らぎの写像は 1 本の式にまとまる。** 測定角 Δ = φ_LO − θ_sq（conv §5）の揺らぎ δ が squeezing 側の分散と独立なら、
   V′ = (V₋ + V₊)/2 − (V₊ − V₋)/2 · cos(2Δ₀) · χ、χ = ⟨cos 2δ⟩。
   ガウス雑音なら χ = e^{−2σ²}（conv §6 の厳密式）、振幅 m の正弦波 dither なら χ = J₀(2m)、独立な成分は χ の積。
   §7-1 の小角近似式との差は σ ≤ 0.1 rad で 0.03 dB 未満（Ref A の光源、[calc]）で、**厳密式と小角近似の違いは実用上は小さい**。効くのは「どの帯域を σ² に入れるか」と「ESA の平均のしかた」のほう。
2. **RBW は位相雑音の帯域を切らない。** 測定周波数 f_s から f だけ離れた位相雑音は、f_s ± f の反 squeezing 雑音を f_s に混ぜ込む（側帯波の混合）。光源の squeezing スペクトルが f_s ± f で平らなら（導波路 OPA は THz 帯域）、効き方は遅い揺らぎと**同じ**。σ² の積分範囲は「1/T（統計を取る時間）から、光源のスペクトルが平らでなくなる周波数まで」で、RBW・VBW・検出器の帯域では決まらない。平らでない場合（OPO）は §3.2 の重み付き畳み込みを使う。
3. **CW の ESA 測定とパルスごとの GBS 入力の違い**は上限ではなく、(a) 下限 f_low（ESA はトレースの長さ、GBS は較正の間隔）、(b) 平均のしかた（ESA の log 平均では RBW より遅い揺らぎが dB 領域で平均され、劣化が小さく見える）にある。Ref A の光源では log 平均による過小評価は 9 mrad で 0.03 dB、30 mrad で 0.8 dB [calc]。
4. **Ref C の dither**: ±4.5° の正弦波 → 0.0555 rad rms、J₀(2m) = 0.99384 はガウス近似 e^{−2σ²} = 0.99385 と一致。残り √(0.1² − 0.0555²) ≈ 0.083 rad が dither 以外。Ref C の光源（損失 21%、反 squeezing 10 dB）では位相雑音なしで −5.63 dB、dither だけで −5.18 dB、全体 0.1 rad で −4.32 dB [calc]。**Ref B の 295 kHz dither は probe にかかるので、ループが 295 kHz に応答しない限り θ には入らない**（survey E65 の記述を修正する提案、§15）。
5. **ESA の放射計雑音**: RBW 1 MHz・VBW 100 Hz のトレースの標準偏差は 0.07–0.10 dB（power / log 平均）[calc, AN 1303 の式]。Ref A が報告している「標準偏差 0.07 dB」はこれとほぼ同じで、**位相雑音の情報はほぼ含まない**。また Ref A Fig. 3 の位相掃引トレースは、1 極の VBW 100 Hz と読み取った掃引速度（π rad / 44 ms）では極小が +1.8 dB（power）/ −3.2 dB（log）まで埋まるはずなのに、図では約 −10 dB に達している。**原因は未解決**（§5.4）。
6. **D16 の等高線**: ポンプを最適にしたときの必要 θ_rms は閉じた式 σ_req = ½·√(−ln[1 − ((V_t − L)/(1 − L))²])。目標 10 dB では全損失 L < 10% が必要条件で、L = 3% で 36 mrad、8% で 10.9 mrad、9% で 5.5 mrad [calc]。将来の構成では、63 m の SMF の損失（0.3%）より固定損失と振り分け素子の損失が効く。指標 (ii) の成分 (3)（時間多重によるレーザー位相の相関の崩れ）は、白色 FM・遅延 31·T_rep に対して σ² = 2πΔν·31T_rep。f_rep = 100 MHz なら、σ ≤ 10 mrad に Δν ≤ 51 Hz、30 mrad に 460 Hz が必要 [calc]。**この線幅への要求はかなり厳しい**。

## 2. 記号（symbols, English）

- `q̂, p̂`: quadratures of the squeezed mode in its own squeezing frame (q̂ = squeezed axis), vacuum variance 1/2 (conv §2).
- `θ_sq`: squeezing angle (conv §3). `φ_LO`: LO phase (conv §4–5). `Δ = φ_LO − θ_sq`: measured-quadrature angle relative to the squeezing angle (conv §5, §10 (i)).
- `Δ_set`: intended set point (0 = squeezing, π/2 = anti-squeezing). `δ0`: static (mean) lock-point offset. `Δ₀ = Δ_set + δ0`.
- `δ(t)`: zero-mean fluctuation of Δ (rad). `σ² = Var(δ)`; `σ = θ_rms` (Ref A: θ̃).
- `S_δ(f)`: one-sided PSD of δ [rad²/Hz] (conv §7). `C_δ(τ) = ⟨δ(t)δ(t+τ)⟩`: autocovariance.
- `V₋, V₊`: squeezed / anti-squeezed quadrature variance **normalized to shot noise** (V_SN = 1/2 → 1) with σ = 0; Ref A's R∓. `V′`: measured (jitter-averaged) normalized variance; Ref A's R′.
- `V_±(f)`: same as function of sideband (Fourier) frequency f; `S_q(f), S_p(f)`: one-sided PSDs of q, p normalized to shot noise (= V₋(f), V₊(f) for a stationary state with frequency-independent squeezing angle).
- `χ = ⟨cos 2δ⟩`: phase-noise "characteristic factor"; χ = e^{−2σ²} for Gaussian δ, J₀(2m) for δ = m sin(2πf_m t).
- `m`: peak amplitude of a sinusoidal phase modulation [rad]; rms = m/√2. `f_m`: its frequency.
- `J_n`: Bessel function of the first kind.
- `f_s`: ESA centre (measurement) frequency [Hz] (conv §1). `RBW`, `VBW`: resolution / video bandwidth (−3 dB) [Hz]. `ENBW = k_RBW·RBW`: equivalent noise bandwidth of the RBW filter (k_RBW = 1.128 for 4-pole synchronous analog, ≈1.05 FFT, per AN 1303).
- `T`: duration over which the statistic is accumulated [s] (ESA trace/average time, pulsed-BHD record, or GBS calibration interval). `f_low ≈ 1/T`.
- `L`: total effective loss (Ref A sense, `R± = L + (1 − L)e^{±2r}`). `r`: squeeze parameter (Ref A: r = √(αP)). `α`: SHG efficiency [1/W]. `P`: pump power [W].
- `S_t` [dB]: target squeezing (positive number); `V_t = 10^{−S_t/10}`.
- `Θ_k = θ_sq,k − θ_sq,0`: metric (ii) (conv §10, D31). `τ_k`: emission-time offset of port k pulse relative to its GBS arrival; `Δτ_k = τ_k − τ_0`.
- `φ_L`: laser (fundamental) optical phase. `Δν`: Lorentzian FWHM linewidth (white FM, S_ν = Δν/π, conv §7).
- `f_N`: Nyquist frequency of the time-loop sampling. `S_open(f)`: open-loop (unsuppressed) phase-noise PSD above the loop bandwidth.
- `w(f)`: spectral weight of phase noise at offset f (= 1 for a flat source spectrum).
- `N_eff`: effective number of independent samples in a variance estimate.

## 3. 支配方程式

### 3.1 準静的な揺らぎ（δ が q, p の相関時間より遅い）

ホモダイン検出で測る直交位相は x̂_Δ = q̂ cos Δ + p̂ sin Δ（squeezing の枠、conv §2–3 を Δ だけ回したもの）。分散は
V(Δ) = V₋ cos²Δ + V₊ sin²Δ = (V₋ + V₊)/2 − (V₊ − V₋)/2 · cos 2Δ。
Δ = Δ₀ + δ で、δ の分布が対称（⟨sin 2δ⟩ = 0）かつ q, p と独立なら、

**V′ = (V₋ + V₊)/2 − (V₊ − V₋)/2 · cos(2Δ₀) · χ、χ = ⟨cos 2δ⟩**　(1)

- ガウス δ: χ = e^{−2σ²}。Δ₀ = 0 で V′ = (V₋ + V₊)/2 + (V₋ − V₊)/2 · e^{−2σ²}（conv §6 の厳密式と同じ）。
- 小角近似（§7-1、Ref A/B の式）: V′ ≈ V₋ cos²σ + V₊ sin²σ。厳密式との差は ⟨sin²δ⟩ − sin²σ = (2/3)σ⁴ + O(σ⁶)。V′ ≥ V₊σ² なので相対誤差は (2/3)σ² 以下、σ = 0.1 rad で 0.027 dB、0.03 rad で 0.0014 dB（Ref A の光源、[calc]）。
- 正弦波 δ = m sin(2πf_m t)（時間平均）: χ = J₀(2m)。σ² = m²/2 のガウスとの差は O(m⁴)。
- 独立な成分の和（ガウス成分＋複数の正弦波）: χ = e^{−2σ_G²} · Π_i J₀(2m_i)。
- 分布が任意（ロック外れ、飽和、ドリフト）: χ = 時間平均 ⟨cos 2δ(t)⟩ を時系列から直接取る。**実装の基本はこれ**（§13）。
- 最適ポンプ（Ref A の光源 R± = L + (1 − L)e^{±2r}、Δ₀ = 0）:
  - 小角近似: e^{2r_opt} = cot σ、P_opt = (ln cot σ)²/(4α)、V′_min = L + (1 − L) sin 2σ。
  - 厳密式: tanh 2r_opt = e^{−2σ²}、V′_min = L + (1 − L)·√(1 − e^{−4σ²})。
- 劣化の目安: 10 log10(V′/V₋) ≈ 10 log10[1 + σ²(V₊ − V₋)/V₋]。損失への換算: ΔL_eq ≈ σ²V₊/(1 − V₋)。Ref A の 10 mrad では ΔL_eq ≈ 1.25%（Ref A §1 の「約 1%」と整合、[calc]）。

### 3.2 一般の場合: 側帯波の混合（δ が速いとき）

δ(t) が q(t), p(t)（squeezing の枠での場の直交位相の時系列）と独立で定常なガウス過程、q と p が無相関（squeezing 角が周波数によらない）とする。測定される時系列 x(t) = q(t) cos δ(t) + p(t) sin δ(t) の自己相関は、ガウス δ について ⟨cos δ₁ cos δ₂⟩ = e^{−σ²} cosh C_δ(τ)、⟨sin δ₁ sin δ₂⟩ = e^{−σ²} sinh C_δ(τ) を使うと

R_x(τ) = e^{−σ²}[R_q(τ) cosh C_δ(τ) + R_p(τ) sinh C_δ(τ)]　(2)

（このノートでの導出。δ が準静的、すなわち C_δ(τ) ≈ σ² なら式 (1) に戻る。）
C_δ について 1 次まで展開すると、Δ_set = 0 のときの片側スペクトルは

**S_x(f_s) ≈ S_q(f_s)(1 − σ²) + ∫₀^∞ S_δ(f) · [S_p(f_s + f) + S_p(|f_s − f|)]/2 df**　(3)

- S_p が f_s ± f で平らなら、第 2 項は σ²S_p(f_s) になり、式 (1) の小角近似と一致する。**RBW も検出器の帯域も入らない**（混合は光の段階・積の段階で起き、そのあとに帯域制限がかかるため）。
- 平らでない場合の重み: w(f) = [S_p(f_s + f) + S_p(|f_s − f|)]/(2 S_p(f_s))。Ref B の OPO（γ/2π = c(T + L)/(2πl) ≈ 12 MHz、[calc]、2π の扱いは T02 で確認）では、f が 10 MHz を超える位相雑音の重みが下がる。Ref C（線幅 1.75 GHz）や導波路 OPA（THz）では、関係する f の範囲で w ≈ 1。
- 正弦波 dither（ガウスでない）は厳密に扱える: cos δ = Σ_{n 偶数} J_n(m) e^{inω_m t}、sin δ = −i Σ_{n 奇数} J_n(m) e^{inω_m t} より、両側スペクトルで
  S_x(f) = Σ_{n 偶数} J_n²(m) S_q(f − n f_m) + Σ_{n 奇数} J_n²(m) S_p(f − n f_m)（周期定常の交差項は VBW ≪ f_m の平均で消える）。
  平らなスペクトルでは Σ_{n 偶数} J_n² = (1 + J₀(2m))/2 になり、式 (1) の χ = J₀(2m) に戻る（数値で確認済み、[calc]）。McKenzie Eq. (6) はこの J₁² の項（= 混合の項）を落とした形。
- 適用外: 周波数によって squeezing 角が回る場合（OPO の離調、分散。T02・T03・T06）は q と p が相関する。その場合は 2×2 のスペクトル共分散行列を回す一般形が必要で、Phase 3 で扱う。

### 3.3 有限時間の窓（下限 f_low）

長さ T の区間で平均（または較正）を取り、区間ごとの平均位相を差し引くときの δ の分散は

σ²(T) = ∫₀^∞ S_δ(f)·[1 − sinc²(πfT)] df（sinc x = sin x / x）　(4)

f ≪ 1/T の成分は区間の平均のずれ（ESA ならトレースごとのオフセット）として現れ、σ² には入らない。複数の区間をまとめて報告する場合は、区間の平均の散らばりが別に入る。

### 3.4 どの帯域が θ_rms に入るか（全トピック共通の写像）

| 測定 | f_low | 上限 | 平均のしかた |
|---|---|---|---|
| CW BHD ＋ ESA ゼロスパン（Ref A/B/C、metric (i)） | 1/T_trace。複数のトレースを平均するなら 1/T_total（式 4 の窓） | 光源のスペクトルが f_s ± f で平らでなくなる周波数（式 3 の w(f)）。RBW・VBW・PD の帯域では決まらない | f < RBW: ESA の平均の種類に従う（log なら dB 領域、§5.2）。f > RBW: スペクトルの混合なので power 領域 |
| パルス BHD（パルスごとの積分、将来の metric (i)） | 1/T_record | 約 1/τ_p（これより速い位相の変化はパルス内の時間モードのゆがみ → T03） | サンプルの分散（power 領域） |
| GBS 入力（metric (ii)） | 1/T_cal（位相の較正の間隔、または 1 組のサンプルを取る時間） | 約 1/τ_p | ショットの集団平均（power 領域） |

- パルスごとのサンプリングでは f_rep/2 より上の成分が折り返すが、分散は変わらない（エイリアスはスペクトルを動かすだけ）。
- 時間ループ段は δ(t) を f_N までしか持たない。f_N より上の位相雑音（GAWBS、レーザーの MHz 帯の雑音など）は、ループ帯域の外なので開ループの PSD のまま解析的に足す: χ_HF = exp(−2∫_{f_N}^∞ S_open(f) w(f) df)。

### 3.5 文献の規約からの変換

- **Ref A**: R± は shot noise で規格化した線形の分散 = このノートの V±。θ̃ は「LO と（反）squeezing の直交位相の相対角の標準偏差、小角、ガウス」= σ。R± = L + (1 − L)exp(±2√(αP)) は r = √(αP)。α = 906 %/W = 9.06 W⁻¹。
- **Ref B** Eq. (1) は Ref A と同じ形。θ̃ は rms（度で書かれている）。Eq. (2) の x・G の定義の問題は T02 で扱う。T04 で使うのは Eq. (1) と、本文の「θ̃ = 3.9° なら損失なしでも −8.66 dB、1.5° なら −12.81 dB」。これは V′_min = sin 2σ（損失 0、ポンプ最適）と一致する（§9.1）。
- **McKenzie 2005**: θ = θ₀ + θ₁ sin Ωt。θ₀ は平均の位相差（このノートの Δ₀ だが、McKenzie は振幅 squeezing を基準にしていて、Ṽ(1)・Ṽ(2) の割り当てがこのノートの V₋, V₊ と入れ替わる場合がある）。θ₁ は変調の深さ（= m）。J₀(θ₁)² の因子は式 (1) の χ の小振幅形に当たる。
- **Ref C**: 「±4.5°」は peak 振幅 m = 0.0785 rad。rms = m/√2 = 0.0555 rad（prompt §7-10 の 0.056 rad と同じ）。「位相雑音 約 0.1 rad」は当てはめのモデル（Ref C の文献 [2,6]）での rms で、ガウスの小角近似と解釈する（**UNVERIFIED**: 当てはめに使った式は本文にない）。
- **符号**: Δ = φ_LO − θ_sq（conv §5）。V′ は δ の符号によらないので、写像そのものは符号の取り方に影響されない。ただし式 (1) の Δ₀（オフセット）と誤差信号の向き（T14）は conv に従う。

## 4. 指標 (i)・(ii) の運用上の定義（D12、D31、conv §10）

位相はすべて conv §9 の区間の和で計算する。θ_sq,k(t) = φ_p(t_gen,k)/2 + π/2 + Σ_{seg ∈ sqz path k} ψ_seg(t_pass,seg)（conv §4。区間を通ると A → A e^{iψ} なので squeezing 角も +ψ 回る）。ポンプがレーザーの SHG なら φ_p = 2φ_L + const なので、θ_sq,k = φ_L(t_gen,k) + Σψ + const。

**(i) BHD で測る値**
- δ_i,k(t) = wrap_π(φ_LO,k(t) − θ_sq,k(t) − Δ_set)。wrap_π は (−π/2, π/2] への折り返し（直交位相の角度は mod π）。
- 出力: δ_i,k(t) の時系列と PSD、σ_i² = σ²(T)（式 4）＋高周波の上乗せ（§3.4）、V′_i(f_s)（式 1 ＋ χ_HF）、ESA の模擬トレース（§5）。
- 使う V±: BHD までの損失をすべて含む値（O の表、T01/T02/T12/T13）。

**(ii) GBS 入力での値**
- Θ_k(t) = wrap_π(θ_sq,k(t) − θ_sq,0(t) − Θ_k,set)。ポート k と 0 の、**同じ GBS ショットに到着するパルス**で評価する（生成時刻は τ_k だけ違う）。
- GBS の光子計数は、全ポートに共通の位相回転では変わらない（受動線形回路は共通の位相と交換する）。したがって意味があるのはポート間の相対角だけで、ポート 0 を基準にした (ii) はこの性質を満たす。
- 分解（D31）: 雑音源のオン・オフを切り替えて走らせ、線形の重ね合わせで分ける（アクチュエータが飽和せず、誤差信号が線形の範囲にあるときに限る。§12 で確認）。
  1. 共通点の間の光路の雑音だけをオン（ループはオン）→ ロックの残差 σ²_(1)
  2. 共通点そのもの（probe を分ける光路、参照光を配る光路）の雑音だけをオン → σ²_(2)
  3. レーザー位相雑音 φ_L だけをオン → σ²_(3)。解析式 σ²_(3) = ∫ 4 sin²(πfΔτ_k) S_φL(f) df（conv §7）。白色 FM で |Δτ_k| ≪ 1/Δν なら σ²_(3) = 2πΔν|Δτ_k|（[calc] で積分を確認）。これを**線幅への要求** Δν ≤ σ²_(3),budget/(2π|Δτ_k|) として出す。
- 実効的な squeezing（GBS 入力の proxy）: V′_eff,k = 式 (1)（χ は Θ_k から）。V± は GBS 入力までの損失（光源、ファイバ、振り分け素子。BHD の QE・モード整合・回路雑音は含めない）で計算する。各ショットの状態はどのモードも純粋な squeezed 状態のままなので、これは「ポート 0 の枠で見た、ショット集団の 2 次モーメント」という proxy である（§15 の未決事項 Q1）。

## 5. ESA の模擬（M 段）

### 5.1 ゼロスパン（Ref A Fig. 3、Ref B Fig. 2、Ref C Fig. 2・4）

ESA のチェーン: 中心 f_s の RBW フィルタ → 包絡線 |z| → （log の場合）対数 → VBW（1 極）→ 検出器（sample / average）→ 表示点 → トレース平均。実装は複素ベースバンドで行う:
1. 準静的な分散の時系列 V′(t) を作る（ループの δ(t) の速い成分は χ_HF にまとめて、V′(t) = V(Δ(t)) を χ_HF で補正する。§3.4）。
2. 複素ガウス白色雑音 → RBW の形のフィルタ（ガウス型、−3 dB 幅 = RBW、ENBW 係数はパラメータ）→ √V′(t) を掛ける → P(t) = |z|²。サンプリングは RBW の 4 倍程度で足りる（Ref A: 4 MS/s × 0.15 s = 6×10⁵ 点）。
3. 平均の種類: power（P に VBW）、log（10 log10 P に VBW、Gaussian 雑音では −2.51 dB 低く出る）、voltage（|z| に VBW、−1.05 dB）（AN 1303）。
4. shot noise のトレースを同じ設定で作り、dB で差を取って規格化する。

### 5.2 log 平均のバイアス

- 定常のガウス雑音では log 平均は power 平均より 2.51 dB 低い（AN 1303）。shot noise と squeezed noise で分布の形が同じなら、規格化で打ち消される。
- 位相の揺らぎが RBW より遅いと、パワーの分布が平均の違う指数分布の混合になる。log 平均は ⟨10 log10 V(δ)⟩、power 平均は 10 log10⟨V(δ)⟩ を返し、Jensen の不等式で log のほうが squeezing を大きく見せる。Ref A の光源（640 mW）、ガウスの遅い δ [calc]:

| σ [mrad] | power 平均 [dB] | log 平均 [dB] | 差 [dB] |
|---|---|---|---|
| 9 | −10.149 | −10.181 | 0.033 |
| 20 | −8.76 | −9.10 | 0.34 |
| 30 | −7.22 | −8.02 | 0.80 |

- RBW より速い揺らぎは式 (3) の混合でガウス雑音の分散そのものを変えるので、平均の種類によらず power 領域で効く。
- Ref A・B・C のどれも ESA の平均の種類（log / power）を書いていない（**UNVERIFIED**、パラメータとして持つ）。

### 5.3 放射計雑音（トレースの揺らぎ）

1 極の VBW では、トレースの相対揺らぎは power 平均で σ_P/P = √(π·VBW/ENBW)、log 平均では分散がその 1.64 倍（AN 1303: log では σ = 9.3 dB·√(VBW/RBW)）。[calc]:

| 構成 | RBW / VBW | power 平均 | log 平均 | 報告値 |
|---|---|---|---|---|
| Ref A | 1 MHz / 100 Hz | 0.072–0.077 dB | 0.093–0.099 dB | 標準偏差 0.07 dB [Ref A §3.2] |
| Ref B | 30 kHz / 300 Hz | 0.73–0.77 dB（20 回平均で 0.17 dB） | 0.93–0.99 dB | ±0.14 dB（20 回平均）[Ref B §4] |
| Ref C | 300 kHz / 300 Hz | 0.23–0.24 dB | 0.29–0.31 dB | 図のみ |

（幅は ENBW 係数 1.0–1.128 の範囲。）Ref A・B の報告値の揺らぎは、ほぼ放射計雑音で説明できる。トレースの標準偏差から θ_rms を推定してはいけない。トレースの平均値の不確かさは σ_trace/√N_eff、N_eff ≈ π·VBW·T_trace（1 極の VBW の ENBW が (π/2)VBW のため。Ref A で約 47 → 0.011 dB）。

### 5.4 位相を掃引したトレース（Ref A Fig. 3）

- Δ(t) = ω_scan·t ＋ 雑音（ロックはオフ、LO 経路のアクチュエータにランプ）。V(Δ) は Δ = 0 で鋭い極小を持つ。V = 2V₋ になる半幅は約 √(V₋/V₊) ≈ 29 mrad、0 dB までが約 89 mrad（Ref A の値、[calc]）。
- 図から読み取ると、極小の間隔は約 44 ms（π rad）なので ω_scan ≈ 71 rad/s。極小の付近（0 dB より下）にいる時間は約 2.5 ms で、1 極の VBW 100 Hz の時定数 1.6 ms と同じ程度になる。
- [calc] では、極小は VBW 100 Hz で +1.76 dB（power）/ −3.19 dB（log）まで埋まり、−10 dB に達するには VBW ≳ 1 kHz が必要。極大（20.5 dB）は VBW に鈍感。
- 一方、Ref A Fig. 3 の極小は約 −10 dB に達している。**食い違いの原因は未解決**（可能性: 掃引トレースだけ別の VBW や検出器、掃引速度の読み違い、Keysight N9010B の VBW の実装が 1 極ではない）。
- 掃引トレースの検証は「ロック時の値と極大」と「VBW に対する極小の深さの単調性」に限り、極小の深さの一致は要求しない（§12 V9）。

### 5.5 スペクトル（Ref A Fig. 4）

- 各周波数ビンで S(f) = ∫ H_RBW(f − f′)·[G_det(f′)·(V′(f′)·S_SN) + S_circ(f′)] df′ ＋ 線スペクトル（probe のビート 1 MHz など、T14）。V′(f) は式 (3) を f ごとに使う。
- 規格化: (a) shot noise で規格化した値、(b) 回路雑音を差し引いた値 N(f) = (S − S_circ)/(S_SN − S_circ)（Ref A Fig. 4(b)）。
- 各ビンの揺らぎは §5.3 の放射計雑音（滞在時間 = 掃引時間/点数 と VBW で決まる）。
- 検出器の応答 G_det と回路雑音 S_circ は T13 から受け取る。

## 6. 仮定と適用範囲

- δ は squeezing 側の量子雑音と独立（古典的な位相雑音）。量子雑音ロックのセンサ雑音（T14）が δ に入るのは、ループを通った後の古典的な信号としてだけ。
- 式 (1) は δ がガウスまたは既知の分布（正弦波、経験分布）で、q, p の相関時間より遅い場合に厳密。速い成分は式 (3)（1 次の展開、σ² ≪ 1）か、正弦波の Bessel 和で扱う。
- q と p が無相関（squeezing 角が f によらない）。そうでない場合は §3.2 の最後を参照。
- 単一モード（CW では各サイドバンド周波数ごと、パルスでは LO の時間モード 1 つ）。多モードの効果は T03。
- 位相に依存しない損失は V± に含まれている（T01/T02/T12/T13）。純損失は位相回転と交換するので、位相雑音が損失の前に入っても後に入っても V′ は同じ。

## 7. 周波数帯

- 効く位相雑音: f_low = 1/T（式 4）から、光源のスペクトルが平らな範囲の上限まで（§3.4）。Ref A: 約 7 Hz（T_trace = 0.15 s）〜 THz（OPA の帯域）。ただし実際の雑音源は DC〜GHz（ファイバ環境 Hz〜kHz、GAWBS MHz〜GHz、dither 1 MHz）。
- dither: Ref C 1 MHz（LO 上、全量が効く）。Ref B 295 kHz（probe 上、ループの 295 kHz での応答分だけ）。
- ESA: f_s（Ref A 3 MHz、B 1 MHz、C 30 MHz–1.21 GHz）、RBW・VBW は測定の統計と log バイアスにだけ効く。

## 8. スケーリング

- ポンプパワー: 劣化は σ²·V₊ ∝ σ²·e^{2r}。Ref A の形では r = √(αP) なので、P に対して指数的に厳しくなる。最適ポンプ P_opt = (ln cot σ)²/(4α) が存在する（Ref A Fig. 5、Ref B Fig. 3 の x_opt）。
- ファイバ長: 写像そのものには V± を通してだけ入る（損失が増えると V₋ が上がり、必要な σ は緩くなるが、達成できる squeezing は下がる。§10）。σ²（入力）は T05・T08 で L とともに増える（環境雑音 ∝ L、遅延干渉 ∝ τ）。
- 温度: 直接は入らない（T08 の PSD を通して）。
- パルス幅: τ_p は帯域の上限 ~1/τ_p を決めるだけ。実際の雑音源はそれよりずっと遅いので、影響はない。パルス内の位相の変化（SPM など）は T03/T06。
- 遅延（指標 (ii) の成分 3）: σ²_(3) ∝ Δν·Δτ_k ∝ k·T_rep（白色 FM の場合）。

## 9. 概算

### 9.1 リファレンスの検算 [calc]

| 条件 | R′₋（小角） | R′₋（厳密） | R′₊ | P_opt | 最適での R′₋ |
|---|---|---|---|---|---|
| Ref A: L = 8%、θ = 9 mrad、640 mW | −10.148 dB | −10.148 dB | +20.556 dB | 612 mW | −10.152 dB |
| L = 12%、θ = 14 mrad | −8.285 | −8.285 | +20.364 | 503 mW | −8.397 |
| L = 8%、θ = 14 mrad | −9.597 | −9.598 | +20.556 | 503 mW | −9.757 |
| L = 8%、θ = 0 | −10.582 | −10.582 | +20.556 | ∞ | → −10.969（= 10 log10 L） |

- 感度（Ref A、640 mW）: θ が +5% → R′₋ は +0.042 dB、+10% → +0.086 dB。L が ±1% → R′₋ は約 ±0.45 dB。P が ±40 mW → R′₊ は ±0.65 dB、R′₋ は 0.02 dB 以内。
- Ref A の測定値の反 squeezing 20.2 ± 0.2 dB に対し、モデルは 20.56 dB（0.36 dB 高い）。論文の当てはめは squeezing 側に強く引かれていると見られ、モデルと測定のずれとして記録する（合わせ込まない）。
- Ref B: 損失 0・ポンプ最適の上限 sin 2σ: 3.9° → −8.67 dB（論文 −8.66）、1.5° → −12.81 dB（論文 −12.81）。厳密式では −8.68 / −12.81 dB。
- Ref C: 損失 21%、反 squeezing 10 dB（30 MHz）から e^{2r} = 12.4 とおくと、σ = 0 で −5.63 dB、dither のみ（0.0555 rad）で −5.18 dB、0.1 rad で −4.32 dB。測定値 −4.7 dB（暗雑音補正後）と同程度（Fig. 5 の当てはめは別のポンプパワーなので、厳密な比較ではない）。

### 9.2 Ref D（D02）

Ref D は squeezed 光を使わない MZI の DC fringe lock なので、T04 の出力は「残留位相 → 仮に Ref A の光源（640 mW、L = 8%）を測ったときの squeezing」の換算表になる。ブラインドの値（D04）は使わない。

| θ_rms | 1 mrad | 3 mrad | 10 mrad | 30 mrad | 100 mrad |
|---|---|---|---|---|---|
| R′₋ [dB] | −10.58 | −10.53 | −10.05 | −7.22 | +0.83 |

Ref D の θ_rms の帯域: f_low = 1/T_record、上限は開ループの干渉計の雑音（T05 の遅延干渉、T08）がある範囲まで。ループ帯域（DSC1 で数十 kHz 以下、T15/T16）より上は抑えられずにそのまま入る。

### 9.3 将来の 32 ポート構成（D07–D13、D30–D32）

値はすべて代用で、f_rep・振り分け素子の損失・PM の損失は Open。
- ポート長 L_k = k·c·T_rep/n_g（survey §3）: f_rep = 100 MHz で L_31 = 63 m、10 MHz で 630 m。
- 固定損失 L₀ の代用: 3%（Ref A の導波路 2% ＋光学部品 1%。BHD の損失は GBS 入力の指標から除く）と 8%（Ref A の全体）。振り分け素子の損失は 0 とした（**placeholder**、実際は効く）。
- ファイバ損失: SMF 0.2 dB/km（Ref C）、PM 1 dB/km（**placeholder、UNVERIFIED**、T09）。

必要な θ_rms [mrad]（ポンプ最適、厳密式、[calc]）:

| L₀ | ファイバ | L_k | 全損失 | 3 dB | 5 dB | 7 dB | 10 dB |
|---|---|---|---|---|---|---|---|
| 3% | SMF | 63 m | 3.3% | 259 | 150 | 87 | 35 |
| 3% | SMF | 630 m | 5.8% | 250 | 140 | 76 | 22 |
| 3% | PM | 63 m | 4.4% | 255 | 145 | 82 | 29 |
| 3% | PM | 630 m | 16.1% | 212 | 93 | 23 | 不可能 |
| 8% | SMF | 63 m | 8.3% | 242 | 129 | 64 | 9.4 |
| 8% | PM | 630 m | 20.4% | 194 | 71 | 不可能 | 不可能 |

- r を固定した場合（GBS の要求で r が決まる場合。Liu 2026 は r = 0.9–1.8）: r = 1.8、L = 5% で 3/5/7/10 dB に 111/84/60/26 mrad。r = 1.2 では 10 dB に届かない（[calc]）。ポンプ最適の表は上限（最も緩い要求）である。
- 成分 (3) の線幅への要求（白色 FM、Δτ = 31·T_rep、[calc]）:

| f_rep | Δτ_31 | σ_(3) ≤ 10 mrad | ≤ 30 mrad | ≤ 100 mrad |
|---|---|---|---|---|
| 10 MHz | 3.1 µs | Δν ≤ 5.1 Hz | 46 Hz | 513 Hz |
| 20 MHz | 1.55 µs | 10 Hz | 92 Hz | 1.0 kHz |
| 100 MHz | 310 ns | 51 Hz | 462 Hz | 5.1 kHz |
| 1 GHz | 31 ns | 513 Hz | 4.6 kHz | 51 kHz |

  ここでの Δν は、パルス間の光の位相の相関を決める量（モード同期レーザーなら、キャリアに近い櫛の線の線幅）である。1/f 雑音や技術的なピークがある場合は、T05 の S_ν を §4 の成分 3 の式で積分する。
- 比較: Liu 2026 は GBS 全体の位相を λ/200（1550 nm）= 31 mrad に安定化したと報告している。各ロックの精度は 1 / 3 / 10 nm（4 / 12 / 41 mrad）[Liu, Methods]。

## 10. D16 の等高線「目標 squeezing × ファイバ長 → 必要な θ_rms」の方法

1. 損失のモデル（O）: L_tot(L_f) = 1 − (1 − L₀)(1 − L_sw)·10^{−a_f L_f/10}（a_f [dB/m]）。L₀・L_sw・a_f は T01・T11・T12・T06/T09 から受け取る。指標 (ii) では BHD の損失を含めず、(i) では含める。
2. 必要な θ_rms（解析式、O）:
   - ポンプ最適: s = (V_t − L_tot)/(1 − L_tot)、σ_req = ½·√(−ln(1 − s²))（s ≤ 0 なら実現不可能。境界 L_tot = 10^{−S_t/10}）。
   - r 固定: e^{−2σ²_req} = [cosh 2r − s]/sinh 2r（右辺が (0, 1] の外なら実現不可能）。
   - 一般形（r 固定、dither あり）: 式 (1) から χ_req = [(V₊ + V₋)/2 − V_t]/[(V₊ − V₋)/2]（V± は L_tot と r から）。e^{−2σ²_G,req} = χ_req/Π_i J₀(2m_i) を解く（σ_G の予算が dither の分だけ減る）。
3. 格子: S_t = 3–10 dB を 0.5 dB 刻み（D16）、L_f = 各ポートの L_k（f_rep の関数）または 0–1 km の連続。
4. 達成できる θ_rms(L_f)（T）: 時間ループで指標 (ii)（と (i)）の σ を L_f ごとに出し、§4 の 3 成分に分ける。
5. 結果: σ_req(S_t, L_f) の等高線に、達成値 σ_ach(L_f) の曲線を重ね、σ_ach ≤ σ_req の範囲を「実現可能」とする。予算の配分 σ²_req = σ²_(1) + σ²_(2) + σ²_(3) から、成分 (3) を線幅への要求に、(1) をループの要求（帯域、利得）に戻す。
6. 検算: 数値で V′ を r について最小化した結果と、閉じた式が 10⁻³ dB で一致すること（§12 V2）。

## 11. prompt §7-1 の許容幅の根拠

- **解析の単一モードモデル**（式の実装の検査）: R′₋ = −10.148 dB、R′₊ = +20.556 dB に 10⁻³ dB 以内で一致すること（数値誤差だけ）。小角近似と厳密式の差は 4×10⁻⁴ dB。
- **全体のシミュレーション**（参照の損失収支、ロックを 9 mrad に調整）: |R′₋ − (−10.15 dB)| ≤ 0.10 dB、|R′₊ − 20.56 dB| ≤ 0.05 dB。内訳:
  1. θ_rms の推定の統計誤差: 有限の時系列から σ² を推定する相対誤差 √(2/N_eff)。R′₋ の感度は σ² の 10% で 0.042 dB。N_eff ≥ 800（σ² で 5%）として 0.02 dB。
  2. 損失のモデル（T01 の分布損失 ↔ 一括の L = 8%）: R′₋ は 1% の損失で約 0.45 dB 動く。T01 が報告する実効の L を 8 ± 0.1% に合わせ、0.045 dB。
  3. ESA の平均の種類（§5.2）: 9 mrad で最大 0.033 dB。
  4. 模擬トレースの平均値の放射計雑音（§5.3）: 0.011 dB。
  直線和で約 0.11 dB、二乗和で約 0.06 dB。直線和を丸めて 0.10 dB とする。これは Ref A の測定誤差 ±0.2 dB（3σ）の半分で、報告された標準偏差 0.07 dB と同程度。R′₊ は θ にほぼ依存せず、L と P だけで決まるので 0.05 dB とする。
- **測定値との比較**（参考、合否には使わない）: squeezing 10.1 ± 0.2 dB とは一致する。反 squeezing はモデル自体が測定値より 0.36 dB 高い（§9.1）。
- **§7-2**（Fig. 5 の形）: P_opt と最適値（612 mW・−10.152 dB、503 mW・−8.397 dB）を 1 mW・0.01 dB で再現すること、θ = 0 で単調に 10 log10 L へ近づくこと。

## 12. 検証テスト（実装が通すべきもの）

| # | テスト | 合格条件 |
|---|---|---|
| V1 | §7-1 の解析式 | R′₋ = −10.148 dB、R′₊ = +20.556 dB（±10⁻³ dB）。小角と厳密の差 < 10⁻³ dB |
| V2 | 最適ポンプの閉じた式（§7-2） | (8%, 9 mrad) → 612 mW、−10.152 dB。(12%, 14 mrad) → 503 mW、−8.397 dB。数値最小化と 10⁻³ dB で一致 |
| V3 | Ref B の上限 | 損失 0: 3.9° → −8.67 dB（小角）/ −8.68 dB（厳密）、1.5° → −12.81 dB |
| V4 | 厳密式 ↔ Monte Carlo | ガウスの δ を 10⁶ 個サンプルし、V(Δ) を平均した値が式 (1) と統計誤差内で一致。σ = 0.1 rad で小角近似との差が 0.027 dB（Ref A の光源） |
| V5 | 正弦波 dither（Ref C、§7-10） | m = 4.5° → rms 0.0555 rad、χ = J₀(2m) = 0.99384。Ref C の光源（L = 21%、V₊ = 10 dB）で −5.63 → −5.18 dB。時間領域で x = q cos δ + p sin δ を作り（q, p 白色、f_m = 1 MHz、f_s = 30 MHz、RBW 300 kHz）、PSD が式 (1) と ±0.05 dB で一致 |
| V6 | 側帯波の混合 | 平らなスペクトルで、f = 20 MHz の位相変調（f_s = 3 MHz、RBW 1 MHz）が、同じ rms の遅い変調と同じ劣化を与える（±0.05 dB）。Lorentz 型の S_p（HWHM 12 MHz）では式 (3) の畳み込みと一致 |
| V7 | ESA の放射計雑音と log バイアス | 一定の V で RBW 1 MHz・VBW 100 Hz のトレースの標準偏差が 0.072–0.077 dB（power）/ 0.093–0.099 dB（log）。log 平均の絶対値は power 平均より 2.51 ± 0.05 dB 低い。δ = 0 なら squeezed と shot の比は平均の種類によらない（差 < 0.01 dB） |
| V8 | 遅い揺らぎの log バイアス | Ref A の光源で σ = 9 / 30 mrad（RBW より遅い）: log 平均と power 平均の差が 0.033 / 0.80 dB（±10%） |
| V9 | 位相掃引トレース | ロック時の値と極大（20.5 dB）を再現。極小の深さは VBW を上げると単調に深くなる（100 Hz → +1.8 dB（power）、1 kHz → −9.9 dB、10 kHz → −10.6 dB）。Ref A Fig. 3 の極小の深さとの一致は要求しない（§5.4） |
| V10 | 窓の式 (4) | random walk の位相で、長さ T の窓ごとの標本分散の平均が ∫S[1 − sinc²(πfT)]df と 5% 以内で一致 |
| V11 | 成分 (3) | 白色 FM Δν = 1 kHz、τ = 310 ns → σ² = 2πΔντ = 1.95×10⁻³ rad²（σ = 44 mrad）。4 sin²(πfτ)·S_φ を数値積分した値、および遅延バッファで作った時系列の分散と 2% 以内で一致 |
| V12 | 指標の不変性 | 全ポートに同じ雑音（共通の区間）を足しても Θ_k は変わらず、(i) は変わる。分解した 3 成分の分散の和が、全部オンの分散と 5% 以内で一致（飽和なし） |
| V13 | 高周波の上乗せ | f_N より上の雑音を (a) 解析の χ_HF で足した場合と、(b) 細かい時間刻みで時系列に含めた場合で、V′ が 0.01 dB で一致 |

## 13. シミュレータへの要求

| 部品 | 段 | 入力 | 出力 | フラグ |
|---|---|---|---|---|
| `phase_to_quadrature`（写像） | M（時間ループの後処理） | δ(t)（指標ごと）、V±(f_s) の表、χ_HF、dither のリスト {m_i, f_i}（時系列に含めない場合）、Δ_set、T、平均の種類 | V′(t)、V′、R′± [dB]、σ_rms（帯域ごと：ループ内、HF、dither） | `phase_noise.enable`、`average: exact\|small_angle`（§7-1 の照合用）、`include_hf_phase_noise`、`dither_mode: timeseries\|analytic` |
| `hf_phase_noise_integral` | O | S_open(f)（T05・T07・T08）、V±(f) の表、f_s、f_N | χ_HF（または式 (3) の ΔV_HF(f_s)） | `sideband_weighting: flat\|spectral` |
| `vpm_table`（V±(f)、損失込み） | O | T01/T02（光源）、T12/T13（損失、検出）、指標 (i)/(ii) で損失の集合を切り替え | V₋(f)、V₊(f) | — |
| `esa_zero_span` | M | V′(t)、f_s、RBW、ENBW 係数、VBW、検出器、平均の種類、掃引時間、点数、トレース平均回数、乱数シード | トレース [dB rel. SN]、平均・標準偏差 | `esa.average_type: log\|power\|voltage`、`esa.radiometer_noise` |
| `esa_swept_spectrum` | M | V′(f)（式 3）、G_det(f)、S_circ(f)、線スペクトル（T14）、RBW/VBW、スパン | S(f)、規格化 (a)、回路雑音を引いた (b) | `subtract_circuit_noise` |
| `phase_scan` の設定 | T（構成） | LO 経路のアクチュエータへのランプ、ロックはオフ | Δ(t) → `esa_zero_span` | — |
| `metric_i`・`metric_ii` | T（State から読むプローブ） | 区間ごとの位相（conv §9）、ポートの到着時刻 | δ_i,k(t)、Θ_k(t)、3 成分（雑音源の切り替えで走らせる） | 雑音源ごとのオン・オフ（保護リスト） |
| `required_theta_contour` | O | S_t の格子、L_f の格子、損失モデル、r の扱い（最適 / 固定）、dither | σ_req(S_t, L_f)、実現可能な範囲 | `pump: optimal\|fixed` |

実装の注意:
- 基本は「時系列からの直接平均」χ = ⟨cos 2δ(t)⟩（分布によらず正しい）× χ_HF。ガウスの閉じた式は検証と等高線の計算に使う。
- ループの速度で δ(t) を持つ（prompt §5.4）。dither（1 MHz）を時系列に含める場合は f_N > f_m が必要。含めない場合は解析の J₀ を使い、二重に数えないこと。
- ESA の模擬は RBW の数倍のサンプリングで時間領域でやっても軽い（Ref A のトレース 1 本で 6×10⁵ 点）。放射計雑音の解析式は検証用。

## 14. 他のトピックとの依存（代用した値）

- V±(f): T01（導波路 OPA の分布損失。代用は Ref A の R±）、T02（OPO の Lorentz 型スペクトル、x・G の定義。代用は Ref C の当てはめ値）。
- 損失: T12・T13（代用は Ref A の 8% / 3%、Ref C の 21%）。振り分け素子 T11（代用 0、**placeholder**）。PM の損失 T09（代用 1 dB/km、**UNVERIFIED**）。
- S_δ の入力: T05（レーザー、遅延干渉）、T07（GAWBS などの高周波）、T08（環境）、T10（相殺の崩れ）、T14（センサ雑音、dither、probe の線）、T16（アクチュエータ、ループ）。
- パルス内の位相とモードの重なり: T03。

## 15. DECISIONS・conventions・survey への提案と、未決事項

調整役への提案（このノートでは変更していない）:
1. **conv §6** の「σ² は測定時間の逆数から、測定帯域で効く上限までの積分」に次を追記する: 「上限は RBW や検出器の帯域ではなく、光源の squeezing スペクトルが f_s ± f で平らな範囲（T04 式 3）。下限は T04 式 (4) の窓。RBW より遅い成分は ESA の平均の種類（log / power）に依存する」。また、一般式 V′ = (V₋ + V₊)/2 − (V₊ − V₋)/2·cos(2Δ₀)·χ を conv §6 に載せる（オフセットと dither を含む）。
2. **survey E65**: 「B（295 kHz）」の dither は probe にかかるので、θ に入るのはループの 295 kHz での応答分だけ。Ref C の dither は LO にかかり、全量が θ に入る。
3. **survey E63**: log 平均のバイアス −2.5 dB は AN 1303 で確認できた（2.51 dB）。UNVERIFIED を外してよい。
4. **survey E64**: 「CW の ESA とパルスごとの GBS で帯域が違う」は、上限ではなく下限と平均のしかたの違い（§3.4）と書き直す。
5. **指標 (ii)** の分解（D31）の実行方法として、§4 の「雑音源の切り替えによる重ね合わせ」と線形性の確認（V12）を PLAN に入れる。

ユーザの決定が必要な点（推測していない）:
- **Q1**: GBS 入力での評価量。このノートの既定は「ポート 0 の枠で見たショット集団の 2 次モーメント」による実効 squeezing（proxy）。GBS では各モードは純粋なままで、位相雑音は干渉計の入力位相の誤差として効くので、出力分布の忠実度（total variation distance など）で評価する必要があるかどうか。
- **Q2**: D16 の θ_rms は、Var(Θ_k)（ポート 0 との差、独立なら 1 ポートの分散の 2 倍）を指すのか、1 ポートあたりの σ_port を指すのか。既定では両方を出す。
- **Q3**: D16 の「目標 squeezing」は GBS 入力での値（BHD の損失を含めない）か、BHD で測る値か。既定は GBS 入力（指標 (ii)）で、(i) の等高線も並べて出す。
- **Q4**: 将来の BHD 測定（指標 (i)）は ESA（CW 的な f_s での測定）か、パルスごとの積分か（T13 とも関係する）。§3.4 の f_low と平均のしかたが変わる。
- 参考（ユーザの決定ではなく調査）: Ref A Fig. 3 の掃引トレースの VBW の食い違い（§5.4）。ESA の平均の種類（Ref A・B・C とも記載なし）。

## 16. 文献

一次資料（本文を読んだもの）:
- [Ref A] K. Hirota, T. Kashiwazaki, G. Ha, T. Yamashima, P. Jaturaphagorn, T. Suzuki, K. Takahashi, A. Kawasaki, A. Inoue, W. Asavanant, M. Endo, T. Umeki, A. Furusawa, "Generation of 10-dB squeezed light from a broadband waveguide optical parametric amplifier with improved phase locking method," arXiv:2511.15082 (2025), DOI 10.48550/arXiv.2511.15082。**掲載誌は UNVERIFIED**（arXiv の記録で確認、雑誌の DOI は見つからなかった）。§3.2・§4.1・Fig. 3–5 を使用。
- [Ref B] Y. Takeno, M. Yukawa, H. Yonezawa, A. Furusawa, "Observation of −9 dB quadrature squeezing with improvement of phase stability in homodyne measurement," Opt. Express 15, 4321–4327 (2007)。DOI 10.1364/OE.15.004321（**UNVERIFIED**: 雑誌・巻・頁は PDF で確認、DOI は OSA の規則から推定）。
- [Ref C] S. Verclas, C. Darsow-Fromm, B. Tohermes, R. Schnabel, "Fiber distribution of phase-stabilized GHz-bandwidth squeezed vacuum states of light between two buildings," Optica Quantum 4, 1 (2026)。DOI 10.1364/OPTICAQ.567418（PDF に記載）。
- K. McKenzie, E. E. Mikhailov, K. Goda, P. K. Lam, N. Grosse, M. B. Gray, N. Mavalvala, D. E. McClelland, "Quantum noise locking," J. Opt. B: Quantum Semiclass. Opt. 7, S421–S428 (2005)。DOI 10.1088/1464-4266/7/10/032（PDF に記載）。Eq. (6)–(7) を使用。
- H.-L. Liu et al., "Gaussian boson sampling with 1,024 squeezed states in 8,176 modes," Nature (2026)。DOI 10.1038/s41586-026-10523-6（PDF に記載）。位相安定度 λ/200、1/3/10 nm、r = 0.9–1.8 を使用。
- Hewlett-Packard, "Spectrum Analyzer Measurements and Noise," Application Note 1303, publication 5966-4008E (1998)。https://www.hpmemoryproject.org/an/pdf/an_1303.pdf（閲覧 2026-10-28）。log 平均 2.51 dB、voltage 平均 1.05 dB、ENBW 係数、σ = 9.3 dB·√(VBW/RBW)、log では分散が 1.64 倍、を使用。

書誌だけ確認したもの（本文は未読）:
- S. Suzuki, H. Yonezawa, F. Kannari, M. Sasaki, A. Furusawa, "7 dB quadrature squeezing at 860 nm with periodically poled KTiOPO4," Appl. Phys. Lett. 89, 061116 (2006)。DOI 10.1063/1.2335806（ウェブ検索で確認）。Ref A の文献 [25]、Ref B の文献 [1]。
- A. Franzen, B. Hage, J. DiGuglielmo, J. Fiurášek, R. Schnabel, "Experimental demonstration of continuous variable purification of squeezed states," Phys. Rev. Lett. 97, 150505 (2006)。DOI 10.1103/PhysRevLett.97.150505（**DOI は UNVERIFIED**、書誌はウェブ検索で確認）。Ref C の文献 [22]（位相拡散した squeezed 状態）。
- T. C. Zhang, K. W. Goh, C. W. Chou, P. Lodahl, H. J. Kimble, "Quantum teleportation of light beams," Phys. Rev. A 67, 033802 (2003)、および T. Aoki, G. Takahashi, A. Furusawa, "Squeezing at 946 nm with periodically poled KTiOPO4," Opt. Express 14, 6930 (2006)。Ref A の文献 [26]・[27]。DOI は **UNVERIFIED**。

式 (2)–(4) と正弦波の Bessel 和は、このノートで導出した（文献からの引用ではない）。§12 の V4–V8・V10 は、その導出を数値で確認するためのテストでもある。
