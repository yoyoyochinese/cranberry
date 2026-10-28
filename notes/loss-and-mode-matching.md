# T12 損失とモード整合（loss-and-mode-matching）

- Phase: 1（deep-dive、トピック T12）
- 日付: 2026-10-28
- 担当する効果（`notes/survey.md`）: E28（後方反射とエタロン）、E52（光学部品の損失）、E53（空間モードの不整合）、E55（squeezed 光の経路上のアクチュエータ・変調器の挿入損失）、E56（conventional tap のトレードオフ）、prompt §7-3・§7-8。E45（probe を入れ出しする素子の損失）は T10 と共有し、損失の値はここで出す。
- 入力: `docs/squeezing_sim_prompt.md`、`DECISIONS.md`（D01–D33）、`notes/conventions.md`、Ref A §2・§4、Ref B §4、Ref C §3–4、Liu et al. 2026（Methods）、`refs/manuals/914-915B-Specifications.pdf`、メーカー資料（§9、閲覧日 2026-10-28）。
- 値の印は survey §0.3 と同じ: **[Ref A]** など（本文から読んだ値）、**[manual]**（`refs/manuals/`）、**[web]**（メーカー資料を開いて読んだ値、URL は §9）、**[AE]**（この担当による概算）、**[UNVERIFIED]**（原典を開いて確かめていない値。検索結果の抜粋だけで見たものを含む）。

---

## 0. 要約

1. 損失はすべて**純損失チャネル** V → ηV + (1−η)V_vac として合成できる。直列の要素は透過率の積 η_tot = Π η_i になり、順序にも依らない（受動的な損失のみの場合）。受動ファイバの分布損失は、同じ総透過率の集中損失と**厳密に**同じ。違いが出るのは、損失と同時に利得や雑音の注入がある場合だけ（導波路 OPA の分布損失は T01、ファイバ内の散乱雑音は T07）。
2. モード不整合は、整合効率 η_mm = |⟨u_s|u_LO⟩|² の純損失。縞の可視度 V（同じパワーのとき）とは η_mm = V² の関係。Ref A の「干渉効率 約 99% → 損失 約 2%」も Ref B の 1 − ρζξ²η も、この関係で読めば合う（§1.3）。
3. 3–10 dB の目標（D16）に対する**全損失の上限**は、位相ジッタ θ = 9 mrad、スクイーズ度 r を最適にとったとき、3 dB で 49%、6 dB で 24%、10 dB で **8.4%**（Ref A の全損失 8% とほぼ同じ）。近似式は L_max ≈ (10^(−S/10) − 2θ)/(1 − 2θ)（§4.1）。
4. 将来の 32 ポート構成で、ポート上流（導波路 2%、光学部品 1%、ファイバ結合 2% の仮の値）だけで約 4.9% を使う。10 dB を GBS 入力で出すには、**ポートごとの部品（振り分け素子、接続、ストレッチャ、probe の入れ出し、遅延ファイバ）の合計を 0.16 dB 以下**に抑える必要がある。FC/APC の嵌合 1 対（0.1–0.3 dB）が 1 個入るだけで超える。ファイバ型 LN 変調器の振り分け（2.7 dB typ.）を 1 段置くだけで、上流と合わせた損失は 49%、squeezing の上限は 3.1 dB になる。コネクタとファイバ型 PBS も含めたすべてファイバ部品の想定では上限 約 1.6 dB（§4.2）。
5. 遅延ファイバの減衰は、f_rep = 100 MHz（L_31 = 63 m）なら PM で 0.73%、SMF で 0.26% と小さい。ポート間の差が問題になるのは f_rep ≲ 20 MHz（L_31 ≳ 300 m）から（§4.3）。PM と SMF の差（D32）は、この長さになるまでは損失の面からは決め手にならない。
6. 後方反射による寄生エタロンは、FC/APC 同士では位相のリップルが 10⁻⁶ rad 程度で無視できる。反射率 −40 dB 級の素子（LN 変調器の ORL）と組み合わせても 10⁻⁵ rad 程度（§5）。
7. conventional tap では、tap 比 r_t の損失と、tap で取り出した probe の SNR の間で最適な点がある。tap 用の PD の電気雑音を入れると、例として P_pr = 10 µW、NEP = 1 pW/√Hz のとき最適な r_t ≈ 3×10⁻³ で −9.94 dB（tap なしの −10.15 dB より 0.2 dB 悪い）。phase-detection OPA では r_t に依存しない（§6）。Liu 型（D30）では、tap 比の代わりに probe を入れ出しする素子の挿入損失が、probe パワーと無関係な固定の損失として残る。

---

## 1. 支配方程式

### 1.1 記号（symbols）

- `η` (eta): power transmission (efficiency) of a passive element or channel, 0 ≤ η ≤ 1. `L = 1 − η`: loss.
- `η_tot`: total transmission of a cascade. `L_tot = 1 − η_tot`.
- `IL_dB`: insertion loss in dB, `η = 10^(−IL_dB/10)`.
- `V`, `V_in`, `V_out`: variance of a quadrature x̂_θ (conventions §2, `V_vac = V_SN = 1/2`).
- `R = V / V_SN`: variance normalized to shot noise (Ref A/B notation). `R±`: anti-squeezed (+) / squeezed (−).
- `σ`: 2×2 quadrature covariance matrix of a single mode, `σ_ij = ⟨{Δx̂_i, Δx̂_j}⟩/2`, `x̂ = (q̂, p̂)`.
- `r`: squeezing parameter (conventions §3). `θ`: RMS phase error between LO and the (anti-)squeezed quadrature (conventions §6).
- `α_dB`: fiber attenuation [dB/km]. `α = α_dB · ln(10)/10 / 1000`: attenuation coefficient [1/m]. `L_f`: fiber length [m].
- `ν(z)`: excess-noise injection density along a fiber [variance units per m], measured at position `z`.
- `u_s(x,y)`, `u_LO(x,y)`: normalized transverse mode functions of the squeezed beam and the LO, `∫|u|² dA = 1`.
- `η_mm = |⟨u_s|u_LO⟩|²`: mode-matching (overlap) efficiency.
- `V_vis`: fringe visibility measured between two beams of powers `P_1`, `P_2`.
- `w_1`, `w_2`: Gaussian beam radii (1/e² intensity) at a common waist plane. `d`: lateral offset. `θ_t`: tilt angle. `λ`: vacuum wavelength, `n`: refractive index of the medium.
- `MFD`: mode-field diameter of a fiber, `w = MFD/2`.
- `r_1`, `r_2`: amplitude reflectances of two reflectors; `R_i = |r_i|²`. `ψ_e`: one-way phase across the parasitic etalon.
- `RL_dB`, `ORL`: return loss / optical return loss in dB, `R = 10^(−RL_dB/10)`.
- `r_t`: tap ratio (fraction of power tapped off the squeezed beam). `P_pr`: probe power exiting the squeezer OPA. `P_tap = r_t · P_pr`.
- `S_φ,sens(f)`: one-sided PSD of the sensor-referred phase noise [rad²/Hz] (conventions §7).
- `η_d`: quantum efficiency of the tap photodiode. `NEP`: noise-equivalent power of the tap detector [W/√Hz], one-sided. `h`, `ν_opt = c/λ`: Planck constant, optical frequency.
- `f_BW`: unity-gain (closed-loop) bandwidth of the lock [Hz]. `κ_sn`, `m`: dimensionless constants of the error-signal scheme (defined in §6.1, value from T14).
- `E_el = 10^(−C_dB/10)`: circuit-noise variance relative to shot noise, `C_dB`: clearance [dB].
- `k = 0…31`: port index (D08). `L_k = k·c·T_rep/n_g`: delay-fiber length of port k (D08). `f_rep = 1/T_rep`.
- `N_j`: number of fiber joints (mated connector pairs or splices) per port. `IL_j`: IL per joint.
- `N_inj`: number of probe injection/extraction elements on the squeezed path per port. `IL_inj`: IL per element.
- `IL_dmx(k)`, `N_sw(k)`: IL of the time-demux switch seen by port k; number of switch stages passed.
- `η_up`: transmission upstream of the demux (waveguide effective loss × optics × free-space-to-fiber coupling).
- `η_port(k)`: transmission of all per-port elements from demux input to GBS input.

### 1.2 純損失チャネル（prompt §7-3、conventions §2）

ビームスプリッタ模型: â_out = √η â_in + √(1−η) v̂、v̂ は真空。共分散行列では

  σ_out = η σ_in + (1 − η) V_vac 𝟙,   V_vac = 1/2.

任意の角度 θ の直交位相の分散は V_out(θ) = η V_in(θ) + (1−η)/2。ショット雑音で規格化すると（Ref A/B の書き方）

  R_out = η R_in + (1 − η).

**規約の変換**: Ref A の R± = L + (1−L)·exp(±2√(αP)) は、R_in = e^{±2r}（r = √(αP)、conventions §3 の Var = e^{∓2r}/2 を V_SN = 1/2 で割ったもの）に η = 1 − L の純損失をかけた式そのもの。Ref B の R± = 1 ± ηξ²ζρ·4x/((1∓x)² + 4Ω²) も、全効率 η_tot = ηξ²ζρ の純損失として同じ形（ここでの η は Ref B の PD の量子効率で、本ノートの η とは別の記号）。

**合成**: 純損失チャネルの直列は純損失チャネルで、η_tot = Π_i η_i（σ に順に作用させれば明らか）。純損失は位相回転と可換（σ → R(φ)σR(φ)ᵀ と 𝟙 が可換）なので、受動的な損失と位相シフトだけの区間では**順序は結果に効かない**。dB では IL_tot = Σ IL_i。

**分布損失**: 減衰係数 α の受動ファイバは、長さ L_f を M 個に切った η_i = e^{−αL_f/M} の直列なので、M に依らず η = e^{−αL_f}（α_dB [dB/km]、L_f [m] で書けば η = 10^{−α_dB L_f/10⁴}）。受動ファイバでは**分布と集中は厳密に同じ**。

分布が効くのは、損失と同時に別の過程があるとき:
- 損失と利得が同時（導波路 OPA の伝搬損失）: 実効損失は利得に依存する。Ref A は 0.1 dB/cm・45 mm・640 mW で実効 2% [Ref A §4.2]。T01 の担当。
- 損失と雑音の注入が同時（自発ラマン、GAWBS など、T07）: 出力の分散は

  V_out = η(0,L_f) V_in + [1 − η(0,L_f)] V_vac + ∫_0^{L_f} ν(z) η(z,L_f) dz,   η(z,L_f) = e^{−α(L_f − z)}.

  ν が一様なら注入雑音の合計は ν(1 − e^{−αL_f})/α ≈ νL_f(1 − αL_f/2)。我々の長さ（αL_f ≤ 0.07、§4.3）では集中モデルとの差は 4% 以下の補正 [AE]。
- 非線形（SPM/XPM）がある区間: squeezed vacuum 自体は弱いので損失との順序は問題にならないが、強い LO・probe のパルスでは損失で非線形位相が変わる（T06）。

### 1.3 モード不整合

BHD は LO のモード u_LO に射影した signal の直交位相を測る。signal 側のモード u_s とのずれを

  â_LO-mode = ⟨u_LO|u_s⟩ â_s + √(1 − η_mm) â_⊥,   η_mm = |⟨u_s|u_LO⟩|²

と書くと、â_⊥ は u_s に直交するモードの演算子。**光源が単一モード（u_s 以外は真空）なら** â_⊥ は真空で、モード不整合は η = η_mm の純損失と厳密に同じ。

可視度との関係: パワー P_1、P_2 の 2 本の干渉で V_vis = 2√(P_1P_2)|⟨u_1|u_2⟩|/(P_1 + P_2)、すなわち

  η_mm = V_vis² (P_1 + P_2)² / (4 P_1 P_2)  →（P_1 = P_2 のとき）η_mm = V_vis².

- Ref A: probe（squeezer OPA で増幅された光）と LO の「干渉効率」約 99% → 約 2% の損失 [Ref A §4.2]。V_vis = 0.99 → η_mm = 0.980 で整合。残り 2% は「squeezed 光と probe の空間モードのわずかな違い」と推定（原因は特定されていない）[Ref A §4.2]。
- Ref B: ξ = 0.988（homodyne visibility）を全効率に ξ² で入れている [Ref B Eq. (2)、§4] → η_mm = 0.976。

ガウスビームの重なり（同じ位置にウエストがある場合、導出は 1 次元の積分の積）:

  η_mm = [2 w_1 w_2/(w_1² + w_2²)]² · exp[−2d²/(w_1² + w_2²)] · exp[−(π n θ_t)² · 2w_1²w_2²/(λ²(w_1² + w_2²))]

同じ径（w_1 = w_2 = w）なら η_mm = exp(−d²/w²) · exp(−(π n w θ_t/λ)²)。

**ファイバ系での読み替え（将来の構成 F、Ref D）**: 単一モードファイバの中では squeezed 光も LO も同じ基本モードなので、ファイバ型の BHD（2×2 カプラ）では横モードの不整合はほぼ 0 になる。代わりに次が効く。
- 自由空間 → ファイバの結合効率（導波路 OPA の出力をファイバに入れるところ）。これは横モードの不整合だが、「結合損失」として squeezed 光にだけかかる純損失。
- コネクタ・スプライスでの軸ずれ（§3.2）。PANDA PM 1550 の MFD = 10.5 ± 0.5 µm [web: Corning PANDA PM] → w ≈ 5.25 µm。軸ずれ d = 0.5 µm で η = exp(−(0.5/5.25)²) = 0.991（0.04 dB）、d = 1 µm で 0.964（0.16 dB）[AE]。
- 偏光の不整合 η_pol（PM の PER、SMF の偏波ドリフト）: T09。
- 時間・スペクトルモードの不整合（パルス）: T03。

**適用範囲の注意**: u_s に直交するモードも squeezing されている場合（準単一モード導波路の高次モード [Ref A §4.2]、パルスの複数の時間モード）、â_⊥ は真空ではなく、不整合は純損失ではなく「別のモードの（反）squeezing の混入」になる。このとき余分な雑音は純損失より悪くなりうる。シミュレータでは、O 段のモード重なりの表（T03・T01）が多モードを扱うまでは純損失として扱い、その仮定をログに残す。

### 1.4 回路雑音の等価損失（T13 と共有）

白色でガウスの回路雑音（ショット雑音比 E_el）を差し引かずに測ると、

  R_meas = (R + E_el)/(1 + E_el) = η_el R + (1 − η_el),   η_el = 1/(1 + E_el)

で、**厳密に**純損失と同じ形。Ref A のクリアランス 28 dB → E_el = 1.58×10⁻³ → L_el = 0.158%（Ref A の「0.2% 未満」と整合）[Ref A §4.2]。損失の収支表ではこの項を必ず分けて書く（T13 が値を出す）。

### 1.5 位相ジッタとの合成

conventions §6 の厳密式で、

  R′_− = (R_− + R_+)/2 + (R_− − R_+)/2 · e^{−2θ²},   R_± = L_tot + (1 − L_tot) e^{±2r}.

r が大きいとき R′_− ≈ L_tot + (1−L_tot)(e^{−2r} + θ² e^{2r})。r について最小にすると e^{2r} = 1/θ（r_opt = ½ ln(1/θ)）で

  min_r R′_− ≈ L_tot + 2θ(1 − L_tot).

ここから全損失の上限（§4.1）が出る。

---

## 2. 仮定と適用範囲

- 損失はすべて線形・受動で、パワーに依らない（例外: 導波路 OPA の実効損失はポンプに依存する（T01）、Ref B の L(x) はポンプの関数（T02））。
- 損失は位相に依らない（偏光依存損失 PDL は偏光が固定された PM 系では 1 つの値に吸収できる。SMF では T09）。
- 単一モード近似（§1.3 の注意）。
- 損失は squeezing 帯域内で周波数に依らない。例外: (i) 検出器の帯域（T13）、(ii) 広帯域パルス（100 fs で FWHM 約 25 nm [AE]: sech² で Δν = 0.315/τ_p）に対する部品の波長依存（PBS・WDM・サーキュレータの帯域は ±20–40 nm 程度 [UNVERIFIED]）。(ii) のとき透過 t(ω) は「純損失」ではなく、モードの形を変える。単一の時間モード u(ω) に対する実効透過は、LO を t·u に合わせた場合で η = ∫|t(ω)|²|u(ω)|² dω。合わせない場合は T03 の重なりの表で扱う。
- エタロン（§5）は反射率が小さい（√(R_1R_2) ≪ 1）場合の 1 次近似。

---

## 3. 損失の値（リファレンスと部品）

### 3.1 リファレンス実験の損失収支

| 構成 | 項目 | 値 | 出典 |
|---|---|---|---|
| Ref A | 導波路の実効損失（0.1 dB/cm、45 mm、640 mW） | 2% | [Ref A §4.2] |
| Ref A | 光学部品（ミラー、レンズ） | 1%（測定値） | [Ref A §4.2] |
| Ref A | PD の光電変換（QE 約 99%） | 約 1% | [Ref A §3.1・§4.2] |
| Ref A | 回路雑音の等価損失（28 dB、3 MHz、LO 15.4 mW） | < 0.2%（§1.4 で 0.16%） | [Ref A §4.2] |
| Ref A | モード不整合（残差から推定） | 約 4%（うち 2% は probe–LO の可視度 0.99、2% は原因不明） | [Ref A §4.2] |
| Ref A | **合計（当てはめ）** | **8%**（積で 1 − 0.98·0.99·0.99·0.998·0.96 = 7.98%） | [Ref A §4.1] |
| Ref A（以前） | 合計 | 12%（θ = 14 mrad） | [Ref A §4.1] |
| Ref B | ρ = T/(T + L) = 0.123/(0.123 + 0.0038)、ζ = 0.99、ξ = 0.988、η = 0.998 | 1 − ρζξ²η = 0.0645（Eq. (2) の当てはめでは 0.0709 ± 0.0045） | [Ref B §4] |
| Ref C | ファイバ前（卓上の伝搬、モード整合、BHD の QE 97% を含む） | 15% | [Ref C §3] |
| Ref C | ファイバへの結合 ＋ 1 km の伝搬 ＋ 取り出し | 8%（0.4 dB） | [Ref C §3–4] |
| Ref C | うち 1 km × 0.2 dB/km の伝搬 | 4.5% [AE] → 結合と取り出しの合計 約 3.7%（1 − 0.92/0.955）[AE] | [Ref C §3] から計算 |
| Ref C | 当てはめの全損失 | 21%（積では 1 − 0.85·0.92 = 21.8%） | [Ref C §3] |
| Liu 2026 | 光源の効率 92%、16 モード干渉計 99%/モード、検出効率 93%、全体 51% | 残り（遅延ループ・結合など）は 0.51/(0.92·0.93·0.99³) ≈ 0.61 [AE] | [Liu 2026 本文・Methods] |

### 3.2 部品の挿入損失（将来の構成 F と Ref D）

| 部品 | 値 | 出典・確度 |
|---|---|---|
| FC/APC の嵌合 1 対（ランダム嵌合、SMF） | IEC 61753-1 の Grade B: 97% 以上の組み合わせで ≤ 0.25 dB。平均 ≤ 0.12 dB は [UNVERIFIED] | SENKO 白書（Grade B ≤ 0.25 dB ≥ 97%）[web]。規格本文は未読 [UNVERIFIED] |
| FC/APC の嵌合 1 対（PM パッチケーブル） | 最大 0.30 dB（Thorlabs P3-1550PMP-2 の検査票） | 検索結果の抜粋のみ [UNVERIFIED] |
| LN 変調器のコネクタ 1 個あたりの追加損失 | 最大 0.25 dB / FC/APC | exail MPX/MPZ データシート [web] |
| **survey E52 の「0.1–0.3 dB / 対」** | 上の 3 つと整合。本ノートでは**名目 0.2 dB**（範囲 0.12–0.3 dB）を使う | — |
| 融着接続（スプライス） | 典型値 SMF 0.03 dB、PMF 0.06 dB（偏光クロストーク典型 −40 dB、PANDA） | Fujikura FSM-100P データシート [web] |
| PM ファイバ PANDA PM 1550（PM15-U25D） | 減衰 最大 0.5 dB/km、MFD 10.5 ± 0.5 µm、クロストーク最大 −30 dB @ 100 m、ビート長 3.0–5.0 mm | Corning PANDA PM（M0500017、2019 年 1 月）[web] |
| SMF（Corning SMF-28 Ultra） | 減衰 ≤ 0.18 dB/km @ 1550 nm、曲げ損失 10 巻 × 半径 15 mm で ≤ 0.05 dB @ 1550 nm | Corning SMF-28 Ultra（PB115、2014 年 9 月）[web]。現行版 PI-1424 は未読 |
| ファイバ EOM（LN 位相変調器、exail MPX-LN-0.1、Ref A の EOM） | IL typ. 2.7 dB / max 3.5 dB（コネクタなし）＋ コネクタ 1 個あたり ≤ 0.25 dB。ORL typ. −45 dB / min −40 dB | exail MPX/MPZ データシート 11_2023_ED16 [web] |
| ファイバ EOM（exail MPZ-LN-01、Ref C の MPZ-LN-10 と同系列） | IL typ. 2.5 dB / max 3.5 dB、低損失オプション typ. 2 / max 2.5 dB | 同上 [web]（MPZ-LN-10 そのものの値は未確認）|
| PM サーキュレータ（Thorlabs CIR1550PM） | IL < 0.9 dB | GoPhotonics の掲載（検索結果の抜粋のみ）[UNVERIFIED] |
| ファイバ型 PBS/PBC（1550 nm） | typ. 0.4–0.5 dB、max 0.6–0.7 dB | 複数メーカーの掲載（検索結果の抜粋のみ）[UNVERIFIED] |
| PM WDM（1550 通過） | 0.03 dB（1 台の検査値）– ≤ 0.8 dB（仕様） | 検索結果の抜粋のみ [UNVERIFIED] |
| 2×2 PM カプラ PNH1550R5A2（Ref D） | 過剰損失: **資料を取得できず** | Thorlabs のサイトが JS で描画され取得できない。PDF の追加が必要 [UNVERIFIED] |
| ストレッチャ 915B（Ref D、D10 の候補） | 38 巻、PM15-U25D、FC/APC（slow 軸にキー）、全電圧範囲の伸び 16 µm（#8669-1）/ 17 µm（#8669-2）。**挿入損失の記載なし** | [manual] |
| 915B のファイバ長 | 約 9 m（「9 m のファイバで 16 µm」） | photonics.com の掲載（抜粋のみ）[UNVERIFIED] |
| 他社のストレッチャ（参考、General Photonics FST-001） | IL < 0.2 dB | 検索結果の抜粋のみ [UNVERIFIED] |
| 自由空間の PBS・ミラーなど（Liu 型の probe の入れ出しを自由空間で行う場合） | 1 素子あたり約 1%（0.044 dB）を仮置き | Ref A の「光学部品 1%」、Liu の干渉計 99%/モードからの類推 [AE] |
| 振り分け素子（T11 の仮の値） | (a) ファイバ型 LN 変調器: 上の EOM の値（2.7 dB）で代用、(b) 低損失（自由空間の Pockels セルなど）: 2%（0.09 dB）を仮置き | (a) [web] の代用、(b) [AE]。T11 の結果で置き換える |

dB と損失の対応（参考）: 0.03 dB = 0.69%、0.06 dB = 1.37%、0.12 dB = 2.73%、0.2 dB = 4.5%、0.3 dB = 6.7%、0.5 dB = 10.9%、0.9 dB = 18.7%、2.7 dB = 46.3%。

### 3.3 915B を squeezed 光の経路に置く場合（D10、D30）

915B の挿入損失は資料にない。中身は (i) 両端の FC/APC ピッグテール（嵌合 2 対、計 約 0.4 dB、または融着で 約 0.12 dB）、(ii) PZT 円筒への巻き付けによる曲げ損失（円筒の径は資料にない）、(iii) 約 9 m のファイバ（0.0045 dB、無視できる）。本ノートでは (ii) を **0.05 dB [AE, UNVERIFIED]** と仮置きし、接続は §4.2 の接続数 N_j に数える。巻き付けによる PER の劣化は T09、伸びとレンジは T16。

---

## 4. 我々のパラメータでの概算

### 4.1 目標 squeezing に対する全損失の上限（D16）

GBS 入力（または BHD）で S [dB] の squeezing を得るのに許される全損失 L_max。§1.5 の厳密式で、r を最適にとる場合と、Ref A の光源（r = √(αP) = 2.408、α = 906 %/W、P = 640 mW）に固定する場合。

**r を最適化（r_opt ≈ ½ ln(1/θ)）**

| θ [mrad] \ S [dB] | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | r_opt |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 50.1% | 39.8% | 31.6% | 25.1% | 20.0% | 15.8% | 12.6% | 10.0% | ∞ |
| 5 | 49.6% | 39.2% | 30.9% | 24.4% | 19.1% | 15.0% | 11.7% | 9.1% | 2.65 |
| 9 | 49.2% | 38.7% | 30.4% | 23.7% | 18.5% | 14.3% | 11.0% | **8.4%** | 2.36 |
| 14 | 48.7% | 38.1% | 29.7% | 23.0% | 17.6% | 13.4% | 10.1% | 7.4% | 2.13 |
| 20 | 48.0% | 37.3% | 28.8% | 22.0% | 16.6% | 12.3% | 8.9% | 6.3% | 1.96 |
| 30 | 46.9% | 36.0% | 27.3% | 20.3% | 14.8% | 10.5% | 7.0% | 4.3% | 1.75 |
| 50 | 44.6% | 33.1% | 24.0% | 16.8% | 11.1% | 6.5% | 2.9% | 不可 | 1.50 |

**r = 2.408（Ref A の光源）に固定**

| θ [mrad] \ S [dB] | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|
| 0 | 49.7% | 39.3% | 31.1% | 24.5% | 19.3% | 15.2% | 11.9% | 9.3% |
| 9 | 49.2% | 38.7% | 30.4% | 23.7% | 18.5% | 14.3% | 11.0% | 8.3% |
| 14 | 48.5% | 37.8% | 29.3% | 22.6% | 17.3% | 13.0% | 9.7% | 7.0% |
| 20 | 47.1% | 36.1% | 27.5% | 20.6% | 15.1% | 10.7% | 7.3% | 4.5% |
| 30 | 43.4% | 31.7% | 22.4% | 15.0% | 9.1% | 4.5% | 0.8% | 不可 |
| 50 | 27.1% | 12.0% | 0% | 不可 | 不可 | 不可 | 不可 | 不可 |

近似式（r 最適化）: L_max ≈ (10^(−S/10) − 2θ)/(1 − 2θ)。10 dB・9 mrad で 0.0835（表 0.084）。θ = 0 では L_max = 10^(−S/10)（3 dB で 50.1%、10 dB で 10.0%）。

θ_rms の帯域の定義（ESA と GBS で違う）は T04 が決める。ここでの θ は §1.5 の式に入る値。

### 4.2 32 ポート構成の損失収支（ポート k の関数）

**モデル**（D08、D10、D30。値はすべて仮で、§3.2 の出典つき）:

  η_k = η_up · η_dmx(k) · 10^{−[N_j IL_j + N_inj IL_inj + IL_str + α_dB L_k/1000]/10},
  η_dmx(k) = 10^{−N_sw(k) IL_sw/10}.

- η_up = 0.98（導波路の実効損失、Ref A の値で代用、T01）× 0.99（光学部品、Ref A）× 0.98（自由空間 → ファイバの結合、Ref C から [AE]）= **0.951**。
- 振り分け素子の構造（T11）: 2 分木（1×2 を 5 段）なら N_sw = 5（k に依らない）、バスから 1 個ずつ取り出す鎖なら N_sw(k) = k + 1。下の表は「1 段（N_sw = 1）」とし、構造による倍数は T11 の結果を待つ。
- 1 ポートの接続数 N_j = 6 [AE]: 上流ファイバ → 振り分け素子、振り分け素子 → probe の取り出し素子、→ ストレッチャ、→ 遅延ファイバ、→ probe の入射素子、→ GBS 入力。
- N_inj = 2（Liu 型: 上流の取り出し素子と、GBS 側の入射素子。D30）。
- L_k = k c/(f_rep n_g)、n_g = 1.468 [AE, UNVERIFIED、T08 で確認]。f_rep は Open（D07）なので 100 MHz と 10 MHz を並べる。

**3 つの想定**

| 記号 | IL_sw | 接続（IL_j） | probe の入れ出し（IL_inj） | ストレッチャ IL_str |
|---|---|---|---|---|
| NF（すべてファイバ部品、コネクタ） | 2.7 dB（ファイバ LN） | FC/APC 0.2 dB × 6 | ファイバ PBS 0.5 dB × 2 | 0.05 dB |
| MF（融着、ファイバ PBS、低損失の振り分け） | 0.09 dB | 融着 0.06 dB × 6 | ファイバ PBS 0.5 dB × 2 | 0.05 dB |
| LL（融着、自由空間の入れ出し、低損失の振り分け） | 0.09 dB | 融着 0.06 dB × 6 | 自由空間 0.044 dB × 2 | 0.05 dB |

**結果（GBS 入力までの全損失 L = 1 − η_k、S_max = 10 log10(L)（r → ∞、θ = 0）、S(A, 9 mrad) は Ref A の光源 r = 2.408・θ = 9 mrad）**

f_rep = 100 MHz、PM（0.5 dB/km）:

| 想定 | k | L_k [m] | ポートの IL [dB] | L | S_max [dB] | S(A, 9 mrad) [dB] |
|---|---|---|---|---|---|---|
| NF | 0 | 0 | 4.95 | 69.6% | −1.57 | −1.54 |
| NF | 31 | 63.3 | 4.98 | 69.8% | −1.56 | −1.53 |
| MF | 0 | 0 | 1.50 | 32.7% | −4.86 | −4.70 |
| MF | 31 | 63.3 | 1.53 | 33.2% | −4.79 | −4.64 |
| LL | 0 | 0 | 0.59 | 17.0% | −7.71 | −7.34 |
| LL | 1 | 2.0 | 0.59 | 17.0% | −7.70 | −7.33 |
| LL | 15 | 30.6 | 0.60 | 17.3% | −7.63 | −7.27 |
| LL | 31 | 63.3 | 0.62 | 17.6% | −7.55 | −7.20 |

f_rep = 10 MHz（L_31 = 633 m）:

| 想定 | ファイバ | k = 0 の L | k = 15 の L | k = 31 の L | k = 31 の S(A, 9 mrad) [dB] |
|---|---|---|---|---|---|
| LL | PM | 17.0% | 19.8% | 22.8% | −6.16 |
| LL | SMF | 17.0% | 18.0% | 19.1% | −6.87 |
| MF | PM | 32.7% | 35.0% | 37.4% | −4.14 |
| MF | SMF | 32.7% | 33.5% | 34.4% | −4.48 |

**ポートごとに許される損失**（η_up = 0.951、θ = 9 mrad、r 最適化、§4.1 の L_max から η_port ≥ (1 − L_max)/η_up）:

| S [dB] | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|
| GBS 入力でのポートの IL 上限 [dB] | 2.72 | 1.91 | 1.36 | 0.96 | 0.67 | 0.45 | 0.29 | **0.16** |
| 同じ値を BHD で測る場合（LO の整合 0.99、QE 0.99、回路雑音 0.998 を追加） [dB] | 2.63 | 1.81 | 1.26 | 0.86 | 0.57 | 0.36 | 0.19 | **0.07** |

読み取り:
- 10 dB は、上流を Ref A 並み（4.9%）にしても、ポートの部品に 0.16 dB しか残らない。FC/APC 1 対で使い切る。融着 0.06 dB × 6 = 0.36 dB でも超える。**10 dB を GBS 入力で出すには、ポートの部品を自由空間の低損失光学系にまとめるか、接続数を 1–2 個に減らす必要がある** [AE]。
- 6 dB なら 0.96 dB まで許され、MF でも 1.5 dB で超えるが、LL（0.6 dB）なら入る。
- ファイバ型 LN 変調器（2.7 dB）1 段だけで、上流と合わせて L = 1 − 0.951 × 0.537 = 49%、squeezing の上限は 3.1 dB。想定 NF 全体では上限 約 1.6 dB で、3 dB の目標にも届かない。
- ファイバ型 PBS（0.5 dB × 2）は Liu 型の probe の入れ出しの代償として大きい。§6.3 も参照。

### 4.3 遅延ファイバの減衰（PM と SMF、D32）

| f_rep | L_31 | PM（0.5 dB/km） | SMF（0.18 dB/km） |
|---|---|---|---|
| 1 GHz | 6.3 m | 0.003 dB（0.07%） | 0.001 dB（0.03%） |
| 100 MHz | 63.3 m | 0.032 dB（0.73%） | 0.011 dB（0.26%） |
| 20 MHz | 316.5 m | 0.158 dB（3.6%） | 0.057 dB（1.3%） |
| 10 MHz | 633 m | 0.317 dB（7.0%） | 0.114 dB（2.6%） |
| 1 MHz | 6.33 km | 3.17 dB（52%） | 1.14 dB（23%） |

PM の 0.5 dB/km は**最大値**（Corning の仕様）で、実際の値はこれより小さい可能性がある。ポート間の損失の差（k = 0 と k = 31）は GBS の入力ごとに squeezing が揃わない原因になる（GBS_REQUIREMENTS は Open）。巻いた遅延ファイバの曲げ損失は、SMF-28 Ultra では半径 15 mm × 10 巻で ≤ 0.05 dB [web] なので、半径の大きいスプールでは小さいと見込む。PM15-U25D の曲げ損失の仕様は Corning の資料にない [UNVERIFIED]。

### 4.4 Ref D（現行系）での損失の役割

Ref D の光は古典的なコヒーレント光で、損失は squeezing ではなく**誤差信号の傾き K_D**（T15）と縞の可視度に効く。
- 透過 η_opt は K_D に比例して入る。接続 約 6 対 × 0.2 dB ＋ カプラの過剰損失 2 個（値なし、0.1–0.3 dB と仮置き [UNVERIFIED]）で η_opt ≈ 0.6–0.75 [AE]。PNH1550R5A2 の資料がないので、D02 の「メーカーの資料の値だけ」を満たせない。
- 腕の損失の差は可視度を下げる: 50:50 の結合で V_vis = 2√(T_a T_b)/(T_a + T_b)。915B の腕だけ 0.3 dB 多くても V_vis = 0.99988 で無視できる [AE]。カプラの分岐比のずれと偏光（PER）のほうが効く（T15・T09）。
- ブラインド検証（D04）に従い、実測の損失・可視度は使っていない。

---

## 5. 後方反射と寄生エタロン（E28）

2 つの反射点（振幅反射率 r_1、r_2、間の片道の位相 ψ_e）がある区間の透過は

  t = t_1 t_2 e^{iψ_e} / (1 − r_1 r_2 e^{2iψ_e}).

√(R_1R_2) ≪ 1 で、位相のリップル δφ ≈ √(R_1R_2) sin(2ψ_e + φ_r) [rad]、パワーのリップル δP/P ≈ 2√(R_1R_2) cos(2ψ_e + φ_r)（φ_r は反射の位相）。

- FC/APC 同士（反射 −60 dB 程度と仮定 [UNVERIFIED]）: δφ ≈ 10⁻⁶ rad。
- LN 変調器（ORL −40 dB min、−45 dB typ. [web]）と FC/APC: δφ ≈ √(10⁻⁴ × 10⁻⁶) = 10⁻⁵ rad。
- 温度で ψ_e がドリフトする（survey E22 の 約 44 rad/(m·K) [AE]）と、このリップルは遅い位相雑音になる。振幅は 10⁻⁵ rad 以下で、θ = 9 mrad に比べて無視できる。
- パルス系（F）では、エタロンの往復 2n_gL_e/c が T_rep の整数倍になると、別スロットのパルスの複製が重なる。パワーで R_1R_2 ≈ 10⁻¹⁰–10⁻¹² で無視できる [AE]。
- Rayleigh 後方散乱を分布反射として扱うと、63 m で R ≈ 3×10⁻⁶ [survey §3、AE]。2 回散乱で前向きに戻る成分は R² のオーダーで無視できる。逆向き probe の 1 回散乱が squeezed 光に混ざる効果（E44）は T10。

結論: 寄生エタロンは squeezing にも位相ロックにも効かない。T ブロックとしては「既定でオフ」の雑音源で十分。ただし、BHD の PD 面で反射した強い LO が 2×2 カプラを逆にたどって squeezed 光の経路を戻り、どこかで再反射して戻ってくる経路は、LO の光子を squeezed モードに入れる。PD の反射率は資料にないので T13・T10 で確認する [UNVERIFIED]。

---

## 6. conventional tap のトレードオフ（E56、prompt §7-8）

### 6.1 モデル

conventional 方式（Ref A Fig. 1(a)）では、squeezer OPA の後で比率 r_t を tap し、同軸の probe（出力パワー P_pr）から誤差信号を作る。

- squeezed 光への損失: L_tot(r_t) = 1 − (1 − L_0)(1 − r_t)。
- センサ起因の位相雑音（片側 PSD、tap の PD で決まる）:

  S_φ,sens = κ_sn h ν_opt / (η_d P_tap) + NEP² / (m P_tap)²,   P_tap = r_t P_pr.

  第 1 項はショット雑音。理想的なホモダイン型の位相の読み出しで κ_sn = 1/2（conventions §2 で Var(φ) = 1/(4N) から）、ヘテロダイン・復調では大きくなる。第 2 項は PD の電気雑音、m = (dP_det/dφ)/P_tap は誤差信号の傾きの係数。κ_sn と m の値は T14 が誤差信号の方式ごとに出す。ここでは κ_sn = 1/2、m = 1 を仮置き。
- 閉ループでセンサ雑音は帯域内で位相に移る。1 次のループ（単位利得 f_BW）で θ_sens² = (π/2) f_BW S_φ,sens。
- θ² = θ_0² + θ_sens²、R′_− は §1.5 の式。

### 6.2 概算（Ref A 相当: L_0 = 8%、r = 2.408、θ_0 = 9 mrad、f_BW = 10 kHz、η_d = 0.9。P_pr と NEP は仮の値）

| P_pr | NEP [pW/√Hz] | 最適な r_t | そのときの S [dB] | θ [mrad] |
|---|---|---|---|---|
| 1 µW | 0（ショット雑音のみ） | 3×10⁻⁴ | −10.12 | 9.2 |
| 1 µW | 1 | 1×10⁻² | −9.09 | 15.4 |
| 1 µW | 10 | 1×10⁻¹ | −6.92 | 15.4 |
| 10 µW | 1 | 3×10⁻³ | −9.94 | 9.9 |
| 10 µW | 10 | 1×10⁻² | −9.09 | 15.4 |
| 100 µW | 10 | 3×10⁻³ | −9.94 | 9.9 |
| tap なし（phase-detection OPA） | — | — | −10.15 | 9.0 |

読み取り: ショット雑音だけなら µW の probe でも tap は 10⁻³ 以下で足り、損失はほぼ問題にならない。トレードオフを作るのは tap 用の PD の**電気雑音**で、Ref A が conventional 方式の probe を「数 µW」に抑えている（probe が BHD に入るのを避けるため [Ref A §2]）こととあわせると、最適点で 0.2–3 dB 悪くなる。phase-detection OPA では tap の損失がなく、probe のパワーを mW まで上げて SNR を稼げる [Ref A §2]。

### 6.3 我々の構成（Liu 型、D30）との関係

Liu 型では squeezed 光を tap しないが、probe を入れ出しする素子（PBS、サーキュレータ、WDM）が squeezed 光の経路に 2 個入る（§4.2 の N_inj = 2）。その損失は probe のパワーに依らない固定の値で、SNR は probe のパワーで上げられる（上限は GBS への混入、T10）。したがって「tap 比 ↔ SNR」のトレードオフは「素子の挿入損失 ↔ probe の混入（PER、Rayleigh）」に置き換わる。直交偏光の probe なら PBS（ファイバ型 0.4–0.7 dB [UNVERIFIED]、自由空間なら約 1% [AE]）、同じ偏光ならサーキュレータ（< 0.9 dB [UNVERIFIED]）か tap カプラ（r_t の損失）。D30 の「Liu 型と phase-detection OPA の統合」でも、ポートごとの位相を測る限り、この入れ出しの素子は残る [AE]。

---

## 7. 周波数帯とスケーリング

| 効果 | 周波数帯 | ファイバ長 | パワー | 温度 | パルス幅 |
|---|---|---|---|---|---|
| 部品の挿入損失（E52、E55、E45） | 静的（DC）、squeezing 帯域内で一定 | 依らない（接続数に比例） | 依らない | コネクタの温度変化は未調査 | 100 fs 級では部品の帯域の端で波長依存（§2） |
| ファイバの減衰 | 静的 | ∝ L_k ∝ k/f_rep | 依らない | 依存は小さいと仮定（資料なし） | 同上 |
| 空間モードの不整合（E53） | 静的 ＋ 遅いドリフト（アライメント、偏光） | ファイバ内では 0（結合点だけ） | 依らない（熱レンズは別） | 遅いドリフト | 時間モードは T03 |
| 回路雑音の等価損失 | 測定周波数で変わる | 依らない | ∝ 1/P_LO | — | パルス BHD は T13 |
| 寄生エタロン（E28） | DC–Hz（温度で ψ_e がドリフト） | L_e に比例して温度感度が上がる | 依らない | ∝ ΔT | — |
| tap のトレードオフ（E56） | 誤差信号の雑音は f_BW 以下 | 依らない | SNR ∝ r_t P_pr | — | — |

---

## 8. シミュレータが必要とするもの

### 8.1 ブロック・表

| 段 | ブロック / 表 | 入力 | 出力 | オン・オフ |
|---|---|---|---|---|
| O | `loss_budget` 表: 光路の区間（D11、conventions §9）ごとの要素の列 | 要素ごとに id、IL_dB または η、区間 id、種類（connector / splice / fiber / component / mode_overlap / electronic）、出典の印、（任意）波長依存 t(λ) | 光路ごとの η_tot と内訳（squeezed → GBS_k、squeezed → BHD_k、LO_k → BHD_k、probe の経路）、L_k の関数 | 要素ごとの `enabled` |
| O | `mode_overlap`: η_mm の計算 | (w_1, w_2, d, θ_t, λ, n) または (V_vis, P_1, P_2) | η_mm | `MODE_MISMATCH_ENABLE` |
| Q / O | `pure_loss` チャネル（共分散行列への作用 σ → ησ + (1−η)/2·𝟙） | η | σ_out | `LOSS_ENABLE` |
| Q | 分布損失 ＋ 雑音注入（T07 と共有）: 区間を分割して損失と注入を交互に | α、ν(z)、L_f、分割数 | σ_out | T07 の雑音項の toggle |
| T | `tap`（conventional トポロジー）: squeezed 光に r_t の損失、P_tap を誤差信号生成（T14）へ | r_t、P_pr | squeezed 光の η、P_tap | `TAP_ENABLE`、`tap_ratio` |
| T | `etalon_noise`（既定でオフ） | R_1、R_2、L_e、温度の時系列 | 位相・強度のリップル | `ETALON_ENABLE`（既定 off） |
| T | （任意）η(t) のドリフト: 偏光（T09）やアライメントによる遅い変化 | η の時系列 | R±(t) | `LOSS_DRIFT_ENABLE`（既定 off） |
| M | 回路雑音の等価損失（T13 と共有。差し引く・差し引かないの両方を出す） | C_dB(f)、P_LO | η_el(f) | `CIRCUIT_NOISE_ENABLE` |

### 8.2 実装上の要件

- squeezed 光の経路に置く変調器・アクチュエータ・probe の入れ出し素子は、必ず `loss_budget` に要素として登録し、その損失を squeezed 状態にかける（prompt §5.5、D10）。
- 損失の値には出典の印（[Ref A] / [manual] / [web] / [AE] / [UNVERIFIED]）を持たせ、実行ログに保存する（設定とともに記録、prompt §6）。
- GBS 入力の squeezing（D16 の目標）と BHD で測る squeezing を別々に出す。BHD の側だけにある損失（LO との整合、QE、回路雑音）を GBS 入力の値に混ぜない。
- 振り分け素子の N_sw(k) と IL_sw は T11 から受け取る。L_k、n_g は T08/T11、η_up の導波路の項は T01、η_pol は T09、η_mm（時間モード）は T03 から受け取る。

---

## 9. 検証テスト

| # | テスト | 合格の基準 |
|---|---|---|
| V1（§7-3） | 純損失: V_in ∈ {e^{−2r}/2, e^{2r}/2, 熱的 (2n̄+1)/2} × η ∈ {0, 0.3, 0.92, 1} で V_out = ηV_in + (1−η)/2 | 相対誤差 < 10⁻¹² |
| V2 | 合成: η_i をランダムに 10 個並べた直列 = 単一の Π η_i。位相回転を間に挟んでも同じ | 相対誤差 < 10⁻¹² |
| V3 | 分布損失: α_dB = 0.5 dB/km、L_f = 633 m を M = 1, 10, 1000 に分割 → η = 10^(−0.3165/10) = 0.9297、M に依らない | 相対誤差 < 10⁻¹² |
| V4（§7-1 と共有） | Ref A の損失収支: (0.98, 0.99, 0.99, 0.998, 0.96) → L = 7.98%。L = 8%、α = 906 %/W、P = 640 mW、θ = 9 mrad で R′_− = −10.15 dB、R_+ = +20.56 dB | ±0.01 dB |
| V5 | Ref B: ρ = 0.123/0.1268、ζ = 0.99、ξ = 0.988、η = 0.998 → 1 − ρζξ²η = 0.0645 | ±0.0001 |
| V6 | Ref C: 1 − 0.85 × 0.92 = 0.218（当てはめ 21%）。1 km × 0.2 dB/km → 4.50% | ±0.001 |
| V7 | モードの重なり: w_2/w_1 = 1.1 → η_mm = 0.99097。同じ径・d/w = 0.0952（d = 0.5 µm、w = 5.25 µm）→ 0.99097。同じ径・θ_t = 8.7 mrad、w = 5.25 µm、n = 1、λ = 1.55 µm → exp(−(πwθ_t/λ)²) = 0.9915。V_vis = 0.99（同じパワー）→ 0.9801。P_1/P_2 = 4、|⟨u_1\|u_2⟩| = 1 → V_vis = 0.8 | 相対誤差 < 10⁻⁶ |
| V8 | 回路雑音: C = 28 dB で差し引かない測定値が η_el = 1/(1 + 10^(−2.8)) の純損失と一致、L_el = 0.158% | 相対誤差 < 10⁻⁹ |
| V9 | 全損失の上限: θ = 0、r → ∞ で L_max = 10^(−S/10)（3 dB → 0.501、10 dB → 0.100）。θ = 9 mrad・r 最適化で 10 dB → 0.084、6 dB → 0.237 | ±0.001 |
| V10（§7-8） | tap: §6.2 の条件（L_0 = 8%、θ_0 = 9 mrad、f_BW = 10 kHz、η_d = 0.9、κ_sn = 1/2、m = 1、P_pr = 10 µW、NEP = 1 pW/√Hz）で R′_−(r_t) が内部に最小を持ち、r_t ≈ 3×10⁻³ で −9.94 dB。phase-detection OPA のトポロジーでは r_t を変えても R′_− = −10.15 dB で変わらない | 最小の位置が r_t の格子 1 点以内、値 ±0.02 dB |
| V11 | 32 ポートの収支: 想定 LL・f_rep = 100 MHz・PM で k = 0 → L = 17.0%、k = 31 → 17.6%。N_sw(k) = k + 1（鎖）なら L は k について単調増加、N_sw = 5（木）なら k の差はファイバだけ | ±0.1% |
| V12 | エタロン: R_1 = R_2 = 10⁻⁴ で、厳密な透過の位相のリップルの振幅が √(R_1R_2) = 10⁻⁴ rad と一致 | 相対誤差 < 2% |
| V13 | GBS 入力と BHD の分離: BHD 側だけの損失（LO 整合、QE、回路雑音）をオフにすると BHD の値が GBS 入力の値に一致 | 完全一致 |

---

## 10. 文献

一次資料（このノートで本文を読んだもの）:
- [Ref A] K. Hirota, T. Kashiwazaki, G. Ha, T. Yamashima, P. Jaturaphagorn, T. Suzuki, K. Takahashi, A. Kawasaki, A. Inoue, W. Asavanant, M. Endo, T. Umeki, A. Furusawa, "Generation of 10-dB squeezed light from a broadband waveguide optical parametric amplifier with improved phase locking method." **掲載誌・年・DOI は PDF から読み取れない [UNVERIFIED]**。§2、§4.1、§4.2 を使用。
- [Ref B] Y. Takeno, M. Yukawa, H. Yonezawa, A. Furusawa, "Observation of −9 dB quadrature squeezing with improvement of phase stability in homodyne measurement," Opt. Express 15, 4321–4327 (2007)。DOI は PDF に記載なし [UNVERIFIED]。Eq. (2)、§4 を使用。
- [Ref C] S. Verclas, C. Darsow-Fromm, B. Tohermes, R. Schnabel, "Fiber distribution of phase-stabilized GHz-bandwidth squeezed vacuum states of light between two buildings," Optica Quantum 4, 1 (2026)。DOI 10.1364/OPTICAQ.567418。§3–4 を使用。
- [Liu 2026] H.-L. Liu et al., "Gaussian boson sampling with 1,024 squeezed states in 8,176 modes," Nature 653, 2026。DOI 10.1038/s41586-026-10523-6。本文と Methods（光源の効率、干渉計の透過、逆向き probe の PBS）を使用。

メーカー資料（`refs/manuals/`）:
- [manual] Evanescent Optics / ハナムラオプティクス、"Model 914-2 Piezo Controller"・"915B Piezo Assembly" 仕様書（`refs/manuals/914-915B-Specifications.pdf`、3 ページのスキャン）。巻数、ファイバの種類、伸び、コネクタを使用。挿入損失の記載はない。

メーカー資料（ウェブ、閲覧日 2026-10-28、D24）。開いて読んだもの [web]:
- Corning, "PANDA PM Specialty Fibers"（M0500017、Issued January 2019）。AMS Technologies が掲載: https://shop.amstechnologies.com/media/b9/6b/7b/1720716991/PANDA-PM-Polarization-Maintaining-Optical-Fibers-Corning-Datasheets.pdf （corning.com 上の版は 403 で取得できず）
- Corning, "SMF-28 Ultra Optical Fiber"（PB115、September 2014）。第三者の掲載: https://www.icsgroup.ru/upload/medialibrary/8a2/SMF-28_UltraFiberBrochure_10-2014.pdf （現行の PI-1424 は未読。検索結果では同じ 0.18 dB/km）
- exail (iXblue), "MPX and MPZ LN Series" modulator datasheet（11_2023_ED16）。Laser Components が掲載: https://lasercomponents.com/fileadmin/user_upload/home/Datasheets/ixblue/mpx-mpz-series-modulator.pdf
- Fujikura, "FSM-100P data sheet": https://www.europe.fujikura.com/wp-content/uploads/2025/08/FSM-100P-data-sheet.pdf
- SENKO, "MPO Connector Random Mating IL versus IL by Master Jumper"（白書）: https://www.senko.com/wp-content/uploads/2024/02/White-Paper_MPO-Connector-Random-Mating-IL-versus-IL-by-Master-Jumper.pdf （IEC 61753-1 の Grade B を引用。規格本文は未読）

検索結果の抜粋だけで見たもの（**UNVERIFIED**、原典を開けていない）:
- Thorlabs P3-1550PMP-2 検査票（IL 最大 0.30 dB）: https://www.thorlabs.de/images/TabImages/P3-1550PMP-2_SampleSheet.pdf （thorlabs.co.jp へ転送され、JS で描画されるため取得できない）
- Thorlabs CIR1550PM（IL < 0.9 dB）: https://www.gophotonics.com/products/fiber-optic-circulators/thorlabs-inc/15-46-cir1550pm-fc （403）
- ファイバ型 PBS/PBC 1550 nm（typ. 0.4–0.5 dB）: Newport https://www.newport.com/f/polarization-beam-combiner-and-splitter ほか
- Evanescent Optics 915B（9 m で 16 µm）: https://www.photonics.com/Buyers_Guide/Products/Model_915B_Piezo-electric_Fiber_Stretcher/psp9424
- General Photonics FST-001（IL < 0.2 dB）: https://lunainc.com/fst-001-data-sheet
- Thorlabs PNH1550R5A2: 製品ページ https://www.thorlabs.com/thorproduct.cfm?partnumber=PNH1550R5A2 は取得できず、値なし。

式について: 純損失チャネル、モード重なり、可視度と効率の関係、エタロンの透過、閉ループのセンサ雑音の積分は、本ノートの中で導出した標準的な結果で、特定の文献の式は引用していない。教科書（H.-A. Bachor & T. C. Ralph, "A Guide to Experiments in Quantum Optics," 2nd ed., Wiley-VCH, 2004）に同じ内容があるはずだが、本文で確認していない [UNVERIFIED]。

---

## 11. 未解決の点（調整役・ユーザへ）

1. **D16 の目標の基準点**: TARGET_SQZ_DB は prompt §2 で「GBS 入力」。将来の構成で BHD（LO_k）を各ポートの監視に使うなら、BHD 側の損失（約 2.2%）の分だけ測定値は小さくなる。検証・最適化でどちらを目的関数にするか。
2. **squeezed 光の経路の接続方法**: D08 は「コネクタは FC/APC」だが、§4.2 のとおり FC/APC 1 対（約 0.2 dB）で 10 dB の目標の残りの余裕（0.16 dB）を使い切る。squeezed 光の経路は融着か自由空間にする前提で比べてよいか。
3. **probe の入れ出しを自由空間で行うか**（D30）: ファイバ型 PBS 2 個（約 1 dB）は 6 dB 以上の目標と両立しない。Liu 2026 は自由空間の PBS で入れている。T10 の比較の前提として、自由空間の選択肢を残すか。
4. 資料の不足（D02、D24）: PNH1550R5A2 と 915B の挿入損失、Thorlabs の PM パッチケーブル・サーキュレータの仕様の PDF を `refs/manuals/` に置いてほしい（Thorlabs のサイトは JS で描画され、エージェントから取得できない）。
5. 振り分け素子の構造（木か鎖か）と種類（T11、D08 Open）で、ポートごとの損失が k に依るかどうかが決まる。

## 12. 調整役への提案（直接は編集しない）

- `notes/survey.md` E52: 「FC/APC 1 対 0.1–0.3 dB（UNVERIFIED）」→ Grade B ≤ 0.25 dB（97%）、LN 変調器のコネクタ ≤ 0.25 dB で裏づけ（§3.2）。
- `notes/survey.md` E55: 「ファイバ EOM 2–3 dB（UNVERIFIED）」→ exail MPX-LN-0.1 で typ. 2.7 / max 3.5 dB ＋ コネクタ 1 個 ≤ 0.25 dB [web]。「ストレッチャは小」→ 915B の資料に値なし、接続 2 個分（0.12–0.4 dB）＋ 曲げ損失（未知）。
- `notes/survey.md` E45: 「サーキュレータ 0.5–1 dB」→ Thorlabs CIR1550PM < 0.9 dB（UNVERIFIED のまま）。
- `DECISIONS.md` Open に追加の候補: §11 の 1–3。
