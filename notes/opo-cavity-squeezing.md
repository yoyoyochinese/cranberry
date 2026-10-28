# T02 OPO キャビティの squeezing（opo-cavity-squeezing）

- Phase: 1（deep dive）／ 日付: 2026-10-28
- 担当する効果（`notes/survey.md`）: E05・E06（OPO 側）、E08、E09、E10、prompt §7-9 の Ref B Eq. (2) の確認
- 入力: `docs/squeezing_sim_prompt.md`、`DECISIONS.md`（D01–D33）、`notes/conventions.md`、`notes/survey.md`、
  Ref B（Takeno 2007）、Ref C（Verclas 2026）、Liu 2026（`refs/`）、Tohermes 2024（arXiv）
- 値の印はサーベイと同じ: [Ref B §n] など＝本文から読んだ値、[AE]＝この担当の概算（unverified）、[UNVERIFIED]＝原典で確認していない値・書誌。
- 数値の検算に使ったスクリプトは scratchpad にだけ置いた（リポジトリには入れていない）。式と結果はこのノートに全部書く。

## 0. 結論（先に）

1. **x の定義の食い違いはない。** Ref B は本文 §2、§4、Eq. (5)、Fig. 3 のキャプションのすべてで **x = √(P/P_th)**。
   サーベイ E10 の「本文は x = P/P_th」は `pdftotext` が根号の記号を落とした抽出の誤り（根号は前の行に浮いた「p」として出る）。
   PDF を画像にして確認した。Fig. 3 の測定点の横軸も、10, 20, 50, 65, 90, 100, 120, 130, 150 mW を √(P/180 mW) にした 0.236 … 0.913 と一致する。
2. **Ω の 2π は必要。** 正しくは Ω_B = 2π f / γ_B（γ_B = c(T + L)/l は**エネルギーの減衰率＝ FWHM [rad/s]**）。
   論文の「Ω = f/γ」をそのまま使うと、x → 1 で反 squeezing が +37 dB に発散し、Fig. 3 の実線（x = 1 で約 +21.5 dB、計算値 +21.4 dB）と合わない。
   2π を入れると、Fig. 3 の実線・点線の形と最小値（−9.35 dB、−7.17 dB）が再現できる。
3. **G = 18.7 は論文の他の数値と合わない（未解決）。** Eq. (4)–(5) で P = 100 mW、P_th = 180 mW なら x = 0.745、G = 15.4。
   G = 18.7 なら x = 0.769（P = 106 mW、または P_th = 169 mW）。論文の生成 squeezing の推定値 R± = −10.12/+15.15 dB を Eq. (2)（2π あり）で逆に解くと
   x = 0.753（P = 102 mW、G = 16.4）、全損失 0.0708（論文は 0.0709）。つまり P/P_th のほうが論文の squeezing データと合い、G = 18.7 は外れている。
   **検証目標では x = √(P/P_th)、P_th = 180 mW を使い、G = 18.7 は使わない。**
4. **§7-9 の検証目標（訂正版、§6 の V-B3）**: f = 1 MHz、L(x) = 0.00249 + 0.00222x、η = 0.998、ξ = 0.988、ζ = 0.99、T = 0.123、l = 0.5 m で、
   θ = 1.5° のとき x_opt = 0.758（P = 104 mW）、R′−,min = −9.35 dB。θ = 3.9° のとき x_opt = 0.607（P = 66 mW）、R′−,min = −7.17 dB。
   論文の「x ≈ 0.82 より上で劣化」は測定点についての記述で、モデルの最小は 0.758。最小の付近は平ら（0.1 dB 以内が x = 0.71–0.81）なので、x_opt だけでなく曲線全体で比べる。
5. **離調はほぼ純粋に squeezing 角の回転として効く。** 小さい離調 Δ（Ω ≪ κ）で δθ_sq ≈ −(Δ/κ)/(1 + x²)。
   Ref B では往復の光路長 1 nm で Δ/κ = 0.115、δθ_sq ≈ −4.2°/nm。**同軸の probe はこの回転を完全には追わない**:
   probe の位相は −Δ/(2κ) だけ回るので、ロック後に残る誤差は (Δ/κ)(1 − x²)/(2(1 + x²))（x = 0.745 で 0.143 Δ/κ ≈ 0.94°/nm）。
   prompt §5.5 の「同軸 = 区間の共有なので共通の雑音は自動的に打ち消し合う」は、OPO キャビティの区間では成り立たない（§8、提案 P3）。
6. Ref C の線幅 1.75 GHz は **FWHM**（Tohermes 2024 が「full-width-half-maximum」と明記。Ref C Fig. 4 の周波数依存も FWHM と仮定したときだけ合う）。
7. Ref D には OPO がないので、T02 は Ref D には関係しない。将来の構成 F（D07）の基本は単一パスの導波路 OPA で、D07 のパルス幅 100 fs–10 ps
   （帯域 4.4 THz–44 GHz）は CW 型 OPO の線幅（Ref B 12 MHz、Ref C 1.75 GHz）よりずっと広い。OPO は Ref B/C 用だが、Liu 2026 型の
   「ns パルスで励起する OPO」を F の光源の候補にするかどうかはユーザが決めること（§9）。

## 1. 記号（symbols, English）

conventions.md の Ω（光のサイドバンド角周波数 [rad/s]）と衝突しないように、Ref B の規格化周波数には `Ω_B` を使う。
Ref B の `L`（共振器内の損失）はファイバ長 `L` と衝突するので `L_ic` と書く。Ref B の `γ` はファイバの非線形係数 `γ` と衝突するので `γ_B` と書く。

- `ω0`: carrier (signal/idler degenerate) angular frequency [rad/s]; pump at `2ω0`. `ν0 = ω0/2π`.
- `Ω`: sideband angular frequency of the squeezed field relative to `ω0` [rad/s]; `f_s = Ω/2π` [Hz] (ESA measurement frequency, Ref B: `f`).
- `T`: power transmittance of the output coupler (Ref B: 0.123; Ref C: 1 − 0.64 = 0.36 at 1550 nm).
- `L_ic`: round-trip intracavity power loss excluding `T` (Ref B: `L`); may be a callable `L_ic(x)` or `L_ic(P)`, and a slow function of time.
- `l`: round-trip optical path length of the cavity [m] (Ref B: ≈ 0.5 m). `FSR = c/l` [Hz].
- `κ_T = c·T/(2l)`, `κ_L = c·L_ic/(2l)`: amplitude decay rates through the coupler and through internal loss [rad/s] (good-cavity limit, `T + L_ic ≪ 1`).
- `κ = κ_T + κ_L`: total amplitude decay rate = cavity HWHM [rad/s]. Cavity linewidth FWHM in Hz: `Δν_c = κ/π`.
- `γ_B = 2κ = c(T + L_ic)/l`: Ref B's "cavity decay rate" (energy decay rate = FWHM in rad/s).
- `Ω_B = Ω/γ_B = 2π f_s/γ_B`: Ref B's normalized frequency (corrected, see §3.2). Then `4 Ω_B² = (Ω/κ)²`.
- `ρ = κ_T/κ = T/(T + L_ic)`: escape efficiency.
- `P`: pump power at the OPO [W]. `P_th`: oscillation threshold [W]. `x = √(P/P_th)` = `ε/κ`: normalized pump amplitude (0 ≤ x < 1).
- `ε`: parametric coupling rate (amplitude) [rad/s], `ε = κ_c|β_p|` in the notation of conventions §4 (`κ_c` there is the coupling constant; renamed here to avoid a clash with `κ`).
- `φ_p`: pump phase (conventions §4).
- `Δ = ω_cav − ω0`: cavity detuning [rad/s] (rotating frame at `ω0`, `H_0 = Δ â†â`). `d = Δ/κ`.
- `G = 1/(1 − x)²`: classical parametric gain of a probe transmitted through the OPO at zero detuning (Ref B Eq. 4).
- `η`: photodiode QE, `ξ`: homodyne visibility (mode-matching efficiency `ξ²`), `ζ`: propagation efficiency (Ref B).
- `η_tot = ρ ζ ξ² η`: total detection efficiency of the squeezed field (for Ref B).
- `R±(Ω)`: anti-squeezed (+) / squeezed (−) noise variance normalized to shot noise, `R = V/V_SN`, `V_SN = 1/2` (conventions §2).
- `R′±`: same, averaged over residual phase jitter with RMS `θ̃` (Ref B notation; `σ` in conventions §6).
- `r(Ω)`: effective squeeze parameter of a lossless OPO sideband pair, defined by `e^{2r(Ω)} = ((1+x)² + (Ω/κ)²)/((1−x)² + (Ω/κ)²)`.
- `θ_sq`: squeezing angle (conventions §3). `δθ_sq`: its shift due to detuning.
- `δl`: change of round-trip optical path length [m]. `T_c`: crystal temperature [K].
- `L_c`: crystal length [m]. `Δk`: phase mismatch [1/m].
- `C_dB`: detector circuit noise relative to shot noise [dB] (Ref B: −21.7 dB at 5.5 mW LO).

## 2. 支配方程式

### 2.1 量子 Langevin 方程式（規約 §4 に合わせた形）

キャビティ内の縮退モード â（ω0 で回る枠）に対して、conventions §4 の H_int = i(κ_c/2)(β_p â†² − β_p* â²) と離調 H_0 = Δ â†â から
（入出力理論、Gardiner & Collett 1985）:

  dâ/dt = −(κ + iΔ) â + ε e^{iφ_p} â† + √(2κ_T) â_in + √(2κ_L) b̂_in,
  â_out = √(2κ_T) â − â_in.

â_in は出力鏡から入る真空、b̂_in は内部損失から入る真空。ε = κ x、発振しきい値は x = 1（Δ = 0）。ε ∝ ポンプの振幅なので x = √(P/P_th)。

φ_p = 0 の枠で直交位相 q̂ = (â + â†)/√2、p̂ = (â − â†)/(i√2) に分けると:

  dq̂/dt = −(κ − ε) q̂ + Δ p̂ + …,  dp̂/dt = −(κ + ε) p̂ − Δ q̂ + …

Fourier 変換（conventions §8、d/dt → −iΩ）で行列 M(Ω) = [[κ − ε − iΩ, −Δ], [Δ, κ + ε − iΩ]] を作ると、

  X_out(Ω) = A(Ω) X_in(Ω) + B(Ω) X_L(Ω),  A = 2κ_T M⁻¹ − I,  B = 2√(κ_T κ_L) M⁻¹,

ショット雑音で規格化したスペクトル行列は S(Ω) = A A† + B B†、角度 θ の直交位相で測るスペクトルは S_θ(Ω) = Re[u_θᵀ S(Ω) u_θ]、u_θ = (cos θ, sin θ)ᵀ。
検出効率 η_d（ρ 以外の損失）を通すと S_det = η_d S_θ + 1 − η_d（conventions §2 の純損失）。
**φ_p ≠ 0 のときは全体を θ → θ − φ_p/2 で回す。** Δ = 0 では増幅される角が θ_amp = φ_p/2、squeezing 角が θ_sq = φ_p/2 + π/2 で、conventions §4 と一致する。
（この 2×2 の行列計算は、実装で numpy の `inv` だけでできる。下の閉じた式はこの行列の Δ = 0 の特別な場合。）

### 2.2 Δ = 0 の閉じた式（Ref B Eq. 2 の正しい形）

  R±(Ω) = 1 ± η_tot · 4x / ((1 ∓ x)² + (Ω/κ)²),  η_tot = ρ ζ ξ² η,  ρ = T/(T + L_ic)

（上の行列計算で Δ = 0 とすると、検出損失 0 で R± = 1 ± ρ·4x/((1∓x)² + (Ω/κ)²) が数値でも 10⁻¹⁵ まで一致することを確認した [AE]。）

同じ式を「純粋な squeezed 状態に純損失をかけたもの」として書ける:

  R±(Ω) = η_tot e^{±2r(Ω)} + (1 − η_tot),  e^{2r(Ω)} = ((1+x)² + (Ω/κ)²) / ((1−x)² + (Ω/κ)²).

損失がなければ R+R− = 1（各サイドバンドで最小不確定）。Ω → 0 では e^{r} = (1+x)/(1−x)、r = 2 artanh x。
**これは Ref A（単一パス OPA）の R± = L + (1 − L) e^{±2r} と同じ構造**なので、OPA と OPO のブロックは「r(Ω) と固有の取り出し効率 ρ を出す」部分だけを差し替え、
損失・位相ジッタ・検出は共通のブロック（T12、T04、T13）にまかせられる。

### 2.3 位相ジッタ（Ref B Eq. 1）

  R′± ≈ R± cos²θ̃ + R∓ sin²θ̃（Ref B Eq. 1）。厳密なガウス平均は conventions §6: R′ = (R− + R+)/2 ± (R− − R+)/2 · e^{−2θ̃²}。
  Ref B の値（θ̃ ≤ 3.9°）では両者の差は 0.01 dB 以下 [AE]。

### 2.4 離調による squeezing 角の回転（E09）

2.1 の M⁻¹ を Δ の 1 次まで展開すると（Ω ≪ κ、損失なし）:

- squeezing の楕円の軸の回転: **δθ_sq ≈ −(Δ/κ)/(1 + x²)**（Δ > 0 で θ_sq が減る向き）。
- 最小の雑音の劣化は Δ の 2 次（Ref B、x = 0.745、Δ/κ = 0.1 で −12.96 → −12.88 dB、検出損失なし）。
- キャビティを透過する probe（同じ周波数 ω0、同じモード）の位相の回転: **δφ_pr ≈ −Δ/(2κ)**（x によらない）。
  これは lock 1（透過 probe の強度の極値にポンプ位相を合わせる）がかかった状態の値で、増幅側・減衰側のどちらにロックしても同じ。
- したがって lock 2（LO を probe に合わせる）の後に残る squeezing 角の誤差は **δθ_res ≈ (Δ/κ)(1 − x²)/(2(1 + x²))**。
  x = 0.3 / 0.745 / 0.9 で 0.417 / 0.143 / 0.053 × Δ/κ。数値（2×2 行列で極値を探す）と 10⁻⁴ 以下の相対誤差で一致 [AE]。
- サイドバンド周波数への依存: Ref B、x = 0.745、Δ/κ = 0.01 で回転は f_s = 1 kHz / 1 MHz / 5 MHz で −0.369° / −0.360° / −0.252°。
  Ω ~ κ では squeezing 角が周波数で回る（周波数依存の squeezing 角）。

離調とキャビティ長の関係（conventions §4 の符号: 光路が伸びると位相が増える）:

  δΔ = −ω0 δl / l（往復の光路長 l が δl 伸びると共振周波数が下がる）。
  モノリシック共振器（Ref C）では長さは温度で決まる: δΔ/(2π) = −ν0 (α_L + (1/n) dn/dT) δT_c（α_L: 線膨張係数、dn/dT: 屈折率の温度係数）。

### 2.5 ポンプに依存する損失（E06）とポンプのパワー揺らぎ

- L_ic は x（または P）の関数: Ref B Eq. (3) L_ic(x) = 0.00249 + 0.00222x。ρ(x) = T/(T + L_ic(x))、κ(x) も L_ic(x) を通して x に依存する。
- Ref B §3: ポンプなしで L_ic = 0.0020、100 mW で 0.0038。波長を 約 1 nm 短くすると、100 mW 注入から 2000 s 後に 0.0073 以上まで時間とともに増えた。
  → L_ic(x, t) として、遅いドリフトを時間ループ段で持たせる。
- ポンプのパワーの揺らぎ: δx/x = δP/(2P)。Ω = 0 で ∂ ln(R+ − 1)/∂x = 1/x + 2/(1 − x)。
  Ref B（x = 0.745）で δ(R+ − 1)/(R+ − 1) ≈ 3.4 δP/P [AE]。しきい値の近くほど敏感（E03、T01 と共通の雑音源）。
- 速さの限界: 増幅される直交位相の緩和率は κ(1 − x)。Ref B（x = 0.745）で κ(1 − x)/2π ≈ 1.5 MHz。
  これより遅いポンプの揺らぎ（時間ループ段の DC–kHz）は準静的に R±(x(t)) で扱ってよい。

### 2.6 位相整合の温度ドリフト（E05、OPO 側）

- 位相整合が Δk だけずれると非線形結合が sinc(Δk L_c/2) 倍になり、しきい値は P_th(T_c) = P_th,0 / sinc²(Δk(T_c) L_c/2)、x = √(P/P_th(T_c))。
  Δk(T_c) ≈ (∂Δk/∂T)(T_c − T_pm)。位相整合点の近くでは x の変化は (T_c − T_pm) の 2 次。
- 結合の位相も e^{iΔk L_c/2} だけずれるが、これはポンプ位相の静的なオフセットと同じで、lock 1（ポンプ–probe）が吸収する。
- ∂Δk/∂T（PPKTP、860/430 nm と 1550/775 nm）の値は今回確認していない [UNVERIFIED]。パラメータとして持つ。
- モノリシック（Ref C）では、同じ温度が共振（離調、2.4）と位相整合の両方を動かす。ポンプも共振するので（Ref C Fig. 1: 775 nm の反射率 98% / 99.93%、
  図中の 111 MHz は 775 nm の線幅と読める [Ref C Fig. 1]）、2 波長の同時共振の条件で温度が決まる（Tohermes 2024 も「crystal temperatures を使って共振を保った」と書いている）。

## 3. Ref B Eq. (2) の確認（prompt §1 と §7-9 が求めている項目）

使った値 [Ref B §2–4]: T = 0.123、l ≈ 0.5 m、L_ic = 0.0038（100 mW）、η = 0.998、ξ = 0.988、ζ = 0.99、θ̃ = 1.5°、f = 1 MHz（ESA、RBW 30 kHz、VBW 300 Hz）、
P_th = 180 ± 2 mW、回路雑音 −21.7 dB（LO 5.5 mW）。報告値: 観測 −9.01/+15.12 dB、回路雑音を引いて −9.22/+15.15 dB、Eq. (1) で θ̃ を除いた生成値 −10.12/+15.15 dB、
Eq. (2) による全損失の推定 0.0709 ± 0.0045、測定値からの 1 − ρζξ²η = 0.0645。

導いた量 [AE]: ρ = 0.9700、η_tot = 0.9355（1 − η_tot = 0.0645、論文と一致）、γ_B = 7.60×10⁷ s⁻¹、FWHM = γ_B/2π = 12.1 MHz、HWHM = 6.05 MHz、FSR = 600 MHz。
f = 1 MHz で 4Ω_B² = (f/HWHM)² = 0.0273（2π あり）、6.9×10⁻⁴（2π なし）。

### 3.1 x の定義

§0-1 のとおり、論文は一貫して x = √(P/P_th)。サーベイの E10・§3・§6-4 の「本文は x = P/P_th」は誤り（テキスト抽出の誤り）。
参考に x = P/P_th（= 0.556）を入れると R′− = −8.03 dB、R′+ = +10.10 dB で、観測値 −9.0/+15.1 dB からまったく外れる。

### 3.2 Ω の 2π

2.1 の導出から、Eq. (2) の 4Ω² は (Ω/κ)² = (2π f_s/κ)²、すなわち Ω_B = 2π f_s/γ_B（γ_B = 2κ = c(T + L)/l）。論文の「Ω = f/γ」は、f を角周波数と読むか 2π が抜けている。
決め手は Fig. 3 の x → 1 の振る舞い:

| x → 0.999、θ̃ = 1.5° | R′− [dB] | R′+ [dB] |
|---|---|---|
| Ω_B = 2πf/γ_B | −7.63 | +21.4 |
| Ω_B = f/γ_B | +5.8 | +37.4 |
| Fig. 3 実線（図から読んだ値） | 約 −7.5 | 約 +21.5 |

点線（θ̃ = 3.9°）も x = 1 で約 −1.5 dB（計算 −1.43 dB）。**2π を入れた形だけが Fig. 3 と合う。** 論文の曲線は実際には 2π を入れて計算されている。
（Ω = 0.08 程度でも、しきい値の近くでは (1 − x)² と同じ大きさになるので効く。）

### 3.3 100 mW の点（G = 18.7 を含む）

| x の決め方 | x | G | R− | R+ | R′− | R′+ |
|---|---|---|---|---|---|---|
| √(P/P_th)、2π あり | 0.7454 | 15.4 | −10.34 | +14.95 | −9.44 | +14.95 |
| G = 18.7 から、2π あり | 0.7688 | 18.7 | −10.54 | +15.64 | −9.45 | +15.63 |
| √(P/P_th)、2π なし | 0.7454 | 15.4 | −10.73 | +16.39 | −9.42 | +16.39 |
| 報告値（回路雑音を引いた値） | — | 18.7 | −10.12 | +15.15 | −9.22 | +15.15 |

（単位 dB。L_ic = 0.0038、θ̃ = 1.5°、f = 1 MHz。）

- 回路雑音 −21.7 dB を足して比べると、モデル（1 行目）は −9.21/+14.92 dB、観測は −9.01/+15.12 dB。差は 0.2 dB。
- 報告された生成値 R± = −10.12/+15.15 dB から (x, η_tot) を逆に解くと（比 (R+ − 1)/(1 − R−) は η_tot によらないので x が決まる）:
  2π あり → x = 0.7532（P = 102 mW、G = 16.4）、全損失 0.0708。2π なし → x = 0.7124（P = 91 mW、G = 12.1）、全損失 0.0708。
  全損失は論文の 0.0709 と一致し（Ω の扱いによらない）、x は 2π ありのときだけ P = 100 mW と合う。
- **G = 18.7 は説明できない。** P_th の不確かさ ±2 mW は G で ±3% にしかならず、15.4 と 18.7 の 21% の差を説明できない。
  G の測り方（probe の透過パワーの比など）は本文に書かれていない。推測で補わない（UNRESOLVED）。G は Eq. (2) に直接は入らず、P_th を決めるのに使われただけなので、
  検証では x = √(P/P_th)、P_th = 180 mW を使えばよい。

### 3.4 L_ic(x) の変数

Eq. (3) の x は、論文の定義どおりなら √(P/P_th)。ただし L_ic(0.745) = 0.00414 で、100 mW での実測 0.0038 と 0.0003 ずれる。
x² = P/P_th を入れると 0.00372 で実測に近い。また、ポンプ誘起の吸収はポンプ強度（∝ x²）に比例すると考えるのが自然 [AE、UNVERIFIED]。
Fig. 4 の実線（x_opt に沿った L_ic）は、θ̃ → 0 で L_ic → 0.0047（= L_ic(x = 1)）となり、どちらの読み方とも矛盾しない。
**検証目標は論文の定義（√）で作る。** L_ic(x²) にすると、R′−,min は 1.5° で −9.46 dB（−0.11 dB）、3.9° で −7.25 dB（−0.08 dB）変わる。x_opt はほとんど変わらない。

### 3.5 論文の理想値の確認

損失なし、Ω → 0 では R′− = R− cos²θ̃ + sin²θ̃/R− を R− について最小にすると R′−,min = sin 2θ̃。
θ̃ = 3.9° → −8.67 dB（論文 −8.66 dB）、θ̃ = 1.5° → −12.81 dB（論文 −12.81 dB）。Eq. (1) の使い方と一致。

## 4. 仮定と適用範囲

- しきい値より下（0 ≤ x < 1）。ポンプの枯渇なし（squeezed vacuum）。線形化した揺らぎ、ガウス状態。
- good-cavity 近似: T + L_ic ≪ 1 で κ_T = cT/(2l)。Ref B（T = 0.123）では −ln(1 − T) = 0.131 との差が 7%。
  **Ref C（T = 0.36）では c(T + L)/l から線幅を出すと 24% ずれる**（−ln(0.64) = 0.446）。Ref C では線幅（κ）を当てはめ値で直接与える。
  確認: Ref C の設計（R = 99.98% と 64%、2.6 mm）は、n = 1.73 なら FWHM = FSR·(−ln R1R2)/2π = 2.37 GHz、n = 1.816 なら 2.26 GHz [AE、n の値は UNVERIFIED]。
- 単一縦モード: Ω ≪ FSR（Ref B: FSR 600 MHz に対し測定 1 MHz、Ref C: FSR 約 32 GHz に対し測定 ≤ 1.3 GHz）。
  隣の縦モードの組（ω0 ± m·FSR）も 2 モード squeezing を作る。Liu 2026 はこれを MZI で除いている。Eq. (2) はこれを含まない。
- 離調は小さい（|Δ| ≪ κ）として 2.4 の 1 次の式を使う。大きい離調は 2.1 の行列で計算する。
- ポンプは共振しない（Ref B、本文に明記はない [UNVERIFIED]）か、ポンプのキャビティの線幅が雑音の周波数よりずっと広い（Ref C: 111 MHz）。
- 時間ループ段の雑音（DC–kHz）は準静的: 雑音の周波数 ≪ κ(1 − x)。
- ESA の RBW は線幅よりずっと狭い（Ref B: 30 kHz ≪ 6 MHz）。Ref C のように線幅と同じ程度のスパンを測るときは、各 f_s で R±(f_s) を評価する。

## 5. 周波数帯とスケーリング

| 効果 | 周波数帯 | スケーリング |
|---|---|---|
| E08 キャビティのスペクトル | サイドバンド 0 – 数 κ（Ref B: HWHM 6.05 MHz、Ref C: HWHM 0.875 GHz） | (Ω/κ)²、x、ρ |
| E09 離調・長さ雑音 | 雑音は DC–kHz（音響、長さロックの残差、温度）。効果は Ω ≪ κ の全域で squeezing 角の回転 | δθ_sq ∝ δl·ω0/(lκ)/(1 + x²)、probe との差 ∝ (1 − x²)/(2(1 + x²)) |
| E06 ポンプ誘起の損失 | 準静的（x の関数）＋遅いドリフト（s–h） | L_ic(x)、ρ = T/(T + L_ic) |
| E05 位相整合の温度 | < Hz | sinc²(Δk L_c/2) を通して x |
| ポンプのパワー雑音 | DC – κ(1 − x)/2π | δx/x = δP/(2P)、しきい値の近くで 1/(1 − x) で増える |

- **ファイバ長**: OPO の squeezing そのものはファイバ長によらない。ファイバに入る V_in(Ω) を決める。
- **パワー**: x = √(P/P_th)。
- **温度**: 位相整合（2.6）と、モノリシックでは離調（2.4）。
- **パルス幅**: CW 型 OPO の時間モードはキャビティの寿命（1/(2κ)、Ref B 約 13 ns、Ref C 約 91 ps）で決まり、パルス幅では決まらない。
  D07 のパルス（100 fs–10 ps、変換限界のガウスで 4.4 THz–44 GHz）は OPO の線幅の外にある。

## 6. 概算（Ref B/C の値と、我々の構成）

### 6.1 Ref B（§3 の表に加えて）

- 離調: ν0 = 349 THz、l = 0.5 m で dν/dl = 0.70 MHz/nm、**Δ/κ = 0.115 /nm**（往復の光路長）。
  δθ_sq ≈ −0.074 rad/nm（−4.2°/nm、x = 0.745）、probe でロックした後の残差 ≈ 0.0165 rad/nm（0.94°/nm）[AE]。
  残差を θ̃ = 1.5° の 1/3（0.5°）以下にするには、OPO の長さロックの残差（lock 1・lock 2 の帯域より遅い成分）が往復の光路長で 約 0.5 nm rms 以下である必要がある [AE]。
  Ref B は OPO の長さロックを本文に書いていない（Fig. 1 に「Lock Beam」の表示があるだけ）。値は時間ループ段のパラメータとする。
- ポンプ誘起の損失の時間ドリフト: L_ic が 0.0038 → 0.0073 になると ρ は 0.970 → 0.944、R′−（x = 0.745、θ̃ = 1.5°）は −9.44 → −8.61 dB。
  ドリフトの速さは ≥ 1.75×10⁻⁶ /s（「2000 s で 0.0073 以上」からの下限）[Ref B §3 から AE]。

### 6.2 Ref C（線幅 1.75 GHz FWHM）

- Fig. 4（建物 B、ファイバの後、ポンプを下げた測定）を、η_tot = 0.79（当てはめの損失 21%）、θ̃ = 0.1 rad、FWHM 1.75 GHz、x は 30 MHz の反 squeezing 10.0 dB から決めて（x = 0.56）計算すると [AE]:

| f_s | 30 MHz | 0.25 GHz | 0.69 GHz | 1.21 GHz |
|---|---|---|---|---|
| 反 squeezing 計算（FWHM と読む） | 10.0（合わせた点） | 8.7 | 5.0 | 2.6 |
| 反 squeezing 計算（HWHM と読む） | 10.0 | 9.6 | 7.8 | 5.6 |
| 反 squeezing Fig. 4 | 10.0 | 8.5 | 5.2 | 2.8 |
| squeezing 計算（FWHM） | −4.3 | −4.35 | −3.5 | −2.2 |
| squeezing Fig. 4 | −4.7 | −4.7 | −3.0 | −2.2 |

  FWHM の読み方だけが合う（差 ≤ 0.5 dB。dither による周波数依存のデコヒーレンス（Ref C §4、文献 [22]）は入れていない。T04/T14 の範囲）。
- ファイバの前（Ref C §2）: 30 MHz で 7 dB / 16 dB、損失 15% → η_tot = 0.85、16 dB から x ≈ 0.75、位相雑音なしで squeezing −7.7 dB [AE]。測定の 7 dB と矛盾しない。
- 離調（温度）: (α_L + (1/n) dn/dT) を 10⁻⁵ /K の桁とすると [AE、KTP の値は UNVERIFIED]、δν_c/δT ~ 2 GHz/K、HWHM 0.875 GHz に対し Δ/κ ~ 2 /K。
  δθ_sq（x = 0.56）~ 1.7 rad/K。10 mrad 以下にするには 約 6 mK 以下の温度安定度が必要。
  Ref C は probe を使わない量子雑音ロックなので、実際の楕円を見ている（2.4 の probe との差は生じない）。

### 6.3 Ref D（D02）

OPO も squeezer もない（coherent 光の MZI）。T02 は関係しない。

### 6.4 将来の構成 F（D07–D13、D30–D32）

- D07 の基本は単一パスの導波路 OPA（T01）。OPO は Ref B/C 用に持つ。将来の値はないので、この節はリファレンスの値での見積もり。
- **CW 型 OPO（Eq. 2）と D07 のパルス幅は両立しない**: 100 fs–10 ps の帯域は 4.4 THz–44 GHz で、Ref C の 1.75 GHz の線幅でも 25 倍以上狭い。
- Liu 2026（D29 の文献）の光源は「1.6 ns の矩形パルスで励起する OPO」4 台、PDH ロック、パルス間隔 50 ns（20 MHz）、
  ポンプが squeezed 光と同じくキャビティを通るのでポンプ位相のロックの probe に使える、OPO の位相ロックの精度は 約 1 nm [Liu 2026 Methods]。
  非縮退の縦モードは 3 段の不平衡 MZI と 0.4 nm のフィルタで除く。F にこの型の光源を採るなら、T02 のモデルに
  （i）パルス励起の時間モード（キャビティ寿命と励起パルス幅の関係）、（ii）縦モードの櫛とフィルタの損失、が必要（§9 の Q1）。
- 時間多重（D08）への目安: CW 励起の OPO を時間ビン T_bin で切り出すなら、時間モードがビンに収まるために 1/(2κ) ≪ T_bin が必要。
  f_rep = 100 MHz（T_bin = 10 ns）では Ref C 型（寿命 91 ps）は使えるが、Ref B 型（寿命 13 ns）は使えない [AE]。
- 1 nm の長さ精度（Liu）は 1550 nm で 2π·1/1551 ≈ 4 mrad。Ref C 型の OPO（l_opt ≈ 9.4 mm、Δ/κ ≈ 2π·ν0·δl/(l κ)）では
  1 nm の往復光路長の変化で Δ/κ ≈ 0.023、δθ_sq ≈ 0.018 rad（x = 0.56）[AE]。

## 7. 検証テスト（実装が通らなければならないもの）

すべて「この担当の計算値」を目標にし、論文の図との比較は読み取りの精度（約 ±0.3 dB）で別に行う。

- **V-B1 閉じた式と行列の一致**: Δ = 0 で、2.1 の行列計算の S_θ（θ = θ_sq, θ_sq + π/2）が 2.2 の R± と相対誤差 10⁻¹² 以内で一致する（Ref B の値、f_s = 0, 1, 5, 20 MHz）。
- **V-B2 純度と極限**: ρ = 1、η_tot = 1 で全 Ω について R+R− = 1（10⁻¹²）。x = 0 で R± = 1（真空は真空のまま、§7-4 の OPO 版）。
  Ω → ∞ で R± → 1。η_tot に対して R± − 1 が線形（§7-3 と同じ純損失）。
- **V-B3 Ref B Fig. 3 の再現（prompt §7-9 の訂正版の目標）**: f = 1 MHz、T = 0.123、l = 0.5 m、η = 0.998、ξ = 0.988、ζ = 0.99、P_th = 180 mW、
  L_ic(x) = 0.00249 + 0.00222x（x = √(P/P_th)）、Ω_B = 2πf/γ_B、Eq. (1)。許容差は 0.02 dB（同じ式なので数値誤差だけ）。

  | P [mW] | x | θ̃ = 1.5°: R′− / R′+ [dB] | θ̃ = 3.9°: R′− / R′+ [dB] |
  |---|---|---|---|
  | 10 | 0.236 | −3.66 / +3.89 | −3.58 / +3.88 |
  | 20 | 0.333 | −5.11 / +5.63 | −4.93 / +5.61 |
  | 50 | 0.527 | −7.72 / +9.47 | −6.92 / +9.46 |
  | 65 | 0.601 | −8.51 / +11.15 | −7.17 / +11.14 |
  | 90 | 0.707 | −9.25 / +13.87 | −6.69 / +13.85 |
  | 100 | 0.745 | −9.35 / +14.94 | −6.22 / +14.93 |
  | 120 | 0.816 | −9.20 / +17.08 | −4.93 / +17.06 |
  | 130 | 0.850 | −8.98 / +18.11 | −4.19 / +18.10 |
  | 150 | 0.913 | −8.35 / +19.96 | −2.70 / +19.94 |
  | — | 0.999 | −7.63 / +21.43 | −1.43 / +21.41 |

  最小: θ̃ = 1.5° で x_opt = 0.758（許容 ±0.01）、R′−,min = −9.35 dB。θ̃ = 3.9° で x_opt = 0.607（±0.01）、R′−,min = −7.17 dB。
  最小から 0.1 dB 以内の範囲は 1.5° で x = 0.708–0.806、3.9° で x = 0.558–0.654。
  （この表の 100 mW の値は L_ic(x) = 0.00414 を使う。実測の L_ic = 0.0038 を使うと −9.44 / +14.95 dB、§3.3。）
  Fig. 3 の実線・点線との比較は ±0.3 dB。論文の「x ≈ 0.82 より上で劣化」は測定点の記述なので、目標にはしない。
- **V-B4 Ω の定義の判別**: 同じ条件で x = 0.999 の R′+ が +21.4 dB（2π あり）であること。+37 dB（2π なし）になったら Ω の定義の誤り。
- **V-B5 100 mW の点と逆算**: L_ic = 0.0038、θ̃ = 1.5° で R′± = −9.44 / +14.95 dB。回路雑音 −21.7 dB を足して −9.21 / +14.92 dB（観測 −9.01 / +15.12 dB に 0.25 dB 以内）。
  生成値 −10.12 / +15.15 dB からの逆算で x = 0.753 ± 0.002、全損失 0.0708（論文 0.0709）。
- **V-B6 理想の極限**: 損失なし、Ω = 0 で min_x R′− = 10 log10(sin 2θ̃): 3.9° → −8.67 dB、1.5° → −12.81 dB（論文 −8.66、−12.81、許容 0.02 dB）。
- **V-B7 離調**: Ref B の値、x = 0.745、Δ/κ = 10⁻³ で δθ_sq = −(Δ/κ)/(1 + x²)（相対 10⁻³）。lock 1 の極値で決めた probe の出力位相の回転が −Δ/(2κ)、
  probe に合わせた LO と楕円の差が (Δ/κ)(1 − x²)/(2(1 + x²))。符号: 往復光路長を伸ばす（δl > 0）と Δ < 0、θ_sq は増える。
  Δ/κ = 0.1、f_s = 1 kHz で最小の雑音 −12.88 dB、回転 −3.69°（検出損失なし、ρ = 0.970）。
- **V-C1 Ref C の線幅**: FWHM 1.75 GHz、η_tot = 0.79、θ̃ = 0.1 rad、30 MHz の反 squeezing 10.0 dB で x を決めたとき、0.25 / 0.69 / 1.21 GHz の反 squeezing が
  8.7 / 5.0 / 2.6 dB（±0.05 dB、計算値）。Ref C Fig. 4 の 8.5 / 5.2 / 2.8 dB とは ±0.5 dB。
- **V-B8 パラメータが関数でもよいこと**: L_ic を定数で渡しても x（または P）の callable で渡しても同じ結果になる（L_ic(x) ≡ 0.0038 のとき V-B5 と一致）。

## 8. シミュレータに必要なもの

| 名前（案） | 段 | 入力 | 出力 | toggle |
|---|---|---|---|---|
| `opo_squeezing_spectrum` | O | T、L_ic（定数 / callable(x, P, t)）、l または κ（線幅を直接）、P、P_th（または x）、Δ、f_s の配列 | r(Ω)、ρ、κ、R±(f_s)（検出前）、2×2 の S(Ω)、δθ_sq(Δ, x, f_s) | `opo.frequency_dependence`（Lorentz / Ω = 0）、`opo.pump_dependent_loss` |
| 共通の損失・検出（T12/T13） | O | η_tot の残り（ζ、ξ²、η）、回路雑音 | R±（検出後） | 既存 |
| `opo_detuning_noise` | T（雑音源） | 長さ雑音 S_δl(f) または温度 T_c(t)、ω0、l、κ、x | Δ(t)、δθ_sq(t)、probe の位相 δφ_pr(t) = −Δ/(2κ) | `opo.detuning_noise` |
| `opo_loss_drift` | T（遅いドリフト） | L_ic(x, t) のモデル（速さ、上限） | ρ(t)、κ(t) | `opo.loss_drift` |
| `opo_phase_matching` | T（遅いドリフト） | T_c(t)、∂Δk/∂T、L_c、T_pm | P_th(t) → x(t) | `opo.phase_matching_drift` |
| ポンプのパワー雑音 → x(t) | T | S_RIN（T01 と共通） | x(t) = √(P(t)/P_th) | T01 と共通 |
| ESA の評価（T04） | M | R±(f_s)、RBW | 測定値 | 既存 |

- **光路の位相のモデル（§5.2、D11）への要求**: OPO キャビティの区間は「共通の ψ を足すだけ」の区間にできない。
  離調 Δ(t) に対して、squeezed 光には δθ_sq = −(Δ/κ)/(1 + x²)、同軸の probe には −Δ/(2κ) を、別の係数で与える必要がある。
  （ポンプ位相 φ_p から θ_sq = φ_p/2 + π/2 は conventions §4 のまま。）
- パラメータは callable を許す（L_ic(x)、P_th(T_c)）。x は P から導くか直接与えるか選べるようにし、P_th = 180 mW を Ref B の設定に書く。G = 18.7 は設定に入れない。
- Ref C の設定では線幅を当てはめ値（1.75 GHz FWHM）で直接与え、c(T + L)/l からは計算しない（§4）。

## 9. 既存の決定・規約・サーベイとの食い違い（調整役への提案）と、ユーザの決定が必要な点

提案（調整役が判断）:

- **P1（survey.md）**: E10 の「本文は x = P/P_th」と §3「Ref B の G」、§6-4 を訂正。x はすべて √(P/P_th)。残る問題は G = 18.7 だけ（UNRESOLVED、§3.3）。
  E08 の「γ_c = c(T + L)/l ≈ 7.6×10⁷ s⁻¹」は FWHM [rad/s] で、Ω = 2πf/γ_c が正しい。
- **P2（prompt §1・§7-9 の Eq. 2）**: 「Ω = f/γ」を「Ω = 2πf/γ（f は Hz、γ = c(T + L)/l は FWHM の角周波数）」として使う。§7-9 の目標は §7 の V-B3 の表で置き換える。
- **P3（prompt §5.5、§7-11、D11）**: 「同軸 = 区間の共有なので共通の雑音は自動で打ち消し合う」は、OPO キャビティの中では離調について成り立たない。
  §7-11 のテストは「区間の位相が同じ係数で効く場合」に限るか、OPO の区間にビームごとの係数を持たせる（§8）ことを計画に入れる。
- **P4（conventions.md）**: OPO の記号（κ = HWHM の振幅減衰率 [rad/s]、線幅 FWHM Δν_c = κ/π [Hz]、x = √(P/P_th)、Δ = ω_cav − ω0）と、
  文献の Ω（規格化周波数）・γ（減衰率）・L（共振器内損失）を書き換える規則を §1 に足す。

ユーザの決定が必要な点（推測で埋めない）:

- **Q1**: 将来の構成 F の光源に、Liu 2026 型の「ns パルスで励起する OPO」を候補として入れるか（D07 は単一パス OPA が基本）。
  入れるなら T02 に、パルス励起の時間モードと縦モードのフィルタを足す必要がある。
- **Q2**: §7-9 の合格基準を、この担当が計算した曲線（V-B3、許容 0.02 dB）にするか、Fig. 3 の図から読んだ値（±0.3 dB）にするか、両方にするか。
  （この担当の推奨は両方: 実装の正しさは V-B3、論文との一致は図で。）

調整役か後の Phase で決めればよい未解決の点:

- L_ic(x) の変数（√ か P/P_th か、§3.4）。影響は R′−,min で約 0.1 dB。
- Ref B の OPO の長さロックの残差（本文に記述なし）、PPKTP の ∂Δk/∂T、KTP の dn/dT と α_L（いずれもパラメータとして掃引する）。

## 10. 文献

読んだもの:

- Ref B: Y. Takeno, M. Yukawa, H. Yonezawa, A. Furusawa, "Observation of −9 dB quadrature squeezing with improvement of phase stability in homodyne measurement,"
  Opt. Express **15**(7), 4321–4327 (2007). DOI: 10.1364/OE.15.004321（Crossref で確認）。PDF を画像にして式・図を確認した。
- Ref C: S. Verclas, C. Darsow-Fromm, B. Tohermes, R. Schnabel, "Fiber distribution of phase-stabilized GHz-bandwidth squeezed vacuum states of light between two buildings,"
  Optica Quantum **4**, 1 (2026). DOI: 10.1364/OPTICAQ.567418（PDF に記載）。
- B. Tohermes, S. Verclas, R. Schnabel, "Directly measured squeeze factors over GHz bandwidth from monolithic ppKTP resonators," arXiv:2412.03221 (2024)
  （arXiv の HTML 版を読んだ: 線幅 γ₂ = 1.750 ± 0.013 GHz は「full-width-half-maximum」、η_esc ≈ 99%、温度で共振を保つ）。
  Ref C は掲載誌を "New J. Phys. 27, 094506 (2024)" と書いているが、巻と年が合わず、確認できなかった **[UNVERIFIED]**。
- Liu et al., "Gaussian boson sampling with 1,024 squeezed states in 8,176 modes," Nature **653** (2026)（PDF のページ表記 "Nature | Vol 653 | 21 May 2026" から）。
  Methods の OPO と位相ロックの節を読んだ。著者の全リスト・ページ・DOI は確認していない **[UNVERIFIED]**。

書誌だけ確認したもの（本文は未読。2.1–2.2 の式はこの担当が入出力理論から導き、Ref B の式と数値で照合した）:

- M. J. Collett, C. W. Gardiner, "Squeezing of intracavity and traveling-wave light fields produced in parametric amplification," Phys. Rev. A **30**, 1386–1391 (1984). DOI: 10.1103/PhysRevA.30.1386（Crossref）。
- C. W. Gardiner, M. J. Collett, "Input and output in damped quantum systems: Quantum stochastic differential equations and the master equation," Phys. Rev. A **31**, 3761–3774 (1985). DOI: 10.1103/PhysRevA.31.3761（Crossref）。
- E. S. Polzik, J. Carri, H. J. Kimble, "Atomic spectroscopy with squeezed light for sensitivity beyond the vacuum-state limit," Appl. Phys. B **55**, 279–290 (1992). DOI: 10.1007/BF00325016（Springer の記録、検索で確認）。
- P. K. Lam, T. C. Ralph, B. C. Buchler, D. E. McClelland, H.-A. Bachor, J. Gao, "Optimization and transfer of vacuum squeezing from an optical parametric oscillator," J. Opt. B **1**, 469–474 (1999). DOI: 10.1088/1464-4266/1/4/319（検索で確認）。
- T. C. Zhang, K. W. Goh, C. W. Chou, P. Lodahl, H. J. Kimble, "Quantum teleportation of light beams," Phys. Rev. A **67**, 033802 (2003). DOI: 10.1103/PhysRevA.67.033802（Crossref）。
- S. Ast, M. Mehmet, R. Schnabel, "High-bandwidth squeezed light at 1550 nm from a compact monolithic PPKTP cavity," Opt. Express **21**(11), 13572 (2013). DOI: 10.1364/OE.21.013572（DOI が opg.optica.org の oe-21-11-13572 に解決することを確認）。

書誌も確認していないもの（Ref B の文献リストからの引き写し、**[UNVERIFIED]**）:

- S. Wang, V. Pasiskevicius, F. Laurell, "Dynamics of green light-induced infrared absorption in KTiOPO4 and periodically poled KTiOPO4," J. Appl. Phys. 96, 2023–2028 (2004)。
- H. Mabuchi, E. S. Polzik, H. J. Kimble, "Blue-light-induced infrared absorption in KNbO3," J. Opt. Soc. Am. B 11, 2023–2029 (1994)。
