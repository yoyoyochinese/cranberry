# T06 fiber-gnlse-classical: ファイバ中の古典パルス伝搬（GNLSE / RK4IP）

- Phase: 1（deep-dive）、トピック T06（`notes/survey.md` §4）
- 日付: 2026-10-28
- 担当する効果: E20（伝搬損失）、E21（β2・β3）、E29（SPM）、E30（XPM）、E31（誘導ラマン・パルス内ラマン）、prompt §7-5（RK4IP の検証）
- 入力: `docs/squeezing_sim_prompt.md`、`DECISIONS.md`（D01–D33）、`notes/conventions.md`（D33 で承認済み）、`notes/survey.md`
- 依存: なし。T03（モードの重なり）、T07（量子雑音層の線形化の背景場）、T10（probe による XPM）、T12（損失収支）に値を渡す。
- 値の印: survey §0.3 と同じ。**[AE]** = この担当の概算（式と入力を併記）。**[sim]** = この担当が scratchpad の試作コード（numpy の RK4IP、§3）で計算した値。コードはコミットしない。**[UNVERIFIED]** = 一次資料で確認できなかった値・文献情報。

## 0. 要約（主な結論）

1. **GNLSE は Agrawal の規約で書けば `notes/conventions.md` §4・§8 とそのまま一致する。** 線形演算子は D̂(Ω) = i[β(ω0+Ω) − β0 − β1Ω] − α/2（伸びると位相が増える向き）、self-steepening は (1 + Ω/ω0)。既存の Python コード（gnlse-python 2.0.0、PyNLO 0.1.2）も Ã ∝ `ifft(A)` の同じ符号だが、**どちらも NumPy 2 系で動かない**（`np.math`、`np.trapz` を使用）。gnlse-python は固定刻みの RK4IP ではなく `solve_ivp`（RK45）なので、§7-5 の刻み半減テストに使えない。→ RK4IP は自前で書き（約 60 行）、既存コードは Raman 応答の式の照合にだけ使うことを勧める（§4）。
2. **survey の SPM の概算（100 fs で約 8 rad）は分散を無視した過大評価。** 100 fs では L_D ≈ 0.15 m なので、63 m の間にパルスが約 400 倍に広がり、非線形位相の大部分は最初の数 L_D でたまる。P_avg = 1 mW・f_rep = 100 MHz・L = 63 m の LO で、非線形位相（オーバーラップの位相）は 0.05–0.3 rad、最大は τ ≈ 1 ps 付近（L_D ~ L のとき）**[sim]**。f_rep = 10 MHz・L = 633 m・1 mW では、パルスエネルギーが 10 倍、長さも 10 倍で、**τ ≳ 0.3 ps は異常分散でソリトンになる**（N = 0.7–4）。非線形位相は数十 rad、スペクトルの形も大きく変わり、線形に伝搬した時間モードとの重なり |ov|² は 0.02–0.4 に落ちる（100 fs でも 0.60）**[sim]**（§8.4）。
3. **分散そのもの（線形）が 100 fs 級の時間モードを強く変える。** β2 L = −1.37 ps²（63 m）、−13.7 ps²（633 m）。ガウス型 100 fs パルスで「重なりの損失 1%」に相当する長さの合わせ残りは δL ≈ 5 cm。1 ps なら約 5 m、10 ps なら約 470 m **[AE]**。ポート間で L_k が最大 61 m 違うので、100 fs–1 ps では**ポートごとの分散補償か、LO_k を同じ種類・同じ長さのファイバで送ること（D09）が必須**になる（T03 へ）。
4. **逆向きの probe による XPM は小さいが、時間をずらしても避けられない。** 平均パワー 1 mW あたり、同じ偏光（b = 2）で 0.16 mrad（63 m）・1.5 mrad（633 m）、直交偏光（b = 2/3）でその 1/3 **[AE]**。逆向きなので、squeezed 光の各スライスは 2n_gL/c の窓の probe パワーの平均を受ける（sinc 型の低域通過、最初の零点は 63 m で 1.6 MHz、633 m で 161 kHz）。窓は T_rep よりずっと長いので、**AOM で probe を時間的にずらしても（Liu 型）XPM は平均パワーで決まり、消えない**。
5. **Ref D（CW、3.33 mW）と Ref C（CW）では GNLSE は事実上不要。** CW では GNLSE が厳密に A(L) = A0 e^{−αL/2} e^{iγP L_eff} に縮退する。Ref D の SPM は 1 腕・10 m あたり約 2×10⁻⁵ rad **[AE]**。量子サイドバンド（Ref C の GHz 帯）に対する β2 の効果は、直交位相の角度を β2Ω²L/2 だけ回すことで、1 km・1 GHz で 4×10⁻⁴ rad **[AE]**。
6. **Raman のモデルの選び方で結果が変わる。** Blow–Wood（f_R = 0.18）では、ソリトン自己周波数シフトが Gordon の式と 0.5 ps で 0.4%、0.2 ps で 1.5% の差で一致する **[sim]**。Lin–Agrawal（f_R = 0.245）では同じ条件で 1.5–2.2 倍大きい。これは h_b の 1 次モーメントが 0 なので、T_R の定義（1 次モーメント）に入らないためである **[sim, 解析]**。また **h_R を時間格子でサンプルすると、Δt ≳ τ1/4 ≈ 3 fs で偽の遅延が出る**（試作で確認）。h̃_R(Ω) は解析式で与えるべき（§3.4）。

## 1. 記号（symbols, English）

- `z` [m]: propagation distance. `T = t − β1 z` [s]: retarded time (conventions §1). `L` [m]: fiber (segment) length.
- `A(z,T)` [√W]: slowly varying envelope, `|A|²` = instantaneous power. Field `E ∝ A e^{i(β0 z − ω0 t)}` (positive-frequency part, conventions §4).
- `Ã(z,Ω) = ∫ A e^{+iΩT} dT` [√W·s]: spectrum (conventions §8). `Ω` [rad/s]: offset from carrier `ω0 = 2πc/λ0`.
- `β(ω)` [rad/m]: propagation constant. `β_k = d^kβ/dω^k|_{ω0}`: `β1 = n_g/c` [s/m], `β2` [s²/m], `β3` [s³/m].
- `D = −(2πc/λ²) β2` [s/m², usually ps/(nm·km)]: dispersion parameter. `S = dD/dλ` [ps/(nm²·km)]: dispersion slope.
- `α` [1/m]: power attenuation, `α = α_dB · ln(10)/10`, `α_dB` [dB/m]. `L_eff = (1 − e^{−αL})/α`.
- `γ = n2 ω0 /(c A_eff)` [1/(W·m)]: nonlinear coefficient (total Kerr + Raman, see §2.3). `n2` [m²/W]: nonlinear index. `A_eff` [m²]: effective mode area, `≈ π (MFD/2)²` for a Gaussian mode.
- `R(t) = (1 − f_R) δ(t) + f_R h_R(t)`: nonlinear response. `f_R`: fractional Raman contribution. `h_R(t)` [1/s]: causal Raman response, `∫ h_R dt = 1`. `h̃_R(Ω) = ∫ h_R e^{+iΩt} dt`.
- `τ1, τ2, τ_b, f_b, f_c`: parameters of the Raman models (§2.3). `T_R = f_R ∫ t h_R(t) dt` [s]: first moment of the Raman response.
- `τ_shock = 1/ω0` [s]: self-steepening time (constant-A_eff approximation).
- `P0` [W]: peak power. `P_avg` [W]: average power. `f_rep = 1/T_rep`: repetition rate. `E_p = P_avg/f_rep` [J]: pulse energy.
- `τ_FWHM` [s]: intensity FWHM. `T0` [s]: `τ_FWHM/1.7627` (sech²) or `τ_FWHM/1.6651` (Gaussian, intensity ∝ e^{−T²/T0²}). sech²: `P0 = 0.8814 E_p/τ_FWHM`.
- `L_D = T0²/|β2|`, `L'_D = T0³/|β3|`, `L_NL = 1/(γ P0)`, `N² = L_D/L_NL`: soliton order. `z0 = (π/2) L_D`: soliton period.
- `φ_NL` [rad]: nonlinear phase. `b`: XPM coefficient (2 same polarization, 2/3 orthogonal linear polarization in a high-birefringence fiber).
- `P_pr(t)` [W]: probe power. `v_g = c/n_g`. `T_w = 2 n_g L / c` [s]: counter-propagation averaging window.
- `N_ph ∝ ∫ |Ã(Ω)|²/(ω0 + Ω) dΩ`: photon number. `E = ∫ |A|² dT` [J]: pulse energy.
- `h` [m]: RK4IP step. `N_t`: number of time-grid points. `Δt` [s]: time step. `T_win = N_t Δt`: time window.
- `L_k ≈ k c T_rep / n_g` (D08), `δL_k`: residual length mismatch between squeezed fiber k and LO fiber k (D09).

## 2. 支配方程式

### 2.1 GNLSE（時間領域、この規約）

スカラーの GNLSE（Agrawal の規約、Dudley et al. 2006 Eq. (3) と同じ形）:

```
∂A/∂z = −(α/2) A + Σ_{k≥2} (i^{k+1} β_k / k!) ∂^k A/∂T^k
        + i γ (1 + i τ_shock ∂/∂T) [ A(z,T) ∫_0^∞ R(T′) |A(z, T − T′)|² dT′ ]        (1)
```

- 各項の意味: 第 1 項は損失（E20）、第 2 項は分散（E21、k = 2, 3 が β2, β3）、第 3 項は非線形の項。R の中の δ 関数の部分が SPM、h_R の部分が誘導ラマンとパルス内ラマン（E31）、(1 + iτ_shock ∂_T) が self-steepening。
- **符号の確認（conventions §4・§8）**: §8 の Fourier 変換では ∂/∂T → −iΩ なので、線形の部分は周波数領域で

```
∂Ã/∂z = D̂(Ω) Ã,   D̂(Ω) = −α/2 + i Σ_{k≥2} β_k Ω^k / k! = −α/2 + i[β(ω0+Ω) − β0 − β1 Ω]           (2)
```

  となる（i^{k+1}(−i)^k = i）。Ω = 0 の成分は位相を変えず（β0 z は場の定義に入っている）、Ω の成分は exp(+iβ(ω0+Ω)z) の位相を受ける。これは conventions §4「長さ L の区間を通ると A → A e^{+iψ}、ψ = βL」と同じ向きである。
- 停留位相で確認: Ã ∝ exp(iβ2Ω²L/2) から A(T) の T は T = β2 Ω L。異常分散（β2 < 0）では Ω > 0（青い側）が T < 0（早く着く）。物理的に正しい。§9 の V6 でテストする。
- self-steepening: (1 + iτ_shock ∂_T) → (1 + Ω/ω0)。周波数領域の非線形の項は iγ(1 + Ω/ω0)·F{…} で、Dudley et al. 2006 と gnlse-python の実装（`W = V + w_0`、`gamma/w_0`）と同じ。
- A_eff の周波数依存（τ_shock の補正）は入れない。帯域 Δω ≪ ω0（100 fs の sech² パルスでも Δω/ω0 ≈ 0.3%）なので無視できる **[AE]**。

### 2.2 相互作用表示（RK4IP で解く形）

A_I(z) = exp(−(z − z′)D̂) A(z) とおくと、∂A_I/∂z = exp(−(z − z′)D̂) N̂(A) A。D̂ は周波数領域で対角なので exp(hD̂/2) は FFT 1 往復で厳密に計算できる。分散と損失の誤差は 0 で、刻みの誤差は非線形の項だけから出る（§3）。

### 2.3 Raman 応答のモデル（E31）

すべて因果的（t < 0 で 0）で、∫h_R = 1 に規格化する。

| モデル | h_R(t) | パラメータ | 出典 |
|---|---|---|---|
| Blow–Wood（BW、単一の減衰振動子） | h_a(t) = (τ1² + τ2²)/(τ1 τ2²) · e^{−t/τ2} sin(t/τ1) | f_R = 0.18、τ1 = 12.2 fs、τ2 = 32 fs | Blow & Wood 1989。数値は gnlse-python 2.0.0 と PyNLO 0.1.2 の 2 つの実装が一致することを確認。原論文と Dudley et al. 2006 の本文は未読 [値は UNVERIFIED against primary] |
| Lin–Agrawal（LA、Boson peak を含む） | h_R = (f_a + f_c) h_a(t) + f_b h_b(t)、h_b(t) = (2τ_b − t)/τ_b² · e^{−t/τ_b} | f_R = 0.245、f_a = 0.75、f_b = 0.21、f_c = 0.04、τ_b = 96 fs、h_a は BW と同じ τ1・τ2 | Lin & Agrawal 2006。f_R = 0.245（n2 = 2.6×10⁻²⁰ m²/W のときピークのラマン利得が実測と合うように選ぶ）は検索で得た本文の抜粋による。f_b・f_c・τ_b は gnlse-python 2.0.0 と PyNLO 0.1.2 の 2 つの実装が一致することだけを確認 [UNVERIFIED against primary] |
| Hollenbeck–Cantrell（13 振動モード） | 中間ガウス型の和 | f_R ≈ 0.2（gnlse-python の実装） | Hollenbeck & Cantrell 2002。必要になるまで実装しない |

解析的なスペクトル応答（§3.4 で使う。conventions §8 の e^{+iΩt} の規約。h̃(0) = 1）:

```
h̃_a(Ω) = [(τ1² + τ2²)/(τ1² τ2²)] / [ (1/τ2 − iΩ)² + 1/τ1² ]
h̃_b(Ω) = (1 − 2iΩτ_b) / (1 − iΩτ_b)²
```

- **Raman 利得の符号の確認**: この規約ではストークス側（Ω < 0）が利得を受ける。§9 の V4 で「Raman をオンにするとスペクトルの重心が Ω < 0 へ動く」ことを確認する。試作では BW・LA とも重心が負の向きに動くことを確認した。
- **T_R**（1 次モーメント、Agrawal の T_R の定義）: BW で f_R⟨t⟩ = 1.46 fs、LA で 1.57 fs **[sim, 数値積分]**。h̃_b = 1 + x² + 2ix³ + …（x = Ωτ_b）なので h_b の 1 次モーメントは厳密に 0 で、LA の利得の傾きの追加分は Ω³ で立ち上がる。このため Gordon の SSFS の式（T_R を使う）は LA ではサブ ps のソリトンで成り立たない（§0 の 6、§9 の V9）。Agrawal の教科書でよく引かれる「1.55 µm で T_R ≈ 3 fs」とは 2 倍の差がある [3 fs の出典は UNVERIFIED]。
- **γ の意味**: (1) では CW・長いパルスの非線形位相は γ[(1 − f_R) + f_R ∫h_R] P = γP。つまり γ は**Kerr とラマンを合わせた全体の値**で、測定値（または n2 から求めた値）をそのまま入れる。Raman をオフにするときは f_R → 0 として γ はそのままにし、CW の非線形位相が変わらないようにすることを勧める（Phase 3 の設計判断、§12）。

### 2.4 保存量

- α = 0 で self-steepening をオンにすると、Raman の有無によらず光子数 N_ph ∝ ∫|Ã|²/(ω0 + Ω) dΩ が保存される（Blow & Wood 1989、Dudley et al. 2006 [保存の主張は UNVERIFIED against primary、試作で数値的に確認]）。エネルギー E は Raman で減る（赤方偏移する分がフォノンに渡る）。
- self-steepening がオフなら、Raman があっても E が局所的に保存される（非線形の項が A × 実数の形になるため）。光子数は厳密には保存されない。
- 試作の結果 **[sim]**: N = 2 ソリトン（T0 = 100 fs、β2 = −21.7 ps²/km、γ = 1.3 /(W·km)）を 3z0 伝搬し、800 ステップ、Raman ＋ shock の場合: |ΔN_ph/N_ph| = 8×10⁻⁷（BW）、7×10⁻⁷（LA）。ΔE/E = −3.6%（BW）、−4.6%（LA）。重心のずれは −6.9 THz（BW）、−8.8 THz（LA）。

### 2.5 XPM（E30）

1. **同じ向きに伝搬し、周波数が近く同じ偏光の光**（例: Ref A のように squeezed 光と 1 MHz ずれた probe が同じ経路を通る場合）: 同じ包絡線 A の中の成分として (1) にそのまま含まれる。別の XPM の項は要らない（入れると二重に数える）。量子雑音層での扱い（Kerr による squeezing の回転・FWM）は T07 の担当。
2. **偏光が直交する光（PM ファイバの slow 軸と fast 軸）**: 軸ごとにスカラーの (1) を解き、相手の軸の強度で iγ b |A_other|² A を加える（b = 2/3）。複屈折が大きい（ビート長 L_B ≤ 5 mm → Δβ = 2π/L_B ≥ 1.3×10³ rad/m）ので、コヒーレントな結合の項（A_other² A*）は位相が合わず、落とす [b = 2/3 と落とせる条件は Agrawal の教科書による、UNVERIFIED]。
3. **逆向きに伝搬する光（Liu 型の probe、D30）**: 定常波のグレーティングの項は位相が合わないので平均で消え、強度の項だけが残る（同じ偏光なら b = 2、PM の直交偏光なら b = 2/3）。probe は squeezed 光の影響をほとんど受けないので、連成した GNLSE は解かず、**外から与えた強度の項**として扱う:

```
∂A_s/∂z ⊃ i b γ P_pr(z, t) A_s
φ_XPM(t_in) = b γ ∫_0^L P_pr,in( t_in + (2z − L)/v_g ) e^{−α(L−z)} dz                         (3)
```

  t_in は squeezed 光のスライスがファイバに入る時刻、P_pr,in は z = L に入る probe のパワーである。α = 0 なら φ_XPM = (bγ v_g/2) ∫_{t_in − L/v_g}^{t_in + L/v_g} P_pr,in(t′) dt′、つまり**長さ T_w = 2n_gL/c の窓での平均**になる。probe のパワーの揺らぎ δP(f) から位相への伝達関数は

```
H_XPM(f) = b γ L · sinc(π f T_w) · e^{iπ f (T_w − 2 n_g L/c)}   (α = 0、sinc(x) = sin x / x)          (4)
```

  （位相の因子は t_in を基準にした遅延。|H| の最初の零点は f = 1/T_w。）
  probe を AOM でパルス化しても、周期 ≪ T_w なら φ_XPM は**平均パワー**で決まる。各スライスは T_w の窓のなかに入ったすべての probe 光と交差するので、「squeezed パルスと時間をずらす」ことでは避けられない（ずらすことの効果は GBS 側の検出時刻での混入光子数、T10）。
4. **SMF（ランダムな複屈折）**: 偏光の関係が長さ方向でランダムに変わるので、b は 2/3 から 2 の間になる。長いファイバで平均すると、SPM の実効係数は 8γ/9（Manakov-PMD、Marcuse, Menyuk & Wai 1997）[8/9 の係数は書誌のみ確認、本文は未読、UNVERIFIED]。この調査では b をパラメータとして持つことを勧める。

### 2.6 CW での縮退（prompt §5.1）

A(z,T) = A(z) が T によらないとき、(1) の ∂_T の項はすべて 0、∫R|A|² = |A|² なので

```
A(L) = A(0) · exp(−αL/2) · exp(i γ P(0) L_eff)                                                 (5)
```

は厳密解である。`LASER_MODE = CW` では FFT を使わず (5) を使う（CW の近道）。同じ入力を時間格子にのせて RK4IP で解いたものと一致すること（§9 の V5）。逆向きの CW probe の XPM は (3)・(4)、量子サイドバンドへの分散の効果は §8.7 に示す。

### 2.7 他の規約からの変換

| 出典・コード | 規約 | この規約への変換 |
|---|---|---|
| Agrawal「Nonlinear Fiber Optics」、Dudley et al. 2006 | E ∝ A e^{i(β0z − ω0t)}、Ã = ∫A e^{+iΩT}dT、(1) と同じ | 変換なし |
| Travers, Frosz & Dudley 2010（MATLAB `gnlse.m`）、gnlse-python 2.0.0 | `AW = ifft(AT)`、`AT = fft(AW)`、`L = 1i*B − alpha/2`、`W = V + w0` | 符号は同じ。規格化は §8 の Ã = NΔt·ifft(A) と比べて 1/(NΔt) 倍（gnlse-python は出力で `N*dt` を掛けて戻す）。gnlse-python の時間格子は `linspace(−T/2, T/2, N)` で端点を両方含むので、Δt = T_win/(N − 1)、周期は NΔt ≠ T_win（小さなずれ）。損失は dB/m で `alpha = ln(10^(loss/10))`、§1 と同じ |
| PyNLO 0.1.2 | `FFT_t = ifftshift(ifft(fftshift(A)))`（e^{+iΩT}）、`exp_D = exp(h/2·(i·betas − alpha/2))` | 符号は同じ。fftshift の扱いと規格化はラッパーで §8 に揃える。単位は ps・km・W |
| 工学の規約（E ∝ A e^{+jω0t}、j = −i） | 位相の符号が逆 | A_ours = A_src*、方程式の i → −i。β_k の値はそのまま |
| 量子光学の文献でよくある Ã = ∫A e^{−iΩT}dT | Ω の符号が逆 | Ω → −Ω（奇数次の β_k と self-steepening、Raman 利得の向きの符号が変わる） |

## 3. 数値法: RK4IP（Hult 2007）と誤差の振る舞い

### 3.1 アルゴリズム

刻み h、D̂ は (2)、N̂(A)A は (1) の非線形の項（周波数領域で iγ(1 + Ω/ω0)F{A·[(1 − f_R)|A|² + f_R h_R ⊛ |A|²]}）。

```
A_I = e^{(h/2)D̂} A(z)
k1  = e^{(h/2)D̂} [h N̂(A(z)) A(z)]
k2  = h N̂(A_I + k1/2) (A_I + k1/2)
k3  = h N̂(A_I + k2/2) (A_I + k2/2)
k4  = h N̂(e^{(h/2)D̂}(A_I + k3)) · e^{(h/2)D̂}(A_I + k3)
A(z + h) = e^{(h/2)D̂} [A_I + k1/6 + k2/3 + k3/3] + k4/6
```

（Hult 2007。PyNLO 0.1.2 の `RK4IP` がこの式どおりに実装していることを確認した。1 ステップあたり FFT は 8 往復程度。）

### 3.2 誤差の振る舞い

- 線形の部分は厳密なので、大域誤差は非線形の項からの O(h⁴) である。刻みを半分にすると、漸近的に誤差は 1/16 になる（prompt §7-5）。
- 試作の結果 **[sim]**（N = 2 ソリトン、z0 まで、基準は 6400 ステップ、2^12 点、T_win = 80 T0）:

| ステップ数 | 50 | 100 | 200 | 400 | 800 |
|---|---|---|---|---|---|
| 相対 L2 誤差 | 1.0×10⁻² | 9.6×10⁻⁴ | 7.3×10⁻⁵ | 5.1×10⁻⁶ | 3.4×10⁻⁷ |
| 前の行との比 | — | 10.7 | 13.2 | 14.2 | 15.2 |

  比が 16 に近づくのは、1 ステップの非線形位相 γP0h が十分小さくなってから（ここでは h ≤ z0/200、つまり 1 ステップあたり約 0.03 rad 以下）。粗い刻みでは比が 16 より小さい。**検証の合格の基準は漸近域で判定する**（§9 の V3）。
- 誤差の下限: 窓の端での打ち切り・エイリアス（時間と周波数の両端のエネルギーの割合で監視）、倍精度の FFT の丸め（1 ステップ約 10⁻¹⁶·√N_t の大きさ、ステップ数とともに増える）。
- 自動刻み: (a) ステップ倍化による局所誤差の推定（Sinkin et al. 2003。PyNLO はこれと Richardson 外挿 16/15·A_fine − 1/15·A_coarse を使う。外挿すると 5 次になり、§7-5 の 16× テストの対象ではなくなる）、(b) 光子数の保存の誤差による刻みの制御（Heidt 2009 の CQE 法）。§7-5 のテストは**固定刻み・外挿なし**で行う。本番の計算では (a) か (b) を使ってよい。
- 刻みの目安: h ≲ 0.03 L_NL(z) = 0.03/(γ max|A|²)。分散は厳密なので L_D による制限はない。試作ではこの規則で刻みを z ごとに変え、§8.4 の表を計算した。

### 3.3 格子の要件

- T_win: 分散で広がったパルス（≈ T0√(1 + (L/L_D)²)）の 10 倍以上。Raman の応答（h_R は約 1 ps で減衰）の周回の畳み込みが巻き戻らないこと。
- Δt: スペクトル全体（SPM・ソリトンの圧縮で広がった分を含む）が ±π/Δt の内側に収まること。§3.4 の Raman の条件も満たすこと。
- 監視: 出力で |A|² と |Ã|² の両端 5% の点に入るエネルギーの割合を記録し、10⁻⁶ を超えたら格子を広げる。

### 3.4 Raman の畳み込みの実装（試作で見つけた落とし穴）

- 畳み込み h_R ⊛ I は FFT で h̃_R(Ω)·Ĩ(Ω)。このとき h_R(t) を時間格子でサンプルしてから FFT すると、Δt ≳ τ1/4（約 3 fs）で h_R の振動（周期 2πτ1 ≈ 77 fs）が解像されない。h_R ≈ δ(t − Δt) の偽の遅延になり、偽の SSFS が生じる（Δt = 470 fs の試作では T_R が約 85 fs 相当に化けた）。ps パルスの計算では Δt を τ1 より細かくしないことが多いので、**h̃_R(Ω) は §2.3 の解析式で与える**（PyNLO も BW の解析式 `R0` を持っている）。
- 同じ理由で、gnlse-python の `raman_blowwood(T)` など、時間格子の上で h_R を作る関数は、Δt ≤ 約 3 fs の格子でだけ使う。

## 4. 既存の Python パッケージと、包み方

| パッケージ（PyPI の最新版、2026-10-28 確認） | 方法 | Raman | 規約 | NumPy 2.5 との互換 | 判断 |
|---|---|---|---|---|---|
| gnlse-python 2.0.0（WUST-FOG、MIT、Redman et al. 2021 arXiv:2110.00298） | 相互作用表示＋`scipy.integrate.solve_ivp`（既定 RK45、適応刻み）。依存に pyfftw、tqdm | BW、LA、Hollenbeck–Cantrell（時間格子で作る） | ifft が前進変換（§2.7） | **不可**: `dispersion.py` が `np.math.factorial` を使う。`np.math` は NumPy 2 で削除された（NumPy 2.5.4 で `hasattr(np, 'math') == False` を確認） | 固定刻みの RK4IP ではないので §7-5 の 16× テストができない。Raman の式の照合だけに使う |
| PyNLO 0.1.2（pyNLO/PyNLO、GPL-3） | SSFM と RK4IP（Hult）、Sinkin の適応刻み＋外挿 | BW（解析式と時間格子の両方）、LA | ifft が前進変換 | **不可**: `np.trapz` を使う（NumPy 2.5.4 では削除済み）。GPL-3 | RK4IP の式の照合だけに使う |

勧め（Phase 3 で判断）: RK4IP の本体は numpy の `fft` で自前に書く（約 60 行）。§8 の規約で直接書けばラッパーは要らない。既存のパッケージは、固定した古い環境（numpy < 2）で同じ入力を流して結果を比べる任意のクロスチェックにとどめる。Julia の Luna.jl など他言語の実装は調べていない。

## 5. 仮定と適用範囲

- スカラー・単一モード・緩やかに変わる包絡線（Δω ≪ ω0）。100 fs（Δω/ω0 ≈ 0.3%）から CW まで成り立つ。
- β(ω) は β2・β3 までの Taylor 展開。ガウス型 100 fs パルスのスペクトル幅（約 ±10 THz の範囲）では β4 の寄与は小さい [AE、β4 は資料がない]。必要なら β(ω) の表で置き換えられる形にしておく（gnlse-python の `DispersionFiberFromInterpolation` に相当）。
- γ は周波数によらない（A_eff の分散は無視、§2.1）。
- 損失は周波数によらない（C バンドの中では平坦）。曲げ損失・コネクタの損失は T12・T16。
- 偏光: PM は軸ごとのスカラー方程式。軸の間の結合（PER）は T09。SMF は Manakov の平均（§2.5）。
- 古典場の方程式。自発ラマン・GAWBS・量子雑音は T07（Q 段）。T06 の解は T07 の線形化の背景場になる。
- ブリルアン（SBS）は (1) に入らない（T07 の閾値チェック）。
- 温度: β2、γ、h_R の温度依存は弱く、ここでは無視する（ラマンの自発雑音の n_th(T) は T07、熱による位相ドリフトは T08）。

## 6. 周波数帯

| 効果 | 働く帯域 | 補足 |
|---|---|---|
| 損失 | DC（パワーの比） | 量子状態には純損失として効く（T12、conventions §2） |
| β2・β3 | 光のスペクトル幅（100 fs で約 THz） | MHz–GHz の量子サイドバンドには β2Ω²L/2 の回転として効く（§8.7）。小さい |
| SPM（LO・probe） | DC の位相の偏り＋強度雑音 → 位相雑音への変換（RIN の帯域全体） | パルスでは、パルスごとの P0 の揺らぎがそのまま位相の揺らぎになる |
| XPM（逆向きの probe） | DC ＋ RIN の f < 1/T_w の成分（63 m で約 1.6 MHz、633 m で約 160 kHz） | (4) の sinc 型の低域通過 |
| パルス内ラマン（SSFS） | ラマン利得の帯域（ピーク約 13 THz） | τ ≲ 1 ps でパルスのスペクトルがかかる |
| 誘導ラマン（CW） | ストークス側 約 13 THz | 閾値が W 級なので CW の mW 光では無関係（§8.5） |
| self-steepening | 時間尺度 1/ω0 ≈ 0.82 fs | 100 fs で 1% 程度の補正 |

## 7. スケーリング

- 損失: 透過率 e^{−αL}。L_eff → L（αL ≪ 1）。
- 分散: パルスの広がり T0√(1 + (L/L_D)²)、L_D ∝ τ²。GDD = β2 L ∝ L。重なりの許容 δL ∝ τ²（§8.3）。
- SPM: 分散を無視すれば φ_NL = γP0L_eff ∝ P_avg/(f_rep·τ)·L。分散が効く（L ≫ L_D、N ≪ 1）と、ガウス型の近似で φ_NL ≈ γP0L_D·asinh(L/L_D)。L への依存が対数になり、∝ P_avg τ/f_rep で、**短いパルスのほうが小さい**。最大になるのは L_D ~ L のとき。
- ソリトンの閾値: sech² の N = 1 のエネルギーは E_1 = 2P0T0 = 2|β2|/(γT0)、平均パワーは P_avg,1 = 2|β2| f_rep/(γT0) ∝ f_rep/τ（T0 = τ_FWHM/1.7627）。
- XPM: φ_XPM = bγP̄_pr L_eff ∝ P̄_pr·L。帯域 1/T_w ∝ 1/L。
- SSFS: Gordon の式 dΩ/dz = −8T_R|β2|/(15T0⁴) ∝ τ⁻⁴（BW。LA は §2.3）。
- 温度: このトピックの量にはほぼ依存しない（§5）。

## 8. 概算（この系のパラメータ）

### 8.1 ファイバのパラメータ（D32: PM と SMF を同じ扱いで並べる）

| 量 | SMF-28 Ultra 型 | PM1550（PANDA） | 出典・状態 |
|---|---|---|---|
| 減衰 α_dB @1550 | ≤ 0.18 dB/km（Ref C は約 0.2 dB/km を使用） | ≤ 0.5 dB/km（Coherent/Nufern PM1550B-XP、Fujikura 1550 nm 帯 PANDA）、≤ 1.0 dB/km（Thorlabs PM1550-XP の表） | SMF: Corning PI-1424 の値（PDF は 403 で取得できず、検索結果の抜粋で確認。0.18 は UNVERIFIED）。PM1550B-XP: https://www.coherent.com/components-accessories/specialty-optical-fibers/single-mode/PM1550B-XP（閲覧 2026-10-28）。Thorlabs・Fujikura: 検索結果の抜粋のみ（UNVERIFIED） |
| MFD @1550 | 10.4 ± 0.5 µm | 10.1 ± 0.4 µm（PM1550B-XP、PM1550-XP）、10.5 ± 0.5 µm（Fujikura） | 同上 |
| A_eff ≈ π(MFD/2)² | 84.9 µm² | 80.1 µm²（10.1 µm）、86.6 µm²（10.5 µm） | [AE] |
| n2 | 2.6×10⁻²⁰ m²/W | 同じ | Lin & Agrawal 2006 が使う値（検索結果の抜粋で確認）。石英の一般的な値として使う [UNVERIFIED] |
| γ = n2ω0/(cA_eff) @1548.8 nm | 1.24 /(W·km) | 1.32 /(W·km)（10.1 µm）、1.22 /(W·km)（10.5 µm） | [AE] |
| D @1550 | ≤ 18 ps/(nm·km)（仕様の上限）。計算には 17 を使う | **メーカーの資料に記載がない**。SMF と同じ値を仮定 | SMF: Corning PI-1424（抜粋）。PM: UNVERIFIED、要測定または要問い合わせ |
| β2（D = 17） | −21.7 ps²/km | −21.7 ps²/km（仮定） | [AE] |
| S @1550 | 約 0.056 ps/(nm²·km)（Corning 型の式 S0/4·(1 + 3λ0⁴/λ⁴)、λ0 = 1314 nm、S0 = 0.086 で 0.055） | 同じと仮定 | [AE, λ0・S0 は UNVERIFIED] |
| β3 | 0.126 ps³/km | 同じと仮定 | [AE] |
| n_g @1550 | 1.4682 | 未確認（1.468 で代用） | SMF: Corning PI-1424（抜粋） |
| ビート長 L_B | — | ≤ 5 mm（PM1550-XP の表） | 抜粋のみ、UNVERIFIED（T09） |

**PM の β2・β3・n_g はメーカーの資料にない。** D23 により SMF の値で代用し、代用を PROCESS_LOG に記録する必要がある（§12）。

### 8.2 損失（E20）

| | L_31 = 63.3 m（f_rep = 100 MHz） | L_31 = 633 m（f_rep = 10 MHz） |
|---|---|---|
| SMF 0.18 dB/km | 透過 0.9974（損失 0.26%） | 0.974（2.6%） |
| PM 0.5 dB/km | 0.9927（0.73%） | 0.930（7.0%） |
| PM 1.0 dB/km | 0.9855（1.45%） | 0.864（13.6%） |

[AE]。633 m・PM 1 dB/km の 13.6% は、Ref A の全損失 8% より大きい。目標 3–10 dB（D16）に対しては、ファイバの種類の選定（D32）で効く量である（T12 へ）。

### 8.3 分散（E21）: L_D と時間モード

| τ_FWHM（sech²） | T0 | L_D | L'_D（β3） | 63 m での広がり（ガウス型） | 重なりの損失 1% の δL（ガウス型） |
|---|---|---|---|---|---|
| 100 fs | 56.7 fs | 0.149 m | 1.4 m | 約 380 倍 | 4.7 cm |
| 300 fs | 170 fs | 1.34 m | 39 m | — | — |
| 1 ps | 567 fs | 14.9 m | 1.4 km | 3.9 倍 | 4.7 m |
| 3 ps | 1.70 ps | 134 m | 39 km | — | — |
| 10 ps | 5.67 ps | 1.49 km | 1400 km | 1.0 倍 | 470 m |

[AE]。重なりの式（ガウス型、GDD の差 δφ2 = β2 δL）: |⟨u1|u2⟩|² = 1/√(1 + (δφ2/(2T0²))²)。
- 100 fs–300 fs では β3 も L ≳ L'_D で効く（63 m ≫ 1.4 m）。
- ポートごとの L_k が違う（最大 61 m）ので、100 fs–1 ps では各ポートの squeezed パルスのチャープがまったく違う。GBS でポート間の干渉が必要なら、ポートごとの分散補償が要る。BHD の指標 (i) は LO_k が同じ長さ・同じ種類のファイバなら打ち消し合い、δL_k の許容は上の表の値になる。→ T03 に渡す。

### 8.4 LO・probe パルスの SPM（E29）と SSFS（E31）

条件: sech² パルス、P_avg = 1 mW（survey の仮定。LO・probe の平均パワーは D07 で Open）、SMF のパラメータ（γ = 1.24 /(W·km)、β2 = −21.7 ps²/km、β3 = 0.126 ps³/km、0.18 dB/km）、BW の Raman、self-steepening をオン。RK4IP で計算 **[sim]**。
φ_NL,peak = γ∫max_T|A|²dz。arg ov と |ov|² は、線形だけで伝搬したパルス（squeezed 光の時間モードの代わり）との重なり ⟨A_lin|A_NL⟩ の位相と大きさ。∂(arg ov)/∂ln P は、P_avg を 1% 変えたときの差から求めた値で、**LO の強度雑音（RIN）が BHD の位相雑音になる係数**である。

| f_rep, L | τ_FWHM | P0 [W] | N | φ_NL,peak [rad] | arg ov [rad] | |ov|² | ∂(arg ov)/∂ln P [rad] | rms スペクトル幅の比 | 重心の移動 [GHz] | 出力 FWHM [ps] |
|---|---|---|---|---|---|---|---|---|---|---|
| 100 MHz, 63.3 m | 0.1 ps | 88 | 0.13 | 0.16 | 0.106 | 0.9989 | 0.108 | 0.98 | −1.2 | 26.2 |
| | 0.3 ps | 29 | 0.22 | 0.30 | 0.208 | 0.9967 | 0.214 | 0.95 | −0.40 | 8.4 |
| | 1 ps | 8.8 | 0.40 | 0.41 | 0.286 | 0.9941 | 0.299 | 0.89 | −0.11 | 3.0 |
| | 3 ps | 2.9 | 0.70 | 0.22 | 0.151 | 0.9959 | 0.153 | 0.98 | −0.01 | 3.3 |
| | 10 ps | 0.88 | 1.28 | 0.069 | 0.046 | 0.9996 | 0.046 | 1.00 | ≈ 0 | 10.0 |
| 10 MHz, 633 m | 0.1 ps | 881 | 0.40 | 2.8 | 1.70 | 0.60 | 2.1 | 0.82 | −14 | 162 |
| | 0.3 ps | 294 | 0.70 | 72 | 2.18 | 0.023 | −17 | 0.53 | −42 | 0.74 |
| | 1 ps | 88 | 1.28 | 99 | 0.17 | 0.15 | 45 | 1.39 | −86 | 0.68 |
| | 3 ps | 29 | 2.21 | 50 | −2.28 | 0.11 | 24 | 3.46 | −21 | 0.79 |
| | 10 ps | 8.8 | 4.03 | 19 | −0.50 | 0.38 | 13 | 7.82 | −3.2 | 6.2 |

（格子: Δt = min(T0/12, 10 fs)、T_win ≥ 10T0√(1 + (L/L_D)²) + 40T0、刻み h = min(0.03 L_NL(z), L/100)。両端のエネルギーの割合は 10 MHz・0.1 ps で 3×10⁻⁶ で、V11 の基準をわずかに超える。他は ≤ 4×10⁻⁷。10 MHz の τ ≥ 0.3 ps の行は強い非線形の領域で、∂(arg ov)/∂ln P は 1% の差分による目安。arg ov の折り返しも含めて、数値は傾向を示すものとして読むこと。）

- 分散を無視した見積もり（survey の 100 fs で 6.9 rad）は、100 MHz・63 m では 40 倍以上の過大評価である。
- 10 MHz・633 m・1 mW: τ ≥ 0.3 ps では N ≥ 0.7 になり、ソリトンができる（パルスが圧縮され、スペクトルが変わり、SSFS で赤方偏移する）。LO の時間モードが squeezed 光と大きくずれ、RIN から位相への変換係数も 1 rad 以上になる。**この領域では LO の平均パワーを下げるか、正常分散のファイバ・分散補償を使う必要がある**（光源の仕様 D07 と LO の配送 D09 への入力）。
- N = 1 になる平均パワー（sech²、SMF）: 100 fs で 61 mW（100 MHz）/ 6.1 mW（10 MHz）、1 ps で 6.1 / 0.61 mW、10 ps で 0.61 / 0.061 mW **[AE]**。
- LO の RIN からの位相雑音: δφ_rms ≈ |∂(arg ov)/∂ln P| · σ_RIN。例: 100 MHz・1 ps で 0.3 × σ_RIN。パルスのエネルギーの揺らぎが 10⁻³ rms なら 0.3 mrad **[AE]**。

### 8.5 誘導ラマン（CW）

CW の閾値（Smith 1972 の 16 A_eff/(g_R L_eff) の形、g_R ≈ 6×10⁻¹⁴ m/W @1.55 µm は UNVERIFIED）: 63 m で約 360 W、633 m で約 36 W、1 km で約 23 W **[AE]**。Ref C・Ref D・CW probe（mW）では無関係。

### 8.6 逆向きの probe による XPM（E30、Liu 型、D30）

| | 63.3 m | 633 m |
|---|---|---|
| T_w = 2n_gL/c | 620 ns | 6.2 µs |
| (4) の最初の零点 1/T_w | 1.61 MHz | 161 kHz |
| φ_XPM / P̄_pr（同じ偏光 b = 2、SMF） | 0.157 mrad/mW | 1.55 mrad/mW |
| φ_XPM / P̄_pr（直交偏光 b = 2/3、PM 1 dB/km） | 0.055 mrad/mW | 0.52 mrad/mW |

[AE]。
- 一定の部分はロックの設定値に吸収される。問題になるのは揺らぎで、δφ_XPM = φ_XPM·δP̄_pr/P̄_pr（f < 1/T_w）。P̄_pr = 1 mW・633 m・パワーのドリフト 10% で 0.16 mrad。小さい。
- **ロックの誤差としての SPM と XPM の差**: Liu 型では probe の位相（probe 自身の SPM γP_pr,peak L_eff を含む）を固定し、squeezed 光は XPM（bγP̄_pr L_eff）を受ける。probe のパワーが変わると、ロックで固定される量と squeezed 光の位相が (b P̄_pr − P_pr,peak)γL_eff の割合でずれる。パルス化した probe（デューティ d）では P_pr,peak = P̄_pr/d なので、d が小さいほど大きい。例: P̄_pr = 1 mW、d = 0.1、633 m で SPM 約 7.7 mrad、その 1% の揺らぎは 0.08 mrad **[AE]**。これは T-段の「強度 → 位相」係数として持つ（§10）。

### 8.7 Ref D と Ref C（CW）

- **Ref D**（D02、3.33 mW、PM ファイバの MZI）: 1 腕あたり 1.67 mW として、SPM は γPL = 1.32×10⁻³ × 1.67×10⁻³ × L → **10 m あたり 2.2×10⁻⁵ rad** **[AE]**。腕の長さの差 ΔL（不明、D02）に比例する差分は、さらに小さい。RIN からの変換も無視できる。**GNLSE の段は CW の近道 (5) で済み、Ref D の結果には実質効かない。** 915B の巻線ファイバの長さは `914-915B-Specifications.pdf` からテキストで取り出せなかった（画像の PDF）。T15・T16 の値を待つ。
- **Ref C**（SMF 1 km、CW の LO）: LO の SPM は 1.2 mrad/mW。一定の部分はロックに吸収される。
- **量子サイドバンドへの分散**: ω0 ± Ω の 2 つのサイドバンドはどちらも位相 β2Ω²L/2（偶数次）を受けるので、Ω で測る直交位相の角度が θ_disp(Ω) = β2Ω²L/2 だけ回る。β3 の項（奇数次）は時間のずれで、角度を変えない。Ref C（1 km）で Ω/2π = 1 GHz なら 4.3×10⁻⁴ rad、2 GHz で 1.7×10⁻³ rad **[AE]**。Ref C の位相雑音の当てはめ値 0.1 rad に比べて無視できる。MHz 帯（ESA 0–100 MHz）ではさらに 4 桁以上小さい。

## 9. 検証テスト（実装が通るべきもの）

試作の数値は 2^12 点、T_win = 80 T0、β2 = −21.7 ps²/km、γ = 1.3 /(W·km)、T0 = 100 fs（L_D = 0.461 m、N = 1 の P0 = 1669 W）。基準値は試作で達成できた値に約 3 倍の余裕を持たせた。

| ID | テスト（prompt §7-5 ほか） | 設定 | 合格の基準 | 試作の値 [sim] |
|---|---|---|---|---|
| V1 | 基本ソリトンの不変性 | N = 1、α = β3 = f_R = 0、shock オフ、5z0、200 ステップ | 解析解 √P0 sech(T/T0) e^{iz/(2L_D)} との相対 L2 誤差 ≤ 1×10⁻⁵（位相も含む） | 2.9×10⁻⁶ |
| V2 | 2 次ソリトンの回帰 | N = 2、同上、z0、400 ステップ | max_T |I(z0) − I(0)| / max I(0) ≤ 1×10⁻⁵ | 2.6×10⁻⁶ |
| V3 | 刻み半減で誤差 約 16 分の 1 | N = 2、z0、固定刻み・外挿なし、ステップ数 200/400/800、基準は 6400 ステップ | 連続する誤差の比が 13–19（漸近域）。h を粗くした場合の比が 16 未満でも不合格にしない | 14.2、15.2 |
| V4 | Raman ありの光子数の保存 | N = 2、BW と LA、shock オン、α = 0、3z0、800 ステップ | |ΔN_ph/N_ph| ≤ 1×10⁻⁵。かつ ΔE < 0、スペクトルの重心の移動 < 0（符号のテスト） | 8×10⁻⁷（BW）、7×10⁻⁷（LA）。ΔE/E = −3.6%、−4.6%。重心 −6.9、−8.8 THz |
| V5 | CW での縮退 | CW 入力を時間格子にのせ、α・γ・f_R・β をすべてオン | RK4IP の結果と (5) の差 ≤ 1×10⁻¹²（相対）。CW の近道の経路と一致 | 未実施（解析的に自明） |
| V6 | 分散の符号と大きさ | 線形のみ（γ = 0）、ガウス型、L = 3L_D | 出力のスペクトル位相が β2Ω²L/2 + β3Ω³L/6 と一致（≤ 1×10⁻¹² rad）。β2 < 0 で青い側が先に着く（周波数分解した重心の時刻 T(Ω) = β2ΩL の符号）。幅が √(1 + (L/L_D)²) 倍 | 未実施 |
| V7 | 損失 | γ = 0、α > 0 | E(L)/E(0) = e^{−αL}（≤ 1×10⁻¹²）。dB/km → 1/m の換算（§1）のテスト | 未実施 |
| V8 | 分散なしの SPM | β = 0、f_R = 0、shock オフ | |A| は不変、位相 = γ|A(0,T)|²L_eff（相対 ≤ 1×10⁻⁸） | 未実施 |
| V9 | SSFS の速さ（Gordon） | BW、shock オフ、β3 = 0、N = 1、T0 = 0.5 ps、L = L_D | スペクトルの重心の移動が −8T_R|β2|L/(15T0⁴)（T_R はモデルの 1 次モーメント 1.46 fs）と 3% 以内。**LA ではこのテストを使わない**（§2.3） | 比 1.004（T0 = 0.5 ps）、1.015（0.2 ps） |
| V10 | 逆向き XPM の伝達関数 | 正弦波で変調した P_pr、長さ L | 位相の応答の振幅が |bγL sinc(πfT_w)| と ≤ 1% で一致、最初の零点 f = c/(2n_gL) | 未実施 |
| V11 | 格子の端の監視 | すべての本番の計算 | 時間・周波数の両端 5% のエネルギーの割合 ≤ 1×10⁻⁶ | §8.4 の計算で記録 |
| V12 | 規約のラッパー | 既存コード（gnlse-python・PyNLO）と比べる場合 | 同じ入力で |A|² が ≤ 1×10⁻⁶、スペクトルの重心の符号が一致 | 任意 |

## 10. シミュレータへの要求

### 10.1 ブロック・表（段の略号は survey §0.2）

| 名前（案） | 段 | 内容 |
|---|---|---|
| `fiber_gnlse`（オフラインのソルバ） | G | (1) を RK4IP で解く。項ごとにオン・オフ（§10.3）。CW では (5) に切り替える |
| `fiber_segment_params` | O | §8.1 の値を持つファイバの種類のプリセット（SMF、PM）。区間ごとの L、n_g、向き、偏光の軸（D11、conventions §9）と同じ区間の定義を共有する |
| `lo_nonlinear_table` | O | P_avg × τ × L_k の格子で: arg ov、|ov|²、∂(arg ov)/∂ln P、φ_NL,peak、スペクトルの重心と rms 幅、出力の時間幅（§8.4 の表を一般化したもの）。T03（モードの重なり）と T04（位相雑音）が読む |
| `fiber_linear_table` | O | 透過率 e^{−αL}、群遅延 β1L、GDD β2L、β3L。T12・T03・T08 が読む |
| `xpm_counterprop_coeff` | O → T | bγL_eff と T_w（(4) の伝達関数）。時間ループでは probe の RIN → squeezed 光の位相の雑音ブロックの係数になる |
| `spm_rin_to_phase` | O → T | LO・probe の RIN → 位相の係数（∂(arg ov)/∂ln P、CW では γPL_eff） |
| 背景場の受け渡し | Q | A(z,T) を z ごとに保存し、T07 の線形化の背景場として渡す |

### 10.2 入力

- 場: A(0,T)（パルスの形と τ、P_avg、f_rep、チャープ）、または CW のパワー。λ0。
- ファイバ: L、α_dB、β_k（または β(ω) の表）、γ（または n2 と A_eff）、n_g、Raman のモデルと f_R、偏光の扱い（PM の軸 / SMF の Manakov 係数）、XPM の b。
- 外部の強度: 逆向き probe の P_pr,in(t)（平均、変調、デューティ）。
- 数値: N_t、T_win（または Δt）、刻み h（固定）か許容誤差（適応）、保存する z の点の数。

### 10.3 オン・オフの旗（prompt §5.1、保護リスト）

`loss`、`beta2`、`beta3`（または `dispersion_order`）、`kerr`（瞬時の SPM）、`raman ∈ {off, blow_wood, lin_agrawal, hollenbeck_cantrell}`、`self_steepening`、`xpm_external`（逆向き probe）、`cw_shortcut ∈ {auto, on, off}`。
`raman = off` のとき γ の全体は保つ（§2.3、§12 の設計判断）。

### 10.4 出力

A(L,T)、Ã(L,Ω)、z ごとの保存点、§10.1 の表、保存量の記録（E、N_ph、格子の端のエネルギー）、設定と乱数の種と commit hash（prompt §6。GNLSE は決定論的だが、同じ記録の仕組みを使う）。

## 11. 他のトピックとの関係

- T03（パルスの時間モード）: §8.3 の GDD と重なりの許容、§8.4 の `lo_nonlinear_table` を渡す。
- T07（ファイバの量子雑音）: 背景場 A(z,T)、h̃_R(Ω) の規約（§2.3）を共有する。自発ラマンの雑音の周波数依存は Im h̃_R(Ω) で決まるので、モデルの選択（BW / LA）は T07 と揃える。
- T08（環境の位相雑音）: 群遅延 β1L と n_g（§8.1）。
- T09（偏光・PM）: ビート長、PM の β2 の値（資料がない）、直交偏光の XPM の係数 b = 2/3。
- T10（probe と squeezed 光の相殺の崩れ）: §8.6 の XPM と、probe の SPM と XPM の差。
- T12（損失）: §8.2 の透過率。

## 12. 規約・DECISIONS との関係、ユーザの判断が必要な点

規約との食い違い: **なし**。conventions §8 の Fourier 変換は Agrawal の規約そのもので、(2) の線形演算子は §4 の「伸びると ψ が増える」と同じ向きである。

調整役への提案（DECISIONS・conventions の追記の案。この担当は編集しない）:
1. conventions §8 に一文を追加する: 「GNLSE の線形演算子は D̂(Ω) = −α/2 + i[β(ω0+Ω) − β0 − β1Ω]。Raman の応答 h̃_R(Ω) は ∫h_R e^{+iΩt}dt で定義し、解析式で与える。」
2. D23 の代用の記録（PROCESS_LOG）: PM ファイバの β2・β3・n_g は資料にないため SMF-28 Ultra 型の値で代用。n2 = 2.6×10⁻²⁰ m²/W。LO・probe の P_avg = 1 mW（survey の仮定のまま）。
3. survey の E29 の概算（100 fs で約 8 rad）を §8.4 の値に置き換える。

ユーザの判断が必要な点（推測で決めない）:
1. **PM ファイバの分散の値**: メーカーの資料に β2 がない。SMF の値で代用してよいか。実際に使うファイバの品番（Ref D の PNH1550R5A2 のファイバを含む）と、メーカーへの問い合わせ・測定の可否。
2. **LO・probe の平均パワーとパルス幅**（D07・D28 が Open）: 10 MHz・数百 m では、1 mW の ps パルスの LO がソリトンになる（§8.4）。LO の平均パワーの上限、または分散補償・正常分散のファイバを候補に加えるか。
3. **Raman のモデルの既定値**: BW（Gordon の式で検証できる）と LA（実測のラマン利得に近いとされる）のどちらを既定にするか。T07 の自発ラマン雑音と揃える必要がある。
4. **Raman をオフにしたときの γ の扱い**（§2.3）: 全体の γ を保つ（勧め）か、瞬時の部分 (1 − f_R)γ だけを残すか。アブレーション（prompt §8）の解釈が変わる。

## 13. 文献

一次資料（書誌は Crossref で確認。本文を読んでいないものはそう書く）:
- J. Hult, "A fourth-order Runge–Kutta in the interaction picture method for simulating supercontinuum generation in optical fibers," J. Lightwave Technol. 25(12), 3770–3775 (2007). doi:10.1109/JLT.2007.909373。書誌は確認、本文は未読。アルゴリズムは PyNLO の実装（Hult を引用）で照合。
- J. M. Dudley, G. Genty, S. Coen, "Supercontinuum generation in photonic crystal fiber," Rev. Mod. Phys. 78(4), 1135–1184 (2006). doi:10.1103/RevModPhys.78.1135。書誌は確認、本文は未読（(1) の形と BW の値は記憶と 2 つの実装による。UNVERIFIED against primary）。
- K. J. Blow, D. Wood, "Theoretical description of transient stimulated Raman scattering in optical fibers," IEEE J. Quantum Electron. 25(12), 2665–2673 (1989). doi:10.1109/3.40655。書誌は確認、本文は未読。
- Q. Lin, G. P. Agrawal, "Raman response function for silica fibers," Opt. Lett. 31(21), 3086 (2006). doi:10.1364/OL.31.003086。書誌は確認、本文は未読（f_R = 0.245 は検索の抜粋、他の値は実装からのみ）。
- R. H. Stolen, J. P. Gordon, W. J. Tomlinson, H. A. Haus, "Raman response function of silica-core fibers," J. Opt. Soc. Am. B 6(6), 1159 (1989). doi:10.1364/JOSAB.6.001159。書誌は確認（f_R ≈ 0.18 の由来とされる。本文は未読、UNVERIFIED）。
- D. Hollenbeck, C. D. Cantrell, "Multiple-vibrational-mode model for fiber-optic Raman gain spectrum and response function," J. Opt. Soc. Am. B 19(12), 2886 (2002). doi:10.1364/JOSAB.19.002886。書誌は確認。
- J. P. Gordon, "Theory of the soliton self-frequency shift," Opt. Lett. 11(10), 662 (1986). doi:10.1364/OL.11.000662。書誌は確認（V9 の式の形 −8T_R|β2|/(15T0⁴) は Agrawal の教科書の形で、本文は未読。試作の数値で BW と 0.4% で一致）。
- F. M. Mitschke, L. F. Mollenauer, "Discovery of the soliton self-frequency shift," Opt. Lett. 11(10), 659 (1986). doi:10.1364/OL.11.000659。書誌は確認。
- O. V. Sinkin, R. Holzlöhner, J. Zweck, C. R. Menyuk, "Optimization of the split-step Fourier method in modeling optical-fiber communications systems," J. Lightwave Technol. 21(1), 61–68 (2003). doi:10.1109/JLT.2003.808628。書誌は確認。
- A. M. Heidt, "Efficient adaptive step size method for the simulation of supercontinuum generation in optical fibers," J. Lightwave Technol. 27(18), 3984–3991 (2009). doi:10.1109/JLT.2009.2021538。書誌は確認。
- R. G. Smith, "Optical power handling capacity of low loss optical fibers as determined by stimulated Raman and Brillouin scattering," Appl. Opt. 11(11), 2489 (1972). doi:10.1364/AO.11.002489。書誌は確認（閾値の係数 16 は UNVERIFIED）。
- D. Marcuse, C. R. Menyuk, P. K. A. Wai, "Application of the Manakov-PMD equation to studies of signal propagation in optical fibers with randomly varying birefringence," J. Lightwave Technol. 15, 1735–1746 (1997). doi:10.1109/50.622902。書誌は確認（8/9 の係数は UNVERIFIED）。

教科書・コード:
- G. P. Agrawal, *Nonlinear Fiber Optics*, 6th ed. (Academic Press, 2019). doi:10.1016/C2018-0-01168-8。書誌は確認、本文は未読。規約 (1)、XPM の b、T_R ≈ 3 fs、g_R の値は UNVERIFIED。
- J. C. Travers, M. H. Frosz, J. M. Dudley, "Nonlinear fibre optics overview," in J. M. Dudley, J. R. Taylor (eds.), *Supercontinuum Generation in Optical Fibers* (Cambridge University Press, 2010), pp. 32–51. doi:10.1017/CBO9780511750465.004。書誌は確認（MATLAB `gnlse.m` の出典。gnlse-python はこれの移植）。
- P. Redman, M. Zatorska, A. Pawłowski, D. Szulc, S. Majchrowska, K. Tarnowski, "gnlse-python: Open source software to simulate nonlinear light propagation in optical fibers," arXiv:2110.00298 (2021)。arXiv の API で確認。コード: gnlse 2.0.0（PyPI、https://github.com/WUST-FOG/gnlse-python、ソースを確認 2026-10-28）。
- PyNLO 0.1.2（PyPI、https://github.com/pyNLO/PyNLO、ソースを確認 2026-10-28）。

メーカーの資料（D24）:
- Coherent (Nufern) PM1550B-XP product page, https://www.coherent.com/components-accessories/specialty-optical-fibers/single-mode/PM1550B-XP（閲覧 2026-10-28）。
- Corning SMF-28 Ultra Product Information Sheet PI-1424, https://www.corning.com/media/worldwide/coc/documents/Fiber/product-information-sheets/PI-1424-AEN.pdf（2026-10-28 に取得を試みたが 403。値は検索結果の抜粋による。要確認）。
- Thorlabs PM1550-XP、Fujikura 1550 nm 帯 PANDA: 2026-10-28 の検索結果の抜粋のみ（製品ページは取得できず、UNVERIFIED）。

参照したリファレンス（`refs/`）: Ref C（Verclas et al. 2026、SMF 1 km・約 0.2 dB/km）、Liu et al. 2026（逆向きで直交偏光の probe、AOM でパルス化して時間をずらす。probe のパワーは本文に記載なし）。
