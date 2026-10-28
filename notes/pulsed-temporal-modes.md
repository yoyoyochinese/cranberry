# T03 パルスの時間・スペクトルモード（pulsed-temporal-modes）

- Phase: 1（deep-dive、トピック T03）
- 日付: 2026-10-28
- 担当する効果（`notes/survey.md`）: E11（パルス励起の多モード squeezing）、E19（タイミングジッタ）、E21 のモードの重なりへの帰結、E26 の帰結（熱によるタイミングドリフト）、E51 のタイミング部分、E54（時間・スペクトルモードの不整合）
- 依存: T06（ファイバの β2・β3・SPM）、T08（群遅延の熱係数）。どちらもまだ値がないので、下の値はプレースホルダ（**[PH]** と書く）。T01（OPA の利得・損失）とも接する。
- 規約: `notes/conventions.md` に従う。記号の定義は英語、説明は日本語（D19）。
- 値の印は survey §0.3 と同じ（[Ref A]、[AE]、[UNVERIFIED] など）。**[PH]** = 他トピックから受け取るまでのプレースホルダ。

## 0. 結論の要約

1. **Ref A 型の 45 mm PPLN（バルク LN で近似）では、1545/773 nm の群速度差（GVM）は 0.304 ps/mm、45 mm で τ_GVM ≈ 13.7 ps。**
   D07 の掃引範囲 100 fs – 10 ps はすべて τ_p ≲ τ_GVM で、出力は強い多モード（低利得の Schmidt 数 K ≈ 2×10² 以上）。
   さらに τ_p ≪ τ_GVM では有効な相互作用長が L_w = τ_p/|Δβ1| に縮むので、同じ r を得るのに必要なピークパワーが (τ_GVM/(√(2π) T0p))² 倍になる（100 fs で約 8×10³ 倍）[AE]。
2. **時間のずれの許容は非常に厳しい。** 時間モードの重なりで 99% を保つには |Δt| ≤ 0.142 T0（100 fs → 8.5 fs、1 ps → 85 fs、10 ps → 0.85 ps）。ファイバ長に直すと 1.7 µm / 17 µm / 174 µm。
3. **熱によるタイミングドリフト**: 39 ps/(km·K) [AE, T08 待ち] なら 63 m で 2.5 ps/K。重なり 99% のための温度差の許容は 3.5 mK（100 fs）/ 35 mK（1 ps）/ 0.35 K（10 ps）（補正なし、2 本の温度が無相関の最悪値）。
   ただし squeezed 光の経路にストレッチャを置いて位相をロックすると（D10）、群遅延の大部分も一緒に補正される。残る割合は (1 − R_τ,st/R_τ,th) で、survey の値では約 21% [AE]。**State は区間ごとに位相だけでなく群遅延も持つ必要がある**（§8、提案）。
4. **分散**: LO_k を同じ長さ・同じ種類のファイバで送れば（D09）BHD 側は相殺される。しかし **ポート間（指標 (ii)、GBS での干渉）では L_k − L_0 の分散が残る。** f_rep = 100 MHz（L_31 = 63 m）、β2 = −21.7 ps²/km [PH] で、補償なしの重なり η は 100 fs で 0.005、1 ps で 0.47、10 ps で 0.9998。
   η ≥ 0.99 には τ_FWHM ≳ 3.7 ps（100 MHz）/ 11.6 ps（10 MHz）/ 1.2 ps（1 GHz）、それより短いならポートごとの分散補償が必要。
5. **CW の極限**: τ_p → ∞ でポンプのスペクトルが δ 関数になり、超モードは ±Ω の側帯波の対に戻り、重なりの係数は 1、タイミングの状態は不要になる。LASER_MODE = CW では本トピックのブロックは恒等写像になり、prompt §7-1 の数値（−10.15 dB / +20.56 dB）を変えない。
6. Ref D（CW）・Ref A・B・C（CW または OPO の CW）では本トピックの効果は作用しない。将来の構成 F だけが対象。

## 1. 記号（symbols）

- `t`, `T`, `z`, `Ω`, `ω0`, `f`: as in `notes/conventions.md` §1. `T = t − β1 z` (retarded time in the frame of the signal group velocity).
- `τ_p`: intensity FWHM of a pulse [s]. `T0 = τ_p / (2√ln2) ≈ τ_p/1.665`: 1/e half-width of the **intensity** of a Gaussian pulse, intensity ∝ exp(−T²/T0²), field ∝ exp(−T²/(2T0²)).
- `T0p`, `τ_p,pump`: same for the pump pulse inside the OPA.
- `u(T)`: normalized temporal mode function, ∫|u|² dT = 1. `ũ(Ω) = ∫ u(T) e^{+iΩT} dT` (conventions §8).
- `u_LO,k(T)`: temporal mode of the LO pulse at the BHD of port k (after all segments). `u_s,k(T)`: mode the squeezed pulse would occupy if it were single-mode.
- `ψ_n(T)`, `r_n ≥ 0`, `φ_n`: n-th output supermode (Bloch–Messiah / Takagi mode), its squeezing parameter and squeeze phase (`φ_ξ` of conventions §3 for that mode). `θ_sq,n = φ_n/2`.
- `κ_n = ∫ u_LO*(T) ψ_n(T) dT`: overlap of the LO mode with supermode n (complex).
- `η_t = |κ|²`: temporal/spectral mode-matching efficiency (single-mode case).
- `Δt_k = τ_s,k − τ_LO,k`: arrival-time difference at the BHD of port k (squeezed later → positive) [s].
- `φ2 = φ2,s − φ2,LO`: net group-delay dispersion (GDD) difference between the squeezed and LO paths [s²]; for a fiber segment `φ2 = β2 L`. Likewise `φ3 = β3 L` [s³].
- `χ = φ2 / (2 T0²)`: normalized chirp mismatch.
- `β1,s`, `β1,p`: inverse group velocities of signal (ω0) and pump (2ω0) in the OPA. `Δβ1 = β1,p − β1,s` (GVM) [s/m]. `τ_GVM = |Δβ1| L_c`. `L_c`: crystal (waveguide) length. `L_w = τ_p/|Δβ1|`: walk-off length.
- `k_s(ω)`, `k_p(ω)`: propagation constants of signal and pump in the OPA. `β2,s`, `β2,p`: their GVD.
- `Δk(Ω1, Ω2) = k_p(2ω0 + Ω1 + Ω2) − k_s(ω0 + Ω1) − k_s(ω0 + Ω2) − K_QPM` (phase mismatch; `K_QPM` chosen so that `Δk(0,0) = 0`).
- `F(Ω1, Ω2)`: joint spectral amplitude (JSA) in the low-gain limit. `λ_n`: its normalized Schmidt weights, Σλ_n = 1. `K = 1/Σλ_n²`: Schmidt number.
- `α`: SHG normalized efficiency of the OPA [1/W] (Ref A: 9.06 /W). `P`, `P_peak`, `P_avg`: CW / peak / average pump power. `g`: parametric coupling per length [W^{-1/2} m^{-1}].
- `A_p(z, T)`: pump complex envelope normalized so |A_p|² = power [W], phase `φ_p` (conventions §4).
- `b̂(z, T)`: signal annihilation operator envelope; `b̃(z, Ω)` its transform; `C`, `S`: input–output kernels, `b̃_out = C ã_in + S ã_in†`.
- `L_k`: fiber length of port k (D08). `δL_k`: residual length mismatch squeezed vs LO_k (D09). `L_D = T0²/|β2|`, `L_D′ = T0³/|β3|`: dispersion lengths.
- `n_g`: group index of the fiber. `c`: speed of light.
- `D_τ ≡ (1/L) dτ_g/dT`: thermal coefficient of group delay per length [s/(m·K)] (T08). `D_ψ ≡ (1/L) dψ/dT` [rad/(m·K)] (T08).
- `R_τ,mech ≡ dτ_g/dψ` for a given mechanism (thermal, PZT stretch) [s/rad].
- `S_t(f)`: one-sided PSD of the pulse-train timing jitter of the source [s²/Hz]. `σ_t`: RMS timing difference.
- `V_−`, `V_+`, `θ_off`: minimum and maximum quadrature variance seen by the LO mode and the angle offset of the minimum (`notes/conventions.md` §2, vacuum = 1/2).

## 2. 支配方程式

### 2.1 パルス励起の導波路 OPA（高利得、ポンプ非枯渇）

結晶の中で信号の群速度で動く座標（T = t − β1,s z）をとると、ポンプ包絡線と信号の演算子は次に従う（ポンプ非枯渇、縮退 type-0、QPM で Δk(0,0) = 0）。

- ポンプ: ∂_z A_p = −Δβ1 ∂_T A_p − (i β2,p/2) ∂_T² A_p（Agrawal の形。conventions §8 の ∂/∂T → −iΩ で、周波数領域では Ã_p(z,Ω) = Ã_p(0,Ω) exp[i(Δβ1 Ω + β2,p Ω²/2) z]）。
- 信号（周波数領域、厳密な分散を残す形）:

  ∂_z b̃(z, Ω) = i D_s(Ω) b̃(z, Ω) + g ∫ (dΩ′/2π) Ã_p(z, Ω + Ω′) b̃(z, Ω′)†

  D_s(Ω) = k_s(ω0 + Ω) − k_s(ω0) − β1,s Ω（GVD 以上の残り）。ポンプ側の Δk の残りは Ã_p の伝搬の中に入っている。
  （A_p(T) b̂†(T) の変換が (1/2π)∫Ã_p(Ω+Ω′) b̃(Ω′)† dΩ′ になることは conventions §8 の規約で確認した。）
- 規格化: CW・分散なし（A_p = √P e^{iφ_p}、D_s = 0）なら ∂_z b̂ = g√P e^{iφ_p} b̂†、すなわち b̂(L) = b̂ cosh(g√P L_c) + b̂† e^{iφ_p} sinh(g√P L_c)。
  conventions §3 の S†âS = â cosh r − â† e^{iφ_ξ} sinh r と比べて **r = g√P L_c、φ_ξ = φ_p + π、θ_sq = φ_p/2 + π/2**（conventions §4 と一致）。
  Ref A の R± = exp(±2√(αP)) と合わせると **g = √α / L_c**。
- 解は線形なので b̃_out = C ã + S ã†。数値では Ω を N 点に離散化し、2N×2N の線形 ODE を z について積分する（区間ごとに行列指数）。

この形は Wasilewski et al. 2006、Christ et al. 2013 と同じ構造（文献の式は本文で未照合、**UNVERIFIED**。ここでの符号と規格化はこのノートで conventions から導いたもの）。

### 2.2 超モード（Bloch–Messiah）

- Bloch–Messiah 分解（Braunstein 2005）: C = U diag(cosh r_n) V†、S = −U diag(e^{iφ_n} sinh r_n) V^T。
  出力モード b̂′_n = (U† b̂)_n、入力モード â′_n = (V† â)_n とおくと b̂′_n = â′_n cosh r_n − â′_n† e^{iφ_n} sinh r_n（各モードが conventions §3 の単一モードの squeezer。この式は代入で確認した）。
- 入力が真空なら出力の状態は U の列（**出力の超モード ψ_n**）と r_n、φ_n だけで決まる。数値では対称行列 M = ⟨b̃ b̃⟩ = C S^T の Takagi 分解
  M = U diag(−e^{iφ_n} sinh r_n cosh r_n) U^T で求める（V は不要）。
- 低利得の極限（r_n ≪ 1）では S ≈ g L_c × F で、超モードは JSA の Schmidt モードに一致し、r_n ∝ √λ_n。

低利得の JSA:

  F(Ω1, Ω2) ∝ Ã_p(0, Ω1 + Ω2) · sinc(Δk(Ω1,Ω2) L_c / 2) · e^{iΔk L_c/2}

1 次まで展開すると Δk ≈ Δβ1 (Ω1 + Ω2) + (β2,p/2)(Ω1+Ω2)² − (β2,s/2)(Ω1² + Ω2²)。
- 和の方向（Ω1 + Ω2）の幅: ポンプの帯域と、GVM による 2π/τ_GVM の小さい方。
- 差の方向（Ω1 = −Ω2）の幅: Δk = −β2,s Ω² から、最初の零点 Ω_PM = √(2π/(|β2,s| L_c))。これが CW の squeezing 帯域。
- ガウス近似（対角方向の 1/e 幅 σ_+、σ_−）では K = (1/2)(σ_−/σ_+ + σ_+/σ_−)（Law, Walmsley, Eberly 2000 の型の結果。係数の定義は本文で未照合、**UNVERIFIED**）。

高利得では超モードの形も重みも低利得の Schmidt モードから変わる（Christ et al. 2013 の主題。具体的な結論は本文未確認、**UNVERIFIED**）。したがって、量の見積もりには低利得の K を使ってよいが、実装では §2.1 の高利得の解を使う。

### 2.3 時間の切り出しで見た利得の低下（GVM が支配的なとき）[AE]

τ_p ≪ τ_GVM で GVD を無視すると、信号の時刻 T の部分はポンプと長さ L_w しか重ならない。

  r(T) ≈ (√α / L_c) ∫_0^{L_c} √(P_p(T − Δβ1 z)) dz

ガウスのポンプ（P_p(t) = P_peak exp(−t²/T0p²)）で τ_p ≪ τ_GVM なら

  r_max ≈ √(α P_peak) · √(2π) T0p / τ_GVM、τ_p ≫ τ_GVM なら r_max → √(α P_peak)（準 CW）。

squeezed 光の時間の広がりは τ_GVM 程度（ポンプが信号を追い越す時間）に伸びる。ポンプの GVD（後述 β2,p ≈ 388 fs²/mm）も 100 fs ではポンプ自身を数 mm で広げるので、この式は桁の目安にとどめる。

### 2.4 BHD で測る分散（多モード）

LO のモードの演算子は â_LO = Σ_n κ_n b̂′_n + (超モードに直交する成分、真空)。conventions §5 の x̂_θ を当てはめると

  V(θ) = Σ_n |κ_n|² · Var_n(θ − arg κ_n) + (1 − Σ_n |κ_n|²) · (1/2)

  Var_n(θ) = (1/2)[e^{−2r_n} cos²(θ − θ_sq,n) + e^{2r_n} sin²(θ − θ_sq,n)]（conventions §3）

各項は θ について a_n + b_n cos(2θ − c_n) の形なので、和も a + b cos(2θ − c) の形になる。したがって **V(θ) は 3 つの数 (V_−, V_+, θ_off) で完全に表せる**:
V_− = a − |b|、V_+ = a + |b|、θ_off = 最小になる角度 − θ_sq,0。これが時間ループへ渡す表の中身（§8）。
V_− V_+ > 1/4 になるのは、多モードによる「実効的な損失」と「モードごとの角度のばらつき」の両方の結果。

単一モード（ψ_0 だけ、r_0）のときは V = η_t Var_0 + (1 − η_t)/2、η_t = |κ_0|²（conventions §2 の純損失と同じ形。prompt §7-3 と整合）。

### 2.5 ガウスパルスの重なり（タイミングと分散の不一致）

LO と squeezed のモードが同じ形 u で、squeezed 側に Δt の遅れと φ2 の GDD の差があるとき（ũ_s = ũ e^{iΩΔt + iφ2Ω²/2}、conventions §8 の符号）:

  κ = (1/2π) ∫ |ũ(Ω)|² e^{iΩΔt + iφ2Ω²/2} dΩ = (1 − iχ)^{−1/2} exp[−Δt² (1 + iχ) / (4 T0² (1 + χ²))]、χ = φ2/(2T0²)

  η_t = |κ|² = (1 + χ²)^{−1/2} · exp[−Δt² / (2 T0² (1 + χ²))]

  arg κ = (1/2) arctan χ − χ Δt² / (4 T0² (1 + χ²))

- 特別な場合: φ2 = 0 で η_t = exp(−Δt²/(2T0²))。Δt = 0 で η_t = (1 + χ²)^{−1/2}、ファイバでは χ = δL/(2 L_D)。
- 幅が違う場合（T1, T2、φ2 = 0）: η_t = [2 T1 T2/(T1² + T2²)] exp[−Δt²/(T1² + T2²)]。
- **arg κ は squeezing 角のずれ（θ_off）として効く。** GDD の不一致があるとタイミングのずれが位相のずれに変わる（Δt² に比例）。ロックが測る位相（probe や CW 成分の位相）はこの項を含まないので、静的なオフセットとその揺らぎとして時間ループに入れる。
- β3 の不一致（φ3 = β3 δL）は解析式を使わず、ũ に e^{iφ3Ω³/6} を掛けて数値積分する（100 fs で L_D′ = 1.7 m [PH] なので 63 m の不一致では無視できない）。

### 2.6 タイミングジッタの平均

Δt がパルスごとに平均 μ、分散 σ_t² のガウスで揺らぐとき（φ2 = 0）、測定は η_t について線形なので平均だけが効く:

  ⟨η_t⟩ = (1 + σ_t²/T0²)^{−1/2} · exp[−μ² / (2 (T0² + σ_t²))]

ソースのジッタが、異なる時刻に出たパルスの間で相関を失う分（遅延差 τ_d）は conventions §7 の遅延干渉と同じ形:

  S_Δt(f) = 4 sin²(π f τ_d) S_t(f)、σ_t² = ∫ S_Δt(f) df

- D09 の配送（LO_k を同じパルスから分け、L_k に合わせて遅らせる）では τ_d = n_g δL_k / c ≈ 0 で、ソースのジッタはほぼ相殺される。
- Liu 型（D30、共通の参照光）でポートごとに τ_d = k T_rep（最大 31 T_rep = 310 ns @ 100 MHz）のときは、約 1/τ_d より上のジッタが効く。

### 2.7 熱によるタイミングドリフトと位相ロックによる共補正

区間の群遅延: τ_seg(t) = n_g L/c + D_τ L ΔT(t) + R_τ,st · ψ_st(t)。ψ_st はストレッチャが加える位相（conventions §4 の δψ = ξ_eo β δL）。

- 熱だけのとき: Δτ_th = D_τ L ΔT、Δψ_th = D_ψ L ΔT → R_τ,th = D_τ / D_ψ。
- ストレッチャ（伸び δL）: δτ = (n_g/c) δL（光弾性による群屈折率の変化は省略 [AE]）、δψ = ξ_eo (n ω0/c) δL → R_τ,st ≈ n_g / (ξ_eo n ω0)。
- 位相ロックがストレッチャで熱の位相を打ち消す（ψ_st = −Δψ_th）と、残るタイミングは
  **Δτ_res = Δτ_th (1 − R_τ,st / R_τ,th)**。
- TEC で温度そのものを戻す（Liu 型の遅い補正）なら同じ機構なので Δτ_res ≈ 0。
- ロックがレンジの端で 2π m 飛ぶとタイミングは 2π m R_τ,st だけ飛ぶ（1 フリンジ約 6.7 fs [AE]）。

## 3. 仮定と適用範囲

- ポンプ非枯渇（squeezed vacuum、T01 の E02）。χ(3) の効果（導波路内の SPM・XPM、光屈折）は入れない。kW 級のピークパワー（§5.2）ではこの仮定が崩れうる → T01 に確認を依頼。
- 結晶の分散はバルクの congruent LN（e 光、Zelmon et al. 1997）で代用。**Ref A の導波路（擬似単一モード PPLN）の導波路分散と組成（ドープ）は未確認**で、GVM・GVD は数十 % 変わりうる。
- 縮退 type-0、単一の空間モード、偏光は 1 つ。空間の多モード性は T01/T12（E07）。
- ファイバは線形（分散と群遅延だけ）。SPM で LO のスペクトルが広がる効果は T06 が φ_NL を出し、その後の ũ_LO を §2.4 の κ_n に入れる。
- ガウスの解析式（§2.5、§2.6）は検証と見積もり用。実装は任意の ũ で数値積分する。
- ジッタの平均（§2.6）は、測定がパルスごとの η_t に線形で、arg κ が揺らがない（φ2 ≈ 0）ときに厳密。φ2 ≠ 0 では表で平均する。
- R_τ の見積もりは survey の [AE] の熱係数を使った。値そのものは T08（熱）と T16（ξ_eo）待ち。

## 4. 作用する周波数帯

| 効果 | 帯域 | 備考 |
|---|---|---|
| 多モード squeezing（E11） | 光の帯域（差方向 約 6 THz、和方向 ≲ 0.07 THz） | 静的。オフライン段 |
| 分散の不一致（E21、E54） | 静的 ＋ 温度による β2 の変化（< 1 Hz） | β2 の温度係数は小さいと予想、T06 に依頼 |
| 熱によるタイミングドリフト（E26） | < 1 Hz | 時間ループでは遅い状態変数 |
| 音響の光路長変化によるタイミング | Hz – kHz | 1 rad あたり約 0.9 fs [AE]。位相雑音として重要でも、タイミングとしては無視できる |
| ソースのタイミングジッタ（E19） | DC – f_rep/2 | 速い成分は平均 ⟨η_t⟩（オフライン）、遅い成分は Δt(t)（時間ループ） |
| L_k、δL_k の静的な誤差（E51） | 静的 | オフラインの表の中心値 |

## 5. スケーリングと、自分たちのパラメータでの桁の見積もり

### 5.1 使った値

| 量 | 値 | 出典 |
|---|---|---|
| λ0, λ_p | 1545.32 nm, 772.66 nm | D07（Ref A の値で代用） |
| L_c | 45 mm | Ref A §3（将来の OPA の長さは未定、代用） |
| α | 9.06 /W | Ref A（CW） |
| LN の分散 | Zelmon e 光、バルク、21 °C | Zelmon et al. 1997（refractiveindex.info で係数を確認、2026-10-28 閲覧） |
| β2（ファイバ） | −21.7 ps²/km | **[PH]**、SMF 型の典型値、UNVERIFIED、T06 待ち |
| β3（ファイバ） | 0.127 ps³/km | **[PH]**、UNVERIFIED、T06 待ち |
| n_g（ファイバ） | 1.468 | survey [AE] |
| D_τ | 39 ps/(km·K) | survey E26 [AE]、T08 待ち |
| D_ψ | 44 rad/(m·K) | survey E22 [AE]、T08 待ち |
| ξ_eo, n | 0.78, 1.444 | conventions §4（T16 で確認）、n は [AE] |
| f_rep | 10 MHz – 1 GHz（主に 100 MHz） | D07 で Open、パラメータ |

計算は Python（numpy）で行った。スクリプトはコミットしていない（Phase 4 で検証テストとして書き直す）。

### 5.2 OPA（Ref A 型の 45 mm を代用）

- n_g(1545 nm) = 2.182、n_g(773 nm) = 2.274 → **GVM Δβ1 = 0.304 ps/mm、τ_GVM = 13.7 ps**（survey E11 の「約 0.3 ps/mm」[AE] を確認。ただしバルクの値）。
- GVD: β2,s = 101 fs²/mm、β2,p = 388 fs²/mm。CW の squeezing 帯域（低利得の最初の零点）Ω_PM/2π = **5.9 THz**（Kashiwazaki 2021 の 6 THz と同程度、Ref A の引用 [19]）。
- 低利得の Schmidt 数（厳密な Sellmeier の Δk、ガウスのポンプ、±12 THz の格子 2400–3200 点で SVD）:

| τ_p (pump) | K | 備考 |
|---|---|---|
| 100 fs | 約 320（±24 THz の格子では 約 640） | sinc の差方向の裾で格子に依存。≥ 10² とだけ言える |
| 1 ps | 約 300 | |
| 3 ps | 約 210 | |
| 10 ps | 約 195 | 格子について収束（195.0 → 196.2） |

  先頭のモードの重みは λ_0 ≈ 0.004–0.009（ほぼ等しい重みの対が続く）。**掃引範囲のすべてで強い多モード**。
- 利得の低下（§2.3、Ref A の CW の r = √(αP) = 2.41 を 640 mW で得るのに必要な量、[AE]）:

| τ_p | 係数 √(2π)T0p/τ_GVM | 必要な P_peak | P_avg @ 100 MHz |
|---|---|---|---|
| 100 fs | 0.011 | 約 5.3 kW | 約 56 mW |
| 1 ps | 0.11 | 約 53 W | 約 5.6 mW |
| 10 ps | ≥ 1（準 CW） | 約 0.64 W | 約 0.7 mW |

  この表は GVD・ポンプ枯渇・χ(3)・光屈折を無視しているので、桁の目安。fs 領域では導波路の GVM を 0 に近づける設計（Nehra et al. 2022 の分散設計した薄膜 LN が該当。Ref A の引用 [21]）か、短い結晶が必要になる。

### 5.3 時間モードの重なりの許容（LO と squeezed、または 2 ポート間）

η_t ≥ 0.99 の条件（ガウス、§2.5）:

| τ_p | T0 | \|Δt\| | ファイバ長に換算 | δL（GDD だけ、Δt = 0） | L_D | L_D′ |
|---|---|---|---|---|---|---|
| 100 fs | 60 fs | 8.5 fs | 1.7 µm | 4.7 cm | 0.17 m | 1.7 m |
| 1 ps | 0.60 ps | 85 fs | 17 µm | 4.7 m | 16.6 m | 1.7 km |
| 10 ps | 6.0 ps | 0.85 ps | 174 µm | 474 m | 1.7 km | 1.7×10³ km |

- **長さの不一致 δL_k は、タイミングとしては µm、分散としては cm – 100 m で効く。** 普通はタイミングを別の遅延（自由空間など）で合わせるので、2 つは別のパラメータとして持つべき（§9 の提案）。
- ジッタ: ⟨η_t⟩ ≥ 0.99 には σ_t ≤ 0.142 T0（同じく 8.5 fs / 85 fs / 0.85 ps）。光源の仕様は Open（D07）なので、これは**光源への要求**として出す。

### 5.4 ポート間の分散（指標 (ii)）

L_k = k c T_rep / n_g（D08）: f_rep = 1 GHz / 100 MHz / 20 MHz / 10 MHz で L_31 = 6.3 m / 63 m / 317 m / 633 m。

補償なしのポート 31 とポート 0 の重なり η_t（Δt = 0、β2 = −21.7 ps²/km [PH]）:

| L_31 − L_0 | φ2 | 100 fs | 1 ps | 10 ps | η_t ≥ 0.99 に要る τ_FWHM |
|---|---|---|---|---|---|
| 6.3 m（1 GHz） | −0.137 ps² | 0.053 | 0.98 | 1.000 | ≥ 1.2 ps |
| 63 m（100 MHz） | −1.37 ps² | 0.005 | 0.47 | 0.9998 | ≥ 3.7 ps |
| 633 m（10 MHz） | −13.7 ps² | 0.0005 | 0.053 | 0.98 | ≥ 11.6 ps |

ここでの τ は「使うモード」（BHD では LO のモード、GBS ではフィルタ後の単一モード）の幅。§5.2 のとおり squeezed 光そのものは THz 幅の多モードなので、GBS に入れる前にスペクトルを絞るか、ポートごとに分散を補償する必要がある（Liu et al. 2026 は 1.6 ns のポンプと MZI・0.4 nm フィルタで単一モードにし、波束の重なり > 99.7% を報告）。

### 5.5 熱によるタイミングドリフト

- 補償なし: D_τ L = 39 ps/(km·K) × 63 m = **2.5 ps/K**。η_t ≥ 0.99 に許される温度差（2 本が無相関の最悪値）: **3.5 mK（100 fs）、35 mK（1 ps）、0.35 K（10 ps）**。
  squeezed 側と LO_k 側の温度の相関係数 ρ（D09）があれば実効的な差は √(2(1 − ρ)) 倍程度（同じ長さのとき）。
- 位相ロックがあるとき（§2.7）: R_τ,th = 0.886 fs/rad、R_τ,st = 1.069 fs/rad、1/ω0 = 0.820 fs/rad → 残る割合 |1 − R_τ,st/R_τ,th| ≈ **0.21** [AE]。許容は約 5 倍に広がる（100 fs で約 17 mK）。
  この割合は dn/dT と dn_g/dT の差、ξ_eo に敏感で、符号も含めて T08・T16 の値が要る。
- 位相に直すと、100 fs でタイミングが問題になるのは補償なしで約 10 rad、ロックありで約 50 rad の熱の位相に相当する。63 m では 2800 rad/K [AE] なので、**fs パルスではタイミングの状態を時間ループで追う必要がある。** 10 ps では 0.35 K 以内なら無視できる。

### 5.6 Ref D

Ref D は CW（D02）なので本トピックの効果はない（η_t = 1、θ_off = 0）。ΔL による遅延干渉は T05。

### 5.7 スケーリングのまとめ

| 量 | ファイバ長 L | パワー | 温度 | パルス幅 T0 |
|---|---|---|---|---|
| タイミングの許容 | — | — | — | ∝ T0 |
| 熱ドリフト Δt | ∝ L ΔT | — | ∝ ΔT | — |
| 許される ΔT | ∝ 1/L | — | — | ∝ T0 |
| GDD の不一致の効果 χ | ∝ ΔL | SPM（T06）で LO が広がる: ∝ P_peak L | β2(T)（小、T06） | ∝ 1/T0² |
| TOD の効果 | ∝ ΔL | — | — | ∝ 1/T0³ |
| OPA の r（GVM 支配） | — | ∝ √P_peak | QPM の温度（T01） | ∝ T0（τ_p ≪ τ_GVM）、一定（τ_p ≫ τ_GVM） |
| Schmidt 数 K | — | 高利得で変わる | — | 弱い依存（表 5.2） |
| ジッタの効果（Liu 型） | 遅延 k T_rep で増える | — | — | ∝ 1/T0² |

## 6. CW の極限（prompt §5.1、D17）

- τ_p → ∞ で Ã_p(Ω1 + Ω2) → 2π √P δ(Ω1 + Ω2) e^{iφ_p}。§2.1 の式は Ω と −Ω の対ごとに閉じ、各 Ω で 2 モード squeezing（側帯波の描像）になる。
  r(Ω) = √(αP)·（位相整合の因子）、検出帯域（MHz）≪ Ω_PM（THz）では r(Ω) は一定で Ref A の R± に戻る。
- LO も CW なら κ はスペクトル上の「LO の周波数の位置」だけで決まり、η_t = 1、θ_off = 0。タイミングの状態は定義されない。
- 実装: `LASER_MODE = CW` のとき
  - オフライン: `opa_supermodes` は呼ばず、単一モード r = √(αP) を返す。`mode_overlap_table` は 1 点 (V_−, V_+, θ_off) = Ref A のモデルの値、θ_off = 0。
  - 時間ループ: `timing_state` と `mode_overlap_lookup` は作らない（または恒等写像）。
  - prompt §7-1 の数値（R′− = −10.15 dB、R+ = +20.56 dB）が**ビット単位で同じ**になることをテストにする（V5）。
- 途中の極限（τ_p ≫ τ_GVM、L_D ≫ 全ファイバ長）では、パルスのモデルが CW の値に連続的に近づくことを確かめる（V6）。

## 7. 検証テスト（実装が通すべきもの）

| # | テスト | 合格条件 |
|---|---|---|
| V1 | ガウスの重なり: 数値積分の κ を §2.5 の解析式と比べる。Δt ∈ [−3T0, 3T0]、χ ∈ [−5, 5] | \|κ_num − κ_ana\| < 1e−6。符号（Δt > 0 で squeezed が遅い、arg κ の符号）も一致 |
| V2 | OPA の解のシンプレクティック性: C C† − S S† = I、C S^T = S C^T | 最大要素の誤差 < 1e−8（N = 512 点の格子） |
| V3 | 低利得: g → 0 で S/(g L_c) が JSA F に、超モードが F の Schmidt モードに一致。ガウス近似の JSA で K が (1/2)(σ_−/σ_+ + σ_+/σ_−) に一致 | K の相対誤差 < 1e−3 |
| V4 | 多モードの BHD（§2.4）: u_LO = ψ_0 なら V_− = e^{−2r_0}/2。u_LO がすべての ψ_n に直交なら V = 1/2（prompt §7-4 の類似） | 相対誤差 < 1e−9 |
| V5 | CW の極限: LASER_MODE = CW で prompt §7-1 の数値が変わらない | 出力が pulsed ブロックなしの単一モードモデルとビット単位で一致 |
| V6 | 準 CW への連続性: τ_p = 10, 30, 100 τ_GVM（ポンプの GVD を 0）で r_0 / √(αP_peak) → 1 に単調に近づく | 100 τ_GVM で \|1 − r_0/√(αP_peak)\| < 1e−2 |
| V7 | 単一モードの不整合 = 損失: 単一の超モードで V = η_t V_sq + (1 − η_t)/2 | prompt §7-3 と同じ精度で一致 |
| V8 | ジッタの平均: Δt ~ N(μ, σ_t²) の Monte Carlo（10⁶ 標本）の ⟨η_t⟩ と §2.6 の式 | 統計誤差（3σ）以内 |
| V9 | 共通の分散の相殺: LO と squeezed が同じ区間の列を通る（同じ β2、β3、長さ）と η_t = 1、θ_off = 0（prompt §7-11 の類似） | \|1 − η_t\| < 1e−12 |
| V10 | 熱と位相ロックの共補正（§2.7）: 熱ドリフトだけ・理想的なロックで、残るタイミングが Δτ_th (1 − R_τ,st/R_τ,th) に一致。TEC で補正する設定では 0 | 相対誤差 < 1e−6 |
| V11 | 表の補間: 時間ループで Δt_k(t) から引いた (V_−, V_+, θ_off) を、同じ Δt でオフラインの直接計算と比べる | 表の格子の半分の点で相対誤差 < 1e−3 |

## 8. シミュレータへの要求

### 8.1 ブロックと表

| 名前 | 段 | 入力 | 出力 | 切り替え |
|---|---|---|---|---|
| `opa_supermodes` | O（Q と接する） | 結晶の分散（Sellmeier か n_g・β2 の表）、L_c、α、ポンプの包絡線（τ_p、P_peak、チャープ）、格子 | ψ_n(T)、r_n、φ_n、K、C、S | `OPA_MULTIMODE`（off → 単一モード r = √(αP_peak)）、`OPA_GVM`、`OPA_GVD` |
| `path_mode_transfer` | O（T06 の出力を使う） | 区間の列（L、β2、β3、T06 の SPM の結果）、LO の初期の ũ_LO | 各ポート k の BHD での ũ_LO,k、squeezed 側の伝達 e^{iφ(Ω)} | `FIBER_DISPERSION_IN_OVERLAP`、`LO_SPM_IN_OVERLAP` |
| `mode_overlap_table` | O | 上の 2 つ、Δt の格子、（必要なら）ΔT や φ2 の格子、ジッタの σ_t,fast | 各ポートの表 (V_−, V_+, θ_off)(Δt, …)、ポート間の η_t,kj(Δt) | `TIMING_JITTER_FAST` |
| `timing_state` | T | 区間ごとの ΔT(t)（T08）、ストレッチャの ψ_st(t)（T16）、ソースのジッタの遅い成分 S_t(f)、R_τ,th、R_τ,st | 区間ごとの τ_seg(t)、ポートごとの Δt_k(t) と Δt_kj(t) | `THERMAL_TIMING`、`TIMING_PHASE_COUPLING`、`TIMING_JITTER_SLOW` |
| `mode_overlap_lookup` | T | Δt_k(t)、表 | V_−(t)、V_+(t)、θ_off(t)（T04 の写像へ） | `PULSED_MODE`（CW では恒等） |
| 指標 (ii) の重なり | M | Δt_kj(t)、ポート間の表 | ポート間の波束の重なり η_t,k0(t) | — |

### 8.2 オフライン段から時間ループへの受け渡し

- 遅い変数は Δt_k だけ（必要なら φ2 の温度変化）。表は 1 次元（Δt）か 2 次元（Δt, φ2）で十分。§2.4 のとおり 3 つの数で V(θ) が完全に決まるので、表は小さい。
- ジッタは帯域で分ける: チャンクのナイキスト周波数より速い成分は σ_t,fast としてオフラインで ⟨·⟩ をとって表に入れ、遅い成分は `timing_state` が Δt_k(t) として生成する（二重に数えない）。
- θ_off(t) は squeezing 角に足す（conventions §10 の Δ_k = φ_LO,k − θ_sq,k に θ_off を加える）。ロックが測る位相には入らないので、ロックで補正されない。
- 表を更新する頻度はチャンクごとで足りる（熱のタイミングは < 1 Hz）。

### 8.3 State への要求（D11 への追加の提案）

- 区間ごとに位相 ψ_seg,a(t) に加えて **群遅延 τ_seg(t)** を持つ。雑音源ごとに R_τ,mech を持たせ、τ_seg を ψ_seg の同じ雑音から計算する（熱とストレッチャで係数が違うので、位相から一律に換算してはいけない）。
- 区間ごとに β2、β3 を持つ（T06 と共通）。
- 2 ポート組の設定（D13）で、ポート間の Δt_kj を出す。

## 9. DECISIONS・conventions との関係（調整役への提案）

1. **D09 の δL_k を 2 つに分ける提案**: 長さの不一致（分散の不一致 φ2 = β2 δL_k を決める）と、タイミングの合わせ残り δτ_k（遅延の微調整後の値）。同じ記号にすると µm と cm のスケールが混ざる（§5.3）。
2. **conventions §9 の区間に群遅延 τ_seg と β2、β3 を追加する提案**（§8.3）。§10 の指標 (ii) に、ポート間の波束の重なり（不可識別性）を並べて出すことの提案。
3. survey の更新の提案: E11 の GVM「約 0.3 ps/mm [AE, UNVERIFIED]」→ 0.304 ps/mm（Zelmon 1997 のバルク e 光から計算、導波路では未確認）。
   §5 の Ref A の書誌（UNVERIFIED）→ Optics Express 34, 7958 (2026)、DOI 10.1364/OE.585323（Crossref で確認。巻・ページは Crossref の値）。
4. conventions §6 の θ_rms の式は、パルスの多モードでは V_−、V_+ を §2.4 の表の値で置き換えて使う（式そのものは変更不要）。

## 10. ユーザの決定が必要な点（推測しない）

1. **将来の OPA の設計**: Ref A 型の 45 mm（τ_GVM ≈ 14 ps）を使うのか、GVM を設計した導波路（Nehra 2022 型）か、短い結晶か。100 fs – 1 ps を使うなら、§5.2 のとおり利得と単一モード性の両方で大きく不利。
2. **GBS に入れる前のスペクトルの絞り込み**（Liu 型のフィルタ）を前提にするか。その帯域が §5.4 の T0 を決め、損失（T12）とのトレードオフになる。
3. **ポートごとの分散補償**（DCF、チャープ FBG など）を候補に入れてよいか。100 MHz なら τ < 約 3.7 ps で必要になる。
4. **タイミングの微調整**（ポートごと・LO ごとの遅延ステージ）と、タイミングの能動的なロックを候補に入れるか。fs パルスでは位相ロックによる共補正だけでは足りない可能性がある（§5.5）。
5. **LO パルスの作り方**: 基本波のパルスをそのまま使うのか、先頭の超モードに合わせて整形するのか。§2.4 の κ_n を決める。
6. LO_k と squeezed 側のファイバの種類（PM と SMF の組み合わせ、D32）。種類が違えば β2 が違い、φ2 の不一致が L に比例して残る。

## 11. 文献

一次資料（書誌は Crossref で確認。本文を読んだものは明記）:

- K. Hirota, T. Kashiwazaki, et al., "Generation of 10-dB squeezed light from a broadband waveguide optical parametric amplifier with improved phase locking method," Optics Express **34**, 7958 (2026). DOI: 10.1364/OE.585323。（Ref A。本文は `refs/` の PDF で読んだ。書誌は Crossref）
- H.-L. Liu, H. Su, Y.-H. Deng, S.-Q. Gong, et al., "Gaussian boson sampling with 1,024 squeezed states in 8,176 modes," Nature **653** (2026)。DOI: 10.1038/s41586-026-10523-6。（本文を `refs/` の PDF で読んだ: ポンプ 1.6 ns、繰り返し 20 MHz（50 ns）、波束の重なり > 99.7%、識別不能性 97%（r = 1.8）、逆向き probe と TEC）
- W. Wasilewski, A. I. Lvovsky, K. Banaszek, C. Radzewicz, "Pulsed squeezed light: Simultaneous squeezing of multiple modes," Phys. Rev. A **73**, 063819 (2006). DOI: 10.1103/PhysRevA.73.063819。（書誌のみ確認、本文未読）
- A. I. Lvovsky, W. Wasilewski, K. Banaszek, "Decomposing a pulsed optical parametric amplifier into independent squeezers," J. Mod. Opt. **54**, 721–733 (2007). DOI: 10.1080/09500340600777805。（書誌のみ確認）
- S. L. Braunstein, "Squeezing as an irreducible resource," Phys. Rev. A **71**, 055801 (2005). DOI: 10.1103/PhysRevA.71.055801。（Bloch–Messiah。書誌のみ確認）
- G. Cariolaro, G. Pierobon, "Reexamination of Bloch-Messiah reduction," Phys. Rev. A **93**, 062115 (2016). DOI: 10.1103/PhysRevA.93.062115。（数値の分解の手順の参考。書誌のみ確認）
- A. Christ, B. Brecht, W. Mauerer, C. Silberhorn, "Theory of quantum frequency conversion and type-II parametric down-conversion in the high-gain regime," New J. Phys. **15**, 053038 (2013). DOI: 10.1088/1367-2630/15/5/053038。（高利得の超モード。書誌のみ確認）
- C. K. Law, I. A. Walmsley, J. H. Eberly, "Continuous frequency entanglement: effective finite Hilbert space and entropy control," Phys. Rev. Lett. **84**, 5304 (2000). DOI: 10.1103/PhysRevLett.84.5304。（Schmidt 数。書誌のみ確認）
- G. Patera, N. Treps, C. Fabre, G. J. de Valcárcel, "Quantum theory of synchronously pumped type I optical parametric oscillators: characterization of the squeezed supermodes," Eur. Phys. J. D **56**, 123–140 (2010; Crossref の issued は 2009). DOI: 10.1140/epjd/e2009-00299-9。（超モードの用語。書誌のみ確認）
- J. Roslund, R. M. de Araújo, S. Jiang, C. Fabre, N. Treps, "Wavelength-multiplexed quantum networks with ultrafast frequency combs," Nat. Photonics **8**, 109–112 (2014). DOI: 10.1038/nphoton.2013.340。（書誌のみ確認）
- M. E. Anderson, J. D. Bierlein, M. Beck, M. G. Raymer, "Quadrature squeezing with ultrashort pulses in nonlinear-optical waveguides," Opt. Lett. **20**, 620 (1995). DOI: 10.1364/OL.20.000620。（Ref A [15]。書誌のみ確認）
- D. K. Serkland, M. M. Fejer, R. L. Byer, Y. Yamamoto, "Squeezing in a quasi-phase-matched LiNbO3 waveguide," Opt. Lett. **20**, 1649 (1995). DOI: 10.1364/OL.20.001649。（Ref A [16]。書誌のみ確認）
- R. Nehra, R. Sekine, L. Ledezma, Q. Guo, R. M. Gray, A. Roy, A. Marandi, "Few-cycle vacuum squeezing in nanophotonics," Science **377**, 1333–1337 (2022). DOI: 10.1126/science.abo6213。（Ref A [21]。「25 THz で 4 dB」は Ref A の本文の記述。書誌のみ確認）
- T. Kashiwazaki et al., "Fabrication of low-loss quasi-single-mode PPLN waveguide and its application to a modularized broadband high-level squeezer," Appl. Phys. Lett. **119**, 251104 (2021). DOI: 10.1063/5.0063118。（Ref A [19]、「6 THz で 6 dB」は Ref A の本文の記述）
- T. Kashiwazaki et al., "Over-8-dB squeezed light generation by a broadband waveguide optical parametric amplifier toward fault-tolerant ultra-fast quantum computers," Appl. Phys. Lett. **122**, 234003 (2023). DOI: 10.1063/5.0144385。（Ref A [23]。書誌のみ確認）
- A. Shaked et al., "Lifting the bandwidth limit of optical homodyne measurement with broadband parametric amplification," Nat. Commun. **9**, 609 (2018). DOI: 10.1038/s41467-018-03083-5。（パルス・広帯域の測定の代替。書誌のみ確認）
- A. I. Lvovsky, M. G. Raymer, "Continuous-variable optical quantum-state tomography," Rev. Mod. Phys. **81**, 299–332 (2009). DOI: 10.1103/RevModPhys.81.299。（パルス BHD のモードの重なりの一般的な扱い。書誌のみ確認）
- D. E. Zelmon, D. L. Small, D. Jundt, "Infrared corrected Sellmeier coefficients for congruently grown lithium niobate and 5 mol.% magnesium oxide-doped lithium niobate," J. Opt. Soc. Am. B **14**, 3319–3322 (1997). DOI: 10.1364/JOSAB.14.003319。（e 光の係数 A = 2.9804, 0.5981, 8.9543、B = 0.02047, 0.0666, 416.08 µm² を https://refractiveindex.info/?shelf=main&book=LiNbO3&page=Zelmon-e から取得、2026-10-28 閲覧。原典の表との照合は未実施 → **UNVERIFIED**）
- O. Gayer, Z. Sacks, E. Galun, A. Arie, "Temperature and wavelength dependent refractive index equations for MgO-doped congruent and stoichiometric LiNbO3," Appl. Phys. B **91**, 343–348 (2008). DOI: 10.1007/s00340-008-2998-2。（温度依存の分散が必要になったときの候補。書誌のみ確認、使っていない）
- A. H. Hartog, A. J. Conduit, D. N. Payne, "Variation of pulse delay with stress and temperature in jacketed and unjacketed optical fibres," Opt. Quantum Electron. **11**, 265–273 (1979). DOI: 10.1007/BF00620112。（群遅延の熱係数・応力係数の一次資料の候補。T08 に渡す。書誌のみ確認、値は未読）

値の出典として **UNVERIFIED** のもの: ファイバの β2、β3（T06）、D_τ、D_ψ（T08）、ξ_eo（T16）、n = 1.444。§2.2 の K の式の係数、§2.3 の利得の低下の式（agent の導出）。
