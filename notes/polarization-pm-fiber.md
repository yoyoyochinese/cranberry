# T09 偏光と PM ファイバ（polarization-pm-fiber）

- Phase: 1（deep-dive）／ 日付: 2026-10-28
- 担当する効果: E36–E40（survey §2.5、D18-3）。E27 の偏光に関わる部分（915B の巻線部）。
- 入力: `DECISIONS.md`（D02、D07–D13、D18、D24、D27、D29–D33）、`notes/conventions.md`、`notes/survey.md`、
  `refs/manuals/914-915B-Specifications.pdf`、Liu et al. 2026（`refs/`）、Ref C（`refs/`）。
- 依存: 温度の時間変化（PSD・ドリフト速度）は T08 の担当。ここでは係数（rad/(m·K)）までを出し、温度の揺らぎは D27 のパラメータ（ドリフト速度 [K/s]）で代用する。
  probe の混入光子数（E43）と時刻ずれ（E42）は T10 の担当。ここでは T10 が使う PER の収支と、偏光による位相の食い違いの式を渡す。
- 値の印は survey §0.3 と同じ: **[manual]**、**[Ref X]**、**[AE]**（この担当の概算、未検証）、**UNVERIFIED**（原典を開いて確認していない）。

## 0. 結論の要約

1. **D11 の概算はおおむね正しい（向きは負）。** 一次資料（Ritari et al. 2004、標準 PANDA、1550 nm、2 方法で測定）では
   K_T = ∂β_b/∂T = −2.06 / −1.8 rad/(K·m)（L_b = 4.2 mm）、規格化すると K_T/β_b = −1.37×10⁻³ /K。
   B = λ/L_b = 3.7×10⁻⁴ に対して dB/dT = −5.1×10⁻⁷ /K。D11 の「6×10⁻⁷ /K、2.4 rad/(m·K)」は大きさで +20% 程度の過大で、範囲内。
   915B の PM15-U25D（Corning PANDA PM 1550、L_b = 3.0–5.0 mm [Corning datasheet]）に規格化した値を当てはめると
   **|K_T| = 1.7–2.9 rad/(m·K)、名目 2.05 rad/(m·K)。60 m で 123 rad/K（範囲 103–172 rad/K）**。符号は負（温度が上がると slow 軸と fast 軸の位相差が減る）。
   ただし別の文献には 3–6 倍小さい値もある（§2.3）。掃引の範囲を広めにとる。
2. **PM ファイバに直交偏光の probe（Liu 型）を通すと、squeezed 光の位相が K_T·L·ΔT だけずれたままロックされる。** 63 m で 10 mrad に抑えるには、
   ファイバ全体の平均温度を **約 80 µK** の範囲に保つ必要がある（630 m なら 8 µK）。
3. **同じ偏光で逆向きに通す probe は、可逆な（非相反の素子がない）光路では、偏光がどう変わっても squeezed 光と同じ位相を測る**（Jones 行列の転置の関係、§1.6）。
   PM でも SMF でもそうなる。残る誤差は、(a) PER によるクロストーク光の DGD（PM で 1.36 ps/m）がパルス幅を超えるときの CW probe とパルスの違い（≈ Σ κ_iκ_j、10⁻³–10⁻² rad）、
   (b) SMF の PMD によるパルスの平均化（100 fs・63 m で ≤ 10⁻² rad [AE]）、(c) Faraday 効果などの非相反の項（≤ 2×10⁻³ rad [AE]）。
   一方で、偏光で probe を squeezed モードから分けられない（混入の抑制は T10）。
4. **SMF + 直交偏光の probe の誤差は 2·arg(a)**（a は偏光の回転 U の対角要素）。複屈折は PM より約 3 桁小さい（巻いた SMF の曲げ複屈折 B ≈ 2×10⁻⁷、R = 5 cm）が、
   偏光が保たれないので誤差が予測できず、|a| → 0 のときは probe の信号そのものが消える。
5. **Ref D:** PDB450C は偏光を区別しないので、両腕の fast 軸の成分どうしの縞が slow 軸の縞に重なり、ロック点が **δφ ≈ ε·sin δ_b(t)** だけ動く。
   ε は 2 つのクロストークの振幅の積（PER 20–25 dB なら ε ≈ 3×10⁻³–10⁻²）、δ_b は 2 腕の複屈折位相の差。915B の腕だけ約 8.4 m 長いとすると、
   δ_b は約 17 rad/K で回る（0.37 K で 1 周）。
6. **915B の巻線部:** 38 巻 × 220 mm/巻 ≈ 8.4 m [manual + 製品紹介ページ、UNVERIFIED]。曲げ半径は資料になく、筐体の寸法から 15–20 mm と推定 [AE]。
   曲げの複屈折による PER の低下は −31〜−36 dB 程度（向きがランダムな場合）[AE]。曲げ損失は R に指数的に効くので、0.01–2 dB の幅で不明（§3.4）。
   D10 により将来は squeezed 光の経路にストレッチャが入るので、R の実測かメーカーへの問い合わせが要る。

## 1. 支配方程式

### 1.1 記号（symbols, English）

- `λ`: vacuum wavelength [m]; Ref D: λ = 1548.8 nm (D02); config F placeholder: 1545.32 nm (D07). `k0 = 2π/λ`.
- `n_s`, `n_f`: effective (phase) index of the slow / fast eigenmode of PM fiber. `B = n_s − n_f > 0`: modal (phase) birefringence [–].
- `L_b = λ/B`: beat length [m]. `β_b = 2π/L_b = k0·B`: birefringent propagation constant [rad/m].
- `Δψ_b = ψ_slow − ψ_fast = β_b·L`: birefringent phase of a segment of length `L` [rad] (> 0 by conventions §4, since ψ = n k0 L).
- `K_T = ∂β_b/∂T`: temperature coefficient of birefringent phase per unit length [rad/(m·K)] (definition of Ritari et al. 2004). `k_T = K_T/β_b = (1/B)·dB/dT` [1/K] (length-normalized, relative).
- `B_g = B − λ·dB/dλ`: group birefringence. `τ_DGD = B_g·L/c`: differential group delay (DGD) between axes [s]. `b_g = B_g/c` [s/m] (Ritari: 1.36 ps/m).
- `T`: fiber temperature [K]; `ΔT(t)`: deviation from the value at lock acquisition. `T_soft`: effective softening temperature of the stress-applying parts.
- `θ_k`: angle between the slow axes of two fibers at joint `k` (connector, splice, component port) [rad].
- `κ_k`: amplitude cross-coupling at joint/point `k` (slow → fast), `κ_k ≈ sin θ_k`; `PER_k = −10·log10(κ_k²/(1−κ_k²)) = −10·log10(tan² θ_k)` [dB].
- `h`: distributed polarization-coupling (h-) parameter of PM fiber [1/m]; power in the orthogonal axis after length L: `P_x/P = tanh(h·L)`.
- `J`: 2×2 Jones matrix of a segment in the fixed lab transverse basis (x, y); `J = e^{iψ}·U`, `U ∈ SU(2)`, `U = [[a, −b*], [b, a*]]`, `|a|² + |b|² = 1`.
- `e_in`, `g`: unit Jones vectors of the polarizer at the source side and at the GBS / far side (real for linear polarizers). `⊥`: the orthogonal state.
- `Θ`: angle between two polarization states on the Poincaré sphere [rad]. `η_pol = |e_1^† e_2|² = cos²(Θ/2)`: power overlap.
- `V`: fringe visibility. `ε`: relative amplitude of the parasitic (fast-axis) fringe in the Ref D MZI. `δ_b`: phase of the parasitic fringe relative to the main one.
- `R`: bend radius [m]. `r`: cladding radius (62.5 µm). `B_bend`: bend-induced birefringence. `δ0 = B_bend/(2B)`: rotation of the principal axes in a bent PM fiber.
- `n`, `p11`, `p12`, `ν`: refractive index, photoelastic (Pockels) coefficients, Poisson ratio of silica.
- `D_PMD`: PMD coefficient of SMF [s/√m]; mean DGD `⟨τ_DGD⟩ = D_PMD·√L`.
- `τ_p`: pulse FWHM (intensity), `σ_t = τ_p/2.355`: rms intensity width (Gaussian assumed).
- `PDL_dB = 10·log10(T_max/T_min)`: polarization-dependent loss.
- `V_F`: Verdet constant [rad/(T·m)]; `B_∥`: magnetic flux density along the fiber [T].
- `L_k`: fiber length of port k (D08). `L`: generic segment length. `ΔL_arm`: PM length difference between the two Ref D arms.

### 1.2 PM ファイバの軸ごとの位相と温度係数（E37）

区間の 2 軸の位相（conventions §4・§9 の ψ の定義に従う）:

  ψ_seg,slow = n_s·k0·L,  ψ_seg,fast = n_f·k0·L = ψ_seg,slow − Δψ_b,  Δψ_b = β_b·L = k0·B·L

温度の変化に対して:

  dΔψ_b/dT = L·K_T + β_b·L·α_L ≈ L·K_T,  K_T = k0·dB/dT = k_T·β_b = k_T·2π/L_b

（α_L ≈ 5×10⁻⁷ /K の熱膨張の項は β_b·α_L ≈ 8×10⁻⁴ rad/(m·K) で、|K_T| ≈ 2 rad/(m·K) に比べて無視できる。Ritari の K_T は測った位相の変化なので両方を含む。）

PANDA の複屈折は応力付与部（B ドープのシリカ）とクラッドの熱膨張係数の差で生じる熱応力によるもので、応力は「凍結温度（軟化点）から室温までの冷却量」に比例する。
このモデル（Rashleigh 1983 など、UNVERIFIED）では B ∝ (T_soft − T) なので

  k_T = (1/B)·dB/dT ≈ −1/(T_soft − T)

Ritari の k_T = −1.37×10⁻³ /K は T_soft − T ≈ 730 K（T_soft ≈ 1020 K ≈ 750 °C）に相当し、このモデルと矛盾しない。
Ritari は「応力付与部が同じ種類なら k_T は同じ値になる」と報告している（PANDA PCF 2 種と標準 PANDA で −1.38〜−1.65×10⁻³ /K）。
**シミュレーションでは k_T を基本のパラメータにし、K_T = k_T·2π/L_b で計算する**（L_b の製造ばらつき 3–5 mm をそのまま反映できる）。

軸ごとの温度係数（D11「温度は区間で共通、応答は軸で異なる」）:

  dψ_seg,slow/dT = L·C_T（T08 の係数、survey E22 では約 44 rad/(m·K) [AE]）
  dψ_seg,fast/dT = L·(C_T − K_T) = L·(C_T + |K_T|)

T08 の C_T がどちらの軸の値かは文献で明記されないことが多い。差は |K_T|/C_T ≈ 5% なので、C_T は slow 軸の値とみなし、fast 軸に K_T の補正をかける。

群遅延の差（DGD）と、その温度係数 [Ritari Table 1、2]:

  τ_DGD = b_g·L, b_g = 1.36 ps/m；dτ_DGD/dT = K_τ·L, K_τ = −1.78 fs/(K·m)

### 1.3 PER とクロストーク（E36）

**接合 1 か所**（PM コネクタ、融着、PM 部品のポート）: 軸のずれ θ_k で、slow 軸の振幅の κ_k = sin θ_k が fast 軸に移る。

  PER_k = −10·log10(tan² θ_k)

θ = 1° → 35.2 dB、2° → 29.1 dB、3° → 25.6 dB、5.7° → 20.0 dB。応力（コネクタのフェルールの接着、締め付け）による部分は角度に換算して同じ式で扱う。

**分布したクロストーク**（ファイバ自体、h パラメータ）:

  P_x(L)/P_total = tanh(h·L) ≈ h·L（h·L ≪ 1）

結合点ごとの位相がランダムなので、出力では fast 軸の成分は互いにほぼインコヒーレントに足し合わされる。

**連結したときの PER**: 1 か所ずつの漏れが小さく位相がランダムなら、パワーで足す: 10^(−PER_tot/10) ≈ Σ_k 10^(−PER_k/10) + h·L。
最悪（全部の振幅が同位相）なら振幅で足す: κ_tot = Σ κ_k。シミュレーションでは位相を区間の Δψ_b(t) から計算する（§1.5）ので、どちらも仮定しない。

**軸ずれが squeezed 光に与える影響**:
- 偏光子・PBS・GBS の入力が slow 軸だけを通すなら、fast 軸に移った成分は**損失**: η_PER = 1 − κ_tot²（conventions §2 の純損失）。
- 移った成分が再び slow 軸に戻る（2 回目の結合）と、主成分と干渉して**位相とパワーのリップル**になる（§1.5）。

### 1.4 偏光の重なりと可視度（E38、E39）

2 つのビーム（偏光 e_1、e_2）の干渉の可視度、BHD の LO と squeezed 光の偏光の重なり:

  V = 2√(P_1 P_2)/(P_1 + P_2)·|e_1^† e_2|,  η_pol = |e_1^† e_2|² = cos²(Θ/2)

BHD では η_pol はモード整合の効率の一部で、conventions §2 の純損失として squeezed 光に入る（T12 のモード整合の収支に入れる）。
1% の損失に相当するのは Θ = 0.20 rad（11.5°）、0.1% なら Θ = 0.063 rad。

**偏光子の前の偏光が揺らぐとき**（Ref D: Cobrite → 減衰器 → Optoquest PC → Opneti 偏光子）: 偏光子が偏光のずれをパワーの変化に変える。

  P_out(t) = P_in·[cos²(Θ(t)/2) + 10^(−ER_pol/10)·sin²(Θ(t)/2)]

Θ は偏光子の通過軸からの角度（Poincaré 球上）。設定時 Θ = 0 なら、δP/P ≈ −Θ²/4（二次）。偏光子の後のパワーのドリフトは DC side-fringe lock の傾きとロック点に効く（T15、D06）。

**PDL** の部品を通るとき: 透過率は偏光で T_min–T_max の間を動く。偏光がドリフトすると、T_mean·(1 ± (T_max−T_min)/(T_max+T_min)) の範囲で損失が揺らぐ。

### 1.5 2 つの結合点の干渉（PM で同じ偏光の probe、Ref D の寄生的な縞）

slow 軸の光が点 i で κ_i だけ fast 軸に移り、点 j（i より下流）で κ_j だけ slow 軸に戻ると、出力の slow 軸の振幅は

  A_slow,out = A·e^{iψ_slow}·[√(1−κ_i²)√(1−κ_j²) + κ_i κ_j·e^{−iΔψ_b,ij}]，Δψ_b,ij = β_b·ℓ_ij（ℓ_ij は i–j 間の長さ）

（符号は結合点の軸の回転の向きで変わる。下の位相は符号まで含めて κ_i κ_j の積に吸収する。）
これによる位相のずれ: δφ_ij ≈ −κ_i κ_j·sin(Δψ_b,ij(t))。Δψ_b,ij は温度で回る（ℓ_ij·K_T rad/K）。

**CW とパルスの違い**: 戻ってきた成分は主成分より τ_ij = b_g·ℓ_ij だけ遅れる。
- CW の probe: 遅れに関係なく干渉し、上の δφ_ij を受ける。
- パルスの squeezed 光（幅 τ_p）: τ_ij ≫ τ_p なら主パルスとは重ならず、別の時間のサテライト（κ_iκ_j の振幅、損失として κ_i²κ_j²）になる。位相のずれは受けない。
  b_g = 1.36 ps/m なので、τ_p = 100 fs なら ℓ_ij > 7 cm、10 ps なら ℓ_ij > 7 m でこの状況になる。
- したがって、**同じ偏光の CW probe とパルスの squeezed 光の位相差は** Σ_{i<j} κ_iκ_j·sin(Δψ_b,ij(t))·[1 − O(τ_ij)]、O(τ) はパルスの自己相関の重なり（τ = 0 で 1）。

### 1.6 逆向きに伝搬する probe と相反性（E41 の偏光の部分、D30・D32 の判断の核心）

**相反性**: 磁気光学効果（Faraday）のない光路では、固定した実験室の横方向の基底 (x, y) で、逆向きの Jones 行列は順方向の転置 J_back = J^T になる
（Martinelli 1989 の Faraday ミラーの議論などで使われる性質、UNVERIFIED。ここでは Lorentz の相反定理から出る結果として使い、§5 のテスト 2 で数値的に確かめる）。

光路: 光源側の偏光子 e_in →（共有区間 J = e^{iψ}U）→ GBS 側の偏光子 g（どちらも直線偏光、実ベクトル）。

- squeezed 光（順方向、e_in から入り g に射影）: t_s = g^T J e_in = e^{iψ}·(g^T U e_in)
- probe（逆方向、g' から入り光源側で r に射影）: t_p = r^T J^T g' = e^{iψ}·(g'^T U r)

**同じ偏光**（g' = g、r = e_in、つまり probe を GBS 側の同じ偏光子から入れ、光源側の同じ偏光子で取り出す）:

  t_p = e^{iψ}·(g^T U e_in) = t_s　→　**arg t_p = arg t_s（U によらない）**

偏光の回転 U が時間とともにどう変わっても、probe が測る位相は squeezed 光の（GBS の偏光子に射影した）位相と一致する。|t_p|² = |t_s|² なので、probe のパワーの低下は squeezed 光の偏光損失そのものを表す。

**直交偏光**（g' = g⊥、r = e_in⊥。Liu et al. 2026 の方式）:
e_in = x、g = x の基底で書くと t_s = e^{iψ}·a、t_p = e^{iψ}·(y^T U y) = e^{iψ}·a*。したがって

  **Θ_err ≡ arg t_s − arg t_p = 2·arg(a)**

- PM（軸を合わせた）: U = diag(e^{iΔψ_b/2}, e^{−iΔψ_b/2}) なので 2·arg(a) = Δψ_b。D11 の結果（2 軸の位相差がそのまま誤差になる）に一致する。
  ロックは probe の位相を固定するので、squeezed 光の位相は Δψ_b(t) − Δψ_b(0) = K_T·L·ΔT(t) だけずれる。
- SMF: U は弱いランダムな複屈折の積で、a(t) の偏角は環境で動く。|a|² も 0 から 1 まで動き、|a| → 0 では probe の信号がなくなり（同時に squeezed 光も g を通らない）、誤差が定義されない。

**この節が成り立つ条件**: (1) 共有区間に非相反の素子（アイソレータ、サーキュレータの内部、Faraday 回転子）がない。probe の出し入れはその外側で行う。
(2) 準静的: 雑音の周波数が通過時間の逆数 c/(n_g·L) より十分低い（高い周波数の逆向き伝搬の相関の崩れは T10、E42）。
(3) probe と squeezed 光の周波数が同じ（周波数をずらすと U(ω) と ψ(ω) の分散が入る。D32 のパラメータ。T10）。
(4) 単一の周波数（CW）の probe と広帯域のパルスの違いは §1.5 と §1.7 で別に扱う。

### 1.7 SMF の偏光のドリフト、PMD（E38、E40）

**曲げの複屈折**（Ulrich, Rashleigh, Eickhoff 1980、UNVERIFIED）:

  B_bend = (1/4)·n³·(p11 − p12)·(1 + ν)·(r/R)² ≈ 0.136·(r/R)²

（n = 1.46、p11 = 0.121、p12 = 0.270、ν = 0.17 で係数 0.136。この係数の値も UNVERIFIED。）fast 軸は曲げの面内。
半径 R = 5 cm のスプールに巻いた SMF では B_bend ≈ 2.1×10⁻⁷（L_b ≈ 7.3 m）。巻いた SMF は「弱い PM ファイバ」として振る舞う。

**偏光のドリフトのモデル**: 区間 j の Jones 行列 U_j(t) = exp(−i·(β⃗_j(t)·σ⃗)·ℓ_j/2)、β⃗_j は 3 次元の複屈折ベクトル（Stokes 空間、Gordon & Kogelnik 2000 の記法）。
環境（温度・振動・曲げ）で β⃗_j(t) が変わる。一次資料の普遍的な値はないので、**出力の偏光の Poincaré 球上の拡散係数 D_pol [rad²/s]** と、温度に対する回転の係数をパラメータにする（§3.6）。

**PMD**: 平均 DGD ⟨τ_DGD⟩ = D_PMD·√L（長いファイバ、ランダム結合の極限）。短いファイバ（L が結合長より短い）では τ_DGD ≈ B_g·L/c に近づく。
パルスの squeezed 光が DGD τ_DGD の 2 つの主偏光状態（PSP）に分かれると、固定した偏光への射影で時間モードがずれる。
ガウスパルスで、振幅が 2 つの PSP に等分された最悪の場合のモード整合の効率 [AE]:

  η_PMD,worst = [(1 + exp(−τ_DGD²/(8σ_t²)))/2]²

（中心周波数での偏光は補正済みとする。CW の probe と、パルスで平均した位相の差は τ_DGD·σ_ω の 2 次、≈ (τ_DGD/(2σ_t))²/8 程度 [AE]。）

### 1.8 非相反の効果（参考）

Faraday 効果: ファイバに沿った磁場 B_∥ で円偏光の位相が ±V_F·B_∥·L だけずれ、向きを逆にすると符号が変わる。§1.6 の相反性を破る。
V_F（シリカ、1550 nm）≈ 0.6 rad/(T·m)（UNVERIFIED）、地磁気 50 µT、L = 63 m がすべて平行という最悪の場合で V_F·B_∥·L ≈ 2×10⁻³ rad [AE]。
スプールに巻けば向きが打ち消しあう。PM ファイバでは大きな直線複屈折が Faraday 回転を抑える。既定ではオフでよい。
（Sagnac 効果（地球の自転）も非相反だが、63 m を半径 5 cm に巻いて 2×10⁻⁶ rad 程度 [AE] で無視できる。）

## 2. PM の複屈折の温度係数の確認（D11 の検証）

### 2.1 一次資料

**Ritari et al., Opt. Express 12, 5931 (2004)**（全文を読んだ）。1550 nm、標準 PANDA（比較用）と PANDA PCF。温度を制御した水槽、出力偏光の回転と POFDR の 2 方法。

| 量 | 標準 PANDA | PANDA5 (PCF) | PANDA9 (PCF) | 出典 |
|---|---|---|---|---|
| L_b [mm] | 4.2 | 16.6 / 18.3 | 14.4 | Table 1 |
| 位相の遅れ β [ps/m] | 1.23 | 0.309 | 0.359 | Table 1 |
| 群遅延の差 b_g [ps/m] | 1.36 | 0.497 | 0.392 | Table 1 |
| K_T [rad/(K·m)]（S_out / POFDR） | −2.06 / −1.8 | −0.624 / −0.6 | −0.601 / −0.6 | Table 2 |
| K_T/β_rad [1/K] | −1.37×10⁻³ | −1.65×10⁻³ | −1.38×10⁻³ | Table 2 |
| K_τ [fs/(K·m)] | −1.78 | −0.731 | −0.564 | Table 2 |

本文: 「標準 PANDA では、応力付与部による複屈折が温度とともに強く減る（応力のない線引き時の状態に近づくため）」。K_T は 1500–1600 nm で波長とともに |K_T| が増える（Fig. 3 では約 −2.0 → −2.4 rad/(K·m)）。

### 2.2 D11 との比較

| | D11（agent の概算） | Ritari（L_b = 4.2 mm） | PM15-U25D に規格化（L_b = 3.0–5.0 mm） |
|---|---|---|---|
| B | 4×10⁻⁴ | 3.7×10⁻⁴ | 3.1–5.2×10⁻⁴ |
| dB/dT [1/K] | 6×10⁻⁷（符号なし） | −5.1×10⁻⁷ | −(4.2–7.1)×10⁻⁷ |
| K_T [rad/(m·K)] | 2.4 | −2.05（−1.8〜−2.06） | −(1.72–2.87) |
| 60 m [rad/K] | 約 150 | −123 | −(103–172) |
| 63 m（L_31、f_rep = 100 MHz） | — | −129 | −(108–181) |
| 8.4 m（915B の巻線部） | — | −17 | −(14–24) |

**判定: D11 の値は大きさの範囲内（Ritari の値より約 20% 大きい）。符号は負。** D11 の結論（直交偏光の probe は PM ファイバでは約 2 rad/(m·K) の差動ドリフトを検出できない）は変わらない。
PM15-U25D の L_b は資料に範囲でしか書かれていない（3.0–5.0 mm）ので、名目値は Ritari の 4.2 mm を使い、そう明記する。

### 2.3 食い違う文献値

- Dvořák et al., Radioengineering 18, 639 (2009) は、別の文献（その論文の [1]）の値として、**633 nm**、L_b = 2 mm の PANDA で 0.757 rad/(m·°C)、
  bow-tie（L_b = 1 mm）で 0.735 rad/(m·°C) を引いている。k_T に直すと PANDA で 2.4×10⁻⁴ /K で、Ritari より約 6 倍小さい。
  可視域用の別のファイバで、元の測定は確認できていない（UNVERIFIED）。
- Zhang & Lit, Appl. Opt. 32, 2213 (1993) は PANDA・bow-tie・楕円コアの温度・歪みの感度を測っている（書誌は確認、本文は取得できず数値は UNVERIFIED）。
  検索結果の要約に「PANDA で 7.57 rad/(m·°C)」とあったが、波長と条件が不明なので使わない。
- PANDA の Sagnac ループ温度センサーの報告（約 1.4 nm/°C、UNVERIFIED）は、dλ/dT = λ·k_T·(β/B_g) から k_T ≈ −0.9×10⁻³ /K に相当し、Ritari と同じ桁。

**シミュレーションでの扱い**: k_T の名目値 −1.37×10⁻³ /K（Ritari）、掃引の範囲 −(0.25–1.6)×10⁻³ /K。
下限は Dvořák の引用値、上限は Ritari の PCF の値に余裕を足したもの。ユーザのファイバで実測するのが最も確実だが、D04 のブラインド検証（Ref D）とは独立の量なので、実測してよいかはユーザが決める（§7）。

## 3. 仮定と適用範囲、周波数帯、スケーリング

### 3.1 仮定と適用範囲

- 単一モード、弱導波。2 軸の固有モードの間の結合は摂動（κ ≪ 1、h·L ≪ 1）。
- K_T は 20–40 °C の範囲で一定（Ritari の水槽の範囲。線形からのずれは文献で確認していない）。
- 相反性（§1.6）は、共有区間に磁気光学の素子がなく、外部の磁場が小さいときに成り立つ。逆向きの伝搬による時刻のずれは準静的な近似の外（T10）。
- 曲げの複屈折の式（§1.7）は R ≫ r のときの弾性論の結果。巻くときの張力による複屈折（tension-coil）と、側圧・接着による複屈折は含まない（値がない、UNVERIFIED）。
- SMF の偏光のドリフトには普遍的な一次資料の値がない。拡散係数と温度係数はパラメータ（§3.6）。
- PMD の式は Maxwell 分布の平均。63 m 程度の短いファイバでは結合長に近く、平均値の式は目安。

### 3.2 周波数帯

| 効果 | 周波数帯 | 備考 |
|---|---|---|
| PM の複屈折の温度ドリフト（E37） | < 0.1 Hz（熱の時定数） | 温度の PSD に K_T·L をかける。温度の PSD は T08 |
| PER（E36） | 静的 ＋ 温度で回る位相 | 2 点の干渉の位相が K_T·ℓ_ij で回るので、リップルの周期は温度で 2π/(|K_T|·ℓ_ij) |
| SMF の偏光のドリフト（E38） | mHz–Hz（熱）、Hz–kHz（振動・音響で曲げとねじれが変わる） | D_pol とそのスペクトルはパラメータ |
| 偏光子の前の偏光の揺らぎ → パワー（E39） | mHz–Hz | 二次の効果（設定点が Θ = 0 のとき） |
| PDL、PMD（E40） | 静的 ＋ 偏光のドリフトと同じ帯域 | |

### 3.3 スケーリング

| 量 | ファイバ長 L | パワー | 温度 | パルス幅 τ_p | 波長 |
|---|---|---|---|---|---|
| 直交 probe の誤差（PM） K_T·L·ΔT | ∝ L | なし | ∝ ΔT | なし | K_T ∝ 1/λ（B 一定） |
| 分布クロストーク h·L | ∝ L | なし | 弱い | なし | — |
| 2 点干渉の位相 δφ_ij | κ_iκ_j（L によらない）、点の数 N で ∝ N² の組 | なし | 周期 2π/(K_T ℓ_ij) | τ_ij ≫ τ_p でパルスには効かない | — |
| PMD の η 低下 | ∝ L（τ_DGD² ∝ L） | なし | 弱い | ∝ τ_p⁻² | — |
| 曲げの PER 低下 | ∝ 巻き数、∝ R⁻⁴ | なし | 区間の間の位相が 1 巻あたり約 0.45 rad/K で回り、PER が温度で揺らぐ | なし | B_bend は λ にほぼよらない |
| 偏光子の後のパワー変化 | 前の SMF の L、環境 | ∝ P | Θ(T) | なし | — |

複屈折の効果はパワーによらない（非線形の偏光回転は T06 の SPM/XPM の範囲。PM ファイバでは軸が固定されるので無視できる [AE]）。

### 3.4 915B の巻線部（E27 の偏光の部分）

`refs/manuals/914-915B-Specifications.pdf` [manual] の値: 915B（serial 8669-1、8669-2）、巻き数 38、ファイバ PM15-U25D、
フル電圧での伸び（光路長）16 µm / 17 µm、ピグテールは FC/APC で **slow 軸にキー合わせ**。914 の入力 ±5 V。

製品紹介（photonics.com の 915B のページ、検索結果の抜粋、2026-10-28 閲覧、ページ本体は 403 で読めず UNVERIFIED）:
「楕円形の枠に弾性的に取り付けた圧電板でファイバを伸ばす。伸ばすのは直線部だけで、PM ファイバの曲げによるクロスカップリングを防ぐ」「1 巻あたり 220 mm」「9 m のファイバで最大 16 µm」「筐体 4.0 × 2.0 × 1.75 in」。

- 巻線部のファイバ長: 38 × 0.22 m = **8.4 m** [AE、manual と紹介ページの積]。9 m の記載と合う。Ref D の ΔL（D02）の下限の目安になる（T15・T16 に渡す）。
- 曲げ半径 R: 資料にない。周長 220 mm が筐体の幅 2.0 in（50.8 mm）に収まることから、両端の半円の R ≈ 15–20 mm と推定 [AE]。
- 曲げの複屈折（§1.7）: R = 15 mm で B_bend = 2.4×10⁻⁶（B の 0.6%）、R = 20 mm で 1.3×10⁻⁶。
- **PER の低下** [AE]: 曲がった部分では主軸が δ0 sin 2α（δ0 = B_bend/(2B)、α は曲げの面と slow 軸の角度）だけ回る。
  直線部と曲線部の境界ごとに振幅 ≈ δ0 sin 2α が移り、1 つの曲線部（入口と出口）の寄与は平均で δ0²（α がランダム）。
  1 巻に半円が 2 つなので、76 区間で P_x ≈ 76·δ0²。R = 15 mm: δ0 = 3.2×10⁻³ → **−31 dB**、R = 20 mm: δ0 = 1.8×10⁻³ → **−36 dB**。
  比較: Corning の曲げに強い PMBI 15 の規格「直径 15 mm、10 巻で曲げのクロストーク ≤ −30 dB」[Corning datasheet]。同じ式で標準 PANDA・R = 7.5 mm・20 区間なら −25 dB で、桁が合う。
  α が全区間で同じ（ファイバをねじらずに巻いた）で、区間の間の位相が揃う場合は振幅で足し合わされ、最悪で数倍〜一桁悪くなる。
  区間の間の位相は 1 巻あたり K_T × 0.22 m ≈ 0.45 rad/K で回るので、**巻線部の PER は温度で変わる**（数 K の周期）。
- **差動のアクチュエータ利得**（直交偏光の probe でストレッチャを使う場合）: 伸び δL で Δψ_b は β_b·δL·(1 + O(1)) 変わる。16 µm で約 0.024 rad、共通の位相（約 65 rad）に対して 4×10⁻⁴。
  2 軸のアクチュエータ利得の差は 10⁻⁴ 程度で無視できる [AE]。
- **曲げ損失**: PM15-U25D は MFD 10.5 µm、遮断波長 1300–1440 nm [Corning datasheet]（MFD/λ_c ≈ 7.7 で SMF-28 と同程度の曲げへの弱さ）。
  曲げ損失は R に指数的に効き、R = 15 mm と 25 mm で桁が違う。SMF-28 型の規格を代わりに使うと、38 巻で 0.01 dB（R = 25 mm）〜 2 dB（R = 16 mm）の幅になる（SMF-28 の規格値は UNVERIFIED）。
  **D10 で将来ストレッチャが squeezed 光の経路に入る**ので、この損失は squeezed 光にかかる。R の実測か問い合わせが必要（§7）。

### 3.5 Ref D の偏光の収支

光路（D02）: Cobrite DX1 → MN9605C 減衰器 → Optoquest PC → [以下 PM、FC/APC] → Opneti 偏光子 → PNH1550R5A2 → 腕 1（915B、約 8.4 m ＋ ピグテール）／腕 2（短い PM、長さ不明）→ PNH1550R5A2 → PDB450C。

部品の値:

| 部品 | 値 | 出典 |
|---|---|---|
| PM15-U25D（915B のファイバ） | L_b 3.0–5.0 mm、クロストーク ≤ −30 dB/100 m（h ≤ 1×10⁻⁵ /m）、典型 −40 dB @ 4 m、損失 ≤ 0.5 dB/km、MFD 10.5 µm | Corning PANDA PM datasheet M0500017（2019-01）、URL は §8 |
| 腕 2 の PM ファイバ | 種類不明。PM15-U25D と同じ値を仮定 | 仮定 |
| FC/APC の PM コネクタ（接続 1 組） | PER 20–30 dB（θ = 1.8–5.7°）| [AE]。メーカーの値がない。参考に、ある in-line 偏光子の資料は「コネクタ 1 つごとに ER が約 2 dB 下がる」と書く（UNVERIFIED） |
| Opneti 偏光子 | ER 不明。同種の製品（1550 nm、PM、in-line）の典型 25–30 dB | メーカー資料が見つからない（UNVERIFIED） |
| PNH1550R5A2 | PER・fast 軸の扱い（通すか止めるか）・結合比・過剰損失が不明 | Thorlabs のページが読めなかった（UNVERIFIED、§7） |
| Optoquest PC | 種類・挿入損失・PDL が不明 | メーカー資料が見つからない（UNVERIFIED） |
| PDB450C | 偏光に依存しない（フォトダイオード）とみなす | [AE] |

**寄生的な縞によるロック点のずれ**: PDB450C は両方の偏光の光を受けるので、差信号は

  i_− ∝ |a_s1||a_s2|·[cos Δφ_s + ε·cos(Δφ_s − δ_b)]，ε = |a_f1 a_f2|/|a_s1 a_s2|

（a_s,f は腕ごとの slow・fast 軸の振幅、Δφ_s は slow 軸の 2 腕の位相差、δ_b = Δψ_b,arm1 − Δψ_b,arm2 ＋（部品での位相）。符号は conventions §5 の cos(φ_LO − φ_pr) に合わせる。）
= Re[e^{iΔφ_s}(1 + ε e^{−iδ_b})] なので、零点は

  Δφ_s,lock = π/2 − arg(1 + ε e^{−iδ_b}) ≈ π/2 + ε·sin δ_b(t)

（厳密には δφ = −arg(1 + ε e^{−iδ_b})、最大値 arcsin ε。）ロックの傾き（K_D）も |1 + ε e^{−iδ_b}| の分（±ε）だけ変わる。

- ε の源: (i) 偏光子を通った後に残る fast 軸の成分（共通、両腕でコヒーレント）: ε ≈ 10^(−ER_pol/10)。
  (ii) 腕ごとに生じた fast 軸の成分の積 κ_1κ_2。(iii) 2 つ目のカプラで fast → slow に戻る成分（§1.5 の 2 点の干渉）。いずれも 2 つの結合の振幅の積で、PER 20–25 dB なら **ε ≈ 3×10⁻³–10⁻²**。
- δ_b(t) = δ_b0 + K_T·ΔL_arm·ΔT(t)。ΔL_arm ≈ 8–10 m（915B の腕が長いと仮定）なら |K_T|·ΔL_arm ≈ 17–20 rad/K、**0.3–0.4 K で 1 周**。
- 結果: **ロック点が ±ε ≈ ±(3–10) mrad の範囲で、温度に対して周期的に動く**（両腕の温度が共通に動く場合。腕ごとに違えば δ_b はそれぞれの ΔT で動く）。
  D06 の「縞の可視度の変化（偏波ずれを含む）」はこの項として入れる。可視度は V(1 ± ε) で揺らぐ。
- 偏光子の前の SMF（Cobrite → PC）の偏光のドリフトは、偏光子の後のパワー変化 δP/P ≈ −Θ²/4 になる（§1.4）。Ref D の 3.33 mW はこの設定点での値。
  Θ の時間変化には値がない（パラメータ、§3.6）。PDL のある部品（減衰器、PC）があると、一次の効果 δP/P ∝ PDL·Θ が加わる（PDL の値は不明）。

### 3.6 SMF の偏光のドリフトのパラメータ

- 曲げの複屈折: 半径 5 cm のスプールで B ≈ 2×10⁻⁷、63 m で Δψ_b ≈ 54 rad（PM の約 1/1800）[AE]。
- その温度係数: 曲げの複屈折の温度係数は光弾性係数・半径・巻きの張力（スプールの材質の熱膨張）で決まり、一次資料の値を見つけられなかった。
  k_T,SMF ~ 10⁻⁴ /K を仮定すると、63 m で 5 mrad/K（PM の約 2.5×10⁴ 分の 1）[AE、根拠は弱い]。被覆の応力などで大きくなりうる。
- 外乱によるドリフト: 文献の値は設置条件に強く依存する（例: 埋設した 7.1 km の SMF-28 で約 0.1°/h という報告がある、arXiv 2209.06920、UNVERIFIED）。
  **D_pol [rad²/s] を掃引のパラメータ**（例: 10⁻⁸–10⁻² rad²/s）とし、ユーザの実験室の値は D04・D27 と同じ扱い（Phase 5 の終わりまで使わない）にする。
- PMD: SMF-28 Ultra の規格は D_PMD ≤ 0.1 ps/√km（1 本）、リンクの設計値 ≤ 0.04 ps/√km（Corning PI1424、検索結果の抜粋で確認、PDF は開けず UNVERIFIED）。

## 4. 概算（自分たちのパラメータ）

### 4.1 Ref D（D02）

| 量 | 値 | 根拠 |
|---|---|---|
| 915B の巻線部の 2 軸の位相差の温度係数 | −17 rad/K（範囲 −14〜−24） | 8.4 m × K_T |
| PER 由来のロック点のずれ | ±3〜±10 mrad、0.3–0.4 K で 1 周 | §3.5、ε = 3×10⁻³–10⁻² [AE] |
| 縞の可視度の揺らぎ | ±ε ≈ ±(0.3–1)% | §3.5 |
| 偏光子の前の偏光の揺らぎ → パワー | δP/P ≈ −Θ²/4（Θ = 0.1 rad で −0.25%） | §1.4、Θ はパラメータ |
| 915B の巻線部の PER | −31〜−36 dB（R = 15–20 mm を仮定） | §3.4 [AE] |
| 915B の曲げ損失 | 0.01–2 dB（R 不明） | §3.4、UNVERIFIED |
| PM ファイバの分布クロストーク（10 m） | ≤ −40 dB | h ≤ 10⁻⁵ /m |

Ref D の DC side-fringe lock の位相の分解能（survey で 0.09 mrad 程度）に比べ、PER によるロック点のずれ（mrad）は大きい。
ただしロック点がずれるだけで、Ref D の 2 腕の位相差は「ロック点の周り」で固定される。D04 の比較では、温度の記録があれば温度と相関した遅いドリフトとして見えるはず。

### 4.2 将来の構成（D07–D13、D30–D32）

波長は Ref A の 1545.32 nm で代用（D07）。K_T の λ 依存（1/λ）による差は 0.2% で無視できる。f_rep が Open なので、L_k は survey の目安（f_rep = 100 MHz なら L_1 = 2.04 m、L_31 = 63 m。10 MHz なら L_31 = 630 m）を使う。

| 量 | L = 2.04 m | L = 63 m | L = 630 m | 根拠 |
|---|---|---|---|---|
| PM の 2 軸の位相差の温度係数 K_T·L [rad/K] | −4.2 | −129（−108〜−181） | −1290 | §2 |
| 直交偏光の probe で Θ_err < 10 mrad に必要な温度の安定度 | 2.4 mK | **77 µK** | 7.7 µK | 0.01/(|K_T|·L) |
| 同 0.1 K/h のドリフト（D27 の例）での誤差の速さ | 0.12 mrad/s | 3.6 mrad/s | 36 mrad/s | K_T·L·dT/dt |
| PM の DGD τ_DGD | 2.8 ps | 86 ps | 0.86 ns | 1.36 ps/m |
| DGD の温度係数 | −3.6 fs/K | −0.11 ps/K | −1.1 ps/K | K_τ |
| PM の分布クロストーク h·L（h = 10⁻⁵ /m） | −47 dB | −32 dB | −22 dB | Corning の規格の上限 |
| PM・同じ偏光・CW probe とパルスの位相差（PER 25 dB の接続 2 か所） | ≤ 3.2 mrad | 同 | 同 | κ_iκ_j、§1.5 |
| SMF の平均 DGD（0.1 / 0.04 ps/√km） | 4.5 / 1.8 fs | 25 / 10 fs | 79 / 32 fs | §1.7 |
| SMF の PMD による η の低下（最悪、τ_p = 100 fs） | < 0.2% | 4.2% / 0.7% | 32% / 6.6% | §1.7 [AE] |
| 同（τ_p = 1 ps） | 無視 | 0.04% / 0.01% | 0.4% / 0.07% | §1.7 [AE] |
| SMF・同じ偏光・CW probe とパルスの位相差（τ_p = 100 fs） | 4×10⁻⁴ | ≤ 1×10⁻² | ≤ 1×10⁻¹ | (τ_DGD/2σ_t)²/8 [AE] |
| SMF の曲げの複屈折による 2 軸の位相差（R = 5 cm） | 1.8 rad | 54 rad | 540 rad | §3.6 [AE] |
| Faraday（最悪、地磁気が全長に平行） | 6×10⁻⁵ rad | 2×10⁻³ rad | 2×10⁻² rad | §1.8 [AE] |

PMD の「パルスで平均した位相」の概算は粗い（最悪の分け方、ガウスパルス）。100 fs・630 m の SMF では PMD がパルスの squeezed 光の主要な問題になりうる。

### 4.3 Liu 型のロック（D30）での PM と SMF の比較（D32 への入力）

ロックが固定するのは probe の位相。表の「誤差」は、squeezed 光（GBS 入力の偏光子に射影）の位相と probe が測る位相の差（指標 (ii) に直接入る）。
L = 63 m、準静的（通過時間の逆数 約 3 MHz より十分低い周波数）での値。逆向きの伝搬による高い周波数の相関の崩れ（E42）と probe の混入（E43、E44）は T10。

| 項目 | PM・同じ偏光 | PM・直交偏光（Liu et al. 2026 の方式） | SMF・同じ偏光 | SMF・直交偏光 |
|---|---|---|---|---|
| 偏光による位相の誤差の主因 | クロストークの 2 点の干渉（CW probe とパルスの違い） | 2 軸の複屈折の温度ドリフト | PMD（パルスの帯域での平均）、非相反の項 | 2·arg(a(t))（偏光の回転） |
| 大きさ（63 m） | ≤ Σκ_iκ_j ≈ 10⁻³–10⁻² rad、温度で周期的 | **129 rad/K × ΔT**（10 mrad には 77 µK） | ≤ 10⁻² rad（100 fs）、≤ 10⁻⁴ rad（1 ps）、Faraday ≤ 2×10⁻³ rad | 予測できない。曲げの複屈折 54 rad × k_T,SMF ΔT ＋ 環境による回転 |
| 温度の制御 | 不要（PER のリップルの位相が回るだけ） | **必須**（µK〜mK、Liu は TEC のフィードバックを併用） | 不要 | 必要（値は不明） |
| squeezed 光の偏光損失 | κ_tot² + h·L（−30 dB 程度、静的） | 同左 | **能動的な偏光の制御が必要**（|a|² が動く）。probe のパワーがそのまま偏光損失の誤差信号になる | 同左、ただし probe のパワー（|a|²）からも分かる |
| probe の信号のフェージング | なし（静的） | なし | |a|² で消えうる（偏光の制御で防ぐ） | 同左 |
| probe と squeezed 光の分離 | 偏光では分けられない（時間ゲート、周波数のずれ、サーキュレータ。T10） | PBS で分けられる（PER の分だけ混入、E43） | 偏光では分けられない（T10） | PBS で分けられるが、偏光が保たれないので混入が時間とともに変わる |
| パルスの時間モード | DGD 86 ps でクロストーク光はサテライトになる（損失） | probe は fast 軸を 86 ps 遅れて走る（AOM のパルス化の時間スケール µs に比べ無視） | PMD による η の低下 0.7–4%（100 fs）| 同左 |
| 結論 [AE] | 偏光の点では最も良い。課題は T10 の混入 | 長い L では温度の制御なしには使えない | 位相は相反性で一致。偏光の能動的な制御と PMD（fs パルス）が課題 | 誤差が予測できず、最も悪い |

注:
- PM・直交偏光で温度を測ってフィードフォワードで補正する場合、残差は |K_T| の不確かさ（10–30%、§2.3）と温度センサーの位置の偏り（ファイバの平均温度との差）で決まる。
  63 m・1 K の変化なら 13–40 rad の補正の誤差が残るので、温度の安定化（Liu 型の TEC）の併用が前提になる。
- 同じ偏光の probe は、相反性で偏光の影響を消せる代わりに、probe と squeezed 光を偏光で分けられない。この長所と短所の大小は T10（混入光子数）と GBS_REQUIREMENTS（Open）で決まる。

## 5. 検証テスト（実装が通すべきもの）

1. **PM の差動ドリフト（D11 の数値の固定）**: PM 区間 L = 63 m、L_b = 4.2 mm、k_T = −1.37×10⁻³ /K、λ = 1545.32 nm。温度を ΔT = +1 mK のステップで変える。
   直交偏光の probe をロックした状態（理想的なロック）で、squeezed 光の位相のずれが K_T·L·ΔT = −0.1291 rad（相対誤差 10⁻⁶ 以内）。同じ偏光の probe（クロストークなし）では 0（10⁻¹² 以内）。
2. **相反性**: 長さ 63 m の SMF を 100 区間のランダムな SU(2)（固定の乱数の種）で作る。順方向の t_s = g^T J e_in と逆方向の t_p = e_in^T J^T g で、arg t_s − arg t_p = 0（10⁻¹² 以内）、|t_s| = |t_p|。
   直交偏光では arg t_s − arg t_p = 2 arg(a)（a = x^T U x、10⁻¹² 以内）。U を PM（対角）にすると 2 arg(a) = Δψ_b。
3. **コネクタの PER**: 軸ずれ θ = 3° の接合 1 か所で、fast 軸のパワーの比が tan² 3° → PER = 25.61 dB（0.01 dB 以内）。
4. **分布クロストーク**: h = 10⁻⁵ /m、L = 100 m で P_x/P = tanh(10⁻³) → −30.00 dB（多数のランダムな小さい結合点で作る場合は、十分な数の実現の平均で ±0.5 dB 以内）。
5. **Ref D のロック点**: ε = 0.01 で δ_b を 0 から 2π まで掃引すると、零点の位置のずれの最大値が arcsin(0.01) = 10.0002 mrad、波形が −arg(1 + ε e^{−iδ_b})（10⁻⁹ 以内）。ε = 0 で 0。
6. **偏光の重なり**: Poincaré 球上で Θ = π/2 離れた 2 つの偏光で η_pol = 0.5、Θ = 0.2 rad で 0.990（10⁻⁶ 以内）。偏光子（ER = ∞）の前で Θ = 0.1 rad なら P_out/P_in = cos²(0.05) = 0.997502。
7. **CW とパルスの違い**: PM の 2 か所の接合（κ_1 = κ_2 = 0.0562、PER 25 dB）の間を ℓ = 1 m とする。CW の probe（同じ偏光）では位相のずれの振幅が κ_1κ_2 = 3.16×10⁻³ rad（Δψ_b を掃引したときの最大値、1% 以内）。
   τ_p = 1 ps のパルス（τ_12 = 1.36 ps ＞ τ_p）では主パルスの位相のずれが 10⁻⁴ 倍以下に小さくなり、遅れ 1.36 ps の位置にパワー比 κ_1²κ_2² = 10⁻⁵ のサテライトが出る。
8. **トグル**: 偏光の効果をすべてオフにすると、軸ごとの位相が一致し（ψ_fast = ψ_slow）、Jones 行列が e^{iψ}·1 になり、既存の単一偏光のテスト（prompt §7-11 の共通光路の打ち消し）がそのまま通る。

## 6. シミュレータが必要とするもの

### 6.1 State の構造について（conventions §9 への提案）

conventions §9 は区間に偏光の「軸」 a ∈ {slow, fast}（PM）または {x, y}（SMF）を持たせる。PM ではこれで足りる（クロストークは摂動として §1.5 の項を足す）。
**SMF では軸が保たれないので、軸のラベルでは §1.6 の相反性も 2·arg(a) の誤差も計算できない。** 提案: 区間ごとに 2×2 の複素 Jones 行列 J_seg(t)（PM は対角 ＋ 小さな非対角）を持たせ、
ビームの複素振幅 = 射影ベクトル^T · (区間の J の積) · 入射の偏光、逆向きは転置の積、とする。軸のラベルは「PM で対角だけを使う速い経路」として残せる。
これは D11（後から State を変えると高くつく）に関わるので、調整役の判断を求める（§7）。

### 6.2 ブロックと表

| 段 | ブロック / 表 | 入力 | 出力 | トグル |
|---|---|---|---|---|
| O | `pm_fiber_params`（静的な表） | L_b（または B）、k_T、b_g、K_τ、h、λ | β_b、K_T、τ_DGD、軸ごとの温度係数 | — |
| O | `per_budget`（部品ごとの κ_k と位置） | 部品の PER（コネクタ、融着、偏光子 ER、カプラ PER、巻線部）、位置 z_k | κ_k、合計の PER、偏光損失 η_PER（squeezed 光に純損失として渡す、T12） | `pol_crosstalk` |
| O | `coil_birefringence`（ストレッチャ・スプール） | R、巻き数、1 巻の長さ、α の分布 | 曲げの PER、巻線部の曲げ損失（R が不明なのでパラメータ） | `pol_bend_coil` |
| O | `pmd_mode_overlap` | D_PMD（または区間の J(ω)）、L、τ_p | η_PMD、CW とパルスの位相差 | `pol_pmd` |
| T | `segment_jones`（区間の雑音源の一部。T08 の温度の雑音源から温度を受ける） | T_seg(t)、区間の種類（PM/SMF）、D_pol | J_seg(t)（PM: 軸ごとの ψ、Δψ_b(t) = Δψ_b0 + K_T L ΔT） | `pol_birefringence_drift`、`pol_smf_drift` |
| T | `polarizer_projection` | 入射の偏光、通過軸、ER | 透過パワー（偏光子の後のパワーのドリフト）と射影した振幅 | `pol_drift_before_polarizer` |
| T | `pdl_element` | PDL_dB、主軸 | 偏光に依存する透過率 | `pol_pdl` |
| T | `faraday_nonreciprocal`（既定オフ） | V_F、B_∥(t)、L | 順・逆で符号が逆の位相 | `pol_faraday` |
| M | 偏光の指標の出力 | | η_pol(t)、Θ_err(t)（probe と squeezed 光の偏光による位相差）、PER の合計 | — |

- 検出器（PDB450C、BHD）の入力では、偏光ごとの振幅を足し合わせて差信号を作る（偏光を区別しない受光）。BHD は LO の偏光への射影。
- 温度の時間変化 T_seg(t)、相関、PSD は T08 のブロックから受け取る。ここでは係数だけを持つ。
- 逆向きの伝搬の時刻のずれ（conventions §9）は T10 の遅延バッファで扱う。ここでの J_seg(t) はその時刻で評価する。

## 7. ユーザの決定が必要な点・未解決

1. **State の偏光の表し方**（§6.1）: SMF を比較対象（D32）にするなら、区間の偏光を「軸のラベル」でなく「Jones 行列」で持つ必要がある。conventions §9 と D11 の変更になる。
2. **915B の曲げ半径と、巻線部の損失・PER**: 資料にない。将来 squeezed 光の経路にストレッチャを入れる（D10）なら、R と挿入損失の実測かメーカーへの問い合わせが必要。実測するなら D04 との関係（Ref D の予測に使ってよいか）も決める。
3. **部品の資料**: PNH1550R5A2（PER、fast 軸を通すか止めるか、結合比）、Opneti 偏光子（ER、型番）、Optoquest PC（型番、種類、PDL）のメーカー資料を `refs/manuals/` に置いてほしい。
   Thorlabs のページは自動では読めなかった。それまでは同種の製品の典型値（§3.5）を [AE] として使う。
4. **K_T の実測**: 文献値に 6 倍の幅がある（§2.3）。PM15-U25D の K_T を実測してよいか（Ref D のブラインド検証の対象そのものではないが、Ref D の予測に入る量）。
5. **Liu et al. 2026 のファイバの種類**: 本文（Methods）には遅延ループの「fibre」の種類（PM か SMF か）が書かれていない。Supplementary を確認したい（`refs/` にない）。直交偏光の probe が PM で働いているなら、TEC の安定度の要求（§4.2）と照らして参考になる。

## 8. 文献

一次資料（本文を読んだもの）:
- T. Ritari, H. Ludvigsen, M. Wegmuller, M. Legré, N. Gisin, J. R. Folkenberg, M. D. Nielsen, "Experimental study of polarization properties of highly birefringent photonic crystal fibers,"
  Opt. Express 12, 5931–5939 (2004)。DOI 10.1364/OPEX.12.005931（PDF には DOI の記載なし。Optica の形式から組み立てたもので UNVERIFIED）。全文は Aalto の公開版
  https://aaltodoc.aalto.fi/server/api/core/bitstreams/c7398e11-178c-40fb-9da0-6521992265ec/content（2026-10-28 閲覧）。
- F. Dvořák, J. Maschke, Č. Vlček, "Utilization of birefringent fiber as sensor of temperature field disturbance," Radioengineering 18(4), 639–643 (2009)。
  https://www.radioeng.cz/fulltexts/2009/09_04_639_643.pdf（2026-10-28 閲覧）。温度係数は同論文の文献 [1] の引用（二次）。
- Liu et al., "Gaussian boson sampling with 1,024 squeezed states in 8,176 modes," Nature 653 (2026)（`refs/`、Methods「Phase stabilization of the system」）。
- S. Verclas, C. Darsow-Fromm, B. Tohermes, R. Schnabel, Optica Quantum 4, 1 (2026)、DOI 10.1364/OPTICAQ.567418（Ref C、Fig. 1 の PC）。

メーカー資料（D24）:
- Evanescent Optics / ハナムラオプティクス「Model 914-2 Piezo Controller / 915B Piezo Assembly」`refs/manuals/914-915B-Specifications.pdf`（巻き数 38、PM15-U25D、16/17 µm、FC/APC slow 軸キー）。
- Corning, "PANDA PM — High Performance Polarization Maintaining Fibers," M0500017, issued January 2019（PANDA PM Bend Insensitive M0500021 を含む）。
  AMS Technologies の転載 https://shop.amstechnologies.com/media/b9/6b/7b/1720716991/PANDA-PM-Polarization-Maintaining-Optical-Fibers-Corning-Datasheets.pdf（2026-10-28 閲覧）。
  Corning の元の URL https://www.corning.com/microsites/coc/oem/documents/specialty-fiber/PANDA-PM.pdf は自動取得できなかった。
- photonics.com「Model 915B Piezo-electric Fiber Stretcher」https://www.photonics.com/Buyers-Guide/Products/Model-915B-Piezo-electric-Fiber-Stretcher/psp9424
  （2026-10-28、検索エンジンの抜粋のみ。ページ本体は 403。**UNVERIFIED**）。
- Corning SMF-28 Ultra Product Information PI1424（PMD ≤ 0.1 ps/√km、LDV ≤ 0.04 ps/√km）。
  https://csmedia.corning.com/media/worldwide/coc/documents/Fiber/product-information-sheets/PI-1424-AEN.pdf（2026-10-28、検索結果の抜粋のみ、PDF は取得できず **UNVERIFIED**）。
- Thorlabs PNH1550R5A2: https://www.thorlabs.com/thorproduct.cfm?partnumber=PNH1550R5A2（2026-10-28、自動取得で内容が読めず、**値はすべて UNVERIFIED**）。
- Opneti 偏光子、Optoquest 偏波コントローラ: メーカーの資料を見つけられなかった（**UNVERIFIED**）。

書誌だけ確認したもの、または未確認のもの:
- F. Zhang, J. W. Y. Lit, "Temperature and strain sensitivity measurements of high-birefringent polarization-maintaining fibers," Appl. Opt. 32, 2213–2218 (1993)（書誌は複数の索引で確認。本文と数値は UNVERIFIED。DOI 10.1364/AO.32.002213 は形式からの推定で UNVERIFIED）。
- J. P. Gordon, H. Kogelnik, "PMD fundamentals: Polarization mode dispersion in optical fibers," Proc. Natl. Acad. Sci. USA 97, 4541–4550 (2000)、DOI 10.1073/pnas.97.9.4541（書誌は確認、本文は未読）。
- R. Ulrich, S. C. Rashleigh, W. Eickhoff, "Bending-induced birefringence in single-mode fibers," Opt. Lett. 5, 273–275 (1980)、DOI 10.1364/OL.5.000273（**UNVERIFIED**: 著者・DOI とも記憶による。式の係数 0.136 も要確認）。
- S. C. Rashleigh, "Origins and control of polarization effects in single-mode fibers," J. Lightwave Technol. 1, 312–331 (1983)（**UNVERIFIED**）。
- M. Martinelli, "A universal compensator for polarization changes induced by birefringence on a retracing beam," Opt. Commun. 72, 341–344 (1989)（**UNVERIFIED**、相反性 J_back = J^T の文献として）。
- SOP のドリフトの例: arXiv:2209.06920（埋設 7.1 km SMF-28 で約 0.1°/h、検索結果の抜粋のみ、**UNVERIFIED**）。
