---
type: discussion
project: gleamworks
date: 2026-10-07
status: plan
tags: [gleamworks, simulation, implementation-plan, python, jax, nlse, quantum-noise, squeezed-light, fiber, phase-stabilization, gawbs, raman]
related:
  - "[[2026-10-07 phase-stab-sim_plan-review]]"
  - "[[Squeezed States]]"
  - "[[Homodyne Detection]]"
  - "[[Coherent Probe Locking]]"
  - "[[Phase Space Representation]]"
---

# 実装計画：ファイバを伝搬するスクイーズド光の量子雑音と位相安定化のシミュレーション（Python）

関連: [[2026-10-07 phase-stab-sim_plan-review]]（計画書のレビュー、制御層のプロトタイプ）

出典の表記: `(paper)` 論文・資料に書いてある、`(agent)` 私（Claude）の推論・概算・導出（式や計算で確かめられる）、`(me)` 自分の判断。

種コード（動作確認済み）:

- `assets/2026-10-07_nlse_linearized_seed.py` — RK4IP による NLSE、jacfwd による線形化、共分散、Bloch–Messiah（10 節）
- `assets/2026-10-07_phase_lock_budget.py` — 制御層の周波数領域・時間領域モデル（前回の Discussion）

---

## 0. 位置づけ

- **方針 (me)**: このシミュレーションは Moku 実験のための道具ではなく、それと関連した **独立した一つの研究** として行う。Moku 実験の参考になればよい。したがって非線形ソルバも含める。
- **この計画の役割 (agent)**: 自分で Python 実装するときの設計図。各段階に「何を作るか」「どう確かめるか（合格基準）」を付けた。
- **前回からの変更点 (agent)**:
  - 前回は「CW なら非線形ソルバは後回しでよい」とした。今回は研究として、**どの条件で非線形効果が効き始めるか** を調べること自体を目的に含める。
  - CW（サイドバンド）とパルス（時間領域）の **2 つの計算レジーム** を用意し、両者が重なる条件で一致を確かめる。

---

## 1. 研究としての問いと成果物

### 1.1 問い

- **RQ1（伝搬による劣化）**: ファイバ長 $L$、パルス幅、LO の送り方（別ファイバ・同一ファイバ・偏光多重）に対して、測定されるスクイージングはどれだけ劣化するか。その劣化は「損失」「モード不整合（分散）」「非線形（SPM/XPM）」「過剰雑音（Raman, GAWBS）」のどれで決まるか。
- **RQ2（広帯域 CW の分散回転）**: 広帯域スクイーズド光を CW の LO で測るとき、分散によってスクイーズ角が周波数ごとに回る（3.5 節）。どの帯域・長さまで使えるか、どう補償するか。
- **RQ3（ロックできる雑音・できない雑音）**: ファイバ由来の位相雑音を周波数帯で分けると、どれがフィードバックで消せて（kHz 以下の熱・音響）、どれが消せないか（MHz 以上の GAWBS など）。Moku の遅延を含むループで、消せる部分をどこまで消せるか。
- **RQ4（設計則）**: 上の結果から、probe/LO のパワー、偏光、周波数オフセット $\Omega$、ファイバ長の上限などの設計則を導けるか。

### 1.2 成果物

- Python パッケージ（仮名 `fibersqz`）: テスト付き、再現可能な設定ファイル
- 計算実験 E1–E5（6 節）の図と表
- 研究ノート（この vault）と、可能なら短い報告書

---

## 2. 全体像

```text
                ┌──────────────────────────────────────────────────────────┐
  レジーム P    │ L1 古典場: GNLSE (分散・SPM/XPM・Raman・自己急峻化・損失) │
 （パルス, 時間 │ L2 量子揺らぎ: 線形化 (jacfwd/jvp/vjp) ＋ 共分散            │
   領域グリッド）│      └ 検算: truncated Wigner (TWA) のアンサンブル        │
                └──────────────────────┬───────────────────────────────────┘
                                       │ η_mm, 過剰雑音, パルスごとの位相
                ┌──────────────────────▼───────────────────────────────────┐
  レジーム C    │ サイドバンド Ω ごとの 2×2 共分散: 損失, 分散回転 θ(Ω),     │
 （CW, 周波数   │ XPM, GAWBS の線スペクトル, Raman 自然散乱               │
   領域）        └──────────────────────┬───────────────────────────────────┘
                                       │ 位相雑音 PSD（ロック可能/不可能に分ける）
                ┌──────────────────────▼───────────────────────────────────┐
  L3 測定と制御 │ ホモダイン射影, Bloch–Messiah, 検出効率                    │
                │ ループ（Moku 遅延, EOM/AOM, オフロード）→ θ_rms → 測定 dB │
                └──────────────────────────────────────────────────────────┘
```

- **L1・L2 を分けない理由 (agent)**: 線形化の行列は古典場の写像のヤコビアンそのもの。古典場ソルバを JAX で書けば、量子揺らぎの伝搬は自動微分で「ただで」得られる。
- **GAWBS を L1 の時間グリッドに入れない理由 (agent)**: GAWBS の線は数十 MHz–GHz で、ps のパルス窓（数十–数百 ps）では「パルスごとにほぼ一定の位相」にしかならない。パルスでは「パルス間の位相揺らぎ」として、CW ではサイドバンドの過剰雑音としてレジーム C で扱う。

---

## 3. 物理モデル

### 3.1 古典場: 一般化非線形シュレーディンガー方程式（GNLSE）

$$
\frac{\partial A}{\partial z} = -\frac{\alpha}{2}A + \sum_{k\ge 2}\frac{i^{k+1}\beta_k}{k!}\frac{\partial^k A}{\partial T^k}
+ i\gamma\left(1+\frac{i}{\omega_0}\frac{\partial}{\partial T}\right)\left[A\int_0^\infty R(t')\,|A(T-t')|^2\,dt'\right]
$$

- $A(z,T)$: パルスの包絡線 [$\sqrt{\mathrm W}$]。$|A|^2$ が瞬時パワー
- $z$: 伝搬距離、$T$: 群速度で動く座標系での時間
- $\alpha$: パワーの損失係数 [1/m]（0.2 dB/km なら $\alpha \approx 4.6\times10^{-5}$ /m）
- $\beta_k$: 分散係数。$\beta_2$ は群速度分散 [ps²/m]（SMF, 1550 nm で約 −22 ps²/km (agent; 典型値)）
- $\gamma$: 非線形係数 [1/(W·m)]（SMF で約 1.3 /(W·km) (agent; 典型値)）
- $\omega_0$: 中心角周波数。$\frac{i}{\omega_0}\partial_T$ は自己急峻化（数 100 fs 以上のパルスではほぼ無視できる）
- $R(t) = (1-f_R)\delta(t) + f_R h_R(t)$: 非線形応答。前半が瞬時の Kerr、後半が Raman
- Raman 応答（Blow & Wood 型）(paper: Blow & Wood 1989):

$$h_R(t) = \frac{\tau_1^2+\tau_2^2}{\tau_1\tau_2^2}\,e^{-t/\tau_2}\sin(t/\tau_1),\qquad \tau_1 = 12.2\ \mathrm{fs},\ \tau_2 = 32\ \mathrm{fs},\ f_R = 0.18$$

  より正確な形として Lin & Agrawal 2006 のモデル（boson peak を含む）も選べるようにする (paper)。

**数値解法**: RK4IP（相互作用描像の 4 次ルンゲ＝クッタ）(paper: Hult 2007)。

- **動作**: 線形部分（分散・損失）は周波数領域で厳密に解き、非線形部分だけを RK4 で積む。
- **意味**: 刻み幅 $h$ に対して誤差が $O(h^4)$。対称スプリットステップ（$O(h^2)$）より同じ精度でずっと少ないステップ数で済む。
- **実用上の目的**: 線形化で得る行列の精度（シンプレクティック性）を $10^{-9}$ 程度まで保つため。種コードで $|S^T\Omega S-\Omega| \approx 5\times10^{-10}$ を確認済み (agent)。

**特徴的な長さ** (agent):

- 分散長 $L_D = T_0^2/|\beta_2|$（$T_0$ はパルスの半幅）: 1 ps → 約 45 m、100 fs → 約 0.45 m
- 非線形長 $L_{NL} = 1/(\gamma P_0)$（$P_0$ はピークパワー）: 1 W → 約 770 m
- ソリトン次数 $N^2 = L_D/L_{NL}$
- **意味**: $L \ll L_D, L_{NL}$ なら線形・無分散の近似でよい。どちらかに近づいたらソルバが要る。この判定を自動で出す関数を最初に作る。

### 3.2 量子揺らぎ: 線形化と共分散行列

**基本の考え方**: 場を「古典場＋小さな量子揺らぎ」に分ける (paper: Haus & Lai, Carter et al. の線形化の考え方)。

$$\hat A(z,T) = A(z,T) + \hat a(z,T)$$

$\hat a$ は 1 次まで残すと線形の方程式に従う。Kerr 項だけ書くと (agent):

$$\frac{\partial \hat a}{\partial z} = \mathcal D\,\hat a + i\gamma\left(2|A|^2\hat a + A^2\hat a^\dagger\right) \;(+\ \text{Raman の項})$$

- $\mathcal D$: 分散・損失の線形演算子
- $2|A|^2\hat a$: 古典場による位相変調（XPM 的な項）
- $A^2\hat a^\dagger$: 位相共役の項。これが揺らぎを **スクイーズ** する

**実ベクトル表現と伝達行列**:

- 時間グリッド $T_j$（$j=1..N$, 刻み $\Delta t$）上で、揺らぎを実ベクトル $y = \sqrt{4\Delta t}\,(\mathrm{Re}\,\delta a_1..\mathrm{Re}\,\delta a_N, \mathrm{Im}\,\delta a_1..\mathrm{Im}\,\delta a_N)$ にまとめる。こう正規化すると **真空の共分散が単位行列** になる（10 節のコードで確認）(agent)。
- 伝搬は線形写像 $y_\mathrm{out} = S\,y_\mathrm{in}$。$S$ は古典場の写像 $A_\mathrm{in}\mapsto A_\mathrm{out}$ の（実表現の）ヤコビアンそのもの。
- 共分散は $V_\mathrm{out} = S\,V_\mathrm{in}\,S^T$。

**JAX での計算法** (agent):

| やりたいこと | JAX の道具 | コスト |
|---|---|---|
| $S$ 全体（Bloch–Messiah 用） | `jax.jacfwd`（＝ $2N$ 本の JVP） | 伝搬 $2N$ 回分 |
| 入力側の少数モード（スクイーズド入力の $K$ モード）だけ | `jax.jvp` を $K$ 本、`vmap` | 伝搬 $2K$ 回分 |
| 出力側の測定モード（LO の時間モード $M$ 個）だけ | `jax.vjp` を $M$ 本 | 伝搬 $\sim M$ 回分 |

- **使い分け**: 測定量は $\mathrm{Var}(X_\theta) = w^T S V_\mathrm{in} S^T w$。$g = S^T w$ を VJP 1 本で得れば、$g^T V_\mathrm{in} g$ で済む。$N$ が大きくても安い。
- **注意**: 明るい場がファイバ中にあると、真空の入力も Kerr でスクイーズされる（$SS^T \ne I$）。入力の少数モードだけを追う JVP 方式は、明るい場がない（$S$ が直交行列）ときにだけ使える。

**損失**: 損失は真空を混ぜる非ユニタリな過程。$z$ を $K$ 区間に分け、区間ごとに (agent)

$$V \leftarrow S_k V S_k^T,\qquad V \leftarrow \eta_k V + (1-\eta_k) I,\qquad \eta_k = e^{-\alpha L/K}$$

- 0.2 dB/km の損失は非線形長に比べてずっとゆっくりなので、$K$ は少なくてよい。$K$ を倍にして結果が変わらないことで確かめる。
- 区間ごとの $S_k$ 全体が必要になるので、$N$ が大きいときは「少数モード基底に射影してから」使う。

**シンプレクティック性**: 損失がなければ $S^T\Omega S = \Omega$（$\Omega$ は $\begin{pmatrix}0&I\\-I&0\end{pmatrix}$）。ソルバとヤコビアンの正しさの一番強い検算になる。

**ホモダイン射影**: LO の時間モード $u(T)$（$\sum|u_j|^2\Delta t = 1$）と位相 $\theta$ に対し

$$X_\theta = \int\left(e^{-i\theta}u^*(T)\,\hat a(T) + \mathrm{h.c.}\right)dT,\qquad \mathrm{Var}(X_\theta) = w_\theta^T V w_\theta$$

- 真空で $\mathrm{Var}=1$（shot noise 単位）。$w_\theta$ は $u$ と $\theta$ から決まる長さ 1 のベクトル。

**Bloch–Messiah 分解** (paper: Braunstein 2005):

$$S = O_1\,\mathrm{diag}(s_1..s_N, 1/s_1..1/s_N)\,O_2$$

- $O_1, O_2$: 直交かつシンプレクティック（受動的な線形光学）
- $s_k$: スーパーモード $k$ のスクイーズ因子（$10\log_{10}s_k^2$ dB）
- **意味**: 理想的な LO 形状で得られる最良のスクイージング。実際の LO との差が「モード不整合による損失」 (agent)。
- **実装**: 損失なしの $S$ なら SVD の特異値がそのまま $s_k$, $1/s_k$ になる（種コードで確認）。損失ありの場合は共分散の固有値で評価する。

### 3.3 確率的な方法: truncated Wigner（TWA）

**動作** (paper: Drummond & Corney 2001 の考え方; agent):

1. 入力場に真空雑音を足す: グリッド点ごとに複素ガウス雑音 $\delta a_j$、$\langle|\delta a_j|^2\rangle = 1/(2\Delta t)$（光子流束の単位。パワー単位なら $\hbar\omega_0$ 倍）
2. スクイーズド入力なら、共分散 $V_\mathrm{in}$ に従う雑音を足す
3. 各ステップで、損失による真空、Raman の熱雑音を足す
4. 非線形の GNLSE をそのまま（線形化せずに）解く
5. $M$ 本の軌道の統計から共分散を推定する

- **意味**: 線形化では捨てた高次の効果（強いスクイージングでの曲がり、長距離での揺らぎの増大）も入る。線形化の検算と、線形化が破れる条件の判定に使う。
- **実用上の目的**: `jax.vmap` で軌道を並列化できる。分散 1 % 精度（約 0.04 dB）には $M \sim 2\times10^4$、0.1 dB 程度なら $M\sim 4\times10^3$（$\sqrt{2/M}$ の相対誤差）(agent)。
- **注意**: Wigner 関数は対称順序なので、正規順序の量（光子数など）を出すときは 1/2 の補正が要る。線形化と TWA は、雑音が小さい極限で一致しなければならない（テスト）。

### 3.4 ファイバ固有の雑音源

#### Raman 自然散乱（熱雑音）

- **動作**: Raman 応答の虚部 $\mathrm{Im}\,\tilde h_R(\Omega)$ に比例する揺らぎ・散逸。熱的な占有数 $n_\mathrm{th}(\Omega) = 1/(e^{\hbar\Omega/k_BT}-1)$ を持つ雑音が、明るい場を通して揺らぎに入る (paper: Drummond & Corney 2001。雑音相関の正確な形は原典で確認すること)。
- **数値の感覚** (agent): 室温で $n_\mathrm{th}$ は 1 THz で約 5.7、100 GHz で約 61、1 GHz で約 6100。
- **意味** (agent): $\Omega\to 0$ で $\mathrm{Im}\,\tilde h_R \propto \Omega$、$n_\mathrm{th} \approx k_BT/\hbar\Omega$ なので積は有限に残る。**明るい場のすぐ近く（LO の周波数付近）にも Raman 雑音は入る**。
- **実装**: L2（線形化）では拡散項 $dV/dz = \ldots + D_R(z)$ として、TWA では各ステップの乱数として入れる。明るい場がなければ寄与はない。

#### GAWBS（導波音響波ブリルアン散乱）

- **動作**: クラッドの熱的な音響振動（R₀ₘ: 径方向、TR₂ₘ: 径方向とねじれの混合）が屈折率を揺らし、通る光に位相変調・偏光変調をかける (paper: Shelby, Levenson, Bayer 1985)。
- **周波数** (agent; Shelby 1985 の境界条件で要確認): R₀ₘ は $f_m = y_m v_L/(2\pi a)$。$a$ はクラッド半径、$v_L$ は縦波音速（シリカで約 5.97 km/s）、$y_m$ は境界条件の根。根の間隔は約 $\pi$ なので線の間隔は $v_L/(2a)$。SMF（直径 125 µm）で約 48 MHz 間隔、〜1 GHz 付近まで。
- **大きさ**: 位相変調の分散は $L$ と温度 $T$ に比例（熱的な励起なので）(agent)。線幅は被覆による減衰で決まる (paper: Elser et al. の観察)。
- **意味**: 明るい場（LO, probe, ポンプ）の GAWBS 側帯波が、スクイーズド光の測定帯に入ると過剰雑音になる (paper: Liu et al. 2024 は GAWBS を「大きな障害」とし、測定帯をずらして回避)。スクイーズド真空だけが通るなら、GAWBS は楕円をわずかに回すだけ（µrad 級）で小さい (agent)。
- **実装**: レジーム C で、Lorentz 線の和として位相雑音 PSD $S_\phi^G(f) = L\,T\sum_m c_m\,\mathcal L(f; f_m, \Gamma_m)$ を作る。係数 $c_m$ は文献値に合わせる（未決, Open questions）。

#### 熱的位相雑音（基本限界）

- **動作**: ファイバの温度の熱力学的な揺らぎが熱光学効果で位相雑音になる。PSD は $L$ に比例 (paper: Wanser 1992; Duan 2012 が一般的な扱い)。
- **意味**: 環境雑音（音響・温度ドリフト）の下限。実験室では環境雑音の方が大きいことが多いので、**下限の参照線** として使い、環境雑音は実測（Moku の Phasemeter）か仮置きモデルで与える (agent)。

### 3.5 レジーム C: CW のサイドバンドモデル

CW のスクイーズド光は、キャリア $\omega_0$ の両側 $\pm\Omega$ のサイドバンド対の相関として測る。$\Omega$ ごとに 2×2 の quadrature 共分散 $V(\Omega)$ を持つ (agent)。

**損失**: $V \leftarrow \eta V + (1-\eta)I$（3.2 節と同じ）

**分散による回転** (agent; 導出):

- ファイバは $\pm\Omega$ のサイドバンドにそれぞれ位相 $\beta(\omega_0\pm\Omega)L$ を与える。
- スクイージングの相関 $\langle a(\Omega)a(-\Omega)\rangle$ には和 $\beta(\omega_0+\Omega)L+\beta(\omega_0-\Omega)L = 2\beta_0L + \beta_2\Omega^2L + \ldots$ が乗る。
- $2\beta_0L$ はキャリア（＝同じファイバを通る LO）と共通。奇数次（$\beta_1$, $\beta_3$）は時間遅れで、quadrature を回さない。
- 残る **スクイーズ角の回転** は

$$\theta(\Omega) = \frac{\beta_2\,\Omega^2 L}{2}$$

| 測定周波数 $\Omega/2\pi$ | 1 km での $\lvert\theta\rvert$ | 1 m での $\lvert\theta\rvert$ |
|---|---|---|
| 10 GHz | 0.043 rad | — |
| 43 GHz | 0.80 rad | 0.8 mrad |
| 1 THz | — | 0.43 rad |

（SMF, $|\beta_2| = 22$ ps²/km で計算）

- **意味**: ロックで合わせられるのは 1 つの $\Omega$ の角度だけ。広帯域（THz 級）の導波路 OPA の光は、**1 m の SMF でも帯域の端でスクイーズ角が大きく回る**。単色の CW LO では補償できないので、分散補償ファイバ（損失が増える）や分散を逆にした光路が必要 (agent)。これは RQ2 の中心。

**XPM（明るい LO/probe と同じファイバ）** (agent):

- 同じ偏光の明るい場（パワー $P$）から、スクイーズド光は $2\gamma PL$ の XPM 位相を受ける。LO 自身は $\gamma PL$ の SPM を受けるので、相対的なずれは $\gamma PL$。直交偏光なら XPM の係数は 2 ではなく 2/3。
- 明るい場の強度雑音 $\delta P$ は、そのまま位相雑音 $\gamma L\,\delta P$ になる。

**GAWBS・Raman**: 3.4 節の PSD を、測定帯の $\Omega$ ごとの過剰雑音として $V(\Omega)$ に足す。

**レジーム P との整合**: 長いパルス（$L_D \gg L$）の極限で、レジーム P の時間グリッドの結果を周波数に直したものが、レジーム C と一致すること（テスト）。

### 3.6 制御層（L3）

前回のプロトタイプ（`assets/2026-10-07_phase_lock_budget.py`）をパッケージに取り込む。

- **入力**: レジーム P/C から得た位相雑音 PSD の成分（熱・音響・GAWBS・XPM 由来）と、誤差信号の傾き・センサ雑音
- **処理**: ループ（Moku の遅延、EOM/AOM、オフロード）で $|S|^2S_d + |T|^2S_n$ を計算
- **出力**: $\theta_\mathrm{rms}$ と、$\langle V\rangle = V_- + (V_+-V_-)(1-e^{-2\sigma^2})/2$ による測定スクイージング
- **ロックできる/できない雑音の分け方** (agent):
  - ループ帯域（〜100 kHz）より十分下の雑音 → $|S|^2$ で抑圧される
  - ループ帯域付近 → サーボバンプで増幅されることもある
  - 測定帯（MHz–GHz）にある雑音（GAWBS など）→ ループでは消せず、過剰雑音として直接スクイージングを劣化させる
  - **パルスごとの位相揺らぎ**: 繰り返し周波数より遅い成分はループで、速い成分はパルスごとのランダムな回転として扱う

---

## 4. ソフトウェア設計

### 4.1 技術選定

| 項目 | 選択 | 理由 (agent) |
|---|---|---|
| 数値計算 | **JAX**（`jit`, `lax.scan`, `vmap`, `jacfwd`/`jvp`/`vjp`） | 自動微分で線形化がただで得られる。TWA の並列化も `vmap` で済む |
| 精度 | `jax_enable_x64 = True`（必須） | float32 だとシンプレクティック誤差が大きく、−15 dB 級のスクイージングが壊れる |
| 代替 | PyTorch（`torch.func.jacfwd`, `vmap`） | BoTorch で慣れているなら選択肢。ただし `scan` 相当の JIT は JAX の方が素直 |
| 古典ソルバの検算 | gnlse-python（Raman モデル、分散入力つき）(paper: arXiv:2110.00298)、PyNLO | 自作 RK4IP と同じ入力でパルス形状・スペクトルを比べる |
| ガウス状態の補助 | thewalrus（共分散、シンプレクティック行列の道具）| GBS 側の道具と共通化できる |
| 制御 | numpy（＋ python-control） | 前回のプロトタイプをそのまま使う |
| 光学テーブル側（CW） | **Finesse 3** | 共振器・BHD・RF 復調・モード不整合（HG モード）・量子雑音スペクトルを周波数領域で。4.4 節 |
| 設定・記録 | dataclass ＋ YAML、結果は `.npz`（大きければ zarr） | 図ごとに設定ファイルと git のコミットを残す |
| テスト | pytest（遅いテストは `@pytest.mark.slow`） | 合格基準をそのままテストにする |

### 4.2 パッケージ構成

```text
fibersqz/
  conventions.py   # 単位（ps, m, W）, FFT の符号, quadrature 規約（真空 = 1）を 1 か所に
  grid.py          # 時間・周波数グリッド, 窓の自動決定（L_D, L_NL, 走り去りから）
  fiber.py         # Fiber dataclass: α, β_k, γ, Raman モデル, 偏光, 長さ
  raman.py         # h_R(t): Blow–Wood, Lin–Agrawal; 周波数応答; 熱雑音の相関
  gnlse.py         # RK4IP（単一場 / 2 場結合: 偏光・周波数多重）, lax.scan, 刻み幅の収束判定
  linearize.py     # 実表現の写像, S（jacfwd）, 少数モード（jvp）, 測定モード（vjp）
  covariance.py    # 共分散の伝搬, 区間ごとの損失, 拡散項（Raman, GAWBS）, シンプレクティック検査
  twa.py           # 真空雑音の注入, Raman 雑音, アンサンブル（vmap, 分割実行）
  sources.py       # 入力状態: 単一モードのスクイーズド真空（時間モード指定）, 多モード, コヒーレント
  sideband.py      # レジーム C: Ω ごとの 2×2, 分散回転, XPM, GAWBS 線, Raman
  optics_finesse.py # Finesse 3 のモデル構築（受信側の光学系, CCSB の誤差信号, 光学プラント）と結果の受け渡し
  noise_psd.py     # 位相雑音 PSD: 熱的基本限界, GAWBS, 環境（仮置き/実測読み込み）
  measure.py       # ホモダイン射影, LO モード, 最適 θ, Bloch–Messiah, 検出効率・電子雑音
  control.py       # 前回プロトタイプ: ループ, 余裕, 残留 PSD, 時間領域, オフロード
  experiments/     # E1–E5 のスクリプト（設定 YAML と対）
tests/
  test_linear.py, test_gnlse.py, test_linearize.py, test_twa.py, test_sideband.py, test_control.py
notebooks/         # 図の作成と探索（ロジックは置かない）
```

### 4.3 主要な関数の形（案）

```python
@dataclass(frozen=True)
class Fiber:
    length: float              # m
    alpha: float               # 1/m（パワー）
    betas: tuple[float, ...]   # (β2, β3, ...) [ps^k/m]
    gamma: float               # 1/(W m)
    fR: float = 0.18
    raman: str = "blow_wood"   # "none" | "blow_wood" | "lin_agrawal"

def propagate(A0: Array, fiber: Fiber, grid: Grid, nsteps: int) -> Array: ...
def transfer_matrix(A0, fiber, grid, nsteps) -> Array:            # 2N×2N, jacfwd
def measure_rows(A0, fiber, grid, nsteps, W: Array) -> Array:     # S^T W, vjp
def propagate_cov(V, A0, fiber, grid, nseg: int) -> Array:        # 損失を区間ごとに
def homodyne_var(V, u, theta, grid) -> float
def bloch_messiah(S) -> tuple[Array, Array, Array]
def twa_ensemble(A0, V_in, fiber, grid, n_traj, key) -> Array    # (n_traj, N) complex
def sideband_cov(V0, Omega, fiber, P_bright, psd_excess) -> Array # (len(Omega), 2, 2)
```

- **設計の約束** (agent): 物理の関数はすべて純粋関数（JAX の `jit` に乗る）。乱数は `jax.random` の key を明示的に渡す。単位の変換は `conventions.py` 以外でしない。

### 4.4 Finesse 3 の位置づけ

Finesse 3 は「CW の光学系を周波数領域で解く」道具。ファイバの伝搬（分散・非線形・パルス）は自作コードで、ファイバの外の光学テーブル側は Finesse で、と役割を分ける (agent)。

**Finesse に任せるもの**:

- **受信側の光学系**: LO とスクイーズド光のビームスプリッタ、BHD（Balanced Homodyne Readout）、検出効率・損失
  - **意味**: 自作の `sideband.py` と同じ量を、独立に作られたコードで計算できる。M8b の相互検算になる。
- **モード不整合**: Gaussian beam（HG モード）を使った LO とスクイーズド光の空間モードの重なり
  - **意味**: ファイバ出力のコリメートや自由空間部分の不整合を、効率の数字ではなく光学系から求められる。
- **CCSB の誤差信号と光学プラント**: AOM で周波数をずらした probe（周波数オフセットした Laser）、EOM（Modulator）、RF 復調（RF photodetector readout）
  - **意味**: 誤差信号の形と傾き、アクチュエータ（EOM の位相）から誤差信号までの伝達関数（signal を使った周波数応答）を得る。これを `control.py` のプラントとして渡す。
- **共振器を使うスクイーザー**: OPO の共振器の応答（サイドバンドの周波数依存）、フィルタ共振器（LO の ASE 除去, [[Hirota2025_squeezing10db]] 型）
- **Lock**: ロックの動作点（DC）を決める。例: 復調位相のずれを掃引したときのスクイージングの変化

**Finesse に任せないもの**（自作コードで扱う）:

- **ファイバの分散・非線形・パルス**: Space は長さと屈折率で位相を与えるだけで、分散や Kerr は扱わない（文書に該当する要素が見当たらない）(web; agent)。したがって θ(Ω) = β₂Ω²L/2 の回転、XPM、Raman、GAWBS は出てこない。
- **パラメトリック利得そのもの**: Squeezer は「指定した dB と角度のスクイーズド光を出す理想的な光源」。ポンプと probe の相対位相で probe が増幅・減衰する過程（CCSB の loop 1 の誤差信号の元）は、要素一覧に非線形結晶・OPA が見当たらないので、そのままでは扱えない可能性が高い（要確認）(web; agent)。
- **ループの動特性**: Lock は動作点を合わせるためのもので、遅延・飽和・時間応答は扱わない。ループは `control.py` で (agent)。

**受け渡しの形** (agent):

```text
自作（L1/L2, sideband.py） ──→ ファイバの実効パラメータ: η_fiber, 回転 θ(Ω), 過剰雑音 PSD, 位相雑音 PSD
                                     │
                                     ▼
Finesse 3（optics_finesse.py）: ファイバ = Space（長さ）＋ 損失（BS/Mirror の透過）＋ 位相雑音（signal 注入）
   └→ BHD スペクトル, 誤差信号の傾き K_d, 光学プラント P(f), モード不整合 η_mm
                                     │
                                     ▼
control.py: L = K_d C A e^{-sτ}, 残留位相 → 測定スクイージング
```

- θ(Ω) のような周波数依存の回転は Finesse にそのまま入らないので、比較は「ロック点 Ω_lock 付近」か「回転なしの条件」で行う。周波数依存の部分は `sideband.py` を正とする (agent)。

**他の光学系ライブラリ** (agent):

- **thewalrus / MrMustard**: 離散モードのガウス回路（ビームスプリッタ、損失、ホモダイン）。多モード（スーパーモード）や GBS 側との接続に。
- **QuTiP**: 非ガウス状態を扱うときに。
- **gnlse-python / PyNLO**: 自作 GNLSE の検算（M2）。
- 光線追跡・回折計算のライブラリ（自由空間のビーム伝搬など）は、この研究の問いには不要。

---

## 5. マイルストーン（合格基準つき）

期間はあくまで目安（週 10–15 時間を想定）(agent)。

| # | 内容 | 合格基準（テスト） | 目安 |
|---|---|---|---|
| M0 | 環境・規約 | `conventions.py` と README に単位・FFT・quadrature の規約。x64 が有効なことを検査 | 0.5 週 |
| M1 | 線形伝搬（分散・損失） | ガウスパルスの広がり $T_1 = T_0\sqrt{1+(z/L_D)^2}$、エネルギー $e^{-\alpha z}$ に一致 | 0.5 週 |
| M2 | GNLSE（SPM, Raman, 自己急峻化） | 1 次ソリトンの形が保存（誤差 < 1e-6）、2 次ソリトンが周期 $z_0 = \pi L_D/2$ で戻る、刻み幅半分で結果が収束、gnlse-python と同条件で一致 | 1–1.5 週 |
| M3 | 線形化・共分散 | 損失なしで $\lvert S^T\Omega S-\Omega\rvert < 10^{-8}$、ヤコビアンが有限差分と一致、無分散の Kerr で解析式（下記）に一致、JVP/VJP 方式と jacfwd 方式が一致 | 1 週（**種コードで一部済**） |
| M4 | 損失・検出 | 区間数を倍にしても結果が変わらない、効率 $\eta$ の式に一致、LO が別経路のときのモード重なり | 0.5 週 |
| M5 | 文献の再現 | ソリトン・スクイージングの長さ依存と最適 LO の傾向を Carter 1987 / Rosenbluh & Shelby 1991 と比べる（定性） | 1 週 |
| M6 | TWA | 雑音が小さい極限で M3 と一致（統計誤差内）、$M$ を増やすと $1/\sqrt M$ で収束 | 1–1.5 週 |
| M7 | 雑音源（Raman, GAWBS, 熱） | Raman 熱雑音を入れたときの過剰雑音が Drummond & Corney の傾向と合う、GAWBS の線の位置が Shelby 1985 と合う | 1 週 |
| M8 | レジーム C | 長パルス極限でレジーム P と一致、$\theta(\Omega) = \beta_2\Omega^2L/2$ を数値で再現 | 1 週 |
| M8b | Finesse 3 の光学テーブルモデル | ファイバを「損失＋位相」に置き換えた条件で、Finesse の BHD スペクトルが `sideband.py` と一致（0.05 dB 以内）、誤差信号の傾きが [[2026-10-07 phase-stab-sim_plan-review]] 4.2 節の式と一致 | 1 週 |
| M9 | 制御層の統合 | 前回のプロトタイプの検算（時間領域と周波数領域が 1 % 以内）をパッケージ内で再現 | 1 週 |
| M10 | 計算実験 E1–E5 とまとめ | 6 節の図がすべて設定ファイルから再生成できる | 2 週 |

**無分散 Kerr の解析式**（M3 の主テスト）(agent; 導出):

- 無分散なら各時刻のスライスが独立に $\Phi = \gamma|A|^2L$ の自己位相変調を受ける。
- 平均場の位相を基準にした振幅・位相 quadrature $(X, Y)$ は $X\to X$, $Y\to Y + 2\Phi X$ と変換される。真空入力なら

$$V = \begin{pmatrix}1 & 2\Phi\\ 2\Phi & 1+4\Phi^2\end{pmatrix},\qquad V_\mathrm{min,max} = 1+2\Phi^2 \mp 2\Phi\sqrt{1+\Phi^2}$$

- $\Phi$: 非線形位相 [rad]。Kitagawa & Yamamoto 1986 の Kerr スクイージングと同じ形。

---

## 6. 計算実験（研究の中身）

### E1: スクイーズド真空のパルスだけがファイバを通る（RQ1）

- **設定**: 入力 = ガウス時間モードの単一モードスクイーズド真空（例: 10 dB）。ファイバは線形（明るい場なし）。LO は (a) 同じ長さの別ファイバ、(b) ファイバなし（分散なし）、(c) 理想（Bloch–Messiah）。
- **掃引**: $L$（1 m–10 km）、パルス幅（100 fs–10 ps）
- **出力**: 測定スクイージング vs $L/L_D$。「LO が同じ分散を受けるなら劣化は損失だけ」「受けないなら $L\sim L_D$ でモード不整合が支配的」を定量化 (agent の予想)。
- **研究の価値**: パルス幅ごとの「使えるファイバ長」の設計則。

### E2: 明るい LO / probe と同じファイバを通る（RQ1, RQ4）

- **設定**: 2 場結合の GNLSE（偏光多重、または周波数多重）。明るい場のパワー、偏光（平行/直交）、周波数オフセットを変える。
- **出力**: XPM によるパルス内で不均一な楕円の回転（単一の LO 位相で直せない分）、LO の強度雑音由来の位相雑音、Raman による過剰雑音を、それぞれ分けて dB で示す。
- **研究の価値**: 「同一ファイバで送ると差動位相雑音は減るが、非線形の代償はどこから効くか」の境界。

### E3: 広帯域 CW スクイージングの分散回転（RQ2）

- **設定**: レジーム C。OPA の出力スペクトル（平坦 or 有限帯域）、$L$、分散補償の有無。
- **出力**: $\Omega$ ごとの測定スクイージング（ロック点 $\Omega_\mathrm{lock}$ を変えて）。分散補償ファイバの損失との兼ね合い。
- **研究の価値**: THz 級の導波路 OPA（[[Hirota2025_squeezing10db]] 型）の光をファイバで送るときの上限。

### E4: ロックできる雑音・できない雑音（RQ3）

- **設定**: 位相雑音 PSD = 熱的基本限界 ＋ 環境（仮置き → 実測）＋ GAWBS ＋ XPM 由来。前回の制御層（Moku の遅延、EOM/AOM）。
- **出力**: 周波数帯ごとの寄与の積み上げ図、$L$ に対する $\theta_\mathrm{rms}$ と測定スクイージング、受信側ロックと送信側 FNC の比較。
- **研究の価値**: Moku 実験に直接返せる結果（必要なループ帯域、アクチュエータの組み合わせ）。

### E5: ファイバ中の Kerr スクイージング（検証を兼ねる）

- **設定**: コヒーレントなパルス（ソリトン含む）を入れ、出力のスクイージングと最適 LO を計算。
- **出力**: 文献の傾向の再現。種コードでは、2 $L_D$ 伝搬した 1 次ソリトンで最良のスーパーモード −11.9 dB に対し、LO をパルス形状にすると −7.1 dB（モード不整合の大きさの例）。
- **研究の価値**: ソルバ全体の信頼性の担保。

---

## 7. 計算コストの見積もり (agent)

- 種コード（コンテナの CPU）: $N=512$, 200 ステップの jacfwd（1024 列）と、$N=256$, 2000 ステップの jacfwd（512 列）、合わせて約 23 秒。
- $S$ 全体のコストは「$2N$ × 伝搬 1 回」。$N = 2048$ なら $N=256$ の約 64 倍以上（FFT の $\log N$ も効く）。メモリは $S$ が 4096×4096 の倍精度で約 134 MB。
- **重くなったら**:
  - 測定モードだけ欲しいなら VJP 方式（$M$ 本）に切り替える
  - 入力が少数モードで明るい場がないなら JVP 方式（$K$ 本）
  - jacfwd は列を分割して `lax.map` で回し、メモリを抑える
- TWA: 1 軌道は伝搬 1 回。$M = 10^4$ 本でも `vmap`（分割）で現実的。

---

## 8. つまずきやすい点 (agent)

- **FFT の符号**: numpy/JAX の順方向 FFT は $e^{-i\omega t}$。$\partial_T \to -i\omega$。分散項の符号を取り違えるとソリトンが崩れる（M2 のテストで検出できる）。
- **真空雑音の正規化**: 複素サンプルあたり $\langle|\delta a_j|^2\rangle = 1/(2\Delta t)$。ホモダイン射影で真空が 1 になることを毎回テストする。
- **平均場の位相**: SPM で平均場の位相が回るので、LO の位相は出力の平均場を基準に取る（種コードでは $u \propto e^{i\arg A_\mathrm{out}}$）。
- **時間窓とエイリアシング**: Raman の畳み込みは FFT で行うので、窓は $\tau_2$ よりずっと長くし、端がゼロに落ちていることを確認する。
- **刻み幅**: 固定刻みなら半分にして収束を確かめる。適応刻みは JAX の `scan` と相性が悪いので、まず固定で。
- **損失の扱い**: 区間ごとの損失は近似。区間数を倍にして確かめる。
- **TWA の適用範囲**: モードあたりの光子数が十分多いときに有効。真空に近いモードが支配的な量は線形化の方が信頼できる。
- **GAWBS を時間グリッドに入れない**: 2 節の理由。
- **dB の基準**: shot noise（真空 = 1）を 0 dB に統一。電子雑音の補正をする場合は別の列で。

---

## 9. 参考文献（モジュール別）

### 古典場ソルバ（L1）

- **Agrawal, *Nonlinear Fiber Optics***（教科書）— GNLSE、分散長・非線形長、ソリトン。
- **Hult, J. Lightwave Technol. 25, 3770 (2007)** — RK4IP。[Cambridge のページ](https://laser.ceb.cam.ac.uk/node/121)
- **Dudley, Genty, Coen, Rev. Mod. Phys. 78, 1135 (2006)** — スーパーコンティニュームの総説。GNLSE の数値解法と入力の量子雑音（1 モードあたり 1 光子のモデル）。
- **Blow & Wood, IEEE J. Quantum Electron. 25, 2665 (1989)** — Raman 応答 $h_R$。
- **Lin & Agrawal, Opt. Lett. 31, 3086 (2006)** — 改良した Raman 応答。
- **gnlse-python, arXiv:2110.00298** — 検算用の既存ソルバ。[ドキュメント](https://gnlse.readthedocs.io/en/stable/gnlse_intro.html)
- **arXiv:2206.07526** — 入力パルスの量子雑音モデル（スペクトル位相雑音、ガウス雑音）を実装した GNLSE の Python モジュール。古典的な背景場として雑音を模擬する (paper; 要旨より)。TWA の入力雑音の実装の参考。

### 量子揺らぎ（L2）

- **Carter, Drummond, Reid, Shelby, PRL 58, 1841 (1987)** — ソリトンのスクイージング（量子化した NLSE）。
- **Lai & Haus, PRA 40, 844 / 854 (1989)**、**Haus & Lai, JOSA B 7, 386 (1990)** — 線形化によるソリトンの量子論（書誌は要確認）。
- **Rosenbluh & Shelby, PRL 66, 153 (1991)** — ソリトン・スクイージングの実験。
- **Drummond & Corney, JOSA B 18, 139 / 153 (2001)**, arXiv:quant-ph/9912095, 9912096 — Raman と GAWBS を含む確率微分方程式。TWA と雑音の相関の式の出典。
- **Kitagawa & Yamamoto, PRA 34, 3974 (1986)** — Kerr 媒質による数・位相スクイージング。
- **Braunstein, PRA 71, 055801 (2005)** — Bloch–Messiah 分解。
- **Wasilewski, Lvovsky et al., PRA 73, 063819 (2006)** — パルスのスクイージングの多モード性。
- **Weedbrook et al., RMP 84, 621 (2012)** — ガウス状態と共分散の規約。

### ファイバの雑音源

- **Shelby, Levenson, Bayer, PRB 31, 5244 (1985)** — GAWBS の理論と測定。モード周波数の式の出典。
- **Elser et al., arXiv:quant-ph/0512044** — フォトニック結晶ファイバによる GAWBS の低減。通常ファイバとの比較。
- **Liu, Iu, Hamdash, Helmy, arXiv:2406.19991** — 全ファイバのスクイーズド光（7.5 dB）。GAWBS と Raman を測定帯の選び方で避ける (paper; 要旨)。
- **Wanser, Electron. Lett. 28, 53 (1992)**（書誌は要確認）、**Duan, PRA 86, 023817 (2012)** — ファイバの熱的位相雑音の基本限界。
- **arXiv:1512.02799** — ファイバの基本熱雑音を低周波まで観測。

### 位相制御・ファイバ伝送（L3）

- 前回の Discussion の文献一覧（[[2026-10-07 phase-stab-sim_plan-review]] 5 節）: Vahlbruch 2006（CCSB）、Williams 2008（FNC の遅延限界）、[[Verclas2026_distribution]]、[[Hirota2025_squeezing10db]] など。

---

## 10. 種コードの結果（動作確認済み）

`assets/2026-10-07_nlse_linearized_seed.py`（JAX 0.11, x64, CPU）。合計約 23 秒。

| テスト | 結果 |
|---|---|
| (a) 無分散 Kerr、$\Phi = 1.3$ rad（1 W, 1 km） | 数値 −9.367 / +9.367 dB、解析 −9.367 / +9.367 dB（一致） |
| (b) シンプレクティック性（分散＋Kerr, $N=256$） | $\max\lvert S^T\Omega S-\Omega\rvert = 5\times10^{-10}$ |
| (c) 1 次ソリトン（2 $L_D$ ≈ 91 m 伝搬） | $\lvert A\rvert$ の最大相対誤差 $4\times10^{-9}$ |
| (d) Bloch–Messiah（同じソリトン） | 最良スーパーモード −11.85 dB / +11.85 dB、LO = 出力パルス形状では −7.14 dB / +11.21 dB |

- **読み取れること** (agent):
  - JAX の jacfwd で得た行列が、解析式・シンプレクティック性の両面で正しい。計画の L1 + L2 の土台は成り立つ。
  - (d) の差（約 4.7 dB）が、パルス・スクイージングで「LO の時間モード」がどれだけ重要かを示す。E1・E5 の出発点になる。
- **次にやること**: この種コードを `gnlse.py` と `linearize.py` に分け、M1–M3 のテストを pytest に移す。

---

## Open questions

- 対象を CW とパルスのどちらに重心を置くか（両方作るが、E1–E5 の優先順位が変わる）。
- パルスの場合の具体的なパラメータ（中心波長、パルス幅、繰り返し、スクイーザーの種類）。
- GAWBS の係数 $c_m$ と線幅 $\Gamma_m$ を、どの文献の測定値に合わせるか（SMF-28 の実測値を探す）。
- Raman 熱雑音の相関の正確な形（Drummond & Corney の式を原典で確認してから実装）。
- 2 場結合（偏光多重）で、複屈折のランダムさ（偏光モード分散）をどこまで入れるか。
- スクイーザー自体（χ(2) の導波路 OPA）も同じ枠組み（線形化）で作るか。作れば「生成から測定まで」を一貫して扱える。

## Related concepts

- [[Squeezed States]]
- [[Homodyne Detection]]
- [[Phase Space Representation]]
- [[Coherent Probe Locking]]
- [[Noise Locking]]
- [[Shot-Noise Limit]]
- [[Active Phase Stabilization in Photonic QC]]

## Related papers

- [[Hirota2025_squeezing10db]]
- [[Verclas2026_distribution]]
- [[Vahlbruch2016_squeezing15db]]

## Related discussions

- [[2026-10-07 phase-stab-sim_plan-review]]
