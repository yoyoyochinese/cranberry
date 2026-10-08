---
type: discussion
project: gleamworks
date: 2026-10-07
status: draft
tags: [gleamworks, discussion, simulation, phase-stabilization, squeezed-vacuum, fiber, Moku, AOM, EOM, homodyne-detection, control]
related:
  - "[[Verclas2026_distribution]]"
  - "[[Hirota2025_squeezing10db]]"
  - "[[Vahlbruch2016_squeezing15db]]"
  - "[[McKenzie2005_noiselocking]]"
  - "[[Coherent Probe Locking]]"
  - "[[Noise Locking]]"
  - "[[Homodyne Detection]]"
  - "[[Squeezed States]]"
  - "[[PID Control]]"
  - "[[Active Phase Stabilization in Photonic QC]]"
---

# Discussion：位相安定化シミュレーション計画のレビューと再設計（Moku・高 QE PD・AOM/EOM）

関連: [[2026-09-24 LiquidInst2025_pll]]（開ループ伝達関数の扱い）, [[2026-10-05 Hirota2025_squeezing10db]]（probe 方式）, [[2026-09-28 Verclas2026_distribution]]（ファイバ伝送＋ noise lock）

出典の表記: `(paper)` 論文に書いてある、`(manual/web)` メーカーの資料・KB に書いてある、`(agent)` 私（Claude）の推論・概算（式や計算で確かめられる）、`(me)` 自分の判断。

生成スクリプト: `assets/2026-10-07_phase_lock_budget.py`（numpy と matplotlib だけで動く）

---

## 0. この Discussion の位置づけ

- **入力**: 「パルススクイーズド光の位相安定化シミュレーション計画」（plan.md）。ファイバ効果（分散・Kerr・Raman・GAWBS）、Finesse 3、python-control の 3 ツール連携を想定したもの。
- **方針の変更 (me)**: 現在の DSC1 / 915B の系はいったん忘れ、**Moku・高量子効率 PD・AOM/EOM** を使った位相安定化実験を前提にシミュレーションを組み直す。
- **この Discussion でやったこと (agent)**:
  1. 計画書の評価（良い点・足りない点）
  2. 必要な理論の整理と、参考にできる論文・資料の調査
  3. 具体的な実装案（モジュール構成・段階的な進め方・検証テスト・BO との接続）
  4. 最小プロトタイプを実装し、ファイバ長に対する残留位相とスクイージングの傾向を確認
- **置いた前提 (agent)**: 特に断らない限り **CW, 1550 nm の squeezed vacuum** を BHD で測る系を考える。パルスの場合は 4.7 節で分けて扱う。

---

## 1. 結論の要約

- 計画書の **物理の押さえどころ（平均場がない → probe/CCSB が必要、EOM のレンジは有限、同一ファイバ伝搬の利点）は正しい** (agent)。
- 一方で、重心が「ファイバ中の非線形な量子雑音伝搬」に寄りすぎていて、**位相ロックの性能を実際に決める要素**（ループ内の遅延、誤差信号の雑音、アクチュエータの種類とレンジ、Moku のデジタル実装）が「zpk で足す」の一行で済まされている (agent)。
- 新しい目的（Moku + AOM/EOM の位相安定化）に合わせるなら、**優先順位を逆転**させるのがよい (agent):
  1. まず「遅延＋雑音＋アクチュエータ」の制御モデル（周波数領域＋時間領域）
  2. 次に CCSB の 2 ループ（probe↔pump、LO↔probe）
  3. ファイバの非線形・パルスの扱いは「η_eff・過剰雑音・位相雑音 PSD の数個の数値を返すモジュール」に縮め、必要になってから作る
- プロトタイプで確認できたこと (agent, 7 節):
  - **ロックをどこで閉じるか** が一番効く。受信側の LO を動かす構成ではファイバの伝搬遅延がループに入らず、送信側 FNC（往復光で検出）より 1 km で 1 桁以上ジッタが小さい（仮置きの雑音モデルで）。
  - EOM 単独ではドリフトでレンジ（±π）を使い切って 2π スリップが起きる。AOM（周波数＝位相の積分）への **オフロード** で解消する。
  - 時間領域シミュレーションと周波数領域の予測が 1 % で一致（27.5 vs 27.2 mrad）。モデルの自己整合は確認済み。

---

## 2. 計画書の評価

### 2.1 良い点

- **平均場がないことへの対処**: シード、CCSB（AOM で周波数をずらした制御光）、ノイズロックの 3 択を挙げ、CCSB では 2 ループが要ると正しく書いている (agent; Vahlbruch 2006 の構成と一致)。
- **EOM のレンジ制限と 2 段ループ**: 長いファイバの熱ドリフトは数 rad 以上になり EOM は飽和する、という指摘は正しい。[[Verclas2026_distribution]] も EOM＋ピエゾ＋ファイバストレッチャの 3 段で帯域を分担している (paper)。
- **同一ファイバ伝搬の利点**: LO と信号が別ファイバだと差動の位相雑音が支配的になる、という指摘は正しい (agent; Williams 2008 でも往路と復路を別ファイバにすると差動雑音が支配的になると述べている (paper))。
- **Finesse の CW 限界の認識**: サーボ帯域は繰り返し周波数よりずっと低いので、ロックはパルス平均量しか見ない、という整理は妥当 (agent)。
- **位相ジッタからスクイージングへの換算式** $V \approx V_- + (V_+-V_-)\theta_\mathrm{rms}^2$ は小角近似として正しい (agent; Aoki 2006 の式の展開)。

### 2.2 問題点・足りない点

#### (1) 目的と労力の配分が合っていない

- **内容**: ステップ 1（NLSE＋JAX の線形化、Raman、GAWBS、truncated Wigner）は、それだけで 1 本の研究になる規模。一方、新しい目的は「Moku と AOM/EOM で位相を安定化する実験」。
- **なぜ問題か**: CW の squeezed vacuum と µW 級の probe なら、ファイバ中の非線形効果は無視できる（4.7 節の概算で Kerr 位相は $10^{-6}$ rad 程度）。重い非線形ソルバーを作っても、ロック性能の予測はほとんど変わらない (agent)。
- **提案**: ファイバは「損失＋位相雑音 PSD＋伝搬遅延」の線形モデルから始める。非線形・パルスのモジュールは、(a) パルス動作に決めた、(b) 明るい LO を同じファイバに通す、のどちらかになった時点で作る (agent)。

#### (2) ループ内の遅延がモデルにない

- **内容**: 計画では「Finesse で P(f) を出し、EOM と回路を zpk で足す」。しかし Finesse の伝達関数には、Moku の処理遅延、AOM の音響遅延、ファイバの伝搬遅延が入らない (agent)。
- **なぜ問題か**: 位相ロックの帯域の上限は、ほとんどの場合この遅延で決まる。$L \propto 1/s$ のループでは位相余裕 45° に対して UGF の上限はおよそ $1/(8\tau)$（4.4 節）。遅延の見積もりがないと、必要な帯域が取れるかを判断できない (agent)。
- **提案**: 遅延は周波数領域で $e^{-s\tau}$ をそのまま掛け、時間領域では整数サンプルの遅れとして入れる（Padé 近似は不要）(agent)。

#### (3) 誤差信号の雑音（センサ雑音）が扱われていない

- **内容**: 残留位相は $|S|^2 S_d + |T|^2 S_n$ で、ループゲインを上げるほど 2 項目のセンサ雑音 $S_n$ がそのまま残る。計画には $S_n$ の項がない (agent)。
- **なぜ問題か**: probe 方式では probe の shot noise、ノイズロック方式ではスクイーズド光そのものの雑音が $S_n$ を決める。「高 QE の PD を使う」効果がどこに効くかも、$S_n$ と効率 $\eta$ を分けないと評価できない（4.2 節）(agent)。

#### (4) Finesse の役割が大きすぎる

- **内容**: Finesse 3 の Lock は動作点（DC）を合わせるためのもので、ループの過渡応答や遅延・飽和は扱わない (agent; Finesse の文書は周波数領域の解析が中心で、時間領域の解析ページは見当たらない (web))。
- **提案**: ループの評価は python-control / numpy に一本化する。Finesse は「OPO 共振器を含む CCSB の誤差信号の形と傾き」を出すときの検算用に限定する。単一パスの導波路 OPA なら解析式で足りる (agent)。

#### (5) 制御光（probe）と信号の「周波数差」による位相のずれ

- **内容**: CCSB の probe はスクイーズド光より $\Omega$（AOM のシフト量）だけ周波数が高い。同じファイバを通っても、probe とキャリアの位相差には $\Omega\,\tau_g$（$\tau_g$ は群遅延）が残る (agent)。
- **意味**: 温度で $\tau_g$ が変わると、probe にロックしていてもスクイーズ角がゆっくりずれる。概算で $\Omega/2\pi = 40$ MHz, 1 km, $\delta\tau_g/\tau_g \approx 7\times10^{-6}$/K とすると、約 9 mrad/K（4.2 節）(agent)。
- **実用上の目的**: $\Omega$ をループ帯域より十分上で、かつ必要以上に大きくしない。長距離では温度依存のオフセットを別途補正するか、測定帯のスクイージング自体でゆっくり補正する（ノイズロックの併用）(agent)。

#### (6) その他

- **換算式の適用範囲**: $\theta_\mathrm{rms}$ が 0.1 rad を超えるあたりから小角近似が崩れる。ガウス揺らぎの厳密な平均を使う（4.1 節）(agent)。
- **QuTiP・Tidy3D・XMDS2**: ガウス状態なら 2×2（多モードなら 2N×2N）の共分散行列で足りる。Tidy3D・XMDS2 は SMF を使う限り不要 (agent)。
- **検証の仕組みがない**: 既知の結果（Williams 2008 の遅延限界、Aoki 2006 の式、時間領域と周波数領域の一致）をテストとして持つべき（6.4 節）(agent)。

---

## 3. シミュレーション対象とする実験系の案

### 3.1 構成 A：受信側ロック（推奨の出発点）

```text
[送信側]
 laser ─┬─ SHG → pump → squeezer (OPA/OPO)
        ├─ AOM (+Ω) → EOM₁ → probe ─┘ (pump と一緒に squeezer へ。loop 1)
        └──────────────── LO 用 ─────────────┐
 squeezed vacuum + probe ──[ファイバ L]──→   │  (LO を別ファイバで送るか、同じファイバに偏光多重)
[受信側]                                      ↓
                                  LO → EOM₂ (速い) → AOM/ストレッチャ (遅い) → BHD
 BHD 出力 → (Ω で復調: Moku LLB / Lock-in) → PID (fast) → EOM₂
                                           └→ PID (slow, offload) → AOM 周波数 or ストレッチャ   (loop 2)
```

- **動作**: probe はスクイーズド光と同じファイバを通り、ファイバの位相雑音を同じように受ける。受信側で probe と LO のビート（周波数 $\Omega$）を復調すると、LO とスクイーズド光の相対位相に比例する誤差信号が出る (agent; [[Coherent Probe Locking]])。
- **意味**: アクチュエータ（EOM₂）が受信側にあるので、**ファイバの伝搬遅延はループに入らない**。ファイバは「外乱」としてだけ現れる (agent)。
- **注意**:
  - LO を別ファイバで送ると、そのファイバの雑音も外乱に加わる（独立なら PSD が約 2 倍、ジッタは約 √2 倍）(agent)。
  - 1 台の Moku で AOM の駆動と復調を両方行うには、送信側と受信側が同じ机にある（ファイバスプールで模擬する）ことが前提。離れた建物なら RF 基準の伝送が要る (agent)。

### 3.2 構成 B：送信側 FNC ＋ 受信側の短いループ

```text
 laser → AOM (FNC 補正) → [ファイバ L] → 受信側で一部を反射 (FRM) → [ファイバ L] → 送信側で検出
```

- **動作**: 往復した光と手元の光を干渉させて往復の位相を測り、入口の AOM で補正する（Ma 1994 型の FNC）(paper: Ma 1994)。
- **意味**: 往復の遅延 $2\tau$ がループに入るので、UGF の上限が下がる。さらに、遅延のせいで原理的に打ち消せない残差 $\approx \tfrac{1}{3}(2\pi f\tau)^2 S_\mathrm{fiber}$ が残る（4.4 節）(paper: Williams 2008)。
- **実用上の目的**: LO を受信側で用意できない（受信側に光源しかない、など）ときに使う。今回のように同じ机で閉じられるなら構成 A の方が有利 (agent)。

### 3.3 構成 C：ノイズロック（probe なし）

- **動作**: BHD 出力の雑音の大きさ（エンベロープ）をディザで変調・復調して誤差信号にする (paper: [[McKenzie2005_noiselocking]], [[Verclas2026_distribution]])。
- **意味**: probe 由来の余分な光がない代わりに、誤差信号の SNR がスクイージングの「見え方」で決まり、ディザがそのまま位相雑音になる（Verclas ではディザ ±4.5° がフィットの位相雑音 約 0.1 rad の主因 (paper)）。
- **実用上の目的**: 構成 A の補助（遅いオフセット補正）か、probe を入れたくない場合の比較対象 (agent)。

### 3.4 Moku の割り当て（構成 A の例）

| 役割 | Moku の機能 | メモ |
|---|---|---|
| probe 用 AOM の RF 生成、復調、fast/slow PID | Laser Lock Box | LO 1 mHz–600 MHz、IIR LPF 2.6 kHz–35 MHz（2 次/4 次）、fast/slow の 2 系統 PID (web: Moku:Pro LLB) |
| 残留位相の **ループ外** 測定 | Phasemeter | 1 kHz–400 MHz、4 ch、tracking BW 10 Hz–10 kHz、最大 6 µrad/√Hz (web: Moku:Pro Phasemeter) |
| アクチュエータ・ループの伝達関数測定 | Frequency Response Analyzer | 実機のプラント同定（[[2026-09-24 LiquidInst2025_pll]] の手順を Moku 上で） |
| BHD スペクトル | Spectrum Analyzer | スクイージングの測定帯 |
| 自動化 | Python API (`moku` パッケージ) | `PIDController.set_by_frequency(...)` などでゲインを設定できる (web: Liquid Instruments KB) |

- Moku:Delta は 8 入力 / 8 出力、2 GHz 帯域、5 GSa/s、同時 8 計器 (web: Moku:Delta 製品ページ)。構成 A の 2 ループ＋測定を 1 台に載せられる見込み (agent)。

---

## 4. 理論

### 4.1 ホモダイン測定とスクイージングの劣化

**損失**（shot noise を 1 とする単位）:

$$V_\pm = \eta\, e^{\pm 2r} + (1-\eta)$$

- $r$: スクイーズパラメータ。生成時のスクイージング量は $10\log_{10} e^{2r}$ dB
- $V_-$, $V_+$: スクイーズド・反スクイーズド quadrature の分散
- $\eta$: 全効率。$\eta = \eta_\mathrm{src}\,\eta_\mathrm{prop}\,\eta_\mathrm{mm}\,\eta_\mathrm{QE}\,\eta_\mathrm{el}$
  - $\eta_\mathrm{src}$: 共振器の脱出効率や導波路の内部損失
  - $\eta_\mathrm{prop} = 10^{-\alpha L/10}\cdot\eta_\mathrm{conn}$: ファイバ損失（$\alpha \approx 0.2$ dB/km, 1550 nm の SMF の典型値 (agent)）とコネクタ・スプライス
  - $\eta_\mathrm{mm}$: LO とのモード整合（ビジビリティ $\mathcal V$ なら $\mathcal V^2$）
  - $\eta_\mathrm{QE}$: PD の量子効率
  - $\eta_\mathrm{el}$: 電子雑音による実効的な効率（クリアランス $C$ dB なら $\approx 1-10^{-C/10}$）(agent)
- **意味**: 1 km のファイバ損失は約 4.5 %。短いファイバではコネクタ損失の方が大きい (agent)。

**位相揺らぎ**:

$$V(\theta) = V_-\cos^2\theta + V_+\sin^2\theta$$

- $\theta$: LO の位相とスクイーズ楕円の向きのずれ [rad]
- Aoki ら は rms 位相揺らぎ $\tilde\theta$ を「オフセット $\tilde\theta$」として扱う形 $R'_\pm \approx R_\pm\cos^2\tilde\theta + R_\mp\sin^2\tilde\theta$ を使っている (paper: Aoki 2006, Eq. 5)。
- $\theta$ が平均 0、rms $\sigma$ のガウス分布なら、平均は (agent):

$$\langle V\rangle = V_- + (V_+ - V_-)\,\frac{1-e^{-2\sigma^2}}{2} \;\approx\; V_- + (V_+-V_-)\,\sigma^2$$

- **数値例** (agent): $V_- = -10$ dB, $V_+ = +20$ dB のとき
  - $\sigma = 10$ mrad → $-9.6$ dB
  - $\sigma = 30$ mrad → $-7.2$ dB
  - 反スクイージングが大きいほど、同じジッタでの劣化が大きい（[[Hirota2025_squeezing10db]] の「10 mrad ≈ 1 % 損失」と同じ話）

**どの周波数の位相揺らぎが効くか** (agent):

- スクイージングはサイドバンド周波数 $\Omega_\mathrm{meas}$（MHz 帯）で測る。位相揺らぎの周波数 $f$ が $\Omega_\mathrm{meas}$ や測定の RBW より十分低ければ、ゆっくりした回転として上の式がそのまま使える。
- したがって $\sigma^2$ は「測定時間の逆数から、ループの UGF より十分上まで」残留位相 PSD を積分したものと考える。プロトタイプでは 1 Hz–20 MHz で積分した。

### 4.2 誤差信号と検出雑音（probe 方式）

**誤差信号**: probe（パワー $P_\mathrm{pr}$、周波数 $\omega+\Omega$）と LO（$P_\mathrm{LO}$）が BHD でビートする (agent):

$$i_-(t) = 2\mathcal R\sqrt{P_\mathrm{LO}P_\mathrm{pr}}\;\cos(\Omega t + \Delta\phi), \qquad \mathcal R = \frac{\eta_\mathrm{QE}\,e}{h\nu}$$

- $i_-$: BHD の差電流 [A]、$\mathcal R$: 応答度 [A/W]、$e$: 電気素量、$h\nu$: 光子エネルギー
- $\Delta\phi$: LO と probe の相対位相。$\Omega$ で I/Q 復調すると $\sin\Delta\phi$ に比例する誤差信号になる
- 傾き $K_d \propto \sqrt{P_\mathrm{LO}P_\mathrm{pr}}$。LO を強くすれば probe を弱くしても傾きは稼げる

**shot noise で決まる位相換算の雑音床** (agent; 強い LO の極限で導出):

$$S_{\phi,\mathrm{sens}} \approx \frac{h\nu}{\eta_\mathrm{QE}\,P_\mathrm{pr}}\quad[\mathrm{rad^2/Hz}]\quad(\text{ヘテロダイン復調。ホモダインならこの半分})$$

- **数値例**: 1550 nm, $P_\mathrm{pr} = 1$ µW, $\eta_\mathrm{QE}=0.99$ → $1.3\times10^{-13}$ rad²/Hz。1 MHz 帯域で積分しても 0.4 mrad (agent)。
- **意味**: probe 方式では、µW の probe でもセンサ雑音はほぼ問題にならない。**高 QE の PD が効くのは主にスクイージングの効率 $\eta$ の方**（1 % の QE 改善 ≒ 1 % の損失減）(agent)。
- **例外**: ノイズロック（構成 C）では誤差信号の SNR がスクイージングの見え方そのものなので、QE や電子雑音が直接ロック性能にも効く (agent)。
- **probe パワーの上限**: probe は BHD に入るので、測定帯に余分な光・雑音を入れないよう µW 程度に抑える必要があった (paper: [[Hirota2025_squeezing10db]] Sec. 2)。

**probe とキャリアの周波数差によるオフセット** (agent):

$$\delta\theta_\mathrm{off} = \Omega\,\delta\tau_g,\qquad \tau_g = \frac{n_g L}{c}$$

- $\Omega$: probe のシフト周波数 [rad/s]、$\tau_g$: 群遅延、$n_g \approx 1.468$（SMF）
- 1 km で $\tau_g \approx 4.9$ µs。$\Omega/2\pi=40$ MHz なら $\Omega\tau_g \approx 1.2\times10^3$ rad。温度係数 $\delta\tau_g/\tau_g \sim 7\times10^{-6}$/K（熱光学係数と熱膨張の和の概算）を使うと約 9 mrad/K

### 4.3 ループのモデル

$$\phi_\mathrm{res} = S\,\phi_d - T\,n,\qquad S = \frac{1}{1+L},\quad T=\frac{L}{1+L}$$

$$L(s) = K_d\,C(s)\,A(s)\,e^{-s\tau_\mathrm{loop}}$$

$$S_\mathrm{res}(f) = |S|^2 S_d(f) + |T|^2 S_n(f),\qquad \theta_\mathrm{rms}^2 = \int S_\mathrm{res}\,df$$

- $\phi_d$: ファイバなどの外乱位相、$n$: センサ雑音（位相換算）
- $K_d$: 誤差信号の傾き [V/rad]、$C$: 制御器（Moku の PID）、$A$: アクチュエータ [rad/V]
- $\tau_\mathrm{loop}$: ループ内の遅延の合計
- $S_d$, $S_n$: それぞれの片側 PSD [rad²/Hz]

**アクチュエータで変わるプラントの形** (agent; [[PID Control]] の「平坦プラント vs 積分プラント」の区別と同じ):

- **EOM（位相アクチュエータ）**: $A = \pi/V_\pi$（平坦）。UGF 付近で $L\propto 1/s$ にするには **制御器側に積分が必要**（I や PI²）。P や PI だけだと高域で $|L|$ が下がらず、遅延で発振する（プロトタイプでも確認）。
- **AOM（周波数アクチュエータ）**: 周波数 $\delta f = K_\mathrm{VCO}V$ を変えると位相は積分される → $A = 2\pi K_\mathrm{VCO}/s$。プラント自体が積分器なので P や PI で $L\propto 1/s$ になる。
- **ファイバストレッチャ / PZT**: 平坦だが機械共振を持つ（915B の経験）。遅い大レンジ用に限定する。

### 4.4 遅延が決める帯域の上限

UGF 付近で $L \propto 1/s$ のとき、位相余裕は (agent):

$$\mathrm{PM} \approx 90^\circ - 360^\circ\, f_c\,\tau_\mathrm{loop} - (\text{ブースト等による遅れ})$$

- PM = 45° なら $f_c \approx 1/(8\tau_\mathrm{loop})$。PI² のブーストの零点を $f_c/4$ に置くとさらに 14° 失うので、$f_c \approx 0.086/\tau_\mathrm{loop}$ 程度

**遅延の内訳**:

| 要素 | 遅延 | 出典 |
|---|---|---|
| Moku:Pro PID Controller | 520 ns | (web: KB "closed-loop bandwidth of the PID") |
| Moku:Lab / Moku:Go PID | 約 800 ns / 約 3.6 µs | (web: 同上) |
| Moku:Pro LLB fast（LPF 1 MHz / 100 kHz） | 約 0.98 µs / 約 4.9 µs | (web: KB "latency for Laser Lock Box") |
| Moku:Pro LLB slow（LPF 1 MHz / 100 kHz） | 約 7.2 µs / 約 11 µs | (web: 同上) |
| Moku:Delta LLB fast / slow（LPF 1 MHz） | 約 0.97 µs / 約 7.13 µs | (web: 同上) |
| AOM の音響遅延 | $d/v_s$。TeO₂ の縦波 $v_s\approx4.2$ km/s, $d=1$ mm なら約 0.24 µs | (agent; 実機の値は要測定) |
| EOM ＋ドライバ | ns 程度（ドライバ帯域次第） | (agent) |
| ファイバ（片道） | $n_gL/c \approx 4.9$ µs/km | (agent) |

- **意味**: 構成 A（Moku＋EOM、約 1 µs）なら UGF は 80–120 kHz 程度まで。構成 B は 1 km で往復約 10 µs が加わり、UGF は約 10 kHz 以下になる (agent)。
- **LLB の LPF**: LPF の角周波数を下げるほど遅延が増える（100 kHz で約 5 µs）。復調のリップル除去と遅延のトレードオフになる (web; agent)。

**FNC の原理的な残差** (paper: Williams, Swann & Newbury 2008; 係数 1/3 は vault の過去の議論とも一致。式番号は未確認):

$$S_\mathrm{res,FNC}(f) \approx \frac{1}{3}\,(2\pi f\tau)^2\,S_\mathrm{fiber}(f)$$

- $\tau$: 片道の伝搬遅延、$S_\mathrm{fiber}$: ファイバ全体の自由走行位相雑音（雑音がファイバ全体に一様に分布すると仮定）
- **意味**: ループゲインをいくら上げても、往復の間に変わってしまう分は消せない。$S_\mathrm{fiber}\propto 1/f^2$ なら残差は周波数によらず一定の床になる (agent)

### 4.5 アクチュエータの選び方とレンジ

| | EOM | AOM（周波数制御） | ファイバストレッチャ |
|---|---|---|---|
| 速さ | 速い（MHz 級） | 中（音響遅延 0.x µs） | 遅い（機械共振） |
| レンジ | 有限: $\pm\pi\,V_\mathrm{max}/V_\pi$ | 実質無制限（周波数 → 位相は積分） | 数十 rad 級（機種による） |
| プラントの形 | 平坦 | 積分 | 平坦＋共振 |
| 注意点 | 駆動電圧。Moku の出力振幅で足りるか要確認 | 周波数を変える。RF ドライバの FM 帯域 | 共振、ヒステリシス |

- **オフロード（多段化）** (agent; Verclas の縦続 PID と同じ形 (paper: [[Verclas2026_distribution]] Fig. 3)):
  - **動作**: 速い PID の出力（EOM の電圧）を、遅い PID の入力にする。遅い側は EOM の電圧を 0 に戻す向きに AOM の周波数（またはストレッチャ）を動かす。
  - **意味**: 低周波の大きなドリフトを遅い大レンジのアクチュエータが受け持ち、EOM はレンジの中央付近で速い揺らぎだけを打ち消す。
  - **条件**: 遅いループのクロスオーバーは、速いループの UGF より十分低くする（目安 1/10 以下）。クロスオーバーで 2 経路の位相差が 180° に近いと、そこでゲインが打ち消し合う。
- **Moku 内で AOM の RF を作る場合**: Waveform Generator の変調源には内部正弦波、アナログ入力、他のアナログ出力を使える (web: KB "another output channel as a modulation source")。PID の出力で FM/PM をかける具体的な配線（Multi-instrument Mode 内のルーティング）は要確認 (agent)。

### 4.6 Moku のデジタル実装として入れるもの (agent)

- 遅延（4.4 節の表）と、LLB の復調 LPF（2 次/4 次 IIR）
- PID の積分器の飽和（`int_saturation`）と出力のクリップ
- サンプリング（時間領域シミュレーションは 1 µs 刻みでも 100 kHz 程度のループなら十分。ZOH とエイリアシングは [[2026-09-24 LiquidInst2025_pll]] 2.5 節と同じ扱い）
- 入力雑音（Moku:Delta は 100 Hz で 10 nV/√Hz 未満 (web)）を $K_d$ で位相に換算してセンサ雑音に加える

### 4.7 ファイバ効果: どこまで入れるべきか

| 状況 | 必要なモデル | 理由 (agent) |
|---|---|---|
| CW, squeezed vacuum＋µW の probe, SMF ≤ 数 km | 損失＋位相雑音 PSD＋遅延（線形） | Kerr 位相 $\gamma P L \approx 1.3\,\mathrm{W^{-1}km^{-1}}\times1\,\mu\mathrm{W}\times1\,\mathrm{km}\approx10^{-6}$ rad（$\gamma$ は SMF の典型値） |
| 明るい LO（mW 級）を同じファイバに通す | 上記＋XPM のオフセット・揺らぎ＋GAWBS | XPM 位相 $\approx 2\gamma P_\mathrm{LO}L$ は 10 mW, 1 km で約 0.03 rad。LO の強度雑音がそのまま位相雑音になる。GAWBS は明るい光に位相変調の側帯波（数十–数百 MHz）を作り、測定帯に入ると過剰雑音になる (paper: Shelby 1985) |
| パルス（ps 以下） | 分散によるモード不整合 → 時間モードの扱いが必要 | 分散長 $L_D = T_0^2/\lvert\beta_2\rvert$。$\beta_2\approx-22\,\mathrm{ps^2/km}$ で $T_0=1$ ps なら約 45 m、100 fs なら約 0.5 m |
| パルス＋明るい場の同伝搬 | 計画書のステップ 1（NLSE の線形化、Raman、GAWBS） | Drummond & Corney の確率微分方程式の枠組みが参考になる (paper) |

- **結論 (agent)**: CW で始めるなら、計画書のステップ 1 は丸ごと後回しにできる。パルスに決めた場合も、LO が同じ分散を受ける（または補償する）なら、まず時間モードの重なり $\eta_\mathrm{mm}$ を線形計算するだけで足りる。

---

## 5. 参考にできる論文・資料

### スクイーズド光の位相制御

- **Vahlbruch et al., "Coherent control of vacuum squeezing in the gravitational-wave detection band", PRL 97, 011101 (2006)** — CCSB（AOM でずらした制御光と 2 ループ）の原典。構成 A の loop 1 / loop 2 の設計に使う。
- **Chelkowski et al., "Coherent control of broadband vacuum squeezing", PRA 75, 043814 (2007)**, arXiv:0704.3796 — 広帯域版の CCSB。
- **"Phase control of squeezed vacuum states of light in gravitational wave detectors"**, arXiv:1411.3454 — 誤差信号の作り方を 2 通り比較し、スクイージングの安定性とミスアラインメントへの感度を論じる (web: arXiv 要旨)。著者・掲載誌は未確認。
- **[[McKenzie2005_noiselocking]]** — ノイズロック。構成 C。
- **Aoki, Takahashi, Furusawa, "Squeezing at 946 nm with periodically poled KTiOPO₄"**, arXiv:quant-ph/0511239 — 位相揺らぎの式 $R'_\pm \approx R_\pm\cos^2\tilde\theta + R_\mp\sin^2\tilde\theta$（Eq. 5）。掲載誌は未確認。
- **[[Hirota2025_squeezing10db]]** — probe を squeezed light から tap しない phase-detection OPA。高 QE（約 99 %）の自作 BHD。
- **[[Vahlbruch2016_squeezing15db]]** — 15 dB。損失・位相雑音の予算の立て方の手本。

### ファイバ伝送と位相雑音

- **[[Verclas2026_distribution]]** — 1 km ファイバでの squeezed vacuum 伝送。Moku:Pro で復調＋3 PID＋加算。構成 A/C に最も近い先行例。
- **Ma, Jungner, Ye, Hall, "Delivering the same optical frequency at two places: accurate cancellation of phase noise introduced by an optical fiber or other time-varying path", Opt. Lett. 19, 1777 (1994)** — AOM を使う FNC の原典。構成 B。
- **Williams, Swann, Newbury, "High-stability transfer of an optical frequency over long fiber-optic links", JOSA B 25, 1284 (2008)** — 遅延で抑圧できない残差、往路と復路を別ファイバにしたときの差動雑音 (paper)。
- **Huo et al., "Deterministic quantum teleportation through fiber channels", Sci. Adv. 4, eaas9401 (2018)** — EPR 状態を 3.0 km、テレポーテーションを 6.0 km のファイバで (web: 書誌)。CV 量子状態のファイバ伝送の先行例。位相ロックの詳細は未確認。

### ファイバ中の量子雑音（パルス・明るい場の場合）

- **Shelby, Levenson, Bayer, "Guided acoustic-wave Brillouin scattering", PRB 31, 5244 (1985)** — GAWBS。
- **Drummond & Corney, "Quantum noise in optical fibers I: stochastic equations" / "II: Raman jitter in soliton communications", JOSA B 18 (2001)**, arXiv:quant-ph/9912095, 9912096 — Raman・GAWBS を含む確率微分方程式。計画書の truncated Wigner の拡張に対応。
- **Wasilewski, Lvovsky et al., "Pulsed squeezed light: simultaneous squeezing of multiple modes", PRA 73, 063819 (2006)** — パルスの多モード性と Bloch–Messiah。
- **Lvovsky, "Squeezed light" (review)**, arXiv:1401.4118 — 損失・位相揺らぎ・ホモダインの基礎。

### 制御・ツール

- **Bechhoefer, "Feedback for physicists", RMP 77, 783 (2005)** — 遅延・感度関数・ループ整形。
- **Linien (arXiv:2203.02947), PyRPL (arXiv:2310.00086)** — FPGA ロックの実装例。
- **Weedbrook et al., "Gaussian quantum information", RMP 84, 621 (2012)** — 共分散行列の扱い。
- Moku 関連 (web):
  - [Laser Lock Box の遅延](https://knowledge.liquidinstruments.com/lase-lock-box/what-is-the-input-to-output-latency-for-laser-lock-box)
  - [PID Controller の閉ループ帯域と遅延](https://knowledge.liquidinstruments.com/141089-pid-controller/what-is-the-closed-loop-bandwidth-of-the-pid)
  - [Moku:Pro Laser Lock Box](https://dev.liquidinstruments.com/products/integrated-instruments/laser-lock-box-mokupro)
  - [Moku:Pro Phasemeter](https://dev.liquidinstruments.com/products/integrated-instruments/phasemeter-mokupro)
  - [Moku:Delta](https://liquidinstruments.com/products/hardware-platforms/Mokudelta/)
  - [Python API: PID controller の例](https://knowledge.liquidinstruments.com/python-api-examples/python:-pid-controller)、[moku-examples](https://github.com/liquidinstruments/moku-examples)
- Finesse 3 (web): [Squeezing and homodyne detection](https://finesse.ifosim.org/docs/latest/examples/quantum_mechanics/homodyne.html)、[Control and filtering](https://finesse.ifosim.org/docs/latest/usage/model_building/elements/control_and_filtering.html)、[Locking actions](https://finesse.ifosim.org/docs/latest/examples/control/lock_actions.html)

---

## 6. 実装案

### 6.1 ツールの見直し

| ツール | 判断 | 理由 (agent) |
|---|---|---|
| numpy / scipy（＋ python-control） | **主役** | 遅延は周波数領域で厳密に掛けられる。python-control は margin 計算や離散化に便利だが必須ではない |
| ガウス状態の共分散行列（自前 numpy） | 主役 | 単一モードは 2×2。多モード・GBS と接続するなら thewalrus / MrMustard |
| Finesse 3 | 検算用 | OPO 共振器を使う場合の CCSB 誤差信号の形と傾き |
| gnlse-python / PyNLO | 保留 | パルス＋明るい場の同伝搬に決まったら |
| QuTiP | 保留 | 非ガウス状態や測定ベースのフィードバックを扱うときに |
| Tidy3D / XMDS2 | 不要 | SMF ならデータシートの値で足りる |

### 6.2 モジュール構成

```text
phasesim/
  params.py      # 全パラメータを dataclass で一元管理（仮定値には出典・仮置きの印をつける）
  noise.py       # 外乱 PSD: ファイバ（長さ依存）, 熱ドリフト, 音響; PSD から時系列を合成
  optics.py      # 共分散行列, 効率の積, 位相ジッタ → 測定スクイージング; 誤差信号の傾きと shot noise
  actuators.py   # EOM（平坦, レンジ）, AOM（積分, 音響遅延）, ストレッチャ（共振モデル or 実測 FRF）
  moku.py        # 遅延, 復調 LPF（IIR）, PID（set_by_frequency と同じパラメータ化）, 飽和
  loop.py        # L, S, T, 余裕, 残留 PSD と θ_rms（周波数領域）
  timesim.py     # 離散時間ループ: sin 型誤差信号, 飽和, オフロード, ロック引き込み, 2π スリップ検出
  ccsb.py        # loop 1（probe↔pump, OPA 利得の位相依存）と loop 2（LO↔probe）の結合
  sweep.py       # L, P_pr, η_QE, Ω, ゲインの掃引
  bo_bridge.py   # 実機と同じ目的関数（θ_rms, IAE）を返す → Ax/BoTorch に接続
tests/           # 6.4 節の検証
```

- **パラメータ化の方針**: PID は Moku の API と同じ形（`prop_gain`, `int_crossover`, `diff_crossover`, `int_saturation`）で持つ。シミュレーションで決めた値をそのまま実機に書ける (agent)。

### 6.3 段階的な進め方

1. **第 0 段（今回のプロトタイプ）**: 周波数領域の雑音バジェット。構成 A/B、EOM/AOM の比較。時間領域でオフロードの必要性を確認。
2. **第 1 段: Moku のデジタルモデル**: LLB の LPF、fast/slow の 2 系統、積分器飽和、出力レンジを入れる。ロック引き込み（スキャン → 捕獲）と、ロックが外れたときの再捕獲（積分器のリセット）を時間領域で。
3. **第 2 段: CCSB の 2 ループ**: loop 1 で OPA の利得の位相依存（probe の増幅・減衰）を入れ、loop 1 の残差が loop 2 の誤差信号にどう漏れるかを見る。スクイーザーを OPO にするなら、ここで Finesse を検算に使う。
4. **第 3 段: 実機データで置き換え**: Moku の Phasemeter で実際のファイバの自由走行位相雑音を測り、`noise.py` の仮置きモデルを置き換える。FRA で EOM・AOM・ストレッチャの伝達関数と遅延を測り、`actuators.py` を実測 FRF に置き換える。
5. **第 4 段（必要なら）**: パルス・非線形のモジュール。出力は「$\eta_\mathrm{mm}$、過剰雑音、位相雑音 PSD」の数個の数値に限定し、上の枠組みに差し込む。

### 6.4 検証テスト

- **時間領域 ↔ 周波数領域**: 同じ離散ループで $\theta_\mathrm{rms}$ が一致すること（プロトタイプで 27.5 vs 27.2 mrad, 1 % 差）。
- **雑音合成の正規化**: 合成した時系列の分散が PSD の積分と一致すること（低周波のビンが少ないと統計誤差が大きいので、複数シードで平均）。
- **FNC の遅延限界**: ゲインを上げた極限で $\tfrac13(2\pi f\tau)^2 S_\mathrm{fiber}$ に漸近すること（Williams 2008）。
- **スクイージング**: $\sigma=0$ で損失だけの式に戻ること、小角で $V_-+(V_+-V_-)\sigma^2$ に一致すること。
- **既知の実験の再現**: [[Vahlbruch2016_squeezing15db]] や [[Hirota2025_squeezing10db]] の損失・位相雑音の予算を入れて、報告値に近いスクイージングが出ること。

### 6.5 BO との接続 (agent)

- シミュレーションに実機と同じ目的関数（$\theta_\mathrm{rms}$ や IAE）を持たせれば、Ax/BoTorch のパイプライン（制約の定義、初期 Sobol 点、ロック失敗の扱い）を **実機なしで** 試せる。
- Moku の Python API でゲインを設定できるので、シミュレーションで決めた探索空間を実機にそのまま移せる。
- 発展案: シミュレーションを低忠実度、実機を高忠実度とする多忠実度 BO や、シミュレーションの最適解を事前分布に使う方法。

---

## 7. プロトタイプの結果

スクリプト: `assets/2026-10-07_phase_lock_budget.py`。数値はすべて **仮置きの雑音モデル** に依存するので、絶対値ではなく **構成どうしの比較と傾向** として読む (agent)。

### 7.1 ファイバ長に対する残留位相と測定スクイージング

![[assets/2026-10-07_budget_vs_length.png]]

- 各構成・各長さで、位相余裕 45° 以上・ゲイン余裕 6 dB 以上を満たす UGF を全探索し、$\theta_\mathrm{rms}$ 最小の設定を選んだ。
- 生成 15 dB、$\eta \approx 0.89$（脱出 0.95, モード整合 0.98, QE 0.99, コネクタ 0.97）＋ファイバ損失。

| L | A-EOM | A-AOM | B-FNC | 位相揺らぎなし |
|---|---|---|---|---|
| 10 m | 1.7 mrad | 3.0 mrad | 3.2 mrad | −8.71 dB |
| 100 m | 5.4 mrad / −8.57 dB | 9.3 mrad / −8.52 dB | 18 mrad / −8.31 dB | −8.59 dB |
| 1 km | 17 mrad / −7.42 dB | 30 mrad / −7.06 dB | 355 mrad（ロック不成立相当） | −7.61 dB |
| 10 km | 54 mrad / −2.96 dB | 93 mrad / −2.17 dB | — | −3.43 dB |

- **読み取れること** (agent):
  - 構成 A では、短いファイバのジッタは小さく、1 km を超えるとファイバ損失が支配的になる。
  - 構成 B は 100 m あたりから往復遅延の効果で急に悪化する。
  - A-EOM と A-AOM の差は AOM の音響遅延（0.5 µs と仮定）による UGF の差。

### 7.2 1 km での PSD の内訳

![[assets/2026-10-07_psd_1km.png]]

- A-EOM は UGF 83 kHz、B-FNC は UGF 7 kHz（往復遅延で制限）。B-FNC の低周波の平らな床が、遅延で打ち消せない残差 $\tfrac13(2\pi f\tau)^2S_\mathrm{fiber}$ (agent)。

### 7.3 時間領域: EOM のレンジ飽和とオフロード

![[assets/2026-10-07_range_offload.png]]

- 100 m 相当の音響雑音（仮置き）＋ 30 rad/s の熱ドリフト（仮置き）。誤差信号は $\sin\Delta\phi$、Moku の遅延 1 µs。
- **EOM のみ**: EOM が ±π に張り付き、2π スリップが起きる（後半の rms 0.72 rad）。
- **EOM ＋ AOM オフロード**（クロスオーバー約 300 Hz）: EOM は中央付近にとどまり、スリップなし（後半の rms 27 mrad）。
- オフロードを 50 Hz にすると、この仮置きの雑音（低周波が大きい）ではまだスリップが残った。オフロードの帯域は実測の低周波雑音で決める必要がある (agent)。

---

## 8. 仮定一覧（プロトタイプ）

| 項目 | 値 | 種別 |
|---|---|---|
| ファイバ位相雑音 | $S_d = (L/1\,\mathrm m)\cdot10^{-4}\,(100\,\mathrm{Hz}/f)^2/(1+(f/5\,\mathrm{kHz})^2)$ rad²/Hz | 仮置き（実測で置き換える） |
| 長さ依存 | PSD ∝ $L$（独立な区間の和） | 仮定 |
| 熱ドリフト | 30 rad/s のランプ | 仮置き |
| Moku の遅延 | 1.0 µs（LLB fast, LPF 1 MHz 相当） | web の値を丸め |
| AOM の音響遅延 | 0.5 µs | 仮定 |
| probe パワー（BHD 到達） | 1 µW | 仮定 |
| PD 量子効率 | 0.99 | 仮定 |
| BHD クリアランス | 15 dB | 仮定 |
| 群屈折率 | 1.468 | 典型値 |
| ファイバ損失 | 0.2 dB/km、コネクタ合計 0.97 | 典型値・仮定 |
| 安定性の要求 | PM ≥ 45°, GM ≥ 6 dB | 設計値 |

---

## Open questions

- **CW かパルスか**: これで計画書のステップ 1 が要るかどうかが決まる。
- **LO の送り方**: 別ファイバ / 同じファイバ（偏光多重）/ 受信側でローカル生成（光 PLL）。差動雑音と XPM の有無が変わる。
- **スクイーザー**: OPO 共振器か、単一パスの導波路 OPA か（loop 1 のプラントと Finesse の要否が変わる）。
- **使う Moku**: Pro か Delta か。LLB と Phasemeter と FRA を Multi-instrument Mode で同時に載せられるか。
- **EOM の駆動**: Moku のアナログ出力の振幅で $V_\pi$（例: 4 V）に足りるか。足りなければアンプが要り、その帯域と雑音もループに入る。
- **AOM の実測遅延と FM 帯域**: 音響遅延、RF ドライバの FM 応答。Moku 内で AOM の RF を作り PID で FM/PM をかけられるか。
- **probe のシフト周波数 $\Omega$**: ループ帯域・測定帯・$\Omega\tau_g$ の温度依存オフセットのバランス。
- **実際のファイバの自由走行位相雑音**: Moku の Phasemeter で probe のビートを測って PSD を得る（仮置きモデルを置き換える最重要の入力）。

## Related concepts

- [[Coherent Probe Locking]]
- [[Noise Locking]]
- [[Homodyne Detection]]
- [[Squeezed States]]
- [[PID Control]]
- [[Lock-in Detection]]
- [[Phase Modulation]]
- [[Shot-Noise Limit]]
- [[Active Phase Stabilization in Photonic QC]]

## Related papers

- [[Verclas2026_distribution]]
- [[Hirota2025_squeezing10db]]
- [[Vahlbruch2016_squeezing15db]]
- [[McKenzie2005_noiselocking]]

## Related discussions

- [[2026-09-24 LiquidInst2025_pll]]
- [[2026-10-05 Hirota2025_squeezing10db]]
- [[2026-09-28 Verclas2026_distribution]]
- [[2026-10-05 stalwart_loop-id-design]]
- [[2026-10-07 phase-stab-sim_implementation-plan]] — 非線形ソルバを含む Python 実装計画（続編）
