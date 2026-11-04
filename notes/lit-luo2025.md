# 追加文献 Luo et al. 2025 — 420 km の量子メモリのもつれと「2 周波数」のロック光（lit-luo2025）

- Phase: 6（Procedure step 1、ユーザが追加した文献。P1–P7 の表にはない）、日付: 2026-11-04
- 原典: X.-Y. Luo, C.-Y. Wang, M.-Y. Zheng, B. Wang, J.-L. Liu, B.-F. Gao, J. Li, Z. Yan, Q.-M. Ke, D. Teng, R.-C. Wang, J. Wu, J. Huang, H. Li,
  L.-X. You, X.-P. Xie, F. Xu, Q. Zhang, X.-H. Bao, J.-W. Pan, "Entangling quantum memories over 420 km in fiber",
  arXiv:2504.05660v1 [quant-ph]（2025-04-08）（`refs/2026-11-04_extra_refs/2504.05660v1.pdf`、本文 pp.1–8・Supplemental Material pp.9–15、
  `pdftotext -layout` で読み、Fig. 2 は p.3 を画像にして確認した）。
- 入力: `docs/phase6_prompt.md`、`DECISIONS.md`（D30・D31・D40・D77・D78・D80）、`notes/conventions.md`、
  `notes/probe-cancellation-breakdown.md`（T10 §3.3）、`notes/lit-pittaluga2021.md`（P5）、`presets/future/F_common.yaml`、`presets/axes/probe/`。
- 値の印: **[L p.x]** は原典の頁（arXiv の頁番号）から読んだ値、**[AE-P6]** はこのノートの概算、**UNVERIFIED** は原典で確認できない値。
- 規約: 本文は日本語、記号は英語（D19）。原典の k は真空でなくファイバ中の波数（式 (S12) で屈折率を含む）。

## 0. 結論の要約

1. **方式**: 2 つの原子集団メモリ（DLCZ）からの write-out 光子を QFC で 1522 nm に変換し、中央の Charlie で単一光子干渉させる。位相は 2 段で固定する [L p.4, Fig. 2b]:
   (a) **常時ループ（full-time loop）**: 1600 nm（QFC のポンプのレーザーから）のロック光を長いファイバに重ね、Charlie で WDM で取り出して干渉させ、APD → EOM。78 nm 離れているので常時動かせ、迷光は < 1 Hz に抑えられる。
   (b) **間欠ループ（signal loop）**: write-out 光子の周波数の**上下対称に ±δν（δν = 675 MHz）の 2 成分**を持つロック光を、メモリの準備（MOT の準備）の時間だけ MEMS スイッチで APD に通し、ピエゾのストレッチャを動かす。もつれ生成の時間は電圧を保持する。
2. **2 周波数（原典の "dual-band"）の要点**: 2 成分の強度が等しいと、干渉の DC 成分は I ≈ S + Ã cos(k0ΔL)、Ã = 2(A1A2 + B1B2) cos(δk ΔL) になり、
   **位相は信号の波数 k0 で決まり、離調 δk は振幅 Ã のゆっくりした変調にしか効かない**（式 (2)、(S11)–(S14)）[L p.4, p.12]。
   単一の離調光だと (k_wo − k_lo)ΔL の偏りが残り、100 km で 1.16 rad/K と見積もっている（式 (S8)）[L p.11]。
3. **数値**: 位相の残差 3.1°（20 km）、3.5°（120 km）、3.6°（220 km）、6.4°（320 km）、6.2°（420 km）、全体で約 7°。可視度 0.98–0.93 [L p.3 Fig. 2d, p.2]。
4. **P8 を提案する**: P1–P7 にない考え方は「**上下対称の 2 周波数で、周波数をずらした参照光の偏りを一次で打ち消す**」こと。
   私たちの Liu 型 probe を ±δν の 2 周波数にすると、周波数をずらした probe の熱の偏り 2πΔν·κ_τ·L（T10 式 (7)）が**2 成分の強度の不均衡 r 倍に減る**（§3.2、[AE-P6]）。
   代わりに振幅 Ã が温度で変わる（ループ利得の変化、0 を横切ればロック外れ）。間欠の部分は P2 の窓と同じ。
5. **シミュレータ**: 2 周波数の probe はビームのオーバーレイ ＋ `liu_probe_interferometer` 2 つ ＋ `adder` で表せる（新しいブロックなし）。
   原典どおりの「常時（離れた波長）＋ 間欠（2 周波数）を保持」には、誤差を時間窓で止めて制御器を保持する NEW ブロック `window_gate` が要る（P2 と共用）。

## 1. 何をしているか

### 1.1 系（Fig. 1、pp.1–3）

- Alice・Bob の 2 つのメモリノードは USTC の同じ実験室、配備ファイバ 10.1 km × 2 本で Charlie（Hefei Software Park）へ [L p.2, Fig. 1a]。
  Fig. 1a の地図に 6.5 km の表示もあるが、本文に説明がない（**UNVERIFIED**、直線距離と思われる）。
- 各ノード: 共振器で増強した DLCZ メモリ（⁸⁷Rb の集団、780 nm）＋ QFC（RPE-PPLN 導波路、780 nm → 1522 nm、ポンプ 1600 nm）[L p.2, Supp. I]。
  1522 nm は 1600 nm のポンプから 78 nm 離れ、反ストークス Raman の雑音を減らすため C 帯のすぐ下を選んだ [L p.3]。
- Charlie: BS_wo で 2 つの write-out 光子を干渉させ、狭帯域のフィルタ（70 MHz のエタロン、VBG、LP・SP・BP）を通して SNSPD 2 台で数える [L p.3, Fig. 1c]。
- 距離を伸ばすため、120 km 以上では ULL のコイル（G.654.E）を write-out 側に足す [L p.5, Table I]。

### 1.2 位相の 2 つの干渉計（Fig. 2a・2b、p.3–4、Supp. II）

- 測る位相 Δφ_wr = Δ(ϕ_w + ϕ_r) + Δ(ϕ_wo + ϕ_ro)（式 (S4)）を 2 つの干渉計で別々に固定する [L p.10]。
  - write–read 干渉計（局所、自由空間）: ロック光を PBS_w のもう一方の入力から入れ、ピエゾで固定 [L p.4, Fig. 2a]。
  - write-out–read-out 干渉計（長いファイバと QFC を含む）: 下の 2 つのループで固定 [L p.4, Fig. 2b]。

### 1.3 2 つのループ（p.4）

| ループ | 光 | 入れ・出し | 検出 | アクチュエータ | 時間 | 出典 |
|---|---|---|---|---|---|---|
| 常時（full-time） | 1600 nm（ポンプのレーザーから） | ノードで QFC の直後に DM で入れ、Charlie で BS_wo の手前に WDM で出す | 別の BS で干渉、APD | EOM | 常時。probe 由来の雑音は SNSPD で < 1 Hz | [L p.4] |
| 間欠（signal loop） | write-out 光子の周波数 ±δν の 2 成分（δν = 675 MHz）、read-out 経路の BS_ro の後の不平衡 BS（T:R ≈ 88:12）から入れる | BS_wo の後で MEMS スイッチで APD へ | APD（帯域は MHz 級） | ピエゾのストレッチャ | **MOT の準備の時間だけロックし、もつれ生成の時間は電圧を保持** | [L p.4, Fig. 1c, Fig. 2b] |

- 以前の実験（ref. [12] = Yu et al., Nature 578, 240 (2020)）では間欠ループのロック光は単一の離調光だった（原子と相互作用しないようリング共振器の共振と write-out の周波数から離した）。
  数百 km では 2 経路の長さの差のドリフトで偏りが出るので、2 周波数に替えた [L p.4]。

### 1.4 単一の離調光の偏り（Supp. II、式 (S5)–(S8)、p.10–11）

- 離調したロック光 k_lo で固定すると、残る位相は φ_rs = (k_ro − k_lo)ΔL_ro + (k_wo − k_lo)ΔL_wol（式 (S8)）。
- 数値: k_ro − k_lo ≈ 215 m⁻¹、k_wo − k_lo ≈ 21.2 m⁻¹、L_ro ≈ 10 m、L_wo ≈ 100 km、線膨張 5.5×10⁻⁷ /K →
  read-out 側 0.0012 rad/K（無視できる）、write-out 側 **1.16 rad/K**（2 本のファイバに 0.2 K の温度差があれば効く）[L p.11]。Fig. S3 で遅延線で ΔL を変えて確かめた [L p.11, p.13]。
- **注意（[AE-P6]）**: この見積もりは線膨張係数だけで、屈折率の温度係数を入れていない。離調した 2 つの光の位相差の温度係数は群遅延の温度係数で決まり、
  T08・T10 の κ_τ = 39 ps/(km·K) を使うと (c/n_g)κ_τ ≈ 8×10⁻⁶ /K で、5.5×10⁻⁷ の約 15 倍。偏りは原典の見積もりより大きく、2 周波数にする理由は強まる。

### 1.5 2 周波数の式（式 (2)、(S9)–(S14)、p.4、p.11–12）

- 2 成分 ν± = ν0 ± δν（振幅 A（青）・B（赤）、添字 1・2 = Alice・Bob の経路）を BS で合わせると、APD の強度は（2δν = Δω/2π ≈ 1.3 GHz のビートは MHz の帯域で平均される）

  I = S + 2(A1A2 + B1B2) cos(k̄ΔL) cos(½ΔkΔL) + 2(A1A2 − B1B2) sin(k̄ΔL) sin(½ΔkΔL)  （式 (S11)）

  S = A1² + A2² + B1² + B2²、k̄ = (k+ + k−)/2、ΔL = L1 − L2。
- 屈折率が周波数に一次で変わるとすると k̄ = k0 + (2π/c)δnδν ≈ k0、½Δk = (2π/c)(n0δν + ν0δn) = δk ≪ k0（式 (S13)。δk は群屈折率で決まる）。
- 2 成分の強度を等しくする（A1A2 − B1B2 ≈ 0）と I ≈ S + Ã cos(k0ΔL)、Ã = 2(A1A2 + B1B2) cos(δkΔL)（式 (2)・(S14)）。
  **干渉の位相は write-out 光子の k0 で決まり、偏りがない**。δkΔL は振幅のゆっくりした変調。
- 確認: 空気の遅延線で ΔL を 0–15 cm 動かし、Vpp の変化を正弦で当てはめて δk = 0.142 cm⁻¹ → δν ≈ 678 MHz（設定 675 MHz と一致）[L p.4, Fig. 2c, p.12]。
  Fig. 2c の Vpp は約 3.9 V から ΔL ≈ 11–12 cm で最小（π/(2δk) = 11.1 cm [AE-P6]）。
- 振幅の変化は 15 または 30 分ごとにファイバの遅延線で補償した [L p.12]。

## 2. 数値

### 2.1 波長・周波数・パワー

| 量 | 値 | 出典 |
|---|---|---|
| メモリの光子 | 780 nm（Rb） | [L p.2] |
| 変換後の信号 | 1522 nm（S 帯） | [L p.2–3] |
| QFC のポンプ・常時ループのロック光 | 1600 nm（信号から 78 nm、周波数で約 9.6 THz [AE-P6]） | [L p.3–4] |
| 間欠ループの 2 成分 | ±δν、δν = 675 MHz（当てはめ 678 MHz） | [L p.4] |
| 2 成分の作り方（EOM の側帯など） | 記載なし | **UNVERIFIED** |
| ロック光のパワー（常時・間欠とも） | 記載なし | **UNVERIFIED** |
| 不平衡 BS | T:R ≈ 88:12 | [L p.4] |
| QFC のポンプ | 300 mW、端から端の効率 約 44%、雑音 約 3.2 kHz（SNR > 40:1） | [L p.3–4, p.9] |
| レーザー | 780 nm の 2 台と 1600 nm の 1 台をノードで共有し、すべて基準共振器にロック | [L p.4, p.11] |

### 2.2 SNSPD と計数率

| 量 | 値 | 出典 |
|---|---|---|
| SNSPD の暗計数 | 約 1 Hz | [L p.4] |
| 検出効率 | 約 60%（D1）、30%（D2） | [L p.4] |
| 常時ループの probe による雑音 | 実験を通して < 1 Hz | [L p.4] |
| 遠方のクリック率 | Hz 級まで下がる | [L p.4] |
| 距離ごとの SNR（Ψ+） | 32.2（0 km）、48.8、48、26、15、3.5（420 km） | [L p.6, Table I] |

### 2.3 ループの帯域・時間の配分

- 常時ループ・間欠ループの帯域、APD の帯域（「MHz 級」以外）、MOT の準備ともつれ生成の時間の比（間欠ループのデューティ）は**記載なし**（**UNVERIFIED**）。
- 振幅の補償: 15 または 30 分ごと [L p.12]。
- 可視度が残差からの見積もりより低い理由として、ロックし直す時間（距離とともに頻度が増える。ストロークを伸ばせば改善）と強度の不均衡を挙げている [L p.4]。

### 2.4 距離に対する残差（Fig. 2d、p.3。弱いプローブパルスを BS_ro から入れ、BS_wo で干渉させて SNSPD で数えた）

| 距離 | 0 km | 20 km | 120 km | 220 km | 320 km | 420 km |
|---|---|---|---|---|---|---|
| ファイバの損失（dB、Table I） | — | 8.0 | 26.7 | 43.3 | 61.2 | 78.7 |
| 位相の残差（Fig. 2d の左下） | 表示なし（局所、QFC なし） | 3.1° | 3.5° | 3.6° | 6.4° | 6.2° |
| 同（mrad）[AE-P6 換算] | — | 54 | 61 | 63 | 112 | 108 |
| 可視度 V_P（Fig. 2d の右下、Table S5） | 0.98(1) | 0.96(1) | 0.96(2) | 0.95(1) | 0.94(2) | 0.93(4) |

- 本文冒頭の「全体の位相の安定度は約 7°」（約 0.12 rad）[L p.2]。
- 残差と可視度の対応は Fig. 2d の画像で確認した（各パネルの左下が残差、右下が可視度。キャプションの "bottom right (left)"）。

### 2.5 迷光・クロストークの扱い

- 1600 nm のロック光は長いファイバで 1522 nm に Raman 雑音を作るので、検出前に除く [L p.3]。
- ノード側の広帯域フィルタ（DM 2 枚: 1522 nm で 99.9% 反射・1600 nm で < 5%、SP 1570 nm、LP 1500 nm）、Charlie 側の狭帯域フィルタ（エタロン FWHM 70 MHz・FSR 19 GHz・透過 92%、
  VBG FWHM 0.1 nm・回折効率 92%、BP 1 nm × 2）[L p.3, p.9]。
  本文 p.3 は「19 GHz の帯域の VBG」、Supp. p.9 は「VBG の FWHM 0.1 nm、エタロンの FSR 19 GHz」で食い違う（**UNVERIFIED**、Supp. の方が詳しい）。
- ±675 MHz の 2 周波数は 70 MHz のエタロンの外だが、もつれ生成の時間は MEMS で APD 側に切り替えないので、SNSPD には時間で入れない（間欠にする理由の一つ。原典は明示していない [AE-P6]）。

## 3. 私たちの系に移せるパラメータ

### 3.1 移せる考え方と、そのままでは移せないもの

| Luo の要素 | 私たちの系 | スケール |
|---|---|---|
| 2 周波数 ±δν で偏りを打ち消す | Liu 型 probe を ±δν の 2 周波数に（参照光も同じ 2 周波数を含む必要がある。DC のビートは同じ周波数どうしで出るため） | 長さに依らない考え方。必要性は δν·L で決まる（§3.2） |
| 2 周波数の強度を等しく | 不均衡 r = (A1A2 − B1B2)/(A1A2 + B1B2)（各成分のパワー比 1 + ε なら r ≈ ε/2） | 長さに依らない |
| 振幅 Ã ∝ cos(δkΔL) | Liu の干渉計の腕の差 ΔL_k ≈ L_k ＋ ポートに依らない部分（probe の分岐・参照の配送） | **ΔL が m 級以上**なので、δν をポートの遅延と合わせる必要がある（§3.3） |
| 常時の離れた波長のループ（78 nm） | WDM の probe（`probe/freq_1THz`・`sep_wdm`、P3） | 偏り 2πΔν·κ_τ·L は長さに比例。1 THz・63 m で 15 rad/K（T10） |
| 間欠ループの保持 | P2 の測定・安定化の窓と同じ | 窓の間のドリフトは長さに比例 |
| 残差 3–6°（20–420 km） | 指標 (ii) の目標 10 mrad 級より大きい | km の値は移せない。方式だけ移す |

### 3.2 2 周波数の偏りの式（[AE-P6]、式 (S11) から）

式 (S11) で側（I = S）にロックすると cos(k0ΔL)cos(δkΔL) + r sin(k0ΔL) sin(δkΔL) = 0 なので、ロック点の偏りは

  b = arctan[ r · tan(δkΔL) ] ≈ r·δkΔL  （|δkΔL| が π の整数倍に近いとき）

- r = 1（単一の離調光）で b = δkΔL（T10 式 (7) の偏り）、r = 0 で b = 0。**2 周波数は偏りを r 倍にする**。
- 温度で変わる部分: db/dT ≈ r·2πδν·κ_τ·ΔL（κ_τ = 3.9×10⁻¹⁴ s/(m·K)、`F_common.yaml`）。
- 振幅: |cos(δkΔL)| が下がるとループ利得が下がり、0 で符号が変わってロックが外れる。
- GVD の二次の項（2 成分の平均の位相 ½(2πδν)²β2L）は定数で、設定点に吸収される（δν = 7 GHz・1.27 km・β2 = −2.17×10⁻²⁶ s²/m で −27 mrad、温度係数は無視できる [AE-P6、dβ2/dT は UNVERIFIED]）。

| 例 | 単一の離調光の偏り | 2 周波数（r = 0.05）の偏り | 振幅 \|cos\| ≥ 0.9 の温度幅 | 10 mrad に要る温度の安定度（単一 / 2 周波数） |
|---|---|---|---|---|
| δν = 7 GHz、ΔL = 63 m | 0.108 rad/K | 5.4 mrad/K | ±4.2 K | 93 mK / 1.9 K |
| δν = 7 GHz、ΔL = 1.27 km（F_liu20 のポート 31） | 2.18 rad/K | 0.11 rad/K | ±0.21 K | 4.6 mK / 92 mK |
| δν = 80 MHz（AOM）、1.27 km | 0.025 rad/K | — | — | 0.4 K（単一で足りる） |

（すべて [AE-P6]。T10 の表の 7 GHz・63 m の 0.11 rad/K と一致。）
→ **2 周波数が効くのは δν を GHz 以上にしたい場合**（迷光を FBG・FSR のフィルタで落とすため。Borealis の +7 GHz、D40・D50）。AOM 程度の離調なら単一で足りる。

### 3.3 ポートごとの振幅の条件（[AE-P6]）

- ポート k の腕の差は ΔL_k = ΔL_0 + (L_k − L_0)、L_k − L_0 = k(N_b + n_guard)·c/(f_rep n_g)（D39）。δk = 2π n_g δν/c なので
  δk(L_k − L_0) = 2π δν k (N_b + n_guard)/f_rep。
- **2δν = m·f_rep/(N_b + n_guard)（m は偶数）にすると、全ポートで cos(δk(L_k − L_0)) = 1**（群屈折率が打ち消し合うので n と n_g の違いに依らない）。m が奇数だとポートごとに符号が交互になる。
  例: F_liu20（f_rep/N_b = 5 MHz）で δν = 7 GHz → 2δν = 14 GHz = 2800 × 5 MHz（偶数で条件を満たす）。
- 残るのはポートに依らない ΔL_0（probe の分岐 PRsrc・PRS_k と参照の REF・REF_k の長さの差）で、|cos| ≥ 0.9 には δk·ΔL_0 を π の整数倍から 0.45 rad 以内に合わせる:
  δν = 7 GHz（δk = 215 m⁻¹）で約 ±2 mm、675 MHz（δk = 20.8 m⁻¹）で約 ±2 cm。Luo は 15–30 分ごとに遅延線で補償した。
- 温度で δk·L_k が変わる分（上の表の「振幅の温度幅」）はポートごとに違う。TEC のあるポート（F のプリセット）では問題になりにくい。

### 3.4 間欠ループ（Luo そのもの）の私たちの版

- 常時: WDM の probe（離れた波長、`probe/freq_1THz`。Luo は 9.6 THz）で速く、迷光は WDM と Raman のフィルタで落とす。
- 間欠: 2 周波数の近い probe で遅く、GBS の測定中は切って制御器を保持（P2 の窓）。近い 2 周波数の迷光（±δν の Rayleigh）は窓の外では出ない。
- P5 との違い: P5 の遅い段は in-band の弱い光をクリックで数える（光子数で精度が決まる）。Luo の間欠ループは古典の明るい光を APD で測る（精度は高いが、参照は in-band ではなく ±δν の平均）。
  2 周波数の偏りの打ち消しが一次なので、Luo 型は P5 の「本物の光」の代わりになりうる（ただし r と二次の項、非共通区間は残る）。

## 4. シミュレータへの対応

### 4.1 既存のブロック・オーバーレイで表せること

1. **2 周波数の probe**（全時間で動かす版、P8-a）:
   - トポロジーのオーバーレイでビーム `pr1p`・`pr1m`（`freq_offset_hz` = ±δν）と、参照の `ref1p`・`ref1m`（同じ ±δν。参照も 2 周波数を含む）を足す。ポート 0 も同じ。
   - 熱の係数は `probe/freq_7GHz` と同じ形で、各ビームに ±δν の項（`F1.dT.pr1p.phase` など）。参照の区間 REF・REF1 にも ±δν の項を足す（長さが違うので ΔL_0 の分がここに入る）。
   - 誤差信号: `liu_probe_interferometer` を 2 つ（`in_pr: beam.pr1p.phi`・`in_ref: beam.ref1p.phi`、もう一つは m）＋ `adder`（`gains: [1, 1]`）→ v = K+ cos X+ + K− cos X− = 2K̄ cos X̄ cos(ΔX/2) ＋ 不均衡の項。
     各ブロックのパワーを成分ごとの値にすればショット雑音（2eℜP）は二重に数えない。電気雑音は片方のブロックだけオンにする。不均衡 r は `p_pr_avg_w`・`p_ref_w` を成分ごとに変えて掃引する。
   - 迷光: `stray_light_table` を成分ごとに 2 つ（`probe_freq_offset_hz` = ±δν）。
2. **常時の離れた波長のループ**: `probe/freq_1THz` ＋ `probe/sep_wdm`（既存）。Luo の 9.6 THz では T10 式 (7) の一次の式（2πΔν·κ_τ·L）が成り立つかが問題（κ_τ の分散、T06・T08。**モデルの限界**）。
3. **ホップの偏り**: `lit-pittaluga2021.md` §4.2(3) と同じ（probe のビームのストレッチャの係数を (1 + Δν/ν0) 倍）。2 周波数では ± で打ち消すので、一次では残らない。

### 4.2 NEW ブロック

**`window_gate`**（族 CTRL。Luo の間欠ループの保持。P2 の測定・安定化の窓と共用）

- 入力スロット: `in`（誤差など）。出力スロット: `out`、`ctrl.<name>.gate`（int、窓の中 1・外 0）。
- パラメータ（D34、既定値なし）: `period_s`、`duty`、`t_offset_s`、`mode: zero | hold`
  （zero: 窓の中は out = in、外は 0。積分だけの PID（F のプリセットの `pid*f`・`pid*s` は kp = 0）なら出力がそのまま保持される。
  hold: 窓の外は最後の窓の中の値を出し続ける）。トグル: `gating`（off で out = in）。
- **オーバーレイで足りない理由**: 時間で変わる掛け算・保持をするブロックがない。`adder` の `gains` は実行中に変わらない定数、
  `pid` の `feedback: false` は出力を `offset_v` に戻す（保持ではない）、`rate_div` はデューティ 1 の周期的な標本化・保持しかできない
  （`detect_after_demux` の 1 フレーム 1 回の標本化はこれ）、`probe_pulser` のゲートは迷光の勘定にしか効かない。

ponytail: 2 周波数を 1 ブロックにまとめる新しいブロックは作らない（2 ブロック ＋ `adder` で同じ式になる）。精度が足りなければ `liu_probe_interferometer` に成分のリストを持たせる拡張で済む（新しいブロックではない。`src/` の変更なので承認が要る）。

## 5. P8 の提案

`docs/phase6_prompt.md` の表と同じ列（英語、表にそのまま足せる形）:

| ID | Paper | Idea to import | Fusion with our system / what to test |
|---|---|---|---|
| P8 | Luo et al. 2025, arXiv:2504.05660 (420 km memory–memory entanglement) | Symmetric two-tone locking light (±δν about the signal, δν = 675 MHz, equal powers): the DC fringe follows k0·ΔL, so the frequency-offset bias cancels to first order and δk·ΔL only modulates the fringe amplitude. Used as an intermittent slow loop (held during data taking) under a full-time far-detuned (78 nm) fast loop | Liu probe and its reference as two tones ±δν (e.g. ±7 GHz, removed by FBG/FSR filters) instead of one offset tone (`probe/freq_7GHz`). Bias becomes r·2πδν·κ_τ·ΔL (r = tone imbalance); choose 2δν = m·f_rep/(N_b + n_guard), m even, so the envelope cos(δk·ΔL_k) is equal for all ports. Sweep r, δν, ΔL_0 trim and spool ΔT: bias and envelope (loop-gain) drift vs the single-tone probe, stray clicks at ±δν. Variant: full-time WDM probe (P3) + intermittent two-tone slow loop with hold windows (as P2) |

**P8 が P1–P7 と違う点**: P4(a)・P5 は周波数・波長の違いによる proxy の偏りを「モデルで見積もる」（P4）か「本物の光で遅く消す」（P5）。
P8 は**参照光の作り方で偏りを一次で消す**。P2（窓）・P3（WDM ＋ ゲート）・P5（2 段）と組み合わせられる。

## 6. UNVERIFIED と未解決の問題

**UNVERIFIED（原典で確認できない）**
1. ロック光のパワー（1600 nm、±δν の 2 成分）、2 成分の作り方、APD の帯域（「MHz 級」以外）。
2. 常時・間欠ループの帯域、EOM・ピエゾのレンジ、MOT の準備ともつれ生成の時間の比（間欠のデューティ）、ロックし直しの頻度。
3. Fig. 1a の 6.5 km の意味。VBG の帯域（本文 19 GHz、Supp. 0.1 nm）。
4. 0 km の残差（Fig. 2d に表示なし）。
5. §1.4 の注意（熱の係数が約 15 倍大きい）と §3.2 の偏りの式、§3.3 の条件はこのノートの導出・概算 [AE-P6]。dβ2/dT の値。

**未解決の問題**
1. ±δν の迷光: 逆向きの probe の Rayleigh・Brillouin は GBS 側へ戻る。±δν の両方を落とすフィルタ（FBG 2 枚、OPO の FSR フィルタの阻止帯に置く）が要る。OPO の FSR（3.7・7.4 GHz、T17）との位置関係。
2. 2 周波数の参照光の配送: 参照 REF も 2 周波数にすると LO_k（REF から分岐、D09）にも ±δν が乗る。LO を別に取るか、LO 側でフィルタするか。
3. 不均衡 r の時間変化（2 成分を EOM の側帯で作るなら変調度のドリフト）。r の揺れは偏りの揺れ（r·δkΔL）になる。
4. Liu の DC fringe（`homodyne_dc`）では振幅 Ã の変化がそのままループ利得の変化になる（E76 と同じ形）。正規化（P4(b)）と組み合わせるか。
5. 間欠版の窓の長さとデューティは P2 の掃引と共通にする（窓の間のドリフト・Dick 効果）。
