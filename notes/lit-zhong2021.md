# P6 Zhong et al. 2021（Jiuzhang 2.0）: 参照信号の遅延によるポートの位相のプログラム（lit-zhong2021）

- Phase: 6（`docs/phase6_prompt.md` の Procedure 1）、日付: 2026-11-04
- 入力: `refs/2026-11-04_extra_refs/PhysRevLett.127.180502.pdf`（本文 9 ページ、`pdftotext -layout` と図のページの画像で読んだ）、
  `DECISIONS.md`（D01–D83）、`notes/conventions.md`、`notes/00_summary.md`、`notes/survey.md`、`notes/probe-cancellation-breakdown.md`（T10）、
  `notes/error-signal-generation.md`（T14）、`notes/actuators-control-digital.md`（T16）、`presets/future/F_common.yaml`、`presets/axes/probe/`、`src/cranberry/`
- 値の印: **[Zhong p.N]** は PRL のページ番号 180502-N（PDF のページと同じ）で本文から読んだ値、**[Zhong Fig.]** は図から読んだ値、
  **[SM]** は本文が Supplemental Material [25] に回している内容（`refs/` に SM はなく**未読**）、**[AE]** はこのノートの概算・解釈（未検証）、
  **[code]** は `src/cranberry/` で確認した事実、**UNVERIFIED** は原典で確認できなかったもの。
- D34 に従い、ここで挙げる構成はすべて切り替えて試せる同格のプリセット・オーバーレイの案で、既存のプリセットは変えない。

## 0. 要点

1. **本文に書かれている位相のプログラムの方法は 1 文だけ**: "To change the phases, we add adjustable electric delay lines to the reference signal to which each TMSS is phase locked." [Zhong p.6]。
   参照信号の周波数・種類、ロックの方式・帯域・残差、遅延線の型は本文にない（[SM]）。
2. **`phase6_prompt.md` の「1450 nm の参照光」は本文に一語も出てこない。** [SM] にある可能性が高いが未確認（UNVERIFIED）。このノートは 1450 nm を前提にしない。
3. 本文の数値: 系全体の位相の揺れは **1 時間で 0.06 rad 以内**（Fig. 1(e) のヒストグラム、横軸 −0.04〜0.02 rad）、それによる光子の干渉の可視度の低下は**約 0.5 %** [Zhong p.3]。
   プログラムしたのは **25 個の TMSS の位相 × 30 組のランダムな設定**（Fig. 4(a)）[Zhong p.6]。
4. **自分たちの系への移し方** [AE]: Liu 型のポートごとのロックで、probe の検出を**ヘテロダイン**（AOM で周波数をずらす）にすると、復調の参照の位相 φ_d,k を動かすだけで
   ロック点 X_k が **−Δφ_d,k だけ 1:1 で動く**（`liu_probe_interferometer` の heterodyne は (a, b) = K(cos X, sin X)、`demodulator` は I = K cos(X + φ_d) [code]）。
   傾き K（probe のパワー・可視度・PD の利得）に依らないので、DC の設定点（X = arccos(V_set/K)、K の揺れ・RIN に敏感）や dither（参照の位相は傾きを変えるだけでロック点は極値のまま）より、プログラムの「目盛り」が正確。
   squeezed 光の経路に変調器を足さない（D10・T16 §3.2 の「EOM は squeezed 光路に置けない」と整合）。
5. **既存の部品で表せる部分**: ヘテロダインの Liu 型（`presets/axes/probe/detect_heterodyne.yaml`）＋ `demodulator` の実行時パラメータ `phi_d_rad` ＋ `metric_ii` の実行時パラメータ `setpoint_rad`（Θ_k,set）で、
   **1 回の実行の中で一定の**プログラム値は今すぐ試せる [code]。
   **実行の途中でプログラムを切り替える（ステップ応答・整定・ポート間の一致）には新しいブロックが要る**（§4.3）。

## 1. 論文がしたこと

| 段落（PRL は節番号なし） | ページ・図 | 内容 |
|---|---|---|
| Abstract | p.1 | 144 モードの干渉計で最大 113 クリックの GBS。入力の squeezed 状態の位相を変えてプログラムできる。Hilbert 空間の次元 ~10⁴³、総当たりの古典計算に比べ ~10²⁴ 倍速いサンプリング |
| Stimulated emission of squeezed photons | p.2、Fig. 1(a)–(c) | 775 nm の変換限界パルスを PPKTP に集光して TMSS を作り、凹面鏡でポンプと TMSS を戻して 2 回目のパラメトリック過程を誘導放出で起こす（2 パス）。石英のくさび板で戻すポンプと TMSS の相対位相を合わせ、干渉の可視度 0.951（Fig. S3 [SM]）。輝度は 1 パスの約 3 倍。ポンプと信号の群速度の一致で "self-phase stable" と書く。図の中の表記で繰り返しは 250 kHz [Zhong Fig. 1(a)] |
| 144-mode phase-locked interferometer | p.3、Fig. 1(d)(e) | 25 個の TMSS（= 50 個の単一モード squeezer 相当、Fig. 1(d) の "50 SMSSs"）を 72 空間モード × 2 偏光の 3 次元干渉計へ。透過 96.5 %、ファイバ結合・PBS・波長板・フィルタ込みで 78 %。ポート間の光路長差の較正の標準偏差 1.6 µm（Fig. S6 [SM]）→ 波束の重なり > 0.999。SNSPD 144 個、平均効率 83 %。**光路全体の能動的な位相ロックと干渉計内部の受動的な安定化**で、系全体の位相の揺れを 1 時間 0.06 rad 以内（Fig. 1(e)）、可視度の低下 約 0.5 % [SM] |
| Validating GBS | p.3–4、Fig. 2 | ポンプ 0.15–1.6 W、集光径 65・125 µm。Bayes 検定・部分系による検証。全体の透過率 48 %（65 µm）・54 %（125 µm）[Zhong p.4] |
| Limits on classical simulability | p.4–6、Fig. 3 | 非古典性の不等式（式 (2)）、切断した高次相関（最大 19 ± 1 次） |
| **Phase-programmable GBS** | **p.6、Fig. 4** | 干渉計は固定し、入力の TMSS の位相を変えて行列を組み替える。**30 組のランダムな位相**（Fig. 4(a)、25 個の TMSS、縦軸 0–6.3 rad を組ごとに 2π ずらして重ねて表示）。**位相は、各 TMSS がロックされている参照信号に「調整できる電気の遅延線」を入れて変える。** 2 点相関 C_ij の全変動距離: 同じ設定どうし 0.068(1)、違う設定どうし 0.32(1)（Fig. 4(b)）で 30 組を区別できた |
| The classical computational cost | p.6–7、Fig. 5 | 1 組の試料の収集時間は 200 s |
| 文献 [25] | p.8 | SM の内容: "stimulated parametric down-conversion sources, three-dimensional interferometer, single photon detectors, subsystem validation method, and **phase locking**" |

## 2. 数値（本文で確認できたもの・できなかったもの）

| 項目 | 値 | 出典 |
|---|---|---|
| 参照光の波長 | **本文に記載なし**（prompt の 1450 nm は未確認） | [SM]、UNVERIFIED |
| 参照信号（電気）の周波数・種類（変調・ビート・dither） | 記載なし | [SM]、UNVERIFIED |
| 位相をプログラムする手段 | 参照信号に入れた調整可能な電気の遅延線 | [Zhong p.6] |
| プログラムした位相の数 | 25（TMSS ごとに 1 つ）× 30 組 | [Zhong p.6、Fig. 4(a)] |
| 位相の範囲 | 0–2π（Fig. 4(a) の縦軸 0.0–6.3 の目盛りと「2π ずつずらして重ねた」という説明から） | [Zhong Fig. 4(a)] |
| ロックの帯域・利得・アクチュエータ | 記載なし | [SM]、UNVERIFIED |
| ロックごとの残差 | 記載なし | [SM]、UNVERIFIED |
| 系全体の位相の揺れ | 1 時間で 0.06 rad 以内（ヒストグラムの横軸 −0.04〜0.02 rad。標準偏差の値は本文にない） | [Zhong p.3、Fig. 1(e)] |
| 位相の揺れによる可視度の低下 | 約 0.5 % | [Zhong p.3]（導出は [SM]） |
| ポート間の光路長の較正 | 標準偏差 1.6 µm、波束の重なり > 0.999 | [Zhong p.3] |
| 設定を切り替えた後の整定時間、切り替えの速さ | 記載なし | UNVERIFIED |
| 設定が狙いどおり入ったかの検証 | 光子のデータの 2 点相関で 30 組を区別（Fig. 4(b)(c)）。位相そのものの誤差は書かれていない | [Zhong p.6] |

- 0.06 rad を自分たちの指標の言葉に直すと: ヒストグラムの幅（全幅）0.06 rad を ±3σ と見れば σ ≈ 10 mrad [AE]。D78 の線幅の表の「10 mrad」の列と同じ水準。
- ヒストグラムが「系全体」のどの量（ポート間の相対位相か、各ロックの誤差信号か、in-loop か out-of-loop か）なのかは本文から読めない（UNVERIFIED）。in-loop の誤差信号なら、T10 §0-1 の「誤差信号に現れない偏り」を含まない。

## 3. 自分たちの系に移せること

### 3.1 将来の構成（単一の squeezer → 振り分け → ポートごとの PM ファイバ遅延 → GBS、Liu 型の逆向き probe が基準）

構成の違い（移す前の注意）[AE]:
- Jiuzhang 2.0 は 25 個の独立な TMSS 源で、**各源の位相**（ポンプと TMSS の相対位相、2 パスの干渉）をロックしている。自分たちは単一の squeezer の時間多重なので、
  ポート k の入力の squeezing 角は θ_sq,k = θ_sq,0(源) + ψ_k（ポート k の光路の位相、conventions §4・§10）。**プログラムの対象は光路の位相 ψ_k** で、Liu 型のロックがそのまま受け口になる。
- TMSS の位相は 2π 周期だが、単一モードの squeezing 角は π 周期（conventions §10 の `wrap_quadrature`）。**ポートごとに要る位相の範囲は (−π/2, π/2]**（光路の位相として最大 π/2 の段差）。

移せる考え方（すべて [AE]）:
1. **プログラムはロックの設定点に入れる。** squeezed 光の経路に変調器を足さない（EOM の挿入損失 2.5–3.5 dB で squeezing の上限が 2.6–3.6 dB になる、T16 §3.2）。
   アクチュエータは今の Liu 型の PZT ストレッチャ ＋ TEC のまま。プログラムの段差 ≤ π/2 はレンジの消費が小さい（915B 級の ±26 rad に比べて 6 %）。
2. **誤差信号はヘテロダインにする。** Liu 型の probe を AOM で f_sh ずらし、ポートごとのビートを復調する（D40 の検出方式、`probe/detect_heterodyne`）。
   I_k = K_k cos(X_k + φ_d,k) なので、ロック点は X_k = π/2 − φ_d,k。**φ_d,k を Δ 動かすと ψ_k が −Δ（ロックの符号で ±）動き、K_k に依らない。**
   - DC の設定点（今の `homodyne_dc`）でプログラムすると X = arccos(V_set/K) で、K（probe パワー、可視度、偏光）の揺れがそのまま位相の誤差になり、傾き sin X も変わる（Ref D・P4 (b) と同じ問題）。
   - dither では復調の位相はロック点を動かさない（I ∝ sin X の極値ロック、`probe/detect_dither`）。プログラムには使えない。
   - Ref A 型の 2 f_sh の誤差（Φ1 = 2φ_pr − φ_p）では φ_d の Δ は位相の Δ/2 になる（T14 §3.2）。目盛りが半分になるだけで、考え方は同じ。
3. **参照の位相の作り方**: Jiuzhang は電気の遅延線。遅延 δτ は位相 2π f_ref δτ になる。f_ref = 1 MHz なら 1 ns で 6.3 mrad、80 MHz（`detect_heterodyne` の仮置きの f_shift）なら 1 ns で 0.50 rad。
   遅延線で π/2 の段差を作るには 1/(4 f_ref) の遅延（1 MHz で 250 ns）が要るので、**低い f_ref では遅延線より DDS（数値で位相を設定する発振器）の位相の方が現実的**。
   Moku:Delta の復調の参照に位相を設定できるか、ポートごとに独立な参照を何本出せるかは `notes/actuators-control-digital.md` §4.5 の表にない（UNVERIFIED、資料で確認する）。
4. **ポート間の一致**: 全ポートの参照を 1 つの RF 源から作れば、源の位相雑音は Θ_k = θ_sq,k − θ_sq,0 で打ち消し、残るのは参照の配線の遅延差 τ_RF の分 S = 4 sin²(πfτ_RF) S_φRF（T14 §3.6、`rf_reference` ＋ `delay_buffer` [code]）と、
   遅延線・配線の温度による遅延のドリフト（2π f_ref δτ(T)、高い f_ref ほど効く）。プログラムの「零点」（φ_d,k = 0 での Θ_k）はポートごとに未知の定数なので、BHD か GBS のデータで 1 回較正する必要がある。
   直交偏光の PM ＋ TEC の構成では、probe と squeezed 光の偏り（T10 §0-1、D78 の M6-04）がこの零点に足される。
5. **プログラムの切り替えはステップではなくランプで** [AE]: ヘテロダインの誤差は sin 型なので、段差 Δ が |Δ| < π なら正しい零点に引き込む。ただし大きな段差ではロックが非線形の範囲を通り、
   アクチュエータが飽和すると整定が遅れる。段差 ≤ π/2 なら引き込みは安全だが、整定時間は約 1/(2π f_u) × 数倍（f_u: ループの交差周波数）。
   φ_d をループの帯域より遅いランプで動かせば、位相は常にロックの線形の範囲にある。どちらが速いかはシミュレーションで比べる（§4.3）。

### 3.2 次のベンチ実験（レーザー ＋ 光路 ＋ BHD、GBS なし。Ref A・Ref C のような位相の安定な squeezing を狙う）

- この実験にはプログラムの必要はない。移せるのは「**LO の位相の設定点を、誤差信号の傾きに依らない形で動かす**」ことだけ [AE]:
  ヘテロダインのロック（Ref A の lock 2 は f_sh = 1 MHz のビート、T14 §3.4）では、復調の位相 φ_d2 を動かすと測る直交位相の角度 Δ が動く。
  **squeezing と反 squeezing の間を ESA のゼロスパンで測る角度の掃引**（Ref A・P7 の Fig. 3 のような位相の走査の代わりに、ロックしたまま角度を段階的に変える）に使える。
  ただし Ref A 型の lock 2 は squeezer を通った probe を使うので、零点は tan ε = (A−/A+) tan δu′ の形で A± に縮められ（T14 §3.4）、φ_d2 と Δ の関係は 1:1 ではない（u′ と v の 2 つの角度が絡む）。1:1 になるのは Liu 型のように probe が squeezer を通らない場合。
- Ref C 型の dither ロックでは、上と同じ理由で復調の位相はロック点を動かさない。

## 4. シミュレータへの対応

### 4.1 既存のブロック・オーバーレイ・プリセットで表せること [code]

| 要素 | 既存のもの | 備考 |
|---|---|---|
| ポートごとの Liu 型ロック（PZT ＋ TEC の積み重ね） | `presets/future/F_common.yaml` の `liu0`/`liu1` → `adc0`/`adc1` → `lock0.stack`/`lock1.stack`（`pid0f` → `sup0` → `pid0s`、アクチュエータ `dac0f` → `dl0f` → `pzt0`、TEC `tec0`） | ポート 0 とポート k の 2 ポート組（D13） |
| ヘテロダインの probe 検出 | `presets/axes/probe/detect_heterodyne.yaml`（`liu_probe_interferometer` の `mode: heterodyne` ＋ `demodulator`（`dem0`/`dem1`）、`liu_f_shift_hz` = 80 MHz [AE]、`liu_demod_lpf_hz` = 300 kHz [AE]） | 名目は `phi_d_rad: 0.0`（I = K cos X、X = π/2 にロック） |
| 設定点（復調の参照の位相） | `demodulator` の実行時パラメータ `phi_d_rad`（`dem1.phi_d_rad`）。`SimController.set_params({"dem1.phi_d_rad": Δ})` で変えられる | **1 回の `acquire` の中では一定**（SimController は実行ごとに仕様を作り直す） |
| 参照の RF 位相の雑音と配線の遅延 | `rf_reference`（白色の位相雑音 ＋ 白色の周波数雑音）＋ `delay_buffer`（τ_RF）→ `demodulator` の `in_phase`（トグル `rf_reference_phase_noise`） | ポート間で共通の源から配るなら、1 つの `rf_reference` を 2 本の `delay_buffer` で分ける |
| 指標 (ii) の目標値 | `metric_ii`（`mii`）の実行時パラメータ `setpoint_rad`（Θ_k,set、ポートごとに `{k: 値}`） | プログラムした値を Θ_k,set に入れれば、指標 (ii) は「プログラムからのずれ」になる |
| ロックの真の誤差の基準 | `lock_error`（`path.lock`）の実行時パラメータ `L1.setpoint_rad`（F_common の名目 π/2） | φ_d をずらしたら、ここも π/2 − Δ にそろえないと `lock.L1.err` と `sup1` のスリップの判定がずれる |
| ロックの安定判定 | `loop_stability_check`（F_common に入っている） | ヘテロダインの復調の LPF の遅れは判定に入っていない（`detect_heterodyne.yaml` の注記: liu_demod_lpf_hz ≫ f_u を保つ） |
| スリップの数え方 | `lock_supervisor` の `slip_center_rad`（D83） | プログラムで零点が動くと、スリップの中心もずらす必要がある（下の 4.3 の 3） |

これで「**一定のプログラム値 Δ_k を入れたときの定常の指標 (ii)**」は、新しいブロックなしで試せる:
`detect_heterodyne` ＋ オーバーレイで `dem1.phi_d_rad = −Δ_k`・`mii.setpoint_rad = {1: Δ_k}`・`path.lock.L1.setpoint_rad = π/2 + Δ_k`（トポロジーの `L1` の設定点。φ_d = −Δ_k ならロック点は X = π/2 + Δ_k）（符号はロックの符号と `coupler_sign` で決まるので、最初に小さい Δ で確かめる）。
Δ_k を `!sweep` で振れば、目盛りの直線性（Θ_k の平均 − Δ_k）と、Δ_k による雑音の変化（ロックの線形の範囲・傾きは変わらないはず）を確認できる。

### 4.2 オーバーレイでは足りないこと

- **時間で変わる設定点がない。** `phi_d_rad`・`setpoint_rad`・PID の `setpoint_v` はどれも実行時パラメータで、1 回の実行の中では一定。`SimController` は `acquire` ごとに初期状態から走らせ直すので、「ロックを保ったまま設定を切り替える」ことは、インターフェースの呼び出しの組み合わせでも表せない [code]。
- 信号の補助ブロック（`ramp`: 鋸歯状波・三角波、`sine_injection`: 正弦）には**ステップ・区分一定の波形がない** [code]。
  `demodulator` の `in_phase` に `ramp` を入れればランプ状の切り替えは作れるが、(a) トグル名が `rf_reference_phase_noise` で意味が違う、(b) 周期的にしか動かない、(c) 指標 (ii) の Θ_k,set は同時に動かない（ずれが全部「誤差」に数えられる）。
  鋸歯状波を半周期遅らせて引き算すれば方形波になるが（`ramp` ＋ `delay_buffer` ＋ `adder`）、遅延バッファの最初の半周期が過渡になり、ポートごとの任意の表（30 組 × K ポート）は作れない。
- 指標 (ii) の Θ_k,set を時間で動かす口がない。`lock_error` の `<id>.setpoint_rad`（`lock.L<k>.err` の基準。`lock_supervisor` の `true_in` が読む）も同じく実行中は一定。
- ランナー（`runner/run.py`）にはチャンクの間のフック・パラメータの予定表・事象の仕組みがない。カーネルは毎サンプル P を読むので、
  Python のスクリプトで `build` → `compile_loop` → `run_compiled` をチャンクごとに呼び、間で `params.set(...)` すれば段差は作れる。
  ただしこれは YAML の設定と記録（D55-8）の外の手順で、プリセットとして残せず、`SimController` 経由で実機と同じ手順にもならない。
- 試験用の `pathtest_table`（`path/testing.py`、表をスロットに流す。インパルス・ステップの予定を持つ）は `ensure_loaded` で読み込まれない試験専用のブロックで、プリセットからは使えない。
  N6-1 はこれを本番用にしたもの、と考えてもよい。

### 4.3 新しいブロックの提案

**N6-1 `phase_program`（CTRL、時間ループ）**
- 役割: ポートごとのプログラム値の表を時刻に沿って出す。区分一定（ステップ）と、傾きの上限付きの遷移（ランプ）を選べる。
- 入力: なし（時刻 `t` だけ）。
- 出力: `rf.<name>.port<k>` [rad]（各ポートの参照の位相。`demodulator` の `in_phase` へ）、`ctrl.<name>.target<k>` [rad]（同じ時刻の Θ_k の目標値。指標 (ii) と監視へ）、
  `ctrl.<name>.index`（今の設定の番号）。
- 構造パラメータ: `ports`（k の並び）、`table_rad`（設定の数 × ポート数の表、`!fn` か YAML の配列）、`dwell_s`（各設定の保持時間）、`transition: step | slew`、`sign`（φ_d と Θ の符号、±1）。
- 実行時パラメータ: `slew_rad_per_s`（`slew` のとき）、`offset_rad`（ポートごとの較正の零点）。
- トグル: `enabled`（False なら全ポート 0）。
- オーバーレイでは足りない理由: 上の 4.2。表を持った区分一定の駆動は既存のどの波形でも作れない。
- 検証の案: 表の値がそのまま出力されること、`slew` の傾きの上限、`enabled` の切り替え。

**N6-2 `metric_ii` の目標値を信号から読む口（既存ブロックへの後方互換の追加。新しい kind ではない）**
- 追加の案: 任意の構造パラメータ `setpoint_slot`（null なら今の実行時パラメータ `setpoint_rad`）。読んだ値を Θ_k,set として Θ_k(t) = wrap_π(θ_sq,k − θ_sq,0 − Θ_k,set(t)) を作る。
- 入力: `ctrl.<program>.target<k>`。出力: 今と同じ `metric.<block>.port<k>`。
- 理由: プログラムを切り替えた瞬間の段差を「誤差」から外し、整定の過程だけを測るため。`lock_supervisor` の `slip_center_rad` にも同じ口が要る（零点が動くとスリップの中心も動く）。
- `src/` の変更なので、このノートでは提案だけ（Procedure 3 の grill-me で決める）。

**N6-3 `step_response`（MEAS、測定段のフック。ESA と同じく時間ループの後に記録から計算）**
- 入力: `ctrl.<program>.index`・`ctrl.<program>.target<k>`、`metric.<mii>.port<k>`（Θ_k(t)）、ロックの真の誤差 `lock.L<k>.err`、アクチュエータの駆動 `ctrl.sup<k>.u`。
- 出力（`RunResult.measurements`）: 切り替えごとの立ち上がり時間、行き過ぎ、整定時間（|Θ_k − 目標| ≤ ε_set に入って出なくなる時刻、ε_set はパラメータ）、
  定常の偏り（目盛りの誤差）、整定後の σ、アクチュエータのレンジの消費、スリップの有無、**ポート間の一致**（同じ時刻に切り替えたポートどうしの整定時間の差と、整定後の Θ_k − Θ_j の偏り）。
- 理由: 既存の測定（`measure.esa`・`measure.pulse_sample`・指標の分解）は定常の統計で、切り替えの事象ごとの時間の指標を出すものがない。

### 4.4 試験の設計（Procedure 3・4 で実行する案）

1. **定常の目盛り（新しいブロックなしで今すぐ）**: F_liu20 ＋ `detect_heterodyne`、Δ ∈ {−π/2, −π/4, 0, π/4, π/2}、10 seed。Θ_k の平均 − Δ（目盛りの誤差）と σ が Δ に依らないことを見る。DC の設定点（`homodyne_dc` の `setpoint_v`）で同じ Δ を入れた場合と並べ、K の揺れ（probe のパワー・可視度の掃引軸）への依存を比べる。
2. **ステップ応答（N6-1〜N6-3 が要る）**: 段差 π/8・π/4・π/2、`step` と `slew`（傾き = 0.1・1・10 × 2π f_u rad/s）。整定時間を f_u（`f_u_hz` の軸）と復調の LPF（`liu_demod_lpf_hz`）に対して出す。
3. **ポート間の一致**: 2 ポート組で同時に切り替え、整定後の Θ_1 の偏りと、参照の配線の遅延差 τ_RF・RF 源の位相雑音（`rf_reference` の掃引）による残差。32 ポートは D13 のとおり必要になったときに。
4. **GBS 側の指標**: 切り替え 1 回あたりの「使えない時間」（整定時間）× 切り替えの頻度を、試料の収集時間（Jiuzhang は 1 組 200 s [Zhong p.7]）と比べる。

## 5. UNVERIFIED と未解決の問い

- 1450 nm の参照光、参照信号の周波数・種類、ロックの帯域・残差・アクチュエータ、遅延線の型と分解能: すべて SM にあると思われるが未読（UNVERIFIED）。**SM の PDF があれば §2 の表の [SM] の行を埋められる。**
- Fig. 1(e) の 0.06 rad がどの量の、どの時間スケールの統計か（in-loop の誤差か、ポート間の真の位相差か）。
- 「参照信号」が光（1450 nm の参照光の位相変調・ビート）なのか電気（RF）なのかは本文からは決められない。「electric delay lines」とあるので、遅延をかけているのは電気信号と読んだ [AE]。
- 2 パスの源の "self-phase stable" は源の内部の話で、源から干渉計までの光路のロックとは別。どの区間を能動的にロックしているかは本文にない。
- 自分たちの側の未確定: Moku:Delta でポートごとに独立に位相を設定できる復調の参照が何本あるか（T16 §4.5 の表になし）、32 ポートで参照を何台の機器から配るか（D15・D13）。
- 符号: `demodulator` の φ_d と Θ_k の符号は `coupler_sign`・ロックの符号・アクチュエータの符号の積で決まる。最初の試験（§4.4 の 1）で確かめてから N6-1 の `sign` を決める。

## 6. 文献

- H.-S. Zhong, Y.-H. Deng, J. Qin, H. Wang, M.-C. Chen, L.-C. Peng, Y.-H. Luo, D. Wu, S.-Q. Gong, H. Su, Y. Hu, P. Hu, X.-Y. Yang, W.-J. Zhang, H. Li, Y. Li, X. Jiang, L. Gan, G. Yang, L. You, Z. Wang, L. Li, N.-L. Liu, J. J. Renema, C.-Y. Lu, J.-W. Pan,
  "Phase-Programmable Gaussian Boson Sampling Using Stimulated Squeezed Light," Phys. Rev. Lett. **127**, 180502 (2021). DOI 10.1103/PhysRevLett.127.180502（PDF の 1 ページ目で確認）。本文 p.1–7、Fig. 1・4 を使用。
- 同論文の Supplemental Material（文献 [25]、http://link.aps.org/supplemental/10.1103/PhysRevLett.127.180502）: **未読**。
