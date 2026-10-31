# Phase 5 物理レビュー: M1 の Ref D 系（新しい文脈のレビュー担当）

- 日付: 2026-10-31 ／ 対象: dev @ 66ccb3a（`git merge --ff-only dev` 後）
- 範囲: WP-B2（`laser_source` の DX1、`segment_thermal`、`segment_mech_noise`、`pd_offset_drift`）、WP-B3 `dc_fringe`、
  WP-B4（`pid` の `dsc1` 離散化、`adc`・`dac` の量子化、`pzt_fiber_stretcher`、`hv_amp_914`、`lock_supervisor`、`thermal_actuator`）、
  WP-B7（`segment_jones`、`polarizer_projection`、`per_fringe`、`per_budget`、`coil_birefringence`）、
  WP-F5（`devices/dsc1.yaml`、`devices/pdb450c.yaml`、`hardware/mapping.py` の DSC1 写像）、WP-I2（`presets/reference/refD*.yaml`）、
  テスト VT-13・14・21・24・38、`tests/integration/*refD*`。
- 読んだもの: AGENTS.md、DECISIONS.md（D02・D04・D06・D34–D36・D49・D51・D53・D55–D57・D61）、`notes/conventions.md`、PLAN §5.2・§6.4・§6.5・§7.4、
  T15（主）、T16、T05・T08・T09・T13 の該当節、マニュアル（`refs/manuals/`: thorlabs-dsc1.pdf、thorlabs-pdb.pdf、914-915B-Specifications.pdf（スキャン画像を目視）、
  pnh1550r5a2.pdf、Opneti の In Line Polarizer）。
- **`blind/` は読んでいない。** 実機の挙動は一切参照せず、資料と物理だけで判定した（D04）。ソースコードは変更していない（`reviews/` だけ）。
- 証拠のスクリプト: `reviews/scripts/m1_refD_*.py`（`uv run python reviews/scripts/<name>.py` で再現）。
- 実行したテスト: VT-13・14・21・24、`tests/integration/*refD*` → 37 passed / 1 skipped（ブラインドの置き場所）。VT-38 → 19 passed。

## 0. 結論（要約）

- **符号・単位・係数 2・機器の値に誤りは見つからなかった。** DSC1 の p.20 の単位の写像（D53）、e = S − v の変換、conv §4（δOPL に ξ_eo を掛けない）・§5（腕と LO/pr の対応、カプラの符号 s）・§9 は、コードとプリセットで正しく実装されている。
- **D61-3 の判定: 観測は正しい（独立のモデルで 1.245 LSB を再現、WP-I2 は 1.243 LSB）。** 原因は「しきい値に乗ったリレー」の記述関数で説明でき、リミットサイクルはループの位相交差周波数（硬い 915B で ≈ 17 kHz）の振動。T15 §2.11・§7 test 10 と D56-7 の「±LSB/2（±6.3 mrad）」は、遅延のない理想のループでも瞬時の誤差の上限としては成り立たない（成り立つのは rms）。修正文案を §2 に示す。
- **D61-2（遅延の勘定）は正しい**（−180° は 60.3 kHz、DSC1 p.2 の > 58 kHz を満たす）。ただし 1.5 µs のうち 0.5 µs は 914 の ZOH 離散化の半サンプルが担っており、プリセットの説明と違う（M1D-05）。
- VT-14 の最低帯域の 9 % 高い値は、窓（Hann ＋ 区間ごとの一次トレンド除去）の偏り +3.5 % と統計のゆらぎ（1.5σ）の和。5 seed の平均は窓の期待値と 0.8 % で一致する。VT-14 の許容幅は帯域平均に対して約 7 倍緩い（M1D-04）。
- f_sim × 2 の変化は −0.3 %（seed のばらつき 0.3 % と同程度）で、PLAN §6.4 の 3.3 % より十分小さい。収束として受け入れてよい。
- refD のどの値も実機の挙動に合わせた形跡はない。ブラインドの置き場所（`test_refD_blind_procedure.py`）に目標値は含まれない。

## 1. 指摘の一覧

| ID | 重大度 | 場所 | 内容 | 提案 |
|---|---|---|---|---|
| M1D-01 | Medium | `DECISIONS.md:370`（D56-7）、`DECISIONS.md:421`（D61-3）、`notes/dc-fringe-lock-ref-d.md:286-287, 342, 380`（T15 §2.11、§6.1、§7 test 10）、`tests/vt/test_vt_24_refD_analytic.py:17-20, 255-257, 275-284` | DAC 量子化のリミットサイクルの上限「±LSB/2 = ±6.3 mrad」と「12.7 mrad は常に成り立つ」は誤り。VT-24 test 6 の説明「±LSB/2 is exact for a static actuator」も誤り（理想のループでも固定点で最大 \|e\| 0.67 LSB、傾斜で帯 1.94 LSB）。test 6 の固定点の判定 `max|e| ≤ 0.6 LSB` は動作点しだいで境界上（独立モデルで最大 0.60 LSB）。 | §2 の文案で D56-7・T15 §2.11・test 10 を置き換える。VT-24 test 6 は (a) 固定点で帯 ≤ 1 LSB、(b) rms ≤ (√2/π)·k·q·\|P(i2πf₁₈₀)\|·1.2（硬い 915B で 5.7 mrad）、(c) 傾斜で帯 ≤ 2 LSB、で判定する（§2）。 |
| M1D-02 | Medium | `presets/reference/refD.yaml:94-95`（`f_res_hz`・`q_res` の掃引）、D61-4 | D53 のゲインのまま 915B の共振を掃引すると、格子 16 点のうち 6–7 点（f_res = 5 kHz の全 Q、11.45 kHz の Q ≥ 9.65 または ≥ 31）が線形に不安定。そのような点は発振（フリンジ上を滑る）を「予測」として出してしまう。 | 掃引点ごとに解析的な安定判定（Nyquist、`m1_refD_stability_grid.py` の方法）を前処理として行い、不安定な点は「D53 のゲインでは不安定」と印を付けて予測帯から分けて記録する。「実機はこのゲインでロックする」という定性的な事実で格子を絞ってよいかは D04 に関わるので、調整役がユーザに確認する（Open）。 |
| M1D-03 | Medium | `src/cranberry/blocks/loss/coil_birefringence.py:55, 123`、`src/cranberry/blocks/loss/per_budget.py:24, 92-99`（D61-5 の未処理分） | YAML から読むと `bend_loss_db` → `bend_loss_lin`、要素の `per_db` → `per_lin` に変換され、両ブロックとも組み立てで失敗する（`m1_refD_lin_names.py`: `missing required parameter(s) ... bend_loss_db`、`give exactly one of theta_deg, per_db, x, ...`）。保護リストの「YAML の読み込みと Python リストの等価性」に反する。 | `polarizer_projection` と同じく `_lin` の名前を受け付ける（`bend_loss_lin` は損失の因子 ≥ 1、要素の `per_lin` は 10^(PER/10)）。YAML 経由の組み立てのテストを 1 つ加える。 |
| M1D-04 | Low | `tests/vt/test_vt_14_refD_free_running.py:11-13, 69-76` | 許容幅 3/√N_avg（= 0.249）は 1 ビンの統計誤差で、帯域平均の誤差（最低帯域で 1σ ≈ 3.5 %、最高帯域で 0.3 %）の 7–80 倍。説明の「per band average, which is the stricter statistic」は逆。さらに最低帯域（200–500 Hz、ビン 4–8）では窓と一次トレンド除去が f⁻² の項に +3.5 % の決まった偏りを与える（`m1_refD_vt14_bands.py`: E[Welch]/model = 1.035）。9 % = 偏り 3.5 % ＋ 統計 1.5σ。 | 期待値を窓の核で畳み込んだ E[Welch]（スクリプトの `expected_welch`）と比べ、許容幅を帯域平均の 3σ = 3·√(1.5/(N_avg·N_bins)) にする（または最低 3 ビンを除く）。説明を直す。 |
| M1D-05 | Low | `presets/reference/refD.yaml:17-20, 85, 170-173, 180` | D61-2 の 1.5 µs（2 点平均 0.5 µs ＋ 500 kHz の ZOH 1 µs）のうち、pid・dac の 2 サンプル保持が担うのは平均 0.5 µs で、残り 0.5 µs は ZOH で離散化した 914 の半サンプル（y[n] = C x[n]）が担う（`m1_refD_delay_accounting.py`: 914 の余分な遅延 0.500 µs、`tustin` で 0.007–0.026 µs、`poles: false` で 0）。説明は「adc / pid / dac が持つ」としており、914 の `poles` を切るか `tustin` にすると、ループ遅延が黙って 0.5 µs 短くなる。180 行目の「input held by the DAC (zoh exact)」も、入力は 100 kHz の LPF と線形補間を通った列なので正確ではない。 | 説明を「0.5 µs（2 点平均）＋ 0.5 µs（2 サンプル保持）＋ 0.5 µs（914 の ZOH の半サンプル）」に直す（VT-24 test 10 の注記と同じ内容）。914 の `poles` を切る・`tustin` にするときは `dsc1_explicit_delay_s` を 1.0 µs にする旨を書くか、`tau_s` の式に 914 の離散化の項を入れる。 |
| M1D-06 | Low | `presets/reference/refD.yaml:64` | Opneti の ER は「30 dB typ、≥ 28 dB、コネクタ付きは −2 dB」。Ref D は FC/APC なので typ 28 dB・最小 26 dB。掃引の上端 30 dB は裸ファイバの typ。名目 28 dB の説明「≥ 28 dB」もコネクタなしの最小値を指している。 | 掃引を [26, 28]（または上端 30 を「コネクタなしの typ、楽観側」と明記）、名目の説明を「コネクタ付きの typ 28 dB」に直す。 |
| M1D-07 | Low | `notes/dc-fringe-lock-ref-d.md:189, 227, 412` | T15 の `ki_entry` の既定「per_second（p.12）」と共振の範囲（10–60 kHz、Q 5–30）は D53・D51-11 で置き換え済みだが、ノートに注記がない。 | ノートに「D53 で p.20（per_sample）が正しいと確定」「D51-11 で 5–60 kHz、Q 3–100」の注記を加える（M1D-01 の T15 の修正と同時に）。 |
| M1D-08 | Info | `presets/reference/refD.yaml:61` | `k_eps_rad_per_m: 4.64e+6` の説明は n = 1.468 と書くが、値は n = 1.46（ξ_eo·2π·1.468/λ0 = 4.667×10⁶）。ひずみの位相係数には位相屈折率を使うのが正しいので値は妥当、説明だけ違う。位相の PSD は直接与えるので結果には効かない（ひずみの群遅延に 0.6 %）。 | 説明を「n_eff = 1.46（位相屈折率）」に直す。 |
| M1D-09 | Info | `presets/reference/refD.yaml:114` | 腕 1 の説明「nominal 8.4 + 0.5 = 8.9 ≈ 9 m」は、腕 2 にも 0.5 m があるので ΔL = 8.4 m（巻線分）。9 m は D53 の掃引の中心なので値は問題ない。 | 説明を「ΔL ≈ 8.4 m（巻線）、掃引の中心は D53 の 9 m」に直す。 |
| M1D-10 | Info | `presets/reference/refD.yaml:333` | γ_coh を `BandLedger`（`gamma_coh_from_ledger`）ではなく白色 FM だけの式で与えている。flicker FM（`h_m1_hz2` の上端 5×10⁸ Hz²）の帯域外の分 ≈ 1.2×10⁻⁴ rad² が抜ける（白色の 0.022 rad² の 0.5 %）。 | 影響は無視できる。将来ほかの高周波の雑音源を入れるときのため、ledger から作る形に揃えるとよい。 |
| M1D-11 | Info | `tests/integration/test_refD_convergence.py:75-85` | f_sim × 2 の確認はレーザー雑音が支配する帯域分散だけを見ている。ループの動特性（位相交差、リミットサイクル）は見ていない。独立に確かめた限り、seed のない量（量子化のリミットサイクルの rms）は 3.171 → 3.165 mrad（硬い）、12.06 → 11.67 mrad（共振、−3.3 %）。 | 決定的な量（例: `sensitivity_z` の 20 kHz の値、またはリミットサイクルの rms）を f_sim × 2 の確認に加える。共振ありのリミットサイクルは 1/3 の許容幅の境目にある。 |

## 2. D61-3（量子化のリミットサイクル）の判定

### 2.1 独立の検証

`reviews/scripts/m1_refD_dac_limit_cycle.py` は cranberry のコードを使わない、資料だけから組んだ準連続のモデル（dt = 0.1 µs）。
v = −V0 sin e → 120 kHz LPF → 1 MHz の ADC（2 点平均、16 bit）→ 2 µs ごとの PID（p.19 の式、e = S − v、D53 → p.20）→ 12 bit の DAC（2 µs 保持）→
100 kHz LPF → 純遅延 4.2 µs（5.7 − 1.5、D61-2）→ 914（29.5 / 112 kHz）→ 915B（硬い、または 20 kHz・Q 10）。アナログ部は 10 MHz の双一次変換（ZOH の半サンプルを入れない）。
V0 = 3.33 V、k_act = 6.491 rad/V、g = 0.99956、1 LSB = 12.67 mrad。

| 場合 | 帯（ptp）[LSB] | 最大 \|e\| [LSB] | rms [mrad] | 主な周波数 | WP-I2（VT-24 test 6） |
|---|---|---|---|---|---|
| 理想（遅延・フィルタなし、硬い）固定点 | 1.000 | 0.674 | 5.94 | — | — |
| 理想、傾斜 0.4 rad/s | 1.940 | 0.996 | 5.16 | 250 kHz（最速の切り替え） | — |
| refD、硬い、固定点（9 点を 1 LSB にわたって） | 0.65–0.97 | 0.49–0.60 | 1.7–4.6 | 17.4 kHz | 0.740 / 0.533 |
| refD、硬い、傾斜 0.4 rad/s | **1.245** | 0.634 | 2.90 | 17–20 kHz | **1.243** / 0.628 |
| refD、共振 20 kHz・Q 10、固定点 | 3.14 | 1.63 | 13.6 | 16.8 kHz | 3.81（帯）、13.96 mrad |
| refD、共振、傾斜 0.4 rad/s | 4.44 | 2.23 | 12.3 | 16.7 kHz | — |

### 2.2 記述関数による説明

DAC の入力 u が 2 つのコードの間（四捨五入のしきい値）にあるとき、量子化器は振幅 q/2 のリレーとして働き、その記述関数は N(A) = 2q/(πA)（A → 0 で無限大）。
したがって積分のあるループでは必ず N(A)·L(iω) = −1 の解があり、リミットサイクルは **ループの位相交差周波数 f₁₈₀** で起き、振幅 A = 2q\|L(i2πf₁₈₀)\|/π（\|L\| が小さいほど小さい）。
光の位相に出る基本波の振幅は (4/π)(q/2)·k·g·\|P(i2πf₁₈₀)\|（P: DAC 出力から光位相までの経路、DC で 1）なので:

| | f₁₈₀ | \|L(f₁₈₀)\|（利得余裕） | \|P(f₁₈₀)\| | 基本波の ptp（予測） | 独立シミュレーション |
|---|---|---|---|---|---|
| 硬い 915B | 17.5 kHz | 0.071（22.9 dB） | 0.83 | 1.06 LSB | 傾斜で 1.245 LSB、PSD のピーク 17.4 kHz |
| 共振 20 kHz・Q 10 | 15.3 kHz | 0.229（12.8 dB） | 2.03 | 2.58 LSB | 3.1–4.4 LSB、16.7 kHz |

傾斜では動作点がしきい値を横切るたびにこの最大の場合が現れるので、帯は固定点（しきい値から離れた動作点）より広くなる。高調波が加わるので、シミュレーションの ptp は基本波だけの予測より 15–70 % 大きい。

### 2.3 雑音があるときにリミットサイクルが残るか（`m1_refD_quant_in_noise.py`、cranberry の refD_nominal、同じ seed で DAC 量子化のオン・オフ）

| 構成 | 残差の rms（オン / オフ） | 量子化による増分（二乗差） |
|---|---|---|
| ΔL = 9 m、硬い | 31.49 / 31.45 mrad | 1.6 mrad（0.13 LSB、T16 §5 の白色の見積もり 1.5 mrad と一致） |
| ΔL = 9 m、共振 | 31.74 / 31.50 mrad | 3.9 mrad（0.31 LSB） |
| ΔL = 0、硬い | 2.81 / 0.19 mrad | 2.8 mrad（0.22 LSB） |
| ΔL = 0、共振 | 13.48 / 0.19 mrad | **13.5 mrad（1.06 LSB）** |

ΔL = 9 m ではレーザー雑音（≫ 1 LSB）が量子化器をディザするので白色の量子化雑音として振る舞う。ΔL → 0 ではリミットサイクルがそのまま残り、共振があると残差を支配する。
**Ref D のブラインド予測では、ΔL の小さい側・共振の強い側で DAC の量子化が最大の項になりうる**ので、アブレーション（`dsc1.dac.quantization`）の結果を予測と一緒に記録すること。

### 2.4 判定

- D61-3 の数値（傾斜で 1.22–1.25 LSB、固定点で 0.74 LSB・最大 \|e\| 0.53 LSB、共振で約 ±2 LSB）は**正しい**。シミュレータの不具合ではなく、遅延のある積分ループと量子化器の組み合わせの物理。
- 「12.7 mrad は常に成り立つ」（D56-7）と「振幅 ±1/2 LSB」（T15 §2.11・test 10）は**誤り**。理想のループでも、2 コードのリミットサイクルの瞬時の誤差は kδ と −k(q − δ)（δ: コード間の位置）なので、固定点で最大 \|e\| は q/2 から q まで、傾斜では帯が 2 LSB に近づく。±LSB/2 で抑えられるのは **rms**（デューティ比 d で k q √(d(1 − d)) ≤ kq/2）だけで、それもアクチュエータに共振があれば成り立たない。

### 2.5 D56-7 の文案（D61-3 と統合して置き換え）

> 7. DAC の量子化（外乱 ≪ 1 LSB のとき）: 積分のあるループでは、隣り合うコードの間のリミットサイクルになる。遅延も動特性もない理想のループでも、瞬時の誤差は ±1 LSB 未満（固定した動作点で帯 ≤ 1 LSB、ゆっくりした傾斜では帯 ≤ 2 LSB）で、±LSB/2 で抑えられるのは rms（k q √(d(1 − d)) ≤ kq/2）だけ。遅延とアクチュエータの動特性があると、リミットサイクルはループの位相交差周波数 f₁₈₀ の振動になり（しきい値に乗ったリレーの記述関数）、真の位相の ptp ≈ (4/π)·k·q·\|P(i2πf₁₈₀)\|、rms ≈ (√2/π)·k·q·\|P(i2πf₁₈₀)\|（P: DAC 出力から光位相まで、DC で 1）。refD（D53 のゲイン）: 硬い 915B で f₁₈₀ ≈ 17 kHz、rms 1.7–4.6 mrad、帯 0.65–0.97 LSB（固定点）・≤ 1.25 LSB（傾斜）。名目の共振（20 kHz、Q 10）で rms 12–14 mrad、帯 3–4.4 LSB。外乱が LSB より十分大きい（ΔL = 9 m のレーザー雑音）ときはディザされ、増分は 1.6 mrad（硬い）/ 3.9 mrad（共振）。（Phase 5 M1 のレビューで確認、`reviews/phase5_M1_refD.md` §2）

T15 §2.11 の 287 行目と §6.1 の表（342 行目）、§7 test 10（380 行目）も同じ内容に直す。VT-24 test 6 は M1D-01 の提案の判定に置き換える。

## 3. 機器の値の照合（マニュアルの印刷頁）

DSC1 = User Guide Rev. B（2025-11-03）、PDB = Operation Manual v2.4（2019）、EO = 914-915B-Specifications.pdf（PDF 頁、スキャン）。

| 機器 | 値（コード・プリセット） | 資料 | 判定 |
|---|---|---|---|
| DSC1 | ADC 1 MHz、2 点平均（実効 500 kHz） | p.1 | 一致 |
| DSC1 | ADC 16 bit、±4 V、入力終端 1 MΩ | p.2 | 一致 |
| DSC1 | 入力帯域 120 kHz typ（次数の記載なし → 1 次は仮定） | p.1 | 一致（次数は UNVERIFIED と明記済み） |
| DSC1 | 入力 RMS 雑音 0.3 mV（100 Hz–250 kHz）、平均 −120 dB V²/Hz、ピーク −105 dB V²/Hz | p.2（注 c） | 一致。0.3 mV は 500 kHz の標本の標準偏差として正しく使われている（ADC は 1 MHz、PID は 1 つおきに読むので標本ごとに独立）。−120 dB V²/Hz × 250 kHz = 0.5 mV で 0.3 mV と食い違うので両方を掃引しているのは妥当 |
| DSC1 | DAC 12 bit、±4 V、出力インピーダンス 220 Ω、負荷で R/(R + 220 Ω) | p.2（注 b） | 一致。1 LSB = 1.953 mV → 12.67 mrad |
| DSC1 | 出力帯域 100 kHz（LPF） | p.1、p.2（注 d） | 一致 |
| DSC1 | PID 更新 500 kHz、h = 2 µs | p.2、p.20 | 一致 |
| DSC1 | e = S − v | p.18 §5.4.1 | 一致。`kp = −p`、`ki = −i/h`、`kd = −d·h` で conv §9 の e = v − S に変換。負帰還の向きを独立に確認（腕 1 を伸ばす → φ_D 減少 → v 増加 → e_DSC1 減少 → u 減少） |
| DSC1 | u = K_p e + h Σ K_i e + (K_d/h)(e[n] − e[n−1])（bumpless、微分のフィルタなし） | p.18–19 | 一致（`first_order(0, 0, ki, 0, h)` = I += K_i h e、e[n] を含む） |
| DSC1 | 入力値 I = K_i·2×10⁻⁶、D = K_d/2×10⁻⁶ | p.20（Ziegler–Nichols の段落）、D53 | 一致。p.12 の「Ki in 1/s、Kd in s」は使わない（D53）。D53 の設定 → **K_p = −2.41×10⁻⁴、K_i = −703 s⁻¹、K_d = −2.216×10⁻⁸ s** を独立に再計算して一致 |
| DSC1 | P・I・D の符号は連動、\|入力値\| 10⁻⁵–10⁵ | p.19、p.6 | 一致 |
| DSC1 | S・O は ±4 V、CV の上下限は設定画面 | p.6、p.9 | 一致（D35: −4〜+4 V） |
| DSC1 | −180° typ 60 kHz、> 58 kHz | p.1、p.2（注 a） | 一致。連続の当てはめ 5.67 µs → 60.00 kHz、5.7 µs → 59.77 kHz、純遅延 8.3 µs → 60.24 kHz、二重に数えると（5.7 + 1.5 µs）50.35 kHz（`m1_refD_delay_accounting.py`）。プリセットの離散の鎖で 60.3 kHz（VT-24 test 10）。D61-2 は正しい |
| DSC1 | ランプ 10–400 ms、PEAK の変調 100 Hz–100 kHz（100 Hz 刻み）、K は V/s、トレースの点 2 µs | p.7、p.2・p.8、p.8、p.6 | 一致 |
| DSC1 | 積分器の anti-windup | 記載なし（p.6・p.9 は CV のクリップだけ） | UNVERIFIED として掃引軸（`none`/`clamp`）にしているのは妥当 |
| PDB450C | 利得 10³ V/A（高 Z 負荷）、50 Ω で 1/2 | p.19（注 1） | 一致（DSC1 の 1 MΩ 入力、D53） |
| PDB450C | RF OUTPUT の振幅 ±4.6 V（PDB450、10³–10⁶） | p.14 | 一致 |
| PDB450C | DC オフセット < ±15 mV、工場で 23 °C に零、暖機 約 5 分 | p.19、p.8 | 一致 |
| PDB450C | CMRR > 25 dB（typ > 30 dB） | p.19 | 一致。`eps_r` の ±0.056 = 10^(−25/20)（振幅比、電気の dB）は妥当 |
| PDB450C | 帯域 DC–150 MHz（10³） | p.19 | 一致 |
| PDB450C | 感度 1.0 A/W（最大、typ）、1549 nm で 1.03 A/W（典型曲線） | p.19、p.15 | 一致（曲線は 1570 nm 付近で約 1.04 の最大。掃引 [1.0, 1.03] は妥当） |
| PDB450C | CW 飽和 4.5 mW @ 1550 nm（差に対して）、損傷 20 mW | p.19、p.14、p.5・p.7 | 一致（P_det ≤ 3.33 mW） |
| PDB450C | 最小 NEP 68.6 pW/√Hz（DC–150 MHz、10³） | p.19 | 一致 |
| PDB450C | モニタ 10 V/mW、10 V（高 Z）/ 1.55 V（50 Ω）、1 MHz、< 180 µV rms、< ±2 mV、RF 出力 50 Ω | p.14 | 一致 |
| 914 | 入力 ±5 V でフルスイング | EO p.1 | 一致 |
| 914 | Fig. 1: 52 V で平坦、2 極 29.5 / 112 kHz | EO p.2 | 2 極の当てはめは 10 / 20 / 50 / 100 kHz で 49.1 / 42.4 / 24.1 / 11.0 V、画像の読み 50 / 42 / 23–24 / 11 V と読み取り誤差の範囲で一致。Fig. 1 は大信号・振幅のみ（位相なし）なので、小信号の 2 極モデルは [AE] のまま |
| 914 | スルーレート（入力換算）6.7×10⁵ V/s | EO p.2 | 2π·100 kHz·11 V / 10.4 = 6.6×10⁵ V/s で一致 |
| 915B | 8669-1: 38 巻、PM15-U25D、FC/APC（slow にキー）、全レンジの光路長変化 16 µm | EO p.1 | 一致。k_act = 2π·16 µm/(1548.8 nm·10 V) = 6.491 rad/V、ξ_eo は掛けない（「optical path length measurement」、conv §4） |
| 915B | 共振・容量・巻線長 | 記載なし | 掃引（D36、D51-11）で妥当 |
| PNH1550R5A2 | 50:50 ± 1.5 %、PER ≥ 25 dB（コネクタ付きで規定）、過剰損失 < 0.3 dB typ、リード 0.8 m | 仕様書 p.1 | 一致（`kappa1/2`、`theta_coupler_rad` = atan(10^(−25/20)) = 0.05617、`coupler_excess_db`、`L_arm2_m`） |
| Opneti 偏光子 | IL 0.3 dB typ、ER 30 dB typ / ≥ 28 dB、コネクタ付きは IL + 0.3 dB、ER − 2 dB | 仕様書 | IL は帳簿だけ（P_in は偏光子の後で測った値、D02）で影響なし。ER の掃引の上端は M1D-06 |
| Cobrite DX1（SC） | 線幅 80 kHz typ → h0 = 2.546×10⁴ Hz²/Hz、RIN −140 / −150 dB/Hz | D35、T05 | 一致（Δν = π h0、conv §7）。メーカーの原本は refs にないので D35 の値の転記として確認 |

## 4. PLAN §7.4 の点検表

| # | 項目 | 結果 |
|---|---|---|
| 1 | 式がノートと一致 | `dc_fringe` = T15 式 (1)–(7)（ε の s の扱い、𝒱 = 4√(κ₂(1−κ₂)P₁P₂)/P_det、μ の積）。`per_fringe` の μ = (a_pr†a_lo)(t_pr t_lo*)/\|…\| と逆変換 x = \|μ−1\|²/(2(1−Re μ)) は厳密（VT-38 の全場の和と一致）。T15 式 (7) の Δχ の向き（腕 2 − 腕 1）の書き違いは per_fringe の説明どおり。`segment_thermal` の 1 次ホールドの厳密な離散化、ランダムウォークの増分 √(2π²bh)、`laser_source` の √(h0 f_b/2)、`segment_mech_noise` の T08 式 2.7、`pzt` の 2 次モード（Tustin、共振で prewarp）、TEC の厳密な 1 次遅れ、すべて一致 |
| 2 | 符号 | conv §4（伸ばすと ψ 増）、§5（腕 2 = LO、腕 1 = pr、φ_D = φ_LO − φ_pr、s）、§9（e = v − S への変換）。VT-24 test 1（φ_D = φ₀ − k·g·u）、test 8（loop_sign 反転で −π/2）で確認。偏光の Δψ_b の K_T < 0（D51-2） |
| 3 | 2 倍の因子 | 片側 PSD → 標本の分散 S·f_s/2（dc_fringe のショット・電子雑音、914 の雑音、レーザー）。ショット雑音は 2 つの PD の合計電流 2eℜP_det。CMRR → ε_R は振幅比 |
| 4 | 単位 | k_act は 914 の入力電圧あたり（rad/V）、δOPL に ξ_eo を掛けない。D53 の入力値 → 1/s・s は p.20。dB → 線形（D56-1） |
| 5 | FFT の符号 | VT-13 の復調 Σx e^{−i2πfn/f_s} と `freq_response` の z = e^{+i2πfh}（Laplace の s = +iω）で整合 |
| 6 | 収束 | f_sim × 2: −0.3 %（seed のばらつき 0.3 %）、チャンク長: ビット単位で一致（テストが通過）。M1D-11 |
| 7 | トグルで恒等 | `pid.feedback` off → u = clip(O)、量子化 off → 素通し、`per_crosscoupling` off → x = 0、`segment_jones` の全トグル off → U = 1（VT-38）、`resonance` off → D = 1 |
| 8 | 損失・帯域の登録 | 915B の IL は `dc_fringe.t_arm1` だけ（`insertion_loss_lin: null`）で二重なし。偏光子の IL は帳簿のみ。γ_coh（帯域外の可視度）と ledger の lockD の帯域外分散は別の物理量で二重ではない（M1D-10 は情報）。電子回路の雑音は `bands()` なし（D56-5） |
| 9 | 出典の印 | ほぼすべての値に出典。説明の食い違いは M1D-06・08・09 |
| 10 | 保護リスト | 削除なし。YAML の等価性に M1D-03 |

## 5. D04（ブラインド）の確認

- `refD.yaml`・`refD_nominal.yaml` のすべての値をたどった。値は (a) マニュアル（§3）、(b) DECISIONS のユーザの構成の事実（D02・D35・D53: 3.33 mW、10³ V/A、DSC1 の設定値、S = O = 0、±4 V、オシロスコープ 1 MΩ、ΔL の掃引範囲）、(c) ノートの名目値・掃引範囲の中心（T05・T08・T09・T15・T16、D49・D51）だけから来ている。測定した応答や PSD に合わせた値はない。`git log` でもプリセットの変更は WP-I2 の作成と M4 のキーの追加の 2 回だけ。
- `refD_nominal` の `dpsi_b0_rad = π/2`（一般的な点）、`r_T_k_per_s = 1e-4`（T08 §6 の例）、`f_res 20 kHz・Q 10`（T15 の中心）も資料由来で、実機に合わせたものではない。
- `tests/integration/test_refD_blind_procedure.py` は skip され、目標値・測定値を含まない。手順は PLAN §6.5 と一致。`src`・`tests`・`presets`・`devices` に `blind/` の内容への参照はない。
- 注意（M1D-02）: 「実機はこのゲインでロックする」を使って共振の格子を絞るのは、実機の挙動を予測に入れることになりうる。判断は調整役とユーザに委ねる。

## 6. 確認して正しかったもの（抜粋）

- DSC1 の p.20 の単位の写像と D53 の数値（K_p −2.41×10⁻⁴、K_i −703 s⁻¹、K_d −2.216×10⁻⁸ s）。`hardware/mapping.py` の `dsc1_to_canonical`・`canonical_to_dsc1`（符号の連動の検査、範囲外の扱い）。
- D61-2 の遅延の勘定（−180° が 60 kHz 付近、文字通りの読み方では 50 kHz で仕様違反）。f_sim = 2 MHz でも 0.5 + 0.75 + 0.25 µs = 1.5 µs が保たれる（`rate_div` 4 の保持 1.5 サンプル＋ 914 の半サンプル）。
- D61-1（駆動値を密度にする）: `DrvTarget.scale = n_sub/(hi − lo)`、refD はスカラーの配置（`test_refD_actuator_placement.py`）。
- VT-13: 多レートの z 領域の式（間引きの像 f + 250 kHz を含む）と時間ループが 10⁻³・0.1° で一致。
- VT-21: ランプでのロック外れの時刻が式と 1 % で一致、外れた後にフリンジを滑る（`n_slip` は巻き戻した真の位相で数える）。
- VT-24 test 9: レーザー × ΔL の全分散（帯域内＋ledger）= 2πΔντ。線形補間の端数遅延での帯域内の分散は 4π²τ²h0 f_N に一致（τ ≪ 1/f_sim）。
- T16 §5 の白色の量子化雑音の見積もり 1.5 mrad は、ΔL = 9 m の雑音下での増分 1.6 mrad と一致。
- 914 の 2 極の当てはめ、915B の 6.491 rad/V、PNH1550R5A2 と Opneti の値（M1D-06 を除く）。
- VT-38（偏光）19 件が通過。PM の U = R(θ_out)diag(1, e^{−iΔψ_b})R(θ_in) と per_budget の x = sin²θ、PER = −10 log tan²θ、曲げ複屈折の係数 ¼n³(p11 − p12)(1 + ν) ≈ 0.136。

## 7. スクリプト

| ファイル | 内容 |
|---|---|
| `reviews/scripts/m1_refD_dac_limit_cycle.py` | 独立の準連続モデルでのリミットサイクル（§2.1）と記述関数（§2.2） |
| `reviews/scripts/m1_refD_quant_in_noise.py` | 名目の雑音下での DAC 量子化の増分（§2.3） |
| `reviews/scripts/m1_refD_vt14_bands.py` | VT-14 の帯域ごとの比、窓の核での期待値、seed のばらつき（M1D-04） |
| `reviews/scripts/m1_refD_fsim_doubling.py` | f_sim × 2 の変化と seed のばらつき、決定的なリミットサイクルの rms（M1D-11） |
| `reviews/scripts/m1_refD_delay_accounting.py` | 914 の ZOH の半サンプルと −180° の周波数（M1D-05、D61-2） |
| `reviews/scripts/m1_refD_stability_grid.py` | D53 のゲインでの共振の掃引格子の安定判定（M1D-02） |
| `reviews/scripts/m1_refD_lin_names.py` | coil_birefringence・per_budget の YAML 読み込み（M1D-03） |
