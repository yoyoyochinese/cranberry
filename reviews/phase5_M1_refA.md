# Phase 5 物理レビュー: M1 の Ref A の系（WP-B1・B3（Ref A 部分）・B5・F3/F4 の測定部分・B2 のレーザーとフィルタキャビティ・I1 の `refA` / `refA_tap`）

- レビュー担当: 新しい文脈の物理レビュー担当（実装には関わっていない）。`blind/` は読んでいない。ponytail は使っていない。
- 基準: `dev`（worktree で `git merge --ff-only dev` した時点）。ソースコードは変更していない。追加したのは `reviews/` の下だけ。
- 読んだもの: `AGENTS.md`、`DECISIONS.md`（D34・D46・D48・D51・D55-4/5/6・D56-8・D58-4・D63 ほか）、`notes/conventions.md`、`PLAN.md` §7.4、T01・T04・T05・T12・T13・T14、`notes/00_summary.md` §4–§5（C23・C24）、Ref A の PDF（本文と Fig. 3–5 をページ画像で確認）、下記のコードとプリセットと VT。
- VT-01・02・04・05・15・16・30・33・42・43 は、この worktree で全件通過した（112 件、98 s）。
- 数値の確認に使ったスクリプトは `reviews/scripts/` に置いた（§8）。いずれも末尾に assert の自己検査がある。

---

## 1. 指摘の一覧

| ID | 重大度 | 場所 | 要旨 |
|---|---|---|---|
| M1-A-01 | **major** | `src/cranberry/blocks/noise/laser.py:123-142, 316-318`（`oob_pairs` の `dnu_fc_hz`）、`src/cranberry/offline/band_ledger.py:12-15`、`presets/reference/refA.yaml:31-32, 200-201` | レーザー × フィルタキャビティの帯域外の分散を台帳で二重に数えている。時間ループの白色 FM はブラウン運動の厳密な標本なので、FC の離散低域フィルタはエイリアスを含めた全分散 πh0/Δν_FC をすでに再現する。そこに台帳が ≈ h0/f_N（20 MHz で FC 項の +10.5 %）を足す。逆に、サブサンプルの遅延の線形補間は遅延差の項を約 10 % 小さくする。この二つは符号が逆で、20 MHz ではたまたま打ち消し合う。そのため D63-1 の「f_sim × 2 で 0.006 dB」は収束の証明になっていない |
| M1-A-02 | minor | `presets/reference/refA.yaml:17-25, 199-201, 257-259`、`presets/axes/refA_curve/previous_work.yaml` | Δν_FC = 3.3 MHz の合わせ込みについて（§6.3）。θ² ≈ πh0/Δν_FC なので、Δν_FC は h0（UNVERIFIED の代用値）と比でしか決まらない。効くのは S_ν の 0.1–3 MHz の帯域で、Lorentz 線幅ではない。この帯域にある緩和振動のピーク（`peaks: false`）が外れている。`dnu_fc_hz` と `in_delay_s` は手で二重に書いていて、整合を検査していない |
| M1-A-03 | minor | `src/cranberry/blocks/err/qnl.py:67, 71`、`tests/vt/test_vt_43_error_signals.py:410-434` | QNL の LO の RAM が暗雑音（回路雑音）にもかかっている。LO の RAM が変調するのは LO に比例する部分だけ。RAM によるずれの正しい目標は 2.54 mrad で、3.06 mrad ではない（D56-8 の判定、§4） |
| M1-A-04 | minor | `src/cranberry/blocks/err/_common.py:22-26`、`notes/error-signal-generation.md:16, 316` | 周期定常の重み 1.23 は 2 乗検出だけの値。ログ検出（Ref C）では雑音の重みは 1 で、傾きが 0.905 倍になる（PSD で 1.222 倍）。T14 §6.3 の「κ_det = 1.32 で 4.2 mrad」は rms に κ を掛けた算術の誤り。定常の値は 3.72 mrad、厳密な値は 4.11 mrad（§4） |
| M1-A-05 | minor | `presets/reference/refA.yaml:169-170, 186-187`、`src/cranberry/blocks/meas/esa.py:9, 114` | `rbw_shape: gaussian` と `enbw_factor: 1.0` を組み合わせると、−3 dB 幅が 0.939 MHz になる。分析器の RBW は −3 dB 幅なので、理想的な Gaussian では enbw_factor = 1.0645 が要る。B_stat は 3 %、ばらつきは 1.5 % 変わる |
| M1-A-06 | minor | `tests/vt/test_vt_33_esa.py:7-13, 52`、`notes/phase-noise-to-squeezing.md:155, 289`、`notes/00_summary.md:491`、D58-4 | log 平均のばらつき「0.083–0.089 dB」は矩形の RBW だけの値。Gaussian の RBW（refA が使う形）では 0.0706–0.0750 dB で、VT では検査していない。oversample 4 の Gaussian の log MC は分散が +1.3 % 偏る（§5） |
| M1-A-07 | minor（記録） | `tests/vt/test_vt_04_refA_fig5.py:11-16`、D63-7 | 1200 mW で許容幅を 0.20 dB に広げた件。θ(P) が増えるのは、ロック 2 の利得が A+(P) に比例する（電子回路を固定した）という [AE] の仮定から出る。独立のモデルで大きさは確かめた（+1.4 % / +3.6 %、§6.1）。ただし利得を P ごとに合わせ直す変種（UGF 一定）も並べないと、物理の発見なのか仮定の帰結なのかを区別できない |

blocker はない。M1-A-01 の Ref A への影響は小さい（真の値 θ = 9.04 mrad に対してシミュレーションは 8.97 mrad、R′− で 0.012 dB）。ただし、帯域外の分散を台帳に積む規則そのものの欠陥なので、M2 以降（Ref C・将来の構成で f_sim と Δν_FC・遅延の比が違う場合）に大きく効きうる。そのため major とした。

---

## 2. M1-A-01（major）: レーザー × FC の帯域外の分散の二重計上と、遅延の補間による不足

### 2.1 何が起きているか

- `laser_source` の白色 FM は 1 ステップごとに std √(h0 f_b/2) の周波数を引き、その累積和を位相にする（`laser.py:13-14, 187, 208`）。Var[φ(t) − φ(t − kh)] = 2π²h0·kh は連続のランダムウォークと一致する。つまり φ_n は**連続のブラウン運動の厳密な標本**で、その離散時間の PSD は f_N より上の成分を折り返して全部含んでいる。
- `filter_cavity` は φ_L(t − τ_b) に y_n = a y_{n−1} + (1 − a)(x_n + x_{n−1})/2（a = e^{−πΔν h}）をかけ、ψ = y − x を書く（`filter_cavity.py:106-118`）。ランダムウォークの増分 w（分散 2π²h0 h）からの伝達関数は G(z) = −(1 + a)/2 / (1 − a z⁻¹) なので、次が厳密に成り立つ。

  Var(y − x) = 2π²h0 h (1 + a)/(4(1 − a)) = πh0/Δν_FC · [1 + ε²/12 + …]、ε = πΔν_FC h。

  これは**連続の全帯域の値**（f_N より上を含む）に O(h²) で一致する。
- それにもかかわらず、`oob_pairs: {target: FC, dnu_fc_hz}` は台帳に psd_outside = |1 − H_FC|² S_ν を宣言する。`phase_to_quadrature` は ∫_{f_N}^{f_max} |1 − H_FC|² h0/f² df ≈ h0/f_N を足す。これは帯域内ですでに数えた分を**もう一度**数えている。
- 一方、LO と squeezed 光の遅延差 δτ = 16.8 ns（f_sim = 20 MHz で 0.34 サンプル）は、beam_assembler と FC の入力の線形補間で作られる。ブラウン運動を線形補間すると、分散がブリッジの分 2π²h0 h μ(1 − μ) だけ減る（D57-7 の白色雑音の話と同じ性質）。この不足は台帳で補われていない。

### 2.2 数値（`reviews/scripts/m1_fc_ledger_double_count.py`、refA で白色 FM だけをオン、mech・センサ雑音・量子化はオフ、0.05 s）

連続の真の値（FC と遅延の両方、全帯域）は ∫|H_FC e^{−iωδτ} − 1|² h0/f² df = 8.173×10⁻⁵ rad²（9.041 mrad）。

| f_sim | ループ内（シミュレーション） | 台帳の σ²_HF | 合計 / 真の値 | σ_total | R′− |
|---|---|---|---|---|---|
| 10 MHz | 7.58×10⁻⁵ | 1.23×10⁻⁵ | 1.078 | 9.385 mrad | −10.112 dB |
| 20 MHz | 7.42×10⁻⁵ | 6.31×10⁻⁶ | 0.985 | 8.973 mrad | −10.150 dB |
| 40 MHz | 7.62×10⁻⁵ | 3.18×10⁻⁶ | 0.971 | 8.907 mrad | −10.157 dB |

離散の伝達関数から解析的に予想したループ内の値（ロックなし）は次のとおりで、シミュレーションと 1.5 % 以内で一致した。

| f_sim | 10 MHz | 20 MHz | 40 MHz | 80 MHz | 160 MHz |
|---|---|---|---|---|---|
| ループ内 [10⁻⁵ rad²] | 7.47 | 7.31 | 7.49 | 7.84 | 7.99 |
| ＋台帳 h0/f_N [10⁻⁵ rad²] | 8.74 | 7.94 | 7.81 | 8.00 | 8.07 |

- 合計は f_sim について単調でない（遅延の端数 μ が f_sim とともに変わるため）。20 MHz と 40 MHz の差が小さい（0.006 dB）のは偶然で、収束先（8.17×10⁻⁵）とは 20 MHz で −2.8 %（R′− で約 0.012 dB）違う。D63-1 の「f_sim × 2 の確認に合格」は、二つの O(h) の誤差（台帳で +h0/f_N、補間で −2π²h0 h μ(1 − μ)）が打ち消し合った結果である。
- Ref D 型の **サブサンプルの遅延だけ** の組（`tau_s`、τ < h）は問題ない。ループ内は μ·2π²h0τ、台帳は (1 − μ)·2π²h0τ（2τf_N = μ）で、和は 1 次まで厳密になる。欠陥が出るのは、IIR のフィルタ（FC）がサンプルの上で働く場合と、遅延が 1 サンプルを超える場合（k ≥ 1 では台帳が 0.2026·2π²h0h を足すのに対し、不足は μ − μ²）である。

### 2.3 修正案

1. `oob_pairs` の意味を「f_N より上の PSD」から「連続の真の値 − 離散のモデルの期待値」に変える。帯域外の分散は、ブロックの離散伝達関数（補間・Tustin/ZOH を含む）から offline で厳密に計算できる。`variance_outside` で宣言すればよい（2.2 の解析式と同じ計算で、スクリプトの `interp.py` 相当の 20 行程度）。
2. 少なくとも `dnu_fc_hz` の組は帯域外を 0 にする（FC 部分はループ内で完結している）。LO と squeezed 光の遅延差は `tau_s` の組で補う。ただし FC との交差項があるので、1. の方法が確実。
3. 収束の確認（PLAN §6.4）は f_sim を ×2 と ×4 で行い、端数 μ を固定した比較（または解析の期待値との比較）にする。D63-1 の記述を差し替える。
4. VT に「白色 FM ＋ FC ＋ 遅延の系で、σ² が連続の解析値に 1 % で一致する」を追加する（VT-34 の FC の項目の拡張）。

---

## 3. その他の minor の詳細

### M1-A-02: Δν_FC = 3.3 MHz の合わせ込み（§6.3 も参照）

- θ² の主項は ∫|1 − H_FC|² S_ν/f² df で、白色 FM なら π h0/Δν_FC = Δν_L/Δν_FC になる。**Δν_FC と h0 は比でしか効かない**。Koheras の仕様（Lorentz 線幅 < 0.1 kHz 程度）を使うと、遅延項を含めて Δν_FC ≈ 1.4 MHz で 9 mrad になる。
- 重み |1 − H_FC|²/f² は f ≲ Δν_FC/2 で平らなので、効くのは S_ν の 0.1–3 MHz の帯域である。ここには EDFA・ファイバレーザーの緩和振動のピーク（T05 E16、RIN の 0.7 MHz）がある。白色 FM の Lorentz 線幅で代表させるのは弱い。
- 推奨: Δν_FC を「h0 との組（または S_ν(f ≲ Δν_FC/2)/Δν_FC）」の 1 軸として記録する。`peaks: true`（緩和振動）の変種で掃引する。`previous_work` の 14 mrad を Δν_FC = 1.175 MHz で作る仮定は、[23] の構成が不明なので「仮定」と明記する（今もそう書いてはあるが、C23 の説明に使わないこと）。
- `oob_pairs.dnu_fc_hz` と `fc.dnu_fc_hz`（refA.yaml:201, 257）、`fc.in_delay_s` と LO の区間の群遅延（:259）は手書きの重複で、検査がない。組み立て時に一致を検査するか、片方から導出する。

### M1-A-03・04: D56-8 の付随（§4 に詳細）

### M1-A-05: Gaussian の RBW と ENBW

- `esa.py` の Gaussian は |H|² = exp(−πf²/ENBW²)（ENBW = enbw_factor·RBW）。−3 dB の全幅は 0.9394·ENBW。分析器の RBW の公称値は −3 dB 幅なので、理想的な Gaussian では enbw_factor = 1.0645 になる（Keysight X シリーズの表記の約 1.06 と同程度、UNVERIFIED）。refA の 1.0 は B_stat を 3 % 大きく、ばらつきを 1.5 % 小さく見積もる。VT の許容幅に比べれば小さい。

### M1-A-06: VT-33 の目標

- §5 の表のとおり、RBW の形ごとに目標を書き分ける。Gaussian の log の値（0.0706–0.0750 dB）を VT に追加する。Gaussian の log の MC は oversample ≥ 8 にするか、分散で +1.3 % を許す（Li₂ の τ² ln τ の尖りを粗い標本化で拾うため）。

### M1-A-07: 1200 mW での θ(P)

- 独立の線形モデル（`m1_lock2_bode.py`: 積分器、遅延 1 µs、復調 LPF 250 kHz、UGF ∝ A+(P)）で、640 → 1200 mW の θ の増加は +1.4 %（9 mrad）、+3.6 %（14 mrad の変種）だった。D63-7 / VT-04 の +1 % / +3.4 % と合う。物理（遅延のあるループの感度の積分 ∫(|S|² − 1)df > 0）は正しい。
- ただしこれは「電子回路の利得を P によらず固定した」という [AE] の仮定から出る。実験で P ごとに利得を合わせ直していれば、Fig. 5 の θ 一定は自然である。UGF 一定の変種を VT-04 に並べて、両方を出すことを勧める。

---

## 4. D56-8 の判定: **確認（物理）。数値は修正**

**主張**: dither ロックのセンサ雑音は dither の周期で変わる（周期定常）ので、復調後の PSD は 2⟨2 sin²θ S_v(θ)⟩ になる。Ref C の残差は 3.24 → 3.62 mrad、RAM のずれは 2.77 → 3.06 mrad に変わる。

**独立の導出**（`reviews/scripts/d56_8_cyclostationary.py`）

- θ は dither の位相 2πf_m t（直交位相の角ではない）。radiometer の出力雑音の相関時間は約 1/B ≈ 3 ns で 1/f_m よりずっと短い。そのため局所的に白色とみなせる: n(t) = √S_v(Δ0 + m sin θ(t))·w(t)。
- 2 sin θ で復調し（誤差は sin 成分にある。コードでは φ_d = −π/2 で I = b）、f_m より下に低域通過させると、片側 PSD は S_e = ⟨4 sin²θ S_v⟩ = 2⟨2 sin²θ S_v⟩ になる。S_v が一定なら 2S_v で、T14 §2.2 の「復調で 2S」と一致する（片側 PSD・2 倍の因子は整合）。
- 重みが 1 を超えるのは、sin²θ が dither の端で最大になり、squeezing 側にロックしているとそこで S_v(Δ) が最大になるからである。

**Monte Carlo**（2 乗検出の radiometer を直接模擬、4000 ブロック、3σ ≈ 0.067）

| 条件 | 測定 / 定常の式 | 測定 / 周期定常の式 |
|---|---|---|
| Ref C、m = 0.0785、squeezing 側 | 1.251 | 1.021 |
| m = 0.3（誇張） | 6.79 | 1.029 |
| 反 squeezing 側 | 0.993 | 1.002 |

**修正した目標**

| 量 | D56-8 / T14 | 修正値 |
|---|---|---|
| Ref C のセンサ雑音による残差、2 乗検出、連続のループ（b_loop 10 kHz） | 3.24 → 3.62 mrad | **3.24 → 3.58 mrad**（重み 1.2257）。3.62 は VT-43 の離散ループ（f_sim 1 MHz）の因子 1/(1 − ε/2) = 1.0204 を含んだ値。同じ条件どうしでは 3.27 → 3.62 |
| 同、ログ検出（κ_det = 1.32、Ref C の構成） | T14: 4.2 mrad | **4.11 mrad**。ログ検出では雑音の重みは 1 で、傾きが 0.9047 倍になる（PSD で 1.222 倍）。T14 の 4.2 は rms に κ を掛けた誤り（定常の値は 3.24 × √1.32 = 3.72 mrad）（M1-A-04） |
| 反 squeezing 側のロック | 78 mrad | 変わらない（重み約 0.99） |
| RAM によるずれ（a = 10⁻²、φ_a = 0） | 2.77 → 3.06 mrad | **2.54 mrad**（2 乗検出。ログ検出は 2.53）。3.06 はコードのモデル（暗雑音にも RAM がかかる、`qnl.py:67, 71`）でなら正しい。しかし LO の RAM は回路雑音を変調しない（T14 §6.4 自身も「回路雑音は除く」と書いている）（M1-A-03） |
| 同、a = 10⁻³ | — | 0.254 mrad |

- シミュレータの `harmonic()`（`_common.py:107-130`）は、検出の法則によらず周期内の求積で正しく扱っている。直すのは、RAM が暗雑音にかかる点（光の部分 p_mean − S_el と p_delta だけに RAM をかける）と、docstring・ノートの数値である。
- 更新先: T14 §6.3・§6.4・V12、00_summary の VT-20（「QNL のセンサ雑音だけで 3.2 mrad」→ 2 乗検出 3.58 / ログ検出 4.11 mrad）、`test_vt_43_error_signals.py:420-434` の 3.06 → 2.54。
- 適用条件: video 帯域 ≫ f_m、f_m ≪ B、ループ帯域 ≪ f_m。

---

## 5. D58-4 の判定: **確認（物理）。数値は RBW の形ごとに修正**

**主張**: ESA の log 平均のトレースのばらつきは π²/6 ではなく ∫Li₂(|ρ|²)/∫|ρ|² ≈ 1.325（矩形）/ 1.341（Gaussian）で決まり、0.083–0.089 dB になる。Gaussian の RBW の統計的な帯域は √2·ENBW。

**独立の導出**（`reviews/scripts/d58_4_log_scatter.py`）

- X = |z₁|²、Y = |z₂|²（相関 ρ の複素ガウス）。結合密度の Laguerre 展開 e^{−x−y} Σ|ρ|^{2k} L_k(x) L_k(y) と ∫ln x·e^{−x}L_k(x)dx = −1/k から、Cov(ln X, ln Y) = Σ|ρ|^{2k}/k² = Li₂(|ρ|²) になる。ρ = 1 では π²/6 で、これは**独立な標本**の分散でしかない。
- VBW ≪ RBW では、平均の分散は log で (1/T)∫Li₂(|ρ|²)dτ、power で (1/T)∫|ρ|²dτ。B_stat = (∫|H|²)²/∫|H|⁴。1 極の video フィルタの両側の雑音帯域は πVBW なので、Var(P̄)/P̄² = πVBW/B_stat。
- Gaussian では閉じた形が得られる: 比 = ζ(5/2) = 1.3415、B_stat = √2·ENBW。

**数値**: Li₂ の恒等式は Monte Carlo と 10⁻³ で一致した（|ρ| = 0.3/0.6/0.9/1.0）。時間領域の MC（RBW → |z|² → 1 極 / boxcar）では、power の分散が πVBW/B_stat に 0.4 % で一致した（Gaussian の √2 を含む）。log/power の比は、矩形 1.330 ± 0.009、Gaussian 1.360 ± 0.004 だった。Gaussian の超過分は標本化の粗さ（fs = 4·ENBW）で説明でき、fs = 16 では 1.3418 になる。

**修正した目標**（RBW 1 MHz、VBW 100 Hz。範囲は ENBW 係数 1.128–1.0）

| RBW の形 | power 平均 [dB] | log 平均 [dB] | 参考: log で π²/6 [dB] |
|---|---|---|---|
| 矩形 | 0.0725–0.0770 | **0.0834–0.0886** | 0.0930–0.0987（T04 の値） |
| Gaussian | **0.0609–0.0647** | **0.0706–0.0750** | — |
| Gaussian（−3 dB 幅 = RBW、ENBW = 1.0645·RBW） | 0.0627 | 0.0727 | — |
| 4 極の同調形 | 0.0593–0.0630 | 0.0683–0.0726 | — |

- D58-4 の「0.083–0.089 dB」は**矩形の値**なので、Gaussian の文の横に置かないこと。refA は Gaussian を使う。
- Ref A の報告値 0.07 dB は、Gaussian の log、矩形の power、4 極のどれとも合う。したがって、平均の種類は判別できない（T04 §5.2 の UNVERIFIED のまま）。
- AN 1303 の「9.3 dB·√(VBW/RBW)」は、矩形で B_stat = ENBW = 1.128·RBW、π²/6 を仮定した場合に当たる。実際の 4 極のフィルタには両方の点で当てはまらない。
- コードの整合: `esa.py:60` の LOG_VAR_RATIO {1.3248108, 1.3414873}、`:62` の _BSTAT {1, √2}、`:19-23` の πVBW、`radiometer_std_db`、掃引のガンマ分布（形状 k = B_stat/(πVBW)）はすべて正しい。WP-B3 の κ_det（ログ）1.32 も 1.325 と合う。
- 更新先: T04 §5.3 の表（:155）と §12 V7（:289）、00_summary の VT-33（:491）、`test_vt_33_esa.py` の docstring。Gaussian の log の目標を VT に追加する（M1-A-06）。

---

## 6. C23・C24・Δν_FC の評価

### 6.1 C23（Ref A の lock 1 の SNR と 14 → 9 mrad）

**結論: シミュレータの説明（lock 1 の SNR では 14 → 9 mrad を説明できない）は物理的に正しい。コードがこの結論を偏らせている点は見つからなかった。ただし「9 mrad はフィルタキャビティを通ったレーザー雑音で決まる」は、合わせ込みで入れた仮説であって、発見ではない。**

- `bhd_probe_beat` の式を、場から独立に導出し直した。probe（+Ω）は e^{i(u − Ωt)}、idler（−Ω）は S†âS の â† e^{iφ_p} の項から e^{i(Ωt − u)}·g_s。LO とのビートは g_c cos(Ωt − (u − v)) + g_s cos(Ωt − (u + v)) で、コード（`coherent.py:150-157`）と一致する。φ_d を合わせると I ≈ K2[A−·δu − A+·Δ] になり、零点は Δ = (A−/A+)·δu。lock 1 の誤差 δu は 1/123 に縮んで Δ に入る。これは構造（LO を squeezer を通った probe にロックする）から出る性質で、プリセットの数値には依らない。
- 必要な大きさ: 14 → 9 mrad の差を lock 1 で説明するには、δu_rms ≈ 123 × √(14² − 9²) mrad ≈ **1.3 rad** が要る。ロックが保たれている状態ではありえない。
- フレームの扱い（`bhdbeat.phi_opa_rad = −π` と `dem2.phi_d_rad = −π/2`）も確かめた。u に一定のずれを入れることは u − v と u + v の両方に同じ位相を足すのと同じなので、復調位相で厳密に打ち消される。誤差の式は物理のフレームの式と同一になる。
- lock 2 のセンサ雑音（BHD のショット雑音、noise_rel_shot 0.1、K2·A+ = 2ℜ√(P_LO·1.6 nW)·G）は約 0.5 mrad で、ヘッダの値と一致した。1.6 nW は Fig. 4(a) の **反 squeezing の** トレースの 1 MHz の線（ショット雑音より約 41 dB 上）から出る値で、P_pr A+² の読み方と整合する。squeezing のトレースの左端（約 −17 dBm）は、回路雑音のトレースにも出ている ESA の 0 Hz の漏れなので、使ってはいけない（T14 はこれを使っていない）。
- モデルにない（＝結論が依存する）仮定:
  - (a) [23] の lock 2 の構成が同じであること。
  - (b) EOM1 の RAM・AOM の駆動位相など、lock 1 のアクチュエータから LO 経路への結合がないこと。
  - (c) 14 → 9 mrad は、Ref A が同時に変えた他の要素（光学部品・フィルタキャビティ・損失）で起きたと考えられること。

  Ref A §4.1 の最後の文（「seed レーザーの雑音を高フィネスのキャビティで減らすべき」）は、レーザー雑音が支配的という仮説と整合する。しかし θ の内訳は論文にない。
- 推奨: C23 は「lock 1 の SNR では説明できない（構造的に A−/A+ で抑えられる）」として閉じる。「9 mrad はレーザー × FC で決まる」は仮説（M1-A-02 の縮退つき）として記録する。

### 6.2 C24（Fig. 3 の極小の深さと VBW）

**結論: T04 の見積もりは正しく、シミュレータも偏っていない。食い違いは論文の Fig. 3 の掃引トレースの取得条件の側にあると判断する（解決はできないが、物理の不具合ではない）。**

- 独立の計算（`reviews/scripts/m1_c24_scan_minimum.py`、ω = 71 rad/s、1 極の video フィルタ）による極小の深さ:

  | VBW | 100 Hz | 300 Hz | 1 kHz | 3 kHz | 10 kHz |
  |---|---|---|---|---|---|
  | power | +1.72 dB | −6.06 dB | −9.92 dB | −10.50 dB | −10.58 dB |
  | log | −3.23 dB | −7.93 dB | −10.07 dB | −10.51 dB | −10.58 dB |

  T04 の +1.76 / −3.19 dB と合う。VBW 100 Hz で −9.5 dB に届くには、ω ≤ 9.3 rad/s（π の周期 ≥ 0.34 s）が要る。図は 0.044 s である。
- Fig. 3 を 110 dpi で描画して確認した。掃引トレースの極小は**左右対称の尖り**で、ロック時の −10 dB と最大の 20 dB に正確に届いている。時定数が尖りの幅（約 2.5 ms）と同程度の 1 極のフィルタなら、極小は浅くなるうえに時間方向に非対称（遅れて裾を引く）になるはずである。このことから、掃引トレースには VBW 100 Hz の 1 極のフィルタがかかっていない（掃引だけ別の VBW、または後処理）と考えるのが自然。
- シミュレータ（`esa.py:142-155`）は video フィルタを検出の領域（power / log / voltage）で 1 極としてかけていて、T04 §5.1 と一致する。ロック時と掃引で同じ処理なので、結論を偏らせない。VT-33 V9 を「単調性だけ」にしているのは妥当。
- 推奨: C24 は「Fig. 3 の掃引トレースは VBW 100 Hz の 1 極とは両立しない（VBW ≳ 1 kHz 相当）。合わせ込まない」として閉じる。

### 6.3 Δν_FC = 3.3 MHz は妥当な代用値か

- **代用値としては擁護できる。** ASE 除去用のフィルタキャビティとして MHz 級の線幅は自然であり、Ref A 自身がレーザー雑音を限界に挙げている。モデルの誤りを隠すために合わせ込んだものでもない。M1-A-01 の誤差を全部直しても、同じ Δν_FC で θ は 8.97 → 約 9.04 mrad（白色 FM だけの場合）になるだけで、合わせ込みの値はほとんど動かない。
- ただし、次の三つは記録として残すべき。
  1. h0 との縮退（実質的に合わせたのは h0/Δν_FC の比）。
  2. 効くのは S_ν(0.1–3 MHz) なのに、緩和振動のピークを外していること。
  3. 「previous work」の 14 mrad を Δν_FC だけで作っていること。

  この三つのため、Δν_FC の値から物理的な結論（例えば「線幅 X kHz のレーザーが必要」）を引き出してはいけない。将来の構成のレーザーへの要求は、S_ν(f) の曲線として出すこと（D49 の「要求の曲線」）。

---

## 7. 確認して正しかったもの（PLAN §7.4 の点検表）

(1) **式がノートと一致**
- `single_mode_reference`: R± = L + (1 − L)e^{±2√(αP)}、厳密（χ = e^{−2θ²}）と小角。独立に計算すると R′− = −10.1479 dB、R′+ = 20.5561 dB で、小角との差は 10⁻⁴ dB 未満。P_opt = 612.3 mW、−10.152 dB。VT-16(b) の T_tap = 0.1 / 1 / 3 % は −10.108 / −9.760 / −9.075 dB。いずれも目標と一致した。
- `opa_waveguide`: 結合波の式 → 直交位相の行列 G = [[Re c, Im c], [Im c, −Re c]] + δ·J − (a_s/2)𝟙、拡散 a_s V_vac、κ = √α/ℓ（r = √(αP)）、分布損失の dR±/dz = ±2gR± − a_s(R± − 1)。A_p(0) = √P e^{i(φ_OPA,0 − π/2)} によって θ_amp = φ_OPA,0/2 になる。
- `parametric_beat_detector`: P = ηP_in[Ḡ + ΔG cos(2Ωt − Φ1)]、Ḡ = cosh 2r、ΔG = sinh 2r。場から再導出して一致した。PD のショット雑音は DC の Ḡ に比例。
- `bhd_probe_beat`: §6.1 のとおり。
- `filter_cavity`: H_FC = 1/(1 + 2if/Δν)、極 πΔν [rad/s] ＝ 半値半幅 Δν/2 [Hz]。`in_delay_s` を使うと LO = H_FC φ_L(t − τ_b) が厳密に成り立つ。
- `pulsed_bhd_kernel`: c_m = e^{−2πBm/f_rep}（1 次、任意の W で厳密）、ρ_adj = c_1。
- `detector_efficiency`: η_e = C/(1 + C)。28 dB で損失 0.158 %。C(f) = |H|²C_0 で DC に換算（refA の frontend の (1 + (3/35)²)）。
- 損失の総和は 0.98 · 0.99 · 0.99 · η_e · η_mm = 0.92（Ref A §4.2 の内訳。η_mm は差し引きで 4.07 %）。

(2) **符号**
- θ_sq = φ_p/2 + π/2（VT-30）。SHG で φ_p = 2φ_L、`sqz` の記録用ビームは φ_L。
- AOM の +Ω（A ∝ e^{i(φ − Ωt)}）、ビートの cos(Ωt − (u ∓ v))、復調の I = a cos φ − b sin φ（2cos の参照）。
- lock 1 の零点 Φ1 = 0（増幅側）。
- φ_L = −2π∫δν dt。

(3) **2 倍の因子**
- 片側 PSD で、白色 FM の 1 ステップの分散 h0 f_N。
- 復調で 2S（基底帯域で a・b 各 2S_v、RF で S_v、f_b/2 の換算）。
- R± の定義（ショット雑音単位、V_SN = 1/2）。
- radiometer の πVBW（両側）。

(4) **単位**
- Hz と rad/s: FC の a = e^{−πΔν h}、LPF の α = 1 − e^{−2πf h}、`ki` の 2π·f_u/K。
- dB と線形: `clearance_lin`、`_lin` への換算。
- ℜ = 0.99eλ/(hc) = 1.2338 A/W。V_π の読み方: lock 1 は 2π/V_π（Φ1 = 2φ_pr）、lock 2 は π/V_π。

(5) **FFT**: `esa._z2` は fft_agrawal / ifft_agrawal と omega_grid を使い、|H|² は偶関数なので符号に依らない。規格化は E|z|² = 1。

(6) **数値の収束**: M1-A-01 を除き問題なし。opa_waveguide は DOP853（rtol 10⁻¹²）。VT-02 の統計誤差は約 0.001 dB。

(7) **トグルでオフにしたときの恒等**: tap（`enabled` が off なら要素は無効で記録され、P_tap = 0）、opa の pump_rin が off なら bands は空、FC の各トグル、detector_efficiency の各トグル。

(8) **損失・帯域の登録の漏れと二重**: 損失は lumped の導波路を opa の `losses()` だけで数え（single_mode_reference は L を R± に含めるので登録しない）、FC の T0 は LO 経路だけ、検出器はオフラインの 1 か所、frontend は `exclude_losses` で電子・不平衡を外す。二重は見つからなかった。帯域の欠陥は M1-A-01 だけ。センサ雑音の帯域はスロットが対象なので、指標には入らない。

(9) **値の出典の印**: refA / refA_tap のすべての値にタグがある。1.6 nW の読み方は §6.1 で確認した。

(10) **保護リスト**: Ref A の系のブロック（single_mode_reference、opa_waveguide、tap、parametric_beat_detector、bhd_probe_beat、demodulator、detector_efficiency、bhd_frontend、ESA の 2 つ、mode_overlap、pulsed_bhd_kernel、laser_source、filter_cavity、loss_budget、metric_i、quantum_noise_power）は、すべて登録されている（`registry.registered_kinds()` で確認）。トグルも残っている。

---

## 8. スクリプト

| スクリプト | 内容 | 実行 |
|---|---|---|
| `reviews/scripts/m1_fc_ledger_double_count.py` | M1-A-01: FC の離散フィルタの閉じた式、refA（白色 FM だけ）を 10/20/40 MHz で実行 | `uv run python reviews/scripts/m1_fc_ledger_double_count.py --sim`（約 3 分） |
| `reviews/scripts/m1_lock2_bode.py` | C23 / D63-7 / M1-A-07: lock 2 の線形モデルで θ(P) | `uv run python reviews/scripts/m1_lock2_bode.py` |
| `reviews/scripts/m1_c24_scan_minimum.py` | C24: 掃引トレースの極小の深さと VBW | `uv run python reviews/scripts/m1_c24_scan_minimum.py` |
| `reviews/scripts/d56_8_cyclostationary.py` | D56-8: 周期定常の radiometer の MC、Ref C の残差、RAM のずれ | `uv run python reviews/scripts/d56_8_cyclostationary.py` |
| `reviews/scripts/d58_4_log_scatter.py` | D58-4: Li₂ の恒等式、RBW の形ごとの比と B_stat、時間領域の MC、Ref A のばらつき | `uv run python reviews/scripts/d58_4_log_scatter.py` |

M1-A-01 の 2.2 の第 2 表（線形補間の遅延を含めた離散モデルの期待値、10–160 MHz）は、`m1_fc_ledger_double_count.py` の Part 1b が出力する（`--sim` なしで数秒）。
