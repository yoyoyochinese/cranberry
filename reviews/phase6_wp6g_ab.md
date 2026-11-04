# Phase 5 レビュー: WP6-G（D90、ゲート）と WP6-AB（D89、実験室の音響・915B の複数モード）

- 日付: 2026-11-04 ／ 対象: dev @ d0ec2cd（`git merge --ff-only dev` 後）
- 読んだもの: AGENTS.md、DECISIONS（D61・D68・D75・D80・D82・D85・D89・D90）、`notes/conventions.md`、`notes/lit-madsen2022.md`・`lit-nardelli2026.md`・`lit-luo2025.md`、
  `predictions/refD/comparison/COMPARISON.md`（読むだけ）、`reviews/phase5_refD_stability_mismatch.md`（形式）、
  `blocks/mod/gate_schedule.py`、`blocks/ctrl/gated_hold.py`、`blocks/chk/gated_loop_stability.py`、`blocks/chk/loop_stability.py`、`blocks/err/coherent.py`（`in_gate`）、
  `blocks/meas/metrics.py`・`runner/decompose.py`・`runner/run.py`（`mask_slot`）、`blocks/noise/_lgs.py`・`mech.py`（`table`）、`presets/axes/gate/*`・`lab/*`・`controller/dsc1.yaml`・`moku_delta.yaml`、
  `presets/future/F_common.yaml`、`presets/reference/refD.yaml`、`devices/dsc1.yaml`、関連のテスト。
- **`blind/` は読んでいない。`predictions/` は何も書き換えていない。プリセット・ソースは変更していない。** DECISIONS.md も変更していない。
- 高速のテスト: `uv run pytest -m "not slow"` → **1550 passed, 8 skipped**（matplotlib・botorch が無いための skip を含む）, 263 deselected（36 分）。
- 根拠のスクリプト（`reviews/scripts/`、すべて読むだけで、プリセットを書き換えない）:
  - `p6_wp6g_floquet_indep.py` … ゲートしたループの独立な再実装（scipy の bilinear の区間、自前の遅延線・PID・ゲート）で Floquet 乗数、周期 Lyapunov 方程式で周期定常の分散、純音で変換（Dick）の項。TEC の rate_div の保持、multimode と dsc1 の組み合わせ（4 秒）
  - `p6_wp6g_timeloop.py` … **実際の F_liu20 ＋ gate/borealis_window の時間ループ**と検査の比較（不安定点、安定点の測定窓の分散）（約 2 分）
  - `p6_wp6g_compose.py` … オーバーレイの組み合わせで、何が大きな声で失敗し、何が黙って物理を変えるか
  - `p6_wp6g_mask_decim.py` … `mask_slot` と記録器のデシメーション、窓に標本の無い格子点
  - `p6_ab_check.py`（既存、WP6-AB の作者のもの）を実行して D89-4 の値を再現
  - `p6_wp6ab_table_psd.py` … 表の PSD の腕ごとの分け方（時間ループ）と外挿の部分の分散
  - `p6_wp6ab_quant.py` … 25・40 kHz の山の DAC 量子化の寄与（量子化の on/off、白色の見積もり）
  - `p6_wp6ab_multimode_alt.py` … 5.3 kHz の副モードの非一意性（遅延・極・直列の biquad、r と Q の走査）
  - `p6_wp6ab_dsc1_h.py` … 項目 6: `controller/dsc1` の h = 1/f_sim と装置の 500 kHz

## 0. 結論

| 対象 | 判定 | 要点 |
|---|---|---|
| A-1 Floquet・周期定常の分散・Dick の項 | **正しい** | 独立な再実装と乗数が 10⁻⁹ まで一致（分数遅延、tustin/dsc1、zero/hold、漏れ、ゲートの遅れ、settle）。分散は周期 Lyapunov と 0.00 %、変換の項は純音の時間領域と 5 桁一致。実際の F_liu20 の時間ループと測定窓の分散が 0.2 %・1.7 %・−7 %（hold、統計の範囲）。TEC の rate_div の保持は遅い乗数を 10⁻⁵ 動かすだけ（無視してよい） |
| A-2 プローブのゲートの物理 | **正しい**（小さな注記のみ） | 縞 ∝ √g・ショット雑音 ∝ g は「強い参照と弱いプローブの DC/ヘテロダインのビート」で正しい。τ_pr,k は AOM → det_k のビームの経路と一致。マスクはカーネル内の積算では正しく同じ時刻 |
| A-2' マスクと記録の系列 | **major ×2** | デシメーションした記録ではマスクした系列が空・全部・混ざる（G-1）。出荷した掃引の格子点に窓の標本が 0 の点があり、run が例外で止まる（G-2） |
| A-3 上書きの組み合わせ | **major ×3** | `controller/dsc1` → gate で ADC が黙って fpga_generic（LSB が 8 倍細かい）、ゲートの検査は tustin のまま（G-3）。`lab/multimode915b_future` ＋ gate で検査が単一モードのまま（108 点中 6 点で楽観的な誤判定、G-4、D90 に記載なし）。ヘテロダインと周波数の上書きで後勝ち（G-5） |
| B-4 表の PSD | **正しい**（片側、腕の分け方も時間ループで確認）／**major 1** | 外挿（3 kHz より上の f⁻² を無限まで、γ = 1 で L/13.2 m 倍）が将来のプリセットの結果を支配する（B-1） |
| B-5 複数モード | 値は再現／**major ×2** | D89-4 の値は再現した。25・40 kHz の山がモデルの DAC 量子化なのは確認したが、決定論的なリミットサイクル（白色の見積もりの 2.6 倍）で装置では未検証、主モードの Q = 25 はそれに条件付き（B-2）。5.3 kHz の構造は \|S\|² からは一意に決まらない（B-3） |
| B-6 dsc1 の h | 項目 6 に回答（§6） | h を 1 → 2 µs にしても 27 点で判定は変わらない。大きいのは F_* ＋ DSC1 に入出力の極（120/100 kHz）が無いこと（B-7） |

blocker はない。新しい 3 ブロックの数値の核（モノドロミー、分散、Dick の項）は正しく、結果に使ってよい。ただし G-1〜G-5・B-7 を直すまで、**ゲートと制御器・複数モードの上書きを組み合わせた点、デシメーションした記録の系列（ゼロスパン、decompose の有効性の検査）、外挿に依存する acoustic_future（γ = 1）の数字は使わない**こと。

## 1. 指摘の一覧

| # | 重さ | 指摘 | 再現 |
|---|---|---|---|
| G-1 | major | `mask_slot` を記録した（デシメーションした）マスクでかけるので、系列が空・全部・窓をまたいで混ざる。decompose の有効性の検査も同じ（黙って通る） | `p6_wp6g_mask_decim.py` A・B |
| G-2 | major | gate/borealis_window の格子点 (T = 10 µs, D_lock = 0.9) は f_sim = 1 MHz で測定窓の標本が 0。`PhaseStats.from_acc` が ValueError、`runner.sweep` は点ごとに捕まえないので掃引全体が止まる（D75 の「全点を残して印を付ける」に反する） | 同 A |
| G-3 | major | `controller/dsc1` → `gate/*` は読み込めてしまい、ADC が dsc1（16 bit ±4 V）から fpga_generic（±0.5 V）に戻る。F_liu20 の名目（K = 6.3 mV/rad）で ADC の LSB が 19 → 2.4 mrad。`chk.gated_stability` は tustin のまま（時間ループは dsc1、27 点中 1 点で楽観的） | `p6_wp6g_compose.py`、`p6_wp6g_floquet_indep.py` (e) |
| G-4 | major | `lab/multimode915b_future` と gate の組み合わせ（順序によらず）で、`chk.gated_stability` は F_common の単一モード（主モード、r = 1）のまま。さらに `_dynamics_d` は r を無視して直列の積にするので、モードの表を写しても和にならない。108 点中 6 点で「安定」と誤判定（逆向きは 0） | `p6_wp6g_floquet_indep.py` (d)、`p6_wp6g_compose.py` |
| G-5 | major | `probe/detect_heterodyne` と `probe/freq_*` は `probe_freq_offset_hz` を取り合い、後に当てたほうが黙って勝つ（1 THz → 80 MHz）。D90-6 に記録済みだが、黙って違う物理になる | `p6_wp6g_compose.py` |
| G-6 | minor | 大きな声で失敗する組み合わせ（gate → controller/*、gate ↔ detect_*）のメッセージが `pipeline: unknown keys ['params'] in seq element` だけで、原因の上書きを指さない | `p6_wp6g_compose.py` |
| G-7 | minor | 周期の軸 `!sweep_log [1e-5, 1e-2, 7]` のうち 3 点（31.6 µs など）は標本数が整数でない。検査は N = round(T f_sim) の周期で解くが、時間ループは準周期。0.65 × 100 µs = 65 µs のように窓の端がちょうど標本に乗ると 1 標本ずれる | 読み（docstring にも注記あり） |
| G-8 | minor | `chk.lock_stability` は 0.5/f_sim の ZOH を足し、`chk.gated_stability`（とその `ungated_*`）は足さない。時間ループに合うのはゲートの検査の方。二つの判定が境界でずれうることを記録すべき | 読み |
| G-9 | minor | チョッパの端（2 µs の線形のランプ、√g で振幅が変わる）を検査は 0/1 で扱う（167 µs の周期で約 1 %） | 読み |
| B-1 | major | 表の外挿（3 kHz より上の f⁻² を無限まで）が、γ = 1 で L/13.2 m 倍（F_liu20 で 96 倍、F_wg100 で 240 倍）されて、ループの帯域の外の分散を支配する: 3 kHz より上 161 mrad、10 kHz より上 88 mrad（F_liu20 ポート 1）。f_N（500 kHz）より上の 12 mrad は、厳密に標本化した OU と台帳の行の両方に入る（二重） | `p6_wp6ab_table_psd.py --no-time-loop` |
| B-2 | major | 25・40 kHz の山がモデルの DAC 量子化であることは確認（量子化を切ると 8.2 → 1.76、5.0 → 2.25）。ただしモデルの寄与は白色の量子化雑音の見積もりの 2.6 倍（25 kHz）で、ディザも出力の雑音も無い決定論的な量子化器のリミットサイクルによる。DSC1 の出力の 3 kHz より上はオシロの床に埋もれていて未検証。主モードの Q = 25 の当てはめはこの仮定に条件付き（白色なら山 ≈ 4 で、Q はもっと高い必要がある） | `p6_wp6ab_quant.py` |
| B-3 | major | 5.3 kHz の構造は \|S\|² からは一意に決まらない。ループのどこに置いた直列の biquad でも同じ \|S\|² になる。遅延・実極は除外できる（山 1.8–2.3、幅 2 kHz、K_d で周波数が動く）。当てはめ（Q 15、r 0.11）は高さは合うが半値幅 125 Hz（実測 ≈ 300 Hz）。Q 5・r 0.34（軸の範囲の外）なら高さ・幅・S = 2 V の値がそろって合う。オーバーレイの「Q ≈ 18 は \|S\|² の半値幅から」は閉ループの山の幅とモードの Q の混同 | `p6_wp6ab_multimode_alt.py` |
| B-4 | minor | D89-4 の「r に非常に敏感（0.08 → 8、0.10 → 21、≥ 0.15 で不安定）」は、当てはめた点の解析値（0.08 → 5.0、0.10 → 8.7、0.15 でも安定で 235）と合わない。時間ループの値か別の点の値なら明記 | `p6_wp6ab_multimode_alt.py` |
| B-5 | minor | 25 kHz の実測は K_d に対して単調でない（7.0、15.4、13.6、9.1）が、モデルは 7.2、6.1、5.7。「25 kHz で 7.2（実測 7–15）」の書き方はこれを隠す。6–20 kHz の \|S\|² はモデル 1.22–1.34、実測 1.48–1.77 で、凍結した名目（1.48–1.75）より悪くなった。実測の \|S\|² は 14 区間の Welch（統計の誤差 ≈ ±27 %）なので、r を 10 % の精度で合わせるのはデータの精度を超える | `p6_ab_check.py`、`p6_wp6ab_multimode_alt.py` |
| B-6 | minor | 腕の分け方 ρ = 0 は、共通の拾い（ρ > 0）があれば腕ごとの PSD の下限（腕ごと = T/[2(1 − ρ)]）。ρ は 0 に固定されていて軸がない | `p6_wp6ab_table_psd.py`（時間ループ部分） |
| B-7 | major | `controller/dsc1` の latency 5.7 µs は、Ref D では DSC1 の入出力の極（120 kHz・100 kHz）と組の値（refD.yaml「5.7 µs with the 120 kHz / 100 kHz poles, 8.3 µs pure delay」）。F_* ＋ DSC1 には極が無い。27 点中 3 点で判定が変わる（極を足すと不安定 17 → 14）。項目 6 の h より効く | `p6_wp6ab_dsc1_h.py` |

## 2. A-1: Floquet・周期定常の分散・Dick の項

```
uv run python reviews/scripts/p6_wp6g_floquet_indep.py
uv run python reviews/scripts/p6_wp6g_timeloop.py --duration 0.1
```

### 2.1 定式化の確認（読み）

- 状態 z、入力 d（ロックする組の自由走行の位相）、e = Cz + d の線形周期系で、ゲートの中は A_on、外は A_off。周期 N = round(T f_sim) 標本の M = Π A_n、**ρ = max\|eig M\| < 1 で安定**。分数遅延は u の履歴の線形補間（`delay_buffer` linear と同じ）、遅延 ≥ 1/f_sim。PID は `pid.first_order` の式（dsc1: I[n] = I[n−1] + K_i h c[n]、tustin: 台形）。zero は c = off_gain·e、hold は最後の有効な e。制御器のゲートは lock 窓を gate_delay_s だけ遅らせ、先頭の ⌈settle_s f_sim⌉ 標本を捨てる — `gated_hold` の「開いてから k·h ≥ settle_s」と同じ。
- 周期定常の応答: d = e^{iωn} に対し x_N = M x_0 + (強制応答) = e^{iωN} x_0 → x_0 = (e^{iωN} − M)⁻¹·(強制応答)。正しい。
- 分散: 実の LPTV 系では e_{−f}[n] = conj(e_f[n]) なので、片側 PSD（規約 §7）で Var e[n] = ∫_0^{f_s/2} S(f)\|e_f[n]\|² df。正しい。
- 変換（Dick）の項: e_f[n] = e^{iωn} p_f[n]（p_f は周期 N）、\|e_f\|² = \|p_f\|²、⟨\|p\|²⟩ − \|⟨p⟩\|² = Σ_{k≠0}\|c_k\|²（Parseval）＝ 入力 f が出力 f + k/T に移る電力。定常な S で積分すると、ゲートによる折り返し（k/T 付近の雑音の直流への折り返しと、その逆）の全体。ゲートなしで 0。正しい。
- 保持中のドリフト: D(τ) = ∫S·4 sin²(πfτ) df は定常過程の構造関数 E[(d(t+τ) − d(t))²]（片側 PSD で正しい）。開ループの近似で、窓の終わりの残差を含まない（docstring に明記）。

### 2.2 独立な再実装との比較

再実装は scipy.signal.bilinear の区間（PZT のモードは fs_eff = ω0/[2 tan(ω0 h/2)] で f_res に prewarp）、自前の遅延線・PID・ゲート。TEC の経路は時間ループの配線どおり（遅延前の u から、PZT のモードを通らない）にした — 検査はこの経路も PZT のモードと遅延に通すが、結果は同じ。

| 点（F_liu20 ＋ borealis_si が基準） | 検査の乗数（上位 3） | 再実装 |
|---|---|---|
| 名目（zero、f_u 1 kHz、τ_pr 6.25 µs、settle 1 µs） | 0.999995, 0.999995, 0.6614254 | 同じ |
| f_u 5.5 kHz、5 kHz Q 100 | 2.8041215, 1.8654847, 0.999995 | 同じ |
| f_u 30 kHz、60 kHz Q 3 | 1.3434422 ×2 | 同じ |
| T 10 µs、D 0.9 | 0.9999995 ×2, 0.9630541 | 同じ |
| T 1 ms、D 0.1 | 0.99995 ×2, 0.5343117 | 同じ |
| hold、T 2 ms、D 0.02、f_u 3 kHz | 10.7484704 | 同じ |
| dsc1、遅延 6.0 µs（分数）、f_u 5.5 kHz | 1.4037299 ×2 | 同じ |
| 漏れ off_gain √10⁻⁵、遅延 2.4 µs | 0.6617076（速い乗数） | 同じ |

- **TEC の段の rate_div 10000（100 Hz の保持）**: lcm(T, 10 ms) の上でモノドロミーを取って ρ^(1/周期数) で比べた。名目は変化なし。最悪（f_12 = 0.3 Hz、τ_th = 1 s、T = 1 ms）で遅い乗数が 0.9995004 → 0.9995053（減衰率の 1 %）。速い乗数は変わらない。**無視してよい**（D90-8 の主張どおり）。
- 周期定常の分散（AR(1) の外乱、周期 Lyapunov で時間領域の厳密値）: var_meas・var_lock・var_period が 4 点すべてで **0.00 %** 一致。
- Dick の項: 純音（2、7、10.05 kHz）を時間領域で 400 周期流して、トーンの bin に残る電力の割合 = 検査の \|c_0\|²/⟨\|p\|²⟩（0.99497、0.98689、0.81945、5 桁一致）。

### 2.3 実際のプリセットの時間ループとの比較

単体テストはおもちゃのループ（利得 1 の植物）と比べている。ここでは F_liu20_nominal ＋ gate/borealis_window の時間ループそのもの（liu の縞、ADC、gated_hold、TEC の段、PZT、遅延したゲート）と比べた。外乱は同じ seed のフィードバックを切った run の `lock.L1.err` の Welch PSD。

| 点 | 時間ループ var(測定窓) | 検査 var_meas | 差 |
|---|---|---|---|
| borealis_si（zero、65/35 µs）、雑音・量子化 off | 7.983e-3 rad² | 7.999e-3 | −0.2 % |
| 同、全部 on | 7.986e-3 | 7.999e-3 | −0.2 % |
| zero、T 1 ms、D 0.1 | 8.646e-3 | 8.499e-3 | +1.7 % |
| hold、T 300 µs、D 0.3 | 1.237e-2 | 1.331e-2 | −7 %（ロック窓も −7 %、100 ms の統計の範囲） |

不安定点（hold、T 1 ms、D 0.1、検査 ρ = 2.42）は、最初の周期から誤差が 1.3 rad（窓の位相が測定窓から始まるので最初の 0.9 ms は保持）で非線形の領域に入り、成長率の定量比較はできなかった（有界に飽和し `valid = False`）。成長率はおもちゃのループの単体テストと、独立な再実装（§2.2）で確認済みとする。

## 3. A-2: プローブのゲートの物理

- **縞 ∝ √g、ショット雑音 ∝ g**: Liu 型の誤差信号は、ゲートしない強い参照（P_ref = 1 mW）と弱いプローブ（gP̄）のビート 2ℜG√(gP̄·P_ref)·V（DC ホモダインでも、ヘテロダインの振幅でも同じ）なので振幅 ∝ √g は正しい。雑音 2eℜ(P_ref + gP̄) も正しく、g = 0 で参照のショット雑音と電気雑音が残る（`in_gate` の docstring どおり）。既存の `gating`（パルスのデューティ）とは掛け算で重なり、二重にならない。
  - 注記: 片側の検出器ならプローブ自身の DC gP̄ の段差が窓の端で誤差に乗るが、`homodyne_dc` の v_off は定数。P̄ = 1 µW・P_ref = 1 mW で縞の振幅の 1.6 % 以下、平衡検出なら不均衡の分だけ。無視できる。
  - 注記: F_common のプローブはもともと AOM（probe_duty）を通るので、窓のゲートで周波数のずれが新しく入ることはない。
- **飛行時間 τ_pr,k**: `tau_prk_s = [n_g(L_prsrc + L_prs + L_k) + 1.0003 L_d]/c` は pr1 のビームの経路（PRsrc → PRS1 → F1（逆向き、fast 軸）→ D1）と一致（fast 軸の群屈折率の差は無視できる）。gd1 で `gate.gw.lock` を遅らせて liu の `in_gate` と `gated_hold` の両方に使うので、検出器に着くパワーと制御器の窓が同じ遅れ。検査も `gate_delay_s = tau_prk_s`。F_liu20 で τ_pr = 6.25 µs、F_wg100 で ≈ 15.6 µs なので、測定窓の最初の τ_pr の間は（飛行中のプローブで）制御器が動き続ける。これは物理的に正しい（その光はロック窓の間に出た）。マスクはその時間も測定窓として数える — GBS が見る真の位相なので正しい。
- **マスクの一貫性（カーネル内）**: `metric_i/ii` の積算は同じ時刻の `gate.gw.meas` を読む（ゲートのブロックが先に書く）。settle の標本は遅らせたロック窓の先頭にあるので測定窓に入らない（T = 10 µs・τ_pr > T のように回り込む場合も、検査と時間ループは同じ規則）。**カーネル内の積算（報告する指標の値）は正しい**。問題は記録した系列の方（G-1）。

## 4. G-1・G-2: マスクと記録の系列、窓に標本の無い点

```
uv run python reviews/scripts/p6_wp6g_mask_decim.py
```

- **A**（F_common の記録そのもの: decim 10、antialias none、T = 10 µs）: 記録したマスクは周期ごとに同じ位置で点標本化されるので、400 値すべて 0。マスクした系列は 0 標本（std NaN、エラーなし）、カーネル内は 2400・1200 標本で正しい。周期と decim が通約できると、窓によらず系列は空か全部になる。decompose の `validity` も記録したマスクを使うので、窓の中の飽和・ロック外れを**黙って見逃す**。
- **A（D_lock = 0.9）**: 測定窓 1 µs − 端 0.1 µs に標本が無く、run が `ValueError: PhaseStats.from_acc: no samples accumulated` で止まる。gate/borealis_window の 63 点の格子に含まれる点で、`runner.sweep` には点ごとの例外処理がない（G-2）。`chk.gated_stability` も `n_meas = 0`（var_meas は NaN）になるだけで CHK を出さない。
- **B**（T = 100 µs、D 0.65）: boxcar は端の値 0.4 が捨てられるので偶然混ざらない（ただし系列は平滑化される）。**fir（161 タップ、1.6 周期）では残した 118 値のうち 78 値のマスクが (0.5, 1) — ロック窓の標本を含む値を測定窓として数えている。**

直し方（案）: マスクの系列を「デシメーションの区間の全標本が窓の中」（最小値・論理積）で記録する。少なくとも `mask_slot` と `decim > 1` かつ `antialias != none` の組み合わせは ConfigError、`antialias none` で T·f_sim/decim がほぼ整数なら警告。窓に標本が 0 の点は、`chk.gated_stability` が無効化する CHK（例 `gate_no_meas_sample`）を出し、`evaluate_metrics` は例外でなく NaN の統計と `valid = False` で点を残す（D75）。

## 5. A-3: 上書きの組み合わせ（D90-4 の限界と、記載のないもの）

```
uv run python reviews/scripts/p6_wp6g_compose.py
```

| 組み合わせ（順） | 結果 | 物理への影響 |
|---|---|---|
| gate → controller/dsc1、moku_delta | ConfigError（大きな声） | — （メッセージが不親切、G-6） |
| detect_heterodyne ↔ gate（両順） | ConfigError | — |
| **controller/dsc1 → gate** | **読み込める** | ADC: dsc1 ±4 V → fpga_generic ±0.5 V（16 bit）。F_liu20 名目の K = liu_k_v = 6.3 mV/rad で、ADC の LSB の位相換算が 19.3 → 2.4 mrad（実機の DSC1 なら ±½ LSB ≈ ±10 mrad の不感帯・リミットサイクルがありうる、ロック経路の予算 5.9 mrad と同じ桁）。DAC・PID は dsc1 のまま。`chk.gated_stability` は tustin（時間ループは dsc1）、5.7 µs の遅延で 27 点中 1 点が楽観的（ρ 1.0000 対 1.0008） |
| controller/moku_delta → gate | 読み込める | ADC の device が moku_delta → fpga_generic（値は同じ。雑音を on にすると device の雑音の値が変わる） |
| **lab/multimode915b_future ＋ gate（両順）** | **読み込める**（D90 に記載なし） | `chk.lock_stability` は 4 モードの和、`chk.gated_stability` は 1 モード（主モード、r = 1）のまま。当てはめた点で f_u = 2 kHz は検査「不安定」・実際は安定、3 kHz で ρ 1.33 対 1.09。走査（f_u × 主の Q × mm2_r × mm4_r、108 点）で **6 点が「安定」と誤判定**（逆は 0） |
| **detect_heterodyne ↔ freq_1THz** | 読み込める | 後勝ち: freq_1THz → detect_heterodyne で 80 MHz、逆で 1 THz。プローブの周波数に依る熱の係数（2πΔν κ_τ L）と迷光の分離が黙って変わる |

`_dynamics_d`（`gated_loop_stability.py`）は各モードを `modes_ss([{..., "r": 1.0}])` で直列につなぐので、r の付いた表を写しても和にならない（`loop_stability._dynamics` の規則と違う）。

**直し方（案、根本で 1 回）**:
1. 制御器の選択を F_common の `offline.params` に出す（`adc_device`・`adc_bits`・`adc_range_v`・`pid_law`・`pid_h_s`）。adc0/adc1・pid*f・`chk.lock_stability`・`chk.gated_stability` はそれを `!expr` で読み、`controller/*` はパイプラインのキーではなくこの params を書く。ゲートの上書きも ADC を同じ `!expr` で書き直せば、順序によらず値が失われない。`pid_discretisation: !expr "offline.params.pid_law"`。
2. 同じく、アクチュエータのモードの表を 1 か所（例 `offline.params.pzt_modes`）に置き、pzt0/pzt1・`chk.lock_stability`・`chk.gated_stability` が読む。`_dynamics_d` に「全モードに r があれば和」の規則を入れる（`loop_stability._dynamics` と同じ）。
3. 枠組みの側の安全装置: `config.load` で、ある上書きの `!replace` が、それより前の上書きが変えた節点を置き換えるとき、または 2 つの上書きが同じ `offline.params` のスカラーに違う値を書くときは ConfigError（両方の上書きの名前を出す）。意図的なら明示のタグ（例 `!override`）。D80 の「上書きが自分のキーを上書きしないかの検査」を上書きの間に広げる形。G-5（D90-6）もこれで大きな声になる。

## 6. B-4: 表の PSD（`mech_shape: table`）

```
uv run python reviews/scripts/p6_wp6ab_table_psd.py            # 時間ループ（2 × 2 s）込み
uv run python reviews/scripts/p6_wp6ab_table_psd.py --no-time-loop
```

- **規格化**: `table_fit` の重み S0_k は `lowpass` の低域の片側 PSD（S0/(1 + (f/f_k)²)）で、`_lgs` の離散化は片側 PSD \|H\|²（規約 §7）。解析の OU の和 / T は 30 Hz–100 kHz で 0.945–1.027（4 個/decade）。refD の時間ループ（フィードバック off、片側 Welch）でも腕の和 / T = 0.92–1.02（100 Hz–10 kHz、Welch の散らばりの範囲）。正しい。
- **腕の分け方**: γ = 1 で S_a1 : S_a2 = 0.77–0.85 : 0.14–0.17（11.1 : 2.1 = 5.3 : 1 どおり）、γ = 0 でそれぞれ ≈ 0.5、組の PSD = 腕の和（ρ = 0、比 0.94–1.05、10 kHz ではレーザーの FM が足されて 1.24）。実測の PSD は φ_lo − φ_pr（両腕の差の位相）なので「無相関なら和」の扱いは正しい。ただし共通の拾い（テーブルの振動）で ρ > 0 なら腕ごとの PSD は T/[2(1 − ρ)] で、ρ = 0 は下限（B-6）。将来のプリセットに 1 本ずつ移すときに効くので、ρ の軸を残すべき。
- **V0²/2 の換算**: 自由走行の記録の中で縞を何度も横切る前提。ただしロック時のループ内の PSD（\|S\|² が分かっている線形の領域）を T が 0.96–1.48 倍で再現する（`p6_ab_check.py`）ので、T の大きさは独立に支えられている。
- **外挿（B-1、major）**: 100 Hz より下の f^−1.5（10 Hz まで）と平坦はループの利得で抑えられるので影響は小さい。問題は 3 kHz より上: T の分散は 3 kHz より上 2.7e-4 rad²、10 kHz より上 8.1e-5 rad²（Ref D の 13.2 m の和）。これがループの帯域（f_u ≈ 1 kHz）の外なので抑えられず、`lab/acoustic_future` の γ = 1 では L/13.2 m 倍される:

  | ファイバ | 倍率 | 3 kHz より上 | 10 kHz より上 | f_N = 500 kHz より上（二重に数える分） | 10 kHz から f⁻⁴ にした場合の 3 kHz より上 |
  |---|---|---|---|---|---|
  | F_liu20 ポート 1（1268 m） | 96 | 161 mrad | 88 mrad | 12.4 mrad | 144 mrad |
  | F_wg100 ポート 31（3170 m） | 240 | 255 mrad | 139 mrad | 19.7 mrad | 228 mrad |
  | F_wg6/F_bor6 ポート 31（1060 m） | 80 | 147 mrad | 81 mrad | 11.4 mrad | 132 mrad |

  3 kHz より上の値は Ref D でも測っていない（50–200 kHz の床 2.2e-9 rad²/Hz はレーザーの FM で、T はその下にあるという制約だけ）。γ = 1（ファイバ 1 m あたりの拾いが Ref D と同じ）は箱に入れた km のスプールでは過大の可能性が高いが、軸としては残してよい（D34）。**ただし結果がこの外挿でほぼ決まるので、(a) 高域の傾き・折れ点を軸にする、(b) γ = 1 の結果には「3 kHz 以上は外挿」と必ず付ける**こと。
- **f_N より上の二重計上（minor、既存の two_pole と共通）**: OU は Van Loan で厳密に標本化されるので、標本の分散は f_N より上の分も（折り返して）含む。台帳の行（`bands`、psd_outside = T）はさらに f_N より上を足す。Ref D では 1.3 mrad で無視できるが、γ = 1 の F_* では上表の 11–20 mrad。table_f_max ≪ f_N の docstring の注意は「合成と台帳が一致する」ことしか保証しない。

## 7. B-5: 複数モード

```
uv run python reviews/scripts/p6_ab_check.py            # 64 s
uv run python reviews/scripts/p6_wp6ab_quant.py
uv run python reviews/scripts/p6_wp6ab_multimode_alt.py
```

- **D89-4 の再現**: 利得余裕 3.3 dB（5412 Hz）、時間ループの \|S\|² 5.3 kHz 11.3/10.6/6.4（S = 0/0.5/1.5 V）、25 kHz 7.2、10 kHz 3.1、40 kHz 5.0、ループ内 300 Hz・1 kHz は実測の 0.96–1.48 倍（S = 2 V は走らせていない）。**主張どおり**。
- **DAC 量子化の主張（B-2）**: 同じ seed で dsc1.dac の量子化を切ると 25 kHz 8.18 → 1.76、40 kHz 5.02 → 2.25（5.3 kHz は 11.2 → 11.1）。ADC の量子化・雑音は寄与しない。主張はモデルの中では正しい。しかし白色の量子化雑音（Δ = 8 V/2¹²、500 kHz、Δ²/12/(f/2)）をアクチュエータの経路 \|S·P\|² に通した見積もりは 25 kHz で +2.45、40 kHz で +1.29 で、モデルの +6.4、+2.8 の 1/2.6・1/2.2。モデルの DAC はディザも出力の雑音も無い決定論的な量子化器で、D61-3 のリミットサイクル（±1–2 LSB）が山を作っている。実機の DSC1 の出力の 3 kHz より上は測れていない（COMPARISON §4）。白色なら 25 kHz の山は ≈ 4 で、実測 7–15 に合わせるには主モードの Q を上げる必要がある（作者の「解析だけの当てはめなら Q ≈ 93」と同じ方向）。**Q = 25 は「DAC の量子化がモデルどおりのリミットサイクル」という仮定に条件付きの値**で、Q の軸（10–300）は残すこと。
- **Σr = 1、r_main = 0.67 の物理的なもっともらしさ**: PZT に巻いたファイバはアクチュエータとセンサが同じ場所（collocated）なので、留数は正で極と零点が交互に並ぶのが典型。当てはめの r がすべて正なのはこれと合う。有限個のモードで Σr = 1 にすると、高次のモードの静的な寄与（残留コンプライアンス）を全部この 4 モードに押し込むことになるので、r_main = 0.67 は「主モードが DC の利得の 2/3」ではなく「この 4 モードの表現での値」と読むべき。5.4 kHz に r 0.11 の低いモードは、915B の筐体・取り付けのモードならありうる。
- **非一意性（B-3）**: \|S\|² はループ全体の L だけで決まる。
  - (a) 遅延（最良 Δτ = 19 µs）・(b) 実極（最良 10.9 kHz）: 5.3 kHz の山は 2.3・1.8 にしかならず、幅 2 kHz、S を変えると 5.35 → 4.9 kHz に動き、3 kHz の \|S\|² を 1.3–1.5 に壊す。**実測の「K_d によらず 5.3–5.45 kHz に固定された狭い山」は共振でしか作れない**（除外できる）。
  - (c) 直列の biquad ((1 − r)(s² + sω/Q) + ω²)/(s² + sω/Q + ω²) を**ループのどこに置いても**（914 の駆動の電気的な共振でも、検出の側でも）\|S\|² は同じ形（17/15.6/8.2/4.9、10 % の差は和と積の違い）。\|S\|² からは共振の場所を決められない。場所はシミュレータでは効く: アクチュエータの側なら DAC の雑音が共振で持ち上がり、真の（ループ外の）位相誤差は S·d。検出の側ならループ内の信号と真の位相が違う。
  - r と Q も一意でない: 当てはめ（Q 15、r 0.11）は高さ 12.4 だが**半値幅 125 Hz**（実測 ≈ 300 Hz、2.4 倍狭い）。Q 5・r 0.34（軸の範囲 Q 8–40・r ±0.15 の外）なら高さ 13.9、幅 337 Hz、S = 2 V で 4.5（実測 3.8）とそろって合う。閉ループの山の幅は閉ループの極の減衰（利得余裕）で決まり、モードの Q ではない — オーバーレイの「Q ≈ 18 は \|S\|² の FWHM から」は混同（B-3）。r の感度の表（B-4）: r 0.06/0.08/0.10/0.11/0.12/0.15 → \|S\|² 3.2/5.0/8.7/12.4/19.0/235（0.15 でも安定）。
- **区別する測定**:
  1. 914 → 915B の開ループの伝達関数（DSC1 か Moku から正弦波を掃引し、ループを 10 Hz 程度の遅いロックにして、縞の傾きで位相の応答を読む）: アクチュエータの経路だけの f_i・Q_i・r_i と遅延・極が直接決まる。共振が駆動・検出の側にあれば現れない。
  2. 閉ループでの注入（DSC1 の入力か出力の加算点に掃引した正弦波、L = −B/A）: コヒーレントな \|S\|² でループ内の雑音を含まない — 25・40 kHz の山が量子化（非コヒーレント）か \|S\|²（コヒーレント）かを分ける。補助: 914 の利得と DSC1 の利得を逆比で変える（L は同じ、位相換算の DAC の LSB だけ変わる）。
  3. 検出の側の確認: 915B を使わずに光の位相を変える（レーザーの周波数変調を ΔL で位相に変える）と、PDB450C → DSC1 の ADC の伝達が平坦か分かる。

## 8. 項目 6: `controller/dsc1` の h = 1/f_sim と装置の 500 kHz（B-7）

```
uv run python reviews/scripts/p6_wp6ab_dsc1_h.py
```

F_* の微調整の段は I だけ（K_p = 0）なので、dsc1 の則は K_i(h/2 + 1/s)、つまり位相の進み arctan(ωh/2)（10 kHz で h = 1 µs 1.8°、2 µs 3.6°；30 kHz で 5.4°・10.7°）と高域の利得 K_i h/2。

| 変種 | 内容 | 27 点（f_u × f_res × Q）の不安定 | 名目 (1 kHz, 20 kHz, Q 10): 余裕 / max\|S\|² | (5.5 kHz, 32.5 kHz, Q 17) | (30 kHz, 60 kHz, Q 3) |
|---|---|---|---|---|---|
| V0 出荷 | h 1 µs、遅延 5.7 + 0.5 µs、極なし | 17 | 7.5 dB / 3.28 | 不安定（−3.5 dB） | 不安定 |
| V1 | rate_div 2、h 2 µs、τ を補償（6.2 µs） | 17（V0 と同じ点） | 7.2 dB / 3.43 | 不安定（−4.7 dB） | 不安定 |
| V2 | rate_div 2、h 2 µs、補償なし（6.7 µs） | 17（同じ） | 7.5 dB / 3.31 | 不安定（−3.5 dB） | 不安定 |
| V0+p | V0 ＋ 入出力の極 120/100 kHz | 14 | 10.2 dB / 2.50 | 安定（4.7 dB、14.1） | 不安定 |
| Vdev | Ref D で確かめた DSC1: h 2 µs、τ_L 5.7 µs、極 120/100 kHz | 15 | 9.2 dB / 2.73 | 安定（2.2 dB、67.7） | 不安定 |

（余裕は 100 Hz より上の −180° の交差での値。TEC の二重積分が作る低域の交差は除いた。）

**回答**:
1. **h だけを 1 → 2 µs（rate_div 2）にしても、27 点の判定は 1 点も変わらない**（V0・V1・V2 は同じ 17 点が不安定）。\|S\|² の山は名目で 3.28 → 3.43（+5 %、V2 では 3.31）、名目の余裕 7.5 → 7.2 dB（V2 では 7.5 dB のまま）。不安定な点の \|S\|² は意味がないので比べていない。
2. τ を補償しない rate_div 2（V2）は今の V0 とほぼ同じ: I だけの制御器では h の増加 Δh による進み ωΔh/2 が保持の遅れ Δh/2 の位相 ωΔh/2 と一次で打ち消す。補償する（V1、Ref D の `dsc1_explicit_delay_s` の流儀）と K_i h/2 の分だけ高域の利得が増えて少し悲観的。
3. **h より効くのは、F_* ＋ DSC1 に入出力の極が無いこと**（B-7）。5.7 µs は Ref D で 120/100 kHz の極と組の τ_L で、極なしなら 8.3 µs の純遅延が組。V0 と Vdev は 2 点で判定が違い（どちらも V0 が悲観的）、名目の余裕は 7.5 対 9.2 dB、(5.5 kHz, 32.5 kHz, Q 17) では V0 が不安定・Vdev が安定で \|S\|² 67.7（境界）。
4. **Ref D**: すでに pid・dac が rate_div 2、`h_s` = 2 µs で装置と一致しているので**影響なし**。仮に h = 1 µs にすると、K_i・K_d（1/s、s）を固定なら余裕 3.3 → 3.2 dB、\|S\|² 5.3 kHz 12.4 → 12.9、25 kHz 1.76 → 1.39 で小さいが、DSC1 に入れる値（D53 の I・D、標本あたり）を固定すると K_i が 2 倍・K_d が 1/2 になって**不安定**（−3.1 dB、10 kHz の \|S\|² 209）。h は装置の更新レートと D53 のゲインの単位に結び付いている。
5. **推奨（ユーザの判断の材料）**: rate_div 2・h = 2 µs を採るなら、Ref D と同じ DSC1 の鎖にそろえる — `tau_s = τ_L − 1/f_sim − (保持の (r − 1)/(2 f_sim))`、入出力の極（`lpf` 120 kHz・100 kHz）を時間ループと `chk.lock_stability` に足す（または極なしで 8.3 µs）。h だけ変えても判定は変わらないので、h だけの変更は優先度が低い。なお `gated_loop_stability_check` は全段を f_sim で回すので、PID の rate_div 2 は表せない（I だけなら 2. の打ち消しで差は小さいが、採るなら検査にも PID の rate_div を入れる必要がある）。

## 9. 確認して問題がなかったこと

- `loop_stability_check` の和の形（D89-1）は `modes_ss` と同じ式で、r の無い表は従来どおり積。単一モード r = 1 は両方で同じ。
- `gated_hold` の zero ＋ I だけの tustin PID は出力を厳密に保持する（窓が閉じてから 2 標本目以降。台形の最後の半ステップ）。検査も c[n−1] を状態に持つので同じ。
- `gate_schedule` の周期の始まりの丸め（1e-9 周期）、ジッタの乱数の消費（毎ステップ 1 個、G11）、`anti_sync` の窓の位置（P3 の 180°）。
- 迷光: ゲートの上書きはオフラインの迷光の表を変えないので、`stray_light_table` の値はプローブが常に入っている場合のまま（保守的）。窓の中の値は `window_stray_light`（D88、別のレビュー）で出すこと — 上書きのコメントに記載あり。

## 10. 直すときの優先順位（提案、決めるのはユーザ）

1. G-3・G-4・G-5: §5 の「1 か所の params ＋ `!replace` の上書きの検出」。Phase 6 の軸 1 の候補（D85-12: 制御器 DSC1・Moku ＋ P2/P3/P8 のゲート、915B の複数モード）はまさにこの組み合わせを使う。
2. G-1・G-2: マスクの記録を論理積に、窓に標本が無い点を CHK にして残す。
3. B-1: 表の高域の傾き・折れ点の軸、acoustic_future（γ = 1）の結果に外挿の注記。
4. B-7・項目 6: F_* の DSC1 の鎖を Ref D と同じにする（極と遅延の勘定）。
5. B-2・B-3: §7 の 3 つの測定（特に 1 と 2）。それまで主モードの Q と副モードの (Q, r) は「当てはめ」ではなく「実測と矛盾しない一例」と書く。
