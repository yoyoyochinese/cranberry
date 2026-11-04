# Phase 5 レビュー: WP6-I（D93、P7 のロックしたままの同定）と WP6-X（D94、軸 1 のベンチ構成とプローブの上書き）

- 日付: 2026-11-04 ／ 対象: dev @ 56e1167（`git merge --ff-only dev` 後）
- 読んだもの: AGENTS.md、DECISIONS（D55・D61・D82・D85・D86・D89・D91・D93・D94）、`notes/conventions.md`、`notes/lit-kashiwazaki2023.md`、
  `notes/fact-p6-thermo-optic-moku.md`、`notes/lit-svarc2023.md` §5、`notes/lit-luo2025.md`、`reviews/phase6_wp6g_ab.md`（形式）、
  `blocks/ctrl/excitation.py`、`blocks/meas/loop_identification.py`、`hardware/interface.py`・`sim.py`・`moku.py`、`runner/run.py`（フック）、
  `presets/axes/ident/*`、`presets/bench/*`、`presets/axes/bench/*`、`presets/axes/probe/{actuator_pathlength,bad_kT_kappa_independent,two_tone,detect_heterodyne_refshift}.yaml`、
  `presets/reference/{refA,refC,refD}.yaml`、`presets/future/F_common.yaml`、`presets/axes/controller/dsc1.yaml`、`devices/dsc1.yaml`、
  `reviews/scripts/p6_wp6i_ident.py`・`p6_bench_smoke.py`・`p6_p4a_check.py`、`tests/integration/{test_phase6_ident,test_p6_bench,refD_support}.py`、
  `tests/unit/{blocks/ctrl/test_excitation,hardware/test_sim_ident,meas/test_loop_identification}.py`、`refs/manuals/MokuDelta-Specifications.pdf` p30。
- **`blind/` は読んでいない。`predictions/` は何も書き換えていない。プリセット・ソース・DECISIONS.md は変更していない。**
- 高速のテスト: `uv run pytest -m "not slow"` → **1673 passed, 8 skipped**, 288 deselected（56 分、§5）。
- 根拠のスクリプト（`reviews/scripts/`、すべて読むだけ）:
  - `p6_wp6i_ident.py`（作者のもの）の part 1・4・5 を実行して D93-5 を再現
  - `p6_wp6i_x_checks.py` … Moku の加算点の換算の感度（係数の誤り、経路の遅延）、p7_literal の TEC 保持の減衰、wire_F ＋ controller/dsc1、量子化と粗い掃引、余裕の「最初の交差」
  - `p6_bench_smoke.py`（作者のもの）を fpga・`bench/dsc1_refA`・`bench/moku_delta_refA` で実行、`p6_p4a_check.py` part 4 を実行
  - `p6_wp6x_smoke_settle.py` … smoke の表からロック獲得の過渡を除いた値
  - `p6_wp6x_fc_delay.py` … LO の遅延の不整合とフィルタキャビティ（FC）の遅延補償の閉じた式（シミュレータと独立）、DX1 の数字
  - `p6_wp6x_refA_demod_sign.py` … lock 2 の復調位相の吸収の符号（50 m で +X／0／−X）
  - `p6_wp6x_two_tone_bg.py` … 直交偏光のプローブ（速い軸）の群複屈折による E のポートごとのずれ

## 0. 結論

| 対象 | 判定 | 要点 |
|---|---|---|
| A-1 推定量（minus_l・plant）、コヒーレンス、余裕 | **正しい**（minor ×3） | L̂ = −resp/ref、S = ref/x、T = −resp/x、L̂ = −P̂·C₁ は道具変数の比として正しい。refD の検証値を再現: 雑音なしで厳密な z 領域と 0.010 %・3×10⁻¹⁰°、判定（dsc1）と 2.03 %・2.19°、雑音ありで 4.33 %・2.26°、f_u 2089 Hz・PM 77.8°、追加掃引で f_180 16.79 kHz・GM 12.92 dB（厳密 16.78・12.79）。粗い 25 点の掃引の余裕は共振の近くで 2 dB 以上ずれる（I-3） |
| A-2 インターフェース、`SimController` | **正しい** | 追加はすべて既定値付き・非抽象で後方互換。`inject_topology` は「無いときだけ足す」なので runner と二重にならない |
| A-3 Moku の換算 L̂ = −H/(g+H) | **代数は正しい**（minor ×2） | e = −L(e + g x) から厳密。ただし「1 dB で 12 %」は \|L\| ≲ 1 だけ（\|L\| = 7 で −16〜−28 %、\|L\| → ∞ で L̂ → 8.2）。FRA のモード名が仕様と違う（I-1、I-2） |
| A-4 p7_literal の TEC 保持 | **一つの整合した読み**（minor） | 符号・配線は整合。ただし保持の帯域を f_12 に結んだので格子の隅で ζ = 0.036・PM 4°、時定数 0.8–83 s ≫ 0.65 s の掃引（掃引中はオフセットだけで動作点を保つ）。0.65 s の検証は保持の符号を試していない（I-5） |
| B-5 refA_fiber | 物理は**正しい**／**major 1** | SQ の係数・損失・復調位相の吸収（符号を 50 m で独立に確認）・肩代わりは正しい。smoke の表は自前で再現（9.37/9.48/10.88、9.03/5.97/20.57）。**ファイバ整合 100 m の 10.9 mrad と DSC1 の「1.5 mrad 悪い」はロック獲得の過渡**（最初の約 10 ms に EOM が ±0.5 V で飽和）。定常は 9.33/9.34/9.41 mrad で L によらない（X-1）。「自由空間は 10 m 付近が最良」は 3 点の中の最良で、閉じた式の最適は 17.8 m（X-2）。DX1 の 400 倍は正しいが、「LO の遅延整合 ≲ 1 ns」は FC があるときの要求として不正確（X-3） |
| B-6 DSC1・Moku（refA/refC） | 遅延の勘定は**正しい**（D61・D91 どおり）／minor | τ_L 8.3 µs − 1/f_sim − (r−1)/(2f_sim)、検査は delay 8.3 µs・dsc1・h 2 µs で時間ループと同じ全遅延。−180° ≈ 32 kHz・GM ≈ 4 dB を確認。「DSC1 は Ref A 型には遅すぎる」は K23 の MHz 級のループに対しては正しいが、refA のモデル（20 kHz 設計）では定常で 0.3–1.0 mrad 悪いだけ（X-6） |
| B-7 actuator_pathlength | **正しい** | g = 1 + f·Δν/ν0、c_r = g(k_T + 2πΔνR)L、ホップ 2π(1/g − 1)、R′ の f = 1.007 はすべて `fact-p6-thermo-optic-moku.md` §1・`lit-svarc2023.md` §5.1 から独立に導いた式と一致（p4a part 4 も通過）。軸の組み合わせの無駄（X-9） |
| B-8 two_tone | 偏りの式・格子は**正しい**／**major 1**（Open の確認） | −atan(r tan E) と 2δν = m f_rep/(N_b+n_guard)（m 偶数）は F_common の L_k = L_0 + k(N_b+n_guard)c/(f_rep n_g) と整合。直交偏光（速い軸）のプローブでは b_g = 1.36 ps/m が E に入らず、F_liu20 の 7 GHz でポート 2 は \|cos E\| = 0.17・偏り 0.85 rad（r = 0.2）。「全ポート共通で打ち消す」は同偏光のプローブでしか成り立たない（X-7） |
| B-9 refC_fiber・pulsed_wg_1port | **整合** | L_deployed = 0 は係数がすべて L_deployed 比例で 0、書き手は 0 の adder で、ゲート・台帳も矛盾なし。L_k の上書きは F1・LO1・τ_k・TEC 利得・迷光・SBS の全部に伝わり、port_k は無効にした demux だけが読む |

blocker はない。**D93-5 の数字はすべて再現し、成り立つ。** D94 の結果は、自由空間の列（9.0/6.0/20.6 mrad）と DX1 の倍率は成り立つ。**ファイバ整合 100 m の値と DSC1 の差は過渡を除いて直す**こと（X-1）。P8 の「格子上で全ポート共通」は、直交偏光のプローブについては群複屈折を入れるまで使わないこと（X-7）。

## 1. 指摘の一覧

| # | 重さ | 指摘 | 再現 |
|---|---|---|---|
| I-1 | minor | Moku の換算の「係数 1 dB の誤りで L̂ に 12 %」は \|L\| ≲ 1 の帯だけ。(1+L̂)/L̂ = (g′/g)(1+L)/L なので \|L\| が大きいほど誤りが増え、\|L\| → ∞ で L̂ → 1/(g′/g − 1) = 8.2。refD のモデルで +1 dB: \|L\| = 7 で −15.8 %、f_u +9 %、PM −7.1°、GM −0.8 dB。複素の g(f) を帯域全体で較正すべき（ループを開いて PID 入力のモニタと FRA 出力の比）。経路の遅延 10–100 ns は無視できる（0.1 % 以下） | `p6_wp6i_x_checks.py 1` |
| I-2 | minor | `loop_gain_plan` の FRA のモード "In1 (input / output normalisation)" は仕様に無い名前。仕様（p30）のモードは In/Out (dB)・In/In1 (dB)・In (dBm…)で、H = e/x に当たるのは **In/Out**。同じ頁の「Normalization（参照掃引）」でケーブル・前段の遅延を補正できることも書いておくとよい | PDF p30、`hardware/moku.py` |
| I-3 | minor | 余裕を粗い 25 点の対数掃引から補間で出すと、共振の近くで大きくずれ、印も付かない: F_wg6 ＋ controller/dsc1 ＋ wire_F で GM 9.85–9.93 dB・f_180 18.19 kHz（判定 7.53 dB・19.20 kHz）、16–22 kHz の追加掃引では 7.39–7.47 dB・19.28 kHz で一致。refD でも雑音ありの 25 点は 16.10 kHz・12.40 dB（厳密 16.78・12.79）。D93-5 は追加掃引の値を書いているので結論は正しいが、ブロックの出力の f_180・GM は「交差の区間が共振をまたぐ／区間の比が大きい」ときに印を付けるか追加掃引を促すべき | `p6_wp6i_x_checks.py 3 4`、`p6_wp6i_ident.py 1 5` |
| I-4 | minor | `loop_margins` は最初の −180° の交差の GM を返し、`loop_stability_check` は全交差の最大 \|L\| で判定する。refD_lab は −180° を 7 回横切り（5.41・5.69・9.9・10.5・21.9・35.9・39.9 kHz）、今は最初が最悪（0.682）なので一致するが、保証はない。全交差の最小の GM を返すべき | `p6_wp6i_x_checks.py 5` |
| I-5 | minor | p7_literal の保持ループ L = ω_h/(s(1 + sτ_th))、ω_h = 2π f_12: F_common の格子（f_12 3×10⁻³–0.3 Hz × τ_th 1–100 s）で ζ = 0.036–3.6、PM 4–89°、時定数 1/(2πf_c) = 0.8–83 s。0.65 s の掃引の間は実質オフセットだけで動作点を保ち、検証（1.1 %・3.6°）は保持の符号も減衰も試していない。読みとしては整合（符号は pid0f と同じで正しい）が、`ident_hold_f_hz` を f_12 から切り離して ζ ≥ 0.5 に選ぶ（ω_h ≤ 1/τ_th）か、論文どおり速いストレッチャで保持する変種を足すべき。ドリフトの大きい点では P̂ の縞の傾きの偏りに注意（ヘッダに記載あり） | `p6_wp6i_x_checks.py 2` |
| I-6 | minor | controller/dsc1 では 12 bit の DAC の LSB（8 V/4096 × 8.125 rad/V = 15.9 mrad）が励起の 20 mrad とほぼ同じで、25 点の掃引の大きさの誤りが 6.9 %・位相 6.3°（量子化を切ると 1.5 %）。ident の振幅は DAC の LSB の数倍にする規則が要る。（「制御器の上書きと組み合わせ可」という D93-4 自体は確認した） | `p6_wp6i_x_checks.py 3 4` |
| I-7 | 記録 | wire_F・p7_literal の励起は rate_div を持たない。D91-1（F 系の DSC1 を 500 kHz で回す）は `controller/dsc1.yaml` にまだ入っておらず（rate_div 1、ヘッダにも「500 kHz は課さない」）、入れるときは wire_refD と同じく励起の rate_div と記録の decim を合わせること | 読み |
| X-1 | **major** | smoke の表（D94 結果）のファイバ整合 100 m（10.88 ± 1.16 mrad）と DSC1 の 100 m（12.33 ± 1.29）は、**ロック獲得の過渡**を含む: 3 seed 中 2 つで最初の約 10 ms に lock 2 の EOM が ±0.5 V で飽和し（100 Hz の肩代わりが初期のファイバの揺れに追いつくまで）、その間の指標 (i) が上がる。t ≥ 15 ms では fpga 9.33/9.34/9.41 mrad、DSC1 9.66/9.64/9.76 mrad で **L によらず平ら**。ファイバ整合で 100 m に伸ばすと悪くなるという読みは成り立たず、DSC1 の差は 0.3–0.4 mrad（自由空間 100 m で 1.0）。smoke は t_settle を捨てて出し直すこと。`test_refA_fiber_100m_runs_and_offload_holds_eom` は既定の seed だけで EOM < 0.45 V を確かめている | `p6_wp6x_smoke_settle.py`（`bench/dsc1_refA` も） |
| X-2 | minor | 「自由空間は 10 m 付近が最良（FC の遅延の補償）」: 機構は正しい（閉じた式で 1/10/100 m を 9.03 → 6.00/20.29 mrad と予測、smoke は 5.97/20.57）が、最適は Δ = −0.69 τ_FC の **L = 17.8 m**（レーザーの項 4.84 mrad）。10 m は 3 点の中の最良 | `p6_wp6x_fc_delay.py` |
| X-3 | minor | ヘッダと D94 結果の DX1 の数字は加法の式 2π²h0\|Δ\| + πh0/Δν_FC を使うが、スクイーズ光の経路が長い（補償の向き）ときは成り立たない（DX1・自由空間 100 m で厳密 0.405 rad、加法 0.511 rad）。「LO の遅延整合 ≲ 1 ns」は FC があるときの要求として不正確: σ² ∝ h0 なので最適の Δ は h0 によらず −66.9 ns（スクイーズ光が長い）で、DX1 でも 156 → 97 mrad に下がり、最適の ±10 ns で 1 % しか変わらない。FC を使わないなら 10 mrad に 0.2 ns（1 ns で 22 mrad）。要るのは FC の再設計（Δν_FC）と、FC に合わせた遅延の置き方 | `p6_wp6x_fc_delay.py` |
| X-4 | minor | refA_fiber のヘッダ「残るのは帯域外の部分と、共通でない LO ファイバ（fiber_matched）」: lock 2 のセンサ `bhdbeat` は `beam.lo`（LO ファイバの後）と `pr_bhd` の差を見るので、LO ファイバの雑音も帯域内では打ち消される。100 m の増加は X-1 の過渡（音響を切ると 3 seed とも 9.39 mrad） | `refA.yaml` の bhdbeat（in_lo: beam.lo.phi）、`p6_wp6x_smoke_settle.py` |
| X-5 | minor | LO のストレッチャ（±32.5 rad、UNVERIFIED）: 915B 型の仮の値・交差 100 Hz・20 kHz Q 10 の共振（lock 2 の 20 kHz の単位利得でストレッチャの経路は 5 %）は妥当。ただし自由空間 100 m では SQ の熱のランプ 44 rad/(m·K) × 100 m × 10⁻⁴ K/s = 0.44 rad/s で、τ_th 300 s の過渡の後は約 1 分で範囲を使い切る。LO の 2π のリセット（BHD には見えない）か温度の段が要る。新しい lock 2 の安定判定は肩代わりの段を含まない（DSC1 の f_180 でストレッチャの経路は 0.2 %、影響なし） | 読み・計算 |
| X-6 | minor | DSC1 の判定: GM ≈ 4 dB（31–33 kHz）、12 bit の LSB が EOM で 1.75 mrad、2 MHz の復調不可はすべて正しい。ただし refA のモデル（20 kHz 設計）では定常で 0.3–1.0 mrad 悪いだけで、「遅すぎる」の根拠は K23 の MHz 級のループ（refA のプリセットにはない）。D94-3 の書き方にその前提を明記すべき | `p6_bench_smoke.py bench/dsc1_refA`、`p6_wp6x_smoke_settle.py bench/dsc1_refA` |
| X-7 | **major**（D94-5 の Open の確認） | two_tone のプローブは F_k の**速い軸**を通る（ビームの経路 `axis: fast`）が、E = 2πδν(τ_pr − τ_ref) は `Topology.timing`（n_g·L）で、群複屈折 b_g = 1.36×10⁻¹² s/m（F_common の pol_params、PM15-U25D）が入らない。実際にはポート k に 2πδν b_g(L_k − L_0) が加わり（F_liu20・7 GHz で 1 ポートあたり 2.4 rad、ポート 31 で 75.7 rad）、格子の条件が崩れる: ポート 2 で \|cos E\| = 0.17（ループ利得 ×0.17）、偏り 0.52 rad（r = 0.1）・0.85 rad（r = 0.2）。b_g の温度係数 K_τ = −1.78×10⁻¹⁵ s/(m·K) も ±δν の行に無い（ポート 31 で dE/dT = −0.10 rad/K、r = 0.2 で 20 mrad/K）。直し方: 格子をプローブの軸の群屈折率 n_g + c·b_g で作る（δν × 0.999722、7 GHz で −1.94 MHz）か、同偏光のプローブにする（今は pol_same と組み合わせられない）。前者は b_g を 1.3×10⁻³ の相対精度で知る必要がある（最長ポートで 0.1 rad） | `p6_wp6x_two_tone_bg.py` |
| X-8 | minor | two_tone のヘッダ「lock.L_k.err = −b at the lock」は符号が逆。テスト（`test_two_tone_bias_matches_formula`）は err ≈ +b を確かめ、ヘッダ自身の式（X̄ = π/2 − atan(r tan E)、誤差 = X̄ − π/2）も +b | 読み、テスト |
| X-9 | minor | actuator_pathlength の `R_s_per_m_K`（5 点）と `R_prime_s_per_m_K`（3 点）は `bias_model` によらず常に掃引の軸になり、3 × 5 × 3 = 45 点のうち異なるのは 9 点。条件付きの軸にするか、掃引の計画で固定すること。f は R で n_g/n = 1.0166、R′・implicit で 1.007（D86 のとおり R は ξ の分散なし）で整合 | 読み |
| X-10 | 記録 | refC_fiber: 机の上の 2 本のスプール（SQ・LO）の音響は refC のまま無相関（安全側）。pulsed_wg_1port: F_common の `port_k` の軸は F_wg6_nominal が 31 に固定しているので重複点はできない。どちらも直す必要はない | 読み |

## 2. A: D93（P7 のロックしたままの同定）

```
uv run python reviews/scripts/p6_wp6i_ident.py 1 5 4     # 36 s
uv run python reviews/scripts/p6_wp6i_x_checks.py        # 約 3 分
```

### 2.1 推定量（項目 1）

- **minus_l**: 加算点の後 u = e + x（reference）、前 e（response）、ループ e = d − L u（d は x と無相関の雑音）。u = (x + d)/(1 + L)、e = −L(x + d)/(1 + L) + …なので、x との相互スペクトル（または x でのロックイン）で L̂ = −S_xe/S_xu = L、S = S_xu/S_xx = 1/(1 + L)、T = −S_xe/S_xx = L/(1 + L)。x を道具変数にしているので d は偏りを作らずコヒーレンスを下げるだけ。コードの `lh = tt/s`（= −h_resp/h_ref）は正しい。
- **plant**: 解析器が C₁ を置き換え、x がアクチュエータの入力。e = P x + d（P は x → 誤差、保持ループの 1/(1 + L_hold) を含む）。閉じたときのループは u_act = C₁ e、e = P u_act なので 1 + L の形では L = −P C₁。コードの `lh = −p·C` と、比較の `model_stages: first`（速い経路だけ）は整合。C は `open_loop` に k0 = 1・遅延 0・極なしで C だけを評価させており、dsc1 の離散化も渡る。
- **ロックイン**: 段ごとの窓 [t_k + settle, t_k + dwell) は整数周期（励起が f_k を 1/(dwell − settle) の格子に寄せる）で、記録の n から窓を切るので rate_div 2（wire_refD）でも正しい。コヒーレンスは整数周期の部分区間で γ² = \|ΣR_m X_m*\|²/(Σ\|R_m\|²Σ\|X_m\|²)。部分区間 5 個では雑音だけのときの期待値が約 1/5 なので、閾値 0.9 は十分に保守的。
- **相互スペクトル**: multisine は 1 周期の矩形窓・重なりなし・周期の境界から（漏れなし）。チャープ・PRBS は S_xx の 10⁻⁶ 以上のビン。正しい。
- **余裕**: 位相を周波数に沿って unwrap し、最低周波数の点を (−360°, 0°]（モデルがあればモデルの枝）に置く。f_u は最初の下向きの \|L̂\| = 1、f_180 は最初の −180°。I-3・I-4 を除いて正しい。

**再現（refD_nominal ＋ wire_refD ＋ in_lock_stepped、1–30 kHz、コヒーレンス ≥ 0.9）**

| | 本レビュー | D93-5 |
|---|---|---|
| 雑音なし 0.02 V: 厳密な z 領域との差 | 0.010 %・3.2×10⁻¹⁰° | 1e−4 |
| 雑音なし: 判定（dsc1）との差 | 2.03 %・2.19° | 2.0 %・2.2° |
| 雑音あり 0.1 V: コヒーレンスの最小 | 0.988 | ≥ 0.988 |
| 雑音あり: 判定との差 | 4.33 %・2.26° | 4.3 %・2.3° |
| f_u・PM（モデル） | 2089 Hz・77.8°（2111・78.7°） | 同じ |
| 追加掃引 14–19.5 kHz: f_180・GM（厳密） | 16.79 kHz・12.92 dB（16.78・12.79） | 同じ |
| F_wg6: wire_F・p7_literal の判定との差 | 2.02 %・3.57°、1.07 %・3.56° | 2.0 %・3.6°、1.1 %・3.6° |

厳密な z 領域（`tests/integration/refD_support.loop_parts`、PID の時刻での多速度の L_d = −gC(BA + B′A′)）は WP6-I より前からある独立なモデルで、雑音なしで 10⁻⁴ の一致は推定量・記録の時刻・励起の配線がすべて正しいことを示す。判定との 2 % は判定（連続＋dsc1 の PID）側の近似。F_wg6 の 3.6° も 3 つの配線で同じ値で、判定側（共振の離散化など）のもの。

### 2.2 インターフェースと `SimController`（項目 2）

- `engage`/`disengage(loop_id=None)`、`acquire(..., injections_on=False)` は既定値付きで、既存の呼び出し（引数なし）の意味は変わらない。既存のサブクラスの `engage(self)` も ABC としては有効のまま（新しい引数で呼ぶときだけ失敗する）。`measure_loop_gain` は非抽象で基底は `NotImplementedError` — 保護リストの「抽象制御インターフェース」を壊していない。`status()["engaged"]` は「全ループが engaged」になり、全ループを一度に切り替える従来の使い方では同じ値。
- `inject_topology(spec.pipeline, spec.topology)` は葉の発見（実行時パラメータの名前の表）のためだけに使い、`_run` は元の spec を runner に渡す。runner も同じ関数で注入し、`_inject` は「topology が無い葉にだけ足す」ので二重にならない。refD・F_* が組み立てられなかった不具合の根本の修正として正しい。
- `measure_loop_gain` は波形のキーを構造パラメータとしてその 1 回のランにだけ入れ、記録のスロットを足し、記録が無ければ KeyError（D34）。期間の計算（t_start + 段数 × dwell、周期ものは settle を周期に切り上げ + n_periods）は `_spectral` の開始と同じ。

### 2.3 Moku の換算（項目 3）

PID の制御行列の中で u = 1·e + g·x、FRA は H = e/x を測る。e = −L u = −L(e + g x) → H = −Lg/(1 + L) → **L = −H/(g + H)**、S = u/(g x) = (g + H)/g、T = −H/g。コードと一致し、真の g で L̂ = L が 10⁻¹² で戻る（`p6_wp6i_x_checks.py 1`）。

g の誤りの効き方（refD のモデル、仮定 g = 1）:

| | \|L\| = 7 | 3 | 1 | 0.3 | 0.1 | f_u | PM | GM |
|---|---|---|---|---|---|---|---|---|
| 真 | | | | | | 2111 Hz | 78.70° | 12.79 dB |
| g_true = +1 dB | −15.8 % | +3.2 % | +8.8 % | +9.6 % | +13.4 % | 2301 | 71.58° | 12.03 |
| g_true = −1 dB | −28.1 % | −13.6 % | −9.5 % | −9.0 % | −11.7 % | 1908 | 86.68° | 13.57 |
| 経路の遅延 100 ns | −0.1 % | −0.1 % | −0.1 % | −0.1 % | +0.1 % | 2108 | 78.65° | 12.85 |

(1 + L̂)/L̂ = (g′/g)(1 + L)/L なので、低い周波数の大きな \|L\| ほど誤りが大きい（I-1）。dB/線形の取り違えでは g = 20 と 10 のような桁の誤りもありうる。FRA のモード名（I-2）は仕様 p30 の In/Out。

### 2.4 p7_literal（項目 4）

- 配線: pid0f は feedback off で出力 = オフセット (π/2)/k_pzt（速いループが保っていた X = +π/2）、x は sup0 の後に足し、pid0s は誤差 v.adc0 を直接読む。K_i = 2π f_hold/(K·k_TEC) は pid0f と同じ符号で、PZT と TEC はともに F0 の上で正の駆動で位相を増やす（F_common の `chk.lock_stability` の段も k_act = k_TEC/k_pzt > 0）ので負帰還になる。整合している。
- 応答 P̂ = v.adc0/x は 1/(1 + L_hold) を含むが、200 Hz 以上では L_hold ≈ 10⁻⁴ 以下で無視できる。
- 物理: TEC の保持ループは I ＋ 熱の遅れの 2 次系で、連続時間では常に安定だが（100 Hz の保持の遅れも無視できる）、格子の隅では強く振動的（I-5 の表）。論文の遅い段はストレッチャ（帯域不明）で、TEC では 1 Hz より上を保てない。0.65 s の掃引中は保持が実質働いていないので、ドリフトの大きい構成・長い掃引では動作点が縞の傾きから外れる。オーバーレイのヘッダはこれを記録しており、F_wg6_nominal では成り立つ（1.1 %）。

```
f_hold  τ_th   ζ      f_c [Hz]  PM    1/(2πf_c) [s]
0.003    1   3.642   0.0030   88.9    53.1
0.03    10   0.364   0.0192   39.7     8.3
0.3    100   0.036   0.0218    4.2     7.3
```

### 2.5 wire_F と制御器の上書き

`F_wg6_nominal + controller/dsc1 + ident/wire_F + ident/in_lock_stepped` は組み立てられ、走る（D93-4 の主張は正しい）。25 点の掃引では大きさ 6.9 %・位相 6.3°・GM 9.9 dB（判定 7.5 dB）とずれるが、原因は (a) 12 bit の DAC の量子化（量子化を切ると 1.5 %）と (b) 共振付近の粗い補間（追加掃引で 7.39–7.47 dB・19.28 kHz、判定 7.53・19.20）で、判定と時間ループの食い違いではない（I-3、I-6）。Moku と fpga では追加掃引で 6.25 dB・19.87 kHz（判定 6.01・19.81）。

## 3. B: D94（軸 1 のベンチ構成とプローブの上書き）

```
uv run python reviews/scripts/p6_bench_smoke.py [bench/dsc1_refA | bench/moku_delta_refA]   # 各約 1 分
uv run python reviews/scripts/p6_wp6x_smoke_settle.py [bench/dsc1_refA]                    # 各約 1 分
uv run python reviews/scripts/p6_wp6x_fc_delay.py
uv run python reviews/scripts/p6_wp6x_refA_demod_sign.py 50
uv run python reviews/scripts/p6_p4a_check.py
uv run python reviews/scripts/p6_wp6x_two_tone_bg.py
```

### 3.1 refA_fiber（項目 5）

**SQ を PM ファイバに**: 駆動は strain・dT、係数は k_ε L（ξ_eo 0.78、refA の k_ε と同じ）、k_T L（44 rad/(m·K)）、群遅延 κ_τ L（39 ps/(km·K)）と strain の n_g L/c（refD と同じ形）。LO ファイバの行も L_lo_m に比例。音響は refA の 2 極のモデル（T08 §2.5、γ = 1 で L に比例）に SQ を足し、熱はドリフトだけ（r_T 10⁻⁴ K/s、τ_th 300 s）。損失 10^(−0.05 L/km)（0.5 dB/km）× 0.98 × η_out は正しい。指摘は X-4（ヘッダの説明）と X-10 相当の注記（fiber_matched で 2 本を並べても音響は無相関、ただし両方とも lock 2 の中なので結果にほぼ効かない）。

**lock 2 の復調位相の吸収 2π f_sh Δτ**: 100 m では 2π f_sh Δτ = 3.06 rad ≈ π なので、符号を逆にしても誤差 ≈ 2π で同じようにロックし、smoke の表では符号を試せない。50 m（1.52 rad）で 3 通りを走らせた:

| 50 m、fiber_matched、20 ms | 指標 (i) の平均・std | lock 2 の誤差（平均・std） |
|---|---|---|
| コードの符号（+X） | +0.17・9.26 mrad | +0.000・0.009 rad |
| 吸収なし（0） | +244・23.1 mrad | +0.244・0.010 rad |
| 逆の符号（−X） | +18.6・658 mrad | −3.143・0.009 rad（隣の縞） |

コードの符号が正しい。吸収なしは傾きが cos(1.52) ≈ 0.05 に落ちるだけでなく 0.24 rad の静的なずれが出る（悪い構成として残す価値がある）。

**EOM → ストレッチャの肩代わり**: pid2s は EOM の駆動 ctrl.pid2.u を積分してストレッチャ（LO 上、正の係数）を動かすので、EOM の出力を 0 に戻す向きで正しい。|C₂k_fs| = k_eom を f_12 = 100 Hz で交差させる標準の形で、交差の位相差は 90°。20 kHz・Q 10 の共振は lock 2 の単位利得で (100/20000)·10 = 5 %。値は 915B 型の仮の値として妥当（X-5 の範囲の注記）。EOM だけでは 100 m で 3 seed 中 2 つが 19.7・44.1 mrad に悪化（肩代わりの必要性を確認）。

**LO の遅延の不整合と FC の補償**: θ に入るレーザーの位相は H_FC φ_L(t − τ_lo) − φ_L(t − τ_sq)。低い周波数で H_FC ≈ e^{−iωτ_FC}（τ_FC = 1/(πΔν_FC) = 96.5 ns）なので、スクイーズ光の経路が LO より長いと補償する。白色 FM の閉じた式 σ²(Δ) = ∫\|H_FC e^{−i2πfΔ} − 1\|² h0/f² df（Δ = τ_lo − τ_sq）で:

| L | 自由空間 Δ | 閉じた式（レーザーのみ） | 閉じた式 ＋ 1 m で合わせた残り 1.33 mrad | smoke |
|---|---|---|---|---|
| 1 m | +15.3 ns | 8.93 mrad | 9.03 | 9.03 |
| 10 m | −28.8 ns | 5.85 | 6.00 | 5.97 |
| 17.8 m（最適） | −66.9 ns | 4.84 | — | — |
| 100 m | −469.5 ns | 20.25 | 20.29 | 20.57 |

時間ループの値は閉じた式と 1.4 % 以内で合い、機構の説明は正しい（X-2 は最適の位置だけ）。ファイバ整合は Δ = 20.1 ns で一定（9.27 mrad）。

**smoke の表と過渡（X-1）**:

| L | 経路 | σ_total（作者の表） | σ_total（t ≥ 15 ms） | 最初の 15 ms の max\|v.dac2\| |
|---|---|---|---|---|
| 1 | fiber_matched | 9.37 ± 0.04 | 9.33 ± 0.07 | 0.068（最大 0.117） |
| 10 | fiber_matched | 9.48 ± 0.12 | 9.34 ± 0.08 | 0.150（0.226） |
| 100 | fiber_matched | 10.88 ± 1.16 | **9.41 ± 0.07** | 0.442（**0.500**） |
| 1 | free_space | 9.03 ± 0.02 | 9.09 ± 0.04 | 0.050 |
| 10 | free_space | 5.97 ± 0.07 | 5.88 ± 0.03 | 0.061 |
| 100 | free_space | 20.57 ± 0.06 | 20.56 ± 0.10 | 0.115 |
| 100 | fiber_matched、DSC1 | 12.33 ± 1.29 | **9.76 ± 0.10** | 0.697（0.863、±4 V） |

（σ_total の「t ≥ 15 ms」は時間ループの分を落ち着いた区間の std に置き換えた値。台帳の分はそのまま。）四分割すると 100 m・seed 3 の指標 (i) の std は 15.96 / 8.85 / 8.99 / 8.86 mrad で、最初の 12.5 ms だけが高い。

**DX1 の 400 倍**: h0 の比 2.55×10⁴ / (200/π) = 400.6、FC の項 πh0/Δν_FC = 0.0243 rad²（156 mrad）、20.1 ns の遅延の項 0.0101 rad²、自由空間 100 m の遅延の項 0.236 rad²（0.486 rad）はすべて正しい（分散で 400 倍、rms で 20 倍）。要求の書き方は X-3。

### 3.2 DSC1・Moku の上書き（項目 6）

- `bench/dsc1_refA`: rate_div = floor(20 MHz × 2 µs + 0.5) = 40、dl の τ = 8.3 µs − 1/f_sim − 39/(2 f_sim) = 7.275 µs。時間ループの全遅延 = τ + 1/f_sim（D55-1）+ 保持 (r − 1)/(2f_sim) = 8.3 µs で、`chk.lock{1,2}_stability` の delay_s 8.3 µs・dsc1・h 2 µs と一致する。ADC の 2 点平均は avg_samples 1 として τ_L の中（純遅延の形、D61-2 の規則で明示的にモデル化したものだけを差し引く）。D91-1 の「保持の遅延を τ_L の中に数え二重に数えない」とも一致。refC（4 MHz、rate_div 8）も同じ形で正しい。
- I ループの −180°: 90° + 360° f·8.3 µs + atan(f/500 kHz) − (dsc1 の進み ≈ 360° f·1 µs) = 180° から ≈ 32 kHz、20 kHz の単位利得で \|L\| ≈ 0.6（GM ≈ 4 dB）。テストの 0.5 < max\|L\| < 0.8・25–40 kHz と合う。
- Moku（refA）: 0.21 µs ＋ 復調器の LPF 500/250 kHz で −180° は 0.52/0.39 MHz（手計算で 46° + 44°、57° + 33°）。検査は 0.5/f_sim（25 ns）を足すが時間ループには無い（G-8 と同じ種類、無視できる）。
- 判定「DSC1 は Ref A に遅すぎる」は X-6 のとおり条件付きで確認。Moku は fpga とほぼ同じ（9.36/9.45/10.77、9.02/5.96/20.52）。

### 3.3 actuator_pathlength（項目 7）

`lit-svarc2023.md` §5.1・`fact-p6-thermo-optic-moku.md` §1 からの独立な導出:

- ストレッチャ・ひずみの効きの比 g = (ξnω)_r/(ξnω)_s = 1 + (Δν/ν0)(n_g/n + ω∂lnξ/∂ω)。R のモデル（ξ の分散なし）で f = n_g/n = 1.0166、R′ で 1.007（ω∂lnξ/∂ω ≈ −0.0096）。コードと一致。
- ストレッチャのロック（g δψ_a + c_r δT = 0）でスクイーズ光に残るのは (c_s − c_r/g)δT。c_r = g(k_T + 2πΔνR)L で −2πΔν L R δT（R の定義そのもの）。物理の c_r = k_T + Δω κ_τ = g k_T + Δω R との差は (g − 1)ΔωR で二次。TEC は環境と同じ係数なので 0。
- R′: c_r/g − c_s = Δω[R − (β + nα)(ω∂ξ/∂ω)/(cξ)] = Δω R′ で、ノートの R′ の定義と一致。implicit は R = κ_τ − f k_T/ω0 で `probe/freq_*` の (k_T + 2πΔνκ_τ)L に一次で戻る（44/39 で 2.56 ps/(km·K)、ヘッダの 2.6）。bad 版の範囲 −2.36…+5.87 も再計算で一致。
- ホップ: 監視器のリセットを u_2π = 2π/(g K_pzt) にするとプローブで 2π、スクイーズ光で 2π/g、偏り 2π(1/g − 1) ≈ −2πfΔν/ν0（Pittaluga の 2πΔν/ν0 × f）。`largerange/fringe_hop` は u_2pi_v を書かないので順序の衝突もない。
- `p6_p4a_check.py` part 4: 4 つの Δν × 3 モデルで全 assert が通り、残差 = −2πΔνLR を 10⁻⁶ で再現、TEC 0、ひずみ ≤ 1.4×10⁻⁶ rad/K（数値の丸め）。

### 3.4 two_tone（項目 8）

- 偏り: v = K₊cos(X̄ + E) + K₋cos(X̄ − E) = (K₊ + K₋)[cos X̄ cos E − r sin X̄ sin E] → cot X̄ = r tan E → X̄ = π/2 − atan(r tan E)。E の符号は規約 §4（+Ω のビームの位相 −Ω(t − τ)）どおり beam_assembler の 2π f_off τ_b。パワー P_± = P̄(1 ± r)²/(2(1 + r²)) で K_± ∝ 1 ± r、和の傾きは r = 0 で単一トーンと同じ、×1/√(1 + r²)。正しい。
- 格子: L_k − L_0 = k(N_b + n_guard)c/(f_rep n_g)（F_common）なので δk(L_k − L_0) = 2πδν k(N_b + n_guard)/f_rep = πmk。m 偶数で 2π の倍数（傾きの符号も保つ）。PRS_k・D_k・REF_k の長さはポートに共通の 1 つのパラメータなので、モデルの中では条件はこれだけ。実機ではこれらの分配ファイバをポート間で µm 級（7 GHz で 1 mrad ↔ 4.6 µm、偏りは r 倍）にそろえる必要がある — 記録しておくべき。
- ±δν の熱の行に κ_τ を使うのは正しい（E のドリフトは群遅延 2δω κ_τ L そのもので、R は平均 X̄ で打ち消す）。
- **開いている点（X-7）**: 本レビューで「解ける問題だが、今の結論は同偏光のプローブにしか使えない」ことを数で確認した。表（F_liu20、7 GHz、格子上、E_0 = 0）:

| ポート | L_k − L_0 | b_g による E のずれ（mod 2π） | cos E | 偏り r = 0.1 | 偏り r = 0.2 |
|---|---|---|---|---|---|
| 1 | 40.8 m | +2.443 | −0.766 | +0.084 | +0.166 |
| 2 | 81.7 m | −1.397 | +0.173 | +0.518 | +0.850 |
| 15 | 612.7 m | −1.053 | +0.495 | +0.174 | +0.337 |
| 31 | 1266 m | +0.338 | +0.943 | −0.035 | −0.070 |

（cos E < 0 は隣の零点でロックするだけで直交位相としては同じだが、\|cos E\| の小さいポートはループ利得が落ちる。）

### 3.5 refC_fiber と pulsed_wg_1port（項目 9）

- refC_fiber: SQF_dep・LOF_dep の係数（dT・strain の位相・群遅延）はすべて L_deployed_m に比例して 0、スプールは L_fiber − L_deployed = L。`env.mech.conduit` は 0 の adder で drv.SQF_dep.strain・drv.LOF_dep.strain を書き（読み手が書き手を要するため）、熱の conduit はドリフトを切る。η_link = 0.973 × 10^(−0.2 dB/km · L)（0.9294/10^(−0.02) = 0.973）。整合。
- pulsed_wg_1port: `L_k_m: !expr L_port_m`（D64）が F1 の長さ、L_lo1_m（LO1 の整合ファイバ）、τ_k、`tec_k_rad_per_v_k`（検査とプラント）、stray1・sbs_probe の L に伝わる。delay1 は port 0・L0_m = L_k で自分の長さを L_k にする（タイミングの静的な不整合 0）。port_k = 31 を読むのは demux/none で無効にした dmx・dmxsw だけ。整合。
- 両者は同じ方針（D64: 軸・式の置き換えで既存プリセットを変えずに長さを変え、依存する量はすべて式で追従）で一貫している。

## 4. D93-5 と D94 の結果が成り立つか

- **D93-5: 成り立つ。** 全数値を再現。I-1〜I-6 は使い方・記述の改善で、検証値は変わらない。
- **D94-1〜6（構成）: 成り立つ。** 物理（係数、損失、復調位相の符号、肩代わり、遅延の勘定、g・R・ホップ、two_tone の偏りと格子、refC/pulsed の上書き）は正しい。
- **D94 結果（Ref A ＋ファイバ）: 一部直す。**
  - 自由空間 9.0 / 6.0 / 20.6 mrad は成り立つ（閉じた式とも一致）。最良は 10 m 付近ではなく 17.8 m（X-2）。
  - ファイバ整合 9.4 / 9.5 / **10.9** → 定常は 9.3 / 9.3 / **9.4** mrad（X-1）。
  - DSC1 の差「0.3–1.5 mrad」→ 定常で 0.3–1.0 mrad（X-1、X-6）。
  - DX1 の 400 倍・FC の項 155 mrad は成り立つ。「LO の遅延整合 ≲ 1 ns」は FC があるときは「スクイーズ光を約 0.7 τ_FC 長く（±10 ns）」が正しく、主な要求は FC の再設計（X-3）。
- **D94-5 の P8 の「全ポート共通で指標 (ii) で打ち消す」**: 同偏光のプローブについては成り立つ。直交偏光（今の two_tone の唯一の形）では群複屈折で崩れる（X-7）。

## 5. テスト

`uv run pytest -m "not slow"` → **1673 passed, 8 skipped**（matplotlib・botorch が無いための skip、Ref D のブラインド手順、PID の Euler 形の既知の skip）, 288 deselected, 22 warnings（56 分）。WP6-I・WP6-X のテスト（`test_phase6_ident.py`、`test_p6_bench.py`、`test_excitation.py`、`test_sim_ident.py`、`test_loop_identification.py`）はすべて通る。X-1 の過渡は既定の seed のテストでは見えない（seed 2・3 で現れる）。
