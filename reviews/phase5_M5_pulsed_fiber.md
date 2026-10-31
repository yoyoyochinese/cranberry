# Phase 5 物理レビュー: M5 のパルス光源とファイバ（WP-B9・B10・B11、WP-B5 の `pulsed_bhd_kernel`）

- レビュー担当: 新しい文脈の物理レビュー担当（実装には関わっていない）。`blind/` は読んでいない。ponytail は使っていない。
- 基準: `dev`（worktree で `git merge --ff-only dev` した時点、`beb30fa`）。ソースコードは変更していない。追加したのは `reviews/` の下だけ。
- 読んだもの: `AGENTS.md`、`DECISIONS.md`（D07・D43・D44・D48 §7-6・D49・D50・D51-6・D52 G03・D55-4・D58-2・D59・D60）、`notes/conventions.md`（§2・§3・§4・§8・§9）、`PLAN.md` §4.3・§6.2–§6.4・§7.4、
  T03・T06・T07・T10（必要な節）・T17（全文）、`reviews/phase5_M1_*.md`（書式）。
  コード: `blocks/fib/{gnlse,quantum,scattering,mode_overlap,tables}.py`、`blocks/source/{pulsed_opo,opo_angle_map,opa_supermodes,temporal_modes}.py`、`blocks/fsr/comb_filter.py`、
  `blocks/spf/{spectral_filter,dispersion_comp}.py`、`blocks/los/lo_shape.py`、`blocks/chk/rep_rate_tail.py`、`blocks/det/pulsed.py`、B8 の `blocks/mod/probe_pulser.py`（迷光のゲート）、`core/fft.py`、`offline/source_output.py`。
  一次資料（Madsen 2022 SI、Liu 2026 SI）は T17 の引用の範囲で、値を T17 経由で照合した（新しい値は読んでいない）。
- VT-03・06・07・08・09・10・11・12・32・35・36・45 は、この worktree で全件通過した（74 件、32 s、`-m ""` で slow を含む）。
- 数値の確認に使ったスクリプトは `reviews/scripts/m5_*.py` に置いた（§7）。いずれも末尾に assert の自己検査がある。

---

## 1. 指摘の一覧

| ID | 重大度 | 場所 | 要旨 |
|---|---|---|---|
| M5-01 | **major** | `src/cranberry/blocks/fib/quantum.py:286-312`（区間ごとの台形の雑音注入 309 行）、D59-5 | Raman をオンにした線形化の共分散が、区間（`steps.n_save`）の数が足りないと**非物理**になる（V + iΩ/2 の最小固有値が負）。D59-5 の「T = 0 では区間を 16 以上」は特定の条件の値で、一般には足りない。1 ps の基本ソリトンを 10 L_D 伝搬させると、T = 0・16 区間で λ_min = −6.8×10⁻³（真空 1/2 の 1.4 %）、295 K でも 4 区間以下で −1.2×10⁻²。コードには物理性の検査も警告もない |
| M5-02 | minor | `src/cranberry/blocks/mod/probe_pulser.py:31-32, 98-110`（WP-B8 だが D59-7 で B10 の迷光と組になる部分） | パルスにした probe の迷光で、「モードあたりの光子数」にもゲートで一様に平均した因子をかけている。conv §2 の N = ∫n|ũ|²df はモードの重みで決まり、ゲートの重みではない。Borealis 型の ns のモード（25 ns のゲート）で、同じゲート因子 0.40 に対し、正しいモードの因子は 0.035–0.47（狭線幅の反射）と 10 倍以上ばらつく |
| M5-03 | minor | `src/cranberry/blocks/fib/scattering.py:233-257, 322-332`（`mode` の指定） | `stray_light_table`・`stray_photons` のモードは `gaussian`（t_fwhm）か `exponential`（τ_eff）のパラメトリックな形だけで、`pulsed_opo_source` の最適モード（τ_eff = 14.5 / 7.5 ns）と結び付ける手段がない。プリセットが 1 ps の既定値（dmx の fixture）のまま OPO に使われると、Rayleigh の光子数を 3.5–6.8×10³ 倍少なく見積もる（T17 の結論が消える）。`exponential` の Lorentz 形は τ_eff = 3.0/κ_c の実際のモードと |ũ(30.5 MHz)|²/τ_eff で 0.34 対 0.46 と違う |
| M5-04 | minor | `notes/pulsed-opo-source.md:244, 246, 333`、`tests/vt/test_vt_45_pulsed_opo.py:228-229`、`src/cranberry/blocks/los/lo_shape.py:19-21`、D50 | T17 の訂正（D60-5）は r₀ と Gaussian の LO だけを直したが、同じ表の R±（§5.2、V-P13 の目標）も誤り。収束値は Borealis −7.72 / +8.82 dB（T17 −7.66 / +8.75）、Liu −8.95 / +9.22 dB（T17 −8.99 / +9.26）。V-P13 のテストは −7.66 ± 0.1 dB で古い目標を持ち続けている。`lo_shape` の docstring と D50 は 0.951 のまま |
| M5-05 | minor | `src/cranberry/blocks/source/opa_supermodes.py:18-19, 109-121` | 超モードの格子: r₀ は収束するが、Schmidt 数 K（と全光子数）は格子のスペクトル幅に比例して増え続ける（1 ps・1 mm の VT-32 の fixture で N = 128 / 256 / 512 → K = 24 / 48 / 88）。位相整合の帯域（この fixture で約 40 THz）を格子が覆っていないことを知らせる検査がない（GNLSE の `grid_edge_monitor` に当たるものがない） |
| M5-06 | minor | `src/cranberry/blocks/chk/rep_rate_tail.py:71-73`、`src/cranberry/blocks/source/pulsed_opo.py:266-267`、`src/cranberry/blocks/source/temporal_modes.py:54-59` | 1/f_rep が光源の計算窓（`window_kappa`/κ_c）の外にあると、`tail_model: source` の裾は 0 になり、警告が出ない。`window_kappa = 3`、60 MHz の Borealis 型で、正しい裾 8.8×10⁻³（閾値 10⁻³ の 9 倍）を見逃す。`tail_into_next_slot` も同じ |
| M5-07 | minor（記録） | `src/cranberry/blocks/source/pulsed_opo.py:18-25, 106-138` | 箱関数への射影による不確定性関係の破れ（D60 の docstring）。δt = 1/(40κ_c) で C の最小のシンプレクティック固有値は 0.99979（2.1×10⁻⁴ の破れ）、1/80 で 0.99995。dB の値には効かない（0.002 dB）が、下流で熱的光子数 (ν − 1)/2 が −10⁻⁴ になり、GBS の標本化などで物理性を要求する処理が失敗しうる。ビン内の雑音を Van Loan の方法で厳密に積めば、破れは丸め誤差まで消せる |
| M5-08 | info | `src/cranberry/offline/source_output.py:64, 77`、`src/cranberry/blocks/source/pulsed_opo.py:257-258, 276-285, 361-366` | `pump_path_phase_rad` を読むコードがまだない（runner・path のどこにもない）。D60-1 の「½ψ で θ_sq に入る」は統合（WP-I5）で初めて効く。`opo_angle_map` と光源の静的な離調の二重計上も規約だけで防いでいる。`pulsed: off` の経路は Δ_p とポンプの共振を捨て、係数表の名前（`.probe`、軸 `x`）が違うので `angle_map_params` が KeyError になる |
| M5-09 | info | `src/cranberry/blocks/los/lo_shape.py:166-175`、`src/cranberry/blocks/fib/mode_overlap.py:149-169` | `projection`（η_t の真空をすでに含む）と `losses()` の `lo_shape` 要素（同じ η_t）の両方を持つ。統合で両方を使うと二重計上になる。どちらを使うかを統合の契約に書くこと |
| M5-10 | info | `notes/fiber-quantum-noise-scattering.md:214` | T07 §3.1 の「⟨T⟩_h = 8.7 fs、g_R′ = 4.1×10⁻¹⁸」は Blow–Wood の値として誤り。正しくは 2τ₁²τ₂/(τ₁² + τ₂²) = 8.122 fs、g_R′ = 3.80×10⁻¹⁸ /(W m)/(rad/s)。コードは正しく、T07 §3.4 の表（0.1 THz で 1.5×10⁻⁴）も正しい値と整合する（本文の数字だけの誤り） |
| M5-11 | info | `src/cranberry/blocks/det/pulsed.py:3-12` | `pulsed_bhd_kernel` は各パルスの電荷を t = jT のδ関数とみなす。ns のモード（τ_eff 14.5 ns、裾 ~20 ns）で f_rep が 45 MHz 付近になると、モード自身の長さが T と同程度になり、c_m を過小に見積もる。OPO のプリセットでは窓の中のモードの重み（`lo_temporal_shape` の `cw_digital` と同じ畳み込み）で c_m を作る変種が要る |

blocker はない。D60-5 の訂正（r₀ ≈ 1.04、Gaussian の LO 0.955）は**確認**（§4）。D59 の選択は 5 を除いて確認、D60 の選択はすべて確認（§5）。

---

## 2. M5-01（major）: Raman をオンにした線形化の共分散の物理性（D59-5）

### 2.1 何が起きているか

- `fiber_quantum_linear` は `steps.n_save` 個の区間ごとに V ← M V Mᵀ + N とし、Raman の雑音を台形で N = (h_seg/2)(M D(z_k) Mᵀ + D(z_{k+1})) と入れる（`quantum.py:309`）。
- M には Raman の利得（Stokes）と吸収（anti-Stokes）が入っている（B9 の接線、Im h̃_R）。対称順序の雑音 S_R = hν₀ g_R (n_th + ½) は、この利得・吸収とちょうど釣り合う最小の雑音（揺動散逸）なので、T = 0 では物理性の境界に余裕がない。台形の注入は区間の中の M(z_{k+1}, s) D(s) M(z_{k+1}, s)ᵀ の変化を 2 点でしか拾わず、利得・長さの積が大きい区間では雑音が足りなくなる。
- D59-5 は「T = 0 では区間を 16 以上」とする。これは VT の fixture（1 ps の基本ソリトン、1 L_D）の値で、区間あたりの Raman の利得・長さの積で決まる量を固定の数で書いている。

### 2.2 数値（`reviews/scripts/m5_fiber_raman_physicality.py`、1 ps 基本ソリトン P = 16.7 W、64 点、400 RK4IP ステップ、Kerr + β₂ + Raman + Raman 雑音）

λ_min(V + iΩ/2)（conv の単位、真空 1/2。物理 ⇔ ≥ 0）:

| L | n_save = 1 | 4 | 16 | 50 | 100 | Raman オフ（RK4IP の分） |
|---|---|---|---|---|---|---|
| 1 L_D, T = 0 | −2.9×10⁻³ | −3.8×10⁻⁴ | −9.0×10⁻⁶ | −6.5×10⁻⁷ | −1.3×10⁻⁷ | −1.1×10⁻⁹ |
| 3 L_D, T = 0 | −3.4×10⁻³ | −4.1×10⁻³ | −1.9×10⁻⁴ | −7.6×10⁻⁶ | −1.5×10⁻⁶ | −5.1×10⁻⁷ |
| 10 L_D, T = 0 | −1.2×10⁻² | −1.2×10⁻² | **−6.8×10⁻³** | −2.0×10⁻⁴ | −8.1×10⁻⁵ | −4.7×10⁻⁴ |
| 10 L_D, 295 K | −1.2×10⁻² | −1.2×10⁻² | +1.4×10⁻⁵ | +1.3×10⁻⁵ | — | — |

- 1 L_D（D59-5 の根拠の条件）では 16 区間で −9×10⁻⁶ と小さいが、10 L_D では 16 区間で真空の 1.4 % の破れが残る。squeezed の入力（r = 1.15）でも同じ値になる（最も破れるのは真空の部分）。
- 295 K でも 4 区間以下では破れる（熱の余裕は区間が粗いときの不足を補えない）。
- Raman オフの列は RK4IP の接線のシンプレクティック性の誤差（O(h⁴)）で、10 L_D・400 ステップでも −4.7×10⁻⁴。Raman の分はその 15 倍。
- VT-36 V3（物理性）は 295 K・1 L_D・4 区間だけを検査していて、T = 0 の Raman の物理性のテストはない。

実際に効く場面: Fw の 100 fs–1 ps の LO・probe（ピーク数十–100 W、63 m）では帯域内の Raman の利得・長さの積が 1 を超えうる。T07 §2.5 はこの LO の振幅雑音 → 位相雑音を線形化で扱う想定である。

### 2.3 修正案

1. `fiber_quantum_linear` の出力で λ_min(V + iΩ/2) を計算し、−tol（例 10⁻⁹·max|V|）を下回れば `ChkLog` に `fiber_quantum_unphysical` を出す（2N のエルミート行列の固有値 1 回、M 全体の計算より安い）。D34 のとおり止めない。
2. 雑音を RK4IP のステップごとに積む選択肢（Q ← M_step Q M_stepᵀ + h·D）を加えるか、区間の長さを「区間あたりの帯域内の Raman 利得・長さの積 g_R,band·P_pk·L_seg ≤ ε」から自動で決める。
3. D59-5 を「区間あたりの g_R P L_seg ≪ 1（T = 0 では 10⁻² 程度で 10⁻⁵ の破れ）」に書き換え、T = 0 の物理性を VT-36 に加える（1 L_D・16 区間と、10 L_D で警告が出ること）。

---

## 3. その他の指摘の詳細

### M5-02（minor）: パルスの probe の迷光で、モードあたりの光子数にゲートの平均をかけている

`gated_stray_photons` は分布散乱・反射のどちらでも、`clicks_per_gate` と `photons_per_mode` に同じ因子 f（ゲートの中の光子束の一様平均）をかける（`probe_pulser.py:101-108`）。クリック数（ゲートの積分）には正しいが、モードの光子数は conv §2 の N = ∫n|ũ|²df、時間領域では狭線幅の光（Rayleigh・反射、コヒーレンス時間 ≫ モード）で N = Φ_pk |∫u(t)·1_pulse(t) dt|²、広帯域で N = Φ_pk ∫|u|²·1_pulse dt になる。

`reviews/scripts/m5_stray_gate_vs_mode.py`（Borealis 型の最適モード、ゲート 25 ns をモードの重心に置く、10 ns の反射パルス）:

| パルスの開始（ゲート中心から） | ゲートの因子 f | モードの因子（狭線幅） | モードの因子（広帯域） |
|---|---|---|---|
| −12.5 ns | 0.400 | 0.035 | 0.190 |
| −7.5 ns | 0.400 | 0.431 | 0.828 |
| −2.5 ns | 0.400 | 0.471 | 0.787 |
| +2.5 ns | 0.400 | 0.090 | 0.169 |

ps のモード（Fw）でも、ゲート（1 ns）の中にパルスの縁があると、モードの位置の光子束（0 か 1）ではなく割合をかけることになる。分布散乱（W = 2n_gL/c ≫ τ_pr の台形）は時間変化が緩いので誤差は小さい。

修正案: モードの強度 |u(t)|² の重み（狭線幅の線は |∫u√Φ|²）で因子を作る。モード関数は `pulsed_opo_source` の最適モードまたは `mode` のパラメトリックな形から作り、ゲートの因子とは別の列にする。

### M5-03（minor）: 迷光の表のモードを光源の最適モードと結べない

`stray_photons` の `mode` は `gaussian`（強度 FWHM）か `exponential`（τ_eff、Lorentz 形）。式（|ũ(f)|² = τ_eff e^{−π τ_eff² f²}、τ_eff/(1 + (π τ_eff f)²)）は正しく、両方とも ∫|ũ|²df = 1 を満たす。しかし:
- OPO の最適モード（τ_eff = 3.0/κ_c）を渡す手段がない。`lo_temporal_shape`・`mode_overlap_table` が D60-2 の `source:` の入れ子で光源を受け取るのに対し、迷光の表は独立の数値を手で書く。dmx の fixture は 1 ps の Gaussian で、OPO のプリセット（未作成、WP-I5）で書き換え忘れると Rayleigh の光子数が 3.5–6.8×10³ 倍小さくなる（T17 の最大の結論が消える）。
- `exponential` の Lorentz 形は片側指数（τ_eff = 2/κ）の形なので、τ_eff = 14.5 ns を入れても帯域の形が違う（30.5 MHz で 0.34 対最適モードの 0.46、T17 §4）。GAWBS の低次のピーク（30.5、82.1 MHz）の重みが 25 % 程度ずれる。
修正案: `mode: {source: <temporal source spec>}` を加え、`temporal_modes.mode_spectrum` の |ũ(f)|² を使う（V-P14 と同じ）。

### M5-04（minor）: T17 の R± の表、V-P13 の目標、0.951 の残り

`reviews/scripts/m5_opo_dt_convergence.py`（δt = 1/(20, 40, 80, 160 κ_c) と Richardson 外挿。差の比は 3.5–4.4 で δt² の収束）:

| 量（x_ss = 2.5） | 外挿値 | T17 | 判定 |
|---|---|---|---|
| Borealis r₀（ρ = 1） | 1.0380 | 1.029（§5.2 本文） | D60-5 の 1.04 を確認 |
| Borealis R_− / R_+（ρ = 0.95） | −7.717 / +8.821 dB | −7.66 / +8.75 dB | T17 §5.2 の表も誤り |
| Liu r₀（ρ = 1） | 1.066 | （−8.99 dB から逆算 1.071） | — |
| Liu R_− / R_+（ρ = 0.99） | −8.948 / +9.221 dB | −8.99 / +9.26 dB | T17 §5.2 の表も誤り（逆向き） |

- 最適モードの R_− は単一モードの ρe^{−2r₀} + 1 − ρ と 10⁻³ dB で一致する（V-P13 の趣旨は正しい）。誤っているのは T17 の r₀ の値だけ。
- `test_v_p13_pure_loss` は −7.66 ± 0.1 dB を検査していて、−7.72 dB でも通る。目標を −7.717 dB（許容 0.01 dB）に直す。
- `lo_shape.py:19-21` の docstring、D50 の本文は 0.951 / 0.864 のまま。矩形の最良値は 0.866（T17 0.864、§4 の表）。

### M5-05（minor）: `opa_supermodes` の K が格子に依存する

`reviews/scripts/m5_opa_supermode_grid.py`（VT-32 の fixture: 1 mm、α = 1 /W、0.25 W、GVM 0.3 ps/mm、β₂,s = 100 fs²/mm）:

| τ_p | W, N, n_z | r₀ | K |
|---|---|---|---|
| 1 ps | 8 ps, 128, 20 | 0.49418 | 23.9 |
| 1 ps | 8 ps, 256, 20 | 0.49440 | 47.8 |
| 1 ps | 16 ps, 256, 20 | 0.49418 | 24.0 |
| 1 ps | 8 ps, 512, 20 | 0.49440 | 87.7 |
| 1 ps | 8 ps, 256, 40 | 0.49439 | 47.8 |
| 100 fs | 8 ps, 128 / 256 / 512 | 0.2447 / 0.2469 / 0.2470 | 4.9 / 9.8 / 19.2 |

- K はスペクトルの幅 N/W に比例する（差の方向の位相整合の帯域 Ω_PM = √(2π/(β₂,s L)) ≈ 2π × 40 THz を格子が覆っていない）。r₀・先頭の超モードは収束している（先頭モードのスペクトルの端のエネルギーは 10⁻⁸ 以下でも K は増え続けるので、先頭モードの端の検査では不十分）。n_z の 2 倍は影響しない（中点の行列指数でシンプレクティック、VT-32 V2 を確認）。
- BHD で先頭の超モードを測る用途（r₀、V_−）には問題ない。GBS の純度（1/K）と全光子数 Σ sinh² r_n は格子で切られた値になる。
修正案: NΔΩ/2 ≥ 数 × Ω_PM（解析式で計算できる）または ⟨b_j†b_j⟩ の端の割合を検査する CHK を加え、K を報告するときに警告する。T03 §5.2 の「100 fs で K ≈ 320（±24 THz で 640）」と同じ性質で、ブロックの docstring の「格子は使う側の責任」に検査を足す。

### M5-06（minor）: 計算窓の外の裾

`tail_energy` は窓の外を 0 とする（`temporal_modes.py:54-59`）。`reviews/scripts/m5_rep_rate_window.py`: Borealis 型 x_ss = 2.4、f_rep = 60 MHz で、`window_kappa = 12` なら裾 9.2×10⁻³（解析式 8.8×10⁻³）で警告、`window_kappa = 3` では裾 0・`ok = True`。最適モード自体も窓で切られて変わる。修正案: 1/f_rep > 窓の終わりなら、解析式 8.8e^{−2κ_c T} で外挿するか、窓が短いことを警告する。

### M5-07（minor・記録）: パルス OPO の離散化の物理性

δt = 1/(40κ_c) で ρ = 1 の固有値の積の誤差 4.2×10⁻⁴（V-P4 の許容 1.3×10⁻³ の内）、最小のシンプレクティック固有値 0.99979（Borealis・Liu とも）。1/80 で 0.99995、1/4 の収束。T17 V-P1 の CW の解析値との一致、Richardson の比も確認した。破れの原因はビンの中の雑音を箱関数の成分だけで駆動していること（ビン内の高い周波数の成分を捨てている）で、各ビンの (X(t_{k+1}), ∫X dt) の雑音共分散を Van Loan の拡大行列の指数で厳密に計算し、独立な列として加えれば、C は丸め誤差まで物理的になる（計算量はほぼ同じ）。dB の値は 0.002 dB しか変わらないので、下流が物理性を要求するときだけ効く。

### M5-08（info）: `pump_path_phase_rad` の消費者と `pulsed: off` の経路

- `reviews/scripts/m5_opo_frame_detuning.py`: Δ_p/κ_c = 0.05、0.2、0.5 で、完全な b(t) で測った角度 = `cov` の角度 + ½·`pump_path_phase_rad` が 10⁻¹⁶ で成り立つ（D60-1 の枠は正しい）。
- ただし runner・path にこのフィールドを読むコードがない。統合で ½ψ を θ_sq に入れる処理と、そのテスト（静的な Δ_p を光源に入れた場合と、ポンプ経路の区間に入れた場合が同じ角度になる）が要る。
- `opo_angle_map` が区間を駆動し、かつ光源に静的な離調があるときの二重計上は、規約（docstring）だけで防いでいる。両方がオンなら警告する検査を統合で入れる。
- `pulsed: off` では `_cw_block` が Δ_p・κ_p を捨て、係数表は `<seg>.detuning.{squeezed,probe}`（軸 `x`）になる。`angle_map_params` は `.pump` と軸 `x_ss` を読むので KeyError になる。CW の極限（VT-03）でも同じ API で使えるように名前を揃えるか、明示的に未対応のエラーを出す。

### M5-09（info）: `projection` と `losses()` の二重計上の危険

`lo_temporal_shape.projection` は (V_−, V_+, θ_off) に η_t（とビンの外の真空）をすでに含み、同時に η_t を `lo_shape` の損失要素として登録する。`mode_overlap_table` の表も V_± に η_t を含み、η_t の列を別に持つ。どちらも docstring に注記があるが、統合（WP-I5）で「表の V_± を使うなら `lo_shape` の要素は無効にする」を契約として書くこと。

### M5-10（info）: T07 の ⟨T⟩_h

`reviews/scripts/m5_scattering_checks.py`: 直接の数値積分 ∫t h_R dt = 8.1220 fs（コード 8.1220 fs）。T07 の 8.7 fs・4.1×10⁻¹⁸ は本文の誤記で、T07 §3.4 の表の n_R(0.1 THz) = 1.5×10⁻⁴ /(W m)（コード 1.468×10⁻⁴）は正しい値から出ている。

### M5-11（info）: `pulsed_bhd_kernel` のδパルスの仮定

c_m = ∫_0^W h(t + mT)dt / ∫_0^W h(t)dt はパルスの電荷が t = jT に集中しているとした式。OPO のモードは重心が立ち上がりから 6.9 ns、裾が約 20 ns あるので、f_rep ≈ 45 MHz（T = 22 ns）では c_1 を過小に見積もる。Borealis の 6 MHz では無視できる。OPO のプリセットで高い繰り返しを試すときは、窓とモードの畳み込みで c_m を作る変種（`cw_digital` と同じ考え方）が要る。

---

## 4. T17 の訂正（D60-5）の判定

| 項目 | T17 | D60-5 | このレビューの収束値 | 判定 |
|---|---|---|---|---|
| Borealis r₀（x_ss = 2.5、ρ = 1） | 1.029 | ≈ 1.04 | 1.0380（δt² の Richardson 外挿、1/40 で 1.0378） | **確認**。T17 自身の §2.6 の表の補間（1.044）とも 0.6 % 以内 |
| Gaussian の LO の最良の重なり | 0.951 | 0.955 | 0.9554（x_ss = 2.5）、0.9552（2.4）、Liu 0.9547。`lo_temporal_shape` のブロックでも 0.9552 | **確認** |
| 矩形の窓 | 0.864 | — | 0.866（幅 2.13/κ_c） | 0.2 % の差。D50 の 0.864 は 0.866 に |
| 片側指数 | 0.841–0.845 | — | 0.840–0.845 | 確認 |
| R± の表（§5.2・V-P13） | −7.66/+8.75、−8.99/+9.26 dB | 記載なし | −7.72/+8.82、−8.95/+9.22 dB | **追加の訂正が必要**（M5-04） |
| K（x_ss = 2.4） | 1.072（§2.4 の表） | — | 1.0713 | 確認 |
| x_ss ↔ r₀ の表（§2.6） | 2.401 ↔ 1.0 など | — | 2.4 → 0.9950 | 表は 0.5 % 程度の精度（テストの 1 % の許容は妥当） |

T17 のほかの数値（V-P1 の CW の解析値、V-P5 の K = 1.10、V-P6 の 1/K = 0.975、V-P9 の離調の係数、V-P10・V-P11）は VT-45 と独立の計算で一致した。

## 5. D59・D60 の選択の判定

| 決定 | 判定 | 根拠 |
|---|---|---|
| D59-1 接線写像で V ← MVMᵀ + N | 確認 | 接線は RK4IP の段に沿った厳密な変分（有限差分と一致、VT-36）。複素線形でない δ|A|² = 2Re(A*δA) を正しく扱っている。損失の雑音 (1 − η)(𝟙 + MMᵀ)/(2(1 + η)) は、一様な損失では M = √η S（S はシンプレクティック）なので、Kerr と交互に働く区間でも物理的であることを証明できた（(𝟙 + iΩ) + η(SSᵀ + iΩ) ≥ 0） |
| D59-2 自己急峻化では光子数の直交位相 | 確認 | 自己急峻化つきの GNLSE は ψ̃ = Ã/√(1 + Ω/ω₀) でハミルトン形（H = (γ/2)∫|A|⁴）になる。雑音 D は場の基底で作って `tq` で変換、真空・損失は光子の基底で一様。VT-36 の 10⁻⁹ を確認。光子数 ∫|Ã|²/(ω₀ + Ω) の保存も Raman・自己急峻化つきで 2×10⁻⁹（`m5_gnlse_checks.py`） |
| D59-3 線スペクトルは周期的でない遅延 | 確認 | C(τ) は非負のスペクトルの Fourier 変換なので、Toeplitz 行列は Bochner の定理で半正定値。周期化すると格子の周波数にない線で正定値が崩れる。Raman（広帯域）は格子の周波数で標本化した循環行列で、固有値 = S_R(Ω_k) ≥ 0 |
| D59-4 減衰は写像の中 | 確認 | `losses()` は空。プリセット側の二重登録の禁止は統合の検査で |
| D59-5 T = 0 の Raman は区間 16 以上 | **否定（一般には不足）** | M5-01 |
| D59-6 TW の適用範囲 | 確認（再計算はしていない） | 対称順序の初期値・損失の増分 αh/2 の式は正しい。自己急峻化つきでは TW は場の基底で雑音を入れる（線形化と O(Ω/ω₀) で違う、検証用なので許容） |
| D59-7 迷光は B10、B8 はゲートだけ | 確認（分担）。ゲートの因子の使い方は M5-02 | |
| D59-8 偏光解消の GAWBS は未実装 | 確認（記録として妥当） | R0,m の式は 7 番目のピーク 323.3 MHz が測定 321.75 MHz と 0.5 % で一致、σ² = 1.32×10⁻¹⁰ rad²/m（T07 1.33×10⁻¹⁰）。「2γ_OM」の因子 2 と `mode_sum_factor` は原典未確認のまま（UNVERIFIED の印を残す）。`mode_sum_factor` は全モードの γ にかかるので、プリセットで R0,m の全 12 本に同じ γ を入れると「R0,7 の 1–5 倍」（T07）を 12 倍に数える。プリセットでは R0,7 だけに γ を置き、合計を `mode_sum_factor` で表すこと（VT の fixture はそうなっている） |
| D60-1 静的なポンプの位相を `pump_path_phase_rad` へ | 確認 | 10⁻¹⁶ で整合（M5-08）。消費者はまだない |
| D60-2 入れ子の `source:` | 確認 | 迷光の表にも広げる（M5-03） |
| D60-3 `mode_overlap_table` の 2 軸 | 確認 | δτ_k（`dt_s`）と δL_k（`dL_m`、φ₂ = β₂δL、φ₃ = β₃δL）が独立。`dispersion_comp` の残差 δτ + delay、φ₂ + GDD とも D51-6 のとおり。U†u_LO で LO を光源の枠に戻す扱い、Gaussian の閉形式 κ = (1 − iχ)^{−1/2} exp[−Δt²(1 + iχ)/(4T₀²(1 + χ²))] を導出し直して一致 |
| D60-4 `pulsed_bhd_kernel` の一本化 | 確認 | c_m の式（tail(mT) − tail(mT + W)）、一次応答で c_m = e^{−2πBmT}、ρ_adj = Σc_m c_{m+1}/Σc_m²、C_pulse = η_q N_LO/σ_el² は正しい。ns のモードでの制限は M5-11 |
| D60-5 T17 の訂正 | 確認（追加あり） | §4 |
| docstring の O((κδt)²) の破れ | 確認（大きさ 2×10⁻⁴、δt² で減る） | M5-07 |

---

## 6. 確認して正しかったもの（PLAN §7.4 の点検表）

1. **式とノート**
   - GNLSE: D̂(Ω) = −α/2 + iΣβ_kΩ^k/k!（Agrawal の ∂_T → −iΩ で導出し直した）、自己急峻化 (1 + Ω/ω₀)、Raman の畳み込みを h̃_R(Ω) の解析式で（conv §8、T06 §3.4）。Blow–Wood の h̃_R = [(τ₁² + τ₂²)/(τ₁²τ₂²)]/[(1/τ₂ − iΩ)² + 1/τ₁²]、Lin–Agrawal の h̃_b = (1 − 2ix)/(1 − ix)²、h̃(0) = 1、h_b の 1 次のモーメント 0 を導出で確認。g_R のピーク 0.678 /(W km)、13.08 THz（T07 0.68、13.1）。
   - RK4IP（Hult 2007）の段、接線の段、step doubling の誤差の制御。Raman を切ったときの γ の合計の保存（D49）。
   - CW の近道 A(L) = √P e^{−αL/2} e^{iγ_eff P L_eff}、側帯波の回転 θ_disp = LΣ_{偶数}β_kΩ^k/k!（奇数次は遅延）、量子の CW の行列 K = [[−α/2, −β_e], [β_e + 2γ_eff P, −α/2]] の符号（キャリアの枠で Kerr は振幅 → 位相、分散は â → âe^{iβ_eL}）。Im h̃_R を CW で捨てても、位相によらない雑音は 1 次で正しい（Stokes の G·½ + gL(n_th + ½) と anti-Stokes の平均 = ½ + gL(n_th + ½)）。
   - Raman の雑音: S_R = hν₀g_R(n_th + ½) の P L S_R/hν₀ が Stokes・anti-Stokes の光子密度の平均に 10⁻¹² で一致。D48 §7-6: T → 0 で熱の部分は 0、Stokes g_R P L は残る。Ω → 0 の極限 g_R′k_BT/ħ、詳細釣り合い aS/S = e^{−ħΩ/kT}（10⁻⁶）。
   - GAWBS: ピークの式、片側の S_φ = 2L S_G、σ_m² = 2πL hν₀ γ (n + ½)Δf、C(τ) の Fourier 対、N_G = Σ(Pσ_m²/2hν₀)(|ũ(f_m)|² + |ũ(−f_m)|²) を導出し直して一致。
   - Brillouin: b_B = (g_B/A_eff)hν₀(2n_th + 1)(πΔν_B/2)、Stokes・anti-Stokes の分け方、Smith の閾値、L_bs = ∫e^{−2αz}dz。
   - 後方散乱の共分散 ⟨â_n†â_m⟩ = ΔtG(T_n − T_m) と [[Re, Im], [−Im, Re]] のブロック、射影が ∫n|ũ|²df（conv §2）になることを導出。
   - パルス OPO: T17 §2.1 の M(t)（Δ_s = ω_cav − ω₀、conv §4 の e^{−iω₀t} の枠で −iΔ_s â）、出力 √(2κ_T)X − X_in、拡大行列の指数のブロック（∫e^{Ms}ds と ∫(δt − u)e^{Mu}du）、初期のキャビティ真空の列、ポンプのキャビティ（定常 κ_p/(κ_p + iΔ_p)、位相 −arctan(Δ_p/κ_p)）、LO のモードの規約（u = Σ(w_p − iw_q)f_n で p_LO = w·X、e^{−iφ/2} で θ_sq,0 = φ/2 + π/2）。
   - `opa_supermodes`: (b, b†) の生成子 [[iD, G], [G*, −iD]]、G_jk = g(ΔΩ/2π)Ã_p(Ω_j + Ω_k)（conv §8 で A_p b† → ∫Ã_p(Ω + Ω′)b̃(Ω′)†dΩ′/2π を確認）、ポンプの伝搬 e^{i(Δβ₁Ω + β₂,pΩ²/2)z}、DFT のモードの規格化。
   - `spectral_filter`: A_mn = ∫f_m*(h∗f_n)、実形式 R、RRᵀ の固有値 = AA† の固有値 ≤ 1 なので cov′ = RcovRᵀ + (𝟙 − RRᵀ)/2 は物理的。Lorentz の t = γ/(γ − iΩ) が因果的（conv §8 の符号）。
   - `fsr_comb_filter`: n_FSR = η_ins Σ[T_f(+m) + T_f(−m)] sinh² r_m（T17 §2.7）、MZI の式（Liu Methods Eq. 2）。
   - `lo_temporal_shape` の `cw_digital`: u_eff(s) = Σ_j k_j w(s + jδt)（因果的な応答で後のサンプルに出る向き）。η_t の正確なノルム（ビン平均の過大評価を避ける）。
   - 迷光のモードの重み: Gaussian の |ũ|² = τ_eff e^{−πτ_eff²f²}（τ_eff = τ√(π/ln 2)）と片側指数の Lorentz 形を導出し直して一致。VT-36 の 435 光子/モード（T17 §5.3）。
2. **符号**: Raman の SSFS は赤（重心 −0.066 THz、5 L_D）、自己急峻化でピークは遅れる（+39 fs）・青い側に裾（歪度 +0.044）、Δt > 0 は遅いパルス（e^{+iΩΔt}）、離調 d > 0 で θ_sq が減る（T02 V-B7 と同じ）、長いポンプで −0.795 d → CW の −0.8 d、長さの変化の 1 次の写像は d = 0.1 まで 6×10⁻⁴ で正確（`m5_opo_frame_detuning.py`）。
3. **2 倍の因子**: 2 側と片側の PSD（S_R、S_G は 2 側、S_φ の表は片側）、対称順序の ½、Brillouin の (2n + 1)、x = √(P/P_th) に当たる x_ss、κ_c = πΔν_c（HWHM）。
4. **単位**: g_R [1/(W m)]、S_R [rad² s/m]、C の格子化（白色の極限で D_nn = 2|A|²S/hν₀ が CW の D と一致）、ρ = κ_T/κ_c、κ_p = πΔν_p。
5. **FFT の符号と規格化**: 全ブロックが `core.fft`（conv §8）を通る。GNLSE の内部は t0 = 0、利用者に返すスペクトルは t0 = T[0] で絶対位相。`stationary_noise_cov` の密度の項、`ifft_agrawal(S_R)` の (1/2π)∫S e^{−iΩτ}dΩ を確認。
6. **数値の収束**: GNLSE の刻み半減で誤差の比 11.9 → 14.6 → 15.3（全項オン、N = 2 ソリトン）。パルス OPO の δt は 1/40 → 1/160 で δt² の収束（R_− の 1/40 と外挿値の差 0.0013 dB）。超モードは r₀ が収束、K は未収束（M5-05）。線形化の区間数は M5-01。
7. **トグルでオフ ⇒ 恒等**: 量子層の全トグルオフ（VT-36）、`spectral_filter`・`dispersion_comp`・`lo_shape`（`shape_mismatch`）・`opo_angle_map`（`map`）・`mode_overlap_table`（`dispersion`）、`pulsed: off` → CW の OPO とビット単位で一致（VT-03）、`multimode: off` → §7-1 の単一モード（VT-03）。
8. **損失・帯域の登録**: FSR フィルタの通過帯、スペクトルフィルタの平らな挿入損失（形の部分は `apply` で、登録しない）、分散補償の挿入損失、`lo_shape`（BHD のみ）。ファイバの減衰は量子写像の中（`losses()` は空、D59-4）。二重計上の危険は M5-09。
9. **値の出典の印**: Raman・Brillouin・GAWBS の定数に UNVERIFIED、fixture の Borealis・Liu の値に頁付きの出典。
10. **保護リスト**: GNLSE の項のトグル（`kerr`・`beta2`・`beta3`・`beta_higher`・`raman`・`self_steepening`・`loss`・`xpm_external`・`cw_shortcut`、D58-2）と量子層の雑音のトグル（7 個）、OPO のトグル（`pulsed`・`pump_resonant`・`detuning`・`include_initial_vacuum`）はすべて残っている。

---

## 7. スクリプト（`reviews/scripts/`、`uv run python reviews/scripts/<name>.py`）

| スクリプト | 内容 | 指摘 |
|---|---|---|
| `m5_opo_dt_convergence.py` | パルス OPO の δt の収束（r₀、K、R_−、Gaussian・矩形・片側指数の LO）、Richardson、最小のシンプレクティック固有値（数分かかる） | M5-04、M5-07、§4 |
| `m5_lo_shape_block.py` | `lo_temporal_shape` のブロックの当てはめの重なり | §4 |
| `m5_opo_frame_detuning.py` | D60-1 の枠（½ψ）、長いポンプの離調の回転、1 次の写像の精度 | M5-08、§6 |
| `m5_fiber_raman_physicality.py` | Raman つきの線形化の物理性と区間数、長さの依存 | M5-01 |
| `m5_gnlse_checks.py` | 刻み半減（全項）、SSFS・自己急峻化の符号、光子数の保存 | §6 |
| `m5_scattering_checks.py` | Blow–Wood の 1 次のモーメント、g_R のピーク、D48 §7-6、詳細釣り合い、GAWBS のピークと σ²、S_R と光子密度の整合 | M5-10、§5、§6 |
| `m5_stray_gate_vs_mode.py` | ゲートの因子とモードの因子の比較（反射の 10 ns パルス、ns のモード） | M5-02 |
| `m5_rep_rate_window.py` | 窓の外の裾で警告が消えること | M5-06 |
| `m5_opa_supermode_grid.py` | 超モードの格子・n_z の収束 | M5-05 |
