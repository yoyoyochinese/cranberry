# Phase 1 効果サーベイ（survey）

- Phase: 1（サーベイ、1 エージェント）
- 日付: 2026-10-28
- 入力: `docs/squeezing_sim_prompt.md`（§4 Phase 1）、`DECISIONS.md`（D01–D23、prompt と食い違う場合はこちらが優先）、`AGENTS.md`、
  `refs/` の Ref A・B・C、McKenzie 2005、`refs/manuals/`（Thorlabs DSC1、Thorlabs PDB4xx、Moku:Delta 仕様書）
- 目的: squeezing の量、または系のいずれかのロックに影響しうる効果をすべて挙げ、deep-dive の分担（§4）を提案する。
  ここでの数値は「当たり」をつけるためのもので、根拠にはしない。根拠は各 deep-dive ノートで一次資料から取る。

## 0. 凡例と記号

### 0.1 構成の略号

| 略号 | 構成 | 出典 |
|---|---|---|
| A | Hirota et al.、PPLN 導波路 OPA、phase-detection OPA によるロック（＋比較用の conventional tap） | Ref A |
| B | Takeno et al. 2007、PPKTP bow-tie OPO、同軸 probe、dither | Ref B |
| C | Verclas et al. 2026、モノリシック ppKTP OPO、squeezed 光と LO を別々の 1 km ファイバ、量子雑音ロック | Ref C（補助: McKenzie 2005） |
| D | 現行系: Cobrite DX1（CW、1548.8 nm）→ PM ファイバ MZI（片腕に 915B）→ PDB450C → DSC1 の DC side-fringe lock | D02 |
| F | 将来の構成: パルス光源 → 単一パス導波路 OPA → 時間で 32 ポートへ振り分け → ポートごとの遅延ファイバ L_k → GBS。LO_k は別ファイバ。アクチュエータは probe と squeezed 光が一緒に通る区間 | D07–D13 |

### 0.2 「段」列の略号（prompt §5.1）

| 略号 | 意味 |
|---|---|
| G | オフライン段の GNLSE の**項**（オン・オフは方程式の項として切り替える） |
| Q | オフライン段の量子雑音層（GNLSE の雑音項、または線形化した揺らぎ） |
| O | オフライン段の GNLSE 以外のモデル・表（OPA/OPO の伝達、モード重なり、損失収支などの静的パラメータ） |
| T | 時間ループ段のブロック（雑音源、センサ、誤差信号、制御器、アクチュエータ、遅延バッファ） |
| M | 時間ループ段の測定模擬（ESA、出力処理） |

### 0.3 値の印

- **[Ref A]** など: そのリファレンスの本文から読んだ値（節を併記）。
- **[manual]**: `refs/manuals/` のメーカー資料から読んだ値。
- **[AE]**: agent estimate（このサーベイ担当による概算）。**unverified**。deep-dive で一次資料により確認すること。
- **[UNVERIFIED]**: 文献値・文献情報を記憶から書いたもので、原典で確認していない。

### 0.4 記号（symbols, English）

- `L_k`: fiber length of port k (D08), `L_k ≈ k·c·T_rep/n_g`, k = 0…31.
- `T_rep`, `f_rep = 1/T_rep`: pulse repetition period / rate (Open, D07).
- `τ_p`: pulse duration (swept 100 fs – 10 ps, D07). `T0`: 1/e half-width used in `L_D`.
- `n`, `n_g`: phase / group index of the fiber. `c`: speed of light in vacuum.
- `λ`: vacuum wavelength. `ν = c/λ`: optical frequency. `h`: Planck constant.
- `β2`, `β3`: group-velocity dispersion and third-order dispersion. `L_D = T0²/|β2|`: dispersion length.
- `γ`: fiber nonlinear coefficient. `P_peak`, `P_avg`: peak / average power. `L_eff`: effective length.
- `φ_NL = γ·P_peak·L_eff`: self-phase-modulation phase.
- `ΔL`: arm length difference of the Ref D MZI. `τ = n_g·ΔL/c`: corresponding delay.
- `δL_k`: residual length mismatch between squeezed fiber k and LO fiber k (D09).
- `S_ν(f)`: one-sided PSD of laser frequency noise [Hz²/Hz]. `Δν`: Lorentzian (white-FM) linewidth, `S_ν = Δν/π`.
- `S_Δφ(f)`: one-sided PSD of the phase difference between two delayed copies [rad²/Hz].
- `θ`, `θ_rms`: phase error between LO and the squeezed (or anti-squeezed) quadrature; RMS value.
- `R±`: anti-squeezed (+) / squeezed (−) noise variance normalized to shot noise. `R'±`: same, including `θ`.
- `L` (Ref A context): total effective optical loss. `α`: SHG efficiency [1/W]. `P`: pump power.
- `x`, `P_th`, `T`, `L_ic`, `l`, `γ_c`, `Ω`, `η`, `ξ`, `ζ`, `ρ` (Ref B context): normalized pump, threshold, output-coupler transmittance, intracavity loss (`L` in the paper), round-trip length, cavity decay rate, normalized frequency, PD QE, homodyne visibility, propagation efficiency, escape efficiency.
- `PER`: polarization extinction ratio. `B`: modal birefringence of PM fiber.
- `ER`: extinction ratio of the time-demultiplexing switch.
- `K_D`: error-signal slope [V/rad]. `G_TIA`: transimpedance gain [V/A]. `ℜ`: PD responsivity [A/W]. `V`: fringe visibility.
- `ε`: power/responsivity imbalance of the balanced detector or 2×2 coupler.

## 1. 構成ごとの要点（リファレンスから読んだ値）

- **Ref A**（§3–4）: 波長 1545.32 nm / 772.66 nm。導波路長 45 mm、伝搬損失 約 0.1 dB/cm → 640 mW で導波路の実効損失 2%。
  光学部品 1%、光電変換 1%、回路雑音の等価損失 <0.2%（3 MHz、LO 15.4 mW でクリアランス 約 28 dB）、モード不整合 約 4%（残差から推定）。合計 L = 8%、θ = 9 mrad、α = 906 %/W。
  640 mW は「squeezer OPA 直後のポンプ出力の推定値」（§3.2）。probe は AOM 2 台で 1 MHz シフト、probe と LO の各経路に EOM。LO 経路に ASE 抑圧用のフィルタキャビティ。
  BHD の PD 帯域 35 MHz、QE 約 99%。conventional 方式の probe は数 µW、新方式では数 mW まで上げられる（§2）。
  検算 [AE]: L = 8%、α = 9.06 /W、P = 0.64 W、θ = 9 mrad で R'− = −10.15 dB、R+ = +20.56 dB（prompt §7-1 と整合）。
- **Ref B**（§2–4）: T = 0.123、往復 約 500 mm、PPKTP 10 mm、L = 0.0038（100 mW）、L(x) = 0.00249 + 0.00222x、G = 18.7、ζ = 0.99、θ = 1.5°、ξ = 0.988、η = 0.998、
  回路雑音 −21.7 dB（LO 5.5 mW）、P_th = 180 ± 2 mW。PZT 共振 14 kHz（PZT1）・22 kHz（PZT3）、変調 295 kHz（PZT2 の共振でもある）、レーザー雑音 89・177 kHz。
  ESA: 1 MHz、RBW 30 kHz、VBW 300 Hz。
- **Ref C**（§2–4）: モノリシック ppKTP 2.6 mm、ポンプ 約 180 mW、線幅は図 1 の設計値 2.37 GHz、Fig. 5 の当てはめ 1.75 GHz。
  SMF 1 km（0.2 dB/km）を squeezed 光用と LO 用に 2 本、うち 約 180 m は敷設、残りは巻いたまま。別々のジャケットで同じ電線管に並べて置く。
  ファイバの結合・伝搬・取り出しで 8%、ファイバ前で 15%（BHD の QE 97% を含む）、当てはめで損失 21%、位相雑音 約 0.1 rad。
  ファイバ EOM（exail MPZ-LN-10、Uπ = 4 V、200 mVpp → ±4.5°）で LO に 1 MHz の位相変調。包絡線検出器（ログアンプ、8–300 MHz）→ 1 MHz で復調。
  3 段の PID（ファイバ EOM: kHz、ピエゾミラー: その下、ファイバ位相シフタ: Hz 以下）。**制御器は Moku:Pro**（本文）。ロック時間は遅いアクチュエータのレンジで数分。
  BHD の干渉ディップ 0.96 GHz（BHD_A）、0.83 GHz（BHD_B）。暗雑音クリアランス 11 / 7.4 / 5 / 5.7 dB（30 MHz / 0.25 / 0.69 / 1.21 GHz）。
- **McKenzie 2005**: 量子雑音ロックの誤差信号 ∝ J0(θ1)J1(θ1) sin2θ0 (V1 − V2)·Δω（Eq. 7）。センサ雑音は分散の分散（kurtosis）で決まり、
  ロック安定度は (1/Δω)^{1/4} で弱く改善（Eq. 13, 14）。squeezed 側のロック点のほうが安定。損失・検出効率で悪化。
- **Ref D**（D02 ＋ manuals）: 入力 3.33 mW。PDB450C 利得 10³ V/A（高インピーダンス負荷時）で帯域 DC–150 MHz、CMRR >25 dB（typ. >30 dB）、
  CW 飽和 4.5 mW @ 1550 nm、RF 出力の最大振幅 ±4.6 V、DC オフセット < ±15 mV [manual]。
  DSC1: 入力 16 bit・1 MHz サンプル（2 点平均で実効 500 kHz）、出力 12 bit、入出力 ±4 V、帯域 DC–100 kHz、入力→出力の −180° 周波数 >58 kHz（typ. 60 kHz）、
  入力 RMS 雑音 0.3 mV（100 Hz–250 kHz）、PID 更新 500 kHz、積分は「積分の中に Ki」を入れる bumpless 形、出力インピーダンス 220 Ω [manual]。
  **915B、Model 914、Cobrite DX1 の資料は `refs/manuals/` にない**（§6-2）。

## 2. 効果の一覧

列: ID / 効果 / 物理機構（1 行） / 構成 / 周波数帯 / スケーリング / 概算 / 段 / 担当トピック（§4 の番号）

### 2.1 光源・スクイーザー（OPA / OPO）

| ID | 効果 | 機構 | 構成 | 周波数帯 | スケーリング | 概算 | 段 | Topic |
|---|---|---|---|---|---|---|---|---|
| E01 | 導波路内の分布損失を伴う OPA squeezing | 生成と減衰が同時に起きるので、実効損失は伝搬損失と利得の両方で決まる | A, F | 位相整合帯域（THz）まで | 伝搬損失 × 導波路長、利得とともに増える | 0.1 dB/cm・45 mm・640 mW で 2% [Ref A §4.2] | O | T01 |
| E02 | ポンプ枯渇 | 信号側へのエネルギー移動でポンプが減る | A（phase-detection OPA の mW probe）、F（高利得パルス） | — | 信号パワー / ポンプパワー | squeezed vacuum では無視できる [AE] | O（パルスでは結合波方程式） | T01 |
| E03 | ポンプ強度雑音 → 利得の揺らぎ | r = √(αP) が揺らぎ R± が揺らぐ | A, B, C, F | DC–検出帯 | δr/r = δP/(2P) | RIN 10⁻³ なら δr/r ≈ 5×10⁻⁴、小さい [AE] | O 表 ＋ T 雑音源 | T01 |
| E04 | ポンプ位相雑音 → squeezing 角の揺らぎ | squeezing 角はポンプ位相の半分に従う（符号は conventions で決める） | A, B, C, F | DC–MHz | ポンプと LO の非共通光路長に比例 | — | T（経路ごとの位相） | T01 |
| E05 | 位相不整合の温度ドリフト | 結晶温度で Δk が変わり、利得が下がり、squeezing 角が周波数に依存して回る | A, B, C, F | < Hz | dΔk/dT × 結晶長 | — | O 表 ＋ T（遅いドリフト） | T01 / T02 |
| E06 | ポンプ誘起の吸収（GRIIRA / BLIIRA / 光屈折） | 強いポンプで結晶の吸収が増え、損失がポンプの関数になる | A, B, F | 遅い（s–h） | ポンプ強度、波長 | L(x) = 0.00249 + 0.00222x、短波長側では時間とともに増加 [Ref B §3–4] | O（パラメータを関数で） ＋ T（遅いドリフト） | T01 / T02 |
| E07 | 空間多モード性（準単一モード導波路） | squeezed 光と probe の空間モードの差で LO との重なりが下がる | A, F | — | コア寸法 | モード不整合 4% のうち原因不明 2% [Ref A §4.2] | O | T01 / T12 |
| E08 | OPO キャビティの伝達（Ref B Eq. 2） | Lorentz 型の squeezing スペクトル、取り出し効率 ρ = T/(T + L) | B, C | 0 – 線幅（C: 1.75 GHz 当てはめ） | ρ、x、Ω | B: γ_c = c(T + L)/l ≈ 7.6×10⁷ s⁻¹。Ω = f/γ_c の 2π の扱いを要確認 [AE] | O | T02 |
| E09 | OPO の離調・キャビティ長雑音 | 離調で squeezing 角が周波数とともに回り、長さ雑音が位相雑音になる | B, C | DC–kHz | フィネス、長さ雑音 | C はモノリシックで長さロックの記述なし | O ＋ T | T02 |
| E10 | パラメトリック利得 G と x の定義 | G = 1/(1 − x)²、x = √(P/P_th)（Eq. 4, 5）と、本文の x = P/P_th が食い違う | B | — | — | P = 100 mW・P_th = 180 mW を Eq. 4–5 に入れると G ≈ 15.4、報告値は 18.7 [AE] | O | T02 |
| E11 | パルス励起の多モード squeezing（supermode） | 結晶内の GVM・GVD で squeezing が複数の時間モードに分かれ、LO の時間モードとの重なりで決まる | F | THz | 結晶長 × GVM と τ_p の比 | PPLN の 775/1550 nm の GVM を 約 0.3 ps/mm とすると 45 mm で 約 10 ps、100 fs ≪ 10 ps [AE, GVM 値 UNVERIFIED] | O | T03 |

### 2.2 レーザー雑音と遅延差のある干渉（D18-1）

| ID | 効果 | 機構 | 構成 | 周波数帯 | スケーリング | 概算 | 段 | Topic |
|---|---|---|---|---|---|---|---|---|
| E12 | レーザーの周波数雑音・位相雑音 | S_ν(f)（白色 FM ＋ 1/f ＋技術的ピーク） | 全部 | DC–MHz | — | Cobrite DX1 の値は資料待ち（§6-2） | T 雑音源 | T05 |
| E13 | 遅延差のある干渉で位相差雑音に変わる | S_Δφ(f) = 4 sin²(πfτ)·S_ν(f)/f²。低周波では (2πτ)²·S_ν | D（ΔL）、F（δL_k、LO_k） | DC – 1/τ | 低周波で τ² に比例 | ΔL = 10 m・Δν = 100 kHz（どちらも仮定）なら全帯域 0.18 rad rms、100 kHz 以下 17 mrad rms [AE] | T（遅延バッファ） | T05 |
| E14 | ポート間のポンプ位相のコヒーレンス（D12 の指標 (ii)） | 32 ポートの squeezed パルスは k·T_rep ずつ違う時刻に生まれ、squeezing 角はその時刻のポンプ位相を引き継ぐ | F | DC – f_rep | 遅延 k·T_rep | 31 × 10 ns・Δν = 100 kHz（仮定）で 0.44 rad rms。モード同期レーザーの櫛線幅しだい [AE] | T | T05 |
| E15 | 強度雑音（RIN） | CMRR を通って BHD に、不平衡を通って DC ロックに、SPM を通って LO 位相に入る | 全部 | DC–MHz | — | — | T 雑音源 | T05 |
| E16 | レーザーの過剰雑音ピーク | 89 kHz、177 kHz のピーク（変調周波数をここから外した） | B | 89、177 kHz | — | [Ref B §3] | T 雑音源 | T05 |
| E17 | EDFA の ASE と LO フィルタキャビティ | 広帯域 ASE と広帯域 squeezed 光のビートが BHD の雑音になる | A, F | 広帯域 | ASE の PSD × LO | フィルタキャビティで抑圧 [Ref A §3.1] | O（スペクトル） | T05 |
| E18 | フィルタキャビティ / MCC のロック残差 | 離調で透過光の位相・振幅が変わり、LO の位相雑音になる | A, B, F | DC–kHz | キャビティ線幅 | — | T | T05 |
| E19 | パルスのタイミングジッタ | LO パルスと squeezed パルスの時間のずれで時間モードの重なりが下がる | F | DC – f_rep/2 | — | — | T ＋ O 表 | T03 |

### 2.3 ファイバ: 線形の伝搬と環境

| ID | 効果 | 機構 | 構成 | 周波数帯 | スケーリング | 概算 | 段 | Topic |
|---|---|---|---|---|---|---|---|---|
| E20 | ファイバの伝搬損失 | Rayleigh 散乱と吸収 | C, D, F | — | L に比例 | C: 結合・取り出し込みで 8% [Ref C §3]。F: f_rep = 10 MHz なら L_31 ≈ 630 m、PM 0.5–1 dB/km で 7–14% [AE, PM の損失 UNVERIFIED] | G（損失項） | T06 |
| E21 | 色分散（β2、β3） | パルスが広がる。LO 側とファイバの種類・長さが違うとモードの重なりが下がる | F | — | L_D = T0²/\|β2\| | β2 ≈ −22 ps²/km（SMF、UNVERIFIED）なら L_D は 100 fs で 0.45 m、1 ps で 45 m、10 ps で 4.5 km [AE] | G | T06（モードの重なりは T03） |
| E22 | 温度による位相ドリフト | dφ/dT = (2π/λ)(dn/dT + n·α_L)·L | C, D, F | < 1 Hz | L·ΔT に比例 | 約 44 rad/(m·K)、63 m で 約 2800 rad/K [AE] | T 雑音源 | T08 |
| E23 | 音響・振動による位相雑音 | 歪みと圧力で光路長が変わる | C, D, F | Hz–kHz | L、巻き方、設置 | 一次資料のモデルか値が必要（§6-6） | T 雑音源 | T08 |
| E24 | 熱力学的な位相雑音（基底） | 温度の熱揺らぎ（Wanser 型） | C, F | Hz–kHz | L に比例 | 環境雑音に比べ小さいと予想 [AE] | T 雑音源 | T08 |
| E25 | ファイバ間の雑音相関 | 同じ電線管・束の中のファイバは同じ環境雑音を受ける。相関係数は周波数に依存 | C, F（squeezed と LO_k、2 ポート組） | DC–kHz | 間隔、束ね方 | C は定性的な記述のみ [Ref C §4] | T（相関係数つき雑音源） | T08 |
| E26 | 群遅延の温度ドリフト（タイミング） | dτ_g/dT でパルスの到着時刻がずれる | F | < 1 Hz | L に比例 | 約 39 ps/(km·K) → 63 m で 2.5 ps/K。100 fs パルスなら 0.04 K で 1 パルス幅ずれる [AE] | T（タイミング状態） ＋ O 表 | T08（係数）/ T03（帰結） |
| E27 | ストレッチャの巻線部の曲げ損失・応力複屈折 | 小径に巻くと損失が増え、PM の PER が下がる | D, F | — | 巻き径、巻き数 | 915B は 38 巻 [D02] | O（静的） | T16 / T09 |
| E28 | 後方反射とエタロン | コネクタなどの反射が寄生干渉計を作り、位相・強度の雑音になる | D, F | DC–kHz | 反射率の積 | FC/APC で十分小さいはず [AE] | T（小） | T12 |

### 2.4 ファイバ: 非線形と量子雑音

| ID | 効果 | 機構 | 構成 | 周波数帯 | スケーリング | 概算 | 段 | Topic |
|---|---|---|---|---|---|---|---|---|
| E29 | SPM（LO、probe パルス） | φ_NL = γ·P_peak·L_eff。RIN が LO の位相雑音になり、スペクトルが広がってモードの重なりが下がる | F（パルス LO）、C/D（CW では小） | — | P_peak、L | LO P_avg = 1 mW・f_rep = 100 MHz・1 ps・63 m・γ = 1.3 /(W·km) なら 約 0.8 rad、100 fs なら 約 8 rad [AE, 入力はすべて仮定] | G | T06 |
| E30 | XPM（probe・LO → squeezed 光） | 同じファイバの強い光が squeezed 光の位相を変える。逆向き伝搬では走り抜けで平均される | F（probe が同じファイバ） | — | P_probe、L | — | G | T06 |
| E31 | 誘導ラマン / パルス内ラマン（SSFS） | fs パルスの LO でスペクトルが赤方偏移する | F（fs LO） | — | P_peak、τ_p | — | G | T06 |
| E32 | 自発ラマン雑音（熱フォノン） | 強い光が squeezed モードへ光子を散乱する。熱フォノン数 n_th(T) に比例 | F（probe が同じファイバ）、C（squeezed と LO が別ファイバなので小） | 広帯域（THz）、パルスの squeezing 帯域に入る | P_probe·L·(n_th + 1) | — | Q | T07 |
| E33 | GAWBS | クラッドの熱的な横方向音響モードが位相・偏光を変調する | C（LO を別ファイバにして避けた [Ref C §1]）、F | 約 20 MHz – 1 GHz（離散ピーク） | 分散は T·L に比例、ピーク周波数はクラッド径で決まる | パルス光ではパルスごとのランダムな位相 [AE] | Q（＋ T の位相雑音） | T07 |
| E34 | SBS の閾値（CW の LO・probe） | 強い CW 光が後方へブリルアン散乱し、パワーが頭打ちになり雑音が増える | C, F（長いファイバに CW） | — | P_th ≈ 21·A_eff/(g_B·L_eff) | 1 km で 約 30 mW、63 m で 約 0.5 W [AE, g_B・A_eff は UNVERIFIED]。ps パルスでは抑えられる | O（チェック） | T07 |
| E35 | 自発ブリルアン後方散乱 | 逆向きの probe から約 11 GHz シフトした光が squeezed 光の向きに戻る | F | 約 11 GHz [UNVERIFIED] | P_probe·L·T | — | Q / O | T07 / T10 |

### 2.5 偏光と PM ファイバ（D18-3）

| ID | 効果 | 機構 | 構成 | 周波数帯 | スケーリング | 概算 | 段 | Topic |
|---|---|---|---|---|---|---|---|---|
| E36 | PM ファイバの PER とクロストーク | 接続・コネクタ・応力で直交軸へ漏れる | D, F | 静的 ＋ 遅いドリフト | 部品ごと | 1 部品あたり 20–30 dB [AE] | O ＋ T | T09 |
| E37 | PM の複屈折の温度依存 | dB/dT で 2 軸の位相差がドリフトする | F（直交偏光の probe）、D（PER が有限のとき） | < 1 Hz | L·ΔT に比例 | 2.4 rad/(m·K)、60 m で 約 150 rad/K（D11 の概算、未検証）。dB/dT ≈ 6×10⁻⁷ /K に相当 [AE] | T | T09 |
| E38 | SMF の偏波ドリフト | ランダムな複屈折の変化で LO との偏光の重なり（可視度）が下がる | C, F（SMF の比較設定） | < 1 Hz – Hz | L、環境 | C は PC で調整 [Ref C Fig. 1] | T | T09 |
| E39 | 偏光子の消光比・偏波コントローラ | MZI の 2 腕の偏光がずれると縞の可視度が下がる | D | < 1 Hz | — | — | T | T09 / T15 |
| E40 | PDL、PMD | 偏光に依存する損失、偏光モード分散 | F（SMF） | — | PMD は √L に比例 | 短いファイバでは小さい [AE] | O | T09 |

### 2.6 probe と squeezed 光の相殺の崩れ（D18-5）

| ID | 効果 | 機構 | 構成 | 周波数帯 | スケーリング | 概算 | 段 | Topic |
|---|---|---|---|---|---|---|---|---|
| E41 | 直交偏光の probe の複屈折ドリフト | probe が測る位相と squeezed 光の位相が E37 の分だけ食い違う | F | < 1 Hz | L·ΔT に比例 | E37 と同じ | T（区間ごとに軸を持つ、D11） | T10 |
| E42 | 逆向き伝搬による時刻のずれ | 位置 z の雑音を probe と squeezed 光は違う時刻に受けるので、高い周波数では相関が消える | F | 約 c/(2·n_g·L) より上 | 1/L に比例 | 63 m で 約 1.6 MHz [AE] | T（区間ごとの遅延バッファ） | T10 |
| E43 | PER による probe → squeezed 光の偏光への漏れ | probe のコヒーレントな光子が squeezed モードに入る（GBS の光子計数を汚す） | F | — | P_probe × 10^(−PER/10) | — | O ＋ T | T10 |
| E44 | Rayleigh 後方散乱による probe の混入 | 逆向きの probe の Rayleigh 散乱光は squeezed 光と同じ向き・同じ周波数 | F | DC –（散乱光の位相は揺らぐ） | 捕捉率 S × α_s × L × P_probe | S ≈ 1.5×10⁻³、α_s ≈ 0.15 dB/km（UNVERIFIED）、63 m・1 mW CW で 約 3 nW、10 ns スロットあたり 約 10² 光子 [AE] | O ＋ T | T10 |
| E45 | probe を入れ出しする素子の損失 | サーキュレータ・PBS・WDM の挿入損失が squeezed 光にかかる | F | — | — | サーキュレータで 0.5–1 dB 程度 [AE, UNVERIFIED] | O（静的） | T10 / T12 |
| E46 | 対策の候補 | 波長をずらした probe（WDM で分離）、周波数シフト、同じ偏光で逆向き、時間ゲート。波長が違うと熱係数が違い差動ドリフトが残る | F | — | — | — | 設計の比較 | T10 |

### 2.7 時間での振り分けと 32 ポートの遅延線（D18-2、D08）

| ID | 効果 | 機構 | 構成 | 周波数帯 | スケーリング | 概算 | 段 | Topic |
|---|---|---|---|---|---|---|---|---|
| E47 | 振り分け素子の挿入損失 | squeezed 光の経路上の部品 | F | — | 素子の種類 | ファイバ型 LN 変調器 3–5 dB、自由空間ポッケルスセルは小（UNVERIFIED） | O（静的） | T11 |
| E48 | 消光比 → 隣のスロットの混入 | 別のスロットの squeezed パルスが混ざる（損失と同時にモード間の相関も生む） | F | — | 10^(−ER/10) | ER 20–30 dB [AE] | O | T11 |
| E49 | スイッチ駆動による位相雑音・バイアスドリフト | 駆動電圧の雑音 × π/Vπ、MZM のバイアスドリフト | F | DC – GHz | Vπ、駆動の雑音 | — | T | T11 |
| E50 | 立ち上がり時間とスロットとの同期 | 遷移中のパルスが一部切られ、損失と位相の誤差になる | F | — | 立ち上がり時間 / T_rep | — | O | T11 |
| E51 | 遅延線長 L_k の誤差と δL_k | 静的な長さの不一致で時間モードの重なりが下がり、レーザー位相雑音の相関も崩れる | F | 静的 | — | f_rep = 100 MHz なら L_1 ≈ 2.04 m、L_31 ≈ 63 m [AE] | O ＋ T | T11 / T03 |

### 2.8 損失とモード整合

| ID | 効果 | 機構 | 構成 | 周波数帯 | スケーリング | 概算 | 段 | Topic |
|---|---|---|---|---|---|---|---|---|
| E52 | 光学部品の損失（自由空間、ファイバ結合、コネクタ、スプライス、カプラ） | — | 全部 | — | 部品数 | A: 光学部品 1% [Ref A]。C: ファイバ前 15%、ファイバ 8% [Ref C]。FC/APC 1 対 0.1–0.3 dB（UNVERIFIED）で 10 dB 目標の損失収支に大きく効く [AE] | O（静的） | T12 |
| E53 | 空間モードの不整合（squeezed 光と LO） | — | 全部 | — | — | A: 約 4% [Ref A §4.2]。B: ξ = 0.988 [Ref B §4] | O | T12 |
| E54 | 時間・スペクトルモードの不整合（パルス） | E11・E19・E21・E26・E51 の結果 | F | — | — | — | O 表 | T03 |
| E55 | squeezed 光の経路上のアクチュエータ・変調器の挿入損失 | D10 と prompt §5.5 | F, C（EOM は LO 側） | — | — | ファイバ EOM 2–3 dB（UNVERIFIED）、ストレッチャは小 [AE] | O（静的） | T12 |
| E56 | conventional tap のトレードオフ | tap 比 r で損失 r が増え、誤差信号の SNR は r に比例して上がる | A（比較用トポロジー） | — | r | prompt §7-8 | O ＋ T | T12 |

### 2.9 検出と測定の模擬

| ID | 効果 | 機構 | 構成 | 周波数帯 | スケーリング | 概算 | 段 | Topic |
|---|---|---|---|---|---|---|---|---|
| E57 | PD の量子効率 | — | 全部 | — | — | A: 約 99%（変換損失 約 1%）。B: 0.998。C: 97% | O（静的） | T13 |
| E58 | 検出器の帯域と周波数応答（干渉ディップ） | 有限の帯域、TIA の応答、2 PD 間の遅延による干渉 | A（35 MHz）、C（0.83・0.96 GHz のディップ）、D（DC–150 MHz） | 帯域の端 | — | [Ref A §3.2、Ref C Fig. 5、manual] | O ＋ M | T13 |
| E59 | 回路雑音（暗雑音）のクリアランス | 等価的な損失 | 全部 | 測定周波数 | 1/P_LO に比例 | A: 28 dB → <0.2%。B: −21.7 dB。C: 11 → 5 dB | O ＋ M | T13 |
| E60 | CMRR と平衡のずれ | LO の古典雑音（RIN、ASE）が漏れる | 全部 | — | ε | PDB450C: >25 dB [manual] | T | T13 |
| E61 | PD の飽和・非線形 | — | D, F | — | — | D: PD 1 個あたり 約 1.7–3.3 mW、飽和 4.5 mW に近い [AE] | O | T13 |
| E62 | パルス BHD 特有の効果 | パルスごとの積分、パルス間のクロストーク（リンギング）、f_rep に対する帯域 | F | f_rep | — | — | O | T13 |
| E63 | ESA の模擬 | RBW・VBW の平均、ゼロスパン、ログ平均のバイアス（ガウス雑音で 約 −2.5 dB、UNVERIFIED）、位相ジッタがあると分布が変わりバイアスが打ち消されない | A, B, C | — | — | — | M | T04 |
| E64 | 位相揺らぎ → 測定される squeezing | R'± ≈ R± cos²θ + R∓ sin²θ（小角）と、分布での厳密な平均。どの帯域の位相雑音を θ_rms に入れるか（CW の ESA と、パルスごとの GBS で違う） | 全部 | — | — | A: 10 mrad で等価損失 約 1%（クリアランス >30 dB）[Ref A §1] | O ＋ M | T04 |
| E65 | 速い位相雑音と測定周波数 | 測定周波数付近の位相変調は直交成分ではなく側帯波を混ぜる（dither の寄与） | B（295 kHz）、C（1 MHz） | 変調周波数 | — | C: dither ±4.5° → 0.056 rad rms、全体 約 0.1 rad [Ref C §3] | O ＋ M | T04 |

### 2.10 誤差信号

| ID | 効果 | 機構 | 構成 | 周波数帯 | スケーリング | 概算 | 段 | Topic |
|---|---|---|---|---|---|---|---|---|
| E66 | ヘテロダイン（AOM 1 MHz）と復調 | probe と LO・ポンプとのビートを復調 | A | 1 MHz | — | BHD のスペクトルに 1 MHz のピーク [Ref A Fig. 4] | T | T14 |
| E67 | phase-detection OPA の信号 | 利得のポンプ・probe 位相依存。probe とポンプのパワーを独立に選べる | A | — | P_probe、P_pump | 新方式は probe 数 mW [Ref A §2] | T | T14 |
| E68 | phase-detection OPA までの非共通光路 | 分岐後の 2 経路の差が誤差信号に入り、squeezer の位相を補正しない | A, F | DC–kHz | 非共通光路長 | 自由空間 cm では小、ファイバで分けると大 [AE] | T（経路ごとの位相） | T14 |
| E69 | dither ＋ DC 極値ロック（Ref B lock 1） | 透過 probe の強度を最大・最小に | B | 295 kHz | 変調深さ | — | T | T14 |
| E70 | probe–LO 干渉 ＋ dither 復調（Ref B lock 2） | BHD 出力を復調 | B | 295 kHz | — | — | T | T14 |
| E71 | 量子雑音ロック | 帯域制限した雑音パワー → 包絡線・ログ検出 → 復調。有限帯域・有限時間の推定の揺らぎ（radiometer 雑音）がセンサ雑音 | C | 8–300 MHz → 1 MHz | (1/Δω)^{1/4} [McKenzie Eq. 13–14] | 全体 約 0.1 rad [Ref C] | T | T14 |
| E72 | probe の BHD への漏れ | 残った probe と LO のビート、probe の古典雑音 | A（conventional）、B、F | 変調・シフト周波数 | P_probe | conventional では数 µW [Ref A §2] | T | T14 |
| E73 | ショット雑音で決まる位相検出の雑音 | S_φ ≈ hν/(η·P_det)（係数は要確認） | probe を使う全部、D | 白色 | 1/P に比例 | 1 µW で 約 1.3×10⁻¹³ rad²/Hz、小さい [AE] | T | T14 |
| E74 | 復調の RF 参照の位相雑音・ミキサのオフセット | ロック点がずれる | A, B, C | DC–kHz | — | — | T | T14 |
| E75 | DC side-fringe lock の誤差の伝わり方（D18-4） | 差信号 ∝ ℜ·G_TIA·P·[ε + V cos φ]。φ = π/2 では強度雑音は ε を通してだけ入り、オフセットドリフト δV_off は δφ = δV_off/K_D | D | DC–100 kHz | P、V、ε | K_D ≈ 1 A/W × 3.33 mW × 10³ V/A ≈ 3.3 V/rad（損失を除く）、DSC1 入力雑音 0.3 mV → 約 0.1 mrad [AE] | T | T15 |
| E76 | 可視度の変化 | 偏光・パワー比のずれで K_D が変わり、ループ利得が変わる | D | < 1 Hz | V | — | T | T15 |

### 2.11 アクチュエータ・制御・デジタル

| ID | 効果 | 機構 | 構成 | 周波数帯 | スケーリング | 概算 | 段 | Topic |
|---|---|---|---|---|---|---|---|---|
| E77 | EOM | Vπ、帯域、範囲、残留振幅変調（RAM）、挿入損失 | A, C, F | DC – GHz | — | C: Uπ = 4 V [Ref C §3] | T | T16 |
| E78 | PZT ファイバストレッチャ（915B） | rad/V、機械共振、レンジ、ヒステリシス、クリープ、偏光の変調 | C, D, F | DC – kHz | — | 915B の値は資料待ち（§6-2） | T | T16 |
| E79 | PZT ミラー | 機械共振 | B（14・22 kHz、単層 PZT では 約 2 kHz）、C | DC – kHz | — | [Ref B §3] | T | T16 |
| E80 | 高圧アンプ（Model 914 など） | 帯域、容量負荷でのスルーレート、出力雑音 → 位相雑音 = rad/V × 雑音電圧 | B, C, D | DC – kHz | PZT 容量 | Model 914 の値は資料待ち | T | T16 |
| E81 | ADC/DAC の量子化とレンジ | — | 全部 | — | ビット数 | DSC1: 入力 LSB 122 µV、出力 LSB 1.95 mV（±4 V、16 / 12 bit）[manual から計算] | T | T16 |
| E82 | ループ遅延とサンプリング | — | 全部 | — | — | DSC1: −180° が 約 60 kHz → 純遅延なら 約 8 µs 相当 [AE] | T（遅延バッファ） | T16 |
| E83 | レンジ飽和、ロック外れ、再引き込み | 遅いドリフトがアクチュエータのレンジを超えるとロックが外れる | C, D, F | < 1 Hz | 熱ドリフト × L / レンジ | C: 数分 [Ref C §4]。F: 63 m・0.1 K で 約 280 rad [AE] | T | T16 |
| E84 | 制御器の構造 | 並列の帯域分割、スタック（Ref C）、PID の離散化、DSC1 の bumpless 積分、加算器、ランプ | 全部 | — | — | — | T | T16 |

## 3. 主な概算の根拠（[AE] の式）

- **Ref A の検算**: R± = L + (1 − L)·exp(±2√(αP))、R'± = R± cos²θ + R∓ sin²θ。L = 0.08、α = 9.06 /W、P = 0.64 W、θ = 0.009 rad → R'− = −10.15 dB、R+ = +20.56 dB。
- **Ref B の G**: x = √(P/P_th) = √(100/180) = 0.745、G = 1/(1 − x)² = 15.4（報告値は 18.7）。逆に G = 18.7 なら x = 0.769、P = 106 mW。
  γ_c = c(T + L)/l = 3×10⁸ × 0.1268 / 0.5 ≈ 7.6×10⁷ s⁻¹。Ω = f/γ_c で f が Hz なら 2π の扱いを確認する必要がある。
- **遅延差のある干渉**: S_Δφ(f) = 4 sin²(πfτ)·S_ν(f)/f²。白色 FM（S_ν = Δν/π）で τ ≪ 1/Δν なら、全帯域の分散 σ² ≈ 2πΔν·τ、
  低周波 S_Δφ ≈ 4πτ²Δν。τ = n_g·ΔL/c = 1.468 × 10 m / c = 49 ns、Δν = 100 kHz（仮定）→ σ² = 0.031 rad²、S_Δφ = 3.0×10⁻⁹ rad²/Hz、100 kHz までで 17 mrad rms。
- **ポート間（指標 (ii)）**: 同じ式で τ = 31 T_rep = 310 ns（f_rep = 100 MHz）→ σ ≈ 0.44 rad（Δν = 100 kHz を仮定）。パルス光源の光学的コヒーレンス（櫛の線幅）で大きく変わる。
- **熱位相**: dφ/dT / L = (2π/λ)(dn/dT + n·α_L)。dn/dT ≈ 1.0×10⁻⁵ /K、α_L ≈ 5.5×10⁻⁷ /K（どちらも UNVERIFIED）→ 約 44 rad/(m·K)。
- **群遅延の熱係数**: (dn_g/dT + n_g·α_L)/c ≈ 39 ps/(km·K)（同じ仮定）。
- **ポートの長さ**: L_k = k·c·T_rep/n_g。f_rep = 100 MHz なら L_1 = 2.04 m、L_31 = 63 m。f_rep = 10 MHz なら L_31 ≈ 630 m。f_rep = 1 GHz なら L_31 ≈ 6.3 m。
- **SPM**: φ_NL = γ·P_avg/(f_rep·τ_p)·L。γ = 1.3 /(W·km)（UNVERIFIED）、1 mW、100 MHz、1 ps、63 m → 0.82 rad。
- **逆向き伝搬の相関周波数**: c/(2·n_g·L) = 3×10⁸ / (2 × 1.468 × 63) ≈ 1.6 MHz（D11 の「通過時間の逆数より上では相関しない」の目安）。
- **Rayleigh**: 後方散乱の割合 ≈ S·α_s·L（αL ≪ 1）。S ≈ 1.5×10⁻³、α_s = 3.5×10⁻⁵ /m → 63 m で 3.3×10⁻⁶。1 mW なら 3.3 nW ≈ 2.6×10¹⁰ 光子/s → 10 ns で 約 260 光子。
  直交偏光にすれば PER の分だけ減るが、GBS（光子計数）で許される量は GBS_REQUIREMENTS が Open なので決められない（§6-7）。
- **DC side-fringe lock**: 差信号 = ℜ·G_TIA·P_in·η_opt·[ε + V cos φ]（P_in は MZI への入力、η_opt は MZI の透過）。K_D = ℜ·G_TIA·P_in·η_opt·V。
  DSC1 の入力雑音 0.3 mV rms → δφ = 0.3 mV / 3.3 V/rad ≈ 0.09 mrad（η_opt = V = 1 とした上限）。

## 4. deep-dive の分担（提案）

方針: 1 トピック 1 エージェント、ノートは `notes/<slug>.md`。依存は「数値を受け取る」関係だけで、並行して始められる（依存先の値が出るまではリファレンス値か仮の値を使い、そう書く）。
すべてのトピックは prompt §4 のノートの形式（式と記号の定義、仮定と適用範囲、周波数帯、スケーリング、概算、検証テスト、文献）に従う。

| # | slug | 担当する効果 | 最初に読む文献 | 依存 |
|---|---|---|---|---|
| T01 | `opa-waveguide-squeezing` | E01–E07（OPA 側） | Ref A、Ref A の文献 [19]、[23]、[28] | — |
| T02 | `opo-cavity-squeezing` | E05・E06（OPO 側）、E08–E10、prompt §7-9 の Eq. (2) の確認 | Ref B、Ref C、Ref B の文献 [8]–[10]（Polzik 1992、Lam 1999、Zhang 2003） | — |
| T03 | `pulsed-temporal-modes` | E11、E19、E21 のモード重なり、E26 の帰結、E51 のタイミング、E54 | Ref A の文献 [15]、[16]、[21]（パルス・広帯域 squeezing）。他は要調査 | T06（分散）、T08（熱係数） |
| T04 | `phase-noise-to-squeezing` | E63–E65、D12 の 2 つの指標の定義、D16 の「目標 × L → 必要な θ_rms」表の方法、prompt §7-1・§7-2 | Ref A Eq. と文献 [25]–[27]、Ref B Eq. (1)、Ref C の文献 [22]（Franzen 2006） | — |
| T05 | `laser-noise-delay-interferometry` | E12–E18（D18-1） | Cobrite DX1 資料（§6-2）、Ref A §4.1、Ref B §3 | — |
| T06 | `fiber-gnlse-classical` | E20、E21、E29–E31、RK4IP と prompt §7-5 | Agrawal「Nonlinear Fiber Optics」、Dudley et al. RMP 2006、Hult JLT 2007（いずれも UNVERIFIED、要確認） | — |
| T07 | `fiber-quantum-noise-scattering` | E32–E35、量子雑音層の方式（線形化 vs truncated Wigner）、prompt §7-6 | Shelby, Levenson, Bayer, PRB 1985（GAWBS）、Drummond & Corney, JOSA B 2001（UNVERIFIED、要確認）、Ref C §1 | — |
| T08 | `fiber-environmental-phase-noise` | E22–E26、E27 の熱的な部分 | Ref C、ファイバ雑音補償の文献（Ma et al. Opt. Lett. 1994、UNVERIFIED） | — |
| T09 | `polarization-pm-fiber` | E36–E40（D18-3） | PM ファイバのメーカー資料、D11 | — |
| T10 | `probe-cancellation-breakdown` | E41–E46（D18-5） | D10、D11、D18-5、T07–T09 の結果 | T07、T08、T09 |
| T11 | `time-demux-switch` | E47–E51（D18-2、D08） | 要調査（時間多重 CV 実験の文献、変調器のメーカー資料） | T03（タイミング） |
| T12 | `loss-and-mode-matching` | E28、E52、E53、E55、E56、prompt §7-3・§7-8 | Ref A §4.2、Ref B §4、Ref C §3 | — |
| T13 | `bhd-detection` | E57–E62、PDB450C の値 | Ref A §3–4、Ref C Fig. 5、`refs/manuals/thorlabs-pdb.pdf` | — |
| T14 | `error-signal-generation` | E66–E74、prompt §5.5 の誤差信号 1–4 | Ref A §2–3、Ref B §2–3、Ref C §2、McKenzie 2005 | T13（C の検出帯域） |
| T15 | `dc-fringe-lock-ref-d` | E39 の D 部分、E75、E76、D05・D06 の DSC1 プリセット、Ref D 全体の解析モデル（prompt §7-7 と D04） | `refs/manuals/thorlabs-dsc1.pdf`、`thorlabs-pdb.pdf`、915B・914・DX1 の資料 | T05、T16 |
| T16 | `actuators-control-digital` | E27 の機械的な部分、E77–E84、機器プリセット（Moku:Delta、DSC1、必要なら Moku:Pro）、prompt §7-7・§7-10 | `refs/manuals/`、Ref B §3、Ref C §2–4、915B・914 の資料 | — |

依存のまとめ: T10 ← T07・T08・T09、T03 ← T06・T08、T11 ← T03、T14 ← T13、T15 ← T05・T16。それ以外は独立。

補足:
- `notes/conventions.md`（prompt §3、全位相の符号の約束）はどの Phase の担当かが決まっていない（§6-5）。T04 か調整役が Phase 2 の前に作るのが自然。
- 「どの帯域の位相雑音が squeezing を下げるか」は T04 が全トピック共通の写像を定める。他のトピックは位相雑音の PSD を出すところまでを担当する。

## 5. 文献

このサーベイで実際に読んだもの:

- Ref A: K. Hirota, T. Kashiwazaki, et al., "Generation of 10-dB squeezed light from a broadband waveguide optical parametric amplifier with improved phase locking method."
  **掲載誌・年・DOI は PDF の本文から読み取れなかった（UNVERIFIED）。**
- Ref B: Y. Takeno, M. Yukawa, H. Yonezawa, A. Furusawa, "Observation of −9 dB quadrature squeezing with improvement of phase stability in homodyne measurement,"
  Opt. Express 15, 4321–4327 (2007)。DOI は PDF に記載がない（UNVERIFIED）。
- Ref C: S. Verclas, C. Darsow-Fromm, B. Tohermes, R. Schnabel, "Fiber distribution of phase-stabilized GHz-bandwidth squeezed vacuum states of light between two buildings,"
  Optica Quantum 4, 1 (2026)。DOI 10.1364/OPTICAQ.567418（PDF に記載）。
- K. McKenzie, E. E. Mikhailov, K. Goda, P. K. Lam, N. Grosse, M. B. Gray, N. Mavalvala, D. E. McClelland, "Quantum noise locking,"
  J. Opt. B: Quantum Semiclass. Opt. 7, S421–S428 (2005)。doi:10.1088/1464-4266/7/10/032（PDF に記載）。
- Thorlabs DSC1 manual（Rev. B, November 3, 2025、DOC-103674）、Thorlabs PDB440x/PDB450x manual（© 2019）、Moku:Delta specifications（`refs/manuals/`）。

deep-dive の出発点として挙げたもので、このサーベイでは確認していないもの（**すべて UNVERIFIED**）:
Agrawal「Nonlinear Fiber Optics」、J. M. Dudley, G. Genty, S. Coen, Rev. Mod. Phys. 78, 1135 (2006)、J. Hult, J. Lightwave Technol. 25, 3770 (2007)、
R. M. Shelby, M. D. Levenson, P. W. Bayer, Phys. Rev. B 31, 5244 (1985)、P. D. Drummond, J. F. Corney, J. Opt. Soc. Am. B 18, 139 (2001)、
L.-S. Ma, P. Jungner, J. Ye, J. L. Hall, Opt. Lett. 19, 1777 (1994)、K. H. Wanser, Electron. Lett. (1992)（ファイバの熱力学的雑音）。
Ref B・C の参考文献リストにある文献（Polzik 1992、Lam 1999、Zhang 2003、Franzen 2006 など）は、リストに書かれている書誌情報までは確認済みで、本文は未読。

## 6. 仕様の食い違いと、ユーザの決定が必要な点

1. **worktree の起点**: この worktree は `main`（init commit）から作られていて、`AGENTS.md`・`DECISIONS.md`・`docs/squeezing_sim_prompt.md` がなかった。
   `dev` へ fast-forward して作業した。他の worktree も `dev` から作る必要がある（運用、D22）。
2. **メーカー資料がない**: D02 は「915B の rad/V・共振・レンジ、DSC1 の仕様、Cobrite DX1 の線幅・周波数雑音・RIN はメーカー資料の値だけを使う」と決めているが、
   `refs/manuals/` には DSC1、PDB4xx、Moku:Delta しかない。915B、Model 914、Cobrite DX1、PNH1550R5A2、Opneti 偏光子、Optoquest PC、MN9605C の資料がない。
   ユーザが PDF を `refs/manuals/` に置くか、deep-dive がメーカーのウェブページを引いてよいか（ウェブの値でも出典を書けば D02 を満たすか）、決めてほしい。
3. **Ref C の制御器は Moku:Pro**（Ref C §2、Fig. 3）。prompt は Moku:Delta をプリセットに挙げている。Ref C の再現に Moku:Pro のプリセットが必要か、汎用値でよいか。
4. **Ref B Eq. (2) の定義**: 本文は x = P/P_th、Eq. (5) は x = √(P/P_th)。Eq. (4)–(5) からは G ≈ 15.4、報告値は 18.7。Ω = f/γ_c の 2π も要確認。
   prompt §1 が research agent に確認を求めている項目そのもので、T02 が調べる。ユーザの決定は不要だが、検証目標（§7-9）の数値に効く。
5. **`notes/conventions.md` の担当**: prompt §3 は「実装の前に」とだけ書いていて、Phase と担当が決まっていない。
6. **環境雑音のデータ**: FIBER_PHASE_NOISE_DATA と TEMPERATURE（環境の温度ドリフト・音響 PSD）にはリファレンス値がない（Ref C は定性的な記述のみ）。D23 では「リファレンス値がなければ止めて聞く」。
   文献のモデル（一般的な実験室環境の値）で代用してよいか。ユーザの実験室で測った値を使う場合、D04 のブラインド検証（Ref D）と衝突しないか。
7. **指標 (ii) の「共通の基準」の定義**（D12）: GBS 入力の位相を何に対して測るか（ポンプの位相？ある時刻のレーザー位相？ポート 0？）。
   E14 の概算どおり、ポート間の k·T_rep の遅延によるレーザー・ポンプ位相雑音が (ii) で支配的になりうるので、定義で結果が大きく変わる。
8. **将来の probe の仕様**（D10、D18-5）: CW かパルスか、波長、偏光（直交か同じか）、パワー。さらに GBS の入力に許されるコヒーレント光子の混入量（GBS_REQUIREMENTS、Open）。
   E44 の概算では、1 mW の CW probe の Rayleigh 散乱だけでスロットあたり 約 10² 光子になりうる。T10 の結論はこの許容量で決まる。
9. 小さな点: prompt §1 は Ref A の光源を「CW seed + EDFA + SHG」と書くが、Ref A の本文は「NKT Koheras Harmonik が基本波と第 2 高調波を出す」とだけ書いている（EDFA は図 2 にあるかもしれないが、図はテキスト抽出で読めなかった）。
   また 640 mW は「squeezer OPA 直後のポンプ出力の推定値」で、入力パワーではない。ASE（E17）とポンプ枯渇（E02）の扱いに効く。T01・T05 で図を確認する。
