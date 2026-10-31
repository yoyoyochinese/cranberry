# T11 時間 → 空間の振り分け素子（time-demux-switch）

> **訂正（2026-11-01、DECISIONS D62-4・D77-2）:** 立ち上がりの縁の条件は ≤ (1 + n_guard)·T_rep（N_b·T_rep ではない）。Borealis の時間配分は T17 の読み方（ロック光 65 µs・squeezed 光 35 µs）。木・鎖のスイッチの偏りのドリフトはポート間で共有される（D77、M4-01）。


- Phase: 1（deep dive、T11）
- 日付: 2026-10-28
- 担当する効果: E47–E51（`notes/survey.md` §2.7）。D08、D09、D10、D18-2、D30–D32。
- 依存: T03（パルスの時間モードの重なり）。T03 の結果が出るまでは、ガウス型パルスの重なりの式を仮に使う（§3.6、**placeholder**）。
  熱による群遅延の係数は T08、n_g・GDD は T06、BHD の帯域は T13、probe の経路は T10/T14 の結果で置き換える。
- 値の印: **[ref:X]** = 文献 X の本文から読んだ値、**[mfr]** = メーカー資料（URL と閲覧日を §9 に記載）、
  **[AE]** = この担当による概算（未検証）、**[UNVERIFIED]** = 原典で確認していない値、**[placeholder]** = 他トピックの値が出るまでの仮の値。

## 0. 要約（結論を先に）

1. **振り分けの「漏れ」（消光比 ER が有限）は、squeezed 光にとっては純粋な損失として働く。** ポートごとの遅延が固定で光源が 1 台のとき、
   スロット j の光がポート k に漏れても、GBS 入力では必ず**ターゲットでない時刻のビン**に着く（§3.3 の証明）。
   したがって測定がパルスごとに時間分解されていれば（BHD のゲート、時間分解できる光子検出）、ポート間にコヒーレントな相関は生まれない。
   相関（クリックの相関）が出るのは、検出器が T_rep より長く積分する場合だけ。survey E48 の「モード間の相関も生む」はこの条件付きで直す（§8）。
2. **挿入損失が支配的。** 32 = 2⁵ なので二分木なら各ポートが S = 5 段を通る。
   ファイバ型 LiNbO3（LN）スイッチは 1 段あたり 2–4 dB [mfr: EOSPACE] → 5 段で 10–20 dB → 入力 −10 dB の squeezing は −0.4 dB 以下になり、**使えない**。
   自由空間の Pockels セル（PC）＋PBS は 1 段あたり約 1–4%（[ref: Wang 2017] PC 透過 >99%、[ref: Kaneda 2019] PC 0.8%、[ref: Madsen 2022] 4 段で約 15%）
   → 5 段で η ≈ 0.82–0.96 → −10 dB が −5.8 〜 −8.7 dB。
3. **低損失の素子（PC）は遅い**（繰り返し 0.5–3 MHz [ref: Kaneda 2019, Madsen 2022, Wang 2017]、立ち上がり 8–10 ns [ref: Wang 2017, Takeda 2019]）。
   D08 のようにパルスごとに振り分ける（N_b = 1）なら、1 段目は f_rep/2 で切り替わるので **f_rep ≲ 6 MHz**（Borealis と同じ: 6 MHz の列を 3 MHz の 1 段目で分けている）。
   このとき L_31 ≈ 1.06 km。f_rep を上げたまま PC を使うには、N_b 個の連続パルスをまとめて振り分ける「ブロック振り分け」が必要で、遅延線は N_b 倍に伸びる
   （f_rep = 100 MHz、N_b = 50 なら L_31 ≈ 3.2 km、0.2 dB/km で 14% の損失）。**素子の損失とファイバの長さ（損失・位相雑音・熱ドリフト）のトレードオフ**がこのトピックの中心。
4. **タイミングの要求は τ_p に比例して厳しい。** ガウスパルスで重なり η_mm ≥ 0.99 には時間ずれ δt ≤ 0.085 τ_p（FWHM）：
   τ_p = 100 fs で 8.5 fs（ファイバ長 1.7 µm 相当）。L_31 = 63 m では熱で 2.5 ps/K [AE/placeholder] → 許される温度差 3.4 mK。
   さらに T_rep の相対変化 δT/T がポート 31 とポート 0 の到着時刻を 31·T_rep·δT/T ずらす（100 MHz、10⁻⁶ で 310 fs）。
5. **スイッチ駆動の位相雑音**は δφ = κ_pm·π·δV/V_π（1 段あたり）。5 段で合計 10 mrad に抑えるには、κ_pm = 1 として各段の駆動電圧の揺らぎ δV/V_π ≤ 1.4×10⁻³ rms。
   スイッチが probe のロックループの外にあると（D30 の干渉点がスイッチより下流のとき）、この雑音は補正されず指標 (ii) に直接入る。
6. **LO を全パルス列のまま配る（D09 の「分岐＋遅延」）なら、BHD はパルスごとに分解して読む必要がある。** 積分型（ESA の連続測定）だと、ポート k では 32 ビン中 1 ビンだけが squeezed で、
   −10 dB が −0.12 dB に薄まる（§3.7）。

## 1. 記号（symbols, English）

- `K`: number of output ports (= 32, D08). `k = 0 … K−1`: port index. `s = 1 … S`: switch stage index; binary tree `S = log2 K = 5`.
- `f_rep`, `T_rep = 1/f_rep`: pulse repetition rate / period of the squeezer output (Open, D07).
- `N_b`: number of consecutive pulses routed to the same port per frame ("block size"); D08 corresponds to `N_b = 1`.
- `T_F = K·N_b·T_rep`: frame period (one GBS shot per frame per `N_b` pulses).
- `j`: input time-slot index within a frame (slot j arrives at the switch at `t_j = t_0 + j·T_rep`).
- `σ(k)`: slot (block) assigned to port k. Default `σ(k) = K−1−k` so that `L_k ∝ k` (see §3.1).
- `L_k`: fiber delay length of port k [m]. `L_0`: common offset length [m]. `n_g`: group index of the delay fiber. `c`: speed of light in vacuum.
- `δL_k`: residual static length mismatch between the squeezed path of port k and LO_k (D09) [m]. `δt_k = n_g·δL_k/c` [s].
- `τ_p`: pulse duration (intensity FWHM) [s]. `σ_t = τ_p/(2√(2 ln 2))`: intensity RMS width [s].
- `η_mm`: temporal mode-overlap efficiency (power) between two pulses (placeholder for T03).
- `t_s`: single-pass power transmission of stage s in its "pass" state (insertion loss `IL_s = −10 log10 t_s` [dB]).
- `ER_s`: extinction ratio of stage s [dB]; `ε_s = 10^(−ER_s/10)`: fraction of power sent to the wrong output.
- `η_sw,k`: total power transmission of the switch network for the intended slot into port k.
- `α_f`: fiber attenuation [dB/m]. `η_f,k = 10^(−α_f·L_k/10)`: fiber transmission of port k.
- `η_k = η_sw,k · η_f,k · η_other`: total transmission of port k's squeezed path (other losses from T12).
- `V_in`, `V_out`: quadrature variance before/after (vacuum `V_vac = 1/2`). `R = V/V_vac`: variance normalised to shot noise.
- `V_s(t)`: drive voltage of stage s [V]. `V_π,s`: voltage giving a π change of the switching retardance/phase difference of stage s [V].
- `Γ_s = π·V_s/V_π,s + Γ_b,s`: switching retardance (Pockels cell) or arm phase difference (MZM) of stage s [rad]; `Γ_b,s`: bias (static + drift) [rad].
- `κ_pm,s`: dimensionless coefficient converting a change of `Γ_s` into a change of the common (output) optical phase of the passing light; 0 for an ideal push-pull MZM.
- `ψ_sw,k(t)`: optical phase added by the switch network to light routed into port k (convention §4: `A → A·e^{+iψ}`) [rad].
- `δV_s`: drive-voltage fluctuation at the instant the pulse passes stage s [V].
- `δΓ_b,s`: bias drift of stage s [rad].
- `t_r`: 10–90% rise time of a switching edge [s]. `σ_j`: RMS timing jitter between the switch edge and the optical pulse [s].
- `f_tog,max`: maximum toggle (repetition) rate of the switch driver [Hz].
- `κ_T`: thermal coefficient of fiber group delay per length, `κ_T = d(n_g L/c)/dT / L` [s/(m·K)] (placeholder 39 ps/(km·K), T08).
- `θ_sq,k`: squeezing angle of the pulse in port k (conventions §3). `Δ_k`, `Θ_k`: metrics (i), (ii) (conventions §10).
- `r`: squeezing parameter. `n̄_leak`: mean photon number per frame leaked into non-target time bins of one port.

## 2. 候補の素子と文献の値

| 方式 | 1 段の損失 | 消光比 | 立ち上がり / 繰り返し | 出典・備考 |
|---|---|---|---|---|
| 自由空間 Pockels セル（RTP など）＋PBS（偏光を回して PBS で分ける） | PC 透過 >99% [ref: Wang 2017]；PC 単体 0.8% [ref: Kaneda 2019]；ファイバ結合込みで 4 段 約 15%（1 段あたり約 4% [AE: 0.85^{1/4}]）[ref: Madsen 2022] | PC >100:1（20 dB）、PBS >1200:1 [ref: Wang 2017]；>200:1（23 dB）[ref: Madsen 2022] | 8 ns、0.76 MHz [ref: Wang 2017]；約 10 ns、時間ビン 66 ns ごとに設定 [ref: Takeda 2019]；500 kHz（ドライバ律速）、「市販のバルク PC は 1 MHz」[ref: Kaneda 2019]；1 段目 3 MHz（6 MHz のパルス列）[ref: Madsen 2022] | GBS の振り分けの実績あり。Wang 2017/2019 は連続するパルスの塊（segment）ごとに振り分け（N_b > 1）。Madsen 2022 はパルスごと（6 MHz を 1→16 の二分木、15 個の EOM、4 段: 3 MHz / 1.5 MHz / 750 kHz / 375 kHz） |
| ファイバ型 LN 導波路 MZM スイッチ（1×2、2×2） | <4 dB（<3 dB オプション、約 2 dB はカスタム）[mfr: EOSPACE] | 高消光比版 >30–60 dB（カスタム）[mfr: EOSPACE]；標準品の値は資料にない | 帯域 10–65 GHz、V_π ≈ 4–5 V @ 1 GHz [mfr: EOSPACE] | 速いが損失が大きすぎる（§3.2）。バイアスのドリフトあり（§3.5、Salvestrini 2011） |
| ファイバ型の固体結晶スイッチ（Agiltron NanoSpeed 1×2） | 0.6 dB typ、1.0 dB max（1260–1650 nm）[mfr: Agiltron] | クロストーク 18/25/35 dB（min/typ/max）、PM 版 ER ≥18 dB [mfr] | 立ち上がり ≤300 ns、繰り返し ≤100 kHz または ≤300 kHz（ドライバ）[mfr]；「素子の共振で周波数応答の悪い領域がある」[mfr] | 損失・速さとも足りない（5 段で 3 dB）。比較用 |
| 受動的な分岐＋全ビンを使う（Liu 型の時空間ハイブリッド） | 干渉計 99%/モード [ref: Liu 2026] | — | スイッチなし。遅延ループ（τ = 50 ns、16τ）を PZT で位相ロック | 振り分け素子ではない。全時間ビンを GBS のモードとして使う設計で、D08 の「32 ポートを同時刻に揃える」とは別のアーキテクチャ。比較の参考 |
| AOM、薄膜 LN（TFLN）スイッチ | — | — | — | 文献値を確認していない（UNVERIFIED）。AOM は回折光が f_AOM だけ周波数シフトするので、ポートごとに LO との周波数がずれる点に注意 [AE] |

関連する位相の安定度の参考値（同じ種類の時間多重系）:
- Madsen 2022（Borealis）: ループの位相雑音の標準偏差 0.02、0.03、0.15 rad（ループ 1、2、3。古典パルスで測定）、ロックは逆向きのレーザー光＋PZT＋lock-in、測定 65 µs と安定化 35 µs を交互。
- Liu 2026: 遅延ループの位相ロック精度 3 nm（短いループ）、10 nm（長いループ）→ 1550 nm で 2π×3/1550 = 12 mrad、41 mrad [AE: 換算]。長いファイバの温度ドリフトは TEC のフィードバックを併用。
- Takeda 2019: ループの位相のドリフト 約 7° がメモリ寿命（約 400 ns）の主な制限。ロック（2 kHz で帰還と測定を交互）。

## 3. 支配方程式

### 3.1 遅延線の長さ L_k（D08 の一般化）

フレームの先頭時刻を t_0 とし、スロット j（j = 0 … K·N_b − 1）は t_j = t_0 + j·T_rep にスイッチに着く。ポート k はブロック σ(k)（N_b 個の連続スロット）を受け取る。
全ポートの光が GBS 入力で同時に着くための条件は、ポートごとの遅延 τ_k = n_g L_k / c が

  τ_k = τ_0 + (K − 1 − σ(k))·N_b·T_rep

を満たすこと。σ(k) = K − 1 − k と番号を付けると（**ポート k は最後から k 番目のブロックを受け取る**）

  L_k = L_0 + k·N_b·c·T_rep / n_g  （k = 0 … K − 1）

で、D08 の L_k ≈ k·c·T_rep/n_g は N_b = 1、L_0 = 0 の場合。ブロック内の i 番目（i = 0 … N_b − 1）のパルスは、ポートを問わず t_0 + τ_0 + (K−1)·N_b·T_rep + i·T_rep に着くので、
1 フレームで N_b 回の GBS ショット（それぞれ K モード）になる。

- 符号: L_k の単位は m、τ_k は群遅延（位相 ψ = β L とは別。位相は conventions §4、ロックで扱う）。
- 群遅延の熱変化: dτ_k/dT = κ_T·L_k。

### 3.2 挿入損失（squeezed 光への作用、prompt §5.5）

スイッチ網はポート k へのターゲットのスロットに対し、透過率 η_sw,k のビームスプリッタ（残りのポートは真空）として働く:

  â_out,k = √η_sw,k · e^{iψ_sw,k} · â_in,σ(k) + √(1 − η_sw,k) · v̂_k

（v̂_k は真空、または §3.3 のようにターゲットでないビンの光。どちらでもターゲットの時間モードとは直交）。二分木（S 段）なら

  η_sw,k = Π_{s=1}^{S} t_s·(1 − ε_s(t))  （ε_s(t) は §3.4・§3.5 で時刻に依存しうる）

で、分散は純損失の式（conventions §2）:

  V_out = η_k V_in + (1 − η_k) V_vac、 η_k = η_sw,k·η_f,k·η_other.

squeezing の上限（r → ∞）は R_min = 1 − η_k、すなわち 10 log10(1 − η_k) dB。

**直列（chain）型**の配置（PC と PBS を 1 列に並べ、ポート k は k + 1 段を通る）では η_sw,k が k に依存する。二分木では全ポートが S 段で同じ（ただしファイバ損失 η_f,k は k に依存）。

### 3.3 消光比と隣のスロットの混入（E48）

ε_s のため、スロット j の光の一部はターゲットでないポート k′ ≠ σ⁻¹(j) に入る。ポート k′ の出力時刻は「入力時刻 + τ_k′」なので、スロット j の漏れはポート k′ の出力で

  t_out = t_j + τ_k′ = （ポート k′ のターゲットの到着時刻）+ (j − j_target(k′))·T_rep

に着く。j ≠ j_target(k′) だから、**漏れた光はポート k′ のターゲットのビンと時間的に重ならない**（ずれは T_rep の整数倍、0 でない）。
他のフレームから来た光も同じ（ずれは K·N_b·T_rep の整数倍だけ増えるだけ）。さらに、ターゲットの時刻にポート k′ にあるのは「その時刻にポート k′ に入ったもの」= スロット j_target(k′) だけ。
したがって、遅延がポートごとに一定（時間不変）で光源が 1 つの振り分けでは、

- ターゲットの時間モード（GBS のショット、LO_k のパルス）から見ると、ε はスロットの**損失**（η_sw,k に入る）であって、他のスロットの squeezed 光は混ざらない。
- 漏れた光は、各ポートのターゲットでない時刻のビンに「弱い squeezed パルス」として残る。1 ポート・1 フレームあたりの漏れ光子数は、1 段目の漏れが支配的なら

    n̄_leak ≈ (Σ_s ε_s)·sinh² r  （ε_s ≪ 1、各段で一度だけ漏れる近似）

- これが GBS の結果に効くのは、検出器がこれらのビンを分けられない場合（応答が T_rep より遅い積分型、ゲートなし）だけ。そのときは「ポート k′ の検出器にスロット j の光子が数えられる」ので、
  ポート間の光子数の相関（ビームスプリッタ比 ε のインコヒーレントな混合）として現れる。BHD で LO_k がターゲットの時間モードだけにあるか、パルスごとに読むなら影響はない。

例外（コヒーレントに混ざる場合）: (a) τ_p が T_rep に近く隣のパルスと時間的に重なる、(b) 遅延が T_rep の整数倍からずれている（|δt| ~ τ_p）。
(a) は D07 の範囲（τ_p ≤ 10 ps、T_rep ≥ 1 ns）では起きない。(b) は §3.6 の時間ずれの問題として T03 で扱う。

### 3.4 立ち上がり時間とスロットとの同期（E50）

スイッチ段 s の状態が時刻 t に Γ_s(t) で、パルス（中心 t_a、幅 τ_p ≪ t_r）が通るとき、PC＋PBS（または MZM）の通過側の電力透過は

  T_pass(t_a) = t_s · cos²((Γ_s(t_a) − Γ_target)/2)、 漏れ = t_s · sin²((Γ_s(t_a) − Γ_target)/2)

（Γ_target = 0 または π、通過させたいポートで決まる）。τ_p が t_r と同程度なら、パルス内で透過と位相が変わり時間モードがゆがむ（T03）。
エッジの間にパルスが来ないための条件（パルスごとの振り分け、1 段目の各スロットの間）:

  t_r + 2·(3σ_j) + τ_p + t_ring ≤ N_b·T_rep、  かつ  f_rep / (2N_b) ≤ f_tog,max

（t_ring は駆動の振動が許容値まで減衰する時間。二分木の段 s は f_rep/(2^s N_b) で切り替わるので、1 段目が最も厳しい）。
条件を満たせない場合は、エッジの間のパルスを捨てる（ガードスロット）ことで squeezing を保ったまま**レートだけ**を下げられる。そのときは L_k の式にガードスロット数を入れる（§3.1 の N_b T_rep を周期に置き換える）。

同期: スイッチの駆動はパルス列から取ったクロック（PD → 2^s 分周 → 各段の HV ドライバ）に位相同期させる。遅延線の長さは T_rep に合わせて固定されるので、T_rep の相対変化 δT_rep/T_rep があると
ポート k と ポート 0 の GBS 入力での時間ずれは

  δt_k^{rep} = k·N_b·T_rep·(δT_rep/T_rep).

（LO_k が squeezed パルスと同じ時刻に生まれたパルスで、同じ長さの遅延を通るなら Δ_k には効かず、Θ_k とポート間の干渉（GBS）にだけ効く。）

### 3.5 スイッチ駆動の位相雑音とバイアスのドリフト（E49）

段 s を通過したパルスの光の位相（conventions §4 の符号: A → A e^{+iψ}）は

  ψ_s = ψ_s^{(0)} + κ_pm,s·(π·δV_s/V_π,s + δΓ_b,s)

で、ポート k の合計は ψ_sw,k = Σ_{s ∈ path(k)} ψ_s。

- κ_pm: 理想的な push-pull MZM では 0（両腕の位相変化が逆符号で打ち消す）。片側駆動の MZM では 1/2（一方の腕だけ π δV/V_π 変わり、出力は平均の位相）。
  PC＋PBS では、結晶の 2 つの主軸の位相変化の和と差の比で決まり、結晶と配置に依存する（RTP の横型で r33、r23 を使うと約 1 になる [AE、未検証]）。**値は素子ごとのパラメータとし、メーカーまたは測定で決める。**
- スクイーズ角への写像（conventions §3–4）: â → â e^{iψ} のとき x̂′_θ = x̂_{θ−ψ} なので **θ_sq,k → θ_sq,k + ψ_sw,k**。
- 静的な部分 ψ_s^{(0)}（パターンに依存するが各ポートで決まった駆動波形を見るので一定）はロックの設定値またはオフセットに吸収される。
- 揺らぎの分散（各段の駆動雑音が独立、パルスが通過する瞬間の値をサンプル）:

    σ_sw,k² = Σ_s κ_pm,s²·π²·⟨δV_s²⟩/V_π,s² + Σ_s κ_pm,s²·⟨δΓ_b,s²⟩.

  ポート k のパルスは周期 T_F = K N_b T_rep で同じ段を通るので、駆動雑音のうち f_F/2 = 1/(2T_F) より上の成分は折り返されて、パルスごとに独立な（白色の）位相誤差になる。
  これはロックの帯域より速いので補正できず、conventions §6 の σ² にそのまま加わる。
- バイアスのドリフト δΓ_b は透過と漏れも変える: 透過 ∝ cos²(δΓ_b/2) ≈ 1 − δΓ_b²/4、漏れ ε_b = sin²(δΓ_b/2)。
  LN MZM のドリフトの原因は焦電効果（温度変化）、印加 DC 電圧による電荷の蓄積（DC ドリフト）、光屈折効果（Salvestrini 2011、Wooten 2000。速さの値は確認していない、UNVERIFIED）。
  バイアス制御にディザを使うと、squeezed 光の経路で振幅・位相の変調になる（squeezed 光への余分な変調として扱う必要がある）。
- 振幅（透過）への駆動雑音の効果は、通過点（Γ = 0 または π）では 1 次の感度が 0 で、2 次 (π δV/(2V_π))² の損失揺らぎになる。

### 3.6 静的な長さの不一致 δL_k と時間モードの重なり（E51、T03 の placeholder）

squeezed パルスと LO_k（またはポート間）の時間ずれ δt に対して、変換限界のガウスパルス（強度の RMS 幅 σ_t）の振幅の重なりは
∫u(t)u(t − δt)dt = exp(−δt²/(8σ_t²)) なので、電力の重なり効率は

  η_mm(δt) = exp(−δt²/(4σ_t²))、 σ_t = τ_p/(2√(2 ln 2)).

η_mm は損失として η_k に掛かる（BHD での効果）。チャープ・分散の違いによる重なりの低下（PC 結晶 S 段分の GDD は squeezed 側にだけ入る）は T03/T06 で扱う。**この式は T03 の結果で置き換える。**
δt の内訳: 製作時の長さの誤差 n_g δL_k/c、温度 κ_T L_k ΔT（squeezed 側と LO 側のファイバの温度差・長さの差に効く）、§3.4 の δt_k^{rep}、スイッチの段ごとの光路長差（自由空間の PC 段）。
レーザー位相雑音の相関の崩れ（遅延 τ = n_g δL_k/c）は T05 の S_Δφ(f) = 4 sin²(πfτ) S_φ(f) を使う（conventions §7）。

### 3.7 測定の模擬への条件（M）

ポート k には 1 フレームあたり N_b 個のターゲットのパルスと K·N_b − N_b 個の空き（真空＋漏れ）のビンがある。LO_k が全パルス列（D09 の分岐＋遅延）で、BHD がパルスを分解せずに積分する（ESA で f ≪ f_rep を見る）と、測定される分散はビンの平均

  R_meas = (1/K)·R_target + (1 − 1/K)·R_empty、 R_empty ≈ 1

になる。M ブロックは「パルスごとに分解（ゲート）」か「積分」かを設定で選べるようにする。

## 4. 仮定と適用範囲

- 光源は 1 台、遅延はポートごとに時間不変（§3.3 の結論はこの仮定に依存）。
- τ_p ≪ T_rep（D07: τ_p ≤ 10 ps、T_rep ≥ 1 ns なので常に成り立つ）。
- スイッチは線形の受動素子として扱う（squeezed 光は弱いので、LN の光屈折効果・非線形は無視。強い probe や LO を同じ素子に通す場合は別途確認）。
- 漏れの式 n̄_leak は ε_s ≪ 1 の 1 次まで。ε_s が大きい（< 15 dB）ときは段ごとに伝搬させる（実装はそうする）。
- ガウスの重なり（§3.6）は placeholder。チャープ、非ガウス形、多モードの squeezed 光（パルスの時間モードが複数）は T03。
- κ_pm、V_π、駆動雑音のスペクトルは素子が決まるまでパラメータ。
- Ref D（CW、MZI）には振り分け素子がないので、このトピックのブロックは Ref D と Ref A/B/C ではオフ。

## 5. 効く周波数帯

| 効果 | 帯域 | 備考 |
|---|---|---|
| E47 挿入損失 | 静的 | 温度依存（例: NanoSpeed の IL 温度依存 0.25 typ / 0.5 max dB [mfr]）は遅いドリフト |
| E48 消光比 | 静的〜バイアスのドリフトの帯域（< 1 Hz、[AE]） | ターゲットのモードには損失として働く |
| E49 駆動の位相雑音 | 駆動回路の帯域（DC – GHz）だが、ポート k ではフレーム周期 T_F ごとにサンプルされる → [0, 1/(2T_F)] に折り返す。f_rep = 6 MHz、K = 32、N_b = 1 なら 1/(2T_F) = 94 kHz | 速い部分はパルスごとに独立な白色の位相誤差（ロック不能）。遅い部分（バイアス、温度）はロックで補正できるのは、スイッチがループの中にあるときだけ |
| E50 立ち上がり・同期 | 静的（条件を満たすかどうか）＋タイミングジッタ σ_j（パルスごと） | |
| E51 δL_k、T_rep のドリフト | 静的＋温度ドリフト（< 1 Hz） | |

## 6. スケーリング

- **ファイバ長**: L_k ∝ k·N_b/f_rep。ファイバ損失 1 − η_f,k ≈ (ln 10/10)·α_f·L_k、熱による時間ずれ κ_T L_k ΔT、環境位相雑音（T08）は L_k とともに増える。
  f_rep を下げる（PC が使える）ほど L_k は長くなる。
- **パワー**: squeezed 光自体は弱い。スイッチの効果はパワーに依存しない（光屈折効果は強い光でのみ）。漏れ光子数 n̄_leak ∝ sinh² r（ポンプパワーとともに増える）。
- **温度**: LN のバイアスドリフト（焦電効果）、PC 結晶の自然複屈折の温度変化（2 個の結晶を対にして補償するのが一般的 [UNVERIFIED]）、ファイバの群遅延 κ_T L_k。
- **パルス幅**: 時間ずれの許容値 δt_max ∝ τ_p（§3.6）。結晶の GDD による広がりは 1/τ_p² で効く（短いパルスほど厳しい）。ガード条件（§3.4）は τ_p にほぼ依存しない。
- **ポート数 K**: 二分木の段数 S = log2 K、η_sw = Π t_s。L_max ∝ (K − 1)。積分型測定の薄まり 1/K。
- **繰り返し f_rep**: 1 段目のトグル f_rep/(2N_b) ≤ f_tog,max、エッジ条件 t_r ≲ N_b T_rep、T_rep のドリフトの効き k N_b T_rep δT/T。

## 7. 概算（Ref D と将来の 32 ポート構成）

**Ref D**: 振り分け素子はない（D02）。このトピックのブロックはオフ。以下は将来の構成（D07–D13、D30–D32）。値の代用:
n_g = 1.468（survey の値、T06 で確認、placeholder）、α_f = 0.2 dB/km（Ref C の SMF の値。PM ファイバの値は T09、placeholder）、
κ_T = 39 ps/(km·K)（survey [AE]、T08 で確認、placeholder）、スイッチ入力での squeezing −10 dB（D16 の範囲の上端）。

### 7.1 遅延線の長さ（N_b = 1、L_0 = 0）

| f_rep | T_rep | 隣のポートとの差 c T_rep/n_g | L_31 | ファイバ透過 η_f,31（0.2 dB/km） | 群遅延の熱係数 κ_T L_31 |
|---|---|---|---|---|---|
| 1 MHz | 1000 ns | 204 m | 6.33 km | 0.747 | 247 ps/K |
| 3 MHz | 333 ns | 68.1 m | 2.11 km | 0.907 | 82 ps/K |
| 6 MHz | 167 ns | 34.0 m | 1.06 km | 0.953 | 41 ps/K |
| 20 MHz | 50 ns | 10.2 m | 317 m | 0.986 | 12 ps/K |
| 100 MHz | 10 ns | 2.04 m | 63.3 m | 0.997 | 2.5 ps/K |
| 1 GHz | 1 ns | 0.204 m | 6.33 m | 0.9997 | 0.25 ps/K |

ポートごとに損失が違う（0 から 1 − η_f,31 まで）ので、GBS 入力の squeezing はポートごとに変わる。

### 7.2 二分木（S = 5）の挿入損失と squeezing

| 1 段の損失 | η_sw（5 段） | −10 dB 入力の出力 | −6 dB 入力の出力 | 上限 10 log10(1 − η) |
|---|---|---|---|---|
| LN MZM 4 dB [mfr 上限] | 0.010（−20 dB） | −0.04 dB | −0.03 dB | −0.04 dB |
| LN MZM 3 dB [mfr オプション] | 0.032 | −0.13 dB | −0.10 dB | −0.14 dB |
| LN MZM 2 dB [mfr カスタム] | 0.100 | −0.41 dB | −0.34 dB | −0.46 dB |
| NanoSpeed 0.6 dB [mfr typ] | 0.50 | −2.60 dB | −2.04 dB | −3.0 dB |
| PC 段 0.96 [AE: Madsen 2022 の 4 段 15% から] | 0.815 | −5.75 dB | −4.10 dB | −7.3 dB |
| PC 段 0.99 [ref: Wang 2017 PC >99%、PBS・結合は除く] | 0.951 | −8.41 dB | −5.41 dB | −13.1 dB |
| PC 段 0.992 [ref: Kaneda 2019 PC 単体] | 0.961 | −8.68 dB | −5.52 dB | −14.1 dB |

→ ファイバ型 LN スイッチの木では squeezing はほぼ残らない。実用的な候補は自由空間の PC＋PBS の木（1 段 1–4%）。PC＋PBS では PBS の消光比（>1200:1、>2000:1 [ref: Wang 2017, 2019]）、
ファイバへの再結合（Wang 2017 は効率の主な制限を「ファイバの結合と伝搬損失」としている）が 1 段の損失に入る。

### 7.3 速さの制約

- PC（f_tog,max ≈ 0.5–3 MHz、t_r ≈ 8–10 ns）でパルスごと（N_b = 1）: f_rep ≤ 2 f_tog,max ≈ 1–6 MHz。6 MHz なら L_31 = 1.06 km（Borealis と同じ f_rep と 1 段目の 3 MHz）。
- PC で f_rep = 100 MHz を保つ: N_b ≥ f_rep/(2 f_tog,max) = 50（f_tog,max = 1 MHz）→ L_31 = 31 × 50 × 2.04 m = 3.17 km、η_f,31 = 0.864。
- LN MZM（t_r ~ 数十 ps）: f_rep 1 GHz まで可能だが損失で不可（§7.2）。

### 7.4 タイミングの許容値（ガウス、η_mm ≥ 0.99、δt ≤ 0.085 τ_p、placeholder）

| τ_p | δt_max | δL_max（n_g = 1.468） | 許される温度差（L_31 = 63 m、κ_T = 39 ps/(km·K)） |
|---|---|---|---|
| 100 fs | 8.5 fs | 1.7 µm | 3.4 mK |
| 1 ps | 85 fs | 17 µm | 34 mK |
| 10 ps | 0.85 ps | 0.17 mm | 0.34 K |

L_31 が 1.06 km（6 MHz）なら許される温度差は 1/17 になる。ただしこの温度差は「ポート間（GBS）」または「squeezed 側と LO_k 側の差」に対するもので、同じ環境に並べたファイバの共通の温度変化は打ち消す（打ち消しの程度は D09 の相関係数）。

T_rep のドリフト（§3.4）: δT/T = 10⁻⁶ で δt_31^{rep} = 310 fs（100 MHz）、5.2 ps（6 MHz）。τ_p = 100 fs で η_mm ≥ 0.99 には δT/T ≤ 2.7×10⁻⁸（100 MHz）。
CW レーザーを変調してパルスを作る方式（Liu 2026 のポンプ、Madsen 2022）なら T_rep は電子回路のクロックで決まり安定。モード同期レーザーなら f_rep の安定化が要る。

### 7.5 駆動の位相雑音とバイアス

- 5 段の合計で σ_sw ≤ 10 mrad（Ref A の θ = 9 mrad と同程度）にするには、κ_pm = 1 で δV/V_π ≤ 10 mrad/(π√5) = 1.4×10⁻³ rms（各段、独立）。
  PC の V_π が kV 級なら 1 V rms 級、LN の V_π = 4 V なら 5.7 mV rms。
- バイアスの誤差 δΓ_b と漏れ: 0.05 rad → 6.2×10⁻⁴（−32 dB）、0.1 rad → 2.5×10⁻³（−26 dB）、0.2 rad → 1.0×10⁻²（−20 dB）。200:1（23 dB）を保つには |δΓ_b| ≤ 0.14 rad。

### 7.6 漏れ光子（GBS へ、§3.3）

ε_s = 1/200（5 段）、r = 1（sinh² r = 1.38）で n̄_leak ≈ 0.035（1 ポート・1 フレーム、ターゲットでないビン）。
時間分解できる検出器（ジッタ ≪ T_rep）ならゲートで除ける。積分型（TES など）では GBS の誤差になる。許容量は GBS_REQUIREMENTS（Open）しだい。

### 7.7 積分型測定の薄まり（§3.7）

K = 32、η = 1、−10 dB: R_meas = (0.1 + 31)/32 = 0.972 → −0.12 dB。パルスごとに分解すれば −10 dB。
LO_k を全パルス列で配るなら BHD の帯域は f_rep より十分高く、パルスごとに積分して読む必要がある（f_rep = 6 MHz なら 35–100 MHz 級の BHD で可能、100 MHz なら GHz 級。T13）。

## 8. survey・DECISIONS への提案（調整役が判断）

1. survey E48: 「損失と同時にモード間の相関も生む」→「ターゲットの時間モードには損失として働く。相関（クリックの相関）は検出がビンを分けられない場合だけ」（§3.3）。
2. survey E47: 「ファイバ型 LN 変調器 3–5 dB」→「1×2/2×2 LN スイッチ <4 dB（<3 dB オプション、約 2 dB カスタム）[EOSPACE 2024]。5 段の木で 10–20 dB」。「自由空間ポッケルスセルは小」→「PC 0.8–1%、PBS・結合込みで 1 段 約 1–4%」。
3. D08 の L_k ≈ k·c·T_rep/n_g は N_b = 1（パルスごとの振り分け）と「ポート k は最後から k 番目のスロット」という番号付けを前提にしている。一般形 L_k = L_0 + k·N_b·c·T_rep/n_g を設定のパラメータにすることを提案（N_b = 1 が既定）。
4. survey E51 の帯域「静的」に「T_rep のドリフト（k N_b T_rep δT/T）」を加える。

## 9. 検証テスト（実装が通すべきもの）

1. **純損失**（prompt §7-3）: 各ポートで V_out = η_k V_in + (1 − η_k)·1/2 が数値誤差（< 1e-12 相対）で成り立つ。
   二分木で t_s = 0.99、ε_s = 0.005（5 段）なら η_sw = (0.99 × 0.995)⁵ = 0.92745、V_in = 0.05（−10 dB）で V_out = 0.082647（−7.817 dB）。
2. **漏れはターゲットのモードと直交**: K = 4 と K = 32 の時間ビンのモードの共分散行列で、有限の ER（10–30 dB）のスイッチ網を通したあと、
   (a) ターゲットのビン同士の共分散（ポート j ≠ k）が厳密に 0、(b) ポート k のターゲットの分散が「η_sw,k の純損失」と一致（差 < 1e-12）、
   (c) ターゲットでないビンの平均光子数の合計が Σ_k (1 − Π(1−ε_s)) sinh² r に一致（1 次近似の n̄_leak との差は O(ε²)）。
3. **薄まり**: K = 32、η = 1、−10 dB 入力、LO が全パルス列、積分型の測定で −0.124 dB、パルスごとの分解で −10.00 dB。
4. **スクイーズ角の符号**: スイッチの位相 ψ_sw,k = 0.3 rad を入れると、分散が最小になる LO の位相が θ_sq + 0.3 rad に移る（conventions §3–4）。
5. **駆動の位相雑音**: 5 段、κ_pm = 1、各段 δV/V_π = 1.4235×10⁻³ rms の白色ガウス雑音（パルスごとに独立）で σ_sw = 10.0 mrad。
   R− = 0.1、R+ = 100（入力）で、測定値が conventions §6 の厳密式 R′ = (R− + R+)/2 + (R− − R+)/2·e^{−2σ²} = 0.10999（−9.587 dB）に一致（統計誤差内）。
6. **バイアスのドリフト**: δΓ_b = 0.1 rad で、その段の透過が cos²(0.05) = 0.997502 倍、漏れが 2.4979×10⁻³（−26.02 dB）。
7. **遅延線の設計**: f_rep = 100 MHz、n_g = 1.468、N_b = 1、L_0 = 0 で L_1 = 2.0422 m、L_31 = 63.31 m。N_b = 50 で L_31 = 3165 m。
   全ポートのターゲットのパルスが同じ時刻に着く（数値誤差内）。T_rep を δT/T = 10⁻⁶ 変えると、ポート 31 と 0 の時間ずれが 310 fs。
8. **重なり（placeholder）**: δt = 0.0849 τ_p で η_mm = 0.990。T03 の式に置き換えたら、ガウスの極限でこの値に戻ること。
9. **立ち上がり**: 誤差関数形のエッジ（t_r = 10 ns）で、パルスの中心をエッジの中央に置くと漏れ 0.5、エッジから 2 t_r 以上離すと漏れ < 1e-6（バイアス 0 のとき）。
   §3.4 の条件を満たさない設定（例: PC、f_tog,max = 1 MHz、f_rep = 100 MHz、N_b = 1）はオフライン段で「実現不能」の警告を出す。
10. **オフで恒等**: 振り分けブロックをオフにすると、全出力が入力と一致（Ref A/B/C/D の結果を変えない）。

## 10. シミュレータへの要求

### 10.1 ブロック

| ブロック | 段 | 入力 | 出力 |
|---|---|---|---|
| `demux_design`（オフライン表） | O | K、f_rep（または T_rep）、N_b、L_0、ガードスロット数、n_g、α_f、κ_T、素子のプリセット（段数 S・配置（木/直列）、t_s、ER_s、V_π,s、κ_pm,s、t_r、f_tog,max、t_ring、段の GDD） | ポートごとの L_k、τ_k、η_sw,k、η_f,k、静的な ψ_sw,k^{(0)}、漏れ行列 ε_{k′j}、GDD_k、実現可能性の判定（§3.4）、n̄_leak |
| `demux_switch`（時間ループ） | T | ポート k、`demux_design` の表、駆動雑音の PSD（または rms とスペクトル形）、バイアスドリフトのモデル（ランダムウォーク・温度係数）、温度の時系列（T08 と共有） | 光路の区間としての位相 ψ_sw,k(t)（遅い部分）、パルスごとの速い位相の分散 σ_sw,fast²（T04 の σ² に加える）、η_sw,k(t)（バイアスで変わる） |
| `delay_line_k` | O ＋ T | L_k、n_g、κ_T、δL_k、温度 | 群遅延 τ_k(t)、squeezed 側と LO_k 側の時間ずれ δt_k(t) → T03 の重なり表 η_mm(δt) を引く |
| 測定の模擬 | M | `gating: pulse_resolved | integrating`、K、N_b | §3.7 の薄まりを含む分散 |

- `demux_switch` は squeezed 光の光路の区間（conventions §9: d = +1、偏光の軸はスイッチの出口の偏光（PBS の出口で H/V が変わるので波長板を含めて設定））として置き、損失を squeezed 状態にかける（prompt §5.5）。
- 区間の属性 `shared_with_probe`（D30 の逆向き probe がこの区間を通るか）を持たせる。通らない（干渉点がスイッチより下流）場合、ψ_sw,k の遅い部分もロックで補正されない。
  通る場合は、probe もスイッチの状態で切られる（ポート k が開いている時間だけ通る）ので、probe のパルスのタイミングをスイッチの窓に合わせる必要がある（T10/T14 に渡す）。
- LO_k を同じ種類のスイッチで振り分ける設定（D09 の Open）も同じブロックで表す。LO 側では損失は LO のパワーにだけ効き、位相雑音は Δ_k にだけ入る（Θ_k には入らない）。

### 10.2 オン・オフ（toggle flags）

`demux.enable`、`demux.insertion_loss`、`demux.crosstalk`、`demux.drive_phase_noise`、`demux.bias_drift`、`demux.rise_time_clipping`、
`demux.delay_thermal`、`demux.rep_rate_drift`、`demux.static_mismatch`（δL_k）。既定: Ref A/B/C/D では `demux.enable = false`、将来の構成では全部オン。

### 10.3 デバイスのプリセット（設定ファイル、値は §2 の出典付き）

`pockels_pbs_tree`（t_s = 0.96–0.99、ER 20–23 dB、t_r = 8–10 ns、f_tog,max = 0.5–3 MHz）、`ln_mzm_fiber_tree`（IL 2–4 dB、ER 20–60 dB、V_π 4–5 V、t_r ~ 数十 ps、帯域 10–65 GHz）、
`crystal_fiber_switch`（NanoSpeed: IL 0.6 dB、ER 18–25 dB、t_r ≤ 300 ns、≤ 300 kHz）。κ_pm と駆動雑音はすべて未定のパラメータ。

## 11. 文献

一次資料（本文を読んで値を確認したもの）:

- H.-L. Liu et al., "Gaussian boson sampling with 1,024 squeezed states in 8,176 modes," Nature 653, 687 (2026). doi:10.1038/s41586-026-10523-6（`refs/` の PDF、Methods を読んだ）。
- L. S. Madsen et al., "Quantum computational advantage with a programmable photonic processor," Nature 606, 75–81 (2022). doi:10.1038/s41586-022-04725-x
  （PMC9159949 の本文から: 1→16 の振り分け、15 個の EOM・4 段、3/1.5/0.75/0.375 MHz、消光比 >200:1、振り分けの損失 約 15%、6 MHz、ループの位相雑音 0.02/0.03/0.15 rad。
  巻・ページは書誌データベースによる。振り分けの値は PMC の本文の抽出（WebFetch）で確認し、原 PDF の表とは照合していない）。
- H. Wang et al., "High-efficiency multiphoton boson sampling," Nat. Photonics 11, 361–365 (2017). doi:10.1038/nphoton.2017.63
  （St Andrews の著者最終稿を読んだ: PC 透過 >99%、消光比 >100:1、PBS >1200:1、0.76 MHz、立ち上がり 8 ns、平均効率 約 84.5%、76 MHz の単一光子列）。
- H. Wang et al., "Boson sampling with 20 input photons and a 60-mode interferometer in a 10^14-dimensional Hilbert space," Phys. Rev. Lett. 123, 250503 (2019). doi:10.1103/PhysRevLett.123.250503
  （arXiv:1910.09930 を読んだ: 19 組の PC と PBS、PBS >2000:1、PC >100:1、平均効率 約 83%、木の形、パルス列を 20 の segment に分ける）。
- S. Takeda, K. Takase, A. Furusawa, "On-demand photonic entanglement synthesizer," Sci. Adv. 5, eaaw4530 (2019). doi:10.1126/sciadv.aaw4530
  （PMC6524974: RTP のバルク EOM、立ち上がり 約 10 ns、時間ビン 66 ns、ループ 19.8 m、往復損失 約 7%、位相のドリフト 約 7°）。
- F. Kaneda, P. G. Kwiat, "High-efficiency single-photon generation via large-scale active time multiplexing," Sci. Adv. 5, eaaw8586 (2019). doi:10.1126/sciadv.aaw8586
  （PMC6777972: RTP 対の PC、PC の損失 0.8%、1 周 1.2%、500 kHz（ドライバ律速）、「市販のバルク PC は 1 MHz」）。

書誌のみ確認（本文未読、値は使っていない。ドリフトの機構の出典）:

- J. P. Salvestrini, L. Guilbert, M. Fontana, M. Abarkan, S. Gille, "Analysis and control of the DC drift in LiNbO3-based Mach–Zehnder modulators," J. Lightwave Technol. 29(10), 1522–1534 (2011). doi:10.1109/JLT.2011.2136322
- E. L. Wooten et al., "A review of lithium niobate modulators for fiber-optic communications systems," IEEE J. Sel. Top. Quantum Electron. 6(1), 69–82 (2000). doi:10.1109/2944.826874

メーカー資料（D24、閲覧日 2026-10-28）:

- EOSPACE, "Product Brochure 2024": https://static1.squarespace.com/static/5b4391a0a2772c7a5b0dae43/t/65fc777677b0500bb240bc2e/1711044471997/EOSPACE_Product_Brochure_2024.pdf
  （1×2/2×2 スイッチの個別の表はなく、同じ LN 素子群の値: 損失 <4 dB・<3 dB・約 2 dB、V_π 4–5 V @ 1 GHz、帯域 10–65 GHz、高消光比版 >30–60 dB）。
- Agiltron, "NanoSpeed 1x2 Series Fiber Optical Switch" datasheet（Rev 06/12/23、Laser Components 経由）:
  https://www.lasercomponents.com/fileadmin/user_upload/home/Datasheets/diverse-fiber-optics/passive-components/ns-switch-1x2.pdf

UNVERIFIED（記憶または概算で、原典で確認していないもの）: PC 結晶の κ_pm の値（§3.5）、PC 結晶の対による温度補償、PC・LN の GDD の値、AOM・TFLN スイッチの値、
κ_T = 39 ps/(km·K) と n_g = 1.468（survey [AE]、T06/T08 で確認）、α_f = 0.2 dB/km を PM ファイバに使うこと（Ref C は SMF）。

## 12. ユーザの決定が必要な点（推測しない）

1. **振り分け素子の種類と (f_rep, N_b) の領域**（D08 Open）: ファイバ型 LN の木は損失で不可（§7.2）。PC の木を使うなら f_rep ≲ 6 MHz（N_b = 1、L_31 ≈ 1 km 以上）か、
   ブロック振り分け（N_b > 1、km 級の遅延）か。f_rep 自体（D07 Open）とあわせて決める必要がある。
2. **ブロック振り分け（N_b > 1）を設定として認めるか**: D08 の L_k の式を一般化することになる（§8-3）。
3. **GBS の検出器が時間ビンを分けられるか**（ゲート、ジッタ）と、ターゲットでないビンの漏れ光子の許容量（GBS_REQUIREMENTS Open）。これで ER の要求が「損失」だけで決まるか、「迷光子」で決まるかが変わる。
4. **D30 の逆向き probe がスイッチを通るか**（干渉点がスイッチの上流か下流か）。通らないなら駆動の遅い位相雑音・バイアスドリフトもロックの外に残り、指標 (ii) の予算に直接入る。
5. **LO_k の配り方**（D09 Open）: 全パルス列を分岐（BHD はパルスごとに分解が必要、T13）か、LO も振り分ける（LO 側のスイッチの位相雑音が Δ_k に入る）か。
6. **D07 のパルス幅の下限（100 fs）**: L_31 = 63 m でも mK 級の温度差、f_rep の相対安定度 ~10⁻⁸ が要る（§7.4）。掃引範囲をそのまま残すか、時間ずれの能動補正（遅延のロック）をブロックとして加えるか。
