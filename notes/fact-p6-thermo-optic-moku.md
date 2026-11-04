# Phase 6 の事実調査: 熱光学係数の波長依存性と Moku:Delta の機能（2026-11-04）

担当の調査（`anthropic-skills:research-papers` を使用）の要約。[AE] は担当の導出・計算。

## 1. 異なる波長の代理信号の偏り係数 R（`lit-svarc2023.md` §5.1 の定義）

- SMF-28、1550 nm、293–295 K: **R ≈ +0.3 ps/(km·K)、不確かさ約 ±0.6（−0.3〜+0.9）。符号は決まらない**。R/κ_τ ≈ +0.008。支配的な不確かさは 1.55 µm での d(dn/dT)/dλ。[AE]
- 現行のプリセットの k_T = 44、κ_τ = 39 から出る R ≈ 2.3 ps/(km·K) は、ガラスのデータが許す範囲の外。k_T と κ_τ は別のファイバ・被覆の値で、被覆のひずみの項（n·ξ·ε）は R では打ち消すが各係数では支配的。独立に掃引すると R が約 10 倍広くなる。
- ストレッチャのロックの偏り: R だけで −1.9×10⁻³ rad/(K·m·THz)（−5.7〜+1.9×10⁻³）。ひずみ光学係数 ξ の分散を入れると R′ = R − (β + nα)(ω∂ξ/∂ω)/(cξ) ≈ +0.56（0.0〜+1.2）ps/(km·K)、−3.5×10⁻³ rad/(K·m·THz)。p11・p12 の分散は未確認。[AE]
- 上書きの比 g_r/g_s = 1 + (Δν/ν)(n_g/n + ω∂lnξ/∂ω) ≈ 1 + 1.007·Δν/ν。Pittaluga 型の 2π リセットをストレッチャで行うと 2πΔν/ν0 より約 0.7 % 大きい。[AE]
- 出典: Leviton & Frey, arXiv:0805.0091（Proc. SPIE 6273、番号・DOI は未確認）Table 3・6; Rego 2023, Sensors 23, 6023, doi:10.3390/s23136023（Table 1・2）; Rego 2024, Sensors 24, 4857, doi:10.3390/s24154857（Table 3）; Slavík et al. 2015, Sci. Rep. 5, 15447, doi:10.1038/srep15447（Fig. 3・4: 37.4 ps/(km·K) at 1535 nm、1510–1635 nm で同じ）; Bousonville et al. BIW2012 MOPG033、DIPAC09 TUPB35（計算値）。ファイバでの熱位相係数の波長依存性の直接測定は見つからない。
- 推奨: R（または R′）を独立の軸に、名目 +0.3（+0.56）ps/(km·K)、範囲 −0.5〜+1.5。k_T・κ_τ を独立に掃引しない。

## 2. Moku:Delta（`refs/manuals/MokuDelta-Specifications.pdf` v26-0817、Liquid Instruments の API 文書）

- 復調位相: Lock-in Amplifier（p48–49）と Laser Lock Box（p43–44）は 0–360°、分解能 10⁻⁶°、API の `set_demodulation(mode, frequency, phase)` で実行中に設定できる。ロック中の変更で乱れないかは記載なし。PID Controller には復調がない（p64–65）。
- ロック中の励起: FRA の In/In1 モード（p30）を Multi-instrument Mode（p79–80）で PID と並べられる。MiM の配線は加算できない（各行き先に 1 つの元）ので、PID Controller の 4 入力の制御行列（係数 −20〜+20、p64–65）で誤差と FRA の出力を足す。API 文書では 2 チャンネル・dB 表記、PDF では線形で、要確認。Laser Lock Box には注入点の記載がない。
- 変換器: 入力 8 ch、ADC 14 bit（広帯域）／20 bit（低帯域）、ENOB 10.8、5 GSa/s（p6）。出力 8 ch、DAC 14 bit、5 GSa/s（p7）。PID Controller は 4 入力・4 出力・4 PID、「最大帯域 400 kHz で位相遅れ 30°」（p64–65、[AE] 約 0.2 µs の遅延）。Laser Lock Box の入出力の遅延 約 186.9 ns（Liquid Instruments の knowledge base）。
