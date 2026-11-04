# P7 Kashiwazaki et al. 2023: ロックしたままの開ループ応答の測定と、周波数を最適化したロック（lit-kashiwazaki2023）

- Phase: 6（`docs/phase6_prompt.md` の Procedure 1）、日付: 2026-11-04
- 入力: `refs/2026-11-04_extra_refs/234003_1_5.0144385.pdf`（7 ページ。1 ページ目は出版社の表紙、本文は 234003-1〜-6。`pdftotext -layout` と、Fig. 2–4 のページを 130–300 dpi の画像にして読んだ）、
  `DECISIONS.md`（D01–D83）、`notes/conventions.md`、`notes/00_summary.md`、`notes/survey.md`、`notes/error-signal-generation.md`（T14）、`notes/actuators-control-digital.md`（T16）、
  `notes/opa-waveguide-squeezing.md`（T01）、`notes/phase-noise-to-squeezing.md`（T04）、`presets/reference/refA.yaml`・`refA_tap.yaml`、`presets/future/F_common.yaml`、`src/cranberry/`
- 値の印: **[K23 p.N]** は本文のページ 234003-N、**[K23 Fig.]** は図から読んだ値（読み取りの誤差つき）、**[Ref A]** は既存のノート経由の Ref A の値、
  **[calc]** はこのノートで式を数値で評価した値、**[AE]** は概算・解釈（未検証）、**[code]** は `src/cranberry/` で確認した事実、**UNVERIFIED** は原典で確認できなかったもの。
- D34 に従い、ここで挙げる構成はすべて切り替えて試せる案で、既存のプリセットは変えない。

## 0. 要点

1. **8.3 ± 0.1 dB の squeezing（11 MHz、損失・回路雑音の補正なし）**、ポンプ 660 mW、LO 15 mW、実効損失 12 %、当てはめた位相の揺れ **0.8°（= 14.0 mrad）** [K23 p.4–5]。
   論文の式 (1)(2) にこの値を入れると R′− = −8.35 dB、R′+ = +19.66 dB [calc]（測定と一致）。
2. **Ref A の「前の研究」の破線（14 mrad、12 %、D48 §7-2）はこの論文の値** [AE: 0.8° = 13.96 mrad と 12 % が一致]。ただしこの論文の α = 820 %/W で 1200 mW の点を計算すると −6.72 dB [calc] で、
   Ref A の図の破線の右端（約 −4.8 dB、T01 §7.3）とはさらに離れる。D48 §7-2 の食い違いはこの論文の値では説明できない。
3. **ロックのトポロジーは Ref A の「従来の tap 方式」（`refA_tap`）と同じ形**: 99/1 の鏡で squeezed 光と probe の 1 % を APD に取り、**2 MHz で復調**（probe が OPA で差周波を作るのでビートは AOM のずれ 1 MHz の 2 倍、[K23 p.4]）。
   これは T14 §3.2 で [AE] としていた「2 f_sh で復調」を、この論文の tap 方式について**原典で確認**したことになる（Ref A 自身の phase-detection OPA の復調周波数は UNVERIFIED のまま）。probe–LO のロックは BHD の信号を 1 MHz で復調。
4. **ロックしたままの応答の測定**: 速いループの「PID ＋ LPF」の箱（Fig. 2(a) の A または B）をネットワークアナライザ（Keysight E5061B）の箱 C に置き換え、**ミキサの出力と EOM の間**に入れる。
   その間、**同じループの 2 つ目の PID（遅い I）とファイバストレッチャがロックを保つ**。これで遅いループより高い周波数の応答が測れる [K23 p.4]。
   位相の遅れが −180° に達するのは squeezed–probe のループで約 4 MHz、probe–LO のループで約 2 MHz、利得は 1 MHz より下でほぼ平ら → AOM のずれを 1 MHz に決めた [K23 p.4]。
   **論文は「位相余裕」という言葉を使っていない。** 書かれている基準は「周波数は高いほどよいが、ロックの系が発振しない範囲で」だけ。
5. **アクチュエータの組**: probe と LO の両方の経路に **EOM（速い、PI）＋ ファイバストレッチャ（遅い、I、広いレンジ）**。遅い PID の入力は速い PID の出力（積み重ね、T16 §4.4、Ref C と同じ）[K23 p.3、Fig. 2(a)]。
   **EOM・ストレッチャの型番・帯域・レンジ、PID の利得、ループの交差周波数、LPF の遮断周波数は本文にない**（UNVERIFIED）。
6. **シミュレータ** [code]: `sine_injection` ＋ `SimController.measure_transfer_function` で、注入点を選べば開ループ L(f) も植物 P(f) も測れる。ただし
   (a) 周波数ごとに初期状態から走らせ直す（ロックを保ったまま掃引しない）、(b) 1 回の実行で 1 つの正弦、(c) `engage`/`disengage` は全ループ一括、(d) `acquire` は注入を必ずオフにする、
   ので、P7・Moku:Delta と**同じ手順**にはならない。新しいブロック 2 つ（N7-1 `excitation`、N7-2 `loop_identification`）とインターフェースへの後方互換の追加（N7-3）を提案する（§5.3）。

## 1. 論文がしたこと

| 段落（APL は節番号なし） | ページ・図 | 内容 |
|---|---|---|
| Abstract・導入 | p.1–2 | 前作（Ref 16 = Kashiwazaki et al. APL 119, 251104 (2021)）の 6 THz 帯域・6.3 dB の PPLN 導波路 OPA を使い、測定系の改良で CW 8.3 dB。改良は 3 つ: (1) 余分な光学部品を除き、全光学の直交位相測定から低損失の BHD に替えて実効損失を約 12 % に、(2) 系の周波数応答を測って「周波数を最適化した位相ロック」を作り、(3) 測定周波数とポンプパワーを依存性から最適化 |
| 前作の構成 | p.2、Fig. 1(a) | ファイバ・ピグテールの OPA モジュール（45 mm の PPLN 導波路、SHG 820 %/W、導波路の損失 7 %）。全光学の測定で実効損失 21 %（squeezer の導波路 約 4 %、squeezer モジュールの出力側 8 %、増幅 OPA の入力側 8 %、増幅 OPA の導波路 4 %）[K23 p.2–3] |
| 光路 | p.3、Fig. 2(a) | NKT Photonics BASIK のファイバレーザー 1545.32 nm → EDFA（Keyopsys CEFA-C-PB-HP）→ ポンプ・probe・LO の 3 本に分け、各経路に 1 nm の BPF（ASE 除去）。ポンプは NTT Electronics WH-0772-000-F-B-C（PPLN 導波路）で倍波にして OPA へ。probe は AOM 2 台（Chongqing Smart Science & Technology SGTF-40-1550-1P、駆動 +40 MHz + f と −40 MHz、Fig. 2(a) の表記）で 1 MHz ずらして OPA に入れ、squeezed 光と同じ向きに進む。OPA の後は自由空間。DM1 で残りのポンプを PD1 へ。99/1 の鏡で 1 % を APD（Thorlabs APD430C、PD2）へ、99 % を BHD（Laser Components IGHQEX0100-1550-10-1.0-SPAR-TH-40、QE 約 98 %）へ。ESA は Keysight N9010B。LO と probe の経路にファイバストレッチャ（FS）と EOM |
| 電気の経路 | p.3–4、Fig. 2(a) | PD2 の信号をミキサ ＋ LPF で 2 MHz から下げ、PID に入れる（probe と squeezed 光のロック）。PID の出力は EOM と、もう 1 つの PID へ。2 つ目の PID は FS への遅い信号を作る。probe–LO のロックは BHD の信号を 1 MHz から下げて、別の 2 つの PID で同様に。Fig. 2(a) の箱は「PID (PI)」（速い）と「PID (I)」（遅い）、その前に LPF |
| 系の応答の測定 | p.4、Fig. 2(b)(c) | §0-4 のとおり。測定の参考文献は Ref 26（D. Hou et al., Opt. Express 24, 17242 (2016)、未読） |
| squeezing の測定 | p.4、Fig. 3 | ESA のゼロスパン、中心 11 MHz、RBW 1 MHz、VBW 1 kHz。LO 15 mW（BHD の直前）。ポンプ 660 mW は PD1 の値から導波路の出口から PD1 までの損失を見込んで推定（ポンプの減衰なしを仮定）。位相をロックした squeezed 雑音はショット雑音より 8.3 ± 0.1 dB 低い |
| 測定周波数の選び方 | p.4–5、Fig. 4(a) | ショット雑音と回路雑音のクリアランスが 11 MHz 付近で約 25 dB（squeezed 光に対して約 0.3 % の損失に相当）。高い周波数では回路雑音が増えて squeezing が下がるので 11 MHz |
| ポンプパワーの選び方 | p.4–5、Fig. 4(b) | squeezing はロックして、反 squeezing は LO の位相を走査して測る。式 (1) R′± ≈ R± cos²θ̃ + R∓ sin²θ̃、式 (2) R± = (1 − η) + η exp(±2√(aP))、a = 820 %/W、損失 12 %（検出 8 % ＋ 回路雑音 0.3 % ＋ 導波路の実効 4 %）。θ̃ = 0.8° の実線が測定に合い、θ̃ = 0 の破線より高いポンプで squeezing が下がる。高いポンプでの劣化は主にロックした位相の揺れによると著者は考える。約 660 mW で squeezing が最大 |
| 結論 | p.5 | 検出の 8 % を差し引くと導波路の直後で 10 dB 超と推定。導波路の損失が次の制限 |

注:
- 式 (2) の本文の説明は「g は squeezed 光の実効的な透過率」だが、値には「12 %」（損失）を入れている [K23 p.5]。式の形 (1 − η) + η e^{±2r} からは η = 0.88（透過率）、損失 12 % と読む。
- Fig. 2(a) のポンプの波長の表記は「775 nm」。1545.32/2 = 772.66 nm なので図の表記は丸めた値と読む [AE]。

## 2. 数値

| 項目 | 値 | 出典 |
|---|---|---|
| squeezing（補正なし） | −8.3 ± 0.1 dB @ 11 MHz | [K23 p.4] |
| 反 squeezing | 約 +19〜+20 dB（Fig. 3 の位相走査の山の高さ、±1 dB の読み取り） | [K23 Fig. 3] |
| モデルの再計算（a = 8.2 /W、P = 0.66 W、損失 12 %、θ̃ = 0.8°） | R′− = −8.35 dB、R′+ = +19.66 dB | [calc] |
| 同じ損失・θ̃ での最適ポンプ | 556 mW で −8.40 dB（a = 8.2 /W）。a = 9.06 /W なら 503 mW（T04 V2 と同じ） | [calc] |
| 残留位相の揺れ | θ̃ = 0.8° = 13.96 mrad（「at most 0.8°」、Fig. 4(b) の当てはめ。直接の測定ではない） | [K23 p.4–5] |
| 実効損失 | 12 % = 検出 8 %（PD の QE 約 98 % → 2 %、probe–LO の可視度 98.5 % → 約 3 %、squeezed 光の経路 約 3 %（1 % の tap を含む））＋ 回路雑音 0.3 % ＋ 導波路 4 % | [K23 p.3、p.5] |
| クリアランス | 約 25 dB @ 11 MHz | [K23 p.4] |
| LO | 15 mW | [K23 p.4] |
| ポンプ | 660 mW（導波路の中、PD1 から推定） | [K23 p.4] |
| SHG 効率 | 820 %/W（前作 Ref 16 の値をそのまま使用） | [K23 p.2、p.5] |
| AOM のずれ | 1 MHz（AOM 2 台の差）。squeezed–probe のビートは 2 MHz（OPA の差周波発生）、probe–LO のビートは 1 MHz | [K23 p.3–4] |
| 復調（ダウンコンバート）の周波数 | squeezed–probe: 2 MHz、probe–LO: 1 MHz | [K23 p.3–4] |
| tap | 1 %（99/1 の鏡）→ APD430C | [K23 p.3] |
| 応答の測定 | ネットワークアナライザ E5061B、注入はミキサの出力と EOM の間（速い PID ＋ LPF の箱を置き換え）、遅い PID ＋ FS でロックを保持 | [K23 p.4] |
| 植物の応答（Fig. 2(b)、squeezed–probe） | 利得 約 −25 dB で 1 kHz から 1 MHz までほぼ平ら（±2 dB の読み取り）、−180° は約 4 MHz、2 MHz にビートの鋭い山 | [K23 p.4、Fig. 2(b)] |
| 植物の応答（Fig. 2(c)、probe–LO） | 利得 約 −25 dB で平ら、−180° は約 2 MHz、1 MHz にビートの山。4–5 MHz で約 −450° まで回る（読み取り） | [K23 p.4、Fig. 2(c)] |
| 等価な純遅延（−180° の周波数を純遅延とみなした値） | 約 125 ns（4 MHz）、約 250 ns（2 MHz） | [AE] |
| 測定した周波数の範囲 | 1 kHz – 約 5 MHz（図の横軸） | [K23 Fig. 2(b)(c)] |
| EOM・FS の型番・V_π・帯域・レンジ | 記載なし | UNVERIFIED |
| PID の利得・交差周波数・LPF の遮断周波数 | 記載なし | UNVERIFIED |
| ロックの時間（どれだけ続くか） | 記載なし | UNVERIFIED |

観察 [AE]: どちらのループでも、復調の周波数（2 MHz と 1 MHz）が植物の −180° の周波数（4 MHz と 2 MHz）の**ちょうど半分**。論文はこの関係を書いていないので、設計の規則なのか偶然なのかは分からない。
ミキサの出力に残る 2 倍の周波数の成分（2 × 復調周波数 = −180° の周波数）が、LPF の後も残るとそこで発振しやすい、という読み方はできるが、確認できない。

## 3. Ref A（Hirota et al. 2026、`presets/reference/refA.yaml`）との違い

| 項目 | P7（Kashiwazaki 2023） | Ref A（既存のノート・プリセットから） |
|---|---|---|
| squeezing | −8.3 dB @ 11 MHz | −10 dB @ 3 MHz（`refA.yaml` の ESA 中心 3 MHz） |
| 実効損失 | 12 % | 8 %（D48） |
| 残留位相 | 0.8° = 14 mrad | 9 mrad（T01・D63） |
| ポンプと α | 660 mW（導波路の中、PD1 から推定）、a = 820 %/W（SHG の測定値） | 640 mW（squeezer の出口）、α = 906 %/W（当てはめ、T01 §0-2） |
| **lock 1（ポンプ–probe の位相）** | **1 % の tap ＋ APD**（squeezed 光に 1 % の損失）。2 MHz で復調 | **別の phase-detection OPA**（tap の損失なし、probe を数 mW にできる）。復調周波数は本文になく、T14 は 2 f_sh と推定（UNVERIFIED） |
| lock 2（probe–LO） | BHD の 1 MHz ビート | 同じ（BHD、f_sh = 1 MHz） |
| アクチュエータ | 各経路に EOM（速い）＋ FS（遅い、広いレンジ）の積み重ね | 各経路に EOM（iXblue MPX-LN-0.1、T16 §3.2）。`refA.yaml` では FS は LO の位相走査（Fig. 3）だけ |
| ロックの設計の方法 | ネットワークアナライザで植物の応答を測り、AOM のずれ（= 復調周波数）を決めた | 記載なし（T14 は Ref A の lock 1 の SNR の議論のみ） |
| BHD | IGHQEX0100（QE 約 98 %）、クリアランス 約 25 dB @ 11 MHz、LO 15 mW | IGHQEX0500（QE 約 99 %）、クリアランス 28 dB @ 3 MHz、LO 15.4 mW（`refA.yaml`） |
| 光源 | NKT BASIK ＋ EDFA ＋ 1 nm の BPF | CW seed ＋ EDFA ＋ SHG（D48）、LO にフィルタキャビティ（`refA.yaml`、D63-6 で Δν_FC を掃引） |

要するに、P7 は **Ref A の一つ前の世代で、`refA_tap` の構成に「遅い FS の段」と「測った植物の応答から決めた復調周波数」を足したもの** [AE]。
Ref A が 14 → 9 mrad に下げた主な変更（Ref A は lock 1 の SNR と説明、T14 §3.4 はそれに疑問）を比べる基準点として使える。

## 4. 自分たちの系に移せること

### 4.1 将来の構成（単一の squeezer → 振り分け → ポートごとの PM ファイバ遅延 → GBS、Liu 型の逆向き probe が基準）

1. **手順そのもの**（いちばん移しやすい）[AE]: 遅い段（TEC、F_common の `pid0s` → `tec0`）でロックを保ったまま、速い段（PZT ストレッチャ）の植物の応答 P(f) を測る。
   PZT の共振（`pzt_f_res_hz` 5–60 kHz・`pzt_q`、D51-11）と純遅延（`latency_s`）は資料から決められない量で、F の掃引軸の中で最も効く（D80-3: L-Liu の掃引の多くの点が線形に不安定）。
   **測った P(f) で `f_u_hz`・`pzt_f_res_hz`・`latency_s` の掃引範囲を狭め、BoTorch の上下限にする**のが Phase 6 の目的（phase6_prompt の P7 の行）に合う。
2. **アクチュエータの組**: P7 の「EOM（速い）＋ FS（遅い）」はそのままは移せない。EOM は squeezed 光の経路に置けない（挿入損失 2.5–3.5 dB、T16 §3.2）ので、
   ポートのロック（D10: probe と squeezed 光が一緒に通る区間にアクチュエータ）では「PZT（速い）＋ TEC（遅い）」が対応する組。EOM ＋ FS の組が使えるのは LO の配送のロック（`L-LOdeliv`）や参照光の側。
3. **復調周波数の選び方**: P7 では MHz のループ（EOM）なので、ビートの周波数がループの帯域に近く、植物の位相との関係で選ぶ必要があった。
   自分たちの Liu 型（PZT、f_u = 1–30 kHz）では、ヘテロダインの f_sh（`detect_heterodyne` は 80 MHz [AE]）は帯域よりはるかに高く、制約は「復調の LPF ≫ f_u、f_sh ≫ LPF」（`detect_heterodyne.yaml` の注記）だけ。
   P7 の基準が効くのは、次のベンチ実験の EOM のループの方。
4. **残差の水準**: P7 の 14 mrad（12 % の損失で −8.3 dB）と Ref A の 9 mrad は、D78 の線幅の表の 10 mrad の列と同じ水準。将来の構成の (ii) の σ は名目で約 100 mrad（D74・D80）なので、まだ 1 桁離れている。

### 4.2 次のベンチ実験（レーザー ＋ 光路 ＋ BHD、GBS なし。Ref A・Ref C のような位相の安定な squeezing を狙う）

P7 は**このベンチ実験のいちばん直接のひな形**になる [AE]:
- 構成: probe を AOM で 1 MHz ずらし、OPA の後で 1 % を取って 2 MHz で復調（lock 1）、BHD のビートを 1 MHz で復調（lock 2）。各経路に EOM（速い PI）＋ FS（遅い I）。
- 手順: (1) 遅いループだけでロック、(2) 速いループの箱を周波数応答の測定に置き換えて植物の応答を測る、(3) −180° の周波数から復調周波数・LPF・PI の利得を決める、
  (4) 測定周波数をクリアランスで選ぶ（P7 は 25 dB @ 11 MHz）、(5) θ̃ を当てはめてポンプパワーを最適化（T04 V2 の閉じた式。P7 の値では 556 mW で −8.40 dB [calc]）。
- 期待値の目安: 損失 12 %・14 mrad で −8.3 dB、損失 8 %・9 mrad（Ref A）で −10.1 dB（T01）。どちらに近づけるかは損失と θ̃ の予算の配分で、`phase-noise-to-squeezing.md` の表で読める。
- Ref C 型（dither ＋ 積み重ね 3 段）とは、ロックの誤差信号の作り方が違うだけで、「積み重ねの遅い段で保ちながら速い段を測る」手順は共通に使える。
- Moku:Delta で同じことをするには、PID の前後に励起を足して 2 つの信号を同時に取る機能が要る。Moku の Frequency Response Analyzer や Multi-Instrument でそれができるかは、`notes/actuators-control-digital.md` §4.5 の表にない（UNVERIFIED、資料で確認）。

## 5. シミュレータへの対応

### 5.1 既存のブロック・プリセットで表せること [code]

| 要素 | 既存のもの |
|---|---|
| P7 の光学系とロック | `presets/reference/refA_tap.yaml`（`parametric_beat_detector` heterodyne → 2 MHz の `demodulator` → `pid` → `eom_phase`、`tap_ratio: 0.01` [AE]、lock 2 は `bhd_probe_beat` → 1 MHz の `dem2` → `pid2` → `eom2`） |
| P7 の値 | オーバーレイで: `opa_waveguide` の `alpha_per_w: 8.2`・`pump_w: 0.66`、損失 12 %、ESA の `f_s_hz: 1.1e7`・クリアランス 25 dB・LO 15 mW・QE 0.98 |
| 遅い FS の段 | `stack`（F_common の `lock0.stack`、Ref C の 3 段と同じ）。`pid1` → `eom1` を `!replace` で `stack: {stages: [pid1, pid1s], actuators: [eom1, pzt_fiber_stretcher]}` に。PR 区間に `stretcher` の駆動を足す（区間のリストは `id` でマージ） |
| 速い段だけ開く | 各 `pid` のトグル `feedback`（`RunSpec.toggles` の `"pid1.feedback": false`） |
| 正弦の注入 | `sine_injection`（`in` ＋ A sin(2πft) → `out`、実行時パラメータ `freq_hz`・`amplitude`、トグル `enabled`）。オーバーレイの `!replace` で `dac1` を `seq: [sine_injection, dac]` にすれば、PID の出力とアクチュエータの入力の間に入る |
| 掃引正弦の伝達関数 | `SimController.configure({"injections": {pt: {block, reference, response, amplitude, settle_s, periods}}})` → `measure_transfer_function(pt, freqs)` = 応答/参照（1 ビンのロックイン） |
| モデルの L(f) | `loop_stability_check`（積み重ねの段、極、共振、遅延、DSC1 の離散の式。`max_l_at_180`・`f_180_hz` を出す）、`pid.freq_response`・`pzt_fiber_stretcher.freq_response`・フィルタ・遅延の `freq_response` |
| 記録 | `measure.record.slots` に書いた任意のスロット（`decim` で間引き） |

**新しいブロックなしで今できる同定** [code、AE]:
- **(B) 閉ループのままの L(f)（Moku でも同じにできる形）**: 誤差の加算点に注入する。`sine_injection` を `in: v.adc0, out: v.inj0` で入れ、`pid0f` の `in` を `v.inj0` に。
  参照 = `v.inj0`（e + d）、応答 = `v.adc0`（e）とすると e = −L(e + d) なので **H = 応答/参照 = −L(f)**。積み重ねの全段（`pid0f` と、その出力を読む `pid0s`）が同じ誤差から駆動されるので、測れるのは**全体の L**。
  アクチュエータの入力の前（`ctrl.sup0.u` の後）に入れると、遅い段は速い段の出力を読むので、測れるのは L_f/(1 + L_s′) の形（遅い段の帯域より上で −L_f に一致）[AE]。
- **(A) P7 そのまま（速い段を開き、遅い段で保持）**: `pid0f.feedback: false`、注入をアクチュエータの入力に、応答を `v.adc0` に。ただし F_common の遅い段は**速い段の出力を読む積み重ね**なので、速い段を開くと遅い段の入力もなくなる。
  P7 のように遅い段だけでロックを保つには、測定のときだけ `pid0s` の `in` を誤差（`v.adc0`）につなぎ替え、利得も付け直す別の設定（オーバーレイ）が要る。
  P7 の論文も、測定中に遅い PID の入力をどこから取ったかは書いていない（Fig. 2(a) では遅い PID は速い PID の出力から LPF 経由。箱 C に置き換えた後の配線は不明、UNVERIFIED）。

### 5.2 今の仕組みで足りないこと [code]

1. **ロックを保ったまま掃引しない。** `measure_transfer_function` は周波数ごとに `settle_s + periods/f` の実行を**初期状態から**行う（`SimController` は状態を持ち越さない）。
   各周波数でロックの引き込みからやり直し、遅い段（F_common の TEC は `rate_div` 10000、mHz 帯）は `settle_s` の間に定常にならない。周波数ごとに遅いドリフトの状態も違う。P7 と実機（Moku）は 1 回のロックの中で掃引する。
2. **1 回の実行で 1 つの正弦。** マルチサイン・チャープ・PRBS・段階的な掃引がない。`sine_injection` を直列に並べれば複数の正弦は足せるが、位相は全部 0 から（波高率を下げられない）、周波数は実行中に変わらず、`SimController` が動かすのは 1 つの注入点だけ。
3. **`engage`/`disengage` は全ループ一括**（`_toggle_state` が全 `pid_block` の `feedback` を同じ値にする）。速い段だけ開く、は `RunSpec.toggles` の書き換え（オーバーレイ）でしかできず、インターフェースの手順にならない。
4. **`acquire` は注入を必ずオフにする**（`_toggle_state` で `enabled: False`）。励起をかけたままの in-loop の時系列（d、e、u）をインターフェースから取れない。`measure_transfer_function` は比しか返さず、コヒーレンスも出さない。
5. **時系列から伝達関数を推定する道具がない。** `core/fft.psd_onesided`（Welch）はあるが、相互スペクトル（H1 = S_xy/S_xx）・コヒーレンス・ロックインの段ごとの平均は `src/` にない（grep で `csd`・chirp・同定のコードは見つからなかった）。
6. ランナーにチャンクの間のフック・パラメータの予定表がない（`runner/run.py`）。

### 5.3 新しいブロックの提案

**N7-1 `excitation`（CTRL、時間ループ）**
- 役割: 実行中に周波数が変わる励起を、任意のスロットに足す。
- 入力: `in`（null 可）。
- 出力: `out` = in + x(t)、`ctrl.<name>.x`（x(t) そのもの。推定の道具変数）、`ctrl.<name>.f_hz`（段階的な掃引のときの今の周波数）、`ctrl.<name>.step`（段の番号）。
- 構造パラメータ: `waveform: stepped_sine | multisine | chirp | prbs | white`、`freqs_hz`（表、または `f_start_hz`・`f_stop_hz`・`n_freqs`・`spacing: log | lin`）、`dwell_s`・`settle_s`（段ごと）、
  `phases: schroeder | zero | random`（`random` は乱数の子 `SeedSequence`、G11）、`amplitude_shape: flat | table`、`unit`。
- 実行時パラメータ: `amplitude`、`t_start_s`（ロックの引き込みの後に始める）。
- トグル: `enabled`。
- オーバーレイでは足りない理由: 5.2 の 1・2。`sine_injection` は 1 回の実行で周波数が固定された 1 つの正弦で、直列に並べても位相・波高率・段階的な掃引を作れない。

**N7-2 `loop_identification`（MEAS、ESA と同じ測定段のフック。記録した系列からループの後で計算）**
- 入力（スロット名）: `x`（`ctrl.<exc>.x`）、`ref`（注入の後の信号、例 `v.inj0`）、`resp`（注入の前の信号、例 `v.adc0`）、任意で `stage_outputs`（`ctrl.pid0f.u`・`ctrl.pid0s.u` など、段ごとの寄与）、`step`（N7-1 の段の番号）。
- 方法: `lockin_per_step`（段ごとに整数周期で復調、`settle_s` を捨てる）| `welch_h1`（x を道具変数にした H = S_x,resp / S_x,ref、`nperseg`、窓）。
- 出力（`RunResult.measurements`）: H(f)、コヒーレンス、`relation: minus_l | plant | sensitivity` に応じた L̂(f)・S = 1/(1 + L)・T、交差周波数 f_u、位相余裕、利得余裕、f_180、
  モデル（同じ実行の `loop_stability_check` の L(f)）との比。
- オーバーレイでは足りない理由: 5.2 の 4・5。既存の測定のフックは ESA・パルスの標本・指標の分解だけで、時系列から多周波数の伝達関数とコヒーレンスを出すものがない。

**N7-3 抽象制御器インターフェースへの後方互換の追加（ブロックではない。D55 で凍結した API なので調整役の判断）**
- `engage(loop_id=None)`・`disengage(loop_id=None)`（None は今と同じ全ループ）。段ごと（`<loop>.<stage>`）も指定できるように。
- `acquire(duration_s, signals, injections_on=())`: 指定した注入（N7-1 を含む）をオンにしたまま記録する。
- `measure_loop_gain(point, excitation, relation)`: N7-1 ＋ N7-2 を 1 回の `acquire` で行い、H(f)・コヒーレンス・余裕を返す。`SimController` は 1 回の実行、実機（Moku）は機器の周波数応答の測定か、出力の加算 ＋ 2 チャンネルの記録に写す。
- 理由: phase6_prompt の「同じ手順を Moku:Delta で走らせる」には、手順がインターフェースの呼び出しだけで書けている必要がある。今は 5.2 の 3・4 のため、オーバーレイ ＋ ランナーの直接呼び出しが要る。

### 5.4 検証と使い方の案（Procedure 3–5）

1. **検証（解析解のある系）**: `refD_nominal`（雑音オフ）で N7-1 ＋ N7-2 の (B) が −L(f) を、`tests/integration/refD_support.sensitivity_z` の厳密な z 領域の式と 1e−3・0.1° 程度で再現すること（VT-13 の (a) と同じ基準）。
   同じ系で今の `measure_transfer_function`（周波数ごとの実行）とも一致すること。
2. **ロックを保ったままの効果**: F_liu20 の名目で、(A) と (B) を TEC の段あり・なしで比べ、遅い段の帯域より上で L_f に一致すること、下で S_s の分だけずれることを示す。
3. **同定したモデルの使い方**: `pzt_f_res_hz`・`pzt_q`・`latency_s` を同定の結果に当てはめ、BoTorch の `f_u_hz`・`f_12_hz` の上下限を「f_180 から位相余裕 45° 以上」「利得余裕 6 dB 以上」で決める（基準の値は [AE]、Phase 6 の決定事項）。
   `loop_stability_check` が「不安定」と判定した点（D80-3）が同定でも不安定に見えるかも確かめる。
4. **P7 の再現**: `refA_tap` ＋ P7 の値 ＋ 遅い FS の段（5.1）で R′− ≈ −8.3 dB になる θ̃ の内訳（どの雑音源が 14 mrad を作るか）。Ref A の 9 mrad との差（D63-7 の「9 mrad はフィルタキャビティを通ったレーザー雑音で決まる」）と並べる。

## 6. UNVERIFIED と未解決の問い

- EOM・FS の型番・V_π・帯域・レンジ、PID の構造と利得、ループの交差周波数、LPF の遮断周波数、ロックの持続時間: 本文にない。
- 箱 C に置き換えたとき、遅い PID の入力をどこから取ったか（Fig. 2(a) の配線のままでは遅い PID は速い PID の出力から駆動される）。ネットワークアナライザのどちらのポートが注入で、どちらが参照か。
- Fig. 2(b)(c) の利得の軸（約 −25 dB）が何の比か（V/V と思われるが、EOM のドライバ・ミキサの変換利得を含む植物全体か）。
- 「復調周波数 = −180° の周波数の半分」（§2 の観察）が設計の規則か偶然か。
- 0.8° は Fig. 4(b) の当てはめで、ロックの残差だけでなく位相の雑音すべてを含む。どのループの何が 14 mrad を作るかは書かれていない。
- Ref 26（Hou et al. 2016、ロック中のループ利得と帯域の直接測定）は未読。in-loop の注入の式（5.1 の (B)）は一般的な形で書いたもので、Ref 26 の方法と同じかは確かめていない。
- D48 §7-2（Ref A の「前の研究」の破線の右端）: P7 の値（a = 8.2 /W、12 %、14 mrad）では 1200 mW で −6.72 dB で、図の約 −4.8 dB とさらに離れる。Ref A が破線を描いたときの α・P の定義は未確認。
- Moku:Delta で、PID の内側に励起を足しながら 2 点を同時に記録できるか（`notes/actuators-control-digital.md` §4.5 の資料の範囲外）。

## 7. 文献

- T. Kashiwazaki, T. Yamashima, K. Enbutsu, T. Kazama, A. Inoue, K. Fukui, M. Endo, T. Umeki, A. Furusawa,
  "Over-8-dB squeezed light generation by a broadband waveguide optical parametric amplifier toward fault-tolerant ultra-fast quantum computers,"
  Appl. Phys. Lett. **122**, 234003 (2023). DOI 10.1063/5.0144385（PDF の表紙と各ページの脚注で確認）。本文 p.1–5、Fig. 1–4、式 (1)(2) を使用。
- [Ref 16 of P7] T. Kashiwazaki, T. Yamashima, N. Takanashi, A. Inoue, T. Umeki, A. Furusawa, Appl. Phys. Lett. **119**, 251104 (2021)（P7 の文献表から。本文は未読）。
- [Ref 26 of P7] D. Hou, J. Tian, F. Sun, X. Huang, "Direct measurement of loop gain and bandwidth of phase-locked loop for mode-locked laser," Opt. Express **24**, 17242 (2016)（P7 の文献表から。**未読**）。
- Ref A: K. Hirota, T. Kashiwazaki et al., Opt. Express **34**, 7958 (2026), DOI 10.1364/OE.585323（D51-1。値は既存のノートとプリセットから）。
