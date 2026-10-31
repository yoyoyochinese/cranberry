# Phase 2 統合とトリアージ（00_summary）

- Phase: 2（統合とトリアージ、1 エージェント）／ 日付: 2026-10-31
- 入力: `docs/squeezing_sim_prompt.md`（§4 Phase 2、§5、§7）、`AGENTS.md`、`DECISIONS.md`（D01–D50、最上位の方針は **D34**）、`notes/conventions.md`、
  `notes/survey.md`（この Phase で改訂）、deep dive 17 本（T01 `opa-waveguide-squeezing` … T17 `pulsed-opo-source`）、調整役の Phase 1 ダイジェスト。
  D37 で追加された `refs/manuals/pnh1550r5a2.pdf` と Opneti 偏光子の仕様書は、この Phase で値を読んだ（§5 の C21）。
- 規約: 本文は日本語、記号は英語（D19）。符号・規格化は `notes/conventions.md` に従う。値の印は deep dive と同じ（[AE] は概算、UNVERIFIED は原典未確認）。
  このノートで新しく出した概算には **[AE-P2]** と書く。
- D34 に従い、このノートの「プリセット」はすべて同格の名前付き設定で、「既定」ではない。悪い構成も比較用のプリセットとして残す。
  効果はひとつも削除しない。すべての構成で無視できる効果は optional（ブロックは実装し、参照プリセットではオフ）にする（prompt §4 Phase 2）。

## 0. 凡例

### 0.1 構成の列（効果の大きさを書く列）

| 列 | 構成 | 大きさを見積もった代表値（ノートの値をそのまま引用。長さはセルに併記） |
|---|---|---|
| A | Ref A（PPLN 導波路 OPA、PD-OPA lock 1 ＋ BHD lock 2）。A-tap は比較用の conventional tap | CW、P = 0.64 W（出口）、L = 8%、θ = 9 mrad |
| B | Ref B（CW bow-tie OPO、同軸 probe、dither） | x = √(P/P_th)、f_s = 1 MHz |
| C | Ref C（CW モノリシック OPO、1 km × 2 本、量子雑音ロック） | 1.75 GHz FWHM、損失 21% |
| D | Ref D（DX1 → PM MZI ＋ 915B、PDB450C、DSC1 の DC side-fringe lock） | 3.33 mW、ΔL 名目 9 m（0–15 m 掃引） |
| Fb | 将来・Borealis 型パルス励起 OPO（D43・D50） | 66 MHz FWHM、T_p = 3 ns、f_rep ≤ 45 MHz（参照点 6 MHz、L_31 ≈ 1.06 km） |
| Fl | 将来・Liu 2026 型パルス励起 OPO（D50） | 125 MHz FWHM、T_p = 1.6 ns、f_rep 20 MHz（L_31 = 317 m（N_b = 1）または 1.27 km（N_b = 4）） |
| Fw | 将来・Ref A 型の単一パス導波路 OPA（D07・D44） | 45 mm PPLN（代用）、τ_p = 100 fs–10 ps（参照点 1 ps）、f_rep 100 MHz（L_31 = 63 m）。630 m は f_rep = 10 MHz の値 |

将来の 3 列で probe・ロック・経路を決めていない値は、ノートの仮定（逆向き probe 1 mW、同じ偏光・同じ周波数など）で書き、そう明記した。どの将来の構成も 32 ポート（D08）。

### 0.2 判定

- **kept**: 少なくとも 1 つの参照プリセット（§2.1 の A・A-tap・B・C・D、§2.2 の将来のプリセット）でオンになる効果。
- **opt**（optional）: すべての参照プリセットで無視できる、またはデータがなくオフにする効果。ブロックは実装し、比較用・アブレーションの設定ではオンにできる（D34、prompt §4）。
- どちらの場合も、オン・オフは設定のフラグで切り替える（保護リストの「雑音源のトグル」）。

### 0.3 段（prompt §5.1、survey §0.2）

G = GNLSE の項、Q = 量子雑音層、O = オフラインのモデル・表、T = 時間ループのブロック、M = 測定の模擬・出力。

### 0.4 ブロック族の略号（prompt §5.5 ＋ DECISIONS で増えた族）

| 略号 | 族 | 出典 |
|---|---|---|
| SRC | 光源（OPA、OPO、パルス OPO、レーザー、ポンプ） | prompt §5.5 |
| PATH | 光路・区間・ビームの位相の組み立て（Jones 行列、区間 × ビームの係数、遅延バッファ） | prompt §5.2、D11、D47 |
| NOISE | 雑音源（レーザー、環境、ファイバ散乱、RF） | prompt §5.5 |
| MOD | 変調器（EOM、dither、AOM、probe のパルス化） | prompt §5.5 |
| ERR | 誤差信号の生成（#1–#6） | prompt §5.5、D06、D30 |
| CTRL | 制御器・ADC/DAC・ロック監視 | prompt §5.5 |
| ACT | アクチュエータ（PZT ストレッチャ、HV アンプ、PZT ミラー、EOM、AOM） | prompt §5.5 |
| TEC | 大レンジ・熱アクチュエータ（TEC、スプールの温度制御、2π リセット） | D41 |
| DET | 検出器（BHD、平衡 PD、パルス BHD） | prompt §5.5 |
| DMX | 時間の振り分け（Pockels 木・ブロック・受動ループ・ファイバ LN、遅延線） | D39 |
| LOSS | 損失収支・モード整合・tap | T12 |
| FSR | OPO の FSR 櫛のフィルタ（FBG の組、MZI の段） | D50 |
| SPF | GBS 前のスペクトルフィルタ、ポートごとの分散補償・遅延の微調整 | D44 |
| LOS | LO の時間波形・配送 | D45、D50 |
| FIB | ファイバのオフライン計算（GNLSE、量子雑音層、散乱表） | prompt §5.1 |
| MEAS | 測定の模擬と指標（ESA、位相 → squeezing、指標 (i)(ii)、迷光） | prompt §6、D12、D46 |
| CHK | オフラインの整合性チェック（SBS、繰り返しの上限、振り分けの実現可能性） | T07、T11、T17 |

### 0.5 集計

- 効果は 103 件（survey の E01–E84 と、deep dive で見つかった N01–N19）。**kept 94 件、opt 9 件**（E02、E17、E24、E28、E30、E33、N06、N18、N19）。
- opt の 9 件のうち、E17（ASE）はデータがないためのオフ（D49）で、他の 8 件はすべての構成で無視できるためのオフ。

## 1. 効果のトリアージ表

各セルはノートの値の要約（記号と数値は出典のノート）。「—」はその構成に該当しない。長さ・パワーの仮定はセルに併記した。

### 1.1 光源・スクイーザー（E01–E11）

| ID | 効果 | 判定 | 段 | 族 | A | B | C | D | Fb | Fl | Fw | 出典 | 理由 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E01 | 導波路 OPA の分布損失 | kept | O | SRC | L_wg 2.03%、同じ α で R′+ −0.27 dB | — | — | — | — | — | 同左（45 mm 代用）、L_wg(P) | T01 | Ref A 再現。§7-1 は集中損失、分布は切り替え（D48） |
| E02 | ポンプ枯渇 | opt | O | SRC | 生成 2×10⁻⁵ W / 0.64 W、PD-OPA の傾きのみ | 無視 | 無視 | — | ≲ 3 光子/パルス | 同左 | 無視 | T01、T17 | すべてで無視。結合波方程式は V-T01-8 用に実装 |
| E03 | ポンプ強度雑音 | kept | O＋T | SRC | R′+ 16.3 dB/W（+1% で +0.10 dB） | 3.4·δP/P | 同形 | — | d ln r/d ln P ≈ 0.51、AOM で安定化 | 同左（補正なし） | A と同じ（ピーク） | T01、T02、T17 | R′+ の揺れ、ロックの傾き |
| E04 | ポンプ位相雑音 | kept | T | PATH | ½ψ_p、lock 1・2 で補正 | lock 1 | QNL が直接 | — | ポンプ位相ロックなし。(ii) には 4 sin²(πfΔτ_k) でのみ | 透過ポンプのホモダイン 1 nm | ポートごとに違う時刻でサンプル | T01、T10、T17 | 経路ごとに計算（D11） |
| E05 | 位相不整合の温度 | kept | O＋T | SRC | 10 mK で 2 mrad（lock 2 で共通） | P_th(T_c)（2 次） | 温度で離調と連動 | — | 二重共振を温度で保つ | 同左 | probe を通さないと補正されない | T01、T02、T17 | 遅いドリフトと角度の回転 |
| E06 | ポンプ誘起吸収 | kept | O＋T | SRC | 1 W まで ΔL ≲ 数% | L(x) ƒ、ドリフトで −9.44 → −8.61 dB | — | — | 値なし（掃引） | 同左 | L_wg(P) ƒ | T01、T02 | §7-9 に必須（callable） |
| E07 | 空間多モード | kept | O | LOSS | 4% = 0.98 × 0.98 | ξ² 0.976 | 15% に含む | — | 値なし | 値なし | 値なし | T01、T12 | 静的な損失 |
| E08 | OPO キャビティの伝達 | kept | O | SRC | — | FWHM 12.1 MHz、ρ 0.970、Ω_B = 2πf/γ_B | 1.75 GHz FWHM | — | 66 MHz、ρ ≈ 0.95 | 125 MHz、ρ ≈ 0.99 | — | T02、T17 | OPO の光源 |
| E09 | OPO の離調・長さ雑音 | kept | O＋T | SRC | — | −4.2°/nm | 約 1.7 rad/K（6 mK 必要） | — | −0.6Δ_s/κ_c、ロック光の Δ(dn/dT) で 0.8–8 rad/K | 同左（ロック光が同じモード） | — | T02、T17 | 区間 × ビームの係数（N02） |
| E10 | Ref B の G と x の定義 | kept | O | SRC | — | x = √(P/P_th)、G = 18.7 は使わない | — | — | — | — | — | T02、D48 | 検証目標の定義 |
| E11 | パルス励起の多モード（GVM） | kept | O | SRC | — | — | — | — | → N13 | → N13 | τ_GVM 13.7 ps、K 195–320、必要 P_peak 増 | T03 | 単一モード近似が破れる |

### 1.2 レーザー雑音と遅延差（E12–E19）

| ID | 効果 | 判定 | 段 | 族 | A | B | C | D | Fb | Fl | Fw | 出典 | 理由 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E12 | レーザー周波数雑音 | kept | T | NOISE | 白色 0.2 kHz（代用、UNVERIFIED） | Ti:S（値なし） | 値なし（掃引） | h0 = 8×10³（D49）–3.2×10⁴ Hz²/Hz、h₋₁ 0–5×10⁸ | 光源未定（要求は E14） | 同左 | 同左 | T05、D35、D49 | すべての構成の基本 |
| E13 | 遅延差の干渉 | kept | T | PATH | LO だけフィルタキャビティ（τ_eff = 1/(πΔν_FC)） | — | δL（Ref C、名目 0） | ΔL 9 m で 2.5×10⁻⁵ rad/√Hz、DC–100 kHz で 7.8 mrad | LO_k の δL_k はモード重なりで自動的に小さい | 同左 | δL_k ≤ 13 cm（CW LO の場合） | T05、T15 | 放出時刻の遅延バッファ |
| E14 | ポート間のコヒーレンス | kept | T | NOISE | — | — | — | — | 10 mrad に Δν ≤ 約 3 Hz [AE-P2]（6 MHz） | ≤ 10 Hz（20 MHz） | ≤ 51 Hz（100 MHz）、DX1 では 0.22 rad | T04、T05、T10 | 帯域外は線幅への要求、帯域内は Liu 型で打ち消し（D47） |
| E15 | RIN | kept | T | NOISE | CMRR 漏れ 1.9%（−145 dB/Hz、ε 0.01） | 同形 | 同形 | (S − V_off)/\|K_d\| 経由、4.5×10⁻¹⁰ rad/√Hz | probe 古典雑音 N01 | 同左 | SPM 経由で 0.3·σ_RIN（1 ps） | T05、T13、T15、T06 | 小さいが経路ごとに必要 |
| E16 | 過剰雑音ピーク | kept | T | NOISE | RIN 0.7 MHz ピーク（Koheras 参考） | 89・177 kHz | — | — | — | — | — | T05 | Ref B の 295 kHz の理由 |
| E17 | ASE | opt | O | NOISE | ε_ASE 10⁻³ で −10.58 → −8.41 dB、値なし | — | — | EDFA なし | EDFA の使い方未定 | 同左 | 同左 | T05、T13、D49 | 文献値がなく、A プリセットではオフ（D49）。掃引用 |
| E18 | フィルタキャビティ・MCC の残差 | kept | T | PATH | Δν_FC 1 MHz・δ 1 kHz で 2 mrad（掃引） | MCC（線幅なし） | — | — | — | — | — | T05 | A の LO 経路の要素（prompt §1） |
| E19 | タイミングジッタ | kept | T＋O | NOISE | — | — | — | — | 許容 436 ps（無視） | 227 ps（無視） | σ_t ≤ 0.142 T0（1 ps で 85 fs） | T03、T17 | Fw で必須 |

### 1.3 ファイバの線形伝搬と環境（E20–E28）

| ID | 効果 | 判定 | 段 | 族 | A | B | C | D | Fb | Fl | Fw | 出典 | 理由 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E20 | 伝搬損失 | kept | G/O | LOSS | — | — | 1 km で 4.5% | 数 m、無視 | 1.06 km: PM 11.5%、SMF 4.3% [AE-P2] | 317 m: PM 3.6%、SMF 1.3% | 63 m: PM 0.73%、SMF 0.26% | T06、T12 | ポートごとに損失が違う |
| E21 | 色分散 | kept | G/O | SPF | — | — | 量子サイドバンドの回転 4×10⁻⁴ rad（無視） | — | L_D ~ 10⁸ m（無視） | 同左 | 100 fs で L_D 0.15 m、1 ps のポート間重なり 0.47 | T03、T06、T17 | Fw で必須 |
| E22 | 温度位相ドリフト | kept | T | NOISE | — | — | 1 km で 4.4×10⁴ rad/K（差で効く） | ΔL 9 m で 約 390 rad/K（レンジは 0.13 K） | 4.7×10⁴ rad/K（1.06 km）[AE-P2] | 1.4×10⁴ rad/K（317 m） | 2.8×10³ rad/K（63 m） | T08 | ロック外れの主因 |
| E23 | 音響・振動 | kept | T | NOISE | — | — | 差の h ≈ 1.5 rad²Hz、f_u 1 kHz で 49 mrad | 9 m spool で 0.1 Hz 以上 34–60 mrad（ロックなし） | 630 m で 89–160 mrad（ロックなし） | 同程度 | 63 m で 28–50 mrad、f_u 1 kHz で 2 mrad | T08 | ループで抑える外乱 |
| E24 | 熱力学的な基底 | opt | T | NOISE | — | — | 0.1 mrad | 10 µrad | 85 µrad（630 m） | 同 | 27 µrad（63 m） | T08、D49 | すべてで無視 |
| E25 | ファイバ間の相関 | kept | T | NOISE | — | — | ρ(f) 不明（掃引、Ref C は定性的） | 2 腕 | squeezed と LO_k、2 ポート組 | 同左 | 同左 | T08 | Ref C と D13 の 2 ポート組 |
| E26 | 群遅延の熱ドリフト | kept | T | PATH | — | — | — | — | 41 ps/K vs 436 ps（約 10 K の余裕） | 同様 | 2.5 ps/K（1 ps で 35 mK、ロック後 21–28% 残る） | T03、T08、T17 | Fw で必須 |
| E27 | ストレッチャ巻線（曲げ） | kept | O | ACT | — | — | — | 8.4 m、PER −31〜−36 dB、損失 0.01–2 dB | squeezed 経路に入る（D10） | 同左 | 同左 | T09、T12、T16 | 将来は損失が直接効く |
| E28 | エタロン | opt | T | NOISE | — | — | — | 10⁻⁶ rad | 10⁻⁶–10⁻⁵ rad | 同 | 同 | T12、D49 | すべてで無視 |

### 1.4 ファイバの非線形と量子雑音（E29–E35）

| ID | 効果 | 判定 | 段 | 族 | A | B | C | D | Fb | Fl | Fw | 出典 | 理由 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E29 | SPM（LO・probe） | kept | G | FIB | — | — | 1.2 mrad/mW（定数、ロックで吸収） | 2×10⁻⁵ rad/10 m | 15 ns LO で 0.9 mrad（63 m） | 同程度 | 100 MHz/63 m/1 mW: 0.05–0.3 rad。10 MHz/633 m はソリトン | T06、T17 | Fw の LO 位相雑音・モード重なり |
| E30 | XPM | opt | G/T | FIB | — | — | — | — | 逆向き probe 0.16 mrad/mW（63 m） | 同 | 同 | T06、T10 | P̄_pr ≲ 100 mW で無視 |
| E31 | 誘導ラマン・SSFS | kept | G | FIB | — | — | 閾値 23 W、無関係 | 無関係 | 無関係 | 無関係 | 10 MHz/633 m で重心 −14〜−86 GHz | T06 | Fw の GNLSE の項 |
| E32 | 自発ラマン | kept | Q/O | FIB | — | — | 無視（別ファイバ） | — | 9×10⁻⁶（1 mW、63 m） | 同 | 同（モードによらない） | T07、T17 | 迷光表の一項（安価） |
| E33 | GAWBS | opt | Q/T | FIB | — | — | 0.4–0.8 mrad（別ファイバ） | 0.1 mrad（帯域外） | 低次ピークだけ帯域内 | 同 | squeezed 自身 0.1–0.2 mrad。同じ向きの明るい光で SN の 10 倍 | T07、T17 | 参照プリセットに同じ向きの明るい光がない。§7-6 用に実装 |
| E34 | SBS 閾値 | kept | O | CHK | — | — | 34 mW（1 km） | 遠い | ns パルスで抑えられない | 同 | 63 m 0.53 W、630 m 54 mW | T07、T17 | 警告のチェック |
| E35 | 自発ブリルアン後方散乱 | kept | Q/O | FIB | — | — | — | — | 帯域外（重み ≤ 10⁻⁵） | 同 | 4.5×10⁻³（1 ps、1 mW、63 m） | T07、T10、T17 | Fw の迷光 |

### 1.5 偏光と PM ファイバ（E36–E40）

| ID | 効果 | 判定 | 段 | 族 | A | B | C | D | Fb | Fl | Fw | 出典 | 理由 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E36 | PER とクロストーク | kept | O＋T | PATH | — | — | — | 偏光子 ER ≥ 26 dB、カプラ PER ≥ 25 dB [manual] | 部品 20–30 dB → ε_pol 2×10⁻³–10⁻² | 同 | 同 | T09、T10、T15 | 寄生縞（N05）と迷光（E43） |
| E37 | PM の複屈折の温度 | kept | T | PATH | — | — | — | 915B 巻線 −17 rad/K | 直交 probe なら 63 m で −129 rad/K | 同 | 同 | T09 | 直交偏光 probe の誤差 |
| E38 | SMF の偏波ドリフト | kept | T | PATH | — | — | PC で調整（値なし） | 偏光子の前の SMF | SMF 変種で D_pol 掃引 | 同 | 同 | T09 | D32 の比較 |
| E39 | 偏光子・PC | kept | T | PATH | — | — | — | δP/P = −Θ²/4 | — | — | — | T09、T15 | Ref D のパワードリフト |
| E40 | PDL・PMD | kept | O | PATH | — | — | — | — | 無視（ns） | 無視 | SMF・100 fs で 63 m 0.7–4.2%、630 m 6.6–32% | T09 | Fw の SMF 変種 |

### 1.6 probe と squeezed 光の相殺の崩れ（E41–E46）

| ID | 効果 | 判定 | 段 | 族 | A | B | C | D | Fb | Fl | Fw | 出典 | 理由 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E41 | 直交偏光 probe の複屈折ドリフト | kept | T | PATH | — | — | — | — | PM: −129 rad/K（63 m）、10 mrad に 77 µK | 同（Liu の方式） | 同 | T09、T10 | in-loop に見えない誤差 |
| E42 | 逆向き伝搬の時刻ずれ | kept | T | PATH | — | — | — | — | f₁ = 0.151c/(n_gL)（1.06 km で 約 29 kHz [AE-P2]） | 317 m で 約 97 kHz [AE-P2] | 63 m で 0.49 MHz、1 kHz で −53 dB | T08、T10 | 音響帯には効かないが長いポートで効く |
| E43 | PER による漏れ（ε_pol） | kept | O | NOISE | — | — | — | — | PM 直交 2×10⁻³–10⁻²、SMF 直交 1/3 | 同 | 同 | T10 | 迷光の係数 |
| E44 | Rayleigh 迷光 | kept | O＋T | NOISE | — | — | — | — | 1 モード 435 光子（1 mW、63 m）→ ≲ 1.3 nW | 226 → ≲ 2.6 nW | 10 ns あたり 3.0×10²、1 ps モード 0.03 | T10、T17 | GBS の迷光と squeezing の劣化 |
| E45 | 入れ出し素子の損失 | kept | O | LOSS | — | — | — | — | 自由空間 PBS 約 1%、CIR1550PM ≤ 0.9 dB | 同 | 同 | T10、T12、D37 | 10 dB 目標の損失収支 |
| E46 | 対策の候補 | kept | — | （プリセットの軸） | — | — | — | — | +7 GHz ＋ FBG が使える | 同 | WDM・時間ゲート・低パワー | T10、T17 | §2.3 の軸として持つ |

### 1.7 振り分けと遅延線（E47–E51）

| ID | 効果 | 判定 | 段 | 族 | A | B | C | D | Fb | Fl | Fw | 出典 | 理由 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E47 | 振り分けの挿入損失 | kept | O | DMX | — | — | — | — | PC 木 5 段 η 0.82–0.96 | 同（受動ループなら別） | 同。LN 木 10–20 dB は不可 | T11 | 損失収支の主要項 |
| E48 | 消光比 | kept | O | DMX | — | — | — | — | 損失として働く。裾は N11 | 同 | 損失、n̄_leak 0.035（r = 1） | T11、T17 | 迷光・損失 |
| E49 | 駆動位相雑音・バイアス | kept | T | DMX | — | — | — | — | 5 段 10 mrad に δV/V_π ≤ 1.4×10⁻³ | 同 | 同 | T11 | (ii) に直接入りうる |
| E50 | 立ち上がり・同期 | kept | O | DMX | — | — | — | — | PC で f_rep ≤ 6 MHz（N_b = 1） | N_b ≥ 4 | N_b ≥ 50 | T11、T17 | 実現可能性の判定 |
| E51 | L_k の誤差・T_rep ドリフト | kept | O＋T | DMX | — | — | — | — | 436 ps に δT/T ≤ 8×10⁻⁵ | 同様 | δT/T 10⁻⁶ で 310 fs | T11、T03、T17 | Fw で厳しい |

### 1.8 損失とモード整合（E52–E56）

| ID | 効果 | 判定 | 段 | 族 | A | B | C | D | Fb | Fl | Fw | 出典 | 理由 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E52 | 光学部品の損失 | kept | O | LOSS | 1% | ζ 0.99 | 15% ＋ 8% | カプラ過剰損失 < 0.3 dB × 2 [manual]、接続 | η_up ≈ 0.76（ρ × FSR フィルタ × 結合） | 0.92–0.94 | η_up 0.951 ＋ ポートの部品 | T12、T17 | 10 dB 目標ではポートの部品 ≤ 0.16 dB |
| E53 | 空間モードの不整合 | kept | O | LOSS | 4% | ξ² 0.976 | 含む | 可視度にだけ効く | ファイバ BHD では結合損失 | 同 | 同 | T12 | 静的 |
| E54 | 時間・スペクトルの不整合 | kept | O | SRC/LOS | — | — | — | — | → N14 | → N14 | Δt ≤ 0.142 T0、GDD、θ_off | T03 | Fw で必須 |
| E55 | squeezed 経路のアクチュエータ・変調器の損失 | kept | O | LOSS | — | — | EOM は LO 側 | — | EOM は置けない（上限 2.6–3.6 dB）、915B 0.12–0.4 dB ＋ 曲げ | 同 | 同 | T12、T16 | D10 の配置規則 |
| E56 | conventional tap | kept | O＋T | ERR | A-tap: NEP 1 pW/√Hz で最適 r_t 3×10⁻³、−9.94 dB | — | — | — | — | — | — | T12、T14 | §7-8 |

### 1.9 検出と測定の模擬（E57–E65）

| ID | 効果 | 判定 | 段 | 族 | A | B | C | D | Fb | Fl | Fw | 出典 | 理由 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E57 | QE | kept | O | DET | 0.99 | 0.998 | 0.97 | PDB450C 0.80（傾きにだけ） | 0.99（代用） | 同 | 同 | T13 | 損失 |
| E58 | 帯域・ディップ | kept | O＋M | DET | 35 MHz | — | 0.96 / 0.83 GHz ディップ | 150 MHz | 200–500 MHz（パルス分解） | 同 | GHz 級（100 MHz） | T13、T17 | クリアランスを通して効く |
| E59 | 回路雑音 | kept | O＋M | DET | 0.16% | 0.67% | 7–24% | SN/電子 +1.2 dB | 28 dB を P_LO で換算（代用） | 同 | 同 | T13 | 損失と等価 |
| E60 | CMRR | kept | T/O | DET | 30 dB で RIN −150 なら 1.5% | 同形 | 同形 | > 25 dB | 同 | 同 | 同 | T13 | LO 古典雑音の漏れ |
| E61 | 飽和 | kept | O/T | DET | — | — | — | DSC1 ±4 V が先（余裕 14–17%） | — | — | — | T13 | Ref D のクリップ |
| E62 | パルス BHD | kept | O | DET | — | — | — | — | ゲート ≥ 25 ns、c_1 | 同 | B ≳ f_rep | T13、T17 | 将来の BHD |
| E63 | ESA の模擬 | kept | M | MEAS | 放射計雑音 0.07–0.10 dB、log バイアス 0.03 dB | 0.73–0.99 dB（20 回平均 0.17） | 0.23–0.31 dB | — | （パルスでは標本分散） | 同 | 同 | T04 | §6 の出力 |
| E64 | 位相 → squeezing の写像 | kept | M | MEAS | 10 mrad で ΔL_eq 1.25% | 同形 | 同形 | 換算表 | 上限はモードの帯域 | 同 | f_low と平均法 | T04 | 全トピック共通の写像 |
| E65 | 速い位相雑音・dither | kept | M | MEAS | — | 295 kHz は probe 上（ループ応答分だけ） | dither 0.0555 rad → −5.63 → −5.18 dB | — | — | — | — | T04 | §7-10 |

### 1.10 誤差信号（E66–E76）

| ID | 効果 | 判定 | 段 | 族 | A | B | C | D | Fb | Fl | Fw | 出典 | 理由 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E66 | ヘテロダイン | kept | T | ERR | lock 1 2 MHz、lock 2 1 MHz | — | — | — | Liu 型の heterodyne 変種 | 同 | 同 | T01、T14 | Ref A |
| E67 | PD-OPA の信号 | kept | T | ERR | センサ雑音 4 µrad、誤差は A−/A+ = 1/123 で Δ へ | — | — | — | 統合案 A/B | 同 | 同 | T01、T14 | Ref A と D42 |
| E68 | PD-OPA までの非共通光路 | kept | T | PATH | 空気 16 rad/m、10 cm・1 K で 5.4 mrad（Φ1） | — | — | — | ファイバで分けると 79 mrad/µm、縮まない | 同 | 同 | T14 | 統合案の要点 |
| E69 | Ref B lock 1 | kept | T | ERR | — | m・P_in 掃引 | — | — | — | — | — | T14 | Ref B |
| E70 | Ref B lock 2 | kept | T | ERR | — | lock 1 誤差 × 0.146 | — | — | — | — | — | T14、T02 | Ref B |
| E71 | 量子雑音ロック | kept | T | ERR | — | — | 3.2 mrad（10 kHz）、log 4.2 mrad | — | ポートごとの QNL 変種（比較） | 同 | 同 | T14 | Ref C |
| E72 | probe の BHD 漏れ | kept | T | ERR/DET | 1 MHz ビート（帯域外） | ±295 kHz（帯域外） | — | — | GBS への混入は E43–E44 | 同 | 同 | T14、T05 | 漏れ情報の出力 |
| E73 | ショット雑音の位相検出 | kept | T | ERR | ヘテロダイン 2 倍 | — | — | 9.6×10⁻¹⁷ rad²/Hz | 1 nW で 6.4×10⁻¹¹ rad²/Hz | 同 | 同 | T14、T10 | 迷光との積の式 |
| E74 | RF 参照位相・ミキサ | kept | T | ERR | φ_d1 1:1（Φ1）、lock 2 は A−/A+ | 同形 | 同形 | — | heterodyne 変種で 1:1 | 同 | 同 | T14 | オフセット |
| E75 | DC fringe の誤差 | kept | T | ERR | — | — | — | V_off 15 mV → 4.5 mrad、ε 0.05 → 50 mrad、RIN は (S − V_off) 経由 | Liu の homodyne_dc で同式 | 同 | 同 | T15 | Ref D |
| E76 | 可視度の変化 | kept | T | ERR | — | — | — | K_d ∝ P𝒱 | 同式 | 同 | 同 | T15 | ループ利得 |

### 1.11 アクチュエータ・制御・デジタル（E77–E84）

| ID | 効果 | 判定 | 段 | 族 | A | B | C | D | Fb | Fl | Fw | 出典 | 理由 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E77 | EOM | kept | T | ACT/MOD | LO・probe（V_π 不明） | — | V_π 4 V、±4.5° | — | squeezed 経路以外 | 同 | 同 | T16 | §5.5 |
| E78 | PZT ストレッチャ | kept | T | ACT | — | — | レンジ不明（数分のロック時間から） | 6.49 rad/V、±26 rad、共振不明 | 915B 級は 63 m で 0.023 K | 同 | 同 | T15、T16 | Ref D と将来 |
| E79 | PZT ミラー | kept | T | ACT | — | 14・22 kHz | 10–20 kHz（掃引） | — | Liu 型の PZT ミラー変種 | 同 | 同 | T16 | Ref B・C |
| E80 | HV アンプ | kept | T | ACT | — | 自作 | HV アンプ | 914: 2 極 29.5/112 kHz、スルー 約 7 V/µs | — | — | — | T15、T16 | Ref D |
| E81 | 量子化 | kept | T | CTRL | — | — | 14–16 bit で差 < 10% | DAC 12.7 mrad/LSB → 1–2 mrad、limit cycle ±6.3 mrad | 細かい段に ≤ 1 mrad/LSB（16 bit 相当） | 同 | 同 | T15、T16、D49 | Ref A の 9 mrad と同じ桁 |
| E82 | ループ遅延 | kept | T | CTRL | 汎用 0.1–2 µs | 同 | 0.2–1 µs | 5.7 µs（上限 8.3） | ＋ probe 飛行時間 5 µs（1.06 km、下流のアクチュエータ）[AE-P2] | ＋1.6 µs（317 m） | ＋0.31 µs（63 m） | T15、T16、T10 | 遅延バッファ |
| E83 | 飽和・ロック外れ | kept | T | CTRL/TEC | — | — | 数分 | 1 mK/s で 約 1 分 | TEC 必須（1.06 km） | 同 | 63 m・10⁻⁴ K/s で 235 s | T08、T15、T16 | §7-10 |
| E84 | 制御器の構造 | kept | T | CTRL | 2 ループ | 2 ループ | 3 段スタック | DSC1 単一 | PZT＋TEC の 2 段 | 同 | 同 | T16 | §5.2 |

### 1.12 deep dive で見つかった効果（N01–N19）

| ID | 効果 | 判定 | 段 | 族 | A | B | C | D | Fb | Fl | Fw | 出典 | 理由 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| N01 | 残留 probe の古典雑音（RIN・位相雑音・ASE）が BHD に入る | kept | O＋M | DET | 残留量は不明（掃引）、1 µW で RIN 0.012、位相雑音は遅延差しだい | 同軸 probe で同形 | — | — | Liu 型の漏れ | 同 | 同 | T05 | survey になかった。ECL では pW–nW に |
| N02 | OPO 区間の区間 × ビーム係数差 | kept | T | PATH/SRC | — | 0.94°/nm | QNL は楕円を直接見る（なし） | — | 透過ポンプロック後 −0.55Δ_s/κ_c ≈ 40 mrad/nm | 同 42 mrad/nm | — | T02、T17 | §7-11 の前提を破る |
| N03 | squeezer 結晶の熱によるポンプ位相 | kept | T | PATH | lock 2 で相殺 | — | — | — | キャビティ側は N02 | 同 | 約 6 mrad/mK、probe を通さないと残る | T01 | D42 の結晶オプション |
| N04 | Liu 型ロックのレーザードリフト書き込み・(i)/(ii) トレードオフ | kept | T | PATH/MEAS | — | — | — | — | (i) に 2πδν(τ_k + τ_p − τ_r) | 同 | 310 ns で 1.95 rad/MHz | T05、T10 | T_a = T_b の設計 |
| N05 | PER の 2 結合点の干渉（寄生縞） | kept | O＋T | PATH | — | — | — | ±2.5 mrad（ER 26 dB）–±10 mrad、0.3–0.4 K 周期 | PM 同偏光で ≤ 3.2 mrad | 同 | 同 | T09、T15 | Ref D のロック点の遅い揺れ |
| N06 | 非相反（Faraday・Sagnac） | opt | T | PATH | — | — | — | — | ≤ 2×10⁻³ rad（63 m 最悪） | 同 | 同 | T09 | 無視 |
| N07 | 迷光の補助経路（反射・AOM 漏れ） | kept | O | NOISE | — | — | — | — | FC/APC 1 個で 78 光子/10 ns（1 mW） | 同 | 同 | T10 | 時間ずらしで避けられる |
| N08 | 周波数をずらした probe の熱ドリフト | kept | T | PATH | — | — | — | — | +7 GHz で 0.11 rad/K（63 m） | 同一周波数（0） | AOM 80 MHz で 1.2 mrad/K | T10、T17 | probe の周波数の選択 |
| N09 | 時間多重による誤差信号の標本化制約 | kept | T | ERR | — | — | — | — | 6 MHz でナイキスト 94 kHz [AE-P2] | 20 MHz で 313 kHz | 100 MHz で 1.56 MHz | T14 | 検出点の選択 |
| N10 | 積分型 BHD の薄まり | kept | M | DET | — | — | — | — | −10 → −0.12 dB | 同 | 同 | T11 | D45 の比較用 |
| N11 | OPO モードの裾による隣スロットの混合 | kept | O | CHK/SRC | — | — | — | — | f_rep ≲ 45 MHz（100 MHz で 13.8%） | ≲ 87 MHz | — | T17 | 繰り返しの上限 |
| N12 | FSR 櫛と FSR フィルタ | kept | O | FSR | — | — | — | — | FBG 30 dB、透過 約 80% | MZI > 95% | — | T17 | OPO の光子計数 |
| N13 | パルス OPO の多モード性 K | kept | O | SRC | — | — | — | — | K = 1.10 | K = 1.026 | — | T17 | GBS の純度 |
| N14 | LO の時間波形の整合 | kept | O | LOS | — | — | — | — | ガウス 0.951、矩形 0.864 | 同 | 超モードに整形 | T17、T03 | D50 |
| N15 | GBS 前の単一モード化フィルタ | kept | O | SPF | — | — | — | — | 不要 | 不要 | 必要（K ≈ 200） | T03、D44 | Fw の前提 |
| N16 | ポートごとの分散補償・遅延の微調整 | kept | O | SPF | — | — | — | — | 不要 | 不要 | 1 ps・100 MHz で必要 | T03、T06、D44 | Fw |
| N17 | レーザー内蔵の dither・SBS 抑圧 | kept | T | NOISE | — | — | — | オフで扱い（D35）、オンの変種は 0.31 rad（1 MHz 偏移）–310 rad | — | — | — | T05、D35 | Ref D の最大の不確かさになりうる |
| N18 | 明るい LO・probe 自身の量子雑音（線形化の破れ） | opt | Q | FIB | — | — | — | — | 無視 | 無視 | 100 MHz/63 m では線形化で十分 | T07 | TW は検証用 |
| N19 | PD 面で反射した LO の戻り | opt | O | DET | — | — | — | — | 値なし | 同 | 同 | T12 | UNVERIFIED、既定オフ |

## 2. プリセット

組み合わせは列挙しない。参照プリセット（Ref A–D）、将来の基本プリセット 4 つ、その上に重ねる「軸」（パラメータ）、比較用の悪い構成に分ける（D34）。
ブロック名は §3 の名前。「ƒ」は callable を許すパラメータ。

### 2.1 参照プリセット（Ref A–D）

| 軸 | `refA` | `refA_tap` | `refB` | `refC` | `refD` |
|---|---|---|---|---|---|
| 光源 | `opa_waveguide`: PPLN 45 mm、α = 9.06 /W、P = 0.64 W（出口、P の基準を明示）、損失は集中 8%（`distributed_loss` は切り替え）、1545.32 nm CW | 同左 | `opo_cw_spectrum`: T = 0.123、l = 0.5 m、L_ic(x) ƒ、P_th = 180 mW、x = √(P/P_th)、Ω_B = 2πf/γ_B | `opo_cw_spectrum`: 線幅 1.75 GHz FWHM を直接、η_tot 0.79 | なし（DX1 SC、1548.8 nm、3.33 mW） |
| レーザー雑音 | 白色 0.2 kHz ＋ RIN ピーク 0.7 MHz（代用、UNVERIFIED） | 同左 | Ti:S、89・177 kHz のピーク（FM/RIN を選ぶ） | 値なし（掃引） | h0 = 8×10³ Hz²/Hz、h₋₁ 0–5×10⁸、RIN −140 dB/Hz 平坦、dither・SBS オフ（N17 はフラグ） |
| probe | AOM 2 台 +1 MHz、EOM、tapping mirror で squeezer と PD-OPA へ（mW） | 同軸 probe 数 µW、tap r_t | 同軸（OPO 区間は係数差 N02）、PZT2 で 295 kHz dither | なし | 腕 1（915B）= pr、腕 2 = LO（conventions §5） |
| 光路 | 自由空間。非共通区間 tapping → 2 つの OPA（E68）、LO に EOM ＋ フィルタキャビティ | ＋ tap | 共有区間 | squeezed・LO が別々の SMF 1 km（180 m 敷設、ρ(f) 掃引） | PM MZI、ΔL 名目 9 m（0–15 m）、PNH1550R5A2（50:50 ± 1.5%、PER ≥ 25 dB）、Opneti 偏光子（ER ≥ 26 dB）、DX1 → 偏光子は SMF |
| 誤差信号 | #1 PD-OPA 2 f_sh（lock 1）、#1c BHD f_sh（lock 2） | #1b tap PD 2 f_sh ＋ BHD f_sh | #2 dither 極値（lock 1）、#3 BHD dither（lock 2） | #4 QNL（8–300 MHz、ログ、κ_det、1 MHz、m = 0.0785 rad） | #5 DC fringe（S = 0、+π/2、loop_sign） |
| 制御器 | `fpga_generic`（Tustin） | 同左 | 同左 | 3 段スタック（EOM → PZT ミラー → ストレッチャ）＋ ランプ、`fpga_generic`（τ_L 0.2–1 µs） | `dsc1`（h = 2 µs、τ 5.7 µs、12/16 bit、e = S − v、ki_entry 切り替え）、K_p・K_i 掃引 |
| アクチュエータ | EOM（probe: lock 1、LO: lock 2） | 同左 | PZT1 14 kHz、PZT3 22 kHz | EOM（V_π 4 V）、PZT ミラー（10–20 kHz）、ストレッチャ（レンジ掃引） | 914（2 極、スルー）＋ 915B（6.49 rad/V、共振 5–60 kHz・Q 3–100 掃引） |
| 検出 | BHD: QE 0.99、35 MHz、C 28 dB @ 15.4 mW | 同左 | QE 0.998、C −21.7 dB @ 5.5 mW | QE 0.97、C(f)、ディップ 0.96/0.83 GHz、暗雑音補正あり | PDB450C: G 10³ V/A、ℜ 1.0、±4.6 V、V_off ≤ 15 mV |
| 測定 | ESA 3 MHz、RBW 1 MHz、VBW 100 Hz、0–100 MHz、平均の種類は掃引 | 同左 | ESA 1 MHz、RBW 30 kHz、VBW 300 Hz | 30 MHz–1.21 GHz、RBW 300 kHz | 真の φ_D と e_DSC1 |
| 主な検証 | §7-1、§7-2 | §7-8 | §7-9、§7-11 | §7-10 | §7-7、ブラインド（D04） |
| オフにする opt | ASE（D49）、エタロン、熱力学、GNLSE（CW 近道） | 同左 | 同左 | GAWBS（別ファイバ） | T07 の散乱（D49）、ASE、熱力学、エタロン |

### 2.2 将来の基本プリセット

| 軸 | `F_wg100` | `F_wg6` | `F_bor6` | `F_liu20` |
|---|---|---|---|---|
| 光源 | `opa_supermodes`: 45 mm PPLN（代用）、τ_p 1 ps（100 fs–10 ps 掃引）、ℓ ƒ | 同左 | `pulsed_opo_source`: 66 MHz、T_p = 3 ns、κ_p/κ_c = 1.924、ρ 0.95、x_ss 1.3–2.4 | `pulsed_opo_source`: 125 MHz、T_p = 1.6 ns、κ_p/κ_c = 2.03、ρ 0.99、x_ss ≤ 4.1 |
| f_rep | 100 MHz | 6 MHz | 6 MHz | 20 MHz |
| 振り分け | Pockels ブロック N_b = 50（L_31 ≈ 3.2 km、E50） | Pockels 木 N_b = 1（L_31 ≈ 1.06 km） | 同左 | Pockels ブロック N_b = 4（L_31 ≈ 1.27 km）。変種: 受動ループ（Liu の構成） |
| FSR・スペクトル | 単一モード化フィルタ（N15）＋ ポートごとの分散補償（N16） | 同左 | FBG の組（透過 0.80） | MZI 3 段（> 0.95） |
| 経路の組み方 | LL（融着 ＋ 自由空間の入れ出し）、PM | 同左 | 同左 | 同左 |
| probe | Liu 型: 逆向き、PM・直交、同じ周波数・同じレーザー、パルス化、P̄ 1 µW（nW–mW 掃引）、自由空間 PBS、検出は振り分けの手前 | 同左 | 同左（同じ周波数なら ≲ 1–3 nW）。変種 +7 GHz ＋ FBG | 同左 |
| ポート光路のロック | `L-Liu` ＋ TEC（同じ区間の温度）＋ PZT、アクチュエータは上流 | 同左 | 同左 | 同左 |
| ポンプ・光源のロック | `L-LiuPD-A`（結晶オプション (a): 温度 1 mK） | 同左 | OPO 長さ: Borealis 型ロック光（+7 GHz、直交、CW）。ポンプ位相ロックなし、パワーは AOM | OPO 長さ: Liu 型（同じ偏光・同じ周波数・時間ずらし）。ポンプ: (d) 透過ポンプのホモダイン |
| LO・BHD | パルス LO（分岐 ＋ 遅延）、超モードに整形、パルス分解 BHD（GHz 級） | 同左（数十 MHz 級） | パルス LO をガウス整形（0.951）、ゲート ≥ 25 ns、BHD 200–500 MHz | 同左 |
| 制御器 | 細かい段 `moku_delta` か `fpga_generic`（14–16 bit）＋ TEC の遅い段 | 同左 | 同左 | 同左 |
| ポート数 | 1 ポートで L_k 掃引、2 ポート組（D13） | 同左 | 同左 | 同左 |
| チェック | E50: N_b = 1 は不可の警告 | — | N11: 6 MHz < 45 MHz で可 | N11: 20 MHz < 87 MHz で可 |

### 2.3 将来のプリセットに重ねる軸（変種）

| 軸 | 値（名前） | 出典 |
|---|---|---|
| 光源 | `waveguide_opa`（ℓ = 10–45 mm、GVM を設計した導波路も）、`borealis_opo`、`liu_opo`、`cw_opo_refB`・`cw_opo_refC` | D07、D43、D44、D50 |
| 振り分け | `pockels_tree`、`pockels_block(N_b)`、`passive_loops`、`ln_fiber_tree`（比較用）、`none`。配置 `tree` / `chain` | D39、T11 |
| 経路の組み方 | `NF`（FC/APC、約 −1.5 dB）、`MF`（融着 ＋ ファイバ PBS、約 −4.7 dB）、`LL`（融着 ＋ 自由空間、約 −7.3 dB。ただし振り分け 1 段の仮定、§5 の C20） | D38、T12 |
| ファイバ | `pm`、`smf`、軸の入れ替え（90° 融着）の位置のリスト | D32、T09、T10 |
| probe | パワー 0.1 nW–1 mW、CW / パルス（デューティ、時間ずらし）、偏光 `same` / `orthogonal`、周波数オフセット 0 / 80 MHz / 7 GHz / 100 GHz / 1 THz、光源 `same_laser` / `other_laser`、分離 `pbs` / `circulator` / `wdm` / `coupler`、検出 `homodyne_dc` / `heterodyne` / `dither`、検出点（振り分けの前 / 後）、アクチュエータ位置（上流 / 下流） | D40、T10、T14 |
| ロック | §2.4 の名前 | D30、D42、D50 |
| 大レンジ | PZT のみ、PZT ＋ TEC、PZT ＋ スプール温度制御、2π リセット（`fringe_hop`）、EOM（squeezed 経路以外） | D41、T16 |
| LO | パルス（分岐 ＋ 遅延）、パルス（振り分け）、CW ＋ デジタル重み。波形 `optimal` / `gauss` / `rect` / `one_sided_exp`。LO 配送のロック on / off | D09、D12、D45、D50 |
| BHD | `pulse_resolved`、`gated`（ゲート ≥ モード長）、`integrating`（比較用） | D45、D46、D50 |
| OPO のロック光 | `borealis`（+7 GHz、直交、CW）、`liu`（同じ偏光・同じ周波数・パルス） | D50 |
| ポンプ | パワーロック on / off、位相ロック `none` / `(b)` / `(d)` | D42、D50 |
| 制御器 | `dsc1`、`moku_delta`、`fpga_generic` | D05、T16 |
| ポート | 1 ポート、2 ポート組（ρ、Δτ）、32 ポート（Phase 6） | D13 |

### 2.4 ロックのトポロジーの名前

| 名前 | 中身 | ロックする位相 | 出典 |
|---|---|---|---|
| `L-A` | Ref A: PD-OPA の 2 f_sh（lock 1）＋ BHD の f_sh（lock 2、probe は squeezer を通る） | Φ1^det、Δ | T01、T14 |
| `L-tap` | 従来: squeezer 出口の tap で lock 1 ＋ BHD lock 2 | Φ1^sq、Δ | T12、T14 |
| `L-B` | Ref B: 透過 probe の dither 極値（lock 1）＋ BHD dither（lock 2） | Φ1、Δ | T14 |
| `L-C` | Ref C: 量子雑音ロック、3 段スタック | Δ | T14、T16 |
| `L-D` | Ref D: DC side-fringe | φ_D | T15 |
| `L-Liu` | ポートごとの逆向き probe と共通参照光の干渉（D30） | φ_pr,k − φ_ref,k | T10、T14 |
| `L-LiuPD-A(a/b/c/d)` | 統合案 A: PD-OPA で共通参照光とポンプをロック ＋ 各ポート `L-Liu`。結晶・キャビティの位相: (a) 結晶温度 約 1 mK、(b) 残ったポンプと 2ω の参照光のロック、(c) squeezer に順方向の probe（GBS に光子、比較用）、(d) 透過ポンプのホモダイン（Liu、OPO） | 共通項と ψ_k | D42、D50、T14 |
| `L-LiuPD-B` | 統合案 B: 逆向き probe をポンプと直接比べる（ポンプの係数が 2 倍、比較用） | 2φ_pr,k − φ_p | T14 |
| `L-A-seed` | Ref A の lock 2 を GBS 構成で使う（seed が squeezer を通る、比較用） | Δ | D42、T10 |
| `L-LOdeliv` | LO 配送のファイバ雑音補償（オン・オフ） | LO_k の配送光路 | D12 |

### 2.5 比較用の悪い構成（D34）

| 名前 | 何を示すか | 出典 |
|---|---|---|
| `bad_all_fiber_NF` | FC/APC とファイバ部品だけの squeezed 経路（上限 約 1.6 dB） | D38、T12 |
| `bad_ln_switch_tree` | ファイバ LN スイッチの木（10–20 dB の損失） | D39、T11 |
| `bad_mw_probe_samepol` | 1 mW・同じ偏光・同じ周波数の逆向き probe（迷光 3×10²/10 ns、OPO では 435/モード） | D34、T10、T17 |
| `bad_orth_pm_no_tec` | 直交偏光 probe を PM で TEC なし（129 rad/K の誤差） | D34、T09、T10 |
| `bad_smf_orth` | SMF・直交偏光 probe（予測できない 2 arg(a)、迷光 1/3） | T09、T10 |
| `bad_bor_100MHz` | Borealis 型 OPO を 100 MHz（裾 13.8% の混合） | D50、T17 |
| `bad_no_fsr_filter` | FSR フィルタなし | T17 |
| `bad_refA_seed_in_gbs` | `L-A-seed`（1 µW で 10 ns に約 8×10⁴ 光子） | D42、T10 |
| `bad_integrating_bhd` | 全パルス列の LO と積分型 BHD（−10 → −0.12 dB） | D45、T11 |
| `bad_eom_on_squeezed` | squeezed 経路にファイバ EOM（上限 2.6–3.6 dB） | T12、T16 |
| `bad_liuPD_B` | `L-LiuPD-B` | T14 |
| `bad_copropagating_probe` | 同じ向きの明るい probe（GAWBS・XPM・ラマンが squeezed に入る） | T07 |
| `bad_moku_pid_05V` | Moku:Delta PID の ±0.5 V で 915B を駆動（レンジ 10–20%） | T16 |
| `bad_wdm_1THz_probe` | 1 THz ずらした probe（15 rad/K の差動ドリフト、ラマン増） | T10、T07 |
| `bad_gauss_rect_lo` | ガウス・矩形の LO（0.951 / 0.864） | T17 |
| `bad_long_pump_opo` | T_p κ_c = 5（K ≈ 2） | T17 |

### 2.6 プリセットでの opt 効果の扱い

opt の 9 件は、すべての参照プリセットでオフ。比較用・アブレーションでオンにする: E17（ASE を掃引）、E33（`bad_copropagating_probe`）、E30（probe パワーの上限の確認）、E02（PD-OPA の高利得の比較）、E24・E28（床の確認）、N06・N18・N19（検証用）。

## 3. Phase 3 のためのブロック一覧

### 3.1 State のフィールド（提案。conventions §9–§10、D47 に沿った名前の案で、Phase 3 で決める）

| フィールド | 中身 | 出典 |
|---|---|---|
| `t`、`dt` | 時刻、時間刻み（多重レートなら刻みごと） | T16 |
| `laser[s].dnu`、`laser[s].phiL`、`laser[s].dp`、`laser[s].hist` | 周波数ゆらぎ、位相 φ_L = −2π∫δν、相対強度、放出時刻で評価するための履歴（端数の遅延を補間） | T05、conventions §7 |
| `seg[i]`（静的） | id、L、n_g、τ_g0、β2、β3、向き d、種類（PM/SMF/自由空間/OPO/結晶/振り分け）、熱ゾーン、N_sub、`shared_with_probe` | conventions §9、D47、T08、T11 |
| `seg[i]`（動的） | ΔT、歪み ε、各小区間の値、2×2 Jones 行列 J（PM は軸ラベルの簡易版）、群遅延 τ_g(t)、ビームごとの位相 ψ[beam] | D47、T03、T08、T09 |
| `coef[i][driver][beam]` | 区間 × 要因 × ビームの係数（k_T,a、k_ε,a、R_τ、OPO の離調の係数・表、Δν_pr による項） | D47、T02、T10、T17 |
| `beam[b]` | 組み立てた位相 φ、パワー P、偏光、放出時刻、周波数オフセット | conventions §9 |
| `pump`、`opa[j]`、`opo[j]` | φ_p、P、x_ss、φ_OPA,0、r、δθ_mm、Δ_s、Δ_p | conventions §4、T01、T17 |
| `v[name]`、`dem[name]` | 電圧（誤差信号・検出器出力）、復調後の I/Q | T14 |
| `ctrl[name]` | 出力 u、積分値、飽和フラグ | T16 |
| `delay[name]` | 遅延バッファ（整数 ＋ 端数） | prompt §5.2、T16 |
| `act[name]` | 機械の状態（共振のモード）、δψ、TEC の温度 | T16 |
| `timing[k]` | Δt_k(t)、Δt_kj(t) | T03 |
| `leak[name]` | probe のパワー・タイミング・偏光・周波数（迷光の計算用） | T14、T10 |
| `metric` | (i) δ_i,k、(ii) Θ_k と 4 成分、ポート間の重なり、迷光（クリック数 ＋ モード内光子数）、η_k（GBS 入力と BHD の両方） | D12、D31、D46、D47 |
| `lock[name]` | ロック状態、滑り・リセットのイベント | T16 |

### 3.2 ブロック（族ごと）

「読む」「書く」は §3.1 のフィールド。ƒ は callable を許すパラメータ。

| 族 | ブロック | 段 | 読む | 書く | 主なパラメータ | 必要な表 | 出典 |
|---|---|---|---|---|---|---|---|
| SRC | `laser_source` | T | `t` | `laser[s].*` | h0、h₋₁、ピーク、ドリフト ν̇、dither、SBS、RIN、seed。プリセット `cobrite_dx1_SC` など | `laser_oob_variance` | T05、T15 |
| SRC | `opa_waveguide` | O | — | Gaussian チャネル（S、N）、r_eff、L_wg、δθ_mm、ψ_p,xtal | α、P（基準 in/out）、ℓ、a_s、a_p、Δk(T)、φ_OPA,0、L_wg(P) ƒ | OPA 表（P、T） | T01 |
| SRC | `single_mode_reference` | O | — | R′± | L、α、P、θ（保護リスト） | — | T01、T04 |
| SRC | `opa_supermodes` | O | — | ψ_n、r_n、φ_n、K、C、S | 結晶の分散、ポンプ包絡線、格子 | → `mode_overlap_table` | T03 |
| SRC | `opo_cw_spectrum` | O | — | R±(f_s)、2×2 S(Ω)、ρ、κ | T、L_ic ƒ(x, P, t)、l または線幅、P_th、Δ | — | T02 |
| SRC | `pulsed_opo_source` | O | — | C（2N×2N）、r_n、K、u_opt、τ_eff、\|ũ(f)\|²、ε_tail | Δν_c、ρ、κ_p、ポンプ波形、x_ss、Δ_s、Δ_p、φ_p,in | `opo_angle_map` の表 | T17 |
| SRC | `opo_detuning_noise` / `opo_length_lock` | T | `seg`（OPO）、温度 | `opo[j].Delta_s/Delta_p` | 長さ雑音、ロック光の種類、Δ(dn/dT) | — | T02、T17 |
| SRC | `opo_angle_map` | T | `opo[j]` | OPO 区間の `seg[i].psi[beam]` | x_ss、ロック状態 | (Δ_s, Δ_p, x_ss) 表 | T17 |
| SRC | `pump_power_stabilization` | T | `pump.P`、RIN | `pump.x_ss` | AOM ループ | — | T17 |
| SRC | `crystal_thermal` | T | 温度 | `seg[xtal].psi[pump]`、Δk | 結晶長、dn/dT、TEC | — | T01 |
| PATH | `segment_phase` | T | `seg[i]` の要因、`coef` | `seg[i].psi[beam]`、`tau_g` | ゾーン、空間構造 lumped/distributed/mixed | `fiber_env_coeffs` | T08、T10 |
| PATH | `segment_jones` | T | `seg[i].dT`、D_pol | `seg[i].J` | 種類（PM 対角 / SMF）、K_T、L_b | `pm_fiber_params` | T09 |
| PATH | `beam_assembler` | T | `laser.hist`、`seg[*]` | `beam[b].phi`、偏光 | 光路の区間リスト、向き | — | conventions §9、T05 |
| PATH | `delay_buffer` | T | 任意 | 遅らせた値 | τ（整数 ＋ 端数、線形補間 / Thiran） | — | T05、T16 |
| PATH | `filter_cavity` | O＋T | 離調 δ(t) | LO の位相・透過、H_FC | Δν_FC、FSR、T0 | — | T05 |
| PATH | `polarizer_projection`、`pdl_element`、`faraday_nonreciprocal` | T | `beam` の偏光 | 透過パワー、射影振幅 | ER、PDL、V_F | — | T09 |
| NOISE | `ambient_temperature`、`segment_thermal` | T | `t` | ゾーンの T、`seg[i].dT` | r_T、b_T、τ_th、周期成分 | — | T08 |
| NOISE | `segment_mech_noise` | T | `t` | `seg[i].strain`（小区間ごと） | h_ref、L_ref、γ、f_c、f_lo、ピーク、ρ_ij(f) | — | T08 |
| NOISE | `segment_thermodynamic`（opt） | T | — | δψ^thd | Dong の物性値 | — | T08 |
| NOISE | `gawbs_phase_jitter`（opt） | T | — | 区間の白色位相、または θ² の定数 | σ_G² | `gawbs_table` | T07 |
| NOISE | `rf_reference` | T | — | φ_RF,j(t) | S_φRF、τ_RF | — | T14 |
| NOISE | `pd_offset_drift` | T | 温度 | V_off(t) | 15 mV、ドリフト | — | T13、T15 |
| MOD | `eom_phase` | T | u | 区間の δψ（squeezed 経路なら `loss_budget` に登録） | V_π,LF/RF、IL、RAM | — | T16、T14 |
| MOD | `dither_modulator` | T | — | 指定光路の位相 m sin ω_m t | m、f_m、a_RAM、φ_a | — | T14 |
| MOD | `aom_shift`、`probe_pulser` | T | δf_RF / φ_RF、ゲート | 位相、周波数、パワー（`leak`） | m 次、ε_AOM、デューティ、t_offset | — | T16、T10 |
| ERR | `parametric_beat_detector` | T | φ_p、φ_pr（相互作用点）、A± | `v`（DC ＋ 2f_sh I/Q） | f_sh、T_tap、ℜ、G、i_n | A±(P) | T14、T01 |
| ERR | `bhd_probe_beat` | T | φ_LO、φ_pr、φ_p、A± | `v`（f_sh または f_m の I/Q） | BHD の R(f) | — | T14 |
| ERR | `quantum_noise_power` | T | Δ | `v`（包絡線・ログ ＋ radiometer 雑音） | f1、f2、κ_det、検出器の法則 | P(Δ)、∫S_x²(Δ) | T14 |
| ERR | `dc_fringe`（`mzi_dc_fringe_detector` と `balanced_pd_dc` を統合） | T | 2 ビームの場（位相・パワー・偏光成分） | `v` | κ_1、κ_2、T_opt、ℜ、ε_R、G、V_off、飽和、μ | — | T13、T14、T15 |
| ERR | `liu_probe_interferometer` | T | φ_pr,k（逆向き）、φ_ref,k | `v` | mode（homodyne_dc/heterodyne/dither）、δf、デューティ、ゲート、標本化（N09） | — | T14、T10 |
| ERR | `demodulator` | T | `v`、φ_d | `dem` | f_ref、LPF、V_off | — | T14 |
| CTRL | `pid` | T | e | `ctrl.u` | 正準係数、離散化（dsc1/moku/tustin）、符号、制限、anti-windup | — | T16、T15 |
| CTRL | `adc`、`dac` | T | 電圧 | 量子化した電圧 | ビット、レンジ、f_s、平均、雑音、LPF、出力抵抗と負荷 | — | T16 |
| CTRL | `iir`、`hpf`、`lpf`、`adder`、`ramp` | T | 任意 | 任意 | SOS 係数、周期・振幅 | — | T16 |
| CTRL | 並列（帯域分割）と積み重ねの合成 | T | e | 各段の u | 段の接続 | — | T16、prompt §5.2 |
| CTRL | `lock_supervisor` | T | e、u、飽和、真の位相 | `lock` | しきい値、待ち時間、再引き込み（reset_center/fringe_hop/ramp_search） | — | T16 |
| ACT | `pzt_fiber_stretcher` | T | `delay` 経由の u | 区間の δψ（軸ごと）、飽和 | K_s（ΔOPL 読み）、σ_act、モード（f_i、Q_i、r_i）、ヒステリシス、クリープ | 914 Fig. 1 の表 | T16、T15 |
| ACT | `hv_amp_914` | T | u | 駆動量 | 2 極（29.5/112 kHz）または 1 極（28 kHz）、スルー、入力レンジ、雑音 | 同 | T15、T16 |
| ACT | `pzt_mirror` | T | u | δψ | d_x、θ_i、モード | — | T16 |
| TEC | `thermal_actuator` | T | u | **区間（小区間）の温度 ΔT_h**（ψ は `segment_phase` が係数から計算、C13） | K_T、τ_th、τ_tr、ΔT_max、L_h | — | T08、T10、T16 |
| DET | `detector_efficiency` | O | — | η_det(f) と内訳 | η_q、C(f)、δ(f)、S_RIN | — | T13 |
| DET | `bhd_frontend` | M | 信号の R(f) / 直交位相 | 差電流の時系列・PSD | H(f)（TIA ＋ dip）、S_el、G、subtract_dark | — | T13 |
| DET | `pulsed_bhd_kernel` | O | — | c_m、C_pulse、隣接相関 | h(t) または B・Q、w(t)、σ_el、N_LO | — | T13 |
| DMX | `demux_design` | O | — | L_k、τ_k、η_sw,k、ψ_sw,k^(0)、漏れ行列、実現可能性 | K、f_rep、N_b、L_0、素子のプリセット | — | T11 |
| DMX | `demux_switch`、`delay_line_k` | T | 温度、駆動雑音 | 区間 ψ_sw,k(t)、σ_sw,fast²、τ_k(t)、δt_k(t) | κ_pm、V_π、バイアスのモデル | `demux_design` | T11 |
| LOSS | `loss_budget` | O | 全要素の登録 | 経路ごとの η（GBS 入力・BHD・LO・probe） | 要素（id、IL、区間、種類、出典の印） | — | T12 |
| LOSS | `pure_loss`、`mode_overlap`、`per_budget`、`coil_birefringence`、`tap` | O / Q / T | — | η、σ_out | η、ビーム径、PER、巻き半径、r_t | — | T12、T09 |
| FSR | `fsr_comb_filter` | O | — | (0,0) の損失、n_FSR | FSR、r_m ƒ、T_f(ν) | — | T17 |
| SPF | `spectral_filter`、`dispersion_comp` | O | — | 帯域・損失、補償後の φ2・φ3 | 帯域、GDD | — | T03、D44 |
| LOS | `lo_temporal_shape` | O | — | η_t、θ_off | 波形、BHD 帯域 | — | T17、T03 |
| FIB | `fiber_gnlse` | G | — | A(L,T)、保存量 | 項のフラグ（loss、beta2、beta3、kerr、raman、self_steepening、xpm_external、cw_shortcut）、刻み | `fiber_segment_params` | T06 |
| FIB | `lo_nonlinear_table`、`fiber_linear_table`、`xpm_counterprop_coeff`、`spm_rin_to_phase` | O | — | 表 | P_avg × τ × L_k | — | T06 |
| FIB | `fiber_quantum_linearized`、`fiber_quantum_tw`（opt） | Q | 古典解、V(0)、LO モード | V(L) または 2 次形式、加わった雑音 N | 雑音のフラグ | 散乱表 | T07 |
| FIB | `raman_noise_table`、`gawbs_table`、`brillouin_table`、`stray_light_table` | O | — | n(f)、N_mode、N_stray(P̄_pr, …) | モデル、L、T、P、ε_pol、b_R | — | T07、T10 |
| FIB | `path_mode_transfer`、`mode_overlap_table`、`timing_state`、`mode_overlap_lookup` | O / T | 区間の列、Δt_k(t) | (V_−, V_+, θ_off)(Δt)、ポート間 η | 格子、ジッタの速い成分 | — | T03 |
| MEAS | `phase_to_quadrature` | M | δ(t)（指標ごと）、V± の表、χ_HF | V′(t)、R′±、帯域ごとの σ | 平均法、dither のモード | `vpm_table`、`hf_phase_noise_integral` | T04 |
| MEAS | `esa_zero_span`、`esa_swept_spectrum` | M | V′(t)、V′(f) | トレース | f_s、RBW、ENBW、VBW、平均の種類、seed | — | T04 |
| MEAS | `metric_i`、`metric_ii`（4 成分）、重なり、迷光 | T/M | `seg`、`beam`、`timing`、`leak` | `metric` | 分解のための雑音源の切り替え | — | T04、T10、D47 |
| MEAS | `required_theta_contour` | O | — | σ_req(S_t, L_f) | ポンプ最適 / 固定、dither | 損失モデル | T04、T12 |
| CHK | `sbs_check`、`rep_rate_tail_check`、`demux_feasibility`、`grid_edge_monitor` | O | — | 警告（計算は続ける、D34） | しきい値 | — | T07、T17、T11、T06 |

### 3.3 オフラインの表（O → T の受け渡し）

| 表 | 作るブロック | 使うブロック | 次元 |
|---|---|---|---|
| (V_−, V_+, θ_off) | `opa_supermodes` / `pulsed_opo_source` / `opo_cw_spectrum` ＋ LO モード | `mode_overlap_lookup`、`phase_to_quadrature` | Δt（＋ φ2）、OPO は (Δ_s, Δ_p, x_ss) |
| V±(f) の損失込み（指標 (i) と (ii) で損失の集合を変える） | `vpm_table` | `phase_to_quadrature`、`quantum_noise_power` | f |
| P(Δ)、∫S_x²(Δ) | T14 の前処理 | `quantum_noise_power` | Δ |
| σ²_oob、χ_HF | `laser_oob_variance`、`hf_phase_noise_integral` | `phase_to_quadrature` | 経路の組 |
| N_stray、n_mode | `stray_light_table` | `metric`（迷光）、squeezing の劣化 | P̄_pr、偏光、Δν_pr、L、T_gate |
| LO の非線形 | `lo_nonlinear_table` | T03、T04 | P_avg × τ × L_k |
| demux の設計 | `demux_design` | `demux_switch`、`loss_budget` | ポート |
| 損失収支 | `loss_budget` | すべての V± | 経路 |

### 3.4 ノート間で食い違うインターフェース（Phase 3 で決める）

| # | 食い違い | 提案（決定ではない） |
|---|---|---|
| C1 | 偏光: T09 は区間の 2×2 Jones 行列、T08・T10 は軸ごとの係数 k_T,a・κ_{j,a}。D47 は Jones ＋ PM の軸ラベル簡易版を採用したが、SMF の J を環境の要因（ΔT、ε）からどう作るかは T08 にない（T09 は拡散係数 D_pol） | 要因（ΔT、ε）は区間で共通に生成し、`segment_jones` が区間の種類ごとのモデル（PM: 対角、SMF: 複屈折ベクトルの拡散 ＋ 温度係数）で J を作る |
| C2 | 帯域の勘定: T05（σ²_oob）、T04（χ_HF）、T07（GAWBS を定数で θ² に）、T11（σ_sw,fast）、T03（σ_t,fast を表に）、T15（帯域外を可視度 γ_coh に）がそれぞれ「ループで表せない高周波」を別の場所に足す | 各雑音源が「時間ループで表した帯域 [0, f_N]」を宣言し、その外は 1 か所（`phase_to_quadrature`）でだけ足す。二重に数えない検査をテストにする（T04 V13） |
| C3 | 損失の集約: T12 `loss_budget`、T11 η_sw、T09 η_PER、T03 η_t、T17 ρ·η_filter、T13 η_det、T10 の入れ出し素子、T11 の重なりの placeholder | 1 つの登録表にまとめ、GBS 入力と BHD を分ける（T12 V13）。モードの重なり（空間 T12、時間 T03・T11、LO 波形 T17）は 1 回だけ数える |
| C4 | 光源の出力: T01 は単一モードの Gaussian チャネル、T02 は R±(Ω)、T03 は超モード、T17 は時間領域の 2N×2N 共分散 | 共通の出力「光源の共分散 ＋ モードの基底 ＋ LO モードへの射影関数（→ (V_−, V_+, θ_off)、CW なら f_s ごと）」にする。`opo.pulsed = off` で T02 に戻す（T17） |
| C5 | OPO 区間の係数: T02 は線形の係数（squeezed −(Δ/κ)/(1 + x²)、probe −Δ/(2κ)）、T17 はパルスでは表（δθ(Δ_s, Δ_p, x_ss)） | `coef` は定数・callable・表のどれでも受ける |
| C6 | 誤差信号の時間モデル: T14 は複素包絡線（解析的な復調 ＋ 復調後の雑音 PSD）、T16 はコンパイルした 1 サンプルごとのループで 1 MHz の RF も直接扱える | ブロックの API はループレートの電圧に統一し、`rf_rate_mode` は検証用（T14 V1・V2・V6） |
| C7 | 914 の小信号モデル: T15 は 1 極 28 kHz、T16 は 2 極 29.5/112 kHz | T16 の 2 極を参照（6 点を 2% で再現）、1 極は代替 |
| C8 | 915B の共振の掃引範囲: T15 は 10–60 kHz・Q 5–30、T16 は 5–50 kHz・Q 3–100 | 和集合 5–60 kHz・Q 3–100 |
| C9 | 熱とストレッチャの群遅延の比: T03 は残る割合 0.21、T08 は 0.28 | 定数にせず (k_T、κ_τ、ξ_eo、n_g/n) から計算する（D47 の「要因ごとの効き方」） |
| C10 | 逆向き伝搬: T08 は「同じ時刻に両端から出る 2 つの光の差」、T10 は「端のアクチュエータでロックした後の残差」（4 倍違う式） | どちらも小区間の遅延バッファで自動的に出る量。アクチュエータの位置を小区間の位置として持たせ、テストは観測量を明示（T08 V6、T10 3） |
| C11 | DC fringe のブロック: T13 `balanced_pd_dc`、T14 `dc_fringe`、T15 `mzi_dc_fringe_detector` | T15 の式 (3)（μ、PER、κ_1・κ_2、ε_R を含む）を上位集合として 1 つにまとめる |
| C12 | レーザーのブロック: T05 `laser_noise_source`、T15 `cw_laser` | 同じブロック |
| C13 | TEC: T16 のモデルは δψ を直接書く（σ_act c_T L_h ΔT）。T10 は TEC が「雑音と同じ区間の温度」を動かすことで軸の差が打ち消されると述べる | TEC は小区間の温度 ΔT_h を書き、位相は `segment_phase` が係数から計算する（軸ごと・群遅延も同時に） |
| C14 | 指標 (ii) の分解: T04 は雑音源の切り替えで 3 成分、D47・T10 は 4 成分（squeezed 光だけの区間を追加） | 4 成分。線形性の検査（T04 V12）と、飽和があるときは分解が無効と出力する |
| C15 | 迷光の受け渡し: T14 は `state.leak[name]`、T10 の表は P̄_pr・タイミング・偏光・周波数・ε_AOM を要求 | `leak` のフィールドを T10 の表の入力に合わせる |
| C16 | 測定: T04 の ESA（CW のゼロスパン）と、パルス BHD の標本分散（T04 §3.4、T13 §1.8） | `MEAS` に `esa` と `pulse_sample` の 2 モードを持つ |

## 4. 検証テストの一覧

各ノートの V テストを重複をまとめて番号 VT-xx にした。§7 の番号への対応、目標、出典のテスト番号。D48 の読み替え（§7-1・§7-2・§7-6・§7-9）を適用した。

### 4.1 prompt §7 の項目

| VT | §7 | 内容 | 目標 | 出典 |
|---|---|---|---|---|
| VT-01 | 1 | 単一モードの解析式 | L 8%、α 9.06 /W、P 0.64 W、θ 9 mrad で R′− = −10.148 dB、R′+ = +20.556 dB（±10⁻³ dB）。小角と厳密の差 < 10⁻³ dB | T04 V1、T01 V-T01-1、T12 V4 |
| VT-02 | 1 | 全体のシミュレーション（参照の損失収支、ロックを 9 mrad に調整） | **R′− = −10.15 ± 0.10 dB、R′+ = +20.56 ± 0.10 dB、導波路は集中損失**（D48。T01 の 20.4 ± 0.3 と T04 の ±0.05 を置き換える）。分布損失モデルは選択肢（R′+ が 0.27 dB 下がる） | D48、T01 §7.3、T04 §11 |
| VT-03 | 1 | CW の極限でパルスのブロックが恒等 | LASER_MODE = CW で VT-01 の値がビット単位で同じ。GNLSE の CW 近道と RK4IP が 10⁻¹² | T03 V5、T06 V5 |
| VT-04 | 2 | Fig. 5 の 4 本の曲線 | 4 本 × 8 点を ±0.01 dB、P_opt 612 mW・−10.152 dB、503 mW・−8.397 dB。データ点との差 ±0.35 dB。**「前の研究」の破線の右端は形の比較から外す**（D48） | T01 V-T01-2・6、T04 V2 |
| VT-05 | 3 | 純損失 | V_out = ηV_in + (1 − η)/2 を squeezed・反 squeezed・熱的な入力で相対 10⁻¹²。直列の順序に依らない、受動の分布損失 = 集中損失、回路雑音 = 損失（C = 28 dB で 0.158%）、単一モードの不整合 = 損失、OPO の R± − 1 が η_tot に線形、パルス OPO の最適モードが ρe^{−2r} + 1 − ρ に ±0.05 dB | T12 V1–V3・V8、T07 V5、T11 1、T13 1・2、T10 12、T03 V7、T02 V-B2、T17 V-P13 |
| VT-06 | 4 | 真空は真空のまま | ファイバの線形化で max\|V(L) − I/2\| < 10⁻¹²（損失の雑音込みでも）。OPA の g → 0、OPO の x = 0、パルス OPO の x_ss = 0 で固有値 1（初期のキャビティ内真空の項を含むこと、V-P3 の回帰テスト）。LO がすべての超モードに直交なら 1/2 | T07 V1・V1b、T01 V-T01-3、T02 V-B2、T17 V-P2・P3、T03 V4 |
| VT-07 | 5 | 基本ソリトンの不変性 | 5z0 で解析解との相対 L2 誤差 ≤ 10⁻⁵ | T06 V1 |
| VT-08 | 5 | 2 次ソリトンの回帰 | z0 で ≤ 10⁻⁵ | T06 V2 |
| VT-09 | 5 | 刻み半減で誤差 約 1/16 | 固定刻み・外挿なしで連続する比 13–19（漸近域で判定）。接線 RK4IP（量子層）でも同じ | T06 V3、T07 V2 |
| VT-10 | 5 | Raman ありの光子数の保存 | \|ΔN_ph/N_ph\| ≤ 10⁻⁵、ΔE < 0、重心の移動 < 0（BW・LA） | T06 V4 |
| VT-11 | 6 | Raman 雑音と温度（**D48: 熱の部分だけ 0 になる**） | T = 1 mK で熱の部分が 295 K の値の 10⁻⁶ 未満、anti-Stokes = 0、Stokes = g_R(Ω)PL（100 GHz で 2.39×10⁻⁶ /(W m)）。詳細釣り合い aS/S = e^{−ħΩ/k_BT}（13.2 THz・295 K で 0.117）。g_R のピーク 0.68 /(W km)・13.1 THz | T07 V6・V6b・V6c |
| VT-12 | 6 | GAWBS のピーク | 125 µm で R0,m = 30.5, 82.1, …, 323.3 MHz（±0.1 MHz）、7 番目が 321.75 MHz と 1% 以内、80 µm で ×125/80。分散が L・T に比例 | T07 V7・V7b |
| VT-13 | 7 | 閉ループ伝達関数 = 解析解 | Ref D: z 領域の解析式と 10⁻³・0.1°、連続の 1/(1 + L) と f ≤ f_c で 5%・3°、P のみの発振 K_p,u = 0.093・40.9 kHz（±2%）。汎用の線形ループで S(f) が 0.5 dB・3° | T15 4、T16 V10 |
| VT-14 | 7 | フィードバックなし → free-running | 真の φ_D の PSD が外乱の和に Welch の誤差内。レーザーの床 (2πτ)²h0 = 7.5×10⁻¹⁰ rad²/Hz（ΔL 10 m、25 kHz）に ±10%。分散比 1 | T15 5、T16 V10、T05 V6 |
| VT-15 | 7 | 積分器ループの残差と lock 2 の縮小 | σ² = πS_0/(2f_c) に 5%。lock 1 のセンサ雑音だけ入れたとき Δ の残差が A−/A+ 倍（10% 以内） | T15 7、T14 V3 |
| VT-16 | 8 | tap のトレードオフ（2 つの定義） | (a) NEP 1 pW/√Hz・10 µW で R′−(r_t) が内部に最小、r_t ≈ 3×10⁻³ で −9.94 ± 0.02 dB。(b) θ を 9 mrad に固定して T_tap = 0.1/1/3% で −10.11/−9.76/−9.07 ± 0.01 dB。PD-OPA の構成では r_t に依らず −10.15 dB | T12 V10、T14 V9 |
| VT-17 | 9 | Ref B Fig. 3（**D48: x = √(P/P_th)、Ω_B = 2πf/γ_B**） | T02 V-B3 の表を ±0.02 dB、x_opt = 0.758（1.5°）・0.607（3.9°）を ±0.01、R′−,min = −9.35・−7.17 dB。Fig. 3 から読んだ値と ±0.3 dB。「x ≈ 0.82 で劣化」は目標にしない | T02 V-B3 |
| VT-18 | 9 | Ref B の定義の判別と逆算 | x = 0.999 で R′+ = +21.4 dB（+37 dB なら 2π の誤り）。100 mW で −9.44/+14.95 dB、逆算 x = 0.753 ± 0.002・全損失 0.0708。理想の極限 sin 2θ（−8.67/−12.81 dB） | T02 V-B4・B5・B6、T04 V3 |
| VT-19 | 10 | dither の寄与 | m = 0.0785 rad → 0.0555 rad rms、J0(2m) = 0.99384、Ref C の光源で −5.63 → −5.18 dB、時間領域の PSD が ±0.05 dB | T04 V5、T14 V8、T16 V15 |
| VT-20 | 10 | 全体 約 0.1 rad | dither 以外の残差 ≤ 0.083 rad。制御器のモデル（τ_L 0.1–2 µs、14–16 bit）で 10% 以内に変わらない。QNL のセンサ雑音だけで 3.2 mrad（10 kHz）| T08 V8、T16 V19、T14 V7 |
| VT-21 | 10 | ロック時間が遅いアクチュエータのレンジで決まる | 線形ドリフト 0.1 rad/s・R = 25.96 rad で 259.6 s ± 1 サンプル。ランダムウォークで R²/(2D_ψ) に 5%。熱ランプ 63 m・10⁻⁴ K/s・64.9 rad で 235 s ± 1% | T16 V12、T15 8、T08 V7 |
| VT-22 | 11 | 同軸（共有区間）の相殺 | **区間 × ビームの係数が等しい区間**で、共有区間の雑音がロックの残差にも Θ_k にも出ない（相対 10⁻⁶）。レーザーの共通光路で 10⁻¹²。分散の共通区間で η_t = 1。偏光をオフにしても通る | T10 1、T05 V3、T03 V9、T14 V13、T09 8 |
| VT-23 | 11（追加） | 係数が違う共有区間（OPO キャビティ）では相殺しない | Ref B、x = 0.745、Δ/κ = 10⁻³ で残差 (Δ/κ)(1 − x²)/(2(1 + x²))、パルス OPO の透過ポンプロックで −0.543d | T02 V-B7、T17 V-P9 |
| VT-24 | D04 | Ref D のブラインド前の解析的な検証 | 静的な縞（±3.33 V、8.26 縞）、ロック点（ε 0.05 で 50.02 mrad、V_off 15 mV で 4.505 mrad）、S = V_off で RIN の結合が消える、PER（x = 0.01 で 10.10 mrad）、DAC のリミットサイクル ±6.3 mrad、符号の反転、915B の利得 6.491 rad/V、914 の 6 点、DSC1 の −180° 60 ± 3 kHz、bumpless、離散 PID、レーザー × ΔL（9 m で 0.0277 rad²）、周波数ステップ 0.3077 rad、dither の線 0.308 rad | T15 1–3・6・9–11、T09 5、T13 5・6、T16 V1・V4・V6・V8・V9・V11、T05 V5・V7 |

### 4.2 追加のテスト（§7 に対応しないもの）

| VT | トピック | 内容と目標（要約） | 出典 |
|---|---|---|---|
| VT-30 | OPA | 不確定性 R+R− ≥ 1、符号（θ_sq = φ_p/2 + π/2、10⁻⁴ rad）、最適ポンプの閉じた式、位相不整合（Δkℓ/2 = 0.01 で 2.04 mrad）、結合波方程式（保存 10⁻⁹）、ポンプ感度（+1% で +0.104 dB）、2 つのロック（PD-OPA の 2Ω_m、lock 2 の縮小、結晶温度） | T01 V-T01-4・5・7–10 |
| VT-31 | CW OPO | 行列と閉じた式の一致 10⁻¹²、離調の回転、callable パラメータ、Ref C の線幅（8.7/5.0/2.6 dB） | T02 V-B1・B7・B8・C1 |
| VT-32 | パルス OPA | ガウスの重なり 10⁻⁶、シンプレクティック性 10⁻⁸、低利得の Schmidt、準 CW への連続性、ジッタの平均、熱と位相ロックの共補正、表の補間 10⁻³ | T03 V1–V3・V6・V8・V10・V11 |
| VT-33 | 位相 → squeezing | Monte Carlo と厳密式、側帯波の混合（20 MHz の変調 = 遅い変調）、放射計雑音（0.072–0.077 / 0.093–0.099 dB）、log バイアス（0.033/0.80 dB）、掃引トレースの単調性、窓の式、成分 (3)（1.95×10⁻³ rad²）、指標の不変性・分解の和、高周波の上乗せ（0.01 dB） | T04 V4・V6–V13 |
| VT-34 | レーザー | 遅延の伝達関数（±10%）、コヒーレンス e^{−πΔντ}、サンプル間隔より短い遅延（±2%）、ポート間、Liu 型の静的補償（Θ_1 変化 0、Δ_1 +0.0628 rad）、ロックありの残差、フィルタキャビティ（arg T = 0.1974 rad）、ASE（−8.41 dB、π/(2F)）、CMRR（1.88×10⁻²）、Ref B のピークの傾向 | T05 V1・V2・V4・V8–V14 |
| VT-35 | GNLSE | 分散の符号と大きさ、損失、分散なしの SPM、Gordon の SSFS（BW のみ、3%）、逆向き XPM の伝達関数、格子の端（≤ 10⁻⁶）、既存コードとの照合（任意） | T06 V6–V12 |
| VT-36 | ファイバ量子 | シンプレクティック性 10⁻⁹、物理性、Kerr squeezing の解析解（−4.18/−7.66/−12.54 dB）、SBS 閾値（0.534 W/54.1 mW/34.4 mW）、線形化 vs TW（任意）、自発ブリルアン（4.3×10⁻⁹ /m）、位相のない雑音 | T07 V2–V4・V8–V11 |
| VT-37 | 環境 | 熱係数と符号（+40.60/+43.84 rad）、群遅延 2.457 ps、機械雑音の PSD、長さのスケーリング（γ = 1/2）、相関（0.2S、0.81）、逆向き（2S[1 − sinc]）、熱力学（Dong の値）、1 次の遅れ、トグル、軸の係数 | T08 V1–V6・V9–V12 |
| VT-38 | 偏光 | PM の差動ドリフト（−0.1291 rad/mK）、相反性（10⁻¹²、2 arg(a)）、コネクタ PER（25.61 dB）、分布クロストーク（−30.00 dB）、偏光の重なり、CW とパルスの違い | T09 1–4・6・7 |
| VT-39 | probe の相殺 | 点の雑音（4 sin²）、分布した雑音（f₁ 0.489 MHz、1 kHz で −53 dB）、軸の差（102.5 rad、TEC で 0、90° 融着で 0）、周波数の差（0.0154 rad）、迷光の総量（302.8 光子/10 ns）、ゲート（k = 15 で 0）、SMF の 1/3、相反性、レーザー位相雑音（T_a = T_b で 0）、ショット雑音（6.43×10⁻¹¹）、迷光の squeezing への影響（−7.96 dB） | T10 2–11・13 |
| VT-40 | 振り分け | 漏れがターゲットと直交（共分散 0）、薄まり（−0.124 dB）、スクイーズ角の符号、駆動位相雑音（10.0 mrad）、バイアス（−26.02 dB）、遅延線（L_31 = 63.31 m、310 fs）、重なり（0.990）、立ち上がり、オフで恒等 | T11 2–10 |
| VT-41 | 損失 | Ref B・C の収支、モードの重なり（0.99097 など）、全損失の上限（10 dB・9 mrad で 0.084）、32 ポートの収支（LL で 17.0/17.6%）、エタロン、GBS 入力と BHD の分離 | T12 V5–V7・V9・V11–V13 |
| VT-42 | 検出 | 不平衡の漏れ（−9.45 dB）、周波数応答とディップ、パルス間クロストーク、ショット雑音のスケーリング | T13 3・4・7・8 |
| VT-43 | 誤差信号 | Ref A lock 1 の形（2 f_sh、10⁻⁶）、基底帯域と RF レートの等価性、Ref B lock 1、QNL の誤差信号、radiometer（1/(BT) に 5%、κ_det 1.02/1.32）、Liu 型と統合案の符号（案 B でポンプ係数 1）、DC fringe、オフセット（2.8 mrad） | T14 V1・V2・V4–V6・V10–V12 |
| VT-44 | アクチュエータ | ξ_eo（0.7837）、機械共振の離散化、遅延バッファ（端数）、Moku:Delta（400 kHz で ≤ 30°）、積み重ね = 並列、符号の検査、再引き込み、AOM、TEC、コンパイル実装と純 Python の一致（10⁻¹²）と f_sim の収束 | T16 V2・V3・V5・V7・V13・V14・V16–V18・V20 |
| VT-45 | パルス OPO | CW の極限（10⁻⁴、δt²）、純粋さ、Borealis の K = 1.10、Liu の純度 0.975、長いポンプ → CW、位相の規約、時間ずれ（0.090/κ_c）、裾（8.8e^{−2κ_cT}）、短いポンプ、狭帯域・広帯域の雑音、FSR の櫛 | T17 V-P1・P4–P8・P10–P12・P14・P15 |

## 5. ノート間の整合性チェック

| # | 項目 | 食い違い | 解決 / 調整役への旗 |
|---|---|---|---|
| C1 | Ref A の書誌 | T03・T05 は Crossref で Opt. Express 34, 7958 (2026)、DOI 10.1364/OE.585323。T01・T04・T14 は arXiv:2511.15082（掲載誌 UNVERIFIED）。survey・T10・T12・T13 は UNVERIFIED | Crossref で確認した T03・T05 の値を採り、arXiv 番号を併記した（survey §5）。**旗**: 調整役が DECISIONS / PROCESS_LOG に記録 |
| C2 | Ref C の年 | Crossref の issued は 2025（T05）、survey と他のノートは 2026 | 号の年（2026）のまま、注記（survey §5）。旗（小） |
| C3 | Ref B の DOI | T04 は UNVERIFIED、T02 が Crossref で確認 | T02 を採用 |
| C4 | Liu のページ | T05・T08・T14 は 687–692、T17 は柱の 692 を未確認 | 687–692 を使い、未確認と注記 |
| C5 | Borealis の時間配分 | T11: 測定 65 µs・安定化 35 µs、T17（SI p.22）: ロック光 65 µs・squeezed 光 35 µs | **T17 を採用**（SI で確認） |
| C6 | Borealis の振り分けの損失 | T11（本文）約 15%、T17（SI）透過 > 90%（結合込み） | 定義の差の可能性。素子の損失はパラメータ（1 段 1–4%）として両方を範囲に含める。旗（小） |
| C7 | §7-1 の許容幅 | T01: R′+ 20.4 ± 0.3、T04: 20.56 ± 0.05 | **D48 で解決**: ±0.10 dB（集中損失） |
| C8 | 915B の巻線長と ΔL | T05: 9 m（掲載）、ΔL 名目 10 m・掃引 5–15 m。T08: 9 m、ΔL 0–10 m。T09・T16: 38 × 0.22 = 8.4 m。T15: 9 m、5–15 m。いずれも出典は UNVERIFIED の掲載のスニペット | 巻線 8.4–9 m。**ΔL 名目 9 m、掃引 0–15 m**（和集合、物理的な事前は 8–15 m）。PNH1550R5A2 のリードは両腕とも 0.8 m で相殺（[manual]）。D36 で縞を数える静的な測定が許されている → §7 (a) |
| C9 | PM の dB/dT | T09（一次資料 Ritari 2004）: K_T = −2.05 rad/(m·K)（PM15-U25D で −1.7〜−2.9）、**符号は負**。D11・T08・T10・T15 の仮の値 +2.4 rad/(m·K)（T10 は 0.8–2.4）。T08 は k_T,fast = k_T,slow − 2.4（符号が逆） | **T09 を採用**: 名目 2.05、大きさの掃引 0.3–3.4 rad/(m·K)（k_T の掃引 −(0.25–1.6)×10⁻³ /K と L_b 3–5 mm から [AE-P2]）、dψ_fast/dT = L(C_T + \|K_T\|)。T10 の 51–154 rad/K（63 m）は 129（108–181）rad/K に、T15 の 22 rad/K（9 m）は 17（14–24）rad/K（8.4 m）に置き換わる。**旗**: conventions に K_T の符号を追記 |
| C10 | 熱位相・群遅延の係数 | survey の 44 rad/(m·K)、39 ps/(km·K) は AE | T08 で確認（40–48、37.4–39）。整合 |
| C11 | 群遅延の共補正の割合 | T03 0.21、T08 0.28 | 定数にせずパラメータから計算（§3.4 の C9） |
| C12 | 逆向き伝搬の周波数 | survey・T14: c/(2n_gL) = 1.6 MHz。T10: f₁ = 0.151c/(n_gL) = 0.49 MHz。T08: 2S[1 − sinc(2πfτ)] | 観測量が違う（§3.4 の C10）。ロックの残差には T10 の式。survey を修正済み、T14 §8.4 の値は古い |
| C13 | Rayleigh / Brillouin の比 | T07: 約 12（survey の S α_s = 5.3×10⁻⁸ /m）、T10: 14（Corning −82 dB から b_R = 6.2×10⁻⁸ /m） | T10 の b_R（データシート）を採用。T07 の Rayleigh の列は 約 15% 低い |
| C14 | probe のショット雑音の係数 | survey E73: hν/(ηP)、T07: hν/P、T10: κ_het hν/(2ηP)、T12: κ_sn = 1/2、T14: 方式ごと | T14 の方式ごとの式を採用（survey 修正済み） |
| C15 | Ref B の G = 18.7 | T02: 使わない（x = √(P/P_th) = 0.745）。T14 は G = 18.7 から x = 0.769 で A± と Ref B の lock 1 の値を出した | **D48・T02 に従う**。T14 V3 の Ref B の値（13.1 mrad）は x = 0.745 で 14.7 mrad になる [AE-P2: tan ε = 0.146 × tan 0.1]。旗: Phase 4 で T14 の Ref B の数値を作り直す |
| C16 | 915B の物理的な伸びと読んだ場合の倍率 | T15: ×1.13（7.3 rad/V）、T16: ×1.16（7.53 rad/V） | conventions §4 は ΔOPL 読み（6.49 rad/V）。代替の倍率は ξ_eo·n_eff をパラメータで計算 |
| C17 | 914 のモデル、915B の共振範囲 | §3.4 の C7・C8 | T16 の 2 極、範囲は和集合 |
| C18 | DSC1 の遅延 | T15 5.67 µs、T16 5.7 µs、survey 約 8 µs | 整合（5.7 µs、上限 8.3 µs） |
| C19 | PM ファイバの損失 | T04 は 1 dB/km の仮置き、T06 は 0.5（Coherent）/ 1.0（Thorlabs の表）、T12 は Corning PM15-U25D の最大 0.5、T11 は SMF 0.2 で代用 | PM15-U25D の ≤ 0.5 dB/km を参照、0.2–1.0 dB/km を掃引。T04 §9.3 の PM 630 m（16.1%）は 0.5 dB/km なら 約 7%（T12 §4.3）で、表の「不可能」の一部は変わる |
| C20 | 経路の組み方の 3 つのシナリオ（D38 の −1.5 / −4.7 / −7.3 dB） | T12 の表は振り分けを 1 段（IL_sw = 0.09 dB）として計算した。T11 の二分木は 5 段で、Pockels 1 段 1–4% | LL を 5 段 × 2%（0.44 dB）にすると、ポートの IL 約 0.94 dB、全損失 約 23%、Ref A の光源・9 mrad で 約 −6.1 dB [AE-P2]。**旗**: D38 の参照値に段数の仮定を明記 |
| C21 | Ref D の部品の資料 | T09・T12・T15 は PNH1550R5A2・Opneti 偏光子の資料を読めず仮の値（PER 20–25 dB、κ = 0.5、T_opt = 1） | D37 で追加された PDF を Phase 2 で読んだ: PNH1550R5A2 は 50:50 ± 1.5%、PER ≥ 25 dB、過剰損失 < 0.3 dB typ、ORL ≥ 60 dB、リード 0.8 m。Opneti は IL 0.3 dB typ（コネクタで +0.3）、ER 30 dB typ・≥ 28 dB（コネクタで −2 dB）。Ref D のプリセットでは入力 PER x ≤ 2.5×10⁻³（±2.5 mrad の寄生縞、T15 式 (7)）、κ_2 = 0.5 ± 0.015 [AE-P2] |
| C22 | DX1 の h0 | D35 は SC 型（instantaneous < 100 kHz → h0 ≤ 3.2×10⁴）、D49 の参照値 h0 = 25 kHz/π（N・GC 型の 25 kHz） | D49 を参照値、8×10³–3.2×10⁴ Hz²/Hz を掃引。旗（小） |
| C23 | Ref A の lock 1 の SNR と 14 → 9 mrad | Ref A の説明と、T01・T14 の導出（lock 1 の誤差は A−/A+ = 1/123 で縮む）が食い違う | 未解決。Phase 5 で VT-15 によって確認 |
| C24 | Ref A Fig. 3 の極小の深さ | VBW 100 Hz なら極小は +1.8 dB まで埋まるはずが、図は −10 dB（T04） | 未解決。掃引トレースは単調性だけを検証（VT-33） |
| C25 | Ref A の反 squeezing | 単一モードモデル 20.56 dB、測定 20.2 ± 0.2 dB | T01 の分布損失モデルで 0.27 dB 説明できる。合わせ込まず記録（D48） |
| C26 | L_ic(x) の変数 | T02: 論文の定義（√）と、実測 0.0038 に近い x² の読み方。R′−,min で 約 0.1 dB の差 | 検証は √（D48）、x² は変種 |
| C27 | 振り分けの漏れ | survey E48「相関も生む」、T11「ターゲットには損失」、T17「OPO の裾はコヒーレントな混合」 | T11 は固定遅延・単一光源・τ_mode ≪ T_rep で正しい。T11 の「τ_p ≪ T_rep は常に成り立つ」は OPO では成り立たない（N11） |
| C28 | 上流の効率 | T12 の η_up = 0.951（導波路 OPA）、T17 の OPO（Fb 約 0.76、Fl 0.92–0.94） | プリセットごとに `loss_budget` の要素を変える。Fb で GBS 入力 10 dB は不可能 |
| C29 | SPM の大きさ | survey 8 rad、T06 0.05–0.3 rad（100 MHz/63 m）、T07 は survey の値を TW が要る理由に挙げた | T06 を採用。TW は Fw の 10 MHz・633 m の領域でだけ要確認（N18） |
| C30 | GAWBS の最低周波数 | survey 20 MHz、T07 30.5 MHz（R0,1）。TR2,m は未計算 | R0,m は T07。TR2,m の式は Phase 4 で原典から |
| C31 | Ref C の文献番号 | Ref C 本文の [27]（LO をチョップ、1.6 km）と文献リストの [27]（20 km）が合わない（T07） | 旗（小、Ref C の再現には影響なし） |
| C32 | Tohermes 2024 の掲載誌 | Ref C は New J. Phys. 27, 094506 (2024) と書くが巻と年が合わない（T02） | UNVERIFIED のまま |
| C33 | SMF・直交偏光の誤差 | T09: 2·arg(a)（予測できない）、T10: 0.03 rad/K（根拠は弱い） | T09 の一般式を採用、T10 の値は 1 つの見積もり |
| C34 | Moku:Pro | T16: プリセット不要、D25: 足りなければ作る | T16 の結論で解決（VT-20）。旗: D25 に追記 |

## 6. Phase 3 で決める性能・アーキテクチャの論点（判断はしない）

| 論点 | 選択肢 | 根拠（ノート） |
|---|---|---|
| 時間ループの実行方式（prompt §5.4、D49） | (c) コンパイルした 1 サンプルごとのループ（numba など、新しい依存）＋ (a) MHz 以上の搬送波は基底帯域・統計モデル ＋ 遅い段の多重レート。(b) チャンクは開ループ（雑音の生成、ESA の後処理）だけ | T16 §5: (b) は「チャンク長 ≤ ループ遅延」が条件で、Moku 系（≲ 0.21 µs）では 1 サンプルになる。Ref C・D のロック時間は数分（f_sim 1 MHz で 6×10⁷ ステップ）。Python で 1 ステップ・1 ブロック 5–20 µs なら 10 ブロックで 1 時間以上、numba で 50–200 ns なら 3–12 s（AE、未計測）。T15 は DSC1 に合わせ Δt = 2 µs、T14 は複素包絡線で RF を解析的に、T05 はレーザーの短い遅延を補間 ＋ 帯域外をオフライン |
| 誤差信号の RF の扱い | 解析的な復調（基底帯域）と RF レート（f_sim 8–10 MHz なら 1 MHz を直接）の両方。互いの検証 | T14 §12.4、T16 §5 |
| 量子雑音層（prompt §5.1） | 線形化（古典解の Jacobian ＋ 共分散、随伴の伝搬 1 回 / LO モード）を基本、truncated Wigner（N_traj ≈ 10⁴）は検証用 | T07 §2.5: squeezed vacuum の経路は φ_NL ~ 10⁻⁷ rad で線形、雑音はすべて加法的なガウス、決定論的で §7-4 を機械精度で検証。TW は明るい LO の φ_NL が大きいときだけ |
| GNLSE のソルバ | 自前の RK4IP（numpy、約 60 行）、h̃_R(Ω) は解析式。既存コードは照合だけ | T06 §4: gnlse-python 2.0.0 と PyNLO 0.1.2 は NumPy 2 で動かない（np.math、np.trapz）、gnlse-python は RK45 で §7-5 の 16× テストに使えない |
| パルス導波路 OPA | 周波数領域の 2N×2N 線形 ODE ＋ Takagi 分解、格子 512–3200 点 | T03 §2.1–2.2。100 fs では K が格子に依存 |
| パルス OPO | 時間領域の量子 Langevin、2N×2N の共分散（伝搬子の積）、δt = 1/(40–160 κ_c)、初期の真空の項を必ず含める | T17 §2.2、V-P1・P3 |
| CW OPO | 側帯波ごとの 2×2 行列（numpy の inv） | T02 §2.1 |
| ESA の模擬 | 時間領域の基底帯域（RBW の 4 倍のサンプル、Ref A のトレース 1 本 6×10⁵ 点）で軽い | T04 §5、§13 |
| 量子雑音ロックのセンサ雑音 | 等価なガウス雑音（S_P̂ = 2∫S_x² df）と RF レートの Monte Carlo の両方 | T14 §6.3、V6 |
| 遅延 | 整数のリングバッファ ＋ 端数（線形補間 / 1 次 Thiran）。放出時刻の遅延は φ_L の補間 | T05 §3.5、T16 §5 |
| 逆向き伝搬 | 区間を N_sub 個の小区間に分け、位置ごとの遅延で読む。N_sub は最高周波数で f τ_L/N_sub ≪ 1 | T08 §2.7、T10 §6.1 |
| 多重レート | 熱（10²–10³ s）は 10–100 Hz で積分して補間、機械（kHz）は FFT で色付けした雑音をチャンクで生成、TEC・遅い PID は間引いたレート | T08 §7.3、T16 §5 |
| 指標 (ii) の分解 | 雑音源を切り替えた再実行（4 回）と線形性の検査。飽和があると無効 | T04 §4、V12 |
| 表の補間 | Δt の 1 次元（＋ φ2）、OPO は (Δ_s, Δ_p, x_ss)、チャンクごとに更新で足りる | T03 §8.2、T17 §8 |
| ポート数 | 1 ポートで L_k を掃引、2 ポート組、32 ポートは Phase 6 | D13 |
| 抽象制御器インターフェース | `configure`・`set_params`・`engage`・`acquire`・`measure_transfer_function`・`status`・記録。機種への写像（DSC1 の K_i 単位と符号、Moku の交差周波数） | T16 §4.6 |
| 決定論と乱数 | GNLSE と線形化は決定論的、雑音源は seed を記録（prompt §6） | T06、T07 |

## 7. 残っている問い

### 7.1 (a) ユーザだけが答えられる事実

1. **D43 の解釈**: 「Xanadu 2022 の Supplementary と同等の OPA」は SI の PPKTP ヘミリシック OPO（パルス励起）を指すか（DECISIONS Open）。
2. **将来の光源の値**: 波長、f_rep、パルス幅または OPO の線幅・T_p、平均パワー、N_b。OPO を複数台使う予定があるか（Liu は 4 台。複数ならポンプ位相の絶対的なロックが要る）。
3. **Ref D の静的な値（D36 で測ってよい量）**: 腕の長さの差 ΔL（縞を数えた値、または測る予定）、腕 2 のファイバ長、915B のピグテールの長さ、915B（8669-1）の挿入損失。曲げ半径が分かるか。
4. **DSC1 の設定と単位**: K_p・K_i・K_d・S・O・CV 制限（D35 で未提供、掃引で代用）。K_i の単位（マニュアル p.12 と p.20 の食い違い）を Thorlabs に問い合わせられるか。
5. **914 の入力インピーダンス、DSC1 の出力に T 分岐したオシロスコープの入力インピーダンス**（D35 では 1 MΩ を仮定）。
6. **未入手の資料**: Optoquest 偏波コントローラ、Thorlabs PM パッチケーブル（D37）。Ref A の導波路（NTT）の 772 nm の損失 a_p と ΔT_FWHM の資料があるか（T01）。
7. **GBS の要求**: 許される迷光（ゲートあたりのクリック数 / モード内の光子数）、検出器が時間ビンを分けられるか（ゲート・ジッタ）、GBS が毎フレーム動くか（D46 の仮定）、モード数 32 以外の要求（D16・D28・D46）。
8. **Moku:Delta の PID 計器で 10 Vpp の出力レンジが使えるか**（使えなければ 915B のレンジが 10–20%、T16）。
9. **PM15-U25D の K_T の実測を許すか**（文献値に 6 倍の幅。Ref D の予測の入力だが、D36 の「静的な量」に入るか）。
10. 将来のハードウェアの具体的な値（EOM、BHD、制御器、D15）。

### 7.2 (b) 設計の選択（D34 により掃引するパラメータ）

| パラメータ | 範囲・値 | 出典 |
|---|---|---|
| 光源の種類 | `waveguide_opa` / `borealis_opo` / `liu_opo` / `cw_opo` | D07、D43、D50 |
| 導波路 OPA の長さ ℓ | 10–45 mm、GVM を設計した導波路 | T01、T03、D44 |
| パルス幅 τ_p（導波路） | 100 fs – 10 ps | D07 |
| OPO の (Δν_c, T_p κ_c, κ_p/κ_c, x_ss, ρ) | 30–300 MHz、0.1–5、{非共振, 1–2.5}、0.5–4.5、0.90–0.995 | D50、T17 |
| f_rep | 1 MHz – 1 GHz（OPO は裾の上限まで） | D07、T11、T17 |
| 振り分け・N_b・配置 | §2.3、N_b = 1–100、木 / 鎖 | D39 |
| L_0、δL_k（分散）、δτ_k（タイミングの合わせ残り） | 0–数 m、0–1 m、0–1 ps（Fw）/ 0–500 ps（OPO） | T03、T11 |
| ファイバの種類 | PM / SMF、軸の入れ替え位置 | D32 |
| 経路の組み方 | NF / MF / LL、接続数 N_j = 1–6、IL_j 0.03–0.3 dB | D38、T12 |
| probe | P̄ 0.1 nW–1 mW、デューティ 10⁻³–1、t_offset、偏光、Δν_pr 0–1 THz、同じ / 別のレーザー、分離、検出方式、検出点、アクチュエータ位置 | D40、T10、T14 |
| ロックのトポロジー | §2.4 | D14、D30、D42、D50 |
| 遅延の合わせ T_a − T_b | 0–100 ns | T10 §3.8 |
| 大レンジ | TEC の L_h 1–10 m、ΔT_max ±1–10 K、τ_th 1–100 s、2π リセット on/off | D41、T16 |
| 制御器 | `dsc1` / `moku_delta` / `fpga_generic`（τ_L 0.1–2 µs、12–16 bit）、PID ゲイン、anti-windup none/clamp/back_calculation | T16 |
| LO | 配り方、波形、LO 配送のロック on/off、P_LO | D09、D12、D45、D50 |
| BHD | パルス分解 / ゲート / 積分、B/f_rep 0.5–10、QE 0.97–0.99、クリアランス 10–30 dB、CMRR 25–40 dB | D45、T13 |
| FSR フィルタ | なし / FBG の組 / MZI の段、透過 0.80–0.99 | D50 |
| スペクトルフィルタ・分散補償 | 帯域、損失、on/off | D44 |
| OPO のロック光・ポンプのロック | §2.3 | D50 |
| ポートのファイバの置き方 | 1 つのスプール / 別々、squeezed と LO_k を同じ束に入れるか（ρ_0 0–0.999、f_ρ 0.01 Hz–1 kHz） | T08、D09 |
| 目標 squeezing | 3–10 dB（GBS 入力と BHD） | D16、D46 |
| 位相の予算の配分 σ_t | Phase 6 | D49 |

### 7.3 (b′) 物理量の不確かさ（掃引）

k_T 40–48 rad/(m·K)、κ_τ 37.4–39 ps/(km·K)、PM の K_T 0.3–3.4 rad/(m·K)、ξ_eo 0.78–0.80、915B の共振 5–60 kHz・Q 3–100、914 のモデル、h0 8×10³–3.2×10⁴ Hz²/Hz、h₋₁ 0–5×10⁸ Hz²、レーザーのドリフト ≤ 0.3 GHz/h、r_T 10⁻⁶–10⁻³ K/s、τ_th 1–3000 s、h_ref 10⁻³–4 rad²·Hz、γ ∈ {1, 2}、f_c 0.1–10 kHz、D_pol 10⁻⁸–10⁻² rad²/s、ε_pol、b_R の 1–3 倍（PM）、KTP の Δ(dn/dT) 10⁻⁶–10⁻⁵ /K、PPKTP の ∂Δk/∂T、Raman のモデル（BW / LA）、ESA の平均の種類、Ref A のフィルタキャビティの線幅と ASE、PD-OPA の利得、Ref C の δL・線幅、915B の曲げ損失 0.01–2 dB。

## 8. 調整役への提案（DECISIONS・conventions の変更案。このノートでは編集していない）

1. **Ref A の書誌**を Opt. Express 34, 7958 (2026)、DOI 10.1364/OE.585323（arXiv:2511.15082）として DECISIONS か PROCESS_LOG に記録する（§5 の C1）。
2. **conventions に PM の 2 軸の位相の温度係数を追記**: Δψ_b = ψ_slow − ψ_fast、K_T = ∂β_b/∂T < 0（名目 −2.05 rad/(m·K)、T09）、dψ_fast/dT = L(C_T + \|K_T\|)。D11 の 2.4 rad/(m·K) は置き換え（C9）。
3. **D38 の参照値**（−1.5 / −4.7 / −7.3 dB）に「振り分け 1 段」の仮定を明記し、5 段の Pockels 木の値（LL で 約 −6.1 dB [AE-P2]）を併記する（C20）。
4. **D48 の §7-11**: 「区間 × ビームの係数が等しい共有区間」に限ることを明記し、係数が違う OPO 区間の追加テスト（VT-23）を加える。§7-8 は 2 つの定義（VT-16 の (a)(b)）を合格条件にする。
5. **D47 の (3)**: 帯域内のレーザー項 16 sin²(πf(T_a − T_b)) sin²(πf(τ_k − τ_0)) S_φ と、遅延の合わせ T_a − T_b をパラメータに加える（N04、T10 §3.8）。
6. **インターフェースの決定**（§3.4）: C1（偏光の作り方）、C2（帯域の勘定の 1 か所化）、C3（損失の登録表）、C4（光源の共通出力）、C13（TEC は温度を書く）を Phase 3 の PLAN の確認（grill）の対象にする。C7・C8 は T16 の 2 極と和集合の範囲を採る。
7. **δL_k を 2 つに分ける**（分散の δL_k とタイミングの δτ_k、T03 §9-1）を DECISIONS に入れる。
8. **D25 に追記**: Moku:Pro のプリセットは不要（T16、VT-20 で確認）。
9. **Ref D のプリセットに D37 の資料の値**（PNH1550R5A2、Opneti）を入れる（C21）。
10. **D49 の DX1 の h0** は D35 の SC 型と合わせて 8×10³–3.2×10⁴ Hz²/Hz を掃引と明記（C22）。
11. **Ref B**: G = 18.7 を使わないことに合わせ、T14 の Ref B の数値（x = 0.769 から出した A±、lock 1 の値、V3 の 13.1 mrad）を Phase 4 で x = 0.745 から作り直す（C15）。
12. **このノートのトリアージ（kept 94 / opt 9）とプリセットの名前**（§2）を Phase 3 の PLAN の出発点として承認するか判断する。
13. Open に残すもの: §7.1 の (a) の 10 件、DSC1 の K_i の単位、Ref A の lock 1 の SNR と 9 mrad の関係（C23）、Ref A Fig. 3 の極小（C24）。
