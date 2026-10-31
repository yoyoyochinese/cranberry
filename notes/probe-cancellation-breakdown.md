# T10 probe と squeezed 光の相殺の崩れ（probe-cancellation-breakdown）

> **訂正（2026-11-01、DECISIONS D77-1）:** §0-6・式 (12)/(13)・§8 の「T_a = T_b なら帯域内のレーザー項は 0」は誤り。正しい伝達は WP-B8 の H（`reviews/phase5_M4_pol_demux.md`）で、ポートに依らない遅延では 1 ポートしか打ち消せない。帯域内の残差は有限のループ利得 1/(1+L) で決まり、線幅の要求は D31-3 のまま（10 mrad に 51 Hz @ 63 m、5.1 Hz @ 630 m）。


- Phase: 1（deep dive）、日付: 2026-10-28
- 担当する効果: E41–E46（`notes/survey.md` §2.6）、E35 の T10 部分、E30 の probe 部分
- 前提: DECISIONS D08–D13、D18-5、D28、D30–D32。規約は `notes/conventions.md`（特に §4、§9、§10）
- 依存: T07（Brillouin・Raman の係数）、T08（環境雑音の PSD、熱係数）、T09（PM の複屈折の温度係数、PER、SMF の偏波）。
  どれもまだ値が出ていないので、このノートでは**仮の値**を使い、そのつど [AE]（agent estimate、未検証）または UNVERIFIED と書く。
  T07–T09 の値が出たら §3 の数値を差し替えること（式はそのまま使える）。
- 値の印: [Liu]・[Ref A] などは原典から読んだ値、[AE] はこの担当の概算、UNVERIFIED は原典で確認できていない値・文献情報。

## 0. 結論の要約

1. **直交偏光の probe を PM ファイバで使うと、ロックは squeezed 光の位相を 2 軸の位相差のドリフトの分だけ間違える**。
   差は (2π/λ0)(dB/dT)·L·ΔT で、dB/dT = 2–6×10⁻⁷ /K [AE/UNVERIFIED] なら 0.8–2.4 rad/(m·K)、63 m で 51–154 rad/K。
   これは probe の誤差信号（in-loop）には一切現れない。Liu et al. の「3 nm、10 nm の精度」は誤差信号からの見積もりで、この差を含まない可能性がある（§2.4）。
   対策: 軸の入れ替え（90° 融着）、温度アクチュエータ（TEC）での補正、SMF、同じ偏光＋別の分離手段（§4）。
2. **逆向き伝搬による相関の崩れ**は、ファイバに一様に分布した雑音・アクチュエータが端にある場合、
   残差 / 補正前 = 2[1 − sinc(4π f n_g L/c)] で、**1 になる周波数は f₁ = 0.151·c/(n_g L)**（63 m で 0.49 MHz。survey の c/(2 n_g L) = 1.6 MHz より約 3 倍低い）。
   1 kHz で −53 dB、10 kHz で −33 dB（63 m）。音響帯の雑音には効かないが、probe がファイバを通る時間（63 m で 0.31 µs、630 m で 3.1 µs）がループ遅延に加わる。
3. **迷光（probe → GBS へ向かう光子）**: 同じ偏光・同じ周波数・CW 1 mW・63 m で、Rayleigh 後方散乱だけで 10 ns 窓あたり **約 3.0×10² 光子**（survey の 2.6×10² を確認。
   後方散乱係数は Corning のデータシートの −82 dB/ns から取った）。自発 Brillouin（±10.8 GHz）がその約 1/14 [AE]。
   PM の直交偏光＋PBS で 20–30 dB 減る。**SMF の直交偏光では平均で 1/3（−4.8 dB）しか減らず、時間で揺らぐ**（往復の Jones 行列 MᵀM による、§3.4）。
   パルス化しても**分布した Rayleigh は平均パワーで決まり、時間をずらしても消えない**。消えるのは往復時間 2n_g L/c がフレームの空き時間より短いポートだけ（§3.4.4）。
4. **ロックに要る probe のパワーはずっと小さい**: ショット雑音で決まる位相雑音から、迷光の数と位相誤差の積は
   N_stray·σ_lock² = ε·b_R·L·B_loop·T_gate/(2η_det T_ret)（probe のパワーによらない）。63 m・10 kHz・1 ns・同じ偏光で 4.9×10⁻¹¹ rad²。
   σ_lock = 1 mrad なら検出 0.8 nW で足り、迷光は 1 ns 窓あたり 5×10⁻⁵ 光子。**1 mW という仮定は 6 桁大きい**。実際の下限は検出器の電気雑音・DC オフセットなどの技術雑音で決まる（T13/T14）。
5. **入れ出しの素子の損失**: 10 dB の squeezing は 0.3 dB の損失で 8.0 dB、0.6 dB で 6.7 dB、1.2 dB で 5.0 dB に落ちる。
   ファイバ型の PBS・サーキュレータを 1 ポートに 2 個入れると致命的。Liu は自由空間の PBS を使っている。
6. **Ref A の phase-detection OPA との統合**: Ref A の lock 1 はポンプ–probe（seed）の位相、Liu 型はポートの光路の位相を見る。
   単一の OPA を時間で分ける構成（D08）では、ポンプ光路の位相雑音もポートごとに違う時刻でサンプルされるので、Ref A の lock 1 は (ii) にも効く。
   残る非共通の光路は §5 の 5 つ。また、Liu 型のロックは帯域内でレーザー位相雑音の時間差の項（D31 の 3）を
   遅延差 τ_k − τ_0 から遅延差 T_a − T_b（probe と squeezed 光の「レーザーからの経路の遅延」の差）に置き換える（§3.8）。D31 の「ロックでは補正できない」は帯域外でのみ正しい。

## 1. 記号（symbols）

- `λ0`, `ν0 = c/λ0`, `ω0`: carrier wavelength / frequency. Future config uses λ0 = 1545.32 nm (D07 placeholder, Ref A). `hν0 = 1.285×10⁻¹⁹ J`.
- `c`: speed of light in vacuum. `n_g`: group index of the fiber (1.468 assumed, survey; UNVERIFIED, T06). `v_g = c/n_g`.
- `L`, `L_k`: length of the common (probe + squeezed) fiber segment; of port k (D08, `L_k = k·c·T_rep/n_g`).
- `x ∈ [0, L]`: position along the common segment, x = 0 at the upstream (switch / probe-extraction) end, x = L at the GBS (probe-injection) end.
- `x_a`: position of the phase actuator (PZT stretcher / PZT mirror) on the common segment.
- `τ_x = n_g x/c`: one-way group delay from x = 0 to x. `τ_L = n_g L/c`.
- `d ∈ {+1, −1}`: propagation direction (conventions §9). `a ∈ {s, f}` (PM slow/fast) or `{x, y}` (SMF).
- `ψ_{j,a}(t)`: optical phase acquired in segment j on axis a (conventions §4, A → A·e^{+iψ}).
- `φ_L(t)`: laser (seed) phase at emission time t. `φ_p`: pump phase. `θ_sq,k`: squeezing angle of port k at the GBS input (`θ_sq = φ_p/2 + π/2`, conventions §4).
- `φ_pr,k`, `φ_ref`: probe phase of port k at its detector, reference-beam phase at the same detector.
- `Θ_k = θ_sq,k − θ_sq,0`: metric (ii) (conventions §10, D31). `Δ_k = φ_LO,k − θ_sq,k`: metric (i).
- `T_a`: group delay from laser emission to the entrance (x = 0) of a port's common segment along the squeezed path (seed → SHG/pump → OPA → switch).
- `T_b`: delay from laser emission to the probe detector, excluding the common segment (probe distribution + injection + extraction + detector path).
- `τ_k = n_g L_k/c`: delay of port k's common segment.
- `B`: modal birefringence of PM fiber (n_s − n_f). `dB/dT`: its temperature coefficient [1/K].
- `κ_{j,a}`: phase response of noise source / actuator j on axis a (e.g. rad/K, rad/V).
- `ρ_act = κ_{A,s}/κ_{A,f}`: axis ratio of the actuator A.
- `Δν_pr`: probe optical frequency offset from the squeezed carrier [Hz]. `dτ_g/dT`: thermal coefficient of the group delay per unit length [s/(m·K)].
- `P_pr(t)`, `P̄_pr`, `P_pk`: probe power launched into the common segment (instantaneous, average, peak). `D = P̄_pr/P_pk`: duty cycle.
- `Φ = P/(hν0)`: photon flux [1/s].
- `b_R = S·α_s`: Rayleigh backscatter coefficient [1/m] (S: capture fraction, α_s: Rayleigh loss coefficient [1/m]).
- `b_B`: spontaneous Brillouin backscatter coefficient [1/m] (Stokes + anti-Stokes). `ν_B`: Brillouin shift (≈10.8 GHz at 1550 nm in SMF, UNVERIFIED). `Δν_B`: Brillouin linewidth (FWHM).
- `ρ_Ram(Δν)`: spontaneous Raman backscatter coefficient per unit length and bandwidth [1/(m·nm)].
- `α`: total power attenuation coefficient [1/m]. `L_bs = (1 − e^{−2αL})/(2α)`: effective backscatter length (≈ L when αL ≪ 1).
- `R_i`, `z_i`: power reflectance and position of discrete reflector i (connector, collimator facet, splice).
- `ε_pol`: fraction of backscattered probe power that ends up in the squeezed polarization at the GBS input (after the injection PBS).
- `ε_xt`: polarization cross-talk (fiber h-parameter × length + components + depolarized scattering). `T_⊥`: leakage of the wrong polarization through the PBS transmission port (1/extinction ratio).
- `ε_AOM`: power extinction of the AOM in the off state (leakage / on-state).
- `F_GBS(Δν)`: power transmission of the GBS-side spectral filter at offset Δν. `I_WDM(Δν)`: isolation (power transmission) of a WDM separating probe and squeezed.
- `T_gate`: detector integration window per time slot. `τ_mode`: duration of the squeezed temporal mode.
- `T_f = 32·T_rep`: GBS frame period (all 32 ports arrive in the same slot, D08).
- `N_stray`: mean number of stray (probe-derived) photons reaching the GBS input in one window T_gate. `n_mode`: mean stray photon number in the squeezed temporal mode.
- `η_det`, `T_ret`: probe detector quantum efficiency; transmission from the probe launch point to the detector.
- `B_loop`: lock bandwidth (integration upper limit for in-loop sensor noise). `σ_lock`: RMS phase error from probe shot noise.
- `γ`: fiber nonlinear coefficient [1/(W·m)]. `S_x(f)`: one-sided PSD (conventions §7).

## 2. 構成と区間のモデル（conventions §9、D30、D31）

### 2.1 Liu et al. 2026 の方式（原典から読んだこと）

Liu et al. Methods「Phase stabilization of the system」と Extended Data Fig. 4b から [Liu]:

- 遅延ループの配列（第 1 と第 2 の干渉計の間、第 2 と第 3 の間）の各空間モードの光路に PZT 付きミラーを置く。
- probe は「周波数は squeezed 光と同じ、偏光は直交」で、逆向きに伝搬する。16 本に分けられ、各空間モードに PBS で入れる（本文は「第 2 の干渉計の入口で」、図 4b では AOM → PBS が第 2 の干渉計の squeezed 光の出口側に描かれている。本文と図の位置の書き方が違って見えるが、逆向きの probe にとっては同じ場所という意味かもしれない。未確認）。
- 迷光の計数を減らすため、probe を AOM でパルスにし、squeezed 光のパルス列から時間をずらす。
- 第 1 の干渉計の出口で、probe を別の PBS で取り出し、参照光（1550 nm reference beam）と干渉させて（ファイバ PBS と PD1・PD2 の差）誤差信号を作り、PZT を動かす。
- 長いファイバの温度ドリフトが PZT 1 台の範囲を超えるので、TEC のフィードバックで遅いドリフトを抑える（図 4b で TEC は遅延ループのファイバについている）。
- 精度の見積もりは、遅延ループの第 1 の配列で 3 nm、第 2 の配列で 10 nm。系全体は 1550 nm で約 λ/200（Supplementary Fig. 18、本文のみ確認、補足資料は未読）。
- 時間のビンの間隔 τ = 50 ns、遅延は第 1 の配列で t … t + 15τ、第 2 の配列で t … t + 15·16τ。
  ファイバで作るなら最長は 750 ns ≈ 153 m と 12 µs ≈ 2.45 km（n_g = 1.468 とした [AE]）。
- **書かれていないこと**: ファイバの種類（PM か SMF か）、probe のパワー、パルス幅、時間のずらし方、TEC の制御信号（温度センサか、PZT の出力の逃がしか）、精度の測り方（in-loop か out-of-loop か）。
  補足資料にあるかもしれない（未読、UNVERIFIED）。

### 2.2 ポート k の区間（Liu 型、D30）

| 区間 id | 中身 | 通る光（向き d、軸 a） | D31 の分類 |
|---|---|---|---|
| `src` | レーザー → SHG → ポンプ光路 → OPA | ポンプ（+1） | ポンプ位相 φ_p（Ref A の lock 1 の対象） |
| `up` | OPA → 振り分け素子 → 取り出し PBS | squeezed 光（+1, a_sq） | squeezed 光だけが通る（ロックされない） |
| `com_k` | 取り出し PBS → アクチュエータ → ファイバ L_k → 入れる PBS | squeezed 光（+1, a_sq）と probe（−1, a_pr） | (1) ロックされる区間 |
| `down_k` | 入れる PBS → GBS 入力 | squeezed 光（+1） | squeezed 光だけ（ロックされない） |
| `pdist_k` | probe の共通の分岐点 → 入れる PBS | probe（−1） | (2) 共通点の揺れ |
| `pext_k` | 取り出し PBS → probe 検出器 | probe（−1） | (2) |
| `rdist_k` | 参照光の配送点 → probe 検出器 | 参照光 | (2) |

`up` と `down_k` は、ロックの誤差信号に現れず、squeezed 光の位相にはそのまま入る。特に `up` は全ポートで共通だが、ポートごとに違う時刻でサンプルされる（ポート k のパルスは τ_k − τ_0 だけ早く出る）ので、(ii) では 4 sin²(π f (τ_k − τ_0))·S_ψup(f) として残る。短い自由空間なら小さいが、振り分け素子の位相雑音（T11）はここに入る。

### 2.3 ロックの式と指標 (ii) の分解

ビームの位相 = 出射時刻での光源の位相 + 通る区間の ψ の和（conventions §9）。ロックの帯域内でループの利得が十分大きいとき、

- probe: φ_pr,k(t) = φ_L(t − T_b − τ_k) + ψ_pdist,k + ψ_pext,k + Σ_x ψ_com,k,a_pr(x, t − τ_x) + a_k(t − τ_{x_a})
- 参照光: φ_ref(t) = φ_L(t − T_ref) + ψ_rdist,k
- ロック: φ_pr,k − φ_ref = c_k（設定値）→ アクチュエータ a_k が決まる。
- squeezed 光（GBS 入力、時刻 t）: θ_sq,k(t) = φ_L(t − T_a − τ_k) + ψ_src/2 + ψ_up + Σ_x ψ_com,k,a_sq(x, t − τ_L + τ_x) + a_k(t − τ_L + τ_{x_a}) + ψ_down,k + π/2

（ポンプは SHG なので φ_p = 2φ_L + ψ_src、θ_sq = φ_p/2 + π/2（conventions §4）。squeezed 光が区間 ψ を通ると、場 A → A e^{iψ} なのでスクイーズ角も ψ だけ回る。
各項はその光がその区間を通った時刻で評価する。記法を軽くするため、ψ_src・ψ_up・ψ_down の時刻の引数は省いた。）

これを Θ_k = θ_sq,k − θ_sq,0 に入れると、D31 の 3 つに分かれる:

1. ロックされる区間の残差 = 制御の残差（T16）+ §3.1 の軸の差 + §3.2 の時刻のずれ + §3.3 の周波数の差 + probe のショット雑音（§3.7）+ XPM（§3.5）
2. 共通点の揺れ = −(ψ_pdist,k + ψ_pext,k − ψ_rdist,k) + (ψ_down,k) − （ポート 0 の同じ量）、および ψ_up の時刻差サンプル
3. レーザー位相雑音: 帯域外では φ_L(t − T_a − τ_k) − φ_L(t − T_a − τ_0)。帯域内では §3.8 のとおり別の形になる。

## 3. 効果ごとの式・範囲・概算

ここでの「63 m」は f_rep = 100 MHz の L_31（survey §3）、「630 m」は 10 MHz の L_31、「6.3 m」は 1 GHz の L_31。Ref D には probe も squeezed 光もないので、このトピックは Ref D には効かない（PM の PER による可視度の変化は T09/T15）。

### 3.1 偏光の軸の違い（E41、D32）

**式**: 雑音源 j（温度、歪み、音響）と、アクチュエータ A の応答を軸ごとに κ_{j,a} とする。probe の軸 a_pr でロックすると、squeezed 光の軸 a_sq に残る位相は

  δθ_res = Σ_j [κ_{j,a_sq} − κ_{j,a_pr}·(κ_{A,a_sq}/κ_{A,a_pr})]·n_j  (1)

（n_j: 雑音源の大きさ。帯域内、ループの利得が十分大きいとき。）

- PM・直交偏光・PZT のストレッチャ（歪み）: 歪みの応答の軸の比 ρ_act ≈ 1（2 軸の差は光弾性の差の分だけ。UNVERIFIED、T09/T16）。温度の雑音 n = ΔT に対しては
  δθ_res ≈ (2π/λ0)(dB/dT)·L·ΔT  (2)
  dB/dT は応力付与型（PANDA、bow-tie）で負、大きさは B/(T_soft − T) 程度（応力が熱膨張の差で生じ、軟化点 T_soft で消えるため）[AE]。
  B = 3–5×10⁻⁴、T_soft − T ≈ 10³ K で |dB/dT| ≈ 3–5×10⁻⁷ /K。D11 の 2.4 rad/(m·K) は 6×10⁻⁷ /K に相当。
  ここでは 2–6×10⁻⁷ /K を範囲とする（UNVERIFIED、T09 の値で差し替え）。

| |dB/dT| [1/K] | rad/(m·K) | 6.3 m | 63 m | 630 m | Liu の 153 m / 2.45 km [AE] |
|---|---|---|---|---|---|
| 2×10⁻⁷ | 0.81 | 5.1 rad/K | 51 rad/K | 510 rad/K | 124 / 2000 rad/K |
| 4×10⁻⁷ | 1.63 | 10 | 102 | 1020 | 249 / 4000 |
| 6×10⁻⁷ | 2.44 | 15 | 154 | 1540 | 373 / 6000 |

  温度ドリフト 1 mK で 63 m なら 0.05–0.15 rad。D16 の目標（θ_rms ~ 10 mrad 級）には温度を 0.1 mK 以下に保つか、別の対策が要る。
  周波数帯は温度の揺らぎの帯域（< 1 Hz、T08）。

- **PM・直交偏光・TEC のアクチュエータ**: TEC が雑音と同じ区間の温度を動かすなら、温度の雑音とアクチュエータの軸の比が同じなので、(1) の括弧は温度の項について 0 になる。
  残るのは「温度が一様でない部分」（TEC が届かない区間、勾配）と、歪み・音響の雑音のうち 2 軸で応答の違う部分。
  Liu の TEC がこの働きをしているかは書かれていない（§2.1）。

- **軸の入れ替え（90° 融着）**: 共通区間の中点で slow と fast を入れ替えると、probe も squeezed 光も前半と後半で違う軸を通るので、一様な温度変化では (2) が打ち消し合う。残りは
  δθ_res = (2π/λ0)(dB/dT)·[(L/2)(ΔT_1 − ΔT_2) + δL_sw·ΔT]  (3)
  （ΔT_1、ΔT_2: 前半・後半の温度変化、δL_sw: 入れ替え点の位置のずれ）。2 つの半分を同じスプールに一緒に巻けば ΔT_1 ≈ ΔT_2 にできる。
  代わりに融着の損失（squeezed 光にかかる、T12）と PER の低下（迷光、§3.4）が増える。

- **SMF・直交偏光**: 相反性から、往路の Jones 行列が M なら復路は Mᵀ。squeezed 光が e_in → M e_in ∝ p（GBS 側の PBS の透過偏光）となるように合わせると、
  GBS 側で s 偏光（p に直交）で入れた probe は、上流側でちょうど e_in に直交する偏光に出てくる（M を SU(2)×位相に分けて確認 [AE]）。
  このとき 2 つの位相の差は M の SU(2) 部分の位相差（往復で 2 倍に見える楕円リターダンス）だけになり、PM の (2) よりずっと小さいはず（SMF の複屈折は 10⁻⁷ 程度で、PM の 10⁻³ 倍 [AE]）。
  63 m の SMF のリターダンスを 25 rad（B ≈ 10⁻⁷）、その温度による相対変化を 10⁻³ /K とすると 0.03 rad/K [AE、根拠は弱い]。
  ただし M そのものが揺らぐ（T09 の偏波ドリフト）と、PBS で squeezed 光が失われ、2 つの位相の差も変わる。数値は T09 待ち。

- **同じ偏光**: (1) の括弧は 0（応答が同じ）。分離に PBS を使えないので、§3.4・§4 の別の方法が要る。

### 3.2 逆向き伝搬による時刻のずれ（E42）

**式**: 共通区間の位置 x にある雑音 δ(t)（両方の軸・向きで同じ応答）を、アクチュエータが x_a にあるロックで補正する（帯域内、利得は十分大きい）。§2.3 の式から、GBS 入力の squeezed 光に残る位相は

  δθ_res(t) = δ(u) − δ(u − 2(τ_x − τ_{x_a})),  u = t − τ_L + τ_x  (4)

  したがって S_res(f) = 4 sin²(2π f n_g |x − x_a|/c)·S_δ(f)  (5)

位置 x のときの残差が補正前の雑音と同じになるのは f = c/(12 n_g |x − x_a|)。

雑音がファイバに一様に分布し（x ごとに独立で、PSD は単位長さあたり一様）、アクチュエータが端（x_a = 0 か L）にあるとき、x で平均して

  S_res(f)/S_free(f) = 2[1 − sin(y)/y],  y = 4π f n_g L/c  (6)

低い周波数では ≈ y²/3 = (16π²/3)(f n_g L/c)²。比が 1 になるのは y = 1.895、すなわち **f₁ = 0.151·c/(n_g L)**。
（ファイバ雑音補償の往復方式で知られる「遅延で抑えきれない雑音」と同じ形。Williams et al. 2008 は往復方式で (1/3)(2π f τ)² を与える〔UNVERIFIED、本文未読〕。
ここは往路と復路を別の光で比べるので、片道の遅延 τ_L に対して (1/3)(4π f τ_L)² となり、4 倍大きい。式 (4)–(6) はこの担当が導いたもの。）

アクチュエータを中点に置くと |x − x_a| ≤ L/2 になり、f₁ は 2 倍になる。

| 構成 | L | 往復 2n_g L/c | c/(2n_g L)（survey の目安） | f₁（式 6） | 1 kHz での比 | 10 kHz | 100 kHz |
|---|---|---|---|---|---|---|---|
| f_rep 1 GHz, k=31 | 6.3 m | 62 ns | 16 MHz | 4.9 MHz | −73 dB | −53 dB | −33 dB |
| f_rep 100 MHz, k=31 | 63 m | 0.62 µs | 1.6 MHz | 0.49 MHz | −53 dB | −33 dB | −13 dB |
| f_rep 10 MHz, k=31 | 630 m | 6.2 µs | 0.16 MHz | 49 kHz | −33 dB | −13 dB | 〜+4 dB |
| Liu 第 2 配列の最長 [AE] | 2.45 km | 24 µs | 42 kHz | 12.6 kHz | −21 dB | −1.7 dB | — |

**帯域**: 音響・振動の雑音（Hz–kHz、T08）では、630 m までは 30 dB 以上抑えられる。問題になるのは、ループの帯域を上げて kHz より上の雑音を抑えたい場合。

**ループ遅延**: probe はアクチュエータを通ってから検出器まで、x_a から取り出し点までの共通区間を通る。遅延 n_g x_a/c がループ遅延（T16 の遅延バッファ）に加わる。
アクチュエータを GBS 側に置くと 63 m で 0.31 µs、630 m で 3.1 µs、2.45 km で 12 µs。**取り出し側（上流）に置けば加わらない**。式 (6) はどちらの端でも同じ。
Liu の図 4b では PZT は probe を入れる側（第 2 の干渉計の近く）にあるように見える（読み取り、未確認）。

### 3.3 周波数をずらした probe（E46 の一部）

probe の周波数を Δν_pr ずらすと、同じ区間で squeezed 光との位相差は ψ(ω0 + 2πΔν_pr) − ψ(ω0) ≈ 2π Δν_pr·τ_g（τ_g = n_g L/c）。その温度ドリフトは

  δθ_res = 2π Δν_pr·L·(dτ_g/dT)·ΔT  (7)

dτ_g/dT = 39 ps/(km·K)（survey の [AE]、UNVERIFIED、T08）で、63 m のとき:

| Δν_pr | 例 | rad/K（63 m） |
|---|---|---|
| 80 MHz | AOM | 1.2×10⁻³ |
| 1 GHz | EOM の側帯 | 0.015 |
| 10.8 GHz | Brillouin シフトと同じ程度 | 0.17 |
| 100 GHz | DWDM の隣のチャネル（0.8 nm） | 1.5 |
| 1 THz | 8 nm 離した WDM | 15 |

AOM 程度のずれなら無視できる。WDM で分離できるほど離すと、PM の直交偏光（§3.1）より 1–2 桁小さいが、温度を mK に保つ必要がある程度には大きい。
さらに大きく離すと色分散（β2 の温度係数）の項が加わる（T06/T08）。

### 3.4 迷光: probe の光子が GBS に入る量（E43、E44、E35、D28）

#### 3.4.1 経路

逆向きの probe が GBS に向かう向き（d = +1）に戻るのは、次の場合だけである。

1. 分布した Rayleigh 後方散乱（周波数は probe と同じ。散乱体が動くので位相はゆっくり揺らぐ）
2. 自発 Brillouin 後方散乱（±ν_B ≈ ±10.8 GHz、幅 Δν_B ≈ 20–30 MHz、UNVERIFIED、T07）
3. 自発 Raman 後方散乱（THz にわたって広い、T07）
4. 不連続な反射（コネクタ、コリメータの面、融着点、ストレッチャ、入れる PBS のすぐ後ろのファイバ端）

survey の E43「P_probe × 10^(−PER/10)」は、PER で squeezed 光の偏光に移った probe も**逆向き**に進むので、そのままでは GBS に届かない。
PER は「後方散乱・反射のうち、squeezed 光の偏光に入る割合」として 1.–4. に掛かる（ε_pol）。

#### 3.4.2 式

GBS 入力での迷光のパワー（時刻 t）は、probe の波形と後方散乱のインパルス応答の畳み込み:

  P_st(t) = ∫_0^L dz b_tot ε_pol(z) e^{−2αz} P_pr(t − 2n_g z/c)
        + Σ_i R_i ε_i P_pr(t − 2 n_g z_i/c)
        + ε_pol ε_AOM P_pk b_tot L_bs（AOM の漏れによる CW 成分。パルスのときも常にある）  (8)

  b_tot = b_R + b_B F_GBS(ν_B) + ∫ρ_Ram(Δν) F_GBS(Δν) dΔν

（z は probe を入れる点から測った距離。最初の項は OTDR の式と同じで、t' = 2n_g z/c に変数を変えると (v_g/2) b_tot ∫_0^{2τ_L} P_pr(t − t') dt'。αL ≪ 1 で e^{−2αz} ≈ 1。
AOM の漏れの項は、漏れた probe も同じように散乱されることを表す（P̄_pr を P̄_pr + ε_AOM P_pk に置き換えるのと同じ）。
周波数をずらした場合は F_GBS(Δν_pr)·I_WDM(Δν_pr) を Rayleigh と反射の項に、F_GBS(Δν_pr ± ν_B) を Brillouin の項に掛ける。）

窓あたりの光子数: N_stray = (1/hν0) ∫_{窓} P_st(t) dt。CW、または窓がパルスの往復の広がりより十分短く平均してよいとき:

  N_stray ≈ (P̄_pr T_gate/hν0)·[ε_pol (b_R + b_B + ρ_Ram B_f) L_bs + Σ_i R_i ε_i]  (9)

squeezed の時間モードの中にいる光子数（BHD の測定・GBS の干渉に効く量）は n_mode ≈ Φ_st·τ_mode（Φ_st = P_st/hν0、スペクトルの重なりが 1 のとき）。

#### 3.4.3 係数

- **b_R**: Corning のデータシートの Rayleigh 後方散乱係数 −82 dB（1550 nm、1 ns パルス、SMF-28e+・Ultra 200 など、typical 値）から、
  b_R·(v_g·1 ns/2) = 10^{−8.2} → **b_R = 6.2×10⁻⁸ /m**（n_g = 1.468 とした）。
  これは S ≈ 1.5–1.7×10⁻³、α_s ≈ 0.15–0.16 dB/km とも合う（survey の値と一致）。
  PM ファイバは Ge の多いコアで損失が大きい（0.5–1 dB/km、UNVERIFIED）ので、b_R は 1–3 倍とみる（T09）。
- **b_B**: 自発 Brillouin を、利得から b_B ≈ (g_B/A_eff)·hν0·(2n_th + 1)·(π Δν_B/2) と見積もる（Stokes と anti-Stokes の和。n_th = k_BT/(hν_B) ≈ 569 at 295 K）。
  g_B = 5×10⁻¹¹ m/W、A_eff = 80 µm²、Δν_B = 30 MHz（いずれも UNVERIFIED）で b_B ≈ 4.3×10⁻⁹ /m、**b_R/b_B ≈ 14** [AE]。
  この比は Landau–Placzek 比（Rayleigh / Brillouin）にあたり、Wait & Newson 1996 などが分布型センサーで使っている量。値は原典で未確認（T07）。
  b_B は温度に比例（n_th ∝ T）。
- **ρ_Ram**: 1550 nm 帯で数 nm 離れたところの後方 Raman は ρ_Ram ~ 5×10⁻⁹ /(km·nm)（量子鍵配送と古典信号の同居の文献の値、UNVERIFIED、T07）。
- **R_i**: FC/APC で 10⁻⁶ 以下（反射減衰量 60 dB 以上、UNVERIFIED）。1 個で 63 m の Rayleigh 全体の 1/4 に相当（b_R L = 3.9×10⁻⁶）。
- **ε_pol**:
  - 同じ偏光: 1。
  - PM・直交偏光: Rayleigh と Brillouin（等方的な密度揺らぎと縦の音響波）は偏光を保つので、散乱光は probe の軸に乗る。
    squeezed 光の偏光に入るのは ε_xt（ファイバの h·L、コネクタ・融着の PER、散乱の偏光解消成分）と、GBS 側の PBS の透過の漏れ T_⊥ の和:
    ε_pol ≈ ε_xt + T_⊥。部品の PER 20–30 dB、PBS の消光比 30 dB（UNVERIFIED、T09/T12）で ε_pol ≈ 2×10⁻³–10⁻²。
    なお、GBS の検出器は偏光をほとんど区別しないので、**数（クリック）としての汚れ**は「squeezed 光の偏光に入る量」ではなく「PBS を抜けて GBS 側に出る量」で決まる。
    PBS の構造上、probe の偏光（s）の後方散乱は probe の入口へ戻り、GBS 側には T_⊥ しか出ない。
  - SMF・直交偏光: 位置 z からの後方散乱の偏光は M_zᵀM_z s（M_z: 入口から z までの Jones 行列）で、z ごとに違う。
    M_z を一様にランダム（Haar）とすると |⟨p|M_zᵀM_z|s⟩|² の平均は **1/3**（数値で確認 [AE]）。長いファイバではこの平均に近づき、短いファイバでは時間とともに 0–1 の間で揺れる。
    したがって ε_pol ≈ 1/3 + (2/3)T_⊥ ≈ 0.33（−4.8 dB）。

#### 3.4.4 時間のずらし（ゲート）がどこまで効くか

1 つの probe パルス（エネルギー E_pr）の Rayleigh 散乱光は、GBS 側の入口に [t_inj, t_inj + 2n_g L/c + τ_pr] の間に戻り、その総エネルギーは E_pr·b_R·L_bs。
したがって**分布した散乱は、パルスでも平均パワーが同じなら窓あたりの平均光子数は CW と同じ**。時間をずらして避けられるのは:

- 不連続な反射（遅延 2n_g z_i/c が決まっている）と、入れる点のすぐそばの反射。Liu の「時間をずらす」の主な効果はこれと考えられる（原典に理由の説明はない）。
- 往復の広がりが GBS のフレームの空き時間に収まるポートの分布散乱。D08 では 32 ポートが同じスロットに着くので、GBS のスロットは T_f = 32 T_rep ごと。
  probe を T_f ごとのパルスにして、条件 2n_g L_k/c + τ_pr + T_gate < T_f、すなわち L_k = k·c·T_rep/n_g なら **2k·T_rep + τ_pr + T_gate < 32 T_rep（k ≤ 15）** を満たすポートは、散乱光をすべてスロットの外に置ける。
  k ≥ 16 のポートでは、往復の広がりがフレームより長いので、平均として CW と同じになる。
  ただし probe 分配の光路（pdist_k）の散乱と、AOM の漏れ ε_AOM P_pk は残る。ε_AOM P_pk ≪ P̄_pr、すなわち ε_AOM ≪ D であることが条件。
  この条件は「GBS が毎フレーム動いている」ことを仮定している（未確認）。

#### 3.4.5 概算（P̄_pr = 1 mW、L = 63 m、b_R = 6.2×10⁻⁸ /m）

Φ_pr = 7.8×10¹⁵ /s。

| 条件 | ε・フィルタ | N_stray（T_gate = 1 ns） | N_stray（10 ns） | n_mode（τ_mode = 1 ps） |
|---|---|---|---|---|
| 同じ偏光・同じ周波数、Rayleigh | 1 | 30 | **3.0×10²** | 0.030 |
| 同上、Brillouin（GBS のフィルタが ±11 GHz を通すとき） | 1 | 2.1 | 21 | （スペクトルの重なりによる） |
| 同上、Raman（B_f = 0.4 nm / 35 nm） | 1 | 0.001 / 0.09 | 0.01 / 0.9 | — |
| PM・直交偏光（ε_pol = 2×10⁻³） | | 0.06 | 0.6 | 6×10⁻⁵ |
| SMF・直交偏光（ε_pol ≈ 1/3、揺らぐ） | | 10 | 1.0×10² | 0.010 |
| 波長をずらして WDM（I_WDM = 10⁻³–10⁻⁵）、Rayleigh のみ | | 0.03–3×10⁻⁴ | 0.3–3×10⁻³ | — |
| FC/APC 1 個（R = 10⁻⁶、同じ偏光、時間をずらさないとき） | | 7.8 | 78 | 0.0078 |

すべて P̄_pr·L·T_gate に比例する（反射の項は L によらない）。温度依存は Brillouin と Raman だけ（n_th ∝ T）。
光子数の許容量は Open（D28）なので、この表と式 (9) を関数として渡す。

**squeezing への影響**: 迷光が squeezed の時間モードに位相がランダムな光子 n_mode を加えると、分散は n_mode だけ増える（真空の分散 1/2 の規約）。
10 dB（V = 0.05）に n_mode = 0.03 を足すと V = 0.08 で **−8.0 dB**。τ_mode = 10 ps なら n_mode = 0.3 で −1.5 dB。
同じ偏光・1 mW では指標 (i) の squeezing 自体も壊れる。PM 直交偏光、または §3.7 の小さい probe なら無視できる。

### 3.5 XPM（E30 の probe 部分）

probe が squeezed 光に与える XPM の位相は、CW または平均として φ_XPM = b_X·γ·P̄_pr·L（b_X = 2 同じ偏光、2/3 直交の直線偏光、Agrawal の式、UNVERIFIED）。
逆向きでもパルス列の平均は同じで、その揺らぎは往復時間で平均される（走り抜け）。probe 自身はこの位相を受けないので、ロックでは補正されない。
γ = 1.3 /(W·km)（UNVERIFIED、T06）、1 mW、63 m で 1.6×10⁻⁴ rad（同じ偏光）、5.5×10⁻⁵ rad（直交）。probe の RIN による揺らぎはさらに 10⁻³ 倍以下。
**無視できる**（P̄_pr が 100 mW を超えない限り）。ブロックとしては G の XPM 項でなく、T の静的なオフセット＋ RIN 比例の雑音として持てば足りる。

### 3.6 入れ出しの素子の損失（E45）

squeezed 光にかかる損失 IL（dB）と、初めの squeezing S0 に対する出力（V_out = η V_in + (1 − η)/2、conventions §2）:

| S0 \ IL | 0.05 dB | 0.1 dB | 0.3 dB | 0.6 dB | 1.2 dB | 1.8 dB |
|---|---|---|---|---|---|---|
| 3 dB | −2.95 | −2.90 | −2.72 | −2.48 | −2.06 | −1.74 |
| 6 dB | −5.85 | −5.71 | −5.21 | −4.59 | −3.65 | −2.96 |
| 10 dB | −9.57 | −9.19 | −7.96 | −6.65 | −4.99 | −3.92 |

1 ポートに入る素子: 取り出し PBS（または WDM、サーキュレータ）と入れる PBS の 2 個、アクチュエータ、（軸の入れ替えなら）融着 1 点。
目安（すべて UNVERIFIED、T12 で一次資料に置き換える）: 自由空間の PBS 1 個 0.5–1%（0.02–0.04 dB）、ファイバ型 PBS 0.3–0.8 dB、
PM サーキュレータ 1 回の通過で 0.6–0.9 dB（第三者のカタログでは Thorlabs CIR1550PM-FC が < 0.9 dB、メーカーの資料で未確認）、
ファイバ型 WDM 0.3–0.6 dB、PM の融着 0.05–0.3 dB。
ファイバ型の素子 2 個（1–2 dB）では 10 dB が 4–5 dB になる。**自由空間の PBS か、低損失の WDM フィルタ（自由空間の二色鏡）にしないと、どの方式でも損失が支配的になる**。

### 3.7 probe のパワーとロックの雑音の関係

probe（検出パワー P_det = T_ret η_det P̄_pr）を強い参照光とホモダインで混ぜると、ショット雑音で決まる位相の PSD は（片側）

  S_φ,sn = κ_het·hν0/(2 η_det T_ret P̄_pr),  κ_het = 1（ホモダイン、Liu の DC 干渉）, 2（ヘテロダイン、Ref A の 1 MHz）  (10)

（導出: N 光子のコヒーレント光の位相の分散は 1/(4N)、時間 T の平均は帯域 1/(2T) に相当。係数は T14 で確認。）
このセンサ雑音はロックで squeezed 光の位相に書き込まれる: σ_lock² = S_φ,sn·B_loop。
(9) と (10) を掛けると P̄_pr が消える:

  N_stray·σ_lock² = κ_het·ε_pol·(b_R + b_B + …)·L·B_loop·T_gate/(2 η_det T_ret)  (11)

63 m、B_loop = 10 kHz、T_gate = 1 ns、η_det = 0.8、T_ret = 0.5、同じ偏光で **4.9×10⁻¹¹ rad²**。
σ_lock = 1 mrad なら P_det ≈ 0.8 nW（P̄_pr ≈ 1.6 nW）、N_stray ≈ 5×10⁻⁵ /ns。1 µW でも 1 ns あたり 0.03。
したがって**迷光の問題は「1 mW の probe」という仮定から来ており、ショット雑音の点では nW–µW で足りる**。
実際の下限は: 検出器の電気雑音・ADC（参照光を強くすればショット雑音で決まる範囲に入れられるか、T13）、DC 干渉の場合はオフセットのドリフト（D06 と同じ問題）、
参照光の強度雑音、probe の偏光が変わったときの可視度の低下。これらは T13/T14/T15 で決まり、その結果で P̄_pr の下限が決まる。

### 3.8 レーザー位相雑音とロック（D31 の 3）

§2.3 の式から、帯域内でロックの利得が十分大きいとき

  Θ_k ⊃ [φ_L(t − T_a − τ_k) − φ_L(t − T_a − τ_0)] − [φ_L(t − T_b − τ_k) − φ_L(t − T_b − τ_0)]  (12)

  S_Θ,L(f) = 16 sin²(π f (T_a − T_b))·sin²(π f (τ_k − τ_0))·S_φL(f)   （帯域内）  (13)
  S_Θ,L(f) = 4 sin²(π f (τ_k − τ_0))·S_φL(f)                        （帯域外、D31 の 3 そのもの）  (14)

probe もポート k の光路 L_k を（逆向きに）通るので、probe にも τ_k が入るためである。
T_a = T_b（「レーザー → OPA → ポートの入口」と「レーザー → probe 分配 → probe 検出器（共通区間を除く）」の遅延をそろえる）なら、帯域内のレーザー位相雑音は消える。
そろえないと、帯域内では遅延差 |T_a − T_b| の干渉計の雑音が加わる。

白色の周波数雑音（S_φL = Δν/(π f²)）では、(14) の分散のうち B_loop より下の割合は約 2(τ_k − τ_0)B_loop（τ_k − τ_0 = 310 ns、B_loop = 10 kHz で 0.6%）なので、この効果は小さい。
ただし実際のレーザーは低い周波数で 1/f の周波数雑音とドリフトが大きいので、(13) と (14) の差は低い周波数で効く。laser の S_ν(f) は T05 から受け取る。
**D31 の「3 はロックでは補正できない」は帯域外では正しいが、帯域内では T_a − T_b で決まる**（DECISIONS への提案、§8）。

## 4. 対策の比較

| 対策 | 軸の差（§3.1） | 周波数の差（§3.3） | 迷光（§3.4、1 mW・63 m・10 ns） | squeezed 光の損失 | 備考 |
|---|---|---|---|---|---|
| Liu 型そのまま（PM、直交、同じ周波数、PBS） | 51–154 rad/K（63 m） | 0 | 〜0.6 | 自由空間 PBS 2 個で < 0.1 dB | 温度の対策が必須 |
| ＋ 90° 融着で軸を入れ替える | 勾配 ΔT_1 − ΔT_2 のみ（式 3） | 0 | PER 低下で増える | ＋融着 1 点 | 同じスプールに巻く |
| ＋ TEC を雑音と同じ区間に（温度のアクチュエータ） | 温度の項は 0（式 1）、勾配が残る | 0 | 同じ | なし | 遅い（< 1 Hz）。速い雑音は PZT |
| SMF、直交、PBS | 0.03 rad/K [AE、弱い] | 0 | 〜1×10²（ε ≈ 1/3、揺らぐ） | 偏波のずれで PBS 損失が揺らぐ（T09） | 偏波の制御が要る |
| 同じ偏光、同じ周波数、サーキュレータ | 0 | 0 | 3×10²（減らない） | 0.6–0.9 dB × 2 回（UNVERIFIED） | 損失で不可 |
| 同じ偏光、AOM で周波数をずらす | 0 | 10⁻³ rad/K | 3×10²（GBS のフィルタで分けられない） | 分離手段がない | 単独では不可 |
| 同じ偏光、波長をずらして WDM | 0 | 1.5 rad/K（100 GHz）〜 15 rad/K（1 THz） | Rayleigh は I_WDM 倍。Raman・Brillouin が残る | WDM 2 個 | パルスの帯域より外へずらす必要。100 fs（約 35 nm）では難しい |
| パルス化＋時間のずらし | — | — | 反射は消える。分布散乱は k ≤ 15 だけ消える（§3.4.4） | なし | AOM の消光 ε_AOM ≪ D が条件 |
| probe のパワーを下げる（§3.7） | — | — | ∝ P̄_pr。nW で 10⁻⁴ 程度 | なし | 技術雑音で下限（T13/T14） |
| 振り分け素子で probe を取り出す（逆向きの probe を、空いている時刻に別の出口へ送る） | 素子による | — | — | 取り出し PBS が不要 | 案。素子の仕様（T11）次第で、未検討 |

組み合わせとしては「PM・直交・自由空間 PBS・nW–µW の probe・パルス化・TEC と PZT の 2 段（TEC は同じ区間の温度）・軸の入れ替え」が、どの項目でも桁で破綻しない。
ただし軸の差の残り（温度の勾配）が θ_rms の要求に収まるかは、T08 の温度の空間分布と T04 の要求で決まる。

## 5. Ref A（phase-detection OPA）と Liu 型の統合（D30）

**Ref A の lock 1** [Ref A §2、§3.1]: ポンプと probe（seed、AOM で 1 MHz シフト）を squeezer OPA の手前で分け、別の PPLN（phase-detection OPA）に入れ、パラメトリック利得からポンプ–probe の相対位相を読む。
「分けた後、2 つの OPA までの光路でポンプと probe の相対位相を乱すものがなければロックできる」。lock 2 は BHD 出力の 1 MHz のビートで LO–probe の位相。
Ref A では probe は squeezer OPA にも入り、BHD まで squeezed 光と一緒に進む（従来方式で µW、新方式で mW まで上げられる）。

**それぞれが何を見るか**:

| ロック | 見る量 | 補正する区間 |
|---|---|---|
| Ref A lock 1（phase-detection OPA） | φ_p − 2φ_seed（PD-OPA の位置で） | `src`（ポンプ光路、SHG、ポンプと seed の分岐後の非共通部分）。θ_sq の「OPA を出る時点の値」を seed に対して固定する |
| Liu 型（ポートごと） | ψ_com,k（＋ pdist_k、pext_k、rdist_k） | `com_k` |
| Ref A lock 2（BHD の probe–LO） | φ_LO,k − φ_seed（BHD で） | 指標 (i) の LO 側。GBS に入る squeezed 光は直さない（D10） |

**統合で必要なこと**:

1. Ref A の probe（seed）を squeezer OPA に入れると、squeezed 光と同じ向き・同じ偏光で GBS まで届き、窓あたり Φ·T_gate の光子になる（1 µW で 7.8×10¹² /s × 10 ns ≈ 7.8×10⁴ 光子）。
   **squeezer OPA には seed を入れない**（phase-detection OPA だけに入れる）か、seed をポートのスロットから外れた時刻だけに入れる必要がある。Ref A の lock 2 は GBS 構成では使えず、BHD（指標 (i)）は別の測定用に限る。
2. 単一の OPA の出力を時間で分ける（D08）ので、`src` の位相雑音はポートごとに違う時刻でサンプルされる: 帯域外で 4 sin²(π f (τ_k − τ_0))·S_ψsrc/4 = sin²(π f (τ_k − τ_0))·S_ψsrc（θ_sq = φ_p/2 + π/2 なので ψ_src は 1/2 倍で入る）。
   Ref A の lock 1 はこれを帯域内で消す。Liu の構成（Extended Data Fig. 4a）でも、ポンプの光路を別のロックで固定している。
3. Liu 型の参照光と Ref A の seed は同じレーザーから分け、その分岐後の光路（参照光の配送 vs PD-OPA への seed の光路）は非共通。

**統合しても残る非共通の光路**（D31 の分類つき）:

- (a) ポンプ・seed の分岐点から squeezer OPA と PD-OPA までの非共通部分（Ref A の前提。OPA を並べて置けば短い）— (2)
- (b) `up`: squeezer OPA → 振り分け素子 → 取り出し PBS。振り分け素子の位相雑音を含む。どのロックにも見えず、ポートごとに時刻のずれたサンプル — (1) にも (2) にも入らない「squeezed 光だけの区間」。D31 の分類に追加が要る（§8）
- (c) `down_k`: 入れる PBS → GBS 入力 — 同上
- (d) `pdist_k`、`pext_k`、`rdist_k`（probe の分配、取り出し、参照光の配送）— (2)
- (e) seed（PD-OPA へ）と参照光（Liu の検出器へ）の分岐後の光路 — (2)
- 加えて、レーザー位相雑音の帯域内の項は T_a − T_b で決まる（§3.8）。T_b を T_a にそろえる遅延の設計が統合の要点になる。

## 6. シミュレータに必要なもの

### 6.1 ブロックと表

| 段 | ブロック / 表 | 入力 | 出力 | 備考 |
|---|---|---|---|---|
| T | 区間（segment）: 長さ、n_g、向き、軸、雑音源への応答 κ_{j,a} | 雑音源の時系列（T08）、温度（共通） | 軸ごとの ψ_seg,a(t) | D11。**軸ごとの κ を必須にする**（(1) を計算で出すため） |
| T | 逆向き伝搬の遅延バッファ（区間を N_sub 個の小区間に分け、x ごとの遅延 τ_x で読む） | ψ の小区間ごとの値 | probe と squeezed 光それぞれの合計 | (5)・(6) を再現できる細かさ。N_sub は f·τ_L ≪ 1 となる最高周波数で決める |
| T | probe センサ: 位相の線形結合 φ_pr − φ_ref ＋ ショット雑音 (10) ＋ 技術雑音 | ビームの位相、P̄_pr、κ_het | 誤差信号 | T14 と共通化 |
| T | アクチュエータの位置 x_a と軸の比 ρ_act（PZT、TEC） | 制御出力 | 軸ごとの位相 | TEC は温度として区間に入る |
| T | ロックの遅延: probe の伝搬 n_g x_a/c をループ遅延に加える | x_a | 遅延バッファの長さ | T16 の遅延バッファを使う |
| O | 迷光の表 N_stray(P̄_pr, CW/パルス, τ_pr, t_offset, Δν_pr, 偏光, L, T, T_gate) | b_R、b_B（T07）、ρ_Ram（T07）、R_i、ε_xt・T_⊥（T09）、ε_AOM、F_GBS、I_WDM | N_stray、n_mode、時間波形 P_st(t) | 式 (8) の畳み込み。D28 の「関数として出す」 |
| O | 入れ出しの素子の損失 | 素子の種類と IL | 静的な η（squeezed 光に掛ける） | T12 と共通 |
| T | XPM の静的オフセット＋RIN 雑音 | P̄_pr、RIN、γ、L | squeezed 光の軸の ψ | 既定でオフ |
| M | 指標 (ii) の分解: (1) ロックされる区間、(2) 共通点、(3) レーザー、＋「squeezed 光だけの区間」 | 全区間の ψ | 4 つの成分の時系列と PSD | D31 の 3 分割に 1 つ加える提案（§8） |
| M | 迷光による squeezing の劣化: V → V + n_mode | n_mode | 指標 (i) の dB | |

### 6.2 切り替えのフラグ（すべて設定で切り替え、保護リストの「雑音源のトグル」）

- `probe.polarization`: `same` | `orthogonal`
- `probe.mode`: `cw` | `pulsed`（`tau_pr`、`period`、`t_offset`、`aom_extinction`）
- `probe.freq_offset_hz`: Δν_pr（0、AOM、WDM の値）
- `probe.power_avg_w`
- `probe.separation`: `pbs` | `circulator` | `wdm` | `coupler`
- `fiber.type`: `pm` | `smf`（D32）、`fiber.axis_swap`: 位置のリスト
- `actuator.position_m`、`actuator.kind`: `pzt_strain` | `tec`（軸の比を持つ）
- 雑音の項ごとのオン・オフ: `birefringence_drift`、`counterprop_delay`、`freq_offset_drift`、`xpm_probe`、`stray_rayleigh`、`stray_brillouin`、`stray_raman`、`stray_reflections`、`insertion_loss_probe_optics`
- `lock.refA_pdopa`: Ref A の lock 1 を加える（§5）、`lock.delay_match_Ta_Tb`: T_b の設定

### 6.3 他のトピックから受け取る値

T07: b_B、ν_B、Δν_B、ρ_Ram(Δν, T)、偏光の保存の度合い。T08: 雑音源の PSD、温度の空間分布（勾配）、dτ_g/dT。T09: dB/dT、h パラメータ、部品の PER、SMF の M(t) のモデル。
T11: 振り分け素子の位相雑音（`up` 区間）と、逆向きの probe を通せるか。T12: 素子の損失。T13/T14: probe 検出の技術雑音、κ_het。T16: アクチュエータの軸の比、ループ遅延。T05: S_ν(f)。T04: 必要な θ_rms。

## 7. 検証テスト（実装が通すべきもの）

1. **同じ区間・同じ向き（prompt §7-11）**: probe と squeezed 光が同じ向き・同じ軸ですべての区間を共有するとき、共有区間の雑音はロックの残差にも Θ_k にも現れない（ループ利得 → ∞、相対誤差 < 10⁻⁶）。
2. **逆向き・点の雑音**: x にだけ白色の位相雑音、アクチュエータを x_a に置いた理想ロック（帯域無限、遅延 0）で、squeezed 光の残差 PSD が 4 sin²(2π f n_g|x − x_a|/c)·S_δ に 2% 以内で一致（f < 0.3 c/(n_g|x − x_a|)）。x = x_a で残差 0。
3. **逆向き・分布した雑音**: L = 63 m、n_g = 1.468、アクチュエータを端に置き、小区間ごとに独立な雑音を一様に与えると、残差/補正前の比が 2[1 − sinc(4π f n_g L/c)] に一致し、比が 1 になる周波数が 0.489 MHz ± 5%。1 kHz で −53 dB ± 1 dB。
4. **軸の差**: PM、直交偏光、温度を一様に 1 K 上げる（ゆっくり）、PZT で理想ロック → squeezed 光の位相の変化 = (2π/λ0)(dB/dT)L。dB/dT = 4×10⁻⁷ /K、L = 63 m、λ0 = 1545.32 nm で 102.5 rad ± 0.1%。
   同じ設定で TEC（温度）をアクチュエータにすると 0。中点で 90° 融着すると 0、前半だけ 1 K 上げると 51.2 rad。
5. **周波数の差**: Δν_pr = 1 GHz、dτ_g/dT = 39 ps/(km·K)、L = 63 m、1 K で 0.0154 rad ± 0.1%。
6. **迷光の総量**: CW 1 mW、b_R = 6.18×10⁻⁸ /m、L = 63 m、α = 0、同じ偏光、λ0 = 1545.32 nm で 10 ns あたり 302.8 光子 ± 0.5%。
   パルス（任意の波形）で、1 パルスの散乱光の時間積分 = E_pr·b_R·L_bs、広がりの幅 = 2n_g L/c（0.617 µs）。
7. **ゲート**: f_rep = 100 MHz、T_f = 320 ns、k = 15（L = 30.6 m）、τ_pr = 10 ns、T_gate = 1 ns、ε_AOM = 0 で、適切な t_offset のとき窓の中の Rayleigh 光子が 0。k = 16 ではどの t_offset でも 0 にならない。
8. **SMF の偏光**: Haar でランダムな M について |⟨p|MᵀM|s⟩|² の平均が 1/3 ± 1%（10⁵ 標本）。M が対角（PM）なら 0。
9. **相反性**: 損失のない任意の M で M e_in ∝ p のとき、Mᵀ s ∝ e_in^⊥ で、squeezed 光と probe の位相の和（共通位相の 2 倍）と差（リターダンス）が式どおりに分かれる。
10. **レーザー位相雑音**: 白色 FM のレーザー、理想ロック（帯域内）で T_a = T_b なら Θ_k のレーザーの項が 0。T_a − T_b = 10 ns、τ_k − τ_0 = 310 ns で PSD が式 (13) に 2% 以内で一致。ロックを切ると式 (14)。
11. **probe のショット雑音**: P_det = 1 nW、η = 1、ホモダインで S_φ = hν0/(2P_det) = 6.43×10⁻¹¹ rad²/Hz ± 2%、ロックを通して squeezed 光の位相に 1:1 で移る（帯域内）。
12. **損失**: 10 dB の squeezing に IL 0.3 dB で −7.96 dB（conventions §2 の純損失の式、prompt §7-3 と同じ）。
13. **迷光の squeezing への影響**: V = 0.05 に n_mode = 0.03 を加えて −7.96 dB。

## 8. 決定が必要な点・DECISIONS への提案（調整役へ）

提案（このノートでは DECISIONS・conventions・survey を編集していない）:

- **D31 への追加**: (1)–(3) のほかに「(4) squeezed 光だけが通る区間（`up`: OPA → 振り分け素子 → 取り出し点、`down_k`: 入れる点 → GBS）の雑音」を別に出す。振り分け素子の位相雑音（T11）はここに入る。
  また (3) は「帯域外では 4 sin²(π f (τ_k − τ_0)) S_φ、帯域内では 16 sin²(π f (T_a − T_b)) sin²(π f (τ_k − τ_0)) S_φ」と書き直すのがよい（§3.8）。
- **survey E42 の目安の訂正**: c/(2n_g L) ではなく、分布した雑音でアクチュエータが端なら 0.151·c/(n_g L)（63 m で 0.49 MHz）。
- **survey E43 の訂正**: PER で squeezed 光の偏光に移った probe はそのままでは GBS に届かない。PER は後方散乱・反射に掛かる係数として扱う。
- **survey E44 の確認**: 10 ns あたり約 3.0×10²（1 mW、63 m、同じ偏光）。b_R はデータシートの値に基づく。

ユーザの決定が要る点（推測では埋めない）:

1. GBS が許す迷光の数（D28、Open のまま）: 1 窓あたりの平均光子数で、クリックとして数える量（PBS 後、両偏光）か、squeezed モードの中の量か。
2. probe のパワーの設計値: §3.7 のとおり nW–µW で足りる見込み。1 mW を基準にし続けるか。
3. GBS は毎フレーム（T_f = 32 T_rep）動くか（§3.4.4 のゲートの条件）。
4. アクチュエータの位置（上流の取り出し側を推奨: ループ遅延が増えない）と、TEC を「同じ区間の温度」のアクチュエータとしてロックに含めるか。
5. 振り分け素子を逆向きの probe が通れるか（取り出し点を振り分け素子の上流に置くか下流に置くか）。上流なら素子の位相雑音もロックされるが、probe の時刻を素子の窓に合わせる必要がある（T11）。
6. Liu の補足資料（ファイバの種類、probe のパワー、TEC の制御、精度の測り方）を `refs/` に置けるか。PM 直交偏光で 2.45 km のループを λ/200 に保てた理由が分かれば、§3.1 の評価を検証できる。

## 9. 文献

一次資料（本文を読んだもの）:

- H.-L. Liu, H. Su, Y.-H. Deng, S.-Q. Gong, Y.-C. Gu, et al., "Gaussian boson sampling with 1,024 squeezed states in 8,176 modes," Nature 653 (21 May 2026). DOI: 10.1038/s41586-026-10523-6（PDF に記載）。
  Methods「Phase stabilization of the system」、Extended Data Fig. 1、3、4。補足資料は未読。
- K. Hirota, T. Kashiwazaki, et al., "Generation of 10-dB squeezed light from a broadband waveguide optical parametric amplifier with improved phase locking method"（Ref A）。§2、§3.1。掲載誌・年・DOI は UNVERIFIED（survey と同じ）。

書誌情報だけ確認したもの（本文は未読。式・数値は本文で確認していない）:

- E. Brinkmeyer, "Analysis of the backscattering method for single-mode optical fibers," J. Opt. Soc. Am. 70, 1010–1012 (1980). DOI: 10.1364/JOSA.70.001010。捕捉率 S の式（S ≈ (3/2)(λ0/(2π n w0))²）はこの担当の記憶で、UNVERIFIED。数値には使わず、b_R はデータシートから取った。
- A. Ourmazd, M. P. Varnham, R. D. Birch, D. N. Payne, "Thermal properties of highly birefringent optical fibers and preforms," Appl. Opt. 22, 2374–2379 (1983). DOI: 10.1364/AO.22.002374。dB/dT の数値は未抽出（T09）。
- M. Wait, T. P. Newson, "Landau Placzek ratio applied to distributed fibre sensing," Opt. Commun. 122, 141–146 (1996)。DOI は未確認。Rayleigh/Brillouin 比の値は未抽出。
- P. A. Williams, W. C. Swann, N. R. Newbury, "High-stability transfer of an optical frequency over long fiber-optic links," J. Opt. Soc. Am. B 25(8) (2008)。往復の補償で抑えきれない遅延の雑音。式 (6) との比較に使ったが、本文の式は未確認（UNVERIFIED）。

メーカー資料（D24）:

- Corning, Ultra 200 / SMF-28e+ / SMF-28 ULL Product Information: Rayleigh backscatter coefficient −82 dB（1550 nm、1 ns パルス、typical）。
  https://www.corning.com/content/dam/corning/media/worldwide/coc/documents/Fiber/PI1450_3-2015.pdf （2026-10-28 閲覧を試みたが 403。値は検索結果の抜粋で確認したのみ。**部分的に確認**）。
- Thorlabs CIR1550PM-FC: IL < 0.9 dB は第三者のカタログ（https://www.gophotonics.com/products/fiber-optic-circulators/thorlabs-inc/15-46-cir1550pm-fc 、2026-10-28）の値で、メーカーの資料では未確認（UNVERIFIED）。

記憶から書いたもの（すべて UNVERIFIED、T06/T07 で確認）:

- G. P. Agrawal, Nonlinear Fiber Optics（XPM の係数 2 と 2/3、b_X）。
- R. G. Smith, Appl. Opt. 11, 2489 (1972)（SBS の閾値、自発 Brillouin の雑音の扱い）。
- K. A. Peters et al., New J. Phys. 11, 045012 (2009)、P. Eraerds et al., New J. Phys. 12, 063027 (2010)（量子鍵配送と古典光の同居での Raman 雑音の係数）。
- A. H. Hartog, M. P. Gold, "On the theory of backscattering in single-mode optical fibers," J. Lightwave Technol. 2, 76 (1984)。
