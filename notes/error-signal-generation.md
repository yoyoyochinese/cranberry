# 誤差信号の生成（T14 `error-signal-generation`）

> **訂正（2026-10-31、DECISIONS D65-1）:** dither ロックのセンサ雑音は周期定常（出力 PSD = 2⟨2 sin²θ S_v⟩）。Ref C の残差は 3.58 mrad（2 乗則、連続ループ）・4.11 mrad（ログ検出器）、RAM のずれは 2.54 mrad（暗雑音には RAM をかけない）。§6.3・§6.4 の 3.24・4.2・2.77 mrad は置き換え。Ref B は x = 0.745 で計算し直す（D51-10、lock 2 のずれ 14.64 mrad）。


- Phase: 1（deep dive）
- 日付: 2026-10-28
- 担当する効果: survey E66–E74（prompt §5.5 の誤差信号 1–4）、D06 の #5（DC fringe lock、汎用ブロックの入出力だけ。Ref D の詳細は T15）、D30 の #6（Liu 型の逆向き probe）と「Liu 型 ＋ Ref A の phase-detection OPA」の統合案
- 依存: T13（検出器の応答・回路雑音）。T13 の値が出るまでは下の「仮置き」の値を使う。T02（OPO の利得）、T04（位相雑音 → squeezing）、T10（probe の混入・相殺の崩れ）、T16（復調器・ADC のプリセット）とは数値の受け渡しだけ。
- 規約: `notes/conventions.md`（承認済み、D33）。文献の規約との変換は各節に書いた。
- 値の印: [Ref A] などは本文から読んだ値、[manual] はメーカー資料、[AE] はこのノートでの導出・概算（agent estimate、Phase 5 で要確認）、[MC] はこのノートで行った数値実験、**UNVERIFIED** は一次資料で確認できなかったもの、**仮置き** は出典のない仮の値（依存トピックの値で置き換える）。

---

## 0. 結論の要約

1. **Ref A の phase-detection OPA の誤差信号は 2 f_sh（= 2 MHz）に出る** [AE]。1 MHz シフトの probe は OPA で ω0 − Ω_sh の idler を作り、PD の強度には probe と idler のビート `P_in [cosh 2r + sinh 2r cos(2Ω_sh t − Φ1)]`、Φ1 = 2φ_pr − φ_p が現れる。Ref A 本文は復調周波数を書いていない。同じ原理の Vahlbruch et al. 2006 は「2 倍のビート周波数で復調」と明記している。BHD での lock 2 は f_sh（1 MHz）。
2. **probe を squeezer の後で使う lock 2（Ref A、Ref B）は、lock 1 の誤差と lock 2 の復調位相の誤差を A−/A+ に縮めて Δ に伝える** [AE]。ロック点のずれは `tan ε = (A−/A+) tan δu′`（A± は 2 つの直交位相の振幅利得）。Ref A は A−/A+ = e^{−2r} = 1/123、Ref B は (1−x)/(1+x) = 0.131。lock 1 に要るのは「|δu′| ≪ π/2 を保ち、傾きの符号を守る」ことだけ。lock 1 のセンサ雑音は 9 mrad の主因になりえない。Ref A は SNR の改善が 14 → 9 mrad に効いたと書いているので、ここは食い違う（§3.4、Phase 5 で確認）。
3. **量子雑音ロック（Ref C）のセンサ雑音は radiometer 型**。ゼロスパンのパワー推定の低周波 PSD は `S_P̂ = 2∫S_x(f)² df`（平坦な帯域なら相対分散 1/(B T)）。ここから `S_Δ = ∫S_x² df / (4 ΔP² J1(2m)²)`、平坦な場合は `S_Δ ≈ R_−² / (B J1(2m)² (R_+ − R_−)²)`（J1(2m) ≈ m）。Ref C の値（B = 292 MHz、+10 / −4.7 dB、暗雑音 11 dB 下、m = 0.0785 rad）では S_Δ = 1.0×10⁻⁹ rad²/Hz。ループ帯域 10 kHz（仮置き）でスクイーズ側のロックは 3.2 mrad rms（ログ検出では 4.2 mrad）、反スクイーズ側のロックは 78 mrad rms。Ref C の「センサ雑音に制限されていない」と整合する。合計 0.1 rad の主因は dither（0.056 rad rms）。
   McKenzie 2005 の安定度の式（Eq. 13–14、(R_−/(R_+−R_−))^{1/2}・(1/Δω)^{1/4}）は、ロック点で分散を 2 次まで展開した静的な見積もりで、閉ループに入るセンサ雑音とはべき指数が違う。シミュレータは radiometer の PSD を使う（§6.4）。
4. **検出器の法則による radiometer 雑音の倍率** [MC]（平坦な 8–300 MHz、BT = 584）: 2 乗検出 0.97、線形包絡線 1.02、ログ（包絡線の 2 乗の対数を平均）1.32（いずれも 1/(BT) に対する比）。
5. **Liu 型 ＋ phase-detection OPA の統合は 2 段で組むのがよい** [AE]。(A)「phase-detection OPA で共通参照光とポンプをロック」＋「各ポートで逆向き probe と順方向の共通参照光の Liu 型ロック」なら、ポート k の光路の位相 ψ_k は打ち消され、共通の項も tap なしで固定できる。
   (B) 逆向き probe を squeezer まで戻してポンプと検出 OPA で直接比べると、ポート間では ψ_k は消えるが、GBS 側の参照に対してポンプの位相が係数 1（自然な値の 2 倍）、参照光の位相が係数 2 で残り、全ポートのアクチュエータがポンプのドリフトを追う。逆向きの光ではポンプに対する位相の符号が逆に入るため（§8.3）。
6. **時間多重の制約**: 誤差信号を振り分け素子の後でポートごとにサンプルすると、標本化レートは f_rep/32。f_rep = 10 MHz ならナイキスト周波数は 156 kHz で、1 MHz のヘテロダインはポートごとに復調できない（§8.4）。逆向き probe の飛行時間 n_g L_k/c（63 m で 0.31 µs、630 m で 3.1 µs）はループ遅延に加わる。

---

## 1. 記号（symbols, English）

General / detection
- `ℜ`: photodiode responsivity [A/W], `ℜ = η e/(hν)`; `η`: detection quantum efficiency (incl. path loss to the PD); `e`: elementary charge; `hν`: photon energy (1.2855×10⁻¹⁹ J at 1545.32 nm).
- `G_TIA`: transimpedance gain [V/A]. `i_n`: input-referred electronic current noise, one-sided [A/√Hz].
- `V_off`: DC offset at the demodulator/ADC input [V]. `K_d`: error-signal slope at the lock point [V/rad] (with respect to the named phase variable).
- `S_I(f)`: one-sided PSD of the photocurrent [A²/Hz]. `S_e(f)`: one-sided PSD of the demodulated error voltage [V²/Hz].
- `S_φ,sens(f) = S_e(f)/K_d²`: sensor noise referred to the locked phase [rad²/Hz] (conventions §7, one-sided).
- `b_loop`: equivalent noise bandwidth of the closed loop for sensor noise [Hz] (placeholder, from T16).

Phases (conventions §3–5, §9)
- `φ_p`, `φ_pr`, `φ_LO`, `φ_ref`: phases of the pump (complex amplitude at 2ω0), probe, LO, Liu-type reference, each evaluated at the named interaction point (e.g. `φ_p^det` = pump phase at the phase-detection OPA input).
- `θ_amp = φ_p/2` (amplified quadrature), `θ_sq = φ_p/2 + π/2` (squeezing angle), conventions §4.
- `Δ = φ_LO − θ_sq`: homodyne angle relative to the squeezing angle (Δ = 0 → squeezing), conventions §5.
- `u = φ_pr − φ_p/2`: probe angle relative to the amplified quadrature. `v = φ_LO − φ_p/2 = Δ + π/2`.
- `Φ1 = 2φ_pr − φ_p = 2u`: pump–probe relative phase seen by a parametric process (lock 1).
- `φ_OPA,i`: fixed + temperature-dependent phase offset of the parametric gain of OPA i (phase mismatch, focusing); generalises `θ_amp → (φ_p + φ_OPA,i)/2`.

Modulation / demodulation
- `f_sh`, `Ω_sh = 2π f_sh`: net optical frequency shift of the probe (Ref A: 1 MHz). Probe complex amplitude in the ω0 frame: `|A_pr| exp[i(φ_pr − Ω_sh t)]` (E ∝ A e^{−iω0 t}, conventions §4).
- `s_AOM = ±1`: diffraction order; an AOM adds `−s_AOM φ_RF(t)` to the optical phase in this sign convention (E ∝ e^{−iωt}; an up-shift multiplies by e^{−i(Ω t + φ_RF)}).
- `f_m`, `ω_m = 2π f_m`: dither frequency (Ref B 295 kHz, Ref C 1 MHz). `m`: phase-modulation amplitude [rad] (phase = φ + m sin ω_m t).
- `φ_d`: demodulation reference phase. Demodulator = multiply by `2 sin(ω_m t + φ_d)` (or `2 cos(...)`) and low-pass.
- `J_n`: Bessel function of the first kind.
- `a_RAM`: residual amplitude modulation depth (fractional power modulation at f_m), `φ_a` its phase.

Parametric gain
- Single-pass OPA (Ref A): `r = √(α P)` (Ref A's R± formula), amplitude gains `A± = e^{±r}`.
- OPO below threshold (Ref B), quasi-static (f_m ≪ cavity linewidth): amplitude gains `A+ = √G_+ = 1/(1−x)`, `A− = √G_− = 1/(1+x)`; `x` as defined by T02.
- `g_c = (A+ + A−)/2`, `g_s = (A+ − A−)/2` (= cosh r, sinh r for the OPA).
- `P_in`: probe power injected into the OPA/OPO that the PD sees (after taps/splits).

Quantum noise lock
- `S_SN(f)`: shot-noise PSD of the BHD output [V²/Hz]; `R±(f)`: anti-squeezed (+)/squeezed (−) noise normalised to shot noise (incl. loss); `S_el(f)`: electronic noise PSD; `H(f)`: detection-chain frequency response (T13).
- `S_x(f; Δ) = |H|² [S_SN (R_− cos²Δ + R_+ sin²Δ) + S_el]`: PSD of the BHD output.
- `f1, f2`, `B = f2 − f1`: detection band of the envelope/log detector (Ref C: 8–300 MHz).
- `P(Δ) = ∫_{f1}^{f2} S_x df`, `P± = P(π/2), P(0)` … see §6.1; `ΔP = (P_+ − P_−)/2`.
- `T`: averaging time; `κ_det`: detector-law factor of the radiometer noise (§6.3).
- `k_det` [V/W-equivalent] (square law) or `k_dB` [V/dB] (log detector).

Fringe / Liu-type
- `V`: fringe visibility, `ε`: imbalance of a balanced detector, `P_pr,det`, `P_ref`: probe/reference power at the interference detector, `D`: probe duty cycle.
- `ψ_f,k`, `ψ_b,k`: phase accumulated by forward (d = +1) / backward (d = −1) light in the port-k segments (conventions §9).

---

## 2. 共通の枠組み

### 2.1 ブロックの入出力

すべての誤差信号ブロックは、State から**ビームの位相**（conventions §9: 光源の位相＋通る区間の ψ の和、相互作用点で評価）とパワー・可視度などを読み、**電圧**を書く。
RF の搬送波（f_sh、2 f_sh、f_m）は prompt §5.4 (a) の複素包絡線で扱う。すなわちブロックは「理想復調後の DC 項」を位相から解析的に計算し、雑音は復調後の等価な基底帯域 PSD で加える。RF レートで明示的に計算するモードも検証用に持つ（§11 V2）。

### 2.2 復調で雑音がどう写るか

白色雑音（片側 PSD S）に `2 cos(ω t + φ_d)` を掛けて低域通過すると、f ± f_lp の 2 つの側帯域が折り返し、出力の片側 PSD は **2S** になる [AE、標準的な結果]。
したがって、PD 電流の雑音 S_I がキャリア付近で平坦なら:

- `S_e = 2 G_TIA² S_I`、`S_φ,sens = 2 S_I / K_I²`（K_I は電流での傾き [A/rad]）。
- DC ロック（復調なし）では `S_φ,sens = S_I / K_I²`。

### 2.3 オフセットとロック点のずれ

- 復調後の DC オフセット V_off: `δφ = V_off / K_d`。
- RAM などで「位相に依らない f_m の成分」が入ると、その同相成分 V_RAM cos φ_a が V_off と同じように働く。
- ヘテロダインでは誤差が `cos(Φ + φ_d)` の形なので、**復調位相のずれ δφ_d がそのままロック点のずれになる**（δΦ = −δφ_d）。dither 方式では、復調位相のずれは傾きを cos φ_d 倍にするだけで、オフセットは作らない（オフセットがあるときは、その大きさを変える）。

### 2.4 捕捉範囲

誤差が `sin(Φ − Φ0)` 型なら、捕捉範囲は安定点のまわり ±π（Φ について）。以下では、ロックする位相変数について書く。

### 2.5 検出器の仮置き値（T13 待ち）

| 量 | 仮置き値 | 根拠 |
|---|---|---|
| ℜ | 1.0 A/W | PDB450C の「typ. max responsivity」[manual]。1545 nm での値は T13 で確認 |
| G_TIA | 10³ V/A（Ref D、PDB450C 最大帯域設定）| [manual]。50 Ω 負荷では 1/2 |
| i_n | 10 pA/√Hz（BHD 以外の PD）、1 pA/√Hz（高利得の tap PD） | **仮置き（出典なし）**。T13 の値で置き換える |
| η（probe 系の PD） | 0.8 | **仮置き** |

---

## 3. #1 Ref A: AOM で周波数をずらした coherent probe ＋ ヘテロダイン復調

構成 [Ref A §3.1、Fig. 2]: 基本波をファイバカプラで LO と probe に分ける。probe は AOM 2 台（net 1 MHz シフト）と EOM（Control signal 1）を通り、コリメータから自由空間に出て、ダイクロイックミラーでポンプ（772.66 nm）と合わさる。
その直後の tapping mirror で、ポンプと probe を一緒に phase-detection OPA へ分ける。phase-detection OPA の出力 probe はファイバ結合されて PD へ行き（Phase-locking signal 1）、BHD の出力の一部が Phase-locking signal 2 になる。LO 経路にも EOM（Control signal 2）がある。

### 3.1 OPA の入出力（この規約で）

conventions §4 の H_int から、ω0 ± Ω の側帯波モードに対して（OPA の位相整合帯域 ≫ Ω_sh なので瞬時応答）:

```
a_out(t) = g_c a(t) + g_s e^{iφ_p} a*(t)        (OPA: g_c = cosh r, g_s = sinh r)
```

縮退の場合に S†âS = â cosh r − â† e^{iφ_ξ} sinh r と比べると φ_ξ = φ_p + π で、conventions §4 と一致する。
周波数シフトした probe `a(t) = |A_pr| e^{i(φ_pr − Ω_sh t)}` を入れると:

```
a_out(t) = |A_pr| [ g_c e^{i(φ_pr − Ω_sh t)} + g_s e^{i(φ_p − φ_pr + Ω_sh t)} ]
```

第 2 項が ω0 − Ω_sh の idler。

### 3.2 lock 1: phase-detection OPA の PD 信号

```
P_PD(t) = η_det P_in [ (g_c² + g_s²) + 2 g_c g_s cos(2Ω_sh t − Φ1^det) ]
        = η_det P_in [ cosh 2r_d + sinh 2r_d cos(2Ω_sh t − Φ1^det) ],     Φ1^det = 2φ_pr^det − φ_p^det − φ_OPA,det
e1 = G_TIA ℜ_d P_in sinh 2r_d · cos(Φ1^det + φ_d1)        (demod with 2cos(2Ω_sh t + φ_d1), LPF)
```

- **復調周波数は 2 f_sh = 2 MHz** [AE]。f_sh の成分は理想的には出ない（実際には AOM のもれや非線形で残りうるので、トグル付きのリーク項として持つ）。
  同じ原理（周波数シフトした control field を OPO に入れ、ポンプ位相の誤差を 2Ω で復調）は Vahlbruch et al. 2006 で確認した（本文: "demodulating the photo current at twice the beat frequency (80 MHz)"）。Ref A 本文には復調周波数の記述がない（UNVERIFIED）。
- 傾き: `K_d1 = G_TIA ℜ_d P_in sinh 2r_d` [V/rad of Φ1]。φ_pr に対しては 2 倍。
- ロック点: Φ1 + φ_d1 = ±π/2。**どの Φ1 にロックするかは φ_d1 で決まる**（2.3）。
- 捕捉範囲: Φ1 について ±π（φ_pr について ±π/2）。
- ショット雑音の床: `S_Φ1,shot = 4 hν cosh 2r_d / (η_d P_in sinh² 2r_d)`。電気雑音: `S_Φ1,el = 2 i_n² / (ℜ_d P_in sinh 2r_d)²`。
- 概算（仮置き: P_in = 1 mW（Ref A §2「数 mW まで」）、検出 OPA の利得 e^{2r_d} = 10、η = 0.8、G_TIA = 10³、i_n = 10 pA/√Hz）:
  K_d1 = 4.95 V/rad、PD の DC 出力 5.05 V（飽和に注意）、S_Φ1,shot = 1.3×10⁻¹⁶ rad²/Hz、S_Φ1,el = 8×10⁻¹⁸ rad²/Hz、100 kHz で積分して 4 µrad rms。**無視できる。**
  検出 OPA のポンプパワー・利得、PD の仕様は Ref A に書かれていない（PD_OPA_PARAMS はプレースホルダ）。

### 3.3 従来の tap 方式（Ref A Fig. 1(a)、prompt §7-8）

squeezer の出力を割合 T_tap だけ PD に分ける。式は 3.2 と同じで、P_in → T_tap P_pr,sq、r_d → r（squeezer の r）。tap は squeezed 光に損失 T_tap を加える: L → 1 − (1 − L)(1 − T_tap)。

Ref A の値（r = √(9.06 × 0.64) = 2.408、e^{2r} = 123.5（20.9 dB）、L = 8%、θ = 9 mrad、P_pr,sq = 5 µW（Ref A §2「数 µW」）、η = 0.8）[AE]:

| T_tap | S_Φ1,shot [rad²/Hz] | i_n = 10 pA/√Hz のとき 100 kHz で積分した値 | R′− [dB]（θ = 9 mrad 固定） |
|---|---|---|---|
| 0.1% | 2.1×10⁻¹² | 14 mrad rms | −10.11 |
| 1% | 2.1×10⁻¹³ | 1.4 mrad rms | −9.76 |
| 3% | 6.9×10⁻¹⁴ | 0.48 mrad rms | −9.07 |

tap を増やすと SNR は上がるが、squeezing は損失で直接下がる（prompt §7-8 の trade-off）。
ただし §3.4 のとおり、lock 1 の雑音は Δ にほとんど伝わらないので、tap 方式で SNR が問題になるのは「lock 1 が外れない程度」の水準だけ [AE]。

### 3.4 lock 2: BHD での probe–LO ビート（f_sh）と、lock 1 の誤差が縮む仕組み

squeezer を通った probe（信号 ＋ idler）が LO と BHD で干渉する。conventions §5 の `⟨n1 − n2⟩ = 2 Re(α_LO* α)` から:

```
i_−(t) ∝ 2|α_LO||A_pr| [ g_c cos(Ω_sh t − (φ_pr − φ_LO)) + g_s cos(Ω_sh t + (φ_p − φ_pr − φ_LO)) ]
e2 (demod 2cos(Ω_sh t + φ_d2)) = K_2 [ g_c cos(u′ − v) + g_s cos(u′ + v) ]
     = K_2 [ A+ cos u′ cos v + A− sin u′ sin v ],      u′ = u + φ_d2,  v = Δ + π/2,
K_2 = G_BHD · 2ℜ √(P_LO P_pr,sq η_pr)
```

（φ_p、u、v は squeezer での値。）

- u′ = 0 に選ぶと e2 = K_2 A+ cos v。零点は v = π/2、すなわち **Δ = 0（squeezing）**で、傾きは A+ = e^{r}（増幅側）。u′ = π/2 なら e2 = K_2 A− sin v で、零点は Δ = −π/2（反 squeezing）、傾きは A−。
- **lock 1 の誤差の伝わり方** [AE]: u′ = δu′ のとき、零点は `tan ε = (A−/A+) tan δu′`、ε = Δ のずれ。
  lock 1 の残差、lock 1 のセンサ雑音、2 つの OPA の非共通光路、lock 2 の復調位相 φ_d2 のドリフトは、**すべて A−/A+ 倍されてから Δ に入る**（Ref A: 1/123.5。δu′ = 0.1 rad → ε = 0.81 mrad）。
  直観的には、利得の大きい OPA を通った probe は入力の位相にかかわらず増幅軸の向きに揃うので、それを基準にした LO のロックはポンプの位相を参照することになる。
- 帰結: lock 1 に必要なのは u′ を ±π/2 から十分遠ざけて傾きの符号と大きさを保つことで、精度ではない。**Ref A は「lock 1 の SNR を上げて 14 → 9 mrad」と説明しているが、この導出では lock 1 のセンサ雑音は Δ の主因になりえない。** 食い違いの理由（Ref A の lock 2 が片側の側帯波だけを使っている、lock 2 の帯域、など）は本文から特定できない。Phase 5 でレビュアーが確認し、シミュレータでは §11 V3 で数値的に確かめる。
- lock 2 のセンサ雑音: ロック点（Δ = 0）での BHD 出力の雑音は f_sh 付近の squeezed 雑音（R_− × ショット雑音）なので、
  `S_Δ,sens = hν R(f_sh) / (η P_pr,sq A+²)`。
  Ref A Fig. 4(a) の 1 MHz のピーク（反 squeezing の軌跡で、ショット雑音より約 41 dB 上、RBW 1 MHz。図からの読み取りで ±3 dB）から `P_pr,sq A+² ≈ 1.6 nW` [AE]。
  R = 0.1 なら S_Δ = 7.9×10⁻¹² rad²/Hz、10 kHz で 0.28 mrad rms。R = 1 なら 0.89 mrad rms。9 mrad より小さい。

### 3.5 phase-detection OPA までの非共通光路（E68）

lock 1 が固定するのは Φ1^det。squeezer で効くのは Φ1^sq。両者の差は:

```
Φ1^sq − Φ1^det = 2(ψ_pr,sq − ψ_pr,det) − (ψ_p,sq − ψ_p,det) − (φ_OPA,sq − φ_OPA,det)
```

ψ は tapping mirror から各 OPA の入力までの区間。ポンプと probe は同じ区間を同軸に通るので、長さが δL 変わったときは分散の分しか残らない:

```
δ(2ψ_pr − ψ_p) = (2ω0/c) (n(ω0) − n(2ω0)) δL
```

- 空気（Peck–Reeder の分散式で計算 [AE]）: n(1545 nm) − n(773 nm) = −1.96×10⁻⁶ → 16 rad/m。δL = 1 µm で 16 µrad。
  10 cm の非共通区間の空気が 1 K 変わると、密度の変化で 5.4 mrad（Φ1）。
- 溶融石英（Malitson の Sellmeier 式 [AE]）: n(ω0) − n(2ω0) = −9.7×10⁻³ → 7.9×10⁴ rad/m、**1 µm あたり 79 mrad**。
  ポンプと probe をファイバで別の場所の検出 OPA に運ぶ構成（F で検討する場合）では、非共通のファイバ長に強く効く。熱光学係数の分散は T08 の値を使う（UNVERIFIED）。
- φ_OPA,det と φ_OPA,sq は、それぞれの導波路の温度・位相不整合で独立にドリフトする（T01 の値を使う）。
- ここに挙げた量はすべて、§3.4 により Δ には A−/A+ 倍されて入る。**ただし F（BHD による lock 2 がなく、GBS に直接入る）では縮まない**（§8.3）。

### 3.6 RF 参照の位相（E74）

- probe の光の位相には AOM の RF の位相が入る: φ_pr ∋ −Σ s_AOM,j φ_RF,j(t)。復調の参照を同じ RF 源から作れば、時刻のずれ τ_RF の分だけ残る: `S = 4 sin²(π f τ_RF) S_φRF(f)`（conventions §7 と同じ形）。
- lock 1 では 2 f_sh の参照の位相 φ_d1 が 1:1 でロック点を動かす（Φ1 について）。lock 2 の φ_d2 は §3.4 により A−/A+ 倍される。
- 参照の作り方（2 台の AOM の駆動の差周波を混合して作るか、共通クロックの DDS か）は Ref A に書かれていない（UNVERIFIED）。パラメータとする。

### 3.7 BHD・GBS への probe の混入（E72）

- Ref A: squeezer 側の probe（f_sh = 1 MHz のビート）は BHD に入り、Fig. 4(a) に 1 MHz のピークとして見える。測定帯（3 MHz、RBW 1 MHz）の外。測定帯に入るのは probe の古典雑音（RIN、位相雑音の f_sh ± 3 MHz 成分）だけ。
- tap 方式の probe は「数 µW」に制限される（Ref A §2）。phase-detection OPA 方式では、検出 OPA 側の probe を数 mW にできる。squeezer 側の probe パワーとの分け方は Ref A に書かれていない。
- GBS（光子計数）では、squeezed モードに入ったコヒーレントな probe は変位として効く。スロットあたりの混入光子数は T10 の式（D28）で出す。このブロックは「squeezer を通る probe のパワーと周波数」を State に出力する。

---

## 4. #2 Ref B lock 1: dither ＋ 透過 probe の DC 強度の極値ロック

[Ref B §2–3] probe を PZT2 で 295 kHz の位相変調。OPO を透過した probe の強度を PD で見て、PZT1 で最大増幅か最大減衰にロックする。
Fig. 1 では OPO の出力の後に「99%R」があり、PD への光が出力の 1% の tap なのか、別の鏡からのもれなのかは本文から判断できない（UNVERIFIED。tap なら squeezed 光に 1% の損失）。

OPO の線幅（γ_c ≈ 7.6×10⁷ s⁻¹、survey E08）≫ 295 kHz なので準静的に扱う:

```
P_PD(t) = η P_in [ Ḡ + ΔG cos Φ1(t) ],   Φ1(t) = Φ1,0 + 2m sin ω_m t      (Φ1 = 0: max gain)
Ḡ = (G_+ + G_−)/2,  ΔG = (G_+ − G_−)/2      (from |g_c e^{iu} + g_s e^{−iu}|² = Ḡ + ΔG cos 2u, Φ1 = 2u)
first harmonic:  −2 ηP_in ΔG J1(2m) sin Φ1,0 · sin ω_m t
e1 = −2 G_TIA ℜ η P_in ΔG J1(2m) sin(Φ1,0) cos φ_d      (demod 2 sin(ω_m t + φ_d))
```

- ロック点: Φ1,0 = 0（最大増幅）または π（最大減衰）。1 つの符号ではどちらか一方が安定。
- 傾き: `K_d = 2 G_TIA ℜ η P_in ΔG J1(2m)`（Φ1 について）。捕捉範囲 ±π（Φ1）。
- dither による 2 倍波 `2J2(2m) cos Φ1,0 cos 2ω_m t` は LPF で落とす。極値でも DC は `J0(2m)` だけ下がる。
- センサ雑音: `S_Φ1 = 2 S_I / K_I²`、`S_I = 2eℜηP_in(Ḡ ± ΔG)`（ロック点の DC で決まる）。
- 値: G = 18.7 [Ref B §4] → x = 1 − 1/√18.7 = 0.769、G_− = 1/(1+x)² = 0.320（Ref B の定義 G = 1/(1−x)² を採用。x の定義の食い違いは T02）。
  m、P_in、PD の値は Ref B に書かれていない（**仮置き**: パラメータ）。
- オフセット: PZT のミラーの傾きによる RAM（a_RAM）→ `δΦ1 = G(Φ_lock) a_RAM cos φ_a / (2 ΔG J1(2m))`。G(Φ_lock) はロック点での利得（最大増幅なら G_+、最大減衰なら G_−）。最大減衰のロック点のほうが RAM に強い。

---

## 5. #3 Ref B lock 2: BHD での probe–LO 干渉 ＋ dither 復調

probe は OPO の中で squeezed 光と同じモードを通り（同軸）、BHD で LO と干渉する。dither は probe だけにかかっている（PZT2 は probe 経路。squeezed 光と LO には dither がかからない）。

```
⟨i_−⟩ ∝ 2|α_LO||A_pr| [ g_c cos(u − v) + g_s cos(u + v) ],   u = u0 + m sin ω_m t
e2 ∝ 2J1(m) · ∂/∂u [ g_c cos(u − v) + g_s cos(u + v) ] = −2J1(m) [ g_c sin(u0 − v) + g_s sin(u0 + v) ]
```

- lock 1 が最大減衰（u0 = π/2）: e2 ∝ −2J1(m) A+ cos v。零点 v = π/2、すなわち **Δ = 0（squeezing）**、傾き ∝ A+ = 1/(1−x)。
- lock 1 が最大増幅（u0 = 0）: e2 ∝ 2J1(m) A− sin v。零点 Δ = −π/2（反 squeezing）、傾き ∝ A−。
- DC 成分（復調しない）は J0(m) [g_c cos(u0−v) + g_s cos(u0+v)]。u0 = 0 で零点 v = π/2。ただし DC ロックは強度・オフセットに弱い（§7）。
- **lock 1 の誤差の伝わり方**（§3.4 と同じ式）: `tan ε = (A−/A+) tan δu` = 0.131 tan δu（Ref B の値）。
- dither 側帯波（±295 kHz）は BHD に入る。測定は 1 MHz（RBW 30 kHz）なので帯域外。Ref B は「測定周波数でショット雑音が影響を受けないことを確認した」と書いている（§3）。
- センサ雑音: `S_Δ = 2 S_I / K_I²`。S_I はロック点での BHD 出力の 295 kHz 付近の雑音で、レーザー雑音（89・177 kHz のピーク）の影響を受ける（そのために 295 kHz を選んだ、Ref B §3）。

---

## 6. #4 Ref C: 量子雑音ロック（probe なし）

構成 [Ref C §2、Fig. 3]: LO に fiber EOM（MPZ-LN-10、Uπ = 4 V）で 1 MHz の位相変調（200 mVpp → ±4.5°）。BHD の高域通過出力 → 利得段（Fig. 3 に「x50」）→ ログアンプを含む包絡線検出器（8–300 MHz の量子雑音レベルに比例する電圧）→ Moku:Pro で 1 MHz で復調 → 3 段の PID。

### 6.1 雑音パワーの位相依存性（帯域で積分）

```
P(Δ) = ∫_{f1}^{f2} |H(f)|² [ S_SN(f) (R_−(f) cos²Δ + R_+(f) sin²Δ) + S_el(f) ] df
     = P̄ − ΔP cos 2Δ,     P_∓ = P(0), P(π/2),  P̄ = (P_+ + P_−)/2,  ΔP = (P_+ − P_−)/2
```

R±(f) は T02（OPO の線幅、Ref C は当てはめで 1.75 GHz）と T12（損失）から、H と S_el は T13 から受け取る。

- squeezing 角が周波数で回る場合（離調、分散）は `R_− cos²(Δ − δ(f)) + R_+ sin²(Δ − δ(f))` とし、ロック点は帯域の重み付き平均 `Δ* = ½ atan2(∫w sin 2δ, ∫w cos 2δ)`、w = |H|² S_SN (R_+ − R_−) になる [AE]。
  ESA の測定周波数での最適角とはずれる。ファイバの GVD によるずれは β2 = −22 ps²/km（UNVERIFIED）・1 km・300 MHz で 4×10⁻⁵ rad で無視できる [AE]。

### 6.2 dither と誤差信号

LO に dither をかけると Δ(t) = Δ0 + m sin ω_m t。conventions §5 のとおり LO の位相が増えると Δ が増える。

```
P(t) = P̄ − ΔP [ J0(2m) cos 2Δ0 − 2 J1(2m) sin 2Δ0 sin ω_m t + 2 J2(2m) cos 2Δ0 cos 2ω_m t − … ]
square-law detector (v = k_det P):    e = 2 k_det ΔP J1(2m) sin(2Δ0) cos φ_d
log detector (v = k_dB · 10 log10 P): e = (10/ln 10) k_dB · 2 ΔP J1(2m) sin(2Δ0) cos φ_d / P(Δ0)
```

- 零点は Δ0 = 0（squeezing）と π/2（反 squeezing）。1 つの符号ではどちらか一方が安定。捕捉範囲は Δ について ±π/2。
- 傾き（Δ0 = 0、2 乗検出）: `K_d = 4 k_det ΔP J1(2m)`。ログ検出では 1/P_− 倍、反 squeezing 側では 1/P_+ 倍で、**ループ利得が 2 つのロック点で P_+/P_− 倍（Ref C で約 14 dB）違う**。SNR はログ検出でも変わらない（信号と雑音が同じ 1/P で変わる）。
- McKenzie Eq. (7) の `J0(θ1) J1(θ1) sin 2θ0` は、ここでの `J1(2m)/2 · sin 2Δ0` を 1 次まで展開したもの（3 次で違う）。McKenzie は真空の分散を 1、θ0 = π/2 を振幅 squeezing の squeezing 側と置いている。この規約では R±（ショット雑音で規格化）と Δ に置き換える。
- dither そのものは測定される直交位相を ±m 揺らす。Ref C の m = 0.0785 rad（±4.5°）は 0.0555 rad rms（prompt §7-10）。squeezing への写像は T04（E65）。

### 6.3 センサ雑音: radiometer 型

BHD 出力を平均 0 のガウス雑音（片側 PSD S_x）とする。2 乗検出 y = x² の低周波の PSD は（Isserlis の定理から C_y(τ) = 2C_x(τ)²）:

```
S_P̂(f ≪ B) = 2 ∫_0^∞ S_x(f′)² df′          [one-sided, (power unit)²/Hz]
flat band (S_x = N on [f1, f2]):  S_P̂ = 2 N² B = (2/B) P²   → boxcar of length T: Var(P̂)/P² = 1/(B T)   (Dicke radiometer equation)
```

復調（§2.2）で 2 倍になり、`S_e = 2 k_det² S_P̂`。位相に換算して:

```
S_Δ,sens(Δ0) = ∫ S_x(f; Δ0)² df / ( 4 ΔP² J1(2m)² cos²φ_d ) · κ_det
flat band:     S_Δ,sens(0)   ≈ κ_det R_−² / ( B J1(2m)² (R_+ − R_−)² )
               S_Δ,sens(π/2) ≈ κ_det R_+² / ( B J1(2m)² (R_+ − R_−)² )
```

ここで R± は回路雑音を含めた値（R± + S_el/S_SN）。

- **検出器の法則の倍率 κ_det** [MC、平坦な 8–300 MHz、T = 2 µs（BT = 584）、約 2000 ブロック]: 2 乗検出 0.97、線形包絡線（4 × 相対分散）1.02、ログ（包絡線の 2 乗の対数を平均）1.32。
  独立な標本の分布から出る値（指数分布の ln の分散 π²/6 = 1.64、Rayleigh 分布で 1.09）は、標本間の相関があるので当てはまらない。シミュレータは κ_det をパラメータにし、§11 V6 で検証する。
- 条件: f_m ≪ B、ループ帯域 ≪ f_m、検出器のビデオ帯域 ≫ f_m、BT ≫ 1（ガウス近似）。
- **Ref C の値で概算** [AE]: R_+ = 10 dB、R_− = −4.7 dB（Fig. 4、30 MHz、暗雑音補正後）、暗雑音はショット雑音の −11 dB（30 MHz の値を帯域全体に使う。実際は 0.25 GHz で −7.4 dB に悪化するので楽観的）、B = 292 MHz、m = 0.0785 rad（J1(2m) = 0.0783）、b_loop = 10 kHz（**仮置き**。Ref C は「fiber EOM が kHz 域を補正」とだけ書いている）:
  S_Δ,sens(0) = 1.05×10⁻⁹ rad²/Hz → **3.2 mrad rms**（κ_det = 1.32 で 4.2 mrad）。反 squeezing 側のロックでは S = 6.1×10⁻⁷ rad²/Hz → 78 mrad rms。
  Ref C の「センサ雑音にも（遅いアクチュエータ以外の）レンジにも制限されていない」（§4）と整合する。合計 0.1 rad の主因は dither の 0.056 rad rms。
- **McKenzie 2005 との違い**: McKenzie Eq. (11)–(14) は「ロック点のまわりで分散を 2 次まで展開し、それが分散の揺らぎ（kurtosis）に等しくなる位相」を安定度とする静的な見積もりで、δθ ∝ (R_−/(R_+ − R_−))^{1/2} (1/Δω)^{1/4}。
  閉ループでは誤差信号（dither の 1 次）に加わる雑音がそのまま位相雑音になるので、上の式（δΔ ∝ R_−/(R_+ − R_−) · (b_loop/B)^{1/2} / J1(2m)）が適切 [AE]。
  定性的な結論（squeezing 側のロックのほうが安定、損失・検出効率で悪化、帯域を広げると良くなる）は一致する。べき指数が違う。シミュレータは radiometer の PSD を使い、McKenzie の式は比較のためにだけ出す。

### 6.4 オフセット

- EOM の RAM（LO のパワーが a_RAM で変調される）: ショット雑音・squeezed 雑音はどちらも P_LO に比例するので、`P(t) ⊃ P̄_lock a_RAM sin(ω_m t + φ_a)` になる。
  `δΔ_off = a_RAM P(Δ0) cos φ_a / (4 ΔP J1(2m)) ≈ a_RAM R_− / (2 (R_+ − R_−) J1(2m))`（Δ0 = 0、回路雑音は除く）。
  Ref C の値で a_RAM = 10⁻³ → 0.28 mrad、10⁻² → 2.8 mrad [AE]。a_RAM は**仮置き**（T16 の EOM の値）。
- 検出器・ログアンプの出力の f_m のピックアップ、復調器の DC オフセット: §2.3。
- 回路雑音・LO パワーのゆっくりした変化はロック点を動かさず、傾き（ループ利得）だけを変える（ログ検出では打ち消される）。

### 6.5 漏れ

probe はない。BHD の出力に入る余分な光はない。dither の位相変調が LO にかかる（6.2）。LO のファイバと squeezed 光のファイバの雑音の相関は T08。

---

## 7. #5 DC fringe lock（汎用ブロック、D06。Ref D の詳細は T15）

conventions §5 の DC 干渉に、平衡のずれ・可視度・オフセットを加える:

```
v = G_TIA ℜ η P_tot [ ε + V cos(φ_a − φ_b) ] + V_off + n(t),       K_d = −G_TIA ℜ η P_tot V sin(φ_set)
```

- 入力: 2 つのビームの位相 φ_a、φ_b（conventions §9 の線形結合）、P_tot(t)（RIN を含む）、V(t)（偏光・モード、T09・T12）、ε、V_off(t)、雑音 i_n。出力: 電圧。
- 設定値 φ_set = π/2（Ref D の初期値、D02）。強度雑音は ε を通してだけ入る（φ_set = π/2 のとき）。
- Ref D の上限の概算 [AE]: ℜ = 1.0 A/W、G = 10³ V/A、P_tot = 3.33 mW、η = V = 1 → K_d = 3.33 V/rad。
  ショット雑音（2 つの PD の合計電流）`S_φ = 2e/(ℜ P_tot V²)` = 9.6×10⁻¹⁷ rad²/Hz。DSC1 の入力雑音（0.3 mV rms、100 Hz–250 kHz [manual]）を白色とみなすと 3.2×10⁻¹⁴ rad²/Hz で、こちらが支配的。詳細は T15。
- 捕捉範囲: ±π（φ_a − φ_b）。逆の傾きの点に落ちないように、ランプによる再引き込み（T16）が要る。

---

## 8. #6 Liu 型: 逆向きの probe を共通の参照光と干渉させる（D30）

### 8.1 Liu et al. 2026 の記述

[Liu 2026 Methods「Phase stabilization of the system」]:
- 3 つのサーボ。(1) ポンプの分岐から第 1 干渉計まで: OPO を通ったポンプをダイクロイックで分け、参照光と干渉させ、**平衡ホモダイン検出**で誤差信号を作り PZT ミラーを駆動。
- (2)(3) 2 つのファイバ遅延ループの配列: 各空間モードの光路に PZT ミラー。**squeezed 光と周波数が同じで偏光が直交する逆向きの probe** を 16 本に分け、第 2 干渉計の入力で PBS により各モードに入れる。迷光の計数を減らすため、probe は **AOM でパルス化し、squeezed 光のパルス列から時間をずらす**。第 1 干渉計の出力で別の PBS により probe を取り出し、**参照光と干渉させて誤差信号を作り PZT を駆動**。長いファイバの温度ドリフト（1 つの PZT のレンジを超えることが多い）には TEC のフィードバックを併用。
- 精度: OPO 光源の配列 1 nm、第 1 遅延ループ 3 nm、第 2 遅延ループ 10 nm（1550 nm で 4.1、12、41 mrad に相当 [AE]。rms か ptp かは書かれていない、UNVERIFIED）。数十時間連続で維持。
- 参照光との干渉の方式（DC ホモダインか、変調・復調か）、probe と参照光のパワー、検出器は**書かれていない**（UNVERIFIED）。
- 繰り返し 20 MHz、遅延の単位 τ = 50 ns。

### 8.2 ブロックの式

ポート k の probe（共通の分岐点から GBS 側で入り、ポート k を逆向きに通る）と、順方向に配られた参照光の位相:

```
φ_pr,k = φ_L(t − τ_pr,k) + ψ_split,k + Σ_{seg∈port k} ψ_seg,a⊥(t at backward crossing)        (d = −1, axis a⊥)
φ_ref,k = φ_L(t − τ_ref,k) + ψ_refdist,k                                                       (d = +1)
```

検出方式はパラメータにする（Liu の論文で書かれていないため）:

1. `homodyne_dc`（δf = 0）: §7 の DC fringe ブロックと同じ式で、`P_tot V → 2√(P_pr,det P_ref) V_pr,ref`（強い参照光の平衡検出）。
   `K_I = 2ℜ√(P_pr,det P_ref) V`、`S_φ,shot = hν/(2η P_pr,det)`（参照光のショット雑音が支配的）。
2. `heterodyne`（probe の周波数を δf ずらす）: §3.4 と同じ I/Q 復調。δf の RF 参照の位相（§3.6）が 1:1 で入る（後ろに増幅器がないので縮まない）。
3. `dither`: 参照光か probe に位相変調し、§5 と同じ形で復調。

パルス化（デューティ D）:
- 検出をゲートしない場合: 傾きは D 倍、参照光のショット雑音は常にある → `S_φ = hν/(2η P̄_pr D)`（P̄_pr は平均パワー）。
- ゲートする（パルスの間だけ標本化）場合: `S_φ = hν/(2η P̄_pr)`。
- 標本化レートは probe パルスの繰り返し。DC ロックのドリフト・オフセットを避けるため、パルスの on/off の差を取る同期検出もパラメータとして持つ。

F での概算（**すべて仮置き**、D28 で probe 仕様は Open）: P_pr,det = 1 µW（平均、ゲートあり）、P_ref = 1 mW、η = 0.8、V = 1、G_TIA = 10⁴、i_n = 10 pA/√Hz →
K_d = 0.63 V/rad、S_φ,shot = 8.0×10⁻¹⁴ rad²/Hz、S_φ,el = 2.5×10⁻¹⁴ rad²/Hz。10 kHz で積分して 30 µrad rms で、Liu の報告値（mrad–数十 mrad）は誤差信号の雑音では決まっていないと考えられる [AE]。

### 8.3 Liu 型の打ち消しと、phase-detection OPA との統合（ユーザの検討案）

逆向き probe と順方向の参照光の干渉で `φ_pr,k − φ_ref,k = C_k` に固定すると、ψ_f,k = ψ_b,k = ψ_k（同じ軸、相反、遅い雑音）のとき:

```
θ_sq,k(GBS) − φ_pr,0(GBS) = [θ_sq(demux) − φ_ref(demux)] + C_k + (common-point terms)
```

ポート k 固有の ψ_k は消え、残るのは全ポートに共通の項 `θ_sq(demux) − φ_ref(demux)` と共通点の揺れ（D31-2）。共通の項は指標 (ii) では時間多重の時刻の差でしか効かない（D31-3）。
ψ_f ≠ ψ_b となる原因（直交偏光の複屈折ドリフト、逆向き伝搬の時刻のずれ c/(2 n_g L) より上の雑音、Rayleigh・PER による混入）は T10。

**統合案の組み方** [AE]:

- (A) 正しい組み方: phase-detection OPA（§3.2）で「共通参照光（probe/LO の元）とポンプ」を squeezer の近くでロック（lock 1）し、各ポートで Liu 型ロック（lock B_k）をかける。
  共通の項 `θ_sq(demux) − φ_ref(demux)` を、squeezed 光を tap せずに固定できる。指標 (i)（LO_k に対する位相）や外部の LO を使う測定で必要になる。
  非共通の区間（tapping 点 → 検出 OPA と tapping 点 → squeezer、squeezer → demux と参照光の配送）は、§3.5 の分散の係数で効く。F では BHD の lock 2 による A−/A+ の縮小が**ない**ので、§3.5 の量がそのまま残る。
- (B) 不利な組み方: 逆向き probe を振り分け素子と squeezer を逆向きに通して入力側まで戻し（単一パスの OPA では逆向きの光は位相整合せず、増幅されない）、ポンプと一緒に phase-detection OPA に入れて `2φ_pr,k^back − φ_p = C_k` を固定する。
  ψ_b,k = φ_p/2 − φ_pr,0 + C_k/2 となり、θ_sq,k − φ_pr,0 = 2(φ_p/2 − φ_pr,0) + C_k/2 + π/2 + (ψ_f,k − ψ_b,k)。
  ポート間の差（指標 (ii)）では ψ_k は消えるが、(1) GBS 側の参照に対しては**ポンプの位相の揺らぎが係数 1（自然な値 1/2 の 2 倍）、参照光の揺らぎが係数 2 で残る**（逆向きの光ではポンプに対する位相の符号が逆に入るため）。(2) 全ポートのアクチュエータがポンプのドリフトに追従するのでレンジを食う。(3) 検出 OPA は 1 台をスロットで共有することになり、§8.4 の標本化の制約（f_rep/32）がかかる。
- シミュレータでは (A)・(B) をどちらも設定で組めるようにし、§11 V10 で符号と係数を検証する（per-path 位相モデルで「仮定」ではなく「計算」する、D11）。

### 8.4 時間多重による制約

- 振り分け素子が時間スイッチで、probe がスイッチを逆向きに通ってから検出される場合、ポート k の誤差信号は k のスロットでだけ得られる。ポートあたりの標本化は f_rep/32、ナイキストは f_rep/64（f_rep = 10 MHz で 156 kHz、100 MHz で 1.56 MHz）。f_sh・f_m はこれより十分低くする必要がある。D30 のとおり**ポートごとに振り分け素子の手前で検出**すれば、この制約はない（probe パルスの繰り返しだけで決まる）。
- 逆向き probe の飛行時間 `n_g L_k / c`（L = 2.04 m: 10 ns、63 m: 0.31 µs、630 m: 3.1 µs）は、アクチュエータ（ポート k のストレッチャ）の位置から検出点までの遅延としてループ遅延に加わる。区間の位置ごとの遅延は conventions §9 の遅延バッファで表す。
  逆向き伝搬で雑音の相関が切れる周波数 c/(2 n_g L)（2.04 m: 50 MHz、63 m: 1.6 MHz、630 m: 0.16 MHz）は T10 の扱い。

### 8.5 GBS への probe の混入

逆向きで偏光が直交、時間をずらしたパルスなので、直接は GBS に入らない。混入の経路は (a) PBS の消光比・ファイバの PER、(b) Rayleigh 後方散乱（順方向に戻る、同じ周波数）、(c) コネクタの反射（FC/APC）、(d) AOM の消光比による probe パルスの裾がスロットに重なる分。
このブロックは probe のパワー・タイミング・偏光・周波数を State に出力し、混入光子数は T10 が計算する。

---

## 9. 方式ごとのまとめ

| # | 方式 | ロックする位相 | K_d（式） | センサ雑音の床 | 主なオフセット | 捕捉範囲 | BHD/GBS への混入 |
|---|---|---|---|---|---|---|---|
| 1a | Ref A lock 1（検出 OPA、2 f_sh） | Φ1^det = 2φ_pr − φ_p − φ_OPA | G ℜ P_in sinh 2r_d | 4hν cosh 2r_d/(ηP_in sinh² 2r_d)、+ 2i_n²/K_I² | φ_d1（1:1）、非共通光路 §3.5（Δ には A−/A+ 倍） | ±π（Φ1） | 検出 OPA 側は混入なし |
| 1b | 従来の tap | Φ1^sq | G ℜ T_tap P_pr sinh 2r | 同上で P_in = T_tap P_pr | 同上 | ±π | tap 損失 T_tap、probe は数 µW |
| 1c | Ref A lock 2（BHD、f_sh） | Δ（v） | G_BHD 2ℜ√(P_LO P_pr) A+ | hν R(f_sh)/(η P_pr A+²) | φ_d2・lock 1 の誤差は A−/A+ 倍 | ±π/2（Δ） | f_sh のビート（帯域外）、probe の古典雑音 |
| 2 | Ref B lock 1（dither、DC 極値） | Φ1 | 2Gℜη P_in ΔG J1(2m) | 2S_I/K_I² | RAM: Ḡ a_RAM/(2ΔG J1) | ±π（Φ1） | PD が tap なら 1%（UNVERIFIED） |
| 3 | Ref B lock 2（BHD、dither） | Δ | ∝ 2J1(m) A+ | 2S_I/K_I² | lock 1 の誤差は A−/A+ 倍 | ±π/2（Δ） | probe ＋ ±295 kHz 側帯波（帯域外） |
| 4 | Ref C 量子雑音ロック | Δ | 4 k ΔP J1(2m)（ログでは /P(Δ0)） | κ_det ∫S_x²/(4ΔP² J1²) | RAM: a_RAM R_−/(2(R_+−R_−)J1) | ±π/2（Δ） | なし（dither が LO に） |
| 5 | DC fringe | φ_a − φ_b | G ℜ η P V sin φ_set | 2e/(ℜPV²) ＋ 回路 | V_off/K_d、ε × RIN | ±π | —（Ref D は coherent のみ） |
| 6 | Liu 型 | φ_pr,k − φ_ref,k | 2Gℜ√(P_pr P_ref) V | hν/(2ηP̄_pr)（ゲートあり） | V_off/K_d、δf の RF 位相 | ±π | PER・Rayleigh・反射・AOM 消光比（T10） |

---

## 10. 周波数帯とスケーリング

周波数帯:
- 誤差信号のキャリア: Ref A 1 MHz（lock 2）・2 MHz（lock 1）、Ref B 295 kHz、Ref C 1 MHz（検出帯 8–300 MHz）、DC fringe と Liu 型（DC）は DC – ループ帯域。
- 誤差信号が効く（センサ雑音・オフセットが位相に入る）のは DC からループ帯域（kHz – 数十 kHz、T16）まで。その上ではセンサ雑音は抑圧されないが、注入もされない。

スケーリング:
- パワー: coherent 方式の K_d ∝ P（または √(P_pr P_ref)）、ショット雑音の床 ∝ 1/P_pr。量子雑音ロックは LO パワーに依らない（ショット雑音が回路雑音より十分大きい限り）。悪化は R_−/(R_+ − R_−) と 1/B で決まる。
- 利得（ポンプ）: lock 1 の傾き ∝ sinh 2r。lock 2 の縮小率 A−/A+ はポンプを上げるほど小さい。量子雑音ロックは R_+ − R_− とともに改善するが、損失で R_− が下がらなくなると飽和する。
- ファイバ長: 誤差信号の生成そのものは L に依らない。ただし (a) probe の受光パワーがファイバ損失で下がる、(b) 量子雑音ロックでは損失で R± が 1 に近づく、(c) Liu 型では飛行時間がループ遅延に加わり（∝ L）、相関の切れる周波数が下がる（∝ 1/L）、(d) 非共通のファイバ区間は §3.5 で ∝ L。
- 温度: 非共通区間の分散による項（空気: 10 cm・1 K で 5.4 mrad、Φ1）、OPA の位相不整合のドリフト（T01）、RF ケーブルの遅延のドリフト（§3.6）、オフセットのドリフト。
- パルス幅: probe・参照光の時間モードの重なりで V が下がる（T03、分散 T06）。パルスのゲートとデューティ D（§8.2）。パルス励起の OPA では r がピークパワーで決まる（T01・T03）。

---

## 11. 検証テスト（実装が通すべきもの）

V1. **Ref A lock 1 の誤差信号の形**: 瞬時応答の OPA（§3.1）で P_PD(t) を RF レートで作り 2 f_sh で復調する。`e1(Φ1) = G ℜ P_in sinh 2r cos(Φ1 + φ_d1)` と相対誤差 10⁻⁶ 以内で一致。f_sh の成分が理想モデルで 0（数値誤差のレベル）。

V2. **基底帯域と RF レートの等価性**: 白色のショット雑音を RF レートで加えて復調した出力の PSD が `2 S_I G²` と 5% 以内で一致（平均 ≥ 100 回）。基底帯域モードの出力と同じ統計。

V3. **lock 2 による縮小**: Ref A の値（e^{2r} = 123.47）で lock 1 の設定値を δu′ = 0.1 rad ずらすと、lock 2 の零点の Δ が `atan(tan(0.1)/123.47) = 0.813 mrad` だけ動く（相対 10⁻³ 以内）。Ref B（G = 18.7、x = 0.769）では δu = 0.1 rad で 13.1 mrad。
  閉ループで lock 1 のセンサ雑音だけを入れ、Δ の残差が lock 1 の残差の A−/A+ 倍になること（lock 2 の帯域内で 10% 以内）。

V4. **Ref B lock 1**: e1 ∝ −2 ΔG J1(2m) sin Φ1。零点が Φ1 = 0、π で、1 つの符号でどちらか一方だけが安定。m を変えたとき傾きが J1(2m) に比例（m = 0.05–0.5 rad で 10⁻⁶ 以内）。

V5. **量子雑音ロックの誤差信号**: 平坦なスペクトルで e(Δ0) = 2k ΔP J1(2m) sin 2Δ0（2 乗検出）。零点が 0 と π/2。ログ検出で、2 つのロック点の傾きの比が P_+/P_−。

V6. **radiometer 雑音**: 8–300 MHz の帯域制限ガウス雑音（f_s ≥ 1 GHz）から 2 乗検出・平均時間 T で作った P̂ の相対分散が 1/(BT) と 5% 以内で一致（BT ≥ 100、ブロック ≥ 1000）。線形包絡線で倍率 1.02 ± 0.1、ログで 1.32 ± 0.15 [MC]。
  時間ループ段の等価ブロック（`S_P̂ = 2∫S_x² df` のガウス雑音）が、同じ設定で RF レートの計算と同じ誤差信号 PSD を 10% 以内で出す。

V7. **センサ雑音だけで決まるロック残差**: Ref C の値（§6.3）、他の雑音なし、1 次のループでノイズ帯域 b_loop のとき、Δ の rms が `√(S_Δ,sens b_loop)`（b = 10 kHz で 3.2 mrad、2 乗検出）と 10% 以内で一致。反 squeezing 側では 78 mrad。

V8. **dither の寄与**（prompt §7-10）: m = 0.0785 rad の dither だけで、LO 位相の Δ の rms が 0.0555 rad。

V9. **tap の trade-off**（prompt §7-8）: Ref A の値で T_tap = 0.1、1、3% のとき、θ = 9 mrad 固定で R′− = −10.11、−9.76、−9.07 dB（±0.01 dB）。phase-detection OPA の構成では T_tap による損失が 0。

V10. **Liu 型の打ち消しと統合の符号**: ポート k の区間だけに遅い位相雑音（ψ_f = ψ_b、例えば 1 rad の正弦、0.1 Hz）を入れる。
  (A) 構成では、指標 (ii) の残差がループの抑圧で決まる値（開ループ利得の逆数倍）以下。
  (B) 構成では、ポート間の差では ψ_k が消える一方、ポンプに入れた位相雑音 δφ_p が θ_sq,k − φ_pr,0 に係数 1（(A) で lock 1 を切ったときの δφ_p/2 ではなく δφ_p）で現れ、各ポートのアクチュエータ出力が δφ_p/2 を追う。
  ψ_f と ψ_b を別にした（直交偏光の差、T10）ときの残差が ψ_f − ψ_b に一致。

V11. **DC fringe**: v(φ) = Gℜη P[ε + V cos φ] + V_off が厳密に一致。K_d が φ_set = π/2 で −GℜηPV。V_off を加えるとロック点が V_off/K_d だけ動く。

V12. **オフセット**: 量子雑音ロックで a_RAM = 10⁻² を入れると、ロック点が `a_RAM R_−/(2(R_+−R_−)J1(2m))` = 2.8 mrad（Ref C の値）だけ動く（10% 以内）。ヘテロダインで φ_d1 を δ ずらすと Φ1 のロック点が δ だけ動く。

V13. **co-propagation**（prompt §7-11）: Ref B の構成で、probe と squeezed 光が共有する区間にだけ雑音を入れると、lock 2 の誤差にも Δ にも現れない。

---

## 12. シミュレータが必要とするもの

### 12.1 ブロック（段 T）

各ブロックは `make_X(params) -> step(state, params) -> state`。入力はビームの位相（相互作用点で評価、conventions §9）とパワー、出力は電圧 `state.v[name]` と、混入光の情報 `state.leak[name]`（T10 用）。

| ブロック | 入力 | 出力 | 主なパラメータ |
|---|---|---|---|
| `parametric_beat_detector`（1a・1b、汎用: 任意の OPA/OPO の出力の PD） | φ_p、φ_pr（相互作用点）、P_in、OPA の A±（O の表から）、φ_OPA | PD 電圧（DC ＋ 2f_sh の I/Q、または dither の 1 次） | f_sh、T_tap、ℜ、G_TIA、i_n、η |
| `bhd_probe_beat`（1c・3） | φ_LO、φ_pr、φ_p（squeezer）、P_LO、P_pr、A± | BHD の f_sh（または f_m）の I/Q | f_sh または m・f_m、BHD の R(f) |
| `quantum_noise_power`（4） | Δ（φ_LO と θ_sq）、O の表 P(Δ) と ∫S_x² の表 | 包絡線・ログ検出器の出力（雑音つき） | f1、f2、κ_det、k_det または k_dB、ビデオ帯域 |
| `dc_fringe`（5、Liu 型の homodyne_dc） | φ_a、φ_b、P_tot、V、ε | 差電圧 | G_TIA、ℜ、V_off、i_n |
| `liu_probe_interferometer`（6） | φ_pr,k（逆向き、軸 a⊥）、φ_ref,k、P_pr,det、P_ref、ゲート | 電圧（DC／ヘテロダイン／dither） | mode、δf、D、ゲートの有無 |
| `demodulator`（共通） | 電圧、参照位相 φ_d（RF 源の状態から） | I、Q（LPF 後） | f_ref（f_sh、2f_sh、f_m）、LPF、φ_d、V_off |
| `dither_modulator`（共通、§5.5 の Modulators と同じ族） | — | 指定した光路の位相に m sin ω_m t を加える（挿入損失は squeezed 光の経路なら損失として） | m、f_m、a_RAM、φ_a |
| `rf_reference` | — | AOM・復調の RF 位相 φ_RF,j(t)（遅延つき） | S_φRF、τ_RF |

### 12.2 オフライン段の表（段 O）

- OPA/OPO の振幅利得 A±（ポンプの関数）: T01・T02 から。
- 量子雑音ロック用の P(Δ) の係数（P̄、ΔP）と `∫S_x(f;Δ)² df`（Δ の関数、または R± の係数）: T02（R±(f)）、T12（損失）、T13（H、S_el）から、帯域 [f1, f2] で積分。
- 周波数で回る squeezing 角 δ(f) がある場合の重み付きロック点 Δ*（§6.1）。

### 12.3 トグル

`shot_noise`、`electronic_noise`、`radiometer_noise`、`ram_offset`、`demod_offset`、`rf_reference_phase_noise`、`probe_leak_output`（T10 へ）、`dither_on_path`、`f_sh_leak`（1a の理想外の f_sh 成分）、`log_detector`（2 乗／線形／ログ）、`gating`（Liu 型）、`rf_rate_mode`（検証用）。

### 12.4 時間刻み

prompt §5.4 (a) の複素包絡線を想定。dither の f_m・f_sh の成分は解析的に復調し、雑音は復調後の PSD で加える。時間刻みはループ遅延と制御器の更新（T16）で決め、f_m より細かくする必要はない。RF レートのモードは検証（V1、V2、V6）用。

---

## 13. DECISIONS・conventions との関係（調整役への提案）

- 食い違いではないが追記の提案: conventions §5 に「周波数シフトした probe の ω0 系での表し方 `|A| e^{i(φ_pr − Ω_sh t)}` と、AOM の RF 位相の符号（E ∝ e^{−iωt} で up-shift は e^{−i(Ωt + φ_RF)}）」を加えると、T14 と T16 の実装で符号を揃えやすい。
- conventions §4 の θ_amp = φ_p/2 は、位相不整合のある実際の OPA では (φ_p + φ_OPA)/2 になる。φ_OPA を OPA ごとのパラメータとして State に持たせる（Ref A の非共通項）ことを提案する。
- survey E73 の「S_φ ≈ hν/(η P)（係数は要確認）」: 強い参照光とのホモダインでは `hν/(2ηP_pr)`、2 つの同程度の強さのビーム（Ref D）では `2e/(ℜ P V²)`、ヘテロダイン（I だけ）ではさらに 2 倍。方式ごとに §9 の式を使う。
- survey §1 の McKenzie の「(1/Δω)^{1/4}」は静的な見積もり。閉ループのセンサ雑音には §6.3 の (b_loop/B)^{1/2} を使う。

## 14. ユーザの決定が必要な点

1. **統合案（Liu 型 ＋ phase-detection OPA）の目的**: §8.3 の (A) は「共通参照光とポンプの関係」を固定するもので、指標 (ii)（ポート間）には直接効かず、指標 (i)（LO_k に対する位相）や外部 LO を使う測定に効く。どちらの指標のために統合を考えているか。(B) の形（逆向き probe をポンプと直接比べる）を想定しているなら、ポンプの位相の係数とアクチュエータのレンジの問題がある。
2. **Liu 型の検出方式**: DC ホモダイン（Ref D と同じ DSC1 で組める）、ヘテロダイン（probe の周波数をずらすと、Liu の「周波数が同じ」とは違う構成になる）、dither のどれを第一候補にするか。D28（probe 仕様）と合わせて決める必要がある。決まるまでは 3 つとも設定で選べるようにする。
3. **検出点の位置**: ポートごとに振り分け素子の手前で検出する（D30 の記述どおり、§8.4 の制約なし）か、振り分け素子を逆向きに通してから検出する（ポートあたり f_rep/32 の標本化）か。

---

## 15. 文献

一次資料（本文を読んで確認したもの）:
- [Ref A] K. Hirota, T. Kashiwazaki, G. Ha, T. Yamashima, P. Jaturaphagorn, T. Suzuki, K. Takahashi, A. Kawasaki, A. Inoue, W. Asavanant, M. Endo, T. Umeki, A. Furusawa, "Generation of 10-dB squeezed light from a broadband waveguide optical parametric amplifier with improved phase locking method," arXiv:2511.15082 (2025), DOI 10.48550/arXiv.2511.15082（arXiv の記録で確認）。**掲載誌は UNVERIFIED**（検索で見つからなかった）。
- [Ref B] Y. Takeno, M. Yukawa, H. Yonezawa, A. Furusawa, "Observation of −9 dB quadrature squeezing with improvement of phase stability in homodyne measurement," Opt. Express 15, 4321–4327 (2007), DOI 10.1364/OE.15.004321（arXiv quant-ph/0702139 の記録で確認）。
- [Ref C] S. Verclas, C. Darsow-Fromm, B. Tohermes, R. Schnabel, "Fiber distribution of phase-stabilized GHz-bandwidth squeezed vacuum states of light between two buildings," Optica Quantum 4, 1 (2026), DOI 10.1364/OPTICAQ.567418（PDF に記載）。
- K. McKenzie, E. E. Mikhailov, K. Goda, P. K. Lam, N. Grosse, M. B. Gray, N. Mavalvala, D. E. McClelland, "Quantum noise locking," J. Opt. B: Quantum Semiclass. Opt. 7, S421–S428 (2005), DOI 10.1088/1464-4266/7/10/032（PDF に記載）。
- H.-L. Liu et al., "Gaussian boson sampling with 1,024 squeezed states in 8,176 modes," Nature 653, 687–692 (2026), DOI 10.1038/s41586-026-10523-6（PDF に記載）。
- H. Vahlbruch, S. Chelkowski, B. Hage, A. Franzen, K. Danzmann, R. Schnabel, "Coherent control of vacuum squeezing in the gravitational-wave detection band," Phys. Rev. Lett. 97, 011101 (2006), DOI 10.1103/PhysRevLett.97.011101（APS の記録と arXiv:0707.0164 の本文で確認。2 倍のビート周波数での復調の記述を本文で確認）。
- R. H. Dicke, "The measurement of thermal radiation at microwave frequencies," Rev. Sci. Instrum. 17, 268–275 (1946), DOI 10.1063/1.1770483（書誌情報を確認。radiometer の式の出典。本文は未読、式は §6.3 で独立に導出）。

メーカー資料:
- Thorlabs PDB440x/PDB450x manual（© 2019、`refs/manuals/thorlabs-pdb.pdf`）: PDB450C の responsivity 1.0 A/W（typ. max）、利得 10³–10⁷ V/A、帯域、CMRR、飽和 4.5 mW。
- Thorlabs DSC1 manual（`refs/manuals/thorlabs-dsc1.pdf`）: 入力雑音 0.3 mV rms（値は survey §1 から引用、T15 で確認）。

このノートで使った式の出典（記憶による、**UNVERIFIED**）:
- 空気の屈折率の分散式: E. R. Peck, K. Reeder, J. Opt. Soc. Am. 62, 958 (1972)（係数は記憶から。§3.5 の概算だけに使用）。
- 溶融石英の Sellmeier 式: I. H. Malitson, J. Opt. Soc. Am. 55, 1205 (1965)（係数は記憶から。§3.5 の概算だけに使用）。
- Ref A の文献 [23]: T. Kashiwazaki et al., Appl. Phys. Lett. 122, 234003 (2023)、DOI 10.1063/5.0144385（arXiv:2301.12658 の記録にある DOI。本文は未読）。前作の AOM 方式の詳細はこの文献にある可能性がある。

計算の記録: 概算は `est.py`、radiometer の数値実験は `radiometer_mc.py`（調査担当の scratchpad、リポジトリには入れていない）。必要なら Phase 4 の検証テスト V6 にそのまま移せる。
