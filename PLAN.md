# PLAN（Phase 3 草案、承認ゲート前）

- Phase: 3（計画）／ 日付: 2026-10-31 ／ 状態: **草案**。調整役がユーザと確認（grill-with-docs）し、§3 の決定を `DECISIONS.md` に入れてから改訂する。承認まで実装しない（prompt §4 Phase 3）。
- 入力: `docs/squeezing_sim_prompt.md`（§4–§8）、`AGENTS.md`（保護リスト）、`DECISIONS.md`（D01–D51、最上位は **D34**）、`notes/conventions.md`、`notes/00_summary.md`（以下「S」と略し、S §3.2 のように節を書く）、deep dive T01–T17。
- 規約: 本文は日本語、コードの識別子・記号は英語（D19）。符号・規格化・PSD は `notes/conventions.md`（以下「conv」）に従う。
- この計画の約束:
  1. **D34**: どの選択肢にも「既定」を置かない。§3 の推奨は「承認されれば参照プリセットで使う値」であり、他の選択肢を消す意味ではない。悪い構成もプリセットとして残す（S §2.5）。
  2. **保護リスト（AGENTS.md）**: 実験ランナー、`State`、光路ごとの位相モデル、遅延バッファ、ブロックの factory インターフェース、入れ子のパイプライン（並列・積み重ね）、YAML とその Python リストとの等価性、機器に依らないパラメータ・機器プリセット・抽象制御器インターフェース、Ref A/B/C/D のためのブロック、GNLSE の項と雑音源のトグル、config・seed・commit の記録、検証テスト。この計画のどの簡略化もこれらには及ばない。
  3. `blind/` は読まない（D04）。Ref D の予測はメーカー資料と物理だけから作る。

---

## 1. 目標と範囲

### 1.1 Phase 4 全体の目標

prompt §0 の 3 つの目的（squeezing の時系列・スペクトルの予測、ファイバ長ごとの必要な位相ロック性能、構成を切り替えて試せるテストベッド）を満たすシミュレータを、検証テスト（§6）付きで作る。Phase 6 の探索（アブレーション、BoTorch）が設定の変更だけで回せる状態を Phase 4–5 の完了条件とする。

### 1.2 マイルストーン（D17: 時間ループ＋単一モード → GNLSE、D03: A → C → B、D は A と並行）

| M | 内容 | 合格条件（§6 の VT） | 根拠 |
|---|---|---|---|
| M0 | フレームワーク: `State`・`Topology`、ブロックの規約、パイプラインの合成、遅延バッファ、多重レート、コンパイル層と純 Python の参照実装、ランナー、設定（YAML・プリセット・オーバーレイ・callable）、乱数、記録、損失の登録表、帯域の台帳、`SourceOutput`、指標の骨組み、抽象制御器インターフェース | VT-05（純損失の部分）、VT-13（汎用の線形ループ）、VT-22（共有区間、合成した位相で）、VT-44（遅延バッファ、積み重ね = 並列、コンパイル = 純 Python）、YAML ↔ Python の往復 | prompt §5.2–5.3、保護リスト |
| M1 | **Ref A**（単一モード、`refA`・`refA_tap`）と **Ref D**（`refD`、解析的な検証のみ、D04）を並行 | VT-01、02、04、15、16、24、13・14（Ref D）、21（Ref D の部分）、30、33、34、42、43（A・D の部分） | D03、D04 |
| M2 | **Ref C**（CW OPO の R±(f)、量子雑音ロック、3 段の積み重ね、相関のある 2 本のファイバ、ランプと再引き込み） | VT-19、20、21、31、37、41（C の部分） | D03 |
| M3 | **Ref B**（L(x) の callable、dither 極値ロック、同軸 probe、OPO 区間の係数差） | VT-17、18、22（B）、23 | D03、D48 |
| M4 | **将来の構成の時間ループ**（単一モードの光源で代用）: Liu 型ロック、逆向き伝搬の小区間、Jones/PM/SMF、TEC、振り分け、ポートごとの損失、基準ポート 0 ＋ ポート k、2 ポート組、指標 (ii) の 4 成分 | VT-38、39、40、41、33（V11・V12） | D08–D13、D30、D31、D47 |
| M5 | **オフラインの重い部分**: GNLSE（RK4IP、項のトグル、CW の近道）、ファイバの量子雑音層、散乱の表、パルス導波路 OPA の超モード、パルス励起 OPO、パルス BHD、モードの重なりの表 | VT-03、06、07–12、32、35、36、45 | D17、prompt §5.1 |
| M6 | **将来のプリセットの統合**（`F_wg100`・`F_wg6`・`F_bor6`・`F_liu20`、S §2.3 の軸、S §2.5 の悪い構成）、性能の確認（32 ポートの規模の見積もり） | 各プリセットが走り、CHK の警告が S §2.2 の「チェック」行どおり出る | D13、D34 |
| （Phase 5 の終わり） | Ref D のブラインド比較（§6.5） | 合わせ込みなし | D04 |

- §5.1 の要件は M0 から設計に入れる（D17）: `SourceOutput`（G07）、GNLSE の項のフラグの設定の形（`fiber_gnlse.terms.*`）、`LASER_MODE = CW` で GNLSE を飛ばす経路（M1 では `fiber_linear_table` だけで CW を通す。M5 で VT-03 により RK4IP と一致を確認）。
- 32 ポートを同時に走らせる設定は Phase 6 で必要になったときに作る（D13）。M6 では規模（メモリ・時間）を見積もるだけ。ただしデータ構造はポート数に上限を置かない。

### 1.3 範囲外（Phase 4 ではやらない）

| 項目 | 時期 | 理由 |
|---|---|---|
| BoTorch の最適化ループ | Phase 6 | prompt §8。Phase 4 では抽象制御器インターフェースと `set_params`・`acquire` を用意するだけ |
| 実機のラッパ（`MokuController`、`DSC1Controller`） | 実機を使うとき | 実機の API が手元にない（D15）。インターフェースと「正準係数 → 機種」の写像（G16-h）は Phase 4 で作る |
| 32 ポートの同時計算 | Phase 6 | D13 |
| GBS の出力分布の忠実度 | Phase 6 で必要なら | D46 |

---

## 2. アーキテクチャ（prompt §5）

### 2.1 全体像

```
 config (YAML or Python list)             presets/ + overlays + devices/
            │  load + overlay merge + resolve callables (topological order)
            ▼
     RunSpec (pure data tree: run / topology / offline / pipeline / toggles / measure)
            │
   ┌────────┴─────────────────────────────────────────────────────────────┐
   │ OFFLINE STAGE (per configuration, cached by hash)                    │
   │  source models → SourceOutput (covariance + mode basis + projector)  │
   │  fiber: GNLSE(RK4IP, term flags) / CW shortcut → linear tables       │
   │  quantum layer (linearized covariance; TW optional)                  │
   │  LossRegistry → η per path & per metric (GBS input / BHD)            │
   │  BandLedger  → out-of-loop-band variance per (metric, path)          │
   │  tables: V±(f), P(Δ), ∫S_x², A±(P), mode overlap, stray light, ...   │
   └────────┬─────────────────────────────────────────────────────────────┘
            │ Tables (read-only arrays)
   ┌────────▼─────────────────────────────────────────────────────────────┐
   │ TIME-LOOP STAGE                                                      │
   │  build: Topology (segments, beams, locks, ports) + SlotLayout        │
   │         blocks = [make_xxx(params) ...] nested Seq/Parallel/Stack    │
   │  outer loop over chunks (Python): open-loop noise for next chunk,    │
   │         slow multirate parts, decimated recording                    │
   │  inner loop over samples (compiled): step(state, params) of each     │
   │         block in order; feedback only through delay buffers          │
   └────────┬─────────────────────────────────────────────────────────────┘
            │ recorded series: true path phases, errors, u, δψ, lock events
   ┌────────▼─────────────────────────────────────────────────────────────┐
   │ MEASUREMENT / METRICS                                                │
   │  metric (i) δ_i,k(t), metric (ii) Θ_k(t) + 4 components, overlap,    │
   │  stray light; phase_to_quadrature (+ BandLedger once); ESA / pulse   │
   └────────┬─────────────────────────────────────────────────────────────┘
            ▼
   runs/<run_id>/ : config (input + resolved), seeds, git commit, env, results, table hashes
```

### 2.2 オフライン段（prompt §5.1）

| 項目 | 設計 |
|---|---|
| 入力 | `RunSpec.offline`（光源、ファイバ区間、損失要素、検出器、帯域の宣言）。時間ループの状態には依らない |
| 出力 | `Tables`: 名前付きの読み取り専用の配列と、補間の格子・単位・出典の印・生成した設定のハッシュ |
| GNLSE | 自前の RK4IP（S §6、T06 §4。gnlse-python・PyNLO は NumPy 2 で動かず、照合用にだけ使う）。項はフラグ `loss`・`beta2`・`beta3`・`kerr`・`raman`・`self_steepening`・`xpm_external`・`cw_shortcut`（保護リストのトグル）。順に並べるパイプラインの段ではなく**方程式の項**として切り替える（prompt §5.1）。FFT は conv §8 のラッパ（`core/fft.py`）だけを通す |
| CW | `LASER_MODE = CW` では GNLSE を解かず、線形の表（損失、静的な位相、SPM の定数位相）を作る（`cw_shortcut`）。M5 で RK4IP の CW 解と 10⁻¹² で一致（VT-03） |
| 量子雑音層 | G03 の推奨では線形化（古典解の Jacobian ＋ 共分散、決定論的）。TW は任意の検証ブロック |
| キャッシュ | `runs/_cache/<hash>.npz` ＋ JSON のメタ（G16-g）。ハッシュは「その表に効く設定の部分木 ＋ コードの版」 |
| 光源 | すべて `SourceOutput`（G07）を返す。時間ループと測定は光源の種類を知らない |
| 損失 | すべてのブロックが自分の損失要素を `LossRegistry` に登録する（G06） |
| 帯域 | すべての雑音源が「時間ループで表す帯域」を `BandLedger` に宣言する（G05） |

### 2.3 State（conv §9–§10、D11、D47、D51）

`State` は 2 層に分ける。**静的な `Topology`**（ビルド時に決まり、ループ中は変わらない）と、**動的な `State`**（毎サンプル更新）。実装の表現は G02。

#### 2.3.1 Topology（静的）

| 要素 | 中身 | 出典 |
|---|---|---|
| `segments[i]` | id、種類（`pm`・`smf`・`free_space`・`opo`・`crystal`・`demux`・`component`）、L、n_g、τ_g0、β2、β3、熱ゾーン、`N_sub`、空間構造（`lumped`・`distributed`・`mixed`）、偏光モデル（G04）、アクチュエータの位置 x_a（小区間の番号） | conv §9、D47、T08、T10 |
| `beams[b]` | 光源（どのレーザー・ポンプ・光源から出たか）、**経路 = (区間 id, 向き d, 入る小区間, 出る小区間) の列**、周波数オフセット、偏光の入射状態（Jones ベクトル）、名目パワー | conv §9、D11 |
| `coef[i][driver][beam]` | 区間 × 要因（`dT`・`strain`・`detuning`・`tec`・`stretcher`）× ビームの係数。定数・callable・表のどれでもよい（G16-a）。位相と群遅延で別の係数（C9 を定数にせず計算） | D47、T02、T10、T17 |
| `locks[l]` | ビームの位相の線形結合 Σ c_b φ_b と設定点、評価する場所（相互作用点・検出器）。例 Ref A lock 2: φ_LO − φ_pr | conv §9 |
| `ports[k]` | ポート k の squeezed ビーム、LO_k ビーム、probe ビーム、GBS 入力の観測点、BHD の観測点。**基準ポート 0 は常に持つ**（G13） | D12、D31 |
| `observers` | 位相を評価する観測点（検出器、GBS 入力、相互作用点）と、その点までの経路上の累積の群遅延 | conv §7 |

#### 2.3.2 State（動的）— フィールドの最終案（S §3.1 を確定する案）

| フィールド | 中身 | 書くブロック | 読むブロック |
|---|---|---|---|
| `t`、`n` | 時刻、サンプル番号 | ランナー | 全部 |
| `laser[s]` | δν、φ_L = −2π∫δν、δp、**履歴のリング**（放出時刻で評価するため、最長の経路の群遅延 ＋ ポート間の遅延 τ_k を覆う長さ） | `laser_source` | `beam_assembler` |
| `drv[i]` | 区間 i の要因（ΔT、ε、Δ_s など）を**小区間ごと**に。履歴のリング付き（逆向きの光が別の時刻に通るため） | NOISE、TEC、ACT | `segment_phase`、`segment_jones` |
| `psi[i][beam]` | 区間 × ビームの位相 ψ（小区間の和を、そのビームが各小区間を通る時刻で評価した値） | `segment_phase` | `beam_assembler` |
| `tau_g[i]` | 区間の群遅延 τ_g(t) | `segment_phase` | `beam_assembler`、`timing` |
| `J[i]` | 区間の 2×2 Jones 行列（複素 4 成分。PM の対角近似では 2 つの位相だけ）（G04） | `segment_jones` | `beam_assembler` |
| `beam[b]` | 観測点での位相 φ_b、パワー、Jones ベクトル、放出時刻 | `beam_assembler` | ERR、DET、MEAS |
| `pump`、`opa[j]`、`opo[j]` | φ_p、P、x_ss、φ_OPA,0、Δ_s、Δ_p | SRC | ERR、MEAS |
| `v[name]`、`dem[name]` | 電圧、復調後の I/Q | ERR、DET | CTRL |
| `ctrl[name]` | 出力 u、積分値、飽和フラグ | CTRL | `delay`、ACT |
| `delay[name]` | 遅延バッファ（整数サンプルのリング ＋ 端数） | `delay_buffer` | ACT |
| `act[name]` | 機械モードの状態、δψ、ΔT_h、飽和 | ACT、TEC | `segment_phase`（`drv` 経由） |
| `timing[k]` | Δt_k(t)、Δt_kj(t) | `timing_state` | `mode_overlap_lookup` |
| `leak[name]` | probe のパワー・タイミング・偏光・周波数・ε_AOM（T10 の迷光の表の入力と同じ名前、C15） | MOD、ERR | MEAS（迷光） |
| `lock[name]` | ロック状態、滑り・リセットのイベント | `lock_supervisor` | 記録 |
| `metric` | δ_i,k、Θ_k、重なり、迷光、η_k（GBS 入力・BHD） | MEAS | 記録 |

**放出時刻の評価**（conv §7、§9）: 観測点に時刻 t に着くビーム b の位相は
φ_b(t) = φ_L(t − τ_b) + Σ_{(i, m) ∈ 経路 b} ψ_{i,m,b}(t − τ_{b,after(i,m)}) + （変調器・アクチュエータの δψ も区間の ψ に含む）。
ここで τ_b は光源から観測点までの群遅延、τ_{b,after(i,m)} は小区間 (i, m) の出口から観測点までの群遅延。逆向き（d = −1）のビームは同じ小区間を別の時刻で読む。これを `beam_assembler` が `laser.hist` と `drv` の履歴から端数の補間で読む。周波数による伝搬位相の変化として二重に数えない（conv §7）。
共有区間で相殺が起きるのは、この和の結果としてだけ（「仮定」ではなく「計算」、D11）。係数が違えば相殺しない（VT-23）。

**指標**（conv §10、D12、D31、D46、D47）:
- (i) δ_i,k(t) = wrap_π(φ_LO,k − θ_sq,k − Δ_set)（BHD の観測点で）。
- (ii) Θ_k(t) = wrap_π(θ_sq,k − θ_sq,0 − Θ_k,set)（GBS 入力で、同じ GBS ショットに着くパルスどうし）。4 成分（ロックされる光路の残差、共通点、レーザーの相関の崩れのうち帯域外、squeezed 光だけの区間）は G13 の方法で分ける。帯域内のレーザー項 16 sin²(πf(T_a − T_b)) sin²(πf(τ_k − τ_0)) S_φ と T_a − T_b はパラメータ（D51-5）。
- ポート間の波束の重なり、迷光（ゲートあたりのクリック数 ＋ モード内の光子数）、η_k（GBS 入力と BHD を分ける）。

### 2.4 ブロックの規約（factory インターフェース、保護リスト）

```python
# core/block.py (sketch; names are the contract, bodies are Phase 4 work)
class Block(Protocol):
    kind: str                    # registry key, e.g. "pid"
    name: str                    # unique instance name, e.g. "lock1.pid"
    spec: BlockSpec              # (kind, name, params) -> exact YAML round-trip
    reads: tuple[str, ...]       # State slot names read
    writes: tuple[str, ...]      # State slot names written
    rate_div: int                # step every rate_div samples (G10)
    def __call__(self, state: State, params: Params) -> State: ...   # reference step (pure Python)
    kernel: KernelFn | None      # compiled step on flat arrays (G01/G02); None = Python only
    def init(self, layout: SlotLayout, rng: np.random.Generator) -> None: ...
    def losses(self) -> list[LossElement]: ...   # G06 (empty if none)
    def bands(self) -> list[BandDecl]: ...       # G05 (empty if none)

def make_pid(params: dict) -> Block: ...         # every block: make_<kind>(params) -> Block
```

- factory は `make_<kind>(params) -> Block` で、`Block` は呼び出し可能（`step(state, params) -> state`）。prompt §5.2 の「factory が step 関数を返す」をそのまま満たし、メタ情報（読む・書くスロット、損失、帯域）を同じ物に持たせる。
- `params` は 2 種類: **構造パラメータ**（ビルド時に固定。区間の数、フィルタの次数など）と**実行時パラメータ**（`Params` の配列。PID の利得など。`set_params` で走行中に変えられる、BoTorch 用）。
- 各ブロックは `@register_block("pid")` で登録し、YAML の `block: pid` と Python の `make_pid` が同じ物を作る。
- 機器に依存する値はブロックに書かない。`device: dsc1` の参照と、機器プリセットの YAML（§2.9）で与える（prompt §5.3）。
- ブロックの内部の符号の約束は conv に揃える。機器が別の符号を使うとき（DSC1 の e = S − v）はブロックの中で変換する（conv §9）。

### 2.5 パイプラインの合成（入れ子、保護リスト）

| 合成 | 意味 | 例 |
|---|---|---|
| `Seq([a, b, c])` | 順に実行。b は a が書いたスロットを同じサンプルで読める | `adc >> pid >> dac` |
| `Parallel([p, q])` | 同じ入力を読む分岐。**書き込みのスロットが互いに素**であることをビルド時に検査し、宣言順に実行（G09） | `[hpf >> pid1 >> eom] ‖ [lpf >> pid2 >> stretcher]` |
| `Stack([c1, c2, c3], actuators=[a1, a2, a3])` | 積み重ね（Ref C）: c_k の出力が a_k を駆動し、同時に c_{k+1} の入力になる。`Seq`・`Parallel` に展開した糖衣（T16 §4.4: `pid1 >> (eom ‖ (pid2 >> (mirror ‖ (pid3 >> stretcher))))`） | Ref C の EOM → PZT ミラー → ストレッチャ |
| 名前による配線 | ブロックのパラメータ `in:`・`out:` でスロット名（`v.bhd`、`ctrl.lock1.u`）を指定 | — |
| ビルド時の検査 | (1) 読むスロットに書き手がいる、(2) 同じサンプル内の循環（代数ループ）がない — **フィードバックは `delay_buffer` を通る経路だけ許す**、(3) 並列の分岐の書き込みが素、(4) 単位のタグが合う（V と rad など） | prompt §5.2 |

- 入れ子は任意の深さ。合成した物もまた `Block`（`reads`・`writes` は和集合）なので、制御器のまとまりを 1 つのブロックとして YAML の別ファイルに置ける。
- 実行順の既定の並び（推奨、G09）: 光源・雑音 → 区間の位相 → ビームの組み立て → 検出・誤差信号 → 制御器 → 遅延バッファ → アクチュエータ（次のサンプルの区間の要因に効く）→ 指標のプローブ。パイプラインはこの順を強制しないが、ビルド時の検査で「読む前に書かれていない」を見つける。

### 2.6 遅延バッファとフィードバック（保護リスト）

- 制御器は `ctrl[name].u` に書き、アクチュエータは `delay[name]` を通して読む。遅延 τ は「整数サンプルのリング ＋ 端数（線形補間 / 1 次 Thiran の全域通過）」（T05 §3.5、T16 §5）。
- ループ遅延の内訳を足し算で持つ: 制御器の純遅延（機器プリセット）＋ ZOH の h/2（離散化で厳密に扱う場合は 0）＋ 光の伝搬（アクチュエータから検出器までの n_g L/c。Liu 型では逆向き probe の飛行時間、E82）＋ 任意の追加。
- 光の伝搬遅延（放出時刻、逆向き伝搬）は `beam_assembler` が `laser.hist`・`drv` の履歴で扱う。これも遅延バッファの同じ実装（`core/delay.py`）を使う。

### 2.7 多重レートと時間刻み（G10）

- `f_sim` は設定の `run.f_sim_hz`。目安（T16 §5-4）: f_sim ≥ max(制御器の更新レート、モデル化した最高の機械共振 × 20、f_u × 10)。参照プリセットの値は G10 の表。
- ブロックは `rate_div` を持ち、n % rate_div == 0 のサンプルだけ更新し、間は値を保持する。熱（10²–10³ s）は 10–100 Hz、TEC・遅い PID は間引いたレート（T08 §7.3）。
- MHz 以上の RF 搬送波（ESA の 3 MHz、Ref C の 8–300 MHz、Ref A の 1 MHz）は基底帯域・統計モデル（G01）。RF レートで直接作るモード（`rf_rate_mode`）は検証用（C6、T14 V1・V2・V6）。

### 2.8 実験ランナー（保護リスト）

| 関数 | 中身 |
|---|---|
| `run(spec) -> RunResult` | オフライン段 → ビルド → 時間ループ（チャンク） → 測定 → 記録 |
| `sweep(spec, axes) -> list[RunResult]` | パラメータの格子・リスト（D16 の「目標 × 長さ → θ_req」、S §7.2） |
| `ablate(spec, toggles)` | トグルを 1 つずつ切る（Phase 6 の最初のアブレーション、prompt §8） |
| `decompose_metric_ii(spec)` | 雑音源の組を切り替えた 4 回の再実行と線形性の検査（G13） |
| `convergence(spec, what)` | f_sim の倍化・チャンク長・seed の集団（§6.4） |
| `SimController` | 抽象制御器インターフェースのシミュレータ実装（§2.10） |

ランナーは純粋な関数で、結果は `RunResult`（時系列は間引いて記録、PSD、指標、イベント、表のハッシュ）。

### 2.9 設定: YAML ↔ Python リスト、プリセット、オーバーレイ、callable（D34、保護リスト）

- **同じ物を作る**: YAML も Python も `RunSpec`（純粋なデータの木）になり、`build(RunSpec)` が同じブロックを作る。Python で組んだパイプラインは `to_spec()` で `RunSpec` に戻り、`dump_yaml` で YAML になる。往復の等価性をテストにする（`spec == load(dump(spec))`、同じ seed で結果がビット単位で同じ）。
- **プリセット**: S §2 の名前をそのまま使う。参照 `presets/reference/{refA, refA_tap, refB, refC, refD}.yaml`、将来 `presets/future/{F_wg100, F_wg6, F_bor6, F_liu20}.yaml`、比較用の悪い構成 `presets/bad/<name>.yaml`、ロックのトポロジー `presets/locks/<name>.yaml`（S §2.4）。**どれも同格**で「既定」のファイルはない。
  - ディレクトリ名に `refs` を使わない（`.gitignore` の `refs` が全階層に効くため）。
- **オーバーレイ（軸）**: `presets/axes/<axis>/<value>.yaml`（S §2.3 の軸: 光源、振り分け、経路の組み方、ファイバ、probe、ロック、大レンジ、LO、BHD、OPO のロック光、ポンプ、制御器、ポート）。合成の規則は G12。
- **callable パラメータ**（prompt §5.5、Ref B の L(x)）: Python では本物の関数を渡せる。YAML では `!expr` と `!fn`（G12）。依存関係は有向グラフにしてトポロジカル順に解き、循環はエラー。
- **機器プリセット**: `devices/{dsc1, moku_delta, fpga_generic}.yaml`。値ごとに出典（マニュアルのページ）を `source:` に持つ（T15 §8.3、T16 §4.5）。制御器のブロックは `device: dsc1` と上書きの差分で設定する。

例（Ref D の一部、形の確認用。値は T15・T16・D35・D51-8 から。最終形は WP-I2 が書く）:

```yaml
preset: refD
run: {duration_s: 60, f_sim_hz: 1.0e6, chunk_samples: 65536, seed: 20261031}
topology:
  segments:
    - {id: smf0, kind: smf, L_m: 2.0, pol: {model: smf_diffusion, D_pol: !sweep [1e-8, 1e-2]}}
    - {id: arm1_915b, kind: pm, L_m: 8.4, N_sub: 8, zone: lab}
    - {id: arm2, kind: pm, L_m: !expr "arm1_915b.L_m - dL", where: {dL: !sweep [0.0, 15.0]}}
  beams:
    - {id: pr, laser: dx1, path: [smf0, coupler1, arm1_915b, coupler2]}
    - {id: lo, laser: dx1, path: [smf0, coupler1, arm2, coupler2]}
  locks:
    - {id: lockD, combo: {lo: +1, pr: -1}, setpoint: !pi_frac 0.5, observer: pdb450c}
pipeline:
  - {block: laser_source, name: dx1, params: {preset: cobrite_dx1_SC}}
  - {block: dc_fringe, name: pdb450c, params: {gain_v_per_a: 1.0e3, sign: !sweep [+1, -1]}}
  - {stack: [{block: adc, params: {device: dsc1}},
             {block: pid, params: {device: dsc1, kp: !sweep_log [1e-3, 1e1]}},
             {block: dac, params: {device: dsc1}},
             {block: delay_buffer, params: {tau_s: 5.7e-6}},
             {block: hv_amp_914}, {block: pzt_fiber_stretcher, params: {segment: arm1_915b}}]}
```

### 2.10 機器に依らないパラメータ、機器プリセット、抽象制御器インターフェース（prompt §5.3、T16 §4.6）

| 部分 | Phase 4 で作るもの |
|---|---|
| 汎用パラメータ | サンプルレート、ビット数、レンジ、帯域、純遅延、出力抵抗と負荷、PID の正準係数（K_p、K_i [1/s]、K_d [s]、ω_IS、ω_DS、二重積分）、離散化の方式（`dsc1`・`moku`・`tustin`）、anti-windup |
| 機器プリセット | `dsc1.yaml`（h = 2 µs、τ 5.7 µs（8.3 µs まで）、16/12 bit、e = S − v、`ki_entry: per_second \| per_sample`、出力 ±4 V・220 Ω）、`moku_delta.yaml`（τ ≲ 0.21 µs、交差周波数による指定、PID 出力 ±0.5 V（50 Ω））、`fpga_generic.yaml`（τ_L 0.1–2 µs、12–16 bit、Tustin）。Moku:Pro は作らない（D51-7） |
| `ControllerInterface`（抽象） | `configure(config)`、`set_params(dict)`・`get_params()`、`engage()`・`disengage()`、`reset_integrators()`、`ramp(...)`、`acquire(duration, signals)`、`measure_transfer_function(injection_point, freqs)`、`status()`、`log_record()` |
| `SimController` | 上を時間ループで実装（`acquire` は走らせて記録を返す、`measure_transfer_function` は正弦の注入、VT-13 にも使う） |
| 写像 | 正準係数 → DSC1（K_i の単位の切り替え、符号の反転、h）、→ Moku（dB と交差周波数）。範囲外は警告して丸め、丸めた値を返す |

### 2.11 乱数、記録、出力（prompt §6）

- 乱数: ルートの seed から `numpy.random.SeedSequence` を作り、**雑音源ごとにブロック名から決まる鍵で子の流れを作る**（G11）。雑音源を 1 つ切っても他の流れは変わらない（指標 (ii) の分解の再実行が同じ実現で比べられる）。
- 記録: `runs/<UTC時刻>_<preset>_<spechash8>/` に `config.input.yaml`、`config.resolved.yaml`（オーバーレイと callable を解いた後）、`seeds.json`（ルートと子の鍵）、`git.json`（commit、dirty の有無、dirty なら差分のハッシュと差分そのもの）、`env.json`（Python・numpy・scipy・numba の版）、`tables.json`（使った表のハッシュ）、`results.npz`、`metrics.json`、`loss_registry.json`（出典の印付き、T12 §8.2）、`band_ledger.json`。`runs/` は gitignore する（WP-F0）。
- 出力（prompt §6）: ESA ゼロスパンの時系列（ショット雑音、squeezing、反 squeezing、位相掃引トレース）、スペクトル（絶対・ショット雑音で規格化、回路雑音を引く・引かない）、ロックごとの残差の時系列と PSD、squeezing・反 squeezing vs ポンプパワー、指標 (i)(ii)、ロックのイベント。図は任意の `plot` 追加依存（G14）で、データは必ず npz/JSON に残す。

---

## 3. 承認ゲートで決めること

> **ゲートの結果（2026-10-31、ユーザ回答、DECISIONS D52）:** G01–G17 はすべて推奨どおり採用。ただし G08 は「位相を直接書く」変種を別のコード経路としては持たず、**係数の設定で再現するプリセット**（2 軸の係数を等しく、群遅延への効きを 0）として残す（D34）。G17 は M5 を M2 と並行して始める。G14 は numba・pyyaml を追加、pytest は dev、matplotlib は extra `plot`、botorch は Phase 6。

各項目は独立に問える形にした。「推奨」は参照プリセットで使う案で、他の選択肢も設定として残す（D34）。ただし G01・G02・G09 のようなフレームワークの形は 1 つに決める必要がある（保護リストの物の形を決めるため）。

### G01 時間ループの実行方式（prompt §5.4、D49 の「numba の可否は Phase 3 で」）

| 選択肢 | 中身 | 根拠 |
|---|---|---|
| (a) 基底帯域だけ | RF をすべて複素包絡線で。ループ自体は純 Python の 1 サンプルごと | 遅延は正確。しかし Python で 6×10⁷ ステップ × 10 ブロックは 1 時間以上（T16 §5-2 の AE） |
| (b) チャンクのベクトル化 | チャンク長 ≤ 最小のループ遅延 | DSC1 で 5 サンプル、Moku 系（≲ 0.21 µs）では 1 サンプル未満になり、(c) と同じになる（T16 §5-1） |
| (c) コンパイルした 1 サンプルごとのループ（numba） | 各ブロックの `kernel` と合成した内側のループを numba でコンパイル | 1 ステップ 50–200 ns なら 60 s が 3–12 s（AE、未計測）。遅延はリング ＋ 端数で任意 |
| (c′) 別のコンパイラ（Cython、C 拡張、JAX） | — | Cython・C はビルドの手間、JAX の `lax.scan` は可変の分岐・飽和・イベントに不向き |
| **推奨: (c) ＋ (a) ＋ 多重レート、(b) は開ループだけ** | 閉ループは numba の 1 サンプルごと。MHz 以上の搬送波は基底帯域・統計モデル（T14 の解析的な復調 ＋ 復調後の雑音 PSD、radiometer 雑音）。開ループの雑音の生成・ESA の後処理はチャンクでベクトル化。**純 Python の参照実装を常に持ち、短い時間で一致を検査（VT-44: 10⁻¹²）** | T16 §5、T14 §12.4、S §6 |
| 条件 | WP-F1 の最初に性能の試作（ブロック 10 個、10⁶ ステップ）で AE を計測し、1 ステップ 1 µs を超えたら G02 の表現を見直す | — |

### G02 State の実装表現とコンパイル層の形（G01 と組。後から変えると高くつく、AGENTS.md）

| 選択肢 | 中身 | 長所 / 短所 |
|---|---|---|
| (A) Python のオブジェクト（dataclass・dict） | `state.ctrl["lock1"].u` | 読みやすいが numba で使えない。(c) と両立しない |
| (B) numba の `jitclass`・`structref` | 型付きの構造体 | 動的なフィールドの集合（構成ごとに違う）を表しにくい。コンパイル時間が長い |
| **(C) 推奨: 平らな配列 ＋ スロットの登録表 ＋ 名前付きのビュー** | ビルド時に `SlotLayout` が各スロット名（`ctrl.lock1.u`、`psi.arm1.pr`、リングバッファの領域）に float64 配列 `S` の添字（と complex128 配列 `Z`、int64 配列 `I`）を割り当てる。カーネルは添字を定数として閉じ込めた `@njit` の閉包。Python の参照実装と記録は名前付きのビュー（`state["ctrl.lock1.u"]`、`state.ctrl.lock1.u`）で同じ配列を読む | 構成ごとに違うフィールドの集合をそのまま表せる。合成は「カーネルを順に呼ぶカーネル」を閉包で作る（numba は閉包が捕まえた jit 関数の呼び出しをコンパイルできる）。`step(state, params) -> state` は同じオブジェクトを返す（中身は在位置で更新） |
| 注意 | コンパイル時間（構成ごとに数秒–数十秒）。`numba` のキャッシュは閉包には効かないので、掃引では同じ構造の構成でカーネルを使い回す（実行時パラメータは `P` 配列で渡す） | — |

### G03 量子雑音層の方式（prompt §5.1、T07 §2.5）

| 選択肢 | 根拠 |
|---|---|
| **推奨: 線形化（古典解の Jacobian ＋ 共分散）を参照の方式、truncated Wigner は任意の検証ブロック（`fiber_quantum_tw`、opt）** | squeezed vacuum の経路は φ_NL ~ 10⁻⁷ rad で線形化は事実上厳密。散乱雑音はすべて加法的なガウス（D）。決定論的なので VT-05・VT-06 を機械精度で検証でき、Phase 6 の目的関数に Monte Carlo の揺らぎが乗らない。TW が要るのは明るい LO・probe の φ_NL が大きいとき（Fw の 10 MHz・633 m の領域、N18）だけ |
| TW を参照にする | N_traj ≈ 10⁴ で重い。§7-4 の検証が統計誤差でしか言えない |
| 線形化だけ（TW なし） | N18 の確認手段がなくなる。保護リストの「雑音のトグル」の趣旨に反する |

### G04 偏光の表し方（S §3.4 C1、D47）

| 選択肢 | 根拠 |
|---|---|
| (a) 軸ごとの係数だけ（T08・T10 の k_T,a、κ_{j,a}） | PM ではほぼ足りるが、SMF の相反性・2 arg(a) の誤差を計算できない（T09 §6.1）。D32 で SMF を同格にしたので不足 |
| (b) 区間ごとの 2×2 複素 Jones 行列だけ | 一般的だが PM の普通の場合に余計な計算 |
| **推奨: (b) を State の表現にし、PM では対角（2 つの軸の位相）に縮めた速い経路を区間の種類で選ぶ（D47 のとおり）。要因（ΔT、ε）は区間で共通に作り、`segment_jones` が区間の種類ごとのモデルで J を作る（PM: 対角 ＋ 小さな非対角（PER）、SMF: 複屈折ベクトルの拡散 D_pol ＋ 温度係数）。J は区間ごと（小区間ごとではない）、逆向きは J^T** | C1 の提案、T09 §1.6・§6.1。小区間は位相のタイミング（逆向き伝搬）のためで、偏光の空間分布まで小区間で持つ根拠はノートにない |
| 補足 | 偏光の効果をすべて切ると J = e^{iψ}·1 になり、単一偏光のテスト（VT-22）がそのまま通ること（T09 テスト 8）を VT に入れる |

### G05 ループの帯域外の雑音を足す場所（S §3.4 C2）

| 選択肢 | 根拠 |
|---|---|
| 各トピックの流儀のまま（T05 σ²_oob、T04 χ_HF、T07 GAWBS の定数、T11 σ_sw,fast、T03 σ_t,fast、T15 γ_coh） | 二重に数える・数え落とす危険 |
| **推奨: `BandLedger`** — 各雑音源のブロックが `bands()` で「時間ループで表した帯域 [f_lo, f_N]」と「その外の PSD（またはその積分の作り方）」を宣言する。帯域外の分散は**`phase_to_quadrature` の 1 か所だけ**で、指標ごと・ビームの組ごとに足す（χ_HF、パルスでは上限 ~1/τ_p、OPO では ~κ_c、conv §6）。T15 の可視度の補正（γ_coh）も台帳の 1 行として DC fringe ブロックが参照する | C2 の提案、T04 V13（台帳の値と細かい刻みの時系列で 0.01 dB 一致）、加えて「同じ雑音源・同じ帯域の 2 重登録はビルド時にエラー」の検査 |

### G06 損失の登録表（S §3.4 C3、T12 §8）

| 選択肢 | 根拠 |
|---|---|
| ブロックごとに η を掛ける | 重なり（空間 T12、時間 T03・T11、LO 波形 T17）を 2 回数える危険。GBS 入力と BHD を分けにくい |
| **推奨: `LossRegistry`（オフライン段の 1 つの表）** — 各ブロックは `losses()` で要素（id、η または IL_dB、区間 id、種類 `connector`・`splice`・`fiber`・`component`・`mode_overlap_spatial`・`mode_overlap_temporal`・`lo_shape`・`electronic`・`detector`、適用先 {`gbs`, `bhd`, `lo`, `probe`}、出典の印、`enabled`）を返す。登録表がビームの経路ごとに η を出す（squeezed → GBS_k、squeezed → BHD_k、LO_k → BHD_k、probe）。重なりの種類は経路の組ごとに 1 回だけ（2 回目はエラー）。squeezed 光の経路上の変調器・アクチュエータ・probe の入れ出し素子は必ず登録（prompt §5.5、D10、E55） | C3、T12 V13、D46（GBS 入力と BHD を分ける）。記録に出典の印付きで保存（T12 §8.2） |

### G07 光源の共通出力（S §3.4 C4）

| 選択肢 | 根拠 |
|---|---|
| 光源ごとに別の出力（T01 の Gaussian チャネル、T02 の R±(Ω)、T03 の超モード、T17 の 2N×2N 共分散） | 測定・指標のブロックが光源ごとに分岐する |
| **推奨: `SourceOutput`** = {モードの基底（`cw_sideband`: f_s の格子 / `temporal`: 時間の格子とモード関数 / `single_mode`）、その基底での共分散（CW は f_s ごとの 2×2、パルスは 2N×2N）、`project(lo_mode) -> (V_−, V_+, θ_off)`（CW なら f_s ごと）、誤差信号用の振幅利得 A±(P) の表、ポンプ経路の区間としての位相の寄与（結晶 ψ_p,xtal、OPO のポンプキャビティ −arctan(Δ_p/κ_p)）、φ_OPA,0、区間 × ビーム係数の表（OPO の離調、C5）、`tail_into_next_slot`（N11）}。`opo.pulsed = off` で T02 の結果に戻ることをテストにする（T17） | C4、T17 V-P1（CW の極限）、VT-03 |
| 単一モード（`single_mode_reference`、§7-1 の式）も同じ型（N = 1）で返す（保護リストの Ref A のブロック） | — |

### G08 TEC が書く量（S §3.4 C13）

| 選択肢 | 根拠 |
|---|---|
| TEC が位相 δψ を直接書く（T16 のモデル σ_act c_T L_h ΔT） | 軸ごとの差（PM の直交偏光 probe）と群遅延が同時に動くことを表せない |
| **推奨: TEC は小区間の温度 ΔT_h を `drv[i]` に書き、ψ（軸ごと）と τ_g は `segment_phase` が係数から計算する** | C13、T10 §6.1（TEC で軸の差が打ち消される、VT-39 の「TEC で 0」）。熱アクチュエータの係数は `coef[i][tec][beam]` に入る |
| 補足 | 比較用に「位相を直接書く」も `thermal_actuator(mode=phase_direct)` の変種として残すか — 推奨は**残さない**（同じことは係数の設定で再現でき、二重の経路が C2・C3 と同じ数え違いを生むため）。D34 に照らして要確認 |

### G09 パイプラインの合成の意味

| 論点 | 推奨 | 他の選択肢 |
|---|---|---|
| 並列の実行 | 書き込みのスロットが素であることをビルド時に検査し、宣言順に逐次実行（結果は順序に依らない） | 入力のスナップショットを各分岐に渡す（コピーのコスト、numba で不利） |
| 積み重ね | `Stack` は `Seq`・`Parallel` に展開する糖衣。線形の範囲で並列と等価なことをテスト（VT-44） | 専用の実行器 |
| 配線 | スロット名による明示の配線 ＋ ビルド時の 4 つの検査（§2.5） | 位置による暗黙の受け渡し（入れ子で壊れやすい） |
| フィードバック | `delay_buffer` を通る経路だけ。遅延 0 の帰還は禁止（代数ループ） | 遅延 0 を許して反復で解く（ループ遅延のモデル化の趣旨に反する） |

### G10 時間刻みと多重レート

| プリセット | f_sim の推奨 | 根拠 |
|---|---|---|
| `refD` | 1 MHz（DSC1 の ADC）、PID は 2 サンプルごと（500 kHz） | T16 §5-4。T15 §8.5 の Δt = 2 µs とは「ADC の 2 点平均を再現するか」の違いで、1 MHz はそれを含む（§8.4 の食い違い 1） |
| `refA`・`refA_tap` | 2–10 MHz（`fpga_generic`、τ_L 0.1–2 µs）。ヘテロダインは基底帯域、`rf_rate_mode` では f_sim 8–10 MHz | T14 §12.4、T16 §5-3 |
| `refB` | 2 MHz（PZT 14/22 kHz、dither 295 kHz を時系列に含める場合） | T16 |
| `refC` | 4 MHz（1 MHz dither を時系列に含めるなら f_N > f_m）。8–300 MHz の雑音電力は統計モデル | T04 §13、T14 §6.3 |
| 将来 | 1–4 MHz ＋ TEC・熱は 10–100 Hz | T08 §7.3、T16 |
| 収束の確認 | どのプリセットも f_sim を 2 倍にして主要な出力の変化が許容内（§6.4） | prompt §8 |

### G11 色付き雑音の生成と乱数の流れ

| 選択肢 | 根拠 |
|---|---|
| チャンクごとの FFT による色付け（T08 §7.3 の案） | チャンク長を変えると実現が変わり、チャンク長に対する収束の確認（prompt §8）が統計的にしか言えない。チャンクの境目で相関が切れる |
| **推奨: 状態を持つフィルタ（IIR・OU 過程の和で 1/f、共振ピークは 2 次の section）で時系列を逐次に作る。ブロックごとに `SeedSequence` の子（鍵 = ブロック名のハッシュ）** | チャンク長を変えても実現がビット単位で同じ（§6.4 の検査が厳密になる）。雑音源を切っても他の流れが変わらない（G13 の再実行の分解に必須）。FFT による色付けはオフラインの表と検証（PSD の照合）にだけ使う |
| 1/f の近似の精度 | OU の数を設定（帯域あたり 1–2 個/桁）。VT-37 V3（1/3 オクターブ帯で ±10%）を満たす数を WP-B2 が決める | T08 §9 |

### G12 プリセット・オーバーレイの合成規則と callable の書き方

| 論点 | 推奨 |
|---|---|
| 合成 | `preset` を土台に `overlays: [...]` を順に深くマージ。辞書は再帰的にマージ、リストは置き換え、ただし `name`・`id` を持つ要素のリスト（区間、ブロック）は名前でマージ。`null` で削除。合成の結果は `config.resolved.yaml` に保存 |
| 掃引 | `!sweep [a, b]`（線形）、`!sweep_log`、`!choice [..]` は「ここが掃引の軸」という印。単発の実行では各プリセットが値を 1 つ指定する（`sweep` が展開） |
| callable（YAML） | `!expr "0.00249 + 0.00222*x"`（`where:` で他のパラメータの点付きパスを束縛）。**評価は `ast` を使った制限付きの評価器**（四則、べき、`math`・`numpy` の初等関数、比較と三項演算だけ。属性アクセス・呼び出しの任意化・import は不可）— YAML は外から来る入力なので任意コードを実行しない |
| callable（登録関数） | `!fn {name: refB_loss_vs_x, args: {...}}`（Python 側で `@register_fn` した関数）。Python のリストでは本物の関数を渡せ、`to_spec()` は登録名のある関数だけ YAML に戻せる（名前のない lambda は往復できないと明示的にエラー） |
| 時間に依る callable | L_ic(x, P, t) のように時間ループの量に依るものは、オフライン段で格子上の表にし、ループでは補間する（`time_varying: true` の印） |

### G13 指標 (ii) の基準ポートと分解の方法

| 論点 | 推奨 | 根拠・他の選択肢 |
|---|---|---|
| 1 ポートの設定での (ii) | **1 ポートの設定も基準ポート 0（L_0 の経路、検出なし）を持つ**。「1 ポートで L_k を掃引」（D13）は「ポート 0 ＋ ポート k」の位相だけを計算する設定と読む。Ref A–D のようにポートの概念がない構成では (ii) は「該当なし」と出力する | D12 は「常に両方出す」、D31 は (ii) をポート間の差と定義。ポート 0 がないと (ii) が定義できない（§8.4 の食い違い 2） |
| 4 成分の分解 | 雑音源の組を切り替えた再実行（全部オン ＋ 4 回、同じ seed の流れ、G11）と線形性の検査（分散の和が全部オンと 5% 以内、T04 V12）。飽和・ロック外れがあれば「分解は無効」と出す | 1 回の実行の中で成分ごとに位相を足し分ける方法は、閉ループ（制御器が全成分の和に応答）では成り立たない |

### G14 新しい依存関係（それぞれユーザの決定）

現状: Python 3.14.8（`.python-version` 3.14、`pyproject.toml` は ≥ 3.12）、numpy 2.5.3、scipy 1.18.1（`uv.lock` はこの 2 つだけ）。2026-10-31 に `uv pip compile --python-version 3.14` で解決を確認（インストールはしていない）。

| 依存 | 用途 | 推奨 | 解決の確認 |
|---|---|---|---|
| `numba` | G01 の閉ループ | **追加**（G01 で (c) を選ぶ場合） | numba 0.68.1 が numpy 2.5.4 と解決（numpy の上限で 2.5.4 に上がる点に注意）。macOS arm64 の cp314 の wheel は WP-F0 でインストールして確認 |
| `pyyaml` | YAML の読み書き（prompt §5.2 の必須）。カスタムタグ `!expr` 等 | **追加**（標準ライブラリに YAML はない。`tomllib` は TOML で prompt の要件を満たさない） | 6.0.3 で解決 |
| `pytest` | テストの実行（dev のグループ） | **追加（dev のみ）**。標準の `unittest` でも書けるが、マーカー（`slow`、`vt`）と parametrize で VT の管理が楽 | 9.2.0 で解決 |
| `matplotlib` | 図（prompt §6 の図、Fig. 3–5 の比較） | **任意の追加（extra `plot`）**。結果は npz/JSON で残し、図は後処理 | 未確認（WP-F4 で） |
| `botorch`（＋ torch） | Phase 6 の最適化 | **Phase 6 まで追加しない**（重い。Phase 4 は抽象制御器インターフェースまで） | — |
| `h5py`・`xarray`・`pydantic` | 記録・スキーマ | **追加しない**（npz ＋ JSON と dataclass で足りる） | — |
| GNLSE の既存コード | 照合のみ | **追加しない**（NumPy 2 で動かない、T06 §4）。照合が要るなら別環境で一度だけ | — |

### G15 パッケージ構成

| 選択肢 | 推奨 |
|---|---|
| 平らなモジュール（`cranberry/*.py`） | 70 余りのブロックで見通しが悪い |
| **推奨: `src/cranberry/` のサブパッケージ（§4.1）、ブロックは族ごとのサブパッケージ、プリセットと機器はリポジトリ直下の YAML、テストは `tests/`** | 族ごとに worktree を分けやすい（ファイルの衝突がない、Phase 1 の運用と同じ） |
| パッケージ名 | `cranberry`（`pyproject.toml` の名前に合わせる）。`main.py` は WP-F0 で CLI（`python -m cranberry run <yaml>`）に置き換える |

### G16 小さなインターフェースの一括確認（S §3.4 の C5–C16 のうち未決のもの、ほか）

| # | 内容 | 推奨 |
|---|---|---|
| a | C5 区間 × ビームの係数の型 | 定数・callable・表のどれでも受ける（`coef` の値の型を 3 つ）。表は (Δ_s, Δ_p, x_ss) などの格子で補間 |
| b | C6 誤差信号の時間モデル | ブロックの API はループレートの電圧に統一。`rf_rate_mode` は検証用（T14 V1・V2・V6） |
| c | C10 逆向き伝搬 | 小区間の遅延の読み出しで自動的に出す。アクチュエータの位置は小区間の番号。テストは観測量を明示（T08 V6 と T10 テスト 2–3 を別の観測量として両方） |
| d | C11・C12 ブロックの統合 | `dc_fringe` は T15 の式 (3)（μ、PER、κ_1・κ_2、ε_R）を上位集合として 1 つに。レーザーは `laser_source` 1 つ（T05 と T15 の名前を統合） |
| e | C15 迷光の受け渡し | `leak[name]` のフィールドを T10 の迷光の表の入力（P̄_pr、タイミング、偏光、Δν_pr、ε_AOM）に合わせる |
| f | C16 測定のモード | MEAS に `esa`（CW のゼロスパン・掃引）と `pulse_sample`（パルス BHD の標本分散）の 2 つ |
| g | オフラインの表の保存 | npz ＋ JSON のメタ（格子、単位、出典の印、生成した設定の部分木のハッシュ、コードの版）。キャッシュは `runs/_cache/` |
| h | 実機のラッパ | Phase 4 は `ControllerInterface`・`SimController`・機種への写像の関数まで。`MokuController`・`DSC1Controller` は実機を使う時点で |
| i | 単位 | 内部は SI。パラメータ名に単位の接尾辞（`_hz`、`_s`、`_m`、`_w`、`_v`、`_rad`、`_db`）。dB の値は名前で区別（`il_db`）し、読み込み時に線形へ。スロットに単位のタグを持たせ、配線の検査に使う |
| j | 警告と停止 | CHK の警告（SBS、繰り返しの上限、振り分けの実現可能性、格子の端）は計算を止めずに記録に残す（D34、S §3.2 CHK） |

### G17 マイルストーンの順番（§1.2）

| 選択肢 | 推奨 |
|---|---|
| **推奨: M0 → M1（A ‖ D）→ M2（C）→ M3（B）→ M4（将来の時間ループ）→ M5（GNLSE・量子・パルス光源）→ M6（将来のプリセット）** | D17・D03 のとおり。ただし M5 の担当はフレームワーク（WP-F3）だけに依存するので、**M2 と並行して始めてよいか**を確認したい（D17 の「時間ループ → GNLSE」を「時間ループの完成を待つ」と読むか「優先順位」と読むか） |

---

## 4. モジュール一覧とインターフェース

### 4.1 パッケージ構成（G15）

```
pyproject.toml                    # deps per G14; [project.optional-dependencies] plot; [dependency-groups] dev
src/cranberry/
  __main__.py                     # CLI: run / sweep / ablate / decompose / convergence
  core/       state.py  layout.py  block.py  registry.py  pipeline.py  compile.py
              delay.py  rates.py  rng.py  units.py  fft.py (conv §8 wrappers)  checks.py
  config/     spec.py (RunSpec dataclasses)  yaml_io.py (tags !expr !fn !sweep ...)  overlay.py  expr.py  resolve.py
  runner/     run.py  sweep.py  ablate.py  decompose.py  convergence.py  record.py (run dir, git, env, seeds)
  offline/    tables.py (cache)  source_output.py  loss_registry.py  band_ledger.py  gaussian.py (covariance ops, pure loss, symplectic checks)
  path/       topology.py  segments.py (segment_phase, subsegments, drivers)  beams.py (beam_assembler, emission time)  locks.py  coef.py
  blocks/
    source/   noise/   pol/   mod/   err/   ctrl/   act/   det/   dmx/   loss/   fsr/   spf/   los/   fib/   meas/   chk/
  hardware/   interface.py (ControllerInterface)  sim.py (SimController)  mapping.py (canonical -> dsc1 / moku / fpga)
presets/      reference/  future/  axes/<axis>/  bad/  locks/
devices/      dsc1.yaml  moku_delta.yaml  fpga_generic.yaml
tests/        unit/<mirrors src>  vt/test_vt_XX_*.py  integration/  conftest.py (markers: slow, vt, compiled)
runs/         (gitignored)
```

### 4.2 フレームワークのモジュール（ponytail 禁止、保護リスト）

| モジュール | 役割 | 公開インターフェース | WP |
|---|---|---|---|
| `core/layout.py`・`state.py` | スロットの登録表、平らな配列、名前付きビュー（G02） | `SlotLayout.add(name, shape, dtype, unit, ring_len=0)`、`State.view(name)`、`state["a.b"]` | F1 |
| `core/block.py`・`registry.py` | `Block` の規約、`@register_block`、`make_<kind>` の解決 | §2.4 | F0 |
| `core/pipeline.py` | `Seq`・`Parallel`・`Stack`、ビルド時の 4 つの検査 | `build(pipeline_spec, layout) -> Block` | F1 |
| `core/compile.py` | カーネルの合成と内側のループ、純 Python の参照ループ、一致の検査 | `compile_loop(block) -> LoopFn`、`run_reference(block, n)` | F1 |
| `core/delay.py` | リング ＋ 端数の遅延（線形 / Thiran） | `DelayLine(tau_s, f_sim, method)`、カーネル用の関数 | F1 |
| `core/rates.py` | `rate_div` の扱い、保持 | — | F1 |
| `core/rng.py` | `SeedSequence` の鍵付きの子 | `stream(root_seed, block_name)` | F1 |
| `core/fft.py`・`units.py` | conv §8 の FFT、片側 PSD、dB 変換、単位のタグ | `fft_agrawal(A, dt)`、`psd_onesided(x, fs)` | F0 |
| `config/*` | `RunSpec`、YAML のタグ、オーバーレイ、制限付きの評価器、依存関係の解決 | `load(path_or_dict, overlays) -> RunSpec`、`dump_yaml(spec)`、`to_spec(blocks)` | F2 |
| `runner/*` | §2.8 の関数、記録 | `run(spec) -> RunResult` ほか | F1（`run`）、F4（`decompose`・`convergence`）、F2（`record`） |
| `offline/tables.py` | 表のキャッシュ、補間 | `Table(name, grid, values, unit, provenance)`、`cache.get_or_build(key, fn)` | F3 |
| `offline/source_output.py` | `SourceOutput`（G07） | §3 G07 | F3 |
| `offline/loss_registry.py` | G06 | `LossRegistry.add(LossElement)`、`.eta(path, metric)`、`.to_json()` | F3 |
| `offline/band_ledger.py` | G05 | `BandLedger.declare(BandDecl)`、`.hf_variance(metric, beams)` | F3 |
| `offline/gaussian.py` | 共分散の演算（純損失、ビームスプリッタ、symplectic の検査、R± ↔ 共分散） | `pure_loss(V, eta)`、`is_symplectic(S)` | F3 |
| `path/*` | 光路ごとの位相モデル（保護リスト）: 区間、小区間、要因の履歴、`segment_phase`、`beam_assembler`、ロックの線形結合、係数 | §2.3 | F6 |
| `blocks/meas/metrics.py` | 指標 (i)(ii)、分解のための雑音源の組の定義 | `metric_i`、`metric_ii` のブロック | F4 |
| `blocks/meas/phase_to_quadrature.py` | 位相 → squeezing の写像（帯域外は `BandLedger` から 1 回だけ） | T04 §13 | F4 |
| `hardware/*` | 抽象制御器インターフェース、`SimController`、写像 | §2.10 | F5 |
| `devices/*.yaml` | 機器プリセット | §2.10 | F5 |
| `presets/**` | プリセットとオーバーレイ | §2.9 | I1–I5 |
| `tests/vt/*` | 検証テスト（保護リスト） | §6 | 各 WP |

### 4.3 ブロック（族ごと）

列: 段（O オフライン、T 時間ループ、M 測定、Q 量子層、G GNLSE）、読む / 書く（§2.3.2 のスロット）、主なパラメータ（全体は S §3.2 の列を契約とする）、表（→ 生成、← 消費）、ノート、プリセット（A = refA、At = refA_tap、B、C、D、Fw = F_wg100/F_wg6、Fb = F_bor6、Fl = F_liu20、bad = S §2.5、すべて = 全部）。WP は §7。

#### SRC（光源・ポンプ）— WP-B1（単一モード・導波路 OPA）、B6（CW OPO）、B11（パルス）、B2（レーザー）

| ブロック | 目的 | 段 | 読む | 書く | 主なパラメータ | 表 | ノート | プリセット |
|---|---|---|---|---|---|---|---|---|
| `laser_source` | レーザーの周波数雑音・位相・RIN（T05 と T15 を統合、G16-d） | T | `t` | `laser[s]` | h0、h₋₁、ピーク、ν̇、dither、SBS、RIN、機器プリセット `cobrite_dx1_{N,SC,GC}`・`koheras_harmonik_refA`・`tisa_mbr110_refB` | → `laser_oob_variance`（帯域の台帳へ） | T05、T15 | すべて |
| `single_mode_reference` | §7-1 の解析式 R± = L + (1 − L)e^{±2√(αP)}、R′±（保護リスト） | O | — | `SourceOutput`（N = 1） | L、α、P、θ | — | T01、T04 | A、At（検証） |
| `opa_waveguide` | 単一パス導波路 OPA（集中 / 分布損失、位相不整合、結合波方程式の枯渇は opt） | O | — | `SourceOutput`、A±(P)、ψ_p,xtal | α、P（in/out の基準）、ℓ、a_s、a_p、Δk(T)、φ_OPA,0、L_wg(P) ƒ、`distributed_loss` | → OPA 表（P、T） | T01 | A、At、Fw（代用） |
| `opa_supermodes` | パルス導波路 OPA の超モード（2N×2N 線形 ODE ＋ Takagi） | O | — | `SourceOutput`（temporal） | 分散、ポンプ包絡線、格子 512–3200 | → `mode_overlap_table` | T03 | Fw |
| `opo_cw_spectrum` | CW OPO の R±(f_s)、2×2 S(Ω)、L_ic(x) ƒ | O | — | `SourceOutput`（cw_sideband）、A± | T、L_ic ƒ(x, P, t)、l または線幅、P_th、Δ、x = √(P/P_th) | → R±(f) | T02 | B、C |
| `pulsed_opo_source` | パルス励起 OPO（時間領域の量子 Langevin、初期の真空の項） | O | — | `SourceOutput`（temporal）、u_opt、K、ε_tail | Δν_c、ρ、κ_p、ポンプ波形、x_ss、Δ_s、Δ_p、φ_p,in、δt = 1/(40–160 κ_c) | → `opo_angle_map` の表 | T17 | Fb、Fl、bad（`bad_bor_100MHz`、`bad_long_pump_opo`） |
| `opo_detuning_noise` / `opo_length_lock` | OPO の長さ雑音とロック（ロック光 `borealis`・`liu`） | T | `drv[opo]`、温度 | `opo[j].Delta_s/Delta_p` | 長さ雑音、ロック光の種類、Δ(dn/dT) | — | T02、T17 | B、C、Fb、Fl |
| `opo_angle_map` | OPO 区間の区間 × ビーム係数の表を引いて ψ を書く | T | `opo[j]` | `psi[opo][beam]` | x_ss、ロック状態 | ← (Δ_s, Δ_p, x_ss) の表 | T17、T02 | B、Fb、Fl |
| `pump_power_stabilization` | AOM によるポンプパワーのロック | T | `pump.P`、RIN | `pump.x_ss` | ループのパラメータ | — | T17 | Fb |
| `crystal_thermal` | 結晶の温度 → ポンプ経路の位相と Δk（D42 (a)） | T | 温度 | `drv[xtal]`（→ `psi[xtal][pump]`） | 結晶長、dn/dT、TEC | — | T01 | A、Fw（`L-LiuPD-A(a)`） |

#### PATH（光路・区間・ビーム）— フレームワーク WP-F6（保護リスト「光路ごとの位相モデル」）、偏光は WP-B7、フィルタキャビティは WP-B2

| ブロック | 目的 | 段 | 読む | 書く | 主なパラメータ | 表 | ノート | プリセット |
|---|---|---|---|---|---|---|---|---|
| `segment_phase` | 要因（ΔT、ε、TEC、ストレッチャ、離調）と係数から ψ（区間 × ビーム、小区間ごと）と τ_g | T | `drv[i]`、`coef` | `psi[i][b]`、`tau_g[i]` | 空間構造、ゾーン、`N_sub` | ← `fiber_env_coeffs`、`pm_fiber_params` | T08、T10、conv §9 | すべて |
| `beam_assembler` | 放出時刻で φ_L を読み、経路の ψ を各小区間の通過時刻で足し、Jones を掛ける | T | `laser.hist`、`psi`、`tau_g`、`J` | `beam[b]` | 経路、向き、観測点 | — | conv §7・§9、T05 | すべて |
| `delay_buffer` | リング ＋ 端数の遅延（フィードバック・ループ遅延） | T | 任意のスロット | 遅らせた値 | τ、補間 | — | T05、T16 | すべて |
| `segment_jones` | 区間の J（PM 対角 / SMF 拡散 ＋ 温度）（G04） | T | `drv[i].dT`、D_pol | `J[i]` | 種類、K_T、L_b、D_pol | ← `pm_fiber_params` | T09 | D、Fw、Fb、Fl、bad（`bad_smf_orth`） |
| `filter_cavity` | LO のフィルタキャビティ（位相 arg T(δ)、透過、H_FC） | O＋T | 離調 δ(t) | `psi[fc][lo]`、透過 | Δν_FC、FSR、T0 | — | T05 | A |
| `polarizer_projection`、`pdl_element`、`faraday_nonreciprocal`（opt） | 偏光子の射影とパワーのドリフト、PDL、非相反 | T | `beam` の Jones | 透過パワー、射影振幅 | ER、PDL、V_F | — | T09 | D（偏光子）、比較 |

#### NOISE（雑音源）— WP-B2

| ブロック | 目的 | 段 | 読む | 書く | 主なパラメータ | 表 | ノート | プリセット |
|---|---|---|---|---|---|---|---|---|
| `ambient_temperature`、`segment_thermal` | ゾーンの温度、区間の温度（1 次の遅れ） | T（10–100 Hz） | `t` | `drv[i].dT` | r_T、b_T、τ_th、周期成分 | — | T08 | C、D、F |
| `segment_mech_noise` | 機械・音響の歪み（区間間の相関 ρ_ij(f)、区間内の空間構造） | T | `t` | `drv[i].strain` | h_ref、L_ref、γ、f_c、f_lo、ピーク、ρ_0、f_ρ | — | T08 | C、D、F |
| `segment_thermodynamic`（opt） | 熱力学の基底 | T | — | `drv[i]` | Dong の物性値 | — | T08 | （比較） |
| `gawbs_phase_jitter`（opt） | GAWBS の位相（ループ帯域外は台帳へ） | T＋台帳 | — | `drv[i]` | σ_G² | ← `gawbs_table` | T07 | `bad_copropagating_probe` |
| `rf_reference` | AOM・復調の RF 位相の雑音（遅延付き） | T | — | `rf[j].phi` | S_φRF、τ_RF | — | T14 | A、B、C |
| `pd_offset_drift` | PD・TIA のオフセットのドリフト | T | 温度 | `v[...].offset` | 15 mV、ドリフト | — | T13、T15 | D |

#### MOD（変調器）— WP-B4

| ブロック | 目的 | 段 | 読む | 書く | 主なパラメータ | 表 | ノート | プリセット |
|---|---|---|---|---|---|---|---|---|
| `eom_phase` | 位相変調（アクチュエータとしても）。squeezed 経路なら損失を登録 | T（損失は O） | `delay[...]`（u） | `drv[i]`（要因 `eom` の δψ） | V_π,LF/RF、IL、RAM | — | T16、T14 | A、C、`bad_eom_on_squeezed` |
| `dither_modulator` | 指定した経路に m sin ω_m t | T | — | 区間の δψ | m、f_m、a_RAM、φ_a、`dither_mode: timeseries \| analytic` | — | T14、T04 | B、C |
| `aom_shift`、`probe_pulser` | 周波数シフト・パルス化（`leak` を書く） | T | δf_RF / φ_RF、ゲート | 位相、周波数、`leak` | 次数、ε_AOM、デューティ、t_offset | — | T16、T10 | A、Fw、Fb、Fl |

#### ERR（誤差信号 #1–#6）— WP-B3

| ブロック | 目的 | 段 | 読む | 書く | 主なパラメータ | 表 | ノート | プリセット |
|---|---|---|---|---|---|---|---|---|
| `parametric_beat_detector` | #1 PD-OPA / tap の PD（DC ＋ 2f_sh I/Q、dither の 1 次） | T | `beam`（φ_p、φ_pr を相互作用点で）、A± | `v[...]` | f_sh、T_tap、ℜ、G、i_n | ← A±(P) | T14、T01 | A、At、Fw（`L-LiuPD-A`）、`bad_liuPD_B` |
| `bhd_probe_beat` | #1c・#3 BHD での probe–LO のビート・干渉 | T | φ_LO、φ_pr、φ_p、A± | `v[...]`（f_sh / f_m の I/Q） | BHD の R(f) | — | T14 | A、At、B、`bad_refA_seed_in_gbs` |
| `quantum_noise_power` | #4 量子雑音ロック（帯域の雑音電力 → 包絡線・ログ ＋ radiometer 雑音） | T | Δ | `v[...]` | f1、f2、κ_det、検出器の法則、ビデオ帯域 | ← P(Δ)、∫S_x²(Δ) | T14 | C |
| `dc_fringe` | #5 DC fringe（T13・T14・T15 を統合、G16-d） | T | 2 ビームの位相・パワー・Jones | `v[...]` | κ_1、κ_2、T_opt、ℜ、ε_R、G、V_off、飽和、μ、符号 s | ← 帯域の台帳（γ_coh） | T13、T14、T15 | D、Fl（`homodyne_dc`） |
| `liu_probe_interferometer` | #6 Liu 型（逆向き probe と参照光、`homodyne_dc`・`heterodyne`・`dither`、ゲート、標本化 N09） | T | φ_pr,k（逆向き、a⊥）、φ_ref,k | `v[...]`、`leak` | mode、δf、デューティ、ゲート | — | T14、T10 | Fw、Fb、Fl、bad（`bad_mw_probe_samepol` など） |
| `demodulator` | 共通の復調（I/Q、LPF、オフセット） | T | `v`、φ_d（`rf`） | `dem[...]` | f_ref、LPF、V_off | — | T14 | A、At、B、C |

#### CTRL（制御器・デジタル）— WP-B4（積み重ね・並列の合成はフレームワーク F1）

| ブロック | 目的 | 段 | 読む | 書く | 主なパラメータ | 表 | ノート | プリセット |
|---|---|---|---|---|---|---|---|---|
| `pid` | 正準係数、離散化 `dsc1`・`moku`・`tustin`、符号、制限、anti-windup | T | e（`v`・`dem`・前段の `ctrl`） | `ctrl[...]` | §2.10 | — | T16、T15 | すべて |
| `adc`、`dac` | 量子化、レンジ、平均、雑音、LPF、出力抵抗と負荷 | T | 電圧 | 量子化した電圧 | bits、range、f_s、平均点数、雑音、R_out | — | T16 | すべて（D49: 量子化オン） |
| `iir`、`hpf`、`lpf`、`adder`、`ramp` | フィルタ（SOS）、加算、ランプ | T | 任意 | 任意 | SOS、周期・振幅 | — | T16 | C（ランプ）、並列の帯域分割 |
| `lock_supervisor` | ロック外れの検出、再引き込み（`reset_center`・`fringe_hop`・`ramp_search`） | T | e、u、飽和、真の位相 | `lock[...]` | しきい値、待ち時間、モード | — | T16 | C、D、F |

#### ACT・TEC（アクチュエータ）— WP-B4

| ブロック | 目的 | 段 | 読む | 書く | 主なパラメータ | 表 | ノート | プリセット |
|---|---|---|---|---|---|---|---|---|
| `pzt_fiber_stretcher` | 915B など（ΔOPL 読み、機械モード、ヒステリシス、クリープ、飽和）。挿入損失を登録 | T（損失は O） | `delay`（u） | `drv[i].stretcher`（小区間 x_a）、`act` | K_s（6.49 rad/V）、σ_act、モード（f_i 5–60 kHz、Q 3–100）、ξ_eo の代替読み | ← 914 Fig. 1 の表 | T16、T15 | C、D、F |
| `hv_amp_914` | 914（2 極 29.5/112 kHz、1 極 28 kHz は代替、スルー） | T | u | 駆動量 | 極、スルー、入力レンジ、雑音 | ← 同上 | T15、T16 | D |
| `pzt_mirror` | PZT ミラー（共振 14・22 kHz、10–20 kHz） | T | u | δψ | d_x、θ_i、モード | — | T16 | B、C、Fl（変種） |
| `thermal_actuator` | TEC・スプール温度（**温度を書く**、G08） | T（間引き） | u | `drv[i].dT`（小区間） | K_T、τ_th、τ_tr、ΔT_max、L_h | — | T08、T10、T16 | F、`bad_orth_pm_no_tec`（TEC を外した比較） |

#### DET（検出器）— WP-B5

| ブロック | 目的 | 段 | 読む | 書く | 主なパラメータ | 表 | ノート | プリセット |
|---|---|---|---|---|---|---|---|---|
| `detector_efficiency` | η_det(f) と内訳（QE、クリアランス C(f)、ディップ、CMRR）。損失を登録 | O | — | 登録表 | η_q、C(f)、δ(f)、S_RIN | — | T13 | すべて |
| `bhd_frontend` | 差電流の時系列・PSD（TIA、ディップ 0.96/0.83 GHz、回路雑音、暗雑音の差し引き） | M | 直交位相・R(f) | 時系列・PSD | H(f)、S_el、G、`subtract_dark` | — | T13 | A、B、C、F |
| `pulsed_bhd_kernel` | パルス BHD の c_m、パルス間クロストーク、ゲート、積分型の薄まり（N10） | O | — | 表 | h(t) または B・Q、w(t)、σ_el、N_LO、`mode: pulse_resolved\|gated\|integrating` | → 表 | T13、T11、T17 | Fw、Fb、Fl、`bad_integrating_bhd` |

#### DMX（振り分け）— WP-B8

| ブロック | 目的 | 段 | 読む | 書く | 主なパラメータ | 表 | ノート | プリセット |
|---|---|---|---|---|---|---|---|---|
| `demux_design` | L_k = L_0 + k·N_b·c·T_rep/n_g、η_sw,k、ψ_sw,k^(0)、漏れ行列、実現可能性 | O | — | 表、登録表 | K、f_rep、N_b、L_0、素子プリセット（`pockels_tree`・`pockels_block`・`passive_loops`・`ln_fiber_tree`）、`tree`/`chain` | → `demux_design` | T11 | Fw、Fb、Fl、`bad_ln_switch_tree` |
| `demux_switch`、`delay_line_k` | 駆動の位相雑音・バイアス、τ_k(t)、δt_k(t) | T | 温度、駆動雑音 | `drv[demux]`、`timing[k]` | κ_pm、V_π、バイアスのモデル | ← `demux_design` | T11 | 同上 |

#### LOSS（損失・モード整合）— 登録表は F3、要素の計算は WP-B5（検出系）・B7（PER・巻線）・B1（tap）

| ブロック | 目的 | 段 | 読む | 書く | 主なパラメータ | 表 | ノート | プリセット |
|---|---|---|---|---|---|---|---|---|
| `loss_budget`（= `LossRegistry` の要素の集合を設定から作る） | 経路ごとの η（GBS / BHD / LO / probe） | O | 全ブロックの `losses()` | 登録表 | 要素の表（S §2.3 の `NF`・`MF`・`LL` を含む） | → η の表 | T12 | すべて、`bad_all_fiber_NF` |
| `pure_loss` | σ → ησ + (1 − η)/2·𝟙 | Q/O | — | 共分散 | η | — | T12、T07 | すべて |
| `mode_overlap` | 空間モードの重なり η_mm | O | — | 登録表（`mode_overlap_spatial`） | ビーム径・ずれ、または可視度 | — | T12 | A、B、C |
| `per_budget`、`coil_birefringence` | PER の合計と偏光の損失、巻線の曲げ | O | — | 登録表、`coef` | 部品の PER、位置、巻き半径 | — | T09 | D、F |
| `tap` | conventional tap（r_t の損失 ＋ P_tap を誤差信号へ） | O＋T | — | 登録表、`v` | r_t、P_pr | — | T12、T14 | At |
| `etalon_noise`（opt） | エタロンの位相・強度リップル | T | 温度 | `drv` | R_1、R_2、L_e | — | T12 | （比較） |

#### FSR・SPF・LOS — WP-B11

| ブロック | 目的 | 段 | 読む | 書く | 主なパラメータ | 表 | ノート | プリセット |
|---|---|---|---|---|---|---|---|---|
| `fsr_comb_filter` | OPO の FSR 櫛のフィルタ（FBG の組、MZI の段） | O | — | 登録表、n_FSR | FSR、r_m ƒ、T_f(ν) | — | T17 | Fb、Fl、`bad_no_fsr_filter` |
| `spectral_filter`、`dispersion_comp` | GBS 前の単一モード化、ポートごとの分散補償・遅延の微調整 | O | — | 登録表、補償後の φ2・φ3 | 帯域、GDD | — | T03、D44 | Fw |
| `lo_temporal_shape` | LO の時間波形（`optimal`・`gauss`・`rect`・`one_sided_exp`）と重なり | O | — | 登録表（`lo_shape`）、θ_off | 波形、BHD 帯域 | — | T17、T03 | Fw、Fb、Fl、`bad_gauss_rect_lo` |

#### FIB（ファイバのオフライン計算）— WP-B9（古典）、B10（量子・散乱）、B11（モードの重なり）

| ブロック | 目的 | 段 | 読む | 書く | 主なパラメータ | 表 | ノート | プリセット |
|---|---|---|---|---|---|---|---|---|
| `fiber_gnlse` | RK4IP、項のトグル（保護リスト）、CW の近道 | G | — | A(L, T)、保存量 | `terms.*`、刻み、格子、Raman モデル（Blow–Wood / Lin–Agrawal） | → 古典解 | T06 | Fw（LO・probe）、C（SPM） |
| `fiber_linear_table`、`lo_nonlinear_table`、`xpm_counterprop_coeff`、`spm_rin_to_phase` | 線形の表、LO の非線形、逆向き XPM、RIN → 位相 | O | — | 表 | P_avg × τ × L_k | → 表 | T06 | すべて（線形）、Fw |
| `fiber_quantum_linearized` | 線形化の量子雑音層（G03） | Q | 古典解、V(0)、LO モード | V(L)、加わった雑音 N | 雑音のフラグ | ← 散乱表 | T07 | C、F |
| `fiber_quantum_tw`（opt） | truncated Wigner（検証用） | Q | 同上 | 同上 | N_traj | — | T07 | （検証） |
| `raman_noise_table`、`gawbs_table`、`brillouin_table`、`stray_light_table` | 散乱雑音と迷光の表 | O | — | n(f)、N_mode、N_stray(P̄_pr, …) | モデル、L、T、P、ε_pol、b_R | → 表（`metric` の迷光へ） | T07、T10 | C、F、bad |
| `path_mode_transfer`、`mode_overlap_table`、`timing_state`、`mode_overlap_lookup` | パルスの時間モードの伝達、ポート間の重なり、タイミング | O / T | 区間の列、Δt_k(t) | (V_−, V_+, θ_off)(Δt)、ポート間 η | 格子、ジッタ | → / ← 重なりの表 | T03 | Fw、Fb、Fl |

#### MEAS（測定と指標）— フレームワーク F4（指標と写像の枠）、WP-B5（ESA の中身）

| ブロック | 目的 | 段 | 読む | 書く | 主なパラメータ | 表 | ノート | プリセット |
|---|---|---|---|---|---|---|---|---|
| `phase_to_quadrature` | δ(t) → V′(t)、R′±（帯域外は台帳から 1 回、G05） | M | δ（指標ごと）、V± の表 | V′、R′±、帯域ごとの σ | `average: exact\|small_angle`、`include_hf_phase_noise`、`dither_mode` | ← `vpm_table`、台帳 | T04 | すべて（D は換算表のみ） |
| `hf_phase_noise_integral`、`vpm_table` | χ_HF、損失込みの V±(f)（指標 (i)(ii) で損失の集合を変える） | O | — | 表 | `sideband_weighting` | → 表 | T04 | すべて |
| `esa_zero_span`、`esa_swept_spectrum`、`pulse_sample` | ESA の模擬（放射計雑音、log バイアス）、パルス BHD の標本分散（G16-f） | M | V′(t)、V′(f) | トレース | f_s、RBW、ENBW、VBW、平均の種類、seed | — | T04、T13 | A、B、C（esa）、F（pulse_sample） |
| `metric_i`、`metric_ii`、重なり、迷光 | 指標（§2.3） | T/M | `beam`、`timing`、`leak` | `metric` | 分解の雑音源の組（G13） | ← 迷光の表 | T04、T10、D47 | すべて（(ii) は F） |
| `required_theta_contour` | 目標 × 長さ → σ_req（D16） | O | — | 等高線 | ポンプ最適 / 固定、dither | ← 損失 | T04、T12 | F |

#### CHK（整合性の警告）— WP-B10（SBS）、B11（繰り返しの上限）、B8（振り分け）、B9（格子の端）

| ブロック | 目的 | 段 | ノート | プリセット |
|---|---|---|---|---|
| `sbs_check` | SBS の閾値の警告 | O | T07 | C、Fw |
| `rep_rate_tail_check` | OPO の裾による繰り返しの上限（Borealis 型 ≲ 45 MHz、Liu 型 ≲ 87 MHz） | O | T17 | Fb、Fl、`bad_bor_100MHz` |
| `demux_feasibility` | 立ち上がり・同期の実現可能性（PC で N_b = 1 は f_rep ≤ 6 MHz） | O | T11 | Fw（N_b = 1 で警告）、Fb、Fl |
| `grid_edge_monitor` | GNLSE の格子の端（≤ 10⁻⁶） | G | T06 | GNLSE を使う構成 |

### 4.4 オフラインの表の受け渡し（S §3.3 を契約にする）

| 表 | 作る | 使う | 次元 |
|---|---|---|---|
| `SourceOutput` | SRC の O ブロック | `vpm_table`、`mode_overlap_lookup`、`quantum_noise_power`、ERR（A±） | 光源ごと（G07） |
| V±(f) の損失込み（(i) 用と (ii) 用） | `vpm_table` | `phase_to_quadrature`、`quantum_noise_power` | f |
| P(Δ)、∫S_x²(Δ) | `quantum_noise_power` の前処理 | `quantum_noise_power` | Δ |
| 帯域外の分散（σ²_oob、χ_HF、γ_coh） | `BandLedger` ＋ `laser_oob_variance`・`hf_phase_noise_integral` | `phase_to_quadrature`、`dc_fringe` | 指標 × ビームの組 |
| N_stray、n_mode | `stray_light_table` | `metric`（迷光）、squeezing の劣化 | P̄_pr、偏光、Δν_pr、L、T_gate |
| LO の非線形 | `lo_nonlinear_table` | `mode_overlap_lookup`、`phase_to_quadrature` | P_avg × τ × L_k |
| 振り分けの設計 | `demux_design` | `demux_switch`、登録表 | ポート |
| η（経路 × 指標） | `LossRegistry` | すべての V± | 経路 |
| (V_−, V_+, θ_off)(Δt, φ2)、OPO は (Δ_s, Δ_p, x_ss) | `mode_overlap_table`・`opo_angle_map` の表 | `mode_overlap_lookup`、`opo_angle_map` | 格子、チャンクごとに更新で足りる（T03 §8.2） |

---

## 5. データフロー

記号: `seg:X` は区間、`→` は信号、`◄` はアクチュエータの書き込み（`drv` 経由）、`[ ]` はブロック、`(τ)` は遅延バッファ。ロックは「ビームの位相の線形結合」で書く（conv §9）。

### 5.1 Ref A（`refA`、ロックのトポロジー `L-A`）

```
 [laser_source] φ_L (1545.32 nm, CW)
      │ seg:L0 (common)
      ├───────────────────────────┬─────────────────────────────────────┐
 seg:SHG (EDFA+SHG, pump path)   seg:PR (AOM1,AOM2: +f_sh; EOM_pr ◄act1) seg:LO (EOM_LO ◄act2; [filter_cavity])
 φ_p = 2φ_L + ψ_SHG              φ_pr                                   φ_LO
      └───────────┬───────────────┘                                     │
            seg:CMB (pump + probe, free space; tapping mirror)          │
        ┌─────────┴───────────┐                                         │
   seg:NC_sq               seg:NC_det      (non-common, coef per seg×beam, E68)
   seg:XTAL_sq (pump)      seg:XTAL_det (pump)                          │
   [opa_waveguide]         [PD-OPA] → [parametric_beat_detector] → v.pd1
   θ_sq=(φ_p+φ_OPA0)/2+π/2       → [demodulator 2f_sh ← rf_reference] → dem.l1
        │                        → [adc]→[pid lock1]→[dac]→(τ_L1)→ act1 [eom_phase on seg:PR]
   seg:SQ (squeezed + residual probe) ───────────► BHD ◄──────────────────┘ (LO)
                                    [bhd_probe_beat] → v.bhd → [demodulator f_sh] → dem.l2
                                    → [adc]→[pid lock2]→[dac]→(τ_L2)→ act2 [eom_phase on seg:LO]
 lock1: 2φ_pr − φ_p − φ_OPA0,det  (at PD-OPA)        lock2: φ_LO − φ_pr  (at BHD)
 METRICS at BHD observer: metric (i) Δ = φ_LO − θ_sq → [phase_to_quadrature] (V± from SourceOutput + LossRegistry[bhd])
                          → [esa_zero_span] (3 MHz, RBW 1 MHz, VBW 100 Hz), [esa_swept_spectrum] (0–100 MHz), Fig.5 sweep over P
 `refA_tap`: PD-OPA branch removed; [tap] on seg:SQ (r_t → LossRegistry) → [parametric_beat_detector] = lock1 (L-tap)
```

### 5.2 Ref D（`refD`、`L-D`、ブラインド前は解析的な検証のみ）

```
 [laser_source: cobrite_dx1_SC] φ_L, δp
      │ seg:SMF0 (SMF, [segment_jones] D_pol drift) → [polarizer_projection] (power drift δP/P = −Θ²/4)
      ▼
   PNH1550R5A2 #1 (κ1 = 0.5 ± 0.015, PER ≥ 25 dB → per_budget)
      ├── arm1 "pr": seg:A1_pig + seg:A1_915B (PM, N_sub, ΔT/strain from [segment_thermal][segment_mech_noise]) ◄ act
      └── arm2 "LO": seg:A2 (PM, ΔL = L_A1 − L_A2 swept 0–15 m)
      ▼
   PNH1550R5A2 #2 (κ2) → [dc_fringe = PDB450C, G = 1e3 V/A, V_off(t), ±4.6 V sat, sign s]
      │ v.pd
      → [adc: device dsc1, 1 MHz, 2-pt avg, 16 bit] → [pid: dsc1 (e = S − v inside), every 2 samples]
      → [dac: dsc1 12 bit ±4 V, 100 kHz LPF, 220 Ω] → (τ = 5.7 µs; 8.3 µs variant) → [hv_amp_914 2-pole]
      → [pzt_fiber_stretcher 915B, 6.49 rad/V, modes 5–60 kHz Q 3–100] ◄ writes drv[A1_915B].stretcher
      [lock_supervisor] (range ±26 rad → lock loss; relock mode)
 lockD: φ_D = φ_LO − φ_pr = [φ_L(t − τ_A2) − φ_L(t − τ_A1)] + ψ_A2 − ψ_A1  → setpoint +π/2 (conv §5)
 METRICS: true φ_D(t) + PSD, e_DSC1(t), CV/u(t), lock events. No squeezing. (VT-13, 14, 21, 24)
```

### 5.3 Ref C（`refC`、`L-C`）

```
 pump ─► [opo_cw_spectrum: monolithic ppKTP, 1.75 GHz FWHM] squeezed vacuum (no probe)  θ_sq
            │ seg:SQF (SMF 1 km; [segment_thermal][segment_mech_noise], ρ(f) with seg:LOF)
            ▼
 [laser_source] → seg:LOF (SMF 1 km) → act1 [eom_phase] (+ [dither_modulator] 1 MHz, m = 0.0785 rad)
                 → act2 [pzt_mirror] → act3 [pzt_fiber_stretcher]  (placement per Ref C Fig. 3; WP-I3 confirms)
            ▼                                                 ▼
          BHD ◄────────────────────────────────────────────── LO
   [bhd_frontend H(f), dips 0.96/0.83 GHz] → [quantum_noise_power: 8–300 MHz, log detector, radiometer noise]
   → [demodulator 1 MHz] → e
   Stack: e → [pid1] ─┬─► (τ) act1 EOM (kHz)
                      └─► [pid2] (+[ramp]) ─┬─► (τ) act2 PZT mirror
                                            └─► [pid3] ─► (τ) act3 stretcher (Hz–sub-Hz, range → lock time)
   [lock_supervisor] (stretcher range → lock loss → ramp_search reacquire)
 lockC: Δ = φ_LO − θ_sq (the BHD noise power itself; no probe)
 METRICS: metric (i) Δ → [phase_to_quadrature] (dither: timeseries or analytic J0, not both) → [esa_swept_spectrum 30 MHz–1.21 GHz]
          residual decomposition: dither 0.0555 rad + rest ≈ 0.1 rad total (VT-19, 20); lock duration (VT-21)
```

### 5.4 将来・32 ポート（`F_liu20`: Liu 型 OPO ＋ `L-Liu` ＋ `L-LiuPD-A(d)`）— 基準ポート 0 とポート k を描く（32 は k の繰り返し）

```
 [laser_source] φ_L ─┬─ seg:SHG ─► seg:PUMP ─► [pulsed_opo_source: Liu-type, 125 MHz, T_p 1.6 ns]  (OPO segment, coef table N02)
                     │                           │ transmitted pump ─► homodyne vs 2ω ref ─► [pid] ─► (τ) PZT on seg:PUMP   [L-LiuPD-A(d)]
                     │                           │ OPO length lock (Liu type: same pol/freq, time-interleaved) [opo_length_lock]
                     │                           ▼ squeezed train
                     │                 seg:UP (squeezed-only, before demux) ─► [demux_switch: pockels_block N_b = 4] (drv[demux], loss → registry)
                     │                     ├─ port 0: seg:D0 ─ seg:F0 (L_0) ─ PBS0 ─ seg:G0 ─► GBS input 0
                     │                     └─ port k: seg:Dk ─ seg:Fk (L_k = L_0 + k N_b c T_rep / n_g; PM; N_sub; zone) ─ PBSk ─ seg:Gk ─► GBS input k
                     │                                            ▲ PZT at x_a (upstream)   ▲ TEC on subsegments (writes ΔT_h)
                     ├─ seg:REF (common reference light delivery, common point 2) ─► ref_k at detectors near demux
                     ├─ seg:PRsrc (same laser, [aom_shift]+[probe_pulser], orth. pol) ─ seg:PRsplit (common point 1, near GBS)
                     │      └─► enters port k at PBSk, travels backwards (d = −1, axis a⊥) through seg:Fk then seg:Dk
                     │            ─► exits near demux ─► [liu_probe_interferometer k] (vs ref_k; gate; sampling N09) ─► v.k
                     │                 → [adc] → [pid_k fine] ─► (τ_L + n_g L_k / c) ─► PZT on seg:Fk
                     │                                └─► [pid_k slow, rate_div] ─► [thermal_actuator] ─► drv[Fk].dT
                     └─ LO_k delivery: seg:LOk (pulsed LO, split + delay, δτ_k, δL_k) ─► BHD_k ([pulsed_bhd_kernel: pulse_resolved])
 lock_k: φ_pr,k − φ_ref,k = C_k   (computed through backward subsegment sampling; cancellation not assumed)
 METRICS:
   (i)  at BHD_k: δ_i,k = φ_LO,k − θ_sq,k
   (ii) at GBS inputs, same GBS shot: Θ_k = θ_sq,k − θ_sq,0 ; components by reruns (G13):
        (1) locked path residual [seg:Fk..]  (2) common points [seg:PRsplit, seg:REF]
        (3) laser decorrelation out of loop band + in-band 16 sin² term (T_a − T_b)  (4) squeezed-only [seg:UP, demux, seg:Gk]
   overlap (mode_overlap_lookup), stray light (leak → stray_light_table: clicks/gate, photons/mode), η_k (gbs / bhd)
 CHK: rep_rate_tail_check (20 MHz < 87 MHz OK), demux_feasibility, sbs_check
```

（図の注: probe は GBS 側の PBS_k から入り、ポート k の `seg:Fk`・`seg:Dk` を逆向きに通り、振り分けの手前で取り出す（S §2.2「検出は振り分けの手前」）。`seg:Gk` を probe が通るかは probe の入れ点の位置のパラメータで、通らない区間は指標 (ii) の成分 (4) に入る。）

---

## 6. テスト計画

### 6.1 テストの種類と置き場所

| 種類 | 中身 | 置き場所 | 実行 |
|---|---|---|---|
| 単体 | ブロック 1 つ・関数 1 つの解析解との照合、符号、単位、トグルで恒等になること | `tests/unit/<族>/` | 毎コミット（数秒） |
| VT | S §4 の VT-01–45（保護リスト）。1 ファイル 1 VT（`tests/vt/test_vt_13_closed_loop_tf.py`）、docstring に出典のテスト番号と許容幅 | `tests/vt/` | マーカー `vt`。重いものは `slow` |
| 統合 | プリセット単位の実行（Ref A–D、将来）、YAML ↔ Python の往復、チャンク長・f_sim の不変性 | `tests/integration/` | `slow` |
| フレームワークの性質 | ビルド時の 4 つの検査、代数ループの拒否、並列の書き込みの衝突、積み重ね = 並列（線形）、コンパイル = 純 Python（10⁻¹²）、seed の再現、乱数の流れの独立性 | `tests/unit/core/` | 毎コミット |

### 6.2 VT と モジュール・WP・マイルストーンの対応

| VT | §7 | 主なモジュール | 種類 | WP | M |
|---|---|---|---|---|---|
| VT-01 | 1 | `single_mode_reference`、`phase_to_quadrature` | 単体 | B1、F4 | M1 |
| VT-02 | 1 | `refA` 全体（参照の損失収支、lock を 9 mrad に） | 統合 | I1 | M1 |
| VT-03 | 1 | CW の極限（パルスのブロックが恒等）、GNLSE の CW 近道 | 単体＋統合 | B9、B11 | M5 |
| VT-04 | 2 | Fig. 5（4 本 × 8 点、P_opt） | 単体＋統合 | B1、I1 | M1 |
| VT-05 | 3 | `offline/gaussian`、`LossRegistry`、`pure_loss`、回路雑音 = 損失、OPO・パルス OPO の部分 | 単体 | F3（M0）、B5、B6、B11 | M0–M5 |
| VT-06 | 4 | 線形化のファイバ、OPA g → 0、OPO x = 0、パルス OPO（初期の真空の項） | 単体 | B10、B1、B6、B11 | M1–M5 |
| VT-07–10 | 5 | `fiber_gnlse`（ソリトン、回帰、刻み半減 13–19、Raman ありの光子数） | 単体 | B9 | M5 |
| VT-11 | 6 | `raman_noise_table`（D48: 熱の部分だけ 0） | 単体 | B10 | M5 |
| VT-12 | 6 | `gawbs_table` | 単体 | B10 | M5 |
| VT-13 | 7 | 汎用の線形ループ（M0、`SimController.measure_transfer_function`）、Ref D の z 領域の解析式（M1） | 統合 | F1、F5、B4、I2 | M0、M1 |
| VT-14 | 7 | フィードバックなし → free-running（Ref D） | 統合 | I2、F6 | M1 |
| VT-15 | 7 | 積分器ループの残差、lock 2 の縮小（A−/A+） | 統合 | B4、I1 | M1 |
| VT-16 | 8 | tap のトレードオフ（(a)(b) の 2 つの定義） | 統合 | I1 | M1 |
| VT-17・18 | 9 | Ref B Fig. 3、定義の判別（x = √(P/P_th)） | 単体＋統合 | B6、I4 | M3 |
| VT-19 | 10 | dither の寄与 | 単体＋統合 | F4、B3、I3 | M2 |
| VT-20 | 10 | 全体 約 0.1 rad、制御器のモデルに鈍感 | 統合 | I3 | M2 |
| VT-21 | 10 | ロック時間が遅いアクチュエータのレンジで決まる | 統合 | B4、I2、I3、I5 | M1、M2、M6 |
| VT-22 | 11 | 共有区間の相殺（係数が等しい場合）、偏光をオフにしても通る | 統合 | F6（M0）、I4（M3）、B7 | M0、M3、M4 |
| VT-23 | 11 追加 | 係数が違う OPO 区間では相殺しない | 統合 | B6、I4、B11 | M3、M5 |
| VT-24 | D04 | Ref D のブラインド前の解析的な検証 | 単体＋統合 | I2、B3、B4、B2、B7 | M1 |
| VT-30 | — | OPA | 単体 | B1 | M1 |
| VT-31 | — | CW OPO | 単体 | B6 | M2 |
| VT-32 | — | パルス OPA | 単体 | B11 | M5 |
| VT-33 | — | 位相 → squeezing、ESA、指標の不変性・分解（V11・V12）、高周波の上乗せ（V13） | 単体＋統合 | F4、B5 | M1（V1–V10、V13）、M4（V11・V12） |
| VT-34 | — | レーザー（遅延、コヒーレンス、端数の遅延、Liu 型の静的補償、フィルタキャビティ、ASE、CMRR） | 単体＋統合 | B2、F6 | M1、M4 |
| VT-35 | — | GNLSE の追加 | 単体 | B9 | M5 |
| VT-36 | — | ファイバ量子の追加 | 単体 | B10 | M5 |
| VT-37 | — | 環境雑音（係数と符号、PSD、相関、逆向き、軸の係数） | 単体 | B2、F6 | M1、M2、M4 |
| VT-38 | — | 偏光 | 単体 | B7 | M4 |
| VT-39 | — | probe の相殺の崩れ | 統合 | B8、F6、B7 | M4 |
| VT-40 | — | 振り分け | 単体＋統合 | B8 | M4 |
| VT-41 | — | 損失（Ref B・C の収支、32 ポートの収支、GBS 入力と BHD の分離） | 単体 | F3、I3、I4、I5 | M2–M6 |
| VT-42 | — | 検出 | 単体 | B5 | M1 |
| VT-43 | — | 誤差信号（基底帯域と RF レートの等価性を含む） | 単体 | B3 | M1–M4 |
| VT-44 | — | アクチュエータ、遅延バッファ、積み重ね = 並列、コンパイル = 純 Python、f_sim の収束 | 単体 | F1、B4 | M0、M1 |
| VT-45 | — | パルス OPO | 単体 | B11 | M5 |
| VT-F1（新） | — | YAML ↔ Python の往復（`spec == load(dump(spec))`、同じ seed でビット単位で同じ結果） | 統合 | F2 | M0 |
| VT-F2（新） | — | チャンク長の不変性（G11: ビット単位で同じ）、f_sim の倍化（§6.4） | 統合 | F1、B2 | M0、各 M |
| VT-F3（新） | — | 帯域の台帳の 2 重登録がエラー、損失の重なりの 2 重登録がエラー | 単体 | F3 | M0 |
| VT-F4（新） | — | 記録（config・seed・commit・env）が揃い、記録から同じ結果を再現できる | 統合 | F2 | M0 |

prompt §7 の 1–11 の読み替えは D48・D51 のとおり（§7-1 は ±0.10 dB・集中損失、§7-2 は「前の研究」の破線の右端を除く、§7-6 は熱の部分だけ、§7-8 は VT-16 の (a)(b) の両方、§7-9 は x = √(P/P_th)・Ω_B = 2πf/γ_B・T02 V-B3 の表、§7-11 は係数が等しい共有区間に限り VT-23 を追加）。

### 6.3 数値の許容幅の方針

| 種類 | 許容幅 | 理由 |
|---|---|---|
| 解析的な恒等式（純損失、真空、symplectic、相殺、コンパイル = 純 Python） | 相対 10⁻¹²（float64）。symplectic・パルスは S §4 のとおり 10⁻⁸–10⁻⁹ | 丸め誤差のみ |
| 解析式の数値（VT-01 など） | S §4 の値（±10⁻³ dB など）をそのまま | ノートで導出済み |
| 表の補間 | 10⁻³（VT-32） | 格子の粗さ |
| 確率的な量（PSD、分散、放射計雑音） | 固定 seed で、統計誤差の 3σ が許容幅より小さくなるサンプル数を各テストが計算して使う（Welch の自由度から）。S §4 の ±5–10% はその上で判定 | seed への偶然の依存を避ける |
| 全体のシミュレーション（VT-02、VT-20） | D48 の ±0.10 dB、VT-20 の 10% | D48 |
| テストの時間 | 単体 < 1 s、VT ＜ 60 s（`slow` は除く）、`slow` は夜間・マイルストーンの判定時 | 開発の速さ |

### 6.4 数値の収束の確認（prompt §8）

| 対象 | 方法 | 合格 |
|---|---|---|
| 時間刻み f_sim | 2 倍・4 倍にして残差 σ、R′±、ロック時間を比べる | 変化が許容幅の 1/3 以下 |
| チャンク長 | 2 つのチャンク長で同じ seed | ビット単位で同じ（G11 の設計で保証。崩れたらバグ） |
| seed | 10–30 seed の集団で平均と幅を出す。結論（例: 効果の有無）は幅を超える差だけ | Phase 6 の主張の条件（prompt §8） |
| GNLSE の刻み・格子 | 刻み半減で誤差 1/16（VT-09）、格子の端（`grid_edge_monitor`） | VT-09 |
| パルス OPO の δt | 1/(40–160 κ_c) で 2 段 | T17 V-P1 |
| 小区間の数 N_sub | 2 倍で逆向きの残差の PSD が 2% 以内（最高周波数で f τ/N_sub ≪ 1） | T10 テスト 3 |
| `runner.convergence` | 上をまとめて実行し、表を `runs/` に残す | — |

### 6.5 Ref D のブラインド検証の手順（D04、Phase 5 の終わり）

1. Phase 4–5 では Ref D を VT-13・14・21・24 の解析的な検証だけで確かめる。**作業担当は `blind/` を読まない。テスト・プリセット・ノートから `blind/` を参照しない。**
2. D36 で許された静的な量（ΔL を縞で数えた値、ファイバ長、挿入損失）だけは、ユーザから受け取ったら Phase 5 の終わりより前にプリセットに入れてよい（D36）。それ以外の不明な量（915B の共振、914 の帯域、DSC1 の K_p・K_i、K_i の単位、h0、D_pol、PM の K_T）は掃引の範囲のまま。
3. **凍結**: Phase 5 の VT がすべて通った時点で、調整役が `refD` のプリセットと予測の設定を固定し、git のタグ `refD-blind-freeze` を付ける。予測（不明な量の掃引範囲から出した帯 — 真の φ_D と e_DSC1 の PSD、残差の rms、ロック時間、ロック点のずれ）を `runs/` に作り、結果ファイルのハッシュを `PROCESS_LOG.md` に記録する。
4. 調整役だけが `blind/ref_d_targets.md` を開き、凍結した予測と比べる。**パラメータを調整しない。**
5. ずれは合わせ込まず、アブレーション（トグルを 1 つずつ）と掃引範囲の中の位置で、原因になりうるモデルの部分を特定して記録する。モデルを直す場合は、新しい凍結として別の比較にし、最初の比較の結果は残す。

---

## 7. Phase 4 の作業分割（worktree）

### 7.1 作業パッケージ

「種類」: **F** = フレームワーク（ponytail 禁止。保護リスト）、**B** = ブロックの内部（ブロックの step 関数の本体・数値の手続き・フィルタに限って ponytail 可、AGENTS.md）、**I** = プリセットと統合の検証（ponytail 禁止。検証テストとプリセットは保護リスト）。B の WP でも、ブロックの規約（`Block` のメタ情報、`losses()`、`bands()`、トグル、`SourceOutput`）の部分はフレームワークとして扱う。

| WP | 種類 | 中身 | 依存 | M |
|---|---|---|---|---|
| WP-F0 | F | 骨組み: `pyproject.toml`（G14）、`src/cranberry/` の構成、`Block`・`BlockSpec`・`RunSpec`・`SlotLayout` の型と名前（契約）、`core/fft.py`・`units.py`、テストの土台（マーカー）、`.gitignore` に `runs/`、CLI の入口、numba の wheel の確認 | なし（**最初に単独で**） | M0 |
| WP-F1 | F | `core/*`: スロットの表現、パイプラインの合成と 4 つの検査、コンパイル層と純 Python の参照、遅延バッファ、多重レート、乱数の流れ、ランナーの `run`。**最初に性能の試作**（G01 の条件） | F0 | M0 |
| WP-F2 | F | `config/*`: YAML のタグ、オーバーレイ、制限付きの評価器、依存関係の解決、`to_spec`・`dump_yaml`、記録（`runner/record.py`）、VT-F1・F4 | F0 | M0 |
| WP-F3 | F | `offline/*`: 表のキャッシュ、`SourceOutput`、`LossRegistry`、`BandLedger`、`gaussian.py`、VT-05 の基本部分、VT-F3 | F0 | M0 |
| WP-F5 | F | `hardware/*`、`devices/*.yaml`（出典のページ付き）、`SimController`、写像 | F0（F1 の `run` が入ったら `SimController` を接続） | M0 |
| WP-F6 | F | `path/*`: Topology、区間・小区間・要因の履歴、`segment_phase`、`beam_assembler`（放出時刻、逆向き）、ロックの線形結合、`coef`、VT-22（M0 部分） | F1 | M0 |
| WP-F4 | F | 指標と写像の枠: `metric_i`・`metric_ii`、`phase_to_quadrature`（台帳を 1 回）、`decompose`、`convergence`、出力の形式 | F1、F3、F6 | M0–M1 |
| WP-B1 | B | `single_mode_reference`、`opa_waveguide`（集中・分布・結合波）、`tap`、A± の表、VT-01・04・30 | F3 | M1 |
| WP-B2 | B | `laser_source`、NOISE の全部（状態を持つフィルタ、G11）、`filter_cavity`、`laser_oob_variance`、VT-34・37 | F1、F3、F6 | M1 |
| WP-B3 | B | ERR の全部と `demodulator`・`rf_reference`、基底帯域と `rf_rate_mode`、VT-43 | F1、F3（A± の型） | M1–M3 |
| WP-B4 | B | CTRL・ACT・MOD・TEC の全部（`pid` の 3 つの離散化、`adc`・`dac`、フィルタ、`lock_supervisor`、915B・914・PZT ミラー・EOM・AOM・dither・TEC）、VT-44（ブロック部分）・15 | F1 | M1 |
| WP-B5 | B | DET の全部、ESA の中身（`esa_zero_span`・`esa_swept_spectrum`・`pulse_sample`）、`mode_overlap`、VT-42・33（ESA 部分） | F3、F4 | M1 |
| WP-I1 | I | `refA`・`refA_tap` のプリセット、VT-02・04・15・16 | F*、B1–B5 | M1 |
| WP-I2 | I | `refD` のプリセット（D35・D51-8 の値）、`dsc1` との結合、VT-13・14・21・24 | F*、B2–B5、B7（偏光子・`segment_jones` の最小部分） | M1 |
| WP-B6 | B | `opo_cw_spectrum`、`opo_detuning_noise`・`opo_length_lock`、L(x) の callable、OPO 区間の係数、VT-17・18・23・31 | F3、F6 | M2–M3 |
| WP-B7 | B | 偏光（`segment_jones`、`polarizer_projection`、`pdl_element`、`faraday_nonreciprocal`、`per_budget`、`coil_birefringence`）、VT-38 | F6 | M1（最小）–M4 |
| WP-I3 | I | `refC`、VT-19・20・21・41 | B2–B6 | M2 |
| WP-I4 | I | `refB`、VT-17・18・22・23 | B6、I3 の後 | M3 |
| WP-B8 | B | DMX、Liu 型の具体（逆向きの小区間の検証、`probe_pulser`）、迷光の簡易な表（`stray_light_table` の Rayleigh・反射）、VT-39・40 | F6、B3、B4、B7 | M4 |
| WP-B9 | B | `fiber_gnlse`（RK4IP、項のトグル、CW の近道）、線形・非線形の表、`grid_edge_monitor`、VT-03（GNLSE）・07–10・35 | F3 | M5 |
| WP-B10 | B | 量子雑音層（線形化、TW）、散乱の表（ラマン、GAWBS、ブリルアン、迷光の残り）、`sbs_check`、VT-06（ファイバ）・11・12・36 | F3、B9（古典解の型） | M5 |
| WP-B11 | B | `opa_supermodes`、`pulsed_opo_source`、`opo_angle_map`、`fsr_comb_filter`、`spectral_filter`・`dispersion_comp`、`lo_temporal_shape`、モードの重なりの表、`pulsed_bhd_kernel`、`rep_rate_tail_check`、VT-32・45・03（パルス） | F3、B6 | M5 |
| WP-I5 | I | 将来のプリセット 4 つ、S §2.3 の軸のオーバーレイ、S §2.5 の悪い構成、2 ポート組、指標 (ii) の分解の統合、32 ポートの規模の見積もり | B8–B11、F4 | M4（単一モードの光源で）–M6 |

### 7.2 順序と並行

```
wave 0 : F0 (alone)
wave 1 : F1 ‖ F2 ‖ F3 ‖ F5                                  (F1 starts with the performance spike)
wave 2 : F6 ‖ B1 ‖ B4 ‖ B3        (after F1 API freeze)       then F4 (needs F1, F3, F6)
wave 3 : B2 ‖ B5 ‖ B7(min)  → I1 ‖ I2                         => M1 (Ref A ‖ Ref D)
wave 4 : B6 → I3 (Ref C) → I4 (Ref B)   ‖  B7(rest) ‖ B8        => M2, M3, M4
wave 5 : B9 → B10 ‖ B11       (may start at wave 4 if G17 allows)  => M5
wave 6 : I5                                                    => M6
```

- 同時に動かす worktree は 4–5 本を目安（Phase 1 の 16 本並行はノートで衝突がなかったが、コードは契約の変更が波及する）。
- **API の凍結**: wave 1 の終わりに F1・F2・F3 の公開インターフェースを調整役が確認して凍結する。凍結後の変更は「後から変えると高くつく決定」なので調整役が grilling の対象にする（AGENTS.md）。

### 7.3 各 WP が守る契約

| 契約 | 内容 |
|---|---|
| 規約 | conv の符号（§4）、BHD（§5）、位相ジッタの写像（§6）、片側 PSD（§7）、FFT（§8、`core/fft.py` 以外で `numpy.fft` を直接呼ばない）、区間とロック（§9）、指標（§10）。文献の式を使うときはコメントに「文献の規約 → conv」の変換を書く |
| ブロックの規約 | `make_<kind>(params) -> Block`、`reads`・`writes` を正しく宣言、`kernel` を持つなら純 Python の `__call__` と一致するテストを付ける |
| 損失 | squeezed・LO・probe の経路に損失のある物は必ず `losses()` で登録（G06）。ブロックの中で η を直接掛けない |
| 帯域 | 雑音源は必ず `bands()` で時間ループの帯域とその外を宣言（G05）。帯域外を自分で足さない |
| 乱数 | `core/rng.stream(root, self.name)` だけを使う。`numpy.random` の大域状態を使わない |
| トグル | 雑音源・効果ごとにトグル（S §3.2・ノートの名前）を持ち、オフで恒等（または 0）になる単体テストを付ける（保護リスト） |
| 値の出典 | プリセットの値には出典の印（[Ref A]、[manual p.N]、[AE]、UNVERIFIED、D 番号）をコメントか `source:` で付ける |
| ファイル | 自分の WP のファイルだけを編集・コミットする。`DECISIONS.md`・`notes/conventions.md`・`notes/*`・`PROCESS_LOG.md` は編集しない（変更は報告に「提案」として書く、Phase 1 の運用） |
| 起点 | worktree は main 起点なので、作業前に `git merge --ff-only dev`（PROCESS_LOG の教訓） |
| `blind/` | 読まない、参照しない |
| ponytail | F・I の WP ではオフ（プラグインのフックで毎ターン有効なら、フレームワークに戻るときに通常モードへ）。B の WP でもブロックの規約の部分と保護リストの物はオフ。ponytail-review は B のブロック内部のテストが通った後だけ任意 |
| 完了の報告 | 通した VT の一覧、未解決の物理（C23・C24 など）、契約への変更提案、代用した値（PROCESS_LOG に記録するため） |

### 7.4 Phase 5 のレビューの手順

| 段 | 中身 |
|---|---|
| 担当 | モジュール（WP）ごとに、**新しい文脈の物理レビュー担当**を 1 人（実装担当とは別の agent。`blind/` は読まない） |
| 入力 | そのモジュールのコード、対応するノート（T01–T17）と `conv`、引用された一次資料（`refs/`、D21 の絶対パス）、VT の結果 |
| 点検の表 | (1) 式がノートと一致、(2) 符号（conv §4・§5・§9、機器の符号の変換）、(3) 2 倍の因子（片側 PSD、復調の 2S、R± の定義、x = √(P/P_th)）、(4) 単位（rad/V の読み、dB と線形、Hz と rad/s）、(5) FFT の符号と規格化（conv §8）、(6) 数値の収束（§6.4）、(7) トグルでオフにしたとき恒等、(8) 損失・帯域の登録の漏れと 2 重、(9) 値の出典の印、(10) 保護リストの物を消していない |
| 流れ | 指摘 → 実装担当（同じ worktree で修正）→ 再レビュー。VT がすべて通るまで繰り返す（prompt §4 Phase 5） |
| フレームワークのレビュー | F1・F6 は別に、遅延・放出時刻・逆向き伝搬・積み重ねの等価性を集中して見るレビュー担当を 1 人 |
| 順番 | M1 のモジュール → M2 → M3 → M4 → M5 → M6。各 M の判定の前に、そのマイルストーンの統合（I）の VT を全部通す |
| 任意 | B のブロック内部に ponytail-review（正しさは見ないので物理レビューの代わりにはしない） |
| 終わり | Ref D のブラインド比較（§6.5）。その後、ユーザが信頼できると確認したら Phase 6 |

---

## 8. リスクと未決事項

### 8.1 ユーザだけが答えられる事実（S §7.1）と、それなしでの進め方

| # | 事実 | なしでの進め方 | 影響する WP |
|---|---|---|---|
| 1 | D43 の解釈（Borealis の Supplementary の PPKTP ヘミリシック OPO か） | `borealis_opo` を光源の軸の 1 つとして作る（D43 のとおり）。解釈が違えば軸の値を足す | B11、I5 |
| 2 | 将来の光源の値（波長、f_rep、パルス幅・線幅・T_p、平均パワー、N_b、OPO の台数） | S §7.2 の範囲で掃引。波長は 1545.32 nm で代用（D07） | B11、I5 |
| 3 | Ref D の静的な値（ΔL、腕 2 の長さ、915B のピグテール・挿入損失、曲げ半径） | ΔL 0–15 m（名目 9 m）、挿入損失 0.01–2 dB を掃引。受け取ったら凍結の前に入れる（§6.5-2） | I2 |
| 4 | DSC1 の K_p・K_i・K_d・S・O・CV の制限、K_i の単位 | 利得を掃引、`ki_entry: per_second \| per_sample` を両方 | I2、F5 |
| 5 | 914・オシロスコープの入力インピーダンス | 1 MΩ を仮定、分圧をパラメータに | B4、I2 |
| 6 | Optoquest PC・PM パッチケーブル・Ref A の導波路の 772 nm の損失の資料 | 類似品の典型値（UNVERIFIED）を掃引 | I1、I2 |
| 7 | GBS の要求（許される迷光、時間ビンの分解、毎フレームか、32 以外） | 迷光は「クリック数」と「モード内の光子数」を両方出し、しきい値は出さない（D46） | F4、I5 |
| 8 | Moku:Delta の PID で 10 Vpp が使えるか | ±0.5 V を `moku_delta`、10 Vpp を別のオーバーレイにし、`bad_moku_pid_05V` を残す | F5、I5 |
| 9 | PM15-U25D の K_T を実測してよいか | −(0.3–3.4) rad/(m·K) を掃引（D51-2） | B2、B7、I2 |
| 10 | 将来のハードウェアの値（EOM、BHD、制御器、D15） | メーカー資料・リファレンスの値、`fpga_generic` | I5 |

### 8.2 DECISIONS の Open とその扱い

| Open | 計画での扱い |
|---|---|
| LOCK_TOPOLOGY（D14・D30・D42） | S §2.4 の名前をすべてプリセット（`presets/locks/`）として作り、Phase 6 で比べる |
| ファイバの種類 PM/SMF（D32） | 軸 `fiber: pm \| smf` の両方。G04 の Jones の表現が SMF を可能にする |
| LO 配送のロックの要否（D12） | `L-LOdeliv` をオン・オフのオーバーレイに |
| 未解決の物理 C23（Ref A の lock 1 の SNR と 14 → 9 mrad）、C24（Fig. 3 の極小） | Phase 5 で VT-15 と VT-33 V9 で確認。合わせ込まない |
| Ref D の実測データ | §6.5 |
| D51-10（Ref B の数値を x = 0.745 で作り直す） | WP-B3・B6 の作業に含める（T14 の Ref B の値を VT の目標に使わない） |

### 8.3 技術的なリスク

| リスク | 影響 | 対策 |
|---|---|---|
| numba の性能が AE（50–200 ns/ステップ）に届かない、またはコンパイル時間が長い | 数分のロック時間（Ref C・D、VT-21）が回らない | WP-F1 の最初に試作して計測。届かなければ G02 の表現を調整（カーネルの融合、ブロックの数を減らすのではなく内側の分岐を減らす）。ブロックの数・規約は保護リストなので減らさない |
| numba と numpy・Python 3.14 の版の組 | インストールできない | WP-F0 でインストールを確認。numba 0.68.1 は numpy を 2.5.4 に上げる。だめなら Python 3.13 か numpy の版を下げる（ユーザに確認） |
| 平らな配列のスロットの添字の誤り | 静かな誤計算 | 名前付きビューを常に使い、純 Python の参照との一致テスト（VT-44）、スロットの単位のタグ |
| 32 ポートの規模（小区間 × ポート × 履歴のリング） | メモリ・時間 | M6 で見積もり。履歴の長さは「最長の経路の遅延 ＋ τ_k」で上限が決まる（L_31 ≈ 3.2 km で 16 µs、f_sim 4 MHz で 64 サンプル） |
| 契約の変更の波及（並行の WP） | 衝突、手戻り | wave 1 の終わりの API 凍結、変更は調整役経由 |
| ブラインドの漏れ | D04 の破綻 | `blind/` は gitignore 済み。worktree の担当への指示に毎回入れる。レビューの点検表にも入れる |
| ノートの仮の値の残り（T14 の Ref B、T15 の 1 極の 914、+2.4 rad/(m·K)） | 誤った目標 | D51 の置き換えを WP の指示に明記。VT の目標は S §4 の値（D48・D51 を適用済み）だけを使う |
| Monte Carlo の目的関数の揺らぎ（Phase 6） | BoTorch が不安定 | 線形化（G03）と固定 seed の流れ（G11）。集団の seed で幅を出す |

### 8.4 見つけた食い違い（DECISIONS・conv・notes の間）

| # | 食い違い | 計画での扱い（決定は調整役・ユーザ） |
|---|---|---|
| 1 | Ref D の時間刻み: T15 §8.5 は Δt = 2 µs（DSC1 の更新）、T16 §5-4 は f_sim = 1 MHz（ADC）で PID を 2 サンプルごと | G10 で T16 の案を推奨（2 点平均を含めて DSC1 の離散則をビット単位で再現できる） |
| 2 | 指標 (ii) の定義と 1 ポートの設定: D12 は「常に (i)(ii) の両方を出す」、D31 は (ii) をポート k とポート 0 の差、D13 は「基本は 1 ポートで L_k を掃引」 | G13 で「1 ポートの設定も基準ポート 0 を持つ」を推奨。ポートのない Ref A–D では (ii) は「該当なし」 |
| 3 | T16 §9 の `thermal_actuator` の出力は δψ と ΔT_h の両方、S §3.2 は ΔT_h だけ（C13） | G08 |
| 4 | T15 §8.3 の `dsc1.yaml` は「既定」の語を使う（τ 5.67 µs）、T16 は 5.7 µs | D34 に従い「参照値」と読み、5.7 µs（5.67 は同じ当てはめの丸めの違い）、5.7–8.3 µs を掃引 |
| 5 | T08 §7.3 は機械雑音を「FFT で色付けしてチャンクごとに生成」、prompt §8 はチャンク長に対する収束の確認を要求 | G11 で状態を持つフィルタを推奨（チャンク長に厳密に不変） |
| 6 | prompt §5.4 は (a)(b)(c) から「1 つ」を選ぶ書き方、T16 と本計画は (c) ＋ (a) の組み合わせ | (c) が閉ループの方式で、(a) は MHz の搬送波の表し方なので両立する、と G01 に明記。prompt の意図に沿うかをゲートで確認 |
| 7 | `.gitignore` の `refs` は全階層に効く | `presets/refs/` のような名前を使わない（§2.9） |
| 8 | T16 §11-6 は「ブロック内部をコンパイル可能な形に制限すること」を Phase 3 の決定事項に挙げる。AGENTS.md の ponytail の規則はブロック内部の最小化を許す | 両立する: ブロック内部は numba で書ける範囲で最小にし、純 Python の参照と一致させる（G01・G02） |
| 9 | DECISIONS D49 の「Ref D では T07 の散乱雑音をオフ」と S §2.1 の `refD` の「オフにする opt」の行は一致。ただし S §0.2 の定義では T07 の散乱（E32 など）は kept | 矛盾ではない（kept = 少なくとも 1 つのプリセットでオン）。記録のみ |

---

## 9. 承認後に DECISIONS に入れるもの（提案。このファイルからは編集しない）

1. G01–G17 の決定（選ばれた案と理由）。
2. VT-F1–F4（フレームワークの性質のテスト）を検証テストの一覧に加える。
3. §8.4 の食い違い 1・2・5・6 の解決。
4. マイルストーン M0–M6 と WP の分割、API の凍結の時点。
5. G14 で追加する依存関係と版（`pyproject.toml`・`uv.lock` は WP-F0 が更新）。
