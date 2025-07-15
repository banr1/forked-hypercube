# Audit Report

実装を追っていくと、以下のような誤り・不厳密な点が見つかりました。優先度順に並べていますので、まずは①②あたりを直すと動作がおおむね正しくなるはずです。

---

### ① `compute_checksum_cost()` の誤り：合計ではなく最後の余りしか返していない

```python
def compute_checksum_cost(sum_value):
    total_cost = 0
    while sum_value != 0:
        global w_value
        total_cost = sum_value % w_value     # ← ここが “=”，本来は累積すべき
        sum_value //= w_value
    return total_cost
```

* **問題**：`total_cost = sum_value % w_value` と書いているため、ループを回すごとに上書きされ、最終的に「最後の一桁分の余り」しか返りません。
* **修正例**：正しくは各桁（`sum_value % w_value`）を足し合わせるべきです。

  ```diff
   def compute_checksum_cost(sum_value):
       total_cost = 0
       global w_value
       while sum_value != 0:
  -        total_cost = sum_value % w_value
  +        total_cost += sum_value % w_value
           sum_value //= w_value
       return total_cost
  ```

---

### ② 二重定義されたチェックサム関数の混在

* 上記の `compute_checksum_cost`（誤）に加えて、末尾に正しく動く

  ```python
  def checksum_cost(layer, max_layer, w):
      inverted = max_layer - layer
      cost = 0
      while inverted > 0:
          cost += inverted % w
          inverted //= w
      return cost
  ```

  が定義されています。
* **問題**：名前が違うだけで「同じ役割」の関数がふたつあり、コード中の用途によって使い分けられているため非常に分かりにくいです。
* **対策**：

  * 「一貫して１つのチェックサム関数」を使う
  * `compute_extrema` 内も `checksum_cost(...)` のほうを呼ぶように統一する

---

### ③ `binom(n, k)` の境界チェックがない

```python
def binom(n, k):
    if binoms[n][k] == 0:
        value = factorials[n] // (factorials[k] * factorials[n - k])
        binoms[n][k] = value
        return value
    else:
        return binoms[n][k]
```

* **問題**：`k < 0` や `k > n` のときに本来は 0 を返すべきですが、`factorials[n - k]` のインデックスが負になり、Python のリストの負インデックス参照につながる恐れがあります。
* **修正例**：

  ```python
  def binom(n, k):
      if k < 0 or k > n:
          return 0
      if binoms[n][k] == 0:
          value = factorials[n] // (factorials[k] * factorials[n - k])
          binoms[n][k] = value
      return binoms[n][k]
  ```

---

### ④ `print_optimal_L2_values()` での `best_w` 更新ミス

```python
if c < best_cost - 0.5:
    best_cost = c
    best_w = min_w    # ← 本来は `w` を入れたい
    print("w=", w, ...)
```

* **問題**：最良の `w` を記録するときに、ループ下限の `min_w` を代入してしまっている。
* **修正例**：

  ```diff
   if c < best_cost - 0.5:
       best_cost = c
  -    best_w = min_w
  +    best_w = w
       print("w=", w, ...)
  ```

---

### ⑤ `compute_max_entropy_lower_bound()` の未使用変数・命名混乱

```python
mu_d0_scaled = (powcc - mu_scaled * sum_ld) / layer_sizes[d + 1]
# … この変数はその後一度も使われていない
```

* **問題**：計算してはいるものの、何も返り値に含めていません。また、引数名が `mu_scaled` ですが、実際には「スケール前の µ」を渡しているケースもあり、混乱を招きます。
* **対策**：

  * 未使用の `mu_d0_scaled` は削除
  * 引数名を `mu` に直すか、ドキュメントで「ここでは未スケールの µ」を明記

---

### ⑥ 浮動小数点による対数・平方根計算の桁落ち・オーバーフロー

* `compute_extrema()` 内で

  ```python
  v1 = (sum_ld_cd * sum_ld_cd) / (sum_ld * sum_ld)
  v4 = v2 / v3
  discr = v1 - v4
  lambda_21 = sum_ld_cd / sum_ld + sqrt(discr)
  ```

  のように、大きな整数を `float` に変換して計算しています。
* **問題**：`w**v` が非常に大きな値になると、`float` の指数部上限（約10^308）を超えて `inf` になったり、精度が落ちたりします。
* **対策**：

  * Python の `fractions.Fraction` や `decimal.Decimal`、あるいは有理数そのままで計算する
  * 対数は `math.log2()`／`int.bit_length()` など整数向きの関数を使う

---

### ⑦ ビルトイン名 `sum` を変数名に使用

```python
def nb(...):
    sum = 0
    ...
    return sum
```

* **問題**：Python 標準の `sum()` 関数を隠蔽してしまいます。
* **対策**：`total = 0` など、別の変数名に変える

---

### ⑧ `global w_value` の配置が不自然

```python
def compute_checksum_cost(sum_value):
    total_cost = 0
    while sum_value != 0:
        global w_value   # ループ内に書かれている
        total_cost += ...
        ...
```

* **問題**：`global` 宣言は関数先頭に置くのが慣例です。ループ内にあると可読性を下げます。
* **対策**：関数冒頭に移動

  ```python
  def compute_checksum_cost(sum_value):
      global w_value
      total_cost = 0
      ...
  ```

---

### ⑨ 毎回 `precompute_real()` を呼び直す非効率

* 各種探索ループの中で、ハイパーキューブの階層サイズや階乗テーブルを何度も再計算しています。
* **対策**：

  * `chains`（次元）や `val`（基数）が同じならキャッシュ
  * あるいは、外側で一度だけ呼び出してから各関数に必要な配列を渡すようにする

---

以上の修正を加えることで、正確性と可読性、パフォーマンスが大きく改善されると思います。もしさらに具体的な箇所の動作例やテストの追加が必要であれば、お知らせください。
