"""
糸球体計数データ（自前測定、4系統×雌雄5匹＝40個体、8週齢）の解析と検証。
ユーザー報告値の再現確認を兼ねる。
"""
import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols

BASE = "/usr/local/jupyter/FVB_B6_glo"
OUT = f"{BASE}/results/phenotype"

df = pd.read_csv(f"{BASE}/data/glomeruli_4strain.csv")
df["glom_per_g_kw"] = df["glomeruli_per_kidney"] / df["kidney_weight_g"]
df["glom_per_g_bw"] = df["glomeruli_per_kidney"] / df["body_weight_g"]
df["kw_bw_pct"] = df["kidney_weight_g"] / df["body_weight_g"] * 100

ORDER = ["BALB", "DBA", "B6", "FVB"]
df["strain"] = pd.Categorical(df["strain"], categories=ORDER, ordered=True)

print("=" * 72)
print("1. 系統別要約（雌雄プール, n=10）— ユーザー報告値との照合")
print("=" * 72)
rep = {"BALB": (12548, 477, 68592), "DBA": (8311, 406, 61792),
       "B6": (7249, 373, 71410), "FVB": (6104, 253, 36071)}
print(f"{'系統':6s} {'糸球体数':>10s} {'(報告)':>9s} {'/g体重':>8s} {'(報告)':>7s} {'/g腎重':>9s} {'(報告)':>8s}")
for s in ORDER:
    d = df[df.strain == s]
    g, bw, kw = d.glomeruli_per_kidney.mean(), d.glom_per_g_bw.mean(), d.glom_per_g_kw.mean()
    r = rep[s]
    print(f"{s:6s} {g:10.0f} {r[0]:>9d} {bw:8.0f} {r[1]:>7d} {kw:9.0f} {r[2]:>8d}"
          f"   {'OK' if abs(g-r[0])<1 and abs(bw-r[1])<1 and abs(kw-r[2])<1 else '★不一致'}")

print("\n" + "=" * 72)
print("2. FVB vs B6（雄, n=5）— ユーザー報告値との照合")
print("=" * 72)
m = df[df.sex == "M"]
comparisons = [
    ("絶対数", "glomeruli_per_kidney", 5888, 6818, -13.6, 0.028),
    ("/g腎重量", "glom_per_g_kw", 29419, 57442, -48.8, 1.2e-5),
    ("/g体重", "glom_per_g_bw", 217, 303, -28.3, 0.0014),
    ("KW/BW比(%)", "kw_bw_pct", 0.742, 0.528, 40.4, 4.0e-5),
]
for label, col, rf, rb, rpct, rp in comparisons:
    f = m[m.strain == "FVB"][col]; b = m[m.strain == "B6"][col]
    t, p = stats.ttest_ind(f, b)
    pct = (f.mean() - b.mean()) / b.mean() * 100
    print(f"{label:12s} FVB={f.mean():9.3f}(報告{rf:>8}) B6={b.mean():9.3f}(報告{rb:>8}) "
          f"{pct:+6.1f}%(報告{rpct:+.1f}%) p={p:.2e}(報告{rp:.1e})")

print("\n" + "=" * 72)
print("3. 内部対照: FVB vs BALB（雄）— 回収率バイアスの排除")
print("=" * 72)
for label, col in [("腎重量(g)", "kidney_weight_g"), ("KW/BW比(%)", "kw_bw_pct"),
                   ("糸球体数", "glomeruli_per_kidney")]:
    f = m[m.strain == "FVB"][col]; b = m[m.strain == "BALB"][col]
    t, p = stats.ttest_ind(f, b)
    print(f"{label:12s} FVB={f.mean():8.3f}  BALB={b.mean():8.3f}  p={p:.4f}"
          f"   {'← 区別できない' if p > 0.05 else '← 有意差あり'}")

print("\n" + "=" * 72)
print("4. 二元配置ANOVA（系統×性）: 糸球体数")
print("=" * 72)
model = ols("glomeruli_per_kidney ~ C(strain) * C(sex)", data=df).fit()
aov = sm.stats.anova_lm(model, typ=2)
print(aov.to_string())
print("\n報告値: 系統 F=48.9 p=4.7e-12 / 性 F=6.3 p=0.018 / 交互作用 F=1.81 p=0.166")

print("\n" + "=" * 72)
print("5. ANCOVA（腎重量を共変量）: 性の主効果は消えるか")
print("=" * 72)
model2 = ols("glomeruli_per_kidney ~ C(strain) + C(sex) + kidney_weight_g", data=df).fit()
aov2 = sm.stats.anova_lm(model2, typ=2)
print(aov2.to_string())
print("\n報告値: 性の主効果 p=0.456")

print("\n" + "=" * 72)
print("6. 系統間の全ペア比較（Tukey HSD, 糸球体数, 雌雄プール）")
print("=" * 72)
from statsmodels.stats.multicomp import pairwise_tukeyhsd
tk = pairwise_tukeyhsd(df["glomeruli_per_kidney"], df["strain"].astype(str))
print(tk)

print("\n" + "=" * 72)
print("7. 表現型順位（rank concordance解析で使用した順位の根拠）")
print("=" * 72)
summ = df.groupby("strain", observed=True).agg(
    n=("glomeruli_per_kidney", "size"),
    glom_mean=("glomeruli_per_kidney", "mean"), glom_sd=("glomeruli_per_kidney", "std"),
    bw_mean=("glom_per_g_bw", "mean"), bw_sd=("glom_per_g_bw", "std"),
    kw_mean=("glom_per_g_kw", "mean"), kw_sd=("glom_per_g_kw", "std"),
).round(1)
print(summ.to_string())
summ.to_csv(f"{OUT}/summary_by_strain.csv")
df.to_csv(f"{OUT}/phenotype_with_derived.csv", index=False)
print(f"\n出力: {OUT}/")
