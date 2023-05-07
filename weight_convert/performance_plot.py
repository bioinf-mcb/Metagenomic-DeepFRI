import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

df = pd.read_csv('inference_times.csv.gz')
df = df.assign(speedup=lambda x: x['tf_time'] / x['onnx_time'])

fig, axes = plt.subplots(1, 4, figsize=(10, 5), sharey=True, sharex=True)

for mode, axis in zip(df["mode"].unique(), axes.flatten()):
    data = df.query(f"mode == \"{mode}\"")
    data = data.sort_values("prot_len")

    # round protein length to highest 25
    data = data.assign(prot_len=lambda x: np.ceil(x['prot_len'] / 50) * 50)
    sns.lineplot(x="prot_len", y="speedup", hue="net_type", data=data, ax=axis)
    axis.set_xlabel("Protein length")
    axis.set_title(f"Mode {mode}")

axes[0].set_ylabel("Speedup ONNX vs TF2 (times)")
fig.tight_layout()
fig.savefig("onnx_vs_tf2.png", dpi=150)
