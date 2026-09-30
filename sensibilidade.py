# Sensibilidade à calibração do ruído do GNSS: o filtro assume 2,5 m; o ruído real varia.
import json, numpy as np
from simulacao_ekf import run_once, truth, N_RUNS, METHODS
tr = truth(); out = {}
for sg in (2.5, 4.0, 6.0):
    runs = [run_once(1000 + i, "nominal", tr, sig_true=sg)[0] for i in range(N_RUNS)]
    out[str(sg)] = {}
    for m in ("GNSS isolado", "EKF sem validação", "EKF com validação (NIS)"):
        out[str(sg)][m] = {k: float(np.mean([r[m][k] for r in runs])) for k in runs[0][m]}
        out[str(sg)][m].update({k + "_sd": float(np.std([r[m][k] for r in runs])) for k in ("rmse",)})
json.dump(out, open("sensibilidade.json", "w"), indent=2, ensure_ascii=False)
for sg, v in out.items():
    print(sg, {m: (round(x["rmse"], 2), round(x.get("cobertura95", 0), 1), round(x.get("rej_pct", 0), 1)) for m, x in v.items()})
