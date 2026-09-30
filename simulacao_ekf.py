#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulação de fusão de sensores (IMU + odometria + GNSS) com Filtro de Kalman Estendido (EKF)
e monitor de consistência (NIS) para navegação planar 2D.

Reprodutibilidade: sementes fixas. Requer: numpy, matplotlib.
Uso: python simulacao_ekf.py
Saídas: resultados.json, fig_trajetoria.png, fig_erro_tempo.png, fig_rmse.png
"""
import json
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------- Parâmetros gerais -------------------------
DT = 0.05            # passo do IMU/EKF (20 Hz)
T_TOTAL = 300.0      # duração (s)
N = int(T_TOTAL / DT)
GPS_DT = 1.0         # GNSS a 1 Hz
WHEEL_DT = 0.1       # odometria a 10 Hz
SIG_GPS = 2.5        # desvio-padrão do GNSS por eixo (m)
SIG_WHEEL = 0.15     # desvio-padrão da velocidade da roda (m/s)
SIG_GYRO = 0.01      # ruído do giroscópio (rad/s)
SIG_ACC = 0.10       # ruído do acelerômetro (m/s^2)
SIG_BIAS0 = 0.005    # desvio-padrão do viés inicial do giroscópio (rad/s)
SIG_SCALE = 0.005    # erro de fator de escala da odometria (0,5 %)
GATE = 9.21          # limiar qui-quadrado (2 gl, 99 %)
R95_LIM = 5.0        # limiar do raio de incerteza 95 % para modo degradado (m)
N_RUNS = 100

SCENARIOS = {
    "nominal": dict(outage=None, outlier_prob=0.0),
    "queda_gnss": dict(outage=(100.0, 140.0), outlier_prob=0.0),
    "outliers_gnss": dict(outage=None, outlier_prob=0.08),
}
METHODS = ["GNSS isolado", "Estima por odometria", "EKF sem validação", "EKF com validação (NIS)"]


def truth():
    t = np.arange(N) * DT
    v = 8.0 + 1.5 * np.sin(2 * np.pi * t / 60.0)
    w = 0.12 * np.sin(2 * np.pi * t / 45.0) + 0.05 * np.sin(2 * np.pi * t / 17.0)
    a = np.gradient(v, DT)
    psi = np.cumsum(w) * DT
    x = np.cumsum(v * np.cos(psi)) * DT
    y = np.cumsum(v * np.sin(psi)) * DT
    return t, x, y, psi, v, w, a


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def run_once(seed, scen, tr, record=False, sig_true=None):
    rng = np.random.default_rng(seed)
    sg = SIG_GPS if sig_true is None else sig_true
    t, X, Y, PSI, V, W, A = tr
    cfg = SCENARIOS[scen]
    bg = rng.normal(0, SIG_BIAS0)
    scale = 1.0 + rng.normal(0, SIG_SCALE)

    # Medições
    wm = W + bg + rng.normal(0, SIG_GYRO, N)
    am = A + rng.normal(0, SIG_ACC, N)
    k_gps = int(GPS_DT / DT)
    k_whl = int(WHEEL_DT / DT)
    gps = {}
    for k in range(0, N, k_gps):
        tk = t[k]
        if cfg["outage"] and cfg["outage"][0] <= tk < cfg["outage"][1]:
            continue
        z = np.array([X[k], Y[k]]) + rng.normal(0, sg, 2)
        # o primeiro fix é considerado validado (partida em condição nominal)
        if len(gps) > 0 and rng.random() < cfg["outlier_prob"]:
            ang = rng.uniform(0, 2 * np.pi)
            mag = rng.uniform(25, 40)
            z = z + mag * np.array([np.cos(ang), np.sin(ang)])
        gps[k] = z
    whl = {k: scale * V[k] + rng.normal(0, SIG_WHEEL) for k in range(0, N, k_whl)}

    k0 = min(gps.keys())
    p0 = gps[k0]
    psi0 = PSI[k0] + rng.normal(0, 0.1)
    v0 = whl[min(kk for kk in whl if kk >= k0)]

    # ---- GNSS isolado (mantém última posição) ----
    est_gps = np.full((N, 2), np.nan)
    last = None
    for k in range(N):
        if k in gps:
            last = gps[k]
        if last is not None:
            est_gps[k] = last

    # ---- Estima por odometria (giroscópio + roda), inicializada no 1º fix ----
    est_dr = np.full((N, 2), np.nan)
    px, py, ps, vv = p0[0], p0[1], psi0, v0
    for k in range(k0, N):
        if k in whl:
            vv = whl[k]
        ps += wm[k] * DT
        px += vv * np.cos(ps) * DT
        py += vv * np.sin(ps) * DT
        est_dr[k] = (px, py)

    # ---- EKF ----
    out = {}
    for name, gated in (("EKF sem validação", False), ("EKF com validação (NIS)", True)):
        x = np.array([p0[0], p0[1], psi0, v0, 0.0])
        P = np.diag([SIG_GPS**2, SIG_GPS**2, 0.1**2, 0.5**2, 0.01**2])
        Qd = np.diag([1e-4, 1e-4, (SIG_GYRO * DT) ** 2, (SIG_ACC * DT) ** 2, (1e-4) ** 2 * DT])
        Rg = np.eye(2) * SIG_GPS**2
        Rw = SIG_WHEEL**2
        Hg = np.zeros((2, 5)); Hg[0, 0] = 1; Hg[1, 1] = 1
        Hw = np.zeros((1, 5)); Hw[0, 3] = 1
        est = np.full((N, 2), np.nan)
        r95 = np.full(N, np.nan)
        rejected = 0
        consec = 0
        n_upd = 0
        t_acc = 0.0
        for k in range(k0, N):
            tic = time.perf_counter()
            # predição
            px, py, ps, vv, b = x
            wc = wm[k] - b
            c, s = np.cos(ps), np.sin(ps)
            x = np.array([px + vv * c * DT, py + vv * s * DT, ps + wc * DT, vv + am[k] * DT, b])
            F = np.eye(5)
            F[0, 2] = -vv * s * DT; F[0, 3] = c * DT
            F[1, 2] = vv * c * DT; F[1, 3] = s * DT
            F[2, 4] = -DT
            P = F @ P @ F.T + Qd
            # atualização odometria
            if k in whl:
                y = whl[k] - x[3]
                S = P[3, 3] + Rw
                K = (P @ Hw.T) / S
                x = x + (K * y).ravel()
                I_KH = np.eye(5) - K @ Hw
                P = I_KH @ P @ I_KH.T + K @ K.T * Rw
            # atualização GNSS
            if k in gps:
                y = gps[k] - x[:2]
                S = Hg @ P @ Hg.T + Rg
                nis = float(y @ np.linalg.solve(S, y))
                n_upd += 1
                if gated and nis > GATE and consec < 4:
                    rejected += 1
                    consec += 1
                elif gated and nis > GATE:
                    # 5 rejeições consecutivas: reinicializa a posição (recuperação de lock-out)
                    x[:2] = gps[k]
                    P[0:2, :] = 0.0; P[:, 0:2] = 0.0
                    P[0, 0] = P[1, 1] = 4 * SIG_GPS**2
                    consec = 0
                else:
                    consec = 0
                    K = P @ Hg.T @ np.linalg.inv(S)
                    x = x + K @ y
                    I_KH = np.eye(5) - K @ Hg
                    P = I_KH @ P @ I_KH.T + K @ Rg @ K.T
            x[2] = wrap(x[2])
            t_acc += time.perf_counter() - tic
            est[k] = x[:2]
            r95[k] = np.sqrt(5.991 * np.max(np.linalg.eigvalsh(P[:2, :2])))
        out[name] = dict(est=est, r95=r95, rejected=rejected, n_upd=n_upd, t_cycle_us=1e6 * t_acc / (N - k0))
    ests = {"GNSS isolado": est_gps, "Estima por odometria": est_dr,
            "EKF sem validação": out["EKF sem validação"]["est"],
            "EKF com validação (NIS)": out["EKF com validação (NIS)"]["est"]}

    truth_xy = np.stack([X, Y], axis=1)
    res = {}
    for m in METHODS:
        e = np.linalg.norm(ests[m][k0:] - truth_xy[k0:], axis=1)
        d = np.linalg.norm(np.diff(ests[m][k0:], axis=0), axis=1).sum()
        d_true = np.linalg.norm(np.diff(truth_xy[k0:], axis=0), axis=1).sum()
        res[m] = dict(rmse=float(np.sqrt(np.mean(e**2))), emax=float(e.max()),
                      path_err_pct=float(100 * (d - d_true) / d_true))
    for m in ("EKF sem validação", "EKF com validação (NIS)"):
        e = np.linalg.norm(ests[m][k0:] - truth_xy[k0:], axis=1)
        r = out[m]["r95"][k0:]
        res[m]["cobertura95"] = float(100 * np.mean(e <= r))
        flag = r > R95_LIM
        res[m]["pct_degradado"] = float(100 * flag.mean())
        res[m]["emax_sem_alerta"] = float(e[~flag].max()) if (~flag).any() else 0.0
        res[m]["t_cycle_us"] = float(out[m]["t_cycle_us"])
        res[m]["rej_pct"] = float(100 * out[m]["rejected"] / max(out[m]["n_upd"], 1))
        if cfg["outage"]:
            idx = np.where(flag & (t[k0:] >= cfg["outage"][0]))[0]
            res[m]["latencia_alerta_s"] = float(t[k0:][idx[0]] - cfg["outage"][0]) if len(idx) else float("nan")
    payload = None
    if record:
        payload = dict(t=t, truth=truth_xy, gps=np.array([[k, *gps[k]] for k in sorted(gps)]),
                       ests=ests, r95=out["EKF com validação (NIS)"]["r95"], k0=k0)
    return res, payload


def main():
    tr = truth()
    agg = {}
    rec = {}
    for scen in SCENARIOS:
        runs = []
        for i in range(N_RUNS):
            res, p = run_once(1000 + i, scen, tr, record=(i == 0))
            runs.append(res)
            if i == 0:
                rec[scen] = p
        agg[scen] = {}
        for m in METHODS:
            agg[scen][m] = {}
            for key in runs[0][m]:
                vals = np.array([r[m][key] for r in runs], dtype=float)
                agg[scen][m][key] = dict(mean=float(np.nanmean(vals)), std=float(np.nanstd(vals)),
                                         median=float(np.nanmedian(vals)))
        print("cenário", scen, "ok")
    with open("resultados.json", "w") as f:
        json.dump(agg, f, indent=2, ensure_ascii=False)

    # ---------------- Figuras ----------------
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9})
    p = rec["queda_gnss"]
    fig, ax = plt.subplots(figsize=(6.2, 4.6))
    ax.plot(p["truth"][:, 0], p["truth"][:, 1], "k-", lw=1.6, label="Trajetória verdadeira")
    ax.plot(p["gps"][:, 1], p["gps"][:, 2], ".", ms=3, color="#999999", label="Fixes GNSS")
    ax.plot(p["ests"]["Estima por odometria"][p["k0"]:, 0], p["ests"]["Estima por odometria"][p["k0"]:, 1],
            "--", color="#d95f02", lw=1.1, label="Estima por odometria")
    ax.plot(p["ests"]["EKF com validação (NIS)"][p["k0"]:, 0], p["ests"]["EKF com validação (NIS)"][p["k0"]:, 1],
            "-", color="#1b9e77", lw=1.1, label="EKF com validação")
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_aspect("equal", adjustable="datalim")
    ax.grid(alpha=0.3); ax.legend(loc="best", fontsize=8, frameon=True)
    fig.tight_layout(); fig.savefig("fig_trajetoria.png", dpi=200); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 3.9))
    t = p["t"]; k0 = p["k0"]
    for m, col, ls in (("GNSS isolado", "#7570b3", "-"), ("Estima por odometria", "#d95f02", "--"),
                       ("EKF com validação (NIS)", "#1b9e77", "-")):
        e = np.linalg.norm(p["ests"][m][k0:] - p["truth"][k0:], axis=1)
        ax.plot(t[k0:], e, ls, color=col, lw=1.0, label=m)
    ax.plot(t[k0:], p["r95"][k0:], ":", color="k", lw=1.0, label="Raio de incerteza 95 % (EKF)")
    ax.axvspan(100, 140, color="#dddddd", alpha=0.6, label="Queda do GNSS")
    ax.set_xlabel("Tempo (s)"); ax.set_ylabel("Erro de posição (m)")
    ax.set_ylim(0, 40); ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.2))
    fig.tight_layout(); fig.savefig("fig_erro_tempo.png", dpi=200); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    scen_names = {"nominal": "Nominal", "queda_gnss": "Queda do GNSS (40 s)", "outliers_gnss": "Outliers no GNSS (8 %)"}
    cols = ["#7570b3", "#d95f02", "#66a61e", "#1b9e77"]
    wbar = 0.2
    for j, m in enumerate(METHODS):
        vals = [agg[s][m]["rmse"]["mean"] for s in SCENARIOS]
        errs = [agg[s][m]["rmse"]["std"] for s in SCENARIOS]
        ax.bar(np.arange(3) + (j - 1.5) * wbar, vals, wbar, yerr=errs, color=cols[j], label=m, capsize=2)
    ax.set_xticks(range(3)); ax.set_xticklabels([scen_names[s] for s in SCENARIOS])
    ax.set_ylabel("RMSE de posição (m)"); ax.set_yscale("log"); ax.grid(alpha=0.3, axis="y")
    ax.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=4)
    fig.tight_layout(); fig.savefig("fig_rmse.png", dpi=200); plt.close(fig)
    print("figuras ok")


if __name__ == "__main__":
    main()
