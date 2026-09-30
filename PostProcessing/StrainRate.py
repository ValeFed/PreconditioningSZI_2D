import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import config as cfg
from scipy.ndimage import gaussian_filter1d

def analyze_baseline_trend(
    time,
    baseline,
    threshold_ratio=0.3,
    min_below_points=1,
    return_diagnostics=False,
    strict_split=False,
    split_mask=None,
):
    """
    Identify the time when the rapid growth phase of the baseline ends, and determine if the trend after that is increasing, decreasing, or stable.
    
    The split time is calculated as the intersection point between the processed derivative
    and the threshold line, using linear interpolation.
    """


    time = np.array(time)
    baseline = np.array(baseline)
    split_mask = None if split_mask is None else np.asarray(split_mask, dtype=bool)
    if split_mask is not None and split_mask.shape != time.shape:
        raise ValueError("split_mask must have the same shape as time/baseline")

    # Seleziona gli indici dove la baseline supera il 5% del massimo (esclude il tratto iniziale quasi piatto).
    valid = np.where(baseline > np.nanmax(baseline) * 0.05)[0]
    # Se ci sono meno di 5 punti validi, non c'e informazione sufficiente per identificare uno split affidabile.
    if len(valid) < 5:
        if return_diagnostics:
            return None, "stable", None
        return None, "stable"

    i_start = valid[0]
    i_end = valid[-1] - 2
    # Taglia il tempo e la baseline sul tratto utile.
    time_crop = time[i_start:i_end]
    baseline_crop = baseline[i_start:i_end]

    # Crea una maschera che mantiene solo i punti finiti simultaneamente in tempo e baseline.
    finite_mask = np.isfinite(time_crop) & np.isfinite(baseline_crop)
    time_crop = time_crop[finite_mask]
    baseline_crop = baseline_crop[finite_mask]
    split_mask_crop = split_mask[i_start:i_end][finite_mask] if split_mask is not None else None
    # Se dopo la pulizia restano meno di 5 punti, interrompe l'analisi.
    if len(time_crop) < 5:
        if return_diagnostics:
            return None, "stable", None
        return None, "stable"

    # Costruisce una maschera per mantenere solo campioni a tempo strettamente crescente (evita dt=0 nel gradiente).
    increasing_mask = np.concatenate(([True], np.diff(time_crop) > 0))
    time_crop = time_crop[increasing_mask]
    baseline_crop = baseline_crop[increasing_mask]
    split_mask_crop = split_mask_crop[increasing_mask] if split_mask_crop is not None else None
    # Ricontrolla che i punti residui siano sufficienti dopo il filtraggio
    if len(time_crop) < 5:
        if return_diagnostics:
            return None, "stable", None
        return None, "stable"

    time_search = time_crop[split_mask_crop] if split_mask_crop is not None else time_crop
    baseline_search = baseline_crop[split_mask_crop] if split_mask_crop is not None else baseline_crop
    if len(time_search) < 5:
        if return_diagnostics:
            return None, "stable", None
        return None, "stable"
    time_diag = time_search
    baseline_diag = baseline_search

    # Calcola la derivata numerica della baseline rispetto al tempo.
    grad = np.gradient(baseline_search, time_search)
    # Il gradiente viene usato senza smoothing.
    grad_proc = grad.copy()

    # Estrae il valore massimo della derivata processata (picco di crescita).
    peak_grad = np.max(grad_proc)
    peak_idx = int(np.argmax(grad_proc))
    # Prende il tratto della derivata dopo il picco per stimare il livello di fondo post-picco.
    grad_post_peak = grad_proc[peak_idx + 1 :]
    # Calcola la media post-picco (primi 5 Myr); se il tratto è vuoto usa peak_grad come fallback.
    mean_post_peak_raw = np.mean(grad_post_peak[3:13]) if grad_post_peak.size > 0 else peak_grad
    # Se la media post-picco è positiva usa quella come riferimento,
    # altrimenti usa 0 come richiesto.
    mean_post_peak = max(mean_post_peak_raw, 0.0)
    threshold_reference = mean_post_peak if mean_post_peak_raw > 0 else 0.0
    # Definisce la threshold come frazione della differenza picco-riferimento.
    threshold = threshold_reference + threshold_ratio * (peak_grad - threshold_reference)
    if threshold_reference == 0.0: threshold = 0.0
    #print(f"Peak gradient: {peak_grad:.3e}, Mean post-peak: {mean_post_peak:.3e}, Threshold reference: {threshold_reference:.3e}, Threshold: {threshold:.3e}")
    # Inizializza l'indice di split (sarà il campione destro del segmento che attraversa la soglia).
    split_idx = None
    # Inizializza il tempo di split (sarà l'intersezione interpolata).
    split_time = None
    # Flag che indica se e stato trovato un crossing reale gradiente-threshold. <- da eliminare
    crossing_found = False

    # Avvia la ricerca dal campione successivo al picco per cercare crossing solo nel ramo discendente.
    search_start = max(peak_idx + 1, 1)
    # Scorre tutti i campioni successivi al picco.
    for i in range(search_start, len(grad_proc)):
        # Verifica il crossing dall'alto verso il basso tra due campioni consecutivi.
        if grad_proc[i - 1] >= threshold and grad_proc[i] <= threshold:
            # Tempo del punto sinistro del segmento di crossing.
            x1, y1 = time_search[i - 1], grad_proc[i - 1]
            # Tempo del punto destro del segmento di crossing.
            x2, y2 = time_search[i], grad_proc[i]
            # Se il segmento e ben definito, calcola l'intersezione con interpolazione lineare.
            if abs(y2 - y1) > 1e-15 and abs(x2 - x1) > 1e-12:
                # Pendenza del segmento locale del gradiente.
                m = (y2 - y1) / (x2 - x1)
                # Tempo di intersezione tra retta locale del gradiente e threshold.
                split_time = x1 + (threshold - y1) / m
            else:
                # In caso degenerato, usa il tempo sinistro come fallback locale.
                split_time = x1
            # Salva l'indice del campione destro del crossing trovato.
            split_idx = i
            crossing_found = True
            # Interrompe al primo crossing valido.
            break

    # Se non e stato trovato alcun crossing reale, attiva il fallback esplicito.
    if split_time is None:
        if strict_split:
            if return_diagnostics:
                diagnostics = {
                    "time_crop": time_diag,
                    "baseline_crop": baseline_diag,
                    "grad": grad,
                    "grad_proc": grad_proc,
                    "threshold": float(threshold),
                    "peak_idx": int(peak_idx),
                    "peak_grad": float(peak_grad),
                    "mean_post_peak": float(mean_post_peak),
                    "threshold_reference": float(threshold_reference),
                    "split_time": None,
                    "split_grad": None,
                    "crossing_found": bool(crossing_found),
                }
                return None, "stable", diagnostics
            return None, "stable"
        # Se non ci sono campioni nel tratto post-picco, usa l'ultimo campione disponibile.
        if search_start >= len(grad_proc):
            split_idx = len(grad_proc) - 1
            split_time = float(time_search[split_idx])
            split_grad = float(grad_proc[split_idx])
        else:
            # Trova, nel tratto di ricerca, il campione con gradiente piu vicino alla threshold.
            rel_idx = int(np.argmin(np.abs(grad_proc[search_start:] - threshold)))
            # Converte l'indice relativo in indice assoluto su grad_proc.
            split_idx = search_start + rel_idx
            # Imposta split_time al tempo del campione piu vicino.
            split_time = float(time_search[split_idx])
            # Imposta split_grad al valore del gradiente nel campione di fallback.
            split_grad = float(grad_proc[split_idx])
    else:
        # Se c'e crossing reale, valuta il gradiente nel tempo di split tramite interpolazione.
        split_grad = float(np.interp(split_time, time_search, grad_proc))

    # Costruisce la maschera dei punti successivi (o uguali) allo split per analizzare il trend post-split.
    post_mask = time_crop >= split_time
    # Estrae il vettore tempo post-split.
    time_post = time_crop[post_mask]
    # Estrae la baseline post-split.
    baseline_post = baseline_crop[post_mask]

    # Se i punti post-split sono meno di 3, non e possibile stimare una regressione stabile.
    if len(time_post) < 3:
        # Se richiesto, prepara il dizionario diagnostico completo.
        if return_diagnostics:
            # Crea diagnostics con tutte le serie intermedie utili al debug.
            diagnostics = {
                "time_crop": time_diag,
                "baseline_crop": baseline_diag,
                "grad": grad,
                "grad_proc": grad_proc,
                "threshold": float(threshold),
                "peak_idx": int(peak_idx),
                "peak_grad": float(peak_grad),
                "mean_post_peak": float(mean_post_peak),
                "threshold_reference": float(threshold_reference),
                "split_time": float(split_time),
                "split_grad": split_grad,
                "crossing_found": bool(crossing_found),
            }
            # Restituisce split, trend stabile e diagnostics.
            return float(split_time), "stable", diagnostics
        # Restituisce split e trend stabile senza diagnostics.
        return float(split_time), "stable"

    # Calcola la pendenza della regressione lineare della baseline nel tratto post-split.
    slope = np.polyfit(time_post, baseline_post, 1)[0]

    # Se la pendenza e positiva oltre soglia numerica, classifica il trend come crescente.
    if slope > 1e-17:
        # Trend crescente.
        trend = "increasing"
    # Se la pendenza e negativa oltre soglia numerica, classifica il trend come decrescente.
    elif slope < -1e-17:
        # Trend decrescente.
        trend = "decreasing"
    else:
        # Altrimenti classifica il trend come stabile.
        trend = "stable"

    # Se richiesto, prepara il dizionario diagnostico finale.
    if return_diagnostics:
        # Crea diagnostics con tutte le variabili principali del calcolo.
        diagnostics = {
            "time_crop": time_diag,
            "baseline_crop": baseline_diag,
            "grad": grad,
            "grad_proc": grad_proc,
            "threshold": float(threshold),
            "peak_idx": int(peak_idx),
            "peak_grad": float(peak_grad),
            "mean_post_peak": float(mean_post_peak),
            "threshold_reference": float(threshold_reference),
            "split_time": float(split_time),
            "split_grad": split_grad,
            "crossing_found": bool(crossing_found),
        }
        # Restituisce split, trend e diagnostics completi.
        return float(split_time), trend, diagnostics

    # Restituisce split e trend quando diagnostics non sono richiesti.
    return float(split_time), trend

def plot_bootstrap_regression(ax, x, y, trend, split_time, colors=None, style=None, tick_fontsize=None, ploteq=False):
    """
    Plotta la retta di regressione e una banda di errore variabile nel tempo.
    Se il trend è 'stable', assume una retta orizzontale e stima la banda con bootstrap.
    """

    x = np.array(x)
    y = np.array(y)

    # Maschera: dati dopo lo split_time
    mask = x >= split_time
    x_post = x[mask]
    y_post = y[mask]

    if len(x_post) < 3:
        return

    # Assicura che il disegno parta esattamente da split_time anche quando
    # split_time cade tra due time step campionati.
    if x_post[0] > split_time:
        x_plot = np.concatenate(([split_time], x_post))
    else:
        x_plot = x_post.copy()

    # obtain color mapping (fallback to global style)
    if colors is None:
        try:
            colors = cfg.style.colors
        except Exception:
            colors = {}

    if style is None:
        style = cfg.style

    if tick_fontsize is None:
        tick_fontsize = getattr(style, "tick_fontsize", None)

    if trend == "stable":
        y_mean = np.nanmean(y_post)

        # Bootstrap per banda di errore
        n_boot = 1000
        boot_means = np.empty(n_boot)
        for i in range(n_boot):
            sample = np.random.choice(y_post, size=len(y_post), replace=True)
            boot_means[i] = np.nanmean(sample)

        lower_const = np.percentile(boot_means, 2.5)
        upper_const = np.percentile(boot_means, 97.5)
        lower = np.full_like(x_plot, lower_const, dtype=float)
        upper = np.full_like(x_plot, upper_const, dtype=float)

        ax.plot(x_plot, [y_mean] * len(x_plot), color=colors["regression"], lw=1.2)
        ax.fill_between(x_plot, lower, upper, color=colors["regression"], alpha=0.2)

        if tick_fontsize is not None:
            ax.tick_params(axis="both", labelsize=tick_fontsize)

        mean_width = np.mean(upper - lower)
        eqn = f"y = {y_mean:.2e} ± {mean_width:.2e}"

    else:
        coeffs = np.polyfit(x_post, y_post, 1)
        slope, intercept = coeffs
        y_fit = np.polyval(coeffs, x_plot)

        # Bootstrap per intervallo di confidenza
        n_boot = 1000
        boot_preds = np.empty((n_boot, len(x_plot)))
        for i in range(n_boot):
            idx = np.random.choice(len(x_post), len(x_post), replace=True)
            x_s, y_s = x_post[idx], y_post[idx]
            c = np.polyfit(x_s, y_s, 1)
            boot_preds[i] = np.polyval(c, x_plot)

        lower = np.percentile(boot_preds, 2.5, axis=0)
        upper = np.percentile(boot_preds, 97.5, axis=0)

        ax.plot(x_plot, y_fit, color=colors["regression"], lw=1.2)
        ax.fill_between(x_plot, lower, upper, color=colors["regression"], alpha=0.2)

        eqn = f"y = {slope:.2e}x + {intercept:.2e}"

    # Inserisce equazione in basso a destra con parametri da config
    if ploteq:
        ax.text(
            0.98, 0.02,
            eqn,
            transform=ax.transAxes,
            fontsize=style.regression_fontsize,
            fontfamily=style.regression_font_family,
            color=style.regression_color,
            ha='right',
            va='bottom',
            bbox=dict(facecolor=style.legend_facecolor, alpha=0.6, edgecolor=style.legend_edgecolor)
        )


# to check the processing procedure 
def plot_SR_processing_sequence(
    base_path,
    plast,
    visc,
    mod_code,
    vel_code,
    signal_name, # "Nemin" or "SRmean"
    output_dir="outputs/processing_steps",
    style=cfg.style,
    threshold_ratio=cfg.threshold_sr,
    pre_smoothing_sigma=cfg.smoothing_sigma,
    min_nemin_max=cfg.nmin_nemin_max,
):
    """
    
    Build the final diagnostic plot for one SR2 trace without touching the
    dataframe.

    Parameters
    ----------
    base_path : str
        Directory containing SR2 files.
    plast, visc, mod_code, vel_code : str
        Identifiers used in SR2 file names.
    signal_name : str, optional
        "Nemin" or "SRmean".
    output_dir : str, optional
        Root directory where diagnostic images are saved.
    threshold_ratio : float, optional
        Fraction of peak gradient used to define the split threshold.
    pre_smoothing_sigma : float, optional
        Sigma for Gaussian smoothing applied to the signal before split
        analysis.
    min_nemin_max : float, optional
        Execute the processing sequence only if max(Nemin) for the selected
        parameter combination is strictly greater than this value.
    Returns
    -------
    dict
        Summary with output path, split time, trend and number of saved frames.
    """
    os.makedirs(output_dir, exist_ok=True)

    colors = style.colors
    linestyles = style.linestyles
    font_family = style.font_family
    fig_facecolor = style.fig_facecolor
    axes_facecolor = style.axes_facecolor
    title_fontsize = style.title_fontsize
    label_fontsize = style.axis_label_fontsize
    tick_fontsize = style.tick_fontsize
    tick_length   = style.tick_length
    legend_fontsize = style.legend_fontsize
    secondary_lw = style.secondary_linewidth
    main_lw = style.main_linewidth
    split_lw = style.split_linewidth
    grid_lw = style.grid_linewidth
    processing_figsize = getattr(style, "figsize", (12, 8))
    processing_small_figsize = getattr(style, "small_figsize", (12, 4))

    signal_key = signal_name.strip().lower()
    if signal_key not in {"nemin", "srmean"}:
        raise ValueError("signal_name must be 'Nemin' or 'SRmean'")

    print(f"Processing sequence for {signal_name} - {mod_code} {vel_code} - plast {plast} visc {visc}")
    if mod_code == "spont":
        filename = f"SR2_{visc}_{plast}_{mod_code}.txt"
    else:
        filename = f"SR2_{visc}_{plast}_{mod_code}_{vel_code}.txt"
    path = os.path.join(base_path, filename)
    tag = f"{signal_key}_{mod_code}_{vel_code}_{visc}_{plast}"
    os.makedirs(output_dir, exist_ok=True)

    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing SR2 file: {path}")

    data = np.loadtxt(path, skiprows=1)
    if data.ndim != 2 or data.shape[1] < 4:
        raise ValueError(f"Invalid SR2 format in {path}")

    tempo = np.asarray(data[:, 0], dtype=float)
    Nemin = np.nan_to_num(np.asarray(data[:, 2], dtype=float), nan=0.0)
    SRmean = np.nan_to_num(np.asarray(data[:, 3], dtype=float), nan=0.0) *1e16

    # Run diagnostic sequence only when the selected signal reaches a meaningful peak.
    nemin_max = float(np.nanmax(Nemin)) if Nemin.size > 0 else 0.0
    srmean_max = float(np.nanmax(SRmean)) if SRmean.size > 0 else 0.0
    signal_max = nemin_max if signal_key == "nemin" else srmean_max
    if signal_max <= float(min_nemin_max):
        print(
            f"Skipping processing sequence for {filename}: "
            f"max {signal_name}={signal_max:.3g} <= {min_nemin_max} "
        )
        return {
            "file": filename,
            "signal": signal_name,
            "output_dir": output_dir,
            "split_time": None,
            "trend": "skipped",
            "n_frames": 0,
            "skipped": True,
            "skip_reason": f"max {signal_name} <= {min_nemin_max}",
            "nemin_max": nemin_max,
            "srmean_max": srmean_max,
        }

    if Nemin.size > 0:
        Nemin[0] = 0.0
    if SRmean.size > 0:
        SRmean[0] = 0.0

    raw_signal = Nemin if signal_key == "nemin" else SRmean
    if signal_key == "nemin":
        signal_label = "A"
        derivative_label = r"$dA/dt$"
    else:
        signal_label = "mSR"
        derivative_label = r"$d(mSR)/dt$"

    frame_idx = 0

    baseline_plot = gaussian_filter1d(raw_signal.copy(), sigma=float(pre_smoothing_sigma))
    if baseline_plot.size > 0:
        baseline_plot[0] = 0.0

    split_mask = (Nemin >= float(cfg.nmin_nemin_max)) if signal_key == "srmean" else None
    split_time, trend_plot, diagnostics_plot = analyze_baseline_trend(
        tempo,
        baseline_plot,
        threshold_ratio=float(threshold_ratio),
        return_diagnostics=True,
        strict_split=True,
        split_mask=split_mask,
    )
    split_time_sigma = 0.0
    


    if diagnostics_plot is not None:
        fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=processing_figsize, sharex=True, gridspec_kw={'hspace': 0.12})
        if fig_facecolor is not None:
            fig.patch.set_facecolor(fig_facecolor)
        if axes_facecolor is not None:
            ax_top.set_facecolor(axes_facecolor)
            ax_bot.set_facecolor(axes_facecolor)
        if font_family is not None:
            plt.rcParams["font.family"] = font_family

        if signal_key == "nemin":
            ax_top.text(
                0.97,
                0.95,
                "a",
                transform=ax_top.transAxes,
                ha="right",
                va="top",
                fontsize=label_fontsize,
                color=style.label_color,
                bbox=dict(facecolor=style.fig_facecolor, edgecolor=style.label_color, boxstyle="round,pad=0.25", alpha=0.9),
            )
            ax_bot.text(
                0.97,
                0.95,
                "b",
                transform=ax_bot.transAxes,
                ha="right",
                va="top",
                fontsize=label_fontsize,
                color=style.label_color,
                bbox=dict(facecolor=style.fig_facecolor, edgecolor=style.label_color, boxstyle="round,pad=0.25", alpha=0.9),
            )   

        if signal_key == "srmean":
            ax_top.text(
                0.97,
                0.95,
                "c",
                transform=ax_top.transAxes,
                ha="right",
                va="top",
                fontsize=label_fontsize,
                color=style.label_color,
                bbox=dict(facecolor=style.fig_facecolor, edgecolor=style.label_color, boxstyle="round,pad=0.25", alpha=0.9),
            )
            ax_bot.text(
                0.97,
                0.95,
                "d",
                transform=ax_bot.transAxes,
                ha="right",
                va="top",
                fontsize=label_fontsize,
                color=style.label_color,
                bbox=dict(facecolor=style.fig_facecolor, edgecolor=style.label_color, boxstyle="round,pad=0.25", alpha=0.9),
            )

        # Top panel: final step.
        ax_top.plot(tempo, raw_signal, color=colors["raw"], linewidth=secondary_lw, label="original")
        ax_top.plot(tempo, baseline_plot, color=colors["baseline"], linewidth=main_lw, label="smoothed")
        if split_time is not None:
            ax_top.axvline(split_time, color=colors["split"], linestyle=linestyles["split"], linewidth=split_lw, label=fr"$t_{{SL}}$ = {split_time:.1f} Myr")
            y_split = float(np.interp(split_time, tempo, baseline_plot))
            if np.isfinite(split_time_sigma) and split_time_sigma > 0.1:
                ax_top.errorbar(
                    split_time,
                    y_split,
                    xerr=split_time_sigma,
                    fmt="o",
                    color=colors["split"],
                    ecolor=colors["split"],
                    elinewidth=split_lw,
                    capsize=3,
                    markersize=4,
                    label=fr"$t_{{SL}}$ = {split_time:.1f} ± {split_time_sigma:.2f} Myr",
                )
            #elif np.isfinite(split_time_sigma):
            #    ax_top.plot(
            #        split_time,
            #        y_split,
            #        "o",
            #                color=colors["split"],
            #        markersize=4,
            #        label=f"split = {split_time:.1f}",
            #    )
            #else:
            #    ax_top.plot([], [], color="none", label=f"split = {split_time:.1f} Myr")
            plot_bootstrap_regression(ax_top, tempo, baseline_plot, trend_plot, split_time, style=style)
        ax_top.set_ylabel(signal_label, fontsize=label_fontsize, color=style.label_color, labelpad=style.label_pad)
        ax_top.set_xlim(float(np.nanmin(tempo)), float(np.nanmax(tempo)))
        ax_top.tick_params(axis="x", labelbottom=False, labelsize=tick_fontsize, length=tick_length, colors=style.label_color)
        ax_top.tick_params(axis="y", labelsize=tick_fontsize, length=tick_length, colors=style.label_color)
        #ax_top.grid(True, linestyle="--", linewidth=grid_lw, alpha=0.6)
        ax_top.legend(loc="best", fontsize=legend_fontsize, facecolor=style.legend_facecolor, edgecolor=style.legend_edgecolor, labelcolor=style.legend_text_color)

        # Bottom panel: derivative
        time_crop = diagnostics_plot["time_crop"]
        grad = diagnostics_plot["grad"]
        grad_proc = diagnostics_plot["grad_proc"]
        threshold = float(diagnostics_plot["threshold"])
        threshold_reference = float(diagnostics_plot.get("threshold_reference", 0.0))
        grad_label = "gradient"

        #ax_bot.plot(time_crop, grad, color=colors["gradient"], linewidth=main_lw, alpha=0.7, label="grad")
        ax_bot.plot(time_crop, grad_proc, color=colors["gradient_processed"], linewidth=main_lw, label=grad_label)

        ax_bot.axhline(threshold, color=colors["threshold"], linestyle=linestyles["threshold"], linewidth=secondary_lw,
                       label=f"threshold={threshold:.2f}")

        if split_time is not None:
            #split_grad = float(diagnostics_plot.get("split_grad", np.interp(split_time, time_crop, grad_proc)))
            ax_bot.axvline(split_time, color=colors["split"], linestyle=linestyles["split"], linewidth=split_lw, label=fr"$t_{{SL}}$ = {split_time:.1f} Myr")
            #ax_bot.plot([split_time], [threshold], marker="o", color=colors["split"], markersize=4,
            #            label="split-threshold")
            #ax_bot.plot([split_time], [split_grad], marker="x", color=colors["gradient_processed"], markersize=6,
            #            label="split on gradient")

        ax_bot.set_xlabel("Time [Myr]", fontsize=label_fontsize, color=style.label_color, labelpad=style.label_pad)
        ax_bot.set_ylabel(derivative_label, fontsize=label_fontsize, color=style.label_color, labelpad=style.label_pad)
        ax_bot.tick_params(axis="x", labelsize=tick_fontsize, length=tick_length, colors=style.label_color)
        ax_bot.tick_params(axis="y", labelsize=tick_fontsize, length=tick_length, colors=style.label_color)
        #ax_bot.grid(True, linestyle="--", linewidth=grid_lw, alpha=0.6)
        ax_bot.legend(loc="best", fontsize=legend_fontsize, facecolor=style.legend_facecolor, edgecolor=style.legend_edgecolor, labelcolor=style.legend_text_color)

        fig.tight_layout()
        save_path = os.path.join(output_dir, f"{signal_label}_{mod_code}_{vel_code}_{visc}_{plast}.png")
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        frame_idx += 1
    else:
        # Final panel with split and regression overlay.
        fig, ax = plt.subplots(figsize=processing_small_figsize)
        if fig_facecolor is not None:
            fig.patch.set_facecolor(fig_facecolor)
        if axes_facecolor is not None:
            ax.set_facecolor(axes_facecolor)
        ax.plot(tempo, raw_signal, color=colors["raw"], linewidth=main_lw, label=f"{signal_label} raw")
        ax.plot(tempo, baseline_plot, color=colors["baseline"], linewidth=main_lw, label="baseline final")
        if split_time is not None:
            ax.axvline(split_time, color=colors["split"], linestyle=linestyles["split"], linewidth=split_lw, label="_nolegend_")
            y_split = float(np.interp(split_time, tempo, baseline_plot))
            if np.isfinite(split_time_sigma) and split_time_sigma > 0:
                ax.errorbar(
                    split_time,
                    y_split,
                    xerr=split_time_sigma,
                    fmt="o",
                    color=colors["split"],
                    ecolor=colors["split"],
                    elinewidth=split_lw,
                    capsize=3,
                    markersize=4,
                    label=fr"$t_{{SL}}$ = {split_time:.2f} ± {split_time_sigma:.2f} Myr",
                )
            elif np.isfinite(split_time_sigma):
                ax.plot(
                    split_time,
                    y_split,
                    "o",
                    color=colors["split"],
                    markersize=4,
                    label=fr"$t_{{SL}}$ = {split_time:.1f} ± {split_time_sigma:.2f} Myr",
                )
            else:
                ax.plot([], [], color="none", label=fr"$t_{{SL}}$ = {split_time:.1f} Myr")
            plot_bootstrap_regression(ax, tempo, baseline_plot, trend_plot, split_time, style=style)
        ax.set_xlabel("Time [Myr]", fontsize=label_fontsize, color=style.label_color, labelpad=style.label_pad)
        ax.set_ylabel(signal_label, fontsize=label_fontsize, color=style.label_color, labelpad=style.label_pad)
        ax.tick_params(axis="x", labelsize=tick_fontsize, colors=style.label_color)
        ax.tick_params(axis="y", labelsize=tick_fontsize, colors=style.label_color)
        ax.grid(True, linestyle="--", linewidth=grid_lw, alpha=0.6)
        ax.legend(loc="best", fontsize=legend_fontsize, facecolor=style.legend_facecolor, edgecolor=style.legend_edgecolor, labelcolor=style.legend_text_color)
        fig.tight_layout()
        save_path = os.path.join(output_dir, f"{signal_label}_{mod_code}_{vel_code}_{visc}_{plast}.png")
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        frame_idx += 1

    return {
        "file": filename,
        "signal": signal_label,
        "output_dir": output_dir,
        "split_time": split_time,
        "split_time_sigma": split_time_sigma,
        "trend": trend_plot,
        "n_frames": frame_idx,
    }


# Obtain the split times and trends for all the SR2 traces, and save the results in the dataframe. 
def plot_SRabs_singleModel(df, base_path, params, maps, style, output_dir="outputs/"
):
    os.makedirs(output_dir, exist_ok=True)

    colors = style.colors
    linestyles = style.linestyles
    title_fontsize = style.title_fontsize
    axis_fontsize = style.axis_label_fontsize
    tick_fontsize = style.tick_fontsize
    main_lw = style.main_linewidth
    secondary_lw = style.secondary_linewidth
    split_lw = style.split_linewidth

    # Add columns to df
    new_cols = [
        "t_SL_A", "trend_A", "t_SL_mSR", "trend_mSR"
    ]

    for col in new_cols:
        if col not in df.columns:
            df[col] = np.nan

    for col in ["trend_A", "trend_mSR"]:
        if col not in df.columns:
            df[col] = pd.Series([None] * len(df), dtype="object")
        else:
            df[col] = df[col].astype(object)

    for mod_code, mod_title in maps.model_names.items():
        velocities = params.vel_order if mod_code != 'spont' else [None]

        for vel_code in velocities:
            # prepare the plot
            fig, axes = plt.subplots(nrows=6, ncols=3, figsize=(32, 54), sharex=True)
            vel_str = f"v = {params.vel_labels[vel_code]} cm/yr, " if vel_code is not None else ""
            fig.suptitle(f"Strain rate evolution ({vel_str}{mod_title})", fontsize=title_fontsize)

            for i, plast in enumerate(params.plasticities):
                for j, visc in enumerate(params.viscosities):
                    ax_Nemin = axes[i * 2, j]
                    ax_SRmean = axes[i * 2 + 1, j]

                    filename = f"SR2_{visc}_{plast}_{mod_code}_{vel_code}.txt" if vel_code is not None else f"SR2_{visc}_{plast}_spont.txt"
                    path = os.path.join(base_path, filename)

                    if os.path.exists(path):
                        try:
                            # Read the SR tracking file
                            data = np.loadtxt(path, skiprows=1)
                            tempo, Nelem, Nemin, SRmean = data[:, 0], data[:, 1], data[:, 2], data[:, 3]*1e16

                            Nemin[0] = 0
                            SRmean[0] = 0
                            Nemin = np.nan_to_num(Nemin)
                            SRmean = np.nan_to_num(SRmean)

                            nemin_max = float(np.nanmax(Nemin)) if Nemin.size > 0 else 0.0
                            do_split_analysis = nemin_max > float(cfg.nmin_nemin_max)

                            # Process Nemin time series (aligned with plot_SR_processing_sequence).
                            baseline_Nemin = gaussian_filter1d(Nemin.copy(), sigma=float(cfg.smoothing_sigma))
                            if baseline_Nemin.size > 0:
                                baseline_Nemin[0] = 0.0
                            split_Nemin, trend_Nemin = None, "stable"
                            if do_split_analysis:
                                split_Nemin, trend_Nemin = analyze_baseline_trend(
                                    tempo,
                                    baseline_Nemin,
                                    threshold_ratio=float(cfg.threshold_sr),
                                    strict_split=True,
                                )

                            # Process SRmean time series (aligned with plot_SR_processing_sequence).
                            baseline_SRmean = gaussian_filter1d(SRmean.copy(), sigma=float(cfg.smoothing_sigma))
                            if baseline_SRmean.size > 0:
                                baseline_SRmean[0] = 0.0
                            split_SRmean, trend_SRmean = None, "stable"
                            if do_split_analysis:
                                split_SRmean, trend_SRmean = analyze_baseline_trend(
                                    tempo,
                                    baseline_SRmean,
                                    threshold_ratio=float(cfg.threshold_sr),
                                    strict_split=True,
                                    split_mask=(Nemin >= float(cfg.nmin_nemin_max)),
                                )

                            
                            # --- Update df with results
                            mask = (
                                (df["plastic_weakening"] == plast) &
                                (df["viscous_weakening"] == visc) &
                                (df["velocity"] == (float(params.vel_labels[vel_code]) if vel_code else 0.0)) &
                                (df["type"] == mod_title)
                            )

                            # Nemin
                            if split_Nemin is not None:
                                df.loc[mask, "t_SL_A"] = split_Nemin
                                df.loc[mask, "trend_A"] = trend_Nemin

                            # SRmean
                            if split_SRmean is not None:
                                df.loc[mask, "t_SL_mSR"] = split_SRmean
                                df.loc[mask, "trend_mSR"] = trend_SRmean


                            # Plot Nemin
                            ax_Nemin.plot(tempo, Nemin, color=colors["raw"], linewidth=secondary_lw, label='original')
                            ax_Nemin.plot(tempo, baseline_Nemin, color=colors["baseline"], linewidth=secondary_lw, label='smoothed')
                            if do_split_analysis and split_Nemin is not None:
                                ax_Nemin.axvline(split_Nemin, color=colors["split"], linestyle=linestyles["split"], linewidth=split_lw, label=fr"$t_{{SL}}$ = {split_Nemin:.2f} Myr")
                                #ax_Nemin.text(split_Nemin, np.max(Nemin) * 0.9, f'Trend: {trend_Nemin}', fontsize=tick_fontsize, rotation=90)
                                ax_Nemin.text(split_Nemin, np.max(Nemin) * 0.1, fr"{split_Nemin:.2f} Myr", fontsize=tick_fontsize, rotation=90)
                                plot_bootstrap_regression(ax_Nemin, tempo, baseline_Nemin, trend_Nemin, split_Nemin)
                            
                            # Plot SRmean
                            ax_SRmean.plot(tempo, SRmean, color=colors["raw"], linewidth=secondary_lw, label='original')
                            ax_SRmean.plot(tempo, baseline_SRmean, color=colors["baseline"], linewidth=secondary_lw, label='smoothed')
                            if do_split_analysis and split_SRmean is not None:
                                ax_SRmean.axvline(split_SRmean, color=colors["split"], linestyle=linestyles["split"], linewidth=split_lw, label=fr"$t_{{SL}}$ = {split_SRmean:.2f} Myr")
                                #ax_SRmean.text(split_SRmean, np.max(SRmean) * 0.9, f'Trend: {trend_SRmean}', fontsize=tick_fontsize, rotation=90)
                                ax_SRmean.text(split_SRmean, np.max(SRmean) * 0.1, fr"{split_SRmean:.2f} Myr", fontsize=tick_fontsize, rotation=90)
                                plot_bootstrap_regression(ax_SRmean, tempo, baseline_SRmean, trend_SRmean, split_SRmean)
                            
                        except Exception as e:
                            print(f"⚠️ Error with file {filename}: {e}")

                    if i == 0:
                        ax_Nemin.set_title(params.viscosities_titles[j], fontsize=axis_fontsize, pad=cfg.style.title_pad)
                    if j == 0:
                        ax_Nemin.set_ylabel(params.plasticities_titles[i], fontsize=axis_fontsize, labelpad=cfg.style.label_pad)
                    if i == 2:
                        ax_SRmean.set_xlabel("Time [Myr]", fontsize=axis_fontsize, labelpad=cfg.style.label_pad)
                    for ax in [ax_Nemin, ax_SRmean]:
                        ax.set_xlim(0, 30)

            ax.tick_params(axis="x", labelsize=tick_fontsize)
            ax.tick_params(axis="y", labelsize=tick_fontsize)
            plt.tight_layout(rect=[0, 0, 1, 0.95])
            filename_tag = vel_code 
            save_path = os.path.join(output_dir, f"SR2_alltraces_annotated_{filename_tag}_{mod_code}.png")
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            plt.close()

    print("   Dataframe updated with the strain rate analysis results")
    return df





# Nice figures for comparing the strain rate localization time with the weakening percentages and areas at that time, for all models.



# Visualization functions for raw Strain Rate

    
