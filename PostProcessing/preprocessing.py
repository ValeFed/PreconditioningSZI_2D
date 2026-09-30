import os

import pandas as pd

from StrainRate import plot_SRabs_singleModel

CSV_COLUMNS = [
    "viscous_weakening",
    "plastic_weakening",
    "type",
    "velocity",
    "gravitational_phase_Myr",
    "convergent_phase_Myr",
    "t_SL_A",
    "trend_A",
    "t_SL_mSR",
    "trend_mSR",
    "tSL",
    "dt_SL",
    "t_SZI",
    "dtSZI",
]

COLUMN_RENAME = {
    "weakening_viscoso": "viscous_weakening",
    "weakening_plastico": "plastic_weakening",
    "tipo_modello": "type",
    "velocità_imposta": "velocity",
    "tempo_quiescenza_Ma": "gravitational_phase_Myr",
    "tempo_convergenza_Ma": "convergent_phase_Myr",
    "split_time_Nemin": "t_SL_A",
    "trend_Nemin": "trend_A",
    "split_time_SRmean": "t_SL_mSR",
    "trend_SRmean": "trend_mSR",
    "split_time_min_nemin_srmean": "tSL",
    "delta_split_min_combined_minus_forced": "dt_SL",
    "SZI": "t_SZI",
    "delta_SZI_combined_minus_forced": "dtSZI",
}


def create_df(directory: str, SI_depth, params, maps, style):
    """Create the database with only the fields used by the reduced workflow."""
    records = []

    for filename in os.listdir(directory):
        if not (filename.startswith("OLMnew_") and filename.endswith(".txt")):
            continue

        parts = filename.replace(".txt", "").split("_")
        visc, plast, mod_code = parts[1], parts[2], parts[3]
        model_type, t_quiescent, t_convergent = maps.model_type_map.get(mod_code, ("Unknown", None, None))
        velocity = maps.vel_map.get(parts[4] if len(parts) > 4 else None, 0.0 if model_type == "Gravitational" else None)

        records.append({
            "viscous_weakening": visc,
            "plastic_weakening": plast,
            "type": model_type,
            "velocity": velocity,
            "gravitational_phase_Myr": t_quiescent,
            "convergent_phase_Myr": t_convergent,
        })

    df_result = pd.DataFrame(records)
    if df_result.empty:
        raise RuntimeError(f"No OLMnew tracking files found in {directory}")

    # This original function calculates and plots split times for all SR2 traces.
    plot_SRabs_singleModel(df_result, directory, params, maps, style)

    df_result["tSL"] = df_result[["t_SL_A", "t_SL_mSR"]].min(axis=1, skipna=True)
    df_result["dt_SL"] = float("nan")
    df_result["dtSZI"] = float("nan")
    df_result["t_SZI"] = pd.NA

    combined_models = {"Combined_10", "Combined_20", "Combined_30"}
    combined_mask = df_result["type"].isin(combined_models)

    for idx, row in df_result.iterrows():
        if df_result.at[idx, "trend_A"] == "increasing":
            df_result.at[idx, "t_SZI"] = df_result.at[idx, "tSL"]

    key_cols = ["plastic_weakening", "viscous_weakening", "velocity"]
    forced_ref = (
        df_result.loc[df_result["type"] == "Forced", key_cols + ["t_SL_A", "t_SL_mSR", "tSL", "t_SZI"]]
        .drop_duplicates(subset=key_cols, keep="first")
        .rename(columns={
            "t_SL_A": "t_SL_A_forced",
            "t_SL_mSR": "t_SL_mSR_forced",
            "tSL": "tSL_forced",
            "t_SZI": "t_SZI_forced",
        })
    )
    forced_ref["t_SZI_forced"] = forced_ref["t_SZI_forced"].fillna(0)
    combined_rows = df_result.loc[combined_mask, key_cols + ["t_SL_A", "t_SL_mSR", "tSL", "t_SZI"]]
    combined_with_forced = combined_rows.merge(forced_ref, on=key_cols, how="left")
    df_result.loc[combined_mask, "dt_SL"] = (combined_with_forced["tSL"] - combined_with_forced["tSL_forced"]).to_numpy()
    df_result.loc[combined_mask, "dtSZI"] = (combined_with_forced["t_SZI"] - combined_with_forced["t_SZI_forced"]).to_numpy()

    df_result = df_result[CSV_COLUMNS]
    df_result.to_csv("tabella_subduzione.csv", index=False)
    print("   Database generated: tabella_subduzione.csv")
    return df_result


def load_df(maps, path="tabella_subduzione.csv"):
    df = pd.read_csv(path).rename(columns=COLUMN_RENAME).loc[:, CSV_COLUMNS]
    numeric_cols = [
        "velocity", "gravitational_phase_Myr", "convergent_phase_Myr",
        "t_SL_A", "t_SL_mSR", "tSL", "dt_SL", "t_SZI", "dtSZI",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for raw, canonical in maps.model_type_aliases:
        df["type"] = df["type"].str.strip().str.replace(raw, canonical, case=False, regex=False)
    return df