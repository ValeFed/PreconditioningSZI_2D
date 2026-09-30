import sys

import config as cfg
from preprocessing import create_df, load_df
from StrainRate import plot_SR_processing_sequence


def main():
    if cfg.CREATE_DATABASE:
        print(" > Generating strain-rate database from SR2 tracking files...")
        create_df(cfg.TRACKING_DIRECTORY, cfg.SI_depth, cfg.params, cfg.maps, cfg.style)

    if not cfg.LOAD_DATABASE:
        print("No database loaded. Exiting.")
        sys.exit(0)

    print(" > Loading database...")
    df = load_df(cfg.maps)

    if cfg.PLOT_SR_DIAGNOSTIC:
        print(" > Plotting strain-rate diagnostic sequences...")
        for p in ["PW2", "PW3", "PW4"]:
            for v in ["B", "C", "D"]:
                for vel in ["v001", "v005", "v01", "v025", "v05", "v1"]:
                    for s in ["Nemin", "SRmean"]:
                        for m in ["f", "p1", "p2", "p3", "spont"]:
                            plot_SR_processing_sequence(
                                cfg.TRACKING_DIRECTORY,
                                p,
                                v,
                                m,
                                vel,
                                signal_name=s,
                                pre_smoothing_sigma=cfg.smoothing_sigma,
                            )

    print("Workflow completed: mSR/A processing figures and tSL/tSZI table.")


if __name__ == "__main__":
    main()