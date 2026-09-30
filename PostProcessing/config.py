from types import SimpleNamespace

import numpy as np

# start TO MODIFIY:
# === General settings ===
TRACKING_DIRECTORY = "../Tracking/"

# === Physical parameters ===
SI_DEPTH_KM = np.asarray([110, 132, 177])  # Reference depths for the subduction initiation [km]
SI_depth = (700 - SI_DEPTH_KM) * 1000  # Convert to meters


# = peak removal 
prominence = 0.05  # Prominence threshold for peak detection SR signals
distance = 3 # Minimum distance (number of data points) between peaks for SR signals
window = 1 # Window size for SR signals
maxpeaks = 10 # Maximum number of peaks to detect for SR signals

# smoothing and threshold of strain rate time series
smoothing_sigma = 2  # Standard deviation for Gaussian smoothing of strain rate time series (adim - number of samples)
threshold_sr = 0.0  # Threshold ratio for identifying significant strain rate peaks (adim)
nmin_nemin_max = 50  # Minimum number of Nemin peaks required for a valid split time estimate
combined_gravity_overlap_myr = 1.5  # Last part of the gravitational phase appended to combined traces [Myr]

# stop TO MODIFIY



# start DO NOT MODIFY:
# === Model parameters + dispay order ===

params = SimpleNamespace(
    viscosities  = ['B', 'D', 'C'],
    plasticities = ['PW2', 'PW4', 'PW3'],
    viscosities_titles  = ["B: $\\epsilon_v=5-10$", "D: $\\epsilon_v=0.5-1.5$", "C: $\\epsilon_v=0.2-0.5$"],
    plasticities_titles = ["PW2", "PW4", "PW3"],
    vel_order  = ['v1', 'v05', 'v025', 'v01', 'v005', 'v001'],
    vel_labels = {'v1'  : '1.0', 
                'v05' : '0.5', 
                'v025': '0.25', 
                'v01' : '0.1', 
                'v005': '0.05', 
                'v001': '0.01'}
)

# === Mappings ===

maps = SimpleNamespace(
    # velocity code → numerical value [cm/yr]
    vel_map = {
        'v001': 0.01,
        'v005': 0.05,
        'v01' : 0.1,
        'v025': 0.25,
        'v05' : 0.5,
        'v1'  : 1.0,
    },
    
    # model-type code → canonical name, t_gravitational [Myr], t_convergence [Myr]
    model_type_map = {
        'spont': ('Gravitational', 30,  0),
        'f'    : ('Forced',         0, 30),
        'p1'   : ('Combined_10',   10, 30),
        'p2'   : ('Combined_20',   20, 30),
        'p3'   : ('Combined_30',   30, 30),
    },
    model_names = {
        'f': "Forced",
        'p1': "Combined_10",
        'p2': "Combined_20",
        'p3': "Combined_30",
        'spont': "Gravitational"
    },
    
    # viscous-weakening code → (epsilon_0, epsilon_1)
    visc_weaken_map = {
        'B': (5,   10),
        'C': (0.2,  0.5),
        'D': (0.5,  1.5),
    },
    
    # plasticity code → general strength order
    plasticity_map = {
        'PW3': {"UCC": 2, "LCC": 2, "LOC": 2, "UOC": 2},
        'PW2': {"UCC": 1, "LCC": 3, "LOC": 3, "UOC": 3},
        'PW4': {"UCC": 1, "LCC": 2, "LOC": 2, "UOC": 1},
    },
    plasticity_default = {"UCC": None, "LCC": None, "LOC": None, "UOC": None},
    
    # canonical model name → numeric encoding - this is indicative of the t_gravitational
    model_type_encoding = {
        "Forced"      :   0,
        "Combined_10" :  10,
        "Combined_20" :  20,
        "Combined_30" :  30,
        "Gravitational": 999,
    },
    # aliases to normalise model-type strings from lowercase to canonical 
    model_type_aliases = [
        ("combined_10",  "Combined_10"),
        ("combined_20",  "Combined_20"),
        ("combined_30",  "Combined_30"),
        ("forced",       "Forced"),
        ("gravitational","Gravitational"),
    ],
    # --- viscous-weakening code → representative scalar (midpoint) ---------------
    visc_weaken_encoding = {
        'B': 7.5,
        'C': 0.35,
        'D': 1.0,
    },
    
    # --- plasticity code → numeric encoding -------------------------------------
    plasticity_encoding = {
        'PW2': 0,
        'PW3': 1,
        'PW4': 2,
    }
)

# stop DO NOT MODIFY




# === Styles and parameters for plotting ===

# ===================================================================
# === LEVEL 2: ACTIVE STYLE NAMESPACE ===
# ===================================================================
# This namespace contains only values used by the active plotting workflow.

style = SimpleNamespace(
    # Font configuration
    font_family = "Liberation Sans",
    fig_facecolor = "white",
    axes_facecolor =  "white",
    
    # Font sizes
    title_fontsize = 32, # column or row title, as PW3 or B:ev=5-10
    axis_label_fontsize = 28, # axis label, as Time[Myr]
    tick_fontsize = 24, # label of the ticks of the axes, as 0,10,20 ecc for Time
    legend_fontsize = 24,
    regression_fontsize = 24,
    
    # Line widths and alpha
    main_linewidth = 3,
    secondary_linewidth = 1,
    split_linewidth = 3.5,
    grid_linewidth = 0.5,

    # Color and linestyle dictionaries
    colors = {
        "raw": "#718096",
        "baseline": "#0d2b52",
        "gradient": "#0d2b52",
        "gradient_processed": "#0d2b52",
        "split": "#D4622A",
        "threshold": "#2A6FA8",
        "regression": "#D4622A",
    },
    linestyles = {
        "raw": "-",
        "baseline": "-",
        "gradient": "-",
        "gradient_processed": "-",
        "split": "--",
        "threshold": ":",
    },

    # Figure sizes used by processing/diagnostic plots
    small_figsize = (12, 4),
    figsize = (18, 12),
    
    # Tick and label styling
    regression_font_family = "DejaVu Sans Mono",
    regression_color = "black",
    label_color = "black", #font color
    title_pad = 20,
    label_pad = 20,
    legend_facecolor = "white",
    legend_edgecolor = "white",
    legend_text_color ="black",
    tick_length = 10,
)


# ==============================================================
# === Workflow ===
CREATE_DATABASE = True # if true, tSL and tSZI are computed from the Tracking files and tabella_subduzione.csv is created
LOAD_DATABASE = True # if true , tabella_subduzione.csv is loaded into a dataframe and the plot functions can read it
PLOT_SR_DIAGNOSTIC = True # if true, produces a plot for A and a plot for mSR for each model, showing the tSL and tSZI assessment process. Useful for understanding the assessment process and analysing the single models. 
