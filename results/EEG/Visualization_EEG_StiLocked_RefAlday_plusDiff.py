#! /usr/bin/env python3
# -*- coding: UTF-8 -*-

"""
Master Analysis & Visualization Script: Stimulus-Locked Components
=============================================================================
v2.1.1 — Bugfixes for cross-experiment robustness.
v2.1   — Aligned with upstream LMM analyses (Sta_EEG_OfferPhase_RefAlday.Rmd v2.1
         and Sta_EEG_OfferPhase_RefTraditional.Rmd v2.1).

[v2.1.1 BUGFIX SUMMARY]:
    1. load_data_csv now standardizes participant_id BEFORE the valid_set
       membership test. Previously this comparison used the raw filename
       prefix (e.g. 'Vp0000'), which mismatched valid_keys entries
       (canonicalized to 'VP0000' by load_trials_csv_filtered). On
       experiments where epoch filenames used non-canonical casing, all
       trials silently failed the membership test and load_data_csv
       returned (None, None, None), causing RegressionAnalysis to halt
       with a generic "No epoch data" message — no figures emitted.
    2. valid_set membership test rewritten as vectorized .isin() (was
       df.apply(lambda r: ..., axis=1), which scaled with rows × subjects
       and could take 10+ minutes on the full dataset).
    3. Bare except block in load_data_csv replaced with explicit error
       reporting (type + filename + message); previously masked legitimate
       file-level errors as silent skips.
    4. [HALT] message in RegressionAnalysis.run now prints diagnostic
       context (valid_keys size, expected ROI channels, common causes)
       to make root-cause analysis tractable.

Description:
    Unified visualization pipeline for both Experiment 1 (E1) and Experiment 2 (E2).
    Generates publication-ready regression-corrected ERP waveforms and QC plots.
    Component definitions are read from config/erp_components.csv (single source of truth).
    
    STATISTICAL MODEL (REGRESSION METHOD) DYNAMICS:
    - Base model: Signal ~ Intercept + Emotion + Offer + Emo*Offer + Baseline
    - Per-component baseline interactions are auto-detected from the upstream
      LMM REPRODUCIBILITY.txt files (see RegressionAnalysis.run for the
      auto-parsing logic). When the upstream LMM retained `emotion:Baseline_c`
      and/or `offer_type:Baseline_c`, the OLS design matrix here is extended
      accordingly so that the reconstructed marginal-mean waveforms match the
      LMM's fixed-effect structure.
    
    [v2.1 BASELINE PARAMETERIZATION]:
    Baseline ('BL') is now mean-centered (NOT Z-scored) before OLS fitting,
    matching the upstream LMM convention (Alday, 2019). Mean-centering keeps
    baseline in microvolt units so that beta_baseline is interpretable as
    'fraction of pre-stimulus drift propagated forward'. When reconstructing
    the marginal means (setting BL=0 in the equation), the waveforms reflect
    the amplitude at the *average* baseline level.
    
    [v2.1 TRIAL INCLUSION]:
    Trial-set selection now flows from trials.csv (the canonical post-
    preprocessing single-trial data frame) rather than independent peak-to-
    peak rejection. Inclusion criteria mirror Sta_Behaviour_RefAlday.Rmd
    and Sta_EEG_OfferPhase_RefAlday.Rmd v2.1: Offers_Other in {5,6,8,9},
    reaction != 0, RT in [300, 3000] ms, non-NaN component scalar, non-NaN
    Baseline_*. This ensures the figure N matches the LMM N exactly.
    
    [v2.1 SANITY CHECKS]:
    Four-layer alignment self-checks run automatically before figure
    generation:
      Sanity 1: participant_id consistency (trials.csv vs epoch CSVs)
      Sanity 2: per-subject trial coverage (LMM-valid keys present in epoch)
      Sanity 3: post-join trial count equals the LMM N
      Sanity 4: component-window mean from epoch matches trials.csv scalar
    A failure on any layer halts the run with a diagnostic message.

    MODULES:
    1. RegressionAnalysis: OLS-corrected waveforms mapping LMM structure.
    2. StandardAnalysis: Traditional subtraction-based baseline (for reference).
    3. QC_Pipeline_Uncorrected: MNE-generated raw averages (no baseline).
=============================================================================
"""

import matplotlib
# Force non-interactive backend for server/headless execution safety
matplotlib.use('Agg') 

import os
import re
import sys
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path
from tqdm import tqdm
import warnings
from datetime import datetime
from scipy.ndimage import gaussian_filter1d

# Suppress runtime warnings for cleaner standard output logs
warnings.filterwarnings('ignore')

# ==============================================================================
# 1. GLOBAL EXPERIMENT TOGGLE & DYNAMIC CONFIGURATION
# ==============================================================================
# Switch between 'E1' and 'E2' to dynamically adjust component lists and windows
EXPERIMENT_VERSION = 'E1'

# ------------------------------------------------------------------------------
# Component definitions are read from the single source of truth shared with the
# pipeline, stats, and topography scripts: config/erp_components.csv
# (t_min/t_max in SECONDS; roi pipe-delimited). Edit ONLY that file to change a
# window/ROI or add/remove a component.
# ------------------------------------------------------------------------------
def _find_component_config():
    p = Path(__file__).resolve()
    for parent in [p.parent] + list(p.parents):
        candidate = parent / 'config' / 'erp_components.csv'
        if candidate.exists():
            return candidate
    raise FileNotFoundError('config/erp_components.csv not found above ' + str(p))

_COMP_CSV = _find_component_config()
_COMP_DF = pd.read_csv(_COMP_CSV)
_COMP_DF = _COMP_DF[_COMP_DF['experiment'] == EXPERIMENT_VERSION].reset_index(drop=True)
BATCH_COMPONENTS = _COMP_DF['comp_name'].tolist()

# ==============================================================================
# 2. VISUAL CONFIGURATION
# ==============================================================================

# [VISUAL STYLE]
LINE_WIDTH = 2.0          
SMOOTHING_SIGMA = 1.2     

# [MANUAL SCALES] Y-Axis Settings (Microvolts): (Y_MIN, Y_MAX, TICK_STEP)
SCALE_CONFIG = {
    row['comp_name']: (float(row['scale_min']), float(row['scale_max']), float(row['scale_step']))
    for _, row in _COMP_DF.iterrows()
}

# ==============================================================================
# 3. SYSTEM CONFIGURATION
# ==============================================================================

# ANALYSIS_METHOD Options: 'All', 'Regression', 'Standard', 'QC'
ANALYSIS_METHOD = 'All' 

REJECT_PTP_THRESHOLD_UV = 200.0

# ------------------------------------------------------------------------------
# Trial-level inclusion criteria — must match upstream LMM scripts.
# ------------------------------------------------------------------------------
# These thresholds define which trials enter visualization. They are applied
# at the analysis layer (in load_trials_csv_filtered below), NOT in the upstream
# HU pipeline preprocessing. Rationale: the canonical trials.csv reflects only
# EEG-artifact rejection (peak-to-peak > 200 uV marked as NaN); response-related
# and RT-range exclusions are analysis-specific decisions.
#
# IMPORTANT: keep these constants synchronized with:
#   - Sta_EEG_OfferPhase_RefAlday.Rmd v2.1 (rt_lower_ms / rt_upper_ms)
#   - Sta_EEG_OfferPhase_RefTraditional.Rmd v2.1 (same)
#   - Sta_Behaviour_RefAlday.Rmd (same)
# Cross-modality and cross-method sample equivalence requires identical logic.
RT_LOWER_MS = 300
RT_UPPER_MS = 3000
INCLUSION_OFFERS = [5, 6, 8, 9]

# Floating-point tolerance for Sanity 4 (component-window scalar match).
# Recomputing the mean over ~25-150 timepoints differs from the precomputed
# scalar by floating-point rounding only (~1e-7 uV); 0.01 uV is conservative.
SCALAR_TOLERANCE_UV = 0.01

COLORS = {
    'neu': "#8491B4", 'aff': "#91D1C2", 'dis': "#E64B35", 
    'dom': "#F39B7F", 'enj': "#3C5488"
}

LABELS_MAP = {
    'neu': 'Neutral', 'aff': 'Affiliative', 'dis': 'Disgust', 
    'dom': 'Dominance', 'enj': 'Reward'
}

OFFER_TYPES = ['fair', 'unfair']

# Component Definitions synchronized with upstream preprocessing architecture
COMPONENT_SPECS = {
    row['comp_name']: {
        'roi': str(row['roi']).split('|'),
        'window': (float(row['t_min']), float(row['t_max'])),
        'plot_xlim': (-200, 1000),
        'folder_name': 'Stimulus_Locked'
    }
    for _, row in _COMP_DF.iterrows()
}

# Full channel configuration for topographic mapping reconstruction
ALL_EEG_CHANNELS = [
    "Fp1", "Fpz", "Fp2", "AF7", "AF3", "AFz", "AF4", "AF8", "F9", "F7", "F5", "F3", "Fz", "F4", "F6", "F8", "F10",
    "FT7", "FC5", "FC3", "FC1", "FC2", "FC4", "FC6", "FT8", "T7", "C5", "C3", "C1", "Cz", "C2", "C4", "C6", "T8",
    "TP9", "TP7", "CP5", "CP3", "CP1", "CPz", "CP2", "CP4", "CP6", "TP8", "TP10", "P7", "P5", "P3", "Pz", "P4",
    "P6", "P8", "PO9", "PO7", "PO3", "POz", "PO4", "PO8", "PO10", "O1", "Oz", "O2"
]

TARGET_COMPONENT = None
CURR_CFG = None
PLOT_WINDOW = None
PATHS = None

# ==============================================================================
# 4. PATH MANAGEMENT
# ==============================================================================

def get_paths():
    """
    Dynamically resolves input/output paths based on the script's hierarchical location.
    Ensures execution portability across local and remote server environments.
    """
    script_path = Path(__file__).resolve()
    root = None
    
    # Locate project root by searching for pipeline data structure fingerprint
    for parent in [script_path.parent] + list(script_path.parents)[:5]:
        if (parent / "data" / "02_Pipeline_Output").exists():
            root = parent
            break
    if root is None:
        # Fallback evaluation logic
        if (Path.cwd() / "data").exists(): root = Path.cwd()
        else: root = script_path.parents[2] 

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    pipe_root = root / "data" / f"02_Pipeline_Output_{EXPERIMENT_VERSION}"
    folder = CURR_CFG['folder_name']
    
    out_dir_base = script_path.parent / f"Results_{TARGET_COMPONENT}_{EXPERIMENT_VERSION}" / f"Session_{ts}"
    
    paths = {
        "results": out_dir_base,
        "project_root": root,
        
        # Standard Subtraction Method Targets
        "std_ave": pipe_root / "Method_Standard" / folder / "ave.csv",
        "std_locs": pipe_root / "Method_Standard" / folder / "channel_locations.csv",
        "std_trials": pipe_root / "Method_Standard" / folder / "trials.csv",
        
        # Regression / Uncorrected Method Targets
        "reg_ave": pipe_root / "Method_Regression" / folder / "ave.csv",
        "reg_epochs": pipe_root / "Method_Regression" / folder / "epochs",
        "reg_trials": pipe_root / "Method_Regression" / folder / "trials.csv",
        
        # Baseline constraint: Values extracted exclusively from Stimulus-Locked phase
        "reg_base": pipe_root / "Baseline_Raw" / "Stimulus_Locked_Values" / "baseline_values.csv",
        "reg_locs": pipe_root / "Method_Regression" / folder / "channel_locations.csv",
        
        # Upstream LMM output: REPRODUCIBILITY.txt files for auto-parsing
        # the per-component fixed-effect specification. Dual paths support
        # both Alday (regression) and Traditional analysis lineages.
        "lmm_alday_stage1": (root / "results" / "EEG"
                             / f"{EXPERIMENT_VERSION}_TwoStage_Bates_Alday"
                             / "Stage1_TrialLevel"),
        "lmm_traditional_stage1": (root / "results" / "EEG"
                                   / f"{EXPERIMENT_VERSION}_TwoStage_Bates_Traditional"
                                   / "Stage1_TrialLevel"),
    }
    return paths

# ==============================================================================
# 5. VISUALIZATION ENGINE
# ==============================================================================

class VisualizationEngine:
    @staticmethod
    def get_scale_settings():
        if TARGET_COMPONENT in SCALE_CONFIG:
            return SCALE_CONFIG[TARGET_COMPONENT]
        return (-5.0, 5.0, 1.0)

    @staticmethod
    def _calculate_dynamic_limits(data_dict, y_min, y_max):
        """
        Internal utility to calculate dynamic Y-axis bounds for QC plots.
        Prevents waveforms with extreme artifact variance from being clipped.
        """
        all_vals = np.concatenate([v for v in data_dict.values() if not np.isnan(v).all()])
        if len(all_vals) > 0:
            data_range = np.percentile(all_vals, 99) - np.percentile(all_vals, 1)
            mid_point = np.median(all_vals)
            if data_range > (y_max - y_min):
                return mid_point - (data_range/2) - 2, mid_point + (data_range/2) + 2
        return y_min, y_max

    @staticmethod
    def _apply_standard_styling(ax, plot_xlim, win_start_ms, win_end_ms, y_min, y_max, tick_step, title_suffix):
        """
        Centralized helper for consistent plot geometry, shading, and axis formatting.
        """
        ax.axvspan(win_start_ms, win_end_ms, color='#f0f0f0', alpha=1.0, lw=0, zorder=0)
        ax.axvline(0, ls='-', color='black', lw=1.0, zorder=1)
        ax.axhline(0, ls=':', color='black', lw=1.0, zorder=1)
        
        for spine in ax.spines.values():
            spine.set_visible(True); spine.set_linewidth(1.2); spine.set_color('black')

        ax.set_xlim(plot_xlim)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
        ticks = ax.get_xticks()
        ax.set_xticklabels([str(int(x)) if x % 200 == 0 else "" for x in ticks])
        
        ax.set_ylim(y_min, y_max)
        if "Uncorrected" not in title_suffix:
            ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_step))
        
        ax.tick_params(axis='both', which='major', labelsize=14, length=6, width=1.2, direction='out')

    @staticmethod
    def plot_unified_panels(data_dict, time_ms, title_suffix, save_path):
        """
        Generates publication-quality 2-panel marginal means plots (Fair vs Unfair).
        """
        fig, axes = plt.subplots(1, 2, figsize=(18, 6), facecolor='white', dpi=300)
        
        win_start_ms = CURR_CFG['window'][0] * 1000
        win_end_ms = CURR_CFG['window'][1] * 1000
        plot_xlim = CURR_CFG['plot_xlim']
        y_min, y_max, tick_step = VisualizationEngine.get_scale_settings()
        
        if "Uncorrected" in title_suffix:
            y_min, y_max = VisualizationEngine._calculate_dynamic_limits(data_dict, y_min, y_max)

        plt.rcParams.update({'font.size': 14, 'font.family': 'sans-serif', 'font.sans-serif': ['Arial']})

        offer_style = {
            'fair': {'linestyle': '--'},
            'unfair': {'linestyle': '-'}
        }

        for idx, offer in enumerate(OFFER_TYPES):
            ax = axes[idx]
            draw_order = ['neu', 'aff', 'dis', 'dom', 'enj']
            for emo in draw_order:
                key = f'{offer}_{emo}'
                if key in data_dict:
                    raw = data_dict[key]
                    if len(raw) == len(time_ms) and not np.all(np.isnan(raw)):
                        smooth = gaussian_filter1d(raw, sigma=SMOOTHING_SIGMA)
                        ax.plot(time_ms, smooth, 
                                color=COLORS[emo], 
                                linestyle=offer_style[offer]['linestyle'], 
                                lw=LINE_WIDTH, 
                                alpha=0.9,
                                label=LABELS_MAP[emo])

            leg = ax.legend(title='Emotion', loc='upper left', bbox_to_anchor=(1.02, 1.0), fontsize=11, frameon=True)
            leg.get_frame().set_edgecolor('#AAAAAA')
            leg.get_frame().set_linewidth(1.0)

            VisualizationEngine._apply_standard_styling(ax, plot_xlim, win_start_ms, win_end_ms, y_min, y_max, tick_step, title_suffix)
            
            ax.set_xlabel('Time (ms)', fontsize=16, weight='bold')
            if idx == 0:
                ax.set_ylabel(f"{TARGET_COMPONENT} (µV)", fontsize=16, weight='bold')
            
            clean_suffix = title_suffix.replace(" (Exp 1 Model)", "").replace(" (Exp 2 Model)", "")
            ax.set_title(f"{TARGET_COMPONENT} ROI: {offer.capitalize()} Offers", fontsize=16, weight='bold', pad=12)
        
        plt.suptitle(f"{TARGET_COMPONENT} Waveforms - {clean_suffix}", fontsize=20, weight='bold', y=0.98)
        plt.tight_layout(rect=[0, 0.03, 0.92, 0.95], w_pad=8.0)
        
        plt.savefig(save_path, dpi=300, format='tiff')
        plt.savefig(save_path.with_suffix('.png'), dpi=300, format='png')
        plt.close()

    @staticmethod
    def plot_emotion_contrasts(data_dict, time_ms, title_suffix, save_path):
        """
        Generates publication-quality N-panel plots (1 column) mapping Fair vs Unfair.
        """
        draw_order = ['neu', 'aff', 'dis', 'dom', 'enj']
        num_panels = len(draw_order)
        
        fig, axes = plt.subplots(num_panels, 1, figsize=(8, 4 * num_panels), facecolor='white', dpi=300)
        
        win_start_ms = CURR_CFG['window'][0] * 1000
        win_end_ms = CURR_CFG['window'][1] * 1000
        plot_xlim = CURR_CFG['plot_xlim']
        y_min, y_max, tick_step = VisualizationEngine.get_scale_settings()
        
        if "Uncorrected" in title_suffix:
            y_min, y_max = VisualizationEngine._calculate_dynamic_limits(data_dict, y_min, y_max)

        plt.rcParams.update({'font.size': 14, 'font.family': 'sans-serif', 'font.sans-serif': ['Arial']})

        offer_style = {
            'fair': {'linestyle': '--', 'label': 'Fair Offer'},
            'unfair': {'linestyle': '-', 'label': 'Unfair Offer'}
        }

        for idx, emo in enumerate(draw_order):
            ax = axes[idx]
            current_color = COLORS[emo]
            
            for offer in OFFER_TYPES:
                key = f'{offer}_{emo}'
                if key in data_dict:
                    raw = data_dict[key]
                    if len(raw) == len(time_ms) and not np.all(np.isnan(raw)):
                        smooth = gaussian_filter1d(raw, sigma=SMOOTHING_SIGMA)
                        ax.plot(time_ms, smooth, 
                                color=current_color, 
                                linestyle=offer_style[offer]['linestyle'],
                                lw=LINE_WIDTH, 
                                alpha=0.9, 
                                label=offer_style[offer]['label'])

            VisualizationEngine._apply_standard_styling(ax, plot_xlim, win_start_ms, win_end_ms, y_min, y_max, tick_step, title_suffix)
            
            if idx == num_panels - 1:
                ax.set_xlabel('Time (ms)', fontsize=16, weight='bold')
            
            ax.set_ylabel(f"{TARGET_COMPONENT} (µV)", fontsize=16, weight='bold')
            emo_title = LABELS_MAP.get(emo, emo).upper()
            ax.set_title(f"Emotion: {emo_title}", fontsize=14, weight='bold', color=current_color, pad=8)
            ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1.0), fontsize=12, frameon=True, edgecolor='#AAAAAA')

        clean_suffix = title_suffix.replace(" (Exp 1 Model)", "").replace(" (Exp 2 Model)", "")
        plt.suptitle(f"{TARGET_COMPONENT} - Fair vs Unfair Contrasts\n{clean_suffix}", fontsize=18, weight='bold', y=0.99)
        plt.tight_layout(rect=[0, 0.0, 0.85, 0.96], h_pad=2.0)
        
        plt.savefig(save_path, dpi=300, format='tiff')
        plt.savefig(save_path.with_suffix('.png'), dpi=300, format='png')
        plt.close()

    @staticmethod
    def plot_fairness_main_effect(data_dict, time_ms, title_suffix, save_path):
        """
        Plots the main effect of fairness (collapsed across emotions) into a single panel.
        Displays fair and unfair waveforms in black with differentiated linestyles.
        """
        fig, ax = plt.subplots(figsize=(10, 6), facecolor='white', dpi=300)
        
        win_start_ms = CURR_CFG['window'][0] * 1000
        win_end_ms = CURR_CFG['window'][1] * 1000
        plot_xlim = CURR_CFG['plot_xlim']
        y_min, y_max, tick_step = VisualizationEngine.get_scale_settings()

        plt.rcParams.update({'font.size': 14, 'font.family': 'sans-serif', 'font.sans-serif': ['Arial']})

        offer_style = {
            'fair': {'linestyle': '--', 'label': 'Fair Offer (Collapsed)'},
            'unfair': {'linestyle': '-', 'label': 'Unfair Offer (Collapsed)'}
        }

        for offer in OFFER_TYPES:
            if offer in data_dict:
                raw = data_dict[offer]
                if len(raw) == len(time_ms) and not np.all(np.isnan(raw)):
                    smooth = gaussian_filter1d(raw, sigma=SMOOTHING_SIGMA)
                    ax.plot(time_ms, smooth, 
                            color='#000000', 
                            linestyle=offer_style[offer]['linestyle'], 
                            lw=LINE_WIDTH, 
                            alpha=0.9,
                            label=offer_style[offer]['label'])

        leg = ax.legend(title='Fairness Main Effect', loc='upper right', fontsize=12, frameon=True)
        leg.get_frame().set_edgecolor('#AAAAAA')

        VisualizationEngine._apply_standard_styling(ax, plot_xlim, win_start_ms, win_end_ms, y_min, y_max, tick_step, title_suffix)
        ax.set_xlabel('Time (ms)', fontsize=16, weight='bold')
        ax.set_ylabel(f"{TARGET_COMPONENT} (µV)", fontsize=16, weight='bold')
        ax.set_title(f"{TARGET_COMPONENT} ROI: Main Effect of Fairness", fontsize=16, weight='bold', pad=12)

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, format='tiff')
        plt.savefig(save_path.with_suffix('.png'), dpi=300, format='png')
        plt.close()

    @staticmethod
    def plot_emotion_main_effect(data_dict, time_ms, title_suffix, save_path):
        """
        Plots the main effect of emotion (collapsed across fairness levels) into a single panel.
        Only executed when contextually relevant (e.g., N400 component).
        """
        fig, ax = plt.subplots(figsize=(10, 6), facecolor='white', dpi=300)
        
        win_start_ms = CURR_CFG['window'][0] * 1000
        win_end_ms = CURR_CFG['window'][1] * 1000
        plot_xlim = CURR_CFG['plot_xlim']
        y_min, y_max, tick_step = VisualizationEngine.get_scale_settings()

        plt.rcParams.update({'font.size': 14, 'font.family': 'sans-serif', 'font.sans-serif': ['Arial']})

        draw_order = ['neu', 'aff', 'dis', 'dom', 'enj']
        
        for emo in draw_order:
            if emo in data_dict:
                raw = data_dict[emo]
                if len(raw) == len(time_ms) and not np.all(np.isnan(raw)):
                    smooth = gaussian_filter1d(raw, sigma=SMOOTHING_SIGMA)
                    ax.plot(time_ms, smooth, 
                            color=COLORS[emo], 
                            linestyle='-', 
                            lw=LINE_WIDTH, 
                            alpha=0.9,
                            label=LABELS_MAP[emo])

        leg = ax.legend(title='Emotion Main Effect', loc='upper left', bbox_to_anchor=(1.02, 1.0), fontsize=11, frameon=True)
        leg.get_frame().set_edgecolor('#AAAAAA')

        VisualizationEngine._apply_standard_styling(ax, plot_xlim, win_start_ms, win_end_ms, y_min, y_max, tick_step, title_suffix)
        ax.set_xlabel('Time (ms)', fontsize=16, weight='bold')
        ax.set_ylabel(f"{TARGET_COMPONENT} (µV)", fontsize=16, weight='bold')
        ax.set_title(f"{TARGET_COMPONENT} ROI: Main Effect of Emotion", fontsize=16, weight='bold', pad=12)

        plt.tight_layout(rect=[0, 0, 0.85, 1.0])
        plt.savefig(save_path, dpi=300, format='tiff')
        plt.savefig(save_path.with_suffix('.png'), dpi=300, format='png')
        plt.close()

    @staticmethod
    def _symmetric_difference_scale(diff_values_list, pad=1.15):
        """
        Symmetric-about-zero y-axis range and a 'nice' tick step for difference waves.

        Difference waves (condition A minus B) center on zero by construction, so a
        symmetric range keeps the zero reference line central and makes polarity
        directly readable. The fixed SCALE_CONFIG ranges are tuned for raw condition
        means and are inappropriate here. The tick step is the smallest candidate
        yielding <= 8 intervals across the full range, so the axis stays legible
        regardless of effect magnitude. Returns (y_min, y_max, tick_step).
        """
        all_vals = np.concatenate([
            np.asarray(v, dtype=float) for v in diff_values_list
            if v is not None and not np.all(np.isnan(v))
        ]) if diff_values_list else np.array([])
        if all_vals.size == 0 or np.all(np.isnan(all_vals)):
            # No plottable difference data: fall back to a neutral +/-1 uV frame.
            return -1.0, 1.0, 0.5
        limit = np.nanmax(np.abs(all_vals)) * pad
        tick_step = 10.0
        for step in (0.1, 0.25, 0.5, 1.0, 2.0, 2.5, 5.0, 10.0):
            if (2 * limit) / step <= 8:
                tick_step = step
                break
        # Round the limit up to a whole number of steps for clean axis ends.
        limit = float(np.ceil(limit / tick_step) * tick_step)
        return -limit, limit, tick_step

    @staticmethod
    def plot_difference_wave_collapsed(diff_wave, time_ms, title_suffix, save_path):
        """
        Single-panel difference wave: Unfair minus Fair, collapsed across emotion.

        Mirrors the upstream functional localizer Grand_Unfair - Grand_Fair
        (Pipline_UG_OfferPhase.Rmd perm_contrasts_list), isolating the economic
        fairness effect. A negative-going deflection within the shaded FRN window
        is the expected fairness-FRN signature for unfair offers. Descriptive only:
        no confidence band, so statistical inference remains with the upstream LMM.
        """
        fig, ax = plt.subplots(figsize=(10, 6), facecolor='white', dpi=300)

        win_start_ms = CURR_CFG['window'][0] * 1000
        win_end_ms = CURR_CFG['window'][1] * 1000
        plot_xlim = CURR_CFG['plot_xlim']
        y_min, y_max, tick_step = VisualizationEngine._symmetric_difference_scale([diff_wave])

        plt.rcParams.update({'font.size': 14, 'font.family': 'sans-serif', 'font.sans-serif': ['Arial']})

        if len(diff_wave) == len(time_ms) and not np.all(np.isnan(diff_wave)):
            smooth = gaussian_filter1d(diff_wave, sigma=SMOOTHING_SIGMA)
            ax.plot(time_ms, smooth, color='#000000', linestyle='-',
                    lw=LINE_WIDTH, alpha=0.9, label='Unfair - Fair (collapsed)')

        leg = ax.legend(title='Difference Wave', loc='upper right', fontsize=12, frameon=True)
        leg.get_frame().set_edgecolor('#AAAAAA')

        VisualizationEngine._apply_standard_styling(
            ax, plot_xlim, win_start_ms, win_end_ms, y_min, y_max, tick_step, title_suffix)
        ax.set_xlabel('Time (ms)', fontsize=16, weight='bold')
        ax.set_ylabel(f"{TARGET_COMPONENT} difference (µV)", fontsize=16, weight='bold')
        ax.set_title(f"{TARGET_COMPONENT}: Unfair - Fair (collapsed across emotion)",
                     fontsize=16, weight='bold', pad=12)

        plt.tight_layout()
        plt.savefig(save_path, dpi=300, format='tiff')
        plt.savefig(save_path.with_suffix('.png'), dpi=300, format='png')
        plt.close()

    @staticmethod
    def plot_difference_wave_by_emotion(diff_dict, time_ms, title_suffix, save_path):
        """
        Single-panel overlay of per-emotion difference waves (Unfair minus Fair).

        Visualizes the Emotion x Offer interaction directly: divergence among the
        five emotion-coded difference waves within the FRN window indicates the
        fairness-FRN effect is modulated by facial emotion. Difference-wave analogue
        of the per-emotion Fair-vs-Unfair contrast panels (plot_emotion_contrasts).
        """
        fig, ax = plt.subplots(figsize=(10, 6), facecolor='white', dpi=300)

        win_start_ms = CURR_CFG['window'][0] * 1000
        win_end_ms = CURR_CFG['window'][1] * 1000
        plot_xlim = CURR_CFG['plot_xlim']
        y_min, y_max, tick_step = VisualizationEngine._symmetric_difference_scale(
            list(diff_dict.values()))

        plt.rcParams.update({'font.size': 14, 'font.family': 'sans-serif', 'font.sans-serif': ['Arial']})

        # Fixed draw order matches the condition-mean figures for color consistency.
        draw_order = ['neu', 'aff', 'dis', 'dom', 'enj']
        for emo in draw_order:
            wave = diff_dict.get(emo)
            if wave is not None and len(wave) == len(time_ms) and not np.all(np.isnan(wave)):
                smooth = gaussian_filter1d(wave, sigma=SMOOTHING_SIGMA)
                ax.plot(time_ms, smooth, color=COLORS[emo], linestyle='-',
                        lw=LINE_WIDTH, alpha=0.9, label=LABELS_MAP[emo])

        leg = ax.legend(title='Emotion (Unfair - Fair)', loc='upper left',
                        bbox_to_anchor=(1.02, 1.0), fontsize=11, frameon=True)
        leg.get_frame().set_edgecolor('#AAAAAA')
        leg.get_frame().set_linewidth(1.0)

        VisualizationEngine._apply_standard_styling(
            ax, plot_xlim, win_start_ms, win_end_ms, y_min, y_max, tick_step, title_suffix)
        ax.set_xlabel('Time (ms)', fontsize=16, weight='bold')
        ax.set_ylabel(f"{TARGET_COMPONENT} difference (µV)", fontsize=16, weight='bold')
        ax.set_title(f"{TARGET_COMPONENT}: Unfair - Fair by emotion",
                     fontsize=16, weight='bold', pad=12)

        plt.tight_layout(rect=[0, 0, 0.85, 1.0])
        plt.savefig(save_path, dpi=300, format='tiff')
        plt.savefig(save_path.with_suffix('.png'), dpi=300, format='png')
        plt.close()


def emit_difference_waves(viz, pdata, t_ms, title_suffix, out_dir, tag):
    """Compute and plot Unfair-minus-Fair difference waves (collapsed + per-emotion).

    Shared by the Regression (Alday, PRIMARY) and Standard (sensitivity) paths so
    BOTH emit difference waves for EVERY component. `pdata` is a time-aligned dict
    keyed '<offer>_<emo>' (offer in OFFER_TYPES, emo in COLORS); `t_ms` is the time
    axis in milliseconds; `tag` labels the output files (e.g. 'Regression',
    'Standard'). Difference is computed on condition means; Gaussian smoothing is
    applied at plot time only, matching the condition-mean figures.

    Replaces the previous `if TARGET_COMPONENT == 'FRN':` gate, which never fired
    because the component names are FRN_pre / FRN_explor / LPP_pre / P3_explor
    (the bare 'FRN' name no longer exists after the SSOT split), so NO difference
    waves were emitted for any component.
    """
    draw_order = ['neu', 'aff', 'dis', 'dom', 'enj']
    t_ms = np.asarray(t_ms, dtype=float)

    diff_by_emotion = {}
    for emo in draw_order:
        u = pdata.get(f'unfair_{emo}')
        f_ = pdata.get(f'fair_{emo}')
        if u is not None and f_ is not None:
            u = np.asarray(u, dtype=float)
            f_ = np.asarray(f_, dtype=float)
            if len(u) == len(t_ms) and len(f_) == len(t_ms):
                diff_by_emotion[emo] = u - f_

    unfair_stack = [np.asarray(pdata[f'unfair_{e}'], dtype=float)
                    for e in draw_order if f'unfair_{e}' in pdata]
    fair_stack = [np.asarray(pdata[f'fair_{e}'], dtype=float)
                  for e in draw_order if f'fair_{e}' in pdata]
    if unfair_stack and fair_stack:
        diff_collapsed = (np.nanmean(np.vstack(unfair_stack), axis=0)
                          - np.nanmean(np.vstack(fair_stack), axis=0))
        viz.plot_difference_wave_collapsed(
            diff_collapsed, t_ms, title_suffix,
            out_dir / f"{TARGET_COMPONENT}_{tag}_Diff_Collapsed.tif")

    if diff_by_emotion:
        viz.plot_difference_wave_by_emotion(
            diff_by_emotion, t_ms, title_suffix,
            out_dir / f"{TARGET_COMPONENT}_{tag}_Diff_ByEmotion.tif")


# ==============================================================================
# 6. DATA INITIALIZATION & PARSING MODULE
# ==============================================================================

def prep_conditions(df):
    if 'label' in df.columns and 'offer_type' not in df.columns:
        def parse(l):
            try:
                p = l.split('_') 
                if len(p) < 4: return None, None
                emo = p[1]
                if p[2] in ['5','6']: off = 'fair'
                elif p[2] in ['8','9']: off = 'unfair'
                else: off = None
                return emo, off
            except: return None, None
        res = df['label'].apply(parse)
        df['emotion'] = [x[0] for x in res]
        df['offer_type'] = [x[1] for x in res]
        return df.dropna(subset=['offer_type', 'emotion'])
    elif 'Offers_Other' in df.columns:
        f = ((df['Offers_Other'].isin([5,6])) & (df['Offers_You'].isin([5,4])))
        u = ((df['Offers_Other'].isin([8,9])) & (df['Offers_You'].isin([2,1])))
        df['offer_type'] = np.where(f, 'fair', np.where(u, 'unfair', None))
        return df.dropna(subset=['offer_type'])
    return df

def load_data_csv(folder_or_file, desc, is_file=False, valid_keys=None):
    """
    Load epoch CSVs (folder mode) or a single ave/labelled CSV (file mode).

    v2.1 changes (folder mode only):
      - Removed unit heuristic `if max < 1.0: x *= 1e6`. The HU pipeline
        always exports voltages in microvolts; the heuristic never fires
        on real data and obscures unit-related upstream bugs.
      - Removed independent peak-to-peak rejection. Bad-epoch decisions are
        now delegated to the upstream HU pipeline (already encoded as NaN
        in trials.csv). To use this, callers pass `valid_keys` — a DataFrame
        of (participant_id, index) tuples that survived LMM-equivalent
        filtering — and load_data_csv does an inner-join on those keys.
      - When `valid_keys` is None (e.g., QC modules consuming ave.csv in
        file mode), no row-level filtering is applied; the function then
        behaves as a passive loader.

    Parameters
    ----------
    folder_or_file : Path
        Either a directory containing Vp*_epo.csv files (is_file=False) or
        a single ave/labelled CSV (is_file=True).
    desc : str
        Progress-bar description.
    is_file : bool
        File-mode toggle (True for ave.csv; False for epoch directory).
    valid_keys : pandas.DataFrame or None
        DataFrame with columns ['participant_id', 'index'] specifying which
        trials to retain. Only used in folder mode. If None, no filtering.

    Returns
    -------
    topo : DataFrame or None  (component-window mean per (subject, trial))
    wave : DataFrame or None  (full timecourse of valid trials)
    eeg  : list of str        (EEG channel names found in the data)
    """
    roi_chs = CURR_CFG['roi']
    topo, wave = [], []
    eeg = []

    if is_file:
        file_path = folder_or_file
        if not file_path.exists(): return None, None, None
        df = pd.read_csv(file_path)
        eeg = [c for c in ALL_EEG_CHANNELS if c in df.columns]
        return None, df, eeg

    folder = folder_or_file
    if not folder.exists(): return None, None, None
    files = sorted(folder.glob("Vp*_epo.csv"))
    if not files: return None, None, None
    
    temp = pd.read_csv(files[0], nrows=1)
    eeg = [c for c in ALL_EEG_CHANNELS if c in temp.columns]
    roi_avail = [c for c in roi_chs if c in temp.columns]

    # Pre-build a fast (participant_id, index) lookup if filtering requested.
    if valid_keys is not None:
        valid_set = set(zip(valid_keys['participant_id'].astype(str),
                            valid_keys['index'].astype(int)))
    else:
        valid_set = None

    for f in tqdm(files, desc=desc):
        try:
            df = pd.read_csv(f, low_memory=False)
            if df.empty: 
                print(f"   [WARN] {f.name}: empty file")
                continue
            
            # v2.1.1: participant_id is derived from filename and standardized
            # immediately. The HU pipeline uses 'VPxxxx' canonically; epoch
            # filenames may or may not match that exact case (Windows is case-
            # insensitive on disk but pandas string comparison is not).
            # Standardizing here ensures the (participant_id, index) tuple
            # matches valid_set entries (which were also standardized when
            # valid_keys was built in load_trials_csv_filtered).
            sub_id_raw = f.name.split('_')[0]
            sub_id = standardize_participant_id(sub_id_raw)
            df['participant_id'] = sub_id
            
            # Filter against valid_keys (LMM-equivalent trial set), if provided.
            # Inner-join semantics: trials present in the epoch but absent from
            # valid_keys are silently dropped. Sanity 2 (run upstream of this
            # function in RegressionAnalysis) verifies the converse direction.
            #
            # v2.1.1: vectorized membership test replaces the row-by-row
            # df.apply(lambda r: ...) which scaled poorly (~30 sec/subject for
            # large epoch files). The merge form completes in <1 sec/subject.
            if valid_set is not None:
                # Build a per-subject expected-index set in vectorized form.
                expected_idx = {idx for (pid, idx) in valid_set if pid == sub_id}
                if not expected_idx:
                    # This subject has no valid trials per LMM filter; either
                    # it has zero rows in trials.csv after filtering, or there
                    # is a participant_id encoding mismatch. The latter would
                    # already have been caught by Sanity 1 / Sanity 2, so a
                    # warning here flags the former (rare but legitimate).
                    print(f"   [WARN] {f.name}: no LMM-valid trials for {sub_id}; skipping")
                    continue
                df = df[df['index'].isin(expected_idx)]
                if df.empty:
                    print(f"   [WARN] {f.name}: filter produced empty frame")
                    continue
            
            win_min, win_max = CURR_CFG['window']
            mask = (df['time'] >= win_min) & (df['time'] <= win_max)
            
            avg = df[mask].groupby(['participant_id','index'])[eeg].mean().reset_index()
            meta_cols = ['label', 'Offers_Other', 'Offers_You', 'emotion']
            avail = [c for c in meta_cols if c in df.columns]
            meta = df[mask].groupby(['participant_id','index'])[avail].first().reset_index()
            topo.append(pd.merge(avg, meta, on=['participant_id','index']))
            
            cols_wave = ['participant_id','index','time'] + roi_avail + avail
            wave.append(df[cols_wave].copy())
        except Exception as e:
            # v2.1.1: surface error type + filename for actionable debugging.
            # Previously bare `except: pass` silently consumed all errors,
            # which made it impossible to distinguish "0 valid trials" from
            # "ROI column missing" from "file corrupt".
            print(f"   [ERROR] skipped {f.name}: {type(e).__name__}: {e}")
        
    if not topo: return None, None, None
    return pd.concat(topo), pd.concat(wave), eeg

def load_bl(path):
    """
    Load baseline scalars from baseline_values.csv (a separate HU pipeline
    export, distinct from trials.csv).

    DEPRECATED in v2.1 for the Regression pipeline:
        v2.1 RegressionAnalysis reads Baseline_<component> directly from
        trials.csv via load_trials_csv_filtered, eliminating a redundant
        merge step and a separate file dependency. This function is retained
        for backward compatibility and for ad hoc analyses that prefer the
        long-format baseline_values.csv export. Future modules should prefer
        load_trials_csv_filtered.
    """
    if not path.exists(): return None
    df = pd.read_csv(path)
    
    mapping = {c: f'Baseline_{c}' for c in BATCH_COMPONENTS}
    
    tgt = mapping.get(TARGET_COMPONENT)
    
    if not tgt or tgt not in df.columns:
        print(f"⚠️ Critical Validation Error: Baseline column '{tgt}' for {TARGET_COMPONENT} absent in parameter matrix!")
        return None

    if 'participant_id' not in df.columns and 'sub_id' in df.columns:
        df.rename(columns={'sub_id':'participant_id'}, inplace=True)
        
    print(f"   -> Parameter matrix aligned: {tgt}")
    return df[['participant_id', 'index', tgt]].rename(columns={tgt:'BL'})


# ==============================================================================
# 6.1 v2.1 HELPERS: trials.csv FILTERING, REPRODUCIBILITY PARSING, SANITY CHECKS
# ==============================================================================

def standardize_participant_id(pid):
    """
    Normalize participant_id to the canonical format 'VPxxxx' (uppercase 'VP',
    4-digit zero-padded numeric tail). Mirrors standardize_id() in the LMM R
    scripts so that trials.csv keys and epoch-CSV filename-derived keys merge
    cleanly. Idempotent on already-canonical inputs.
    """
    if pd.isna(pid):
        return pid
    s = str(pid).strip()
    m = re.search(r'(\d+)', s)
    if not m:
        return s
    return f"VP{int(m.group(1)):04d}"


def load_trials_csv_filtered(trials_path, component_name, require_baseline=True):
    """
    Load trials.csv and apply LMM-equivalent inclusion criteria.

    Filter chain (must match Sta_EEG_OfferPhase_RefAlday.Rmd v2.1
    load_and_prep, which in turn matches Sta_Behaviour_RefAlday.Rmd):
      (1) Offers_Other in INCLUSION_OFFERS  — drops 7:3 medium-fairness filler
      (2) reaction != 0                     — drops timeout / no-response
      (3) RT in [RT_LOWER_MS, RT_UPPER_MS]  — drops premature / overlong
      (4) component scalar not NaN          — drops HU bad-epoch trials
      (5) Baseline_<component> not NaN      — only when require_baseline=True

    Parameters
    ----------
    trials_path : Path
        Path to trials.csv (Method_Regression for Alday, Method_Standard for
        Traditional).
    component_name : str
        Column name of the target ERP component (e.g. 'FRN').
    require_baseline : bool
        Apply filter (5). True for Alday-framework analyses (Baseline_<comp>
        must be present), False for Traditional (no baseline column).

    Returns
    -------
    valid_keys : DataFrame  with columns ['participant_id', 'index']
    full_filtered : DataFrame  the full filtered trials.csv (for baseline use,
                                Sanity 4 cross-check, etc.)
    """
    if not trials_path.exists():
        raise FileNotFoundError(f"trials.csv not found at {trials_path}")
    df = pd.read_csv(trials_path)
    n_orig = len(df)
    
    required = ['participant_id', 'Offers_Other', 'reaction', 'RT', 'index',
                component_name]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            f"trials.csv missing required columns: {missing}. "
            f"Run upstream HU pipeline with the appropriate components spec."
        )
    
    df['participant_id'] = df['participant_id'].apply(standardize_participant_id)
    
    # Apply inclusion criteria with logging.
    df = df[df['Offers_Other'].isin(INCLUSION_OFFERS)]
    n_after_offer = len(df)
    df = df[df['reaction'] != 0]
    n_after_reaction = len(df)
    df = df[(df['RT'] >= RT_LOWER_MS) & (df['RT'] <= RT_UPPER_MS)]
    n_after_rt = len(df)
    df = df[df[component_name].notna()]
    n_after_amp = len(df)
    
    if require_baseline:
        bl_col = f'Baseline_{component_name}'
        if bl_col not in df.columns:
            raise ValueError(
                f"Baseline column '{bl_col}' missing from trials.csv. "
                f"Re-run upstream HU pipeline with baseline extraction enabled."
            )
        df = df[df[bl_col].notna()]
    n_after_bl = len(df)
    
    print(f"   [trials.csv filter] {component_name}:")
    print(f"     Original           : {n_orig}")
    print(f"     After Offers filter: {n_after_offer}  (-{n_orig - n_after_offer})")
    print(f"     After reaction !=0 : {n_after_reaction}  (-{n_after_offer - n_after_reaction})")
    print(f"     After RT in range  : {n_after_rt}  (-{n_after_reaction - n_after_rt})")
    print(f"     After amp not NaN  : {n_after_amp}  (-{n_after_rt - n_after_amp})")
    if require_baseline:
        print(f"     After BL not NaN   : {n_after_bl}  (-{n_after_amp - n_after_bl})")
    print(f"     Final N            : {n_after_bl}, Subjects: {df['participant_id'].nunique()}")
    
    valid_keys = df[['participant_id', 'index']].drop_duplicates().reset_index(drop=True)
    return valid_keys, df


def parse_reproducibility(repro_dir, component_name):
    """
    Parse <repro_dir>/<component_name>/REPRODUCIBILITY.txt and extract the
    final fixed-effect specification chosen by the upstream LMM (buildmer
    LRT-based selection in Sta_EEG_OfferPhase_RefAlday.Rmd).

    Detects baseline-condition interactions semantically (regex-based on
    'Baseline_c:emotion' or 'emotion:Baseline_c' patterns). Term ordering
    in the formula string is irrelevant.

    Returns
    -------
    spec : dict with keys
        'has_emo_BL'    : bool — whether emotion:Baseline_c was retained
        'has_offer_BL'  : bool — whether offer_type:Baseline_c was retained
        'has_emo_offer' : bool — whether emotion:offer_type was retained (sanity)
        'final_formula' : str  — the parsed final formula (for figure caption)
        'n_obs'         : int  — Observations reported by LMM
        'singular'      : bool — Singular fit flag
        'generated'     : str  — Generation timestamp from REPRODUCIBILITY.txt
        'source_path'   : str  — Path to the parsed file
    """
    repro_file = Path(repro_dir) / component_name / "REPRODUCIBILITY.txt"
    if not repro_file.exists():
        raise FileNotFoundError(
            f"REPRODUCIBILITY.txt not found at {repro_file}. "
            f"Has the LMM analysis been run for {component_name}?"
        )
    text = repro_file.read_text(encoding='utf-8')
    
    # Extract final formula block.
    m = re.search(r"## Final formula:\s*\n(.+?)(?=\n##|\Z)", text, re.DOTALL)
    if not m:
        raise ValueError(f"Cannot parse final formula in {repro_file}")
    formula_str = m.group(1).strip().replace('\n', ' ').strip()
    
    # Semantic detection of interaction terms (order-independent).
    has_emo_BL = bool(re.search(
        r"emotion\s*:\s*Baseline_c|Baseline_c\s*:\s*emotion", formula_str
    ))
    has_offer_BL = bool(re.search(
        r"offer_type\s*:\s*Baseline_c|Baseline_c\s*:\s*offer_type", formula_str
    ))
    has_emo_offer = bool(re.search(
        r"emotion\s*:\s*offer_type|offer_type\s*:\s*emotion", formula_str
    ))
    
    # Diagnostics.
    m_obs = re.search(r"Observations:\s*(\d+)", text)
    m_sing = re.search(r"Singular fit:\s*(TRUE|FALSE)", text)
    m_gen = re.search(r"Generated:\s*([^\n]+)", text)
    
    return {
        'component':      component_name,
        'final_formula':  formula_str,
        'has_emo_BL':     has_emo_BL,
        'has_offer_BL':   has_offer_BL,
        'has_emo_offer':  has_emo_offer,
        'n_obs':          int(m_obs.group(1)) if m_obs else None,
        'singular':       (m_sing.group(1) == 'TRUE') if m_sing else None,
        'generated':      m_gen.group(1).strip() if m_gen else None,
        'source_path':    str(repro_file),
    }


def run_sanity_checks(trials_path, epoch_dir, component_name, valid_keys,
                      lmm_n_expected=None, full_alignment=True):
    """
    Run alignment sanity checks before figure generation.

    Sanity 1 — participant_id consistency:
        trials.csv subjects == epoch CSV subjects (set equality).
        Catches encoding mismatches, missing subjects, extraneous files.
    Sanity 2 — per-subject trial coverage (full_alignment only):
        every (participant_id, index) in valid_keys has corresponding rows
        in the epoch CSVs. Catches trials missing from epoch data.
    Sanity 3 — post-join trial count (full_alignment only):
        len(valid_keys) == LMM Observations from REPRODUCIBILITY.txt.
        Catches drift between viz inclusion criteria and LMM filter chain.
    Sanity 4 — component scalar consistency (full_alignment only):
        deferred until after the actual epoch+window mean is computed in
        RegressionAnalysis.run; verified there before OLS fitting.

    Parameters
    ----------
    full_alignment : bool
        True for RegressionAnalysis (full 4-layer check).
        False for StandardAnalysis / QC modules (Sanity 1 only — they consume
        ave.csv directly and have no per-trial join to verify).

    Raises
    ------
    AssertionError on any failure (with a diagnostic message naming the layer).
    """
    print(f"\n   [Sanity Checks] component = {component_name}")
    
    # ---- Sanity 1 ----
    if not trials_path.exists():
        raise FileNotFoundError(f"trials.csv missing at {trials_path}")
    trials_pids = set(
        pd.read_csv(trials_path, usecols=['participant_id'])['participant_id']
        .apply(standardize_participant_id).unique()
    )
    if epoch_dir is not None and epoch_dir.exists():
        epoch_pids = set()
        for f in sorted(epoch_dir.glob("Vp*_epo.csv")):
            sub_id = standardize_participant_id(f.name.split('_')[0])
            epoch_pids.add(sub_id)
        only_in_trials = trials_pids - epoch_pids
        only_in_epoch = epoch_pids - trials_pids
        if only_in_trials or only_in_epoch:
            raise AssertionError(
                f"Sanity 1 FAIL: participant_id mismatch.\n"
                f"  In trials.csv only: {sorted(only_in_trials)}\n"
                f"  In epoch dir only:  {sorted(only_in_epoch)}"
            )
    print(f"     Sanity 1 PASS: {len(trials_pids)} subjects consistent")
    
    if not full_alignment:
        return  # Standard / QC modules stop here.
    
    # ---- Sanity 2 ----
    valid_set = set(zip(valid_keys['participant_id'].astype(str),
                        valid_keys['index'].astype(int)))
    per_sub_missing = []
    for f in sorted(epoch_dir.glob("Vp*_epo.csv")):
        sub_id = standardize_participant_id(f.name.split('_')[0])
        expected_idx = {idx for (pid, idx) in valid_set if pid == sub_id}
        if not expected_idx:
            continue
        epoch_idx = set(pd.read_csv(f, usecols=['index'])['index'].unique())
        missing = expected_idx - epoch_idx
        if missing:
            per_sub_missing.append((sub_id, len(missing), sorted(missing)[:5]))
    if per_sub_missing:
        msg = "\n".join(
            f"     {sid}: {n} missing, e.g. {sample}"
            for sid, n, sample in per_sub_missing
        )
        raise AssertionError(
            f"Sanity 2 FAIL: trials in valid_keys absent from epoch CSV.\n{msg}"
        )
    print(f"     Sanity 2 PASS: every valid trial has epoch data")
    
    # ---- Sanity 3 ----
    n_valid = len(valid_keys)
    if lmm_n_expected is not None:
        if n_valid != lmm_n_expected:
            raise AssertionError(
                f"Sanity 3 FAIL: viz N ({n_valid}) != LMM N ({lmm_n_expected}). "
                f"Inclusion criteria have drifted between viz and LMM. "
                f"Check RT_LOWER_MS, RT_UPPER_MS, INCLUSION_OFFERS, and "
                f"upstream Sta_EEG_OfferPhase_RefAlday.Rmd v2.1 load_and_prep."
            )
        print(f"     Sanity 3 PASS: viz N ({n_valid}) == LMM N ({lmm_n_expected})")
    else:
        print(f"     Sanity 3 SKIP: no LMM N reference; viz N = {n_valid}")
    # Sanity 4 is verified inside RegressionAnalysis.run after epoch loading.

# ==============================================================================
# 7. REGRESSION ANALYSIS KERNEL (ORDINARY LEAST SQUARES LMM EMULATION)
# ==============================================================================
# v2.1 changes:
#   - Trial-set selection now flows from trials.csv (load_trials_csv_filtered),
#     replacing the prior independent peak-to-peak rejection in load_data_csv.
#   - Per-component fixed-effect structure is auto-parsed from REPRODUCIBILITY.txt
#     (parse_reproducibility) so the OLS design matrix matches the upstream LMM's
#     final formula (with or without emotion:Baseline_c, offer_type:Baseline_c).
#   - Baseline is mean-centered (Alday 2019), not Z-scored.
#   - Four-layer alignment self-check runs before OLS fitting.
# ==============================================================================

class RegressionAnalysis:
    def __init__(self):
        self.out = PATHS["results"] / "Regression_Results"
        self.out.mkdir(parents=True, exist_ok=True)
        self.viz = VisualizationEngine()
        self.topo_dir = self.out / "topomap_data_for_R"
        self.topo_dir.mkdir(exist_ok=True)
        self.lmm_spec = None  # Populated by run() via parse_reproducibility.

    def run(self):
        print(f"\n🚀 Regression Engine: {TARGET_COMPONENT} (Specification: {EXPERIMENT_VERSION})")
        if not PATHS["reg_epochs"].exists():
            print(f"   [SKIP] Epoch directory missing: {PATHS['reg_epochs']}")
            return

        # ----------------------------------------------------------------------
        # Step 1: parse upstream LMM final formula for this component.
        #         Fail-soft on missing files (with prominent warning) so that
        #         viz can still produce figures with a fallback main-effect-only
        #         design matrix.
        # ----------------------------------------------------------------------
        try:
            self.lmm_spec = parse_reproducibility(
                PATHS["lmm_alday_stage1"], TARGET_COMPONENT
            )
            print(f"   [LMM auto-parse] {TARGET_COMPONENT} | "
                  f"emo:BL={self.lmm_spec['has_emo_BL']}, "
                  f"offer:BL={self.lmm_spec['has_offer_BL']}, "
                  f"N_LMM={self.lmm_spec['n_obs']}, "
                  f"generated={self.lmm_spec['generated']}")
        except FileNotFoundError as e:
            print(f"   ⚠️  [LMM auto-parse FAILED] {e}")
            print(f"   ⚠️  Falling back to main-effect-only OLS (no baseline interactions).")
            print(f"   ⚠️  Figure caption will be flagged 'NOT VALIDATED AGAINST LMM'.")
            self.lmm_spec = {
                'component':      TARGET_COMPONENT,
                'has_emo_BL':     False,
                'has_offer_BL':   False,
                'has_emo_offer':  True,
                'final_formula':  '[FALLBACK: REPRODUCIBILITY.txt unavailable]',
                'n_obs':          None,
                'singular':       None,
                'generated':      None,
                'source_path':    None,
            }
        
        # ----------------------------------------------------------------------
        # Step 2: derive valid trial keys from trials.csv via LMM-equivalent
        #         filter chain.
        # ----------------------------------------------------------------------
        try:
            valid_keys, df_trials = load_trials_csv_filtered(
                PATHS["reg_trials"], TARGET_COMPONENT, require_baseline=True
            )
        except (FileNotFoundError, ValueError) as e:
            print(f"   [HALT] Cannot load valid trial keys: {e}")
            return
        
        # ----------------------------------------------------------------------
        # Step 3: Sanity 1-3 alignment self-check.
        # ----------------------------------------------------------------------
        try:
            run_sanity_checks(
                trials_path=PATHS["reg_trials"],
                epoch_dir=PATHS["reg_epochs"],
                component_name=TARGET_COMPONENT,
                valid_keys=valid_keys,
                lmm_n_expected=self.lmm_spec['n_obs'],
                full_alignment=True,
            )
        except AssertionError as e:
            print(f"   [HALT] {e}")
            return
        
        # ----------------------------------------------------------------------
        # Step 4: load epoch data, restricted to LMM-valid trials.
        # ----------------------------------------------------------------------
        df_topo, df_wave, chans = load_data_csv(
            PATHS["reg_epochs"], "Loading epochs", is_file=False, valid_keys=valid_keys
        )
        if df_topo is None:
            # v2.1.1: surface upstream context for diagnosis. Possible causes
            # (in rough order of historical frequency):
            #   1) participant_id encoding mismatch between trials.csv and
            #      epoch filenames (now defended in load_data_csv via
            #      standardize_participant_id, but printed here in case the
            #      defense itself fails or the directory is wrong).
            #   2) Epoch files missing the expected ROI columns for this
            #      component (load_data_csv reports per-file warnings).
            #   3) trials.csv filter chain produced an empty valid set
            #      (Sanity checks should have caught this earlier).
            print(f"   [HALT] No epoch data after filtering for {TARGET_COMPONENT}.")
            print(f"          valid_keys had {len(valid_keys)} rows across "
                  f"{valid_keys['participant_id'].nunique()} subjects.")
            print(f"          Re-check: (a) epoch filename casing matches VPxxxx, "
                  f"(b) ROI channels {CURR_CFG['roi']} present in epoch CSV, "
                  f"(c) inspect [WARN]/[ERROR] lines above.")
            return
        
        # v2.1.1: standardize_participant_id is now applied inside load_data_csv
        # (so valid_set membership tests succeed). The defensive re-application
        # here is retained as a safeguard in case load_data_csv changes; it is
        # idempotent on already-canonical participant_ids.
        df_topo['participant_id']  = df_topo['participant_id'].apply(standardize_participant_id)
        df_wave['participant_id']  = df_wave['participant_id'].apply(standardize_participant_id)
        
        # ----------------------------------------------------------------------
        # Step 5: attach baseline scalar from trials.csv (NOT from baseline_values.csv).
        #         trials.csv is canonical (already filtered by HU NaN logic);
        #         mean-centering matches LMM convention (Alday 2019).
        # ----------------------------------------------------------------------
        bl_col = f'Baseline_{TARGET_COMPONENT}'
        df_bl = df_trials[['participant_id', 'index', bl_col]].rename(columns={bl_col: 'BL'})
        df_bl['participant_id'] = df_bl['participant_id'].apply(standardize_participant_id)
        
        df_fin_topo = pd.merge(df_topo, df_bl, on=['participant_id', 'index'], how='inner')
        df_fin_topo = prep_conditions(df_fin_topo)
        
        df_fin_wave = pd.merge(df_wave, df_bl, on=['participant_id', 'index'], how='inner')
        df_fin_wave = prep_conditions(df_fin_wave)
        
        # Mean-center baseline (replaces v2.0 Z-scoring). Centering preserves
        # microvolt scale and matches the upstream LMM's Baseline_c definition.
        mean_bl = df_fin_wave['BL'].mean()
        df_fin_topo['BL'] = df_fin_topo['BL'] - mean_bl
        df_fin_wave['BL'] = df_fin_wave['BL'] - mean_bl
        print(f"   [Baseline] Mean-centered: original mean = {mean_bl:.4f} uV (now 0)")
        
        # ----------------------------------------------------------------------
        # Step 6: Sanity 4 — verify component-window mean from epochs matches
        #         the precomputed scalar in trials.csv. Strongest evidence that
        #         epoch data and trials.csv reference the same preprocessing run.
        # ----------------------------------------------------------------------
        self._sanity_4_scalar_match(df_fin_topo, df_trials, chans)
        
        # ----------------------------------------------------------------------
        # Step 7: produce topographic data and ROI waveforms.
        # ----------------------------------------------------------------------
        self._process_topomaps(df_fin_topo, chans)
        
        print("   Computing Adjusted Marginal Waveforms...")
        roi_chs = [c for c in CURR_CFG['roi'] if c in df_fin_wave.columns]
        pure, times = self._perform_regression_roi(df_fin_wave, roi_chs)
        
        # ----------------------------------------------------------------------
        # Step 8: figures (unchanged formatting; caption now reflects LMM spec).
        # ----------------------------------------------------------------------
        self.viz.plot_unified_panels(
            pure, np.array(times)*1000, 
            f"Regression-Corrected", 
            self.out/f"{TARGET_COMPONENT}_Regression_Unified.tif"
        )
        self.viz.plot_emotion_contrasts(
            pure, np.array(times)*1000, 
            f"Regression-Corrected", 
            self.out/f"{TARGET_COMPONENT}_Regression_Contrasts.tif"
        )
        
        self._save_waveform_csv(pure, times)
        self._save_lmm_provenance()  # New in v2.1: caption-source provenance.

        # Difference waves on the Alday-corrected marginal means (PRIMARY path).
        # `pure` is keyed '<offer>_<emo>'; times are in seconds.
        emit_difference_waves(self.viz, pure, np.array(times) * 1000,
                              "Regression-Corrected", self.out, "Regression")

        # =======================================================================
        # POSTER VISUALIZATIONS (Main Effects Calculation & Plotting)
        # =======================================================================
        poster_dir = self.out / "teap2026poster"
        poster_dir.mkdir(parents=True, exist_ok=True)

        # 1. Calculate and Plot Fairness Main Effect (Applicable to all components)
        fairness_me_data = {'fair': [], 'unfair': []}
        for i in range(len(times)):
            for off in OFFER_TYPES:
                # Average across emotions per timepoint, safely ignoring NaNs
                vals = [pure[f'{off}_{emo}'][i] for emo in COLORS if not np.isnan(pure[f'{off}_{emo}'][i])]
                fairness_me_data[off].append(np.mean(vals) if vals else np.nan)

        self.viz.plot_fairness_main_effect(
            fairness_me_data, np.array(times)*1000,
            f"Regression-Corrected",
            poster_dir / f"{TARGET_COMPONENT}_MainEffect_Fairness.tif"
        )

        # 2. Calculate and Plot Emotion Main Effect (for EVERY component, mirroring
        #    the Fairness main-effect plot above). Previously gated by
        #    `if TARGET_COMPONENT == 'N400':` -- a stale name absent from the current
        #    component set -- so this poster plot was never produced for any component.
        emotion_me_data = {emo: [] for emo in COLORS}
        for i in range(len(times)):
            for emo in COLORS:
                # Average across offer types per timepoint, safely ignoring NaNs
                vals = [pure[f'{off}_{emo}'][i] for off in OFFER_TYPES if not np.isnan(pure[f'{off}_{emo}'][i])]
                emotion_me_data[emo].append(np.mean(vals) if vals else np.nan)

        self.viz.plot_emotion_main_effect(
            emotion_me_data, np.array(times)*1000,
            f"Regression-Corrected",
            poster_dir / f"{TARGET_COMPONENT}_MainEffect_Emotion.tif"
        )

    # --------------------------------------------------------------------------
    # Sanity 4 (Regression-mode only): epoch-window mean = trials.csv scalar
    # --------------------------------------------------------------------------
    def _sanity_4_scalar_match(self, df_topo, df_trials, channels):
        """
        Verify that recomputing the component-window mean across ROI channels
        from epoch data equals the precomputed scalar in trials.csv (within
        SCALAR_TOLERANCE_UV). Strongest evidence that epoch CSVs and trials.csv
        come from the same HU pipeline run with consistent ROI/window settings.
        """
        roi_chs = [c for c in CURR_CFG['roi'] if c in df_topo.columns]
        if not roi_chs:
            print(f"     Sanity 4 SKIP: no ROI channels in epoch frame")
            return
        df_topo = df_topo.copy()
        df_topo['epoch_window_mean'] = df_topo[roi_chs].mean(axis=1)
        check = pd.merge(
            df_topo[['participant_id', 'index', 'epoch_window_mean']],
            df_trials[['participant_id', 'index', TARGET_COMPONENT]],
            on=['participant_id', 'index'], how='inner'
        )
        if check.empty:
            print(f"     Sanity 4 SKIP: no overlapping (participant_id, index) keys")
            return
        check['diff'] = (check['epoch_window_mean'] - check[TARGET_COMPONENT]).abs()
        max_diff = check['diff'].max()
        mean_diff = check['diff'].mean()
        n_above = (check['diff'] > SCALAR_TOLERANCE_UV).sum()
        if max_diff > SCALAR_TOLERANCE_UV:
            worst = check.nlargest(3, 'diff')[
                ['participant_id', 'index', TARGET_COMPONENT, 'epoch_window_mean', 'diff']
            ]
            raise AssertionError(
                f"Sanity 4 FAIL: epoch-recomputed mean diverges from trials.csv "
                f"scalar (max diff = {max_diff:.4f} uV > tol = {SCALAR_TOLERANCE_UV} uV). "
                f"Likely cause: ROI/window mismatch between viz and HU pipeline, "
                f"or epoch CSV and trials.csv are from different runs. "
                f"Worst-aligned trials:\n{worst.to_string(index=False)}"
            )
        print(f"     Sanity 4 PASS: max diff = {max_diff:.6f} uV "
              f"(checked {len(check)} trials, {n_above} above tol)")

    # --------------------------------------------------------------------------
    # Per-component design matrix builder (v2.1)
    # --------------------------------------------------------------------------
    def _build_design_matrix(self, dt):
        """
        Construct the OLS design matrix for the current component, mirroring
        the LMM's final fixed-effect formula auto-parsed in run().

        Coding scheme (matches Sta_EEG_OfferPhase_RefAlday.Rmd contr.sum):
          - Intercept:        ones
          - Emotion (5 lvls): 4 dummy contrasts (aff/dis/dom/enj vs grand mean,
                              with neu absorbed into the intercept under the
                              -1 row of contr.sum). Implementation here uses
                              treatment dummies (emo == lvl ? 1 : 0) for the
                              4 non-reference levels; this is mathematically
                              equivalent up to a reparameterization of the
                              intercept and produces identical fitted values
                              for the marginal-mean reconstructions used here.
          - Offer (2 lvls):   sum-coded as fair=-1, unfair=+1 (oc).
          - Baseline (BL):    mean-centered scalar.
          - emo:offer:        4 columns = emo_dummy * oc, always included.
          - emo:BL:           4 columns = emo_dummy * BL, only if has_emo_BL.
          - offer:BL:         1 column  = oc * BL,        only if has_offer_BL.

        Returns
        -------
        X : 2D numpy array
        col_names : list of str
        coef_index : dict mapping logical name -> slice or int into b
        """
        spec = self.lmm_spec
        emo_levels = ['aff', 'dis', 'dom', 'enj']  # neu absorbed in intercept
        n = len(dt)
        
        cols = [np.ones(n)]
        names = ['intercept']
        idx = {'intercept': 0}
        ptr = 1
        
        # Emotion main effects.
        idx['emo_main'] = slice(ptr, ptr + 4)
        for e in emo_levels:
            cols.append((dt['emotion'] == e).astype(int).values)
            names.append(f'emo_{e}')
        ptr += 4
        
        # Offer main effect.
        oc = dt['offer_type'].map({'fair': -1, 'unfair': 1}).values.astype(float)
        cols.append(oc)
        names.append('offer')
        idx['offer'] = ptr
        ptr += 1
        
        # Baseline main effect.
        bl = dt['BL'].values.astype(float)
        cols.append(bl)
        names.append('BL')
        idx['BL'] = ptr
        ptr += 1
        
        # emotion:offer interaction (always present per LMM).
        idx['emo_offer'] = slice(ptr, ptr + 4)
        for e in emo_levels:
            cols.append((dt['emotion'] == e).astype(int).values * oc)
            names.append(f'emo_{e}_x_offer')
        ptr += 4
        
        # emotion:Baseline interaction (conditional).
        if spec['has_emo_BL']:
            idx['emo_BL'] = slice(ptr, ptr + 4)
            for e in emo_levels:
                cols.append((dt['emotion'] == e).astype(int).values * bl)
                names.append(f'emo_{e}_x_BL')
            ptr += 4
        
        # offer:Baseline interaction (conditional).
        if spec['has_offer_BL']:
            idx['offer_BL'] = ptr
            cols.append(oc * bl)
            names.append('offer_x_BL')
            ptr += 1
        
        X = np.column_stack(cols)
        return X, names, idx
    
    def _reconstruct_marginal_mean(self, b, idx, ofr, emo, emo_levels):
        """
        Compute the marginal mean for (offer, emotion) by setting BL=0 and
        evaluating the design equation with the chosen contrast values.
        Marginal means are reported at the average baseline level (BL=0 after
        mean-centering), matching emmeans(..., at = list(Baseline_c = 0)) in
        the upstream LMM.
        """
        ov = -1 if ofr == 'fair' else 1
        v = b[idx['intercept']] + b[idx['offer']] * ov
        if emo != 'neu':
            ei = emo_levels.index(emo)
            v += b[idx['emo_main']][ei]
            v += b[idx['emo_offer']][ei] * ov
        # BL terms vanish at BL=0 regardless of has_emo_BL / has_offer_BL.
        return v

    def _process_topomaps(self, df, channels):
        """Build per-channel topographic marginal means via per-component OLS."""
        src = PATHS["reg_locs"] if PATHS["reg_locs"].exists() else PATHS["std_locs"]
        if src.exists(): shutil.copy(src, self.out / "channel_locations.csv")
        
        dt = df.copy()
        X, names, idx = self._build_design_matrix(dt)
        emo_levels = ['aff', 'dis', 'dom', 'enj']
        
        rows = []
        for ch in tqdm(channels, desc="   Topographic Solvers"):
            if ch not in dt.columns: continue
            try:
                y = dt[ch].values
                # Drop rows with NaN in y for this channel.
                mask = ~np.isnan(y)
                if mask.sum() < 20: continue
                b = np.linalg.lstsq(X[mask], y[mask], rcond=None)[0]
                row = {'channel': ch}
                for ofr in OFFER_TYPES:
                    for emo in COLORS:
                        row[f'raw_{emo}_{ofr}'] = self._reconstruct_marginal_mean(
                            b, idx, ofr, emo, emo_levels
                        )
                rows.append(row)
            except Exception as e:
                print(f"   [WARN] topomap solve failed for channel {ch}: {e}")
        pd.DataFrame(rows).to_csv(self.topo_dir / f"{TARGET_COMPONENT}_topomap_data.csv", index=False)

    def _perform_regression_roi(self, df, roi_chs):
        """Per-timepoint OLS over ROI-averaged amplitude with per-component design."""
        df = df.copy()
        df['ROI'] = df[roi_chs].mean(axis=1)
        times = sorted([t for t in df['time'].unique() if PLOT_WINDOW[0] <= t <= PLOT_WINDOW[1]])
        pure = {f'{o}_{e}': [] for o in OFFER_TYPES for e in COLORS}
        emo_levels = ['aff', 'dis', 'dom', 'enj']
        
        for t in tqdm(times, desc="   Temporal Fitting"):
            dt = df[df['time'] == t].dropna(subset=['ROI', 'BL'])
            if len(dt) < 20: 
                for k in pure: pure[k].append(np.nan)
                continue
            X, names, idx = self._build_design_matrix(dt)
            y = dt['ROI'].values
            try:
                b = np.linalg.lstsq(X, y, rcond=None)[0]
                for ofr in OFFER_TYPES:
                    for emo in COLORS:
                        pure[f'{ofr}_{emo}'].append(
                            self._reconstruct_marginal_mean(b, idx, ofr, emo, emo_levels)
                        )
            except Exception:
                for k in pure: pure[k].append(np.nan)
        return pure, times

    def _save_waveform_csv(self, pure_data, times):
        df = pd.DataFrame({'time': times})
        for k, v in pure_data.items():
            if len(v) == len(times): df[k] = v
        df.to_csv(self.out / f"{TARGET_COMPONENT}_reg_waveforms.csv", index=False)
    
    def _save_lmm_provenance(self):
        """
        Write a small text file documenting which LMM REPRODUCIBILITY.txt was
        used and the parsed fixed-effect spec. Provenance for reviewers and
        for downstream scripts (e.g., topographic plotting) that consume
        topomap_data_for_R/*.csv.
        """
        spec = self.lmm_spec
        with open(self.out / f"{TARGET_COMPONENT}_lmm_provenance.txt", 'w', encoding='utf-8') as f:
            f.write(f"# LMM provenance for {TARGET_COMPONENT} ({EXPERIMENT_VERSION})\n")
            f.write(f"# Generated by Visualization_EEG_StiLocked.py v2.1 at "
                    f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"Source REPRODUCIBILITY.txt: {spec['source_path']}\n")
            f.write(f"LMM generated:              {spec['generated']}\n")
            f.write(f"LMM Observations:           {spec['n_obs']}\n")
            f.write(f"LMM Singular fit:           {spec['singular']}\n\n")
            f.write(f"Final formula (LMM):\n  {spec['final_formula']}\n\n")
            f.write(f"OLS design-matrix toggles in this viz run:\n")
            f.write(f"  has_emo_offer  : {spec['has_emo_offer']}\n")
            f.write(f"  has_emo_BL     : {spec['has_emo_BL']}\n")
            f.write(f"  has_offer_BL   : {spec['has_offer_BL']}\n")

# ==============================================================================
# 8. QC MODULE: UNCORRECTED PIPELINE REFERENCE
# ==============================================================================
# v2.1 changes:
#   - Sanity 1 (participant_id consistency) runs at the start of run().
#     The QC module consumes ave.csv directly (already condition-averaged
#     by the HU pipeline), so Sanity 2-4 are not applicable.
#   - QC_Manual_Uncorrected (the prior epoch-level recompute that re-derived
#     condition averages from epochs/Vp*_epo.csv) was REMOVED in v2.1: the
#     new Sanity 1-4 framework subsumes its sanity-check role, and viz
#     should not duplicate the HU pipeline's by-condition averaging.
# ==============================================================================

class QC_Pipeline_Uncorrected:
    def __init__(self):
        self.out = PATHS["results"] / "QC_Uncorrected"
        self.out.mkdir(parents=True, exist_ok=True)
        self.viz = VisualizationEngine()
        self.input = PATHS["reg_ave"]

    def run(self):
        print(f"\n📉 QC Pipeline Validator: {TARGET_COMPONENT}")
        if not self.input.exists():
            print(f"   [SKIP] ave.csv missing: {self.input}")
            return
        
        # Lightweight Sanity 1 only: participant_id consistency between
        # trials.csv and the epoch dir (the same dir HU averaged into ave.csv).
        # ave.csv has no per-trial keys to merge on, so Sanity 2-4 are N/A here.
        try:
            run_sanity_checks(
                trials_path=PATHS["reg_trials"],
                epoch_dir=PATHS["reg_epochs"],
                component_name=TARGET_COMPONENT,
                valid_keys=pd.DataFrame(columns=['participant_id', 'index']),
                lmm_n_expected=None,
                full_alignment=False,
            )
        except (AssertionError, FileNotFoundError) as e:
            print(f"   ⚠️  [Sanity 1 warning] {e}")
            # Do not halt QC: it's a reference visualization, not the main figure.
        
        _, df, _ = load_data_csv(self.input, "", is_file=True)
        df = prep_conditions(df)
        
        roi_chs = [c for c in CURR_CFG['roi'] if c in df.columns]
        df['mean'] = df[roi_chs].mean(axis=1)
        res = df.groupby(['offer_type','emotion','time'])['mean'].mean().reset_index()
        
        self._plot(res, "Pipeline")

    def _plot(self, res, source_name):
        plot_win = CURR_CFG['plot_xlim']
        res_plot = res[(res['time']*1000 >= plot_win[0]) & (res['time']*1000 <= plot_win[1])]
        if res_plot.empty: return
        
        pdata = {}
        t_ms = sorted(res_plot['time'].unique() * 1000)
        
        for ot in OFFER_TYPES:
            for e in COLORS:
                sub = res_plot[(res_plot['offer_type']==ot)&(res_plot['emotion']==e)].sort_values('time')
                ref_df = pd.DataFrame({'time': np.array(t_ms)/1000})
                merged = pd.merge(ref_df, sub, on='time', how='left')
                pdata[f'{ot}_{e}'] = merged['mean'].values
        
        self.viz.plot_unified_panels(
            pdata, np.array(t_ms), 
            f"Uncorrected QC ({source_name})", 
            self.out/f"{TARGET_COMPONENT}_Uncorrected_{source_name}.tif"
        )
        
        self.viz.plot_emotion_contrasts(
            pdata, np.array(t_ms), 
            f"Uncorrected QC ({source_name})", 
            self.out/f"{TARGET_COMPONENT}_Uncorrected_{source_name}_Contrasts.tif"
        )

# ==============================================================================
# 9. STANDARD ANALYSIS (SUBTRACTION BASELINE METHOD)
# ==============================================================================
# v2.1 changes:
#   - Sanity 1 (participant_id consistency) runs at the start. StandardAnalysis
#     consumes Method_Standard/ave.csv (already condition-averaged), so it
#     inherits trial-set inclusion from the upstream HU `average_by` query
#     rather than from trials.csv directly. Sanity 2-4 are N/A here.
#   - The trial set is still aligned with the Traditional LMM
#     (Sta_EEG_OfferPhase_RefTraditional.Rmd v2.1) because both consume the
#     same Method_Standard outputs filtered by the same upstream query.
# ==============================================================================

class StandardAnalysis:
    def __init__(self):
        self.out = PATHS["results"] / "Standard_Baseline"
        self.out.mkdir(parents=True, exist_ok=True)
        self.viz = VisualizationEngine()
        self.input = PATHS["std_ave"]

    def run(self):
        print(f"\n📊 Standard Processor: {TARGET_COMPONENT}")
        if not self.input.exists():
            print(f"   [SKIP] Standard ave.csv missing: {self.input}")
            return
        
        # Sanity 1 against the Standard trials.csv (the canonical post-
        # preprocessing single-trial frame for the Traditional method).
        try:
            run_sanity_checks(
                trials_path=PATHS["std_trials"],
                epoch_dir=None,  # Standard mode has no per-trial epoch consumption
                component_name=TARGET_COMPONENT,
                valid_keys=pd.DataFrame(columns=['participant_id', 'index']),
                lmm_n_expected=None,
                full_alignment=False,
            )
        except (AssertionError, FileNotFoundError) as e:
            print(f"   ⚠️  [Sanity 1 warning] {e}")
            # Continue: ave.csv may still be valid even if trials.csv inspection fails.
        
        if PATHS["std_locs"].exists(): 
            try: shutil.copy(PATHS["std_locs"], self.out / "channel_locations.csv")
            except Exception: pass
        
        _, df, _ = load_data_csv(self.input, "", is_file=True)
        df = prep_conditions(df)
        df.to_csv(self.out / f"{TARGET_COMPONENT}_ave_labelled.csv", index=False)
        
        roi_chs = [c for c in CURR_CFG['roi'] if c in df.columns]
        df['mean'] = df[roi_chs].mean(axis=1)
        res = df.groupby(['offer_type','emotion','time'])['mean'].mean().reset_index()
        
        plot_win = CURR_CFG['plot_xlim']
        res_plot = res[(res['time']*1000 >= plot_win[0]) & (res['time']*1000 <= plot_win[1])]
        
        if res_plot.empty: return
        pdata = {}
        t_ms = sorted(res_plot['time'].unique() * 1000)
        
        for ot in OFFER_TYPES:
            for e in COLORS:
                sub = res_plot[(res_plot['offer_type']==ot)&(res_plot['emotion']==e)].sort_values('time')
                ref_df = pd.DataFrame({'time': np.array(t_ms)/1000})
                merged = pd.merge(ref_df, sub, on='time', how='left')
                pdata[f'{ot}_{e}'] = merged['mean'].values

        self.viz.plot_unified_panels(
            pdata, np.array(t_ms), 
            "Standard Baseline", 
            self.out / f"{TARGET_COMPONENT}_Standard_Unified.tif"
        )
        
        self.viz.plot_emotion_contrasts(
            pdata, np.array(t_ms), 
            "Standard Baseline", 
            self.out / f"{TARGET_COMPONENT}_Standard_Contrasts.tif"
        )

        # ----------------------------------------------------------------------
        # Difference waves (Unfair minus Fair) for EVERY component, via the shared
        # helper. (Previously gated by `if TARGET_COMPONENT == 'FRN':`, which never
        # fired after the component SSOT was split into FRN_pre / FRN_explor etc.,
        # so no difference waves were ever produced.)
        # ----------------------------------------------------------------------
        emit_difference_waves(self.viz, pdata, np.array(t_ms),
                              "Standard Baseline", self.out, "Standard")


# ==============================================================================
# 10. PIPELINE ORCHESTRATION KERNEL
# ==============================================================================

def update_global_config(comp_name):
    global TARGET_COMPONENT, CURR_CFG, PLOT_WINDOW, PATHS
    TARGET_COMPONENT = comp_name
    CURR_CFG = COMPONENT_SPECS[TARGET_COMPONENT]
    PLOT_WINDOW = (CURR_CFG['plot_xlim'][0]/1000, CURR_CFG['plot_xlim'][1]/1000)
    PATHS = get_paths()
    if not PATHS["results"].exists(): 
        PATHS["results"].mkdir(parents=True, exist_ok=True)

if __name__ == "__main__":
    print(f"Visualization Execution Initiated at {datetime.now().strftime('%H:%M:%S')}")
    print(f"  EXPERIMENT_VERSION: {EXPERIMENT_VERSION}")
    print(f"  ANALYSIS_METHOD:    {ANALYSIS_METHOD}")
    print(f"  Inclusion criteria: Offers in {INCLUSION_OFFERS}, "
          f"reaction != 0, RT in [{RT_LOWER_MS}, {RT_UPPER_MS}] ms")
    
    for comp in BATCH_COMPONENTS:
        update_global_config(comp)
        
        if ANALYSIS_METHOD == 'All':
            RegressionAnalysis().run()
            StandardAnalysis().run()
            QC_Pipeline_Uncorrected().run()
            
        elif ANALYSIS_METHOD == 'Regression': 
            RegressionAnalysis().run()
        elif ANALYSIS_METHOD == 'Standard': 
            StandardAnalysis().run()
        elif ANALYSIS_METHOD == 'QC': 
            QC_Pipeline_Uncorrected().run()

    print(f"Visualization Execution Finalized at {datetime.now().strftime('%H:%M:%S')}")
