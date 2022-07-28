#scil_visualize_connectivity.py avg_cm.npy CMavg_cm --display_legend

#scil_visualize_connectivity.py avg_cm.npy CMLavg_cm --display_legend

#scil_visualize_connectivity.py avg_cm.npy CMdiff_cm --display_legend

scil_visualize_connectivity.py avg_cm.npy CMavg_cm --display_legend --chord_chart CMavg_cm_chord --percentile_threshold 75 

scil_visualize_connectivity.py Lavg_cm.npy CMLavg_cm --display_legend --chord_chart CMLavg_cm_chord --percentile_threshold 75 

scil_visualize_connectivity.py diff_cm.npy CMdiff_cm --display_legend --chord_chart CMdiff_cm_chord --percentile_threshold 75
