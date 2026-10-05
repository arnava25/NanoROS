"""Calibrated base parameters (fit to Ludke 2017, early-phase window)."""
CALIBRATED = {
    "k_in": 0.5,
    "k_cl": 0.026702308310478907,
    "k_ros": 0.5,
    "k_sod": 5.0,
    "k_gsh": 2.9999999956821464,
    "k_syn": 0.02,
    "k_deg": 0.001,
    "k_ox": 0.001,
    "k_ox_sustained": 0.0,
    "k_nrf2": 0.1,
    "k_nrf2deg": 0.01,
    "k_gshsyn": 0.1,
    "k_gshdeg": 0.02,
    "k_gshsyn_basal": 5.0,
    "k_gsh_turn": 0.001,
    "Km_gsh": 500.0,
    "k_damage_acute": 0.1,
    "k_repair_acute": 0.001,
    "k_damage_chronic": 0.0005,
    "k_chronic_gsh": 0.0002
}

# fitted: ['k_cl', 'k_gsh'] = [np.float64(0.0267), np.float64(3.0)]
# k_ros, k_in retained at literature base (not identifiable from normalized H2O2)
