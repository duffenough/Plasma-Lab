#!/usr/bin/env python3
"""
Compute a unitless similarity score between a surrogate salt and a reference salt (default FLiNaK).

Score:
  S = exp( - sqrt( w_Re * [ln R_Re]^2 + w_mu * [ln R_mu]^2 + w_rho * [ln R_rho]^2 + w_mp * [ln R_mp]^2 ) )

Where:
  R_Re  = Re_surrogate / Re_reference (or 1 if --match_Re is used)
  R_mu  = mu_surrogate / mu_reference
  R_rho = rho_surrogate / rho_reference
  R_mp  = (T_melt_surrogate[K]) / (T_melt_reference[K])

Inputs:
  - Data are imported from salt_properties.py (pure Python dicts).
  - Temperature T for property evaluation is in Kelvin.
  - U [m/s], L [m] define the Reynolds number scale (nozzle or droplet convention).

Usage examples:
  python salt_similarity.py --T 900
  python salt_similarity.py --surrogate NaCl_CaCl2_eutectic --match_Re
  python salt_similarity.py --weights '{"Re":0.5,"mu":0.2,"rho":0.2,"mp":0.1}'
"""

import argparse, json, math, importlib
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from typing import Dict, Any, List, Tuple

props = importlib.import_module("salt_properties")

def _to_celsius_if_needed(basis: str, T_K: float) -> float:
    """Return temperature in the basis requested ('K' or 'C')."""
    b = (basis or "K").upper()
    if b == "K":
        return T_K
    if b == "C":
        return T_K - 273.15
    raise ValueError(f"Unknown temperature basis: {basis}")

def _density_to_SI(value: float, units: str | None) -> float:
    """
    Convert density to kg/m^3.
    Supported unit strings (case-insensitive, spaces/dots/slashes flexible):
      - 'kg/m^3' (default if None) -> kg/m^3
      - 'g/cm^3', 'g/cc', 'g/mL'   -> x1000 to kg/m^3
    """
    if units is None:
        return value
    u = units.lower().replace(" ", "").replace("·", "").replace("\\", "/")
    if u in ("kg/m^3", "kg/m3"):
        return value
    if u in ("g/cm^3", "g/cc", "g/ml", "g/mL".lower()):
        return value * 1000.0
    # Add other aliases as needed
    raise ValueError(f"Unknown density units: {units}")

def _viscosity_to_SI(value: float, units: str | None) -> float:
    """
    Convert viscosity to Pa*s.
    Supported:
      - 'pa*s' (default if None) -> Pa*s
      - 'mpa*s'                  -> x1e-3 Pa*s
      - 'cp', 'centipoise'       -> x1e-3 Pa*s
      - 'p', 'poise'             -> x0.1   Pa*s
    """
    if units is None:
        return value
    u = units.lower().replace(" ", "").replace("·", "")
    if u in ("pa*s", "pas"):
        return value
    if u in ("mpa*s", "mpas"):
        return value * 1e-3
    if u in ("cp", "centipoise"):
        return value * 1e-3
    if u in ("p", "poise"):
        return value * 0.1
    raise ValueError(f"Unknown viscosity units: {units}")

def _interp_table(table: List[Tuple[float, float]], T: float) -> float:
    if not table:
        raise ValueError("Empty table for property")
    table = sorted(table, key=lambda x: x[0])
    if T <= table[0][0]:
        return table[0][1]
    if T >= table[-1][0]:
        return table[-1][1]
    for i in range(len(table)-1):
        T0, v0 = table[i]
        T1, v1 = table[i+1]
        if T0 <= T <= T1:
            a = (T - T0) / (T1 - T0)
            return v0 + a * (v1 - v0)
    return table[-1][1]

def _poly_value(T: float, coeffs):
    """Evaluate polynomial c0 + c1*T + c2*T^2 + ..."""
    v = 0.0
    p = 1.0
    for c in coeffs:
        v += c * p
        p *= T
    return v

def _check_validity(prop_dict: Dict[str, Any], T_K: float) -> bool:
    """Check if temperature is within valid range for a property."""
    if "valid_range" not in prop_dict:
        return True  # No range specified = assume valid

    T_min, T_max = prop_dict["valid_range"]
    basis = prop_dict.get("basis", "K")

    # Convert stored range to Kelvin for comparison
    if basis.upper() == "C":
        T_min_K = T_min + 273.15
        T_max_K = T_max + 273.15
    else:
        T_min_K = T_min
        T_max_K = T_max

    # T_K is already in Kelvin
    return (T_K >= T_min_K) and (T_K <= T_max_K)

def _density(salt: Dict[str, Any], T_K: float, warn: bool = True) -> float:
    """Return density in kg/m^3."""
    d = salt["density"]
    units = d.get("units", "kg/m^3")
    ttype = d["type"].lower()

    # Check validity range
    is_valid = _check_validity(d, T_K)
    if not is_valid and warn:
        # Warn but continue (user requested not to stop execution)
        basis = d.get("basis", "K")
        T_display = _to_celsius_if_needed(basis, T_K)
        vr = d.get("valid_range")
        if vr is not None:
            warnings.warn(
                f"Temperature {T_display:.1f} {basis} is outside the valid range {vr} ({basis}) for density; proceeding with calculation.",
                stacklevel=2,
            )
        else:
            warnings.warn(
                f"Temperature {T_display:.1f} {basis} may be outside documented valid range for density; proceeding with calculation.",
                stacklevel=2,
            )

    if ttype == "linear":
        # rho = a + b*T  (T in K or C per 'basis')
        T = _to_celsius_if_needed(d.get("basis", "K"), T_K)
        raw = d["a"] + d["b"] * T
        return _density_to_SI(raw, units)

    elif ttype == "poly":
        # rho = poly(T)  (T in K or C per 'basis')
        T = _to_celsius_if_needed(d.get("basis", "K"), T_K)
        raw = _poly_value(T, d["coeffs"])
        return _density_to_SI(raw, units)

    elif ttype == "table":
        # table expects temperatures in K for the first column;
        # values are in 'units' for density-level dict
        raw = _interp_table(d["data"], T_K)
        return _density_to_SI(raw, units)

    else:
        raise ValueError(f"Unknown density type: {d['type']}")

def _viscosity(salt: Dict[str, Any], T_K: float, warn: bool = True) -> float:
    """Return viscosity in Pa*s."""
    v = salt["viscosity"]
    units = v.get("units", "Pa*s")
    ttype = v["type"].lower()

    # Check validity range
    is_valid = _check_validity(v, T_K)
    if not is_valid and warn:
        # Warn but continue
        basis = v.get("basis", "K")
        T_display = _to_celsius_if_needed(basis, T_K)
        vr = v.get("valid_range")
        if vr is not None:
            warnings.warn(
                f"Temperature {T_display:.1f} {basis} is outside the valid range {vr} ({basis}) for viscosity; proceeding with calculation.",
                stacklevel=2,
            )
        else:
            warnings.warn(
                f"Temperature {T_display:.1f} {basis} may be outside documented valid range for viscosity; proceeding with calculation.",
                stacklevel=2,
            )

    if ttype == "arrhenius":
        # mu = A * exp(E_over_R / T)  (T in K)
        raw = v["A"] * math.exp(v["E_over_R"] / T_K)
        return _viscosity_to_SI(raw, units)
    
    elif ttype == "mod_arrhenius":
        # mu = A * exp(D*T0/T-T0)  (T in K)
        raw = v["A"] * math.exp((v["D"] * v["T0"]) / (T_K - v["T0"]))
        return _viscosity_to_SI(raw, units)

    elif ttype == "poly":
        # mu = poly(T)  (T in K or C per 'basis'), polynomial returns in 'units'
        T = _to_celsius_if_needed(v.get("basis", "K"), T_K)
        raw = _poly_value(T, v["coeffs"])
        return _viscosity_to_SI(raw, units)

    elif ttype == "table":
        # table temperatures are in K; values in 'units' for viscosity-level dict
        raw = _interp_table(v["data"], T_K)
        return _viscosity_to_SI(raw, units)

    else:
        raise ValueError(f"Unknown viscosity type: {v['type']}")

def _Re(rho: float, mu: float, U: float, L: float) -> float:
    return rho * U * L / mu

def _ratio(a: float, b: float) -> float:
    return a / b

def _similarity(R_Re: float, R_mu: float, R_rho: float, R_mp: float,
                w_Re=0.5, w_mu=0.2, w_rho=0.2, w_mp=0.1) -> float:
    terms = [
        w_Re  * (math.log(R_Re))**2,
        w_mu  * (math.log(R_mu))**2,
        w_rho * (math.log(R_rho))**2,
        w_mp  * (math.log(R_mp))**2,
    ]
    return math.exp(-math.sqrt(sum(terms)))

def get_validity_ranges(salt: Dict[str, Any]) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Get validity ranges for density and viscosity in Kelvin."""
    density_range = None
    viscosity_range = None
    
    # Get density range
    if "density" in salt and "valid_range" in salt["density"]:
        try:
            T_min, T_max = salt["density"]["valid_range"]
            basis = salt["density"].get("basis", "K").upper()
            
            # Convert to Kelvin if needed
            if basis == "C":
                density_range = (T_min + 273.15, T_max + 273.15)
            else:
                density_range = (T_min, T_max)
            
            # Print debug info
            print(f"DEBUG: density range for {salt.get('composition', 'unknown')}: {T_min}-{T_max} {basis} -> {density_range[0]:.1f}-{density_range[1]:.1f}K")
        except Exception as e:
            density_range = None
            print(f"Error processing density range: {e}")
    else:
        print(f"DEBUG: No density range found for {salt.get('composition', 'unknown')}")
    
    # Get viscosity range
    if "viscosity" in salt and "valid_range" in salt["viscosity"]:
        try:
            T_min, T_max = salt["viscosity"]["valid_range"]
            basis = salt["viscosity"].get("basis", "K").upper()
            
            # Convert to Kelvin if needed
            if basis == "C":
                viscosity_range = (T_min + 273.15, T_max + 273.15)
            else:
                viscosity_range = (T_min, T_max)
            
            # Print debug info
            print(f"DEBUG: viscosity range for {salt.get('composition', 'unknown')}: {T_min}-{T_max} {basis} -> {viscosity_range[0]:.1f}-{viscosity_range[1]:.1f}K")
        except Exception as e:
            viscosity_range = None
            print(f"Error processing viscosity range: {e}")
    else:
        print(f"DEBUG: No viscosity range found for {salt.get('composition', 'unknown')}")
    
    # Additional debug output for verification
    print(f"DEBUG: Final ranges for {salt.get('composition', 'unknown')}:")
    print(f"       Density range (K): {density_range}")
    print(f"       Viscosity range (K): {viscosity_range}")
    
    return density_range, viscosity_range

def compute_similarity_vs_temperature(surrogate: str,
                                  reference: str = None,
                                  T_min: float = 600.0,
                                  T_max: float = 1150.0,
                                  n_points: int = 100,
                                  U: float = 1.0,
                                  L: float = 1e-3,
                                  match_Re: bool = False,
                                  weights: Tuple[float, float, float, float] = (0.5, 0.2, 0.2, 0.1)) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """Calculate similarity score S over a range of temperatures.
    
    Returns:
        Tuple containing:
        - Temperature points
        - Similarity scores
        - Dictionary with validity ranges for reference and surrogate salts
    """
    T_range = np.linspace(T_min, T_max, n_points)
    S_values = []
    
    # Get salt data
    if reference is None:
        reference = props.reference_salt
    ref_salt = props.salts[reference]
    surr_salt = props.salts[surrogate]
    
    # Get validity ranges
    print(f"\nProcessing reference salt: {reference}")
    ref_density_range, ref_visc_range = get_validity_ranges(ref_salt)
    
    print(f"\nProcessing surrogate salt: {surrogate}")
    surr_density_range, surr_visc_range = get_validity_ranges(surr_salt)
    
    print("\nValidity ranges (K):")
    print(f"Reference: density {ref_density_range}, viscosity {ref_visc_range}")
    print(f"Surrogate: density {surr_density_range}, viscosity {surr_visc_range}")
    
    validity_info = {
        "reference": {
            "density": ref_density_range,
            "viscosity": ref_visc_range
        },
        "surrogate": {
            "density": surr_density_range,
            "viscosity": surr_visc_range
        }
    }
    
    # Calculate similarity scores
    for T in T_range:
        # suppress warnings during automated sweep; only the single user T should warn
        result = compute_similarity(surrogate, reference, T, U, L, match_Re, weights, warn=False)
        S_values.append(result["S"])
    
    return T_range, np.array(S_values), validity_info

def compute_similarity(surrogate: str,
                       reference: str = None,
                       T_K: float = 900.0,
                       U: float = 1.0,
                       L: float = 1e-3,
                       match_Re: bool = False,
                       weights: Tuple[float, float, float, float] = (0.5, 0.2, 0.2, 0.1),
                       warn: bool = True) -> Dict[str, Any]:

    if reference is None:
        reference = props.reference_salt

    salts = props.salts
    if reference not in salts:
        raise KeyError(f"Reference salt '{reference}' not found in salt_properties.salts.")
    if surrogate not in salts:
        raise KeyError(f"Surrogate salt '{surrogate}' not found in salt_properties.salts.")

    ref = salts[reference]
    surr = salts[surrogate]

    rho_ref = _density(ref, T_K, warn=warn)
    mu_ref  = _viscosity(ref, T_K, warn=warn)
    rho_s   = _density(surr, T_K, warn=warn)
    mu_s    = _viscosity(surr, T_K, warn=warn)

    Re_ref = _Re(rho_ref, mu_ref, U, L)
    Re_s   = _Re(rho_s,   mu_s,   U, L)

    R_rho = _ratio(rho_s, rho_ref)
    R_mu  = _ratio(mu_s,  mu_ref)
    R_mp  = _ratio(surr["melting_point_C"] + 273.15, ref["melting_point_C"] + 273.15)
    R_Re  = 1.0 if match_Re else _ratio(Re_s, Re_ref)

    w_Re, w_mu, w_rho, w_mp = weights
    S = _similarity(R_Re, R_mu, R_rho, R_mp, w_Re=w_Re, w_mu=w_mu, w_rho=w_rho, w_mp=w_mp)

    return {
        "reference": reference,
        "surrogate": surrogate,
        "T_K": T_K,
        "U": U,
        "L": L,
        "match_Re": match_Re,
        "weights": {"Re": w_Re, "mu": w_mu, "rho": w_rho, "mp": w_mp},
        "rho_ref": rho_ref,
        "mu_ref": mu_ref,
        "rho_surr": rho_s,
        "mu_surr": mu_s,
        "Re_ref": Re_ref,
        "Re_surr": Re_s,
        "R_rho": R_rho,
        "R_mu": R_mu,
        "R_Re": R_Re,
        "R_mp": R_mp,
        "S": S
    }

def main():
    ap = argparse.ArgumentParser(description="Similarity between surrogate and reference salts.")
    ap.add_argument("--surrogate", default="None")
    ap.add_argument("--reference", default=None)
    ap.add_argument("--T", type=float, default=900.0, help="Temperature in K")
    ap.add_argument("--T_min", type=float, default=600.0, help="Minimum temperature for plot in K")
    ap.add_argument("--T_max", type=float, default=1150, help="Maximum temperature for plot in K")
    ap.add_argument("--U", type=float, default=1.0, help="Characteristic velocity [m/s]")
    ap.add_argument("--L", type=float, default=1e-3, help="Characteristic length [m]")
    ap.add_argument("--match_Re", action="store_true", help="Force R_Re=1 (assume you retune U to match Re)")
    ap.add_argument("--weights", type=str, default=None, help='JSON like {"Re":0.5,"mu":0.2,"rho":0.2,"mp":0.1}')
    ap.add_argument("--no-plot", action="store_true", help="Skip plotting temperature dependence")
    args = ap.parse_args()

    weights = (0.5, 0.2, 0.2, 0.1)
    if args.weights:
        jw = json.loads(args.weights)
        weights = (jw.get("Re",0.5), jw.get("mu",0.2), jw.get("rho",0.2), jw.get("mp",0.1))

    # Calculate similarity at the specified temperature
    res = compute_similarity(args.surrogate, args.reference, args.T, args.U, args.L, args.match_Re, weights)
    print(json.dumps(res, indent=2))

    # Plot temperature dependence unless --no-plot is specified
    if not args.no_plot:
        T_range, S_values, validity = compute_similarity_vs_temperature(
            args.surrogate, args.reference,
            T_min=args.T_min, T_max=args.T_max,
            U=args.U, L=args.L,
            match_Re=args.match_Re,
            weights=weights
        )
        
        plt.figure(figsize=(10, 8))
        
        # Create main plot
        plt.subplot(2, 1, 1)
        # Convert temperature axis from K -> C for plotting
        T_plot = T_range - 273.15
        plt.plot(T_plot, S_values, 'b-', linewidth=2)
        plt.grid(True)
        plt.xlabel('Temperature (°C)')
        plt.ylabel('Similarity Score S')
        
        ref_name = args.reference if args.reference else props.reference_salt
        plt.title(f'Temperature Dependence of Similarity Score\n{args.surrogate} vs {ref_name}')
        
    # Note: we no longer plot the fixed user T point here. Instead
    # we'll identify the best temperature within the intersection
    # of validity ranges (or fall back to the global best) and
    # plot that marker below.
        
        # Print ranges to debug
        print("\nDebug - Validity Ranges (K):")
        for salt_name, ranges in [(ref_name, validity["reference"]), (args.surrogate, validity["surrogate"])]:
            print(f"{salt_name}:")
            print(f"  density: {ranges['density']}")
            print(f"  viscosity: {ranges['viscosity']}")
        
        # Add validity ranges as full-height, low-alpha rectangles with contrasting
        # edge colors and hatch patterns so overlaps remain distinguishable.
        # Draw them behind the curve (low zorder) so the curve stays readable.
        alpha = 0.18
        handles = []

        styles = [
            ("tab:green", "darkgreen", "///", f"{ref_name} density valid"),
            ("tab:blue", "navy", "\\\\", f"{ref_name} viscosity valid"),
            ("tab:red", "darkred", "xxx", f"{args.surrogate} density valid"),
            ("tab:orange", "darkorange", "ooo", f"{args.surrogate} viscosity valid"),
        ]

        # reference density (convert K->C for plotting)
        if validity["reference"]["density"]:
            T1_k, T2_k = validity["reference"]["density"]
            T1, T2 = T1_k - 273.15, T2_k - 273.15
            face, edge, hatch, label = styles[0]
            plt.axvspan(T1, T2, color=face, alpha=alpha, edgecolor=edge, linewidth=1.2, hatch=hatch, zorder=0)
            handles.append(mpatches.Patch(facecolor=face, edgecolor=edge, hatch=hatch, label=label))

        # reference viscosity
        if validity["reference"]["viscosity"]:
            T1_k, T2_k = validity["reference"]["viscosity"]
            T1, T2 = T1_k - 273.15, T2_k - 273.15
            face, edge, hatch, label = styles[1]
            plt.axvspan(T1, T2, color=face, alpha=alpha, edgecolor=edge, linewidth=1.2, hatch=hatch, zorder=0)
            handles.append(mpatches.Patch(facecolor=face, edgecolor=edge, hatch=hatch, label=label))

        # surrogate density
        if validity["surrogate"]["density"]:
            T1_k, T2_k = validity["surrogate"]["density"]
            T1, T2 = T1_k - 273.15, T2_k - 273.15
            face, edge, hatch, label = styles[2]
            plt.axvspan(T1, T2, color=face, alpha=alpha, edgecolor=edge, linewidth=1.2, hatch=hatch, zorder=0)
            handles.append(mpatches.Patch(facecolor=face, edgecolor=edge, hatch=hatch, label=label))

        # surrogate viscosity
        if validity["surrogate"]["viscosity"]:
            T1_k, T2_k = validity["surrogate"]["viscosity"]
            T1, T2 = T1_k - 273.15, T2_k - 273.15
            face, edge, hatch, label = styles[3]
            plt.axvspan(T1, T2, color=face, alpha=alpha, edgecolor=edge, linewidth=1.2, hatch=hatch, zorder=0)
            handles.append(mpatches.Patch(facecolor=face, edgecolor=edge, hatch=hatch, label=label))

        # Compute intersection of all valid ranges (if all present) and display it
        all_ranges = [
            validity["reference"]["density"],
            validity["reference"]["viscosity"],
            validity["surrogate"]["density"],
            validity["surrogate"]["viscosity"],
        ]
        present = [r for r in all_ranges if r is not None]
        intersection_range = None
        if present:
            low_k = max(r[0] for r in present)
            high_k = min(r[1] for r in present)
            if low_k <= high_k:
                intersection_range = (low_k, high_k)
                # convert to C for plotting
                low_c = low_k - 273.15
                high_c = high_k - 273.15
                # vertical dashed lines at the boundaries and a faint shaded band
                plt.axvline(low_c, color='k', linewidth=1.5, linestyle='--', zorder=1)
                plt.axvline(high_c, color='k', linewidth=1.5, linestyle='--', zorder=1)
                plt.axvspan(low_c, high_c, color='k', alpha=0.06, zorder=0)
                handles.append(Line2D([0], [0], color='k', lw=1.5, linestyle='--', label='All models valid'))

        # Find the temperature with highest similarity within the intersection
        # If there is no intersection or no sampled points inside it, fall
        # back to the global best over the plotted T_range.
        global_best_idx = int(np.argmax(S_values))
        global_best_T = T_range[global_best_idx]
        global_best_S = S_values[global_best_idx]

        best_idx = global_best_idx
        best_T = global_best_T
        best_S = global_best_S
        best_label = f'Best overall: {best_T:.0f} K'

        if intersection_range is not None:
            low, high = intersection_range
            mask = (T_range >= low) & (T_range <= high)
            if mask.any():
                idxs = np.where(mask)[0]
                sub = S_values[idxs]
                sub_best = int(np.argmax(sub))
                best_idx = idxs[sub_best]
                best_T = T_range[best_idx]
                best_S = S_values[best_idx]
                best_label = f'Best within intersection: {best_T-273.15:.0f} °C'

        # Plot the best-point marker (convert to C) and add it to the legend handles
        best_T_plot = best_T - 273.15
        plt.plot(best_T_plot, best_S, 'ro', markersize=8, zorder=5)
        handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor='r', markersize=8, label=best_label))

        # Annotate the main axes with intersection summary and best S
        ax = plt.gca()
        if intersection_range is not None:
            ax.text(0.02, 0.95, f'Intersection: {intersection_range[0]-273.15:.0f} - {intersection_range[1]-273.15:.0f} °C', transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7))
            ax.text(0.02, 0.88, f'Best S = {best_S:.3f} at {best_T-273.15:.0f} °C', transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7))
        else:
            ax.text(0.02, 0.95, 'Intersection: None', transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7))
            ax.text(0.02, 0.88, f'Best S = {best_S:.3f} at {best_T-273.15:.0f} °C', transform=ax.transAxes, fontsize=10, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7))

        plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Add validity range summary
        plt.subplot(2, 1, 2)
        plt.axis('off')
        summary = "Validity Ranges (K):\n\n"
        
        for salt_name, ranges in [
            (ref_name, validity["reference"]),
            (args.surrogate, validity["surrogate"])
        ]:
            summary += f"{salt_name}:\n"
            if ranges["density"]:
                T1, T2 = ranges["density"]
                summary += f"  Density: {T1:.0f} - {T2:.0f} K\n"
            if ranges["viscosity"]:
                T1, T2 = ranges["viscosity"]
                summary += f"  Viscosity: {T1:.0f} - {T2:.0f} K\n"
            summary += "\n"
        
        plt.text(0, 0.5, summary, fontsize=10, fontfamily='monospace')
        
        plt.tight_layout()
        plt.show()

if __name__ == "__main__":
    main()
