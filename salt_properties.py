"""
Units:
- Density [kg/m^3], Viscosity [Pa*s]. (can also be converted)
- Melting points are stored in deg C.
- Temperature ranges can be specified in either K or C (must match the basis)
Supported forms:
- density:  {
    "type": "linear", 
    "basis": "K" or "C",
    "units": [units],
    "a": float,
    "b": float,
    "valid_range": [T_min, T_max]  # in same basis as property
  }  => rho = a + b*T
- density:  {
    "type": "poly",
    "basis": "K" or "C",
    "units": [units],
    "coeffs": [c0, c1, c2, ...],
    "valid_range": [T_min, T_max]
  } => rho = c0 + c1*T + c2*T^2 + ...
- density:  {
    "type": "table",
    "data": [(T_K, rho), ...],
    "valid_range": [T_min, T_max]  # always in K
  } (piecewise-linear)
- viscosity: Similar to density, with same valid_range format
"""

reference_salt = "FLiNaK"

salts = {
    "FLiNaK": {
        "composition": "LiF-NaF-KF (46.5-11.5-42.0 mol%)",

        "melting_point_C": 454, 
        
        "density": {
            "type": "linear", 
            "basis": "C", 
            "units": "g/cm^3",
            "a": 2.492, 
            "b": -6.846e-4,
            "valid_range": [470, 800]  # in C
        },

        "viscosity": {
            "type": "arrhenius", 
            "basis": "K", 
            "units": "Pa*s",
            "A": 2.487e-5, 
            "E_over_R": 4478.62,
            "valid_range": [770, 970]  # in K
        },

         "sources": [
            "Engineering Database of Liquid Salt Thermophysical and Thermochemical Properties",
        ],

    },

    "NaCl_CaCl2": {
        "composition": "NaCl-CaCl2 ~49:51 mol% (near-eutectic)",

        "melting_point_C": 499.2,

        "density": {
            "type": "linear", 
            "basis": "C", 
            "units": "g/cm^3",
            "a": 2.17124, 
            "b": -0.000402,
            "valid_range": [515, 600]  # in C
        },

        "viscosity": {
            "type": "poly",
            "basis": "C",           
            "units": "cP",          
            "coeffs": [60.45986, -0.1902, 0.000158],
            "valid_range": [520, 600]  # in C
        },

        "sources": [
            "TWO NEW CHLORIDE EUTECTIC MIXTURES AND THEIR THERMO-PHYSICAL PROPERTIES FOR HIGH TEMPERATURE THERMAL ENERGY STORAGE, Pei Xie et al. (2019)",
        ],
    },

    "NaCl_KCl_CaCl2": {
        "composition": "NaCl-KCl-CaCl2 (41.7-6.1-52.2 mol%)",

        "melting_point_C": 504, 
    
        "density": {
            "type": "poly",
            "basis": "C",                 
            "units": "g/cm^3",            
            "coeffs": [2.20483, -4.41e-4], # c0 + c1*T
            "valid_range": [515, 600]  # in C
        },
        
        "viscosity": {
            "type": "poly",
            "basis": "C",                 
            "units": "cP",                
            "coeffs": [67.97221, -0.20979, 1.76e-4],
            "valid_range": [520, 600]  # in C
        },

        "sources": [
        "Energy Proceedings (Xie et al.): DSC mp, density & viscosity fits; valid ~515-600 °C",
        "Yin et al., 2022: confirms composition near 41.72-6.12-52.16 mol% and mp ≈503.8 °C"
        ],
    },
    "LiCl_KCl": {
        "composition": "LiCl-KCl (59-41 mol%)",

        "melting_point_C": 343, 
        
        "density": {
            "type": "poly", 
            "basis": "K", 
            "units": "g/cm^3",
            "coeffs": [2.0049, -0.5148e-3], # rho = c0 + c1*T
            "valid_range": [626, 1300]  # in K
        },

        "viscosity": { # valid 675-1075K
            "type": "mod_arrhenius", 
            "basis": "K", 
            "units": "cP",
            "A": 0.04, 
            "D": 21.2, 
            "T0": 185, # mu = A * exp(D*T0/T-T0)
            "valid_range": [675, 1075]  # in K
        },

         "sources": [
            "Thermophysical Properties of the LiCl-KCl Eutectic Mixture, Janz et al. (1968)",
        ],

    },
        "KCl_MgCl2": {
        "composition": "KCl-MgCl2 (68-32 mol%)",
        "melting_point_C": 424.4,  # averaged experimental mp
        "sources": [
            "Xu et al., ASME J. Solar Energy Eng. (2018): experimental Cp, vapor pressure, viscosity, density; mp 424.4 °C",
        ],
        
        # Density: q(kg/m^3) = 1903.7 - 0.552*T[°C]  (valid ~450–800 °C)
        "density": {
            "type": "poly",
            "basis": "C",
            "units": "kg/m^3",
            "coeffs": [1903.7, -0.552],
            "valid_range": [450, 800]  # in C
        },
        
        # Viscosity: μ(cP) = 14.965 - 0.0291*T + 1.784e-5*T^2   (T in °C; valid ~450–800 °C)
        "viscosity": {
            "type": "poly",
            "basis": "C",
            "units": "cP",
            "coeffs": [14.965, -0.0291, 1.784e-5],
            "valid_range": [450, 800]  # in C
        },

    },
        "NaF_ZrF4": {
        "composition": "NaF_ZrF4 (59.5-40.5 mol%)",
        "melting_point_C": 500,  # averaged experimental mp
        "sources": [
            "https://www.sciencedirect.com/science/article/pii/S0306454917301391#s0005",
        ],
        
        "density": {
            "type": "poly",
            "basis": "C",
            "units": "g/cm^3",
            "coeffs": [3.584, -0.000889],
            "valid_range": [600, 800]  # in C
        },
        
        "viscosity": {
            "type": "arrhenius", 
            "basis": "K", 
            "units": "cP",
            "A": 0.0767, 
            "E_over_R": 3977,
            "valid_range": [770, 970]  # in K
        },  
    }   

}