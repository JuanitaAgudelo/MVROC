#!/usr/bin/env python3
"""
Script to extract asteroid rotation periods, amplitudes, and colors from
Greenstreet et al. (Table 2 in Asteroid_Light_Curve_Rotation.pdf) and
save to CSV format.
"""

import os
import re
import pandas as pd
import pymupdf

def extract_asteroids(pdf_path):
    doc = pymupdf.open(pdf_path)
    
    # Table 2 blocks across the two pages
    lines = (
        [l.strip() for l in doc[0].get_text('blocks')[11][4].splitlines() if l.strip()] +
        [l.strip() for l in doc[1].get_text('blocks')[14][4].splitlines() if l.strip()]
    )
    
    asteroids = []
    current = []
    for line in lines:
        if re.match(r'^2025\s+[A-Z0-9]+', line):
            if current:
                asteroids.append(current)
            current = [line]
        else:
            current.append(line)
    if current:
        asteroids.append(current)
    
    # Footnote a: associations by MPC with known objects
    mpc_known = {
        '2025 MA45': '2015 KR298',
        '2025 MD38': '2014 HR48',
        '2025 MD67': '2015 FU236',
        '2025 MD76': '2021 PA196',
        '2025 ME24': '2014 HT395',
        '2025 MH69': '2016 PW289',
        '2025 MJ71': '2022 CJ2',
        '2025 MK88': '2022 SZ331',
        '2025 ML10': '2015 KP348',
        '2025 MP47': '2015 UC70',
        '2025 MX44': '2016 UO174',
    }
    
    rows = []
    for a in asteroids:
        desig = a[0]
        items = a[1:]
        notes = ''
        if len(items) == 11:
            notes = items[0]
            items = items[1:]
        
        h_mag = float(items[0])
        mean_mag_r = float(items[1])
        n_obs = int(items[2])
        arc_days = int(items[3])
        
        p_lsm, p_four = [float(x.strip()) for x in items[4].split('/')]
        a_lsm, a_four = [float(x.strip()) for x in items[5].split('/')]
        gr_lsm, gr_four = [float(x.strip()) for x in items[6].split('/')]
        gi_lsm, gi_four = [float(x.strip()) for x in items[7].split('/')]
        ri_lsm, ri_four = [float(x.strip()) for x in items[8].split('/')]
        axial = float(items[9])
        
        is_superfast = 'b' in notes
        alt_desig = mpc_known.get(desig, '')
        
        rows.append({
            'designation': desig,
            'provid': desig,
            'alt_designation': alt_desig,
            'notes': notes,
            'superfast_rotator': is_superfast,
            'H_mag': h_mag,
            'mean_mag_r': mean_mag_r,
            'n_obs': n_obs,
            'arc_days': arc_days,
            'period_lsm_hr': p_lsm,
            'period_fourier_hr': p_four,
            'period_str': items[4],
            'amplitude_lsm': a_lsm,
            'amplitude_fourier': a_four,
            'amplitude_str': items[5],
            'g_minus_r_lsm': gr_lsm,
            'g_minus_r_fourier': gr_four,
            'g_minus_r_str': items[6],
            'g_minus_i_lsm': gi_lsm,
            'g_minus_i_fourier': gi_four,
            'g_minus_i_str': items[7],
            'r_minus_i_lsm': ri_lsm,
            'r_minus_i_fourier': ri_four,
            'r_minus_i_str': items[8],
            'axial_elongation': axial
        })
    
    return pd.DataFrame(rows)

if __name__ == '__main__':
    base_dir = os.path.dirname(os.path.abspath(__file__))
    pdf_path = os.path.join(base_dir, 'data', 'Asteroid_Light_Curve_Rotation.pdf')
    df = extract_asteroids(pdf_path)
    
    out_paths = [
        os.path.join(base_dir, 'data', 'asteroids-with-rotation.csv'),
        os.path.join(base_dir, '..', 'data', 'asteroids-with-rotation.csv'),
        os.path.join(base_dir, '..', '..', 'data', 'asteroids-with-rotation.csv'),
    ]
    for p in out_paths:
        os.makedirs(os.path.dirname(p), exist_ok=True)
        df.to_csv(p, index=False)
        print(f"Saved {len(df)} records to {p}")
