"""Descriptions of the columns of the parquet files written by the notebooks.

Each column gets its description in the parquet schema (field metadata 'description'), and each file a short
description of its content (schema metadata 'description'). They are written with the data, so they travel with
the file; pandas does not show them, read them with `column_descriptions(path)`.

Writing a file with descriptions:
    write_parquet(df, 'data/obs_elements_state.parquet')             # a DataFrame
    writer.write_table(describe(table, 'obs_elements_state_at_obstime.parquet'))   # an Arrow table (chunked writing)
"""
import os

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

BANDS = ['u', 'g', 'r', 'i', 'z', 'y']

# ------------------------------------------------------------------------------------------------ object identifiers
IDENTIFIERS = {
    'mpc_designation': "Designation of the object in the MPC orbit catalogue (MPCORB readable designation, e.g. "
                       "'(1) Ceres', '(831800) 2009 WC150' or '2021 AQ'). It is the same for all the detections of an "
                       "object, so use it to group detections by object: one object can have several provid "
                       "(see n_provids). Null if the object was not found (comets, or permid and provid both null).",
    'provid': "Provisional designation with which Rubin reported the detection (e.g. '2025 JQ7'). It is NOT a unique "
              "object id: the MPC may later link two provids to the same object (e.g. '2025 JQ7' and '2025 NF311' are "
              "both '2025 NF311'), so to group by object use mpc_designation. See n_provids.",
    'permid': "Permanent number of the object (numbered minor planets, e.g. '831800'; comets e.g. '228P'), as "
              "reported with the detection. Null for unnumbered objects.",
    'n_provids': "Number of different Rubin provids of the object (same mpc_designation) in the Rubin catalogue: "
                 "1 = one provid; 2 = Rubin reported the object as two separate objects, with two provids that the "
                 "MPC later linked (e.g. '2025 JQ7' and '2025 NF311' are both '2025 NF311'); 0 = reported only by "
                 "its number (permid). The list of provids of each object is in data/mpc_designations.csv. "
                 "Null if mpc_designation is null.",
    'match_by': "How mpc_designation was found (add_mpc_designation.ipynb): 'permid' = by number in MPCORB; "
                "'provid' = provid is the main designation in MPCORB; 'other_desig' = provid is one of the other "
                "designations of the object in MPCORB; 'mpc_api' = with the MPC Designation Identifier API "
                "(objects linked after the MPCORB file was downloaded); null = not found.",
    'mpc_row': "Row (0-based) of the object in MPCORB (data/mpcorb_extended.dat.gz) as read in exploration_data.ipynb. "
               "One value per object, like mpc_designation.",
}

# ----------------------------------------------- columns of the MPC observations table (obs_sbn, ADES format)
# Columns that are always null for the Rubin (X05) observations
SBN_EMPTY = ['artsat', 'trx', 'rcv', 'sys', 'ctr', 'pos1', 'pos2', 'pos3', 'poscov11', 'poscov12', 'poscov13',
             'poscov22', 'poscov23', 'poscov33', 'rastar', 'decstar', 'obscenter', 'deltara', 'deltadec', 'dist', 'pa',
             'rmsdist', 'rmspa', 'delay', 'rmsdelay', 'doppler', 'rmsdoppler', 'photap', 'nucmag', 'logsnr', 'seeing',
             'exp', 'rmsfit', 'com', 'frq', 'disc', 'subfrm', 'prectime', 'precra', 'precdec', 'unctime', 'notes',
             'remarks', 'deprecated', 'localuse', 'nstars', 'prev_desig', 'prev_ref', 'rmstime', 'orbit_id',
             'designation_asterisk', 'shapeocc', 'replacesobsid', 'vel1', 'vel2', 'vel3', 'fltr']

SBN = {
    'id': "Row id of the observation in the MPC observations table (obs_sbn).",
    'obsid': "Unique observation id assigned by the MPC.",
    'obssubid': "Observation id given by the submitter (Rubin source id, e.g. 'LSST-DP2-DS-...').",
    'trksub': "Tracklet id given by the submitter (Rubin): detections of the same object in one night share it.",
    'trkid': "Tracklet id assigned by the MPC.",
    'trkmpc': "ADES trkMPC field (MPC tracklet id); almost always null for Rubin.",
    'group_id': "Group id of the observation in the MPC database.",
    'submission_id': "Id of the MPC submission (batch of observations) that contained the observation.",
    'submission_block_id': "Id of the block, inside the MPC submission, that contained the observation.",
    'obs80': "The observation in the old MPC 80-column format.",
    'status': "MPC processing status code of the observation.",
    'ref': "Publication reference of the observation (e.g. 'MPEC 2026-Q94'); null if not published yet.",
    'healpix': "HEALPix pixel index of the sky position (computed by the MPC).",
    'mode': "Observation mode ('CCD').",
    'stn': "MPC observatory code ('X05' = Vera C. Rubin Observatory).",
    'prog': "MPC program code of the observation.",
    'obstime': "Time of the observation, UTC.",
    'obstime_text': "Time of the observation as ISO text, UTC (same as obstime).",
    'ra': "Right ascension [deg], ICRS, as text.",
    'dec': "Declination [deg], ICRS, as text.",
    'rmsra': "Uncertainty of the right ascension (times cos(dec)) [arcsec], as text.",
    'rmsdec': "Uncertainty of the declination [arcsec], as text.",
    'rmscorr': "Correlation between the right ascension and declination uncertainties, as text.",
    'astcat': "Astrometric reference catalogue ('Gaia3').",
    'mag': "Apparent magnitude reported by Rubin, in the band of the detection (column band), as text.",
    'rmsmag': "Uncertainty of mag [mag], as text.",
    'band': "Photometric band of mag as reported (e.g. 'Lr' or 'r' = LSST r).",
    'photcat': "Photometric reference catalogue.",
    'subfmt': "Format in which the observation was submitted to the MPC.",
    'created_at': "Time when the observation was added to the MPC database.",
    'updated_at': "Time when the observation was last updated in the MPC database.",
    **{c: "ADES field not used for the Rubin observations: always null in this file." for c in SBN_EMPTY},
}

# ------------------------------------------------------------------------- MPC orbit (MPCORB), at the MPC epoch
MPCORB = {
    'mpc_a': "Semimajor axis from MPCORB [AU], at the epoch of the MPC elements (mpc_epoch_mjd).",
    'mpc_e': "Eccentricity from MPCORB, at the epoch of the MPC elements.",
    'mpc_incl': "Inclination from MPCORB [deg], ecliptic J2000, at the epoch of the MPC elements.",
    'mpc_node': "Longitude of the ascending node from MPCORB [deg], ecliptic J2000, at the epoch of the MPC elements.",
    'mpc_peri': "Argument of perihelion from MPCORB [deg], ecliptic J2000, at the epoch of the MPC elements.",
    'mpc_M': "Mean anomaly from MPCORB [deg], at the epoch of the MPC elements.",
    'mpc_q': "Perihelion distance [AU], mpc_a * (1 - mpc_e).",
    'mpc_epoch_mjd': "Epoch of the MPC elements [MJD, TT].",
    'mpc_H': "Absolute magnitude H from MPCORB [mag].",
    'mpc_tp_jd': "Time of perihelion from MPCORB extended [JD].",
    'epoch_et': "Epoch of the MPC elements as ephemeris time [TDB seconds past J2000].",
}

# ------------------------------------------------------------ propagated orbit and geometry, at the observation time
STATE = ['x', 'y', 'z', 'vx', 'vy', 'vz']
STATE_UNITS = {'x': 'position x [km]', 'y': 'position y [km]', 'z': 'position z [km]',
               'vx': 'velocity x [km/s]', 'vy': 'velocity y [km/s]', 'vz': 'velocity z [km/s]'}
STATE_AT_OBSTIME = {c: f"Heliocentric ecliptic J2000 {u} at the observation time (two-body propagation of the MPC orbit)."
                    for c, u in STATE_UNITS.items()}
STATE_AT_EPOCH = {c: f"Heliocentric ecliptic J2000 {u} at the epoch of the MPC elements (epoch_et), from the MPCORB elements."
                  for c, u in STATE_UNITS.items()}

PROPAGATED = {
    'obs_et': "Time of the observation as ephemeris time [TDB seconds past J2000].",
    'dt_from_epoch_days': "Time from the epoch of the MPC elements to the observation [days] (negative = before).",
    'a': "Semimajor axis [AU]: osculating two-body heliocentric element at the observation time (null if e >= 1).",
    'e': "Eccentricity: osculating two-body heliocentric element at the observation time.",
    'q': "Perihelion distance [AU]: osculating two-body heliocentric element at the observation time.",
    'incl': "Inclination [deg], ecliptic J2000: osculating two-body element at the observation time.",
    'node': "Longitude of the ascending node [deg], ecliptic J2000: osculating two-body element at the observation time.",
    'peri': "Argument of perihelion [deg], ecliptic J2000: osculating two-body element at the observation time.",
    'M': "Mean anomaly [deg]: osculating two-body element at the observation time.",
    'r': "Distance from the Sun [AU] at the observation time, light-time corrected (position when the light left the object).",
    'delta': "Distance from the Earth (geocenter) [AU] at the observation time, light-time corrected.",
    'phase_angle': "Phase angle Sun-object-Earth [deg] at the observation time, light-time corrected, geocentric.",
    'H': "Absolute magnitude [mag] from mag with the HG model (G = 0.15, V = mag, no band correction), r, delta and phase_angle.",
    'mag': "Apparent magnitude reported by Rubin, in the band of the detection (not converted to V).",
    'rmsmag': "Uncertainty of mag [mag], as text.",
}

# --------------------------------------------------------------------------------------- colors (one row per object)
COLORS = {
    'provid': "One of the Rubin provids of the object (the first one found); if n_provids = 2 the object has another "
              "one. To identify the object use mpc_designation.",
    'permid': "Permanent number of the object (numbered minor planets). Null for unnumbered objects.",
    'n_provids': IDENTIFIERS['n_provids'],
    'obstime_first': "Time of the first detection used for the object (UTC).",
    'obstime_last': "Time of the last detection used for the object (UTC).",
    **{f'H_{b}': f"Absolute magnitude in band {b} [mag]: mean of H over the object's {b} detections, weighted with 1/rmsmag^2."
       for b in BANDS},
    **{f'Herr_{b}': f"Uncertainty of H_{b} [mag]: 1/sqrt(sum of the weights 1/rmsmag^2)." for b in BANDS},
    **{f'Hstd_{b}': f"Weighted standard deviation of H in band {b} [mag] (light-curve scatter); null with one detection."
       for b in BANDS},
    **{f'N_{b}': f"Number of detections in band {b}." for b in BANDS},
    **{f'{b1}-{b2}': f"Color {b1}-{b2} [mag] = H_{b1} - H_{b2}." for b1, b2 in zip(BANDS[:-1], BANDS[1:])},
}

DESCRIPTIONS = {**SBN, **MPCORB, **PROPAGATED, **IDENTIFIERS}

# Per file: description of the file, and the columns that mean something different in that file
FILES = {
    'obs_sbn_X05_full_mpc_desig.parquet': dict(
        description="Rubin (X05) observations submitted to the Minor Planet Center (MPC obs_sbn table, one row per "
                    "detection), with the MPC designation of the object of each detection (mpc_designation, match_by, "
                    "n_provids). Written by add_mpc_designation.ipynb.",
        columns={}),
    'obs_elements_state.parquet': dict(
        description="Rubin detections with an MPC orbit (one row per detection): MPCORB elements at the MPC epoch "
                    "(mpc_*) and the corresponding state vector (x..vz, epoch_et). Written by exploration_data.ipynb.",
        columns={**STATE_AT_EPOCH, 'mag': "Apparent magnitude reported by Rubin, in the band of the detection."}),
    'obs_elements_state_at_obstime.parquet': dict(
        description="Two-body propagation of the MPC orbit to the time of each Rubin detection (one row per detection, "
                    "same order as obs_elements_state.parquet): osculating elements and state vector. "
                    "Written by propagate_orbit_1_elements.ipynb.",
        columns=STATE_AT_OBSTIME),
    'obs_elements_geometry_at_obstime.parquet': dict(
        description="Osculating elements, distance from the Earth and phase angle at the time of each Rubin detection "
                    "(one row per detection, same order as obs_elements_state.parquet). "
                    "Written by propagate_orbit_2_geometry_H.ipynb.",
        columns={}),
    'obs_final_elements_geometry_H.parquet': dict(
        description="Final table, one row per Rubin detection with an MPC orbit: osculating elements, distances from "
                    "the Sun and the Earth, phase angle, Rubin magnitude and absolute magnitude H. "
                    "Written by propagate_orbit_2_geometry_H.ipynb.",
        columns={}),
    'colors.parquet': dict(
        description="One row per object (mpc_designation): weighted mean H in each LSST band, its uncertainty, "
                    "scatter and number of detections, and the colors between adjacent bands. Written by colors.ipynb.",
        columns=COLORS),
}


def describe(table, file_name):
    """Arrow table with the description of each column and of the file in its schema (data unchanged)."""
    file = FILES[os.path.basename(file_name)]
    descriptions = {**DESCRIPTIONS, **file['columns']}
    missing = [f.name for f in table.schema if f.name not in descriptions]
    if missing:
        print(f'Warning: no description for {missing} in {file_name}: add them to parquet_docs.py')
    fields = [f.with_metadata({**(f.metadata or {}), b'description': descriptions[f.name].encode()})
              if f.name in descriptions else f
              for f in table.schema]
    metadata = {**(table.schema.metadata or {}), b'description': file['description'].encode()}
    return table.cast(pa.schema(fields, metadata=metadata))


def write_parquet(df, path):
    """Write a DataFrame to parquet with the column descriptions (like df.to_parquet(path, index=False))."""
    pq.write_table(describe(pa.Table.from_pandas(df, preserve_index=False), path), path)


def column_descriptions(path):
    """Description of the file and of each of its columns (name, type, description), without reading the data."""
    schema = pq.read_schema(path)
    print((schema.metadata or {}).get(b'description', b'').decode())
    return pd.DataFrame({'column': schema.names,
                         'type': [str(t) for t in schema.types],
                         'description': [(f.metadata or {}).get(b'description', b'').decode() for f in schema]})
