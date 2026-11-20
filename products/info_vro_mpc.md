# Rubin–MPC cross-match report

Built from `code/data/mpc_designations.csv`, `colors.parquet`, and
`rubin_objects_two_provids.csv` (pipeline step 2).

## Summary

- Unique `provid_vro` with an MPC designation: **300,952**
- Unique `provid_mpc`: **297,966**
- Numbered MPC designations: **180,379**
- Provisional / unnumbered MPC designations: **120,573**
- `provid_vro` that share an MPC designation with another Rubin id: **5,972**
- MPC objects known to gather several Rubin provids (`rubin_objects_two_provids.csv`): **2,989**

## Many-to-one associations

A single `provid_mpc` can collect several `provid_vro`. The equivalence
table `db_vro_mpc.csv` records every pair; `n_vro_for_mpc` counts how
many Rubin ids map to that MPC designation.

Objects with ≥2 Rubin provids (from the curated list): 2989.

| mpc_designation | provids | n_obs_total |
|---|---|---:|
| (831800) 2009 WC150 | 2009 WC150, 2025 MN83 | 409 |
| (831814) 2009 WY177 | 2009 WY177, 2025 ML51 | 66 |
| (835088) 2011 BY148 | 2025 MS60, 2011 BY148 | 100 |
| (840358) 2015 AF100 | 2015 AF100, 2025 MH16 | 251 |
| (881568) 2015 KL461 | 2025 MR86, 2015 KL461 | 425 |
| (887157) 2001 FK241 | 2001 FK241, 2025 OH216 | 23 |
| (887234) 2003 SH392 | 2025 OT185, 2003 SH392 | 34 |
| (887333) 2004 RT278 | 2025 NK138, 2004 RT278 | 17 |
| (887532) 2005 WW172 | 2005 WW172, 2025 NW41 | 29 |
| (887538) 2005 XC56 | 2025 MA248, 2005 XC56 | 32 |
| (887595) 2006 CN22 | 2025 OJ274, 2006 CN22 | 21 |
| (887623) 2006 HR137 | 2006 HR137, 2025 OE210 | 22 |
| (887712) 2006 VU71 | 2025 MG103, 2006 VU71 | 27 |
| (887756) 2007 BB33 | 2025 NS198, 2007 BB33 | 17 |
| (888607) 2009 HM64 | 2009 HM64, 2025 MV241 | 68 |
| (888680) 2009 UX196 | 2009 UX196, 2025 MJ202 | 26 |
| (889204) 2011 QP81 | 2011 QP81, 2025 PV158 | 40 |
| (889682) 2012 UF22 | 2025 MD125, 2012 UF22 | 42 |
| (891660) 2015 BB238 | 2015 BB238, 2025 OP156 | 38 |
| (892516) 2015 PQ164 | 2015 PQ164, 2025 NL182 | 42 |
| (895882) 2025 NQ153 | 2025 NQ153, 2025 PC98 | 32 |
| 1995 SN14 | 1995 SN14, 2025 NZ84 | 12 |
| 1999 CL129 | 2025 OA65, 1999 CL129 | 12 |
| 2000 OG62 | 2000 OG62, 2025 PX119 | 34 |
| 2000 OO66 | 2000 OO66, 2025 MO64 | 273 |
| … | (2964 more rows) | |
