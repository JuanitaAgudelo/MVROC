Un pipeline para el análisis fotométrico y dinámico de SB en la era de VRO
===

1. **Obtener datos**. Descargar los datos reportados al MPC del Asteroid Institute: https://b612.ai/rubin-mpc-downloads/, en particular https://storage.googleapis.com/asteroid-institute-public/production/rubin/mpc/obs_sbn/full/parquet/obs_sbn_X05_full.parquet. Con las siguientes columnas:

      ['id', 'trksub', 'trkid', 'obsid', 'submission_id',
            'submission_block_id', 'obs80', 'status', 'ref', 'healpix', 'permid',
            'provid', 'artsat', 'mode', 'stn', 'trx', 'rcv', 'sys', 'ctr', 'pos1',
            'pos2', 'pos3', 'poscov11', 'poscov12', 'poscov13', 'poscov22',
            'poscov23', 'poscov33', 'prog', 'obstime', 'ra', 'dec', 'rastar',
            'decstar', 'obscenter', 'deltara', 'deltadec', 'dist', 'pa', 'rmsra',
            'rmsdec', 'rmsdist', 'rmspa', 'rmscorr', 'delay', 'rmsdelay', 'doppler',
            'rmsdoppler', 'astcat', 'mag', 'rmsmag', 'band', 'photcat', 'photap',
            'nucmag', 'logsnr', 'seeing', 'exp', 'rmsfit', 'com', 'frq', 'disc',
            'subfrm', 'subfmt', 'prectime', 'precra', 'precdec', 'unctime', 'notes',
            'remarks', 'deprecated', 'localuse', 'nstars', 'prev_desig', 'prev_ref',
            'rmstime', 'created_at', 'updated_at', 'trkmpc', 'orbit_id',
            'designation_asterisk', 'shapeocc', 'obssubid', 'replacesobsid',
            'group_id', 'vel1', 'fltr', 'vel3', 'vel2', 'obstime_text']

   1.1 De aquí nos sirven: provid y permid. 

2. **Cross match con MPC**. Se descargar la información orbital de todos los objetos con provisional ID en el MPC: https://www.minorplanetcenter.net/data

   Cada objeto en la base del MPC tiene un atributo que se llama: designation, other_desigs.

   Las designations pueden ser:
   | Format | Example | Objects |
   |---|---|---|
   | ⁠ (number) provisional designation ⁠: numbered, unnamed | ⁠ (4596) 1981 QB ⁠, ⁠ (800001) 2013 YW88 ⁠ | 865,766 |
   | ⁠ provisional designation ⁠: unnumbered | ⁠ 2021 AQ ⁠, ⁠ 1960 SB1 ⁠ | 680,535 |
   | ⁠ (number) Name ⁠: numbered and named | ⁠ (1) Ceres ⁠, ⁠ (22754) Olympus ⁠ | 26,521 |
   | ⁠ (number) NNNN P-L ⁠ / ⁠ T-1 ⁠, ⁠ T-2 ⁠, ⁠ T-3 ⁠: numbered, survey designation | ⁠ (12957) 2258 T-2 ⁠ | 3,623 |
   | ⁠ NNNN P-L ⁠ / ⁠ T-n ⁠: unnumbered, survey designation | ⁠ 4847 P-L ⁠ | 3 |


   2.1. Crossmatch:
        - Busca el provid|permid en designation. Si los encuentra se trae el designation como está en la base de datos del MPC al dataframe de MVROC.

        - Los que no encuentra los bisca en other_desigs. Si los encuentra hace lo mismo asignando el designation principal.

        - Si no los encuentra en ninguna de las dos se hace una búsqueda por API con provid|permiod. Devuelve JSON con toda la información te traes el equivalente a unpacked_primary_provisional_designation. Busco todos los unpacked_primary_provisional_designation en MPC si lo encuentra entonces se va para MVROC el designation de MPC. 

          * Parar la búsqueda en el primer hallazgo
          * Buscarlos todos

      Cada objeto del MPC tiene un solo provid. En el pipeline se llamara provid_mpc.

   2.2 Puede pasar que varios provid de VRO correspondan al mismo objeto en el MPC. En principio por cada provid_vro debería haber un provid_mpc. Pero no, un mismo provid_mpc puede corresponder a varios provid_vro. Construímos una tabla que por cada provid_vro identifca el provid_mpc que le corresponde.

   2.3 Generar un reporte: identificaciones directas, identificaciones secundarias, tabla de provid_mpc y sus provid_vro asociados. 

   Archivos: db_vro_mpc y info_vro_mpc.md

   Objects not indentified may be the MPC are processing and not solution.

   **Cuando llega uno nuevo**: hay que hacer otra vez todo el cross-match. Los archivos db_vro_mpc pueden cambiar.

3. **Obtención de órbitas**. Tomamos todos los provid_mpc de la tabla en 2.2 y extraemos de MPC data los elementos orbitales y la época reportada. Se hace una nueva tabla.

   Archivos: 'db_mvroc_orbits'.

   **Cuando llega un nuevo**: Se debe repetir este procedimiento (cambios de mpc designation y las órbitas pueden cambiar). 

4. **Vector de estado absoluto y relativo de cada observación**. Por cada una de las observaciones, cada 'id', leemos el provid_vro, consultamos la tabla de equivalencias y obtenemos el provid_mpc. Tomamos los elementos orbitales y los propagamos con SPICE.prop2b a la época de observación. Con el vector de estado resultante calculamos: r (distancia heliocéntrica), Delta (distancia geocéntrica) y alpha (ángulo de fase). Se genera una tabla con todos estos data. La table usa como primary field el id. 

   La posición de la Tierra de la Tierra en cada fecha de observación se obtiene a partir de un kernel de SPICE. Atención hay que descargar un kernel preciso para calcular posición de la Tierra respecto del Sol.

   Manejo del tiempo: la fecha en ISO UTC que es provista por VRO (obstime) genera: JD TDB. Se calcula el número de días desde la época de referencia hasta el día de observación en TDB (deltat).

   Archivos: 'db_mvroc_state'

   **Cuando llega uno nuevo**: Se hace un chequeo de objetos que no les haya cambiado los elementos y época (usando un hash) y a esos no se les recalcula nada. Crear una columna de 'flag_updated' para lo que sigue en el pipeline.

   Check for flags in the MPC data pointing out errors in orbit: uncertainty parameter. Study! condition_code (integer telling the quality of the orbits, u; number of opositions). 

5. **Cálculo de H para cada observación**. Por cada observación tomando la magnitud aparente, r, Delta, alpha y asumiendo el mismo valor de G para todos (G=0.15) calculamos H. La columna se guarda en el archivo de estado.

   Archivos: 'db_mvroc_state'

   **Cuando llega uno nuevo**: Se hace un chequeo obstime-mag (hash_obstime_mag). Si flag_updated=0 y el hash_obstime_mag no cambio no se hace nada. Si cualquiera de los dos cambio se recalcula.

6. **Filtro por número de observaciones**. Se filtra por número de observaciones (N_obs_min = 100) a este conjunto lo llamamos MOI, MVROC Object of Interest.

   Archivos: 'db_mvroc_moi' (tiene las columnas del original)

7. **Análisis fotométrico**. 

   7.1. Se hace periodograma de los 3 filtros en los que hayan más observaciones. Producen período y potencia.

   7.2. Se hace periodograma combinados, ej. Lr+Lg, Lr+Li, Lg+Li, Lr+Lg+Li. Producen período y potencia.

   7.3. Votación: se mira cuáles períodos son aliases diarios. Se deja el conjunto de periodos únicos y se vota para ver cuál es el más popular.

   7.4. Se hace un ajuste de todas las bandas con una función armónica con cada período, el popular y los demás y se evaluúa el chi-2 reducido de cada ajuste. El período con el menor chi-2 pasa a ser el preferido. 

   7.5. Se identifican los periodos que son aliases de la time-window (span, campaign duration). 

   Archivos: completamos 'db_mcrov_moi' con: periodo mejor, chi-2, potencia, periodos secundarios, chi-2, potencia, tags de aliases; coeficientes c_x que ajustan las diferencias entre filtos; índices de color (los disponibles).

8. **Selección de las curvas de luz que revelan rotación**. Aquí es donde puede usarse la OI y la IA. 

   8.1. OI: Mirar las curvas folded de todos los períodos por objeto. Se marcan los objetos cuyas curvas muestran rotación y pasan al MOC.

   Archivos: 'db_mvroc_moc': MVROC Objects with Candidate light-curves

   8.2. IA: Entrenar un algoritmo para que de alguna manera saque las mejores curvas.

9. **Generación de los períodos detectados**. Después del MOC se eligen las soluciones, también por inspección visual. Para cada objeto se revisan los períodos candidatos y se conservan los que pliegan bien los datos. Cada período elegido lleva una multiclasificación: una o más etiquetas que dicen cómo pliega.

   - **clean-fold** (c): poca dispersión; los puntos siguen bien el modelo de dos picos.
   - **scattered** (s): errores grandes o mucha dispersión; la forma de la curva se reconoce.
   - **contaminated** (x): unos pocos puntos quedan lejos del pliegue principal.
   - **phase-incomplete** (p): solo se cubre un tramo de la rotación.
   - **low-amplitude** (a): amplitud pico a pico pequeña; la rotación se ve igual (a menudo en rotadores rápidos).
   - **misaligned** (m): problemas de color o de punto cero; la curva de rotación sigue clara.
   - **undersampled** (u): pocos puntos y/o huecos grandes en fase.
   - **excluded** (e): se reclasifica después; sale del conjunto científico.

   Un mismo período puede llevar varias etiquetas (por ejemplo `c` y `p`).

   Archivos: 'db_mvroc_sol': MVROC Solutions. Una fila por período retenido.

10. **Clasificación de curvas de luz**. Con las etiquetas del paso 9 se asigna una calidad Q de 1 a 9: 1–2 es alta confianza y 7–9 pide revisión profunda. El esquema es sencillo. Si es `c` puro, Q = 1. Si es `c-p`, Q = 2. Si es `c-x`, Q = 1, etc. Las bandas de uso son A (Q = 1–2), B (3–4), C (5–6) y D (7–9).

   NOTA: Podemos entrenar algoritmos para clasificación múltiple.

11. Colores a taxonomía.

12. Dinámica.

13. Shape. Solo los de clase A y B.

-----
[MVROC Team](https://github.com/seap-udea/MVROC) (2026-present)
