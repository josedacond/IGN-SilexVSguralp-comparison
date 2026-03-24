# =============================================================================
#         INFORME AUTOMÁTICO: TODAS LAS TANDAS → PDF CON REPORTLAB
# =============================================================================

import sys
import os
import numpy as np
from scipy.signal import hilbert
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from obspy import read, UTCDateTime
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import cm
from reportlab.pdfgen import canvas
import subprocess

# =============================================================================
#                         PARÁMETROS
# =============================================================================

RUTA_CONFIG = "/Users/joseda_cond/Desktop/IGN/rvbecarios/Ensayo_20260303/config.txt"

VISTAS = [
    ("Señal Completa", 1, 32.0),
    ("Primeros 1.5 segundos", 0.2, 1.5),
    ("Primeros 5 segundos", 0.5, 5.0),
]

COMPONENTES = ["Z", "N", "E"]


# =============================================================================
#                           FICHERO DE CONFIGURACIÓN
# =============================================================================

def leer_config(ruta_config):
    with open(ruta_config, 'r', encoding='utf-8') as f:
        lineas = f.readlines()

    ruta_base = ""
    duracion = 60
    tandas = []
    tanda_actual = None

    for linea in lineas:
        linea = linea.strip()

        if linea.startswith("RUTA_BASE:"):
            ruta_base = linea.split(":", 1)[1].strip()
            continue
        if linea.startswith("DURACION:"):
            duracion = int(linea.split(":", 1)[1].strip())
            continue
        if not linea:
            continue

        partes = [p.strip() for p in linea.split(";")]

        if len(partes) == 2:
            if tanda_actual is not None:
                tandas.append(tanda_actual)
            tanda_actual = {
                'cmg_ruta': os.path.join(ruta_base, partes[0]),
                'cmg_hora': partes[1],
                'cmg_nombre': os.path.basename(partes[0]),
                'silex': []
            }
        elif len(partes) == 3:
            tanda_actual['silex'].append({
                'ruta': os.path.join(ruta_base, partes[0]),
                'hora': partes[1],
                'orientacion': partes[2].upper(),
                'nombre': os.path.basename(partes[0])
            })

    if tanda_actual is not None:
        tandas.append(tanda_actual)

    return ruta_base, duracion, tandas


# =============================================================================
#                    FUNCIÓN DE CARGA Y PREPROCESADO
# =============================================================================

def cargar_y_preprocesar(ruta, es_cmg=False, orientacion=None):
    print(f"    Cargando: {os.path.basename(ruta)}")
    st = read(ruta)

    if not es_cmg:
        st.merge(method=1, fill_value="interpolate")

    if es_cmg:
        for tr in st:
            if tr.stats.channel.endswith("Z") or tr.stats.channel.endswith("E"):
                tr.data = -1 * tr.data

    if orientacion == "S":
        for tr in st:
            if tr.stats.channel.endswith("Z") or tr.stats.channel.endswith("N"):
                tr.data = -1 * tr.data
    elif orientacion == "N":
        for tr in st:
            if tr.stats.channel.endswith("Z") or tr.stats.channel.endswith("E"):
                tr.data = -1 * tr.data

    for tr in st:
        tr.detrend("demean")

    if es_cmg:
        factor = 2e-4
        for tr in st:
            tr.data = tr.data * factor

    st.taper(max_percentage=0.1, type='hamming', max_length=5.0)
    st.filter("bandpass", freqmin=1, freqmax=25.0, corners=4, zerophase=True)

    return st


# =============================================================================
#     FUNCIÓN: ALINEACIÓN CON CORRELACIÓN CRUZADA DE ENVOLVENTES (3 COMP)
#              + AJUSTE FINO CON FORMA DE ONDA
# =============================================================================
# NOTA: El ajuste de inicio (umbral) se hace FUERA, una sola vez por tanda.

def alinear_por_envolvente(st_cmg, st_silex, hora_cmg_str, hora_silex_str):
    fs = st_cmg[0].stats.sampling_rate

    fecha_cmg = str(st_cmg[0].stats.starttime.date)
    fecha_silex = str(st_silex[0].stats.starttime.date)
    t_ref_cmg = UTCDateTime(f"{fecha_cmg}T{hora_cmg_str}")
    t_ref_silex = UTCDateTime(f"{fecha_silex}T{hora_silex_str}")

    # --- PASO 1: Energía combinada del CMG (15s desde la hora ref) ---
    ventana_ref = 15.0
    trazas_cmg = {}
    for comp in COMPONENTES:
        tr = st_cmg.select(component=comp)[0].copy()
        tr.trim(t_ref_cmg, t_ref_cmg + ventana_ref)
        trazas_cmg[comp] = tr.data.astype(float)

    n_min = min(len(trazas_cmg[c]) for c in COMPONENTES)
    energia_cmg = np.sqrt(sum(trazas_cmg[c][:n_min]**2 for c in COMPONENTES))

    # --- Energía combinada del Silex (±30s alrededor de la hora ref) ---
    margen = 30.0
    trazas_silex = {}
    for comp in COMPONENTES:
        tr = st_silex.select(component=comp)[0].copy()
        tr.trim(t_ref_silex - margen, t_ref_silex + margen)
        trazas_silex[comp] = tr.data.astype(float)

    n_min_s = min(len(trazas_silex[c]) for c in COMPONENTES)
    energia_silex = np.sqrt(sum(trazas_silex[c][:n_min_s]**2 for c in COMPONENTES))

    # --- Envolventes suavizadas ---
    env_cmg = np.abs(hilbert(energia_cmg))
    env_silex = np.abs(hilbert(energia_silex))

    ventana_suav = int(0.5 * fs)
    kernel = np.ones(ventana_suav) / ventana_suav
    env_cmg_suave = np.convolve(env_cmg, kernel, mode='same')
    env_silex_suave = np.convolve(env_silex, kernel, mode='same')

    # --- Correlación cruzada normalizada ---
    env_cmg_norm = (env_cmg_suave - np.mean(env_cmg_suave)) / (np.std(env_cmg_suave) * len(env_cmg_suave))
    env_silex_norm = (env_silex_suave - np.mean(env_silex_suave)) / np.std(env_silex_suave)

    cc = np.correlate(env_silex_norm, env_cmg_norm, mode='full')

    pico_idx = np.argmax(cc)
    coeficiente = cc[pico_idx]
    posicion_en_ventana = pico_idx - len(env_cmg_norm) + 1

    t_inicio_ventana = t_ref_silex - margen
    t_evento_silex = t_inicio_ventana + posicion_en_ventana / fs
    t_evento_cmg = t_ref_cmg

    desfase = t_evento_silex - t_evento_cmg
    print(f"    Coef. envolventes: {coeficiente:.4f} | Desfase grueso: {desfase*1000:.1f} ms")

    if coeficiente < 0.3:
        print(f"    AVISO: Correlación muy baja ({coeficiente:.4f})")

    # --- PASO 2: Ajuste fino con correlación de forma de onda ---
    ventana_fina = 10.0
    margen_fino = 2.0

    tr_ref = st_cmg.select(component="Z")[0].copy()
    tr_ref.trim(t_evento_cmg, t_evento_cmg + ventana_fina)

    tr_sx = st_silex.select(component="Z")[0].copy()
    tr_sx.trim(t_evento_silex - margen_fino, t_evento_silex + ventana_fina + margen_fino)

    ref = tr_ref.data.astype(float)
    vent = tr_sx.data.astype(float)

    ref_norm = (ref - np.mean(ref)) / (np.std(ref) * len(ref))
    vent_norm = (vent - np.mean(vent)) / np.std(vent)

    cc_fina = np.correlate(vent_norm, ref_norm, mode='full')
    pico_fino = np.argmax(cc_fina)
    coef_fino = cc_fina[pico_fino]
    pos_fina = pico_fino - len(ref_norm) + 1

    t_evento_silex_fino = tr_sx.stats.starttime + pos_fina / fs
    correccion = t_evento_silex_fino - t_evento_silex

    print(f"    Coef. fino: {coef_fino:.4f} | Corrección: {correccion*1000:.1f} ms")

    t_evento_silex = t_evento_silex_fino
    desfase_final = t_evento_silex - t_evento_cmg
    print(f"    Desfase final: {desfase_final*1000:.1f} ms")

    return t_evento_cmg, t_evento_silex, desfase_final


# =============================================================================
#     FUNCIÓN: AJUSTE DE INICIO POR ENERGÍA SOSTENIDA (UNA VEZ POR TANDA)
# =============================================================================

def calcular_avance_inicio(st_cmg, t_ref_cmg):
    fs = st_cmg[0].stats.sampling_rate
    tr_cmg_z = st_cmg.select(component="Z")[0].copy()
    tr_cmg_z.trim(t_ref_cmg, t_ref_cmg + 30)

    env = np.abs(hilbert(tr_cmg_z.data.astype(float)))
    tam_ventana = int(1.0 * fs)
    n_ventanas = len(env) // tam_ventana

    energias = []
    for v in range(n_ventanas):
        trozo = env[v * tam_ventana : (v + 1) * tam_ventana]
        energias.append(np.mean(trozo))

    energias = np.array(energias)
    max_energia = np.max(energias)

    umbral = max_energia * 0.10
    idx_ventana = 0
    for v in range(n_ventanas - 2):
        if energias[v] > umbral and energias[v+1] > umbral and energias[v+2] > umbral:
            idx_ventana = v
            break

    idx = max(0, idx_ventana * tam_ventana - int(0.5 * fs))
    avance = idx / fs
    return avance


# =============================================================================
#           FUNCIÓN: GENERAR PLOT 6x1 Y GUARDAR COMO PNG
# =============================================================================

def generar_plot_6x1(st_cmg, st_silex, t_evento_c, t_evento_s,
                     dur_antes, dur_despues, titulo, nombre_cmg, nombre_silex,
                     ruta_png, mostrar_stats=False):
    st_c = st_cmg.copy()
    st_s = st_silex.copy()
    st_c.trim(t_evento_c - dur_antes, t_evento_c + dur_despues)
    st_s.trim(t_evento_s - dur_antes, t_evento_s + dur_despues)

    fig, axes = plt.subplots(6, 1, figsize=(8, 11.5), sharex=True)
    fig.suptitle(f"{titulo}\n{nombre_cmg} vs {nombre_silex}", fontsize=12, fontweight='bold')

    for j, comp in enumerate(COMPONENTES):
        tr_c = st_c.select(component=comp)[0]
        tr_s = st_s.select(component=comp)[0]

        t_rel_c = tr_c.times() - dur_antes
        t_rel_s = tr_s.times() - dur_antes

        max_val_cmg = max(abs(tr_c.data.min()), abs(tr_c.data.max()))
        max_val_silex = max(abs(tr_s.data.min()), abs(tr_s.data.max()))
        max_val = max(max_val_cmg, max_val_silex)
        limite_y = int(max_val * 1.05) + 1

        idx_cmg = j * 2
        idx_silex = j * 2 + 1

        axes[idx_cmg].plot(t_rel_c, tr_c.data, color='blue', lw=0.6)
        axes[idx_cmg].set_title(f"CMG - Componente {comp}", fontsize=10, loc='left')
        axes[idx_cmg].set_ylabel("mg", fontsize=8)
        axes[idx_cmg].set_ylim(-limite_y, limite_y)
        axes[idx_cmg].grid(True, alpha=0.3)

        axes[idx_silex].plot(t_rel_s, tr_s.data, color='black', lw=0.6)
        axes[idx_silex].set_title(f"Silex ({nombre_silex}) - Componente {comp}", fontsize=10, loc='left')
        axes[idx_silex].set_ylabel("mg", fontsize=8)
        axes[idx_silex].set_ylim(-limite_y, limite_y)
        axes[idx_silex].grid(True, alpha=0.3)

        if mostrar_stats:
            max_c = tr_c.data.max()
            min_c = tr_c.data.min()
            max_s = tr_s.data.max()
            min_s = tr_s.data.min()

            porc_max = (max_s / max_c * 100) if max_c != 0 else 0
            porc_min = (min_s / min_c * 100) if min_c != 0 else 0

            # Texto de estadísticas
            texto_c = f"Max: {max_c:.2f} mg | Min: {min_c:.2f} mg"
            axes[idx_cmg].text(0.98, 0.92, texto_c, transform=axes[idx_cmg].transAxes,
                               fontsize=7, ha='right', va='top',
                               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

            texto_s = f"Max: {max_s:.2f} mg ({porc_max:.1f}%) | Min: {min_s:.2f} mg ({porc_min:.1f}%)"
            axes[idx_silex].text(0.98, 0.92, texto_s, transform=axes[idx_silex].transAxes,
                                 fontsize=7, ha='right', va='top',
                                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

            # Líneas rojas del máximo y mínimo del CMG (en ambos plots)
            for ax_idx in [idx_cmg, idx_silex]:
                axes[ax_idx].axhline(y=max_c, color='#8B0000', linestyle='--', lw=0.8, alpha=0.7)
                axes[ax_idx].axhline(y=min_c, color='#8B0000', linestyle='--', lw=0.8, alpha=0.7)

            # Etiquetas en las líneas rojas del CMG
            axes[idx_cmg].text(0.01, max_c, f" {max_c:.1f}", transform=axes[idx_cmg].get_yaxis_transform(),
                               fontsize=6, color='#8B0000', va='bottom')
            axes[idx_cmg].text(0.01, min_c, f" {min_c:.1f}", transform=axes[idx_cmg].get_yaxis_transform(),
                               fontsize=6, color='#8B0000', va='top')

            # Etiquetas en las líneas rojas del Silex (valores del CMG como referencia)
            axes[idx_silex].text(0.01, max_c, f" {max_c:.1f}", transform=axes[idx_silex].get_yaxis_transform(),
                                 fontsize=6, color='#8B0000', va='bottom')
            axes[idx_silex].text(0.01, min_c, f" {min_c:.1f}", transform=axes[idx_silex].get_yaxis_transform(),
                                 fontsize=6, color='#8B0000', va='top')

    axes[-1].set_xlabel("Tiempo desde el inicio (s)")
    plt.tight_layout()
    plt.savefig(ruta_png, dpi=150, bbox_inches='tight')
    plt.close(fig)


# =============================================================================
#                    FUNCIÓN: GENERAR PDF CON REPORTLAB
# =============================================================================

def generar_pdf(lista_pngs, ruta_pdf):
    ancho, alto = A4
    c = canvas.Canvas(ruta_pdf, pagesize=A4)

    for png_path in lista_pngs:
        img_ancho = ancho - 2 * cm
        img_alto = alto - 2 * cm
        c.drawImage(png_path, 1 * cm, 1 * cm, width=img_ancho, height=img_alto,
                    preserveAspectRatio=True, anchor='c')
        c.showPage()

    c.save()
    print(f"\nPDF generado: {ruta_pdf}")


# =============================================================================
#                         EJECUCIÓN PRINCIPAL
# =============================================================================

print("\n" + "=" * 60)
print("     INFORME AUTOMÁTICO: TODAS LAS TANDAS → PDF")
print("=" * 60)

# 1. Leer configuración
print("\nLeyendo configuración...")
ruta_base, duracion, tandas = leer_config(RUTA_CONFIG)
print(f"  Ruta base: {ruta_base}")
print(f"  Duración: {duracion}s")
print(f"  Tandas: {len(tandas)}")

total_silex = sum(len(t['silex']) for t in tandas)
total_paginas = total_silex * len(VISTAS)
print(f"  Silex totales: {total_silex}")
print(f"  Páginas a generar: {total_paginas}")

# 2. Carpeta temporal para PNGs
carpeta_informe = os.path.join(ruta_base, "INFORME")
if not os.path.exists(carpeta_informe):
    os.makedirs(carpeta_informe)

tmp_dir = os.path.join(carpeta_informe, "PNGs")
# Limpiamos la carpeta si ya existe, si no la creamos
if os.path.exists(tmp_dir):
    for archivo in os.listdir(tmp_dir):
        if archivo.endswith('.png'):
            os.remove(os.path.join(tmp_dir, archivo))
else:
    os.makedirs(tmp_dir)
lista_pngs = []
contador_par = 0

# 3. Bucle principal: recorrer todas las tandas
for idx_tanda, tanda in enumerate(tandas):
    num_tanda = idx_tanda + 1
    num_silex_tanda = len(tanda['silex'])

    print(f"\n{'='*60}")
    print(f"  TANDA {num_tanda}/{len(tandas)}: {tanda['cmg_nombre']} ({num_silex_tanda} Silex)")
    print(f"{'='*60}")

    # Cargar y preprocesar CMG (una vez por tanda)
    print(f"\n  Cargando CMG...")
    st_cmg = cargar_y_preprocesar(tanda['cmg_ruta'], es_cmg=True)
    fs_cmg = st_cmg[0].stats.sampling_rate

    # --- FASE 1: Alinear TODOS los Silex de esta tanda ---
    print(f"\n  FASE 1: Alineando todos los Silex de la tanda...")
    resultados_tanda = []  # Lista de {silex_info, st_silex, st_cmg_par, t_cmg, t_silex, desfase}

    for idx_silex, silex_info in enumerate(tanda['silex']):
        contador_par += 1
        nombre_silex_limpio = os.path.splitext(silex_info['nombre'])[0]
        nombre_cmg_limpio = os.path.splitext(tanda['cmg_nombre'])[0]

        print(f"\n  --- Par {contador_par}/{total_silex}: {nombre_cmg_limpio} vs {nombre_silex_limpio} ---")

        try:
            st_silex = cargar_y_preprocesar(
                silex_info['ruta'], es_cmg=False, orientacion=silex_info['orientacion']
            )
        except Exception as e:
            print(f"    ERROR cargando {silex_info['nombre']}: {e}")
            continue

        # Igualar frecuencias
        fs_silex = st_silex[0].stats.sampling_rate
        st_cmg_par = st_cmg.copy()
        if fs_cmg != fs_silex:
            fs_comun = min(fs_cmg, fs_silex)
            if fs_cmg > fs_comun:
                st_cmg_par.resample(fs_comun)
            if fs_silex > fs_comun:
                st_silex.resample(fs_comun)

        # Alinear
        print(f"    Alineando...")
        try:
            t_evento_c, t_evento_s, desfase = alinear_por_envolvente(
                st_cmg_par, st_silex, tanda['cmg_hora'], silex_info['hora']
            )
        except Exception as e:
            print(f"    ERROR alineando: {e}")
            continue

        resultados_tanda.append({
            'silex_info': silex_info,
            'st_silex': st_silex,
            'st_cmg_par': st_cmg_par,
            't_cmg': t_evento_c,
            't_silex': t_evento_s,
            'desfase': desfase,
            'nombre_silex': nombre_silex_limpio,
            'nombre_cmg': nombre_cmg_limpio,
        })

    if len(resultados_tanda) == 0:
        print(f"  No se alinearon Silex en esta tanda. Saltando.")
        continue

    # --- FASE 2: Corregir outliers de desfase ---
    print(f"\n  FASE 2: Verificando desfases de la tanda...")
    desfases = [r['desfase'] for r in resultados_tanda]
    mediana = np.median(desfases)
    print(f"    Desfase mediano de la tanda: {mediana*1000:.1f} ms")

    for r in resultados_tanda:
        diferencia = abs(r['desfase'] - mediana)
        if diferencia > 0.5:  # Más de 500ms de la mediana = outlier
            correccion = mediana - r['desfase']
            print(f"    OUTLIER: {r['nombre_silex']} (desfase: {r['desfase']*1000:.1f} ms, "
                  f"mediana: {mediana*1000:.1f} ms) -> Corrigiendo {correccion*1000:.1f} ms")
            r['t_silex'] = r['t_silex'] + correccion
            r['desfase'] = mediana

    # --- FASE 3: Ajuste de inicio (una vez por tanda, usando el CMG) ---
    print(f"\n  FASE 3: Ajustando inicio del evento...")
    # Usamos el t_cmg del primer resultado (todos comparten el mismo CMG)
    t_ref_cmg = UTCDateTime(f"{str(st_cmg[0].stats.starttime.date)}T{tanda['cmg_hora']}")
    avance = calcular_avance_inicio(resultados_tanda[0]['st_cmg_par'], t_ref_cmg)
    print(f"    Inicio real del evento: +{avance:.2f}s desde la hora del .txt")

    # Aplicar el avance a todos los resultados
    for r in resultados_tanda:
        r['t_cmg'] = r['t_cmg'] + avance
        r['t_silex'] = r['t_silex'] + avance

    # --- FASE 4: Generar los PNGs ---
    print(f"\n  FASE 4: Generando gráficas...")
    for r in resultados_tanda:
        for i, (vista_nombre, dur_antes, dur_despues) in enumerate(VISTAS):
            titulo = f"Tanda {num_tanda} - {vista_nombre}"
            ruta_png = os.path.join(tmp_dir, f"tanda{num_tanda}_{r['nombre_silex']}_{i}.png")

            try:
                generar_plot_6x1(
                    r['st_cmg_par'], r['st_silex'], r['t_cmg'], r['t_silex'],
                    dur_antes, dur_despues, titulo,
                    r['nombre_cmg'], r['nombre_silex'],
                    ruta_png, mostrar_stats=(i == 0)
                )
                lista_pngs.append(ruta_png)
            except Exception as e:
                print(f"    ERROR generando gráfica '{vista_nombre}' para {r['nombre_silex']}: {e}")

        print(f"    ✓ {r['nombre_silex']} completado ({len(VISTAS)} páginas)")

# 4. Generar PDF final
print(f"\n{'='*60}")
print(f"  Generando PDF final ({len(lista_pngs)} páginas)...")
ruta_pdf = os.path.join(carpeta_informe, "Informe_Ensayo_20260303.pdf")
generar_pdf(lista_pngs, ruta_pdf)

# 5. Abrir el PDF automáticamente (multiplataforma)
if sys.platform == "darwin":
    subprocess.Popen(["open", ruta_pdf])
elif sys.platform == "win32":
    os.startfile(ruta_pdf)
else:
    subprocess.Popen(["xdg-open", ruta_pdf])


print(f"\n{'='*60}")
print(f"  INFORME COMPLETADO")
print(f"  {len(lista_pngs)} páginas generadas para {contador_par} pares CMG-Silex")
print(f"  PDF: {ruta_pdf}")
print(f"{'='*60}")







