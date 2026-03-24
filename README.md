# IGN-SilexVSguralp-comparison

## Herramienta de comparación automatizada de sensores sísmicos Silex vs Guralp CMG-5T

Herramienta desarrollada en Python para el **Instituto Geográfico Nacional (IGN)** como parte del programa de prácticas externas de la Universidad Complutense de Madrid (UCM). Su objetivo es automatizar la comparación de señales sísmicas registradas por sensores **Silex** (en pruebas) frente al sensor patrón **Guralp CMG-5T**, ambos montados sobre una mesa vibrante en las instalaciones del IGN.

El programa procesa los datos sísmicos en formato MiniSEED/SEED, alinea temporalmente las señales de forma automática, y genera un **informe en PDF** con gráficas comparativas de las tres componentes espaciales (Z, N, E) en distintas ventanas temporales.

---

## Índice

- [Requisitos e instalación](#requisitos-e-instalación)
- [Estructura del proyecto](#estructura-del-proyecto)
- [Configuración: archivo config.txt](#configuración-archivo-configtxt)
- [Ejecución](#ejecución)
- [Procesado de señales](#procesado-de-señales)
- [Ejemplo de salida](#ejemplo-de-salida)

---

## Requisitos e instalación

### Python

El script requiere **Python 3.8 o superior**. Se puede verificar la versión instalada con:

```bash
python3 --version
```

### Dependencias

Instalar las siguientes bibliotecas mediante pip:

```bash
pip install obspy numpy scipy matplotlib reportlab
```

**Descripción de cada biblioteca:**

| Biblioteca | Uso en el proyecto |
|---|---|
| `obspy` | Lectura de archivos MiniSEED/SEED, procesado de señales sísmicas |
| `numpy` | Operaciones numéricas y cálculo de correlaciones |
| `scipy` | Transformada de Hilbert para el cálculo de envolventes |
| `matplotlib` | Generación de gráficas comparativas |
| `reportlab` | Exportación del informe final a PDF |

### Verificar instalación

Para comprobar que todo está instalado correctamente:

```bash
python3 -c "import obspy, numpy, scipy, matplotlib, reportlab; print('Todas las dependencias instaladas correctamente.')"
```

---

## Estructura del proyecto

```
IGN-SilexVSguralp-comparison/
├── script_informe.py       # Script principal (genera el informe PDF)
├── config.txt              # Archivo de configuración (rutas, horas, orientaciones)
├── README.md               # Este documento
└── PNGs/                   # Se crea automáticamente dentro de RUTA_BASE al ejecutar el script
```

---

## Configuración: archivo config.txt

El archivo `config.txt` es el **único archivo que el usuario necesita modificar**. Contiene toda la información necesaria para que el script procese los datos automáticamente.

### Formato general

```
RUTA_BASE: /ruta/absoluta/a/la/carpeta/del/ensayo
DURACION: 60

ruta/relativa/al/CMG.seed; HH:MM:SS
ruta/relativa/al/Silex1.mseed; HH:MM:SS; N
ruta/relativa/al/Silex2.mseed; HH:MM:SS; S
ruta/relativa/al/Silex3.mseed; HH:MM:SS; N

ruta/relativa/al/CMG2.seed; HH:MM:SS
ruta/relativa/al/Silex4.mseed; HH:MM:SS; S
...
```

### Parámetros globales

| Parámetro | Descripción |
|---|---|
| `RUTA_BASE` | Ruta absoluta a la carpeta raíz del ensayo. Todas las rutas de archivos son relativas a esta. |
| `DURACION` | Duración en segundos del evento sísmico (no se usa actualmente en el corte, pero se mantiene como referencia). |

### Estructura de tandas

Los datos se organizan en **tandas**. Cada tanda tiene:
- **Una línea de CMG** (sensor patrón): `ruta_relativa; hora_referencia`
- **Varias líneas de Silex** (sensores a comparar): `ruta_relativa; hora_referencia; orientación`

Las tandas se separan entre sí por **una o más líneas vacías**.

### Campos de cada línea

**Línea de CMG (2 campos separados por `;`):**

| Campo | Descripción | Ejemplo |
|---|---|---|
| Ruta relativa | Ruta al archivo .seed/.mseed desde la RUTA_BASE | `Cmg5T/Cmg5tTan1ZNE.seed` |
| Hora de referencia | Hora aproximada del inicio del evento (HH:MM:SS) | `09:20:50` |

**Línea de Silex (3 campos separados por `;`):**

| Campo | Descripción | Ejemplo |
|---|---|---|
| Ruta relativa | Ruta al archivo .mseed/.seed desde la RUTA_BASE | `Silex/carpeta_silex/ES.SX244.mseed` |
| Hora de referencia | Hora aproximada del inicio del evento (HH:MM:SS) | `09:20:50` |
| Orientación | Orientación del sensor en la mesa vibrante: `N` (Norte) o `S` (Sur) | `N` |

### Ejemplo real

```
RUTA_BASE: /home/usuario/Ensayo_20260303
DURACION: 60

Cmg5T/Cmg5tTan1ZNE.seed; 09:20:50
Silex/Silex LA244 - LA245 - LA246 - LA247 - LA248 - LA249 - LA250/ES.SX244.mseed; 09:20:50; N
Silex/Silex LA244 - LA245 - LA246 - LA247 - LA248 - LA249 - LA250/ES.SX245.mseed; 09:20:50; N
Silex/Silex LA244 - LA245 - LA246 - LA247 - LA248 - LA249 - LA250/ES.SX246.mseed; 09:20:50; N
Silex/Silex LA244 - LA245 - LA246 - LA247 - LA248 - LA249 - LA250/ES.SX247.mseed; 09:20:50; S
Silex/Silex LA244 - LA245 - LA246 - LA247 - LA248 - LA249 - LA250/ES.SX248.mseed; 09:20:50; S
Silex/Silex LA244 - LA245 - LA246 - LA247 - LA248 - LA249 - LA250/ES.SX249.mseed; 09:20:50; S
Silex/Silex LA244 - LA245 - LA246 - LA247 - LA248 - LA249 - LA250/ES.SX250.mseed; 09:20:50; S

Cmg5T/Cmg5tTan2ZNE.seed; 09:45:05
Silex/Silex LA251 - LA252 - LA253 - LA254 - LA255 - LA256 - LA257/ES.SX251.mseed; 09:45:05; N
...
```

### Notas importantes sobre el config.txt

- La **hora de referencia** no necesita ser exacta. Es una hora aproximada cercana al inicio del evento. El programa detecta automáticamente el inicio real.
- La **orientación** (`N` o `S`) determina qué corrección de polaridad se aplica al Silex. Es fundamental que sea correcta para que la comparación sea válida.
- Si un archivo no se encuentra en la ruta indicada, el programa lo salta y continúa con el siguiente, mostrando un mensaje de error en la consola.
- Los nombres de carpetas con espacios y caracteres especiales funcionan correctamente siempre que se escriban exactamente igual que en el sistema de archivos.

---

## Ejecución

### Antes de ejecutar

1. Asegurarse de que todas las dependencias están instaladas.
2. Verificar que el archivo `config.txt` está correctamente configurado con las rutas, horas y orientaciones.
3. Verificar que la variable `RUTA_CONFIG` al inicio del script apunta al archivo `config.txt`.

### Ejecutar el script

```bash
python3 script_informe.py
```

### ¿Qué hace el programa al ejecutarse?

1. Lee el archivo `config.txt`.
2. Procesa cada tanda secuencialmente:
   - Carga y preprocesa el CMG (una vez por tanda).
   - Carga y preprocesa cada Silex.
   - Alinea temporalmente cada par CMG–Silex.
   - Corrige automáticamente alineaciones atípicas (outliers).
   - Detecta el inicio real del evento sísmico.
   - Genera 3 gráficas por cada par (señal completa, 1.5 s, 5 s).
3. Compila todas las gráficas en un **único PDF**.
4. Guarda los PNGs individuales en la carpeta `PNGs/`.
5. Abre automáticamente el PDF generado.

### Salida

- **PDF**: Se genera en la carpeta definida en `RUTA_BASE` con el nombre `Informe_Ensayo_XXXXXXXX.pdf`.
- **PNGs**: Se guardan automáticamente en una carpeta `PNGs/` dentro de la `RUTA_BASE` del ensayo. Se sobreescriben cada vez que se ejecuta el script.

### Progreso en consola

El script muestra el progreso detallado en tiempo real:

```
============================================================
  TANDA 1/6: Cmg5tTan1ZNE.seed (7 Silex)
============================================================

  FASE 1: Alineando todos los Silex de la tanda...

  --- Par 1/40: Cmg5tTan1ZNE vs ES.SX244 ---
    Coef. envolventes: 0.7763 | Desfase grueso: 1800.0 ms
    Coef. fino: 0.6572 | Corrección: 807.1 ms
    Desfase final: 2607.1 ms

  FASE 2: Verificando desfases de la tanda...
    Desfase mediano de la tanda: 2622.0 ms

  FASE 3: Ajustando inicio del evento...
    Inicio real del evento: +8.16s desde la hora del .txt

  FASE 4: Generando gráficas...
    ✓ ES.SX244 completado (3 páginas)
```

---

## Procesado de señales

### 1. Carga y preprocesado

Cada señal (CMG y Silex) pasa por la siguiente cadena de procesado:

1. **Merge** de fragmentos del Silex (interpolación de gaps entre chunks).
2. **Igualación de frecuencias** de muestreo: se baja la señal de mayor frecuencia a la menor (downsampling, más limpio que inventar muestras).
3. **Corrección de polaridad**:
   - *CMG (fija)*: se invierten las componentes Z y E.
   - *Silex (dinámica)*: según la orientación indicada en el config.txt (N o S) se invierten las componentes correspondientes.
4. **Eliminación del offset** (demean) para centrar la señal en cero.
5. **Factor de conversión** a miligramos (mg) en el CMG (factor 2×10⁻⁴).
6. **Suavizado de bordes** (taper tipo Hamming, máximo 5 segundos).
7. **Filtro pasabanda 1–25 Hz** con fase cero (zerophase, 4 esquinas).

### 2. Alineación temporal automática

La alineación se realiza en varias fases para lograr precisión milimétrica:

**Fase 1 — Alineación gruesa (envolventes de energía combinada):**
- Se calcula la energía combinada de las 3 componentes: `√(Z² + N² + E²)`.
- Se obtiene la envolvente mediante la transformada de Hilbert.
- Se suaviza con media móvil (0.5 s).
- Se realiza una correlación cruzada normalizada entre las envolventes del CMG y del Silex en una ventana de ±30 s.

**Fase 2 — Ajuste fino (forma de onda):**
- Correlación cruzada de la componente Z en una ventana estrecha (±2 s) alrededor del resultado de la fase anterior.

**Fase 3 — Corrección de outliers:**
- Se calcula la mediana de los desfases de todos los Silex de la misma tanda.
- Si algún Silex se desvía más de 500 ms de la mediana, su desfase se corrige automáticamente al valor mediano.

### 3. Detección del inicio del evento

Se analiza la envolvente de la componente Z del CMG en ventanas de 1 segundo. Se busca la primera ventana donde la energía media supera el 10% del máximo **y las dos ventanas siguientes también lo superan** (energía sostenida, no picos de ruido aislados). Se retrocede 0.5 s desde el punto detectado para no perder el arranque del evento.

### 4. Generación de gráficas

Por cada par CMG–Silex se generan **3 gráficas** en formato 6×1 vertical (tamaño A4):

| Vista | Descripción |
|---|---|
| Señal completa | Evento completo con estadísticas (max, min, porcentaje de similitud) y líneas de referencia del CMG |
| Primeros 1.5 s | Detalle del inicio del evento |
| Primeros 5 s | Fase inicial del evento |

Cada gráfica muestra las 3 componentes (Z, N, E) con el CMG en azul y el Silex en negro, una debajo de la otra para facilitar la comparación visual.

---

## Ejemplo de salida

### Señal completa

<!-- Insertar captura de ejemplo de señal completa aquí -->
![Señal Completa](capturas/ejemplo_senal_completa.png)

### Primeros 1.5 segundos

<!-- Insertar captura de ejemplo de 1.5 segundos aquí -->
![Primeros 1.5s](capturas/ejemplo_1_5s.png)

### Primeros 5 segundos

<!-- Insertar captura de ejemplo de 5 segundos aquí -->
![Primeros 5s](capturas/ejemplo_5s.png)

---

## Parámetros modificables en el script

Si se desea ajustar el comportamiento del programa sin modificar la lógica interna, estos son los parámetros principales al inicio del archivo `script_informe.py`:

| Parámetro | Ubicación | Descripción | Valor por defecto |
|---|---|---|---|
| `RUTA_CONFIG` | Línea ~23 | Ruta al archivo de configuración | (debe configurarse) |
| `VISTAS` | Línea ~25 | Ventanas temporales a generar (nombre, segundos antes, segundos después) | Completa, 1.5s, 5s |

---

## Autor

**José David Conde Quispe**
Prácticas externas — Grado en Ingeniería Electrónica de Comunicaciones (GIEC)
Universidad Complutense de Madrid (UCM)
Instituto Geográfico Nacional (IGN) — 2026
