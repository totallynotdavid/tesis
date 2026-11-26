#set page(
  paper: "a4",
  margin: (left: 25mm, right: 25mm, top: 40mm, bottom: 25mm),
)
#set text(lang: "es", size: 11pt)
#set par(leading: 0.75em, first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")
#show link: set text(fill: blue)

// UNMSM logo, absolute from the page corner like the textpos block did.
#place(
  top + left,
  dx: -14mm,
  dy: -24mm,
  image("unmsm-logo.png", width: 20mm),
)

#let unmsm-line = text.with(weight: "bold", size: 12pt)
#align(center)[
  #unmsm-line[UNIVERSIDAD NACIONAL MAYOR DE SAN MARCOS]
  #v(8pt)
  #unmsm-line[FACULTAD DE CIENCIAS FÍSICAS]
  #v(8pt)
  #unmsm-line[Escuela Profesional de Física]
  #v(8pt)
  #unmsm-line[PROYECTO DE TESIS DE LICENCIATURA]
]
#v(10pt)

= Responsabilidad de Ejecución
Tesista: David Duran \
Asesor: Cesar Jimenez

= Título del Trabajo:
Aspectos físicos de la fuente del terremoto de Yauca-Arequipa de 2013 (7.1 Mw)

= Resumen
En este trabajo se calculan las ecuaciones...

= Contenido

== Introducción

#heading(level: 3, numbering: none)[Estructura sismotectónica del Perú]

La interacción entre la placa oceánica de Nazca y la placa continental sudamericana define la estructura geológica del sur del Perú. La subducción continua de la placa de Nazca, con un espesor promedio de 6.4 km#footnote[*Nota*: Este promedio es para el área de estudio del artículo de Krabbenhöft et al. (2004)] @Krabbenhoft2004, bajo la placa sudamericana a lo largo de la Fosa de Perú-Chile @Krabbenhoft2004 ocurre a una velocidad aproximada de 63 mm/año @VillegasLanza2016.

Las velocidades sísmicas de la corteza continental del sur del Perú varían de 4.5 a 6.5 km/s hasta una profundidad de 25 km, según los modelos de Krabbenhöft et al. (2004), lo que refleja un cambio en el gradiente de velocidad del basamento cristalino asociado al proceso de subducción @Krabbenhoft2004. La deformación cortical resulta en la elevación de la Cordillera Occidental y la altiplanicie, estructuras que son producto de procesos tectónicos milenarios @VillegasLanza2016.

#heading(level: 3, numbering: none)[Antecedentes de trabajos previos]

Arequipa, conocida por su actividad sísmica, proporciona un contexto esencial para comprender el terremoto de Yauca en 2013. Este análisis se basa en el trabajo de E. Silgado, "Historia de los sismos más notables ocurridos en el Perú (1513–1974)" @silgado_ferro_historia_1978. Desde la antigüedad, Arequipa ha sido escenario de terremotos significativos.

Por ejemplo, en la época de Túpac Yupanqui, un gran terremoto, posiblemente relacionado con la erupción del volcán Misti, arrasó los primeros asentamientos. En 1555, ya bajo el dominio español, los registros indican que los temblores eran frecuentes y destructivos. Avanzando en el tiempo, el terremoto de 1604, que generó un tsunami, causó estragos desde Arica hasta Ica, afectando gravemente localidades como el puerto de Ilo y la Villa de Camaná. En 1784, otro terremoto, esta vez de una magnitud estimada en 8.0, duró cerca de cinco minutos y cobró la vida de 154 personas.

El siglo XX también fue testigo de terremotos devastadores. Uno notable ocurrió el 6 de agosto de 1913, afectando a ciudades como Chala, Atico y Ocoña, entre otras; este evento fue seguido por un maremoto que inundó partes de la costa. En 1958, un sismo de magnitud 7.3 resultó en 28 muertes y daños significativos en estructuras antiguas. En 1960, un terremoto y sus réplicas causaron aún más daño y pérdidas humanas en Arequipa. Más recientemente, el 25 de septiembre de 2015, un sismo de magnitud 7.0 Mw sacudió la región, con su epicentro reportado cerca de la ciudad de Yauca.

#heading(level: 3, numbering: none)[Planteamiento del problema]

El terremoto de Yauca-Arequipa de 2018 presenta un proceso de ruptura complejo que requiere un análisis detallado para entender las características físicas de la fuente sísmica y cómo contribuyeron a la generación de un tsunami local @CEN13.

#heading(level: 3, numbering: none)[Planteamiento de la hipótesis]

Mediante la inversión de las formas de onda telesísmicas, se puede obtener una comprensión detallada de las características de la fuente sísmica y el proceso de ruptura del terremoto.

#heading(level: 3, numbering: none)[Objetivos]

Principales:

+ Analizar las características físicas de la fuente sísmica del terremoto de Yauca-Arequipa de 2018 y entender cómo contribuyeron al tsunami local generado.

Secundarios:

+ Automatizar el proceso de inversión de formas de onda telesísmicas.
+ Comparar estos hallazgos con los eventos sísmicos históricos en la región para proporcionar una mejor comprensión de la actividad sísmica en esta área.

== Fundamento teórico

La inversión de formas de onda telesísmicas es una técnica crucial utilizada para analizar el proceso de ruptura de terremotos. Esta técnica permite deducir información sobre el mecanismo focal, que está relacionado con la orientación del plano de ruptura y es coherente con el patrón sismotectónico de la región sur del Perú.

== Metodología

Utilizando el método de inversión de formas de onda telesísmicas de Kikuchi y Kanamori @Kik2003, se analizarán las señales sísmicas de estaciones de la red IRIS. Este método se basa en la deconvolución iterativa de sismogramas para analizar el proceso de ruptura de los terremotos. El método inicia con la modelación de un terremoto como una secuencia de subeventos, cada uno caracterizado por su propio tensor de momento sísmico $M_(i j)$, descompuesto en componentes esféricas y deviatorias:

$
  (1/3) "Tr"[M] & = (M_(1 1) + M_(2 2) + M_(3 3)) / 3 \
        M_(i j) & = M_0 (n_i d_j + n_j d_i)
$

donde $arrow(n)$ es el vector normal al plano de falla y $arrow(d)$ es el vector de deslizamiento.

El tensor de momento sísmico se representa mediante una combinación lineal de seis tensores elementales, que forman la base para representar la fuente sísmica.

La deconvolución se realiza ajustando los sismogramas observados $x_j (t)$ a los sismogramas sintéticos mediante la minimización del error $Delta$:

$
  Delta = sum_(j=1)^(N_s) integral [x_j (t) - sum_(n=1)^(N_b) a_n omega_(j n) (t; p)]^2 dif t
$

Los coeficientes $a_n$ y los parámetros $p$ se determinan minimizando $Delta$. Las ecuaciones normales para la solución son:

$ a_n = sum_(m=1)^(N_b) R_(n m)^(-1) G_m $

El modelo se evalúa a través de la correlación $psi_M (p)$ entre las formas de onda observadas y sintéticas:

$ psi_M (p) = (sum_n G_n a_n^0) / R_x $

== Cronograma

El cronograma del proyecto de tesis se muestra en la Tabla.

#let gantt-bar = rgb("#7fa0d9")
#let gantt-group = rgb("#b8c9e8")
#let months = (
  "Ene",
  "Feb",
  "Mar",
  "Abr",
  "May",
  "Jun",
  "Jul",
  "Ago",
  "Sep",
  "Oct",
  "Nov",
  "Dic",
)
#let cell(months-active) = (
  months
    .enumerate()
    .map(((i, _)) => table.cell(
      fill: if months-active.contains(i + 1) { gantt-bar } else { none },
    )[])
)
#let group-row(name, months-active) = (
  table.cell(fill: gantt-group)[*#name*],
  ..months
    .enumerate()
    .map(((i, _)) => table.cell(
      fill: if months-active.contains(i + 1) { gantt-group } else { none },
    )[]),
)
#let bar-row(name, months-active) = (table.cell[#name], ..cell(months-active))

#figure(
  table(
    columns: (auto,) + (1fr,) * 12,
    align: (left,) + (center,) * 12,
    stroke: 0.4pt,
    table.header([], ..months.map(m => [*#m*])),
    ..group-row("Objetivo Principal", (1, 2, 3, 4, 5, 6, 7, 8)),
    ..bar-row("Revisión de la literatura", (1, 2)),
    ..bar-row("Recolección y preparación de datos", (3, 4)),
    ..bar-row("Implementación del modelo", (5, 6)),
    ..bar-row("Análisis de datos", (7, 8)),
    ..group-row("Objetivos secundarios", (3, 4, 5, 6, 7, 8, 9, 10)),
    ..bar-row("Automatización de inversión", (3, 4, 5)),
    ..bar-row("Comparación con eventos históricos", (6, 7, 8, 9, 10)),
    ..group-row("Redacción y publicación", (9, 10, 11, 12)),
    ..bar-row("Redacción de la tesis", (9, 10)),
    ..bar-row("Preparación de publicación", (11, 12)),
  ),
  caption: [Cronograma del proyecto de tesis, 2023.],
)

= Nivel del Trabajo de Investigación:
Trabajo de Tesis para obtener el Título Profesional en Física.

#bibliography("refs.bib", style: "ieee")

#align(center, line(length: 100mm, stroke: 0.2mm))

#let today-es = {
  let d = datetime.today()
  let month = (
    "enero",
    "febrero",
    "marzo",
    "abril",
    "mayo",
    "junio",
    "julio",
    "agosto",
    "septiembre",
    "octubre",
    "noviembre",
    "diciembre",
  ).at(d.month() - 1)
  [#d.day() de #month de #d.year()]
}

#v(40pt)
#pad(left: 280pt)[Lima, #today-es]

#v(100pt)
#align(center)[... #h(100pt) ...]

#align(center)[(Tesista) #h(140pt) (Asesor)]
