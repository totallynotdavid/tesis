# Event files

```bash
kk fetch events/yauca-2013/event.toml   # writes stations.xml and waveforms.mseed
kk prep events/yauca-2013/event.toml    # writes fort.1
kk green events/yauca-2013/event.toml   # writes fort.3
kk invert events/yauca-2013/event.toml  # writes fort.16 and fort.4
kk plot events/yauca-2013/event.toml    # writes solution.svg and waveforms.svg
kk sweep events/yauca-2013/event.toml structures/*.toml  # one run per structure
```

Each command reads the event file and keeps its files in the event file's
directory. Git tracks only the event file. `kk prep` prints one line per
channel:

```text
II.ASCN.00.BHZ  skipped: no data returned
II.CMLA.00.BHZ  converted
II.EFI.00.BHZ  converted
...
IU.RCBR.00.BHZ  skipped: amplitude 22 times the network median
...
34 of 39 channels written to events/yauca-2013/fort.1
```

A station is skipped when its peak, corrected for geometrical spreading, is
more than 5 times the median of the stations on its component. P radiation
cannot put a working channel that far above the rest.

## Fields

[events/](../events/) holds the three M7.1–7.2 thrusts off Atiquipa:
`yauca-2013`, `yauca-2018` and `yauca-2024`. The 2013 file starts:

```toml
title = "Yauca 2013"
origin = "2013-09-25T16:42:43.5"
latitude = -15.84
longitude = -74.51
depth = 40
structure = "../../structures/arequipa-2001.toml"

[p]
lead = 10
duration = 120
dt = 0.5
low = 0.002
high = 1.0
```

| Field | Meaning | Default |
| --- | --- | --- |
| `title` | The `fort.1` title. The programs keep 20 characters. | required |
| `origin` | Origin time, UTC, ISO 8601. | required |
| `latitude`, `longitude` | Epicentre in degrees. | required |
| `depth` | Source depth in km. It selects the J-B travel-time column. | required |
| `structure` | The structure file, relative to the event file. `kk green` needs it. | none |
| `output` | `displacement`, `velocity` or `acceleration`. | `displacement` |
| `distances` | Epicentral distance range in degrees. | `[30, 90]` |
| `p.lead` | Seconds kept before the P onset. | required |
| `p.duration` | Seconds kept from the window start. | required |
| `p.dt` | Output sample interval in seconds, a whole multiple of the record's. | required |
| `p.low`, `p.high` | Band-pass corners in Hz; 0 leaves that side open. | `0` |

The `[stations]` table selects the channels, in FDSN query syntax:

| Field | Default |
| --- | --- |
| `service` | `https://service.earthscope.org` |
| `networks` | `IU,II` |
| `location` | `00` |
| `channel` | `BHZ` |

## Green's functions

`kk green` needs a structure file, which is the `fort.2` velocity structure.
The event file names it with a path relative to itself:

```toml
structure = "../../structures/arequipa-2001.toml"
```

[structures/](../structures/) holds four:

| File | Structure |
| --- | --- |
| `arequipa-2001.toml` | The reference's 2001 Arequipa example: 3 km of water over one 15 km crustal layer. |
| `prem-ocean.toml` | PREM's crust under 3 km of water at the source, PREM's continental crust under the receivers. |
| `ak135.toml` | ak135's 35 km continental crust at the source and under the receivers. |
| `crust1-atiquipa.toml` | CRUST1.0's cell holding the Atiquipa epicentres at the source, 40 km of crust; ak135 under the receivers. |

A structure file holds these fields:

```toml
title = "Arequipa 2001: 2 layers+3km water"
tq_p = 1.0
tq_s = 4.0
source = [[1.50, 0.00, 1.00, 3.0], [6.50, 3.74, 2.87, 15.0], [7.80, 4.40, 3.30, 0.0]]
receiver = [[5.57, 3.36, 2.65, 15.0], [6.50, 3.74, 2.87, 18.0], [7.80, 4.40, 3.30, 0.0]]
bounce = [[6.50, 3.74, 2.87, 0.0]]
```

| Field | Meaning |
| --- | --- |
| `title` | The `fort.3` title. 40 characters are kept. |
| `tq_p`, `tq_s` | Attenuation t* for P and S, in seconds. |
| `source`, `receiver`, `bounce` | Layers at the source, under the receivers and at the PP bounce point, top down. Each is `[vp, vs, density, thickness]` in km/s, g/cm3 and km. The last layer is the half-space; its thickness is ignored. A source layer with `vs = 0` is water. |

The `[green]` table sets the sampling and the line of source points. Every
field is optional:

| Field | Meaning | Default |
| --- | --- | --- |
| `nt` | Samples per trace, a power of two. | `256` |
| `dt` | Sample interval in seconds. | `p.dt` |
| `h0` | Depth of the reference point in km. | `depth` |
| `spacing` | Distance between points along the line, in km. | `5` |
| `points` | Number of points. | `1` |
| `reference` | The point at `h0`, counted from 1. | `1` |
| `dip` | Dip of the line in degrees; 90 is vertical. | `90` |

`kk green` prints one line per `fort.1` station. Like the reference, it skips
a station whose onset is past its last sample:

```text
CMLA  computed
EFI  computed
FFC  computed
...
34 stations, 1 depths, written to events/yauca-2013/fort.3
```

## Inversion

`kk invert` needs the `[inversion]` table, which holds the `i_inversion`
parameters:

```toml
[inversion]
length = 60
max_velocity = 4.0
time_function = { shape = "trapezoid", t1 = 2, t2 = 4 }
spans = [20, 20]
basis = "double-couple"
shift = -10
grid = { line = { azimuth = 310, spacing = 10, count = 5, reference = 3 } }
```

| Field | Meaning | Default |
| --- | --- | --- |
| `title` | The `fort.16` title. 40 characters are kept. | `title` |
| `length` | Seconds of each record fitted. | required |
| `max_velocity` | A subevent starts no earlier than its distance from the hypocentre over this velocity, in km/s. | required |
| `time_function` | The source time function: `shape`, and `t1`, `t2` in seconds. | required |
| `spans` | The latest onset searched for each subevent, in seconds; one value per subevent. | required |
| `basis` | The moment-tensor elements each subevent combines. | required |
| `shift` | Samples added to every station's onset. | `0` |
| `weight` | The weight of each station `kk green` kept. | `1` |
| `weights` | Weights by station code, such as `{ CMLA = 0.5, RPN = 0 }`. Weight 0 leaves the station out. A code `kk green` did not keep is an error. | `{}` |
| `grid` | The source points searched at each depth of `fort.3`. | required |

The time function `shape` is one of:

| `shape` | Function |
| --- | --- |
| `impulse` | One sample of area 1. |
| `trapezoid` | Rise `t1`, top `t2 - t1`, fall `t1`. A ramp of rise `t1` if `t2 < t1`. |
| `cosine` | The trapezoid with cosine flanks. |
| `triangle` | A triangle peaking at `t1` and ending at `t2`. |

The `basis` is one of:

| `basis` | Elements |
| --- | --- |
| `deviatoric` | The five deviatoric elements. |
| `double-couple` | The five deviatoric elements, constrained to a double couple. |
| `strike-slip` | Mxy and Mxx-Myy: vertical strike-slip. |
| `moment-tensor` | All six elements, with the isotropic part. |
| `vertical-double-couple` | Four elements, constrained to a double couple with a vertical nodal plane. |

The `grid` is a list of points or a line:

```toml
grid = { points = [[0, 0], [30, 120]] }   # [distance km, azimuth deg] from the hypocentre
grid = { line = { azimuth = 310, spacing = 10, count = 5, reference = 3 } }
```

A line has `count` points `spacing` km apart along `azimuth`. Point
`reference`, counted from 1, is at the hypocentre; it defaults to `1`.

`kk invert` writes the subevents to `fort.16` and the plot data to `fort.4`.
It prints the moment magnitude and the residual variance:

```text
 Mw, Var. =    6.88  0.4956
```

## Sweep

`kk sweep` runs `green`, `invert` and `plot` once per structure file. Each
run writes its files into a directory beside the event file, named for the
structure file, and prints one line:

```bash
kk sweep events/yauca-2013/event.toml structures/*.toml
```

```text
USGS Mww: M0 5.00e19 N m, Mw 7.07, 125/67/87, 45 km
structure               Mw M0/e19 ratio    var. strike/dip/rk  kagan
ak135                 6.91   2.96  0.59  0.3941     132/55/96     14
arequipa-2001         6.88   2.59  0.52  0.4956     128/53/98     17
crust1-atiquipa       6.92   2.99  0.60  0.4138    120/53/100     22
prem-ocean            6.99   3.82  0.76  0.5262    129/52/101     19
```

`ratio` is the total moment over the catalogue's. `kagan` is the smallest
rotation, in degrees, between the total double couple and the catalogue's.
Both need the `[catalogue]` table:

```toml
[catalogue]
source = "USGS Mww"
moment = 5.0e19
depth = 45
strike = 125
dip = 67
rake = 87
```

| Field | Meaning |
| --- | --- |
| `source` | The catalogue and solution type, printed on the first line. |
| `moment` | Scalar moment in N m. |
| `depth` | Centroid depth in km. |
| `strike`, `dip`, `rake` | Either nodal plane, in degrees. |

## Plot

`kk plot` reads `fort.4` and writes two figures:

- `solution.svg`: the moment rate, the grid with each subevent's point, and
  one beachball per subevent and for the total. A ball's area follows its
  moment. Open circles are stations on the lower hemisphere.
- `waveforms.svg`: each station's observed record in black and the synthetic
  in red, eight to a column. Every record of one phase type shares one
  amplitude scale, so a station that dwarfs the rest stands out. The number
  under the name is the observed peak-to-peak.
