# Event files

```bash
kk fetch events/yauca-2013/event.toml   # writes stations.xml and waveforms.mseed
kk prep events/yauca-2013/event.toml    # writes fort.1
kk green events/yauca-2013/event.toml   # writes fort.3
```

Each command reads the event file and keeps its files in the event file's
directory. Git tracks only the event file. `kk prep` prints one line per
channel:

```text
II.ASCN.00.BHZ  skipped: no data returned
II.CMLA.00.BHZ  converted
II.EFI.00.BHZ  converted
...
35 of 39 channels written to events/yauca-2013/fort.1
```

## Fields

[events/yauca-2013/event.toml](../events/yauca-2013/event.toml) is the
thesis event:

```toml
title = "Yauca 2013"
origin = "2013-09-25T16:42:43.5"
latitude = -15.84
longitude = -74.51
depth = 40

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

`kk green` needs the `[structure]` table, which is the `fort.2` velocity
structure. The Yauca file uses the structure of the reference's 2001 Arequipa
example:

```toml
[structure]
title = "2 layers+3km water"
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
35 stations, 1 depths, written to events/yauca-2013/fort.3
```
