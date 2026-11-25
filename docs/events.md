# Event files

```bash
kk fetch events/yauca-2013/event.toml   # writes stations.xml and waveforms.mseed
kk prep events/yauca-2013/event.toml    # writes fort.1
```

Both commands read the event file and keep their files in its directory. Git
tracks only the event file. `kk prep` prints one line per channel:

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
