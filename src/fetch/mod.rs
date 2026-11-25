//! The FDSN requests that collect an event's records, and the join of their
//! answers into [`Record`]s. The HTTP calls themselves live in `main`.

pub mod miniseed;
pub mod stationxml;

use libmseed::NsTime;
use libmseed::time::{Subseconds, TimeFormat};

use crate::event::Config;
use crate::prep::Record;
use crate::prep::earth::{Earth, Phase};
use miniseed::Trace;
use stationxml::Channel;

/// Seconds of record requested on each side of the window, so a late onset
/// or a picking change still finds samples.
const MARGIN: f32 = 60.0;

/// The fdsnws-station query for every channel the event's settings select
/// within its distance range, with responses.
pub fn station_url(config: &Config) -> String {
    let s = &config.stations;
    let origin = format_time(config.origin);
    format!(
        "{}/fdsnws/station/1/query?level=response&format=xml\
         &network={}&location={}&channel={}\
         &startbefore={origin}&endafter={origin}\
         &latitude={}&longitude={}&minradius={}&maxradius={}",
        s.service,
        s.networks,
        location(&s.location),
        s.channel,
        config.latitude,
        config.longitude,
        config.distances[0],
        config.distances[1],
    )
}

pub fn dataselect_url(config: &Config) -> String {
    format!("{}/fdsnws/dataselect/1/query", config.stations.service)
}

/// The fdsnws-dataselect POST body: one line per channel, covering its P
/// window with [`MARGIN`] on each side.
pub fn bulk(config: &Config, earth: &Earth, channels: &[Channel]) -> String {
    let mut body = String::new();
    for c in channels {
        let distance = earth
            .geometry(
                (config.latitude, config.longitude),
                (c.latitude, c.longitude),
            )
            .distance;
        let travel = earth.travel(Phase::P, distance);
        let start = travel - config.p.lead - MARGIN;
        let end = start + config.p.duration + 2.0 * MARGIN;
        body += &format!(
            "{} {} {} {} {} {}\n",
            c.network,
            c.station,
            location(&c.location),
            c.code,
            format_time(offset(config.origin, start)),
            format_time(offset(config.origin, end)),
        );
    }
    body
}

/// The channel's record, from the one trace that carries its id.
pub fn record(channel: &Channel, traces: &[Trace]) -> Result<Record, String> {
    let response = channel.response.clone()?;
    let id = channel.id();
    let trace = traces
        .iter()
        .find(|t| t.id == id)
        .ok_or("no data returned")?;
    let [segment] = trace.segments.as_slice() else {
        return Err(format!(
            "{} segments: the record has gaps",
            trace.segments.len()
        ));
    };
    Ok(Record {
        code: channel.station.clone(),
        channel: channel.code.clone(),
        latitude: channel.latitude,
        longitude: channel.longitude,
        start: segment.start,
        dt: segment.dt,
        samples: segment.samples.clone(),
        response,
    })
}

/// FDSN writes the empty location code as `--`.
fn location(code: &str) -> &str {
    if code.is_empty() { "--" } else { code }
}

fn offset(time: NsTime, seconds: f32) -> NsTime {
    time.add_nanos((f64::from(seconds) * 1e9).round() as i64)
}

fn format_time(time: NsTime) -> String {
    time.format(TimeFormat::IsoMonthDay, Subseconds::Micro)
}
