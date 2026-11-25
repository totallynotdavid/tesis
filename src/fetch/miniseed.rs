//! miniSEED records merged into one trace per channel.

use libmseed::{ParseFlags, ReadOptions, Samples, TraceList, sid};

use crate::prep::Time;

pub struct Trace {
    /// `NET.STA.LOC.CHA`, as [`super::stationxml::Channel::id`] gives it.
    pub id: String,
    /// The continuous runs of samples, in time order. More than one means the
    /// record has gaps.
    pub segments: Vec<Segment>,
}

pub struct Segment {
    pub start: Time,
    pub dt: f32,
    pub samples: Vec<f32>,
}

pub fn read(bytes: &[u8]) -> Result<Vec<Trace>, String> {
    let mut list = TraceList::new();
    list.read_buffer(bytes, &ReadOptions::new(ParseFlags::UNPACK_DATA))
        .map_err(|e| e.to_string())?;
    list.ids()
        .iter()
        .map(|id| {
            let nslc = sid::sid_to_nslc(&id.sid).map_err(|e| e.to_string())?;
            let segments = id
                .segments
                .iter()
                .map(|s| {
                    Ok(Segment {
                        start: Time::from_nanos(s.start_time.as_nanos()),
                        dt: (1.0 / s.sample_rate) as f32,
                        samples: samples(s.samples.as_ref())
                            .ok_or_else(|| format!("{}: no numeric samples", id.sid))?,
                    })
                })
                .collect::<Result<_, String>>()?;
            Ok(Trace {
                id: format!(
                    "{}.{}.{}.{}",
                    nslc.network, nslc.station, nslc.location, nslc.channel
                ),
                segments,
            })
        })
        .collect()
}

fn samples(samples: Option<&Samples>) -> Option<Vec<f32>> {
    match samples? {
        Samples::Int32(v) => Some(v.iter().map(|&x| x as f32).collect()),
        Samples::Float32(v) => Some(v.clone()),
        Samples::Float64(v) => Some(v.iter().map(|&x| x as f32).collect()),
        Samples::Text(_) => None,
    }
}
