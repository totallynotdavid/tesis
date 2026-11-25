//! `kk prep` on real FDSN answers for the 2013-09-25 Yauca event: IU.ANMO
//! (poles and zeros in rad/s), II.PFO (in Hz) and IU.KOWA (no data returned),
//! cut from a `kk fetch` of EarthScope.

use std::fs;
use std::path::Path;
use std::process::Command;

use kk::fetch::stationxml;
use kk::formats::waveforms;

const FIXTURE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/yauca-2013");

/// Peak displacement in micrometres of each P window, from ObsPy 1.5.2:
/// `remove_response(output="DISP", pre_filt=(0.001, 0.002, 1, 2),
/// water_level=None, taper=False)` on the same demeaned 120 s window. ObsPy
/// uses every stage and a cosine band-pass, kk the first stage, the overall
/// sensitivity and its own band-pass, so they agree to within 1%.
const OBSPY: [(&str, f32); 2] = [("ANMO", 53.089), ("PFO", 54.888)];

#[test]
fn prep_converts_to_micrometres_of_displacement() {
    let dir = tempfile::tempdir().unwrap();
    for file in ["event.toml", "stations.xml", "waveforms.mseed"] {
        fs::copy(Path::new(FIXTURE).join(file), dir.path().join(file)).unwrap();
    }
    let output = Command::new(env!("CARGO_BIN_EXE_kk"))
        .arg("prep")
        .arg(dir.path().join("event.toml"))
        .output()
        .unwrap();
    assert!(output.status.success(), "{output:?}");
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(
        stdout.contains("IU.KOWA.00.BHZ  skipped: no data returned"),
        "{stdout}"
    );

    let text = fs::read_to_string(dir.path().join("fort.1")).unwrap();
    let fort1 = waveforms::parse(&text).unwrap();
    assert!(fort1.title.starts_with("Yauca 2013          disp."));
    assert_eq!(fort1.stations.len(), OBSPY.len());
    for station in &fort1.stations {
        let code = station.name[..12].trim();
        let (_, expected) = OBSPY.iter().find(|(c, _)| *c == code).unwrap();
        let ratio = station.peak / expected;
        assert!((ratio - 1.0).abs() < 0.01, "{code}: peak ratio {ratio}");
        assert_eq!(station.samples.len(), 240);
        assert_eq!(station.onset, 20);
    }
}

/// A poles-and-zeros stage has gain 1 at its normalisation frequency, in
/// whichever variable it was written; after conversion to rad/s that holds
/// at s = 2 pi i f.
#[test]
fn responses_have_unit_gain_at_their_normalisation_frequency() {
    let xml = fs::read_to_string(Path::new(FIXTURE).join("stations.xml")).unwrap();
    let doc = roxmltree::Document::parse(&xml).unwrap();
    let frequencies = doc
        .descendants()
        .filter(|n| n.has_tag_name("Channel"))
        .map(|channel| {
            let stage = channel
                .descendants()
                .find(|n| n.has_tag_name("PolesZeros"))
                .unwrap();
            let f = stage
                .children()
                .find(|n| n.has_tag_name("NormalizationFrequency"))
                .unwrap();
            f.text().unwrap().trim().parse::<f64>().unwrap()
        });
    let channels = stationxml::parse(&xml).unwrap();
    assert_eq!(channels.len(), 3);
    for (channel, f) in channels.iter().zip(frequencies) {
        let response = channel.response.as_ref().unwrap();
        let w = 2.0 * std::f64::consts::PI * f;
        let distance = |r: &kk::prep::fft::Complex| f64::from(r.re).hypot(w - f64::from(r.im));
        let gain = f64::from(response.normalisation)
            * response.zeros.iter().map(distance).product::<f64>()
            / response.poles.iter().map(distance).product::<f64>();
        assert!((gain - 1.0).abs() < 1e-3, "{}: gain {gain}", channel.id());
    }
}
