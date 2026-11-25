use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{Parser, Subcommand};

use kk::event::{Config, Files};
use kk::fetch::{self, miniseed, stationxml};
use kk::formats::waveforms::{self, Waveforms};
use kk::prep::{self, earth::Tables};

/// Kikuchi-Kanamori teleseismic body-wave inversion.
#[derive(Parser)]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Download the stations and records an event file selects, beside it.
    Fetch { event: PathBuf },
    /// Convert the fetched records into fort.1, beside the event file.
    Prep { event: PathBuf },
}

fn main() -> ExitCode {
    let result = match Cli::parse().command {
        Command::Fetch { event } => fetch(&event),
        Command::Prep { event } => prep(&event),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("kk: {e}");
            ExitCode::FAILURE
        }
    }
}

fn fetch(path: &Path) -> Result<(), String> {
    let config = config(path)?;
    let files = Files::beside(path);
    let tables = Tables::jb();
    let earth = tables.at(config.depth).map_err(|e| e.to_string())?;

    let xml = get(&fetch::station_url(&config))?.ok_or("no station matches the event file")?;
    let xml = String::from_utf8(xml).map_err(|e| format!("station response: {e}"))?;
    let channels = stationxml::parse(&xml).map_err(|e| format!("station response: {e}"))?;
    write(&files.stations, xml.as_bytes())?;

    let body = fetch::bulk(&config, &earth, &channels);
    let data = post(&fetch::dataselect_url(&config), &body)?.unwrap_or_default();
    write(&files.waveforms, &data)?;
    println!(
        "{} channels, {} bytes of miniSEED",
        channels.len(),
        data.len()
    );
    Ok(())
}

fn prep(path: &Path) -> Result<(), String> {
    let config = config(path)?;
    let files = Files::beside(path);
    let tables = Tables::jb();
    let earth = tables.at(config.depth).map_err(|e| e.to_string())?;
    let channels = stationxml::parse(&read_to_string(&files.stations)?)
        .map_err(|e| format!("{}: {e}", files.stations.display()))?;
    let traces = miniseed::read(&read(&files.waveforms)?)
        .map_err(|e| format!("{}: {e}", files.waveforms.display()))?;

    let event = config.event();
    let window = config.window();
    let mut stations = Vec::new();
    for channel in &channels {
        let station = fetch::record(channel, &traces)
            .and_then(|r| prep::convert(&event, &earth, &r, &window).map_err(|s| s.to_string()));
        match station {
            Ok(s) => {
                println!("{}  converted", channel.id());
                stations.push(s);
            }
            Err(e) => println!("{}  skipped: {e}", channel.id()),
        }
    }
    let waveforms = Waveforms {
        title: prep::title(&event, &earth),
        stations,
    };
    write(&files.fort1, waveforms::write(&waveforms).as_bytes())?;
    println!(
        "{} of {} channels written to {}",
        waveforms.stations.len(),
        channels.len(),
        files.fort1.display()
    );
    Ok(())
}

fn config(path: &Path) -> Result<Config, String> {
    Config::parse(&read_to_string(path)?).map_err(|e| format!("{}: {e}", path.display()))
}

/// The body of a successful answer, or `None` for 204 No Content, which FDSN
/// services send when nothing matches.
fn get(url: &str) -> Result<Option<Vec<u8>>, String> {
    body(ureq::get(url).call(), url)
}

fn post(url: &str, body_text: &str) -> Result<Option<Vec<u8>>, String> {
    body(ureq::post(url).send(body_text), url)
}

fn body(
    response: Result<ureq::http::Response<ureq::Body>, ureq::Error>,
    url: &str,
) -> Result<Option<Vec<u8>>, String> {
    const LIMIT: u64 = 1 << 30;
    let mut response = response.map_err(|e| format!("{url}: {e}"))?;
    if response.status() == 204 {
        return Ok(None);
    }
    response
        .body_mut()
        .with_config()
        .limit(LIMIT)
        .read_to_vec()
        .map(Some)
        .map_err(|e| format!("{url}: {e}"))
}

fn read(path: &Path) -> Result<Vec<u8>, String> {
    fs::read(path).map_err(|e| format!("{}: {e}", path.display()))
}

fn read_to_string(path: &Path) -> Result<String, String> {
    fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))
}

fn write(path: &Path, bytes: &[u8]) -> Result<(), String> {
    fs::write(path, bytes).map_err(|e| format!("{}: {e}", path.display()))
}
