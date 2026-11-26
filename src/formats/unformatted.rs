//! gfortran unformatted sequential files: each record is its bytes between
//! two copies of their length as a 4-byte little-endian integer.

use std::fmt;

#[derive(Debug, PartialEq)]
pub struct Error {
    /// The byte offset of the record that could not be read.
    pub offset: usize,
    pub message: String,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "byte {}: {}", self.offset, self.message)
    }
}

impl std::error::Error for Error {}

/// A `CHARACTER*width` value.
pub(super) fn padded(text: &str, width: usize) -> Vec<u8> {
    format!("{text:<width$.width$}").into_bytes()
}

pub(super) fn record(out: &mut Vec<u8>, bytes: &[u8]) {
    let length = (bytes.len() as u32).to_le_bytes();
    out.extend(length);
    out.extend(bytes);
    out.extend(length);
}

pub(super) fn reals(values: &[f32]) -> Vec<u8> {
    values.iter().flat_map(|v| v.to_le_bytes()).collect()
}

pub(super) struct Reader<'a> {
    bytes: &'a [u8],
    offset: usize,
}

impl<'a> Reader<'a> {
    pub(super) fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, offset: 0 }
    }

    pub(super) fn at_end(&self) -> bool {
        self.offset >= self.bytes.len()
    }

    /// The next record, which must hold `length` bytes.
    pub(super) fn next(&mut self, length: usize) -> Result<Vec<u8>, Error> {
        let rest = &self.bytes[self.offset..];
        let marker = |at: usize| rest.get(at..at + 4).map(|b| u32::from_le_bytes(word(b)));
        let found = marker(0).ok_or_else(|| self.error("truncated record".into()))?;
        if found as usize != length {
            return Err(self.error(format!("record of {found} bytes, expected {length}")));
        }
        if marker(4 + length) != Some(found) {
            return Err(self.error("record length markers differ".into()));
        }
        let data = rest[4..4 + length].to_vec();
        self.offset += 8 + length;
        Ok(data)
    }

    pub(super) fn next_text(&mut self, width: usize) -> Result<String, Error> {
        Ok(String::from_utf8_lossy(&self.next(width)?).into_owned())
    }

    pub(super) fn next_reals(&mut self, count: usize) -> Result<Vec<f32>, Error> {
        let bytes = self.next(4 * count)?;
        Ok(bytes
            .as_chunks::<4>()
            .0
            .iter()
            .map(|b| f32::from_le_bytes(*b))
            .collect())
    }

    /// A count read from a record, which must not be negative.
    pub(super) fn count(&self, value: i32) -> Result<usize, Error> {
        usize::try_from(value).map_err(|_| self.error(format!("count {value} is negative")))
    }

    pub(super) fn error(&self, message: String) -> Error {
        Error {
            offset: self.offset,
            message,
        }
    }
}

/// The 4-byte word at the start of `bytes`.
pub(super) fn word(bytes: &[u8]) -> [u8; 4] {
    bytes[..4].try_into().unwrap()
}
