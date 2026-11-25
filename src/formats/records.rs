//! Fortran list-directed input, as `READ(u,*)` and `READ(u,'(A)')` see it.
//!
//! A list-directed read starts at the next record, takes values across as many
//! records as it needs, and discards the rest of the last record it touched.
//! Programs depend on that: `inversion` skips one record per pole, and a
//! `fort.2` read takes a layer count and then that many layers in one
//! statement.

use std::fmt;
use std::str::FromStr;

#[derive(Debug, PartialEq)]
pub struct Error {
    pub line: usize,
    pub message: String,
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "line {}: {}", self.line, self.message)
    }
}

impl std::error::Error for Error {}

pub struct Records<'a> {
    lines: std::str::Lines<'a>,
    line: usize,
}

impl<'a> Records<'a> {
    pub fn new(text: &'a str) -> Self {
        Self {
            lines: text.lines(),
            line: 0,
        }
    }

    /// The next record as text, or `None` at end of file.
    pub fn raw(&mut self) -> Option<&'a str> {
        let record = self.lines.next()?;
        self.line += 1;
        Some(record)
    }

    /// Starts one list-directed read. Dropping it discards what is left of
    /// the record it stopped in.
    pub fn read(&mut self) -> Read<'_, 'a> {
        Read {
            records: self,
            pending: Vec::new(),
        }
    }

    fn error(&self, message: String) -> Error {
        Error {
            line: self.line,
            message,
        }
    }
}

pub struct Read<'r, 'a> {
    records: &'r mut Records<'a>,
    /// Unread tokens of the current record, last token first.
    pending: Vec<&'a str>,
}

impl<'a> Read<'_, 'a> {
    pub fn f32(&mut self) -> Result<f32, Error> {
        self.value(|t| f32::from_str(&t.replace(['d', 'D'], "e")).ok())
    }

    pub fn i32(&mut self) -> Result<i32, Error> {
        self.value(|t| i32::from_str(t).ok())
    }

    pub fn count(&mut self) -> Result<usize, Error> {
        self.value(|t| usize::from_str(t).ok())
    }

    pub fn f32s(&mut self, n: usize) -> Result<Vec<f32>, Error> {
        (0..n).map(|_| self.f32()).collect()
    }

    fn value<T>(&mut self, parse: impl Fn(&str) -> Option<T>) -> Result<T, Error> {
        let token = self.token()?;
        parse(token).ok_or_else(|| self.records.error(format!("bad value {token:?}")))
    }

    fn token(&mut self) -> Result<&'a str, Error> {
        loop {
            if let Some(token) = self.pending.pop() {
                if token == "/" || token.contains('*') {
                    let message = format!("unsupported list-directed token {token:?}");
                    return Err(self.records.error(message));
                }
                return Ok(token);
            }
            let record = self
                .records
                .raw()
                .ok_or_else(|| self.records.error("unexpected end of file".to_owned()))?;
            self.pending = record
                .split([' ', ',', '\t'])
                .filter(|t| !t.is_empty())
                .rev()
                .collect();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_read_spans_records_and_drops_the_rest_of_the_last() {
        let mut records = Records::new("1 2\n3 4 5\n6\n");
        let mut read = records.read();
        let first: Vec<i32> = (0..3).map(|_| read.i32().unwrap()).collect();
        drop(read);
        assert_eq!(first, [1, 2, 3]);
        assert_eq!(records.read().i32().unwrap(), 6);
    }

    #[test]
    fn reals_take_fortran_spellings() {
        let mut records = Records::new("1. .5 0.1548e+03 2.5D-1 -20\n");
        assert_eq!(
            records.read().f32s(5).unwrap(),
            [1.0, 0.5, 154.8, 0.25, -20.0]
        );
    }

    #[test]
    fn a_repeat_count_is_refused_with_its_line() {
        let mut records = Records::new("1\n3*0.0\n");
        records.raw();
        assert_eq!(records.read().f32().unwrap_err().line, 2);
    }
}
