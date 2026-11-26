//! gfortran's `F` and `I` edit descriptors, for the formatted reports.

/// `Fw.d`. gfortran rounds the exact binary value to `d` decimals, half to
/// even, as Rust does. It drops the leading zero of a value under 1 only when
/// the field has no room for it, keeps the sign of a negative value that
/// rounds to zero, and fills a field too narrow with asterisks. It spells an
/// infinity `Infinity` when the field holds it and `Inf` otherwise.
pub fn f(value: f32, width: usize, decimals: usize) -> String {
    let text = if value.is_nan() {
        "NaN".to_owned()
    } else if value.is_infinite() {
        let sign = if value < 0.0 { "-" } else { "" };
        let long = format!("{sign}Infinity");
        if long.len() <= width {
            long
        } else {
            format!("{sign}Inf")
        }
    } else {
        let text = format!("{value:.decimals$}");
        if text.len() > width {
            text.replacen("0.", ".", 1)
        } else {
            text
        }
    };
    if text.len() > width {
        return "*".repeat(width);
    }
    format!("{text:>width$}")
}

/// `Iw`.
pub fn i(value: i32, width: usize) -> String {
    let text = value.to_string();
    if text.len() > width {
        return "*".repeat(width);
    }
    format!("{text:>width$}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f_drops_the_leading_zero_only_without_room() {
        assert_eq!(f(0.5, 6, 3), " 0.500");
        assert_eq!(f(0.5, 5, 4), ".5000");
        assert_eq!(f(-0.5, 6, 4), "-.5000");
        assert_eq!(f(-0.04, 6, 1), "  -0.0");
        assert_eq!(f(1234.5, 5, 1), "*****");
    }

    #[test]
    fn f_spells_values_that_are_not_finite_as_gfortran_does() {
        assert_eq!(f(f32::NAN, 8, 4), "     NaN");
        assert_eq!(f(f32::INFINITY, 9, 3), " Infinity");
        assert_eq!(f(f32::NEG_INFINITY, 8, 4), "    -Inf");
        assert_eq!(f(f32::INFINITY, 2, 1), "**");
    }

    #[test]
    fn i_fills_a_narrow_field_with_asterisks() {
        assert_eq!(i(-30, 5), "  -30");
        assert_eq!(i(123_456, 5), "*****");
    }
}
