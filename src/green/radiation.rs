//! The far-field radiation pattern of a double couple, `sub.radp.f`.

use super::RAD;

/// Amplitudes of the down- and up-going P, SV and SH waves. SH is positive
/// clockwise and SV along SH × P.
pub struct Radiation {
    pub pd: f32,
    pub pu: f32,
    pub svd: f32,
    pub svu: f32,
    pub shd: f32,
    pub shu: f32,
}

/// The pattern at azimuth `az` from strike, for a fault of `dip` and `rake`
/// in degrees, and a ray of slowness `p` leaving a layer of velocities `vp`,
/// `vs`. A dip above 360 selects the isotropic source.
pub fn radp(az: f32, dip: f32, rake: f32, p: f32, vp: f32, vs: f32) -> Radiation {
    if dip > 360.0 {
        return Radiation {
            pd: 1.0,
            pu: 1.0,
            svd: 0.0,
            svu: 0.0,
            shd: 0.0,
            shu: 0.0,
        };
    }
    let dr0 = dip * RAD;
    let ar0 = rake * RAD;
    let thet = az * RAD;
    let sih = p * vp;
    let cih = if sih >= 1.0 {
        0.0
    } else {
        (1.0 - sih * sih).sqrt()
    };
    let sih2 = p * vs;
    let cih2 = (1.0 - sih2 * sih2).sqrt();
    let c1 = sih * sih;
    let c2 = -((2.0 * sih) * cih);
    let c3 = 2.0 - 3.0 * (sih * sih);
    let sv1 = sih2 * cih2;
    let sv2 = 2.0 * (sih2 * sih2) - 1.0;
    let sv3 = -((3.0 * sih2) * cih2);
    let (sh1, sh2) = (sih2, cih2);
    let (sth, cth) = (thet.sin(), thet.cos());
    let (sth2, cth2) = ((thet * 2.0).sin(), (thet * 2.0).cos());
    let (sd0, cd0) = (dr0.sin(), dr0.cos());
    let (sd02, cd02) = ((dr0 * 2.0).sin(), (dr0 * 2.0).cos());
    let (sa0, ca0) = (ar0.sin(), ar0.cos());
    let a1 = ((sth2 * ca0) * sd0) + (((0.5 * cth2) * sa0) * sd02);
    let a2 = ((cth * ca0) * cd0) - ((sth * sa0) * cd02);
    let a3 = (0.5 * sa0) * sd02;
    let a4 = ((cth2 * ca0) * sd0) - (((0.5 * sth2) * sa0) * sd02);
    let a5 = ((sth * ca0) * cd0) + ((cth * sa0) * cd02);
    Radiation {
        pd: (a1 * c1 + a2 * c2) + a3 * c3,
        pu: (a1 * c1 - a2 * c2) + a3 * c3,
        svd: (a1 * sv1 + a2 * sv2) + a3 * sv3,
        svu: ((-(a1 * sv1)) + a2 * sv2) - a3 * sv3,
        shd: a4 * sh1 + a5 * sh2,
        shu: a4 * sh1 - a5 * sh2,
    }
}
