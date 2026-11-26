//! The beachballs read the reference's tensor and mechanism conventions: on
//! the golden Denali `fort.4`, each subevent's nodal planes fall where its
//! tensor radiates no P wave, and its T axis is compressional.

use std::fs;
use std::path::Path;

use kk::formats::plot;
use kk::plot::sphere;

#[test]
fn nodal_planes_and_t_axes_agree_with_the_denali_tensors() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("reference/golden/denali/fort.4");
    let p = plot::read(&fs::read(path).unwrap()).unwrap();
    let mechanisms = p
        .subevents
        .iter()
        .map(|s| (s.mechanism, s.tensor))
        .chain([(p.total, p.tensor)]);
    for (mechanism, tensor) in mechanisms {
        let m = sphere::tensor(&tensor);
        let (normal, slip) = sphere::normal_and_slip(mechanism);
        let t_axis = std::array::from_fn(|i| normal[i] + slip[i]);
        let t_amplitude = sphere::p_amplitude(&m, unit(t_axis));
        assert!(t_amplitude > 0.0, "{mechanism:?}: T axis at {t_amplitude}");
        for pole in [normal, slip] {
            for (right, up) in sphere::great_circle(pole, 36) {
                let on_plane = sphere::p_amplitude(&m, sphere::unproject(right, up));
                assert!(
                    on_plane.abs() < 0.25 * t_amplitude,
                    "{mechanism:?}: {on_plane} on a nodal plane, {t_amplitude} on T"
                );
            }
        }
    }
}

fn unit(v: sphere::Vector) -> sphere::Vector {
    let norm = v.iter().map(|c| c * c).sum::<f64>().sqrt();
    v.map(|c| c / norm)
}
