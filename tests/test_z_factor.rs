#[cfg(test)]
mod tests {
    // use super::*;
    use fluid_props::z_factor::{z_factor, ZCorrelations};

    #[test]
    fn test_z_factor() {
        let test_cases = [
            (2.0, 0.147, 0.9950397850047862, 0.994616021885175),
            (1.6, 8.82, 1.0548103327841676, 1.0521119421909435),
        ];

        for (t_pr, p_pr, hall_z, dranchuck_z) in test_cases {
            assert_eq!(
                hall_z,
                z_factor(t_pr, p_pr, ZCorrelations::HallAndYarborough, 1e-6)
            );
            assert_eq!(
                dranchuck_z,
                z_factor(t_pr, p_pr, ZCorrelations::DranchkAndAbouKassem, 1e-6)
            );
        }
    }
}
