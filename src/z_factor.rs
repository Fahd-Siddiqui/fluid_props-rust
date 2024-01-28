#[derive(PartialEq)]
pub enum ZCorrelations {
    HallAndYarborough,
    DranchkAndAbouKassem,
}

fn hall_and_yarborough(t_pr: f64, p_pr: f64, tolerance: f64) -> f64 {
    let mut fy: f64 = 1.0;
    let mut fdy: f64;
    let mut yn: f64 = 1e-3;
    let mut y: f64 = 1.0;

    let t: f64 = 1.0 / t_pr;
    let t_square: f64 = t.powi(2);
    let t_cube: f64 = t.powi(3);
    let a: f64 = 0.06125 * t * (-1.2 * (1.0 - t).powi(2)).exp();

    let big_a: f64 = 14.76 * t - 9.76 * t_square + 4.58 * t_cube;
    let big_b: f64 = 90.7 * t - 242.2 * t_square + 42.4 * t_cube;
    let big_c: f64 = 2.18 + 2.82 * t;

    for _ in 0..1000 {
        if fy.abs() < tolerance {
            break;
        }

        y = yn;
        let y_square: f64 = y.powi(2);
        let y_cube: f64 = y.powi(3);
        let y_quad: f64 = y.powi(4);
        let y_pow_big_c: f64 = y.powf(big_c);

        fy = -a * p_pr + (y + y_square + y_cube - y_quad) / (1.0 - y).powi(3) - big_a * y.powf(2.0)
            + big_b * y_pow_big_c;
        fdy = (1.0 + 4.0 * y + 4.0 * y_square - 4.0 * y_cube + y_quad) / (1.0 - y).powi(4)
            - 2.0 * big_a * y
            + big_b * big_c * y_pow_big_c;
        yn = y - fy / fdy;
    }
    return a * p_pr / y;
}

fn dranchuck_and_aboukassem(t_pr: f64, p_pr: f64, tolerance: f64) -> f64 {
    let mut yn: f64 = 0.5;
    let mut fy: f64 = 1.0;
    let mut fdy: f64;
    let mut y: f64;
    let mut big_r_r: f64;
    // Double Rr = 0.0, C4 = 0.0, A1 = 0.3265, A2 = -1.0700, A3 = -0.5339, A4 = 0.01569, A5 = -0.05165, A6 = 0.5475, A7 = -0.7361, A8 = 0.1844, A9 = 0.1056, A10 = 0.6134, A11 = 0.7210;

    let big_a: [f64; 11] = [
        0.3265, -1.0700, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361, 0.1844, 0.1056, 0.6134,
        0.7210,
    ];

    let t_pr_cube: f64 = t_pr.powi(3);

    let mut big_r_r_square: f64;
    let mut big_a_11_times_big_r_r_square: f64;
    let mut big_r_r_pent: f64;

    let big_c2_big_c3_common: f64 = big_a[6] / t_pr + big_a[7] / t_pr.powi(2);

    let mut big_c: [f64; 4] = [
        big_a[0]
            + big_a[1] / t_pr
            + big_a[2] / t_pr_cube
            + big_a[3] / t_pr.powi(4)
            + big_a[4] / t_pr.powi(5),
        big_a[5] + big_c2_big_c3_common,
        big_a[8] * big_c2_big_c3_common,
        0.0,
    ];

    for _ in 0..1000 {
        if fy.abs() < tolerance {
            break;
        }
        y = yn;
        big_r_r = 0.27 * p_pr / (y * t_pr);
        big_r_r_square = big_r_r.powi(2);
        big_a_11_times_big_r_r_square = big_a[10] * big_r_r_square;

        big_r_r_pent = big_r_r.powi(5);

        big_c[3] = big_a[9]
            * (1.0 + big_a_11_times_big_r_r_square)
            * (big_r_r_square / t_pr_cube)
            * (-big_a_11_times_big_r_r_square).exp();
        fy = y
            - (1.0 + big_c[0] * big_r_r + big_c[1] * big_r_r_square - big_c[2] * big_r_r_pent
                + big_c[3]);
        fdy = 1.0 + big_c[0] * big_r_r / y + 2.0 * big_c[1] * big_r_r_square / y
            - 5.0 * big_c[2] * big_r_r_pent / y
            + 2.0
                * big_a[9]
                * (big_r_r_square / t_pr_cube / y)
                * (1.0 + big_a_11_times_big_r_r_square + (big_a_11_times_big_r_r_square.powi(2)))
                * (-big_a_11_times_big_r_r_square).exp();
        yn = y - fy / fdy;
    }
    return yn;
}

pub fn z_factor(t_pr: f64, p_pr: f64, correlation: ZCorrelations, tolerance: f64) -> f64 {
    let mut z: f64 = 1.0;
    if correlation == ZCorrelations::HallAndYarborough {
        z = hall_and_yarborough(t_pr, p_pr, tolerance);
    } else if correlation == ZCorrelations::DranchkAndAbouKassem {
        z = dranchuck_and_aboukassem(t_pr, p_pr, tolerance);
    }

    return z.max(0.01);
}
