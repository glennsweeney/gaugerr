const DWARF: f64 = f64::MIN_POSITIVE * 10.0;
const GIANT: f64 = f64::MAX / 10.0;
const PI: f64 = std::f64::consts::PI;
const MACHTOL: f64 = f64::EPSILON;

fn ratfun<const M: usize, const N: usize>(x: f64, ak: [f64; M], bk: [f64; N]) -> f64 {
    let mut num = 0.0;
    for i in 0..M {
        num = num * x + ak[i]
    }
    let mut den = 0.0;
    for i in 0..N {
        den = den * x + bk[i]
    }
    return num / den;
}

fn auxln(x: f64) -> f64 {
    const A: [f64; 5] = [
        3.899341537646e-5,
        -8.310525299547e-4,
        -1.423751838241e-1,
        -5.717084236157e-1,
        -4.999999994526e-1,
    ];
    const B: [f64; 4] = [
        1.575899184525e-1,
        9.914744762863e-1,
        1.810083408290e0,
        1.000000000000e0,
    ];
    if x <= -1.0 {
        return -GIANT;
    } else if x < -0.7 || x > 1.36 {
        return (1.0 + x).ln() - x;
    } else if x.abs() < MACHTOL {
        return -0.5 * x.powi(2);
    }
    if x > 0.0 {
        return x.powi(2) * ratfun(x, A, B);
    } else {
        let z = -x / (1.0 + x);
        if z > 1.36 {
            return -((1.0 + z).ln() - z) + x * z;
        } else {
            return -z.powi(2) * ratfun(z, A, B) + x * z;
        }
    }
}

fn gamma_pos(x: f64) -> f64 {
    const A: [f64; 5] = [
        9.308302710346e-3,
        -4.880928874015e-2,
        2.546766167439e-1,
        -3.965937302325e-1,
        1.000000000000e0,
    ];
    const B: [f64; 5] = [
        -5.024949667262e-3,
        9.766752854610e-2,
        -6.508685450017e-1,
        1.510518912977e0,
        -1.345271397926e-1,
    ];
    if x < DWARF {
        return 1.0 / DWARF;
    }
    let k = x.round();
    let m = x.trunc();
    let dw = if k == 0.0 { DWARF } else { 1.0 + x } * MACHTOL;
    if (k - x).abs() < dw && x <= 15.0 {
        let ki = k as i64;
        assert_eq!(k, ki as f64); // Failing here means k isn't representable in i64...
        let mut g = 1.0;
        for j in 2..ki {
            g *= j as f64;
        }
        return g;
    } else if ((x - m) - 0.5).abs() < (1.0 + x) * MACHTOL && x <= 15.0 {
        let mi = m as i64;
        assert_eq!(m, mi as f64); // Failing here means m isn't representable in i64...
        let mut g = PI.sqrt();
        for j in 1..(mi + 1) {
            g *= j as f64 - 0.5;
        }
        return g;
    }
    if x < 1.0 {
        return ratfun(x + 2.0, A, B) / (x * (x + 1.0));
    } else if x < 2.0 {
        return ratfun(x + 1.0, A, B) / x;
    } else if x < 3.0 {
        return ratfun(x, A, B);
    } else if x < 10.0 {
        let mut g = 1.0;
        let mut a = x;
        while a >= 3.0 {
            a -= 1.0;
            g = a * g;
        }
        return g * ratfun(a, A, B);
    } else if x < GIANT.ln() {
        let a = 1.0 / x.powi(2);
        let g = (1.0 + a * (-3.333333333333e-2 + a * 9.52380952381e-3)) / (12.0 * x);
        let a = -x + (x - 0.5) * x.ln() + g + (2.0 * PI).sqrt().ln();
        if a < GIANT.ln() {
            return a.exp();
        } else {
            return GIANT;
        }
    } else {
        return GIANT;
    }
}

pub fn gamma(x: f64) -> f64 {
    if x == 0.0 {
        return std::f64::INFINITY;
    } else if x > 0.0 {
        return gamma_pos(x);
    } else {
        let xp = x.abs();
        let k = xp.round();
        let dw = if k == 0.0 { DWARF } else { 1.0 + xp } * MACHTOL;
        if (k - xp).abs() < dw {
            return std::f64::INFINITY;
        } else {
            return PI / (PI * x).sin() / gamma_pos(1.0 - x);
        }
    }
}

fn gammastar(x: f64) -> f64 {
    const A1: [f64; 3] = [7.806359425652e-2, 9.781658613041e-1, 1.000000000949e0];
    const B1: [f64; 2] = [8.948328926305e-1, 1.000000000000e0];
    const A2: [f64; 4] = [
        9.999999625957e-1,
        9.404953102900e-1,
        4.990196893575e-1,
        5.115471897484e-2,
    ];
    const B2: [f64; 4] = [
        1.000000000000e0,
        8.571609363101e-1,
        4.241288251916e-1,
        1.544892866413e-2,
    ];
    if x > 1.0e10 {
        if x > 1.0 / (12.0 * MACHTOL) {
            return 1.0;
        } else {
            return 1.0 + 1.0 / (12.0 * x);
        }
    } else if x > 12.0 {
        let a = 1.0 / x;
        return ratfun(a, A1, B1);
    } else if x > 1.0 {
        return ratfun(x, A2, B2);
    } else if x > DWARF {
        let a = 1.0 + 1.0 / x;
        return gammastar(x + 1.0) * a.sqrt() * (-1.0 + x * a.ln()).exp();
    } else {
        return 1.0 / ((PI * 2.0).sqrt() * DWARF.sqrt());
    }
}

fn alfa(x: f64) -> f64 {
    if x > 0.25 {
        return x + 0.25;
    } else {
        let lnx = if x <= DWARF { DWARF.ln() } else { x.ln() };
        return -0.6931 / lnx;
    }
}

fn exmin1(x: f64) -> f64 {
    const A: [f64; 4] = [
        1.107965764952e-3,
        2.331217139081e-2,
        6.652950247674e-2,
        9.999999998390e-1,
    ];
    const B: [f64; 4] = [
        -5.003986850699e-3,
        7.338073943202e-2,
        -4.334704979491e-1,
        1.000000000000e0,
    ];
    if x < MACHTOL.ln() {
        return -1.0;
    } else if x > GIANT.ln() {
        return GIANT;
    } else if x < -0.69 || x > 0.41 {
        return x.exp() - 1.0;
    } else if x.abs() < MACHTOL {
        return x;
    } else {
        return x * ratfun(x, A, B);
    }
}

fn auxgam(x: f64) -> f64 {
    const A: [f64; 4] = [
        -6.127046810372e-3,
        4.369287357367e-2,
        -1.087824060619e-1,
        -5.772156647338e-1,
    ];
    const B: [f64; 5] = [
        8.148654046054e-3,
        2.322361333467e-2,
        1.776068284106e-1,
        3.247396119172e-1,
        1.000000000000e0,
    ];
    if x <= -1.0 {
        return -0.5;
    } else if x < 0.0 {
        return -(1.0 + (x + 1.0).powi(2) * ratfun(x + 1.0, A, B)) / (1.0 - x);
    } else if x <= 1.0 {
        return ratfun(x, A, B);
    } else if x <= 2.0 {
        return ((x - 2.0) * ratfun(x - 1.0, A, B) - 1.0) / x.powi(2);
    } else {
        return (1.0 / gamma(x + 1.0) - 1.0) / (x * (x - 1.0));
    }
}

fn qtaylor(a: f64, x: f64, eps: f64) -> f64 {
    let lnx = if x <= DWARF { DWARF.ln() } else { x.ln() };
    let r = a * lnx;

    let mut q = if r < -0.69 || r > 0.41 {
        r.exp() - 1.0
    } else {
        exmin1(r)
    };
    let s = -a * (a - 1.0) * auxgam(a);
    let u = s - q * (1.0 - s);
    let mut p = a * x;
    q = a + 1.0;
    let mut r = a + 3.0;
    let mut t: f64 = 1.0;
    let mut v: f64 = 1.0;
    while (t / v).abs() >= eps {
        p = p + x;
        q = q + r;
        r += 2.0;
        t = -p * t / q;
        v = v + t;
    }
    v = a * (1.0 - s) * ((a + 1.0) * lnx).exp() * v / (a + 1.0);
    return u + v;
}
fn ptaylor(a: f64, x: f64, eps: f64) -> f64 {
    let mut p = 1.0;
    let mut c = 1.0;
    let mut r = a;
    while c / p >= eps {
        r += 1.0;
        c = x * c / r;
        p = p + c;
    }
    return p;
}
fn qfraction() {}
fn pqasymp() {}

pub fn incomplete_gamma(a: f64, x: f64, eps: f64) -> (f64, f64) {
    // Handle special cases of inputs
    if a == 0.0 && x == 0.0 {
        return (0.5, 0.5);
    } else if x == 0.0 {
        return (0.0, 1.0);
    } else if a == 0.0 {
        return (1.0, 0.0);
    }

    // begin dax
    let mu = (x - a) / a;
    let auxlnmu = auxln(mu);
    let dpp = a * auxlnmu - 0.5 * (2.0 * PI * a);
    if dpp < DWARF.ln() {
        // TODO; is this redundant with the later dp>=DWARF test?
        if mu < 0.0 {
            return (0.0, 1.0);
        } else {
            return (1.0, 0.0);
        }
    }
    let dp = dpp.exp() / gammastar(a);
    // end dax

    if dp >= DWARF {
        if a > 25.0 && mu.abs() < 0.2 {
            pqasymp();
            return (0.0, 0.0); // TBD
        } else if a > alfa(x) {
            let p = ptaylor(a, x, eps) * dp;
            let q = 1.0 - p;
            return (p, q);
        } else {
            if x < 1.0 {
                let q = qtaylor(a, x, eps);
                let p = 1.0 - q;
                return (p, q);
            } else {
                qfraction();
                return (0.0, 0.0); // TBD
            }
        }
    } else {
        if a > x {
            return (0.0, 1.0);
        } else {
            return (1.0, 0.0);
        }
    }
}

#[cfg(test)]
mod test_gamma {
    // https://keisan.casio.com/exec/system/1180573444

    use approx;

    #[test]
    fn special_cases() {
        assert_eq!(super::gamma(0.0), std::f64::INFINITY);
        assert_eq!(super::gamma(std::f64::MIN_POSITIVE), 1.0 / super::DWARF);
        assert_eq!(super::gamma(super::GIANT.ln()), super::GIANT);
        assert_eq!(super::gamma(super::GIANT), super::GIANT);
        assert_eq!(super::gamma(std::f64::MAX), super::GIANT);
    }

    #[test]
    fn negative_integers() {
        assert_eq!(super::gamma(-1.0), std::f64::INFINITY);
        assert_eq!(super::gamma(-2.0), std::f64::INFINITY);
        assert_eq!(super::gamma(-10000.0), std::f64::INFINITY);
    }

    #[test]
    fn negatives() {
        approx::assert_relative_eq!(super::gamma(-0.001), -1000.57820562, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(-0.01), -100.5871979644, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(-0.1), -10.686287021193, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(-0.5), -3.5449077018110, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(-0.9), -10.570564109632, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(-0.99), -100.4369546658, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(-0.999), -1000.42419668, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(-5.5), 0.01091265478191, max_relative = 1.0e-10);
    }

    #[test]
    fn positive_integers() {
        assert_eq!(super::gamma(1.0), 1.0);
        assert_eq!(super::gamma(2.0), 1.0);
        assert_eq!(super::gamma(3.0), 2.0);
        assert_eq!(super::gamma(4.0), 6.0);
        assert_eq!(super::gamma(5.0), 24.0);
        assert_eq!(super::gamma(10.0), 362880.0);
        assert_eq!(super::gamma(14.0), 6227020800.0);
    }

    #[test]
    fn positive_half_integers() {
        approx::assert_relative_eq!(super::gamma(0.5), 1.77245385090551, max_relative = 1.0e-12);
        approx::assert_relative_eq!(super::gamma(1.5), 0.88622692545275, max_relative = 1.0e-12);
        approx::assert_relative_eq!(super::gamma(2.5), 1.32934038817913, max_relative = 1.0e-12);
        approx::assert_relative_eq!(super::gamma(3.5), 3.32335097044784, max_relative = 1.0e-12);
        approx::assert_relative_eq!(super::gamma(4.5), 11.6317283965674, max_relative = 1.0e-12);
        approx::assert_relative_eq!(super::gamma(5.5), 52.3427777845535, max_relative = 1.0e-12);
        approx::assert_relative_eq!(super::gamma(10.5), 1133278.3889487, max_relative = 1.0e-12);
        approx::assert_relative_eq!(super::gamma(14.5), 23092317922.314, max_relative = 1.0e-12);
    }

    #[test]
    fn smallest_reals() {
        approx::assert_relative_eq!(super::gamma(0.001), 999.42377248459, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(0.01), 99.4325851191506, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(0.25), 3.62560990822191, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(0.33), 2.70720622261519, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(0.66), 1.36616419875147, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(0.75), 1.22541670246517, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(0.99), 1.00587197964411, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(0.999), 1.0005782056293, max_relative = 1.0e-10);
    }

    #[test]
    fn smaller_reals() {
        approx::assert_relative_eq!(super::gamma(1.001), 0.9994237724846, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(1.01), 0.99432585119151, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(1.25), 0.90640247705548, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(1.33), 0.89337805346301, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(1.66), 0.90166837117597, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(1.75), 0.91906252684888, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(1.99), 0.99581325984767, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(1.999), 0.9995776274237, max_relative = 1.0e-10);
    }

    #[test]
    fn small_reals() {
        approx::assert_relative_eq!(super::gamma(2.001), 1.0004231962571, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(2.01), 1.00426910970342, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(2.25), 1.13300309631934, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(2.33), 1.18819281110580, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(2.66), 1.49676949615211, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(2.75), 1.60835942198554, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(2.99), 1.98166838709685, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(2.999), 1.9981556772200, max_relative = 1.0e-10);
    }

    #[test]
    fn reals() {
        approx::assert_relative_eq!(super::gamma(3.001), 2.0018468157104, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(4.01), 6.07592854061667, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(5.25), 35.2116118527997, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(6.33), 212.765976208088, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(7.66), 2559.73037972556, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(8.75), 23698.1257017427, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(9.99), 354802.017019831, max_relative = 1.0e-10);
    }

    #[test]
    fn big_reals() {
        approx::assert_relative_eq!(super::gamma(16.0), 1.30767436800e12, max_relative = 1.0e-10);
        approx::assert_relative_eq!(super::gamma(22.5), 2.38280159446e20, max_relative = 1.0e-10);
        approx::assert_relative_eq!(
            super::gamma(104.56),
            1.3328667189E+165,
            max_relative = 1.0e-10
        );
        approx::assert_relative_eq!(
            super::gamma(133.33),
            5.6111035648E+224,
            max_relative = 1.0e-10
        );
        approx::assert_relative_eq!(super::gamma(707.0), super::GIANT, max_relative = 1.0e-10);
    }
}

#[cfg(test)]
mod test_helpers {

    use approx;

    #[test]
    fn test_auxln() {
        assert_eq!(super::auxln(-1.1), -super::GIANT);
        approx::assert_relative_eq!(super::auxln(1.0), -0.306852819440, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::auxln(-0.7), -0.50397280433, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::auxln(0.0), 0.0000000000000, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::auxln(-2.0e-20), -2.000e-40, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::auxln(2.0e-20), -2.0000e-40, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::auxln(6.4), -4.398519999790, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::auxln(-0.8), -0.80943791243, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::auxln(-0.2), -0.02314355131421, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::auxln(-0.1), -0.00536051565782, max_relative = 1.0e-9);
    }

    #[test]
    fn test_exmin1() {
        assert_eq!(super::exmin1(-40.0), -1.0);
        assert_eq!(super::exmin1(1.234e-30), 1.234e-30);
        assert_eq!(super::exmin1(0.0), 0.0);
        assert_eq!(super::exmin1(800.0), super::GIANT);
        approx::assert_relative_eq!(super::exmin1(-0.8), -0.5506710358828, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::exmin1(-0.7), -0.5034146962086, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::exmin1(2.4), 10.02317638064160, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::exmin1(-0.1), -0.0951625819640, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::exmin1(-0.01), -0.009950166251, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::exmin1(0.01), 0.01005016708417, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::exmin1(0.1), 0.105170918075648, max_relative = 1.0e-9);
    }

    #[test]
    fn test_gammastar() {
        assert_eq!(super::gammastar(1e20), 1.0);
        approx::assert_relative_eq!(super::gammastar(1e11), 1.00000000000, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::gammastar(20.0), 1.00417501087, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::gammastar(8.0), 1.010465651062, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::gammastar(6.0), 1.013972849149, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::gammastar(2.5), 1.033718890965, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::gammastar(0.6), 1.139258813362, max_relative = 1.0e-9);
        approx::assert_relative_eq!(super::gammastar(0.1), 1.669860485516, max_relative = 1.0e-9);
        approx::assert_relative_eq!(
            super::gammastar(1e-300),
            3.989422804e149,
            max_relative = 1.0e-9
        );
        approx::assert_relative_eq!(
            super::gammastar(1e-307),
            8.457419059e152,
            max_relative = 1.0e-9
        );
    }
}
