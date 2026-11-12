//! Unit tests: the Python- and numpy-compatible pieces against values taken from Python 3.15 / numpy 2.5
//! (the whole tool is held to locusmasker.py by validation/compare_ports.py and validation/fuzz_ports.py).
use crate::npmath::*;
use crate::pyfmt::*;

#[test]
fn float_repr_is_pythons() {
    let cases = [(0.1, "0.1"), (1e16, "1e+16"), (1e15, "1000000000000000.0"), (1e-5, "1e-05"), (0.0001, "0.0001"), (123456789.123, "123456789.123"), (-0.0, "-0.0"),
                 (5e-324, "5e-324"), (1.7976931348623157e308, "1.7976931348623157e+308"), (0.30000000000000004, "0.30000000000000004"), (100.0, "100.0"), (2.5e-7, "2.5e-07")];
    for (x, want) in cases { assert_eq!(repr_f64(x), want); }
}

#[test]
fn round_and_formats_are_pythons() {
    assert_eq!(round(2.675, 2), 2.67); assert_eq!(round(1.0000005, 6), 1.000001); assert_eq!(round(123.456789, 4), 123.4568);
    assert!(round(-0.00049, 3) == 0.0 && round(-0.00049, 3).is_sign_negative());
    assert_eq!(fnum(-0.00049, 3), "0.000"); assert_eq!(fnum(f64::NAN, 2), ""); assert_eq!(pfmt(f64::NEG_INFINITY), ""); assert_eq!(pfmt(0.5), "1.000000");
    for (x, want) in [(0.05, "0.05"), (1e-6, "1e-06"), (123456789.0, "1.23457e+08"), (0.0001, "0.0001"), (1.5, "1.5"), (100000.0, "100000"), (1e6, "1e+06")] { assert_eq!(g(x), want); }
    assert_eq!(commas(1234567i64), "1,234,567"); assert_eq!(commas(-1000i64), "-1,000"); assert_eq!(commas(999i64), "999");
    assert_eq!(repr_str("it's"), "\"it's\""); assert_eq!(repr_str("a\tb"), "'a\\tb'");
}

#[test]
fn json_and_csv_are_pythons() {
    let j = crate::jdict! { "a" => J::List(vec![J::Int(1), J::Float(2.5), J::Null, J::Bool(true), J::s("x\"y\u{e9}\n")]), "b" => J::Dict(vec![]), "c" => J::List(vec![]), "d" => J::Float(f64::NAN) };
    assert_eq!(j.indent1().replace('\n', "|"), "{| \"a\": [|  1,|  2.5,|  null,|  true,|  \"x\\\"y\\u00e9\\n\"| ],| \"b\": {},| \"c\": [],| \"d\": NaN|}");
    assert_eq!(j.compact(), "{\"a\":[1,2.5,null,true,\"x\\\"y\\u00e9\\n\"],\"b\":{},\"c\":[],\"d\":NaN}");
    assert!(json_parse(&j.compact()).is_some());
    let row: Vec<String> = ["1", "a,b", "q\"r", "x y", "", " lead", "l\nf"].iter().map(|s| s.to_string()).collect();
    assert_eq!(csv_row(&row), "1,\"a,b\",\"q\"\"r\",x y,, lead,\"l\nf\"\r\n");
}

#[test]
fn npy_header_is_numpys() {
    let b = crate::npy::npy(&[3, 4], &[0f32; 12]);
    assert_eq!(b.len(), 176); assert_eq!(&b[..10], b"\x93NUMPY\x01\x00v\x00"); assert!(b[10..].starts_with(b"{'descr': '<f4', 'fortran_order': False, 'shape': (3, 4), } "));
    assert!(crate::npy::npy1(&[1i64, 2]).windows(13).any(|w| w == b"'shape': (2,)")); assert_eq!(crate::npy::read_npy(&b).unwrap().1, vec![3, 4]);
}

#[test]
fn choice_is_numpys() {
    let first6 = |pop, size| crate::rng::choice_sorted(0, pop, size)[..6].to_vec();
    assert_eq!(first6(50, 7), [2, 3, 12, 14, 23, 28]); assert_eq!(first6(20000, 500), [54, 106, 125, 158, 165, 196]);          // Floyd's algorithm; the tail shuffle
    assert_eq!(first6(20000, 100), [54, 106, 165, 329, 440, 564]); assert_eq!(first6(2_000_000, 30), [5476, 33054, 67171, 81946, 150478, 350530]);
    assert_eq!(crate::rng::choice_sorted(0, 9, 9), (0..9).collect::<Vec<_>>());
}

#[test]
fn hashes_and_draws_are_the_tools() {
    assert_eq!([sm64(0), sm64(1), sm64(12345)], [16294208416658607535, 10451216379200822465, 2454886589211414944]);
    use crate::prepare::{bs_keys_cell, bs_keys_taxa};
    assert_eq!(hash_draw(bs_keys_taxa(0), 30, 4), [28, 11, 29, 24]); assert_eq!(hash_draw(bs_keys_taxa(1), 30, 4), [17, 27, 28, 22]);
    assert_eq!(hash_draw(bs_keys_cell(3, 2), 500, 50)[..5], [483, 342, 434, 237, 350]); assert_eq!(hash_draw(bs_keys_cell(9, 2), 500, 50)[..5], [216, 334, 189, 431, 340]);
}

#[test]
fn special_functions_match_to_the_bit() {
    let lg = [(0.3, 1.0957979948180763), (1.0, 0.0), (7.5, 7.534364236758733), (123.25, 468.6144829505166), (4000.0, 29172.970494954596)];
    for (x, want) in lg { assert_eq!(lgamma_v(x).to_bits(), f64::to_bits(want), "lgamma_v({x})"); }
    let lf = log_fact(1002);
    let (o, n, e, phi) = ([30., 2., 0., 55., 10., 3.], [200., 300., 50., 60., 1000., 40.], [10., 12., 3., 20., 25.5, 0.0], [0.0, 0.02, 0.5, 0.001, 1e-12, 0.1]);
    let up = [-16.256292659923375, -0.00023410280570779202, 0.0, -31.97441295682957, -0.00013788904382700174, -690.7755278982137];
    let dn = [-2.478249783166575e-08, -6.6671309609495415, -1.7868390023785707, 0.0, -7.881158662602448, 0.0];
    for i in 0..6 {
        assert_eq!(bb_tail(o[i], n[i], e[i], phi[i], true, &lf).to_bits(), f64::to_bits(up[i]), "upper tail {i}");
        assert_eq!(bb_tail(o[i], n[i], e[i], phi[i], false, &lf).to_bits(), f64::to_bits(dn[i]), "lower tail {i}");
    }
    assert_eq!(norm_quantile(1.0 - 0.05 / 1234.0), 3.9413115842417996);
}

#[test]
fn numpy_reductions() {
    let v = [3.0, 1.0, 4.0, 1.5, 9.0, 2.6, 5.3];
    assert_eq!(percentile(&v, 10.0), 1.3); assert_eq!(percentile(&v, 90.0), 6.780000000000001); assert_eq!(interp(2.2, &[0.0, 1.0, 3.0], &[0.0, 10.0, 11.0]), 10.6);
    assert_eq!(median(&v), 3.0); assert_eq!(median(&v[..6]), 2.8); assert!(median(&[1.0, f64::NAN]).is_nan()); assert_eq!(nanmedian(&[1.0, f64::NAN, 3.0]), 2.0);
    let a: Vec<f64> = (0..300).map(|i| 0.1 * i as f64 + 1e-9 * (i * i) as f64).collect();
    assert_ne!(pairwise_sum(&a).to_bits(), a.iter().sum::<f64>().to_bits());          // numpy's pairwise order is not the running sum
    let mut vals = vec![5.0, 1.0, 3.0, 2.0, 4.0]; let mut w = vec![1.0, 1.0, 1.0, 1.0, 10.0];
    assert_eq!(wmed_select(&mut vals, &mut w, 5), 4.0); assert_eq!(wmedian1(&[5.0, 1.0, 3.0], &[1.0, 1.0, 1.0]), 3.0);
    assert!(npmax(f64::NAN, 1.0).is_nan() && pymax(f64::NAN, 1.0).is_nan() && pymax(1.0, f64::NAN) == 1.0);
}

#[test]
fn nexus_helpers_write_as_biopython() {
    use crate::input::{compact4nexus, safename};
    assert_eq!([compact4nexus(&[0, 1, 2, 3]), compact4nexus(&[0, 2, 4, 6, 9]), compact4nexus(&[5]), compact4nexus(&[0, 1, 2, 10, 11, 20]), compact4nexus(&[])], ["1-4", "1-7\\2 10", "6", "1-3 11-12 21", ""]);
    assert_eq!([safename("abc"), safename("a b"), safename("it's"), safename("x-y")], ["abc", "'a b'", "'it''s'", "'x-y'"]);
    let r = crate::report::ranges;
    assert_eq!(r(&[0, 1, 2, 3]), [[0, 3, 1]]); assert_eq!(r(&[0, 2, 4, 6, 9]), [[0, 6, 2], [9, 9, 1]]); assert_eq!(r(&[5]), [[5, 5, 1]]);
    assert_eq!(r(&[0, 3, 4, 5]), [[0, 0, 1], [3, 5, 1]]); assert_eq!(r(&[1, 4]), [[1, 1, 1], [4, 4, 1]]);
}

#[test]
fn readers() {
    let (names, seqs) = crate::input::read_fasta(b">a desc\nAC\nGT\n>b\r\nAC GT\r\n");
    assert_eq!(names, ["a", "b"]); assert_eq!(seqs, [b"ACGT".to_vec(), b"ACGT".to_vec()]);
    let (names, seqs) = crate::input::read_phylip("2 4\na ACGT\nb AC\nGT\n"); assert_eq!(names, ["a", "b"]); assert_eq!(seqs, ["ACGT", "ACGT"]);
    let (names, seqs) = crate::input::read_phylip("2 8\na ACGT\nb ACGA\n\nTTTT\nCCCC\n"); assert_eq!(names, ["a", "b"]); assert_eq!(seqs, ["ACGTTTTT", "ACGACCCC"]);
}
