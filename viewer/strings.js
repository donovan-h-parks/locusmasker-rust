/* Every human-visible word in the viewer lives in this file, so it can be reviewed and edited in one place.
 *
 * How to read an entry:
 *   - the section comment says WHERE the text appears and WHEN it shows;
 *   - {name} is a placeholder the viewer fills in (the comment names each one);
 *   - edit the words freely: the tests check structure (classes and data-state attributes), never these words.
 *
 * Two bootstrap strings cannot live here, because they must exist before any script runs:
 *   - viewer.html, the static line inside <div id="app">: shown only when the browser runs no scripts at all;
 *   - viewer.html, <title>: the tab title before the report loads (the loader retitles the tab from tab_title).
 *
 * The strings below are being gathered view by view (the loading screen is done); until the sweep finishes,
 * older views still carry their words inline in viewer.js.
 */
const LM_STR = {

  // ---- the loading screen: shown from the first moment until the report is read, in #app.
  // Big logo, the dataset's name, a progress bar, and the live figure under it. No other words while all is well.
  loading: {
    tab_title: "{dataset} · LocusMasker",   // the browser tab, as soon as the dataset's name is known
    no_title: "LocusMasker report",         // in place of the dataset name when the page carries none (template without data)
    progress: "{x} / {y} {u}",              // under the bar, ticking live while the file arrives; x = loaded, y = total, u = MB or KB
    line: "Loading data into browser…",     // the one line under the figures, for the whole journey (a short file keeps
                                            // loading and the page re-reads it by itself; no other words)

    // only when the page's own scripts never ran (an IDE pane or preview that blocks or truncates scripts):
    // reloading cannot help there, so say where to open the file instead.
    stalled: "This page's scripts have not run here. Open report.html in a regular browser (Chrome or Firefox).",
  },

  // ---- the overview's summary (the landing view): the hero figures, then one card per test.
  // The cards run in the tests' order: close-relative, closest-match, broad-sample. {n} is each card's mask count.
  overview: {
    hero_sub: "of {tested} sequences tested ({pct})",      // under the big masked count
    hero_taxa: "{with} / {all} taxa hold a mask",          // the small stats beside the hero
    hero_loci: "{with} / {all} loci hold a mask",
    hero_excluded: "{n} taxa excluded",                    // only when the run excluded taxa
    hero_skipped: "{n} loci skipped",                      // only when the run skipped loci

    // one card per test ARM (four cards): its name, its mask count, what it asks, and when it is used.
    // The words "divergent" and "similar" are tinted the design's red and blue wherever they appear on a card.
    // These cards show on the summary (click: the table, filtered) and again on the dataset card, below the ribbons.
    card_close_relative: "close-relative test",
    card_close_relative_line: "too divergent from its close relatives",
    card_cm_rel: "closest-match test",
    card_cm_rel_line: "too similar to its closest match at this locus",
    card_cm_rel_when: "used when the closest match is one of the close relatives",
    card_cm_out: "closest-match test",
    card_cm_out_line: "too similar to its closest match at this locus",
    card_cm_out_when: "used when the closest match is not one of the close relatives",
    card_broad_sample: "broad-sample test",
    card_broad_close: "{n} sequences too similar to one or more of the broad sample of {k} taxa",
    card_broad_divergent: "{n} sequences too divergent to one or more of the broad sample of {k} taxa",

    worst_taxa: "most masked taxa (share of tested loci)",
    worst_loci: "most masked loci (sequences)",
    hist_taxa: "taxa by their masked sequences",           // the distribution plots under the lists
    hist_loci: "loci by their masked sequences",
    hist_x: "masked sequences in the {unit}",              // {unit} = taxon | locus
  },

  // ---- the matrix view: the zoom controls' tooltips and the one-time hint chip (gone on first zoom)
  matrix: {
    zoom_in: "zoom in", zoom_out: "zoom out", zoom_fit: "fit",
    hint: "scroll to move · ctrl+scroll or pinch to zoom · click a cell",
  },

  // ---- the residual ribbons (overview > loci/taxon ribbons): one row per unit, one mark per sequence, placed by
  // its residual from the SELECTED test's own expectation. The axis caption names each test's units.
  ribbons: {
    test_label: "test",
    test_close_relative: "close-relative",
    test_cm_relative: "closest-match · relative",   // the arm whose partner is one of the focal's relatives
    test_cm_outside: "closest-match · outside",     // the arm whose partner is not
    test_broad_sample: "broad-sample",
    axis_close_relative: "◀ less divergent · excess divergence from its close relatives, in units of m · more divergent ▶",
    axis_closest: "similarity to its closest match, p · stronger evidence ▶",
    axis_broad: "◀ steadier than expected · observed / expected variance over the broad sample, log2 · wilder ▶",
  },

  // ---- the overview table: one row per cell. Headers; the family filter chips; the footer.
  table: {
    h_taxon: "taxon", h_locus: "locus",
    h_test: "test",                                        // which test masked it (the strongest, as masked.csv's mask_test)
    h_p: "p",                                              // the mask's one significance, with its stars
    h_obs: "observed %", h_exp: "expected %",              // divergence to the close relatives, median
    h_excess: "excess (m)", h_partner: "partner", h_best: "best match",
    chips: { all: "all", close_relative: "close-relative", closest_match: "closest-match", broad_sample: "broad-sample" },
    masked_only: "masked cells only",
    foot: "{n} rows · click a row for the alignment · click a header to sort, shift-click for the tie-breaking key",
  },

  // ---- the matrix legend (cell colours) and its names for the three mask kinds
  legend: {
    not_masked: "not masked",
    divergent: "too divergent",
    similar: "too similar",
    both: "both directions",
    untested: "not tested",
    absent: "absent",
  },

  // ---- the alignment (sequence) panel: the gutter names the tests, in the order closest-match, broad-sample,
  // close-relative; the broad row is always present, showing the sample's biggest-residual sequence with the
  // variance test's p.
  seq: {
    gut_closest: "closest-match",
    gut_broad: "broad-sample",
    gut_relatives: "close-relative",
    broad_none: "none drawn",                        // the broad row's name when the report carries no representative
    // the "sequence" column: what each row is to the focal
    role_focal: "focal",
    role_relative: "relative",                       // the closest match, when it is one of the close relatives
    role_nonrelative: "non-relative",                // ... when it is not
    role_largest: "largest residual",                // the broad-sample row
    role_relative_n: "relative {n}",                 // the close relatives, in genome-wide rank order
    // the two trees under the alignment: their titles
    tree_left: "distances across the alignment",     // the expected tree: genome-wide divergences, scaled to this locus
    tree_right: "distances at this locus",           // the observed tree: this locus's own divergences
  },

  // ---- the navigation stack: the deepest open level fills the window; the levels above it collapse to one
  // band each (click: back to it). These are the only words the stack adds.
  nav: {
    up: "back up · Esc",                             // the bands' and the ✕ button's tooltip
    copies: "copies",                                // the unit band's name when the open unit is the copies panel
  },

  // ---- the copy flow (the overview's summary and the dataset card): donors on the left, recipients on the
  // right, one band per flow, its thickness the attributed copies; {n} in the title is the summary's own picker.
  flow: {
    title: "Top {n} donor and recipient taxa, plus others",
    sub: "Number of events per taxon and percentage of all sequences for each taxon (or all others) shown. {n} events from others to others not shown.",
    others: "others",
    hover: "{don} → {rec} · {n}",                    // the band's hover; {n} = its copies
    hover_click: "click: the copies' site ribbons",
  },

  // ---- the cards: one printable page per dataset / locus / taxon / sequence, one design language.
  // The verdict band is the "fine / bad / discuss" line; the spine is the three per-test ribbons on dataset-wide
  // axes, so cards compare side by side; the table names what is involved; the footer makes the page traceable.
  card: {
    button: "card",                                  // on each panel: opens the printable card
    print: "print", close: "close",                  // the preview's two buttons
    type_dataset: "dataset card", type_locus: "locus card", type_taxon: "taxon card", type_sequence: "sequence card",
    verdict_clean: "No problems found",              // the calm verdict when nothing is masked
    verdict_unit: "{n} of {m} sequences masked ({pct})",   // locus and taxon cards
    verdict_dataset: "{n} of {m} sequences masked ({pct})",
    t_taxon: "taxon", t_locus: "locus", t_test: "test", t_p: "p", t_obs: "observed %", t_exp: "expected %", t_partner: "partner",
    rank_taxa: "% masked per taxon",                 // the dataset card's two rank curves (little titles above each)
    rank_loci: "% masked per locus",
    rank_x: "ranked",                                // the curves' x caption
    // the unit cards' body: per test, the top masked sequences as site ribbons, each with its p; the same
    // sequences are numbered on the spine's ribbon plots above
    top_ribbons: "top {n} · {test}",
    align_title: "alignment",                        // the sequence card's alignment: the sequence tab's rows (focal,
                                                     // closest, broad, relatives), whole locus, sized to the page
    trees_title: "expected · observed",              // the sequence card's face-to-face trees
    // the sequence card's results block: one row per test, no interpretation. {p} is the close-relative row's
    // median over the comparator tests; {call} is empty or one of the two call strings below it.
    r_cr: "{nf} of {nc} failed · observed {obs}% · expected {exp}% · p {p}",
    r_cm: "closest match {name} · {arm} · p {p}",
    r_cm_arm_rel: "a relative", r_cm_arm_out: "outside the relatives",
    r_bs: "variance {v}× expected · p {p}{call}",
    r_bs_divergent: " · sequence too divergent",
    r_bs_close: " · sequence too close",
    footer: "{dataset} · LocusMasker {version} · {date} · every test held to p < {bar}",
  },

  // ---- fatal states: shown in #app in place of the report.
  fatal: {
    no_data: "This page carries no data: it is the viewer template without a report. Run locus_outliers.py; the report.html it writes has the data inside it.",
    unreadable: "The report could not be read.",   // followed by the error detail, for a bug report
    // a report written by a NEWER tool than this viewer; {found} and {max} are format numbers, {version} the tool's
    newer_format: "this report was written in bundle format {found}; this viewer reads formats 1 to {max}. Please use the viewer that came with the tool version {version}.",
  },

};
// p-values for people: plain p, two decimals down to 0.01, scientific below (never log10 in a display).
// The inputs everywhere are log10 p (the lossless machine form); this is the one place they become words.
function lmP(lp) { if (lp === null || lp === undefined || Number.isNaN(lp)) return ""; const p = Math.pow(10, lp); return p >= 0.01 ? p.toFixed(2) : p.toExponential(1); }
// the viewer and the loader fill {placeholders} with this
function lmFmt(s, vars) { return s.replace(/\{(\w+)\}/g, (_, k) => vars[k] !== undefined ? vars[k] : "{" + k + "}"); }
