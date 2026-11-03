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
    reading: "reading…",                    // under the full bar, while the loaded data is decoded (big reports: a few seconds)

    // the same screen, bar held at {x} / {y} MB, when the file on disk is short (copied or downloaded before it
    // was whole) or its data is damaged. The page retries by itself ({s} counts down); no action needed.
    incomplete: "The file is incomplete — retrying in {s}s.",

    // the same screen when the page's own scripts never ran (an IDE pane or preview that blocks or truncates
    // scripts). Retrying cannot help there, so say where to open the file instead.
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

    // one card per test: its name, its mask count, and ONE line on what the test asks. Click: the table, filtered.
    card_close_relative: "close-relative test",
    card_close_relative_line: "too divergent from its close relatives",
    card_closest_match: "closest-match test",
    card_closest_match_line: "too similar to another sequence at the locus",
    card_broad_sample: "broad-sample test",
    card_broad_sample_line: "out of line with a broad sample of taxa",
    card_divergent_labels: "{n} divergent cells labelled also-too-close",   // the broad card's second line, when any

    worst_taxa: "most masked taxa (share of tested loci)",
    worst_loci: "most masked loci (sequences)",
    hist_taxa: "taxa by their masked sequences",           // the distribution plots under the lists
    hist_loci: "loci by their masked sequences",
    hist_x: "masked sequences in the {unit}",              // {unit} = taxon | locus
    hist_sig: "masks by significance",                     // the significance histogram
    hist_sig_x: "mask log10 p",
    bar_note: "the run's bar",                             // the dashed line's label on the significance histogram
  },

  // ---- the overview table: one row per cell. Headers; the family filter chips; the footer.
  table: {
    h_taxon: "taxon", h_locus: "locus",
    h_test: "test",                                        // which test masked it (the strongest, as masked.csv's mask_test)
    h_p: "log10 p",                                        // the mask's one significance, with its stars
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

  // ---- fatal states: shown in #app in place of the report.
  fatal: {
    no_data: "This page carries no data: it is the viewer template without a report. Run locus_outliers.py; the report.html it writes has the data inside it.",
    unreadable: "The report could not be read.",   // followed by the error detail, for a bug report
    // a report written by a NEWER tool than this viewer; {found} and {max} are format numbers, {version} the tool's
    newer_format: "this report was written in bundle format {found}; this viewer reads formats 1 to {max}. Please use the viewer that came with the tool version {version}.",
  },

};
// the viewer and the loader fill {placeholders} with this
function lmFmt(s, vars) { return s.replace(/\{(\w+)\}/g, (_, k) => vars[k] !== undefined ? vars[k] : "{" + k + "}"); }
