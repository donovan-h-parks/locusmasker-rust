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
 * Only the loader's own error details (a bad ZIP member, an unsupported dtype) stay in viewer.js: they are shown
 * after fatal.unreadable for a bug report, never as prose.
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
    hero_taxa: "{with} / {all} taxa hold a mask",          // the small stats beside the hero
    hero_loci: "{with} / {all} loci hold a mask",
    hero_excluded: "{n} taxa excluded",                    // only when the run excluded taxa
    hero_skipped: "{n} loci skipped",                      // only when the run skipped loci

    // one card per test ARM (four cards): its name, its mask count, what it asks, and when it is used.
    // The words "divergent" and "similar" are tinted the design's red and blue wherever they appear on a card.
    // These cards show on the summary (click: the table, filtered) and again on the dataset card, below the ribbons.
    card_close_relative: "close-relative test",
    card_close_relative_line: "too divergent from its {k} close relatives",   // {k}: the comparators each sequence gets
    card_cm_rel: "closest-match · relative",
    card_cm_rel_line: "too similar to its closest match at this locus",
    card_cm_rel_when: "used when the closest match is one of the close relatives",
    card_cm_out: "closest-match · outside",
    card_cm_out_line: "too similar to its closest match at this locus",
    card_cm_out_when: "used when the closest match is not one of the close relatives",
    card_broad_sample: "broad-sample test",
    card_broad_close: "{n} sequences too similar to part of the {k} broad samples",
    card_broad_divergent: "{n} sequences too divergent to part of the {k} broad samples",

    worst_taxa: "most masked taxa (share of tested loci)",
    worst_loci: "most masked loci (share of tested taxa)",
    worst_label: "{pct} ({n} mask{s})",                    // after each name in the two lists
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
    axis_broad: "evidence from the broad sample, p (the stronger side) · stronger evidence ▶",
  },

  // ---- the masks bar, under the topbar at every level: the one filter and order behind the table and the ← → paging.
  filter: {
    label: "filter masks",                         // the row's name; a click on it folds the row away
    fold: "fold the filter row away", unfold: "show the filter row",
    all: "all tests",                                // the chip that widens back to every test
    arm_close_relative: "close-relative",            // the five arms, one chip each (a click keeps that arm; more clicks add or drop arms)
    arm_cm_relative: "closest-match · relative",
    arm_cm_outside: "closest-match · outside",
    arm_bs_close: "broad-sample · too close",
    arm_bs_divergent: "broad-sample · too divergent",
    any_taxon: "any taxon", any_locus: "any locus",  // the name boxes' placeholders (type to search, pick to narrow)
    clear: "any",                                    // the ✕ beside a chosen name
    order: "order", order_p: "most significant first", order_locus: "by locus", order_taxon: "by taxon",
    count: "{n} masks",                              // how many the filter leaves
    position: "{k} of {n}", none: "– of {n}",        // where the open sequence sits in the list; a dash when it is not in it
    prev: "previous mask (←)", next: "next mask (→)",
    this_unit: "page its masks",                     // on a taxon or locus panel: narrow the list to this unit
  },

  // ---- the overview table: one row per mask of the bar's list. Headers; the footer.
  table: {
    h_id: "id", h_taxon: "taxon", h_locus: "locus",      // id: masked.csv's row number (rows by locus, then taxon)
    h_reason: "reason",                                    // divergent, similar or both, as masked.csv's reason
    h_test: "test",                                        // which test masked it (the strongest, as masked.csv's mask_test)
    h_p: "p",                                              // the mask's one significance, with its stars
    h_cr: "close-relative p", h_cm: "closest-match p", h_bs: "broad-sample p",   // each test's own p for the cell
    h_partner: "partner",
    r_divergent: "divergent", r_similar: "similar", r_both: "both",
    source: "a summary of masked.csv, which carries every column",
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
    // the alignment's key (one rule for every row; the focal is the reference): focal · red · blue
    key_focal: "focal: the reference",
    key_red: "differs from the focal",
    key_blue: "agrees with the focal where most relatives differ",
    // a tree tip's hover: its name, then what its colour means (the row's p<α colour)
    tip_focal: "focal",
    tip_red: "significantly too divergent from the focal",
    tip_blue: "significantly too similar to the focal",
    tip_plain: "neither test fired",
    // the ruler button over the alignment: what the numbers count. It shows the current choice; click to switch.
    ruler_alignment: "alignment positions",          // the site's column in the whole alignment (the default)
    ruler_locus: "locus positions",                  // the site's number from the locus's start
    ruler_hint: "the ruler: positions in the whole alignment, or from the locus's first site · click to switch",
    // the stage's other words: loading and empty states, the rows' hovers, the verdict line, the controls, the footer
    loading: "loading the sequences…", loading_trees: "loading the trees…", no_trees: "trees were not carried for this cell",
    no_seqs: "sequences were not carried for this locus", nothing: "neither sequences nor site codes were carried for this cell in this report",
    row_no_closest: "no closest match at this locus", row_no_closest_old: "closest match not recorded in this report", row_missing: "row not carried in this report",
    name_none: "none", name_not_recorded: "not recorded",            // the closest row's name when there is none
    // each row's hover: what it is to the focal; {ord} is one of ord below ({n}th beyond the fifth)
    ord: ["closest", "2nd closest", "3rd closest", "4th closest", "5th closest"], ord_n: "{n}th closest",
    tip_focal_row: "the focal sequence",
    tip_no_closest: "no taxon shares enough sites with the focal at this locus",
    tip_no_closest_old: "this report does not record the closest match of a sequence whose closest-match test did not run (it is carried from report format 7)",
    tip_closest: "the focal's closest match at this locus{rel}", tip_closest_rel: " (also relative {n})", tip_closest_nonrel: " (not one of its relatives)",
    tip_no_broad: "the report carries no broad-sample representative for this cell", tip_broad: "the broad comparator with the most extreme standardised residual",
    tip_relative: "relative {n}: the {ord} taxon genome-wide",
    tip_differs: "; differs from the focal at {pct}% of {n} shared sites",
    // the verdict line over the alignment
    ctx_taxon_note: "(points are other loci for this taxon)", ctx_locus_note: "(points are other taxa for this locus)",   // after the two context columns' names
    v_not_tested: "Not tested: {reason}", v_no_comparators: "no comparators",
    v_masked: "Masked", v_not_masked: "Not masked", v_off: "Not masked (masking was off)", v_fewer: " (fewer comparators than requested)",
    says_divergent: "significantly too divergent", says_bs_divergent: "significantly too divergent from {name}, a broad-sample comparator",
    donor_uncertain: " (donor uncertain: both masked)", donor_likely: " ({name} looks more likely to be the donor)",
    says_bs_close: "significantly too close to {name}, a broad-sample comparator", says_similar: "significantly too similar to {name}",
    says_bs_also: "the broad sample also calls it too close to {name}'s side",
    // the tests in words, in the rows' hovers
    t_not_tested: "{what}: not tested", t_p: "{what}: p {p}{stars}", t_bar: " · bar {bar}",
    t_divergent: "too divergent", t_similar_to: "too similar to {name}", t_similar_focal: "too similar to the focal", t_divergent_focal: "too divergent from the focal",
    t_bs_close_focal: "too close to the focal (the broad-sample test, stage 1)", t_bs_div_focal: "too divergent from the focal (the broad-sample test, stage 1)",
    // the controls
    dots: "dots where a sequence matches the focal residue", whole: "whole locus", letters: "letters", trees: "trees",
    trees_hint: "the expected and the observed distances among these sequences as two neighbour-joining trees, face to face",
    // the column heads
    h_test: "test type", h_role: "sequence", h_name: "Taxon name", h_p: "p<α",
    h_p_hint: "significant: * below alpha / N ({bar}), ** below alpha / 5N, *** below alpha / 50N", resize: "drag to resize",
    // the footer and the site hover
    foot: "{n} sites · pinch or ctrl+wheel: zoom · drag: move · letters from 7 px per site",
    k_site: "site", site: "{s} of {len} · alignment position {pos}",
    k_diff: "difference", diff: "{d} of {n} relatives differ from the focal{note}", diff_most: " (most)", diff_few: " (fewer than 3 present)",
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
    verdict_seq_masked: "Masked by {n} test{s}", verdict_seq_clean: "Not masked",   // the sequence card's verdict line
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
    r_cr: "{nf} of {nc} too divergent · observed {obs}% · expected {exp}% · p {p}",
    r_cm: "closest match {name} · {arm} · p {p}",
    r_cm_arm_rel: "a relative", r_cm_arm_out: "outside the relatives",
    r_bs_closest: "closest member {name} · {mis} sites differ · {exp} expected · p {p}",   // the broad test's two rows (0.14)
    r_bs_majority: "{nfail} of {nt} too divergent · observed {obs}% · expected {exp}% · p {p}",
    r_bs: "variance {v}× expected · p {p}{call}",                                        // older reports' one row
    r_bs_divergent: " · sequence too divergent",
    r_bs_close: " · sequence too close",
    footer: "{dataset} · LocusMasker {version} · {date} · every test held to p < {bar}",
  },

  // ---- the topbar: the one-line count under the title, and the "about" fold with the run's provenance.
  // Every {x} is a field of the report's manifest, written by the tool.
  about: {
    sub: "{masked} of {cells} sequences masked ({pct}) · {taxa} taxa × {loci} loci",
    summary: "{tool} {version} · {date}",                                            // the fold's closed line
    run: "commit {commit} · python {python} · numpy {numpy} · {system}",
    input: "input {file} ({mb} MB, sha256 {sha}…){partitions}",
    input_partitions: ", partitions {file}",
    options: "options: datatype {datatype}, n-comp {n_comp} (K {K}, majority {M}), fit-overlap {fit_overlap}, alpha {alpha}{flags}",
    opt_similar_only: ", --mask-similar-only", opt_different_only: ", --mask-different-only",
    test: "test: {summary}, dataset dispersion {phi}; carried: {seqs} sequences, codes for {codes} cells, best matches for {bm} cells",
    // {summary}: the run's tests and their one bar
    tests: "{n} tests ({cr} close-relative, {cmr} closest-match against a relative, {cmo} against an outside match{bs}), one Bonferroni bar for all of them, p < {bar}",
    tests_bs: ", {n} broad-sample",
    tests_old: "{n} tests ({div} too divergent, {sim} too similar, {dir} direction), ",   // reports before the test names
    tests_bar: "one Bonferroni bar for all of them, p < {bar}", tests_k: "Bonferroni k {k}",
    tests_oldest: "{n} distinct comparisons, Bonferroni k {k}",
    command: "command: {cmd}",
    licences: "LocusMasker is MIT-licensed, © Robert Lanfear. This page embeds D3 (ISC), d3-sankey (BSD), Fuse.js (Apache 2.0) and the Inter font (OFL 1.1): the licence texts",   // the fold inside the fold; opens on the full texts
  },

  // ---- a cell's state, in words: the matrix legend's hovers, the cell hover's first row, the histograms' hovers.
  // {k} of {K}: comparators failed; {plus} is "+" on the top grade (3 or more); {s} is "s" when plural.
  state: {
    absent: "absent", untested: "present, not tested", not_masked: "not masked",
    red: "masked, too divergent · fails {k} of {K}",
    blue: "masked, copied from another taxon · {k}{plus} too-similar partner{s}",
    blue_old: "masked, too similar · {k}{plus} pair{s}",                                 // reports before format 4
    both: "masked by both tests: too divergent from its relatives and copied from another taxon · {k}{plus} too-similar pair{s}",
    blue_many: "masked, copied from another taxon · {n} too-similar partners",         // the exact count where the grade says 3+
    bs_div_and_close: "masked, too divergent · also significantly too close to part of the broad sample",
    bs_close: "masked, too close to the broad sample",
    bs_divergent: "masked, too divergent from the broad sample",
    bs_copy_divergent: "masked: copied from another taxon · also too divergent from the broad sample",
    // the three mask families (the histograms' segments)
    fam_red: "masked, too divergent", fam_blue: "masked, copied from another taxon", fam_blue_old: "masked, too similar", fam_both: "masked by both tests",
    // why a cell was not tested; {status} is the tool's word for the taxon or the locus (excluded, skipped, …)
    reason_no_residue: "the taxon has no residue at this locus",
    reason_taxon: "taxon {status}", reason_locus: "locus {status}",
    reason_no_comparator: "no comparator shares enough sites with the taxon here",
  },

  // ---- the cell's partner (its closest match, its copy or source, or the broad-sample comparator that decided it):
  // the TAG names the row, the TEXT is the cell hover's sentence, the NOTE is the trees' explanation. {name} = the partner.
  partner: {
    bs_donor_tag: "close broad comparator",
    bs_donor_text: "{name}: a far broad-sample comparator this divergent sequence is significantly too CLOSE to — the donor side; the divergent mask stands, this names where it points",
    bs_donor_note: "the most extreme too-close broad-sample comparator (the donor side)",
    bs_far_tag: "far broad comparator",
    bs_far_text: "{name}: significantly too DIVERGENT from this sequence for their genome-wide divergence — contamination or a chimera confined to sites its close relatives do not cover (the broad-sample test)",
    bs_far_note: "the most extreme too-divergent broad-sample comparator",
    bs_close_tag: "close broad comparator",
    bs_close_text: "{name}: significantly too close to this sequence for their genome-wide divergence — a coalescence younger than the species split allows (the broad-sample test)",
    bs_close_note: "the most extreme too-close broad-sample comparator",
    pair_tag: "partner", source_tag: "source", copy_tag: "copy", closest_tag: "closest match",
    uncertain_text: "{name}: the two are significantly too similar, and neither fits its own relatives clearly worse than the other (within the donor threshold), so the direction is uncertain and both are masked{off}",
    uncertain_note: "the other member of a too-similar pair whose direction is uncertain; both masked",
    copy_text: "{name}: this sequence is closer to it than its own relatives are, and fits its relatives worse than {name} does, so it is the copy{masked}",
    copy_note: "the source the focal sequence was copied from",
    source_text: "{name} is a copy of this sequence: {name} is closer to it than its relatives are and fits its own relatives worse; this sequence is the source, not masked for it",
    source_note: "a copy of the focal sequence",
    closest_text: "{name}, its closest match at the locus, is not closer than its own relatives are",
    closest_note: "the closest match at the locus; not surprisingly close",
    masked_sfx: ", masked", off_sfx: "; similar masking was off",      // {masked} / {off} above
    // reports before format 5 (a pair's direction from the LRT)
    old_recipient: "{name} (identified as the recipient)", old_undetermined: "{name} (direction not identified)", old_pair_note: "the partner of the too-similar pair",
    old_source_text: "{name}, its source: this sequence is the identified recipient (the copy){masked}",
    old_copy_text: "{name}, its copy: {name} is the identified recipient; this sequence is the source, not masked for it",
    old_copy_note: "the identified recipient, a copy of the focal sequence",
    old_undet_text: "{name}, direction not identified: neither sequence is masked for it",
    old_undet_note: "too similar to the focal sequence, direction not identified",
  },

  // ---- the overview's controls (the view picker, the matrix's orders and finder, the exports) and its small hovers
  controls: {
    view: "view", taxa: "taxa", loci: "loci", sort: "sort",
    view_summary: "summary", view_matrix: "matrix", view_loci: "loci ribbons", view_taxa: "taxon ribbons", view_table: "table",
    order_input: "input order", order_masked: "most masked first", order_name: "name", order_length: "length", order_grouped: "grouped by shared masks",
    search: "find a taxon or locus  /",                                              // the finder's placeholder; / focuses it
    help: "pinch or ctrl+wheel: zoom (shift: loci only, alt: taxa only) · two-finger scroll or drag: move · shift+drag: zoom to a box · double-click: fit · click a label: open its panel · click a cell: open the alignment · arrows move, + and − zoom, / finds, Esc closes",
    tile_taxon: "{masked} of {tested} tested loci masked", tile_locus: "{masked} of {tested} tested taxa masked",   // the worst lists' hovers
    rib_name: "{name}: {masked} masked of {tested} tested",                          // a ribbon row's name hover
    rib_count: "{n} {unit}", unit_loci: "loci", unit_taxa: "taxa",
    identical: "{n} sequences identical to their comparators are not drawn",
    // the histograms under the lists
    hist_hover: "{n} {unit} with {range} masked sequence{s}", hist_range: "{lo} to {hi}",
    hist_fam: "{label}: {n} sequences ({pct})", hist_more: " · click to highlight them in the matrix",
    hist_title: "zero is its own bar; the other bars are coloured by the kinds of mask they hold. The y axis is linear near zero and logarithmic above. Click a bar to highlight its units in the matrix.",
  },

  // ---- the hovers: a cell (the matrix, the scatter, the ribbons' points) and a unit (a row or column label).
  // k_* are the left-column labels; the others are the values. {bar} is the run's bar where the report carries one.
  hover: {
    many: "{n} cells under the pointer, worst first:", more: "… and {n} more; zoom in to see them",
    also_copied: " · also copied from another taxon", also_similar: " · also masked as too similar", fewer_comps: " · fewer comparators than requested",
    flagged_off: " (flagged; divergent masking was off)", copy_off: " · identified as a copy (similar masking was off)", similar_off: " · a similar-side fail (similar masking was off)",
    k_state: "state", state_reason: "{label}: {reason}",
    k_mask_p: "mask significance", mask_p: "p = {p}{stars} (the {test} test)",
    k_residues: "residues",
    k_pairs: "too-similar pairs", pairs: "{n}{plus}{more}", pairs_more: "; the strongest below",
    k_partners: "too-similar partners", partners: "{n}{plus} significant, among every taxon at the locus{more}", partners_source: "; its source below",
    k_copied_from: "copied from", k_too_similar_to: "too similar to", k_too_divergent_from: "too divergent from", k_too_close_to: "too close to", k_copied_by: "copied by", k_closest_match: "closest match",
    k_cm_test: "closest-match test", cm_test: "p = {p}{bar}",
    k_which_copy: "which is the copy", which_copy: "{who}, by the close-relative excess: {x0} m here against {xp} m for {name}{thr}",
    which_uncertain: "uncertain", which_worse: "the member that fits its own relatives worse", which_thr: " (the larger must be more than {thr} times the smaller to call it)",
    // reports before format 5: the direction test and its two factors
    k_dir_lrt: "direction LRT", dir_lrt: "|2 ln LR| = {lam}, {dir} (approximate χ², 1 df{bar})", dir_lrt_old: "|2 ln LR| = {lam}, p = {p} (χ², 1 df)", dir_bar: "; bar {bar}",
    k_row_factor: "row factor", row_factor: "{f} (the pair was tested against {f} × its expectation: {who} sits below expectation against its other comparators too)", row_one_of_two: "one of the two", row_the_sequence: "the sequence",
    k_depth_factor: "depth factor", depth_factor: "{f} (pairs this far apart sit at {f} × expectation at this locus; the pair was tested against row factor × depth factor × its expectation)",
    k_best_match: "best match", best_match: "{name} ({id}% identity over shared sites{close})", best_match_close: ", one of the 5 closest relatives",
    k_broad: "broad sample", broad: "{n} comparators drawn from the taxon's pool{p}{call}", broad_p: "; stage 1 p = {p}", broad_past: " — past the run's bar{bar}", broad_nothing: ", nothing unusual",
    k_rho: "own rate ρ", rho: "{rho} × the locus rate (the median over its comparators)",
    k_bs_close: "broad: too close", bs_close: "the far shortfall, p = {p} at α{covered}",
    k_bs_div: "broad: too divergent", bs_div: "half the far sample leans high, p = {p} at α{covered}", covered: " — another mask already covers the cell",
    k_bs_donor: "broad: the donor side", bs_donor: "the divergent mask stands; without ρ the far shortfall names the donor side, p = {p}{bar}",
    k_bs_lean: "broad: the lean", bs_lean: "{dir}, {name} the most extreme — not significant at α, nothing masked", lean_close: "toward too close", lean_div: "toward too divergent",
    k_comparators: "comparators", comparators: "{nc} tested, {nf} failed as too divergent{sim}", comparators_sim: ", {ns} as too similar",
    k_divergence: "divergence, median over comparators", divergence: "{obs}% observed, {exp}% expected",
    k_closest_comp: "closest comparator", closest_comp: "{name}: {obs}% observed vs {exp}% expected over {n} shared sites",
    k_excess: "excess beyond expectation", excess: "{x} m",
    k_cr_test: "close-relative test", cr_test: "p = {p}, the comparison that decides the majority{bar}",
    k_bs_role: "broad-sample role", bs_role: "excluded as a comparator at this locus (too divergent, or a masked copy)",
    click_cell: "click: open the alignment",
    bar: " (bar {bar})",
    // a taxon or a locus
    k_status: "status", k_sequences: "sequences", sequences: "{present} present, {tested} tested",
    k_masked: "masked", masked: "{n} ({pct}): {div} too divergent, {sim} {simword}; {near} near misses", sim_copied: "copied from another taxon", sim_similar: "too similar",
    k_length: "length", length: "{n} sites", k_rate: "rate", rate: "{r} × the typical locus", k_plateau: "plateau", k_phi: "dispersion φ", k_closest_by_m: "closest by m",
    click_unit: "click: open the panel",
  },

  // ---- the unit panel (one taxon or one locus; or the copies): its two views, their notes behind the "?", the
  // copies' blocks, and the empty states. {blue} picks a format variant; _old entries serve reports before format 4.
  unit: {
    view: "view", help: "what this plot shows", crumb_sub: "{n} of {m} sequences masked ({pct})",
    view_results: "test results", view_ribbons: "site ribbons",   // the two views of a taxon or a locus; the results are the unit's card, live
    // the copies
    copies: "copies", copies_n: "{n} copies", copies_of: "{k} of {n} copies", blocks: "blocks", help_copies: "what this shows",
    by_donor: "by donor", by_recipient: "by recipient", by_pair: "by donor-recipient pair",
    no_copies: "no copies here",
    head_recipient: "recipient", head_locus: "locus", head_donor: "donor", head_taxon: "taxon", head_best: "best match",
    pair_rule: "donor → recipient", block_one: "{label} · 1 copy", block_many: "{label} · {n} copies",
    copy_row: "{rec} at {locus}, copied from {don} · click: the alignment",
    note_copies_head: "Copies.",
    note_copies: " Each row is a copied sequence (the recipient, masked) at one locus, drawn along its sites: blue where it shares its donor's residue while two or more of its comparators differ, grey to red as more of its comparators differ from it. The donor is the member of the too-similar pair that fits its own relatives better. Blocks run from the most copies to the fewest. Click a row for the alignment.",
    // test results
    note_results_head: "Test results.",
    note_results: " Every tested sequence of this taxon or locus on each test's axis: its own masks coloured, the ten strongest numbered. Hover a point for its numbers, click it for the alignment.",
    // site ribbons
    no_masked: "no masked sequence here",
    head_best_hint: "the taxon closest to the sequence at this locus (per-site divergence compared on the same sites); in the outside colour when it is not one of the sequence's comparators",
    best_cell: "{name}: {id}% identity over shared sites{in}", best_in: ", one of the comparators", best_out: ", not one of the comparators",
    note_ribbons_head: "Site ribbons.",
    note_ribbons: " {n} masked sequences, sites along x: red where a majority of the comparators differ, blue {blue}. {rows} The right-hand column names the taxon closest to the sequence at that locus, in the outside colour when it is not one of its comparators. Hover a row for the site; click it for the alignment.",
    note_ribbons_blue: "where a copy shares its source's state while two or more of its comparators differ", note_ribbons_blue_old: "where a pair shares a state against the others",
    note_ribbons_locus: "Rows share the alignment's columns: a vertical band is a column problem (misalignment), a horizontal block in one row is a taxon problem (a wrong exon, a frameshift).",
    note_ribbons_taxon: "Each row is a different locus drawn to its own length on a shared axis.",
    no_codes: "codes not carried for this cell",
    site: "site {s} of {len} (alignment position {pos}): {code}",       // a ribbon's hover
    // the site codes (the ribbons' colours) in words; {n} comparators
    code_gap: "no residue in the focal sequence", code_differs: "the copy differs from its source here", code_differs_old: "the pair differs here",
    code_agree: "the partner agrees with the focal here while {n} comparators differ", code_agree_v4: "the copy shares its source's state while {n} of its comparators differ", code_agree_old: "the pair shares this state while {n} other comparators differ",
    code_differ: "{n} of the comparators present differ{maj}", code_majority: " (a majority)",
  },

  // ---- the two trees: the branches' hovers, the scale bar, and the explanation behind the "?"
  tree: {
    split: "this split is not in the other tree", og_branch: "the outgroup's branch, {len}%, drawn shortened", branch: "branch length {len}%",
    scale: "{bar}% corrected divergence", scale_old: "{bar}% divergence",
    head: "Two neighbour-joining trees",
    over: " over the focal sequence (orange), its comparators (black; red when significantly too divergent) and, when it has one, its too-similar partner or its closest match at the locus (blue when significantly too similar; a hollow dot when not{partner}); the colours are the alignment's, other taxa are grey. Hover a tip for its role.",
    over_v4: " over the focal sequence (orange), its comparators (black; red when significantly too divergent) and up to five of the taxa it is significantly too similar to (blue), its reported partner first ({n} here{partner}). Hover a tip for its role.",
    over_old: " over the focal sequence (orange), its comparators (black dots), the partner of a too-similar pair (blue label) and, when it is not one of the comparators, the best match (crimson dot and label; the taxon closest to the focal by per-site divergence compared on the same sites{id}).",
    partner_is: "; {name} is {note}", none: "none", identity: ", identity {id}%",
    bold: " In each tree the taxon closest to the focal by that tree's own distances is in bold: the closest relative genome-wide on the left, the closest at this locus on the right.",
    left: "Left:", right: "Right:",
    left_text: " the distances the test expects for every pair: each pair's genome-wide divergence, F81-corrected with this locus's plateau and scaled by its rate. ",
    left_text_old: " the distances the test expects for every pair, each pair's genome-wide divergence through this locus's expectation. ",
    right_text: " the divergence observed at this locus over shared sites, corrected the same way (capped at 0.98 of the plateau). ",
    right_og: "Both rooted on the same outgroup (hollow square, muted): the taxon at this locus most distant from these taxa genome-wide; its branch is drawn shortened, with a break. One scale over the ingroup, in corrected divergence; negative branches drawn as zero.",
    right_mid: "Both midpoint-rooted, one scale in corrected divergence; negative branches drawn as zero.",
    right_text_old: " the mismatch fractions observed at this locus over shared sites. Both midpoint-rooted, one scale; negative branches drawn as zero.",
    splits: "Thick branches are splits absent from the other tree; dotted lines join the same taxon on both sides. A black closest taxon on the right means the sequence has moved away from everything; a blue one means it has moved towards a taxon it is significantly too similar to.",
    reading: "A long focal branch with its neighbours unchanged is a damaged copy of the right sequence. A focal that moves across to a too-similar partner points to contamination or a paralogue; the direction test says which of the two is the copy. A whole clade pushed away from the rest on one long branch is a paralogue shared by the clade.",
    reading_old: "A long focal branch with its neighbours unchanged is a damaged copy of the right sequence. A focal that moves across to its best match points to contamination or a paralogue. A whole clade pushed away from the rest on one long branch is a paralogue shared by the clade. A best match on its own long branch is merely the least different taxon, not a relative.",
  },

  // ---- fatal states: shown in #app in place of the report.
  fatal: {
    no_data: "This page carries no data: it is the viewer template without a report. Run locusmasker.py; the report.html it writes has the data inside it.",
    unreadable: "The report could not be read.",   // followed by the error detail, for a bug report
    // the page's data folder is missing (report.html was moved without report_files beside it)
    no_folder: "report_files is not beside report.html. Keep the two together: the page reads its data from that folder.",
    // a report written by a NEWER tool than this viewer; {found} and {max} are format numbers, {version} the tool's
    newer_format: "this report was written in bundle format {found}; this viewer reads formats 1 to {max}. Please use the viewer that came with the tool version {version}.",
  },

};
// p-values for people: plain p, two decimals down to 0.01, scientific below (never log10 in a display).
// The inputs everywhere are log10 p (the lossless machine form); this is the one place they become words.
function lmP(lp) { if (lp === null || lp === undefined || Number.isNaN(lp)) return ""; const p = Math.pow(10, lp); return p >= 0.01 ? p.toFixed(2) : p.toExponential(1); }
// the viewer and the loader fill {placeholders} with this
function lmFmt(s, vars) { return s.replace(/\{(\w+)\}/g, (_, k) => vars[k] !== undefined ? vars[k] : "{" + k + "}"); }
