# Handoff: H2 Chemistry Reaction Drill (session 2 checkpoint)

Paste this whole file as your first message in the new conversation, with
this project folder open as the working directory. This picks up from
`HANDOFF_PROMPT.md` (the original handoff, kept in the repo root for
history) — everything in this document happened in the session(s) after
that one.

## What this is

A single-file web app (`app.html`, ~2200 lines, no dependencies, no build
step) that drills H2 Chemistry organic reaction questions. Molecules are
procedurally generated (never memorizable), and you answer by visually
building the product molecule in an embedded builder — not typing text.
Grading is now real (auto-graded against the reaction engine), not
self-reported.

`engine/` (the reaction/molecule-generation engine, unchanged since the
original handoff, still 17/17 tests passing via `node engine/test.js`) is
the source of truth for chemistry; `app.html` embeds it verbatim plus the
visual builder and quiz UI on top.

## Current status — what's built and verified

All of this has been tested directly (see "How this has been tested"
below — screenshot/visual testing was largely abandoned as unreliable for
this kind of precision work):

1. **Molecule builder** (embedded in `app.html`, a full rewrite of the
   original hex-grid prototype from `prior_outputs/reaction-memoriser.html`):
   - Free-form growth from one starting atom — no fixed grid. Bond angles
     are always 120° (≤3 bonds on an atom) or 90° (exactly 4), computed by
     a `layoutChildren` recursion that re-flows the whole affected area
     automatically when a 4th bond is added.
   - **Ghost-hover-to-grow**: hovering an already-placed atom shows small
     translucent "ghost" atoms at its open bond slots (sized to fit
     entirely inside the hover zone — this took two iterations to get
     right, see bugs below); clicking one grows there. Direct-clicking an
     already-active atom no longer auto-adds a branch (that was the
     original, confusing behavior) — it just selects it.
   - **Connect mode**: a third mode (alongside Place/Erase) to bond two
     *already-placed* atoms together with whatever bond order is currently
     selected — this is how you close a ring by hand (e.g. build a 6-chain,
     connect the ends) instead of only via the dedicated ring tool.
   - Dedicated **Benzene ring** tool still exists as a one-click aromatic
     hexagon shortcut, Kekulé-alternating.
   - Atoms/bonds/functional groups/undo/erase all carry over from the
     original prototype's design, adapted to the new free-form model.
   - Panning (drag empty canvas), a persistent "current bond order"
     selector (always exactly one of Single/Double/Triple selected, used
     for every new bond and every bond-order edit).
2. **Visual question rendering**: the question's starting molecule renders
   as an SVG (same builder engine, read-only) with its Hill-system formula
   below it — replacing the old plain condensed-text display.
3. **Seeded builder**: clicking "+ Build Molecule" opens the builder
   *pre-populated* with the question's starting material, so you edit it
   into your answer instead of rebuilding it from scratch.
4. **Real auto-grading**: `window.__normalizeMoleculeForGrading(builderGraph)`
   translates whatever you build back into `engine.js`'s own molecule
   format (or a clear rejection reason — see rules below), and grading
   compares your submission against the engine's actual computed answer via
   `engine.js`'s existing `answerEqual`/`productSetEqual`. The self-report
   "Got it right/wrong" buttons are gone — there's a computed ✓/✗ verdict
   now. A "+ CO2 + H2O byproduct" toggle exists next to "No reaction" for
   oxidative-cleavage questions whose answer includes that byproduct
   (building CO2 atom-by-atom isn't supported, and wasn't asked for).
5. **Theme**: the whole app (question panel + builder modal) uses one
   consistent light "lab notebook" palette (dark ink on light paper) —
   there used to be a dark quiz-panel theme that made molecules nearly
   invisible against it; that's fixed by unifying to one palette everywhere
   rather than patching colors in two places that can drift apart again.

## Grading rules (what normalizes vs. gets rejected)

`__normalizeMoleculeForGrading` deliberately **rejects rather than
guesses** on anything outside what `engine.js`'s molecule model can
represent — a wrong best-effort mapping risks exactly the "flukes" this
was built to avoid.

- **Supported**: plain carbon chains/branches; halogen (Cl/Br/F/I) and
  OH substituents (a raw single-bonded `O` atom is treated the same as the
  `OH` group chip — formula-equivalent); a double-bonded terminal `O`
  → `oxo:true`; `COOH`/`CHO` groups (the labeled node keeps its own
  identity with `group` set — it does **not** fold into its neighbor,
  that was a real bug, see below); a single plain 6-ring (built via the
  dedicated tool **or** by hand via Connect mode — cycle detection is
  generic, not tag-based) with proper Kekulé alternation and at most one
  external substituent, collapsing to `phenyl:true` on the *outside*
  atom the ring attaches to (not on any ring member).
- **Rejected outright**: `NH2`/`CH3`/`NO2` groups (no engine equivalent),
  triple bonds (engine model is single/double only), disconnected
  submissions, multi-substituted or fused/multiple rings, an in-chain
  heteroatom built by directly labeling a skeletal position.
- **"Which carbon" (radical substitution) questions require naming every
  distinct product** to count as correct — confirmed by explicit user
  answer, not assumed. Partial submissions are marked incorrect.
- If **any** submitted molecule fails to normalize, the whole answer is
  shown as "ungradable" with the specific reason, not silently marked
  wrong with no explanation.

## Real bugs found and fixed this session (why the tests below exist)

These were all found by *actually running* the code (Node simulations
and/or live browser JS-driving), not by inspection — worth knowing since
they're the kind of thing that could resurface if this code gets touched
again:

1. **Chain-spiral bug**: the layout's rule for continuing a plain
   2-bonded chain always picked the same rotational slot, so any chain
   past ~6 atoms retraced a closed hexagon and exactly overlapped itself
   (proven by simulation: node 0 and node 6 landed on the identical
   point). Fixed by alternating the slot choice by depth-parity.
2. **Cross-branch collisions**: separate branches (e.g. a long chain and
   a phenyl ring) could occasionally grow into each other's space — a
   different, global problem the local zigzag fix can't solve. Fixed
   pragmatically: detect it (`hasLayoutOverlap`) and have question
   generation just try a different random molecule, rather than building
   real collision-avoidance layout.
3. **COOH/CHO fold direction was backwards**: originally merged the
   labeled node into its *neighbor* instead of keeping its own identity —
   this also meant "branched COOH must be rejected" was a wrong
   assumption (it's actually representable in engine.js's model; now
   accepted).
4. **Phenyl ring-collapse kept the wrong node**: `phenyl:true` belongs on
   the *outside* carbon a ring attaches to; the ring's own 6 members are
   never represented at all (fully opaque, like the original spec's
   `phenyl` flag). Had this backwards initially.
5. **Ring validation checked `symbol` truthiness directly** instead of
   the `isPlainCarbon` helper, so a hand-built ring's root atom (always
   explicitly labeled 'C') was wrongly rejected as "not a plain carbon."
6. **Ghost hover gap**: ghosts originally rendered at full bond length
   outside the hover-detection zone, so reaching one meant leaving the
   zone that revealed it, causing it to vanish mid-click. Fixed by
   shrinking + repositioning ghosts to fit entirely inside the hover zone.

## Deliberately out of scope

- CO2/H2O byproduct is a toggle, not something you build atom-by-atom.
- No stereochemistry anywhere (not asked for beyond confirming enantiomers
  are already a non-issue — the graph model has no 3D representation at
  all, so it's structurally incapable of distinguishing them).
- Halogen derivatives / arenes topics still not wired into the generator
  (only Alkanes + Alkenes, per the original handoff's scope).
- No persistence beyond the existing manual JSON export/import of
  progress stats (unchanged from the original handoff).

## How this has been tested

**Screenshot-based visual testing was largely abandoned as unreliable**
for this precision work (the sandboxed browser's `file://` preview
renders as a static, non-interactive snapshot; coordinate-based clicking
drifted from real element positions). What actually worked, consistently:

1. Serve the app over a real local HTTP server (`python3 -m http.server`
   in the project root, then `preview_start`/`navigate` to
   `http://localhost:PORT/app.html`) rather than opening the file
   directly — `file://` behaves differently in the sandboxed browser tool.
2. Drive everything via `javascript_tool` dispatching real DOM events
   (`dispatchEvent(new MouseEvent('click'/'mousemove', ...))`) directly on
   elements found via `document.querySelector('[data-type="..."]')`,
   rather than coordinate-based `computer` clicks.
3. Cross-check computed results against `engine.js`'s own
   `canonicalForm`/`productSetEqual` (all exposed via
   `window.__registry['engine']`, `['operators']`, `['generator']`,
   `['render']`) rather than eyeballing — e.g. build the same molecule two
   different ways and confirm identical canonical output; build a known
   structure and diff its canonical form against a hand-constructed engine
   mol.
4. Two exposed debug entry points make this possible from outside the
   builder's closure: `window.__openMoleculeBuilder(onSubmit, seedMol)`
   and `window.__normalizeMoleculeForGrading(builderGraph)` — both were
   used extensively for isolated testing, not just by the real UI.
5. **A gotcha worth remembering**: `render()` inside the builder rebuilds
   the *entire* SVG on every click, so a `querySelectorAll` snapshot taken
   before a loop of clicks goes stale after the first one — re-query
   fresh inside the loop, don't cache element references across clicks.
6. **Another gotcha**: `window.__openMoleculeBuilder` takes the seed
   molecule as its *second* argument — a test that opens it with only a
   callback silently falls back to a blank canvas (this looks like normal,
   correct behavior and can waste time before you notice).
7. **Another gotcha**: since `+ Build Molecule` now seeds with the
   question's starting material (feature described above), any test
   script that assumes a *blank* canvas after clicking it will double-grow
   whatever chain it tries to build on top of the seed — check what's
   already there before growing more.

## What to do first in the new conversation

Read `app.html` (it's long — ~2200 lines — but it's the only file that
matters for anything UI/builder/grading-related) and `engine/*.js` if the
task touches chemistry. Run `node engine/test.js` to confirm the baseline
still passes before changing anything. Ask me what the next task actually
is — this document is a snapshot of where things stand, not a task list.
