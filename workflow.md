# Handoff: H2 Chemistry Reaction Drill — workflow.md

Paste this whole file as your first message in a new conversation, with
this project folder open as the working directory. This supersedes
`HANDOFF_PROMPT.md` and `HANDOFF_2.md` (kept in the repo root for history
only — both describe an ~2200-line early prototype; `app.html` is now
~12,300 lines and has been through a full v2.0 UI rebuild since then).
Read this file, skim `app.html`'s structure per the map below, then ask
what the next task is — this is a snapshot of where things stand, not a
task list.

## What this is

A single-file web app (`app.html`, no build step, no dependencies) that
drills H2 Chemistry organic reaction questions. Molecules are procedurally
generated (never memorizable) across most of the H2 syllabus (Alkanes,
Alkenes, Arenes, Halogen Derivatives, Amines, Hydroxy Compounds, Carbonyl
Compounds, Carboxylic Acids & Derivatives, Nitrogen Compounds). Four
question modes: guess the **product**, guess the **reactant**, guess the
**reagents & conditions** (typed tags), or a **distinguishing-test** table
mode. Answers are built visually in an embedded molecule builder (not
typed structures) and auto-graded against the reaction engine.

Served for local dev via `python3 -m http.server 8791` from the repo
root, then browse to `http://localhost:8791/app.html` — do **not** open
via `file://`, several things behave differently there.

Persistence: `localStorage` only by default (stats/attemptLog/mistakes/
settings), scoped per browser+origin — no account, no server. As of this
session there's an **optional** Google Sign-In + Firestore cloud sync
layered on top (see its own section below) — localStorage stays the
source of truth either way; syncing is opt-in and additive, never
required to use the app.

## File map

- **`app.html`** — everything: markup, CSS, and 8 `<script>` blocks plus
  one `<script type="module">` (the optional Firebase sync, see below;
  the syntax-check snippet further down needs to skip it — an ES module's
  top-level `import` isn't valid inside a plain `Function()` body, that's
  a checker limitation, not a real error). The
  only file that matters for anything UI/builder/grading/stats-related.
- **`engine/`** — `engine.js` (core molecule model + canonical-form
  comparison), `operators.js` (the ~50 named reaction transforms),
  `generator.js` (procedural starting-material generation), `render.js`
  (Hill-formula text), `test.js` (144 assertions, run with
  `node engine/test.js`), `demo.js`/`gen_batch.js` (scratch tools).
  **These are mirrored, not imported** — `app.html` embeds a verbatim copy
  of `engine.js`/`operators.js`/`generator.js`/`render.js`'s logic inside
  its own script blocks (each wrapped as a tiny fake CommonJS module
  writing into a shared `__registry` object, since the app has no build
  step to actually `require()` anything). **Any fix to chemistry logic
  needs to land in both places** — `app.html`'s embedded copy (what
  actually ships) and the matching function in `engine/*.js` (what
  `node engine/test.js` exercises) — or they silently drift apart. This
  bit this project directly this session (see "Ring atom selection..."
  below): a fix landed in `app.html` first, then had to be manually
  ported into `engine/operators.js` plus a new regression test added to
  `engine/test.js` to actually cover it going forward.
- **`reagents_conditions_reference.md`** — generated reference doc (every
  POOL reaction's reagent text + accepted tag combos), already delivered
  to the user; regenerate via the scratch scripts described in an earlier
  session if the POOL data changes meaningfully.
- **`HANDOFF_PROMPT.md`, `HANDOFF_2.md`** — historical, superseded by this
  file. Not worth reading unless curious about the original prototype.
- **`vercel.json`** — deploy config, not touched recently.

## `app.html`'s script-block structure

Roughly, in source order:
1. **Molecule builder engine** (~lines 1900–6700) — the interactive
   builder's own node/edge graph model, rendering, hover-to-grow ghosts,
   selection, and the embedded `engine.js`/`operators.js`/`generator.js`/
   `render.js` copies (each closes with its own `module.exports =
   {...}; __registry['name'] = module.exports;`). This block is its own
   closure — nothing in it is visible to other blocks except via
   `window.__` exports (see below).
2. **App logic** (~lines 6900–12300+) — quiz flow, grading, stats,
   feedback, all-time dashboard, background blobs. A separate IIFE; pulls
   in engine functions via `__registry['operators']`/`['engine']`/
   `['generator']`/`['render']` (aliased locally as `op`, `engineMod`,
   etc.).

**Cross-script-block communication is exclusively via `window.__name`**
(e.g. `window.__openMoleculeBuilder`, `window.__renderStaticMolecule`,
`window.__normalizeMoleculeForGrading`, `window.__showFloatingTooltip` /
`__hideFloatingTooltip`). If a new helper needs to be called from the
other script block, it has to be exported this way — there is no shared
scope otherwise.

## Molecule model: two shapes for "a plain ring with one attachment"

This is the single trickiest recurring source of bugs in this codebase —
worth understanding before touching anything ring-related.

- **`ring:true` node**: a real 6-atom aromatic ring (6 nodes, Kekulé-
  alternating S/D edges tagged `.ring=true`, a `subs[6]` flat-substituent
  array). This is what ring-*reactive* operators need
  (`ringElectrophilicSubstitution`, `ringFriedelCraftsAlkylation`,
  `ringFlatSubSwap`, `ringTribromination`, ...) — they all do
  `mol.nodes.find(n=>n.ring)` and fail/refuse if it's missing.
- **`phenyl:true` flag**: an opaque marker on a *plain carbon* node
  meaning "an unspecified, unsubstituted benzene ring is attached here,
  don't ask about its internal structure." Original, much older
  representation from when only Alkane/Alkene side-chain questions
  existed (`radicalSubstitution`, `sideChainOxidationToBenzoicAcid`) —
  those operators explicitly key off `n.phenyl` and would break if the
  ring were expanded into real atoms instead.
- **`downgradeRingIfMono(mol)`** (in `operators.js`/its embedded copy):
  folds a `ring:true` node down to `phenyl:true` whenever it ends up with
  exactly one occupied position and one outside attachment — used both by
  ring-reactive operators on their own OUTPUT (so a newly-mono-substituted
  product compares equal to the old phenyl-flag shape) and by
  `__normalizeMoleculeForGrading`'s own step 5.5, unconditionally, on
  every hand-built "ring + one plain chain" shape a student submits.
- **`upgradePhenylToRing(mol)`** (added this session, in both
  `app.html` and `engine/operators.js`) — the exact inverse. Needed
  because the normalizer's unconditional fold above is *correct* for
  comparing a submitted answer directly against a precomputed correct
  answer (product/reagent mode), but *wrong* when the submitted molecule
  is about to be **re-simulated** (reactant mode — see next section):
  re-running a ring-reactive operator on a phenyl-folded candidate finds
  no `ring:true` node and silently reports `occurs:false`, so a genuinely
  correct hand-built reactant (e.g. plain ethylbenzene, answering a
  ring-chlorination question) read as wrong. Fixed by trying re-simulation
  against both the as-normalized candidate and its phenyl-upgraded
  version. See `reactantCandidateMatchesGiven` in `app.html` and the
  matching regression test at the bottom of `engine/test.js`.

## Grading chokepoints

- **Product mode**: `gradeAnswer(spec, result, ua)` →
  `computeSubmittedAnswerSet` (normalizes every submitted molecule via
  `window.__normalizeMoleculeForGrading`, "ungradable" if any one fails)
  vs. `getGuessTargetAnswerSet` → `computeCorrectAnswerSet`, compared via
  `engineMod.productSetEqual`.
- **Reagent mode**: `revealTagAnswer` → `getReagentAcceptableCombos()`.
  This is the mechanism behind "this reagent is chemically equivalent, why
  wasn't it accepted?" — it re-runs *every other POOL spec* against the
  current question's own reactant, and if a sibling spec's `run()`
  produces the exact same given product, that sibling's tag combo(s) get
  merged in as additional accepted answers (shown in the reveal as "Also
  acceptable: ..."). Confirmed working this session for Na(s)/NaOH(aq)
  cross-acceptance on phenol reactions; the same mechanism should cover
  PCl5/PCl3/SOCl2 for acyl chloride formation and Na(s)/NaOH(aq)/
  Na2CO3(aq)/NaHCO3(aq) for carboxylic acid salt formation, since they all
  route through the identical underlying operator call.
- **Reactant mode**: `revealReactantAnswer` → `reactantCandidateMatchesGiven`
  — no precomputed correct answer; instead re-simulates the submitted
  candidate through the current spec *and* every sibling spec sharing the
  exact same reagent/condition tag combo, accepting if any reproduces the
  given product(s) exactly (via `productSetEqual`). See the ring-model
  section above for the one real bug found in this path.
- **Distinguish mode**: separate subsystem, `DISTINGUISH_FAMILIES` (~line
  11280+) — each family has 2–4 members sharing a reagent but differing
  in observed *result*; the table/free-text grading wasn't deeply
  audited this session beyond confirming the one reported "wrong answer"
  case was actually correctly designed (two members sharing a reagent
  name but with different expected observations, which IS what
  distinguishes them).

## UI conventions worth knowing before editing CSS/JS together

- **`--ui-scale: 0.75` transform on `<body>`**: `getBoundingClientRect()`
  returns POST-transform (visual) pixels; CSS length values you *write*
  (height, translateY, grid-template track sizes) are PRE-transform. Any
  JS that measures a real rect and feeds it back into a CSS length must
  divide by `uiScale` first, or it double-scales. Also matters when
  reasoning about `scrollHeight`/`clientHeight` (both raw DOM properties,
  both consistently PRE-transform, safe to compare against each other —
  the trap is comparing one of those against a `getBoundingClientRect()`
  value from the same element).
- **Play-pane strict 2-colour palette**: `#playPane` scopes `--text`/
  `--bg`/etc. to a black/cream duotone; "paper" surfaces nested inside it
  (`.eqBox`, `.molVisualBox`, `.submittedCard`, `.answerRect`, `.tagBox`,
  `.distinguishMolCard`) re-scope those same variables to the *opposite*
  pairing, so anything nested inside a paper surface renders black-on-
  cream even though the pane itself is cream-on-black. `--good`/`--bad`
  are the one deliberate exception (real red/green, not black/cream) for
  correctness feedback. **Specificity gotcha**: `#playPane button` (ID +
  tag) beats a plain `.someClass` rule — any new button-like element
  needs its override written as `#playPane .someClass` to actually win.
- **Floating tooltips**: native `title` attributes don't render at all in
  this environment (confirmed by explicit user testing) — use
  `window.__showFloatingTooltip(anchorRect, html, opts)` /
  `__hideFloatingTooltip()` instead, or the `wireDataTipHover(el)` helper
  which wires an element's `data-tip` attribute to that mechanism with
  `autoFlip:true` (measures off-screen first, flips above/below and
  clamps horizontally to avoid edge-clipping).
- **`sizeEquationBoxes()`** shrink-cascade (runs after every render and
  after reveal): shrinks the **answer/submitted-card area first** (down
  to a 140px floor), only touches the **given equation boxes** as a last
  resort — this was flipped this session (used to be the other way
  around, which shrank the given molecules every time an answer was
  revealed).
- **Radar chart** (`buildRadarChartSvg`, All-time stats): non-square,
  asymmetrically-cropped viewBox — top crop assumes axis 0 always lands
  at 12 o'clock (deterministically true, safe for any topic count);
  bottom crop uses the worst-case label (currently "Carboxylic Acids &
  Derivatives", 3 wrapped lines) at the worst-case angle (sin ≤ 1 bound,
  safe for any topic count too). **If you change the crop math, you must
  also hand-update `.radarChartSvg`'s CSS `width: min(100cqw, Ncqh,
  900px)` constant** (`N = 100 * viewBoxWidth / viewBoxHeight`) — there's
  no way for CSS to read the SVG's own computed aspect ratio, so this is
  a manual sync point, same category of gotcha as the engine.js/
  operators.js file duplication above.
- **Molecule builder ring selection**: clicking a single ring atom/bond
  now selects *just that one item* (`selectedItems`) — ring atomicity
  ("delete/resize the whole ring as one unit") only gets applied via
  `expandForRingAtomicity()` at **action time** (Backspace, +/-, not at
  selection time), specifically so pressing **C** to arm a Connect
  gesture works from a single selected ring carbon (that handler requires
  `selectedItems.length===1`). The ring's decorative hexagon/circle glyph
  is `pointer-events:none` — it used to be a click target for "select the
  whole ring," which sat on top of (and silently blocked hover-to-grow
  on) every individual ring atom underneath it.
- **Background blobs** (`#bgBlobs`, top of `<body>`): purely decorative
  drifting gradient layer. Currently 3 blurred circles (`BLOB_COUNT`),
  `will-change: transform, opacity` (pre-promotes the GPU layer so
  movement doesn't stutter on the first transition), background-color
  re-rolls ~1-in-3 retargets (that's the one property here that forces a
  real repaint of the blurred layer — kept deliberately rare). This has
  been tuned down twice now for reported lag/glitching (9→4→3 blobs); if
  it comes up again, the color-reroll frequency and blob count are the
  two levers, in that order of likely impact.

## Testing / verification workflow

1. **Syntax check** (run after any `app.html` edit, catches typos fast) —
   note this skips the one `type="module"` script (the Firebase sync,
   see its own section below): its top-level `import` isn't valid inside
   a plain `Function()` body, so it would always "fail" a naive version of
   this check without that actually meaning anything:
   ```
   node -e "
   const fs = require('fs');
   const html = fs.readFileSync('app.html', 'utf8');
   const re = /<script(\s[^>]*)?>([\s\S]*?)<\/script>/g;
   let m, count=0, errors=0, skipped=0;
   while((m = re.exec(html))){
     const attrs = m[1] || '';
     const src = m[2];
     if(!src.trim()) continue;
     if(/type\s*=\s*[\"']module[\"']/.test(attrs)){ skipped++; continue; }
     count++;
     try { new Function(src); } catch(e){ errors++; console.log('Error in block', count, ':', e.message); }
   }
   console.log('Checked', count, 'classic script blocks,', errors, 'errors,', skipped, 'module block(s) skipped');
   "
   ```
   Expect "Checked 8 classic script blocks, 0 errors, 1 module block(s)
   skipped".
2. **Chemistry regression**: `node engine/test.js` — expect "144 passed, 0
   failed" (as of this session; will grow). Only exercises `engine/*.js`,
   **not** `app.html`'s embedded copy — remember to port any operator fix
   both ways (see file map above).
3. **Live verification**: serve via `python3 -m http.server 8791` from
   the repo root, then the Browser pane's `preview_start`/`navigate` to
   `http://localhost:8791/app.html`. Seed `localStorage` directly to force
   specific question types/topics for testing, e.g.:
   ```js
   var s = JSON.parse(localStorage.getItem('h2ChemDrill.settings.v1'));
   Object.keys(s.selected).forEach(function(k){ s.selected[k] = (k==='someSpecId'); });
   s.questionTypes = {product:true, reactant:false, reagent:false, distinguish:false};
   localStorage.setItem('h2ChemDrill.settings.v1', JSON.stringify(s));
   ```
   then `navigate` to the same URL again to pick it up.
4. **Driving the DOM**: dispatch real events on the actual element under
   the cursor (`document.elementFromPoint(x,y)`, then
   `el.dispatchEvent(new MouseEvent(...))` / `PointerEvent(...)`) rather
   than the `computer` tool's coordinate clicks — this app's SVG hit-
   testing and hover zones are precise enough that coordinate drift
   causes false negatives. `svg.dispatchEvent(...)` (dispatching on the
   container instead of the actual hit element) silently no-ops, since
   `ev.target` becomes the container, not whatever's visually there — a
   real repeat mistake worth remembering.
   Re-query elements fresh after every `render()` call (the builder
   rebuilds its whole SVG on every interaction) — cached references from
   before a re-render are detached and won't receive events. A plain
   `.click()`/`MouseEvent('click')` sometimes doesn't trigger listeners
   bound as `pointerdown` — dispatch `pointerdown`+`pointerup`+`click` in
   sequence if a click seems to silently no-op.
   The Browser pane's screenshot pixel dimensions are **not** the real
   viewport size (screenshots are compressed thumbnails) — always compute
   click coordinates from `getBoundingClientRect()` in the real DOM, never
   by eyeballing a screenshot's pixel grid.

## Cloud sync (Firebase) -- optional, off until configured

`localStorage` stays authoritative on this device either way; this is a
strictly additive layer for "the same progress on another device/browser."

- **Not yet configured on the deployed app** — `FIREBASE_CONFIG` in the
  `<script type="module">` block near the end of `app.html` still has
  placeholder values (`"YOUR_API_KEY"` etc). Until a real Firebase
  project's config is pasted in, that whole script throws inside its own
  try/catch, logs one `console.warn`, and `window.__fb` is simply never
  defined — every call site elsewhere checks for that before touching it,
  so the app behaves exactly as it did before this feature existed. The
  setup steps (create a free Firebase project, enable Google sign-in,
  create a Firestore database, the exact security rules to paste in,
  where to get the config object) are all written out as a comment
  directly above that script block in `app.html` — follow those, not a
  paraphrase, since they're the ones actually kept in sync with the code.
- **Architecture**: the module script is a thin wrapper only — auth +
  Firestore calls, nothing app-specific — exposed as `window.__fb`
  (`signIn`/`signOut`/`onAuthChange`/`fetchDoc`/`saveDoc`). All the merge/
  business logic (`pushToCloudSync`, `mergeCloudIntoLocal`,
  `renderAccountRow`, the `window.__onFirebaseReady` entry point) lives in
  the app-logic script instead, right after the whole-app-tour code,
  since it needs direct access to `stats`/`attemptLog`/`mistakes`/
  settings and their existing `save*()` functions. One Firestore document
  per signed-in user (`users/{uid}`), holding all four pieces of state
  together — simple atomic writes, no partial-sync states to worry about.
- **Merge policy** (runs every time sign-in resolves, including a reload
  while already signed in — deliberately idempotent, not a one-time
  migration): `stats` takes the max of each field per reaction (can't
  merge two independent counters exactly, but max can never erase a real
  attempt); `attemptLog`/`mistakes` union by their own unique key
  (`ts`/`id`) then re-sort/re-cap same as a normal local append; `settings`
  is cloud-wins outright (a preference, not accumulated progress). Verified
  this session with a mocked `window.__fb` (real Firebase SDK not
  reachable from this sandbox) standing in for the actual network calls —
  sign-in/out, first-sync-push, and a synthetic "device B" cloud doc
  merging correctly into local `stats`/`attemptLog`/`mistakes` were all
  exercised end-to-end this way. **Not yet tested against a real Firebase
  project** — do that before considering this fully proven, particularly
  the Google sign-in popup flow itself and Firestore security rules
  actually blocking cross-user access.
- **UI**: `#accountRow` is now the **last child of `.sideNavList` itself**
  (not a separate sibling row with its own `margin-top:auto`) — it moves as
  one unit with the nav buttons via the list's existing `margin-top:auto`,
  with a `border-top` divider above it. This was a deliberate fix this
  session: a standalone `#accountRow` with its own `margin-top:auto` caused
  the whole sidebar to visibly shift vertical position depending on which
  pane was active (different `paneInfoDesc` text heights pushed the nav
  differently) — nesting it inside the same flex list as the nav buttons
  makes it move with them as a single block instead. Hidden by default
  (`.accountRow` / `.accountRow.ready`) and only revealed once
  `window.__onFirebaseReady` actually fires. Status (idle/syncing/synced/
  error) renders as a small dot (`#accountStatusDot`, color + `title`
  tooltip) next to the email, not a separate text line — avoids the
  previous "Synced" line's own height changing sidebar layout.

## Interactive tutorial engine (`window.__runTutorial`)

Replaces the old single whole-app tour. One generic, reusable engine drives
**every** tutorial in the app — the first-visit intro, each pane's own "?"
tutorial, and the molecule builder's tutorial — all through the same code
path (`tutorialResolveTarget`/`tutorialReposition`/`tutorialEnterStep`/
`tutorialAdvance`/`tutorialFinish`, defined in the app-logic script,
right after `switchPane`/`PANE_ORDER`, where the old whole-app-tour code
used to live).

- **Step schema**: `{ title, body, target: fn()=>Element|null,
  waitForClick, actionHint, onEnter }`. `target` is a *function*, resolved
  fresh every time the step is entered or the window resizes (the element
  it points at often doesn't exist yet, or is currently hidden). A target
  that resolves to a real element with a zero-size rect (e.g. `#nextBtn`
  before anything's been revealed) is treated the same as "no target" —
  falls back to a centered card with a full-page dim, no floating spotlight
  glued to a stale 0×0 point.
- **`waitForClick: true`** is what makes a step *actually interactive*
  rather than narrated: it hides the engine's own Next button, adds a real
  (non-`preventDefault`, non-`stopPropagation`) click listener straight
  onto the target element, and lets the page's own handler for that element
  fire first — genuine app state changes — before advancing the tour
  second. A disabled button (e.g. Submit before an answer is entered)
  simply can't be clicked yet, which is exactly the gate you want for "type
  a reagent, then click Submit" style steps — no extra state-watching
  needed, the browser's own disabled-button behavior does it for free.
- **Spotlight mechanics**: one `#tutorialSpotlight` div using the classic
  `box-shadow: 0 0 0 9999px <dim>` trick — dims everything except a rounded
  cutout around the current target. `#tutorialOverlay` (the full-viewport
  container) and `#tutorialSpotlight` both need `pointer-events: none` —
  **this was a real bug found live this session**: without it, the overlay
  (an empty, transparent, but fully click-intercepting div sitting at
  z-index 1001 over the whole page) silently ate every click anywhere on
  the page, so `waitForClick` steps looked right visually but a real click
  on the spotlighted button never reached it. `.tutorialCard` then needs
  its own explicit `pointer-events: auto` so the card and its Skip/Back/
  Next buttons stay clickable.
- **`--ui-scale` double-scaling — the other real bug found live this
  session**: `#tutorialSpotlight`/`.tutorialCard` are `position:fixed`
  descendants of the scaled `<body>` (see the `--ui-scale` note above), so
  a raw `getBoundingClientRect()` value (already POST-scale/visual px)
  written straight back as a CSS `top`/`left`/`width`/`height` gets shrunk
  by `--ui-scale` a *second* time — the spotlight rendered visibly smaller
  and higher/left of the real target. Fixed by dividing every measured
  value by `uiScale` before writing it back (`tutorialUiScale()`,
  used throughout `tutorialReposition`/`tutorialPlaceCard`) — the exact
  same divide-back-out pattern `switchPane`/`switchAllTimePage` already
  use for the identical reason; worth grepping for `tutorialUiScale` as a
  precedent the next time a new `position:fixed` overlay needs to measure
  and reposition against a live element.
- **z-index**: `#mb-overlay` (molecule builder modal) = 1000,
  `#tutorialOverlay` = 1001, `.tutorialCard` = 1002 — specifically so the
  *same* engine can spotlight elements inside the open molecule builder
  too (confirmed working live: `MB_TUTORIAL_STEPS` spotlights
  `#mbAtomCatalog`/`#mbSvg`/`#mbAddBtn` correctly on top of the modal).
- **Per-pane dispatch** (`runPaneTutorial()`, wired to `#tutorialHelpBtn`'s
  click): checks `mistakesModeActive` first (Mistakes is a mode-switch on
  `#playPane`, not a distinct `activePane` value), then dispatches on
  `activePane` to `PLAY_TUTORIAL_STEPS` / `TOPICS_TUTORIAL_STEPS` /
  `STATS_TUTORIAL_STEPS` / `FEEDBACK_TUTORIAL_STEPS` /
  `MISTAKES_TUTORIAL_STEPS`.
- **`PLAY_TUTORIAL_STEPS`** (22 steps) is the one substantial script: intro
  → basic play cycle (spotlighting whatever real question happened to
  already be on screen — forcing a fresh real reaction-mode question first
  if that happened to be a Distinguishing Test) → all 4 question types
  demoed via **fixed, deterministic molecules** (`generateAlcohol([2])`,
  i.e. ethanol — never the random generator, so the walkthrough is
  reproducible), including a real correct-then-wrong Submit cycle for
  product mode to show the green/red verdict colors and ghost cards → a
  brief nav-button-only tour of the other panes (no real pane-switching —
  deliberately simpler and avoids racing `switchPane`'s own 350ms slide
  animation) → close. Demo answers are seeded directly via
  `window.__buildAnswerFromEngineMol`/`tutorialSeedCorrectMolecule`/
  `tutorialSeedWrong` rather than driven through the interactive builder,
  since the spec explicitly asked for the builder itself to only be
  mentioned in passing here — its own deep-dive is `MB_TUTORIAL_STEPS`.
  `tutorialRestorePlayState()` runs as the tour's `onFinish` (fires on
  natural completion *and* Skip alike) to replace whatever demo
  question/answer was left on screen with a real one.
- **First-visit auto-trigger**: `TUTORIAL_INTRO_SEEN_KEY =
  'h2ChemDrill.tutorialIntroSeen.v1'` (new key — the old whole-app-tour's
  `APP_TUTORIAL_STORAGE_KEY` no longer exists), same "once, then only via
  '?'" pattern the molecule builder's own `MB_TUTORIAL_STORAGE_KEY` already
  used. Runs `PLAY_TUTORIAL_STEPS` verbatim, per the spec's explicit
  requirement that the intro be identical to the Play pane's own tutorial.

## Recent work (2026-09-01 session)

- **Tutorial redesign**: the full-page dimming overlay used to spotlight
  targets is gone, replaced by a precise pulsing ring outline
  (`box-shadow` + `.tutorialRingPulse` keyframe) drawn directly around the
  target element — no more "everything looks unclickable" impression, and
  no more overlay/target misalignment since there's no separate dim layer
  to keep in sync. `MB_TUTORIAL_STEPS` copy rewritten to be explicit and
  literal ("press 'Backspace'", "press '1', '2', or '3'") with a new
  `.tutorialKey` keycap-badge CSS class used for every key mention.
- **Back-button soft-lock fixed**: `tutorialEnterStep` used to always
  re-arm `waitForClick` gating on entry, so navigating Back into an
  already-completed step (whose target was now covered by an open modal)
  left the tutorial waiting forever on an unreachable click. Fixed with
  `tutorialState.maxIndexReached` + an `isNewGround` flag computed once
  per step entry, threaded into both card rendering and the click-arming
  condition — Back now re-shows completed steps as already-satisfied
  instead of re-blocking on them.
- **Related occlusion bug found while fixing the above**: even after the
  gating fix, the ring itself could still be drawn around a target that
  was hidden behind a modal. Fixed via `tutorialResolveVisibleTarget`
  (uses `document.elementFromPoint` at the target's center, checks
  ancestor/descendant containment, treats a hit inside the tutorial's own
  card/overlay as non-occluding, and fails open — shows the ring — on an
  inconclusive/null check).
- **Grading-step reorder**: the "how grading works" explanation step now
  appears *after* the Submit step's real click (not before), matching the
  actual flow of building a molecule then submitting it.
- **Content-coverage audit against the syllabus PDFs**: user provided all
  8 topic lecture-note PDFs (`Chem notes.zip`) and asked for a rigorous
  learning-outcomes-vs-POOL audit. Ran via 8 parallel research agents (one
  per topic), then independently re-verified every claim before touching
  code. Root cause of the "CN and amines are missing" impression: most of
  those reactions exist but are filed under their *starting-material*
  topic (e.g. nitrile formation lives under Halogen Derivatives, not
  Nitrogen Compounds) — that's the app's deliberate filing convention, not
  a bug. On top of that, found and fixed genuine gaps:
  - New operators + POOL entries: alkene hydrohalogenation (HBr, HCl),
    alkene hydration (steam), alkene halohydrin formation (Br2(aq)),
    substituted-ring bromination (`methylbenzeneRingBromination`), alkane
    combustion, phenol trinitration (picric acid), amide reduction
    (`reduceAmide`, scoped to primary amide only — matches
    `hydrolyzeAmide`'s existing scope), and the iodoform test for ethanal
    specifically (`iodoformCleavage`'s find condition only matched open
    oxo carbons, not `generateAldehyde`'s closed `group:'CHO'` node shape)
    and for CH3CH(OH)-R alcohols (`iodoformCleavageFromAlcohol`, a new
    operator).
  - New `substitutionWeight`/`moreSubstitutedEnd` helper pair implements
    Markovnikov regiochemistry for the three new alkene operators — counts
    real bonded neighbors plus `subs.length` plus `+1` for `node.phenyl`
    (a boolean "unspecified ring attached here" flag, not a graph edge —
    easy to undercount if you forget it). Halohydrin formation is the odd
    one out: OH goes on the *more* substituted carbon and Br on the
    *less* substituted one (water attacks the bromonium ion's more
    electrophilic carbon), opposite of plain hydrohalogenation/hydration.
  - Real bug fix along the way: `nitrobenzeneReduction` listed an
    incorrect LiAlH4 alternate reagent.
  - Deliberately scoped out (documented in-code near the Nitrogen
    Compounds POOL section): amino acids/proteins — no amino-acid
    molecule representation exists anywhere in the engine, so this is a
    new molecule type from scratch, not a new operator. Also skipped: HI
    hydrohalogenation, anti-Markovnikov minor product as an alternate
    answer, and a couple of minor no-reaction/reverse-reaction gaps
    (I2+alkane, ammonium salt + NaOH regenerating the amine) that are
    trivial inverses of reactions already covered.
  - `engine/operators.js` and `engine/test.js` mirrors updated to match
    (new operators, new/fixed tests, corrected a pre-existing test that
    had asserted incorrect chemistry — ethanal was wrongly asserted to
    fail the iodoform test). `node engine/test.js`: 162/162 passing.
  - POOL now has 105 entries total. `reagents_conditions_reference.md`
    was NOT regenerated (no generator script exists in the repo) and is
    now stale relative to these additions plus the `nitrobenzeneReduction`
    fix — regenerate by hand or write a script next time it matters.

## Recent work (prior session)

- **Feedback pane**: added a second "Mistake in the question?" link next
  to the verdict badge (in addition to the top-right "Report question"
  button), matching `.linkBtn`'s simple opacity-fade hover (an earlier
  dual-scope hover rule was visually broken — fixed by simplifying, and
  separately fixed a real `border` leak from `#playPane button`'s higher-
  specificity rule that was drawing an unwanted box around it). Feedback
  context line simplified to a fixed "Question is attached automatically"
  (the full topic/reaction/mode/difficulty detail still gets sent, just
  isn't spelled out on screen).
- **Reveal/answer-box layout**: Submit/Next buttons moved to live *inside*
  the answer box itself (bottom-right, solid black pill), persisting
  through reveal instead of a separate row below; `sizeEquationBoxes`'s
  shrink priority flipped (see above) so revealing an answer no longer
  shrinks the given equation.
- **Molecule builder**: CHO now draws as a real C=O + explicit H (not one
  opaque box) matching how COOH already did; a genuine grading bug fixed
  where the azo-coupling shape (`Ar-N=N-Ar'`) was incorrectly rejected as
  an invalid amine nitrogen; charged atoms with implicit hydrogens now
  spell out the H count (`NH3+`, not the ambiguous `N+`) in both the live
  builder and static/given molecule rendering; ring atom selection fully
  reworked (see UI conventions above) — this was the main piece of work
  this session, fixing a real "hand-built correct reactant marked wrong"
  grading bug for ring-substitution reactant-mode questions.
- **All-time stats dashboard**: radar chart widened/enlarged via tighter
  top+bottom viewBox cropping (see UI conventions above); grid columns
  rebalanced; Strongest/Weakest topic cards reverted to single-column
  stacked layout; height floors raised to their true safe minimums after
  discovering the old ones only "worked" by accident (see below);
  scroll eliminated at 1280×800/1366×768/1440×900 (was previously
  overflowing even at 1366×768).
- **A subtle bug found while doing the above**: `grid-template-rows`
  using a bare `1fr`/`1.3fr` resolves to its row's own MAX-CONTENT size
  whenever the grid is measured under an indefinite/intrinsic height
  (which is what happens here, since the page is a `flex:1` child with no
  `min-height:0`) — this made two rows balloon to a fixed size regardless
  of viewport, which is what was actually causing the "needs to scroll at
  1366×768" complaint. Fixed by capping the max side in `vh` instead
  (behaves as a normal length, no max-content special case) — but that
  then legitimately started hitting the row floors, which turned out to
  be too low for their content (an unnoticed latent bug, previously
  masked by the fr-inflation always over-provisioning). Raised the floors
  to their real minimums and added `flex-shrink: 0` to the affected
  cards' text children so a floor that's ever hit again fails visibly
  (overflow) rather than silently blanking text.
- **Background blobs**: perf-tuned (see UI conventions above) in response
  to a reported glitching complaint.
- **Sidebar nav reorder + sign-in row fix**: nav order changed to Play →
  Mistakes → Topic and Question Selection → Stats → Feedback (was Play →
  Topics → Mistakes → Stats → Feedback). `#accountRow` moved to be the last
  child of `.sideNavList` itself (see its own note above) to stop the
  sidebar visibly shifting position between panes.
- **Whole-app tutorial replaced with a generic per-pane engine**: see its
  own "Interactive tutorial engine" section above. The single whole-app
  tour and the molecule-builder-specific tutorial UI (two separate,
  non-interactive, plain-text modal implementations) are both gone,
  replaced by one `window.__runTutorial(steps, onFinish)` engine driving
  five step lists (`PLAY_TUTORIAL_STEPS`, `TOPICS_TUTORIAL_STEPS`,
  `STATS_TUTORIAL_STEPS`, `FEEDBACK_TUTORIAL_STEPS`,
  `MISTAKES_TUTORIAL_STEPS`, plus the pre-existing `MB_TUTORIAL_STEPS`
  retrofitted onto the same engine). Every step now spotlights a real DOM
  element, and several steps require an actual click on the real,
  fully-functional UI (Submit, Reveal answers, place an atom on the
  builder canvas, ...) rather than just reading text. Live-verified this
  session end-to-end: first-visit auto-trigger, the "?" button dispatching
  correctly on every pane and in Mistakes mode, all 22 steps of the Play
  tour (including real Submit clicks producing real green/red verdicts and
  a real typed-tag reagent answer), the molecule builder's retrofitted
  tutorial spotlighting elements inside its own modal, and Skip cleanly
  restoring real play state mid-tour.

## Feedback backlog reviewed this session (not all fixed)

The user pasted several weeks of feedback-form submissions; each was
checked against current behavior. Confirmed-fixed items are folded into
"Recent work" above. Still-open, lower-confidence, or by-design items —
worth revisiting if they resurface:

- **Tollens'/Fehling's product is shown as neutral RCOOH, not the
  carboxylate (RCOO⁻/RCOONa) the alkaline reagent would actually give** —
  chemically correct feedback, currently a deliberate simplification
  (`oxidizeAldehyde` doesn't track which specific oxidant/pH condition
  was used, only whether RCHO→RCOOH occurs at all). Fixing this properly
  means threading "was this an alkaline or acidic oxidant" through to the
  operator, which touches every question in the Carbonyl topic that uses
  it — worth a deliberate scoped task, not a quick patch.
- **Distinguish-test "answers are wrong" (Carbonyl and acid derivatives
  family)** — traced and believed correctly designed (two members share a
  *reagent* but differ in observed *result*, which is what actually
  distinguishes them) but not exhaustively verified against the actual
  grading/comparison logic for the table-entry mode.
- **"Amide hydrolysis, can[not?] provide the correct answer"** — feedback
  text was too garbled to act on; needs a repro (screenshot or exact
  molecule) if it comes up again.
- Several "can't tell which specific reagent was used without side
  info" reports (Na metal vs NaOH, PCl5 vs SOCl2, etc.) turned out to
  already be handled by `getReagentAcceptableCombos`'s cross-spec
  acceptance (confirmed live for one case, reasoned-but-not-individually-
  reproduced for the others, since they all route through the identical
  underlying operator call).

## What to do first in a new conversation

1. `node engine/test.js` — confirm the baseline still passes (144/144).
2. Skim this file's "grading chokepoints" and "molecule model" sections
   before touching anything ring- or reactant-mode-related — that's where
   the sharpest edges are.
3. Ask what the actual next task is — this file is a snapshot, not a plan.
