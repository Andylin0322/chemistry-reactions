# Handoff: H2 Chemistry Organic Reaction Practice App

Paste this whole message as your first message to Claude Code, with this
entire folder attached/available in the working directory.

## Who this is for

I (the user) am a pre-university student in Singapore studying H2 Chemistry
(A-Level syllabus). I struggle with memorising organic reactions and mixing
up reagents/conditions/products. I'm building a web app to drill this —
NOT flashcards, because flashcards don't work for me.

**Important working style:** I want strict adherence to my actual syllabus
content — no invented reactions or rules that aren't either (a) explicitly
in my reaction-rule JSON files, or (b) a standard, clearly-flagged extension
of them. If you're inferring or extrapolating beyond the provided data
(e.g. a substitution-pattern rule, a regiochemistry rule), say so explicitly
so I can check it against my lecture notes. Don't just present it as fact.

## The core idea (this took a while to arrive at, don't relitigate it)

Static question banks — even hand-written "hard" ones — get memorised after
a few passes. The fix: a **procedural generation engine** that computes
reaction products algorithmically from the actual mechanism, applied to
randomly generated molecules, so questions are never identical twice and
can't be memorised as question→answer pairs.

The elegant part, already validated: five distinct "trap categories" I
originally thought I'd need to hand-design (chemoselectivity, benzene as
an unreactive spectator, multi-site tracking, branching complexity, "which
carbon" ambiguity) **emerge automatically** from applying a small number of
correct reaction operators to randomly generated molecules. No special-case
trap logic needed — just correct chemistry.

## What's built and tested (in `engine/`)

- `engine.js` — molecule data model (see below), tree canonicalization for
  order/direction-independent equality checking, answer-equality functions.
- `generator.js` — random molecule generator with a **difficulty parameter
  (1-5)** controlling chain length, branching, number of double bonds,
  phenyl-group probability, and pre-existing-halogen probability. Accepts
  constraints (e.g. `minDoubleBonds`, `forcePhenyl`, `forceExistingHalogen`,
  `forceNoDoubleBonds`) so a specific reaction type can force a molecule
  that's actually applicable to it.
- `operators.js` — the actual reaction logic, as pure functions on the
  molecule graph, not stored answers:
  - `radicalSubstitution(mol, X)` — finds every C–H site, returns all
    constitutionally distinct monosubstitution products (deduped via
    canonical form). This is what naturally produces "which carbon"
    questions and further-substitution questions, for free.
  - `additionX2(mol, X)` — electrophilic addition of X2 across every C=C.
  - `hydrogenation(mol)` — H2/Ni, saturates every C=C.
  - `mildOxidationDiol(mol)` — cold dilute alkaline KMnO4, diol across
    every C=C.
  - `oxidativeCleavage(mol)` — hot acidified KMnO4. Breaks every C=C;
    each resulting carbon becomes CO2+H2O (if it had nothing else attached),
    -COOH (if one other chain neighbour), or a ketone via `oxo` (if two
    other chain neighbours/branch point). Handles multi-fragment splitting
    for polyenes correctly — verified against hexa-1,4-diene AND
    independently reproduced the textbook answer for hexa-2,4-diene
    (2x ethanoic acid + oxalic acid) without being told that answer in
    advance.
- `render.js` — condensed structural formula text renderer (e.g.
  `CH3-CH2-CHBr-CH(-CH3)-CH2-CH3`), recursive over branches. **I explicitly
  do not want aesthetics/graphics right now** — this text format is fine.
- `test.js` — 17 unit tests, all passing: hand-verified answer
  reproduction, symmetry deduplication, phenyl-as-spectator check, and a
  1500-molecule valence-validity stress test across all 5 difficulty
  levels. **Run this after any change to the engine before trusting output.**
- `demo.js`, `gen_batch.js` — example scripts showing the engine in action
  and generating batches of question/answer pairs.

## Molecule data model

```
node: { id, subs:[...], oxo:bool, group:null|'COOH'|'CHO', phenyl:bool }
edge: { a:id, b:id, type:'S'|'D' }
molecule: { nodes:[...], edges:[...] }
```

It's a **tree** (branching allowed, no rings in the graph itself). A benzene
ring is deliberately modelled as an opaque `phenyl: true` flag on a node,
NOT as a real ring subgraph — this was a deliberate simplification, and it's
why the ring correctly never reacts under any operator (nothing touches the
`phenyl` flag). Every carbon has 4 valence slots; unused slots are implicit
H, computed from bond orders + substituent count (+2 for oxo, +1 for
phenyl — this was a bug I already found and fixed once, watch for it if you
extend the model further).

## Known bugs I already found and fixed (context for why tests exist)

1. `subSlotsUsed` originally didn't count the phenyl substituent as using a
   valence slot — caused both generation and rendering to be wrong whenever
   a phenyl group was present. Fixed.
2. The renderer originally only printed the *first* atom of a branch and
   silently dropped the rest if a branch was more than one carbon long
   (e.g. printed `CH(-CH2)` for a two-carbon ethyl branch, hiding the second
   carbon entirely). This produced structurally wrong/misleading molecules
   in output. Rewrote as a proper recursive renderer. Watch for this class
   of bug generally — silently-wrong chemistry output is much worse than a
   crash, since I might not catch it by eye.

## Deliberately out of scope / not yet built

- **HX addition and steam addition to unsymmetrical alkenes are not
  implemented as operators** — these require Markovnikov's rule, which
  isn't stated in my reaction JSON, and I didn't want you inventing/assuming
  phrasing of that rule without checking my notes. If we add these, ask me
  first and flag it as an extrapolation beyond the JSON.
- The oxidative cleavage ketone-vs-acid substitution-pattern rule (already
  implemented) is a reasonable but *slightly* extrapolated reading of the
  JSON's generic "carbonyl compound(s) and/or carboxylic acid(s)" line — I
  approved it once already, but flag it again if you touch that logic.
- Inorganic/combustion reactions (from the alkanes topic — NO/NO2/CO2/etc.)
  are NOT procedurally generated, since there's no molecule structure to
  vary. They're simple fixed equations; keep them as a separate static set
  if/when we wire the full topic list in, not part of the generator.
- **Answer input UI is unresolved.** I built and explicitly rejected a
  click-to-build carbon-chain editor (see `prior_outputs/reaction-
  memoriser.html` for the rejected attempt — it's there for historical
  context, not to be reused as-is) because it can't scale to rings and
  larger molecules cleanly. We have NOT decided how I'll actually submit
  answers in the final app. Don't assume a UI approach — ask me, or propose
  options, before building one.
- No game loop / app shell wiring the generator up yet (weighted
  spaced-repetition-style selection toward reactions I get wrong existed
  conceptually in the earlier rejected HTML prototype's storage logic, but
  hasn't been reconnected to this engine).
- Currently only covers Alkanes + Alkenes topics (see
  `reference_data/*.json`, my actual syllabus reaction lists, kept
  unedited). My course is also currently covering halogen derivatives and
  arenes — extending the engine to those topics using the same
  operator-based pattern is a reasonable next step, but needs the same
  "ask before assuming a rule not in the source" discipline.

## Folder contents

- `engine/` — all engine code described above, runnable with plain Node
  (no dependencies).
- `reference_data/` — my original, unedited alkanes/alkenes reaction JSON
  (source of truth for what reactions exist and their exact
  reagents/conditions/mechanism text).
- `prior_outputs/` — artifacts from earlier in this project for context:
  a hand-written complex-problem set (superseded by the generator, but
  useful as a sanity-check reference for expected answer style), a sample
  generated batch, and the rejected HTML prototype (UI approach only,
  don't reuse the builder).

## What I want you to do first

Don't start building immediately. Read through `engine/` and run
`node test.js` inside it to confirm everything still passes in this new
environment. Then let's talk through: (1) how I want to submit answers now
that the click-builder is rejected, and (2) how to wire the engine into an
actual playable app (difficulty slider, weighted practice on my weak spots,
persistence). Ask me questions rather than assuming.
