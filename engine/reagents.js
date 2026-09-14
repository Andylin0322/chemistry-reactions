'use strict';
/* =========================================================================
   REAGENT / CONDITION FACTS
   A ReagentCondition is the structured shape every rule's precondition is
   checked against -- built from a controlled tag vocabulary (never
   freeform text; see TAG_VOCAB below), so both the question generator and
   a student's own tag-based answer produce the exact same shape.

     { reagents:['Br2'], solvent:'aq', catalyst:null, temperature:'room',
       light:false, quantity:'controlled', pH:null }

   `reagents` is always present (possibly empty); every other field is
   either a real value or `null` (explicitly absent/not applicable -- e.g.
   `catalyst:null` on a plain Br2(aq) bromination). There is deliberately
   no "unknown" state for a well-formed condition -- callers building one
   from a tag list (see conditionFromTags) fill every field in.

   A RULE's requirement (see rules.js) is a PARTIAL condition: any field
   it omits means "don't care", any field it sets means "must be
   satisfied" -- `reagents` by subset (extra reagents in the given
   condition never block a match), everything else by the match rules in
   `satisfies` below. This is what gives "minimal satisfaction" its
   precise meaning: a rule fires on the SMALLEST condition that still
   contains everything it asked for, regardless of what else is present.
========================================================================= */

// Ordered so a rule can require a MINIMUM ("needs at least heat") or a
// MAXIMUM ("must stay cold or it decomposes") via {min} / {max} rather
// than an exact string, without hand-listing every value on each side.
// Deliberately no separate 'reflux' level above 'heat' -- reflux is lab
// technique (how the heating is done so volatiles don't escape), not a
// hotter condition or its own markable point in H2 exams, so "heat" and
// "heat under reflux" are the exact same condition here, always (see
// app.html's own tagSetsEqual/normalizeConditionTag for the matching
// equivalence on the display/grading side).
const TEMP_ORDER = ['cold', 'room', 'warm', 'heat'];
function tempIndex(t){ const i = TEMP_ORDER.indexOf(t); if(i===-1) throw new Error('unknown temperature: '+t); return i; }

// The full controlled vocabulary a ReagentCondition/requirement is built
// from -- deliberately small and closed (no freeform reagent text ever
// reaches the engine). Extend this list, not ad-hoc strings, when a new
// rule needs a reagent/catalyst/solvent this doesn't already cover.
const REAGENT_SPECIES = [
  'Cl2','Br2','I2','HBr','HCl','HI','H2','H2O','O2','KMnO4','K2Cr2O7',
  'NaOH','Na','Na2CO3','NaHCO3','HNO3','H2SO4','PCl5','PCl3','SOCl2',
  'LiAlH4','NaBH4','HCN','KCN','NH3','NaNO2','Sn','I2/NaOH',
  'NH2R-primary','NH2R-secondary','R-COCl','R-COOH','R-OH','ArOH','ArNH2',
];
const SOLVENTS = ['aq', 'CCl4', 'ethanolic', 'dry ether', 'aqueous KOH'];
const CATALYSTS = ['FeBr3', 'AlCl3', 'Ni', 'Pt', 'V2O5', 'H2SO4']; // conc H2SO4 acts catalytically in Fischer esterification, not as a stoichiometric reagent -- modelled here, not in `reagents`
const QUANTITIES = ['controlled', 'excess'];
const PHS = ['acidic', 'alkaline'];
const CONCENTRATIONS = ['dilute', 'concentrated']; // solution STRENGTH (dilute/hot-concentrated KMnO4) -- distinct from `quantity`, which is stoichiometric AMOUNT (controlled/excess Br2 water)

function emptyCondition(){
  return { reagents:[], solvent:null, catalyst:null, temperature:'room', light:false, quantity:null, pH:null, concentration:null };
}

// `given` is a full ReagentCondition; `req` is a PARTIAL one (a rule's
// requirement) -- every key req sets must be satisfied, every key it
// omits is ignored. `req.temperature` may be a plain string (exact) or
// `{min}`/`{max}`/`{min,max}` (range).
function satisfies(given, req){
  if(req.reagents && req.reagents.length){
    if(!req.reagents.every(r=>given.reagents.includes(r))) return false;
  }
  // Priority/exclusivity between two rules that would otherwise both
  // match the same site (e.g. a primary aromatic amine's -NH2 reacts
  // with NaNO2/HCl by diazotisation, not simple protonation, even
  // though HCl alone is also present) -- expressed as a requirement,
  // not a runtime tie-break, so resolve.js's clash detector still
  // throws if two rules are ever genuinely under-specified against each
  // other instead of silently picking a winner.
  if(req.reagentsExclude && req.reagentsExclude.some(r=>given.reagents.includes(r))) return false;
  if('solvent' in req && req.solvent !== given.solvent) return false;
  if('catalyst' in req && req.catalyst !== given.catalyst) return false;
  if('light' in req && req.light !== given.light) return false;
  if('quantity' in req && req.quantity !== given.quantity) return false;
  if('pH' in req && req.pH !== given.pH) return false;
  if('concentration' in req && req.concentration !== given.concentration) return false;
  if(req.temperature != null){
    const gi = tempIndex(given.temperature);
    if(typeof req.temperature === 'string'){
      if(gi !== tempIndex(req.temperature)) return false;
    } else {
      if(req.temperature.min != null && gi < tempIndex(req.temperature.min)) return false;
      if(req.temperature.max != null && gi > tempIndex(req.temperature.max)) return false;
    }
  }
  return true;
}

module.exports = { TEMP_ORDER, tempIndex, REAGENT_SPECIES, SOLVENTS, CATALYSTS, QUANTITIES, PHS, CONCENTRATIONS, emptyCondition, satisfies };
