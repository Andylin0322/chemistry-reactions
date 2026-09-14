'use strict';
/* =========================================================================
   RESOLVER
   The one function question generation AND all three grading modes call
   (see rules.js's own header): given a molecule and a structured
   ReagentCondition, scan its functional groups, find every rule whose
   precondition (functional-group shape + reagent condition) is
   satisfied, and apply each one's transform -- ONCE per matched rule,
   never once per fact (see facts.js's header for why). Two independent
   functional groups both reacting under the same reagent (e.g. an
   alkene AND an alcohol, both touched by the same reagent) is the
   expected, common case and is handled by simply chaining both
   transforms. Two DIFFERENT rules genuinely competing for the SAME site
   is a rule-design bug, not a real ambiguity to silently pick a winner
   for -- resolve() throws loudly so it's caught during authoring/vetting
   rather than silently mis-grading a live question.
========================================================================= */
const { scanFunctionalGroups } = require('./facts');
const { satisfies } = require('./reagents');
const { canonicalForm } = require('./engine');

class RuleClashError extends Error {
  constructor(fact, ruleIds){
    super(`Rule clash: [${ruleIds.join(', ')}] all matched the same ${fact.kind} site (fact #${fact.id}) under the same reagent condition -- exactly one rule should have matched. Fix the rules' reagentRequirement/fgRequirement specificity, this is never correct behaviour.`);
    this.fact = fact; this.ruleIds = ruleIds;
  }
}

function resolve(mol, reagentCondition, rules){
  const facts = scanFunctionalGroups(mol);
  const candidates = []; // {rule, facts: matching facts for this rule}
  const claimedBy = new Map(); // factId -> [ruleId, ...]

  for(const rule of rules){
    if(rule.type==='two-reactant') continue; // handled only by resolveTwoReactant, below
    const matchingFacts = facts.filter(f=>rule.fgRequirement(f));
    if(matchingFacts.length===0) continue;
    if(!satisfies(reagentCondition, rule.reagentRequirement || {})) continue;
    candidates.push({ rule, facts: matchingFacts });
    matchingFacts.forEach(f=>{
      const list = claimedBy.get(f.id) || [];
      list.push(rule.id);
      claimedBy.set(f.id, list);
    });
  }

  for(const [factId, ruleIds] of claimedBy){
    if(ruleIds.length>1){
      throw new RuleClashError(facts.find(f=>f.id===factId), [...new Set(ruleIds)]);
    }
  }

  if(candidates.length===0) return { occurs:false, facts, candidateRuleIds:[] };

  let mol2 = mol;
  const applied = [];
  for(const {rule, facts:mf} of candidates){
    const result = rule.transform(mol2, mf);
    // A rule's fgRequirement/reagentRequirement matching is a NECESSARY
    // but not always SUFFICIENT condition (some real chemistry shape
    // checks -- e.g. "does this ring already have 3 substituents" --
    // live inside the operator itself, not duplicated into the fact
    // scanner) -- such a rule is dropped here rather than treated as an
    // engine error, so an unrelated qualifying rule can still apply.
    if(!result || !result.occurs) continue;
    applied.push(rule.id);
    if(result.product !== undefined){ mol2 = result.product; continue; }
    if(result.products){ return { occurs:true, appliedRules:[...applied], products: result.products, facts }; }
    if(result.variants){ return { occurs:true, appliedRules:[...applied], variants: result.variants, facts }; }
  }
  if(applied.length===0) return { occurs:false, facts, candidateRuleIds: candidates.map(c=>c.rule.id) };
  return { occurs:true, appliedRules: applied, product: mol2, facts };
}

/* =========================================================================
   TWO-REACTANT RESOLVER
   For rules shaped `{ type:'two-reactant', roleA, roleB, reagentRequirement,
   transform(first, second) }` (esterification, Friedel-Crafts, acyl
   chloride + amine, azo coupling, peptide bond FORMATION -- see rules.js's
   "Two-reactant rules" section). `roleA`/`roleB` are each just
   `{fgRequirement}`, exactly like a single-reactant rule's own
   fgRequirement -- a two-reactant rule is really just a single-reactant
   rule with a second precondition on a second molecule.

   Order is NOT significant to the caller -- molA/molB are tried in BOTH
   role assignments, since a student's submitted pair (or a generator
   that hasn't committed to which slot is which) has no fixed order. Two
   assignments matching the SAME rule (e.g. two amino acids, each of
   which independently has both a free -NH2 and a free -COOH, so either
   one can supply either role) is NOT a clash -- it's two genuinely
   different, both-correct products (Gly-Ala vs Ala-Gly), collected as
   `variants`. Two assignments matching DIFFERENT rules IS a clash, same
   policy as the single-reactant resolver and for the same reason.
========================================================================= */
class TwoReactantRuleClashError extends Error {
  constructor(ruleIds){
    super(`Two-reactant rule clash: [${ruleIds.join(', ')}] all matched this pair of molecules under the same reagent condition -- fix the rules' role/reagent requirement specificity, this is never correct behaviour.`);
    this.ruleIds = ruleIds;
  }
}

function resolveTwoReactant(molA, molB, reagentCondition, rules){
  const twoReactantRules = rules.filter(r=>r.type==='two-reactant');
  const factsA = scanFunctionalGroups(molA);
  const factsB = scanFunctionalGroups(molB);

  const assignments = []; // {rule, first, second, firstFacts, secondFacts}
  for(const rule of twoReactantRules){
    if(!satisfies(reagentCondition, rule.reagentRequirement || {})) continue;
    const aInRoleA = factsA.filter(f=>rule.roleA.fgRequirement(f));
    const bInRoleB = factsB.filter(f=>rule.roleB.fgRequirement(f));
    if(aInRoleA.length && bInRoleB.length){
      assignments.push({ rule, first:molA, second:molB, firstFacts:aInRoleA, secondFacts:bInRoleB });
    }
    const bInRoleA = factsB.filter(f=>rule.roleA.fgRequirement(f));
    const aInRoleB = factsA.filter(f=>rule.roleB.fgRequirement(f));
    if(bInRoleA.length && aInRoleB.length){
      assignments.push({ rule, first:molB, second:molA, firstFacts:bInRoleA, secondFacts:aInRoleB });
    }
  }

  if(assignments.length===0) return { occurs:false };

  const distinctRuleIds = [...new Set(assignments.map(a=>a.rule.id))];
  if(distinctRuleIds.length>1) throw new TwoReactantRuleClashError(distinctRuleIds);

  const rule = assignments[0].rule;
  const seen = new Set();
  const variants = [];
  let smallResult = null;
  for(const {first, second} of assignments){
    const result = rule.transform(first, second);
    if(!result || !result.occurs) continue;
    if(result.products || result.variants){ smallResult = result; break; } // fragment/variant-producing transform -- return as-is, no cross-assignment dedup needed
    if(result.product){
      const key = canonicalForm(result.product);
      if(!seen.has(key)){ seen.add(key); variants.push(result.product); }
    }
  }
  if(smallResult) return { occurs:true, appliedRules:[rule.id], ...(smallResult.products ? {products:smallResult.products} : {variants:smallResult.variants}) };
  if(variants.length===0) return { occurs:false };
  return variants.length===1
    ? { occurs:true, appliedRules:[rule.id], product: variants[0] }
    : { occurs:true, appliedRules:[rule.id], variants };
}

module.exports = { resolve, RuleClashError, resolveTwoReactant, TwoReactantRuleClashError };
