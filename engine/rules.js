'use strict';
/* =========================================================================
   RULE TABLE
   Each rule: { id, topic, fgRequirement(fact)->bool, reagentRequirement
   (partial ReagentCondition, see reagents.js), transform(mol, matchedFacts)
   -> operator result }. resolve.js (the generic engine) is the only
   consumer of this file's shape; it never special-cases a rule id.

   Two-reactant rules (see that section, near the end of this file) are
   shaped differently -- `{ type:'two-reactant', roleA:{fgRequirement},
   roleB:{fgRequirement}, reagentRequirement, transform(first, second) }`
   -- and are resolved by resolve.js's separate resolveTwoReactant(),
   never mixed into the single-molecule candidate pool.

   Transforms are simply the SAME operators.js functions used by the old
   POOL system -- this table only changes how a reaction gets SELECTED,
   not how it's carried out once selected.
========================================================================= */
const op = require('./operators');
const { neighborsOf } = require('./engine');
const { scanFunctionalGroups } = require('./facts');

const RULES = [];

/* ---------------------------------------------------------------- Alkanes */
RULES.push({
  id:'radicalSubstitutionCl', topic:'Alkanes',
  fgRequirement: f=>f.kind==='alkylCH',
  reagentRequirement: { reagents:['Cl2'], light:true },
  transform: mol=>op.radicalSubstitution(mol, 'Cl'),
});
RULES.push({
  id:'radicalSubstitutionBr', topic:'Alkanes',
  fgRequirement: f=>f.kind==='alkylCH',
  reagentRequirement: { reagents:['Br2'], light:true },
  transform: mol=>op.radicalSubstitution(mol, 'Br'),
});
RULES.push({
  id:'combustion', topic:'Alkanes/Hydroxy',
  fgRequirement: f=>f.kind==='organic',
  reagentRequirement: { reagents:['O2'], quantity:'excess' },
  transform: mol=>op.combustion(mol),
});

/* ---------------------------------------------------------------- Alkenes */
RULES.push({
  id:'alkeneX2Addition', topic:'Alkenes',
  fgRequirement: f=>f.kind==='alkene',
  reagentRequirement: { reagents:['Br2'], solvent:'CCl4' },
  transform: mol=>op.additionX2(mol, 'Br'),
});
RULES.push({
  id:'alkeneCl2Addition', topic:'Alkenes',
  fgRequirement: f=>f.kind==='alkene',
  reagentRequirement: { reagents:['Cl2'] }, // ionic addition -- fast at room temp regardless of light; UV is only ever a distinguishing factor for a SEPARATE alkylCH site (radical substitution elsewhere in the same molecule), not for whether this addition itself occurs
  transform: mol=>op.additionX2(mol, 'Cl'),
});
RULES.push({
  id:'alkeneHalohydrin', topic:'Alkenes',
  fgRequirement: f=>f.kind==='alkene',
  reagentRequirement: { reagents:['Br2'], solvent:'aq' },
  transform: mol=>op.halohydrinFormation(mol, 'Br'),
});
RULES.push({
  id:'alkeneHydrogenation', topic:'Alkenes',
  fgRequirement: f=>f.kind==='alkene',
  reagentRequirement: { reagents:['H2'], catalyst:'Ni' },
  transform: mol=>op.hydrogenation(mol),
});
RULES.push({
  id:'alkeneMildOxidation', topic:'Alkenes',
  fgRequirement: f=>f.kind==='alkene',
  reagentRequirement: { reagents:['KMnO4'], temperature:{max:'cold'}, concentration:'dilute' },
  transform: mol=>op.mildOxidationDiol(mol),
});
RULES.push({
  id:'alkeneOxidativeCleavage', topic:'Alkenes',
  fgRequirement: f=>f.kind==='alkene',
  reagentRequirement: { reagents:['KMnO4'], temperature:{min:'heat'}, concentration:'concentrated', pH:'acidic' },
  transform: mol=>op.oxidativeCleavage(mol),
});
RULES.push({
  id:'alkeneHBr', topic:'Alkenes',
  fgRequirement: f=>f.kind==='alkene',
  reagentRequirement: { reagents:['HBr'] },
  transform: mol=>op.hydrohalogenation(mol, 'Br'),
});
RULES.push({
  id:'alkeneHCl', topic:'Alkenes',
  fgRequirement: f=>f.kind==='alkene',
  reagentRequirement: { reagents:['HCl'] },
  transform: mol=>op.hydrohalogenation(mol, 'Cl'),
});
RULES.push({
  id:'alkeneHydration', topic:'Alkenes',
  fgRequirement: f=>f.kind==='alkene',
  reagentRequirement: { reagents:['H2O'], temperature:{min:'heat'} },
  transform: mol=>op.hydration(mol),
});

/* ----------------------------------------------------------------- Arenes
   The three catalysed EAS rules below all require `!hasActivatingSub` --
   an activated ring (phenol/aniline) is reactive enough that it never
   needs a Lewis-acid catalyst (real chemistry: phenol brominates readily
   with plain aqueous Br2, no FeBr3), so those get their own,
   catalyst-free rules further down. Without this exclusion, a test
   condition on an activated ring that happened to also carry a catalyst
   tag would satisfy BOTH families at once -- a real clash, not a
   legitimate double reaction, since they're the same site.
*/
RULES.push({
  id:'areneChlorination', topic:'Arenes',
  fgRequirement: f=>f.kind==='arene' && f.freePositions.length>0 && !f.hasActivatingSub,
  reagentRequirement: { reagents:['Cl2'], catalyst:'AlCl3' },
  transform: mol=>op.ringElectrophilicSubstitution(mol, 'Cl'),
});
RULES.push({
  id:'areneBromination', topic:'Arenes',
  fgRequirement: f=>f.kind==='arene' && f.freePositions.length>0 && !f.hasActivatingSub,
  reagentRequirement: { reagents:['Br2'], catalyst:'FeBr3' },
  transform: mol=>op.ringElectrophilicSubstitution(mol, 'Br'),
});
RULES.push({
  id:'areneNitration', topic:'Arenes',
  fgRequirement: f=>f.kind==='arene' && f.freePositions.length>0 && !f.hasActivatingSub,
  reagentRequirement: { reagents:['HNO3','H2SO4'] },
  transform: mol=>op.ringElectrophilicSubstitution(mol, 'NO2'),
});
RULES.push({
  id:'activatedRingMonobromination', topic:'Arenes/Hydroxy/Nitrogen',
  fgRequirement: f=>f.kind==='arene' && f.hasActivatingSub && f.freePositions.length>0,
  reagentRequirement: { reagents:['Br2'], solvent:'CCl4', quantity:'controlled' },
  transform: mol=>op.ringElectrophilicSubstitution(mol, 'Br'),
});
RULES.push({
  id:'activatedRingTribromination', topic:'Arenes/Hydroxy/Nitrogen',
  fgRequirement: f=>f.kind==='arene' && f.hasActivatingSub,
  reagentRequirement: { reagents:['Br2'], solvent:'aq', quantity:'excess' },
  transform: (mol, facts)=>{
    const ring = mol.nodes.find(n=>n.id===facts[0].ringId);
    const sub = ring.subs.find(s=>s==='OH'||s==='NH2');
    return op.ringTribromination(mol, sub);
  },
});
RULES.push({
  id:'phenolTrinitration', topic:'Hydroxy Compounds',
  fgRequirement: f=>f.kind==='arene' && f.hasActivatingSub,
  reagentRequirement: { reagents:['HNO3','H2SO4'], quantity:'excess', temperature:{min:'heat'} },
  transform: (mol, facts)=>{
    const ring = mol.nodes.find(n=>n.id===facts[0].ringId);
    const sub = ring.subs.find(s=>s==='OH');
    if(!sub) return { occurs:false };
    return op.ringTrinitration(mol, sub);
  },
});
RULES.push({
  id:'benzylicOxidation', topic:'Arenes',
  fgRequirement: f=>f.kind==='benzylicCH',
  reagentRequirement: { reagents:['KMnO4'], pH:'acidic', temperature:{min:'heat'} },
  transform: mol=>op.sideChainOxidationToBenzoicAcid(mol),
});

/* ------------------------------------------------------- Halogen Derivatives */
const ALIPHATIC_HALS = new Set(['Cl','Br','I']); // F excluded -- C-F bond too strong to hydrolyse/substitute at this level
RULES.push({
  id:'haloalkaneHydrolysis', topic:'Halogen Derivatives',
  fgRequirement: f=>f.kind==='haloC' && !f.aromatic && ALIPHATIC_HALS.has(f.halogen),
  reagentRequirement: { reagents:['NaOH'], solvent:'aq', temperature:{min:'heat'} },
  transform: (mol, facts)=>op.nucleophilicSubstitutionFlat(mol, facts[0].halogen, 'OH'),
});
RULES.push({
  id:'haloalkaneToNitrile', topic:'Halogen Derivatives',
  fgRequirement: f=>f.kind==='haloC' && !f.aromatic && ALIPHATIC_HALS.has(f.halogen),
  reagentRequirement: { reagents:['KCN'], solvent:'ethanolic', temperature:{min:'heat'} },
  transform: (mol, facts)=>op.nitrileFormation(mol, facts[0].halogen),
});
RULES.push({
  id:'haloalkaneToPrimaryAmine', topic:'Halogen Derivatives',
  fgRequirement: f=>f.kind==='haloC' && !f.aromatic && ALIPHATIC_HALS.has(f.halogen),
  reagentRequirement: { reagents:['NH3'], solvent:'ethanolic', temperature:{min:'heat'}, quantity:'excess' },
  transform: (mol, facts)=>op.nucleophilicSubstitutionFlat(mol, facts[0].halogen, 'NH2'),
});
RULES.push({
  id:'haloalkaneElimination', topic:'Halogen Derivatives',
  fgRequirement: f=>f.kind==='haloC' && !f.aromatic && ALIPHATIC_HALS.has(f.halogen),
  reagentRequirement: { reagents:['NaOH'], solvent:'ethanolic', temperature:{min:'heat'} },
  transform: (mol, facts)=>op.eliminationHX(mol, facts[0].halogen),
});

/* ------------------------------------------------------------------ Amines */
RULES.push({
  id:'aminePlusAcid', topic:'Amines/Nitrogen Compounds',
  fgRequirement: f=>f.kind==='amine',
  // NaNO2 excluded: a primary aromatic amine with NaNO2+HCl present is
  // the diazotisation rule's job (below), not simple protonation, even
  // though HCl alone would be enough to protonate it in isolation.
  reagentRequirement: { reagents:['HCl'], reagentsExclude:['NaNO2'] },
  transform: mol=>op.protonateAmine(mol),
});

/* ----------------------------------------------------------- Hydroxy Compounds */
RULES.push({
  id:'alcoholPlusNa', topic:'Hydroxy Compounds',
  fgRequirement: f=>f.kind==='hydroxyl' && !f.aromatic,
  reagentRequirement: { reagents:['Na'] },
  transform: mol=>op.nucleophilicSubstitutionFlat(mol, 'OH', 'ONa'),
});
RULES.push({
  id:'phenolPlusNa', topic:'Hydroxy Compounds',
  fgRequirement: f=>f.kind==='hydroxyl' && f.aromatic,
  reagentRequirement: { reagents:['Na'] },
  transform: mol=>op.ringFlatSubSwap(mol, 'OH', 'ONa'),
});
RULES.push({
  id:'phenolPlusNaOH', topic:'Hydroxy Compounds',
  fgRequirement: f=>f.kind==='hydroxyl' && f.aromatic,
  reagentRequirement: { reagents:['NaOH'], solvent:'aq' },
  transform: mol=>op.ringFlatSubSwap(mol, 'OH', 'ONa'),
});
RULES.push({
  id:'dehydrationOfAlcohol', topic:'Hydroxy Compounds',
  fgRequirement: f=>f.kind==='hydroxyl' && !f.aromatic,
  reagentRequirement: { reagents:['H2SO4'], concentration:'concentrated', temperature:{min:'heat'} },
  transform: mol=>op.eliminationHX(mol, 'OH'),
});
RULES.push({
  id:'oxidationPrimaryAlcoholToAldehyde', topic:'Hydroxy Compounds',
  fgRequirement: f=>f.kind==='hydroxyl' && !f.aromatic && f.class==='1°',
  reagentRequirement: { reagents:['K2Cr2O7'], pH:'acidic', quantity:'controlled' },
  transform: mol=>op.oxidizeAlcohol(mol, 'CHO'),
});
RULES.push({
  id:'oxidationPrimaryAlcoholToAcid', topic:'Hydroxy Compounds',
  fgRequirement: f=>f.kind==='hydroxyl' && !f.aromatic && f.class==='1°',
  reagentRequirement: { reagents:['K2Cr2O7'], pH:'acidic', quantity:'excess', temperature:{min:'heat'} },
  transform: mol=>op.oxidizeAlcohol(mol, 'COOH'),
});
RULES.push({
  id:'oxidationSecondaryAlcoholToKetone', topic:'Hydroxy Compounds',
  fgRequirement: f=>f.kind==='hydroxyl' && !f.aromatic && f.class==='2°',
  reagentRequirement: { reagents:['K2Cr2O7'], pH:'acidic' },
  transform: mol=>op.oxidizeAlcohol(mol, 'ketone'),
});
RULES.push({
  id:'alcoholToBromoalkane', topic:'Hydroxy Compounds',
  fgRequirement: f=>f.kind==='hydroxyl' && !f.aromatic,
  reagentRequirement: { reagents:['HBr'] },
  transform: mol=>op.nucleophilicSubstitutionFlat(mol, 'OH', 'Br'),
});
RULES.push({
  id:'nitrationOfPhenol', topic:'Hydroxy Compounds',
  fgRequirement: f=>f.kind==='arene' && f.hasActivatingSub && f.freePositions.length>0,
  reagentRequirement: { reagents:['HNO3'], concentration:'dilute', temperature:{max:'room'} },
  transform: mol=>op.ringElectrophilicSubstitution(mol, 'NO2'),
});
RULES.push({
  id:'iodoformMethylKetone', topic:'Hydroxy/Carbonyl Compounds',
  fgRequirement: f=>f.kind==='carbonyl' && f.class==='ketone' && f.hasAlphaMethyl,
  reagentRequirement: { reagents:['I2','NaOH'], temperature:{min:'warm'} },
  transform: mol=>op.iodoformCleavage(mol),
});
// Alcohol-side iodoform (CH3-CH(OH)-R): the fgRequirement below is
// deliberately broad (any non-tertiary aliphatic alcohol) -- the real
// "does it have a methyl neighbour on the carbinol carbon" shape check
// already lives inside iodoformCleavageFromAlcohol itself, which
// correctly returns occurs:false (silently dropped by resolve.js) for
// e.g. propan-1-ol, matching how every other operator here is the actual
// source of truth for its own shape, not the fact scanner.
RULES.push({
  id:'iodoformAlcohol', topic:'Hydroxy Compounds',
  fgRequirement: f=>f.kind==='hydroxyl' && !f.aromatic && f.class!=='3°',
  reagentRequirement: { reagents:['I2','NaOH'], temperature:{min:'warm'} },
  transform: mol=>op.iodoformCleavageFromAlcohol(mol),
});

/* --------------------------------------------------------- Carbonyl Compounds */
RULES.push({
  id:'hcnAddition', topic:'Carbonyl Compounds',
  fgRequirement: f=>f.kind==='carbonyl',
  reagentRequirement: { reagents:['HCN'] },
  transform: mol=>op.cyanohydrinFormation(mol),
});
RULES.push({
  id:'carbonylReductionNaBH4', topic:'Carbonyl Compounds',
  fgRequirement: f=>f.kind==='carbonyl',
  reagentRequirement: { reagents:['NaBH4'] },
  transform: mol=>op.reduceCarbonyl(mol),
});
RULES.push({
  id:'carbonylReductionLiAlH4', topic:'Carbonyl Compounds',
  fgRequirement: f=>f.kind==='carbonyl',
  reagentRequirement: { reagents:['LiAlH4'], solvent:'dry ether' },
  transform: mol=>op.reduceCarbonyl(mol),
});
RULES.push({
  id:'ethanalIodoform', topic:'Carbonyl Compounds',
  fgRequirement: f=>f.kind==='carbonyl' && f.class==='aldehyde' && !f.aromatic,
  reagentRequirement: { reagents:['I2','NaOH'], temperature:{min:'warm'} },
  transform: mol=>op.iodoformCleavage(mol),
});
RULES.push({
  id:'oxidationAldehydeToAcid', topic:'Carbonyl Compounds',
  fgRequirement: f=>f.kind==='carbonyl' && f.class==='aldehyde',
  reagentRequirement: { reagents:['K2Cr2O7'], pH:'acidic' },
  transform: mol=>op.oxidizeAldehyde(mol),
});

/* --------------------------------------------------- Carboxylic Acids & Derivatives */
RULES.push({
  id:'carboxylicAcidPlusBase', topic:'Carboxylic Acids & Derivatives',
  fgRequirement: f=>f.kind==='carboxyl' && f.group==='COOH',
  reagentRequirement: { reagents:['NaOH'] },
  transform: mol=>op.carboxylicAcidSaltFormation(mol),
});
// Methanoic acid (still has an aldehyde-like H on its carbonyl carbon)
// and ethanedioic acid (two carboxyls on a weak C-C bond) are the two
// carboxylic acids hot acidified KMnO4 oxidises all the way to CO2 --
// every other carboxylic acid is already at its most oxidised state and
// this rule correctly never fires for them (specialOxidisable is only
// ever true for those two shapes, see facts.js).
RULES.push({
  id:'specialAcidOxidation', topic:'Carboxylic Acids & Derivatives',
  fgRequirement: f=>f.kind==='carboxyl' && f.group==='COOH' && f.specialOxidisable,
  reagentRequirement: { reagents:['KMnO4'], pH:'acidic', temperature:{min:'heat'} },
  transform: mol=>op.combustion(mol),
});
RULES.push({
  id:'acidToAcylChloride', topic:'Carboxylic Acids & Derivatives',
  fgRequirement: f=>f.kind==='carboxyl' && f.group==='COOH',
  reagentRequirement: { reagents:['SOCl2'] },
  transform: mol=>op.acidToAcylChloride(mol),
});
RULES.push({
  id:'acidReductionLiAlH4', topic:'Carboxylic Acids & Derivatives',
  fgRequirement: f=>f.kind==='carboxyl' && f.group==='COOH',
  reagentRequirement: { reagents:['LiAlH4'], solvent:'dry ether' },
  transform: mol=>op.reduceCarboxylicAcid(mol),
});
RULES.push({
  id:'acylChlorideHydrolysis', topic:'Carboxylic Acids & Derivatives',
  fgRequirement: f=>f.kind==='acylHalide',
  reagentRequirement: { reagents:['H2O'] },
  transform: mol=>op.hydrolyzeAcylChloride(mol),
});
RULES.push({
  id:'esterHydrolysisAcidic', topic:'Carboxylic Acids & Derivatives',
  fgRequirement: f=>f.kind==='ester',
  reagentRequirement: { reagents:['H2O'], pH:'acidic', temperature:{min:'heat'} },
  transform: mol=>op.hydrolyzeEster(mol, 'COOH'),
});
RULES.push({
  id:'esterHydrolysisAlkaline', topic:'Carboxylic Acids & Derivatives',
  fgRequirement: f=>f.kind==='ester',
  reagentRequirement: { reagents:['NaOH'], solvent:'aq', temperature:{min:'heat'} },
  transform: mol=>op.hydrolyzeEster(mol, 'COONa'),
});
RULES.push({
  id:'anhydrideHydrolysis', topic:'Carboxylic Acids & Derivatives',
  fgRequirement: f=>f.kind==='anhydride',
  reagentRequirement: { reagents:['H2O'] },
  transform: mol=>op.hydrolyzeAcidAnhydride(mol),
});
RULES.push({
  id:'primaryAmideHydrolysisAcidic', topic:'Carboxylic Acids & Derivatives',
  fgRequirement: f=>f.kind==='amide' && f.class==='primary',
  reagentRequirement: { reagents:['H2O'], pH:'acidic', temperature:{min:'heat'} },
  transform: mol=>op.hydrolyzeAmide(mol, 'COOH'),
});
RULES.push({
  id:'primaryAmideHydrolysisAlkaline', topic:'Carboxylic Acids & Derivatives',
  fgRequirement: f=>f.kind==='amide' && f.class==='primary',
  reagentRequirement: { reagents:['NaOH'], solvent:'aq', temperature:{min:'heat'} },
  transform: mol=>op.hydrolyzeAmide(mol, 'COONa'),
});
RULES.push({
  id:'amideReduction', topic:'Carboxylic Acids & Derivatives',
  fgRequirement: f=>f.kind==='amide' && f.class==='primary',
  reagentRequirement: { reagents:['LiAlH4'], solvent:'dry ether' },
  transform: mol=>op.reduceAmide(mol),
});
RULES.push({
  id:'nitrileHydrolysisAcidic', topic:'Halogen Derivatives',
  fgRequirement: f=>f.kind==='nitrile',
  reagentRequirement: { reagents:['H2O'], pH:'acidic', temperature:{min:'heat'} },
  transform: mol=>op.nitrileHydrolysis(mol, 'COOH'),
});
RULES.push({
  id:'nitrileHydrolysisAlkaline', topic:'Halogen Derivatives',
  fgRequirement: f=>f.kind==='nitrile',
  reagentRequirement: { reagents:['NaOH'], solvent:'aq', temperature:{min:'heat'} },
  transform: mol=>op.nitrileHydrolysis(mol, 'COONa'),
});
RULES.push({
  id:'nitrileReduction', topic:'Halogen Derivatives',
  fgRequirement: f=>f.kind==='nitrile',
  reagentRequirement: { reagents:['LiAlH4'], solvent:'dry ether' },
  transform: mol=>op.nitrileReduction(mol),
});

/* -------------------------------------------------------- Nitrogen Compounds */
RULES.push({
  id:'diazotisation', topic:'Nitrogen Compounds',
  fgRequirement: f=>f.kind==='amine' && f.aromatic,
  reagentRequirement: { reagents:['NaNO2','HCl'], temperature:{max:'cold'} },
  transform: mol=>op.diazotisation(mol),
});
RULES.push({
  id:'nitrobenzeneReduction', topic:'Nitrogen Compounds',
  fgRequirement: f=>f.kind==='nitro' && f.aromatic,
  reagentRequirement: { reagents:['Sn','HCl'], temperature:{min:'heat'} },
  transform: mol=>op.ringFlatSubSwap(mol, 'NO2', 'NH2'),
});
RULES.push({
  id:'substitutedAmideHydrolysisAcidic', topic:'Nitrogen Compounds', // covers peptide-bond hydrolysis too -- a peptide bond IS a secondary amide, no separate rule needed
  fgRequirement: f=>f.kind==='amide' && f.class==='secondary',
  reagentRequirement: { reagents:['H2O'], pH:'acidic', temperature:{min:'heat'} },
  transform: mol=>op.hydrolyzePeptideBond(mol, 'acidic'),
});
RULES.push({
  id:'substitutedAmideHydrolysisAlkaline', topic:'Nitrogen Compounds',
  fgRequirement: f=>f.kind==='amide' && f.class==='secondary',
  reagentRequirement: { reagents:['NaOH'], solvent:'aq', temperature:{min:'heat'} },
  transform: mol=>op.hydrolyzePeptideBond(mol, 'alkaline'),
});

/* =========================================================================
   TWO-REACTANT RULES
   Resolved by resolve.js's resolveTwoReactant(), never mixed into the
   single-molecule RULES pool above. roleA/roleB are role labels only
   (which molecule ends up as the operator's FIRST vs SECOND argument) --
   resolveTwoReactant tries both physical orderings of the two given
   molecules against them, so the caller never has to pre-decide which
   one is "the acid" etc.
========================================================================= */

// Walks a plain unbranched chain from a given start node, counting every
// reachable carbon -- bridges this engine's real molecule graph to
// ringFriedelCraftsAlkylation's own chainLen PARAMETER (the operator
// itself was written for the old POOL system, which passed chain length
// as a hand-carried integer alongside the reactant molecule rather than
// deriving it from the molecule's own structure; this derives the same
// number from the actual haloalkane molecule instead).
function chainLengthFrom(mol, startId){
  const visited = new Set(); let count = 0;
  (function dfs(id){
    if(visited.has(id)) return;
    visited.add(id); count++;
    neighborsOf(mol, id).forEach(nb=>{ if(!visited.has(nb.to)) dfs(nb.to); });
  })(startId);
  return count;
}

RULES.push({
  id:'esterificationAlcoholAcid', topic:'Hydroxy Compounds', type:'two-reactant',
  roleA:{ fgRequirement: f=>f.kind==='carboxyl' && f.group==='COOH' },
  roleB:{ fgRequirement: f=>f.kind==='hydroxyl' && !f.aromatic },
  reagentRequirement: { catalyst:'H2SO4', concentration:'concentrated', temperature:{min:'heat'} },
  transform: (acid, alcohol)=>op.esterifyAcid(acid, alcohol),
});
RULES.push({
  id:'esterificationAcylChloride', topic:'Carboxylic Acids & Derivatives', type:'two-reactant',
  roleA:{ fgRequirement: f=>f.kind==='acylHalide' },
  roleB:{ fgRequirement: f=>f.kind==='hydroxyl' }, // accepts BOTH an aliphatic alcohol and a phenol -- unlike esterifyAcid above, esterifyAcylChloride doesn't refuse a ring
  reagentRequirement: {},
  transform: (acyl, nucleophile)=>op.esterifyAcylChloride(acyl, nucleophile),
});
RULES.push({
  id:'friedelCraftsAlkylation', topic:'Arenes', type:'two-reactant',
  roleA:{ fgRequirement: f=>f.kind==='arene' && f.freePositions.length>0 && !f.hasActivatingSub },
  roleB:{ fgRequirement: f=>f.kind==='haloC' && !f.aromatic },
  reagentRequirement: { catalyst:'AlCl3' },
  transform: (arene, haloalkane)=>{
    const hf = scanFunctionalGroups(haloalkane).find(f=>f.kind==='haloC');
    if(!hf) return { occurs:false };
    return op.ringFriedelCraftsAlkylation(arene, chainLengthFrom(haloalkane, hf.nodeId));
  },
});
RULES.push({
  id:'acylChloridePlusAmine', topic:'Carboxylic Acids & Derivatives', type:'two-reactant',
  roleA:{ fgRequirement: f=>f.kind==='acylHalide' },
  roleB:{ fgRequirement: f=>f.kind==='amine' },
  reagentRequirement: {},
  transform: (acyl, amine)=>op.acylChlorideToAmide(acyl, amine),
});
RULES.push({
  id:'azoCoupling', topic:'Nitrogen Compounds', type:'two-reactant',
  roleA:{ fgRequirement: f=>f.kind==='diazonium' },
  roleB:{ fgRequirement: f=>f.kind==='arene' && f.hasActivatingSub },
  reagentRequirement: { temperature:{max:'cold'} },
  transform: (diazonium, arene)=>op.azoCoupling(diazonium, arene),
});
RULES.push({
  id:'peptideBondFormation', topic:'Nitrogen Compounds', type:'two-reactant',
  roleA:{ fgRequirement: f=>f.kind==='aminoAcid' },
  roleB:{ fgRequirement: f=>f.kind==='aminoAcid' }, // symmetric -- either amino acid can supply either end, see resolveTwoReactant's own header for how that becomes `variants` (Gly-Ala vs Ala-Gly) rather than a clash
  reagentRequirement: {},
  transform: (acid1, acid2)=>op.formPeptideBond(acid1, acid2),
});

module.exports = { RULES };
