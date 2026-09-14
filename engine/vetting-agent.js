'use strict';
/* =========================================================================
   VETTING AGENT
   An independent correctness check on resolve.js + rules.js + facts.js,
   NOT on the individual operators (those already have their own coverage
   in engine/test.js). Every expected answer below was worked out by hand
   from A-level (H2) organic chemistry knowledge, not derived from the
   engine itself -- the point is to catch cases where the RULE TABLE
   mis-selects, over-fires, under-fires, or throws a clash, which
   engine/test.js's direct operator calls can never exercise (they never
   go through the "which reaction, if any, applies here" decision at
   all).

   Deliberately includes: molecules with 2-3 independent functional
   groups (only the intended one(s) should react), functional groups that
   look reactive but structurally aren't (tertiary alcohol, ketone vs
   aldehyde, fluoroalkane, an unactivated ring without a catalyst),
   conditions engineered to distinguish between two real reactions that
   share most of their reagents (acid vs alkaline hydrolysis, dilute-cold
   vs hot-concentrated-acidified KMnO4, HCl alone vs NaNO2/HCl on an
   aromatic amine), and a few full multi-site "everything that can react,
   does" cases.

   Run: node engine/vetting-agent.js
========================================================================= */
const { resolve, RuleClashError, resolveTwoReactant, TwoReactantRuleClashError } = require('./resolve');
const { RULES } = require('./rules');
const { emptyCondition } = require('./reagents');
const { canonicalForm, findNode, neighborsOf } = require('./engine');

/* ------------------------------------------------------- molecule builders */
// Starts high, deliberately -- engine.js's own newNode/newAmineNode/
// newEtherOxygenNode draw from a SEPARATE, module-level `_id` counter
// (starting at 0) that several operators (protonateAmine, alkylateAmine,
// esterifyCommon, azoCoupling, diazotisation, formPeptideBond, ...) call
// internally to allocate a fresh node mid-transform. A hand-built test
// molecule using low ids of its own risks a genuine id COLLISION with
// whatever that internal counter allocates next -- silent corruption
// (two different nodes sharing one id), not always an obvious crash.
// generator.js sidesteps this the same way, starting its own `_gid` at
// 100000; this uses a different high offset so a test that mixes a
// generator-built molecule with a hand-built one (none currently do)
// still wouldn't collide either.
let _id = 900000;
function node(extra){ return Object.assign({ id:_id++, subs:[], oxo:false, group:null, phenyl:false, ring:false }, extra); }
function ringNode(){ return Object.assign(node({ ring:true }), { subs:new Array(6).fill(null) }); }
// A plain unbranched chain of n carbons; edges[i] joins chain[i]-chain[i+1]
// (all single by default -- mutate .type for a double bond). `sub0`, if
// given, is pushed onto the FIRST (terminal) carbon's subs.
function chain(n, sub0){
  const nodes=[]; for(let i=0;i<n;i++) nodes.push(node());
  const edges=[]; for(let i=0;i<n-1;i++) edges.push({a:nodes[i].id, b:nodes[i+1].id, type:'S'});
  if(sub0) nodes[0].subs.push(sub0);
  return { nodes, edges };
}
function merge(...mols){ return { nodes: mols.flatMap(m=>m.nodes), edges: mols.flatMap(m=>m.edges) }; }
function bond(mol, aId, bId, type){ mol.edges.push({a:aId, b:bId, type: type||'S'}); return mol; }
function cond(overrides){ return Object.assign(emptyCondition(), overrides); }

function hasGroup(mol, group){ return mol.nodes.some(n=>n.group===group); }
function hasSub(mol, sub){ return mol.nodes.some(n=>n.subs && n.subs.includes(sub)); }
function subOnNode(mol, id, sub){ const n=findNode(mol,id); return !!(n && n.subs && n.subs.includes(sub)); }
function ringSubIncludes(mol, sub){ const r=mol.nodes.find(n=>n.ring); return !!(r && r.subs.includes(sub)); }
function fragmentCount(products){ return products.filter(p=>p.kind==='chain').length; }
// A resolved result is EITHER a single `.product` OR a `.variants` list
// (when the underlying operator has a real regiochemistry/position
// choice, e.g. radical substitution, elimination, unsubstituted-ring EAS)
// -- this normalises to "the list of molecules to check", 1 or more.
function productsOf(r){ return r.product ? [r.product] : (r.variants || []); }
function countNodesWithOxo(mol){ return mol.nodes.filter(n=>n.oxo).length; }
function countDoubleBonds(mol){ return mol.edges.filter(e=>e.type==='D').length; }

/* ------------------------------------------------------------------- cases */
const cases = [];
function T(id, description, build, condOverrides, expect){
  cases.push({ id, description, build, condOverrides, expect });
}
const twoCases = [];
function T2(id, description, buildA, buildB, condOverrides, expect){
  twoCases.push({ id, description, buildA, buildB, condOverrides, expect });
}

/* ---- Alkanes ---- */
T('A1', 'propane + Cl2/UV -> chloropropane (radical substitution occurs)',
  ()=>chain(3), {reagents:['Cl2'], light:true},
  r=>r.occurs===true && r.variants && r.variants.every(v=>hasSub(v,'Cl')));
T('A2', 'propane + Cl2, NO UV -> no reaction (radical mechanism needs light)',
  ()=>chain(3), {reagents:['Cl2'], light:false},
  r=>r.occurs===false);
T('A3', 'propane + O2 (excess) -> combustion, no organic product',
  ()=>chain(3), {reagents:['O2'], quantity:'excess'},
  r=>r.occurs===true && r.products && r.products.length===1 && r.products[0].kind==='small');

/* ---- Alkenes (regiochemistry is the real test here) ---- */
T('K1', 'but-1-ene + HBr -> 2-bromobutane (Markovnikov: Br on C2, the more-substituted end)',
  ()=>{ const m=chain(4); m.edges[0].type='D'; return m; }, {reagents:['HBr']},
  r=>{
    if(!r.occurs) return false;
    const alphaC = r.product.nodes[1]; // chain[1] = but-1-ene's C2, the more-substituted alkene carbon
    return subOnNode(r.product, alphaC.id, 'Br') && !subOnNode(r.product, r.product.nodes[0].id, 'Br');
  });
T('K2', 'but-1-ene + Br2 (CCl4) -> 1,2-dibromobutane (simple addition, both carbons get Br, no regiochemistry choice)',
  ()=>{ const m=chain(4); m.edges[0].type='D'; return m; }, {reagents:['Br2'], solvent:'CCl4'},
  r=>r.occurs===true && subOnNode(r.product, r.product.nodes[0].id,'Br') && subOnNode(r.product, r.product.nodes[1].id,'Br'));
T('K3', 'but-1-ene + Br2(aq) -> bromohydrin, OH on the MORE substituted carbon (opposite of the HX case), Br on the less substituted',
  ()=>{ const m=chain(4); m.edges[0].type='D'; return m; }, {reagents:['Br2'], solvent:'aq'},
  r=>r.occurs===true && subOnNode(r.product, r.product.nodes[1].id,'OH') && subOnNode(r.product, r.product.nodes[0].id,'Br'));
T('K4', 'but-2-ene + cold dilute KMnO4 -> butane-2,3-diol (both new OH groups, no cleavage)',
  ()=>{ const m=chain(4); m.edges[1].type='D'; return m; }, {reagents:['KMnO4'], temperature:'cold', concentration:'dilute'},
  r=>r.occurs===true && countDoubleBonds(r.product)===0 && subOnNode(r.product,r.product.nodes[1].id,'OH') && subOnNode(r.product,r.product.nodes[2].id,'OH'));
T('K5', 'the SAME but-2-ene + hot concentrated acidified KMnO4 -> oxidative cleavage into two separate ethanoic acid fragments (trick: same molecule, same base oxidant, different conditions -> completely different product shape)',
  ()=>{ const m=chain(4); m.edges[1].type='D'; return m; }, {reagents:['KMnO4'], temperature:'heat', concentration:'concentrated', pH:'acidic'},
  r=>r.occurs===true && r.products && fragmentCount(r.products)===2);
T('K6', 'propene + hot steam -> propan-2-ol (Markovnikov hydration)',
  ()=>{ const m=chain(3); m.edges[0].type='D'; return m; }, {reagents:['H2O'], temperature:'heat'},
  r=>r.occurs===true && subOnNode(r.product, r.product.nodes[1].id,'OH'));
T('K7', 'ethene (symmetric, tied substitution) + HBr -> bromoethane, still occurs, no regiochemistry ambiguity to break',
  ()=>{ const m=chain(2); m.edges[0].type='D'; return m; }, {reagents:['HBr']},
  r=>r.occurs===true && hasSub(r.product,'Br'));

/* ---- Arenes ---- */
function benzene(){ return { nodes:[ringNode()], edges:[] }; }
// Matches generateAlkylbenzene's own shape exactly (a phenyl:true flag
// on the chain's first carbon, not a real ring:true node) -- that's what
// sideChainOxidationToBenzoicAcid and every other side-chain operator
// actually key off (a ring:true chain, EAS-reactive on the ring itself,
// is a structurally different, richer shape this app only builds when a
// question is specifically about ring substitution).
function toluene(){ const me=node(); me.phenyl=true; return { nodes:[me], edges:[] }; }
T('R1', 'benzene + Cl2/AlCl3 -> chlorobenzene',
  benzene, {reagents:['Cl2'], catalyst:'AlCl3'},
  r=>r.occurs===true && (r.product ? ringSubIncludes(r.product,'Cl') : r.variants.every(v=>ringSubIncludes(v,'Cl'))));
T('R2', 'benzene + Cl2, NO catalyst -> no reaction (unactivated ring needs a Lewis acid)',
  benzene, {reagents:['Cl2'], catalyst:null},
  r=>r.occurs===false);
T('R3', 'benzene + HNO3 alone (no H2SO4) -> no reaction',
  benzene, {reagents:['HNO3']},
  r=>r.occurs===false);
T('R4', 'phenol + Br2(aq), EXCESS -> 2,4,6-tribromophenol (all three activated positions, no catalyst needed)',
  ()=>{ const r=ringNode(); r.subs[0]='OH'; return { nodes:[r], edges:[] }; }, {reagents:['Br2'], solvent:'aq', quantity:'excess'},
  r=>r.occurs===true && r.product.nodes[0].subs.filter(s=>s==='Br').length===3);
T('R5', 'phenol + Br2(CCl4), controlled -> monobromophenol only (same reagent species as R4, different solvent/quantity -> very different outcome)',
  ()=>{ const r=ringNode(); r.subs[0]='OH'; return { nodes:[r], edges:[] }; }, {reagents:['Br2'], solvent:'CCl4', quantity:'controlled'},
  r=>{
    if(!r.occurs) return false;
    const prods = r.product ? [r.product] : r.variants;
    return prods.every(p=>p.nodes[0].subs.filter(s=>s==='Br').length===1);
  });
T('R6', 'toluene (methylbenzene) side chain has a benzylic H -> + hot acidified KMnO4 gives benzoic acid',
  toluene, {reagents:['KMnO4'], pH:'acidic', temperature:'heat'},
  r=>r.occurs===true && hasGroup(r.product,'COOH'));
T('R7', 'tert-butylbenzene has NO benzylic H (fully substituted) -> hot acidified KMnO4 does nothing to the side chain',
  ()=>{ const c=node(); c.phenyl=true; const m1=node(); const m2=node(); const m3=node();
        return { nodes:[c,m1,m2,m3], edges:[{a:c.id,b:m1.id,type:'S'},{a:c.id,b:m2.id,type:'S'},{a:c.id,b:m3.id,type:'S'}] }; },
  {reagents:['KMnO4'], pH:'acidic', temperature:'heat'},
  r=>r.occurs===false);

/* ---- Halogen Derivatives ---- */
function haloPropane(hal){ const m=chain(3); m.nodes[0].subs.push(hal); return m; }
T('H1', '1-bromopropane + NaOH(aq), heat -> propan-1-ol (hydrolysis)',
  ()=>haloPropane('Br'), {reagents:['NaOH'], solvent:'aq', temperature:'heat'},
  r=>r.occurs===true && subOnNode(r.product, r.product.nodes[0].id,'OH'));
T('H2', 'the SAME 1-bromopropane + NaOH, ethanolic, heat -> propene (elimination -- trick: same NaOH, different solvent, completely different reaction)',
  ()=>haloPropane('Br'), {reagents:['NaOH'], solvent:'ethanolic', temperature:'heat'},
  r=>r.occurs===true && productsOf(r).every(p=>countDoubleBonds(p)===1));
T('H3', '1-fluoropropane + NaOH(aq), heat -> no reaction (C-F too strong at this level)',
  ()=>haloPropane('F'), {reagents:['NaOH'], solvent:'aq', temperature:'heat'},
  r=>r.occurs===false);
T('H4', '1-bromopropane + ethanolic NH3, excess, heat -> propan-1-amine',
  ()=>haloPropane('Br'), {reagents:['NH3'], solvent:'ethanolic', temperature:'heat', quantity:'excess'},
  r=>r.occurs===true && hasSub(r.product,'NH2'));
T('H5', '1-bromopropane + ethanolic KCN, heat -> butanenitrile (one carbon added)',
  ()=>haloPropane('Br'), {reagents:['KCN'], solvent:'ethanolic', temperature:'heat'},
  r=>r.occurs===true && hasGroup(r.product,'CN'));

/* ---- Hydroxy Compounds ---- */
function alcohol(nArms){ // 1=primary, 2=secondary, 3=tertiary
  const central = node(); const nodes=[central]; const edges=[];
  for(let i=0;i<nArms;i++){ const c=node(); nodes.push(c); edges.push({a:central.id,b:c.id,type:'S'}); }
  central.subs.push('OH');
  return { nodes, edges };
}
T('O1', 'propan-1-ol + K2Cr2O7/H+, CONTROLLED -> propanal (stops at aldehyde)',
  ()=>alcohol(1), {reagents:['K2Cr2O7'], pH:'acidic', quantity:'controlled'},
  r=>r.occurs===true && hasGroup(r.product,'CHO'));
T('O2', 'the SAME propan-1-ol + K2Cr2O7/H+, EXCESS + heat -> propanoic acid (trick: same reagent species as O1, different quantity/heat -> further oxidation)',
  ()=>alcohol(1), {reagents:['K2Cr2O7'], pH:'acidic', quantity:'excess', temperature:'heat'},
  r=>r.occurs===true && hasGroup(r.product,'COOH'));
T('O3', 'propan-2-ol (secondary) + K2Cr2O7/H+ -> propan-2-one (ketone, cannot go further)',
  ()=>alcohol(2), {reagents:['K2Cr2O7'], pH:'acidic'},
  r=>r.occurs===true && countNodesWithOxo(r.product)===1 && !hasGroup(r.product,'COOH') && !hasGroup(r.product,'CHO'));
T('O4', '2-methylpropan-2-ol (tertiary) + K2Cr2O7/H+ -> no reaction (no H on the carbinol carbon)',
  ()=>alcohol(3), {reagents:['K2Cr2O7'], pH:'acidic'},
  r=>r.occurs===false);
T('O5', 'propan-2-ol + I2/NaOH, warm -> positive iodoform (CH3-CH(OH)- shape)',
  ()=>alcohol(2), {reagents:['I2','NaOH'], temperature:'warm'},
  r=>r.occurs===true);
T('O6', 'propan-1-ol (HO-CH2-CH2-CH3) + I2/NaOH, warm -> NO iodoform (the carbinol carbon\'s only neighbour is a CH2, not a methyl -- note alcohol(1) is ethanol, HO-CH2-CH3, which IS iodoform-positive, so this needs its own 3-carbon builder)',
  ()=>chain(3, 'OH'), {reagents:['I2','NaOH'], temperature:'warm'},
  r=>r.occurs===false);
T('O7', 'ethanol + conc H2SO4, heat -> ethene (dehydration)',
  ()=>alcohol(1), {reagents:['H2SO4'], concentration:'concentrated', temperature:'heat'},
  r=>r.occurs===true && productsOf(r).every(p=>countDoubleBonds(p)===1));
T('O8', 'phenol + NaOH(aq) -> sodium phenoxide (acidic enough to react with a base, unlike a plain alcohol)',
  ()=>{ const r=ringNode(); r.subs[0]='OH'; return { nodes:[r], edges:[] }; }, {reagents:['NaOH'], solvent:'aq'},
  r=>r.occurs===true && ringSubIncludes(r.product,'ONa'));

/* ---- Carbonyl Compounds ---- */
function aldehyde(){ const c=node(); const a=node(); a.group='CHO'; return { nodes:[c,a], edges:[{a:c.id,b:a.id,type:'S'}] }; }
function methylKetone(){ // CH3-CO-CH2CH3 (butan-2-one)
  const me=node(); const co=node(); co.oxo=true; const c1=node(); const c2=node();
  return { nodes:[me,co,c1,c2], edges:[{a:me.id,b:co.id,type:'S'},{a:co.id,b:c1.id,type:'S'},{a:c1.id,b:c2.id,type:'S'}] };
}
function symmetricKetone(){ // pentan-3-one, CH3CH2-CO-CH2CH3 -- no methyl directly on the carbonyl, so iodoform-negative
  const co=node(); co.oxo=true;
  const a1=node(), a2=node(), b1=node(), b2=node();
  return { nodes:[co,a1,a2,b1,b2], edges:[{a:co.id,b:a1.id,type:'S'},{a:a1.id,b:a2.id,type:'S'},{a:co.id,b:b1.id,type:'S'},{a:b1.id,b:b2.id,type:'S'}] };
}
T('C1', 'propanal + Tollens-equivalent (K2Cr2O7/H+) -> propanoic acid',
  aldehyde, {reagents:['K2Cr2O7'], pH:'acidic'},
  r=>r.occurs===true && hasGroup(r.product,'COOH'));
T('C2', 'butan-2-one (ketone) + K2Cr2O7/H+ -> no reaction (ketones do not oxidise this way)',
  methylKetone, {reagents:['K2Cr2O7'], pH:'acidic'},
  r=>r.occurs===false);
T('C3', 'butan-2-one + I2/NaOH -> positive iodoform (has a methyl directly on the carbonyl)',
  methylKetone, {reagents:['I2','NaOH'], temperature:'warm'},
  r=>r.occurs===true);
T('C4', 'pentan-3-one + I2/NaOH -> NO iodoform (no methyl directly on the carbonyl carbon)',
  symmetricKetone, {reagents:['I2','NaOH'], temperature:'warm'},
  r=>r.occurs===false);
T('C5', 'propanal + HCN -> cyanohydrin (nucleophilic addition, not oxidation, even though a carbonyl is present)',
  aldehyde, {reagents:['HCN']},
  r=>r.occurs===true && hasGroup(r.product,'CN') && hasSub(r.product,'OH'));

/* ---- Carboxylic Acids & Derivatives ---- */
function acid(n){ const m = (function(){ const nodes=[]; for(let i=0;i<n-1;i++) nodes.push(node()); const edges=[]; for(let i=0;i<nodes.length-1;i++) edges.push({a:nodes[i].id,b:nodes[i+1].id,type:'S'}); const g=node(); g.group='COOH'; nodes.push(g); if(nodes.length>1) edges.push({a:nodes[nodes.length-2].id,b:g.id,type:'S'}); return {nodes,edges}; })(); return m; }
function ester(){ // methyl propanoate: CH3CH2-C(=O)-O-CH3
  const a1=node(), a2=node(), carbonyl=node(); carbonyl.oxo=true;
  const o=node(); o.element='O';
  const me=node();
  return { nodes:[a1,a2,carbonyl,o,me], edges:[{a:a1.id,b:a2.id,type:'S'},{a:a2.id,b:carbonyl.id,type:'S'},{a:carbonyl.id,b:o.id,type:'S'},{a:o.id,b:me.id,type:'S'}] };
}
function primaryAmide(){ const c=node(); const a=node(); a.oxo=true; a.subs.push('NH2'); return { nodes:[c,a], edges:[{a:c.id,b:a.id,type:'S'}] }; }
T('X1', 'propanoic acid + NaOH -> sodium propanoate',
  ()=>acid(3), {reagents:['NaOH']},
  r=>r.occurs===true && hasGroup(r.product,'COONa'));
T('X2', 'propanoic acid + LiAlH4/dry ether -> propan-1-ol',
  ()=>acid(3), {reagents:['LiAlH4'], solvent:'dry ether'},
  r=>r.occurs===true && hasSub(r.product,'OH') && !hasGroup(r.product,'COOH'));
T('X3', 'propanoic acid + NaBH4 -> no reaction (NaBH4 too weak for a carboxylic acid)',
  ()=>acid(3), {reagents:['NaBH4']},
  r=>r.occurs===false);
T('X4', 'methyl propanoate + dilute acid, heat -> propanoic acid + methanol (acid hydrolysis, equilibrium)',
  ester, {reagents:['H2O'], pH:'acidic', temperature:'heat'},
  r=>r.occurs===true && fragmentCount(r.products)===2 && r.products.some(p=>hasGroup(p.mol,'COOH')));
T('X5', 'the SAME methyl propanoate + NaOH(aq), heat -> sodium propanoate + methanol (trick: same ester, alkaline instead of acidic -> salt not free acid)',
  ester, {reagents:['NaOH'], solvent:'aq', temperature:'heat'},
  r=>r.occurs===true && fragmentCount(r.products)===2 && r.products.some(p=>hasGroup(p.mol,'COONa')));
T('X6', 'propanamide (primary) + NaOH(aq), heat -> sodium propanoate',
  primaryAmide, {reagents:['NaOH'], solvent:'aq', temperature:'heat'},
  r=>r.occurs===true && hasGroup(r.product,'COONa'));
T('X7', 'propanamide + LiAlH4/dry ether -> propan-1-amine (C=O reduced all the way to CH2)',
  primaryAmide, {reagents:['LiAlH4'], solvent:'dry ether'},
  r=>r.occurs===true && countNodesWithOxo(r.product)===0 && hasSub(r.product,'NH2'));

/* ---- Nitrogen Compounds ---- */
function phenylamine(){ const r=ringNode(); r.subs[0]='NH2'; return { nodes:[r], edges:[] }; }
function nitrobenzene(){ const r=ringNode(); r.subs[0]='NO2'; return { nodes:[r], edges:[] }; }
T('N1', 'phenylamine + NaNO2/HCl, COLD -> diazonium salt (not protonation, even though HCl alone is present)',
  phenylamine, {reagents:['NaNO2','HCl'], temperature:'cold'},
  r=>r.occurs===true && hasGroup(r.product,'N2Cl'));
T('N2', 'phenylamine + HCl ALONE (no NaNO2) -> phenylammonium chloride (simple protonation)',
  phenylamine, {reagents:['HCl']},
  r=>r.occurs===true && !hasGroup(r.product,'N2Cl') && r.product.nodes.some(n=>n.element==='N' && n.charge===1));
T('N3', 'phenylamine + NaNO2/HCl at ROOM temperature -> no reaction (diazonium salt needs cold; this app does not model a lesser side-reaction, so this should cleanly fail rather than silently mis-fire)',
  phenylamine, {reagents:['NaNO2','HCl'], temperature:'room'},
  r=>r.occurs===false);
T('N4', 'nitrobenzene + Sn/conc HCl, heat -> phenylamine (full reduction)',
  nitrobenzene, {reagents:['Sn','HCl'], temperature:'heat'},
  r=>r.occurs===true && ringSubIncludes(r.product,'NH2'));

/* ---- Multi-functional-group molecules (the real point of this harness) ---- */
function aminoAcidLike(){ // H2N-CH2-COOH (glycine)
  const alpha=node(); alpha.subs.push('NH2'); const acidC=node(); acidC.group='COOH';
  return { nodes:[alpha,acidC], edges:[{a:alpha.id,b:acidC.id,type:'S'}] };
}
T('M1', 'glycine (has BOTH -NH2 and -COOH) + HCl -> only the amine protonates, the -COOH is untouched',
  aminoAcidLike, {reagents:['HCl']},
  r=>r.occurs===true && hasGroup(r.product,'COOH') && r.product.nodes.some(n=>n.element==='N'&&n.charge===1));
T('M2', 'glycine + NaOH -> only the -COOH deprotonates to -COONa, the -NH2 is untouched',
  aminoAcidLike, {reagents:['NaOH']},
  r=>r.occurs===true && hasGroup(r.product,'COONa') && hasSub(r.product,'NH2'));
T('M3', 'an unsaturated primary alcohol (HO-CH2-CH=CH-CH3) + H2/Ni -> ONLY the C=C reduces, the -OH survives untouched',
  ()=>{ const m=chain(4); m.nodes[0].subs.push('OH'); m.edges[1].type='D'; return m; }, {reagents:['H2'], catalyst:'Ni'},
  r=>r.occurs===true && countDoubleBonds(r.product)===0 && hasSub(r.product,'OH'));
T('M4', 'the SAME unsaturated alcohol + hot concentrated acidified KMnO4 -> BOTH groups react: the C=C cleaves AND the primary alcohol end oxidises (a real two-site question)',
  ()=>{ const m=chain(4); m.nodes[0].subs.push('OH'); m.edges[1].type='D'; return m; }, {reagents:['KMnO4'], concentration:'concentrated', pH:'acidic', temperature:'heat'},
  r=>{
    if(!r.occurs) return false;
    // cleavage always returns a fragment set; among the fragments, the one
    // that WAS the alcohol-bearing half should now carry a carboxyl group
    // (oxidised past the cleavage, not merely cleaved)
    return r.products && r.products.some(p=>p.kind==='chain' && (hasGroup(p.mol,'COOH')||hasGroup(p.mol,'COONa')));
  });
T('M5', 'a bromoalkane with a spare C=C elsewhere (Br-CH2-CH2-CH=CH2) + NaOH(aq), heat -> only the C-Br hydrolyses to -OH, the alkene is untouched by aqueous NaOH',
  ()=>{ const m=chain(4); m.nodes[0].subs.push('Br'); m.edges[2].type='D'; return m; }, {reagents:['NaOH'], solvent:'aq', temperature:'heat'},
  r=>r.occurs===true && subOnNode(r.product, r.product.nodes[0].id,'OH') && countDoubleBonds(r.product)===1);
T('M6', 'a ketone that ALSO happens to be a methyl ketone AND has its own separate alkene (CH2=CH-CH2-CO-CH3) + Br2(CCl4) -> only the alkene adds Br2; the ketone is untouched by Br2/CCl4 at this level',
  ()=>{ const c1=node(),c2=node(),c3=node(),co=node(),me=node(); co.oxo=true;
        return { nodes:[c1,c2,c3,co,me], edges:[{a:c1.id,b:c2.id,type:'D'},{a:c2.id,b:c3.id,type:'S'},{a:c3.id,b:co.id,type:'S'},{a:co.id,b:me.id,type:'S'}] }; },
  {reagents:['Br2'], solvent:'CCl4'},
  r=>r.occurs===true && countDoubleBonds(r.product)===0 && countNodesWithOxo(r.product)===1);
T('M7', 'that SAME molecule + I2/NaOH -> the methyl ketone gives a positive iodoform, the alkene is untouched by this reagent',
  ()=>{ const c1=node(),c2=node(),c3=node(),co=node(),me=node(); co.oxo=true;
        return { nodes:[c1,c2,c3,co,me], edges:[{a:c1.id,b:c2.id,type:'D'},{a:c2.id,b:c3.id,type:'S'},{a:c3.id,b:co.id,type:'S'},{a:co.id,b:me.id,type:'S'}] }; },
  {reagents:['I2','NaOH'], temperature:'warm'},
  r=>r.occurs===true);
T('M8', 'a molecule with both a phenolic -OH and a separate benzylic CH3 on the SAME ring (4-methylphenol) + excess Br2(aq) -> ring tribrominates (activated by -OH); the methyl side chain is untouched by aqueous Br2',
  ()=>{ const r=ringNode(); r.subs[0]='OH'; const me=node(); return { nodes:[r,me], edges:[{a:r.id,b:me.id,type:'S',ringPos:3}] }; },
  {reagents:['Br2'], solvent:'aq', quantity:'excess'},
  r=>r.occurs===true && r.product.nodes.find(n=>n.ring).subs.filter(s=>s==='Br').length>=2);

/* ---- More adversarial cases ---- */
function chlorobenzene(){ const r=ringNode(); r.subs[0]='Cl'; return { nodes:[r], edges:[] }; }
T('P1', 'chlorobenzene + NaOH(aq), heat -> no reaction (an aryl C-X bond does not undergo simple nucleophilic substitution at this level, unlike an aliphatic haloalkane)',
  chlorobenzene, {reagents:['NaOH'], solvent:'aq', temperature:'heat'},
  r=>r.occurs===false);
function benzaldehyde(){ const r=ringNode(); const c=node(); c.group='CHO'; return { nodes:[r,c], edges:[{a:r.id,b:c.id,type:'S',ringPos:0}] }; }
T('P2', 'benzaldehyde + K2Cr2O7/H+ -> benzoic acid (an aromatic aldehyde still oxidises with dichromate, even though it fails Fehling\'s -- this app only models the dichromate pathway)',
  benzaldehyde, {reagents:['K2Cr2O7'], pH:'acidic'},
  r=>r.occurs===true && hasGroup(r.product,'COOH'));
function methanoicAcid(){ const c=node(); c.group='COOH'; return { nodes:[c], edges:[] }; }
function ethanedioicAcid(){ const c1=node(); c1.group='COOH'; const c2=node(); c2.group='COOH'; return { nodes:[c1,c2], edges:[{a:c1.id,b:c2.id,type:'S'}] }; }
T('P3', 'methanoic acid + hot acidified KMnO4 -> fully oxidised (CO2+H2O), unlike every other carboxylic acid',
  methanoicAcid, {reagents:['KMnO4'], pH:'acidic', temperature:'heat'},
  r=>r.occurs===true && r.products && r.products.length===1 && r.products[0].kind==='small');
T('P4', 'ethanedioic acid + hot acidified KMnO4 -> also fully oxidised (weak C-C bond between the two carboxyls)',
  ethanedioicAcid, {reagents:['KMnO4'], pH:'acidic', temperature:'heat'},
  r=>r.occurs===true && r.products && r.products.length===1 && r.products[0].kind==='small');
T('P5', 'propanoic acid (a perfectly ordinary carboxylic acid) + hot acidified KMnO4 -> no reaction (already at its most oxidised state -- specialOxidisable must NOT over-fire on a plain acid)',
  ()=>acid(3), {reagents:['KMnO4'], pH:'acidic', temperature:'heat'},
  r=>r.occurs===false);
T('P6', 'an alpha,beta-unsaturated aldehyde (CH2=CH-CHO) + HCN -> adds to the carbonyl; the conjugated C=C does not confuse the site selection',
  ()=>{ const c1=node(),c2=node(),cho=node(); cho.group='CHO'; return { nodes:[c1,c2,cho], edges:[{a:c1.id,b:c2.id,type:'D'},{a:c2.id,b:cho.id,type:'S'}] }; },
  {reagents:['HCN']},
  r=>r.occurs===true && hasGroup(r.product,'CN') && countDoubleBonds(r.product)===1);
T('P7', 'a full dipeptide (built via formPeptideBond itself, not hand-drawn) + dilute acid, heat -> hydrolyses back to its two amino acids, run entirely through the RULE ENGINE (fact scan correctly recognises the peptide bond as a secondary amide)',
  ()=>{
    const op2 = require('./operators');
    const gly = (()=>{ const a=node(); a.subs.push('NH2'); const g=node(); g.group='COOH'; return { nodes:[a,g], edges:[{a:a.id,b:g.id,type:'S'}] }; })();
    const ala = (()=>{ const a=node(); a.subs.push('NH2'); const me=node(); const g=node(); g.group='COOH';
      return { nodes:[a,me,g], edges:[{a:a.id,b:me.id,type:'S'},{a:a.id,b:g.id,type:'S'}] }; })();
    return op2.formPeptideBond(gly, ala).product;
  },
  {reagents:['H2O'], pH:'acidic', temperature:'heat'},
  r=>r.occurs===true && r.products && fragmentCount(r.products)===2 && r.products.every(p=>{
    const n = p.mol.nodes.find(nn=>nn.element==='N'); return n && n.charge===1; // acidic hydrolysis -> both freed amines protonated to -NH3+
  }));

/* ---- Two-reactant reactions ---- */
function acidChain(n){ const nodes=[]; for(let i=0;i<n-1;i++) nodes.push(node()); const edges=[]; for(let i=0;i<nodes.length-1;i++) edges.push({a:nodes[i].id,b:nodes[i+1].id,type:'S'}); const g=node(); g.group='COOH'; nodes.push(g); if(nodes.length>1) edges.push({a:nodes[nodes.length-2].id,b:g.id,type:'S'}); return {nodes,edges}; }
function alcoholSimple(nArms){ const central=node(); const nodes=[central]; const edges=[]; for(let i=0;i<nArms;i++){ const c=node(); nodes.push(c); edges.push({a:central.id,b:c.id,type:'S'}); } central.subs.push('OH'); return {nodes,edges}; }
function phenolSimple(){ const r=ringNode(); r.subs[0]='OH'; return { nodes:[r], edges:[] }; }
function acylChlorideSimple(n){ const nodes=[]; for(let i=0;i<n-1;i++) nodes.push(node()); const edges=[]; for(let i=0;i<nodes.length-1;i++) edges.push({a:nodes[i].id,b:nodes[i+1].id,type:'S'}); const c=node(); c.oxo=true; c.subs.push('Cl'); nodes.push(c); if(nodes.length>1) edges.push({a:nodes[nodes.length-2].id,b:c.id,type:'S'}); return {nodes,edges}; }
function benzeneSimple(){ return { nodes:[ringNode()], edges:[] }; }
function chloroalkane(n){ const m=chain(n); m.nodes[0].subs.push('Cl'); return m; }
function primaryAmineSimple(nArms){ const central=node(); const nodes=[central]; const edges=[]; for(let i=0;i<nArms;i++){ const c=node(); nodes.push(c); edges.push({a:central.id,b:c.id,type:'S'}); } central.subs.push('NH2'); return {nodes,edges}; }
function benzenediazonium(){ const r=ringNode(); const g=node(); g.group='N2Cl'; return { nodes:[r,g], edges:[{a:r.id,b:g.id,type:'S',ringPos:0}] }; }
function glycineLike(){ const a=node(); a.subs.push('NH2'); const g=node(); g.group='COOH'; return { nodes:[a,g], edges:[{a:a.id,b:g.id,type:'S'}] }; }
function alanineLike(){ const a=node(); a.subs.push('NH2'); const me=node(); const g=node(); g.group='COOH'; return { nodes:[a,me,g], edges:[{a:a.id,b:me.id,type:'S'},{a:a.id,b:g.id,type:'S'}] }; }

T2('E1', 'propanoic acid + ethanol + conc H2SO4, heat -> ethyl propanoate (Fischer esterification)',
  ()=>acidChain(3), ()=>alcoholSimple(1), {catalyst:'H2SO4', concentration:'concentrated', temperature:'heat'},
  r=>r.occurs===true && hasGroup(r.product,'CN')===false && r.product.nodes.some(n=>n.element==='O'));
T2('E2', 'ethanoic acid + phenol + conc H2SO4, heat -> NO reaction (Fischer esterification does not work directly on a phenol)',
  ()=>acidChain(2), phenolSimple, {catalyst:'H2SO4', concentration:'concentrated', temperature:'heat'},
  r=>r.occurs===false);
T2('E3', 'the SAME phenol + propanoyl CHLORIDE instead (not the acid), room temperature -> DOES esterify (trick: same alcohol partner as E2, different acyl partner reactivity)',
  ()=>acylChlorideSimple(3), phenolSimple, {},
  r=>r.occurs===true && r.product.nodes.some(n=>n.element==='O'));
T2('E4', 'propanoic acid + propan-2-ol (secondary) + conc H2SO4, heat -> still esterifies (alcohol class does not matter for esterification)',
  ()=>acidChain(3), ()=>alcoholSimple(2), {catalyst:'H2SO4', concentration:'concentrated', temperature:'heat'},
  r=>r.occurs===true);
T2('F1', 'benzene + 1-chloropropane + AlCl3 -> propylbenzene (Friedel-Crafts alkylation)',
  benzeneSimple, ()=>chloroalkane(3), {catalyst:'AlCl3'},
  r=>r.occurs===true && productsOf(r).every(p=>p.nodes.some(n=>n.phenyl)));
T2('F2', 'the SAME pair, NO AlCl3 -> no reaction (needs the Lewis acid to generate the carbocation)',
  benzeneSimple, ()=>chloroalkane(3), {},
  r=>r.occurs===false);
T2('G1', 'propanoyl chloride + propan-1-amine, room temperature -> N-propylpropanamide',
  ()=>acylChlorideSimple(3), ()=>primaryAmineSimple(3), {},
  r=>r.occurs===true && r.product.nodes.some(n=>n.element==='N') && countNodesWithOxo(r.product)===1);
T2('Z1', 'benzenediazonium chloride + phenol, cold, alkaline -> an azo dye (ring coupling)',
  benzenediazonium, phenolSimple, {temperature:'cold'},
  r=>r.occurs===true && productsOf(r).every(p=>p.nodes.some(n=>n.element==='N')));
T2('PEP1', 'glycine + alanine, condensation conditions -> a dipeptide forms; BOTH orderings are chemically valid (Gly-Ala and Ala-Gly are different molecules), so this should come back as 2 distinct variants, not a clash',
  glycineLike, alanineLike, {},
  r=>{
    if(!r.occurs) return false;
    const prods = productsOf(r);
    if(prods.length!==2) return false;
    const forms = prods.map(canonicalForm);
    return forms[0] !== forms[1]; // genuinely two different dipeptides, not one product duplicated
  });

function run(){
  let pass=0, fail=0, errors=0;
  for(const c of cases){
    _id = 900000; // deterministic per-case reset, but still clear of engine.js's own internal _id counter -- see the top-of-file comment
    let mol, condition, result, ok, detail='';
    try{
      mol = c.build();
      condition = cond(c.condOverrides);
      result = resolve(mol, condition, RULES);
      ok = !!c.expect(result);
    } catch(e){
      if(e instanceof RuleClashError){ errors++; console.log(`CLASH  ${c.id} - ${c.description}\n       ${e.message}`); continue; }
      errors++; console.log(`ERROR  ${c.id} - ${c.description}\n       ${e.stack.split('\n').slice(0,3).join('\n       ')}`); continue;
    }
    if(ok){ pass++; console.log(`PASS   ${c.id} - ${c.description}`); }
    else {
      fail++;
      detail = result.occurs
        ? `occurs=true, appliedRules=[${(result.appliedRules||[]).join(',')}]`
        : `occurs=false, candidateRuleIds=[${(result.candidateRuleIds||[]).join(',')}]`;
      console.log(`FAIL   ${c.id} - ${c.description}\n       got: ${detail}`);
    }
  }
  for(const c of twoCases){
    _id = 900000;
    let molA, molB, condition, result, ok;
    try{
      molA = c.buildA();
      molB = c.buildB();
      condition = cond(c.condOverrides);
      result = resolveTwoReactant(molA, molB, condition, RULES);
      ok = !!c.expect(result);
    } catch(e){
      if(e instanceof TwoReactantRuleClashError){ errors++; console.log(`CLASH  ${c.id} - ${c.description}\n       ${e.message}`); continue; }
      errors++; console.log(`ERROR  ${c.id} - ${c.description}\n       ${e.stack.split('\n').slice(0,3).join('\n       ')}`); continue;
    }
    if(ok){ pass++; console.log(`PASS   ${c.id} - ${c.description}`); }
    else {
      fail++;
      console.log(`FAIL   ${c.id} - ${c.description}\n       got: occurs=${result.occurs}, appliedRules=[${(result.appliedRules||[]).join(',')}]`);
    }
  }
  const total = cases.length + twoCases.length;
  console.log(`\n${pass} passed, ${fail} failed, ${errors} errored (clash/exception) -- ${total} total`);
  process.exit(fail>0||errors>0 ? 1 : 0);
}

run();
