'use strict';
const { bondSlotsUsed, subSlotsUsed } = require('./engine');

function ri(min,max){ return min+Math.floor(Math.random()*(max-min+1)); }
function pick(arr){ return arr[Math.floor(Math.random()*arr.length)]; }
function freeSlotsRaw(mol,id){
  const n = mol.nodes.find(x=>x.id===id);
  return 4 - bondSlotsUsed(mol,id) - subSlotsUsed(n);
}

let _gid = 100000;
function mkNode(){ return { id:_gid++, subs:[], oxo:false, group:null, phenyl:false }; }

const DIFF = {
  chainRange:      [[2,3],[2,4],[3,5],[4,7],[5,9]],
  branchProb:      [0,   0.15, 0.3,  0.45, 0.6],
  maxBranches:     [0,   1,    1,    2,    3],
  doubleBondMax:   [1,   1,    1,    2,    2],
  phenylProb:      [0,   0.1,  0.2,  0.3,  0.35],
  existingHalProb: [0,   0.1,  0.2,  0.3,  0.35],
};

/* constraints:
   minDoubleBonds, forceNoDoubleBonds, forcePhenyl, forceExistingHalogen,
   minChainLen
*/
function generateMolecule(difficulty, constraints={}){
  const d = Math.min(5, Math.max(1, difficulty));
  const [lo,hi] = DIFF.chainRange[d-1];
  let chainLen = ri(lo,hi);
  if(constraints.minChainLen) chainLen = Math.max(chainLen, constraints.minChainLen);

  const nodes=[]; const edges=[];
  const chainIds=[];
  for(let i=0;i<chainLen;i++){ const n=mkNode(); nodes.push(n); chainIds.push(n.id); }
  for(let i=0;i<chainLen-1;i++) edges.push({a:chainIds[i], b:chainIds[i+1], type:'S'});
  const mol = { nodes, edges };

  // ---- double bonds (spaced, non-adjacent edges) ----
  let dbCount;
  if(constraints.forceNoDoubleBonds) dbCount = 0;
  else if(constraints.minDoubleBonds!=null){
    const maxPossible = Math.max(constraints.minDoubleBonds, DIFF.doubleBondMax[d-1]);
    dbCount = ri(constraints.minDoubleBonds, maxPossible);
  } else {
    dbCount = Math.random()<0.55 ? 0 : ri(1, DIFF.doubleBondMax[d-1]);
  }
  const availableEdgeIdx = [...Array(chainLen-1).keys()];
  const chosen=[];
  const shuffled = availableEdgeIdx.slice().sort(()=>Math.random()-0.5);
  for(const idx of shuffled){
    if(chosen.length>=dbCount) break;
    if(chosen.some(c=>Math.abs(c-idx)<2)) continue; // no adjacent double bonds (avoid conjugation)
    chosen.push(idx);
  }
  chosen.forEach(idx=>{ edges[idx].type='D'; });
  const sp2Ids = new Set();
  chosen.forEach(idx=>{ sp2Ids.add(edges[idx].a); sp2Ids.add(edges[idx].b); });

  // ---- branches (sp3 chain carbons only, keep sp2 carbons unsubstituted for clean cleavage logic) ----
  const maxBranches = DIFF.maxBranches[d-1];
  let branchesAdded = 0;
  const branchCandidates = chainIds.filter(id=>!sp2Ids.has(id)).sort(()=>Math.random()-0.5);
  for(const cid of branchCandidates){
    if(branchesAdded>=maxBranches) break;
    if(Math.random()>=DIFF.branchProb[d-1]) continue;
    if(freeSlotsRaw(mol,cid) < 1) continue;
    const branchLen = ri(1,2);
    let prev = cid;
    for(let i=0;i<branchLen;i++){
      if(freeSlotsRaw(mol,prev) < 1) break;
      const bn = mkNode(); nodes.push(bn);
      edges.push({a:prev, b:bn.id, type:'S'});
      prev = bn.id;
    }
    branchesAdded++;
  }

  // ---- phenyl (one, attached to an sp3 carbon only -- keep sp2 carbons clean) ----
  const wantsPhenyl = constraints.forcePhenyl===true || (constraints.forcePhenyl!==false && Math.random()<DIFF.phenylProb[d-1]);
  if(wantsPhenyl){
    const cands = mol.nodes.filter(n=>!sp2Ids.has(n.id) && freeSlotsRaw(mol,n.id)>=1);
    if(cands.length){ pick(cands).phenyl = true; }
  }

  // ---- pre-existing halogen (for "further substitution" scenarios; sp3 only) ----
  const wantsHalogen = constraints.forceExistingHalogen===true || (constraints.forceExistingHalogen!==false && Math.random()<DIFF.existingHalProb[d-1]);
  if(wantsHalogen){
    const cands = mol.nodes.filter(n=>!sp2Ids.has(n.id) && !n.phenyl && freeSlotsRaw(mol,n.id)>=1);
    if(cands.length){ pick(cands).subs.push(constraints.forceHalogenType || pick(['Cl','Br','F'])); }
  }

  return mol;
}

// Arene starting materials -- ring:true nodes here (not the opaque
// phenyl:true flag Alkane/Alkene questions use), since the ring itself is
// what the arene-topic reactions act on and its existing substituent (if
// any) needs to be inspectable for ortho/para directing. Uses its own
// mkNode()-shaped node builder throughout (never calls engine.js's own
// node/ring constructors) to stay in this module's own _gid id space,
// matching how generateMolecule already works.
function mkRingNode(){ return { id:_gid++, subs:new Array(6).fill(null), oxo:false, group:null, phenyl:false, ring:true }; }

// chainLen: 0/falsy for plain benzene, or N for a plain unbranched N-carbon
// alkyl chain at position 0 (1=methylbenzene, 2=ethylbenzene, 3=
// propylbenzene, ...) -- ringSubstitutionPositions (operators.js) only
// ever looks at WHICH position is occupied, never what's actually there,
// so any chain length is exactly as valid a "one existing substituent"
// starting shape as methylbenzene was. A ring with 2+ existing
// substituents would need real directing-group conflict resolution, out
// of scope here (ringSubstitutionPositions deliberately refuses that case
// rather than guessing).
function generateArene(chainLen){
  const ring = mkRingNode();
  const nodes = [ring];
  const edges = [];
  if(chainLen){
    const chainIds = [];
    for(let i=0;i<chainLen;i++){ const n=mkNode(); nodes.push(n); chainIds.push(n.id); }
    for(let i=0;i<chainLen-1;i++) edges.push({a:chainIds[i], b:chainIds[i+1], type:'S'});
    edges.push({a:ring.id, b:chainIds[0], type:'S', ringPos:0});
  }
  return { nodes, edges };
}

// Ar-(CH2)(n-1)-CH3 -- a plain, unbranched side chain of chainLen carbons,
// the first one bonded directly to the ring via the OLD opaque phenyl:true
// flag (this reaction family never touches the ring itself, so there's no
// need for the richer ring:true representation -- and using phenyl:true
// keeps the product directly comparable to what a graded submission folds
// down to, with no extra downgrade step required). chainLen=1 is toluene
// methylbenzene); longer chains are generic "alkylbenzene".
function generateAlkylbenzene(chainLen){
  const nodes = []; const edges = [];
  const chainIds = [];
  for(let i=0;i<chainLen;i++){ const n=mkNode(); nodes.push(n); chainIds.push(n.id); }
  for(let i=0;i<chainLen-1;i++) edges.push({a:chainIds[i], b:chainIds[i+1], type:'S'});
  nodes[0].phenyl = true;
  return { nodes, edges };
}

// Ph-C(CH3)3 -- the benzylic carbon is quaternary (bonded to the ring plus
// 3 methyls, filling all 4 valence slots), so it has NO benzylic H and
// can't be oxidised -- the one molecule the "no reaction" side-chain
// oxidation question needs.
function generateTertButylbenzene(){
  const nodes = []; const edges = [];
  const central = mkNode(); central.phenyl = true; nodes.push(central);
  for(let i=0;i<3;i++){ const m=mkNode(); nodes.push(m); edges.push({a:central.id, b:m.id, type:'S'}); }
  return { nodes, edges };
}

// A plain, unbranched 1-chloroalkane (Cl on the terminal carbon) -- the R-Cl
// reactant for Friedel-Crafts alkylation. Shown to the student as its own
// molecule (not just described as "a chloroalkane" in the reagent text),
// since which one is used determines the ring's new substituent -- the
// SAME chainLen this returns has to be threaded through to whichever
// operator consumes it, so the question and the correct answer never
// disagree about which alkyl chloride was actually used.
function generateAlkylChloride(chainLen){
  const nodes = []; const edges = [];
  const chainIds = [];
  for(let i=0;i<chainLen;i++){ const n=mkNode(); nodes.push(n); chainIds.push(n.id); }
  for(let i=0;i<chainLen-1;i++) edges.push({a:chainIds[i], b:chainIds[i+1], type:'S'});
  nodes[chainLen-1].subs.push('Cl');
  return { nodes, edges };
}

// A halogenoarene -- ring with a single halogen X directly on one position,
// nothing else. Used as the starting material for both "haloarenes don't
// undergo nucleophilic substitution/elimination" (always no-reaction) and
// "haloarenes DO still undergo electrophilic aromatic substitution, just
// less readily than benzene" -- the latter reuses ringElectrophilicSubstitution
// unchanged, since a halogen is (unusually) a deactivating but still
// ortho/para-directing group, so the existing ortho+para position logic
// already gives the chemically correct answer.
function generateHaloarene(X){
  const ring = mkRingNode();
  ring.subs[0] = X;
  return { nodes:[ring], edges:[] };
}

// A ring bearing 1 or 2 existing substituents drawn from a broader palette
// than plain alkyl -- flat op-directors (OH, NH2, Cl, Br) alongside alkyl
// chains, or meta-directors (NO2 flat, or a COOH/CHO/CN group node on a
// ringPos edge) -- for richer electrophilic-substitution starting
// materials. When two substituents are placed they're always drawn from the
// SAME directing family (both op-director or both meta-director): operators.js's
// ringSubstitutionPositions only models "favoured by at least one existing
// substituent", which stays chemically well-posed for two groups that agree
// on direction (e.g. p-/m-xylene-style pairs) but has no notion of one
// group's directing effect dominating a genuinely opposed one (e.g. methyl
// vs nitro), so mixed-family pairs are deliberately never generated.
// difficulty gates how adventurous the substituent choice gets: low
// difficulty stays a single short alkyl chain (methylbenzene-style,
// unchanged from before), difficulty 3+ can swap in a non-alkyl op-director,
// difficulty 4+ can add a second substituent, and only difficulty 4+ ever
// risks the (rarer, less commonly drilled) meta-director family at all --
// UNLESS opts.allowMeta is explicitly false, for reactions (Friedel-Crafts
// alkylation) that genuinely fail outright on a strongly deactivated ring in
// real chemistry, not just direct to a different position.
function generateSubstitutedArene(d, opts){
  opts = opts || {};
  const ring = mkRingNode();
  const nodes = [ring]; const edges = [];
  function chainLen(){ const maxLen = d<=2?1:(d<=4?2:3); return 1+ri(0,maxLen-1); }
  function placeAlkyl(pos){
    const len = chainLen();
    const chainIds = [];
    for(let i=0;i<len;i++){ const c=mkNode(); nodes.push(c); chainIds.push(c.id); }
    for(let i=0;i<len-1;i++) edges.push({a:chainIds[i], b:chainIds[i+1], type:'S'});
    edges.push({a:ring.id, b:chainIds[0], type:'S', ringPos:pos});
  }
  function placeFlat(pos, sym){ ring.subs[pos] = sym; }
  function placeGroup(pos, sym){
    const g = mkNode(); g.group = sym; nodes.push(g);
    edges.push({a:ring.id, b:g.id, type:'S', ringPos:pos});
  }
  const family = (opts.allowMeta!==false && d>=4 && Math.random()<0.15) ? 'm' : 'op';
  function placeOne(pos){
    if(family==='m'){
      const sym = pick(['NO2','NO2','COOH','CHO','CN']); // NO2 weighted heavier -- the far more commonly drilled meta-director at this level
      if(sym==='NO2') placeFlat(pos,'NO2'); else placeGroup(pos, sym);
    } else if(d>=3 && Math.random()<0.35){
      placeFlat(pos, pick(['OH','NH2','Cl','Br']));
    } else {
      placeAlkyl(pos);
    }
  }
  const positions = [0,1,2,3,4,5].sort(()=>Math.random()-0.5);
  const count = (d>=4 && Math.random()<0.45) ? 2 : 1;
  for(let i=0;i<count;i++) placeOne(positions[i]);
  return { nodes, edges };
}

// A branching amine-nitrogen starting material with `chainLens.length`
// existing alkyl arms already attached (2 = secondary, 3 = tertiary) -- the
// starting point for "add one more alkyl group" questions that step a
// primary amine up towards a quaternary ammonium salt. A 1-arm (primary)
// amine is deliberately NOT built this way; that stays the flat 'NH2'
// substituent every other amine-formation question already produces (see
// operators.js's alkylateAmine, which transparently promotes that flat form
// to a real node when it needs to extend it further).
function mkAmineNode(){ return { id:_gid++, subs:[], oxo:false, group:null, phenyl:false, ring:false, element:'N', charge:0 }; }
function generateAmine(chainLens){
  const n = mkAmineNode();
  const nodes = [n]; const edges = [];
  chainLens.forEach(function(len){
    const chainIds = [];
    for(let i=0;i<len;i++){ const c=mkNode(); nodes.push(c); chainIds.push(c.id); }
    for(let i=0;i<len-1;i++) edges.push({a:chainIds[i], b:chainIds[i+1], type:'S'});
    edges.push({a:n.id, b:chainIds[0], type:'S'});
  });
  return { nodes, edges };
}

// An alcohol with armChainLens.length alkyl arms on the carbinol carbon (the
// one bearing -OH) -- 1 arm = primary, 2 = secondary, 3 = tertiary, exactly
// mirroring generateAmine's arm-based construction. The carbinol carbon
// itself is a plain carbon (unlike generateAmine's central node), since -OH
// is just an ordinary flat substituent here, not a whole new node kind.
function generateAlcohol(armChainLens){
  const central = mkNode();
  const nodes = [central]; const edges = [];
  armChainLens.forEach(function(len){
    const chainIds = [];
    for(let i=0;i<len;i++){ const c=mkNode(); nodes.push(c); chainIds.push(c.id); }
    for(let i=0;i<len-1;i++) edges.push({a:chainIds[i], b:chainIds[i+1], type:'S'});
    edges.push({a:central.id, b:chainIds[0], type:'S'});
  });
  central.subs.push('OH');
  return { nodes, edges };
}
// A phenol -- ring with a single -OH flat substituent, nothing else.
// Mirrors generateHaloarene exactly (same shape, different substituent).
function generatePhenol(){
  const ring = mkRingNode();
  ring.subs[0] = 'OH';
  return { nodes:[ring], edges:[] };
}
// A plain unbranched carboxylic acid, chainLen carbons total INCLUDING the
// carboxyl carbon itself (chainLen=1 is HCOOH, methanoic acid) -- the
// carboxyl carbon IS its own group:'COOH' node, per engine.js's model, not
// a separate node with a flat sub.
function generateCarboxylicAcid(chainLen){
  const nodes = []; const edges = [];
  const chainIds = [];
  for(let i=0;i<chainLen-1;i++){ const c=mkNode(); nodes.push(c); chainIds.push(c.id); }
  for(let i=0;i<chainIds.length-1;i++) edges.push({a:chainIds[i], b:chainIds[i+1], type:'S'});
  const acidC = mkNode(); acidC.group = 'COOH'; nodes.push(acidC);
  if(chainIds.length) edges.push({a:chainIds[chainIds.length-1], b:acidC.id, type:'S'});
  return { nodes, edges };
}
// A plain unbranched acyl chloride (R-COCl), chainLen carbons total
// including the acyl carbon -- unlike the carboxylic acid above, this
// carbon isn't a `group` node at all, just a plain carbon with oxo:true (a
// real C=O) and a flat 'Cl' sub, exactly like any other terminal
// oxo-bearing carbon already representable in this model.
function generateAcylChloride(chainLen){
  const nodes = []; const edges = [];
  const chainIds = [];
  for(let i=0;i<chainLen-1;i++){ const c=mkNode(); nodes.push(c); chainIds.push(c.id); }
  for(let i=0;i<chainIds.length-1;i++) edges.push({a:chainIds[i], b:chainIds[i+1], type:'S'});
  const acylC = mkNode(); acylC.oxo = true; acylC.subs.push('Cl'); nodes.push(acylC);
  if(chainIds.length) edges.push({a:chainIds[chainIds.length-1], b:acylC.id, type:'S'});
  return { nodes, edges };
}

// A plain unbranched aliphatic aldehyde, chainLen carbons total including
// the CHO carbon itself -- exact mirror of generateCarboxylicAcid, just
// group:'CHO' instead of 'COOH'.
function generateAldehyde(chainLen){
  const nodes = []; const edges = [];
  const chainIds = [];
  for(let i=0;i<chainLen-1;i++){ const c=mkNode(); nodes.push(c); chainIds.push(c.id); }
  for(let i=0;i<chainIds.length-1;i++) edges.push({a:chainIds[i], b:chainIds[i+1], type:'S'});
  const choC = mkNode(); choC.group = 'CHO'; nodes.push(choC);
  if(chainIds.length) edges.push({a:chainIds[chainIds.length-1], b:choC.id, type:'S'});
  return { nodes, edges };
}
// A ketone with two alkyl arms on the carbonyl carbon (oxo:true, exactly 2
// real neighbours to fill valence 4 with no free H) -- arm1/arm2 each >=1
// carbon (a 0-length arm would leave the carbonyl carbon with only 1 real
// neighbour, an aldehyde shape, not a ketone). Mirrors generateAlcohol's
// arm-based construction, just oxo instead of a flat -OH sub.
function generateKetone(arm1, arm2){
  const central = mkNode(); central.oxo = true;
  const nodes = [central]; const edges = [];
  [arm1, arm2].forEach(function(len){
    const chainIds = [];
    for(let i=0;i<len;i++){ const c=mkNode(); nodes.push(c); chainIds.push(c.id); }
    for(let i=0;i<len-1;i++) edges.push({a:chainIds[i], b:chainIds[i+1], type:'S'});
    edges.push({a:central.id, b:chainIds[0], type:'S'});
  });
  return { nodes, edges };
}
// Benzaldehyde -- ring with a group:'CHO' node attached via a ringPos edge,
// same shape ring-attached COOH/CHO/CN groups already use (established
// since the Arenes topic).
function generateBenzaldehyde(){
  const ring = mkRingNode();
  const choC = mkNode(); choC.group = 'CHO';
  return { nodes:[ring, choC], edges:[{a:ring.id, b:choC.id, type:'S', ringPos:0}] };
}

// A 2-connected ether/ester oxygen in generator.js's own id space -- mirrors
// engine.js's newEtherOxygenNode exactly (see generateAmine's mkAmineNode
// for why generator.js keeps its own parallel node builders instead of
// importing engine.js's: staying in this module's own _gid space avoids id
// collisions when a generated molecule is later merged with another one).
function mkEtherOxygenNode(){ return { id:_gid++, subs:[], oxo:false, group:null, phenyl:false, ring:false, element:'O' }; }

// Ethanedioic acid (oxalic acid, HOOC-COOH) -- two group:'COOH' carbons
// bonded directly to each other, no chain in between. Along with methanoic
// acid (generateCarboxylicAcid(1)), this is one of the two acids that
// oxidise all the way to CO2 instead of just resisting oxidation like every
// other carboxylic acid.
function generateEthanedioicAcid(){
  const c1 = mkNode(); c1.group = 'COOH';
  const c2 = mkNode(); c2.group = 'COOH';
  return { nodes:[c1,c2], edges:[{a:c1.id, b:c2.id, type:'S'}] };
}

// An acid anhydride (RCO-O-COR') -- two oxo carbon chains sharing one
// element:'O' link, both sides oxo (unlike an ester, where only one side
// is). Only ever used as a HYDROLYSIS starting material (this app has no
// anhydride-formation question), so there's no matching operator that
// builds this shape -- it's constructed directly here instead.
function generateAcidAnhydride(len1, len2){
  const nodes = []; const edges = [];
  function acidArm(len){
    const chainIds = [];
    for(let i=0;i<len-1;i++){ const c=mkNode(); nodes.push(c); chainIds.push(c.id); }
    for(let i=0;i<chainIds.length-1;i++) edges.push({a:chainIds[i], b:chainIds[i+1], type:'S'});
    const acidC = mkNode(); acidC.oxo = true; nodes.push(acidC);
    if(chainIds.length) edges.push({a:chainIds[chainIds.length-1], b:acidC.id, type:'S'});
    return acidC;
  }
  const c1 = acidArm(len1);
  const c2 = acidArm(len2);
  const o = mkEtherOxygenNode(); nodes.push(o);
  edges.push({a:c1.id, b:o.id, type:'S'});
  edges.push({a:c2.id, b:o.id, type:'S'});
  return { nodes, edges };
}

// Nitrobenzene / phenylamine -- ring with a single flat -NO2/-NH2
// substituent, nothing else. Exact mirror of generatePhenol.
function generateNitrobenzene(){
  const ring = mkRingNode();
  ring.subs[0] = 'NO2';
  return { nodes:[ring], edges:[] };
}
function generatePhenylamine(){
  const ring = mkRingNode();
  ring.subs[0] = 'NH2';
  return { nodes:[ring], edges:[] };
}
// Benzenediazonium chloride -- ring with a closed group:'N2Cl' node
// attached via a ringPos edge, exact mirror of generateBenzaldehyde.
function generateBenzenediazonium(){
  const ring = mkRingNode();
  const dGroup = mkNode(); dGroup.group = 'N2Cl';
  return { nodes:[ring, dGroup], edges:[{a:ring.id, b:dGroup.id, type:'S', ringPos:0}] };
}

module.exports = { generateMolecule, generateArene, generateAlkylbenzene, generateTertButylbenzene, generateAlkylChloride, generateHaloarene, generateSubstitutedArene, generateAmine, generateAlcohol, generatePhenol, generateCarboxylicAcid, generateAcylChloride, generateAldehyde, generateKetone, generateBenzaldehyde, generateEthanedioicAcid, generateAcidAnhydride, generateNitrobenzene, generatePhenylamine, generateBenzenediazonium, ri, pick };
