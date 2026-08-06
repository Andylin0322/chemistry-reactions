'use strict';
const { cloneMol, findNode, neighborsOf, implicitH, canonicalForm, newNode, newAmineNode, newEtherOxygenNode } = require('./engine');

function radicalSubstitution(mol, X){
  const variants = [];
  const seen = new Set();
  mol.nodes.forEach(n=>{
    if(implicitH(mol, n.id) > 0){
      const v = cloneMol(mol);
      findNode(v, n.id).subs.push(X);
      const key = canonicalForm(v);
      if(!seen.has(key)){ seen.add(key); variants.push(v); }
    }
  });
  return { occurs: variants.length>0, variants };
}

function doubleBondEdges(mol){ return mol.edges.filter(e=>e.type==='D'); }

function additionX2(mol, X){
  const dbs = doubleBondEdges(mol);
  if(dbs.length===0) return { occurs:false };
  const v = cloneMol(mol);
  v.edges.filter(e=>e.type==='D').forEach(e=>{
    e.type='S';
    findNode(v,e.a).subs.push(X);
    findNode(v,e.b).subs.push(X);
  });
  return { occurs:true, product: v };
}

function hydrogenation(mol){
  const dbs = doubleBondEdges(mol);
  if(dbs.length===0) return { occurs:false };
  const v = cloneMol(mol);
  v.edges.filter(e=>e.type==='D').forEach(e=>{ e.type='S'; });
  return { occurs:true, product: v };
}

function mildOxidationDiol(mol){
  const dbs = doubleBondEdges(mol);
  if(dbs.length===0) return { occurs:false };
  const v = cloneMol(mol);
  v.edges.filter(e=>e.type==='D').forEach(e=>{
    e.type='S';
    findNode(v,e.a).subs.push('OH');
    findNode(v,e.b).subs.push('OH');
  });
  return { occurs:true, product: v };
}

// Hot acidified KMnO4: break every C=C. Each former alkene carbon's fate depends
// on how much else it was attached to (beyond the double bond itself):
//   0 other attachments -> fully oxidised away as CO2 + H2O
//   1 other attachment  -> becomes -COOH (chain keeps going through that neighbour)
//   2 other attachments -> becomes a ketone carbon (branch point, no H to remove further)
function oxidativeCleavage(mol){
  const dbs = doubleBondEdges(mol);
  if(dbs.length===0) return { occurs:false };
  const v = cloneMol(mol);
  const dbEdgesInV = v.edges.filter(e=>e.type==='D');

  const fate = {}; // nodeId -> 'CO2H2O' | 'COOH' | 'KETONE'
  dbEdgesInV.forEach(e=>{
    [e.a, e.b].forEach(id=>{
      const otherCount = neighborsOf(v, id).filter(nb=> nb.edge!==e ).length; // other bonded neighbours (chain), excluding this double bond
      const node = findNode(v, id);
      const subCount = node.subs.length; // should be 0 by generation constraint, but handle generally
      const k = otherCount + subCount;
      if(k<=0) fate[id]='CO2H2O';
      else if(k===1) fate[id]='COOH';
      else fate[id]='KETONE';
    });
  });

  Object.keys(fate).forEach(idStr=>{
    const id = Number(idStr);
    const node = findNode(v, id);
    if(fate[id]==='COOH') node.group='COOH';
    else if(fate[id]==='KETONE') node.oxo = true;
    // CO2H2O nodes are removed below
  });

  // remove double bond edges
  v.edges = v.edges.filter(e=>e.type!=='D');
  // remove CO2H2O nodes (they have no other edges, safe to drop)
  const removedIds = new Set(Object.keys(fate).filter(id=>fate[id]==='CO2H2O').map(Number));
  v.nodes = v.nodes.filter(n=>!removedIds.has(n.id));
  v.edges = v.edges.filter(e=>!removedIds.has(e.a) && !removedIds.has(e.b));

  // connected components -> fragments
  const idToNode = {}; v.nodes.forEach(n=>idToNode[n.id]=n);
  const adj = {}; v.nodes.forEach(n=>adj[n.id]=[]);
  v.edges.forEach(e=>{ adj[e.a].push(e.b); adj[e.b].push(e.a); });
  const visited = new Set();
  const fragments = [];
  v.nodes.forEach(n=>{
    if(visited.has(n.id)) return;
    const compIds=[]; const stack=[n.id];
    while(stack.length){
      const cur=stack.pop();
      if(visited.has(cur)) continue;
      visited.add(cur); compIds.push(cur);
      adj[cur].forEach(nb=>{ if(!visited.has(nb)) stack.push(nb); });
    }
    const compIdSet = new Set(compIds);
    const compEdges = v.edges.filter(e=>compIdSet.has(e.a) && compIdSet.has(e.b));
    fragments.push({ kind:'chain', nodes: compIds.map(id=>idToNode[id]), edges: compEdges });
  });

  const products = fragments.map(f=>({kind:'chain', mol:{nodes:f.nodes, edges:f.edges}}));
  if(removedIds.size>0) products.push({kind:'small', species:['CO2','H2O']});

  return { occurs:true, products };
}

/* =========================================================================
   ARENE OPERATORS
   These act on a ring:true node's own positions (see engine.js) rather
   than on ordinary carbon valence, so they're a genuinely separate family
   from everything above -- addition/hydrogenation/cleavage all key off
   C=C, which a ring node never has.
========================================================================= */

// A submitted molecule with a ring that ends up having exactly ONE
// occupied position and nothing else folds down to the old opaque
// phenyl:true flag during grading (see app.html's
// __normalizeMoleculeForGrading) -- a correct answer computed here has to
// match that exact shape to compare equal, or a genuinely correct
// submission would silently fail to match. Only a flat single substituent
// directly on the ring (no equivalent in the old model at all) is left as
// a real ring:true node; that shape was never producible by any existing
// Alkane/Alkene operator anyway, so there's nothing it needs to match.
function downgradeRingIfMono(mol){
  const ring = mol.nodes.find(n=>n.ring);
  if(!ring) return mol;
  if(ring.subs.some(s=>s)) return mol;
  const ringEdges = mol.edges.filter(e=>e.a===ring.id || e.b===ring.id);
  if(ringEdges.length!==1) return mol;
  const v = cloneMol(mol);
  const vEdge = v.edges.find(e=>e.a===ring.id || e.b===ring.id);
  const outsideId = vEdge.a===ring.id ? vEdge.b : vEdge.a;
  v.nodes = v.nodes.filter(n=>n.id!==ring.id);
  v.edges = v.edges.filter(e=>e!==vEdge);
  findNode(v, outsideId).phenyl = true;
  return v;
}

// What kind of director sits at ring position `pos` (null if unoccupied):
// 'op' (ortho/para director -- alkyl, OH, NH2, or a halogen, whether given
// as a flat sub or a real chain attached via a ringPos edge) or 'm' (meta
// director -- NO2 flat, or a COOH/CHO/CN group node attached via a ringPos
// edge). Every op-director in this model is ortho/para-directing regardless
// of whether it's activating (alkyl/OH/NH2) or deactivating-but-still-op
// (halogens, matching the Halogen Derivatives topic's haloarene EAS) --
// this app doesn't model relative directing STRENGTH, only which positions
// are favoured at all, so that distinction doesn't need to be tracked here.
const META_DIRECTOR_FLAT = new Set(['NO2']);
const META_DIRECTOR_GROUPS = new Set(['COOH','CHO','CN']);
function ringPositionDirectingType(mol, ring, pos){
  const flat = ring.subs[pos];
  if(flat) return META_DIRECTOR_FLAT.has(flat) ? 'm' : 'op';
  const edge = mol.edges.find(e=>(e.a===ring.id||e.b===ring.id) && e.ringPos===pos);
  if(!edge) return null;
  const outside = findNode(mol, edge.a===ring.id ? edge.b : edge.a);
  if(outside && outside.group && META_DIRECTOR_GROUPS.has(outside.group)) return 'm';
  return 'op';
}
function ringPositionRelation(p, sub){
  const d = Math.min((p-sub+6)%6, (sub-p+6)%6);
  return d===1 ? 'ortho' : d===2 ? 'meta' : 'para';
}
// A position is "favoured" by an existing substituent if the substituent's
// directing type actually points electron density there (ortho/para for an
// op-director, meta for a meta-director).
function isFavoredBy(mol, ring, pos, subPos){
  const rel = ringPositionRelation(pos, subPos);
  const type = ringPositionDirectingType(mol, ring, subPos);
  return type==='op' ? (rel==='ortho'||rel==='para') : rel==='meta';
}

// Candidate ring position(s) for a NEW substituent: every position if the
// ring is currently unsubstituted (all equivalent by symmetry -- any one
// choice gives the same molecule); with 1 or 2 existing substituents, every
// open position favoured by AT LEAST ONE of them (matching how e.g.
// p-xylene brominates cleanly at the position ortho to one methyl and meta
// to the other -- real EAS doesn't require every existing group to agree,
// just that the position isn't disfavoured by all of them at once). Raw
// candidates may include symmetry-duplicates (canonicalForm-deduped by the
// callers below), so there's no need to hand-pick a minimal set here. Three
// or more existing substituents is a real directing-conflict-resolution
// problem this app doesn't model, so that case is refused (null) rather
// than guessed, same as before.
function ringSubstitutionPositions(mol){
  const ring = mol.nodes.find(n=>n.ring);
  if(!ring) return null;
  const occupied = [];
  for(let i=0;i<6;i++){ if(ringPositionDirectingType(mol, ring, i)) occupied.push(i); }
  if(occupied.length===0) return [0];
  if(occupied.length>2) return null;
  const open = [0,1,2,3,4,5].filter(p=>occupied.indexOf(p)===-1);
  const candidates = open.filter(p=>occupied.some(sub=>isFavoredBy(mol, ring, p, sub)));
  return candidates.length ? candidates : null;
}

// Dedupes a list of candidate-position products down to their distinct
// canonicalForm keys, mirroring radicalSubstitution/eliminationHX's own
// enumerate-then-dedup pattern -- needed here because ringSubstitutionPositions
// can return raw positions that are related by the ring's own symmetry (e.g.
// both ortho positions relative to a single existing substituent).
function dedupeByCanonicalForm(products){
  const seen = new Set(); const out = [];
  products.forEach(v=>{ const key=canonicalForm(v); if(!seen.has(key)){ seen.add(key); out.push(v); } });
  return out;
}

// Electrophilic aromatic substitution with a flat single-atom/group
// electrophile (Cl, Br, NO2, ...) -- covers benzene halogenation,
// nitration, and the same reactions on a substituted ring (ortho/para or
// meta to whatever's already there, per ringSubstitutionPositions above).
function ringElectrophilicSubstitution(mol, sub){
  const positions = ringSubstitutionPositions(mol);
  if(!positions) return { occurs:false };
  const results = dedupeByCanonicalForm(positions.map(pos=>{
    const v = cloneMol(mol);
    v.nodes.find(n=>n.ring).subs[pos] = sub;
    return downgradeRingIfMono(v);
  }));
  return results.length===1 ? { occurs:true, product:results[0] } : { occurs:true, variants:results };
}

// Friedel-Crafts alkylation -- like ringElectrophilicSubstitution, but the
// electrophile is a whole fresh plain alkyl chain (chainLen carbons)
// rather than a flat label, so it needs a real subtree + a ringPos-tagged
// edge instead of a subs[] entry.
function ringFriedelCraftsAlkylation(mol, chainLen){
  const positions = ringSubstitutionPositions(mol);
  if(!positions) return { occurs:false };
  const results = dedupeByCanonicalForm(positions.map(pos=>{
    const v = cloneMol(mol);
    const vRing = v.nodes.find(n=>n.ring);
    const chainIds = [];
    for(let i=0;i<chainLen;i++){ const n=newNode(); v.nodes.push(n); chainIds.push(n.id); }
    for(let i=0;i<chainLen-1;i++) v.edges.push({a:chainIds[i], b:chainIds[i+1], type:'S'});
    v.edges.push({a:vRing.id, b:chainIds[0], type:'S', ringPos:pos});
    return downgradeRingIfMono(v);
  }));
  return results.length===1 ? { occurs:true, product:results[0] } : { occurs:true, variants:results };
}

// Hot KMnO4/H2SO4 (or KMnO4/NaOH then acidify) oxidises an alkylbenzene's
// WHOLE side chain down to a single -COOH directly on the ring, regardless
// of how long the chain is -- but only if the benzylic carbon (the one
// bonded to the ring) still has at least one H to lose; a quaternary
// benzylic carbon (tert-butylbenzene) has none and doesn't react. Only
// ever needs the old phenyl:true representation (the ring itself is a
// pure spectator here, never inspected or modified), so there's no
// downgrade step to worry about.
function sideChainOxidationToBenzoicAcid(mol){
  const benzylic = mol.nodes.find(n=>n.phenyl);
  if(!benzylic) return { occurs:false };
  if(implicitH(mol, benzylic.id)<=0) return { occurs:false };
  const v = cloneMol(mol);
  const vBenzylic = findNode(v, benzylic.id);
  const toRemove = new Set();
  (function dfs(id, cameFrom){
    neighborsOf(v, id).forEach(nb=>{
      if(nb.to===cameFrom) return;
      toRemove.add(nb.to);
      dfs(nb.to, id);
    });
  })(benzylic.id, null);
  v.nodes = v.nodes.filter(n=>!toRemove.has(n.id));
  v.edges = v.edges.filter(e=>!toRemove.has(e.a) && !toRemove.has(e.b));
  // vBenzylic.phenyl stays true -- the ring itself is untouched by this
  // reaction, only the side chain is; a node's group and phenyl flags are
  // independent (the carboxyl carbon is directly ring-attached in the
  // product, same as it was before oxidation).
  vBenzylic.subs = [];
  vBenzylic.oxo = false;
  vBenzylic.group = 'COOH';
  return { occurs:true, product:v };
}

/* =========================================================================
   HALOGEN DERIVATIVE OPERATORS
   All of these key off a halogen sitting in some node's `subs` array on a
   PLAIN CARBON -- never a ring node's flat substituent, which is what
   correctly makes a haloarene (Cl/Br/I directly on a ring position) not
   react here at all, matching real chemistry (aryl halides don't undergo
   nucleophilic substitution or elimination under these conditions) without
   any special-casing: a ring's `subs` array only ever gets inspected by
   the ring-specific arene operators above, never by these.
========================================================================= */

// R-X + Nu- -> R-Nu + X- for a simple FLAT substituent nucleophile (OH-,
// NH3) -- the halogen on some plain-carbon node is just swapped for the
// new substituent in place. Used for hydrolysis (X->OH) and formation of
// a primary amine (X->NH2).
function nucleophilicSubstitutionFlat(mol, leavingX, newSub){
  const v = cloneMol(mol);
  let found = false;
  for(const n of v.nodes){
    if(n.ring) continue;
    const idx = n.subs.indexOf(leavingX);
    if(idx!==-1){ n.subs.splice(idx, 1, newSub); found = true; break; }
  }
  if(!found) return { occurs:false };
  return { occurs:true, product:v };
}

// R-X + CN- -> R-CN + X- -- unlike the flat case above, the nitrile carbon
// is a NEW carbon (not a relabeling of the halogen-bearing one), so this
// removes the halogen and attaches a fresh group:'CN' node in its place.
function nitrileFormation(mol, leavingX){
  const v = cloneMol(mol);
  let carbonId = null;
  for(const n of v.nodes){
    if(n.ring) continue;
    const idx = n.subs.indexOf(leavingX);
    if(idx!==-1){ n.subs.splice(idx, 1); carbonId = n.id; break; }
  }
  if(carbonId==null) return { occurs:false };
  const cn = newNode();
  cn.group = 'CN';
  v.nodes.push(cn);
  v.edges.push({ a:carbonId, b:cn.id, type:'S' });
  return { occurs:true, product:v };
}

// R-CN + H2O -> R-COOH (acidic or alkaline conditions both end up at the
// same organic product for grading purposes here, matching how this app
// never requires inorganic byproducts like HX/H2 to be built for any
// other substitution/addition reaction -- only CO2+H2O gets its own
// toggle, and only because oxidative cleavage can produce it as the WHOLE
// answer with no organic product at all, which never happens here).
function nitrileHydrolysis(mol){
  const v = cloneMol(mol);
  const cn = v.nodes.find(n=>n.group==='CN');
  if(!cn) return { occurs:false };
  cn.group = 'COOH';
  return { occurs:true, product:v };
}

// R-CN + 4[H] -> R-CH2-NH2 -- the nitrile carbon survives as a real carbon
// (no longer a closed group node) carrying an NH2 substituent; its other 2
// valence slots are implicit H, matching -CH2-NH2 exactly.
function nitrileReduction(mol){
  const v = cloneMol(mol);
  const cn = v.nodes.find(n=>n.group==='CN');
  if(!cn) return { occurs:false };
  cn.group = null;
  cn.subs = ['NH2'];
  return { occurs:true, product:v };
}

// RCH2-CHX-R' + NaOH(ethanol) -> RCH=CH-R' + HX -- one variant per distinct
// beta-carbon (a real neighbour of the halogen-bearing carbon with at
// least one H to lose), mirroring radicalSubstitution's own
// enumerate-every-position-then-dedup-by-canonicalForm pattern. A halogen
// with no beta-H anywhere (no adjacent carbon has a free H) simply doesn't
// eliminate, same as the real E2 requirement.
function eliminationHX(mol, leavingX){
  const variants = [];
  const seen = new Set();
  mol.nodes.forEach(n=>{
    if(n.ring) return;
    const idx = n.subs.indexOf(leavingX);
    if(idx===-1) return;
    neighborsOf(mol, n.id).forEach(nb=>{
      if(implicitH(mol, nb.to) <= 0) return;
      const v = cloneMol(mol);
      const vn = findNode(v, n.id);
      vn.subs.splice(vn.subs.indexOf(leavingX), 1);
      const ve = v.edges.find(e => (e.a===n.id && e.b===nb.to) || (e.b===n.id && e.a===nb.to));
      ve.type = 'D';
      const key = canonicalForm(v);
      if(!seen.has(key)){ seen.add(key); variants.push(v); }
    });
  });
  return { occurs: variants.length>0, variants };
}

// R2NH/R3N + R'-X -> the next amine up (secondary -> tertiary -> quaternary
// ammonium) by grafting a fresh alkyl group from `reactantMol` onto `mol`'s
// nitrogen centre. `mol` may carry its existing amine as EITHER a flat 'NH2'
// substituent (the primary-amine shape every other operator in this topic
// already produces, e.g. nucleophilicSubstitutionFlat(...,'NH2')) or a real
// newAmineNode -- a flat NH2 is transparently "promoted" to a real node
// first (same molecule, just given a real bondable centre), so a primary
// amine built either way reacts identically here. Degree 3->4 sets
// charge:1 (quaternary ammonium); degree>=4 already (nothing left to
// alkylate) refuses, same as a halogen with no beta-H refusing elimination.
function alkylateAmine(mol, reactantMol, leavingX){
  const v = cloneMol(mol);
  let amine = v.nodes.find(n=>n.element==='N');
  let degree;
  if(amine){
    degree = neighborsOf(v, amine.id).length;
  } else {
    const host = v.nodes.find(n=>!n.ring && n.subs && n.subs.includes('NH2'));
    if(!host) return { occurs:false };
    host.subs = host.subs.filter(s=>s!=='NH2');
    amine = newAmineNode(false);
    v.nodes.push(amine);
    v.edges.push({ a:host.id, b:amine.id, type:'S' });
    degree = 1;
  }
  if(degree>=4) return { occurs:false };
  const halNode = reactantMol.nodes.find(n=>!n.ring && n.subs && n.subs.includes(leavingX));
  if(!halNode) return { occurs:false };
  const rClone = cloneMol(reactantMol);
  const rHal = findNode(rClone, halNode.id);
  rHal.subs = rHal.subs.filter(s=>s!==leavingX);
  v.nodes.push(...rClone.nodes);
  v.edges.push(...rClone.edges);
  v.edges.push({ a:amine.id, b:rHal.id, type:'S' });
  amine.charge = (degree+1===4) ? 1 : 0;
  return { occurs:true, product:v };
}

/* =========================================================================
   HYDROXY COMPOUND OPERATORS
   nucleophilicSubstitutionFlat (X->OH, X->NH2) and eliminationHX (X-beta-H
   -> alkene) are already fully generic on the leaving-group string, so
   ROH -> RCl/RBr/RI/RONa and dehydration (ROH -> alkene) are just direct
   calls with 'OH' as the leaving group -- no new operator code needed for
   those. Only the shapes below (oxidation, esterification/acylation, the
   ring-attached-OH cases, phenol tribromination, combustion) are new.
========================================================================= */

// R-CH2OH (primary) -> R-CHO or R-COOH, R2CH-OH (secondary) -> R2C=O
// (ketone), R3C-OH (tertiary) -> no reaction (no H left on the carbinol
// carbon to remove). target: 'CHO' | 'COOH' | 'ketone' -- callers generate
// the matching alcohol CLASS already, so the class checks below are mostly
// defensive, but they're what makes e.g. calling this with target:'COOH' on
// a tertiary alcohol correctly refuse (used for the
// "tertiary alcohols don't oxidise" no-reaction question, reusing this
// same operator rather than a separate one).
function oxidizeAlcohol(mol, target){
  const v = cloneMol(mol);
  const c = v.nodes.find(n=>!n.ring && n.subs && n.subs.includes('OH'));
  if(!c) return { occurs:false };
  const otherCarbons = neighborsOf(v, c.id).length;
  if(target==='ketone'){
    if(otherCarbons!==2) return { occurs:false }; // must be secondary
    c.subs = c.subs.filter(s=>s!=='OH');
    c.oxo = true;
  } else {
    if(otherCarbons>=2) return { occurs:false }; // must be primary (0 or 1 other carbon)
    c.subs = c.subs.filter(s=>s!=='OH');
    c.group = target; // 'CHO' or 'COOH'
  }
  return { occurs:true, product:v };
}

// Shared merge step for both esterification (from a group:'COOH' acid) and
// acylation (from an oxo+Cl acyl chloride carbon): the carbonyl carbon gets
// a new single-bonded edge to a fresh element:'O' node (newEtherOxygenNode),
// which in turn bonds to whatever supplies the nucleophile's oxygen -- a
// plain carbon's -OH (alcohol) or a ring position's flat OH (phenol),
// either way stripped of its OH first. Mirrors alkylateAmine's "graft a
// whole second molecule onto one specific atom" pattern, generalized to
// also accept a ring nucleophile the way Friedel-Crafts already does.
function esterifyCommon(carbonylMol, carbonylId, fromGroup, nucleophileMol){
  const v = cloneMol(carbonylMol);
  const vC = findNode(v, carbonylId);
  if(fromGroup){ vC.group = null; vC.oxo = true; }
  else { vC.subs = vC.subs.filter(s=>s!=='Cl'); }
  const o = newEtherOxygenNode();
  v.nodes.push(o);
  v.edges.push({ a:vC.id, b:o.id, type:'S' });

  const nClone = cloneMol(nucleophileMol);
  const ring = nClone.nodes.find(n=>n.ring);
  if(ring){
    const pos = ring.subs.findIndex(s=>s==='OH');
    if(pos===-1) return { occurs:false };
    ring.subs[pos] = null;
    v.nodes.push(...nClone.nodes);
    v.edges.push(...nClone.edges);
    v.edges.push({ a:o.id, b:ring.id, type:'S', ringPos:pos });
  } else {
    const c = nClone.nodes.find(n=>!n.ring && n.subs && n.subs.includes('OH'));
    if(!c) return { occurs:false };
    c.subs = c.subs.filter(s=>s!=='OH');
    v.nodes.push(...nClone.nodes);
    v.edges.push(...nClone.edges);
    v.edges.push({ a:o.id, b:c.id, type:'S' });
  }
  return { occurs:true, product:v };
}
// RCOOH + R'OH -> RCOOR' + H2O -- phenols deliberately refused here (a
// carboxylic acid doesn't esterify a phenol directly under these
// conditions, unlike an acyl chloride below), matching the source data.
function esterifyAcid(acidMol, alcoholMol){
  if(alcoholMol.nodes.some(n=>n.ring)) return { occurs:false };
  const acidC = acidMol.nodes.find(n=>n.group==='COOH');
  if(!acidC) return { occurs:false };
  return esterifyCommon(acidMol, acidC.id, true, alcoholMol);
}
// RCOCl + R'OH -> RCOOR' + HCl, or RCOCl + ArOH -> RCOOAr + HCl -- unlike
// esterifyAcid, phenols DO react here.
function esterifyAcylChloride(acylMol, nucleophileMol){
  const acylC = acylMol.nodes.find(n=>!n.ring && n.oxo && n.subs && n.subs.includes('Cl'));
  if(!acylC) return { occurs:false };
  return esterifyCommon(acylMol, acylC.id, false, nucleophileMol);
}

// Swaps a substituent sitting DIRECTLY on a ring position (ArOH -> ArONa,
// i.e. phenol + Na/NaOH). Deliberately separate from
// nucleophilicSubstitutionFlat, which skips ring nodes on purpose (a
// ring-attached HALOGEN never reacts that way) -- a ring-attached OH
// genuinely does react here, so this is the ring-scoped equivalent.
function ringFlatSubSwap(mol, leavingSub, newSub){
  const ring = mol.nodes.find(n=>n.ring);
  if(!ring) return { occurs:false };
  const pos = ring.subs.findIndex(s=>s===leavingSub);
  if(pos===-1) return { occurs:false };
  const v = cloneMol(mol);
  v.nodes.find(n=>n.ring).subs[pos] = newSub;
  return { occurs:true, product:v };
}

// ArOH + excess Br2(aq) -> 2,4,6-tribromophenol -- unlike controlled
// Br2(l)/Br2-in-CCl4 conditions (which give a single monobromination
// product via the ordinary ringElectrophilicSubstitution, ortho/para like
// any other op-director), excess aqueous bromine is strongly activated
// enough by -OH to substitute ALL THREE non-meta positions at once, with no
// "pick one position" choice for the student -- a fixed single answer, not
// a variants list.
// Generalized from phenol specifically -- ArOH + excess Br2(aq) and
// ArNH2 + excess Br2(aq) are the same shape reaction (a strong op-director
// activating all three non-meta positions at once), just keyed on a
// different existing flat substituent.
function ringTribromination(mol, existingSub){
  const ring = mol.nodes.find(n=>n.ring);
  if(!ring) return { occurs:false };
  const subPos = ring.subs.findIndex(s=>s===existingSub);
  if(subPos===-1) return { occurs:false };
  const v = cloneMol(mol);
  const vRing = v.nodes.find(n=>n.ring);
  [1,3,5].forEach(offset=>{ vRing.subs[(subPos+offset)%6] = 'Br'; }); // both ortho + para -- everything but the two meta positions
  return { occurs:true, product:v };
}

// ROH + O2 -> CO2 + H2O -- no organic product at all, same "products" shape
// oxidativeCleavage already uses for its own CO2+H2O byproduct, just as the
// WHOLE answer instead of one fragment among several.
function combustion(mol){
  return { occurs:true, products:[{ kind:'small', species:['CO2','H2O'] }] };
}

/* =========================================================================
   CARBONYL COMPOUND OPERATORS
   Every product here is a shape the grading pipeline already understands
   from earlier topics -- a flat -OH sub, a group:'CN' node on a real edge,
   or a plain oxo/group reassignment -- so none of these needed any new
   engine.js node kind or app.html grading branch, unlike the amine/ester
   work.
========================================================================= */

// Finds the carbonyl carbon in a molecule built by generateAldehyde,
// generateKetone, or generateBenzaldehyde -- a group:'CHO' node or a plain
// carbon with oxo:true, whichever is present. Both aldehyde and ketone
// starting materials share this one lookup, since the operators below all
// act on "whatever the carbonyl carbon is" identically regardless of class.
function findCarbonylCarbon(mol){
  return mol.nodes.find(n=>n.group==='CHO' || n.oxo);
}

// RCHO/RCOR' + HCN -> a cyanohydrin (R2C(OH)CN) -- the carbonyl carbon
// loses its group/oxo closure, gains a flat -OH sub, and gains a real edge
// to a fresh group:'CN' node (mirrors nitrileFormation's own "attach a
// fresh CN node" step). Works identically for an aldehyde or a ketone
// carbonyl carbon -- whichever findCarbonylCarbon finds.
function cyanohydrinFormation(mol){
  const v = cloneMol(mol);
  const c = findCarbonylCarbon(v);
  if(!c) return { occurs:false };
  c.group = null; c.oxo = false;
  c.subs.push('OH');
  const cn = newNode(); cn.group = 'CN';
  v.nodes.push(cn);
  v.edges.push({ a:c.id, b:cn.id, type:'S' });
  return { occurs:true, product:v };
}

// RCHO + 2[H] -> RCH2OH, or RCOR' + 2[H] -> RCH(OH)R' -- the carbonyl
// carbon just swaps its group/oxo closure for a flat -OH sub, same
// reassignment-only pattern as nitrileReduction. One operator covers both
// classes since the shape of the swap is identical either way.
function reduceCarbonyl(mol){
  const v = cloneMol(mol);
  const c = findCarbonylCarbon(v);
  if(!c) return { occurs:false };
  c.group = null; c.oxo = false;
  c.subs.push('OH');
  return { occurs:true, product:v };
}

// RCHO + [O] -> RCOOH -- also doubles as every "ketones/benzaldehyde don't
// give this test" no-reaction question in this topic: a ketone has no
// group:'CHO' node to find at all (oxidation only ever touches the CHO
// shape, unlike reduceCarbonyl/cyanohydrinFormation above which accept
// oxo too), so calling this on a ketone-only molecule correctly refuses
// with no special-casing needed. Covers the shared oxidising-agent
// chemistry behind K2Cr2O7/KMnO4, Tollens', and Fehling's alike -- this
// app doesn't track which inorganic oxidant/observation goes with which,
// only whether the underlying RCHO->RCOOH transform occurs.
function oxidizeAldehyde(mol){
  const v = cloneMol(mol);
  const c = v.nodes.find(n=>n.group==='CHO');
  if(!c) return { occurs:false };
  c.group = 'COOH';
  return { occurs:true, product:v };
}

/* =========================================================================
   CARBOXYLIC ACID & DERIVATIVE OPERATORS
   'COONa' is a new group value alongside 'COOH'/'CHO'/'CN' -- just another
   string in engine.js's existing group field, so it needed NO engine.js
   changes at all. Secondary/tertiary amide formation reuses newAmineNode
   exactly as alkylateAmine's amine centre does -- an amide nitrogen and an
   amine nitrogen are structurally identical in this graph-only model (both
   just a trivalent N with some real carbon neighbours); the resonance/
   planarity distinction real chemistry cares about isn't tracked by this
   engine for ordinary amines either, so there's nothing new to represent.
========================================================================= */

// RCOOH + Na/NaOH/Na2CO3/NaHCO3 -> RCOONa -- one reassignment covers every
// "carboxylic acid meets a base" question in this topic; they only differ
// in reagentText, not in what happens to the molecule.
function carboxylicAcidSaltFormation(mol){
  const v = cloneMol(mol);
  const c = v.nodes.find(n=>n.group==='COOH');
  if(!c) return { occurs:false };
  c.group = 'COONa';
  return { occurs:true, product:v };
}

// RCOOH -> RCOCl (PCl5, PCl3, or SOCl2 -- all three give the same organic
// product, only the reagentText differs) -- the acid carbon stops being a
// closed group node and becomes a plain oxo carbon with a flat Cl sub,
// exactly matching generateAcylChloride's own shape.
function acidToAcylChloride(mol){
  const v = cloneMol(mol);
  const c = v.nodes.find(n=>n.group==='COOH');
  if(!c) return { occurs:false };
  c.group = null;
  c.oxo = true;
  c.subs.push('Cl');
  return { occurs:true, product:v };
}

// RCOOH + LiAlH4 -> RCH2OH -- LiAlH4 is the only reagent in this topic that
// actually reduces a carboxylic acid (NaBH4 and H2/Ni both refuse, covered
// by fixed no-reaction POOL entries rather than this operator).
function reduceCarboxylicAcid(mol){
  const v = cloneMol(mol);
  const c = v.nodes.find(n=>n.group==='COOH');
  if(!c) return { occurs:false };
  c.group = null;
  c.subs.push('OH');
  return { occurs:true, product:v };
}

// RCOCl + H2O -> RCOOH -- reverse of acidToAcylChloride.
function hydrolyzeAcylChloride(mol){
  const v = cloneMol(mol);
  const c = v.nodes.find(n=>!n.ring && n.oxo && n.subs && n.subs.includes('Cl'));
  if(!c) return { occurs:false };
  c.subs = c.subs.filter(s=>s!=='Cl');
  c.oxo = false;
  c.group = 'COOH';
  return { occurs:true, product:v };
}

// RCOCl + amine -> amide + HCl (untracked, same as every other HX
// byproduct in this app). amineMol null means NH3 -> a flat 'NH2' sub
// (primary amide, RCONH2); a real amine molecule (flat-NH2 primary amine,
// promoted exactly like alkylateAmine does, or a real secondary-amine N
// node) grafts its nitrogen directly onto the carbonyl carbon instead,
// giving a secondary (RCONHR') or tertiary (RCONR'2) amide.
function acylChlorideToAmide(acylMol, amineMol){
  const v = cloneMol(acylMol);
  const c = v.nodes.find(n=>!n.ring && n.oxo && n.subs && n.subs.includes('Cl'));
  if(!c) return { occurs:false };
  c.subs = c.subs.filter(s=>s!=='Cl');
  if(!amineMol){
    c.subs.push('NH2');
    return { occurs:true, product:v };
  }
  const aClone = cloneMol(amineMol);
  let amine = aClone.nodes.find(n=>n.element==='N');
  if(!amine){
    const host = aClone.nodes.find(n=>!n.ring && n.subs && n.subs.includes('NH2'));
    if(!host) return { occurs:false };
    host.subs = host.subs.filter(s=>s!=='NH2');
    amine = newAmineNode(false);
    aClone.nodes.push(amine);
    aClone.edges.push({ a:host.id, b:amine.id, type:'S' });
  }
  v.nodes.push(...aClone.nodes);
  v.edges.push(...aClone.edges);
  v.edges.push({ a:c.id, b:amine.id, type:'S' });
  return { occurs:true, product:v };
}

// Splits a molecule into its connected-component fragments, in the same
// {nodes,edges} shape oxidativeCleavage's own fragment-finder produces --
// extracted here as a shared helper since ester and acid-anhydride
// hydrolysis both need it too now.
function connectedComponents(mol){
  const idToNode = {}; mol.nodes.forEach(n=>idToNode[n.id]=n);
  const adj = {}; mol.nodes.forEach(n=>adj[n.id]=[]);
  mol.edges.forEach(e=>{ adj[e.a].push(e.b); adj[e.b].push(e.a); });
  const visited = new Set(); const fragments = [];
  mol.nodes.forEach(n=>{
    if(visited.has(n.id)) return;
    const compIds=[]; const stack=[n.id];
    while(stack.length){
      const cur = stack.pop();
      if(visited.has(cur)) continue;
      visited.add(cur); compIds.push(cur);
      adj[cur].forEach(nb=>{ if(!visited.has(nb)) stack.push(nb); });
    }
    const compIdSet = new Set(compIds);
    const compEdges = mol.edges.filter(e=>compIdSet.has(e.a) && compIdSet.has(e.b));
    fragments.push({ nodes: compIds.map(id=>idToNode[id]), edges: compEdges });
  });
  return fragments;
}

// RCOOR' + H2O -> RCOOH + R'OH (acidic, target:'COOH') or RCOO- + R'OH
// (alkaline, target:'COONa') -- reverse of esterifyAcid/esterifyAcylChloride:
// removes the element:'O' ester link, restores the carbonyl side to an acid
// (or its salt) and the alcohol side to a flat -OH (a plain carbon's sub,
// or a ring position's flat sub for an aryl ester), then splits into two
// separate fragments the same way oxidativeCleavage already does.
function hydrolyzeEster(mol, target){
  const v = cloneMol(mol);
  const o = v.nodes.find(n=>n.element==='O');
  if(!o) return { occurs:false };
  const oNbs = neighborsOf(v, o.id);
  if(oNbs.length!==2) return { occurs:false };
  const carbonylNb = oNbs.find(nb=>{ const n=findNode(v, nb.to); return !n.ring && n.oxo; });
  if(!carbonylNb) return { occurs:false };
  const alcoholNb = oNbs.find(nb=>nb!==carbonylNb);
  const carbonylC = findNode(v, carbonylNb.to);
  carbonylC.oxo = false;
  carbonylC.group = target;
  const alcoholTarget = findNode(v, alcoholNb.to);
  if(alcoholTarget.ring){
    const linkEdge = v.edges.find(e=>(e.a===o.id && e.b===alcoholTarget.id) || (e.b===o.id && e.a===alcoholTarget.id));
    alcoholTarget.subs[linkEdge.ringPos] = 'OH';
  } else {
    alcoholTarget.subs.push('OH');
  }
  v.nodes = v.nodes.filter(n=>n.id!==o.id);
  v.edges = v.edges.filter(e=>e.a!==o.id && e.b!==o.id);
  const products = connectedComponents(v).map(f=>({ kind:'chain', mol:{nodes:f.nodes, edges:f.edges} }));
  return { occurs:true, products };
}

// (RCO)2O + H2O -> 2RCOOH -- both carbonyl carbons flanking the shared
// element:'O' node close back up into COOH groups, then the molecule
// splits into its two fragments exactly like ester hydrolysis above.
function hydrolyzeAcidAnhydride(mol){
  const v = cloneMol(mol);
  const o = v.nodes.find(n=>n.element==='O');
  if(!o) return { occurs:false };
  const oNbs = neighborsOf(v, o.id);
  if(oNbs.length!==2 || !oNbs.every(nb=>{ const n=findNode(v, nb.to); return !n.ring && n.oxo; })) return { occurs:false };
  oNbs.forEach(nb=>{ const c = findNode(v, nb.to); c.oxo = false; c.group = 'COOH'; });
  v.nodes = v.nodes.filter(n=>n.id!==o.id);
  v.edges = v.edges.filter(e=>e.a!==o.id && e.b!==o.id);
  const products = connectedComponents(v).map(f=>({ kind:'chain', mol:{nodes:f.nodes, edges:f.edges} }));
  return { occurs:true, products };
}

// RCONH2 + H2O -> RCOOH (acidic, target:'COOH') or RCOO- (alkaline,
// target:'COONa') -- primary amide only (a flat 'NH2' sub on an oxo
// carbon); secondary/tertiary (N-substituted) amide hydrolysis would need
// the same fragment-splitting treatment as ester hydrolysis and isn't
// covered by a separate entry in the source data, so it's out of scope
// here.
function hydrolyzeAmide(mol, target){
  const v = cloneMol(mol);
  const c = v.nodes.find(n=>!n.ring && n.oxo && n.subs && n.subs.includes('NH2'));
  if(!c) return { occurs:false };
  c.subs = c.subs.filter(s=>s!=='NH2');
  c.oxo = false;
  c.group = target;
  return { occurs:true, product:v };
}

/* =========================================================================
   NITROGEN COMPOUND OPERATORS
   The diazonium group (-N2Cl) is a CLOSED group node, same family as
   COOH/CHO/CN/ONa/COONa -- Ar-N=N+ technically needs a triple bond, which
   this engine has never needed to represent (nitriles sidestep the same
   issue the same way), so rather than add triple-bond support engine-wide
   for one topic, the diazonium salt is opaque until azo coupling actually
   needs to react at it -- at which point it's unpacked into two real
   element:'N' nodes joined by a plain double bond (already fully
   supported), which is everything the PRODUCT side of this chemistry
   actually needs.
========================================================================= */

// ArNH2 + HNO2/HCl -> ArN2+Cl- -- the ring's flat -NH2 is replaced by a
// closed -N2Cl group node, attached the same way a ring-position COOH/CHO
// already is.
function diazotisation(mol){
  const ring = mol.nodes.find(n=>n.ring);
  if(!ring) return { occurs:false };
  const pos = ring.subs.findIndex(s=>s==='NH2');
  if(pos===-1) return { occurs:false };
  const v = cloneMol(mol);
  const vRing = v.nodes.find(n=>n.ring);
  vRing.subs[pos] = null;
  const dGroup = newNode(); dGroup.group = 'N2Cl';
  v.nodes.push(dGroup);
  v.edges.push({ a:vRing.id, b:dGroup.id, type:'S', ringPos:pos });
  return { occurs:true, product:v };
}

// ArN2+Cl- + phenol/phenylamine -> an azo dye (Ar-N=N-Ar'), coupling at an
// ortho/para position on the SECOND ring the same way any other EAS
// substitution in this app picks positions (ringSubstitutionPositions,
// keyed off the nucleophile's own existing -OH/-NH2). Unpacks the closed
// -N2Cl group into two real element:'N' nodes (inner one keeps the
// original ring attachment, outer one gets a fresh edge into the new
// ring) joined by a plain double bond -- both nitrogens end up neutral,
// matching the real product (the cationic diazonium nitrogen is quenched
// by the coupling).
function azoCoupling(diazoniumMol, nucleophileMol){
  const dGroup = diazoniumMol.nodes.find(n=>n.group==='N2Cl');
  if(!dGroup) return { occurs:false };
  const positions = ringSubstitutionPositions(nucleophileMol);
  if(!positions) return { occurs:false };
  const results = dedupeByCanonicalForm(positions.map(pos=>{
    const v = cloneMol(diazoniumMol);
    const vGroup = findNode(v, dGroup.id);
    const groupEdge = v.edges.find(e=>e.a===vGroup.id || e.b===vGroup.id);
    const ringSideId = groupEdge.a===vGroup.id ? groupEdge.b : groupEdge.a;
    const innerN = newAmineNode(false);
    const outerN = newAmineNode(false);
    v.nodes = v.nodes.filter(n=>n.id!==vGroup.id);
    v.edges = v.edges.filter(e=>e!==groupEdge);
    v.nodes.push(innerN, outerN);
    const ringEdge = { a:ringSideId, b:innerN.id, type:'S' };
    if(groupEdge.ringPos!=null) ringEdge.ringPos = groupEdge.ringPos;
    v.edges.push(ringEdge);
    v.edges.push({ a:innerN.id, b:outerN.id, type:'D' });

    const nClone = cloneMol(nucleophileMol);
    v.nodes.push(...nClone.nodes);
    v.edges.push(...nClone.edges);
    const nRing = nClone.nodes.find(n=>n.ring);
    v.edges.push({ a:outerN.id, b:nRing.id, type:'S', ringPos:pos });
    return v;
  }));
  return results.length===1 ? { occurs:true, product:results[0] } : { occurs:true, variants:results };
}

// RNH2 + H+ -> RNH3+ (also ArNH2 -> ArNH3+) -- promotes a flat 'NH2' sub
// (plain-carbon OR ring-attached) to a real element:'N' node exactly like
// alkylateAmine's own promotion step, then just sets charge:1. No new
// bonds are added -- nodeValence/freeSlots/implicitH already fill in the
// resulting extra implicit H generically now that charge is a universal
// per-node field, the same machinery quaternary ammonium already relies on.
function protonateAmine(mol){
  const v = cloneMol(mol);
  let amine = v.nodes.find(n=>n.element==='N');
  if(!amine){
    const ring = v.nodes.find(n=>n.ring && n.subs.includes('NH2'));
    if(ring){
      const pos = ring.subs.indexOf('NH2');
      ring.subs[pos] = null;
      amine = newAmineNode(false);
      v.nodes.push(amine);
      v.edges.push({ a:ring.id, b:amine.id, type:'S', ringPos:pos });
    } else {
      const host = v.nodes.find(n=>!n.ring && n.subs && n.subs.includes('NH2'));
      if(!host) return { occurs:false };
      host.subs = host.subs.filter(s=>s!=='NH2');
      amine = newAmineNode(false);
      v.nodes.push(amine);
      v.edges.push({ a:host.id, b:amine.id, type:'S' });
    }
  }
  const degree = neighborsOf(v, amine.id).length;
  if(degree>=4) return { occurs:false };
  amine.charge = 1;
  return { occurs:true, product:v };
}

module.exports = {
  radicalSubstitution, additionX2, hydrogenation, mildOxidationDiol, oxidativeCleavage, doubleBondEdges,
  ringElectrophilicSubstitution, ringFriedelCraftsAlkylation, sideChainOxidationToBenzoicAcid, downgradeRingIfMono,
  nucleophilicSubstitutionFlat, nitrileFormation, nitrileHydrolysis, nitrileReduction, eliminationHX, alkylateAmine,
  oxidizeAlcohol, esterifyAcid, esterifyAcylChloride, ringFlatSubSwap, ringTribromination, combustion,
  cyanohydrinFormation, reduceCarbonyl, oxidizeAldehyde,
  carboxylicAcidSaltFormation, acidToAcylChloride, reduceCarboxylicAcid, hydrolyzeAcylChloride, acylChlorideToAmide,
  hydrolyzeEster, hydrolyzeAcidAnhydride, hydrolyzeAmide, connectedComponents,
  diazotisation, azoCoupling, protonateAmine
};
