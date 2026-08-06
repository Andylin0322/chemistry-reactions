'use strict';
/* =========================================================================
   MOLECULE MODEL
   A molecule is a TREE of carbon nodes (branching allowed, no rings in the
   editable graph, EXCEPT for the dedicated aromatic-ring node kind below).
   The original spectator-only benzene representation is an opaque 'phenyl'
   flag on a node -- never touched by any reaction operator, which is what
   makes it correctly inert without any special-case code. That stays
   exactly as-is (Alkane/Alkene topics never react AT the ring, only ever
   treat it as a spectator substituent).

   Arene-topic reactions, introduced later, DO react at the ring itself, so
   a second, richer representation exists alongside it: a node with
   ring:true is an aromatic ring with 6 fixed, individually addressable
   positions (0-5, aromatic C-H by default) instead of ordinary valence --
   see canonicalizeRingPositions below for how those positions compare for
   canonical-form purposes. A ring node ignores oxo/group entirely
   (meaningless for it, same as how a 'group' node ignores oxo/phenyl), but
   REUSES subs -- for a ring node specifically, subs[pos] (0-5) is a flat
   substituent string (Cl/Br/F/I/OH/NO2 -- the same simple-substituent
   vocabulary a plain carbon already uses) directly at that position,
   exactly parallel to how subs already works for an ordinary carbon.
   Anything that needs real STRUCTURE at a position instead of a flat label
   (an alkyl chain from Friedel-Crafts, or the rest of the molecule the ring
   is attached to, e.g. toluene's CH3) is an ordinary tree edge tagged with
   which position it occupies (edge.ringPos, 0-5) -- reusing the whole
   existing tree/edge machinery instead of a parallel one. A position is
   never both at once (enforced by callers, not here).

   node: { id, subs:[...substituent strings, or 6 ring-position slots for a ring node], oxo:bool, group:null|'COOH'|'CHO'|'CN', phenyl:bool, ring:bool }
   edge: { a:id, b:id, type:'S'|'D', ringPos:null|0-5 }
   molecule: { nodes:[...], edges:[...] }

   Valence bookkeeping: every plain carbon has 4 slots.
     used = sum(bond order to each neighbour) + subs.length + (oxo?2:0) + (group? treated as fully closing the node)
     implicit H = 4 - used   (a 'group' node's H's are implied by the group itself, not tracked separately)
   A ring node has no valence in this sense -- see computeRingFormula-style
   logic in callers instead (each of its 6 positions is either aromatic C-H
   or a real bonded substituent, never a "free slot" to fill).
========================================================================= */

let _id = 0;
function newNode(){ return { id:_id++, subs:[], oxo:false, group:null, phenyl:false, ring:false, charge:0 }; }
function newRingNode(){ return { id:_id++, subs:new Array(6).fill(null), oxo:false, group:null, phenyl:false, ring:true }; }
// A branching/charged nitrogen centre -- secondary/tertiary amine (charge 0,
// 2 or 3 real carbon neighbours) or quaternary ammonium (charge 1, exactly 4
// real carbon neighbours). Deliberately NOT a new architecture: it's a plain
// tree node like any other carbon, just with element:'N' overriding its
// label and valence (see nodeLabel/nodeValence below) -- canonicalForm's
// existing generic branching-tree recursion (rec() in rootedCanonical)
// already handles any number of real children per node with no changes at
// all, so 2-4 carbon arms "just work" the same way a branched alkyl carbon
// already does. A 1-neighbour amine (primary, R-NH2) is NOT represented
// this way -- that stays the flat 'NH2' substituent string on the carbon,
// exactly as before, since it can't bridge two chains and folding it to a
// flat sub keeps every existing primary-amine answer unchanged.
function newAmineNode(charged){ return { id:_id++, subs:[], oxo:false, group:null, phenyl:false, ring:false, element:'N', charge:charged?1:0 }; }
// A 2-connected ether/ester oxygen (R-O-R') -- same "reuse the plain node,
// override label/valence via element" trick as newAmineNode above, not a
// new node kind. Only ever built with exactly 2 real neighbours (an ester's
// carbonyl carbon on one side, the alcohol's R' chain on the other) --
// nothing here enforces that shape structurally, it's just what every
// caller that builds one actually does.
function newEtherOxygenNode(){ return { id:_id++, subs:[], oxo:false, group:null, phenyl:false, ring:false, element:'O', charge:0 }; }
function cloneMol(mol){
  return { nodes: mol.nodes.map(n=>({...n, subs:n.subs.slice()})), edges: mol.edges.map(e=>({...e})) };
}
function findNode(mol, id){ return mol.nodes.find(n=>n.id===id); }
function neighborsOf(mol, id){
  return mol.edges.filter(e=>e.a===id||e.b===id).map(e=>({to: e.a===id?e.b:e.a, type:e.type, edge:e}));
}
function bondSlotsUsed(mol, id){
  return neighborsOf(mol,id).reduce((s,e)=> s + (e.type==='D'?2:1), 0);
}
function subSlotsUsed(node){
  // A ring node has no ordinary valence at all -- it's 6 fixed positions,
  // not a 4-slot atom, and node.subs.length is always 6 for it regardless
  // of how many positions are actually occupied (see newRingNode), which
  // would otherwise silently feed a meaningless number into the classic
  // valence math below. Explicit here rather than relying on that number
  // happening to come out too large to matter (freeSlots clamps to 0
  // either way, but only by accident, not by design).
  if(node.ring) return 0;
  if(node.group) return 4; // terminal group closes the node entirely
  return node.subs.length + (node.oxo?2:0) + (node.phenyl?1:0);
}
// Charge shifts an atom's bonding valence by +-1 per unit, but which
// direction depends on whether the NEUTRAL atom has a spare lone pair to
// promote into a bond: N/O/S/P/halogens do (ammonium NH4+ =3+1=4, alkoxide
// RO- =2-1=1, oxonium H3O+ =2+1=3, halide Cl- =1-1=0 -- all match this
// directly), so charge simply adds. Carbon (and hydrogen) have NO spare
// lone pair in their normal bonding state -- a carbocation (3 bonds, empty
// orbital) and a carbanion (3 bonds, lone pair) both cost a bonding slot
// regardless of sign, so |charge| always subtracts for them. This is a
// completely generic per-node formula now (not just the amine centre) --
// only ever restricted to charge in {-1,0,1} by the builder UI, since
// that's the only chemically meaningful range at this level for every one
// of these elements (each has at most one lone pair worth promoting).
const LONE_PAIR_ELEMENTS = { N:true, O:true, S:true, P:true, F:true, Cl:true, Br:true, I:true };
function neutralValence(node){
  if(node.element==='N') return 3;
  if(node.element==='O') return 2;
  return 4; // plain carbon (the default -- no `element` field set)
}
function nodeValence(node){
  const charge = node.charge||0;
  if(charge===0) return neutralValence(node);
  const el = node.element || 'C';
  return LONE_PAIR_ELEMENTS[el] ? neutralValence(node)+charge : neutralValence(node)-Math.abs(charge);
}
function usedSlots(mol, id){ return bondSlotsUsed(mol,id) + subSlotsUsed(findNode(mol,id)); }
function freeSlots(mol, id){ return nodeValence(findNode(mol,id)) - usedSlots(mol,id); }
function implicitH(mol, id){
  const n=findNode(mol,id);
  if(n.group || n.ring) return 0; // a ring's aromatic C-H count is per-position, not this generic per-atom count; operators that shouldn't touch the ring (radicalSubstitution etc.) rely on this being 0, not a leftover valence number
  return Math.max(0, freeSlots(mol,id));
}

/* ============================= CANONICAL FORM =============================
   Order-independent equality for branching trees: root at the tree's
   centroid(s) (standard leaf-trimming algorithm), compute a canonical
   recursive string from each candidate root, take the lexicographically
   smallest. This generalises the old "compare forward vs reversed chain"
   trick to arbitrary trees.
========================================================================= */
function nodeLabel(n){
  if(n.ring) return 'RING'; // see rec() below -- a ring node's real label is computed via ringPositionLabels/canonicalizeRingPositions instead, this is only a defensive fallback
  if(n.group) return 'G:'+n.group;
  const parts=[]; if(n.oxo) parts.push('O'); if(n.phenyl) parts.push('Ph');
  if(n.element==='N') parts.push('N');
  if(n.element==='O') parts.push('Oe'); // distinct token from the plain oxo '=O' flag's 'O' -- never the same node kind, but kept unambiguous anyway
  // Generic charge marker, independent of element -- covers a charged plain
  // carbon (carbocation/carbanion) exactly the same way as a charged N/O
  // centre, since canonicalForm needs +1 and -1 (and neutral) to compare as
  // genuinely different molecules regardless of what kind of node carries it.
  if(n.charge) parts.push('chg'+n.charge);
  parts.push(...n.subs.slice().sort());
  return parts.join(',');
}
function buildAdj(mol){
  const adj={}; mol.nodes.forEach(n=>adj[n.id]=[]);
  mol.edges.forEach(e=>{
    const ringPos = e.ringPos!=null ? e.ringPos : null;
    adj[e.a].push({to:e.b,type:e.type,ringPos:ringPos});
    adj[e.b].push({to:e.a,type:e.type,ringPos:ringPos});
  });
  return adj;
}
// Canonicalizes a ring's 6 position labels over its full physical symmetry
// -- D6, 12 elements (6 rotations x 2 reflections) -- ALWAYS, whether or not
// one of the positions is the outside attachment point ('@', see rec()
// below). A uniform rotation/reflection of the whole position array (anchor
// included) never changes any position's arrangement RELATIVE to any
// other, which is the only thing that's ever chemically meaningful (ortho/
// meta/para is a statement about relative distance around the ring, not
// about literal array index) -- so the anchor is just one more label to
// canonicalize along with everything else, not a fixed pivot. Only ever
// treating it as a fixed pivot (an earlier version of this function) was a
// real bug: it made a monosubstituted ring compare unequal to itself
// depending purely on which position happened to record the substituent.
function canonicalizeRingPositions(labels){
  const variants = [];
  for(let refl=0; refl<2; refl++){
    const base = refl ? labels.slice().reverse() : labels.slice();
    for(let rot=0; rot<6; rot++) variants.push(base.slice(rot).concat(base.slice(0,rot)).join('|'));
  }
  return variants.sort()[0];
}
function centroids(mol, adj){
  const ids = mol.nodes.map(n=>n.id);
  if(ids.length<=2) return ids.slice();
  const deg={}; ids.forEach(id=>deg[id]=adj[id].length);
  let remaining = new Set(ids);
  let leaves = ids.filter(id=>deg[id]<=1);
  while(remaining.size>2){
    const next=[];
    leaves.forEach(leaf=>{
      if(!remaining.has(leaf)) return;
      remaining.delete(leaf);
      adj[leaf].forEach(e=>{
        if(remaining.has(e.to)){ deg[e.to]--; if(deg[e.to]===1) next.push(e.to); }
      });
    });
    leaves = next;
    if(leaves.length===0) break;
  }
  return [...remaining];
}
function rootedCanonical(mol, rootId, adj){
  const nodeMap={}; mol.nodes.forEach(n=>nodeMap[n.id]=n);
  function rec(id, parent){
    const node = nodeMap[id];
    if(node.ring){
      const labels = new Array(6).fill('H');
      (node.subs||[]).forEach((s,pos)=>{ if(s) labels[pos] = s; }); // flat simple substituents (Cl/Br/OH/NO2/...)
      adj[id].forEach(e=>{
        if(e.ringPos==null) return; // defensive; every real ring edge carries ringPos
        labels[e.ringPos] = (e.to===parent) ? '@' : (e.type + rec(e.to, id)); // real structure (alkyl chain / outside attachment)
      });
      return 'RING['+canonicalizeRingPositions(labels)+']';
    }
    const label = nodeLabel(node);
    const kids = adj[id].filter(e=>e.to!==parent).map(e=> e.type + rec(e.to, id));
    kids.sort();
    return label+'('+kids.join(',')+')';
  }
  return rec(rootId, null);
}
function canonicalForm(mol){
  if(mol.nodes.length===0) return '';
  const adj = buildAdj(mol);
  const cs = centroids(mol, adj);
  const forms = cs.map(c=>rootedCanonical(mol, c, adj));
  return forms.sort()[0];
}

/* ============================= EQUALITY FOR ANSWERS ============================= */
function normSpecies(s){ return s.replace(/^\d+/,'').trim().toLowerCase(); }
function setsEqual(a,b){
  const sa=new Set(a.map(x=>x.toLowerCase().trim())), sb=new Set(b.map(x=>x.toLowerCase().trim()));
  if(sa.size!==sb.size) return false;
  for(const x of sa) if(!sb.has(x)) return false;
  return true;
}
function answerEqual(a,b){
  if(a.kind!==b.kind) return false;
  if(a.kind==='text') return a.value.trim().toLowerCase()===b.value.trim().toLowerCase();
  if(a.kind==='small') return setsEqual(a.species, b.species);
  if(a.kind==='chain') return canonicalForm(a.mol||a) === canonicalForm(b.mol||b);
  return false;
}
// A "product set" answer = an unordered collection of fragments (chain and/or small).
// Two product sets are equal if there's a 1:1 matching with answerEqual.
function productSetEqual(setA, setB){
  if(setA.length!==setB.length) return false;
  const used = new Array(setB.length).fill(false);
  for(const a of setA){
    let found=-1;
    for(let i=0;i<setB.length;i++){ if(!used[i] && answerEqual(a,setB[i])){ found=i; break; } }
    if(found===-1) return false;
    used[found]=true;
  }
  return true;
}

module.exports = {
  newNode, newRingNode, newAmineNode, newEtherOxygenNode, cloneMol, findNode, neighborsOf, bondSlotsUsed, subSlotsUsed, nodeValence, usedSlots, freeSlots, implicitH,
  canonicalForm, answerEqual, productSetEqual, setsEqual
};
