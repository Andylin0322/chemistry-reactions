'use strict';
/* =========================================================================
   FUNCTIONAL GROUP SCANNER
   scanFunctionalGroups(mol) reads a molecule's raw node/edge graph (see
   engine.js's own header for that shape) and returns an array of
   FunctionalGroupFact objects -- one per distinct reactive SITE found (so
   a diene reports two 'alkene' facts, a diol reports two 'hydroxyl'
   facts). This is the "what is this molecule" half of the rule engine
   (see rules.js for "what happens to it"); it never mutates `mol` and
   never decides what reacts -- it only reports what's there.

   Each fact carries a `kind` discriminator plus enough classification
   detail for a rule's predicate to filter on (1°/2°/3°, aromatic vs
   aliphatic, etc.), and the node id(s) it's anchored to for diagnostics.
   Rules do NOT use a fact's node id to drive the actual transform --
   every existing operator in operators.js already re-locates its own
   target site(s) inside `mol` and (where the real chemistry allows more
   than one qualifying site, e.g. hydrohalogenation looping over every
   C=C) already handles all of them in one call. So a rule's transform is
   invoked at most once per matched rule, never once per fact -- see
   resolve.js.
========================================================================= */
const { findNode, neighborsOf, implicitH } = require('./engine');
const { doubleBondEdges, substitutionWeight, ringSubstitutionPositions, findCarbonylCarbon } = require('./operators');

// Real bonded-carbon count at a plain carbon that carries a flat
// substituent (OH, a halogen, ...) -- the same "otherCarbons" shape
// oxidizeAlcohol/nucleophilicSubstitutionFlat's callers already use to
// tell 1°/2°/3° apart. A ring attachment or a phenyl flag both count as
// one more real neighbour, matching how substitutionWeight treats them.
function carbonClass(mol, node){
  const n = neighborsOf(mol, node.id).length + (node.phenyl ? 1 : 0);
  return n<=1 ? '1°' : n===2 ? '2°' : '3°';
}

const META_DIRECTOR_FLAT = new Set(['NO2']);
const ACTIVATING_FLAT = new Set(['OH','NH2']);

function scanRing(mol, ring, facts){
  const occupiedPositions = [];
  for(let i=0;i<6;i++){
    const flat = ring.subs[i];
    const edge = mol.edges.find(e=>(e.a===ring.id||e.b===ring.id) && e.ringPos===i);
    if(flat || edge) occupiedPositions.push(i);
  }
  const freePositions = ringSubstitutionPositions(mol); // null if >2 substituents (this app's own directing-conflict boundary)
  facts.push({
    kind:'arene', ringId:ring.id,
    substitutionLevel: occupiedPositions.length,
    freePositions: freePositions || [],
    hasActivatingSub: ring.subs.some(s=>s && ACTIVATING_FLAT.has(s)),
  });
  // Ring-attached flat substituents are themselves functional groups (a
  // phenol's -OH, an aniline's -NH2, a haloarene's halogen, nitrobenzene's
  // -NO2) -- reported here as their own facts, `aromatic:true`, so the
  // SAME hydroxyl/amine/haloC rules used for the aliphatic case can match
  // them too (their reagentRequirement/transform pair can special-case
  // `aromatic` where the real chemistry actually differs, e.g. phenol
  // esterification refusing outright).
  ring.subs.forEach((s, pos)=>{
    if(!s) return;
    if(s==='OH') facts.push({ kind:'hydroxyl', nodeId:ring.id, ringPos:pos, class:null, aromatic:true });
    else if(s==='NH2') facts.push({ kind:'amine', nodeId:ring.id, ringPos:pos, class:'1°', aromatic:true });
    else if(s==='NO2') facts.push({ kind:'nitro', nodeId:ring.id, ringPos:pos, aromatic:true });
    else if(['F','Cl','Br','I'].includes(s)) facts.push({ kind:'haloC', nodeId:ring.id, ringPos:pos, halogen:s, class:null, aromatic:true });
  });
  // A real chain hanging off the ring via a ringPos edge, with at least
  // one H on the ring-attached (benzylic) carbon -- side-chain oxidation/
  // halogenation eligibility. tert-butylbenzene's fully-substituted
  // benzylic carbon (implicitH===0) correctly reports no benzylicCH fact,
  // same refusal sideChainOxidationToBenzoicAcid's own implicitH check
  // already relies on.
  mol.edges.forEach(e=>{
    if(e.ringPos==null) return;
    if(e.a!==ring.id && e.b!==ring.id) return;
    const outsideId = e.a===ring.id ? e.b : e.a;
    const outside = findNode(mol, outsideId);
    if(outside && !outside.ring && !outside.group && !outside.oxo && implicitH(mol, outsideId)>0){
      facts.push({ kind:'benzylicCH', nodeId:outsideId, ringId:ring.id });
    }
  });
}

function scanFunctionalGroups(mol){
  const facts = [];
  if(mol.nodes.length>0) facts.push({ kind:'organic' }); // universal -- any carbon skeleton at all; combustion keys off this rather than a specific group

  doubleBondEdges(mol).forEach(e=>{
    const a = findNode(mol, e.a), b = findNode(mol, e.b);
    const wa = substitutionWeight(mol, e.a, e), wb = substitutionWeight(mol, e.b, e);
    const total = wa + wb;
    facts.push({
      kind:'alkene', edge:{a:e.a, b:e.b},
      substitution: total===0?'mono':total===1?'di':total===2?'tri/gem-di':'tetra/tri',
      tied: wa===wb,
      phenylConjugated: a.phenyl || b.phenyl,
    });
  });

  mol.nodes.forEach(node=>{
    if(node.ring){ scanRing(mol, node, facts); return; }

    // A `phenyl:true` spectator flag (an unsubstituted benzene ring
    // attached at exactly one point, see engine.js's own header) is the
    // OTHER shape "attached to a ring" can take alongside a real
    // ring:true node -- generateAlkylbenzene and every side-chain
    // operator (sideChainOxidationToBenzoicAcid, ...) use exactly this
    // shape, not a ringPos edge, so it needs its own benzylicCH check
    // here rather than being folded into scanRing above.
    if(node.phenyl && !node.oxo && !node.group && implicitH(mol, node.id)>0){
      facts.push({ kind:'benzylicCH', nodeId:node.id, ringId:null });
    }

    // -- terminal `group` nodes (COOH/CHO/CN/COONa/N2Cl) --
    if(node.group==='COOH'){
      const nbs = neighborsOf(mol, node.id);
      const isMethanoic = nbs.length===0;
      let isEthanedioic = false;
      if(nbs.length===1){
        const other = findNode(mol, nbs[0].to);
        isEthanedioic = other.group==='COOH';
      }
      facts.push({ kind:'carboxyl', group:'COOH', nodeId:node.id, specialOxidisable: isMethanoic || isEthanedioic });
      return;
    }
    if(node.group==='COONa'){ facts.push({ kind:'carboxyl', group:'COONa', nodeId:node.id }); return; }
    if(node.group==='CN'){ facts.push({ kind:'nitrile', nodeId:node.id }); return; }
    if(node.group==='N2Cl'){ facts.push({ kind:'diazonium', nodeId:node.id }); return; }
    if(node.group==='CHO'){
      const ringEdge = mol.edges.find(e=>(e.a===node.id||e.b===node.id) && e.ringPos!=null);
      facts.push({ kind:'carbonyl', class:'aldehyde', nodeId:node.id, aromatic: !!ringEdge, hasAlphaMethyl:false });
      return;
    }

    // -- oxo (C=O) plain carbons: ketone / acyl chloride / ester / amide --
    if(node.oxo){
      if(node.subs.includes('Cl')){ facts.push({ kind:'acylHalide', nodeId:node.id }); return; }
      if(node.subs.includes('NH2')){ facts.push({ kind:'amide', class:'primary', nodeId:node.id }); return; }
      const etherO = neighborsOf(mol, node.id).find(nb=>findNode(mol, nb.to).element==='O');
      if(etherO){
        const oNode = findNode(mol, etherO.to);
        const oOtherNbs = neighborsOf(mol, oNode.id).filter(nb=>nb.to!==node.id);
        const other = oOtherNbs[0] ? findNode(mol, oOtherNbs[0].to) : null;
        if(other && other.oxo){ facts.push({ kind:'anhydride', nodeId:node.id }); return; }
        facts.push({ kind:'ester', nodeId:node.id }); return;
      }
      const amideN = neighborsOf(mol, node.id).find(nb=>findNode(mol, nb.to).element==='N');
      if(amideN){
        const nNode = findNode(mol, amideN.to);
        const nDegree = neighborsOf(mol, nNode.id).length;
        facts.push({ kind:'amide', class: nDegree>=3 ? 'tertiary' : 'secondary', nodeId:node.id, amideNodeId:nNode.id });
        return;
      }
      // Plain carbonyl carbon, no leaving group / heteroatom neighbour ->
      // ketone (2 real carbon neighbours) or a mid-build aldehyde-shaped
      // node (1 neighbour, oxo, no group -- e.g. right after
      // oxidizeAlcohol's 'CHO' branch is skipped for some hand-built
      // case); either way this app's own iodoformCleavage already treats
      // "oxo && !group" as CHO-equivalent, so mirror that here.
      const realNbs = neighborsOf(mol, node.id).length;
      if(realNbs>=2){
        const hasAlphaMethyl = neighborsOf(mol, node.id).some(nb=>{
          const n = findNode(mol, nb.to);
          return !n.ring && !n.oxo && !n.group && n.subs.length===0 && neighborsOf(mol, n.id).length===1;
        });
        facts.push({ kind:'carbonyl', class:'ketone', nodeId:node.id, aromatic:false, hasAlphaMethyl });
      } else {
        facts.push({ kind:'carbonyl', class:'aldehyde', nodeId:node.id, aromatic:false, hasAlphaMethyl:false });
      }
      return;
    }

    // -- real branching nitrogen (secondary/tertiary amine, or a charged
    //    ammonium/amine-salt centre -- charge>0 ones are already
    //    protonated and reported inert:true so a rule can refuse to
    //    re-protonate an already-protonated amine). --
    if(node.element==='N'){
      const degree = neighborsOf(mol, node.id).length;
      if(node.charge>0){ facts.push({ kind:'amineSalt', nodeId:node.id, degree }); return; }
      facts.push({ kind:'amine', nodeId:node.id, class: degree>=3?'3°':'2°', aromatic:false });
      return;
    }
    if(node.element==='O'){ return; } // ether/ester bridge oxygen -- not a site of its own, covered by the ester/anhydride facts above

    // -- flat substituents on an otherwise-plain carbon --
    if(node.subs.includes('OH')) facts.push({ kind:'hydroxyl', nodeId:node.id, class:carbonClass(mol,node), aromatic:false });
    if(node.subs.includes('NH2')) facts.push({ kind:'amine', nodeId:node.id, class:'1°', aromatic:false });
    if(node.subs.includes('NO2')) facts.push({ kind:'nitro', nodeId:node.id, aromatic:false });
    ['F','Cl','Br','I'].forEach(hal=>{
      if(node.subs.includes(hal)) facts.push({ kind:'haloC', nodeId:node.id, halogen:hal, class:carbonClass(mol,node), aromatic:false });
    });
  });

  // A composite fact -- present only when the SAME molecule carries both
  // a free -COOH and a free (non-aromatic, primary) -NH2, i.e. is itself
  // amino-acid-shaped. Peptide bond FORMATION (see rules.js's
  // two-reactant section) keys off this rather than the plain
  // 'carboxyl'/'amine' facts individually, since real syllabus scope is
  // "two alpha-amino acids condense", not "any carboxylic acid condenses
  // with any primary amine" -- formPeptideBond the OPERATOR doesn't
  // itself enforce that (it would happily amide-couple a plain acid with
  // a plain amine too), so the narrower scope has to live here instead.
  if(facts.some(f=>f.kind==='carboxyl'&&f.group==='COOH') && facts.some(f=>f.kind==='amine'&&!f.aromatic&&f.class==='1°')){
    facts.push({ kind:'aminoAcid' });
  }

  // Alkane C-H -- always structurally present wherever a plain carbon has
  // implicit H left; radicalSubstitution itself enumerates every valid
  // position (and app.html's variants UI lets the student pick among
  // them), so this is deliberately ONE aggregate fact, not one per site.
  if(mol.nodes.some(n=>!n.ring && !n.group && implicitH(mol, n.id)>0)){
    facts.push({ kind:'alkylCH' });
  }

  // Stable per-scan identity for each site -- used by resolve.js to
  // detect two DIFFERENT rules genuinely competing for the same site
  // (as opposed to two independent sites, e.g. an alkene and an alcohol
  // in the same molecule, which should both react).
  facts.forEach((f, i)=>{ f.id = i; });
  return facts;
}

module.exports = { scanFunctionalGroups, carbonClass };
