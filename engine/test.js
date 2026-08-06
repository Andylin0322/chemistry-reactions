'use strict';
const { canonicalForm, answerEqual, productSetEqual, newNode, newRingNode, newAmineNode, cloneMol, implicitH, nodeValence } = require('./engine');
const { radicalSubstitution, additionX2, hydrogenation, mildOxidationDiol, oxidativeCleavage,
  ringElectrophilicSubstitution, ringFriedelCraftsAlkylation, sideChainOxidationToBenzoicAcid,
  nucleophilicSubstitutionFlat, nitrileFormation, nitrileHydrolysis, nitrileReduction, eliminationHX, alkylateAmine,
  oxidizeAlcohol, esterifyAcid, esterifyAcylChloride, ringFlatSubSwap, ringTribromination, combustion,
  diazotisation, azoCoupling, protonateAmine,
  cyanohydrinFormation, reduceCarbonyl, oxidizeAldehyde,
  carboxylicAcidSaltFormation, acidToAcylChloride, reduceCarboxylicAcid, hydrolyzeAcylChloride, acylChlorideToAmide,
  hydrolyzeEster, hydrolyzeAcidAnhydride, hydrolyzeAmide } = require('./operators');
const { generateMolecule, generateArene, generateAlkylbenzene, generateTertButylbenzene, generateHaloarene, generateAmine, generateAlkylChloride, generateSubstitutedArene,
  generateAlcohol, generatePhenol, generateCarboxylicAcid, generateAcylChloride, generateAldehyde, generateKetone, generateBenzaldehyde,
  generateEthanedioicAcid, generateAcidAnhydride, generateNitrobenzene, generatePhenylamine, generateBenzenediazonium } = require('./generator');

let pass=0, fail=0;
function t(name, cond){ if(cond){pass++; console.log('PASS -', name);} else {fail++; console.log('FAIL -', name);} }

function mkMol(n){ // helper: build a plain chain of n carbons, ids 0..n-1
  const nodes=[]; for(let i=0;i<n;i++) nodes.push({id:i,subs:[],oxo:false,group:null,phenyl:false});
  const edges=[]; for(let i=0;i<n-1;i++) edges.push({a:i,b:i+1,type:'S'});
  return {nodes,edges};
}
function setBond(mol,i,type){ mol.edges[i].type=type; return mol; }
function setSubs(mol,i,subs){ mol.nodes[i].subs=subs.slice(); return mol; }
function setGroup(mol,i,g){ mol.nodes[i].group=g; return mol; }

/* ---------- TEST 1: hexa-1,4-diene oxidative cleavage ----------
   CH2=CH-CH2-CH=CH-CH3  (C0=C1-C2-C3=C4-C5)
   expect: CO2+H2O, HOOC-CH2-COOH, CH3COOH  */
{
  const mol = mkMol(6);
  setBond(mol,0,'D'); // C0=C1
  setBond(mol,3,'D'); // C3=C4
  const res = oxidativeCleavage(mol);
  const expectedFrag1 = { kind:'small', species:['CO2','H2O'] };
  const propanedioic = (()=>{ const m=mkMol(3); setGroup(m,0,'COOH'); setGroup(m,2,'COOH'); return {kind:'chain', mol:m}; })();
  const ethanoic = (()=>{ const m=mkMol(2); setGroup(m,1,'COOH'); return {kind:'chain', mol:m}; })();
  const expected = [expectedFrag1, propanedioic, ethanoic];
  t('hexa-1,4-diene cleavage occurs', res.occurs===true);
  t('hexa-1,4-diene cleavage gives 3 fragments', res.products.length===3);
  t('hexa-1,4-diene cleavage matches hand-verified answer', productSetEqual(res.products, expected));
}

/* ---------- TEST 2: 2-methylbut-2-ene oxidative cleavage ----------
   CH3-C(CH3)=CH-CH3 : C0-C1(=C2)(branch C4=CH3)-C2=... wait build directly:
   main chain: C0-C1=C2-C3 ; branch: C4 attached to C1
   C0(CH3)-C1(=C2, branch C4 methyl)-C2(=C1, H)-C3(CH3)
   expect: propanone (CH3-C(=O)-CH3) + ethanoic acid (CH3-COOH) */
{
  const mol = mkMol(4); // C0-C1-C2-C3 chain
  setBond(mol,1,'D'); // C1=C2
  // branch methyl on C1
  mol.nodes.push({id:4,subs:[],oxo:false,group:null,phenyl:false});
  mol.edges.push({a:1,b:4,type:'S'});
  const res = oxidativeCleavage(mol);
  const propanone = (()=>{ const m={nodes:[{id:0,subs:[],oxo:false,group:null,phenyl:false},{id:1,subs:[],oxo:true,group:null,phenyl:false},{id:4,subs:[],oxo:false,group:null,phenyl:false}], edges:[{a:0,b:1,type:'S'},{a:1,b:4,type:'S'}]}; return {kind:'chain', mol:m}; })();
  const ethanoicAcid = (()=>{ const m=mkMol(2); setGroup(m,1,'COOH'); return {kind:'chain', mol:m}; })();
  t('2-methylbut-2-ene cleavage occurs', res.occurs===true);
  t('2-methylbut-2-ene cleavage gives 2 fragments (ketone + acid)', res.products.length===2);
  t('2-methylbut-2-ene cleavage matches hand-verified answer (propanone + ethanoic acid)',
    productSetEqual(res.products, [propanone, ethanoicAcid]));
}

/* ---------- TEST 3: octa-1,6-diene + Br2/CCl4 addition at both double bonds ---------- */
{
  const mol = mkMol(8); // C0..C7
  setBond(mol,0,'D'); // C0=C1
  setBond(mol,5,'D'); // C5=C6 (i.e. positions 6-7 in 1-indexed naming)
  const res = additionX2(mol, 'Br');
  const expected = mkMol(8);
  setSubs(expected,0,['Br']); setSubs(expected,1,['Br']);
  setSubs(expected,5,['Br']); setSubs(expected,6,['Br']);
  t('octa-1,6-diene dibromination at both alkenes occurs', res.occurs===true);
  t('octa-1,6-diene product matches (tetrabromooctane)', canonicalForm(res.product)===canonicalForm(expected));
}

/* ---------- TEST 4: 2-chloropropane further chlorination -> multiple valid products ---------- */
{
  const mol = mkMol(3); // C0-C1-C2, Cl on C1
  setSubs(mol,1,['Cl']);
  const res = radicalSubstitution(mol,'Cl');
  t('further chlorination occurs', res.occurs===true);
  // Expect exactly 2 distinct products: 1,2-dichloropropane and 2,2-dichloropropane
  t('further chlorination gives exactly 2 distinct products', res.variants.length===2);
  const dichloro_1_2 = mkMol(3); setSubs(dichloro_1_2,0,['Cl']); setSubs(dichloro_1_2,1,['Cl']);
  const dichloro_2_2 = mkMol(3); setSubs(dichloro_2_2,1,['Cl','Cl']);
  const canonSet = new Set(res.variants.map(canonicalForm));
  t('includes 1,2-dichloropropane', canonSet.has(canonicalForm(dichloro_1_2)));
  t('includes 2,2-dichloropropane', canonSet.has(canonicalForm(dichloro_2_2)));
}

/* ---------- TEST 5: butane monobromination -> 2 distinct products (1- and 2-bromobutane) ---------- */
{
  const mol = mkMol(4);
  const res = radicalSubstitution(mol,'Br');
  t('butane bromination occurs', res.occurs===true);
  t('butane bromination gives exactly 2 distinct products (1- and 2-bromobutane)', res.variants.length===2);
}

/* ---------- TEST 6: symmetry dedup -- propane bromination should give only 2 products not 3 ----------
   (C1 and C3 methyls are equivalent -> same product) */
{
  const mol = mkMol(3);
  const res = radicalSubstitution(mol,'Br');
  t('propane bromination dedups symmetric terminal carbons (2 distinct products)', res.variants.length===2);
}

/* ---------- TEST 7: phenyl group is never touched by any operator ---------- */
{
  const mol = mkMol(4); // C0-C1-C2-C3
  setBond(mol,2,'D'); // C2=C3
  mol.nodes[0].phenyl = true; // phenyl on C0 (sp3, far from the alkene)
  const res = additionX2(mol,'Br');
  const stillHasPhenyl = res.product.nodes.find(n=>n.id===mol.nodes[0].id).phenyl===true;
  t('phenyl substituent survives unchanged through addition operator (ring is a spectator)', stillHasPhenyl);
}

/* ---------- TEST 8: generator stress test across all difficulties, no crashes, valid valence ---------- */
{
  let ok=true;
  for(let d=1; d<=5; d++){
    for(let i=0;i<300;i++){
      const mol = generateMolecule(d, {});
      // valence check: every node's used slots <= 4
      for(const n of mol.nodes){
        const bonds = mol.edges.filter(e=>e.a===n.id||e.b===n.id).reduce((s,e)=>s+(e.type==='D'?2:1),0);
        const subs = n.subs.length + (n.oxo?2:0) + (n.phenyl?1:0);
        const total = n.group ? 4 : bonds+subs;
        if(total>4){ ok=false; console.log('VALENCE VIOLATION', d, JSON.stringify(n), bonds, subs); }
      }
      canonicalForm(mol); // should not throw
    }
  }
  t('generator produces valid-valence molecules across 1500 random trials, difficulties 1-5', ok);
}

/* ---------- TEST 9: aromatic ring D6 (rotation+reflection) symmetry ----------
   Ring positions (0-5) must compare equal/unequal purely by RELATIVE
   arrangement, never by literal array index -- whether the ring is
   free-standing, has an outside attachment, or sits buried several bonds
   into a larger tree. */
{
  function ringLeaf(sub){ const n=newNode(); n.subs=[sub]; return n; }
  function freeRing(subsAt){ // subsAt: {pos:substituentString}
    const nodes=[]; const edges=[];
    const r = newRingNode(); nodes.push(r);
    Object.keys(subsAt).forEach(posStr=>{
      const leaf = ringLeaf(subsAt[posStr]); nodes.push(leaf);
      edges.push({a:r.id, b:leaf.id, type:'S', ringPos:Number(posStr)});
    });
    return {nodes,edges};
  }

  t('monosubstituted ring: rotation does not change canonical form',
    canonicalForm(freeRing({0:'Cl'}))===canonicalForm(freeRing({3:'Cl'})));

  const ortho=canonicalForm(freeRing({0:'Cl',1:'Cl'}));
  const meta=canonicalForm(freeRing({0:'Cl',2:'Cl'}));
  const para=canonicalForm(freeRing({0:'Cl',3:'Cl'}));
  t('ortho/meta/para disubstituted ring are pairwise distinct', ortho!==meta && meta!==para && ortho!==para);
  t('ortho pattern is rotation-invariant', canonicalForm(freeRing({4:'Cl',5:'Cl'}))===ortho);
  t('meta pattern is rotation-invariant', canonicalForm(freeRing({1:'Cl',5:'Cl'}))===meta);

  t('reflection symmetry: mirrored 2-different-substituent ring is the same molecule',
    canonicalForm(freeRing({0:'Cl',1:'Br'}))===canonicalForm(freeRing({0:'Cl',5:'Br'})));
  t('non-mirror-equivalent placement of 2 different substituents is a different molecule',
    canonicalForm(freeRing({0:'Cl',1:'Br'}))!==canonicalForm(freeRing({0:'Cl',2:'Br'})));

  // 1,3,5-trisubstituted (e.g. 2,4,6-tribromophenol's substitution pattern)
  function triRing(basePos, sub){
    const nodes=[]; const edges=[];
    const r=newRingNode(); nodes.push(r);
    [0,1,2].forEach(k=>{
      const pos=(basePos+2*k)%6;
      const leaf=ringLeaf(sub); nodes.push(leaf);
      edges.push({a:r.id,b:leaf.id,type:'S',ringPos:pos});
    });
    return {nodes,edges};
  }
  t('1,3,5-trisubstituted ring is rotation-invariant regardless of starting position',
    canonicalForm(triRing(0,'Br'))===canonicalForm(triRing(1,'Br')) &&
    canonicalForm(triRing(1,'Br'))===canonicalForm(triRing(4,'Br')));

  // anchored ring (real outside attachment, e.g. toluene's ring) -- ortho/
  // meta/para relative to the attachment must still be correctly
  // distinguished, AND rotating the substituent to the mirror-image side
  // must still compare equal.
  function anchoredRing(subPos){
    const nodes=[]; const edges=[];
    const r=newRingNode(); nodes.push(r);
    const outside=newNode(); nodes.push(outside); // e.g. the CH3 of toluene
    edges.push({a:r.id,b:outside.id,type:'S',ringPos:0});
    if(subPos!=null){
      const leaf=ringLeaf('Cl'); nodes.push(leaf);
      edges.push({a:r.id,b:leaf.id,type:'S',ringPos:subPos});
    }
    return {nodes,edges};
  }
  const aOrtho=canonicalForm(anchoredRing(1)), aMeta=canonicalForm(anchoredRing(2)), aPara=canonicalForm(anchoredRing(3));
  t('anchored ring: ortho/meta/para relative to attachment are pairwise distinct',
    aOrtho!==aMeta && aMeta!==aPara && aOrtho!==aPara);
  t('anchored ring: substituent on either side of the attachment (mirror image) compares equal',
    canonicalForm(anchoredRing(5))===aOrtho && canonicalForm(anchoredRing(4))===aMeta);

  // ring buried several bonds deep in a larger tree, nowhere near the
  // canonicalForm centroid -- must still resolve correctly via whichever
  // root ends up being tried.
  function chainWithRing(pos){
    const nodes=[]; const edges=[];
    const c1=newNode(); c1.subs=['Cl']; nodes.push(c1);
    const c2=newNode(); nodes.push(c2);
    const c3=newNode(); nodes.push(c3);
    const c4=newNode(); nodes.push(c4);
    const c5=newNode(); nodes.push(c5);
    const r=newRingNode(); nodes.push(r);
    edges.push({a:c1.id,b:c2.id,type:'S'},{a:c2.id,b:c3.id,type:'S'},{a:c3.id,b:c4.id,type:'S'},{a:c4.id,b:c5.id,type:'S'},{a:c5.id,b:r.id,type:'S',ringPos:0});
    if(pos!=null){ const leaf=ringLeaf('Br'); nodes.push(leaf); edges.push({a:r.id,b:leaf.id,type:'S',ringPos:pos}); }
    return {nodes,edges};
  }
  t('ring buried deep in a larger tree: mirror-image ortho positions compare equal',
    canonicalForm(chainWithRing(1))===canonicalForm(chainWithRing(5)));
  t('ring buried deep in a larger tree: ortho vs para are distinct',
    canonicalForm(chainWithRing(1))!==canonicalForm(chainWithRing(3)));

  // flat ring.subs[pos] labels (simple substituent directly on a ring
  // position, e.g. chlorobenzene) -- parallel mechanism to the edge-based
  // one above, must obey the exact same D6 symmetry, alone or mixed with a
  // real edge substituent at another position.
  function flatSubRing(pos){ const r=newRingNode(); r.subs[pos]='Cl'; return {nodes:[r], edges:[]}; }
  t('flat ring.subs substituent: rotation-invariant', canonicalForm(flatSubRing(0))===canonicalForm(flatSubRing(4)));

  function toluenePlusFlatCl(pos){
    const r=newRingNode(); const ch3=newNode();
    const edges=[{a:r.id,b:ch3.id,type:'S',ringPos:0}];
    r.subs[pos]='Cl';
    return {nodes:[r,ch3], edges};
  }
  const fo=canonicalForm(toluenePlusFlatCl(1)), fm=canonicalForm(toluenePlusFlatCl(2)), fp=canonicalForm(toluenePlusFlatCl(3));
  t('mixed flat-sub + real edge: ortho/meta/para distinct', fo!==fm && fm!==fp && fo!==fp);
  t('mixed flat-sub + real edge: mirror side compares equal', canonicalForm(toluenePlusFlatCl(5))===fo);
}

/* ---------- TEST 10: arene operators (electrophilic substitution, Friedel-Crafts, side-chain oxidation) ---------- */
{
  const benzeneCl = ringElectrophilicSubstitution(generateArene(null), 'Cl');
  t('benzene chlorination occurs with a single product (symmetric ring)', benzeneCl.occurs===true && !!benzeneCl.product && !benzeneCl.variants);

  const benzeneEthyl = ringFriedelCraftsAlkylation(generateArene(null), 2);
  t('Friedel-Crafts alkylation of benzene occurs', benzeneEthyl.occurs===true && !!benzeneEthyl.product);
  t('a lone alkyl-chain ring product downgrades to phenyl:true (matches what a hand-built submission folds to)',
    !benzeneEthyl.product.nodes.some(n=>n.ring) && benzeneEthyl.product.nodes.some(n=>n.phenyl));

  const chloroBenzene = ringElectrophilicSubstitution(generateArene(null), 'Cl').product;
  t('a lone FLAT ring substituent (nothing to attach a phenyl flag to) stays a real ring node',
    chloroBenzene.nodes.some(n=>n.ring && n.subs.includes('Cl')));

  const tolueneCl = ringElectrophilicSubstitution(generateArene(1), 'Cl');
  t('methylbenzene ring chlorination gives exactly 2 variants (ortho+para, meta excluded)',
    tolueneCl.occurs===true && !!tolueneCl.variants && tolueneCl.variants.length===2);
  t('the ortho and para variants are structurally distinct',
    canonicalForm(tolueneCl.variants[0])!==canonicalForm(tolueneCl.variants[1]));

  // CH3 and Cl are both op-directors (same directing family), so a THIRD
  // EAS substitution on this now-disubstituted ring is chemically
  // well-posed (e.g. real p-/m-xylene-style chemistry) and is now supported
  // -- ringSubstitutionPositions was generalized this round to allow up to
  // 2 existing substituents as long as a candidate position is favoured by
  // at least one of them. Only 3+ existing substituents (a genuine
  // multi-way directing conflict this app doesn't model) still refuses.
  const disubstituted = tolueneCl.variants[0];
  const furtherSub = ringElectrophilicSubstitution(disubstituted, 'Br');
  t('a ring with 2 same-family existing substituents (CH3 + Cl) DOES undergo further EAS', furtherSub.occurs===true);

  const tripleSubstituted = cloneMol(disubstituted);
  const triRing = tripleSubstituted.nodes.find(n=>n.ring);
  const triOpenPos = [0,1,2,3,4,5].find(p=>!triRing.subs[p] && !tripleSubstituted.edges.some(e=>e.ringPos===p));
  triRing.subs[triOpenPos] = 'Br';
  const furtherSub2 = ringElectrophilicSubstitution(tripleSubstituted, 'NO2');
  t('a ring with 3 existing substituents still refuses EAS (unmodeled conflict)', furtherSub2.occurs===false);

  const benzoicAcid = { nodes:[ {id:0, subs:[], oxo:false, group:'COOH', phenyl:true, ring:false} ], edges:[] };
  const tolueneOx = sideChainOxidationToBenzoicAcid(generateAlkylbenzene(1));
  t('toluene side-chain oxidation gives benzoic acid exactly', tolueneOx.occurs===true && answerEqual({kind:'chain',mol:tolueneOx.product},{kind:'chain',mol:benzoicAcid}));
  // canonicalForm alone can't catch this one -- nodeLabel short-circuits on
  // `group` before ever looking at `phenyl`, so a real bug (the ring being
  // silently dropped by clearing phenyl:true on the new COOH node) still
  // passed the canonicalForm-only check above. Check the actual field too.
  t('the benzoic acid product actually keeps its ring attached (phenyl:true on the COOH carbon, not silently dropped)',
    tolueneOx.product.nodes[0].phenyl===true);
  const propylOx = sideChainOxidationToBenzoicAcid(generateAlkylbenzene(3));
  t('a longer alkylbenzene oxidizes down to the SAME benzoic acid regardless of chain length',
    propylOx.occurs===true && answerEqual({kind:'chain',mol:propylOx.product},{kind:'chain',mol:benzoicAcid}));

  const tbbOx = sideChainOxidationToBenzoicAcid(generateTertButylbenzene());
  t('tert-butylbenzene (no benzylic H) does not undergo side-chain oxidation', tbbOx.occurs===false);

  const tolueneSideChainCl = radicalSubstitution(generateAlkylbenzene(1), 'Cl');
  t('methylbenzene side-chain chlorination (reuses radicalSubstitution) occurs with exactly 1 distinct product',
    tolueneSideChainCl.occurs===true && tolueneSideChainCl.variants.length===1);
}

/* ---------- TEST 11: halogen derivative operators ---------- */
{
  function chain(n){
    const nodes=[]; const edges=[];
    for(let i=0;i<n;i++) nodes.push(newNode());
    for(let i=0;i<n-1;i++) edges.push({a:nodes[i].id,b:nodes[i+1].id,type:'S'});
    return {nodes,edges};
  }

  const propylChloride = chain(3); propylChloride.nodes[0].subs.push('Cl');
  const hydrolysis = nucleophilicSubstitutionFlat(propylChloride, 'Cl', 'OH');
  const propanol = chain(3); propanol.nodes[0].subs.push('OH');
  t('hydrolysis (X->OH) gives propan-1-ol exactly', hydrolysis.occurs===true && answerEqual({kind:'chain',mol:hydrolysis.product},{kind:'chain',mol:propanol}));

  const amineForm = nucleophilicSubstitutionFlat(propylChloride, 'Cl', 'NH2');
  const propanamine = chain(3); propanamine.nodes[0].subs.push('NH2');
  t('amine formation (X->NH2) gives propan-1-amine exactly', amineForm.occurs===true && answerEqual({kind:'chain',mol:amineForm.product},{kind:'chain',mol:propanamine}));

  const haloarene = generateHaloarene('Cl');
  t('haloarene hydrolysis does not occur (ring halogen never inspected)', nucleophilicSubstitutionFlat(haloarene,'Cl','OH').occurs===false);
  t('haloarene nitrile formation does not occur', nitrileFormation(haloarene,'Cl').occurs===false);
  t('haloarene elimination does not occur', eliminationHX(haloarene,'Cl').occurs===false);

  const nitrile = nitrileFormation(propylChloride, 'Cl');
  t('nitrile formation occurs and adds exactly one new (CN) node', nitrile.occurs===true && nitrile.product.nodes.length===4 && nitrile.product.nodes.some(n=>n.group==='CN'));

  const acid = nitrileHydrolysis(nitrile.product);
  t('nitrile hydrolysis gives a COOH group node, no CN left', acid.occurs===true && acid.product.nodes.some(n=>n.group==='COOH') && !acid.product.nodes.some(n=>n.group==='CN'));

  const reducedAmine = nitrileReduction(nitrile.product);
  const butanamine = chain(4); butanamine.nodes[3].subs.push('NH2');
  t('nitrile reduction gives the expected primary amine one carbon longer than the nitrile chain',
    reducedAmine.occurs===true && answerEqual({kind:'chain',mol:reducedAmine.product},{kind:'chain',mol:butanamine}));

  const symmetricHalide = chain(3); symmetricHalide.nodes[1].subs.push('Cl');
  const symElim = eliminationHX(symmetricHalide, 'Cl');
  t('elimination on a symmetric halide gives exactly 1 distinct alkene', symElim.occurs===true && symElim.variants.length===1);

  const asymmetricHalide = chain(4); asymmetricHalide.nodes[1].subs.push('Cl');
  const asymElim = eliminationHX(asymmetricHalide, 'Cl');
  t('elimination on an asymmetric halide gives exactly 2 distinct alkenes (but-1-ene + but-2-ene)', asymElim.occurs===true && asymElim.variants.length===2);

  // C(CH3)3-CH2-Cl -- the only beta-carbon (the quaternary center) has no free H
  const central = newNode(), ch2 = newNode(); ch2.subs.push('Cl');
  const m1 = newNode(), m2 = newNode(), m3 = newNode();
  const noBetaH = { nodes:[central, ch2, m1, m2, m3], edges:[
    {a:central.id, b:ch2.id, type:'S'}, {a:central.id, b:m1.id, type:'S'},
    {a:central.id, b:m2.id, type:'S'}, {a:central.id, b:m3.id, type:'S'}
  ]};
  t('elimination does not occur when the only beta-carbon has no H', eliminationHX(noBetaH, 'Cl').occurs===false);

  let sawF = false;
  for(let i=0;i<20;i++){
    const mol = generateMolecule(3, {forceExistingHalogen:true, forceHalogenType:'F', forceNoDoubleBonds:true});
    if(mol.nodes.some(n=>n.subs.includes('F'))) sawF = true;
  }
  t('forceHalogenType:"F" reliably produces a fluorine substituent', sawF);
}

/* ---------- TEST 12: branching/charged amine nitrogen + multi-substituent arene EAS ---------- */
{
  function chain(n){
    const nodes=[]; const edges=[];
    for(let i=0;i<n;i++) nodes.push(newNode());
    for(let i=0;i<n-1;i++) edges.push({a:nodes[i].id,b:nodes[i+1].id,type:'S'});
    return {nodes,edges};
  }

  const secondary = generateAmine([1,1]); // (CH3)2NH
  const secN = secondary.nodes.find(n=>n.element==='N');
  t('secondary amine N has degree 2', secondary.edges.filter(e=>e.a===secN.id||e.b===secN.id).length===2);
  t('secondary amine N has exactly 1 implicit H', implicitH(secondary, secN.id)===1);

  const tertiary = generateAmine([1,1,1]); // (CH3)3N
  const tertN = tertiary.nodes.find(n=>n.element==='N');
  t('tertiary amine N has 0 implicit H', implicitH(tertiary, tertN.id)===0);
  t('secondary and tertiary amines have distinct canonical forms', canonicalForm(secondary)!==canonicalForm(tertiary));

  const arm12 = generateAmine([1,2]), arm21 = generateAmine([2,1]);
  t('amine canonical form is order-independent across arms', canonicalForm(arm12)===canonicalForm(arm21));

  const primaryFlat = chain(2); primaryFlat.nodes[0].subs.push('NH2'); // ethylamine, flat-sub form
  t('flat-sub primary amine and a structural N node are NOT canonical-form-equal', canonicalForm(primaryFlat)!==canonicalForm(secondary));

  // primary (flat NH2) + CH3Cl -> secondary amine, via transparent promotion
  const promoted = alkylateAmine(primaryFlat, generateAlkylChloride(1), 'Cl');
  t('alkylating a flat-NH2 primary amine occurs', promoted.occurs===true);
  const promotedN = promoted.product.nodes.find(n=>n.element==='N');
  t('promoted product has a real N node of degree 2', !!promotedN && promoted.product.edges.filter(e=>e.a===promotedN.id||e.b===promotedN.id).length===2);
  t('promoted product matches the directly-generated secondary amine (ethyl+methyl arms)', canonicalForm(promoted.product)===canonicalForm(generateAmine([2,1])));

  const toTertiary = alkylateAmine(secondary, generateAlkylChloride(1), 'Cl');
  t('secondary + R-X -> tertiary amine occurs', toTertiary.occurs===true);
  const toTertN = toTertiary.product.nodes.find(n=>n.element==='N');
  t('tertiary product has degree 3 and charge 0', toTertiary.product.edges.filter(e=>e.a===toTertN.id||e.b===toTertN.id).length===3 && toTertN.charge===0);

  const toQuat = alkylateAmine(tertiary, generateAlkylChloride(1), 'Cl');
  t('tertiary + R-X -> quaternary ammonium occurs', toQuat.occurs===true);
  const quatN = toQuat.product.nodes.find(n=>n.element==='N');
  t('quaternary ammonium N has degree 4 and charge 1', toQuat.product.edges.filter(e=>e.a===quatN.id||e.b===quatN.id).length===4 && quatN.charge===1);

  const overAlkylate = alkylateAmine(toQuat.product, generateAlkylChloride(1), 'Cl');
  t('a quaternary ammonium centre refuses further alkylation', overAlkylate.occurs===false);

  // multi-substituent arene EAS: nitrobenzene (meta-director) must direct to
  // meta, not ortho/para -- this was WRONG under the old hardcoded
  // ortho/para-only ringSubstitutionPositions, so it's worth pinning down
  // permanently now that directing type is modeled.
  function monoRing(sub){ const r={nodes:[{id:0,subs:new Array(6).fill(null),oxo:false,group:null,phenyl:false,ring:true}],edges:[]}; r.nodes[0].subs[0]=sub; return r; }
  const nitroBromination = ringElectrophilicSubstitution(monoRing('NO2'), 'Br');
  const metaProduct = monoRing('NO2'); metaProduct.nodes[0].subs[2]='Br';
  const orthoProduct = monoRing('NO2'); orthoProduct.nodes[0].subs[1]='Br';
  t('nitrobenzene bromination occurs and gives a single (symmetric) meta product', nitroBromination.occurs===true && (nitroBromination.variants?nitroBromination.variants.length:1)===1);
  const nitroActual = nitroBromination.product || nitroBromination.variants[0];
  t('nitrobenzene bromination product is the META isomer', canonicalForm(nitroActual)===canonicalForm(metaProduct));
  t('nitrobenzene bromination product is NOT the ortho isomer', canonicalForm(nitroActual)!==canonicalForm(orthoProduct));

  let sweepErrors = 0;
  for(let i=0;i<200;i++){
    try { canonicalForm(generateSubstitutedArene(2+Math.floor(Math.random()*4))); }
    catch(e){ sweepErrors++; }
  }
  t('generateSubstitutedArene produces valid molecules across 200 iterations', sweepErrors===0);
}

/* ---------- TEST 13: Hydroxy Compounds (alcohols, phenols, esters) ---------- */
{
  const primary = generateAlcohol([2]);
  const carbinol = primary.nodes.find(n=>n.subs.includes('OH'));
  t('primary alcohol carbinol carbon has 1 real neighbour (2 implicit H)', implicitH(primary, carbinol.id)===2);

  const secondary = generateAlcohol([1,1]);
  const secCarbinol = secondary.nodes.find(n=>n.subs.includes('OH'));
  t('secondary alcohol carbinol carbon has 1 implicit H', implicitH(secondary, secCarbinol.id)===1);

  const tertiary = generateAlcohol([1,1,1]);
  const tertCarbinol = tertiary.nodes.find(n=>n.subs.includes('OH'));
  t('tertiary alcohol carbinol carbon has 0 implicit H', implicitH(tertiary, tertCarbinol.id)===0);

  t('ROH -> RCl via nucleophilicSubstitutionFlat', nucleophilicSubstitutionFlat(primary,'OH','Cl').occurs===true);
  t('ROH -> RONa via nucleophilicSubstitutionFlat', nucleophilicSubstitutionFlat(primary,'OH','ONa').occurs===true);
  t('dehydration (eliminationHX with OH) gives an alkene', eliminationHX(primary,'OH').variants[0].edges.some(e=>e.type==='D'));
  t('dehydration of methanol (no beta-carbon) refuses', eliminationHX(generateAlcohol([]),'OH').occurs===false);

  const phenol = generatePhenol();
  t('nucleophilicSubstitutionFlat ignores ring-attached OH (phenol)', nucleophilicSubstitutionFlat(phenol,'OH','ONa').occurs===false);
  const phenolSwap = ringFlatSubSwap(phenol,'OH','ONa');
  t('ringFlatSubSwap converts phenol -> ArONa', phenolSwap.occurs===true && phenolSwap.product.nodes.find(n=>n.ring).subs.includes('ONa'));

  t('primary alcohol oxidises to CHO', oxidizeAlcohol(primary,'CHO').occurs===true);
  t('primary alcohol oxidises to COOH', oxidizeAlcohol(primary,'COOH').occurs===true);
  const ketoneRes = oxidizeAlcohol(secondary,'ketone');
  t('secondary alcohol oxidises to a ketone (oxo, no group)', ketoneRes.occurs===true && ketoneRes.product.nodes.some(n=>n.oxo && !n.group));
  t('secondary alcohol refuses CHO/COOH target', oxidizeAlcohol(secondary,'CHO').occurs===false);
  t('tertiary alcohol oxidation refuses entirely', oxidizeAlcohol(tertiary,'COOH').occurs===false && oxidizeAlcohol(tertiary,'ketone').occurs===false);

  const acid = generateCarboxylicAcid(2);
  const ester = esterifyAcid(acid, generateAlcohol([1]));
  t('carboxylic acid + alcohol esterifies', ester.occurs===true);
  const esterO = ester.product.nodes.find(n=>n.element==='O');
  t('ester product has one element:O node with 2 real neighbours', !!esterO && ester.product.edges.filter(e=>e.a===esterO.id||e.b===esterO.id).length===2);
  t('carboxylic acid + phenol does NOT esterify', esterifyAcid(acid, generatePhenol()).occurs===false);

  const acyl = generateAcylChloride(2);
  t('acyl chloride + alcohol acylates', esterifyAcylChloride(acyl, generateAlcohol([1])).occurs===true);
  t('acyl chloride + phenol DOES acylate (aryl ester)', esterifyAcylChloride(acyl, generatePhenol()).occurs===true);

  const nitration = ringElectrophilicSubstitution(generatePhenol(), 'NO2');
  t('phenol nitration gives 2 variants (2-/4-nitrophenol)', nitration.occurs===true && nitration.variants.length===2);
  const tribromination = ringTribromination(generatePhenol(), 'OH');
  const triRing = tribromination.product.nodes.find(n=>n.ring);
  t('phenol tribromination gives exactly 3 Br + 1 OH, single fixed product', tribromination.occurs===true && triRing.subs.filter(s=>s==='Br').length===3 && triRing.subs.includes('OH'));

  const combustionRes = combustion(generateAlcohol([2]));
  t('combustion always gives CO2+H2O-only products', combustionRes.occurs===true && combustionRes.products.length===1 &&
    JSON.stringify(combustionRes.products[0].species.slice().sort())===JSON.stringify(['CO2','H2O']));
}

/* ---------- TEST 14: Carbonyl Compounds ---------- */
{
  const ald = generateAldehyde(3);
  const ket = generateKetone(1,1);
  const benz = generateBenzaldehyde();

  const cyanoAld = cyanohydrinFormation(ald);
  t('cyanohydrin formation from aldehyde occurs', cyanoAld.occurs===true);
  const cyanoAldC = cyanoAld.product.nodes.find(n=>n.subs.includes('OH'));
  t('aldehyde cyanohydrin carbon keeps 1 implicit H, has a CN group node', implicitH(cyanoAld.product, cyanoAldC.id)===1 && cyanoAld.product.nodes.some(n=>n.group==='CN'));

  const cyanoKet = cyanohydrinFormation(ket);
  t('cyanohydrin formation from ketone occurs', cyanoKet.occurs===true);
  const cyanoKetC = cyanoKet.product.nodes.find(n=>n.subs.includes('OH'));
  t('ketone cyanohydrin carbon is fully substituted (0 implicit H)', implicitH(cyanoKet.product, cyanoKetC.id)===0);

  const redAld = reduceCarbonyl(ald);
  t('aldehyde reduction to primary alcohol occurs', redAld.occurs===true);
  const redAldC = redAld.product.nodes.find(n=>n.subs.includes('OH'));
  t('primary alcohol product has 2 implicit H', implicitH(redAld.product, redAldC.id)===2);

  const redKet = reduceCarbonyl(ket);
  t('ketone reduction to secondary alcohol occurs', redKet.occurs===true);
  const redKetC = redKet.product.nodes.find(n=>n.subs.includes('OH'));
  t('secondary alcohol product has 1 implicit H', implicitH(redKet.product, redKetC.id)===1);

  t('aldehyde oxidises to carboxylic acid', oxidizeAldehyde(ald).occurs===true && oxidizeAldehyde(ald).product.nodes.some(n=>n.group==='COOH'));
  t('ketone oxidation refuses (Tollens/Fehlings/K2Cr2O7 all negative)', oxidizeAldehyde(ket).occurs===false);
  t('benzaldehyde still oxidises (Tollens is positive for it)', oxidizeAldehyde(benz).occurs===true);
}

/* ---------- TEST 15: Carboxylic Acids & Derivatives ---------- */
{
  const acid = generateCarboxylicAcid(3);
  t('COOH -> COONa salt formation occurs', carboxylicAcidSaltFormation(acid).occurs===true);

  const acyl = acidToAcylChloride(acid);
  t('COOH -> acyl chloride occurs', acyl.occurs===true);
  const backToAcid = hydrolyzeAcylChloride(acyl.product);
  t('acyl chloride hydrolysis round-trips to the original acid', backToAcid.occurs===true && canonicalForm(backToAcid.product)===canonicalForm(acid));

  const reduced = reduceCarboxylicAcid(acid);
  const reducedC = reduced.product.nodes.find(n=>n.subs.includes('OH'));
  t('LiAlH4 reduces acid to primary alcohol (2 implicit H)', reduced.occurs===true && implicitH(reduced.product, reducedC.id)===2);

  t('methanoic acid is a lone COOH node', generateCarboxylicAcid(1).nodes.length===1);
  const oxalic = generateEthanedioicAcid();
  t('ethanedioic acid is two bonded COOH nodes', oxalic.nodes.length===2 && oxalic.nodes.every(n=>n.group==='COOH'));

  const primaryAmide = acylChlorideToAmide(acyl.product, null);
  t('acyl chloride + NH3 -> primary amide occurs', primaryAmide.occurs===true && primaryAmide.product.nodes.some(n=>n.oxo && n.subs.includes('NH2')));

  const primAmineMol = nucleophilicSubstitutionFlat(generateAlkylChloride(2), 'Cl', 'NH2').product;
  const secondaryAmide = acylChlorideToAmide(acyl.product, primAmineMol);
  const secN = secondaryAmide.product.nodes.find(n=>n.element==='N');
  t('acyl chloride + primary amine -> secondary amide (N degree 2, 1 implicit H)', secondaryAmide.occurs===true &&
    secondaryAmide.product.edges.filter(e=>e.a===secN.id||e.b===secN.id).length===2 && implicitH(secondaryAmide.product, secN.id)===1);

  const tertiaryAmide = acylChlorideToAmide(acyl.product, generateAmine([1,1]));
  const tertN = tertiaryAmide.product.nodes.find(n=>n.element==='N');
  t('acyl chloride + secondary amine -> tertiary amide (N degree 3, 0 implicit H)', tertiaryAmide.occurs===true &&
    tertiaryAmide.product.edges.filter(e=>e.a===tertN.id||e.b===tertN.id).length===3 && implicitH(tertiaryAmide.product, tertN.id)===0);

  const alcohol = generateAlcohol([2]);
  const ester = esterifyAcid(acid, alcohol).product;
  const acidicHydrolysis = hydrolyzeEster(ester, 'COOH');
  t('acidic ester hydrolysis gives 2 fragments matching the original acid+alcohol', acidicHydrolysis.occurs===true &&
    productSetEqual(acidicHydrolysis.products, [{kind:'chain',mol:acid},{kind:'chain',mol:alcohol}]));
  const alkalineHydrolysis = hydrolyzeEster(ester, 'COONa');
  t('alkaline ester hydrolysis gives a COONa fragment', alkalineHydrolysis.occurs===true && alkalineHydrolysis.products.some(p=>p.mol.nodes.some(n=>n.group==='COONa')));

  const anhydride = generateAcidAnhydride(2,3);
  const anhydrideHydrolysis = hydrolyzeAcidAnhydride(anhydride);
  t('acid anhydride hydrolysis gives 2 COOH fragments', anhydrideHydrolysis.occurs===true && anhydrideHydrolysis.products.length===2 &&
    anhydrideHydrolysis.products.every(p=>p.mol.nodes.some(n=>n.group==='COOH')));

  const amide = acylChlorideToAmide(acyl.product, null).product;
  t('acidic amide hydrolysis gives COOH', hydrolyzeAmide(amide,'COOH').occurs===true && hydrolyzeAmide(amide,'COOH').product.nodes.some(n=>n.group==='COOH'));
  t('alkaline amide hydrolysis gives COONa', hydrolyzeAmide(amide,'COONa').occurs===true && hydrolyzeAmide(amide,'COONa').product.nodes.some(n=>n.group==='COONa'));
}

/* ---------- TEST 16: generic charge + Nitrogen Compounds ---------- */
{
  function chain(n){
    const nodes=[]; const edges=[];
    for(let i=0;i<n;i++) nodes.push(newNode());
    for(let i=0;i<n-1;i++) edges.push({a:nodes[i].id,b:nodes[i+1].id,type:'S'});
    return {nodes,edges};
  }

  t('neutral carbon valence 4', nodeValence({charge:0})===4);
  t('carbocation (C+) valence 3', nodeValence({charge:1})===3);
  t('carbanion (C-) valence 3', nodeValence({charge:-1})===3);
  t('N+ valence 4 (ammonium, backward-compat)', nodeValence({element:'N',charge:1})===4);
  t('N- valence 2 (amide ion)', nodeValence({element:'N',charge:-1})===2);
  t('O+ valence 3 (oxonium)', nodeValence({element:'O',charge:1})===3);
  t('O- valence 1 (alkoxide)', nodeValence({element:'O',charge:-1})===1);

  const c3 = chain(3);
  c3.nodes[1].charge = 1;
  t('a plain-carbon carbocation carries 1 implicit H (isopropyl-cation shape)', implicitH(c3, c3.nodes[1].id)===1);

  const phenylamine = generatePhenylamine();
  const diazo = diazotisation(phenylamine);
  t('diazotisation occurs, gives a ring-attached N2Cl group', diazo.occurs===true && diazo.product.nodes.some(n=>n.group==='N2Cl'));
  t('diazotisation refuses phenol (no ring -NH2)', diazotisation(generatePhenol()).occurs===false);

  const primaryAmine = chain(2); primaryAmine.nodes[0].subs.push('NH2');
  const salt = protonateAmine(primaryAmine);
  const saltN = salt.product.nodes.find(n=>n.element==='N');
  t('protonateAmine gives RNH3+ (charge 1, degree 1, 3 implicit H)', salt.occurs===true && saltN.charge===1 &&
    salt.product.edges.filter(e=>e.a===saltN.id||e.b===saltN.id).length===1 && implicitH(salt.product, saltN.id)===3);
  const quatBase = generateAmine([1,1,1]);
  const quat = alkylateAmine(quatBase, generateAlkylChloride(1), 'Cl').product;
  t('protonateAmine refuses an already-quaternary N', protonateAmine(quat).occurs===false);

  const phenylamineTri = ringTribromination(generatePhenylamine(), 'NH2');
  const triRing = phenylamineTri.product.nodes.find(n=>n.ring);
  t('phenylamine tribromination gives 3 Br + 1 NH2', phenylamineTri.occurs===true && triRing.subs.filter(s=>s==='Br').length===3 && triRing.subs.includes('NH2'));

  const nitroReduced = ringFlatSubSwap(generateNitrobenzene(), 'NO2', 'NH2');
  t('nitrobenzene -> phenylamine via ringFlatSubSwap', nitroReduced.occurs===true && nitroReduced.product.nodes.find(n=>n.ring).subs.includes('NH2'));

  const diazonium = diazotisation(generatePhenylamine()).product;
  const coupled = azoCoupling(diazonium, generatePhenol());
  const coupledProducts = coupled.variants || [coupled.product];
  const coupledN = coupledProducts[0].nodes.filter(n=>n.element==='N');
  t('azo coupling with phenol occurs, gives 2 neutral N joined by a double bond, 2 rings',
    coupled.occurs===true && coupledN.length===2 && coupledN.every(n=>(n.charge||0)===0) &&
    coupledProducts[0].nodes.filter(n=>n.ring).length===2 &&
    coupledProducts[0].edges.some(e=>e.type==='D' && [e.a,e.b].includes(coupledN[0].id) && [e.a,e.b].includes(coupledN[1].id)));
  t('azo coupling refuses a non-diazonium input', azoCoupling(generateAldehyde(2), generatePhenol()).occurs===false);
}

console.log('\\n'+pass+' passed, '+fail+' failed');
process.exit(fail>0?1:0);
