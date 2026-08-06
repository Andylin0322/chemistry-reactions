'use strict';
const { generateMolecule } = require('./generator');
const { radicalSubstitution, additionX2, hydrogenation, mildOxidationDiol, oxidativeCleavage } = require('./operators');
const { render } = require('./render');

// For "which carbon" style questions: regenerate with a smaller/simpler
// molecule if the site count would exceed maxVariants (default 4), so
// "name all products" stays a reasonable ask rather than an enumeration slog.
function generateCappedSubstitution(difficulty, X, constraints, maxVariants){
  maxVariants = maxVariants || 4;
  let d = difficulty;
  for(let attempt=0; attempt<30; attempt++){
    const mol = generateMolecule(d, constraints);
    const res = radicalSubstitution(mol, X);
    if(res.occurs && res.variants.length<=maxVariants) return {mol, res};
    if(attempt>0 && attempt%5===0 && d>1) d--; // step difficulty down if repeatedly too big
  }
  // fallback: force a small molecule directly
  const mol = generateMolecule(1, {...constraints, minChainLen:2});
  return {mol, res: radicalSubstitution(mol, X)};
}

function fmtProduct(p){
  if(p.kind==='small') return p.species.join('+');
  return render(p.mol || p);
}

console.log('=== CHEMOSELECTIVITY: alkane vs alkene, same reagent Br2 ===');
for(let i=0;i<2;i++){
  const mol = generateMolecule(3, {}); // random: may or may not have a double bond
  const hasDB = mol.edges.some(e=>e.type==='D');
  console.log('\\nGenerated:', render(mol), hasDB?'(has C=C)':'(alkane only)');
  console.log('Reagents: Br2, CCl4, room temperature, dark');
  if(hasDB){
    const r = additionX2(mol,'Br');
    console.log('-> Electrophilic addition occurs:', fmtProduct(r.product));
  } else {
    console.log('-> No reaction (no C=C for addition; no UV/heat for substitution)');
  }
}

console.log('\\n=== BENZENE AS SPECTATOR: forced phenyl + forced double bond ===');
for(let i=0;i<2;i++){
  const mol = generateMolecule(3, {forcePhenyl:true, minDoubleBonds:1});
  console.log('\\nGenerated:', render(mol));
  const r = additionX2(mol,'Br');
  console.log('Reagents: Br2, CCl4, room temperature, dark');
  console.log('-> ', r.occurs? fmtProduct(r.product) : 'No reaction');
}

console.log('\\n=== MULTI-SITE TRACKING: forced 2 double bonds, oxidative cleavage ===');
for(let i=0;i<2;i++){
  const mol = generateMolecule(4, {minDoubleBonds:2});
  console.log('\\nGenerated:', render(mol));
  console.log('Reagents: KMnO4, dilute H2SO4, heat');
  const r = oxidativeCleavage(mol);
  if(r.occurs) console.log('-> Fragments:', r.products.map(fmtProduct).join('  +  '));
  else console.log('-> No reaction');
}

console.log('\\n=== BRANCHING: forced branch, mild oxidation (diol) ===');
for(let i=0;i<2;i++){
  const mol = generateMolecule(4, {minDoubleBonds:1});
  console.log('\\nGenerated:', render(mol));
  console.log('Reagents: KMnO4, dilute NaOH, cold');
  const r = mildOxidationDiol(mol);
  console.log('-> ', r.occurs? fmtProduct(r.product) : 'No reaction');
}

console.log('\\n=== WHICH CARBON (capped at <=4 valid products): further substitution ===');
for(let i=0;i<3;i++){
  const {mol, res} = generateCappedSubstitution(4, 'Cl', {forceExistingHalogen:true, forceNoDoubleBonds:true}, 4);
  console.log('\\nGenerated:', render(mol));
  console.log('Reagents: Cl2, UV light');
  console.log('-> ', res.variants.length, 'distinct valid product(s):');
  res.variants.forEach(v=>console.log('     ', render(v)));
}

