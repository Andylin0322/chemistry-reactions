'use strict';
const { generateMolecule } = require('./generator');
const { radicalSubstitution, additionX2, hydrogenation, mildOxidationDiol, oxidativeCleavage } = require('./operators');
const { render } = require('./render');
const { canonicalForm } = require('./engine');

function fmtProduct(p){
  if(p.kind==='small') return p.species.join(' + ');
  return render(p.mol || p);
}

const POOL = [
  {id:'radSubCl', reagentText:'Cl2, UV light', mode:'variants', X:'Cl',
   gen:(d)=>generateMolecule(d, Math.random()<0.5 ? {forceExistingHalogen:true, forceNoDoubleBonds:true} : {forceNoDoubleBonds:true}),
   run:(mol)=>radicalSubstitution(mol,'Cl')},
  {id:'radSubBr', reagentText:'Br2, UV light', mode:'variants', X:'Br',
   gen:(d)=>generateMolecule(d, Math.random()<0.5 ? {forceExistingHalogen:true, forceNoDoubleBonds:true} : {forceNoDoubleBonds:true}),
   run:(mol)=>radicalSubstitution(mol,'Br')},
  {id:'addBr2', reagentText:'Br2, CCl4, room temperature, in the dark', mode:'single',
   gen:(d)=>generateMolecule(d, Math.random()<0.75 ? {minDoubleBonds:1, forcePhenyl: Math.random()<0.4?true:undefined} : {forceNoDoubleBonds:true}),
   run:(mol)=>additionX2(mol,'Br')},
  {id:'addCl2', reagentText:'Cl2, CCl4, room temperature, in the dark', mode:'single',
   gen:(d)=>generateMolecule(d, {minDoubleBonds:1, forcePhenyl: Math.random()<0.4?true:undefined}),
   run:(mol)=>additionX2(mol,'Cl')},
  {id:'hydrog', reagentText:'H2, Ni catalyst', mode:'single',
   gen:(d)=>generateMolecule(d, {minDoubleBonds:1, forcePhenyl: Math.random()<0.4?true:undefined}),
   run:(mol)=>hydrogenation(mol)},
  {id:'diol', reagentText:'KMnO4, dilute NaOH, cold', mode:'single',
   gen:(d)=>generateMolecule(d, {minDoubleBonds:1, forcePhenyl: Math.random()<0.3?true:undefined}),
   run:(mol)=>mildOxidationDiol(mol)},
  {id:'cleave', reagentText:'KMnO4, dilute H2SO4, heat', mode:'fragments',
   gen:(d)=>generateMolecule(d, {minDoubleBonds: Math.random()<0.5?2:1, forcePhenyl: Math.random()<0.3?true:undefined}),
   run:(mol)=>oxidativeCleavage(mol)},
];

function generateBatch(n){
  const problems=[]; const seenKeys=new Set();
  let guard=0;
  while(problems.length<n && guard<500){
    guard++;
    const spec = POOL[Math.floor(Math.random()*POOL.length)];
    const difficulty = 2 + Math.floor(Math.random()*4); // 2..5
    let mol;
    try{ mol = spec.gen(difficulty); } catch(e){ continue; }

    let variantsIfBig = null;
    if(spec.mode==='variants'){
      const res = spec.run(mol);
      if(!res.occurs) continue;
      if(res.variants.length>4) continue; // keep "which carbon" answers reasonable, skip oversized ones for this batch
      variantsIfBig = res.variants;
    }

    const key = spec.id+'::'+canonicalForm(mol);
    if(seenKeys.has(key)) continue;
    seenKeys.add(key);

    const result = spec.mode==='variants' ? {occurs:true, variants:variantsIfBig} : spec.run(mol);
    problems.push({ spec, mol, difficulty, result });
  }
  return problems;
}

const problems = generateBatch(15);

let md = '# Generated Problem Set\n\n';
md += problems.length+' freshly generated questions (different every run). Predict the product for each, then check the answer key.\n\n---\n\n';
problems.forEach((p,i)=>{
  const askAll = p.spec.mode==='variants' ? Math.random()<0.5 : null;
  md += '**Q'+(i+1)+'.** ('+ 'difficulty '+p.difficulty+') '+render(p.mol)+'\n\n';
  md += 'Reagents/conditions: '+p.spec.reagentText+'\n\n';
  if(askAll===true) md += '*Name all constitutionally distinct products.*\n\n';
  else if(askAll===false) md += '*Name one valid product.*\n\n';
  md += '\n';
});

md += '---\n\n# Answer Key\n\n';
problems.forEach((p,i)=>{
  md += '**A'+(i+1)+'.** ';
  if(p.spec.mode==='single'){
    md += p.result.occurs ? fmtProduct(p.result.product) : 'No reaction.';
  } else if(p.spec.mode==='fragments'){
    md += p.result.occurs ? p.result.products.map(fmtProduct).join('  +  ') : 'No reaction.';
  } else if(p.spec.mode==='variants'){
    md += p.result.variants.map(v=>render(v)).join('   OR   ');
  }
  md += '\n\n';
});

require('fs').writeFileSync('/home/claude/engine/generated_problem_set.md', md);
console.log(md);
