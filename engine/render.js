'use strict';
const { bondSlotsUsed, subSlotsUsed } = require('./engine');

function subsStr(subs){
  const counts={}; subs.forEach(s=>counts[s]=(counts[s]||0)+1);
  return Object.keys(counts).sort().map(k=> counts[k]>1 ? k+counts[k] : k).join('');
}
function formatCarbon(mol, node){
  if(node.group) return node.group; // COOH / CHO closes the atom entirely
  const used = bondSlotsUsed(mol, node.id) + subSlotsUsed(node);
  const h = Math.max(0, 4-used);
  let s = 'C';
  if(h>0) s += 'H'+(h>1?h:'');
  if(node.oxo) s += '(=O)';
  s += subsStr(node.subs);
  if(node.phenyl) s += '(C6H5)';
  return s;
}

function render(mol){
  if(mol.nodes.length===0) return '(empty)';
  const nodeMap={}; mol.nodes.forEach(n=>nodeMap[n.id]=n);
  const adj={}; mol.nodes.forEach(n=>adj[n.id]=[]);
  mol.edges.forEach(e=>{ adj[e.a].push({to:e.b,type:e.type}); adj[e.b].push({to:e.a,type:e.type}); });

  if(mol.nodes.length===1) return formatCarbon(mol, mol.nodes[0]);

  // find a good starting leaf: one end of the tree's longest path (diameter),
  // so the printed "main chain" reads as the longest continuous run.
  function bfsFarthest(start){
    const dist={[start]:0}; const parent={}; const q=[start];
    while(q.length){ const cur=q.shift(); adj[cur].forEach(e=>{ if(!(e.to in dist)){ dist[e.to]=dist[cur]+1; parent[e.to]=cur; q.push(e.to);} }); }
    let far=start; Object.keys(dist).forEach(k=>{ if(dist[k]>dist[far]) far=Number(k); });
    return far;
  }
  const start = bfsFarthest(bfsFarthest(mol.nodes[0].id));

  // subtree size (used to decide which branch continues the "main line" unparenthesized)
  function subtreeSize(id, cameFrom){
    let size=1;
    adj[id].forEach(e=>{ if(e.to!==cameFrom) size += subtreeSize(e.to, id); });
    return size;
  }

  function renderFrom(id, cameFrom){
    let s = formatCarbon(mol, nodeMap[id]);
    const children = adj[id]
      .filter(e=>e.to!==cameFrom)
      .map(e=>({...e, size: subtreeSize(e.to, id)}))
      .sort((a,b)=>a.size-b.size); // smallest first, largest (main continuation) last
    children.forEach((e,idx)=>{
      const sub = renderFrom(e.to, id);
      const bond = e.type==='D' ? '=' : '-';
      if(idx < children.length-1) s += '(' + bond + sub + ')';
      else s += bond + sub;
    });
    return s;
  }

  return renderFrom(start, null);
}
module.exports = { render };
