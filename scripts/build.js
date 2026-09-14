#!/usr/bin/env node
'use strict';
/* =========================================================================
   Assembles engine/*.js (the real, `node engine/test.js`-tested source)
   into app.html's embedded copies, so those two are never hand-mirrored
   again. This is NOT a real bundler -- app.html already has a tiny
   require()/module shim (see the comment above the first <script> in its
   "Engine" section) that lets each unmodified Node-style CommonJS file
   load in the browser as its own IIFE. All this script does is: wrap each
   file the same way that shim expects, and splice it between that file's
   `<!-- BUILD:name --> ... <!-- /BUILD:name -->` marker pair in app.html.

   Run after ANY change to a file in FILES below:
     node scripts/build.js
   Run in CI / before a commit to catch a forgotten build:
     node scripts/build.js --check   (exits 1 if app.html is out of sync)

   To add a new engine module: add its filename (in registry-name order --
   later files may `require()` earlier ones) to FILES below, and add a
   matching `<!-- BUILD:yourfile.js --><script>(function(){\n'use strict';\nvar module = { exports: {} };\n\n})();</script><!-- /BUILD:yourfile.js -->`
   marker pair in app.html's "Engine" section (see the existing four for
   the exact shape) before running this script.
========================================================================= */
const fs = require('fs');
const path = require('path');

const ROOT = path.join(__dirname, '..');
const APP_HTML = path.join(ROOT, 'app.html');
const ENGINE_DIR = path.join(ROOT, 'engine');

// Order matters: a later file's require('./x') must resolve to an earlier
// one's registry entry, so this must be a valid topological order.
const FILES = [
  'engine.js', 'operators.js', 'generator.js', 'render.js',
  // The rule-based reaction engine (see workflow.md's "rule-based reaction
  // engine" section): facts.js needs engine+operators, reagents.js is
  // standalone, resolve.js needs facts+reagents+engine, rules.js needs
  // operators+engine+facts.
  'facts.js', 'reagents.js', 'resolve.js', 'rules.js',
];

function registryName(filename){ return filename.replace(/\.js$/, ''); }

function wrap(filename, source){
  const name = registryName(filename);
  // Every engine/*.js file starts with a bare `'use strict';` line (as
  // required by node engine/test.js's own CommonJS semantics) -- reuse
  // that exact line rather than assuming its exact whitespace, then
  // insert the fake `module` shim right after it, matching what was
  // previously maintained by hand.
  const firstNewline = source.indexOf('\n');
  const firstLine = source.slice(0, firstNewline);
  if(firstLine.trim() !== "'use strict';"){
    throw new Error(`${filename}: expected first line to be 'use strict'; got: ${firstLine}`);
  }
  const rest = source.slice(firstNewline + 1);
  return `(function(){\n'use strict';\nvar module = { exports: {} };\n${rest}\n__registry['${name}'] = module.exports;\n})();\n`;
}

function main(){
  const check = process.argv.includes('--check');
  let html = fs.readFileSync(APP_HTML, 'utf8');
  let changedAny = false;

  for(const filename of FILES){
    const srcPath = path.join(ENGINE_DIR, filename);
    const source = fs.readFileSync(srcPath, 'utf8');
    const wrapped = wrap(filename, source);

    const startMarker = `<!-- BUILD:${filename} -->\n<script>\n`;
    const endMarker = `\n</script>\n<!-- /BUILD:${filename} -->`;
    const startIdx = html.indexOf(startMarker);
    if(startIdx === -1) throw new Error(`Marker not found in app.html: <!-- BUILD:${filename} -->`);
    const contentStart = startIdx + startMarker.length;
    const endIdx = html.indexOf(endMarker, contentStart);
    if(endIdx === -1) throw new Error(`Marker not found in app.html: <!-- /BUILD:${filename} -->`);

    const existing = html.slice(contentStart, endIdx);
    if(existing !== wrapped){
      changedAny = true;
      html = html.slice(0, contentStart) + wrapped + html.slice(endIdx);
    }
  }

  if(check){
    if(changedAny){
      console.error('app.html is OUT OF SYNC with engine/*.js -- run `node scripts/build.js`.');
      process.exit(1);
    }
    console.log('app.html is in sync with engine/*.js.');
    return;
  }

  if(changedAny){
    fs.writeFileSync(APP_HTML, html);
    console.log('app.html updated from: ' + FILES.join(', '));
  } else {
    console.log('app.html already in sync -- nothing to do.');
  }
}

main();
