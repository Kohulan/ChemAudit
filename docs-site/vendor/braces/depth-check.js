// ChemAudit patched copy of braces 3.0.3 (GHSA-vfj7-8cjw-p6xm / CVE-2026-93687).
// Upstream has no fixed release; lib/parse.js now rejects nesting deeper than MAX_DEPTH.
// Drop this folder and the "braces" override in package.json once upstream ships a fix.
// Run: node vendor/braces/depth-check.js
'use strict';
const assert = require('assert');
const braces = require('.');

const deep = '{'.repeat(4900) + 'a,b' + '}'.repeat(4900);
assert.throws(() => braces(deep), /max depth/);
assert.throws(() => braces(deep, { expand: true }), /max depth/);
assert.throws(() => braces('('.repeat(5000) + 'x'), /max depth/);
assert.deepStrictEqual(braces('a/{b,c}.js', { expand: true }), ['a/b.js', 'a/c.js']);
console.log('braces depth guard ok');
