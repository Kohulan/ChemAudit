#!/usr/bin/env node
/**
 * Hardens the image-size parsers against CVE-2025-71330 (ICNS) and
 * CVE-2025-71329 (JXL/HEIF): three loops advance by a length/size field read
 * straight from the input, so a zero-valued field leaves the offset unmoved and
 * spins the event loop forever.
 *
 * image-size 2.0.2 is the newest release and upstream has shipped no fix, so the
 * guards are applied locally on postinstall. image-size reaches us as a
 * transitive dependency of @docusaurus/mdx-loader, which sizes local images at
 * docs build time.
 *
 * ponytail: string-match patching of a bundled dist (esbuild emits the same
 * source into five files). Upgrade path is deleting this script and the
 * postinstall hook as soon as upstream releases above 2.0.2.
 *
 * Usage:
 *   node scripts/patch-image-size.mjs           apply the guards (idempotent)
 *   node scripts/patch-image-size.mjs --check   verify the guards are present
 */
import { readFileSync, writeFileSync, existsSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const PKG_DIR = join(dirname(dirname(fileURLToPath(import.meta.url))), 'node_modules', 'image-size');

// esbuild bundles the same parser source into each of these.
const FILES = ['index.cjs', 'index.mjs', 'lookup.cjs', 'lookup.mjs', 'detector.cjs'];

// Each bundle renames locals independently, so match a single line plus its
// indentation rather than surrounding context. `marker` proves a file is done.
const PATCHES = [
  {
    name: 'ICNS entry length (CVE-2025-71330)',
    marker: 'imageHeader[1] > 0',
    vulnerable: /( *)imageOffset \+= imageHeader\[1\];/g,
    patched: '$1if (!(imageHeader[1] > 0)) break;\n$1imageOffset += imageHeader[1];',
  },
  {
    name: 'HEIF ispe box size (CVE-2025-71329)',
    marker: 'nextIspeOffset',
    vulnerable: /( *)currentOffset = ispeBox\.offset \+ ispeBox\.size;/g,
    patched:
      '$1const nextIspeOffset = ispeBox.offset + ispeBox.size;\n' +
      '$1if (nextIspeOffset <= currentOffset) break;\n' +
      '$1currentOffset = nextIspeOffset;',
  },
  {
    name: 'JXL jxlp box size (CVE-2025-71329)',
    marker: 'nextJxlpOffset',
    vulnerable: /( *)offset = jxlpBox\.offset \+ jxlpBox\.size;/g,
    patched:
      '$1const nextJxlpOffset = jxlpBox.offset + jxlpBox.size;\n' +
      '$1if (nextJxlpOffset <= offset) break;\n' +
      '$1offset = nextJxlpOffset;',
  },
];

const checkOnly = process.argv.includes('--check');

if (!existsSync(PKG_DIR)) {
  // Installed with --ignore-scripts, or image-size dropped from the tree.
  console.log(`patch-image-size: ${PKG_DIR} not found, nothing to do`);
  process.exit(checkOnly ? 1 : 0);
}

const failures = [];
const landed = new Set();
let changed = 0;

for (const file of FILES) {
  const path = join(PKG_DIR, 'dist', file);
  if (!existsSync(path)) {
    failures.push(`${file}: missing`);
    continue;
  }

  const original = readFileSync(path, 'utf8');
  let source = original;

  for (const patch of PATCHES) {
    if (source.includes(patch.marker)) {
      landed.add(patch.name); // already guarded
      continue;
    }
    patch.vulnerable.lastIndex = 0; // regexes are global, so test() is stateful
    if (!patch.vulnerable.test(source)) continue; // parser not in this bundle
    if (checkOnly) {
      failures.push(`${file}: ${patch.name} — guard missing`);
      continue;
    }
    source = source.replace(patch.vulnerable, patch.patched);
    landed.add(patch.name);
  }

  if (source !== original) {
    writeFileSync(path, source);
    changed++;
  }
}

// A guard that landed nowhere means the vulnerable code moved or was renamed.
for (const patch of PATCHES) {
  if (!landed.has(patch.name)) failures.push(`${patch.name} — not found in any bundle`);
}

if (failures.length > 0) {
  console.error('patch-image-size: FAILED');
  for (const failure of failures) console.error(`  - ${failure}`);
  console.error('  Verify whether image-size still needs patching, then update or delete this script.');
  process.exit(1);
}

// Regression proof for the ICNS loop: a zero-valued entry length used to spin
// until the heap died. Only reached once the guards above are confirmed present,
// so this can never be the thing that hangs.
const { imageSize } = await import(join(PKG_DIR, 'dist', 'index.mjs'));
const icns = new Uint8Array(24);
icns.set(new TextEncoder().encode('icns'), 0);
new DataView(icns.buffer).setUint32(4, icns.length); // file length
icns.set(new TextEncoder().encode('ic09'), 8); // entry type -> 512 x 512
new DataView(icns.buffer).setUint32(12, 0); // entry length 0
const { width, height } = imageSize(icns);
if (width !== 512 || height !== 512) {
  console.error(`patch-image-size: FAILED — crafted ICNS returned ${width}x${height}, expected 512x512`);
  process.exit(1);
}

console.log(
  checkOnly
    ? `patch-image-size: all ${PATCHES.length} guards present in ${FILES.length} files, ICNS regression passes`
    : `patch-image-size: guards applied (${changed} of ${FILES.length} files rewritten), ICNS regression passes`
);
