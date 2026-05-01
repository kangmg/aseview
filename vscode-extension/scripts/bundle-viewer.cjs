// scripts/bundle-viewer.cjs
// Builds a self-contained molecular_viewer_bundle.html by inlining all CDN dependencies.
// Run via: node scripts/bundle-viewer.cjs

const fs = require('fs');
const path = require('path');

const mediaDir = path.join(__dirname, '..', 'media');

function readMedia(name) {
  return fs.readFileSync(path.join(mediaDir, name), 'utf8');
}

// Read base template (already downloaded)
let html = readMedia('molecular_viewer.html');

// Helper: replace a <script src="URL"> tag with inline <script>code</script>
function inlineScript(html, urlPattern, code) {
  // Match <script src="URL"></script> (with possible attributes before src)
  const escaped = urlPattern.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
  const re = new RegExp(`<script[^>]*src="${escaped}"[^>]*>\\s*</script>`, 'i');
  const tag = `<script>\n${code}\n</script>`;
  const result = html.replace(re, tag);
  if (result === html) {
    console.warn(`  WARNING: pattern not found: ${urlPattern}`);
  }
  return result;
}

console.log('Bundling molecular_viewer.html...');

// 1. Inline three.js
html = inlineScript(html,
  'https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js',
  readMedia('three.min.js')
);
console.log('  + three.min.js');

// 2. Inline OrbitControls
html = inlineScript(html,
  'https://unpkg.com/three@0.128.0/examples/js/controls/OrbitControls.js',
  readMedia('OrbitControls.js')
);
console.log('  + OrbitControls.js');

// 3. Inline TrackballControls
html = inlineScript(html,
  'https://unpkg.com/three@0.128.0/examples/js/controls/TrackballControls.js',
  readMedia('TrackballControls.js')
);
console.log('  + TrackballControls.js');

// 4. Inline styles.js (handle both CDN URL variants)
const stylesCode = readMedia('styles.js');
const stylesUrls = [
  'https://cdn.jsdelivr.net/gh/kangmg/aseview@main/aseview/static/js/styles.js',
  'https://raw.githack.com/kangmg/aseview/main/aseview/static/js/styles.js',
];
let stylesInlined = false;
for (const url of stylesUrls) {
  const next = inlineScript(html, url, stylesCode);
  if (next !== html) { html = next; stylesInlined = true; break; }
}
if (!stylesInlined) console.warn('  WARNING: styles.js CDN tag not found');
console.log('  + styles.js');

// 5. Inline gifshot (handle direct CDN tag or __GIFSHOT_SRC__ placeholder)
const gifshotCode = readMedia('gifshot.min.js');
const gifshotCdn = 'https://cdnjs.cloudflare.com/ajax/libs/gifshot/0.3.2/gifshot.min.js';
const gifshotPlaceholder = '<script src="__GIFSHOT_SRC__"></script>';
if (html.includes(gifshotPlaceholder)) {
  html = html.replace(gifshotPlaceholder, `<script>\n${gifshotCode}\n</script>`);
  console.log('  + gifshot.min.js (placeholder)');
} else {
  html = inlineScript(html, gifshotCdn, gifshotCode);
  console.log('  + gifshot.min.js');
}

// 6. Also inject download relay (postMessage to parent for GIF/PNG saving)
const downloadRelay = `
(function() {
  document.addEventListener('click', function(e) {
    var a = e.target && e.target.closest ? e.target.closest('a[download]') : null;
    if (!a) return;
    var href = a.getAttribute('href'); if (!href) return;
    e.preventDefault(); e.stopPropagation();
    var fname = a.getAttribute('download') || 'download';
    function relay(dataUrl) { window.parent.postMessage({type:'saveFile', dataUrl:dataUrl, filename:fname}, '*'); }
    if (href.startsWith('data:')) { relay(href); }
    else if (href.startsWith('blob:')) {
      fetch(href).then(function(r){return r.blob();}).then(function(b){
        var fr = new FileReader(); fr.onload = function(){ relay(fr.result); }; fr.readAsDataURL(b);
      });
    }
  }, true);
})();
`.trim();

html = html.replace('</body>', `<script>\n${downloadRelay}\n</script>\n</body>`);
console.log('  + download relay script');

const outPath = path.join(mediaDir, 'molecular_viewer_bundle.html');
fs.writeFileSync(outPath, html, 'utf8');
const size = fs.statSync(outPath).size;
console.log(`Done! molecular_viewer_bundle.html: ${(size / 1024).toFixed(1)} KB`);
