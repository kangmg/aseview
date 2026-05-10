# Styles

ASEView ships multiple molecular rendering styles. The previews below use the low-level iframe API directly: each card loads `molecular_viewer.html`, sends `setSettings`, then sends `setData`.

The molecule is acetamide, matching `ase.build.molecule("CH3CONH2")`.

<div class="aseview-style-grid" id="aseview-style-grid"></div>

<script>
(function () {
  const STYLES = [
    'default',
    'cartoon',
    'glossy',
    'metallic',
    'rowan',
    'grey',
    'bubble',
    'cinematic',
    'neon',
    '2d'
  ];

  const DARK_BACKGROUNDS = new Set(['bubble', 'cinematic', 'neon', '2d']);
  const TEMPLATE_URL = 'https://cdn.jsdelivr.net/gh/kangmg/aseview@main/aseview/templates/molecular_viewer.html?aseview_styles_docs=20260510';

  const MOLECULE = {
    symbols: ['O', 'C', 'N', 'C', 'H', 'H', 'H', 'H', 'H'],
    positions: [
      [ 0.424546,  1.327024,  0.008034],
      [ 0.077158,  0.149789, -0.004249],
      [ 0.985518, -0.878537, -0.048910],
      [-1.371475, -0.288665, -0.000144],
      [ 0.707952, -1.824249,  0.169942],
      [-1.997229,  0.584922, -0.175477],
      [-1.560842, -1.039270, -0.771686],
      [-1.632113, -0.723007,  0.969814],
      [ 1.953133, -0.631574,  0.111866]
    ]
  };

  function titleCase(value) {
    if (value === '2d') return '2D';
    return value.slice(0, 1).toUpperCase() + value.slice(1);
  }

  function withPreviewChrome(html) {
    const base = `<base href="${TEMPLATE_URL}">`;
    const previewCss = `
      <style>
        .sidebar,
        .sidebar-toggle,
        .toolbar {
          display: none !important;
        }
        .app,
        .main-content,
        .viewer-container,
        .viewer {
          width: 100vw !important;
          height: 100vh !important;
        }
      </style>`;

    let next = html.replace(
      /https:\/\/cdn\.jsdelivr\.net\/gh\/kangmg\/aseview@main\/aseview\/static\/js\/styles\.js(?:\?[^"']*)?/g,
      'https://cdn.jsdelivr.net/gh/kangmg/aseview@main/aseview/static/js/styles.js?aseview_styles_docs=20260510'
    );

    if (!/<base\s/i.test(next)) {
      next = next.replace(/<head([^>]*)>/i, `<head$1>${base}`);
    }
    return next.replace('</head>', `${previewCss}</head>`);
  }

  function postPreviewData(iframe, style) {
    const dark = DARK_BACKGROUNDS.has(style);
    const backgroundColor = dark ? '#05070d' : '#ffffff';
    const cellColor = dark ? '#8a93a6' : '#6b7280';

    iframe.contentWindow.postMessage({
      type: 'setSettings',
      settings: {
        style,
        backgroundColor,
        cellColor,
        showCell: false,
        showBond: true,
        showAxis: false,
        atomSize: 0.42,
        bondThickness: 0.09,
        bondThreshold: 1.2,
        radiusContrast: 0.85,
        radiusContrastMode: 'log',
        showEnergyPlot: false,
        showForceMaxPlot: false
      }
    }, '*');

    iframe.contentWindow.postMessage({
      type: 'setData',
      data: [MOLECULE]
    }, '*');
  }

  async function loadTemplate() {
    const response = await fetch(TEMPLATE_URL, { cache: 'force-cache' });
    if (!response.ok) {
      throw new Error(`Failed to load molecular viewer template: ${response.status}`);
    }
    return withPreviewChrome(await response.text());
  }

  async function init() {
    const grid = document.getElementById('aseview-style-grid');
    if (!grid) return;

    grid.innerHTML = STYLES.map(style => {
      const dark = DARK_BACKGROUNDS.has(style);
      return `
        <section class="aseview-style-card ${dark ? 'is-dark' : 'is-light'}" data-style="${style}">
          <div class="aseview-style-title">${titleCase(style)}</div>
          <iframe title="${titleCase(style)} style preview" loading="lazy"></iframe>
        </section>`;
    }).join('');

    const srcdoc = await loadTemplate();
    const frames = Array.from(grid.querySelectorAll('iframe'));
    const pending = new Map();

    window.addEventListener('message', event => {
      if (!event.data || event.data.type !== 'viewerLoaded') return;
      const record = pending.get(event.source);
      if (record) postPreviewData(record.frame, record.style);
    });

    frames.forEach(frame => {
      const style = frame.closest('.aseview-style-card').dataset.style;
      pending.set(frame.contentWindow, { frame, style });
      frame.addEventListener('load', () => {
        setTimeout(() => postPreviewData(frame, style), 350);
      }, { once: true });
      frame.srcdoc = srcdoc;
    });
  }

  init().catch(err => {
    const grid = document.getElementById('aseview-style-grid');
    if (grid) grid.innerHTML = `<div class="aseview-style-error">${String(err.message || err)}</div>`;
  });
})();
</script>

<style>
.aseview-style-grid {
  display: grid;
  grid-template-columns: repeat(4, minmax(0, 1fr));
  gap: 0.85rem;
  margin: 1.25rem 0 2rem;
}

.aseview-style-card {
  border: 1px solid var(--md-default-fg-color--lightest);
  border-radius: 8px;
  overflow: hidden;
  background: var(--md-default-bg-color);
  box-shadow: 0 10px 28px rgba(0, 0, 0, 0.08);
}

.aseview-style-card.is-light {
  background: #ffffff;
  border-color: #e5e7eb;
}

.aseview-style-card.is-dark {
  background: #05070d;
  border-color: #1f2937;
}

.aseview-style-title {
  padding: 0.5rem 0.7rem;
  font-size: 0.78rem;
  font-weight: 700;
  letter-spacing: 0;
  border-bottom: 1px solid var(--md-default-fg-color--lightest);
}

.aseview-style-card.is-light .aseview-style-title {
  color: #111827;
  border-bottom-color: #e5e7eb;
}

.aseview-style-card.is-dark .aseview-style-title {
  color: #f8fafc;
  border-bottom-color: #1f2937;
}

.aseview-style-card iframe {
  display: block;
  width: 100%;
  height: 280px;
  border: 0;
  background: transparent;
}

.aseview-style-error {
  grid-column: 1 / -1;
  padding: 1rem;
  border: 1px solid #ef4444;
  border-radius: 8px;
  color: #991b1b;
  background: #fee2e2;
}

@media (max-width: 1200px) {
  .aseview-style-grid {
    grid-template-columns: repeat(2, minmax(0, 1fr));
  }
}

@media (max-width: 720px) {
  .aseview-style-grid {
    grid-template-columns: 1fr;
  }
}
</style>
