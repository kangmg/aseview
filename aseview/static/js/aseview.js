/**
 * ASEView.js - JavaScript API for molecular visualization
 *
 * Uses the same templates as the Python package via iframe embedding.
 * This ensures consistent UI and behavior between Python and web usage.
 *
 * Usage:
 *   const viewer = new ASEView.MolecularViewer('#container');
 *   viewer.setData({ symbols: ['O','H','H'], positions: [[0,0,0],[0.96,0,0],[-0.24,0.93,0]] });
 */

(function(global) {
    'use strict';

    // CDN base for themes. Templates are fetched as text and injected via srcdoc
    // so they are not blocked by CDN iframe headers or text/plain content types.
    // Themes live at: aseview/themes/{theme}/{template}.html
    // Fallback to legacy templates/ path for backward compat
    const REMOTE_THEMES_BASE = 'https://cdn.jsdelivr.net/gh/kangmg/aseview@main/aseview/themes';
    const REMOTE_FALLBACK_BASE = 'https://cdn.jsdelivr.net/gh/kangmg/aseview@main/aseview/templates';

    function inferLocalTemplateBases() {
        const currentScript = document.currentScript;
        if (!currentScript || !currentScript.src) return null;
        try {
            const scriptUrl = new URL(currentScript.src, global.location && global.location.href ? global.location.href : undefined);
            if (scriptUrl.protocol === 'file:') return null;
            return {
                themes: new URL('../../themes', scriptUrl).href.replace(/\/$/, ''),
                templates: new URL('../../templates', scriptUrl).href.replace(/\/$/, '')
            };
        } catch (err) {
            return null;
        }
    }

    const localTemplateBases = inferLocalTemplateBases();
    const THEMES_BASE = (global.ASEVIEW_THEMES_BASE || (localTemplateBases && localTemplateBases.themes) || REMOTE_THEMES_BASE);
    const FALLBACK_BASE = (global.ASEVIEW_TEMPLATES_BASE || (localTemplateBases && localTemplateBases.templates) || REMOTE_FALLBACK_BASE);
    const CDN_BASE       = THEMES_BASE;   // kept for external access via ASEView.CDN_BASE
    const DEFAULT_CACHE_BUST = 'auto';
    const ASEVIEW_MAIN_ASSET_RE = /https:\/\/cdn\.jsdelivr\.net\/gh\/kangmg\/aseview@main\/aseview\/static\/js\/styles\.js(?:\?[^"']*)?/g;
    const DEFAULT_VIEWER_OPTIONS = {
        bondThreshold: 1.2,
    };

    let _defaultTheme = 'dark';

    function setTheme(name) { _defaultTheme = name; }
    function getTheme()     { return _defaultTheme; }

    /**
     * Base class for all viewers
     */
    class BaseViewer {
        constructor(container, templateName, options = {}) {
            // Get container element
            if (typeof container === 'string') {
                this.container = document.querySelector(container);
            } else {
                this.container = container;
            }

            if (!this.container) {
                throw new Error('ASEView: Container element not found');
            }

            this.options = Object.assign({}, DEFAULT_VIEWER_OPTIONS, options);
            this.templateName = templateName;
            this._theme = this.options.theme || _defaultTheme;
            this.iframe = null;
            this.isReady = false;
            this.pendingMessages = [];
            this._loadId = 0;
            this._cacheBustToken = this._resolveCacheBustToken(
                Object.prototype.hasOwnProperty.call(this.options, 'cacheBust')
                    ? this.options.cacheBust
                    : DEFAULT_CACHE_BUST
            );

            this._createIframe();
        }

        _createIframe() {
            // Clear container
            this.container.innerHTML = '';

            // Create iframe
            this.iframe = document.createElement('iframe');
            this.iframe.style.cssText = 'width:100%;height:100%;border:none;display:block;';
            this.iframe.setAttribute('allowfullscreen', 'true');

            // Listen for messages from iframe
            this._messageHandler = (event) => this._onMessage(event);
            window.addEventListener('message', this._messageHandler);

            // Add to container
            this.container.appendChild(this.iframe);

            this._loadTemplate();
        }

        _templateCandidates() {
            return [
                `${THEMES_BASE}/${this._theme}/${this.templateName}`,
                `${FALLBACK_BASE}/${this.templateName}`,
            ];
        }

        _withBaseHref(html, href) {
            const base = `<base href="${href}">`;
            if (/<base\s/i.test(html)) return html;
            if (/<head[^>]*>/i.test(html)) {
                return html.replace(/<head([^>]*)>/i, `<head$1>${base}`);
            }
            return `${base}${html}`;
        }

        _resolveCacheBustToken(value) {
            if (!value) return null;
            if (value === true || value === 'auto') {
                return `${Date.now().toString(36)}-${Math.random().toString(36).slice(2, 8)}`;
            }
            return String(value);
        }

        _withCacheBust(url) {
            if (!this._cacheBustToken) return url;
            try {
                const next = new URL(url, global.location && global.location.href ? global.location.href : undefined);
                next.searchParams.set('_aseview', this._cacheBustToken);
                return next.href;
            } catch (err) {
                const separator = url.includes('?') ? '&' : '?';
                return `${url}${separator}_aseview=${encodeURIComponent(this._cacheBustToken)}`;
            }
        }

        _withAssetCacheBust(html) {
            if (!this._cacheBustToken) return html;
            return html.replace(ASEVIEW_MAIN_ASSET_RE, (url) => this._withCacheBust(url));
        }

        async _loadTemplate() {
            const loadId = ++this._loadId;
            const candidates = this._templateCandidates();
            let lastError = null;

            for (const url of candidates) {
                try {
                    const requestUrl = this._withCacheBust(url);
                    const response = await fetch(requestUrl, { cache: this._cacheBustToken ? 'reload' : 'no-cache' });
                    if (!response.ok) {
                        throw new Error(`HTTP ${response.status}`);
                    }
                    const html = this._withAssetCacheBust(await response.text());
                    if (loadId !== this._loadId || !this.iframe) return;
                    this.iframe.srcdoc = this._withBaseHref(html, url);
                    return;
                } catch (err) {
                    lastError = err;
                }
            }

            console.error(
                `ASEView: failed to load ${this.templateName}`,
                lastError && lastError.message ? lastError.message : lastError
            );
        }

        _onMessage(event) {
            const msg = event.data;
            if (!msg || !msg.type) return;
            if (!this.iframe || event.source !== this.iframe.contentWindow) return;

            if (msg.type === 'viewerLoaded') {
                this.isReady = true;
                // Apply initial settings before queued data messages.
                if (Object.keys(this.options).length > 0) {
                    this._postMessage({ type: 'setSettings', settings: this.options });
                }
                this.pendingMessages.forEach(m => this._postMessage(m));
                this.pendingMessages = [];
            }
        }

        _postMessage(message) {
            if (this.isReady && this.iframe && this.iframe.contentWindow) {
                this.iframe.contentWindow.postMessage(message, '*');
            } else {
                this.pendingMessages.push(message);
            }
        }

        /**
         * Update viewer settings
         * @param {Object} settings - Settings to update
         */
        setSettings(settings) {
            Object.assign(this.options, settings);
            this._postMessage({ type: 'setSettings', settings });
        }

        /**
         * Dispose the viewer and clean up resources
         */
        dispose() {
            window.removeEventListener('message', this._messageHandler);
            this._loadId++;
            if (this.iframe) {
                this.iframe.remove();
                this.iframe = null;
            }
            this.container.innerHTML = '';
        }
    }

    /**
     * MolecularViewer - Full-featured molecular structure viewer
     *
     * Features: Style selector, bond/cell toggles, atom size controls,
     * trajectory animation, energy plot, and more.
     */
    class MolecularViewer extends BaseViewer {
        constructor(container, options = {}) {
            super(container, 'molecular_viewer.html', options);
        }

        /**
         * Set molecular data
         * @param {Object|Array} data - Molecular data (single structure or trajectory)
         * @param {string[]} data.symbols - Atom symbols
         * @param {number[][]} data.positions - Atom positions
         * @param {number[][]} [data.cell] - Unit cell vectors (optional)
         */
        setData(data) {
            // Ensure data is array format (for trajectory support)
            const dataArray = Array.isArray(data) ? data : [data];
            this._postMessage({ type: 'setData', data: dataArray });
        }
    }

    /**
     * NormalModeViewer - Viewer for molecular vibration animations
     *
     * Features: Mode selector, amplitude control, animation speed,
     * play/pause, and vibration visualization.
     */
    class NormalModeViewer extends BaseViewer {
        constructor(container, options = {}) {
            super(container, 'normal_viewer.html', options);
        }

        /**
         * Set equilibrium structure (for trajectory mode)
         * @param {Object|Array} data - Molecular data
         */
        setData(data) {
            const dataArray = Array.isArray(data) ? data : [data];
            this._postMessage({ type: 'setData', data: dataArray });
        }

        /**
         * Initialize normal mode viewer with vibration data
         * @param {Object} atoms - Equilibrium structure { symbols, positions }
         * @param {Object} vibrationData - Vibration data
         * @param {number[][]} vibrationData.modeVectors - Displacement vectors for each mode
         * @param {string[]} vibrationData.frequencies - Frequency labels (e.g., "1234.56" or "123.45i")
         * @param {boolean[]} vibrationData.isImaginary - Whether each mode is imaginary
         */
        setVibrationData(atoms, vibrationData) {
            this._postMessage({
                type: 'initNormalMode',
                atoms: atoms,
                vibrationData: vibrationData
            });
        }
    }

    /**
     * OverlayViewer - Viewer for comparing multiple molecular structures
     *
     * Features: Structure opacity controls, visibility toggles,
     * per-structure colors, and multi-structure visualization.
     */
    class OverlayViewer extends BaseViewer {
        constructor(container, options = {}) {
            super(container, 'overlay_viewer.html', options);
        }

        /**
         * Set structures for overlay comparison
         * @param {...Object} structures - Two or more structures to compare
         */
        setStructures(...structures) {
            this._postMessage({ type: 'setData', data: structures });
        }

        /**
         * Set molecular data (array of structures)
         * @param {Array} data - Array of molecular structures
         */
        setData(data) {
            const dataArray = Array.isArray(data) ? data : [data];
            this._postMessage({ type: 'setData', data: dataArray });
        }
    }

    /**
     * FragSelector - Interactive fragment selector with synchronized 2D/3D views
     *
     * Features: Click/rect/lasso atom selection, synchronized SVG 2D and THREE.js 3D
     * panels, structure copy in multiple formats, and clipboard copy of selected/unselected indices.
     */
    class FragSelector extends BaseViewer {
        constructor(container, options = {}) {
            super(container, 'frag_selector.html', options);
        }

        /**
         * Set molecular data for fragment selection
         * @param {Object} data - Single molecular structure { symbols, positions }
         * @param {string[]} data.symbols - Atom symbols
         * @param {number[][]} data.positions - Atom positions [[x,y,z], ...]
         */
        setData(data) {
            const dataArray = Array.isArray(data) ? data : [data];
            this._postMessage({ type: 'setData', data: dataArray });
        }

        /**
         * Get current selection (requires iframe to post back)
         * @returns {Promise<number[]>} Selected atom indices
         */
        getSelection() {
            return new Promise((resolve) => {
                const handler = (event) => {
                    const msg = event.data;
                    if (event.source !== (this.iframe && this.iframe.contentWindow)) return;
                    if (msg && msg.type === 'selectionResponse') {
                        window.removeEventListener('message', handler);
                        resolve(msg.selected || []);
                    }
                };
                window.addEventListener('message', handler);
                this._postMessage({ type: 'getSelection' });
            });
        }

        /**
         * Set selection programmatically
         * @param {number[]} indices - Atom indices to select
         */
        setSelection(indices) {
            this._postMessage({ type: 'setSelection', indices });
        }

        /**
         * Clear all selection
         */
        clearSelection() {
            this._postMessage({ type: 'clearSelection' });
        }
    }

    // Export to global namespace
    global.ASEView = {
        MolecularViewer,
        NormalModeViewer,
        OverlayViewer,
        FragSelector,
        setTheme,
        getTheme,
        DEFAULT_OPTIONS: Object.assign({}, DEFAULT_VIEWER_OPTIONS),
        version: '0.0.6',
        CDN_BASE,
    };

})(typeof window !== 'undefined' ? window : this);
