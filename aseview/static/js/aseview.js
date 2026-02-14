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

    // CDN base URL for templates (using raw.githack.com which properly serves HTML with correct content-type)
    // raw.githack.com serves files with correct MIME types, unlike jsDelivr which serves HTML as text
    const CDN_BASE = 'https://raw.githack.com/kangmg/aseview_v2_dev/main/aseview/templates';

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

            this.options = options;
            this.templateName = templateName;
            this.iframe = null;
            this.isReady = false;
            this.pendingMessages = [];

            this._createIframe();
        }

        _createIframe() {
            // Clear container
            this.container.innerHTML = '';

            // Create iframe
            this.iframe = document.createElement('iframe');
            this.iframe.style.cssText = 'width:100%;height:100%;border:none;display:block;';
            this.iframe.setAttribute('allowfullscreen', 'true');

            // Build template URL
            const templateUrl = `${CDN_BASE}/${this.templateName}`;
            this.iframe.src = templateUrl;

            // Listen for messages from iframe
            this._messageHandler = (event) => this._onMessage(event);
            window.addEventListener('message', this._messageHandler);

            // Add to container
            this.container.appendChild(this.iframe);
        }

        _onMessage(event) {
            const msg = event.data;
            if (!msg || !msg.type) return;

            if (msg.type === 'viewerLoaded') {
                this.isReady = true;
                // Send any pending messages
                this.pendingMessages.forEach(m => this._postMessage(m));
                this.pendingMessages = [];
                // Apply initial settings if any
                if (Object.keys(this.options).length > 0) {
                    this._postMessage({ type: 'setSettings', settings: this.options });
                }
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
     * Features: Structure opacity controls, alignment algorithms (Kabsch, Hungarian),
     * RMSD calculation, and multi-structure visualization.
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

    // Export to global namespace
    global.ASEView = {
        MolecularViewer,
        NormalModeViewer,
        OverlayViewer,
        version: '1.0.0',
        CDN_BASE
    };

})(typeof window !== 'undefined' ? window : this);
