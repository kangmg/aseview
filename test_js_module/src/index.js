/**
 * ASEView.js - JavaScript library for molecular visualization
 *
 * Usage:
 * <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
 * <script src="https://unpkg.com/three@0.128.0/examples/js/controls/TrackballControls.js"></script>
 * <script src="aseview.umd.js"></script>
 * <script>
 *   const viewer = new ASEView.MolecularViewer('#container', {
 *     style: 'cartoon',
 *     colorBy: 'Element'
 *   });
 *   viewer.setData({
 *     symbols: ['O', 'H', 'H'],
 *     positions: [[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]]
 *   });
 * </script>
 */

// Core
export { SceneManager } from './core/Scene.js';

// Viewers
export { MolecularViewer } from './viewers/MolecularViewer.js';

// Renderers
export { createAtom, getChargeColor } from './renderers/AtomRenderer.js';
export { createBond } from './renderers/BondRenderer.js';

// Utils
export { atomInfo, getAtomInfo, detectBonds } from './utils/atomInfo.js';

// Version
export const version = '0.1.0';
