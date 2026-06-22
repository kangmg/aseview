import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
ASEVIEW_JS = ROOT / "aseview" / "static" / "js" / "aseview.js"


NODE_HARNESS = r"""
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');

function createHarness(scriptPath) {
  const posted = [];
  const frames = [];
  const listeners = { message: [] };

  const fakeWindow = {
    location: { href: 'https://example.test/app/index.html' },
    addEventListener(type, handler) {
      if (!listeners[type]) listeners[type] = [];
      listeners[type].push(handler);
    },
    removeEventListener(type, handler) {
      if (!listeners[type]) return;
      listeners[type] = listeners[type].filter((item) => item !== handler);
    },
  };
  fakeWindow.window = fakeWindow;

  const container = {
    innerHTML: '',
    child: null,
    appendChild(element) {
      this.child = element;
    },
  };

  const fakeDocument = {
    currentScript: { src: 'https://example.test/aseview/static/js/aseview.js' },
    querySelector(selector) {
      assert.equal(selector, '#container');
      return container;
    },
    createElement(tagName) {
      assert.equal(tagName, 'iframe');
      const frame = {
        style: {},
        attributes: {},
        listeners: {},
        srcdoc: '',
        contentWindow: {
          postMessage(message, targetOrigin) {
            posted.push({ frame, message, targetOrigin });
          },
        },
        setAttribute(name, value) {
          this.attributes[name] = value;
        },
        addEventListener(type, handler) {
          this.listeners[type] = handler;
        },
        removeEventListener(type, handler) {
          if (this.listeners[type] === handler) {
            delete this.listeners[type];
          }
        },
        remove() {
          this.removed = true;
        },
      };
      frames.push(frame);
      return frame;
    },
  };

  const context = {
    console,
    document: fakeDocument,
    fetch: async () => ({
      ok: true,
      text: async () => '<!doctype html><html><head></head><body></body></html>',
    }),
    setTimeout,
    clearTimeout,
    URL,
    window: fakeWindow,
  };
  context.globalThis = context;
  vm.createContext(context);
  vm.runInContext(fs.readFileSync(scriptPath, 'utf8'), context, { filename: scriptPath });

  function dispatchFrom(viewer, data) {
    const event = { source: viewer.iframe.contentWindow, data };
    for (const handler of listeners.message.slice()) {
      handler(event);
    }
  }

  function markReady(viewer) {
    dispatchFrom(viewer, { type: 'viewerLoaded' });
  }

  function drainPosted() {
    return posted.splice(0).map((entry) => JSON.parse(JSON.stringify(entry.message)));
  }

  return {
    ASEView: fakeWindow.ASEView,
    container,
    dispatchFrom,
    drainPosted,
    frames,
    markReady,
    posted,
  };
}
"""


def _run_node(tmp_path, js_source):
    if not shutil.which("node"):
        pytest.skip("Node.js not available")

    script_path = tmp_path / "aseview_api_harness.js"
    script_path.write_text(NODE_HARNESS + "\n" + textwrap.dedent(js_source), encoding="utf-8")
    result = subprocess.run(
        ["node", str(script_path), str(ASEVIEW_JS)],
        cwd=ROOT,
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def test_existing_ready_settings_and_setdata_flow(tmp_path):
    _run_node(
        tmp_path,
        """
        const harness = createHarness(process.argv[2]);
        const viewer = new harness.ASEView.MolecularViewer('#container', {
          theme: 'light',
          showBonds: true,
        });

        viewer.setData({ symbols: ['He'], positions: [[0, 0, 0]] });
        assert.equal(harness.posted.length, 0, 'data should queue before viewerLoaded');

        harness.markReady(viewer);
        const readyMessages = harness.drainPosted();
        assert.equal(readyMessages[0].type, 'setSettings');
        assert.equal(readyMessages[0].settings.showBonds, true);
        assert.equal(readyMessages[1].type, 'setData');
        assert.deepEqual(readyMessages[1].data[0].symbols, ['He']);

        viewer.setSettings({ showBonds: false });
        const settingsMessages = harness.drainPosted();
        assert.deepEqual(settingsMessages, [
          { type: 'setSettings', settings: { showBonds: false } },
        ]);
        """,
    )


def test_camera_export_request_response_contract(tmp_path):
    _run_node(
        tmp_path,
        """
        (async () => {
          const harness = createHarness(process.argv[2]);
          const molecular = new harness.ASEView.MolecularViewer('#container', {
            requestTimeout: 25,
          });
          harness.markReady(molecular);
          harness.drainPosted();

          const overlay = new harness.ASEView.OverlayViewer('#container', {
            requestTimeout: 25,
          });
          harness.markReady(overlay);
          harness.drainPosted();

          const normal = new harness.ASEView.NormalModeViewer('#container', {
            requestTimeout: 25,
          });
          harness.markReady(normal);
          harness.drainPosted();

          for (const viewer of [molecular, overlay, normal]) {
            assert.equal(typeof viewer.setView, 'function', 'setView should be public');
            assert.equal(typeof viewer.resetView, 'function', 'resetView should be public');
            assert.equal(typeof viewer.savePNG, 'function', 'savePNG should be public');
            assert.equal(typeof viewer.saveGIF, 'function', 'saveGIF should be public');
          }

          const viewPromise = molecular.setView({ preset: 'top' });
          const setViewMessage = harness.drainPosted().find((message) => message.type === 'setView');
          assert.ok(setViewMessage.requestId, 'setView should include a requestId');
          assert.deepEqual(setViewMessage.view, { preset: 'top' });
          assert.equal(molecular._pendingRequests.size, 1);
          harness.dispatchFrom(molecular, {
            type: 'viewResult',
            requestId: setViewMessage.requestId,
            ok: true,
            view: { preset: 'top', applied: true },
          });
          assert.deepEqual(JSON.parse(JSON.stringify(await viewPromise)), {
            ok: true,
            type: 'view',
            view: { preset: 'top', applied: true },
          });
          assert.equal(molecular._pendingRequests.size, 0);

          const directionPromise = molecular.setView({ direction: [1, 0, 0], up: [0, 0, 1] });
          const directionMessage = harness.drainPosted().find((message) => message.type === 'setView');
          harness.dispatchFrom(molecular, {
            type: 'viewResult',
            requestId: directionMessage.requestId,
            ok: true,
            view: { direction: [1, 0, 0], up: [0, 0, 1] },
          });
          await directionPromise;

          const presetPromise = molecular.setView({ preset: 'side-a' });
          const presetMessage = harness.drainPosted().find((message) => message.type === 'setView');
          harness.dispatchFrom(molecular, {
            type: 'viewResult',
            requestId: presetMessage.requestId,
            ok: true,
            view: { preset: 'side-a' },
          });
          await presetPromise;
          assert.equal(molecular.options.viewPreset, 'side-a');
          assert.equal(
            Object.prototype.hasOwnProperty.call(molecular.options, 'viewDirection'),
            false,
            'setView preset should clear stale direction options',
          );

          const pngPromise = molecular.savePNG({
            download: false,
            returnDataUrl: true,
            scale: 2,
          });
          const pngMessage = harness.drainPosted().find((message) => message.type === 'exportImage');
          assert.equal(pngMessage.format, 'png');
          assert.ok(pngMessage.requestId, 'savePNG should include a requestId');
          assert.equal(pngMessage.options.download, false);
          assert.equal(pngMessage.options.returnDataUrl, true);
          assert.equal(pngMessage.options.scale, 2);
          harness.dispatchFrom(molecular, {
            type: 'exportResult',
            requestId: pngMessage.requestId,
            ok: true,
            format: 'png',
            filename: 'molecule.png',
            dataUrl: 'data:image/png;base64,abc',
            width: 640,
            height: 480,
          });
          assert.deepEqual(JSON.parse(JSON.stringify(await pngPromise)), {
            ok: true,
            type: 'png',
            filename: 'molecule.png',
            dataUrl: 'data:image/png;base64,abc',
            width: 640,
            height: 480,
          });
          assert.equal(molecular._pendingRequests.size, 0);

          const failedPngPromise = molecular.savePNG({ download: false });
          const failedPngMessage = harness.drainPosted().find((message) => message.type === 'exportImage');
          harness.dispatchFrom(molecular, {
            type: 'exportResult',
            requestId: failedPngMessage.requestId,
            ok: false,
            code: 'render_failed',
            message: 'No canvas available',
          });
          await assert.rejects(
            failedPngPromise,
            (error) => {
              assert.equal(error.code, 'render_failed');
              assert.equal(error.type, 'png');
              assert.equal(error.requestId, failedPngMessage.requestId);
              assert.match(error.message, /No canvas/);
              return true;
            },
          );
          assert.equal(molecular._pendingRequests.size, 0);

          const gifPromise = normal.saveGIF({
            returnDataUrl: true,
            frames: 3,
            delay: 25,
            sampleInterval: 5,
          });
          const gifMessage = harness.drainPosted().find((message) => message.type === 'exportImage');
          assert.equal(gifMessage.format, 'gif');
          assert.equal(gifMessage.options.frames, 3);
          assert.equal(gifMessage.options.delay, 25);
          assert.equal(gifMessage.options.sampleInterval, 5);
          harness.dispatchFrom(normal, {
            type: 'exportResult',
            requestId: gifMessage.requestId,
            ok: true,
            format: 'gif',
            filename: 'normal.gif',
            dataUrl: 'data:image/gif;base64,abc',
          });
          assert.deepEqual(JSON.parse(JSON.stringify(await gifPromise)), {
            ok: true,
            type: 'gif',
            filename: 'normal.gif',
            dataUrl: 'data:image/gif;base64,abc',
          });
          assert.equal(normal._pendingRequests.size, 0);

          await assert.rejects(
            overlay.saveGIF({ returnDataUrl: true }),
            (error) => {
              assert.equal(error.code, 'unsupported_export');
              assert.equal(error.type, 'gif');
              assert.match(error.message, /OverlayViewer/);
              return true;
            },
          );
          assert.equal(
            harness.drainPosted().some((message) => message.type === 'exportImage' && message.format === 'gif'),
            false,
            'Overlay GIF rejection should not post an iframe export command',
          );

          await assert.rejects(
            molecular.setView(['top']),
            (error) => {
              assert.equal(error.code, 'invalid_view');
              assert.equal(error.type, 'view');
              return true;
            },
          );
          await assert.rejects(
            molecular.savePNG({ width: 10.5 }),
            (error) => {
              assert.equal(error.code, 'invalid_options');
              assert.equal(error.type, 'png');
              assert.match(error.message, /width/);
              return true;
            },
          );
          assert.equal(molecular._pendingRequests.size, 0);

          const timeoutPromise = molecular.resetView();
          const resetMessage = harness.drainPosted().find((message) => message.type === 'resetView');
          assert.ok(resetMessage.requestId, 'resetView should include a requestId');
          await assert.rejects(
            timeoutPromise,
            (error) => {
              assert.equal(error.code, 'request_timeout');
              assert.equal(error.type, 'view');
              assert.equal(error.requestId, resetMessage.requestId);
              return true;
            },
          );
          assert.equal(molecular._pendingRequests.size, 0);

          harness.dispatchFrom(molecular, {
            type: 'viewResult',
            requestId: resetMessage.requestId,
            ok: true,
            view: { preset: 'stale' },
          });
          assert.equal(molecular._pendingRequests.size, 0, 'stale timeout response should stay ignored');

          const iframeErrorPromise = molecular.savePNG({ returnDataUrl: true });
          const iframeErrorMessage = harness.drainPosted().find((message) => message.type === 'exportImage');
          harness.dispatchFrom(molecular, {
            type: 'error',
            requestId: iframeErrorMessage.requestId,
            code: 'iframe_error',
            message: 'Renderer crashed',
          });
          await assert.rejects(
            iframeErrorPromise,
            (error) => {
              assert.equal(error.code, 'iframe_error');
              assert.equal(error.type, 'png');
              assert.equal(error.requestId, iframeErrorMessage.requestId);
              assert.match(error.message, /Renderer crashed/);
              return true;
            },
          );
          assert.equal(molecular._pendingRequests.size, 0);
        })().catch((error) => {
          console.error(error && error.stack ? error.stack : error);
          process.exitCode = 1;
        });
        """,
    )


def test_set_view_options_commit_only_after_iframe_success(tmp_path):
    _run_node(
        tmp_path,
        """
        (async () => {
          const harness = createHarness(process.argv[2]);

          const rejectedViewer = new harness.ASEView.MolecularViewer('#container', {
            requestTimeout: 25,
            viewPreset: 'top',
          });
          harness.markReady(rejectedViewer);
          harness.drainPosted();

          const malformedPromise = rejectedViewer.setView(['top']);
          await assert.rejects(
            malformedPromise,
            (error) => {
              assert.equal(error.code, 'invalid_view');
              assert.equal(error.type, 'view');
              return true;
            },
          );
          assert.equal(rejectedViewer.options.viewPreset, 'top');
          assert.equal(
            Object.prototype.hasOwnProperty.call(rejectedViewer.options, 'viewDirection'),
            false,
            'malformed setView input must not mutate local view options',
          );

          const rejectedPromise = rejectedViewer.setView({ direction: [0, 0, 0] });
          const rejectedMessage = harness.drainPosted().find((message) => message.type === 'setView');
          assert.ok(rejectedMessage.requestId);
          assert.equal(
            rejectedViewer.options.viewPreset,
            'top',
            'setView should not mutate local options before iframe success',
          );
          assert.equal(
            Object.prototype.hasOwnProperty.call(rejectedViewer.options, 'viewDirection'),
            false,
          );

          harness.dispatchFrom(rejectedViewer, {
            type: 'viewResult',
            requestId: rejectedMessage.requestId,
            ok: false,
            code: 'invalid_view',
            message: 'zero direction',
          });
          await assert.rejects(
            rejectedPromise,
            (error) => {
              assert.equal(error.code, 'invalid_view');
              assert.equal(error.type, 'view');
              assert.equal(error.requestId, rejectedMessage.requestId);
              return true;
            },
          );
          assert.equal(rejectedViewer.options.viewPreset, 'top');
          assert.equal(
            Object.prototype.hasOwnProperty.call(rejectedViewer.options, 'viewDirection'),
            false,
            'rejected setView must not leave stale local direction state',
          );
          assert.equal(rejectedViewer._pendingRequests.size, 0);

          const timeoutViewer = new harness.ASEView.MolecularViewer('#container', {
            requestTimeout: 5,
          });
          const timeoutPromise = timeoutViewer.setView({ preset: 'side-a' });
          assert.equal(
            Object.prototype.hasOwnProperty.call(timeoutViewer.options, 'viewPreset'),
            false,
            'pre-ready setView should not mutate local options before iframe success',
          );
          await assert.rejects(
            timeoutPromise,
            (error) => {
              assert.equal(error.code, 'request_timeout');
              assert.equal(error.type, 'view');
              return true;
            },
          );
          assert.equal(timeoutViewer._pendingRequests.size, 0);
          assert.equal(timeoutViewer.pendingMessages.length, 0);
          assert.equal(
            Object.prototype.hasOwnProperty.call(timeoutViewer.options, 'viewPreset'),
            false,
            'timed-out pre-ready setView must not leave stale local preset state',
          );

          harness.markReady(timeoutViewer);
          const readyMessages = harness.drainPosted();
          const settingsMessage = readyMessages.find((message) => message.type === 'setSettings');
          assert.ok(settingsMessage, 'viewerLoaded should still send base settings');
          assert.equal(
            Object.prototype.hasOwnProperty.call(settingsMessage.settings, 'viewPreset'),
            false,
            'viewerLoaded settings must not include a stale timed-out viewPreset',
          );
          assert.equal(
            readyMessages.some((message) => message.type === 'setView'),
            false,
            'timed-out pre-ready setView command should be removed from the queue',
          );

          const successViewer = new harness.ASEView.MolecularViewer('#container', {
            requestTimeout: 25,
            viewPreset: 'front',
          });
          harness.markReady(successViewer);
          harness.drainPosted();

          const successPromise = successViewer.setView({ direction: [1, 0, 0], up: [0, 0, 1] });
          const successMessage = harness.drainPosted().find((message) => message.type === 'setView');
          assert.equal(successViewer.options.viewPreset, 'front');
          assert.equal(
            Object.prototype.hasOwnProperty.call(successViewer.options, 'viewDirection'),
            false,
            'successful setView should update local options only after iframe response',
          );
          harness.dispatchFrom(successViewer, {
            type: 'viewResult',
            requestId: successMessage.requestId,
            ok: true,
            view: { direction: [1, 0, 0], up: [0, 0, 1] },
          });
          await successPromise;
          assert.equal(
            Object.prototype.hasOwnProperty.call(successViewer.options, 'viewPreset'),
            false,
            'successful direction setView should clear stale local preset state',
          );
          assert.deepEqual(JSON.parse(JSON.stringify(successViewer.options.viewDirection)), [1, 0, 0]);
          assert.deepEqual(JSON.parse(JSON.stringify(successViewer.options.viewUp)), [0, 0, 1]);
          assert.equal(successViewer._pendingRequests.size, 0);
        })().catch((error) => {
          console.error(error && error.stack ? error.stack : error);
          process.exitCode = 1;
        });
        """,
    )


def test_reset_view_options_commit_only_after_iframe_success(tmp_path):
    _run_node(
        tmp_path,
        """
        (async () => {
          const harness = createHarness(process.argv[2]);

          const rejectedViewer = new harness.ASEView.MolecularViewer('#container', {
            requestTimeout: 25,
            viewPreset: 'side-a',
            viewUp: [0, 0, 1],
          });
          harness.markReady(rejectedViewer);
          harness.drainPosted();

          const rejectedPromise = rejectedViewer.resetView();
          const rejectedMessage = harness.drainPosted().find((message) => message.type === 'resetView');
          assert.ok(rejectedMessage.requestId);
          assert.equal(
            rejectedViewer.options.viewPreset,
            'side-a',
            'resetView should not clear local options before iframe success',
          );
          assert.deepEqual(JSON.parse(JSON.stringify(rejectedViewer.options.viewUp)), [0, 0, 1]);
          harness.dispatchFrom(rejectedViewer, {
            type: 'viewResult',
            requestId: rejectedMessage.requestId,
            ok: false,
            code: 'render_failed',
            message: 'reset failed',
          });
          await assert.rejects(
            rejectedPromise,
            (error) => {
              assert.equal(error.code, 'render_failed');
              assert.equal(error.type, 'view');
              assert.equal(error.requestId, rejectedMessage.requestId);
              return true;
            },
          );
          assert.equal(rejectedViewer.options.viewPreset, 'side-a');
          assert.deepEqual(JSON.parse(JSON.stringify(rejectedViewer.options.viewUp)), [0, 0, 1]);

          const timeoutViewer = new harness.ASEView.MolecularViewer('#container', {
            requestTimeout: 5,
            viewDirection: [1, 0, 0],
            viewEuler: [0, 90, 0],
            viewUp: [0, 0, 1],
          });
          harness.markReady(timeoutViewer);
          harness.drainPosted();

          const timeoutPromise = timeoutViewer.resetView();
          const timeoutMessage = harness.drainPosted().find((message) => message.type === 'resetView');
          assert.ok(timeoutMessage.requestId);
          await assert.rejects(
            timeoutPromise,
            (error) => {
              assert.equal(error.code, 'request_timeout');
              assert.equal(error.type, 'view');
              assert.equal(error.requestId, timeoutMessage.requestId);
              return true;
            },
          );
          assert.deepEqual(JSON.parse(JSON.stringify(timeoutViewer.options.viewDirection)), [1, 0, 0]);
          assert.deepEqual(JSON.parse(JSON.stringify(timeoutViewer.options.viewEuler)), [0, 90, 0]);
          assert.deepEqual(JSON.parse(JSON.stringify(timeoutViewer.options.viewUp)), [0, 0, 1]);
          assert.equal(timeoutViewer._pendingRequests.size, 0);

          harness.dispatchFrom(timeoutViewer, {
            type: 'viewResult',
            requestId: timeoutMessage.requestId,
            ok: true,
            view: {},
          });
          assert.deepEqual(
            JSON.parse(JSON.stringify(timeoutViewer.options.viewDirection)),
            [1, 0, 0],
            'late success response after timeout must stay ignored',
          );

          const successViewer = new harness.ASEView.MolecularViewer('#container', {
            requestTimeout: 25,
            viewPreset: 'front',
            viewDirection: [0, 1, 0],
            viewEuler: [10, 20, 30],
            viewUp: [0, 0, 1],
          });
          harness.markReady(successViewer);
          harness.drainPosted();

          const successPromise = successViewer.resetView();
          const successMessage = harness.drainPosted().find((message) => message.type === 'resetView');
          assert.equal(successViewer.options.viewPreset, 'front');
          assert.deepEqual(JSON.parse(JSON.stringify(successViewer.options.viewDirection)), [0, 1, 0]);
          harness.dispatchFrom(successViewer, {
            type: 'viewResult',
            requestId: successMessage.requestId,
            ok: true,
            view: {},
          });
          await successPromise;
          for (const key of ['viewPreset', 'viewDirection', 'viewEuler', 'viewUp']) {
            assert.equal(
              Object.prototype.hasOwnProperty.call(successViewer.options, key),
              false,
              `successful resetView should clear ${key}`,
            );
          }
          assert.equal(successViewer._pendingRequests.size, 0);
        })().catch((error) => {
          console.error(error && error.stack ? error.stack : error);
          process.exitCode = 1;
        });
        """,
    )
