import * as path from "path";
import * as fs from "fs/promises";
import { spawn } from "child_process";
import * as vscode from "vscode";

const VIEW_TYPE = "aseview.viewer";

type ViewerFrame = {
  symbols: string[];
  positions: number[][];
  cell?: number[][];
  forces?: number[][];
  charges?: number[];
  energy?: number;
  name?: string;
};

class AseviewEditorProvider implements vscode.CustomTextEditorProvider {
  private readonly sessions = new Map<
    string,
    { document: vscode.TextDocument; panel: vscode.WebviewPanel }
  >();

  constructor(private readonly context: vscode.ExtensionContext) {}

  public static register(context: vscode.ExtensionContext): vscode.Disposable {
    const provider = new AseviewEditorProvider(context);
    const providerRegistration = vscode.window.registerCustomEditorProvider(
      VIEW_TYPE,
      provider,
      {
        webviewOptions: {
          retainContextWhenHidden: true,
        },
        supportsMultipleEditorsPerDocument: false,
      }
    );

    const openWith = vscode.commands.registerCommand(
      "aseview.openWith",
      async (arg?: unknown) => {
        const resource = coerceResourceUri(arg) ?? getActiveResourceUri();
        if (!resource) {
          void vscode.window.showWarningMessage(
            "No file selected. Open a structure file first."
          );
          return;
        }

        await vscode.commands.executeCommand("vscode.openWith", resource, VIEW_TYPE);
      }
    );

    const openActiveFile = vscode.commands.registerCommand(
      "aseview.openActiveFile",
      async () => {
        const resource = getActiveResourceUri();
        if (!resource) {
          void vscode.window.showWarningMessage(
            "No active file found. Open a structure file first."
          );
          return;
        }
        await vscode.commands.executeCommand("vscode.openWith", resource, VIEW_TYPE);
      }
    );

    const refresh = vscode.commands.registerCommand("aseview.refresh", async () => {
      const resource = getActiveResourceUri();
      if (!resource) {
        void vscode.window.showWarningMessage(
          "No active ASEView tab found to refresh."
        );
        return;
      }

      const refreshed = await provider.refresh(resource);
      if (!refreshed) {
        void vscode.window.showWarningMessage(
          "Current file is not open in ASEView. Use 'Open With ASEView' first."
        );
      }
    });

    return vscode.Disposable.from(providerRegistration, openWith, openActiveFile, refresh);
  }

  public async resolveCustomTextEditor(
    document: vscode.TextDocument,
    webviewPanel: vscode.WebviewPanel,
    _token: vscode.CancellationToken
  ): Promise<void> {
    webviewPanel.webview.options = {
      enableScripts: true,
      localResourceRoots: [
        this.context.extensionUri,
        vscode.Uri.joinPath(this.context.extensionUri, "media"),
      ],
    };

    const sessionKey = document.uri.toString();
    this.sessions.set(sessionKey, { document, panel: webviewPanel });

    const refreshPanel = async (): Promise<void> => {
      await this.updateWebview(document, webviewPanel);
    };

    const saveListener = vscode.workspace.onDidSaveTextDocument((saved) => {
      if (saved.uri.toString() === sessionKey) {
        void refreshPanel();
      }
    });

    const fileWatcher = vscode.workspace.createFileSystemWatcher(
      new vscode.RelativePattern(path.dirname(document.uri.fsPath), path.basename(document.uri.fsPath))
    );
    fileWatcher.onDidChange(() => void refreshPanel());
    fileWatcher.onDidCreate(() => void refreshPanel());

    webviewPanel.onDidDispose(() => {
      saveListener.dispose();
      fileWatcher.dispose();
      const current = this.sessions.get(sessionKey);
      if (current?.panel === webviewPanel) {
        this.sessions.delete(sessionKey);
      }
    });

    webviewPanel.webview.onDidReceiveMessage(async (msg: unknown) => {
      if (!msg || typeof msg !== "object") { return; }
      const m = msg as Record<string, unknown>;
      if (m["type"] !== "saveFile") { return; }

      const dataUrl = m["dataUrl"] as string | undefined;
      const suggestedName = typeof m["filename"] === "string" ? m["filename"] : "download";
      if (!dataUrl) { return; }

      const saveUri = await vscode.window.showSaveDialog({
        defaultUri: vscode.Uri.file(suggestedName),
        filters: suggestedName.endsWith(".gif")
          ? { "GIF Image": ["gif"] }
          : { "PNG Image": ["png"] },
      });
      if (!saveUri) { return; }

      try {
        const base64 = dataUrl.replace(/^data:[^;]+;base64,/, "");
        const buffer = Buffer.from(base64, "base64");
        await fs.writeFile(saveUri.fsPath, buffer);
        void vscode.window.showInformationMessage(`Saved: ${path.basename(saveUri.fsPath)}`);
      } catch (err) {
        void vscode.window.showErrorMessage(`Failed to save file: ${asError(err).message}`);
      }
    });

    await refreshPanel();
  }

  public async refresh(uri: vscode.Uri): Promise<boolean> {
    const session = this.sessions.get(uri.toString());
    if (!session) {
      return false;
    }
    await this.updateWebview(session.document, session.panel);
    return true;
  }

  private async updateWebview(
    document: vscode.TextDocument,
    panel: vscode.WebviewPanel
  ): Promise<void> {
    panel.webview.html = this.renderLoadingHtml(document.uri.fsPath);

    try {
      const { html, bundleHtml } = await this.renderWithAseTs(document);
      panel.webview.html = html;
      // Send the self-contained template bundle via postMessage to avoid
      // embedding a ~900KB string inline (which causes </script> escaping bugs).
      if (bundleHtml) {
        // Small delay to ensure the WebView script has set up its message listener.
        setTimeout(() => {
          void panel.webview.postMessage({ type: "initTemplate", template: bundleHtml });
        }, 100);
      }
    } catch (error) {
      panel.webview.html = this.renderErrorHtml(document.uri.fsPath, error);
    }
  }

  private async renderWithAseTs(document: vscode.TextDocument): Promise<{ html: string; bundleHtml: string }> {
    const config = vscode.workspace.getConfiguration("aseview");
    const style = config.get<string>("defaultStyle", "cartoon");
    const viewer = config.get<string>("defaultViewer", "auto");
    const indexExpr = config.get<string>("readIndex", ":").trim() || ":";
    const formatSetting = config.get<string>("readFormat", "").trim();

    const format = formatSetting || inferFormatFromFileName(document.uri.fsPath);
    const workerPath = vscode.Uri.joinPath(
      this.context.extensionUri,
      "node",
      "parse_with_ase_ts.bundle.cjs"
    ).fsPath;

    const allFrames = await this.parseWithAseTsWorker(
      workerPath,
      document.uri.fsPath,
      format
    );
    if (allFrames.length === 0) {
      throw new Error("No structure frames were parsed from this file.");
    }

    const selectedFrames = selectFramesByIndex(allFrames, indexExpr);
    if (selectedFrames.length === 0) {
      throw new Error(`No frames selected by index expression: ${indexExpr}`);
    }

    // Read the self-contained bundle (three.js + gifshot + download relay all inlined)
    const bundlePath = vscode.Uri.joinPath(
      this.context.extensionUri, "media", "molecular_viewer_bundle.html"
    ).fsPath;
    let bundleHtml = "";
    try {
      bundleHtml = await fs.readFile(bundlePath, "utf8");
    } catch {
      // Bundle not available; fallback to CDN (GIF saving won't work)
    }

    return {
      html: this.renderAseviewHtml(document.uri.fsPath, selectedFrames, viewer, style),
      bundleHtml,
    };
  }

  private parseWithAseTsWorker(
    workerScriptPath: string,
    filePath: string,
    format?: string
  ): Promise<ViewerFrame[]> {
    return new Promise((resolve, reject) => {
      const args = [workerScriptPath, "--file", filePath];
      if (format) {
        args.push("--format", format);
      }

      const child = spawn(process.execPath, args, {
        windowsHide: true,
      });

      let stdout = "";
      let stderr = "";
      const timeoutMs = 45000;
      const timeout = setTimeout(() => {
        child.kill();
        reject(new Error(`ase-ts worker timed out after ${timeoutMs / 1000}s.`));
      }, timeoutMs);

      child.stdout.setEncoding("utf8");
      child.stderr.setEncoding("utf8");

      child.stdout.on("data", (chunk: string) => {
        stdout += chunk;
      });
      child.stderr.on("data", (chunk: string) => {
        stderr += chunk;
      });

      child.on("error", (err) => {
        clearTimeout(timeout);
        reject(new Error(`Failed to start ase-ts worker: ${err.message}`));
      });

      child.on("close", (code) => {
        clearTimeout(timeout);
        if (code !== 0) {
          const details = [stderr.trim(), stdout.trim()].filter(Boolean).join("\n");
          reject(new Error(details || `ase-ts worker failed with exit code ${code}.`));
          return;
        }

        try {
          const parsed = JSON.parse(stdout) as { frames?: ViewerFrame[] };
          const frames = parsed.frames;
          if (!Array.isArray(frames)) {
            reject(new Error("ase-ts worker returned invalid payload."));
            return;
          }
          resolve(frames);
        } catch (err) {
          reject(new Error(`Failed to parse ase-ts worker output: ${asError(err).message}`));
        }
      });
    });
  }

  private renderAseviewHtml(
    filePath: string,
    frames: ViewerFrame[],
    viewerSetting: string,
    style: string
  ): string {
    const viewerType = normalizeViewerType(viewerSetting);
    const initData =
      viewerType === "frag" || viewerType === "normal" ? frames[0] : frames;
    const dataJson = JSON.stringify(initData);
    const optionsJson = JSON.stringify({ style });
    const baseName = path.basename(filePath, path.extname(filePath));

    return `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>ASEView</title>
  <style>
    html, body { margin:0; padding:0; width:100%; height:100%; overflow:hidden; background:#0b1220; color:#e5e7eb; font-family:Arial,sans-serif; }
    #root { width:100%; height:100%; position:relative; }
    .loading { position:absolute; inset:0; display:flex; align-items:center; justify-content:center; background:#0b1220; color:#93c5fd; z-index:5; font-size:13px; }
  </style>
</head>
<body>
  <div id="root"><div class="loading" id="loading">ASEView is loading...</div></div>
  <script src="https://raw.githack.com/kangmg/aseview/main/aseview/static/js/aseview.js"></script>
  <script>
    (function () {
      const vscodeApi = acquireVsCodeApi();
      const data = ${dataJson};
      const options = ${optionsJson};
      const viewerType = ${JSON.stringify(viewerType)};
      const filePath = ${JSON.stringify(filePath)};
      const baseName = ${JSON.stringify(baseName)};

      let _blobUrl = null;

      // Set up iframe.src override BEFORE aseview.js can create the iframe.
      // When molecular_viewer.html is requested, redirect to our self-contained blob.
      const _srcDesc = Object.getOwnPropertyDescriptor(HTMLIFrameElement.prototype, 'src');
      if (_srcDesc && _srcDesc.set) {
        Object.defineProperty(HTMLIFrameElement.prototype, 'src', {
          configurable: true,
          get: function() { return _srcDesc.get.call(this); },
          set: function(url) {
            if (_blobUrl && url && url.includes('molecular_viewer.html')) {
              _srcDesc.set.call(this, _blobUrl);
            } else {
              _srcDesc.set.call(this, url);
            }
          }
        });
      }

      // Receive bundle HTML or save-file relay from iframe
      window.addEventListener('message', function(e) {
        const msg = e.data;
        if (!msg) return;

        // Template bundle from extension host
        if (msg.type === 'initTemplate' && msg.template) {
          const blob = new Blob([msg.template], { type: 'text/html' });
          _blobUrl = URL.createObjectURL(blob);
          _initViewer();
          return;
        }

        // Save-file relay from iframe
        if (msg.type === 'saveFile') {
          vscodeApi.postMessage(msg);
        }
      });

      let _viewerCreated = false;
      function _initViewer() {
        if (_viewerCreated) return;
        _viewerCreated = true;

        const root = document.getElementById('root');
        const loading = document.getElementById('loading');

        if (!window.ASEView) {
          if (loading) loading.textContent = 'Failed to load ASEView module.';
          return;
        }

        let viewer = null;
        if (viewerType === 'overlay') {
          viewer = new window.ASEView.OverlayViewer(root, options);
        } else if (viewerType === 'frag') {
          viewer = new window.ASEView.FragSelector(root, options);
        } else if (viewerType === 'normal') {
          viewer = new window.ASEView.NormalModeViewer(root, options);
        } else {
          viewer = new window.ASEView.MolecularViewer(root, options);
        }

        viewer.setData(data);
        if (loading) loading.remove();
        window.__aseview_meta = { filePath, viewerType };
      }
    })();
  </script>
</body>
</html>`;
  }


  private renderLoadingHtml(filePath: string): string {
    return `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>ASEView Loading</title>
  <style>
    body {
      font-family: Arial, sans-serif;
      background: #0f172a;
      color: #e2e8f0;
      margin: 0;
      height: 100vh;
      display: flex;
      align-items: center;
      justify-content: center;
    }
    .box {
      max-width: 720px;
      width: calc(100% - 32px);
      padding: 20px;
      border: 1px solid #334155;
      border-radius: 10px;
      background: #111827;
    }
    .title { font-size: 16px; margin-bottom: 10px; color: #93c5fd; }
    .path { font-family: Consolas, monospace; font-size: 12px; color: #cbd5e1; }
  </style>
</head>
<body>
  <div class="box">
    <div class="title">ASEView is rendering the structure...</div>
    <div class="path">${escapeHtml(filePath)}</div>
  </div>
</body>
</html>`;
  }

  private renderErrorHtml(filePath: string, error: unknown): string {
    const err = asError(error);
    const message = escapeHtml(err.message);
    const hint = escapeHtml(
      [
        "This extension parses files with ase-ts (no Python runtime required).",
        "If parsing fails, try setting 'aseview.readFormat' explicitly (e.g. vasp, cif, extxyz, vasprun-xml).",
        "You can also adjust 'aseview.readIndex' (for trajectory files).",
      ].join("\n")
    );

    return `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>ASEView Error</title>
  <style>
    body {
      font-family: Arial, sans-serif;
      background: #1f1f1f;
      color: #f3f4f6;
      margin: 0;
      padding: 16px;
    }
    .card {
      max-width: 980px;
      margin: 0 auto;
      background: #111827;
      border: 1px solid #374151;
      border-radius: 10px;
      overflow: hidden;
    }
    .header {
      padding: 12px 14px;
      background: #7f1d1d;
      border-bottom: 1px solid #991b1b;
      font-weight: 700;
    }
    .body { padding: 14px; }
    .path, .block {
      font-family: Consolas, monospace;
      font-size: 12px;
      white-space: pre-wrap;
      word-break: break-word;
      background: #0b1220;
      border: 1px solid #334155;
      border-radius: 6px;
      padding: 10px;
      margin-top: 10px;
    }
    .hint-title { margin-top: 14px; color: #93c5fd; font-size: 13px; }
  </style>
</head>
<body>
  <div class="card">
    <div class="header">ASEView failed to render this file</div>
    <div class="body">
      <div><strong>File</strong></div>
      <div class="path">${escapeHtml(filePath)}</div>
      <div class="hint-title"><strong>Error</strong></div>
      <div class="block">${message}</div>
      <div class="hint-title"><strong>Quick Fix</strong></div>
      <div class="block">${hint}</div>
    </div>
  </div>
</body>
</html>`;
  }
}

function coerceResourceUri(arg: unknown): vscode.Uri | undefined {
  if (arg instanceof vscode.Uri) {
    return arg;
  }
  if (Array.isArray(arg) && arg.length > 0 && arg[0] instanceof vscode.Uri) {
    return arg[0];
  }
  return undefined;
}

function getActiveResourceUri(): vscode.Uri | undefined {
  const fromEditor = vscode.window.activeTextEditor?.document.uri;
  if (fromEditor) {
    return fromEditor;
  }

  const activeTab = vscode.window.tabGroups.activeTabGroup.activeTab;
  if (!activeTab) {
    return undefined;
  }

  const input = activeTab.input;
  if (input instanceof vscode.TabInputText) {
    return input.uri;
  }
  if (input instanceof vscode.TabInputCustom) {
    return input.uri;
  }
  if (input instanceof vscode.TabInputNotebook) {
    return input.uri;
  }

  return undefined;
}

function escapeHtml(value: string): string {
  return value
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&#39;");
}

function asError(error: unknown): Error {
  if (error instanceof Error) {
    return error;
  }
  return new Error(String(error));
}

function inferFormatFromFileName(filePath: string): string | undefined {
  const base = path.basename(filePath).toLowerCase();
  if (base === "poscar") {
    return "vasp";
  }
  if (base === "contcar") {
    return "vasp";
  }
  if (base === "vasprun.xml") {
    return "vasprun-xml";
  }

  const ext = path.extname(filePath).toLowerCase();
  if (!ext) {
    return undefined;
  }
  if (ext === ".xml") {
    return "vasprun-xml";
  }
  return ext.slice(1);
}

function selectFramesByIndex(frames: ViewerFrame[], indexExpr: string): ViewerFrame[] {
  const expr = indexExpr.trim();
  if (!expr || expr === ":" || expr === "all") {
    return frames;
  }

  if (!expr.includes(":")) {
    const index = toSignedIndex(parseInteger(expr, "index"), frames.length);
    return [frames[index]];
  }

  const [startRaw, stopRaw, stepRaw] = expr.split(":");
  const step = stepRaw === undefined || stepRaw === "" ? 1 : parseInteger(stepRaw, "step");
  if (step === 0) {
    throw new Error("Index slice step cannot be zero.");
  }

  const startDefault = step > 0 ? 0 : frames.length - 1;
  const stopDefault = step > 0 ? frames.length : -1;
  let start = startRaw === "" || startRaw === undefined ? startDefault : parseInteger(startRaw, "start");
  let stop = stopRaw === "" || stopRaw === undefined ? stopDefault : parseInteger(stopRaw, "stop");

  start = normalizeSliceBound(start, frames.length);
  stop = normalizeSliceBound(stop, frames.length);

  const selected: ViewerFrame[] = [];
  if (step > 0) {
    for (let i = start; i < stop; i += step) {
      if (i >= 0 && i < frames.length) {
        selected.push(frames[i]);
      }
    }
  } else {
    for (let i = start; i > stop; i += step) {
      if (i >= 0 && i < frames.length) {
        selected.push(frames[i]);
      }
    }
  }
  return selected;
}

function parseInteger(input: string, label: string): number {
  const value = Number(input);
  if (!Number.isInteger(value)) {
    throw new Error(`Invalid ${label} in index expression: ${input}`);
  }
  return value;
}

function toSignedIndex(index: number, length: number): number {
  const resolved = index < 0 ? length + index : index;
  if (resolved < 0 || resolved >= length) {
    throw new Error(`Index out of range: ${index}`);
  }
  return resolved;
}

function normalizeSliceBound(value: number, length: number): number {
  if (value < 0) {
    return Math.max(-1, length + value);
  }
  return Math.min(length, value);
}

function normalizeViewerType(input: string): "molecular" | "overlay" | "frag" | "normal" {
  if (input === "overlay" || input === "frag" || input === "normal" || input === "molecular") {
    return input;
  }
  return "molecular";
}

export function activate(context: vscode.ExtensionContext): void {
  context.subscriptions.push(AseviewEditorProvider.register(context));
}

export function deactivate(): void {}
