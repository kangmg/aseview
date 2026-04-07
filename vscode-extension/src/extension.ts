import * as path from "path";
import * as os from "os";
import * as fs from "fs/promises";
import { spawn } from "child_process";
import * as vscode from "vscode";

const VIEW_TYPE = "aseview.viewer";

type CommandCandidate = {
  command: string;
  prefixArgs: string[];
  label: string;
};

class RendererLaunchError extends Error {
  public readonly notFound: boolean;

  constructor(message: string, notFound = false) {
    super(message);
    this.name = "RendererLaunchError";
    this.notFound = notFound;
  }
}

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
      localResourceRoots: [this.context.extensionUri],
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
      const renderedHtml = await this.renderWithAseTs(document);
      panel.webview.html = renderedHtml;
    } catch (error) {
      panel.webview.html = this.renderErrorHtml(document.uri.fsPath, error);
    }
  }

  private async renderWithAseTs(document: vscode.TextDocument): Promise<string> {
    const config = vscode.workspace.getConfiguration("aseview");
    const style = config.get<string>("defaultStyle", "cartoon");
    const viewer = config.get<string>("defaultViewer", "auto");
    const indexExpr = config.get<string>("readIndex", ":").trim() || ":";
    const formatSetting = config.get<string>("readFormat", "").trim();

    const format = formatSetting || inferFormatFromFileName(document.uri.fsPath);
    const contentBytes = await vscode.workspace.fs.readFile(document.uri);
    const text = new TextDecoder("utf-8").decode(contentBytes);

    const aseTs = await importAseTs();
    const allFramesRaw = aseTs.readAll(text, {
      data: true,
      ...(format ? { format } : {}),
    });

    const allFrames = allFramesRaw.map((atoms) => atomsToViewerFrame(atoms));
    if (allFrames.length === 0) {
      throw new Error("No structure frames were parsed from this file.");
    }

    const selectedFrames = selectFramesByIndex(allFrames, indexExpr);
    if (selectedFrames.length === 0) {
      throw new Error(`No frames selected by index expression: ${indexExpr}`);
    }

    return this.renderAseviewHtml(document.uri.fsPath, selectedFrames, viewer, style);
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

    return `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>ASEView</title>
  <style>
    html, body {
      margin: 0;
      padding: 0;
      width: 100%;
      height: 100%;
      overflow: hidden;
      background: #0b1220;
      color: #e5e7eb;
      font-family: Arial, sans-serif;
    }
    #root {
      width: 100%;
      height: 100%;
      position: relative;
    }
    .loading {
      position: absolute;
      inset: 0;
      display: flex;
      align-items: center;
      justify-content: center;
      background: #0b1220;
      color: #93c5fd;
      z-index: 5;
      font-size: 13px;
    }
  </style>
</head>
<body>
  <div id="root">
    <div class="loading" id="loading">ASEView is loading...</div>
  </div>
  <script src="https://raw.githack.com/kangmg/aseview/main/aseview/static/js/aseview.js"></script>
  <script>
    (function () {
      const data = ${dataJson};
      const options = ${optionsJson};
      const viewerType = ${JSON.stringify(viewerType)};
      const filePath = ${JSON.stringify(filePath)};

      const root = document.getElementById("root");
      const loading = document.getElementById("loading");

      if (!window.ASEView) {
        if (loading) loading.textContent = "Failed to load ASEView JavaScript module.";
        return;
      }

      let viewer = null;
      if (viewerType === "overlay") {
        viewer = new window.ASEView.OverlayViewer(root, options);
      } else if (viewerType === "frag") {
        viewer = new window.ASEView.FragSelector(root, options);
      } else if (viewerType === "normal") {
        viewer = new window.ASEView.NormalModeViewer(root, options);
      } else {
        viewer = new window.ASEView.MolecularViewer(root, options);
      }

      viewer.setData(data);
      if (loading) loading.remove();
      window.__aseview_meta = { filePath, viewerType };
    })();
  </script>
</body>
</html>`;
  }

  private async renderWithPython(document: vscode.TextDocument): Promise<string> {
    const config = vscode.workspace.getConfiguration("aseview");
    const style = config.get<string>("defaultStyle", "cartoon");
    const viewer = config.get<string>("defaultViewer", "auto");
    const index = config.get<string>("readIndex", ":");
    const format = config.get<string>("readFormat", "").trim();

    const scriptPath = vscode.Uri.joinPath(
      this.context.extensionUri,
      "python",
      "render_viewer.py"
    ).fsPath;

    const workspaceFolder = vscode.workspace.getWorkspaceFolder(document.uri);
    const cwd = workspaceFolder?.uri.fsPath ?? path.dirname(document.uri.fsPath);
    const baseArgs = [
      scriptPath,
      "--file",
      document.uri.fsPath,
      "--viewer",
      viewer,
      "--style",
      style,
      "--index",
      index,
    ];
    if (format) {
      baseArgs.push("--format", format);
    }

    const candidates = this.getPythonCandidates();
    let lastError: Error | undefined;

    for (const candidate of candidates) {
      const args = [...candidate.prefixArgs, ...baseArgs];
      try {
        const output = await this.runRenderer(candidate, args, cwd);
        if (!output.trim()) {
          throw new RendererLaunchError(
            `Renderer returned empty output (${candidate.label}).`
          );
        }
        return output;
      } catch (error) {
        const err = asError(error);
        lastError = err;
        if (err instanceof RendererLaunchError && err.notFound) {
          continue;
        }
      }
    }

    const pythonError = lastError ?? new Error("Failed to run ASEView renderer with Python.");

    try {
      return await this.renderWithAseviewCli(document, style, viewer, index, format, cwd);
    } catch (cliError) {
      const cliErr = asError(cliError);
      throw new Error(
        [
          "Python renderer failed.",
          pythonError.message,
          "",
          "ASEView CLI fallback failed.",
          cliErr.message,
        ].join("\n")
      );
    }
  }

  private getPythonCandidates(): CommandCandidate[] {
    const configured = vscode.workspace
      .getConfiguration("aseview")
      .get<string>("pythonPath", "")
      .trim();

    const rawCandidates: CommandCandidate[] = [];
    if (configured) {
      rawCandidates.push({
        command: configured,
        prefixArgs: [],
        label: configured,
      });
    }

    if (process.platform === "win32") {
      rawCandidates.push(
        { command: "python", prefixArgs: [], label: "python" },
        { command: "py", prefixArgs: ["-3"], label: "py -3" }
      );
    } else {
      rawCandidates.push(
        { command: "python3", prefixArgs: [], label: "python3" },
        { command: "python", prefixArgs: [], label: "python" }
      );
    }

    const dedup = new Map<string, CommandCandidate>();
    for (const candidate of rawCandidates) {
      const key = `${candidate.command}|${candidate.prefixArgs.join(" ")}`;
      if (!dedup.has(key)) {
        dedup.set(key, candidate);
      }
    }
    return [...dedup.values()];
  }

  private getAseviewCandidates(): CommandCandidate[] {
    const launcher = vscode.workspace
      .getConfiguration("aseview")
      .get<string>("launcher", "")
      .trim();
    const configured = vscode.workspace
      .getConfiguration("aseview")
      .get<string>("commandPath", "")
      .trim();

    const rawCandidates: CommandCandidate[] = [];
    if (launcher) {
      const parsed = parseCommandLine(launcher);
      if (parsed) {
        rawCandidates.push({
          command: parsed.command,
          prefixArgs: parsed.args,
          label: launcher,
        });
      }
    }

    if (configured) {
      rawCandidates.push({
        command: configured,
        prefixArgs: [],
        label: configured,
      });
    }

    rawCandidates.push({ command: "aseview", prefixArgs: [], label: "aseview" });
    rawCandidates.push({
      command: "uv",
      prefixArgs: ["run", "aseview"],
      label: "uv run aseview",
    });

    const dedup = new Map<string, CommandCandidate>();
    for (const candidate of rawCandidates) {
      const key = `${candidate.command}|${candidate.prefixArgs.join(" ")}`;
      if (!dedup.has(key)) {
        dedup.set(key, candidate);
      }
    }
    return [...dedup.values()];
  }

  private async renderWithAseviewCli(
    document: vscode.TextDocument,
    style: string,
    viewer: string,
    index: string,
    format: string,
    cwd: string
  ): Promise<string> {
    const candidates = this.getAseviewCandidates();
    let lastError: Error | undefined;

    for (const candidate of candidates) {
      const tempPath = path.join(
        os.tmpdir(),
        `aseview-vscode-${Date.now()}-${Math.random().toString(36).slice(2)}.html`
      );
      const args = [
        ...candidate.prefixArgs,
        document.uri.fsPath,
        "--no-browser",
        "-o",
        tempPath,
        "--style",
        style,
        "-i",
        index,
      ];

      if (viewer !== "auto") {
        args.push("-v", viewer);
      }
      if (format) {
        args.push("-f", format);
      }

      try {
        await this.runRenderer(candidate, args, cwd);
        const html = await fs.readFile(tempPath, "utf8");
        if (!html.trim()) {
          throw new RendererLaunchError(
            `ASEView CLI returned empty output (${candidate.label}).`
          );
        }
        return html;
      } catch (error) {
        const err = asError(error);
        lastError = err;
        if (err instanceof RendererLaunchError) {
          if (err.notFound) {
            continue;
          }
          // For commands that exist but fail, surface the first actionable error.
          throw err;
        }
        throw err;
      } finally {
        try {
          await fs.unlink(tempPath);
        } catch {
          // Ignore temp cleanup failures.
        }
      }
    }

    throw lastError ?? new Error("Failed to run ASEView CLI fallback.");
  }

  private runRenderer(
    candidate: CommandCandidate,
    args: string[],
    cwd: string
  ): Promise<string> {
    return new Promise((resolve, reject) => {
      const env: NodeJS.ProcessEnv = { ...process.env };
      if (isUvCommand(candidate.command)) {
        const uvRoot = inferUvRoot(cwd, args);
        if (!env.UV_CACHE_DIR) {
          env.UV_CACHE_DIR = path.join(uvRoot, ".uv-cache");
        }
        if (!env.UV_PYTHON_INSTALL_DIR) {
          env.UV_PYTHON_INSTALL_DIR = path.join(uvRoot, ".uv-python");
        }
      }

      const child = spawn(candidate.command, args, {
        cwd,
        windowsHide: true,
        env,
      });

      let stdout = "";
      let stderr = "";
      const timeoutMs = 45000;

      const timeout = setTimeout(() => {
        child.kill();
        reject(
          new RendererLaunchError(
            `Renderer timed out after ${timeoutMs / 1000}s (${candidate.label}).`
          )
        );
      }, timeoutMs);

      child.stdout.setEncoding("utf8");
      child.stderr.setEncoding("utf8");

      child.stdout.on("data", (chunk: string) => {
        stdout += chunk;
      });

      child.stderr.on("data", (chunk: string) => {
        stderr += chunk;
      });

      child.on("error", (err: NodeJS.ErrnoException) => {
        clearTimeout(timeout);
        if (err.code === "ENOENT") {
          reject(new RendererLaunchError(`Command not found: ${candidate.label}`, true));
          return;
        }
        reject(
          new RendererLaunchError(
            `Failed to launch renderer (${candidate.label}): ${err.message}`
          )
        );
      });

      child.on("close", (code, signal) => {
        clearTimeout(timeout);
        if (code === 0) {
          resolve(stdout);
          return;
        }

        const details = [stderr.trim(), stdout.trim()].filter(Boolean).join("\n");
        const signalText = signal ? `, signal ${signal}` : "";
        reject(
          new RendererLaunchError(
            `Renderer failed (${candidate.label}, exit ${code ?? "null"}${signalText}).\n${details}`
          )
        );
      });
    });
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
        "Check Python environment and install dependencies:",
        "pip install aseview ase",
        "If needed, set 'aseview.pythonPath' in VS Code settings.",
        "If needed, set 'aseview.launcher' (example: uv run --project ... aseview).",
        "If ASEView CLI is installed, set 'aseview.commandPath' (or keep default 'aseview').",
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

async function importAseTs(): Promise<{
  readAll: (source: string, options?: { data?: boolean; format?: string }) => unknown[];
}> {
  const dynamicImport = new Function(
    "m",
    "return import(m);"
  ) as (moduleName: string) => Promise<{
    readAll: (source: string, options?: { data?: boolean; format?: string }) => unknown[];
  }>;
  return dynamicImport("ase-ts");
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

function atomsToViewerFrame(atoms: unknown): ViewerFrame {
  const a = atoms as {
    getChemicalSymbols?: () => string[];
    getPositions?: () => number[][];
    getCell?: () => number[][];
    getPbc?: () => [boolean, boolean, boolean];
    has?: (name: string) => boolean;
    getArray?: (name: string, copy?: boolean) => unknown;
    getInitialCharges?: () => number[];
    info?: Record<string, unknown>;
  };

  const symbols = ensureStringArray(a.getChemicalSymbols?.(), "symbols");
  const positions = ensureVec3Array(a.getPositions?.(), "positions");
  const frame: ViewerFrame = { symbols, positions };

  const pbc = a.getPbc?.() ?? [false, false, false];
  if (Array.isArray(pbc) && pbc.some(Boolean)) {
    const cell = a.getCell?.();
    if (cell) {
      frame.cell = ensureVec3Array(cell, "cell");
    }
  }

  const forcesRaw = a.has?.("forces") ? a.getArray?.("forces", true) : undefined;
  if (forcesRaw) {
    frame.forces = ensureVec3Array(forcesRaw, "forces");
  }

  const chargesRaw = a.has?.("charges") ? a.getArray?.("charges", true) : a.getInitialCharges?.();
  if (chargesRaw) {
    frame.charges = ensureNumberArray(chargesRaw, "charges");
  }

  const info = a.info ?? {};
  const energy = info["energy"];
  if (typeof energy === "number" && Number.isFinite(energy)) {
    frame.energy = energy;
  }
  const name = info["name"];
  if (typeof name === "string" && name.trim()) {
    frame.name = name;
  }

  return frame;
}

function ensureStringArray(value: unknown, label: string): string[] {
  if (!Array.isArray(value)) {
    throw new Error(`Expected ${label} to be an array.`);
  }
  return value.map((item) => String(item));
}

function ensureNumberArray(value: unknown, label: string): number[] {
  if (!Array.isArray(value)) {
    throw new Error(`Expected ${label} to be an array.`);
  }
  return value.map((item) => {
    const n = Number(item);
    if (!Number.isFinite(n)) {
      throw new Error(`Expected ${label} to contain numeric values.`);
    }
    return n;
  });
}

function ensureVec3Array(value: unknown, label: string): number[][] {
  if (!Array.isArray(value)) {
    throw new Error(`Expected ${label} to be an array.`);
  }
  return value.map((vec) => {
    if (!Array.isArray(vec) || vec.length < 3) {
      throw new Error(`Expected ${label} to contain 3D vectors.`);
    }
    const xyz = vec.slice(0, 3).map((x) => Number(x));
    if (xyz.some((x) => !Number.isFinite(x))) {
      throw new Error(`Expected ${label} vectors to be numeric.`);
    }
    return xyz as number[];
  });
}

function isUvCommand(command: string): boolean {
  const name = path.basename(command).toLowerCase();
  return name === "uv" || name === "uv.exe";
}

function inferUvRoot(cwd: string, args: string[]): string {
  const projectFlagIndex = args.findIndex((arg) => arg === "--project");
  if (projectFlagIndex >= 0 && projectFlagIndex + 1 < args.length) {
    const projectArg = args[projectFlagIndex + 1];
    if (projectArg.trim()) {
      return path.isAbsolute(projectArg)
        ? projectArg
        : path.resolve(cwd, projectArg);
    }
  }
  return cwd;
}

function parseCommandLine(input: string): { command: string; args: string[] } | undefined {
  const parts: string[] = [];
  let current = "";
  let quote: "'" | '"' | null = null;
  let i = 0;

  while (i < input.length) {
    const ch = input[i];
    if (quote) {
      if (ch === quote) {
        quote = null;
      } else {
        current += ch;
      }
      i += 1;
      continue;
    }

    if (ch === "'" || ch === '"') {
      quote = ch;
      i += 1;
      continue;
    }

    if (/\s/.test(ch)) {
      if (current) {
        parts.push(current);
        current = "";
      }
      i += 1;
      continue;
    }

    current += ch;
    i += 1;
  }

  if (current) {
    parts.push(current);
  }

  if (parts.length === 0) {
    return undefined;
  }

  return {
    command: parts[0],
    args: parts.slice(1),
  };
}

export function activate(context: vscode.ExtensionContext): void {
  context.subscriptions.push(AseviewEditorProvider.register(context));
}

export function deactivate(): void {}
