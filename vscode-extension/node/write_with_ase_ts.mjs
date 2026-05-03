#!/usr/bin/env node

import { Atoms } from "ase-ts";
import { writeAll } from "ase-ts/browser";

function parseArgs(argv) {
  const out = {};
  for (let i = 0; i < argv.length; i += 1) {
    const token = argv[i];
    if (token === "--format" && i + 1 < argv.length) {
      out.format = argv[i + 1];
      i += 1;
    }
  }
  return out;
}

function readStdin() {
  return new Promise((resolve, reject) => {
    let text = "";
    process.stdin.setEncoding("utf8");
    process.stdin.on("data", (chunk) => {
      text += chunk;
    });
    process.stdin.on("error", reject);
    process.stdin.on("end", () => resolve(text));
  });
}

function toVec3Array(value, label) {
  if (!Array.isArray(value)) {
    throw new Error(`Expected ${label} to be an array.`);
  }
  return value.map((vec, index) => {
    if (!Array.isArray(vec) || vec.length < 3) {
      throw new Error(`Expected ${label}[${index}] to be a 3D vector.`);
    }
    const xyz = vec.slice(0, 3).map((x) => Number(x));
    if (xyz.some((x) => !Number.isFinite(x))) {
      throw new Error(`Expected ${label}[${index}] to contain numeric values.`);
    }
    return xyz;
  });
}

function normalizeBool(value) {
  if (typeof value === "boolean") return value;
  if (typeof value === "number") return value !== 0;
  if (typeof value === "string") {
    const normalized = value.trim().toLowerCase();
    if (["t", "true", "1"].includes(normalized)) return true;
    if (["f", "false", "0"].includes(normalized)) return false;
  }
  return false;
}

function fixedIndicesFromFrame(frame, natoms) {
  const fixed = new Set();

  if (Array.isArray(frame.fixed)) {
    for (const value of frame.fixed) {
      const index = Number(value);
      if (Number.isInteger(index) && index >= 0 && index < natoms) {
        fixed.add(index);
      }
    }
  }

  if (frame.arrays && Array.isArray(frame.arrays.fixed)) {
    frame.arrays.fixed.forEach((value, index) => {
      if (index < natoms && normalizeBool(value)) {
        fixed.add(index);
      }
    });
  }

  if (frame.arrays && Array.isArray(frame.arrays.move_mask)) {
    frame.arrays.move_mask.forEach((value, index) => {
      if (index >= natoms) return;
      const mask = Array.isArray(value)
        ? value.slice(0, 3).map(normalizeBool)
        : [normalizeBool(value), normalizeBool(value), normalizeBool(value)];
      if (mask.length >= 3 && mask.every((axis) => axis === false)) {
        fixed.add(index);
      }
    });
  }

  return [...fixed].sort((a, b) => a - b);
}

function frameToAtoms(frame) {
  if (!frame || typeof frame !== "object") {
    throw new Error("Expected frame to be an object.");
  }

  const symbols = Array.isArray(frame.symbols) ? frame.symbols.map((x) => String(x)) : [];
  const positions = toVec3Array(frame.positions, "positions");
  if (symbols.length === 0 || symbols.length !== positions.length) {
    throw new Error("Frame symbols and positions must have the same non-zero length.");
  }

  const init = { symbols, positions };

  if (frame.cell) {
    init.cell = toVec3Array(frame.cell, "cell");
  }
  if (Array.isArray(frame.pbc)) {
    init.pbc = frame.pbc.slice(0, 3).map(Boolean);
  } else if (frame.cell) {
    init.pbc = [true, true, true];
  }

  const fixedIndices = fixedIndicesFromFrame(frame, symbols.length);
  if (fixedIndices.length > 0) {
    init.fixedIndices = fixedIndices;
  }

  if (typeof frame.name === "string" && frame.name.trim()) {
    init.info = { comment: frame.name.trim() };
  }

  return new Atoms(init);
}

function normalizeWriteFormat(format) {
  const normalized = String(format || "").trim().toLowerCase();
  if (normalized === "vasp-poscar" || normalized === "poscar") return "vasp";
  if (["xyz", "extxyz", "cif", "vasp"].includes(normalized)) return normalized;
  throw new Error(`Unsupported ase-ts write format: ${format}`);
}

function writeFrames(frames, format) {
  const atoms = frames.map(frameToAtoms);
  if (atoms.length === 0) {
    throw new Error("No frames to write.");
  }

  if (format === "vasp") {
    return atoms.map((single) => writeAll([single], format, { vasp5: true })).join("\n");
  }

  const output = writeAll(atoms, format);
  if (typeof output !== "string") {
    throw new Error(`Expected text output for format: ${format}`);
  }
  return output;
}

async function main() {
  const args = parseArgs(process.argv.slice(2));
  const format = normalizeWriteFormat(args.format);
  const input = JSON.parse(await readStdin());
  if (!Array.isArray(input.frames)) {
    throw new Error("Expected JSON payload with frames array.");
  }
  process.stdout.write(writeFrames(input.frames, format));
}

main().catch((err) => {
  const msg = err && err.message ? err.message : String(err);
  process.stderr.write(`${msg}\n`);
  process.exit(1);
});
