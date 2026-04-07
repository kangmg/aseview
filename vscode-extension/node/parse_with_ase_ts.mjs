#!/usr/bin/env node

import fs from "node:fs";
import { readAll } from "ase-ts";

function parseArgs(argv) {
  const out = {};
  for (let i = 0; i < argv.length; i += 1) {
    const token = argv[i];
    if (token === "--file" && i + 1 < argv.length) {
      out.file = argv[i + 1];
      i += 1;
      continue;
    }
    if (token === "--format" && i + 1 < argv.length) {
      out.format = argv[i + 1];
      i += 1;
      continue;
    }
  }
  return out;
}

function readAllWithFallback(text, format) {
  try {
    return readAll(text, {
      data: true,
      ...(format ? { format } : {}),
    });
  } catch (err) {
    const message = String(err && err.message ? err.message : err).toLowerCase();
    const normalized = (format ?? "").toLowerCase();
    const extxyzLike = normalized === "extxyz" || normalized === "xyz";
    const tagsConflict = message.includes("array 'tags' already exists");
    if (extxyzLike && tagsConflict) {
      return readAll(text, { data: true, format: "xyz" });
    }
    throw err;
  }
}

function toNumberArray(value, label) {
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

function toVec3Array(value, label) {
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
    return xyz;
  });
}

function atomsToViewerFrame(atoms) {
  const symbols = (atoms.getChemicalSymbols?.() ?? []).map((x) => String(x));
  if (!Array.isArray(symbols) || symbols.length === 0) {
    throw new Error("Parsed frame has no symbols.");
  }

  const positions = toVec3Array(atoms.getPositions?.() ?? [], "positions");
  const frame = { symbols, positions };

  const pbc = atoms.getPbc?.() ?? [false, false, false];
  if (Array.isArray(pbc) && pbc.some(Boolean)) {
    const cell = atoms.getCell?.();
    if (cell) {
      frame.cell = toVec3Array(cell, "cell");
    }
  }

  const has = typeof atoms.has === "function" ? atoms.has.bind(atoms) : () => false;
  const getArray =
    typeof atoms.getArray === "function" ? atoms.getArray.bind(atoms) : () => undefined;

  const forcesRaw = has("forces") ? getArray("forces", true) : undefined;
  if (forcesRaw) {
    frame.forces = toVec3Array(forcesRaw, "forces");
  }

  const chargesRaw = has("charges")
    ? getArray("charges", true)
    : atoms.getInitialCharges?.();
  if (chargesRaw) {
    frame.charges = toNumberArray(chargesRaw, "charges");
  }

  const info = atoms.info ?? {};
  if (typeof info.energy === "number" && Number.isFinite(info.energy)) {
    frame.energy = info.energy;
  }
  if (typeof info.name === "string" && info.name.trim()) {
    frame.name = info.name;
  }

  return frame;
}

function main() {
  const args = parseArgs(process.argv.slice(2));
  if (!args.file) {
    throw new Error("Missing required argument: --file");
  }

  const text = fs.readFileSync(args.file, "utf8");
  const frames = readAllWithFallback(text, args.format).map((atoms) =>
    atomsToViewerFrame(atoms)
  );

  process.stdout.write(JSON.stringify({ frames }));
}

try {
  main();
} catch (err) {
  const msg = err && err.message ? err.message : String(err);
  process.stderr.write(`${msg}\n`);
  process.exit(1);
}
