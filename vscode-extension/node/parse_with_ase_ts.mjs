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

function toPerAtomScalarArray(value, nAtoms, label, { allowVectors = false } = {}) {
  if (!Array.isArray(value)) {
    throw new Error(`Expected ${label} to be an array.`);
  }
  if (value.length !== nAtoms) {
    throw new Error(`Expected ${label} to contain ${nAtoms} per-atom values.`);
  }
  return value.map((item, index) => {
    if (allowVectors && Array.isArray(item)) {
      if (item.length < 3) {
        throw new Error(`Expected ${label}[${index}] to contain 3 values.`);
      }
      const xyz = item.slice(0, 3).map((x) => Number(x));
      if (xyz.some((x) => !Number.isFinite(x))) {
        throw new Error(`Expected ${label}[${index}] vector values to be numeric.`);
      }
      return Math.hypot(...xyz);
    }
    const n = Number(item);
    if (!Number.isFinite(n)) {
      throw new Error(`Expected ${label}[${index}] to be numeric.`);
    }
    return n;
  });
}

function toBooleanValue(value, label) {
  if (typeof value === "boolean") {
    return value;
  }
  if (typeof value === "number") {
    return value !== 0;
  }
  if (typeof value === "string") {
    const normalized = value.trim().toLowerCase();
    if (["t", "true", "1"].includes(normalized)) {
      return true;
    }
    if (["f", "false", "0"].includes(normalized)) {
      return false;
    }
  }
  throw new Error(`Expected ${label} to contain boolean-like values.`);
}

function toBooleanArray(value, label) {
  if (!Array.isArray(value)) {
    throw new Error(`Expected ${label} to be an array.`);
  }
  return value.map((item, index) => toBooleanValue(item, `${label}[${index}]`));
}

function toMoveMaskArray(value, label) {
  if (!Array.isArray(value)) {
    throw new Error(`Expected ${label} to be an array.`);
  }
  return value.map((item, index) => {
    if (Array.isArray(item)) {
      if (item.length < 3) {
        throw new Error(`Expected ${label}[${index}] to contain 3 values.`);
      }
      return item.slice(0, 3).map((x, axis) => toBooleanValue(x, `${label}[${index}][${axis}]`));
    }
    return toBooleanValue(item, `${label}[${index}]`);
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
      frame.pbc = pbc.map((value) => Boolean(value));
    }
  }

  const has = typeof atoms.has === "function" ? atoms.has.bind(atoms) : () => false;
  const getArray =
    typeof atoms.getArray === "function" ? atoms.getArray.bind(atoms) : () => undefined;

  let forcesRaw;
  if (typeof atoms.getForces === "function") {
    try {
      forcesRaw = atoms.getForces();
    } catch {
      // Some formats do not carry forces; fall back to arrays below.
    }
  }
  if (!forcesRaw && has("forces")) {
    forcesRaw = getArray("forces", true);
  }
  if (forcesRaw) {
    frame.forces = toVec3Array(forcesRaw, "forces");
  }

  const chargesRaw = has("charges")
    ? getArray("charges", true)
    : atoms.getInitialCharges?.();
  if (chargesRaw) {
    frame.charges = toPerAtomScalarArray(chargesRaw, symbols.length, "charges");
  }

  const magmomsRaw = has("magmoms")
    ? getArray("magmoms", true)
    : has("magmom")
      ? getArray("magmom", true)
      : has("initial_magmoms")
        ? getArray("initial_magmoms", true)
        : atoms.getInitialMagneticMoments?.();
  if (magmomsRaw) {
    frame.magmoms = toPerAtomScalarArray(magmomsRaw, symbols.length, "magmoms", {
      allowVectors: true,
    });
  }

  const fixedRaw = has("fixed") ? getArray("fixed", true) : undefined;
  if (fixedRaw) {
    const fixedMask = toBooleanArray(fixedRaw, "fixed");
    const fixed = [];
    fixedMask.forEach((value, index) => {
      if (value) {
        fixed.push(index);
      }
    });
    if (fixed.length > 0) {
      frame.fixed = fixed;
      frame.arrays = { ...(frame.arrays ?? {}), fixed: fixedMask };
    }
  }

  const moveMaskRaw = has("move_mask") ? getArray("move_mask", true) : undefined;
  if (moveMaskRaw) {
    const moveMask = toMoveMaskArray(moveMaskRaw, "move_mask");
    frame.arrays = { ...(frame.arrays ?? {}), move_mask: moveMask };
    const fixed = [];
    moveMask.forEach((value, index) => {
      const mask = Array.isArray(value) ? value : [value, value, value];
      if (mask.every((axis) => axis === false)) {
        fixed.push(index);
      }
    });
    if (fixed.length > 0 && !frame.fixed) {
      frame.fixed = fixed;
    }
  }

  const info = atoms.info ?? {};
  let energy;
  if (typeof atoms.getPotentialEnergy === "function") {
    try {
      energy = atoms.getPotentialEnergy();
    } catch {
      // Some formats do not carry energy; fall back to info below.
    }
  }
  if (typeof energy !== "number" && atoms.calc?.results && typeof atoms.calc.results.energy === "number") {
    energy = atoms.calc.results.energy;
  }
  if (typeof energy !== "number" && typeof info.energy === "number") {
    energy = info.energy;
  }
  if (typeof energy === "number" && Number.isFinite(energy)) {
    frame.energy = energy;
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
