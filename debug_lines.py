#!/usr/bin/env python
"""Debug script to check for line artifacts on spheres"""

from ase.build import molecule
from aseview import MolecularViewer

water = molecule('H2O')

print("Creating test files with different settings...\n")

# Test 1: Default glossy
v1 = MolecularViewer(water, style='glossy')
v1.save_html('debug_glossy_default.html')
print("✓ debug_glossy_default.html")

# Test 2: Cartoon style
v2 = MolecularViewer(water, style='cartoon')
v2.save_html('debug_cartoon.html')
print("✓ debug_cartoon.html")

# Test 3: Default style
v3 = MolecularViewer(water, style='default')
v3.save_html('debug_default.html')
print("✓ debug_default.html")

# Test 4: Metallic style
v4 = MolecularViewer(water, style='metallic')
v4.save_html('debug_metallic.html')
print("✓ debug_metallic.html")

print("\nOpen these files and check:")
print("1. Are there visible lines on the spheres?")
print("2. Which style shows the lines most clearly?")
print("3. Do the lines move when you rotate?")
print("\nIf lines are visible:")
print("- They might be from low polygon count (but we have 64 segments)")
print("- They might be from lighting angle")
print("- They might be from material shininess/specular")
print("- They might be from shadow acne")
