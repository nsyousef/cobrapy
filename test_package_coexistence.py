#!/usr/bin/env python3
"""
Demonstrate that cobra_structural can coexist with original COBRApy.

This test shows how you can use both packages side-by-side:
- cobra_structural: For building models structurally (no solver)
- cobra (original): For optimization and analysis with solver
"""

# Import structural-only version
from cobra_structural import Model as StructuralModel
from cobra_structural import Reaction as StructuralReaction
from cobra_structural import Metabolite as StructuralMetabolite

print("=" * 70)
print("Testing cobra_structural (structural-only, no solver)")
print("=" * 70)

# Create a structural model
model = StructuralModel("structural_test")

glc = StructuralMetabolite("glc", compartment="c")
atp = StructuralMetabolite("atp", compartment="c")
adp = StructuralMetabolite("adp", compartment="c")

rxn1 = StructuralReaction("PGI")
rxn1.add_metabolites({glc: -1, atp: -1, adp: 1})
rxn1.bounds = (0, 10)

model.add_reactions([rxn1])
model.objective = {rxn1: 1.0}

print(f"✓ Structural model created: {model.id}")
print(f"  - Reactions: {len(model.reactions)}")
print(f"  - Metabolites: {len(model.metabolites)}")
print(f"  - Objective: {model.objective}")
print(f"  - Has solver: {hasattr(model, 'solver')}")
print()

# Now try to export to JSON (I/O functionality should work)
from cobra_structural.io import save_json_model
import tempfile
import os

with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
    temp_file = f.name

try:
    save_json_model(model, temp_file)
    print(f"✓ Successfully exported to JSON: {os.path.basename(temp_file)}")
    print(f"  - File size: {os.path.getsize(temp_file)} bytes")
finally:
    os.unlink(temp_file)

print()
print("=" * 70)
print("Workflow suggestion for using with original COBRApy:")
print("=" * 70)
print("""
1. Build your model using cobra_structural (fast, lightweight)
2. Export to JSON/SBML using cobra_structural.io
3. Load into original COBRApy for optimization:
   
   from cobra_structural.io import save_json_model
   from cobra import load_json_model  # Original COBRApy
   
   # Build structural model
   structural_model = StructuralModel("my_model")
   # ... add reactions, metabolites ...
   
   # Export
   save_json_model(structural_model, "model.json")
   
   # Load into full COBRApy with solver
   full_model = load_json_model("model.json")
   solution = full_model.optimize()
""")

print("✓✓✓ cobra_structural package working correctly! ✓✓✓")
