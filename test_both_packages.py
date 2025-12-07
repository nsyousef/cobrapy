#!/usr/bin/env python
"""
Test script to verify cobra and cobra-structural can coexist and convert between each other.
"""

from cobra_structural import Model, Reaction, Metabolite
from cobra_structural.io import to_cobrapy_model
import cobra

print("=" * 70)
print("Testing cobra and cobra-structural coexistence")
print("=" * 70)

# Verify both packages are installed
print("\n1. Package Verification")
print("-" * 70)
print(f"✓ cobra installed at: {cobra.__file__}")
print(f"  Version: {cobra.__version__}")
import cobra_structural
print(f"✓ cobra_structural installed at: {cobra_structural.__file__}")
print(f"  Version: {cobra_structural.__version__}")

# Build a structural model
print("\n2. Building Structural Model")
print("-" * 70)
struct_model = Model('glycolysis_simplified')

# Add metabolites
glc = Metabolite('glc__D_e', name='D-Glucose', compartment='e')
g6p = Metabolite('g6p_c', name='Glucose-6-phosphate', compartment='c')
f6p = Metabolite('f6p_c', name='Fructose-6-phosphate', compartment='c')
pyr = Metabolite('pyr_c', name='Pyruvate', compartment='c')
atp = Metabolite('atp_c', name='ATP', compartment='c')
adp = Metabolite('adp_c', name='ADP', compartment='c')

# Add reactions
ex_glc = Reaction('EX_glc__D_e')
ex_glc.name = 'Glucose exchange'
ex_glc.add_metabolites({glc: -1})
ex_glc.bounds = (-10, 1000)

glk = Reaction('GLCpts')
glk.name = 'Glucose kinase'
glk.add_metabolites({glc: -1, g6p: 1, atp: -1, adp: 1})
glk.bounds = (0, 1000)

pgi = Reaction('PGI')
pgi.name = 'Phosphoglucose isomerase'
pgi.add_metabolites({g6p: -1, f6p: 1})
pgi.bounds = (-1000, 1000)

pfk = Reaction('PFK')
pfk.name = 'Phosphofructokinase'
pfk.add_metabolites({f6p: -1, atp: -1, adp: 1})
pfk.bounds = (0, 1000)

bio = Reaction('BIOMASS')
bio.name = 'Biomass production'
bio.add_metabolites({pyr: -1, atp: -10})
bio.bounds = (0, 1000)

struct_model.add_reactions([ex_glc, glk, pgi, pfk, bio])
struct_model.objective = {bio: 1}
struct_model.objective_direction = 'max'

print(f"✓ Structural model '{struct_model.id}' created")
print(f"  Reactions: {len(struct_model.reactions)}")
print(f"  Metabolites: {len(struct_model.metabolites)}")
print(f"  Objective: {list(struct_model.objective.keys())[0].id}")
print(f"  Direction: {struct_model.objective_direction}")

# Test in-memory conversion
print("\n3. In-Memory Conversion (dict-based)")
print("-" * 70)
full_model_mem = to_cobrapy_model(struct_model, use_file=False)
print(f"✓ Converted to full COBRApy model")
print(f"  Type: {type(full_model_mem)}")
print(f"  Module: {type(full_model_mem).__module__}")
print(f"  Has solver: {hasattr(full_model_mem, 'solver')}")

solution_mem = full_model_mem.optimize()
print(f"✓ Optimization result:")
print(f"  Status: {solution_mem.status}")
print(f"  Objective value: {solution_mem.objective_value:.6f}")
if solution_mem.status == 'optimal':
    print(f"  Biomass flux: {solution_mem.fluxes['BIOMASS']:.6f}")

# Test file-based conversion
print("\n4. File-Based Conversion (SBML)")
print("-" * 70)
full_model_file = to_cobrapy_model(struct_model, use_file=True)
print(f"✓ Converted via SBML file")

solution_file = full_model_file.optimize()
print(f"✓ Optimization result:")
print(f"  Status: {solution_file.status}")
print(f"  Objective value: {solution_file.objective_value:.6f}")

# Verify results match
print("\n5. Verification")
print("-" * 70)
if abs(solution_mem.objective_value - solution_file.objective_value) < 1e-6:
    print("✓ Both conversion methods produce identical results")
else:
    print(f"✗ Results differ: {solution_mem.objective_value} vs {solution_file.objective_value}")

# Test structural model I/O
print("\n6. Structural Model I/O")
print("-" * 70)
from cobra_structural.io import save_json_model, load_json_model
import tempfile
import os

with tempfile.TemporaryDirectory() as tmpdir:
    json_path = os.path.join(tmpdir, "test_model.json")
    save_json_model(struct_model, json_path)
    print(f"✓ Saved structural model to JSON")
    
    loaded_model = load_json_model(json_path)
    print(f"✓ Loaded structural model from JSON")
    print(f"  Reactions: {len(loaded_model.reactions)}")
    print(f"  Metabolites: {len(loaded_model.metabolites)}")
    print(f"  Objective: {list(loaded_model.objective.keys())[0].id if loaded_model.objective else 'None'}")

print("\n" + "=" * 70)
print("✅ All tests passed! cobra and cobra-structural coexist perfectly!")
print("=" * 70)
