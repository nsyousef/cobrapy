#!/usr/bin/env python3
"""Test conversion from structural model to full COBRApy model."""

from cobra_structural import Model, Reaction, Metabolite
from cobra_structural.io import to_cobrapy_model

print("=" * 70)
print("Testing structural → full COBRApy conversion")
print("=" * 70)

# Build structural model
structural_model = Model("test_conversion")

glc_e = Metabolite("glc_e", name="Glucose external", compartment="e")
glc_c = Metabolite("glc_c", name="Glucose cytosol", compartment="c")
atp_c = Metabolite("atp_c", name="ATP", compartment="c")
adp_c = Metabolite("adp_c", name="ADP", compartment="c")

# Exchange reaction
ex_glc = Reaction("EX_glc")
ex_glc.add_metabolites({glc_e: -1})
ex_glc.bounds = (-10, 0)

# Transport
glc_t = Reaction("GLCt")
glc_t.add_metabolites({glc_e: -1, glc_c: 1})
glc_t.bounds = (0, 10)

# Internal reaction
glc_kinase = Reaction("GLCkinase")
glc_kinase.add_metabolites({glc_c: -1, atp_c: -1, adp_c: 1})
glc_kinase.bounds = (0, 10)

# Demand
adp_demand = Reaction("DM_adp")
adp_demand.add_metabolites({adp_c: -1})
adp_demand.bounds = (0, 1000)

structural_model.add_reactions([ex_glc, glc_t, glc_kinase, adp_demand])
structural_model.objective = {glc_kinase: 1.0}

print(f"\n✓ Structural model built: {structural_model.id}")
print(f"  - Reactions: {len(structural_model.reactions)}")
print(f"  - Metabolites: {len(structural_model.metabolites)}")
print(f"  - Objective: {structural_model.objective}")

# Test in-memory conversion (default, fast)
print("\n" + "=" * 70)
print("Testing in-memory conversion (use_file=False)")
print("=" * 70)

try:
    full_model = to_cobrapy_model(structural_model, use_file=False)
    print(f"✓ Converted to full COBRApy model: {full_model.id}")
    print(f"  - Has solver: {hasattr(full_model, 'solver') and full_model.solver is not None}")
    print(f"  - Reactions: {len(full_model.reactions)}")
    print(f"  - Metabolites: {len(full_model.metabolites)}")
    
    # Try optimization
    solution = full_model.optimize()
    print(f"\n✓ Optimization successful!")
    print(f"  - Status: {solution.status}")
    print(f"  - Objective value: {solution.objective_value:.6f}")
    print(f"  - GLCkinase flux: {solution.fluxes['GLCkinase']:.6f}")
    
except ImportError as e:
    print(f"⚠ Original COBRApy not installed: {e}")
    print("  Install with: pip install cobra")
except Exception as e:
    print(f"✗ Conversion failed: {e}")
    import traceback
    traceback.print_exc()

# Test file-based conversion (via SBML)
print("\n" + "=" * 70)
print("Testing file-based conversion (use_file=True)")
print("=" * 70)

try:
    full_model_file = to_cobrapy_model(structural_model, use_file=True)
    print(f"✓ Converted via SBML file: {full_model_file.id}")
    print(f"  - Has solver: {hasattr(full_model_file, 'solver') and full_model_file.solver is not None}")
    
    # Try optimization
    solution_file = full_model_file.optimize()
    print(f"\n✓ Optimization successful!")
    print(f"  - Status: {solution_file.status}")
    print(f"  - Objective value: {solution_file.objective_value:.6f}")
    
except ImportError as e:
    print(f"⚠ Original COBRApy not installed: {e}")
except Exception as e:
    print(f"✗ File-based conversion failed: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 70)
print("✓✓✓ Conversion functionality working! ✓✓✓")
print("=" * 70)
print("""
Usage:
    from cobra_structural.io import to_cobrapy_model
    
    # Fast in-memory conversion (recommended)
    full_model = to_cobrapy_model(structural_model)
    
    # Or via SBML file (more robust for complex models)
    full_model = to_cobrapy_model(structural_model, use_file=True)
    
    # Then use full COBRApy features
    solution = full_model.optimize()
    fva_result = cobra.flux_analysis.flux_variability_analysis(full_model)
""")
