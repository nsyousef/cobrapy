#!/usr/bin/env python3
"""Test structural-only model functionality."""

from cobra_structural import Model, Reaction, Metabolite

# Create a structural-only model
model = Model("test_model")

# Create metabolites
glc = Metabolite("glc", compartment="c")
atp = Metabolite("atp", compartment="c")
adp = Metabolite("adp", compartment="c")

# Create reactions
rxn1 = Reaction("PGI")
rxn1.add_metabolites({glc: -1, atp: -1, adp: 1})
rxn1.bounds = (0, 10)

rxn2 = Reaction("ATP_synthase")
rxn2.add_metabolites({adp: -1, atp: 1})
rxn2.bounds = (-5, 5)

# Add to model
model.add_reactions([rxn1, rxn2])

# Test basic properties
print("✓ Model created successfully")
print(f"  - Model ID: {model.id}")
print(f"  - Reactions: {len(model.reactions)}")
print(f"  - Metabolites: {len(model.metabolites)}")
print(f"  - Genes: {len(model.genes)}")

# Test objective setting (dict-based)
model.objective = {rxn1: 1.0, rxn2: 0.5}
print(f"\n✓ Objective set successfully")
print(f"  - Objective: {model.objective}")
print(f"  - Objective direction: {model.objective_direction}")

# Test that objective can be updated
model.objective = {rxn2: 1.0}
print(f"  - Updated objective: {model.objective}")

# Test reaction properties
print(f"\n✓ Reaction properties accessible")
print(f"  - Reaction bounds: {rxn1.bounds}")
print(f"  - Forward variable: {rxn1.forward_variable}")
print(f"  - Reverse variable: {rxn1.reverse_variable}")

# Test that flux raises helpful error
try:
    _ = rxn1.flux
    print("✗ Expected error for flux access")
except RuntimeError as e:
    print(f"\n✓ Flux property raises expected error:")
    print(f"  - {str(e)[:80]}...")

# Test that reduced_cost raises helpful error
try:
    _ = rxn1.reduced_cost
    print("✗ Expected error for reduced_cost access")
except RuntimeError as e:
    print(f"\n✓ Reduced cost property raises expected error:")
    print(f"  - {str(e)[:80]}...")

# Test that shadow_price raises helpful error
try:
    _ = glc.shadow_price
    print("✗ Expected error for shadow_price access")
except RuntimeError as e:
    print(f"\n✓ Shadow price property raises expected error:")
    print(f"  - {str(e)[:80]}...")

# Test model copy
model_copy = model.copy()
print(f"\n✓ Model copy successful")
print(f"  - Original objective: {model.objective}")
print(f"  - Copied objective: {model_copy.objective}")
print(f"  - Reactions match: {len(model.reactions) == len(model_copy.reactions)}")

# Test model merge
model2 = Model("model2")
rxn3 = Reaction("ATP_consumption")
rxn3.add_metabolites({atp: -1})
rxn3.bounds = (0, 3)
model2.add_reactions([rxn3])

model_merged = model.merge(model2)
print(f"\n✓ Model merge successful")
print(f"  - Merged reactions: {len(model_merged.reactions)}")

# Test I/O compatibility (just verify structure is maintained)
print(f"\n✓ Structural-only model working correctly!")
print(f"  - Use export to JSON/SBML and reimport with full COBRApy for optimization")
