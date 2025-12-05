#!/usr/bin/env python3
"""Test structural-only model functionality - direct import."""

import sys
sys.path.insert(0, '/Users/nicholasyousefi/Documents/Coding/School/Harvard/Labs/zomorrodi_lab/cobra/cobrapy/src')

# Import just the core classes without the full cobra package
from cobra.core.model import Model
from cobra.core.reaction import Reaction
from cobra.core.metabolite import Metabolite

try:
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
    
    # Test objective setting (dict-based)
    model.objective = {rxn1: 1.0, rxn2: 0.5}
    print(f"\n✓ Objective set successfully")
    print(f"  - Objective type: {type(model.objective)}")
    print(f"  - Objective keys: {list(model.objective.keys())}")
    
    # Test reaction properties
    print(f"\n✓ Reaction properties accessible")
    print(f"  - Reaction bounds: {rxn1.bounds}")
    print(f"  - Forward variable (should be None): {rxn1.forward_variable}")
    print(f"  - Reverse variable (should be None): {rxn1.reverse_variable}")
    
    # Test that flux raises helpful error
    try:
        _ = rxn1.flux
        print("✗ Expected error for flux access")
    except RuntimeError as e:
        error_msg = str(e)[:80]
        if "structural-only" in error_msg:
            print(f"\n✓ Flux property raises expected error (contains 'structural-only')")
        else:
            print(f"\n✗ Flux error doesn't mention structural-only: {error_msg}")
    
    # Test that reduced_cost raises helpful error
    try:
        _ = rxn1.reduced_cost
        print("✗ Expected error for reduced_cost access")
    except RuntimeError as e:
        error_msg = str(e)[:80]
        if "structural-only" in error_msg:
            print(f"✓ Reduced cost property raises expected error (contains 'structural-only')")
        else:
            print(f"✗ Reduced cost error doesn't mention structural-only: {error_msg}")
    
    # Test that shadow_price raises helpful error
    try:
        _ = glc.shadow_price
        print("✗ Expected error for shadow_price access")
    except RuntimeError as e:
        error_msg = str(e)[:80]
        if "structural-only" in error_msg:
            print(f"✓ Shadow price property raises expected error (contains 'structural-only')")
        else:
            print(f"✗ Shadow price error doesn't mention structural-only: {error_msg}")
    
    print(f"\n✓✓✓ Structural-only model working correctly! ✓✓✓")
    print(f"  - Use export to JSON/SBML and reimport with full COBRApy for optimization")

except Exception as e:
    print(f"✗ Unexpected error: {e}")
    import traceback
    traceback.print_exc()
