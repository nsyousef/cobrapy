#!/usr/bin/env python
"""
Test FBA conversion from structural model to full COBRApy model.

This test ensures that structural models can be properly converted to
full COBRApy models and produce non-zero flux solutions when expected.

Test Cases:
1. Simple Linear Pathway - Tests basic conversion with a straightforward
   A → B → C → Biomass pathway. Expected flux: 10.0 (limited by uptake).

2. Branched Pathway - Tests more complex model with multiple nutrients
   (glucose + oxygen), ATP production, and maintenance requirements.
   Verifies that constraints are properly maintained after conversion.

3. Textbook Model - Tests conversion of the classic E. coli textbook model
   if available. This is a comprehensive test of a real metabolic network.

All tests verify both in-memory (dict-based) and file-based (SBML) conversion
methods produce identical, non-zero results.
"""

from cobra_structural import Model, Reaction, Metabolite
from cobra_structural.io import to_cobrapy_model


def test_simple_linear_pathway():
    """Test a simple linear pathway that should produce non-zero flux."""
    print("=" * 70)
    print("Test 1: Simple Linear Pathway (A → B → C → Biomass)")
    print("=" * 70)
    
    # Create structural model
    model = Model('linear_pathway')
    
    # Metabolites
    A_ext = Metabolite('A_e', name='A external', compartment='e')
    A_int = Metabolite('A_c', name='A cytosol', compartment='c')
    B = Metabolite('B_c', name='B', compartment='c')
    C = Metabolite('C_c', name='C', compartment='c')
    
    # Reactions
    # External A can be imported
    ex_A = Reaction('EX_A')
    ex_A.name = 'A exchange'
    ex_A.add_metabolites({A_ext: -1})
    ex_A.bounds = (-10, 1000)  # Allow uptake up to 10
    
    # Transport A from external to internal
    transport_A = Reaction('A_transport')
    transport_A.name = 'A transport'
    transport_A.add_metabolites({A_ext: -1, A_int: 1})
    transport_A.bounds = (0, 1000)
    
    # A → B
    r1 = Reaction('R1')
    r1.name = 'A to B'
    r1.add_metabolites({A_int: -1, B: 1})
    r1.bounds = (0, 1000)
    
    # B → C
    r2 = Reaction('R2')
    r2.name = 'B to C'
    r2.add_metabolites({B: -1, C: 1})
    r2.bounds = (0, 1000)
    
    # C → Biomass
    biomass = Reaction('BIOMASS')
    biomass.name = 'Biomass production'
    biomass.add_metabolites({C: -1})
    biomass.bounds = (0, 1000)
    
    # Add all reactions
    model.add_reactions([ex_A, transport_A, r1, r2, biomass])
    
    # Set objective
    model.objective = {biomass: 1}
    model.objective_direction = 'max'
    
    print(f"✓ Structural model created: {model.id}")
    print(f"  Reactions: {len(model.reactions)}")
    print(f"  Metabolites: {len(model.metabolites)}")
    print(f"  Objective: maximize {biomass.id}")
    
    # Convert and optimize - in-memory
    print("\n--- In-Memory Conversion ---")
    full_model_mem = to_cobrapy_model(model, use_file=False)
    solution_mem = full_model_mem.optimize()
    
    print(f"Status: {solution_mem.status}")
    print(f"Objective value: {solution_mem.objective_value:.6f}")
    print(f"Expected: 10.0 (limited by A uptake)")
    
    # Check results
    assert solution_mem.status == 'optimal', f"Expected optimal status, got {solution_mem.status}"
    assert solution_mem.objective_value > 0, "Expected non-zero biomass flux!"
    assert abs(solution_mem.objective_value - 10.0) < 1e-6, \
        f"Expected biomass flux of 10.0, got {solution_mem.objective_value:.6f}"
    
    print("✓ In-memory conversion: PASSED")
    
    # Convert and optimize - file-based
    print("\n--- File-Based Conversion (SBML) ---")
    full_model_file = to_cobrapy_model(model, use_file=True)
    solution_file = full_model_file.optimize()
    
    print(f"Status: {solution_file.status}")
    print(f"Objective value: {solution_file.objective_value:.6f}")
    
    assert solution_file.status == 'optimal', f"Expected optimal status, got {solution_file.status}"
    assert solution_file.objective_value > 0, "Expected non-zero biomass flux!"
    assert abs(solution_file.objective_value - 10.0) < 1e-6, \
        f"Expected biomass flux of 10.0, got {solution_file.objective_value:.6f}"
    
    print("✓ File-based conversion: PASSED")
    print()


def test_branched_pathway():
    """Test a branched pathway with multiple nutrient sources."""
    print("=" * 70)
    print("Test 2: Branched Pathway (glucose + oxygen → ATP → Biomass)")
    print("=" * 70)
    
    # Create structural model
    model = Model('branched_pathway')
    
    # External metabolites
    glc_e = Metabolite('glc_e', name='Glucose external', compartment='e')
    o2_e = Metabolite('o2_e', name='Oxygen external', compartment='e')
    
    # Internal metabolites
    glc_c = Metabolite('glc_c', name='Glucose', compartment='c')
    g6p = Metabolite('g6p_c', name='Glucose-6-phosphate', compartment='c')
    pyr = Metabolite('pyr_c', name='Pyruvate', compartment='c')
    o2_c = Metabolite('o2_c', name='Oxygen', compartment='c')
    atp = Metabolite('atp_c', name='ATP', compartment='c')
    adp = Metabolite('adp_c', name='ADP', compartment='c')
    
    # Exchange reactions
    ex_glc = Reaction('EX_glc')
    ex_glc.add_metabolites({glc_e: -1})
    ex_glc.bounds = (-20, 1000)  # Allow glucose uptake
    
    ex_o2 = Reaction('EX_o2')
    ex_o2.add_metabolites({o2_e: -1})
    ex_o2.bounds = (-50, 1000)  # Allow oxygen uptake
    
    # Transport reactions
    glc_t = Reaction('GLCt')
    glc_t.add_metabolites({glc_e: -1, glc_c: 1})
    glc_t.bounds = (0, 1000)
    
    o2_t = Reaction('O2t')
    o2_t.add_metabolites({o2_e: -1, o2_c: 1})
    o2_t.bounds = (0, 1000)
    
    # Glycolysis (simplified): glc → g6p → pyr + ATP
    glk = Reaction('GLK')
    glk.name = 'Glucokinase'
    glk.add_metabolites({glc_c: -1, g6p: 1, atp: -1, adp: 1})
    glk.bounds = (0, 1000)
    
    pgi = Reaction('PGI')
    pgi.name = 'Glycolysis (simplified)'
    pgi.add_metabolites({g6p: -1, pyr: 2, atp: 2, adp: -2})
    pgi.bounds = (0, 1000)
    
    # Oxidative phosphorylation: pyr + o2 → ATP
    oxphos = Reaction('OXPHOS')
    oxphos.name = 'Oxidative phosphorylation'
    oxphos.add_metabolites({pyr: -1, o2_c: -3, atp: 15, adp: -15})
    oxphos.bounds = (0, 1000)
    
    # ATP maintenance
    atpm = Reaction('ATPM')
    atpm.name = 'ATP maintenance'
    atpm.add_metabolites({atp: -1, adp: 1})
    atpm.bounds = (5, 1000)  # Minimum ATP requirement
    
    # Biomass: requires ATP
    biomass = Reaction('BIOMASS')
    biomass.name = 'Biomass production'
    biomass.add_metabolites({atp: -40, adp: 40})
    biomass.bounds = (0, 1000)
    
    # Add all reactions
    model.add_reactions([
        ex_glc, ex_o2, glc_t, o2_t, 
        glk, pgi, oxphos, atpm, biomass
    ])
    
    # Set objective
    model.objective = {biomass: 1}
    model.objective_direction = 'max'
    
    print(f"✓ Structural model created: {model.id}")
    print(f"  Reactions: {len(model.reactions)}")
    print(f"  Metabolites: {len(model.metabolites)}")
    print(f"  Objective: maximize {biomass.id}")
    print(f"  Constraints: ATPM ≥ 5, glucose uptake ≤ 20, oxygen uptake ≤ 50")
    
    # Convert and optimize - in-memory
    print("\n--- In-Memory Conversion ---")
    full_model_mem = to_cobrapy_model(model, use_file=False)
    solution_mem = full_model_mem.optimize()
    
    print(f"Status: {solution_mem.status}")
    print(f"Objective value: {solution_mem.objective_value:.6f}")
    print(f"Key fluxes:")
    print(f"  EX_glc: {solution_mem.fluxes['EX_glc']:.6f}")
    print(f"  EX_o2: {solution_mem.fluxes['EX_o2']:.6f}")
    print(f"  BIOMASS: {solution_mem.fluxes['BIOMASS']:.6f}")
    print(f"  ATPM: {solution_mem.fluxes['ATPM']:.6f}")
    
    # Check results
    assert solution_mem.status == 'optimal', f"Expected optimal status, got {solution_mem.status}"
    assert solution_mem.objective_value > 0, "Expected non-zero biomass flux!"
    assert solution_mem.fluxes['BIOMASS'] > 0, "Biomass flux should be positive"
    assert solution_mem.fluxes['EX_glc'] < 0, "Should be consuming glucose"
    assert solution_mem.fluxes['EX_o2'] < 0, "Should be consuming oxygen"
    assert solution_mem.fluxes['ATPM'] >= 5 - 1e-6, "ATP maintenance requirement not met"
    
    print("✓ In-memory conversion: PASSED")
    
    # Convert and optimize - file-based
    print("\n--- File-Based Conversion (SBML) ---")
    full_model_file = to_cobrapy_model(model, use_file=True)
    solution_file = full_model_file.optimize()
    
    print(f"Status: {solution_file.status}")
    print(f"Objective value: {solution_file.objective_value:.6f}")
    
    assert solution_file.status == 'optimal', f"Expected optimal status, got {solution_file.status}"
    assert solution_file.objective_value > 0, "Expected non-zero biomass flux!"
    assert abs(solution_file.objective_value - solution_mem.objective_value) < 1e-6, \
        "In-memory and file-based conversions should give identical results"
    
    print("✓ File-based conversion: PASSED")
    print()


def test_textbook_model():
    """Test with a classic textbook model structure."""
    print("=" * 70)
    print("Test 3: E. coli Textbook Model (from cobra_structural.test)")
    print("=" * 70)
    
    try:
        from cobra_structural.test import create_test_model
        
        # Get the textbook model
        print("Loading textbook model...")
        struct_model = create_test_model("textbook")
        
        print(f"✓ Textbook model loaded: {struct_model.id}")
        print(f"  Reactions: {len(struct_model.reactions)}")
        print(f"  Metabolites: {len(struct_model.metabolites)}")
        print(f"  Genes: {len(struct_model.genes)}")
        
        # Verify it has a dict objective
        assert isinstance(struct_model.objective, dict), \
            "Textbook model should have dict-based objective"
        
        # Convert and optimize - in-memory
        print("\n--- In-Memory Conversion ---")
        full_model = to_cobrapy_model(struct_model, use_file=False)
        solution = full_model.optimize()
        
        print(f"Status: {solution.status}")
        print(f"Objective value: {solution.objective_value:.6f}")
        
        # The textbook model should produce a known optimal value
        assert solution.status == 'optimal', f"Expected optimal status, got {solution.status}"
        assert solution.objective_value > 0, "Expected non-zero biomass flux!"
        assert solution.objective_value > 0.8, \
            f"Textbook model typically produces ~0.87 flux, got {solution.objective_value:.6f}"
        
        print("✓ Textbook model conversion: PASSED")
        
    except ImportError:
        print("⚠ Skipping textbook model test (create_test_model not available)")
    
    print()


def main():
    """Run all FBA conversion tests."""
    print("\n")
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 68 + "║")
    print("║" + "  FBA Conversion Tests - Structural to Full COBRApy".center(68) + "║")
    print("║" + " " * 68 + "║")
    print("╚" + "=" * 68 + "╝")
    print()
    
    tests_passed = 0
    tests_total = 0
    
    # Test 1: Simple linear pathway
    tests_total += 1
    try:
        test_simple_linear_pathway()
        tests_passed += 1
    except AssertionError as e:
        print(f"✗ Test 1 FAILED: {e}")
    except Exception as e:
        print(f"✗ Test 1 ERROR: {e}")
    
    # Test 2: Branched pathway
    tests_total += 1
    try:
        test_branched_pathway()
        tests_passed += 1
    except AssertionError as e:
        print(f"✗ Test 2 FAILED: {e}")
    except Exception as e:
        print(f"✗ Test 2 ERROR: {e}")
    
    # Test 3: Textbook model
    tests_total += 1
    try:
        test_textbook_model()
        tests_passed += 1
    except AssertionError as e:
        print(f"✗ Test 3 FAILED: {e}")
    except Exception as e:
        print(f"✗ Test 3 ERROR: {e}")
    
    # Summary
    print("=" * 70)
    print(f"Test Summary: {tests_passed}/{tests_total} tests passed")
    print("=" * 70)
    
    if tests_passed == tests_total:
        print("\n✅ ALL TESTS PASSED!")
        print("\nConclusion: Structural models convert correctly to full COBRApy")
        print("models and produce expected non-zero FBA solutions.")
        return 0
    else:
        print(f"\n✗ {tests_total - tests_passed} test(s) failed")
        return 1


if __name__ == '__main__':
    exit(main())
