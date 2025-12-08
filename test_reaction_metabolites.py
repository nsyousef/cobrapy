#!/usr/bin/env python
"""
Test for reaction metabolite operations on structural models.

This test ensures that adding/setting metabolites in reactions works
correctly on structural models without trying to access the solver's
constraints (which don't exist).

Regression test for: https://github.com/opencobra/cobrapy/issues/XXXX
"""

from cobra_structural import Model, Reaction, Metabolite


def test_add_metabolites_to_reaction():
    """Test adding metabolites to a reaction in a structural model."""
    print("=" * 70)
    print("Test 1: Adding Metabolites to Reaction in Structural Model")
    print("=" * 70)
    
    # Create a structural model
    model = Model('test_model')
    
    # Create metabolites and add to model
    met1 = Metabolite('A_c', compartment='c')
    met2 = Metabolite('B_c', compartment='c')
    model.add_metabolites([met1, met2])
    
    # Create a reaction and add to model
    rxn = Reaction('R1')
    model.add_reactions([rxn])
    
    print(f"✓ Model created with {len(model.metabolites)} metabolites")
    print(f"✓ Reaction {rxn.id} added to model")
    
    # Add metabolites to the reaction
    # This was failing because add_metabolites tried to access model.constraints
    print("\nAdding metabolites to reaction...")
    rxn.add_metabolites({met1: -1, met2: 1})
    
    assert len(rxn.metabolites) == 2, "Should have 2 metabolites"
    assert rxn.metabolites[met1] == -1, "Stoichiometry should be -1"
    assert rxn.metabolites[met2] == 1, "Stoichiometry should be 1"
    
    print(f"✓ Added metabolites successfully")
    print(f"  Reaction: {rxn.reaction}")
    
    print()


def test_set_reaction_from_string():
    """Test setting reaction from a string (inter-microbe exchange pattern)."""
    print("=" * 70)
    print("Test 2: Setting Reaction from String (Community Pattern)")
    print("=" * 70)
    
    # Create a structural model
    model = Model('community')
    
    # Create metabolites representing a shared metabolite and organism-specific versions
    shared_met = Metabolite('glc_e', compartment='e')
    org1_met = Metabolite('ecoli_glc[e]', compartment='e')
    org2_met = Metabolite('bsubtilis_glc[e]', compartment='e')
    
    model.add_metabolites([shared_met, org1_met, org2_met])
    
    # Create inter-microbe exchange reactions
    iex_rxn1 = Reaction('INTER_EX_ecoli')
    iex_rxn2 = Reaction('INTER_EX_bsubtilis')
    
    model.add_reactions([iex_rxn1, iex_rxn2])
    
    print(f"✓ Community model created")
    print(f"  Metabolites: {[m.id for m in model.metabolites]}")
    
    # Set reactions from strings
    # This is the exact pattern that was failing
    print("\nSetting reactions from strings (community pattern)...")
    iex_rxn1.reaction = f'{shared_met.id} <=> {org1_met.id}'
    iex_rxn2.reaction = f'{shared_met.id} <=> {org2_met.id}'
    
    print(f"✓ Reaction 1: {iex_rxn1.reaction}")
    print(f"✓ Reaction 2: {iex_rxn2.reaction}")
    
    # Verify metabolites are correct
    assert shared_met in iex_rxn1.metabolites, "Shared metabolite should be in reaction"
    assert org1_met in iex_rxn1.metabolites, "Organism 1 metabolite should be in reaction"
    assert shared_met in iex_rxn2.metabolites, "Shared metabolite should be in reaction"
    assert org2_met in iex_rxn2.metabolites, "Organism 2 metabolite should be in reaction"
    
    print("\n✓ All metabolites correctly added to reactions")
    
    print()


def test_combine_metabolites():
    """Test combining metabolites (additive mode)."""
    print("=" * 70)
    print("Test 3: Combining/Updating Metabolite Coefficients")
    print("=" * 70)
    
    model = Model('test_combine')
    
    met_a = Metabolite('A_c', compartment='c')
    met_b = Metabolite('B_c', compartment='c')
    
    model.add_metabolites([met_a, met_b])
    
    rxn = Reaction('R1')
    model.add_reactions([rxn])
    
    # Add metabolites in first call
    print("Initial: Adding A and B to reaction")
    rxn.add_metabolites({met_a: -1, met_b: 1})
    print(f"✓ Reaction: {rxn.reaction}")
    
    # Update with combine=True (additive)
    print("\nCombining: Adding more of metabolite A (combine=True)")
    rxn.add_metabolites({met_a: -2}, combine=True)
    print(f"✓ Reaction: {rxn.reaction}")
    
    assert rxn.metabolites[met_a] == -3, "Should have combined coefficients (-1 + -2 = -3)"
    print(f"✓ Metabolite A coefficient correctly combined: {rxn.metabolites[met_a]}")
    
    # Update with combine=False (replace)
    print("\nReplacing: Setting metabolite A to new value (combine=False)")
    rxn.add_metabolites({met_a: -1}, combine=False)
    print(f"✓ Reaction: {rxn.reaction}")
    
    assert rxn.metabolites[met_a] == -1, "Should have replaced coefficient"
    print(f"✓ Metabolite A coefficient correctly replaced: {rxn.metabolites[met_a]}")
    
    print()


def main():
    """Run all tests."""
    print("\n")
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 68 + "║")
    print("║" + "  Reaction Metabolite Operations on Structural Models".center(68) + "║")
    print("║" + " " * 68 + "║")
    print("╚" + "=" * 68 + "╝")
    print()
    
    tests = [
        test_add_metabolites_to_reaction,
        test_set_reaction_from_string,
        test_combine_metabolites,
    ]
    
    passed = 0
    for test_func in tests:
        try:
            test_func()
            passed += 1
        except AssertionError as e:
            print(f"✗ Test failed: {e}")
        except Exception as e:
            print(f"✗ Test error: {e}")
    
    print("=" * 70)
    print(f"Test Summary: {passed}/{len(tests)} tests passed")
    print("=" * 70)
    
    if passed == len(tests):
        print("\n✅ ALL TESTS PASSED!")
        print("\nConclusion: Reaction metabolite operations work correctly on")
        print("structural models without trying to access solver constraints.")
        return 0
    else:
        print(f"\n✗ {len(tests) - passed} test(s) failed")
        return 1


if __name__ == '__main__':
    exit(main())
