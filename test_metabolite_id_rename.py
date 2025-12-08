#!/usr/bin/env python
"""
Test for metabolite ID renaming on structural models.

This test ensures that metabolite IDs can be changed on structural models
without trying to access the solver's constraints (which don't exist).

Regression test for: https://github.com/opencobra/cobrapy/issues/XXXX
"""

from cobra_structural import Model, Reaction, Metabolite


def test_metabolite_id_rename_structural_model():
    """Test that metabolite IDs can be renamed on structural models."""
    print("=" * 70)
    print("Test: Metabolite ID Renaming on Structural Model")
    print("=" * 70)
    
    # Create a structural model
    model = Model('test_model')
    
    # Add metabolites
    met1 = Metabolite('A_c', compartment='c')
    met2 = Metabolite('B_c', compartment='c')
    model.add_metabolites([met1, met2])
    
    print(f"✓ Initial metabolites: {[m.id for m in model.metabolites]}")
    
    # Test simple ID change
    print("\nTest 1: Simple ID change")
    original_id = met1.id
    new_id = 'X_c'
    met1.id = new_id
    
    assert met1.id == new_id, f"Expected {new_id}, got {met1.id}"
    assert original_id not in [m.id for m in model.metabolites], \
        f"Old ID {original_id} should not be in model"
    assert new_id in [m.id for m in model.metabolites], \
        f"New ID {new_id} should be in model"
    print(f"✓ ID changed from {original_id} to {new_id}")
    
    # Test community gem builder pattern (the original bug)
    print("\nTest 2: Community gem builder tagging pattern")
    model2 = Model('microbe')
    
    # Create metabolites with standard naming
    metabolites = [
        Metabolite('glc_c', compartment='c'),
        Metabolite('atp_c', compartment='c'),
        Metabolite('h2o_e', compartment='e'),
    ]
    model2.add_metabolites(metabolites)
    
    # Add reactions to ensure metabolite references are maintained
    rxn = Reaction('GLK')
    rxn.add_metabolites({metabolites[0]: -1, metabolites[1]: 1})
    model2.add_reactions([rxn])
    
    print(f"Original metabolites: {[m.id for m in model2.metabolites]}")
    
    # Apply community gem builder tagging (the pattern that was failing)
    microbe_name = 'ecoli'
    for met in model2.metabolites:
        compartment = met.compartment
        # Remove compartment notation from the name
        no_c_name = (met.id
                    .replace(f'[{compartment}]', '')
                    .replace(f'_{compartment}', ''))
        # Add microbe prefix - THIS WAS FAILING WITH AttributeError
        met.id = f'{microbe_name}_{no_c_name}[{compartment}]'
    
    print(f"Tagged metabolites: {[m.id for m in model2.metabolites]}")
    
    # Verify all metabolites were renamed
    expected_ids = {
        'ecoli_glc[c]',
        'ecoli_atp[c]',
        'ecoli_h2o[e]',
    }
    actual_ids = {m.id for m in model2.metabolites}
    assert expected_ids == actual_ids, \
        f"Expected {expected_ids}, got {actual_ids}"
    
    # Verify reaction still references metabolites correctly
    assert metabolites[0].id == 'ecoli_glc[c]', \
        "Reaction should still reference updated metabolite ID"
    assert rxn.metabolites[metabolites[0]] == -1, \
        "Reaction stoichiometry should be preserved"
    
    print("✓ All metabolites tagged successfully")
    print("✓ Reaction metabolite references updated correctly")
    
    print("\n" + "=" * 70)
    print("✅ All tests passed!")
    print("=" * 70)


if __name__ == '__main__':
    test_metabolite_id_rename_structural_model()
