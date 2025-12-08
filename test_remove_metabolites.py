#!/usr/bin/env python
"""
Test for removing metabolites and reactions from structural models.

This test ensures that removing metabolites with destructive=True works
correctly without trying to access None model references.

Regression test for: AttributeError when removing metabolites destructively
"""

from cobra_structural import Model, Reaction, Metabolite


def test_remove_metabolites_non_destructive():
    """Test removing metabolites while keeping reactions (default mode)."""
    print("=" * 70)
    print("Test 1: Remove Metabolites Non-Destructively")
    print("=" * 70)
    
    model = Model('test_model')
    
    # Create metabolites
    met_a = Metabolite('A_c', compartment='c')
    met_b = Metabolite('B_c', compartment='c')
    met_c = Metabolite('C_c', compartment='c')
    
    model.add_metabolites([met_a, met_b, met_c])
    
    # Create reactions
    r1 = Reaction('R1')
    r1.add_metabolites({met_a: -1, met_b: 1})
    
    r2 = Reaction('R2')
    r2.add_metabolites({met_b: -1, met_c: 1})
    
    model.add_reactions([r1, r2])
    
    print(f"Initial: {len(model.metabolites)} metabolites, {len(model.reactions)} reactions")
    print(f"  R1: {r1.reaction}")
    print(f"  R2: {r2.reaction}")
    
    # Remove metabolite non-destructively
    # This should remove met_b from reactions, but keep the reactions
    print(f"\nRemoving met_b non-destructively...")
    model.remove_metabolites([met_b], destructive=False)
    
    print(f"After: {len(model.metabolites)} metabolites, {len(model.reactions)} reactions")
    print(f"  R1: {r1.reaction}")
    print(f"  R2: {r2.reaction}")
    
    # Verify metabolite was removed
    assert met_b not in model.metabolites, "met_b should be removed from model"
    assert len(model.metabolites) == 2, "Should have 2 metabolites left"
    
    # Verify reactions still exist
    assert len(model.reactions) == 2, "Both reactions should still exist"
    assert met_b not in r1.metabolites, "met_b should be removed from R1"
    assert met_b not in r2.metabolites, "met_b should be removed from R2"
    
    print("✓ Metabolite removed, reactions preserved")
    print()


def test_remove_metabolites_destructive():
    """Test removing metabolites and their associated reactions (destructive mode)."""
    print("=" * 70)
    print("Test 2: Remove Metabolites Destructively")
    print("=" * 70)
    
    model = Model('test_model')
    
    # Create metabolites
    met_x = Metabolite('X_c', compartment='c')
    met_y = Metabolite('Y_c', compartment='c')
    met_z = Metabolite('Z_c', compartment='c')
    
    model.add_metabolites([met_x, met_y, met_z])
    
    # Create reactions
    # R1 uses X and Y
    r1 = Reaction('R1')
    r1.add_metabolites({met_x: -1, met_y: 1})
    
    # R2 uses only Y (will be removed with Y)
    r2 = Reaction('R2')
    r2.add_metabolites({met_y: -1, met_z: 1})
    
    # R3 uses only Z (should remain)
    r3 = Reaction('R3')
    r3.add_metabolites({met_z: -1})
    
    model.add_reactions([r1, r2, r3])
    
    print(f"Initial: {len(model.metabolites)} metabolites, {len(model.reactions)} reactions")
    print(f"  R1: {r1.reaction}")
    print(f"  R2: {r2.reaction}")
    print(f"  R3: {r3.reaction}")
    
    # Remove metabolite destructively
    # This should remove met_y and all reactions that use met_y (R1, R2)
    print(f"\nRemoving met_y destructively...")
    model.remove_metabolites([met_y], destructive=True)
    
    print(f"After: {len(model.metabolites)} metabolites, {len(model.reactions)} reactions")
    print(f"  Remaining reactions: {[r.id for r in model.reactions]}")
    
    # Verify metabolite was removed
    assert met_y not in model.metabolites, "met_y should be removed"
    assert len(model.metabolites) == 2, "Should have 2 metabolites left"
    
    # Verify reactions using met_y were removed
    assert r1 not in model.reactions, "R1 (uses met_y) should be removed"
    assert r2 not in model.reactions, "R2 (uses met_y) should be removed"
    assert r3 in model.reactions, "R3 (doesn't use met_y) should remain"
    
    print("✓ Metabolite and associated reactions removed correctly")
    print()


def test_remove_multiple_metabolites_destructive():
    """Test removing multiple metabolites destructively (mgPipe/migemox pattern)."""
    print("=" * 70)
    print("Test 3: Remove Multiple Metabolites Destructively")
    print("=" * 70)
    
    model = Model('community')
    
    # Create metabolites representing a community model
    # Some might have zero abundance (to be pruned)
    metabolites = []
    for i in range(5):
        m = Metabolite(f'M{i}_c', compartment='c')
        metabolites.append(m)
    
    model.add_metabolites(metabolites)
    
    # Create reactions
    r1 = Reaction('R1')
    r1.add_metabolites({metabolites[0]: -1, metabolites[1]: 1})
    
    r2 = Reaction('R2')
    r2.add_metabolites({metabolites[1]: -1, metabolites[2]: 1})
    
    r3 = Reaction('R3')
    r3.add_metabolites({metabolites[3]: -1, metabolites[4]: 1})
    
    model.add_reactions([r1, r2, r3])
    
    print(f"Initial: {len(model.metabolites)} metabolites, {len(model.reactions)} reactions")
    
    # Simulate pruning zero-abundance metabolites
    # This is the exact pattern that was failing: prune_zero_abundance_microbe
    zero_abundance_metabolites = [metabolites[3], metabolites[4]]
    
    print(f"\nRemoving {len(zero_abundance_metabolites)} metabolites destructively...")
    model.remove_metabolites(zero_abundance_metabolites, destructive=True)
    
    print(f"After: {len(model.metabolites)} metabolites, {len(model.reactions)} reactions")
    print(f"  Remaining metabolites: {[m.id for m in model.metabolites]}")
    print(f"  Remaining reactions: {[r.id for r in model.reactions]}")
    
    # Verify
    assert len(model.metabolites) == 3, "Should have 3 metabolites"
    assert len(model.reactions) == 2, "Should have 2 reactions (R3 removed)"
    assert r3 not in model.reactions, "R3 should be removed"
    assert r1 in model.reactions, "R1 should remain"
    assert r2 in model.reactions, "R2 should remain"
    
    print("✓ Multiple metabolites removed correctly (migemox pattern)")
    print()


def main():
    """Run all tests."""
    print("\n")
    print("╔" + "=" * 68 + "╗")
    print("║" + " " * 68 + "║")
    print("║" + "  Removing Metabolites and Reactions from Structural Models".center(68) + "║")
    print("║" + " " * 68 + "║")
    print("╚" + "=" * 68 + "╝")
    print()
    
    tests = [
        test_remove_metabolites_non_destructive,
        test_remove_metabolites_destructive,
        test_remove_multiple_metabolites_destructive,
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
            import traceback
            traceback.print_exc()
    
    print("=" * 70)
    print(f"Test Summary: {passed}/{len(tests)} tests passed")
    print("=" * 70)
    
    if passed == len(tests):
        print("\n✅ ALL TESTS PASSED!")
        print("\nConclusion: Metabolite removal works correctly on structural models")
        print("in both destructive and non-destructive modes.")
        return 0
    else:
        print(f"\n✗ {len(tests) - passed} test(s) failed")
        return 1


if __name__ == '__main__':
    exit(main())
