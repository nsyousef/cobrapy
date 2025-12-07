# COBRApy Structural - Installation and Usage Guide

## Overview

This repository contains **cobra-structural**, a lightweight fork of COBRApy that removes the optlang solver dependency and maintains models purely as data structures. It can coexist with the original **cobra** package and provides bidirectional conversion.

## Installation

### 1. Create a Conda Environment

```bash
conda create -n cobra-structural python=3.10
conda activate cobra-structural
```

### 2. Install cobra-structural (Editable)

```bash
cd /path/to/cobrapy
pip install -e .
```

This installs the `cobra-structural` package from this repository.

### 3. Install Original COBRApy (Optional, for Conversion)

To enable conversion from structural models to solver-enabled models:

```bash
pip install cobra
```

After this, you'll have both packages installed:
- `cobra` (from PyPI, with solver) 
- `cobra-structural` (this repo, structural only)

## Key Differences

| Feature | cobra-structural | cobra (original) |
|---------|-----------------|------------------|
| **Import** | `import cobra_structural` | `import cobra` |
| **Solver** | ❌ None | ✅ optlang-based |
| **Optimization** | ❌ Not available | ✅ `model.optimize()` |
| **Model Building** | ✅ Full support | ✅ Full support |
| **I/O (JSON, SBML, MAT)** | ✅ Full support | ✅ Full support |
| **Objective** | `dict` ({Reaction: coef}) | optlang.Objective |
| **Purpose** | Structural manipulation | Analysis & optimization |

## Usage Examples

### Basic Structural Model

```python
from cobra_structural import Model, Reaction, Metabolite

# Build a model
model = Model('my_model')

# Add metabolites
glc = Metabolite('glc_e', compartment='e')
g6p = Metabolite('g6p_c', compartment='c')

# Add reactions
ex_glc = Reaction('EX_glc')
ex_glc.add_metabolites({glc: -1})
ex_glc.bounds = (-10, 1000)

glk = Reaction('GLK')
glk.add_metabolites({glc: -1, g6p: 1})
glk.bounds = (0, 1000)

model.add_reactions([ex_glc, glk])

# Set objective (dict-based)
model.objective = {glk: 1}
model.objective_direction = 'max'

print(f"Model: {model.id}")
print(f"Reactions: {len(model.reactions)}")
print(f"Objective: {list(model.objective.keys())[0].id}")
```

### Convert to Full COBRApy (for Optimization)

```python
from cobra_structural.io import to_cobrapy_model

# Option 1: In-memory conversion (fast)
full_model = to_cobrapy_model(model, use_file=False)
solution = full_model.optimize()
print(f"Objective value: {solution.objective_value}")

# Option 2: File-based conversion (via SBML)
full_model = to_cobrapy_model(model, use_file=True)
solution = full_model.optimize()
```

### Save and Load Structural Models

```python
from cobra_structural.io import save_json_model, load_json_model

# Save
save_json_model(model, 'my_model.json')

# Load
loaded_model = load_json_model('my_model.json')
```

### Export for mgPipe (MATLAB)

```python
from cobra_structural.io import save_community_mat_model
import scipy.io

# For mgPipe community modeling
save_community_mat_model(
    model, 
    'community_model.mat',
    C=coupling_matrix,  # Optional coupling constraints
    d=coupling_bounds,  # Optional
    dsense=coupling_sense,  # Optional
    ctrs=constraint_names  # Optional
)
```

## Typical Workflow

### Scenario 1: Build Structural Model → Convert → Optimize

```python
import cobra_structural
from cobra_structural.io import to_cobrapy_model

# 1. Build structural model (no solver needed)
struct_model = cobra_structural.Model('test')
# ... add reactions ...
struct_model.objective = {reaction: 1}

# 2. Convert to full COBRApy when needed
full_model = to_cobrapy_model(struct_model)

# 3. Optimize
solution = full_model.optimize()
```

### Scenario 2: Use Only Structural Operations

```python
import cobra_structural
from cobra_structural.io import save_json_model, save_community_mat_model

# Build model
model = cobra_structural.Model('community')
# ... add reactions ...

# Export to file (no conversion needed)
save_json_model(model, 'output.json')
save_community_mat_model(model, 'output.mat')  # For MATLAB/mgPipe
```

### Scenario 3: Two Separate Environments

If you want to avoid installing both packages in the same environment:

**Environment 1 (cobra-structural only):**
```python
from cobra_structural import Model
from cobra_structural.io import save_json_model

model = Model('test')
# ... build model ...
save_json_model(model, 'model.json')
```

**Environment 2 (cobra only):**
```python
import cobra

model = cobra.io.load_json_model('model.json')
solution = model.optimize()
```

## Troubleshooting

### Issue: ImportError when importing cobra

**Problem:** After installing cobra-structural, `import cobra` fails.

**Cause:** Old `cobra.egg-info` directory from before renaming.

**Solution:**
```bash
# Remove old metadata
rm -rf src/cobra.egg-info

# Reinstall
pip uninstall -y cobra cobra-structural
pip install -e .
pip install cobra  # If you want both
```

### Issue: AttributeError 'dict' object has no attribute 'direction'

**Problem:** SBML export fails with objective direction error.

**Cause:** You're using an outdated version that doesn't handle dict objectives.

**Solution:** Pull the latest changes - SBML writer now handles both dict and solver-based objectives.

## Package Structure

```
src/cobra_structural/
├── __init__.py
├── core/
│   ├── model.py          # Model class (no solver)
│   ├── reaction.py       # Reaction class
│   ├── metabolite.py     # Metabolite class
│   └── ...
├── io/
│   ├── convert.py        # to_cobrapy_model()
│   ├── json.py           # JSON I/O
│   ├── sbml.py           # SBML I/O
│   ├── mat.py            # MATLAB/mgPipe export
│   └── ...
├── util/
│   └── solver.py         # Handles dict objectives
└── ...
```

## Testing

Run the comprehensive test:

```bash
python test_both_packages.py
```

This verifies:
- Both packages are installed
- Structural model building
- In-memory conversion
- File-based conversion
- I/O operations

## Important Notes

1. **No optimization in structural models:** You must convert to full COBRApy to run `optimize()`.

2. **Objective format:** Structural models use `{Reaction: coefficient}` dict format, not optlang.Objective.

3. **Package names:**
   - PyPI package: `cobra-structural` (with hyphen)
   - Import name: `cobra_structural` (with underscore)
   - Original: `cobra` (both package and import)

4. **Solver dependency:** cobra-structural has no optlang dependency. The original cobra package requires optlang.

## Questions?

- Check test scripts: `test_both_packages.py`, `test_conversion.py`
- Review function docstrings: `help(to_cobrapy_model)`
- See `.github/copilot-instructions.md` for architecture details
