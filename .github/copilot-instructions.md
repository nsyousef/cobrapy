# COBRApy AI Coding Agent Instructions

## Overview
COBRApy is a constraint-based modeling (COBRA) package for genome-scale metabolic network analysis. It provides tools for flux balance analysis (FBA), flux variability analysis (FVA), gene deletion analyses, and other metabolic modeling methods using mathematical optimization via the optlang solver interface.

## Architecture & Components

### Core Model Design (`src/cobra/core/`)
The model follows an **object-graph pattern** with four main entity types:
- **Model**: Container holding reactions, metabolites, genes, and groups. Uses context managers (`__enter__`/`__exit__`) for reversible modifications.
- **Reaction**: Links metabolites (stoichiometric coefficients) and genes (via GPR strings). Maintains forward/reverse variable pairs in the solver.
- **Metabolite**: Represents chemical species with compartment information. Tracks associated reactions via `_reaction` set.
- **Gene**: Encodes gene-to-reaction relationships. Updated from reaction GPR expressions via `update_genes_from_gpr()`.
- **Group**: Optional associations of reactions/metabolites/genes for functional groupings.

All entities inherit from `Object` (base class with id, name, notes, annotation).

### DictList Data Structure
Custom hybrid list+dictionary (`DictList`) provides **O(1) lookups** by object ID while maintaining order:
- Used for `model.reactions`, `model.metabolites`, `model.genes`, `model.groups`
- Always indexed in memory via internal `_dict` mapping
- Call `_generate_index()` after bulk modifications outside the standard API
- Access patterns: `model.reactions.get_by_id(id)`, `model.reactions.query(lambda rxn: ...)`, `model.reactions[index]`

### Solver Integration (`util/solver.py`)
- **Abstraction layer** over optlang (solver interface abstraction)
- Models map to optlang objects: reactions → forward/reverse variables, metabolites → mass-balance constraints
- Objective setting via `set_objective()` (linear only; use custom constraints for nonlinear)
- Key functions: `linear_reaction_coefficients()`, `add_cons_vars_to_problem()`, `remove_cons_vars_from_problem()`
- Status codes: `OPTIMAL`, `FEASIBLE`, `INFEASIBLE` (accessible via `model.solver.status`)

### Context Management Pattern
Models support Python context managers for **atomic reversible modifications**:
```python
with model:
    model.reactions.remove(rxn)
    # Automatically reverted on exit
```
Implementation: `HistoryManager` in `util/context.py` stores undo operations (LIFO stack).

### Flux Analysis Submodule (`flux_analysis/`)
High-level methods for model analysis:
- **FBA** (`optimize()`): Linear optimization
- **FVA** (`variability.flux_variability_analysis()`): Reaction flux ranges
- **Gene deletion** (`deletion.single_gene_deletion()`)
- **Loopless models** (`loopless.py`): Eliminates thermodynamically invalid cycles

### I/O Module (`io/`)
Unified interface for loading/saving models across multiple formats. All formats use a **common dict intermediate** representation via `model_to_dict()` and `model_from_dict()`:

#### Supported Formats

| Format | Module | Functions | Use Case |
|--------|--------|-----------|----------|
| **SBML** | `sbml.py` | `read_sbml_model()`, `write_sbml_model()` | Standard format (FBC extension), most flexible, largest files |
| **JSON** | `json.py` | `load_json_model()`, `save_json_model()` | Human-readable, easy to version control, fast parsing |
| **YAML** | `yaml.py` | `load_yaml_model()`, `save_yaml_model()` | Human-editable, smallest files, readability over speed |
| **MATLAB** | `mat.py` | `load_matlab_model()`, `save_matlab_model()` | Legacy `.mat` files, requires scipy; maps MATLAB field names to annotations |
| **Python Dict** | `dict.py` | `model_to_dict()`, `model_from_dict()` | Internal representation, building custom serialization |
| **Web Repositories** | `web/` | `load_model()` | BiGG, BioModels API integration for remote model loading |

#### Serialization Architecture

**Model → Dict → Format** (all formats converge on `model_to_dict`):
```python
# Unified pipeline
model = load_json_model("model.json")  # → dict → model
model_to_dict(model, sort=False)       # Custom sorting for reproducibility
save_json_model(model, "out.json")     # model → dict → format
```

Key serialization details:
- **Reactions**: id, name, metabolites (dict with stoich), bounds, gene_reaction_rule, subsystem, notes, annotation
- **Metabolites**: id, name, compartment, charge, formula, notes, annotation  
- **Genes**: id, name, notes, annotation
- **Model**: id, name, compartments (dict), objective_coefficient (deprecated), notes, annotation
- **Annotation**: Stored as dict with namespace keys (e.g., `"annotation": {"sbo": "SBO:0000627"}`)

#### SBML-Specific Details

SBML is the most complex format with three parsing modes:
1. **FBC (Flux Balance Constraints) Package**: Modern standard, most efficient parsing
2. **Fallback Mode**: Older SBML files without FBC extension (slower but functional)
3. **Validation**: Use `validate_sbml_model(filename)` to check for specification violations

Special handling:
- **Notes**: Structured as `<p>key: value</p>` HTML tags → converted to `Object.notes` dict
- **Bounds**: Default bounds stored as named parameters (e.g., `cobra_default_ub`); custom bounds as `Parameter` objects
- **ID Escaping**: Non-alphanumeric characters in IDs escaped as `__ASCII_CODE__` (e.g., `__45__` for `-`)
- **Compartment Descriptions**: Not first-class in cobrapy; mapped to model compartments dict

#### MATLAB Format Mapping

MATLAB field names map to annotation namespaces via `MET_MATLAB_TO_PROVIDERS` and `RXN_MATLAB_TO_PROVIDERS`:
```python
# Example mapping in annotations
"metKEGGID" → {"kegg.compound": "C00001"}
"rxnECNumbers" → {"ec-code": ["1.1.1.1"]}
```
MATLAB also stores confidence scores and references in reaction notes.

#### Loading from Remote Repositories
```python
from cobra.io import load_model
model = load_model("iJO1366")  # BiGG Models
model = load_model("BIOMD0000000393", cache=True)  # BioModels
```
Automatic caching and format detection based on repository.

## Developer Workflows

### Testing
Run tests via `tox` (orchestrates linting, type checking, and pytest):
```bash
tox                           # Full test suite (all Python versions + linting)
tox -e py311                  # Specific Python version
tox -e coverage               # Coverage report
pytest tests/test_core/       # Direct pytest on subset
pytest tests/test_core/test_model.py::TestModel::test_add_reactions -xvs
```

**Test fixtures** in `tests/conftest.py`:
- `model`: Function-scoped textbook model (fast default test model)
- `large_model`: E. coli iJO1366 (session-scoped, reused across tests)
- `empty_model`: Empty model for isolation tests
- Always `.copy()` session-level fixtures to avoid state pollution

**Test data**: Located in `tests/data/` (pickled models, XML files) and `cobra/data/` (shipped test models).

### Code Style & Quality
- **Formatter**: Black (line length 88)
- **Import sorting**: isort (config in `pyproject.toml`)
- **Linting**: flake8
- **Pre-commit hooks**: `.pre-commit-config.yaml` runs formatters automatically

Semantic commit messages required (see CONTRIBUTING.rst): `fix:`, `feat:`, `refactor:`, etc.

### Build & Installation
- Setup: `pip install -e ".[development]"` (editable with dev dependencies)
- Version in `setup.py` (single source of truth)
- Package config in `setup.cfg`
- Distribution via PyPI and conda

## Critical Patterns & Conventions

### Metabolite-Reaction Coupling
Reactions don't directly store metabolites; instead, stoichiometry is maintained in `Reaction._metabolites` (dict: Metabolite → coefficient). **Metabolites track reverse relationships** via `_reaction` set:
```python
rxn.add_metabolites({met: -1})  # Updates met._reaction automatically
met._model = model              # Maintains bidirectional reference
```
Always ensure metabolites are in `model.metabolites` before adding to reactions (automatic in `add_reactions()`).

### GPR String Parsing
Gene-protein-reaction associations encoded as boolean expressions (e.g., `"(gene1 AND gene2) OR gene3"`):
- Parsed by `Gene.update_gpr()` to build Gene-Reaction graph
- Called automatically via `Reaction.update_genes_from_gpr()` during model modifications
- **Do not manually modify gene relations**—modify GPR strings instead

### Model State & Serialization
- Use `copy()` for fast shallow copies (faster than `deepcopy`)
- Context stack (`_contexts`) cleared during pickling (`__getstate__`) to avoid serialization issues
- `repair()` method rebuilds all indices and relationships (use after manual manipulation outside API)

### Boundary Reactions
Pre-defined boundary types accessed via `model.exchanges`, `model.demands`, `model.sinks`:
- **Exchange**: Reversible, external metabolite (e.g., nutrient import/export)
- **Demand**: Irreversible sink for internal metabolites
- **Sink**: Reversible sink for internal metabolites
- Create via `model.add_boundary(metabolite, type="exchange")`

### Solver Abstractions
- Always check `solver.status` after optimization (not just return values)
- Don't assume linear objectives—use `linear_reaction_coefficients()` to extract coefficients
- Custom constraints/variables added via `model.add_cons_vars()` (respects context)
- Solver-specific tolerances set via `model.tolerance = value`

## Common Pitfalls

1. **Modifying collections outside the API** (e.g., `model.reactions.pop(0)`): Use `remove_reactions()` to maintain consistency
2. **Assuming solver status**: Always check `status` enum; don't rely on exception raising alone
3. **Forgetting metabolite compartments**: Boundary reactions must match metabolite compartment (external/internal)
4. **Accessing private attributes** (`_metabolites`, `_model`): Use public APIs unless absolutely necessary
5. **Not reverting context changes**: Incomplete `with model:` blocks can leak state

## Key Files to Reference

| File | Purpose |
|------|---------|
| `src/cobra/core/model.py` | Model class definition, core API |
| `src/cobra/core/reaction.py` | Reaction stoichiometry, GPR, bounds |
| `src/cobra/core/dictlist.py` | Fast lookups by ID |
| `src/cobra/util/context.py` | Reversible modifications pattern |
| `src/cobra/util/solver.py` | Solver abstraction, optimization |
| `tests/conftest.py` | Test fixtures and data loading |
| `.github/CONTRIBUTING.rst` | Contribution guidelines |

## Tips for Productivity

- Use `model.summary()` to inspect exchange fluxes and key reactions quickly
- `model.repair()` after bulk manual operations to rebuild indices
- Always use `model.copy()` for test model isolation (not separate objects)
- Reference test files in `tests/test_core/test_*.py` for API usage examples
- Check solver availability with `from cobra.util.solver import solvers`; not all may be installed
