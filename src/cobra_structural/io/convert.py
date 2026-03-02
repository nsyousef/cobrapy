"""Provide conversion between structural models and full COBRApy models."""

import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Literal

if TYPE_CHECKING:
    from ..core import Model as StructuralModel


def to_cobrapy_model(
    structural_model: "StructuralModel",
    method: Literal['direct', 'dict', 'file'] = "direct",
    temp_dir: Optional[str] = None,
):
    """Convert a structural model to a full COBRApy model with solver.

    This function enables conversion from a lightweight structural-only model
    to a full-featured COBRApy model that includes solver capabilities for
    optimization and analysis.

    Parameters
    ----------
    structural_model : cobra_structural.Model
        The structural-only model to convert.
    method : {"direct", "dict", "file"}, optional
        Conversion method (default is "direct"):
        - "direct": Fastest method, directly instantiates cobra objects from
          structural model data without intermediate serialization. Recommended
          for most use cases, especially large community models.
        - "dict": Uses model_to_dict() intermediate representation. Fallback for
          compatibility if direct method encounters issues.
        - "file": Uses SBML file as intermediate (slowest but most robust for
          models with extensive SBML-specific annotations).
    temp_dir : str or pathlib.Path, optional
        Directory for temporary SBML file if method="file". Defaults to system
        temp directory (/tmp on Unix). Ignored for other methods.

    Returns
    -------
    cobra.Model
        A full COBRApy model with solver, ready for optimization.

    Raises
    ------
    ImportError
        If the original COBRApy package is not installed.
    ValueError
        If an invalid method is specified.

    Examples
    --------
    >>> from cobra_structural import Model, Reaction, Metabolite
    >>> from cobra_structural.io import to_cobrapy_model
    >>> # Build structural model
    >>> structural_model = Model("my_model")
    >>> met = Metabolite("glc_c", compartment="c")
    >>> rxn = Reaction("GLCt")
    >>> rxn.add_metabolites({met: -1})
    >>> structural_model.add_reactions([rxn])
    >>> structural_model.objective = {rxn: 1.0}
    >>> # Convert to full COBRApy (fast direct method)
    >>> full_model = to_cobrapy_model(structural_model)
    >>> solution = full_model.optimize()
    >>> print(solution.objective_value)

    Notes
    -----
    - Direct conversion (method="direct") is fastest and recommended for most cases,
      especially for large microbiome community models.
    - Dict conversion (method="dict") provides compatibility fallback.
    - File conversion (method="file") may be more robust for complex models
      with extensive SBML annotations but is significantly slower.
    - The returned model requires the original ``cobra`` package to be installed.
    - After conversion, you can use all standard COBRApy optimization methods:
      optimize(), flux_variability_analysis(), etc.
    """
    try:
        import cobra
    except ImportError:
        raise ImportError(
            "The original 'cobra' package must be installed to convert to a full "
            "COBRApy model. Install it with: pip install cobra"
        )

    if method == "direct":
        return _direct_conversion(structural_model, cobra)
    elif method == "dict":
        return _dict_conversion(structural_model, cobra)
    elif method == "file":
        return _file_conversion(structural_model, temp_dir, cobra)
    else:
        raise ValueError(
            f"Invalid method '{method}'. Must be 'direct', 'dict', or 'file'."
        )


def _direct_conversion(structural_model: "StructuralModel", cobra):
    """Convert via direct object instantiation (fastest method).

    This bypasses all serialization layers and directly constructs
    COBRApy objects from structural model data structures.

    Parameters
    ----------
    structural_model : cobra_structural.Model
        The structural-only model to convert.
    cobra : module
        The imported cobra module.

    Returns
    -------
    cobra.Model
        A full COBRApy model with solver.
    """
    # Create empty full model with same metadata
    full_model = cobra.Model(structural_model.id)
    full_model.name = structural_model.name
    if hasattr(structural_model, "compartments") and structural_model.compartments:
        full_model.compartments = structural_model.compartments.copy()

    # Map structural metabolites → full metabolites (bulk create)
    met_map = {}
    full_metabolites = []
    for s_met in structural_model.metabolites:
        f_met = cobra.Metabolite(
            id=s_met.id,
            name=s_met.name,
            compartment=s_met.compartment,
        )
        # Optional attributes
        if hasattr(s_met, "charge") and s_met.charge is not None:
            f_met.charge = s_met.charge
        if hasattr(s_met, "formula") and s_met.formula:
            f_met.formula = s_met.formula
        if hasattr(s_met, "notes") and s_met.notes:
            f_met.notes = s_met.notes.copy()
        if hasattr(s_met, "annotation") and s_met.annotation:
            f_met.annotation = s_met.annotation.copy()

        met_map[s_met.id] = f_met
        full_metabolites.append(f_met)

    full_model.add_metabolites(full_metabolites)

    # Create reactions with stoichiometry + GPR (bulk add)
    full_reactions = []
    for s_rxn in structural_model.reactions:
        f_rxn = cobra.Reaction(
            id=s_rxn.id,
            name=s_rxn.name,
            lower_bound=s_rxn.lower_bound,
            upper_bound=s_rxn.upper_bound,
        )

        # Optional attributes
        if hasattr(s_rxn, "subsystem") and s_rxn.subsystem:
            f_rxn.subsystem = s_rxn.subsystem

        # Map stoichiometry using pre-built metabolite map
        if hasattr(s_rxn, "_metabolites") and s_rxn._metabolites:
            f_rxn.add_metabolites(
                {met_map[s_met.id]: coeff for s_met, coeff in s_rxn._metabolites.items()}
            )

        # Transfer GPR (this auto-creates genes in full model)
        if hasattr(s_rxn, "gene_reaction_rule") and s_rxn.gene_reaction_rule:
            f_rxn.gene_reaction_rule = s_rxn.gene_reaction_rule

        if hasattr(s_rxn, "notes") and s_rxn.notes:
            f_rxn.notes = s_rxn.notes.copy()
        if hasattr(s_rxn, "annotation") and s_rxn.annotation:
            f_rxn.annotation = s_rxn.annotation.copy()

        full_reactions.append(f_rxn)

    full_model.add_reactions(full_reactions)

    # Set objective (direct solver access)
    if hasattr(structural_model, "objective") and structural_model.objective:
        obj_dict = {}
        for s_rxn, coeff in structural_model.objective.items():
            # Get reaction by ID from full model
            f_rxn = full_model.reactions.get_by_id(s_rxn.id)
            obj_dict[f_rxn] = coeff
        full_model.objective = obj_dict

    return full_model


def _dict_conversion(structural_model: "StructuralModel", cobra):
    """Convert via dict intermediate representation.

    Parameters
    ----------
    structural_model : cobra_structural.Model
        The structural-only model to convert.
    cobra : module
        The imported cobra module.

    Returns
    -------
    cobra.Model
        A full COBRApy model with solver.
    """
    from .dict import model_to_dict

    # Convert structural model to dict
    model_dict = model_to_dict(structural_model, sort=False)

    # Use full COBRApy's dict loader
    full_model = cobra.io.model_from_dict(model_dict)

    return full_model


def _file_conversion(structural_model: "StructuralModel", temp_dir: Optional[str], cobra):
    """Convert via SBML file intermediate.

    Parameters
    ----------
    structural_model : cobra_structural.Model
        The structural-only model to convert.
    temp_dir : str or None
        Directory for temporary file.
    cobra : module
        The imported cobra module.

    Returns
    -------
    cobra.Model
        A full COBRApy model with solver.
    """
    from .sbml import write_sbml_model

    if temp_dir is None:
        temp_dir = tempfile.gettempdir()

    temp_dir = Path(temp_dir)
    temp_file = temp_dir / f"{structural_model.id}_temp.xml"

    try:
        # Write structural model to SBML
        write_sbml_model(structural_model, str(temp_file))

        # Load into full COBRApy
        full_model = cobra.io.read_sbml_model(str(temp_file))

        return full_model
    finally:
        # Clean up temp file
        if temp_file.exists():
            temp_file.unlink()
