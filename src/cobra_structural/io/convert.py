"""Provide conversion between structural models and full COBRApy models."""

import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from ..core import Model as StructuralModel


def to_cobrapy_model(
    structural_model: "StructuralModel",
    use_file: bool = False,
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
    use_file : bool, optional
        If True, uses SBML file as intermediate (writes to temp dir then loads).
        If False (default), uses in-memory dict-based conversion which is faster
        but may not preserve all SBML-specific annotations.
    temp_dir : str or pathlib.Path, optional
        Directory for temporary SBML file if use_file=True. Defaults to system
        temp directory (/tmp on Unix). Ignored if use_file=False.

    Returns
    -------
    cobra.Model
        A full COBRApy model with solver, ready for optimization.

    Raises
    ------
    ImportError
        If the original COBRApy package is not installed.

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
    >>> # Convert to full COBRApy
    >>> full_model = to_cobrapy_model(structural_model)
    >>> solution = full_model.optimize()
    >>> print(solution.objective_value)

    Notes
    -----
    - In-memory conversion (use_file=False) is faster and recommended for most cases.
    - File-based conversion (use_file=True) may be more robust for complex models
      with extensive SBML annotations.
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

    if use_file:
        # File-based conversion via SBML
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
    else:
        # In-memory conversion via dict
        from .dict import model_to_dict

        # Convert structural model to dict
        model_dict = model_to_dict(structural_model, sort=False)

        # Use full COBRApy's dict loader
        full_model = cobra.io.model_from_dict(model_dict)

        return full_model
