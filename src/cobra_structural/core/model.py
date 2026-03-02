"""Define the Model class."""

import logging
from copy import copy, deepcopy
from functools import partial
from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Union
from warnings import warn

from ..medium import find_boundary_types, find_external_compartment, sbo_terms
from ..util.context import HistoryManager, get_context
from .configuration import Configuration
from .dictlist import DictList
from .gene import Gene
from .group import Group
from .metabolite import Metabolite
from .object import Object
from .reaction import Reaction


if TYPE_CHECKING:
    import pandas as pd

    from cobra import Solution
    from cobra.summary import ModelSummary

logger = logging.getLogger(__name__)
configuration = Configuration()


class Model(Object):
    """Class representation for a cobra model.

    Parameters
    ----------
    id_or_model: str or Model, optional
        String to use as model id, or actual model to base new model one.
        If string, it is used as id. If model, a new model object is
        instantiated with the same properties as the original model (default None).
    name: str, optional
        Human readable string to be model description (default None).

    Attributes
    ----------
    reactions : DictList
        A DictList where the key is the reaction identifier and the value a
        Reaction
    metabolites : DictList
        A DictList where the key is the metabolite identifier and the value a
        Metabolite
    genes : DictList
        A DictList where the key is the gene identifier and the value a
        Gene
    groups : DictList
        A DictList where the key is the group identifier and the value a
        Group
    """

    def __init__(
        self, id_or_model: Union[str, "Model", None] = None, name: Optional[str] = None
    ) -> None:
        """Initialize the Model.

        This creates a structural-only model without an optlang solver.
        Model structure is tracked via reactions, metabolites, and genes.
        """
        if isinstance(id_or_model, Model):
            Object.__init__(self, name=name)
            self.__setstate__(id_or_model.__dict__)
        else:
            Object.__init__(self, id_or_model, name=name)
            self.genes = DictList()
            self.reactions = DictList()  # A list of cobra.Reactions
            self.metabolites = DictList()  # A list of cobra.Metabolites
            self.groups = DictList()  # A list of cobra.Groups
            self._compartments = {}
            self._contexts = []

            # Store objective as simple dict: {reaction: coefficient}
            self._objective = {}
            self._objective_direction = "max"

    def __setstate__(self, state: Dict) -> None:
        """Make sure all cobra.Objects in the model point to the model.

        Parameters
        ----------
        state: dict
        """
        self.__dict__.update(state)
        for y in ["reactions", "genes", "metabolites"]:
            for x in getattr(self, y):
                x._model = self
        if not hasattr(self, "name"):
            self.name = None

    def __getstate__(self) -> Dict:
        """Get state for serialization.

        Ensures that the context stack is cleared prior to serialization,
        since partial functions cannot be pickled reliably.

        Returns
        -------
        odict: Dict
            A dictionary of state, based on self.__dict__.
        """
        odict = self.__dict__.copy()
        odict["_contexts"] = []
        return odict

    @property
    def objective(self) -> Dict[Reaction, float]:
        """Get the objective coefficients.

        Returns
        -------
        dict
            Dictionary mapping reactions to their objective coefficients.

        Examples
        --------
        >>> model.objective = {model.reactions.PGI: 1.0}
        >>> model.objective
        {<Reaction PGI at 0x...>: 1.0}
        """
        return self._objective

    @objective.setter
    def objective(
        self, value: Union[Dict[Reaction, float], str, Reaction, List[Reaction]]
    ) -> None:
        """Set objective coefficients.

        Parameters
        ----------
        value : dict or str or Reaction or list
            If dict: {reaction: coefficient} pairs
            If str: reaction ID to maximize (coeff=1.0)
            If Reaction: reaction object to maximize
            If list: reactions to maximize equally

        Raises
        ------
        ValueError
            If reaction is not in the model.
        TypeError
            If value type is not supported.
        """
        if isinstance(value, dict):
            for rxn in value.keys():
                if rxn not in self.reactions:
                    raise ValueError(f"Reaction {rxn.id} not in model")
            self._objective = value
        elif isinstance(value, str):
            rxn = self.reactions.get_by_id(value)
            self._objective = {rxn: 1.0}
        elif isinstance(value, Reaction):
            if value not in self.reactions:
                raise ValueError(f"Reaction {value.id} not in model")
            self._objective = {value: 1.0}
        elif isinstance(value, list):
            obj_dict = {}
            for item in value:
                if isinstance(item, str):
                    rxn = self.reactions.get_by_id(item)
                else:
                    rxn = item
                obj_dict[rxn] = 1.0
            self._objective = obj_dict
        else:
            raise TypeError(f"Cannot set objective from {type(value)}")

    @property
    def objective_direction(self) -> str:
        """Get objective direction.

        Returns
        -------
        str
            "max" or "min"
        """
        return self._objective_direction

    @objective_direction.setter
    def objective_direction(self, value: str) -> None:
        """Set objective direction.

        Parameters
        ----------
        value : {"max", "min"}
            Direction to maximize or minimize.

        Raises
        ------
        ValueError
            If direction is not 'max' or 'min'.
        """
        value = value.lower()
        if value.startswith("max"):
            self._objective_direction = "max"
        elif value.startswith("min"):
            self._objective_direction = "min"
        else:
            raise ValueError(f"Unknown objective direction '{value}'")

    @property
    def compartments(self) -> Dict:
        """Return all metabolites' compartments.

        Returns
        -------
        dict
            A dictionary of metabolite compartments, where the keys are the short
            version (one letter version) of the compartments, and the values are the
            full names (if they exist).
        """
        return {
            met.compartment: self._compartments.get(met.compartment, "")
            for met in self.metabolites
            if met.compartment is not None
        }

    @compartments.setter
    def compartments(self, value: Dict) -> None:
        """Get or set the dictionary of current compartment descriptions.

        Assigning a dictionary to this property updates the model's
        dictionary of compartment descriptions with the new values.

        Parameters
        ----------
        value : dict
            Dictionary mapping compartments abbreviations to full names.

        Examples
        --------
        >>> from cobra.io import load_model
        >>> model = load_model("textbook")
        >>> model.compartments = {'c': 'the cytosol'}
        >>> model.compartments
        {'c': 'the cytosol', 'e': 'extracellular'}
        """
        self._compartments.update(value)

    @property
    def medium(self) -> Dict[str, float]:
        """Get the constraints on the model exchanges.

        `model.medium` returns a dictionary of the bounds for each of the
        boundary reactions, in the form of `{rxn_id: bound}`, where `bound`
        specifies the absolute value of the bound in direction of metabolite
        creation (i.e., lower_bound for `met <--`, upper_bound for `met -->`)

        Returns
        -------
        Dict[str, float]
            A dictionary with rxn.id (str) as key, bound (float) as value.
        """

        def is_active(reaction: Reaction) -> bool:
            """Determine if boundary reaction permits flux towards creating metabolites.

            Parameters
            ----------
            reaction: cobra.Reaction

            Returns
            -------
            bool
                True if reaction produces metaoblites and has upper_bound above 0
                or if reaction consumes metabolites and has lower_bound below 0 (so
                could be reversed).
            """
            return (bool(reaction.products) and (reaction.upper_bound > 0)) or (
                bool(reaction.reactants) and (reaction.lower_bound < 0)
            )

        def get_active_bound(reaction: Reaction) -> float:
            """For an active boundary reaction, return the relevant bound.

            Parameters
            ----------
            reaction: cobra.Reaction

            Returns
            -------
            float:
                upper or minus lower bound, depenending if the reaction produces or
                consumes metaoblties.
            """
            if reaction.reactants:
                return -reaction.lower_bound
            elif reaction.products:
                return reaction.upper_bound

        return {
            rxn.id: get_active_bound(rxn) for rxn in self.exchanges if is_active(rxn)
        }

    @medium.setter
    def medium(self, medium: Dict[str, float]) -> None:
        """Set the constraints on the model exchanges.

        `model.medium` returns a dictionary of the bounds for each of the
        boundary reactions, in the form of `{rxn_id: rxn_bound}`, where `rxn_bound`
        specifies the absolute value of the bound in direction of metabolite
        creation (i.e., lower_bound for `met <--`, upper_bound for `met -->`)

        Parameters
        ----------
        medium: dict
            The medium to initialize. medium should be a dictionary defining
            `{rxn_id: bound}` pairs.
        """

        def set_active_bound(reaction: Reaction, bound: float) -> None:
            """Set active bound.

            Parameters
            ----------
            reaction: cobra.Reaction
                Reaction to set
            bound: float
                Value to set bound to. The bound is reversed and set as lower bound
                if reaction has reactants (metabolites that are consumed). If reaction
                has reactants, it seems the upper bound won't be set.
            """
            if reaction.reactants:
                reaction.lower_bound = -bound
            elif reaction.products:
                reaction.upper_bound = bound

        # Set the given media bounds
        media_rxns = []
        exchange_rxns = frozenset(self.exchanges)
        for rxn_id, rxn_bound in medium.items():
            rxn = self.reactions.get_by_id(rxn_id)
            if rxn not in exchange_rxns:
                logger.warning(
                    f"{rxn.id} does not seem to be an an exchange reaction. "
                    f"Applying bounds anyway."
                )
            media_rxns.append(rxn)
            # noinspection PyTypeChecker
            set_active_bound(rxn, rxn_bound)

        frozen_media_rxns = frozenset(media_rxns)

        # Turn off reactions not present in media
        for rxn in exchange_rxns - frozen_media_rxns:
            is_export = rxn.reactants and not rxn.products
            set_active_bound(
                rxn, min(0.0, -rxn.lower_bound if is_export else rxn.upper_bound)
            )

    def copy(self) -> "Model":
        """Provide a partial 'deepcopy' of the Model.

        All the Metabolite, Gene, and Reaction objects are created anew but
        in a faster fashion than deepcopy.

        Returns
        -------
        cobra.Model: new model copy

        Notes
        -----
        This model is structural-only and does not include solver state.
        """
        new = self.__class__()
        do_not_copy_by_ref = {
            "metabolites",
            "reactions",
            "genes",
            "notes",
            "annotation",
            "groups",
        }
        for attr in self.__dict__:
            if attr not in do_not_copy_by_ref:
                new.__dict__[attr] = self.__dict__[attr]
        new.notes = deepcopy(self.notes)
        new.annotation = deepcopy(self.annotation)

        new.metabolites = DictList()
        do_not_copy_by_ref = {"_reaction", "_model"}
        for metabolite in self.metabolites:
            new_met = metabolite.__class__()
            for attr, value in metabolite.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_met.__dict__[attr] = (
                        copy(value) if attr == "formula" else value
                    )
            new_met._model = new
            new.metabolites.append(new_met)

        new.genes = DictList()
        for gene in self.genes:
            new_gene = gene.__class__(None)
            for attr, value in gene.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_gene.__dict__[attr] = (
                        copy(value) if attr == "formula" else value
                    )
            new_gene._model = new
            new.genes.append(new_gene)

        new.reactions = DictList()
        do_not_copy_by_ref = {"_model", "_metabolites", "_genes"}
        for reaction in self.reactions:
            new_reaction = reaction.__class__()
            for attr, value in reaction.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_reaction.__dict__[attr] = copy(value)
            new_reaction._model = new
            new.reactions.append(new_reaction)
            for metabolite, stoic in reaction._metabolites.items():
                new_met = new.metabolites.get_by_id(metabolite.id)
                new_reaction._metabolites[new_met] = stoic
                new_met._reaction.add(new_reaction)
            new_reaction.update_genes_from_gpr()

        new.groups = DictList()
        do_not_copy_by_ref = {"_model", "_members"}
        for group in self.groups:
            new_group: Group = group.__class__(group.id)
            for attr, value in group.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_group.__dict__[attr] = copy(value)
            new_group._model = new
            new.groups.append(new_group)
        for group in self.groups:
            new_group = new.groups.get_by_id(group.id)
            new_objects = []
            for member in group.members:
                if isinstance(member, Metabolite):
                    new_object = new.metabolites.get_by_id(member.id)
                elif isinstance(member, Reaction):
                    new_object = new.reactions.get_by_id(member.id)
                elif isinstance(member, Gene):
                    new_object = new.genes.get_by_id(member.id)
                elif isinstance(member, Group):
                    new_object = new.groups.get_by_id(member.id)
                else:
                    raise TypeError(
                        f"The group member {member!r} is unexpectedly not a "
                        f"metabolite, reaction, gene, nor another group."
                    )
                new_objects.append(new_object)
            new_group.add_members(new_objects)

        new._contexts = []
        return new

    def add_metabolites(self, metabolite_list: Union[List, Metabolite]) -> None:
        """Add new metabolites to a model.

        Will add a list of metabolites to the model object.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        metabolite_list : list or Metabolite.
            A list of `cobra.core.Metabolite` objects. If it isn't an iterable
            container, the metabolite will be placed into a list.

        Notes
        -----
        This model is structural-only and does not maintain solver constraints.
        """
        if not hasattr(metabolite_list, "__iter__"):
            metabolite_list = [metabolite_list]
        if len(metabolite_list) == 0:
            return None

        metabolite_list = [x for x in metabolite_list if x.id not in self.metabolites]

        bad_ids = [
            m for m in metabolite_list if not isinstance(m.id, str) or len(m.id) < 1
        ]
        if len(bad_ids) != 0:
            raise ValueError(f"invalid identifiers in {repr(bad_ids)}")

        for x in metabolite_list:
            x._model = self
        self.metabolites += metabolite_list

        context = get_context(self)
        if context:
            context(partial(self.metabolites.__isub__, metabolite_list))
            for x in metabolite_list:
                context(partial(setattr, x, "_model", None))

    def remove_metabolites(
        self, metabolite_list: Union[List, Metabolite], destructive: bool = False
    ) -> None:
        """Remove a list of metabolites from the the object.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        metabolite_list : list or Metaoblite
            A list of `cobra.core.Metabolite` objects. If it isn't an iterable
            container, the metabolite will be placed into a list.

        destructive : bool, optional
            If False then the metabolite is removed from all
            associated reactions.  If True then all associated
            reactions are removed from the Model (default False).

        Notes
        -----
        This model is structural-only and does not maintain solver constraints.
        """
        if not hasattr(metabolite_list, "__iter__"):
            metabolite_list = [metabolite_list]
        metabolite_list = [x for x in metabolite_list if x.id in self.metabolites]
        for x in metabolite_list:
            x._model = None

            associated_groups = self.get_associated_groups(x)
            for group in associated_groups:
                group.remove_members(x)
        if not destructive:
            reaction_coefficients: Dict[Reaction, Dict[Metabolite, float]] = {}
            for metabolite in metabolite_list:
                for reaction in list(metabolite._reaction):
                    coefficients = reaction_coefficients.setdefault(reaction, {})
                    coefficients[metabolite] = reaction._metabolites[metabolite]
            for reaction, coefficients in reaction_coefficients.items():
                reaction.subtract_metabolites(coefficients)
        else:
            reactions_to_remove = set()
            for metabolite in metabolite_list:
                reactions_to_remove.update(list(metabolite._reaction))
            for reaction in reactions_to_remove:
                reaction.remove_from_model()

        self.metabolites -= metabolite_list

        context = get_context(self)
        if context:
            context(partial(self.metabolites.__iadd__, metabolite_list))
            for x in metabolite_list:
                context(partial(setattr, x, "_model", self))

    def add_boundary(
        self,
        metabolite: Metabolite,
        type: str = "exchange",
        reaction_id: Optional[str] = None,
        lb: Optional[float] = None,
        ub: Optional[float] = None,
        sbo_term: Optional[str] = None,
    ) -> Reaction:
        """
        Add a boundary reaction for a given metabolite.

        There are three different types of pre-defined boundary reactions:
        exchange, demand, and sink reactions.
        An exchange reaction is a reversible, unbalanced reaction that adds
        to or removes an extracellular metabolite from the extracellular
        compartment.
        A demand reaction is an irreversible reaction that consumes an
        intracellular metabolite.
        A sink is similar to an exchange but specifically for intracellular
        metabolites, i.e., a reversible reaction that adds or removes an
        intracellular metabolite.

        If you set the reaction `type` to something else, you must specify the
        desired identifier of the created reaction along with its upper and
        lower bound. The name will be given by the metabolite name and the
        given `type`.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        metabolite : cobra.Metabolite
            Any given metabolite. The compartment is not checked but you are
            encouraged to stick to the definition of exchanges and sinks.
        type : {"exchange", "demand", "sink"}
            Using one of the pre-defined reaction types is easiest. If you
            want to create your own kind of boundary reaction choose
            any other string, e.g., 'my-boundary' (default "exchange").
        reaction_id : str, optional
            The ID of the resulting reaction. This takes precedence over the
            auto-generated identifiers but beware that it might make boundary
            reactions harder to identify afterwards when using `model.boundary`
            or specifically `model.exchanges` etc. (default None).
        lb : float, optional
            The lower bound of the resulting reaction (default None).
        ub : float, optional
            The upper bound of the resulting reaction (default None).
        sbo_term : str, optional
            A correct SBO term is set for the available types. If a custom
            type is chosen, a suitable SBO term should also be set (default None).

        Returns
        -------
        cobra.Reaction
            The created boundary reaction.

        Examples
        --------
        >>> from cobra.io load_model
        >>> model = load_model("textbook")
        >>> demand = model.add_boundary(model.metabolites.atp_c, type="demand")
        >>> demand.id
        'DM_atp_c'
        >>> demand.name
        'ATP demand'
        >>> demand.bounds
        (0, 1000.0)
        >>> demand.build_reaction_string()
        'atp_c --> '

        """
        ub = configuration.upper_bound if ub is None else ub
        lb = configuration.lower_bound if lb is None else lb
        types = {
            "exchange": ("EX", lb, ub, sbo_terms["exchange"]),
            "demand": ("DM", 0, ub, sbo_terms["demand"]),
            "sink": ("SK", lb, ub, sbo_terms["sink"]),
        }
        if type == "exchange":
            external = find_external_compartment(self)
            if metabolite.compartment != external:
                raise ValueError(
                    f"The metabolite is not an external metabolite (compartment is "
                    f"`{metabolite.compartment}` but should be `{external}`). "
                    f"Did you mean to add a demand or sink? If not, either change"
                    f" its compartment or rename the model compartments to fix this."
                )
        if type in types:
            prefix, lb, ub, default_term = types[type]
            if reaction_id is None:
                reaction_id = f"{prefix}_{metabolite.id}"
            if sbo_term is None:
                sbo_term = default_term
        if reaction_id is None:
            raise ValueError(
                "Custom types of boundary reactions require a custom "
                "identifier. Please set the `reaction_id`."
            )
        if reaction_id in self.reactions:
            # It already exists so just retrieve it.
            logger.info(f"Boundary reaction '{reaction_id}' already exists.")
            return self.reactions.get_by_id(reaction_id)
        name = f"{metabolite.name} {type}"
        rxn = Reaction(id=reaction_id, name=name, lower_bound=lb, upper_bound=ub)
        rxn.add_metabolites({metabolite: -1})
        if sbo_term:
            rxn.annotation["sbo"] = sbo_term
        self.add_reactions([rxn])
        return rxn

    def add_reactions(self, reaction_list: Iterable[Reaction]) -> None:
        """Add reactions to the model.

        Reactions with identifiers identical to a reaction already in the
        model are ignored.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        reaction_list : list
            A list of `cobra.Reaction` objects

        Notes
        -----
        This model is structural-only and does not maintain solver variables.
        """

        def existing_filter(rxn: Reaction) -> bool:
            """Check if the reaction does not exists in the model.

            Parameters
            ----------
            rxn: cobra.Reaction

            Returns
            -------
            bool
                False if reaction exists, True if it doesn't.
                If the reaction exists, will log a warning.
            """
            if rxn.id in self.reactions:
                logger.warning(f"Ignoring reaction '{rxn.id}' since it already exists.")
                return False
            return True

        pruned = DictList(filter(existing_filter, reaction_list))

        context = get_context(self)

        for reaction in pruned:
            reaction._model = self
            if context:
                context(partial(setattr, reaction, "_model", None))

            for metabolite in list(reaction.metabolites):
                if metabolite not in self.metabolites:
                    self.add_metabolites(metabolite)
                else:
                    stoichiometry = reaction._metabolites.pop(metabolite)
                    model_metabolite = self.metabolites.get_by_id(metabolite.id)
                    reaction._metabolites[model_metabolite] = stoichiometry
                    model_metabolite._reaction.add(reaction)
                    if context:
                        context(partial(model_metabolite._reaction.remove, reaction))

            reaction.update_genes_from_gpr()

        self.reactions += pruned

        if context:
            context(partial(self.reactions.__isub__, pruned))

    def remove_reactions(
        self,
        reactions: Union[str, Reaction, List[Union[str, Reaction]]],
        remove_orphans: bool = False,
    ) -> None:
        """Remove reactions from the model.

        The change is reverted upon exit when using the model as a context.

        Parameters
        ----------
        reactions : list or reaction or str
            A list with reactions (`cobra.Reaction`), or their id's, to remove.
            Reaction will be placed in a list. Str will be placed in a list and used to
            find the reaction in the model.
        remove_orphans : bool, optional
            Remove orphaned genes and metabolites from the model as
            well (default False).

        Notes
        -----
        This model is structural-only and does not maintain solver variables.
        """
        if isinstance(reactions, (str, Reaction)) or not hasattr(reactions, "__iter__"):
            warn("need to pass in a list")
            reactions = [reactions]

        context = get_context(self)
        normalized: List[Reaction] = []
        for reaction in reactions:
            reaction_id = reaction.id if hasattr(reaction, "id") else reaction
            try:
                normalized.append(self.reactions.get_by_id(reaction_id))
            except KeyError:
                warn(f"{reaction} not in {self}")

        if not normalized:
            return None

        unique_reactions: List[Reaction] = []
        seen_ids = set()
        for reaction in normalized:
            if reaction.id in seen_ids:
                continue
            seen_ids.add(reaction.id)
            unique_reactions.append(reaction)

        metabolites_to_check: Dict[str, Metabolite] = {}
        genes_to_check: Dict[str, Gene] = {}

        for reaction in unique_reactions:
            if context:
                context(partial(setattr, reaction, "_model", self))
                context(partial(self.reactions.add, reaction))

            self.reactions.remove(reaction)
            reaction._model = None

            for met in reaction._metabolites:
                if reaction in met._reaction:
                    met._reaction.remove(reaction)
                    if context:
                        context(partial(met._reaction.add, reaction))
                metabolites_to_check[met.id] = met

            for gene in reaction._genes:
                if reaction in gene._reaction:
                    gene._reaction.remove(reaction)
                    if context:
                        context(partial(gene._reaction.add, reaction))
                genes_to_check[gene.id] = gene

            associated_groups = self.get_associated_groups(reaction)
            for group in associated_groups:
                group.remove_members(reaction)

        if remove_orphans:
            orphan_metabolites = [
                metabolite
                for metabolite in metabolites_to_check.values()
                if len(metabolite._reaction) == 0
            ]
            if orphan_metabolites:
                self.remove_metabolites(orphan_metabolites)

            for gene in list(genes_to_check.values()):
                if len(gene._reaction) == 0 and gene in self.genes:
                    self.genes.remove(gene)
                    if context:
                        context(partial(self.genes.add, gene))

    def add_groups(self, group_list: Union[str, Group, List[Group]]) -> None:
        """Add groups to the model.

        Groups with identifiers identical to a group already in the model are
        ignored.

        If any group contains members that are not in the model, these members
        are added to the model as well. Only metabolites, reactions, and genes
        can have groups.

        Parameters
        ----------
        group_list : list or str or Group
            A list of `cobra.Group` objects to add to the model. Can also be a single
            group or a string representing group id. If the input is not a list, a
            warning is raised.
        """

        def existing_filter(new_group: Group) -> bool:
            """Check if the group does not exist.

            Parameters
            ----------
            new_group: cobra.Group
                Group to check.

            Returns
            -------
            bool
                False if the group already exists, True if it doesn't.
            """
            if new_group.id in self.groups:
                logger.warning(
                    f"Ignoring group '{new_group.id}'" f" since it already exists."
                )
                return False
            return True

        if isinstance(group_list, str) or hasattr(group_list, "id"):
            warn("need to pass in a list")
            group_list = [group_list]

        pruned = DictList(filter(existing_filter, group_list))

        for group in pruned:
            group._model = self
            for member in group.members:
                # If the member is not associated with the model, add it
                if isinstance(member, Metabolite) and member not in self.metabolites:
                    self.add_metabolites([member])
                if isinstance(member, Reaction) and member not in self.reactions:
                    self.add_reactions([member])
                # TODO(midnighter): `add_genes` method does not exist.
                # if isinstance(member, Gene):
                #     if member not in self.genes:
                #         self.add_genes([member])

            self.groups += [group]

    def remove_groups(self, group_list: Union[str, Group, List[Group]]) -> None:
        """Remove groups from the model.

        Members of each group are not removed
        from the model (i.e. metabolites, reactions, and genes in the group
        stay in the model after any groups containing them are removed).

        Parameters
        ----------
        group_list : list or str or Group
            A list of `cobra.Group` objects to remove from the model. Can also be a
            single group or a string representing group id. If the input is not a list,
             a warning is raised.
        """
        if isinstance(group_list, str) or hasattr(group_list, "id"):
            warn("need to pass in a list")
            group_list = [group_list]

        for group in group_list:
            # make sure the group is in the model
            if group.id not in self.groups:
                logger.warning(f"{group!r} not in {self!r}. Ignored.")
            else:
                self.groups.remove(group)
                group._model = None

    def get_associated_groups(
        self, element: Union[Reaction, Gene, Metabolite]
    ) -> List[Group]:
        """Get list of groups for element.

        Returns a list of groups that an element (reaction, metabolite, gene)
        is associated with.

        Parameters
        ----------
        element: `cobra.Reaction`, `cobra.Metabolite`, or `cobra.Gene`

        Returns
        -------
        list of `cobra.Group`
            All groups that the provided object is a member of
        """
        # check whether the element is associated with the model
        return [g for g in self.groups if element in g.members]

    @property
    def boundary(self) -> List[Reaction]:
        """Boundary reactions in the model.

        Reactions that either have no substrate or product.

        Returns
        -------
        list
            A list of reactions that either have no substrate or product and
            only one metabolite overall.
        """
        return [rxn for rxn in self.reactions if rxn.boundary]

    @property
    def exchanges(self) -> List[Reaction]:
        """Exchange reactions in model.

        Reactions that exchange mass with the exterior. Uses annotations
        and heuristics to exclude non-exchanges such as sink reactions.

        Returns
        -------
        list
            A list of reactions that satisfy the conditions for exchange reactions.

        See Also
        --------
        cobra.medium.find_boundary_types
        """
        return find_boundary_types(self, "exchange", None)

    @property
    def demands(self) -> List[Reaction]:
        """Demand reactions in model.

        Irreversible reactions that accumulate or consume a metabolite in
        the inside of the model.

        Returns
        -------
        list
            A list of reactions that are demand reactions (reactions that
            accumulate/consume a metabolite irreversibly).

        See Also
        --------
        cobra.medium.find_boundary_types
        """
        return find_boundary_types(self, "demand", None)

    @property
    def sinks(self) -> List[Reaction]:
        """Sink reactions in model.

        Reversible reactions that accumulate or consume a metabolite in
        the inside of the model.

        Returns
        -------
        list
            A list of reactions that are demand reactions (reactions that
            accumulate/consume a metabolite reversibly).

        See Also
        --------
        cobra.medium.find_boundary_types
        """
        return find_boundary_types(self, "sink", None)

    def repair(
        self, rebuild_index: bool = True, rebuild_relationships: bool = True
    ) -> None:
        """Update all indexes and pointers in a model.

        Parameters
        ----------
        rebuild_index : bool
            rebuild the indices kept in reactions, metabolites and genes
        rebuild_relationships : bool
             reset all associations between genes, metabolites, model and
             then re-add them.
        """
        if rebuild_index:  # DictList indexes
            self.reactions._generate_index()
            self.metabolites._generate_index()
            self.genes._generate_index()
            self.groups._generate_index()
        if rebuild_relationships:
            for met in self.metabolites:
                met._reaction.clear()
            for gene in self.genes:
                gene._reaction.clear()
            for rxn in self.reactions:
                rxn.update_genes_from_gpr()
                for met in rxn._metabolites:
                    met._reaction.add(rxn)

        for dict_list in (self.reactions, self.genes, self.metabolites, self.groups):
            for entity in dict_list:
                entity._model = self

    def summary(
        self,
        solution: Optional["Solution"] = None,
        fva: Union["pd.DataFrame", float, None] = None,
    ) -> "ModelSummary":
        """
        Create a summary of the exchange fluxes of the model.

        Parameters
        ----------
        solution : cobra.Solution, optional
            A previous model solution to use for generating the summary. If
            ``None``, the summary method will generate a parsimonious flux
            distribution (default None).
        fva : pd.DataFrame or float, optional
            Whether or not to include flux variability analysis in the output.
            If given, `fva` should either be a previous FVA solution matching the
            model or a float between 0 and 1 representing the fraction of the
            optimum objective to be searched (default None).

        Returns
        -------
        cobra.ModelSummary

        See Also
        --------
        Reaction.summary
        Metabolite.summary

        """
        from cobra.summary import ModelSummary

        return ModelSummary(model=self, solution=solution, fva=fva)

    def __enter__(self) -> "Model":
        """Record future changes to the model.

        Record all future changes to the model, undoing them when a call to
        __exit__ is received. Creates a new context and adds it to the stack.

        Returns
        -------
        cobra.Model
            Returns the model with context added.
        """
        try:
            self._contexts.append(HistoryManager())
        except AttributeError:
            self._contexts = [HistoryManager()]

        return self

    def __exit__(self, type, value, traceback) -> None:
        """Pop the top context manager and trigger the undo functions."""
        context = self._contexts.pop()
        context.reset()

    def merge(
        self,
        right: "Model",
        prefix_existing: Optional[str] = None,
        inplace: bool = True,
        objective: str = "left",
    ) -> "Model":
        """Merge two models to create a model with the reactions from both models.

        Parameters
        ----------
        right : cobra.Model
            The model to add reactions from
        prefix_existing : string or optional
            Prefix the reaction identifier in the right that already exist
            in the left model with this string (default None).
        inplace : bool
            Add reactions from right directly to left model object.
            Otherwise, create a new model leaving the left model untouched
            (default True).
        objective : {"left", "right", "sum"}
            One of "left", "right" or "sum" for setting the objective of the
            resulting model to that of the corresponding model or the sum of
            both (default "left").

        Returns
        -------
        cobra.Model
            The merged model.

        Notes
        -----
        This model is structural-only and does not merge solver state.
        """
        if inplace:
            new_model = self
        else:
            new_model = self.copy()
            new_model.id = f"{self.id}_{right.id}"
        new_reactions = deepcopy(right.reactions)
        if prefix_existing is not None:
            existing = new_reactions.query(lambda rxn: rxn.id in self.reactions)
            for reaction in existing:
                reaction.id = f"{prefix_existing}{reaction.id}"
        new_model.add_reactions(new_reactions)

        # Set objective based on requested combination
        if objective == "left":
            new_model._objective = dict(self._objective)
        elif objective == "right":
            new_model._objective = dict(right._objective)
        elif objective == "sum":
            combined_obj = {}
            for rxn, coef in self._objective.items():
                combined_obj[rxn] = combined_obj.get(rxn, 0) + coef
            for rxn, coef in right._objective.items():
                if rxn in new_model.reactions:
                    target_rxn = new_model.reactions.get_by_id(rxn.id)
                    combined_obj[target_rxn] = combined_obj.get(target_rxn, 0) + coef
            new_model._objective = combined_obj
        else:
            raise ValueError(f"Invalid objective option: {objective}")

        return new_model

    def _repr_html_(self) -> str:
        """Get HTML represenation of the model.

        Returns
        -------
        str
            Model representation as HTML string.
        """
        objective_str = (
            ", ".join(f"{rxn.id}: {coef}" for rxn, coef in self._objective.items())
            if self._objective
            else "None"
        )
        return f"""
        <table>
            <tr>
                <td><strong>Name</strong></td>
                <td>{self.id}</td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>{f"{id(self):x}"}</td>
            </tr><tr>
                <td><strong>Number of metabolites</strong></td>
                <td>{len(self.metabolites)}</td>
            </tr><tr>
                <td><strong>Number of reactions</strong></td>
                <td>{len(self.reactions)}</td>
            </tr><tr>
                <td><strong>Number of genes</strong></td>
                <td>{len(self.genes)}</td>
            </tr><tr>
                <td><strong>Number of groups</strong></td>
                <td>{len(self.groups)}</td>
            </tr><tr>
                <td><strong>Objective</strong></td>
                <td>{objective_str}</td>
            </tr><tr>
                <td><strong>Objective direction</strong></td>
                <td>{self._objective_direction}</td>
            </tr><tr>
                <td><strong>Compartments</strong></td>
                <td>{", ".join(v if v else k for k, v in
                               self.compartments.items())}</td>
            </tr>
          </table>"""
