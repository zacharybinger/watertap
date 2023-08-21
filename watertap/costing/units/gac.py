#################################################################################
# WaterTAP Copyright (c) 2020-2023, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory, Oak Ridge National Laboratory,
# National Renewable Energy Laboratory, and National Energy Technology
# Laboratory (subject to receipt of any required approvals from the U.S. Dept.
# of Energy). All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and license
# information, respectively. These files are also available online at the URL
# "https://github.com/watertap-org/watertap/"
#################################################################################

import pyomo.environ as pyo
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.misc import StrEnum
from idaes.core.util.math import smooth_min
from ..util import (
    register_costing_parameter_block,
    make_capital_cost_var,
    make_fixed_operating_cost_var,
)


class ContactorType(StrEnum):
    pressure = "pressure"
    gravity = "gravity"


def build_gac_cost_param_block(blk):

    # ---------------------------------------------------------------------
    # design options

    blk.num_contactors_op = pyo.Var(
        initialize=1,
        units=pyo.units.dimensionless,
        doc="number of GAC contactors in operation in parallel",
    )
    blk.num_contactors_redundant = pyo.Var(
        initialize=1,
        units=pyo.units.dimensionless,
        doc="number of off-line redundant GAC contactors in parallel",
    )
    blk.regen_frac = pyo.Var(
        initialize=0.70,
        units=pyo.units.dimensionless,
        doc="fraction of spent GAC adsorbent that can be regenerated for reuse",
    )

    # ---------------------------------------------------------------------
    # correlation reference points

    blk.bed_mass_max_ref = pyo.Var(
        initialize=18143.7,
        units=pyo.units.kg,
        doc="reference maximum value of GAC mass needed for initial charge where "
        "economy of scale no longer discounts the unit price",
    )

    # ---------------------------------------------------------------------
    # correlation parameter data

    # dummy data is used to initialize, fixed in cost_gac based on the ContactorType
    # USD_2020 embedded in equation
    contactor_cost_coeff_dummy = {0: 10000, 1: 1000, 2: -10, 3: 0.1}
    blk.contactor_cost_coeff = pyo.Var(
        contactor_cost_coeff_dummy.keys(),
        initialize=contactor_cost_coeff_dummy,
        units=pyo.units.dimensionless,
        doc="contactor polynomial cost coefficients",
    )
    # USD_2020 * kg**-1 embedded in equation adsorbent_unit_cost_constraint
    adsorbent_unit_cost_coeff_dummy = {0: 1, 1: -1e-5}
    blk.adsorbent_unit_cost_coeff = pyo.Var(
        adsorbent_unit_cost_coeff_dummy.keys(),
        initialize=adsorbent_unit_cost_coeff_dummy,
        units=pyo.units.dimensionless,
        doc="GAC adsorbent cost exponential function parameters",
    )
    # USD_2020 embedded in equation other_process_cost_constraint
    other_cost_param_dummy = {0: 10000, 1: 0.1}
    blk.other_cost_param = pyo.Var(
        other_cost_param_dummy.keys(),
        initialize=other_cost_param_dummy,
        units=pyo.units.dimensionless,
        doc="other process cost power law parameters",
    )
    blk.regen_unit_cost = pyo.Var(
        initialize=4.28352,
        units=pyo.units.USD_2020 * pyo.units.kg**-1,
        doc="unit cost to regenerate spent GAC adsorbent by an offsite regeneration facility",
    )
    blk.makeup_unit_cost = pyo.Var(
        initialize=4.58223,
        units=pyo.units.USD_2020 * pyo.units.kg**-1,
        doc="unit cost to makeup spent GAC adsorbent with fresh adsorbent",
    )
    # kW embedded in equation energy_consumption_constraint
    energy_consumption_coeff_dummy = {0: 1e-3, 1: 1e-3, 2: 0}
    blk.energy_consumption_coeff = pyo.Var(
        energy_consumption_coeff_dummy.keys(),
        initialize=energy_consumption_coeff_dummy,
        units=pyo.units.dimensionless,
        doc="energy consumption polynomial coefficients",
    )


@register_costing_parameter_block(
    build_rule=build_gac_cost_param_block,
    parameter_block_name="gac",
)
def cost_gac(blk, contactor_type=ContactorType.pressure):
    """
    3 equation capital cost estimation for GAC systems with: (i), contactor/pressure vessel cost by polynomial as a
    function of individual contactor volume; (ii), initial charge of GAC adsorbent cost by exponential as a function of
    required mass of GAC adsorbent; and (iii), other process costs (vessels, pipes, instrumentation, and controls)
    calculated by power law as a function of total contactor(s) volume. Operating costs calculated as the required
    makeup and regeneration of GAC adsorbent. Energy consumption is estimated from that required for booster, backwash,
    and residual pumps as a function of total contactor(s) volume.

    Args:
        contactor_type: ContactorType Enum indicating whether to cost based on steel pressure vessels or concrete,
            default = ContactorType.pressure
    """

    # ---------------------------------------------------------------------
    # with ContactorType not assigned in build_gac_cost_param_block and blk.costing_package.gac variables
    # fixed when register_costing_parameter_block, refix parameters based on contactor type here

    # costing data parameters based on contactor type
    adsorbent_unit_cost_coeff_data = {0: 4.58342, 1: -1.25311e-5}
    if contactor_type == ContactorType.pressure:
        contactor_cost_coeff_data = {0: 10010.9, 1: 2204.95, 2: -15.9378, 3: 0.110592}
        other_cost_param_data = {0: 16660.7, 1: 0.552207}
        energy_consumption_coeff_data = {0: 8.09926e-4, 1: 8.70577e-4, 2: 0}
    elif contactor_type == ContactorType.gravity:
        contactor_cost_coeff_data = {0: 75131.3, 1: 735.550, 2: -1.01827, 3: 0.000000}
        other_cost_param_data = {0: 38846.9, 1: 0.490571}
        energy_consumption_coeff_data = {0: 0.123782, 1: 0.132403, 2: -1.41512e-5}
    else:
        raise ConfigurationError(
            f"{blk.unit_model.name} received invalid argument for contactor_type:"
            f" {contactor_type}. Argument must be a member of the ContactorType Enum."
        )

    # iterable matching coeff_data to vars
    gac_cost = blk.costing_package.gac
    cost_params = (
        (adsorbent_unit_cost_coeff_data, gac_cost.adsorbent_unit_cost_coeff),
        (contactor_cost_coeff_data, gac_cost.contactor_cost_coeff),
        (other_cost_param_data, gac_cost.other_cost_param),
        (energy_consumption_coeff_data, gac_cost.energy_consumption_coeff),
    )

    # refix variables to appropriate costing parameters
    for indexed_data, indexed_var in cost_params:
        for index, var in indexed_var.items():
            var.fix(indexed_data[index])

    # ---------------------------------------------------------------------

    make_capital_cost_var(blk)
    blk.contactor_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Unit contactor(s) capital cost",
    )
    blk.bed_mass_gac_ref = pyo.Var(
        initialize=4,
        domain=pyo.NonNegativeReals,
        units=pyo.units.kg,
        doc="Reference value of GAC mass needed for initial charge where "
        "economy of scale no longer discounts the unit price",
    )
    blk.adsorbent_unit_cost = pyo.Var(
        initialize=2,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency * pyo.units.kg**-1,
        doc="GAC adsorbent cost per unit mass",
    )
    blk.adsorbent_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Unit adsorbent capital cost",
    )
    blk.other_process_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency,
        doc="Unit other process capital cost",
    )
    blk.energy_consumption = pyo.Var(
        initialize=100,
        domain=pyo.NonNegativeReals,
        units=pyo.units.kW,
        doc="Approximate GAC system energy consumption",
    )

    # intermediate variables to shorten constraint code
    num_contactors = (
        blk.costing_package.gac.num_contactors_op
        + blk.costing_package.gac.num_contactors_redundant
    )
    unit_contactor_volume = (
        blk.unit_model.bed_volume / blk.costing_package.gac.num_contactors_op
    )
    total_bed_volume = num_contactors * unit_contactor_volume

    blk.contactor_cost_constraint = pyo.Constraint(
        expr=blk.contactor_cost
        == num_contactors
        * pyo.units.convert(
            (
                blk.costing_package.gac.contactor_cost_coeff[3]
                * (pyo.units.m**3) ** -3
                * unit_contactor_volume**3
                + blk.costing_package.gac.contactor_cost_coeff[2]
                * (pyo.units.m**3) ** -2
                * unit_contactor_volume**2
                + blk.costing_package.gac.contactor_cost_coeff[1]
                * (pyo.units.m**3) ** -1
                * unit_contactor_volume**1
                + blk.costing_package.gac.contactor_cost_coeff[0]
            )
            * pyo.units.USD_2020,
            to_units=blk.costing_package.base_currency,
        )
    )

    blk.bed_mass_gac_ref_constraint = pyo.Constraint(
        expr=blk.bed_mass_gac_ref
        == smooth_min(
            blk.costing_package.gac.bed_mass_max_ref / pyo.units.kg,
            pyo.units.convert(blk.unit_model.bed_mass_gac, to_units=pyo.units.kg)
            / pyo.units.kg,
        )
        * pyo.units.kg
    )
    blk.adsorbent_unit_cost_constraint = pyo.Constraint(
        expr=blk.adsorbent_unit_cost
        == pyo.units.convert(
            blk.costing_package.gac.adsorbent_unit_cost_coeff[0]
            * pyo.exp(
                blk.bed_mass_gac_ref
                * pyo.units.kg**-1
                * blk.costing_package.gac.adsorbent_unit_cost_coeff[1]
            )
            * pyo.units.USD_2020
            * pyo.units.kg**-1,
            to_units=blk.costing_package.base_currency * pyo.units.kg**-1,
        )
    )
    blk.adsorbent_cost_constraint = pyo.Constraint(
        expr=blk.adsorbent_cost == blk.adsorbent_unit_cost * blk.unit_model.bed_mass_gac
    )

    blk.other_process_cost_constraint = pyo.Constraint(
        expr=blk.other_process_cost
        == pyo.units.convert(
            (
                blk.costing_package.gac.other_cost_param[0]
                * (total_bed_volume * pyo.units.m**-3)
                ** blk.costing_package.gac.other_cost_param[1]
            )
            * pyo.units.USD_2020,
            to_units=blk.costing_package.base_currency,
        )
    )

    blk.capital_cost_constraint = pyo.Constraint(
        expr=blk.capital_cost
        == blk.contactor_cost + blk.adsorbent_cost + blk.other_process_cost
    )

    make_fixed_operating_cost_var(blk)
    blk.gac_regen_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency / blk.costing_package.base_period,
        doc="Cost to regenerate spent GAC adsorbent by an offsite regeneration facility",
    )
    blk.gac_makeup_cost = pyo.Var(
        initialize=1e5,
        domain=pyo.NonNegativeReals,
        units=blk.costing_package.base_currency / blk.costing_package.base_period,
        doc="Cost to makeup spent GAC adsorbent with fresh adsorbent",
    )

    blk.gac_regen_cost_constraint = pyo.Constraint(
        expr=blk.gac_regen_cost
        == pyo.units.convert(
            (
                blk.costing_package.gac.regen_unit_cost
                * (blk.costing_package.gac.regen_frac * blk.unit_model.gac_usage_rate)
            ),
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )
    blk.gac_makeup_cost_constraint = pyo.Constraint(
        expr=blk.gac_makeup_cost
        == pyo.units.convert(
            (
                blk.costing_package.gac.makeup_unit_cost
                * (
                    (1 - blk.costing_package.gac.regen_frac)
                    * blk.unit_model.gac_usage_rate
                )
            ),
            to_units=blk.costing_package.base_currency
            / blk.costing_package.base_period,
        )
    )
    blk.fixed_operating_cost_constraint = pyo.Constraint(
        expr=blk.fixed_operating_cost == blk.gac_regen_cost + blk.gac_makeup_cost
    )

    blk.energy_consumption_constraint = pyo.Constraint(
        expr=pyo.units.convert(blk.energy_consumption, to_units=pyo.units.kW)
        == pyo.units.kW
        * (
            blk.costing_package.gac.energy_consumption_coeff[2]
            * total_bed_volume**2
            * (pyo.units.m**3) ** -2
            + blk.costing_package.gac.energy_consumption_coeff[1]
            * total_bed_volume
            * (pyo.units.m**3) ** -1
            + blk.costing_package.gac.energy_consumption_coeff[0]
        )
    )

    blk.costing_package.cost_flow(
        pyo.units.convert(blk.energy_consumption, to_units=pyo.units.kW),
        "electricity",
    )
