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

import pytest
from pyomo.environ import (
    ConcreteModel,
    value,
    Constraint,
    Var,
    assert_optimal_termination,
)
from pyomo.util.check_units import assert_units_consistent

from pyomo.network import Port
from idaes.core import (
    FlowsheetBlock,
    MaterialBalanceType,
    MomentumBalanceType,
    StateBlock,
)
from watertap.unit_models.reverse_osmosis_1D import (
    ReverseOsmosis1D,
    ConcentrationPolarizationType,
    MassTransferCoefficient,
    PressureChangeType,
)
import watertap.property_models.NaCl_prop_pack as props

from idaes.core.solvers import get_solver
from idaes.core.util.model_statistics import (
    number_variables,
    number_unused_variables,
    number_total_constraints,
    degrees_of_freedom,
)

from idaes.core.util.testing import initialization_tester
from idaes.core.util.scaling import (
    calculate_scaling_factors,
    unscaled_variables_generator,
    badly_scaled_var_generator,
)

from watertap.core import (
    MembraneChannel1DBlock,
    FrictionFactor,
)
import idaes.logger as idaeslog

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()
# -----------------------------------------------------------------------------


@pytest.mark.unit
def test_config():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis1D(property_package=m.fs.properties)

    assert len(m.fs.unit.config) == 19
    assert not m.fs.unit.config.dynamic
    assert not m.fs.unit.config.has_holdup
    assert m.fs.unit.config.material_balance_type == MaterialBalanceType.useDefault
    assert m.fs.unit.config.momentum_balance_type == MomentumBalanceType.pressureTotal
    assert not m.fs.unit.config.has_pressure_change
    assert m.fs.unit.config.property_package is m.fs.properties
    assert m.fs.unit.config.pressure_change_type is PressureChangeType.fixed_per_stage
    assert (
        m.fs.unit.config.concentration_polarization_type
        is ConcentrationPolarizationType.calculated
    )
    assert (
        m.fs.unit.config.mass_transfer_coefficient is MassTransferCoefficient.calculated
    )
    assert not m.fs.unit.config.has_full_reporting


@pytest.mark.unit
def test_option_has_pressure_change():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis1D(
        property_package=m.fs.properties, has_pressure_change=True
    )

    assert isinstance(m.fs.unit.feed_side.deltaP, Var)
    assert isinstance(m.fs.unit.feed_side.dP_dx, Var)
    assert isinstance(m.fs.unit.deltaP, Var)


@pytest.mark.unit
def test_option_concentration_polarization_type_fixed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=False,
        concentration_polarization_type=ConcentrationPolarizationType.fixed,
        mass_transfer_coefficient=MassTransferCoefficient.none,
    )

    assert (
        m.fs.unit.config.concentration_polarization_type
        == ConcentrationPolarizationType.fixed
    )
    assert isinstance(m.fs.unit.feed_side.cp_modulus, Var)


@pytest.mark.unit
def test_option_concentration_polarization_type_calculated_kf_fixed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=False,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.fixed,
    )

    assert (
        m.fs.unit.config.concentration_polarization_type
        == ConcentrationPolarizationType.calculated
    )
    assert m.fs.unit.config.mass_transfer_coefficient == MassTransferCoefficient.fixed
    assert isinstance(m.fs.unit.feed_side.K, Var)


@pytest.mark.unit
def test_option_concentration_polarization_type_calculated_kf_calculated():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=False,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
    )

    assert (
        m.fs.unit.config.concentration_polarization_type
        == ConcentrationPolarizationType.calculated
    )
    assert (
        m.fs.unit.config.mass_transfer_coefficient == MassTransferCoefficient.calculated
    )
    assert isinstance(m.fs.unit.feed_side.K, Var)
    assert isinstance(m.fs.unit.feed_side.channel_height, Var)
    assert isinstance(m.fs.unit.feed_side.dh, Var)
    assert isinstance(m.fs.unit.feed_side.spacer_porosity, Var)
    assert isinstance(m.fs.unit.feed_side.N_Sc_comp, Var)
    assert isinstance(m.fs.unit.feed_side.N_Sh_comp, Var)
    assert isinstance(m.fs.unit.feed_side.N_Re, Var)


@pytest.mark.unit
def test_option_pressure_change_calculated():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.none,
        mass_transfer_coefficient=MassTransferCoefficient.none,
        pressure_change_type=PressureChangeType.calculated,
    )

    assert (
        m.fs.unit.config.concentration_polarization_type
        == ConcentrationPolarizationType.none
    )
    assert m.fs.unit.config.mass_transfer_coefficient == MassTransferCoefficient.none
    assert m.fs.unit.config.pressure_change_type == PressureChangeType.calculated
    assert isinstance(m.fs.unit.feed_side.deltaP, Var)
    assert isinstance(m.fs.unit.feed_side.dP_dx, Var)
    assert isinstance(m.fs.unit.deltaP, Var)
    assert isinstance(m.fs.unit.feed_side.channel_height, Var)
    assert isinstance(m.fs.unit.feed_side.dh, Var)
    assert isinstance(m.fs.unit.feed_side.spacer_porosity, Var)
    assert isinstance(m.fs.unit.feed_side.N_Re, Var)


@pytest.mark.unit
def test_option_friction_factor_spiral_wound():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = props.NaClParameterBlock()
    m.fs.unit = ReverseOsmosis1D(
        property_package=m.fs.properties,
        has_pressure_change=True,
        concentration_polarization_type=ConcentrationPolarizationType.calculated,
        mass_transfer_coefficient=MassTransferCoefficient.calculated,
        pressure_change_type=PressureChangeType.calculated,
        friction_factor=FrictionFactor.spiral_wound,
    )

    assert m.fs.unit.config.friction_factor == FrictionFactor.spiral_wound
    assert isinstance(m.fs.unit.feed_side.velocity, Var)
    assert isinstance(m.fs.unit.feed_side.eq_friction_factor, Constraint)


class TestReverseOsmosis:
    @pytest.fixture(scope="class")
    def RO_frame(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            pressure_change_type=PressureChangeType.calculated,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=3,
            has_full_reporting=True,
        )

        # fully specify system
        feed_flow_mass = 1000 / 3600
        feed_mass_frac_NaCl = 0.034283
        feed_pressure = 70e5

        feed_temperature = 273.15 + 25
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 1e5
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )

        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.feed_side.N_Re[0, 0].fix(400)
        m.fs.unit.recovery_mass_phase_comp[0, "Liq", "H2O"].fix(0.5)
        m.fs.unit.feed_side.spacer_porosity.fix(0.97)
        m.fs.unit.feed_side.channel_height.fix(0.001)

        return m

    @pytest.mark.unit
    def test_build(self, RO_frame):
        m = RO_frame

        # test ports
        port_lst = ["inlet", "retentate", "permeate"]
        for port_str in port_lst:
            port = getattr(m.fs.unit, port_str)
            assert isinstance(port, Port)
            # number of state variables for NaCl property package
            assert len(port.vars) == 3

        # test feed-side control volume and associated stateblocks
        assert isinstance(m.fs.unit.feed_side, MembraneChannel1DBlock)
        cv_stateblock_lst = ["properties"]
        for sb_str in cv_stateblock_lst:
            sb = getattr(m.fs.unit.feed_side, sb_str)
            assert isinstance(sb, props.NaClStateBlock)

        stateblock_lst = ["permeate_side", "mixed_permeate"]
        for sb_str in stateblock_lst:
            sb = getattr(m.fs.unit, sb_str)
            assert isinstance(sb, StateBlock)
            assert isinstance(sb, props.NaClStateBlock)

        # test statistics
        assert number_variables(m) == 242
        assert number_total_constraints(m) == 204
        assert number_unused_variables(m) == 14

    @pytest.mark.integration
    def test_units(self, RO_frame):
        m = RO_frame
        assert_units_consistent(m.fs.unit)

    @pytest.mark.unit
    def test_dof(self, RO_frame):
        m = RO_frame
        assert degrees_of_freedom(m) == 0

    @pytest.mark.unit
    def test_calculate_scaling(self, RO_frame):
        m = RO_frame

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

    @pytest.mark.component
    def test_initialize(self, RO_frame):
        initialization_tester(RO_frame, outlvl=idaeslog.DEBUG)

    @pytest.mark.component
    def test_var_scaling(self, RO_frame):
        m = RO_frame
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

    @pytest.mark.component
    def test_solve(self, RO_frame):
        m = RO_frame

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

    @pytest.mark.component
    def test_conservation(self, RO_frame):
        m = RO_frame
        b = m.fs.unit
        comp_lst = ["NaCl", "H2O"]

        flow_mass_inlet = sum(
            b.feed_side.properties[0, 0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            b.feed_side.properties[0, 1].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        assert value(flow_mass_inlet) == pytest.approx(1 / 3.6, rel=1e-3)
        assert value(flow_mass_retentate) == pytest.approx(0.1437, rel=1e-3)
        assert value(flow_mass_permeate) == pytest.approx(0.1341, rel=1e-3)

        assert (
            abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
            <= 1e-2
        )

    @pytest.mark.component
    def test_solution(self, RO_frame):
        m = RO_frame
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)
        assert pytest.approx(-1.755e5, rel=1e-3) == value(m.fs.unit.deltaP[0])
        assert pytest.approx(7.992e-3, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        )
        assert pytest.approx(2.114e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        )
        assert pytest.approx(2.486e-3, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        )
        assert pytest.approx(2.593e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        )
        assert pytest.approx(18.131, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "H2O"] * 3.6e3
        )
        assert pytest.approx(8.543, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp_avg[0, "Liq", "NaCl"] * 3.6e6
        )
        assert pytest.approx(0.1341, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(6.3195e-5, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        )
        assert pytest.approx(371.01, rel=1e-3) == value(m.fs.unit.feed_side.N_Re_avg[0])
        assert pytest.approx(107.48, rel=1e-3) == value(
            m.fs.unit.feed_side.K_avg[0, "NaCl"] * 3.6e6
        )
        assert pytest.approx(26.63, rel=1e-3) == value(m.fs.unit.area)

    @pytest.mark.component
    def testReverseOsmosis_basic(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=False,
            concentration_polarization_type=ConcentrationPolarizationType.none,
            mass_transfer_coefficient=MassTransferCoefficient.none,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=3,
        )

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )

        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.length.fix(8)
        m.fs.unit.recovery_vol_phase[0, "Liq"].fix(0.4)

        # test statistics
        assert number_variables(m) == 198
        assert number_total_constraints(m) == 163
        assert number_unused_variables(m) == 20

        # Test units
        assert_units_consistent(m.fs.unit)

        # Check degrees of freedom = 0
        assert degrees_of_freedom(m) == 0

        # Scaling
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        # Test initialization
        initialization_tester(m, outlvl=idaeslog.DEBUG)

        # Test variable scaling
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        # Solve
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Check mass conservation
        b = m.fs.unit
        comp_lst = ["NaCl", "H2O"]

        flow_mass_inlet = sum(
            b.feed_side.properties[0, 0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            b.feed_side.properties[0, 1].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        assert value(flow_mass_inlet) == pytest.approx(1.0, rel=1e-3)
        assert value(flow_mass_retentate) == pytest.approx(0.6103, rel=1e-3)
        assert value(flow_mass_permeate) == pytest.approx(0.3898, rel=1e-3)

        assert (
            abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
            <= 1e-2
        )

        # Test solution
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)
        assert pytest.approx(4.841e-3, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        )
        assert pytest.approx(1.629e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        )
        assert pytest.approx(144.31, rel=1e-3) == value(m.fs.unit.area)
        assert pytest.approx(1.019e-3, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        )
        assert pytest.approx(2.000e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        )
        assert pytest.approx(0.3896, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(2.652e-4, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        )

    @pytest.mark.component
    def testReverseOsmosis_cp_mod_fixed(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=False,
            concentration_polarization_type=ConcentrationPolarizationType.fixed,
            mass_transfer_coefficient=MassTransferCoefficient.none,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=3,
        )

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )

        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.length.fix(8)
        m.fs.unit.recovery_vol_phase[0, "Liq"].fix(0.4)
        m.fs.unit.feed_side.cp_modulus.fix(1.1)

        # test statistics
        assert number_variables(m) == 202
        assert number_total_constraints(m) == 163
        assert number_unused_variables(m) == 21

        assert_units_consistent(m.fs.unit)

        assert degrees_of_freedom(m) == 0

        # Scaling
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )
        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        # Test initialization
        initialization_tester(m, outlvl=idaeslog.DEBUG)
        # Check for poorly scaled variables
        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        # Solve
        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Check mass conservation
        b = m.fs.unit
        comp_lst = ["NaCl", "H2O"]

        flow_mass_inlet = sum(
            b.feed_side.properties[0, 0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            b.feed_side.properties[0, 1].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        assert value(flow_mass_inlet) == pytest.approx(1.0, rel=1e-3)
        assert value(flow_mass_retentate) == pytest.approx(0.6103, rel=1e-3)
        assert value(flow_mass_permeate) == pytest.approx(0.3898, rel=1e-3)

        assert (
            abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
            <= 1e-2
        )

        # Test solution
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)
        assert pytest.approx(2.449e-3, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        )
        assert pytest.approx(1.8639e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        )
        assert pytest.approx(329.708, rel=1e-3) == value(m.fs.unit.area)
        assert pytest.approx(3.600e-4, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        )
        assert pytest.approx(2.052e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        )
        assert pytest.approx(0.3895, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(6.529e-4, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        )

    @pytest.mark.component
    def testReverseOsmosis_cp_calculated_kf_fixed(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=False,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.fixed,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=3,
        )

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )

        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        # m.fs.unit.dP_dx.fix(0)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.length.fix(8)
        m.fs.unit.recovery_vol_phase[0, "Liq"].fix(0.4)
        m.fs.unit.feed_side.K.fix(2e-5)

        # test statistics
        assert number_variables(m) == 206
        assert number_total_constraints(m) == 166
        assert number_unused_variables(m) == 22

        assert_units_consistent(m.fs.unit)

        assert degrees_of_freedom(m) == 0

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )
        calculate_scaling_factors(m)

        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        initialization_tester(m, outlvl=idaeslog.DEBUG)

        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Check mass conservation
        b = m.fs.unit
        comp_lst = ["NaCl", "H2O"]

        flow_mass_inlet = sum(
            b.feed_side.properties[0, 0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            b.feed_side.properties[0, 1].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        assert value(flow_mass_inlet) == pytest.approx(1.0, rel=1e-3)
        assert value(flow_mass_retentate) == pytest.approx(0.6103, rel=1e-3)
        assert value(flow_mass_permeate) == pytest.approx(0.3898, rel=1e-3)

        assert (
            abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
            <= 1e-2
        )

        # Test solution
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)
        assert pytest.approx(2.792e-3, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        )
        assert pytest.approx(1.831e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        )
        assert pytest.approx(239.28, rel=1e-3) == value(m.fs.unit.area)
        assert pytest.approx(7.082e-4, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        )
        assert pytest.approx(2.027e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        )
        assert pytest.approx(0.3895, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(4.645e-4, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        )

    @pytest.mark.component
    def testReverseOsmosis_cp_calculated_kf_calculated(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=False,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=3,
        )

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )

        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.length.fix(8)
        m.fs.unit.recovery_vol_phase[0, "Liq"].fix(0.4)
        m.fs.unit.feed_side.spacer_porosity.fix(0.75)
        m.fs.unit.feed_side.channel_height.fix(0.002)

        # test statistics
        assert number_variables(m) == 229
        assert number_total_constraints(m) == 191
        assert number_unused_variables(m) == 14

        assert_units_consistent(m.fs.unit)

        assert degrees_of_freedom(m) == 0

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        initialization_tester(m, outlvl=idaeslog.DEBUG)

        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Check mass conservation
        b = m.fs.unit
        comp_lst = ["NaCl", "H2O"]

        flow_mass_inlet = sum(
            b.feed_side.properties[0, 0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            b.feed_side.properties[0, 1].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        assert value(flow_mass_inlet) == pytest.approx(1.0, rel=1e-3)
        assert value(flow_mass_retentate) == pytest.approx(0.6103, rel=1e-3)
        assert value(flow_mass_permeate) == pytest.approx(0.3898, rel=1e-3)

        assert (
            abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
            <= 1e-2
        )

        # Test solution
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)
        assert pytest.approx(2.383e-3, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        )
        assert pytest.approx(1.870e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        )
        assert pytest.approx(278.49, rel=1e-3) == value(m.fs.unit.area)
        assert pytest.approx(6.270e-4, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        )
        assert pytest.approx(2.034e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        )
        assert pytest.approx(0.3895, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(5.467e-4, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        )

    @pytest.mark.component
    def test_friction_factor_spiral_wound(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            pressure_change_type=PressureChangeType.calculated,
            friction_factor=FrictionFactor.spiral_wound,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=3,
        )

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        membrane_area = 19
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )

        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.area.fix(membrane_area)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.length.fix(8)
        m.fs.unit.feed_side.spacer_porosity.fix(0.75)
        m.fs.unit.feed_side.channel_height.fix(0.002)

        # test statistics
        assert number_variables(m) == 242
        assert number_total_constraints(m) == 204
        assert number_unused_variables(m) == 14

        assert_units_consistent(m.fs.unit)

        assert degrees_of_freedom(m) == 0

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        initialization_tester(m, outlvl=idaeslog.DEBUG)

        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Check mass conservation
        b = m.fs.unit
        comp_lst = ["NaCl", "H2O"]

        flow_mass_inlet = sum(
            b.feed_side.properties[0, 0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            b.feed_side.properties[0, 1].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        assert value(flow_mass_inlet) == pytest.approx(1.0, rel=1e-3)
        assert value(flow_mass_retentate) == pytest.approx(0.8981, rel=1e-3)
        assert value(flow_mass_permeate) == pytest.approx(0.1011, rel=1e-3)

        assert (
            abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
            <= 1e-2
        )

        # Test solution
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)
        assert pytest.approx(-1.879e5, rel=1e-3) == value(m.fs.unit.deltaP[0])
        assert pytest.approx(5.908e-3, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        )
        assert pytest.approx(1.495e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        )
        assert pytest.approx(4.745e-3, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        )
        assert pytest.approx(1.559e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        )
        assert pytest.approx(0.1011, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(2.902e-5, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        )

    @pytest.mark.component
    def testRO_cp_calculated_kf_calculated_pdrop_fixed_by_dx(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            pressure_change_type=PressureChangeType.fixed_per_unit_length,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=3,
        )

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )

        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.length.fix(8)
        m.fs.unit.recovery_vol_phase[0, "Liq"].fix(0.4)
        m.fs.unit.feed_side.spacer_porosity.fix(0.75)
        m.fs.unit.feed_side.channel_height.fix(0.002)
        m.fs.unit.feed_side.dP_dx.fix(-0.1e5)

        # test statistics
        assert number_variables(m) == 234
        assert number_total_constraints(m) == 192
        assert number_unused_variables(m) == 15

        assert_units_consistent(m.fs.unit)

        assert degrees_of_freedom(m) == 0

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        initialization_tester(m, outlvl=idaeslog.DEBUG)

        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Check mass conservation
        b = m.fs.unit
        comp_lst = ["NaCl", "H2O"]

        flow_mass_inlet = sum(
            b.feed_side.properties[0, 0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            b.feed_side.properties[0, 1].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        assert value(flow_mass_inlet) == pytest.approx(1.0, rel=1e-3)
        assert value(flow_mass_retentate) == pytest.approx(0.6103, rel=1e-3)
        assert value(flow_mass_permeate) == pytest.approx(0.3898, rel=1e-3)

        assert (
            abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
            <= 1e-2
        )

        # Test solution
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)
        assert pytest.approx(-8.000e4, rel=1e-3) == value(m.fs.unit.deltaP[0])
        assert pytest.approx(2.249e-3, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        )
        assert pytest.approx(1.872e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        )
        assert pytest.approx(301.76, rel=1e-3) == value(m.fs.unit.area)
        assert pytest.approx(5.511e-4, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        )
        assert pytest.approx(2.008e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        )
        assert pytest.approx(0.3895, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(5.888e-4, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        )

    @pytest.mark.component
    def testRO_cp_calculated_kf_calculated_pdrop_fixed_by_stage(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)

        m.fs.properties = props.NaClParameterBlock()

        m.fs.unit = ReverseOsmosis1D(
            property_package=m.fs.properties,
            has_pressure_change=True,
            concentration_polarization_type=ConcentrationPolarizationType.calculated,
            mass_transfer_coefficient=MassTransferCoefficient.calculated,
            pressure_change_type=PressureChangeType.fixed_per_stage,
            transformation_scheme="BACKWARD",
            transformation_method="dae.finite_difference",
            finite_elements=3,
        )

        # fully specify system
        feed_flow_mass = 1
        feed_mass_frac_NaCl = 0.035
        feed_pressure = 50e5
        feed_temperature = 273.15 + 25
        A = 4.2e-12
        B = 3.5e-8
        pressure_atmospheric = 101325
        feed_mass_frac_H2O = 1 - feed_mass_frac_NaCl

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "NaCl"].fix(
            feed_flow_mass * feed_mass_frac_NaCl
        )

        m.fs.unit.inlet.flow_mass_phase_comp[0, "Liq", "H2O"].fix(
            feed_flow_mass * feed_mass_frac_H2O
        )

        m.fs.unit.inlet.pressure[0].fix(feed_pressure)
        m.fs.unit.inlet.temperature[0].fix(feed_temperature)
        m.fs.unit.A_comp.fix(A)
        m.fs.unit.B_comp.fix(B)
        m.fs.unit.permeate.pressure[0].fix(pressure_atmospheric)
        m.fs.unit.length.fix(8)
        m.fs.unit.recovery_vol_phase[0, "Liq"].fix(0.4)
        m.fs.unit.feed_side.spacer_porosity.fix(0.75)
        m.fs.unit.feed_side.channel_height.fix(0.002)
        m.fs.unit.deltaP.fix(-62435.6)

        # test statistics
        assert number_variables(m) == 234
        assert number_total_constraints(m) == 195
        assert number_unused_variables(m) == 14

        assert_units_consistent(m.fs.unit)

        assert degrees_of_freedom(m) == 0

        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e1, index=("Liq", "H2O")
        )
        m.fs.properties.set_default_scaling(
            "flow_mass_phase_comp", 1e3, index=("Liq", "NaCl")
        )

        calculate_scaling_factors(m)

        # check that all variables have scaling factors
        unscaled_var_list = list(unscaled_variables_generator(m))
        assert len(unscaled_var_list) == 0

        initialization_tester(m, outlvl=idaeslog.DEBUG)

        badly_scaled_var_lst = list(badly_scaled_var_generator(m))
        assert badly_scaled_var_lst == []

        results = solver.solve(m)

        # Check for optimal solution
        assert_optimal_termination(results)

        # Check mass conservation
        b = m.fs.unit
        comp_lst = ["NaCl", "H2O"]

        flow_mass_inlet = sum(
            b.feed_side.properties[0, 0].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_retentate = sum(
            b.feed_side.properties[0, 1].flow_mass_phase_comp["Liq", j]
            for j in comp_lst
        )
        flow_mass_permeate = sum(
            b.mixed_permeate[0].flow_mass_phase_comp["Liq", j] for j in comp_lst
        )

        assert value(flow_mass_inlet) == pytest.approx(1.0, rel=1e-3)
        assert value(flow_mass_retentate) == pytest.approx(0.6103, rel=1e-3)
        assert value(flow_mass_permeate) == pytest.approx(0.3898, rel=1e-3)

        assert (
            abs(value(flow_mass_inlet - flow_mass_retentate - flow_mass_permeate))
            <= 1e-2
        )

        # Test solution
        x_interface_in = m.fs.unit.feed_side.length_domain.at(2)
        assert pytest.approx(-6.2436e4, rel=1e-3) == value(m.fs.unit.deltaP[0])
        assert pytest.approx(2.278e-3, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "H2O"]
        )
        assert pytest.approx(1.872e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, x_interface_in, "Liq", "NaCl"]
        )
        assert pytest.approx(296.458, rel=1e-3) == value(m.fs.unit.area)
        assert pytest.approx(5.670e-4, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "H2O"]
        )
        assert pytest.approx(2.013e-6, rel=1e-3) == value(
            m.fs.unit.flux_mass_phase_comp[0, 1, "Liq", "NaCl"]
        )
        assert pytest.approx(0.3895, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "H2O"]
        )
        assert pytest.approx(5.792e-4, rel=1e-3) == value(
            m.fs.unit.mixed_permeate[0].flow_mass_phase_comp["Liq", "NaCl"]
        )

    @pytest.mark.unit
    def test_report(self, RO_frame):
        RO_frame.fs.unit.report()
