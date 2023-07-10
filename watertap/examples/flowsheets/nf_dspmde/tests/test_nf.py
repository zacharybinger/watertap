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
from pyomo.environ import value
from watertap.examples.flowsheets.nf_dspmde.nf import main


@pytest.mark.requires_idaes_solver
@pytest.mark.component
def test_main():
    m = main()
    test_dict = {
        "lcow": [m.fs.costing.LCOW, 0.16811587158493219],
        "pressure": [m.fs.NF.pump.outlet.pressure[0] / 1e5, 6.56],
        "area": [m.fs.NF.nfUnit.area, 285.6900547389303],
        "recovery": [
            m.fs.NF.nfUnit.recovery_vol_phase[0.0, "Liq"] * 100,
            73.47934090302432,
        ],
    }
    for (model_result, testval) in test_dict.values():
        assert pytest.approx(testval, rel=1e-3) == value(model_result)
