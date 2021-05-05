import pytest

@pytest.mark.serial
@pytest.mark.forces
def test_VO2_forces(execute_fleur, grep_number, grep_exists):
    """Check the computation of normal Fleur forces for the starting density.

    Simple test of Fleur with one step:
    1.Generate a starting density.
    2.Calculate the resulting forces.
    """
    test_file_folder = './inputfiles/VO2_forces/'

    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    fx_a8 = grep_number(res_files['out'], "FX_A8", res_index=None)
    assert abs(fx_a8[0] - -0.000153) <= 0.0001
    assert abs(fx_a8[-1] - -0.106870) <= 0.0001
    fy_a8 = grep_number(res_files['out'], "FY_A8", res_index=None)
    assert abs(fy_a8[0] - -0.124257) <= 0.0001
    assert abs(fy_a8[-1] - 0.052459) <= 0.0001
    fz_a8 = grep_number(res_files['out'], "FZ_A8", res_index=None)
    assert abs(fz_a8[0] - -0.060538) <= 0.0001
    assert abs(fz_a8[-1] - 0.120169) <= 0.0001

    fx_a3 = grep_number(res_files['out'], "FX_A3", res_index=None)
    assert abs(fx_a3[0] - -0.043394) <= 0.0001
    assert abs(fx_a3[-1] - -0.181073) <= 0.0001
    fy_a3 = grep_number(res_files['out'], "FY_A3", res_index=None)
    assert abs(fy_a3[0] - -0.161893) <= 0.0001
    assert abs(fy_a3[-1] - 0.085353) <= 0.0001
    fz_a3 = grep_number(res_files['out'], "FZ_A3", res_index=None)
    assert abs(fz_a3[0] - -0.042752) <= 0.0001
    assert abs(fz_a3[-1] - 0.200983) <= 0.0001

    fx_a4 = grep_number(res_files['out'], "FX_A4", res_index=None)
    assert abs(fx_a4[0] - 0.018906) <= 0.0001
    assert abs(fx_a4[-1] - 0.045934) <= 0.0001
    fy_a4 = grep_number(res_files['out'], "FY_A4", res_index=None)
    assert abs(fy_a4[0] - 0.074452) <= 0.0001
    assert abs(fy_a4[-1] - -0.021553) <= 0.0001
    fz_a4 = grep_number(res_files['out'], "FZ_A4", res_index=None)
    assert abs(fz_a4[0] - 0.020155) <= 0.0001
    assert abs(fz_a4[-1] - -0.050975) <= 0.0001

    fx_tot = grep_number(res_files['out'], "FX_TOT", res_index=None)
    assert abs(fx_tot[0] - -0.024332) <= 0.0001
    assert abs(fx_tot[-1] - -0.117370) <= 0.0001
    fy_tot = grep_number(res_files['out'], "FY_TOT", res_index=None)
    assert abs(fy_tot[0] - -0.065113) <= 0.0001
    assert abs(fy_tot[-1] - 0.056519) <= 0.0001
    fz_tot = grep_number(res_files['out'], "FZ_TOT", res_index=None)
    assert abs(fz_tot[0] - -0.012405) <= 0.0001
    assert abs(fz_tot[-1] - 0.130208) <= 0.0001

@pytest.mark.serial
@pytest.mark.forces
def test_VO2_force_levels(execute_fleur, grep_number, grep_exists):
    """Check the computation of Aaron's force levels for the starting density.

    Simple test of Fleur with one steps:
    1.Generate a starting density.
    2.Calculate the resulting forces on all 3 force levels as implemented by Aaron K.
    """
    test_file_folder = './inputfiles/VO2_force_levels/'

    res_files = execute_fleur(test_file_folder)
    res_file_names = list(res_files.keys())
    should_files = ['out']
    for file1 in should_files:
        assert (file1 in res_file_names), f'{file1} missing'

    assert grep_exists(res_files['out'], "FORCES: SURFACE CORRECTION")

    fy_is = grep_number(res_files['out'], "FY_IS")
    assert abs(fy_is - -0.000021) <= 0.00001
    fz_is = grep_number(res_files['out'], "FZ_IS")
    assert abs(fz_is - -0.000003) <= 0.00001

    fx_mt = grep_number(res_files['out'], "FX_MT", res_index=None)
    assert abs(fx_mt[0] - -0.000090) <= 0.00001
    assert abs(fx_mt[-1] - 0.005375) <= 0.00001
    fy_mt = grep_number(res_files['out'], "FY_MT", res_index=None)
    assert abs(fy_mt[0] - -0.000250) <= 0.00001
    assert abs(fy_mt[-1] - -0.020187) <= 0.00001
    fz_mt = grep_number(res_files['out'], "FZ_MT", res_index=None)
    assert abs(fz_mt[0] - -0.000308) <= 0.00001
    assert abs(fz_mt[-1] - 0.001556) <= 0.00001

    fx_sf = grep_number(res_files['out'], "FX_SF", res_index=None)
    assert abs(fx_sf[0] - -0.007899) <= 0.0001
    assert abs(fx_sf[-1] - 0.040664) <= 0.0001
    fy_sf = grep_number(res_files['out'], "FY_SF", res_index=None)
    assert abs(fy_sf[0] - 0.027425) <= 0.0001
    assert abs(fy_sf[-1] - -0.021777) <= 0.0001
    fz_sf = grep_number(res_files['out'], "FZ_SF", res_index=None)
    assert abs(fz_sf[0] - 0.015781) <= 0.0001
    assert abs(fz_sf[-1] - -0.046126) <= 0.0001

    fx_tot = grep_number(res_files['out'], "FX_TOT", res_index=None)
    assert abs(fx_tot[0] - -0.033440) <= 0.0001
    assert abs(fx_tot[-1] - -0.071340) <= 0.0001
    fy_tot = grep_number(res_files['out'], "FY_TOT", res_index=None)
    assert abs(fy_tot[0] - -0.035890) <= 0.0001
    assert abs(fy_tot[-1] - 0.014537) <= 0.0001
    fz_tot = grep_number(res_files['out'], "FZ_TOT", res_index=None)
    assert abs(fz_tot[0] - 0.004741) <= 0.0001
    assert abs(fz_tot[-1] - 0.085641) <= 0.0001
