# type: ignore
from biobb_common.tools import test_fixtures as fx
from biobb_gromacs.gromacs.solvate import solvate
import pytest


class TestSolvateSingularity():
    def setup_class(self):
        fx.test_setup(self, 'solvate_singularity')

    def teardown_class(self):
        # pass
        fx.test_teardown(self)

    @pytest.mark.skip(reason="singularity currently not available")
    def test_solvate_singularity(self):
        returncode = solvate(properties=self.properties, **self.paths)
        assert fx.not_empty(self.paths['output_gro_path'])
        assert fx.equal(self.paths['output_gro_path'], self.paths['ref_output_gro_path'])
        assert fx.not_empty(self.paths['output_top_zip_path'])
        assert fx.equal(self.paths['output_top_zip_path'], self.paths['ref_output_top_zip_path'])
        assert fx.exe_success(returncode)
