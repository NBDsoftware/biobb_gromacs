from biobb_common.tools import test_fixtures as fx
from biobb_gromacs.gromacs.make_ndx import make_ndx


class TestMakeNdx:
    def setup_class(self):
        fx.test_setup(self, 'make_ndx')

    def teardown_class(self):
        fx.test_teardown(self)

    def test_make_ndx(self):
        returncode = make_ndx(properties=self.properties, **self.paths)
        assert fx.not_empty(self.paths['output_ndx_path'])
        assert fx.equal(self.paths['output_ndx_path'], self.paths['ref_output_ndx_path'])
        assert fx.exe_success(returncode)
