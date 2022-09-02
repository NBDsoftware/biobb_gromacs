from biobb_common.tools import test_fixtures as fx
from biobb_gromacs.gromacs.genion import genion


class TestGenionDocker:
    def setUp(self):
        fx.test_setup(self, 'genion_docker')

    def tearDown(self):
        #pass
        fx.test_teardown(self)

    def test_genion_docker(self):
        returncode = genion(properties=self.properties, **self.paths)
        assert fx.not_empty(self.paths['output_gro_path'])
        assert fx.equal(self.paths['output_gro_path'], self.paths['ref_output_gro_path'])
        assert fx.not_empty(self.paths['output_top_zip_path'])
        assert fx.equal(self.paths['output_top_zip_path'], self.paths['ref_output_top_zip_path'])
        assert fx.exe_success(returncode)
