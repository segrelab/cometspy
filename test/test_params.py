import unittest
import os

print(os.getcwd())

import pandas as pd
from cometspy import params


class TestParams(unittest.TestCase):
    def setUp(self):
        self.params = params()

    def test_set_get_param(self):
        # Test setting and getting a parameter
        param_name = 'timeStep'
        expected_value = 0.01
        self.params.set_param(param_name, expected_value)
        actual_value = self.params.get_param(param_name)
        self.assertEqual(actual_value, expected_value)

        # Test setting a non-existing parameter
        param_name = 'nonExistentParam'
        with self.assertWarns(UserWarning):
            self.params.set_param(param_name, 42)

    def test_show_params(self):
        # Test if show_params() returns a DataFrame
        result = self.params.show_params()
        self.assertIsInstance(result, pd.DataFrame)

    def test_write_params(self):
        # Test writing and reading the parameters to/from files
        out_glb_file = 'test_global_params.txt'
        out_pkg_file = 'test_package_params.txt'

        try:
            # Write params to files
            self.params.write_params(out_glb_file, out_pkg_file)

            # Verify that the files are created
            self.assertTrue(os.path.isfile(out_glb_file))
            self.assertTrue(os.path.isfile(out_pkg_file))

            # Read params from files and compare
            with open(out_glb_file) as f:
                global_params_lines = f.readlines()
            with open(out_pkg_file) as f:
                package_params_lines = f.readlines()

            # Here, you could check if the lines contain the expected parameters and values.

        finally:
            # Clean up: remove the temporary files
            if os.path.isfile(out_glb_file):
                os.remove(out_glb_file)
            if os.path.isfile(out_pkg_file):
                os.remove(out_pkg_file)


if __name__ == '__main__':
    unittest.main()
