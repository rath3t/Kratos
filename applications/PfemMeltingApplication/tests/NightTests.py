# Definition of the classes for the NIGHTLY TESTS

# Import TestFactory
import TestFactory as TF

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

class constant_laser_on_box(TF.TestFactory):
    file_parameters = "element_tests/constant_laser_on_box/ProjectParameters.json"


def SetTestSuite(suites):
    night_suite = suites['nightly']

    night_suite.addTests(
        KratosUnittest.TestLoader().loadTestsFromTestCases([
            constant_laser_on_box,
            # test2
        ])
    )

    return night_suite
