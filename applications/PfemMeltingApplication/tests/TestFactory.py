import os
import json

# Importing the Kratos Library
import KratosMultiphysics

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.PfemMeltingApplication.pfem_melting_analysis import PfemMeltingAnalysis

# This utility will control the execution scope
class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

# General test factory
class TestFactory(KratosUnittest.TestCase):
    debug_mode = False

    def setUp(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):

            # Reading the ProjectParameters
            with open(self.file_parameters,'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            if not self.debug_mode:
                parameters['output_processes'] = KratosMultiphysics.Parameters("""{}""")
            
            # Create Model
            model = KratosMultiphysics.Model()

            # To avoid too many prints
            if (parameters["problem_data"]["echo_level"].GetInt() == 0):
                KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)

            self.test = PfemMeltingAnalysis(model,parameters)

    def test_execution(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Run()

    def tearDown(self):
        pass
