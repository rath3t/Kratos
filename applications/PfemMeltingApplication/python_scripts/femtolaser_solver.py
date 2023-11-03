from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import math
# Importing the Kratos Library
import KratosMultiphysics

# Import applications
# import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
# import KratosMultiphysics.ConvectionDiffusionApplication as ConvDiff
# import KratosMultiphysics.MeshingApplication as MeshApp
import KratosMultiphysics.PfemMeltingApplication as PfemM

import time as timer

# Importing the base class
import KratosMultiphysics.PfemMeltingApplication.coupled_fluid_thermal_solverwithoutmeshgeneration
BaseClass = KratosMultiphysics.PfemMeltingApplication.coupled_fluid_thermal_solverwithoutmeshgeneration.PfemCoupledFluidThermalSolver

def CreateSolver(main_model_part, custom_settings):

    return FemtolaserSolver(main_model_part, custom_settings)

class FemtolaserSolver(BaseClass):
    def _SetLaser(self, laser):
        self.laser = laser
        self.fluid_solver._SetLaser(laser)
        # self.thermal_solver._SetLaser(laser)
        self.laser_is_set = True

    def Initialize(self):

        if not self.laser_is_set:
            raise('Laser must be set before Initialize()')

        super(FemtolaserSolver, self).Initialize()

        nodal_area_process = KratosMultiphysics.CalculateNodalAreaProcess(
            self.fluid_solver.main_model_part,
            self.domain_size)
        nodal_area_process.Execute()

        nodal_area = self.fluid_solver.main_model_part.Nodes[1214].GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)
        materials_filename = self.settings["thermal_solver_settings"]["material_import_settings"]["materials_filename"].GetString()

        with open(materials_filename, 'r') as parameter_file:
                materials = KratosMultiphysics.Parameters(parameter_file.read())

        material_settings = materials["properties"][0]["Material"]

        cp = material_settings['Variables']['SPECIFIC_HEAT'].GetDouble()
        rho = material_settings['Variables']['DENSITY'].GetDouble()
        Q = 25e-6
        T0 = self.settings['environment_settings']['ambient_temperature'].GetDouble()
        energy_to_temperature = 1.0 / (nodal_area * cp * rho)
        laser_settings_file_name = self.settings['laser_import_settings']['laser_filename'].GetString()
        with open(laser_settings_file_name, 'r') as parameter_file:
            laser_settings = KratosMultiphysics.Parameters(parameter_file.read())

        for node in self.fluid_solver.main_model_part.Nodes:
            if node.Id==1214:
                initial_temp = T0 + energy_to_temperature * Q
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, initial_temp)
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, initial_temp)



    def SolveSolutionStep(self):

        for node in self.fluid_solver.main_model_part.Nodes:
            T = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            if False and node.Z > 0.0049:
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 400)
                node.Fix(KratosMultiphysics.TEMPERATURE)

        t7=timer.time()

        # fluid_is_converged = self.fluid_solver.SolveSolutionStep()
        #self.Streamline.RungeKutta4ElementbasedSI(self.fluid_solver.main_model_part,100)

        for node in self.fluid_solver.main_model_part.Nodes:
            velocity = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY)
            node.SetSolutionStepValue(KratosMultiphysics.MESH_VELOCITY, velocity)

        self.DecompositionUtility.CalculateDecomposition(self.fluid_solver.main_model_part)
        #self.DecompositionUtility.ElementDeactivation(self.fluid_solver.main_model_part)

        thermal_is_converged = self.thermal_solver.SolveSolutionStep()
        # self.CalculateViscosityaux()
        t8=timer.time()
        self.problemsolution=self.problemsolution + t8 - t7
        step=self.fluid_solver.main_model_part.ProcessInfo[KratosMultiphysics.STEP]
        self.outputfile6.write(str(step)+" "+ str(self.streamlineintegration) +" "+ str(self.meshingprocedure) +" "+ str(self.fillingsubmodelparts) +" "+ str(self.initializeSolutionStep)+" "+ str(self.problemsolution) +"\n")

        return thermal_is_converged