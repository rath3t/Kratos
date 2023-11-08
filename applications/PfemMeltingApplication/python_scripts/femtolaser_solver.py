from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
import numpy as np
import h5py

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

    def FindCenterNodeId(self):
        for node in self.fluid_solver.main_model_part.Nodes:
            if node.X**2 + node.Y**2 + node.Z**2 <= 1e-12:
                center_Id = node.Id
                break
        return center_Id

    def MonitorEnergy(self):
        energy = 0.0
        T0 = self.settings['environment_settings']['ambient_temperature'].GetDouble()

        for node in self.fluid_solver.main_model_part.Nodes:
            nodal_measure = node.GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)
            energy += (node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - T0) * nodal_measure

        return energy

    def CreateResultsFile(self, filename):
        if os.path.exists(self.results_filename):
            os.remove(self.results_filename)
        with h5py.File(filename, 'a') as f:
            f.attrs['ambient_temperature'] = self.T0
            f.attrs['pulse_energy'] = self.Q
            f.attrs['specific_heat_capacity'] = self.cp
            f.attrs['density'] = self.rho
            f.attrs['conductivity'] = self.conductivity

            # Create a dataset to store the radii
            dataset = f.create_dataset('radii', (self.radii.shape), dtype=self.radii.dtype)
            dataset[:] = self.radii[:]
            f.create_group('temperature_increments')

    def WriteResults(self, filename, process_info):

        step = process_info[KratosMultiphysics.STEP]
        time = step = process_info[KratosMultiphysics.TIME]

        # Open the HDF5 file.
        with h5py.File(filename, 'a') as f:
            assert self.radii.shape  == self.temperature_increments.shape

            # Create a dataset to store the radii and temperatures data.
            dataset = f['/temperature_increments'].create_dataset(str(step), self.temperature_increments.shape, dtype=self.temperature_increments.dtype)

            # Write the radii and temperatures data to the dataset.
            dataset[:] = self.temperature_increments

            # Add a time label to the dataset.
            dataset.attrs["time"] = time

    def print_vectors_to_hdf5_at_every_time_step(self, filename, vec1, vec2, time_step_number):
        # Prints two vectors into an hdf5 file at every time step, labeling the data for each time step with the time step number.

        # Args:
        #     filename: The name of the hdf5 file to print the vectors to.
        #     vec1: A numpy array containing the first vector.
        #     vec2: A numpy array containing the second vector.
        #     time_step_number: The time step number.

        # Open the hdf5 file.
        f = h5py.File(filename, 'a')

        # Create a dataset for each vector at the current time step.
        vec1_dset = f.create_dataset("vec1_" + str(time_step_number), (len(vec1),), dtype=vec1.dtype)
        vec2_dset = f.create_dataset("vec2_" + str(time_step_number), (len(vec2),), dtype=vec2.dtype)

        # Write the vectors to the hdf5 file.
        vec1_dset[:] = vec1
        vec2_dset[:] = vec2

        # Close the hdf5 file.
        f.close()

    def ImposeTemperatureDueToLaser(self):
        def bell_curve(radius_squared, R_far, Tmax):
            # Calculate the z-score of the radius.
            z = radius_squared / R_far**2

            # Calculate the value of the bell curve at the z-score.
            bell_curve_value = Tmax * np.exp(-0.5 * z)

            return bell_curve_value

        T0 = self.settings['environment_settings']['ambient_temperature'].GetDouble()

        # for node in self.fluid_solver.main_model_part.Nodes:
        #     r_2 = node.X**2 + node.Y**2 + node.Z**2
        #     if r_2 < R_far**2:
        #         temp = T0 + bell_curve(r_2, R_far, 0.5 * T0)
        #         node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, temp)
        #         node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, temp)
        center_Id = self.FindCenterNodeId()
        center_node_nodal_area = self.fluid_solver.main_model_part.Nodes[center_Id].GetSolutionStepValue(KratosMultiphysics.NODAL_AREA)

        energy_to_temperature_change = 1.0 / (center_node_nodal_area * self.cp * self.rho)
        for node in self.fluid_solver.main_model_part.Nodes:
            if node.Id == center_Id:
                initial_temp = self.T0 + energy_to_temperature_change * self.Q
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, initial_temp)
                node.SetSolutionStepValue(KratosMultiphysics.TEMPERATURE, 1, initial_temp)

    def Initialize(self):

        if not self.laser_is_set:
            raise('Laser must be set before Initialize()')

        super(FemtolaserSolver, self).Initialize()

        nodal_area_process = KratosMultiphysics.CalculateNodalAreaProcess(
            self.fluid_solver.main_model_part,
            self.domain_size)
        nodal_area_process.Execute()

        materials_filename = self.settings["thermal_solver_settings"]["material_import_settings"]["materials_filename"].GetString()

        with open(materials_filename, 'r') as parameter_file:
                materials = KratosMultiphysics.Parameters(parameter_file.read())

        material_settings = materials["properties"][0]["Material"]

        self.Q = 25e-6
        self.R_far = 0.04
        self.cp = material_settings['Variables']['SPECIFIC_HEAT'].GetDouble()
        self.conductivity = material_settings['Variables']['CONDUCTIVITY'].GetDouble()
        self.rho = material_settings['Variables']['DENSITY'].GetDouble()
        self.T0 = self.settings['environment_settings']['ambient_temperature'].GetDouble()

        laser_settings_file_name = self.settings['laser_import_settings']['laser_filename'].GetString()
        with open(laser_settings_file_name, 'r') as parameter_file:
            laser_settings = KratosMultiphysics.Parameters(parameter_file.read())

        self.ImposeTemperatureDueToLaser()

        self.initial_energy = self.MonitorEnergy()

        radius_2 = lambda node: node.X**2 + node.Y**2 + node.Z**2

        self.near_field_nodes = [node for node in self.fluid_solver.main_model_part.Nodes if radius_2(node) < self.R_far**2]
        self.radii = np.sqrt(np.array([radius_2(node) for node in self.near_field_nodes]))
        self.results_filename = 'results.h5'
        self.CreateResultsFile(self.results_filename)

    def SolveSolutionStep(self):
        # Checking whether the energy is conserved
        energy = self.MonitorEnergy()
        print('*'*100, '\nrelative_energy_change = ', (energy - self.initial_energy) / self.initial_energy)

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

        self.temperature_increments = np.array([node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) - self.T0 for node in self.near_field_nodes])
        self.WriteResults(self.results_filename, self.fluid_solver.main_model_part.ProcessInfo)
        return thermal_is_converged