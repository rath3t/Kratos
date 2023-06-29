# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Import base class file
from KratosMultiphysics.PfemMeltingApplication.fluid_solver_customized import FluidSolverCustomized

def CreateSolver(model, custom_settings):
    return FemtolaserNavierStokesSolverMonolithic(model, custom_settings)

class FemtolaserNavierStokesSolverMonolithic(FluidSolverCustomized):

    def __init__(self, model, custom_settings):
        super(FemtolaserNavierStokesSolverMonolithic,self).__init__(model,custom_settings)

        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction of FemtolaserNavierStokesSolverMonolithic finished.")

    def _SetLaser(self, laser):
        self.laser = laser
        self.laser_is_set = True

    def Initialize(self):
        if not self.laser_is_set:
            raise('Laser must be set before Initialize()')

        super(FemtolaserNavierStokesSolverMonolithic, self).Initialize()

    def AdvanceInTime(self, current_time):

        dt = self.settings["time_stepping"]["time_step"].GetDouble()
        dt_reduction_coefficient = 0.1

        if self.laser.IsOn(current_time):
            dt *= dt_reduction_coefficient

        new_time = current_time + dt

        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time
