//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if defined(KRATOS_PYTHON)

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "linear_solvers/linear_solver.h"
#include "processes/process.h"
#include "processes/variational_distance_calculation_process.h"
#include "trilinos_space.h"
#include "spaces/ublas_space.h"
#include "custom_python/add_trilinos_processes_to_python.h"
#include "custom_processes/trilinos_levelset_convection_process.h"
#include "custom_strategies/builder_and_solvers/trilinos_block_builder_and_solver.h"

namespace Kratos::Python
{
namespace py = pybind11;

typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

// Helpers to define Trilinos VariationalDistanceCalculatorProcess
template<unsigned int TDim> using TrilinosVariationalDistanceCalculation = VariationalDistanceCalculationProcess<TDim,TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType>;

template< class TBinder, unsigned int TDim > void DistanceCalculatorConstructionHelper(TBinder& rBinder)
{
    rBinder.def(py::init([](
        Epetra_MpiComm& rCommunicator, Model& rModel,TrilinosLinearSolverType::Pointer pLinearSolver, Parameters ThisParameters)
        {
            constexpr int row_size_guess = TDim == 2 ? 15 : 40;
            auto p_builder_solver = Kratos::make_shared<TrilinosBlockBuilderAndSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > >(
                rCommunicator, row_size_guess, pLinearSolver);
            return Kratos::make_shared<TrilinosVariationalDistanceCalculation<TDim>>(rModel, pLinearSolver, p_builder_solver, ThisParameters);
        }));
    rBinder.def(py::init([](
        Epetra_MpiComm& rComm,ModelPart& rModelPart,TrilinosLinearSolverType::Pointer pLinearSolver,
        unsigned int MaxIter,Flags TheFlags)
        {
            constexpr int RowSizeGuess = (TDim == 2 ? 15 : 40);
            auto p_builder_solver = Kratos::make_shared<TrilinosBlockBuilderAndSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > >(
                rComm, RowSizeGuess, pLinearSolver);
            return Kratos::make_shared<TrilinosVariationalDistanceCalculation<TDim>>(rModelPart, pLinearSolver, p_builder_solver, MaxIter, TheFlags);
        }));
    rBinder.def(py::init([](
        Epetra_MpiComm& rComm,ModelPart& rModelPart,TrilinosLinearSolverType::Pointer pLinearSolver,
        unsigned int MaxIter,Flags TheFlags,std::string& rAuxName)
        {
            constexpr int RowSizeGuess = (TDim == 2 ? 15 : 40);
            auto p_builder_solver = Kratos::make_shared<TrilinosBlockBuilderAndSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > >(
                rComm, RowSizeGuess, pLinearSolver);
            return Kratos::make_shared<TrilinosVariationalDistanceCalculation<TDim>>(rModelPart, pLinearSolver, p_builder_solver, MaxIter, TheFlags, rAuxName);
        }));
    rBinder.def(py::init([](
        Epetra_MpiComm& rComm,ModelPart& rModelPart,TrilinosLinearSolverType::Pointer pLinearSolver,
        unsigned int MaxIter,Flags TheFlags,std::string& rAuxName, double Coefficient1)
        {
            constexpr int RowSizeGuess = (TDim == 2 ? 15 : 40);
            auto p_builder_solver = Kratos::make_shared<TrilinosBlockBuilderAndSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > >(
                rComm, RowSizeGuess, pLinearSolver);
            return Kratos::make_shared<TrilinosVariationalDistanceCalculation<TDim>>(rModelPart, pLinearSolver, p_builder_solver, MaxIter, TheFlags, rAuxName, Coefficient1);
        }));
    rBinder.def(py::init([](
        Epetra_MpiComm& rComm,ModelPart& rModelPart,TrilinosLinearSolverType::Pointer pLinearSolver,
        unsigned int MaxIter,Flags TheFlags,std::string& rAuxName, double Coefficient1, double Coefficient2)
        {
            constexpr int RowSizeGuess = (TDim == 2 ? 15 : 40);
            auto p_builder_solver = Kratos::make_shared<TrilinosBlockBuilderAndSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType > >(
                rComm, RowSizeGuess, pLinearSolver);
            return Kratos::make_shared<TrilinosVariationalDistanceCalculation<TDim>>(rModelPart, pLinearSolver, p_builder_solver, MaxIter, TheFlags, rAuxName, Coefficient1, Coefficient2);
        }));
}

void AddProcesses(pybind11::module& m)
{
    // Variational distance calculation processes
    using DistanceCalculator2DBinderType = py::class_<TrilinosVariationalDistanceCalculation<2>, typename TrilinosVariationalDistanceCalculation<2>::Pointer, Process >;
    using DistanceCalculator3DBinderType = py::class_<TrilinosVariationalDistanceCalculation<3>, typename TrilinosVariationalDistanceCalculation<3>::Pointer, Process >;

    auto distance_calculator_2d_binder = DistanceCalculator2DBinderType(m,"TrilinosVariationalDistanceCalculationProcess2D");
    auto distance_calculator_3d_binder = DistanceCalculator3DBinderType(m,"TrilinosVariationalDistanceCalculationProcess3D");

    DistanceCalculatorConstructionHelper<DistanceCalculator2DBinderType,2>(distance_calculator_2d_binder);
    DistanceCalculatorConstructionHelper<DistanceCalculator3DBinderType,3>(distance_calculator_3d_binder);

    // Level set convection processes
    typedef LevelSetConvectionProcess<2, TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> BaseLevelSetConvectionProcess2D;
    typedef LevelSetConvectionProcess<3, TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> BaseLevelSetConvectionProcess3D;

    py::class_<BaseLevelSetConvectionProcess2D, BaseLevelSetConvectionProcess2D::Pointer, Process>(m,"BaseTrilinosLevelSetConvectionProcess2D");
    py::class_<BaseLevelSetConvectionProcess3D, BaseLevelSetConvectionProcess3D::Pointer, Process>(m,"BaseTrilinosLevelSetConvectionProcess3D");

    typedef TrilinosLevelSetConvectionProcess<2, TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosLevelSetConvectionProcess2D;
    py::class_<TrilinosLevelSetConvectionProcess2D, TrilinosLevelSetConvectionProcess2D::Pointer, BaseLevelSetConvectionProcess2D>(m, "TrilinosLevelSetConvectionProcess2D")
        .def(py::init<Epetra_MpiComm&, Model&, TrilinosLinearSolverType::Pointer, Parameters>())
        .def(py::init<Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, Parameters>())
        ;

    typedef TrilinosLevelSetConvectionProcess<3, TrilinosSparseSpaceType, TrilinosLocalSpaceType, TrilinosLinearSolverType> TrilinosLevelSetConvectionProcess3D;
    py::class_<TrilinosLevelSetConvectionProcess3D, TrilinosLevelSetConvectionProcess3D::Pointer, BaseLevelSetConvectionProcess3D>(m, "TrilinosLevelSetConvectionProcess3D")
        .def(py::init<Epetra_MpiComm&, Model&, TrilinosLinearSolverType::Pointer, Parameters>())
        .def(py::init<Epetra_MpiComm&, ModelPart&, TrilinosLinearSolverType::Pointer, Parameters>())
        ;

}

}

#endif
