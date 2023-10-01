// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    A. Cornejo
//
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include "includes/model_part.h"

// Include kratos definitions

// Project includes

// Configures

// External includes

namespace Kratos {
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class ParallelComputeCLVariablesUtility
 * @ingroup StructuralMechanicsApplication
 * @brief Node Search
 * @details This class provides several methods to perform paralelized loops in c++ to compute stresses and store them in a matrix
 * within a given radius.
 * @author Manuel Messmer
 */

class ParallelComputeCLVariablesUtility
{
    public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ParallelComputeCLVariablesUtility
    KRATOS_CLASS_POINTER_DEFINITION(ParallelComputeCLVariablesUtility);


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ParallelComputeCLVariablesUtility() {}

    /// Destructor.
    ~ParallelComputeCLVariablesUtility(){}

    Matrix ComputeStresses(ModelPart& rModelPart, const Matrix& rStrainList, const Flags& rOptions)
    {

        // #here set the data needed (note that N and DN_DX are not really needed so we set them to dymmy values)
        // N = KratosMultiphysics.Vector(3)
        // DN_DX = KratosMultiphysics.Matrix(3,2)
        // F = KratosMultiphysics.Matrix(np.eye(2)) #F is not needed for small deformation CLs
        // detF = 1.0
        // strain_vector = KratosMultiphysics.Vector(np.zeros(3))
        // stress_vector = KratosMultiphysics.Vector(np.zeros(3))
        // constitutive_matrix = KratosMultiphysics.Matrix(np.zeros((3,3)))

        // cl_params = KratosMultiphysics.ConstitutiveLawParameters()
        // cl_params.SetOptions(cl_options)
        // cl_params.SetDeformationGradientF(F)
        // cl_params.SetDeterminantF(detF)
        // cl_params.SetStrainVector(strain_vector)
        // cl_params.SetStressVector(stress_vector)
        // cl_params.SetConstitutiveMatrix(constitutive_matrix)
        // cl_params.SetShapeFunctionsValues(N)
        // cl_params.SetShapeFunctionsDerivatives(DN_DX)
        // cl_params.SetProcessInfo(mp.ProcessInfo)
        // cl_params.SetMaterialProperties(prop)
        // i = 0
        // for elem in mp.Elements:            
        //     cl = self.cl_list[i]
        //     strain_vector[:] = real_strains[i,:]
        //     cl_params.SetElementGeometry(elem.GetGeometry())

        //     cl.CalculateMaterialResponseCauchy(cl_params)
        //     stress_matrix[i,:] = stress_vector
        //     cl.FinalizeMaterialResponseCauchy(cl_params)            
        //     i = i+1
        // return stress_matrix

        Vector N(3);
        N.clear();
        Matrix DN_DX(3, 2);
        DN_DX.clear();
        Matrix F(2, 2);
        F.clear();
        F(0, 0) = 1.0;
        F(1, 1) = 1.0;
        double detF = 1.0;
        Vector strain_vector(3);
        Vector stress_vector(3);
        strain_vector.clear();
        stress_vector.clear();
        Matrix C(3, 3);
        C.clear();

        auto cl_parameters = ConstitutiveLaw::Parameters();

        cl_parameters.SetOptions(rOptions);
        cl_parameters.SetDeformationGradientF(F);
        cl_parameters.SetDeterminantF(detF);
        cl_parameters.SetStrainVector(strain_vector);
        cl_parameters.SetStressVector(stress_vector);
        cl_parameters.SetConstitutiveMatrix(C);
        cl_parameters.SetShapeFunctionsValues(N);
        cl_parameters.SetShapeFunctionsDerivatives(DN_DX);
        cl_parameters.SetProcessInfo(rModelPart.GetProcessInfo());
        cl_parameters.SetMaterialProperties(rModelPart.GetProperties(0));

        return Matrix();
    }

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}

    private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

    }; // Class ParallelComputeCLVariablesUtility

///@}

///@} addtogroup block

}  // namespace Kratos.