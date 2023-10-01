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
        auto &r_process_info = rModelPart.GetProcessInfo();

        auto cl_parameters = ConstitutiveLaw::Parameters();

        cl_parameters.SetOptions(rOptions);
        cl_parameters.SetDeformationGradientF(F);
        cl_parameters.SetDeterminantF(detF);
        cl_parameters.SetStrainVector(strain_vector);
        cl_parameters.SetStressVector(stress_vector);
        cl_parameters.SetConstitutiveMatrix(C);
        cl_parameters.SetShapeFunctionsValues(N);
        cl_parameters.SetShapeFunctionsDerivatives(DN_DX);
        cl_parameters.SetProcessInfo(r_process_info);
        cl_parameters.SetMaterialProperties(rModelPart.GetProperties(0));

        Matrix stress_list = rStrainList;
        stress_list.clear();
        int elem_iterator = 0;

        block_for_each(rModelPart.Elements(), [&](Element& rElement) {
            std::vector<ConstitutiveLaw::Pointer> cl_list;
            rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, cl_list, r_process_info);
            strain_vector[0] = rStrainList(elem_iterator, 0);
            strain_vector[1] = rStrainList(elem_iterator, 1);
            strain_vector[2] = rStrainList(elem_iterator, 2);
            cl_parameters.SetElementGeometry(rElement.GetGeometry());

            cl_list[0]->CalculateMaterialResponseCauchy(cl_parameters); // Here I assume one IP
            stress_list(elem_iterator, 0) = stress_vector[0];
            stress_list(elem_iterator, 1) = stress_vector[1];
            stress_list(elem_iterator, 2) = stress_vector[2];

            cl_list[0]->FinalizeMaterialResponseCauchy(cl_parameters);

            elem_iterator++;
        });

        return stress_list;
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