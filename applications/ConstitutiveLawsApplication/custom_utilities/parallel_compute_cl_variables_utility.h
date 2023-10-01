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