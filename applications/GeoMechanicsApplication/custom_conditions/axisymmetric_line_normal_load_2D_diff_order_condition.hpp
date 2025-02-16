// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Vahid Galavi
//

#pragma once

// Project includes
#include "custom_conditions/line_normal_load_2D_diff_order_condition.hpp"
#include "includes/serializer.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) AxisymmetricLineNormalLoad2DDiffOrderCondition : public LineNormalLoad2DDiffOrderCondition
{
public:
    using IndexType      = std::size_t;
    using PropertiesType = Properties;
    using NodeType       = Node;
    using GeometryType   = Geometry<NodeType>;
    using NodesArrayType = GeometryType::PointsArrayType;
    using VectorType     = Vector;
    using MatrixType     = Matrix;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AxisymmetricLineNormalLoad2DDiffOrderCondition);

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default constructor
    AxisymmetricLineNormalLoad2DDiffOrderCondition();

    // Constructor 1
    AxisymmetricLineNormalLoad2DDiffOrderCondition(IndexType NewId, GeometryType::Pointer pGeometry);

    // Constructor 2
    AxisymmetricLineNormalLoad2DDiffOrderCondition(IndexType               NewId,
                                                   GeometryType::Pointer   pGeometry,
                                                   PropertiesType::Pointer pProperties);

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    Condition::Pointer Create(IndexType               NewId,
                              NodesArrayType const&   ThisNodes,
                              PropertiesType::Pointer pProperties) const override;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:
    // Member Variables

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    double CalculateIntegrationCoefficient(const IndexType                    PointNumber,
                                           const GeometryType::JacobiansType& JContainer,
                                           const GeometryType::IntegrationPointsArrayType& IntegrationPoints) const override;

    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, LineNormalLoad2DDiffOrderCondition)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, LineNormalLoad2DDiffOrderCondition)
    }

}; // class AxisymmetricLineNormalLoad2DDiffOrderCondition.

} // namespace Kratos.