/*

// Project Includes
#include "connector.h" // Connector
#include "geometries/line_3d_2.h" // Line3D2

// STL Includes
#include <limits> // std::numeric_limits::max


namespace Kratos {


Connector::Connector()
    : Connector(std::numeric_limits<IndexType>::max(),
                GeometryType::PointType::Pointer(new GeometryType::PointType(1, -1.0, 0.0, 0.0)),
                GeometryType::PointType::Pointer(new GeometryType::PointType(2,  1.0, 0.0, 0.0)),
                Properties::Pointer(new Properties))
{
}


Connector::Connector(IndexType Id,
                     const Node::Pointer& rpNode0,
                     const Node::Pointer& rpNode1,
                     const Properties::Pointer& rpProperties)
    : Connector(Id,
                GeometryType::Pointer(new Line3D2<GeometryType::PointType>(rpNode0, rpNode1)),
                rpProperties)
{
}


Connector::Connector(IndexType Id,
                     const GeometryType::Pointer& rpGeometry,
                     const Properties::Pointer& rpProperties)
    : Element(Id, rpGeometry, rpProperties)
{
}


Element::Pointer Connector::Create(IndexType Id,
                                   GeometryType::Pointer pGeometry,
                                   Properties::Pointer pProperties) const
{
    KRATOS_ERROR_IF_NOT(pGeometry->size() == 2);
    return Element::Pointer(new Connector(Id, pGeometry, pProperties));
}


Element::Pointer Connector::Create(IndexType Id,
                                   const NodesArrayType& rNodes,
                                   Properties::Pointer pProperties) const
{
    KRATOS_ERROR << "factory not implemented\n";
}


void Connector::EquationIdVector(EquationIdVectorType& rOutput,
                                 const ProcessInfo& rCurrentProcessInfo) const
{
    rOutput.resize(6);

    unsigned i_dof = 0u;
    for (unsigned i_node=0u; i_node<2; ++i_node) {
        rOutput[i_dof++] = this->GetGeometry()[i_node].GetDof(DISPLACEMENT_X).EquationId();
        rOutput[i_dof++] = this->GetGeometry()[i_node].GetDof(DISPLACEMENT_Y).EquationId();
        rOutput[i_dof++] = this->GetGeometry()[i_node].GetDof(DISPLACEMENT_Z).EquationId();
        rOutput[i_dof++] = this->GetGeometry()[i_node].GetDof(ROTATION_X).EquationId();
        rOutput[i_dof++] = this->GetGeometry()[i_node].GetDof(ROTATION_Y).EquationId();
        rOutput[i_dof++] = this->GetGeometry()[i_node].GetDof(ROTATION_Z).EquationId();
    }
}


} // namespace Kratos
*/
