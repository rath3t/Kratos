#pragma once

// Project Includes
#include "includes/element.h"


namespace Kratos {


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) Connector : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Connector);

    Connector();

    Connector(IndexType Id,
              const Node::Pointer& rpNode0,
              const Node::Pointer& rpNode1,
              const Properties::Pointer& rpProperties);

    Element::Pointer Create(IndexType Id,
                            GeometryType::Pointer pGeometry,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType Id,
                            const NodesArrayType& rNodes,
                            PropertiesType::Pointer pProperties) const override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rProcessInfo) const override;

    void GetDofList(DofsVectorType& rElementalDofList,
                    const ProcessInfo& rProcessInfo) const override;

    void Initialize(const ProcessInfo& rProcessInfo) override;

    void Calculate(const Variable<Matrix>& rVariable, Matrix& rOutput, const ProcessInfo& rProcessInfo) override;

    void Calculate(const Variable<double>& rVariable, double& rOutput, const ProcessInfo& rProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>& rOutput,
                                      const ProcessInfo& rProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable,
                                      std::vector< array_1d<double,3>>& rOutput,
                                      const ProcessInfo& rProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>& rOutput,
                                      const ProcessInfo& rProcessInfo) override;

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                              VectorType& rRightHandSideVector,
                              const ProcessInfo& rProcessInfo) override;


    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                const ProcessInfo& rProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               const ProcessInfo& rProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix,
                             const ProcessInfo& rProcessInfo) override;

    void GetValuesVector(Vector& rValues,
                         int Step = 0) const override;

    void GetSecondDerivativesVector(Vector& rValues,
                                    int Step = 0) const override;

    void GetFirstDerivativesVector(Vector& rValues,
                                   int Step = 0) const override;

    int Check(const ProcessInfo& rProcessInfo) const override;

    void FinalizeSolutionStep(const ProcessInfo& rProcessInfo) override;

    const Parameters GetSpecifications() const override;

private:
    Connector(IndexType Id,
              const GeometryType::Pointer& rpGeometry,
              const Properties::Pointer& rpProperties);

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;
};


} // namespace Kratos
