#include "Element_Test.h"
#include "ae108/elements/CoreElement.h"
#include "ae108/elements/TwoNodeCorotationalBeamElement.h"
#include <gmock/gmock.h>

// #include "Element_Test.h"
// #include "ae108/elements/CoreElement.h"
#include "ae108/elements/compute_consistent_mass_matrix.h"
#include "ae108/elements/compute_lumped_mass_matrix.h"
#include "ae108/elements/embedding/IsoparametricEmbedding.h"
#include "ae108/elements/integrator/IsoparametricIntegrator.h"
#include "ae108/elements/materialmodels/AutomaticStressTrait.h"
#include "ae108/elements/materialmodels/AutomaticTangentMatrixTrait.h"
#include "ae108/elements/materialmodels/Hookean.h"
#include "ae108/elements/quadrature/Quadrature.h"
#include "ae108/elements/shape/Hexa8.h"
#include "ae108/elements/shape/Quad4.h"
#include "ae108/elements/shape/Seg2.h"
// #include <gmock/gmock.h>

using testing::DoubleEq;
using testing::Eq;
using testing::Test;
using testing::Types;


namespace ae108 {
namespace elements {
namespace {
/*
template <class Element_> struct Configuration_1Dxx {
  using Element = Element_;
  static Element create_element() noexcept {
    using Integrator = typename Element::Integrator;
    using MaterialModel = typename Element::MaterialModel;
    using Embedding = typename Integrator::Embedding;

    return Element(MaterialModel(1., 0.), Integrator(Embedding({{
                                              {{0.}},
                                              {{1.}},
                                          }})));
  }

  static typename Element::Time create_time() noexcept {
    return typename Element::Time{0.};
  }
};

using Element_Seg2 = CoreElement<
    materialmodels::Hookean<1>,
    integrator::IsoparametricIntegrator<
        shape::Seg2,
        quadrature::Quadrature<quadrature::QuadratureType::Cube, 1, 3>>>;

using Configurations_1D = Types<Configuration_1Dxx<Element_Seg2>>;
INSTANTIATE_TYPED_TEST_CASE_P(CoreElement_Seg2_Test, Element_Test,
                              Configurations_1D);

struct CoreElement_Seg2_Test : Test {
  using Element = Element_Seg2;
  const Element element = Configuration_1Dxx<Element>::create_element();
};

TEST_F(CoreElement_Seg2_Test, physical_dimension_is_1) {
  EXPECT_THAT(element.dimension(), Eq(1));
}
*/

}
} // namespace elements
} // namespace ae108