// © 2020 ETH Zurich, Mechanics and Materials Lab
//
// This file is part of ae108.
//
// ae108 is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or any
// later version.
//
// ae108 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ae108. If not, see <https://www.gnu.org/licenses/>.

#include "Embedding_Test.h"
#include "ae108/elements/embedding/IsoparametricEmbedding.h"
#include "ae108/elements/embedding/compute_jacobian.h"
#include "ae108/elements/embedding/embed_point.h"
#include "ae108/elements/shape/Quad4.h"
#include "ae108/elements/shape/Seg2.h"
#include <gmock/gmock.h>

using testing::DoubleEq;
using testing::ElementsAre;
using testing::Eq;
using testing::Test;
using testing::Types;

namespace ae108 {
namespace elements {
namespace embedding {
namespace {

struct Seg2_Configuration {
  using Shape = shape::Seg2;
  using Embedding = IsoparametricEmbedding<Shape>;

  static Embedding create_embedding() {
    return Embedding{{{
        {{1.}},
        {{2.}},
    }}};
  }
};

struct Quad4_Configuration {
  using Shape = shape::Quad4;
  using Embedding = IsoparametricEmbedding<Shape>;

  static Embedding create_embedding() {
    return Embedding{{{
        {{0., 0.}},
        {{1., 0.}},
        {{1., 1.}},
        {{0., 1.}},
    }}};
  }
};

using Configurations = Types<Seg2_Configuration, Quad4_Configuration>;

INSTANTIATE_TYPED_TEST_CASE_P(IsoparametricEmbedding_Test, Embedding_Test,
                              Configurations);

struct IsoparametricEmbedding_Seg2_Test : Test {
  using Embedding = Seg2_Configuration::Embedding;
  const Embedding embedding = Seg2_Configuration::create_embedding();
};

TEST_F(IsoparametricEmbedding_Seg2_Test, reference_dimension_is_correct) {
  static_assert(Embedding::reference_dimension() == 1,
                "1 is the dimension of the shape.");
}

TEST_F(IsoparametricEmbedding_Seg2_Test, physical_dimension_is_correct) {
  static_assert(Embedding::physical_dimension() == 1,
                "1 is the dimension of the shape.");
}

TEST_F(IsoparametricEmbedding_Seg2_Test, embedding_point_0_works) {
  const auto result = embed_point(embedding, {{0.}});

  EXPECT_THAT(result, ElementsAre(DoubleEq(1.5)));
}

TEST_F(IsoparametricEmbedding_Seg2_Test, jacobian_is_1_for_id) {
  const Embedding embedding = Embedding{{{
      {{-1.}},
      {{+1.}},
  }}};
  const auto result = compute_jacobian(embedding, {{0.}});

  EXPECT_THAT(result, ElementsAre(ElementsAre(DoubleEq(1.))));
}

struct IsoparametricEmbedding_Quad4_Test : Test {
  using Embedding = Quad4_Configuration::Embedding;
  const Embedding embedding = Quad4_Configuration::create_embedding();
};

TEST_F(IsoparametricEmbedding_Quad4_Test, reference_dimension_is_correct) {
  static_assert(Embedding::reference_dimension() == 2,
                "2 is the dimension of the shape.");
}

TEST_F(IsoparametricEmbedding_Quad4_Test, physical_dimension_is_correct) {
  static_assert(Embedding::physical_dimension() == 2,
                "2 is the dimension of the shape.");
}

TEST_F(IsoparametricEmbedding_Quad4_Test, embedding_point_0_works) {
  const auto result = embed_point(embedding, {{0., 0.}});

  EXPECT_THAT(result, ElementsAre(DoubleEq(.5), DoubleEq(.5)));
}

TEST_F(IsoparametricEmbedding_Quad4_Test, jacobian_is_1_for_id) {
  const Embedding embedding = Embedding{{{
      {{-1., -1.}},
      {{1., -1.}},
      {{1., 1.}},
      {{-1., 1.}},
  }}};
  const auto result = compute_jacobian(embedding, {{0., 0.}});

  EXPECT_THAT(result, ElementsAre(ElementsAre(DoubleEq(1.), DoubleEq(0.)),
                                  ElementsAre(DoubleEq(0.), DoubleEq(1.))));
}

} // namespace
} // namespace embedding
} // namespace elements
} // namespace ae108
