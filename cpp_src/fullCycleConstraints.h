#pragma once

#include "words.h"

#include <utility>
#include <vector>

namespace low_index {

using FullCycleConstraint =
    std::pair<Relator, unsigned int>;

using FullCycleConstraints =
    std::vector<FullCycleConstraint>;

} // namespace low_index
