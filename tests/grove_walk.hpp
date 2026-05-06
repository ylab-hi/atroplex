/*
 * SPDX-License-Identifier: GPLv3
 *
 * Copyright (c) 2025 Richard A. Schäfer
 *
 * This file is part of atroplex and is licensed under the terms of the GPLv3
 * license. See the LICENSE file in the root of the repository for more
 * information.
 *
 * Header-only B+ tree leaf-walk helpers shared across all test binaries.
 *
 * Without this, every test that needs to inspect grove contents copy-pastes
 * a get_root_nodes / get_is_leaf / get_next loop. The audit (#61) called
 * this out as duplicated across 4+ test files; this header is the
 * extraction point.
 */

#ifndef ATROPLEX_TESTS_GROVE_WALK_HPP
#define ATROPLEX_TESTS_GROVE_WALK_HPP

#include <vector>

#include "genomic_feature.hpp"

namespace atroplex::test {

/**
 * Visit every segment key in the grove's spatial index (B+ tree leaves)
 * in per-chromosome leaf-link order, invoking
 * `fn(segment_feature&, key_ptr)` for each. Both live and absorbed
 * (tombstoned) segments are visited; callers that want only live
 * segments should check `seg.absorbed` inside the callback.
 *
 * `Fn` must be callable as `void(segment_feature&, key_ptr)` or
 * `void(const segment_feature&, key_ptr)`. The reference is `const` if
 * the grove is `const`, but every existing caller passes a mutable
 * grove so the non-const overload is sufficient.
 */
template <typename Fn>
void for_each_segment(grove_type& grove, Fn&& fn) {
    for (auto& [seqid, root] : grove.get_root_nodes()) {
        if (!root) continue;
        auto* node = root;
        while (!node->get_is_leaf()) {
            auto& children = node->get_children();
            if (children.empty()) break;
            node = children[0];
        }
        while (node) {
            for (auto* key : node->get_keys()) {
                auto& feature = key->get_data();
                if (!is_segment(feature)) continue;
                fn(get_segment(feature), key);
            }
            node = node->get_next();
        }
    }
}

/**
 * Convenience: collect every live (non-absorbed) segment as
 * `{segment_feature*, key_ptr}` pairs. Built on `for_each_segment`.
 */
struct live_segment_entry {
    const segment_feature* seg;
    key_ptr key;
};

inline std::vector<live_segment_entry> collect_live_segments(grove_type& grove) {
    std::vector<live_segment_entry> out;
    for_each_segment(grove, [&](segment_feature& seg, key_ptr key) {
        if (!seg.absorbed) out.push_back({&seg, key});
    });
    return out;
}

} // namespace atroplex::test

#endif // ATROPLEX_TESTS_GROVE_WALK_HPP
