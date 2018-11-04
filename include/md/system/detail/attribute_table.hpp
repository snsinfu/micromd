// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SYSTEM_DETAIL_ATTRIBUTE_TABLE_HPP
#define MD_SYSTEM_DETAIL_ATTRIBUTE_TABLE_HPP

#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../../basic_types.hpp"
#include "../attribute.hpp"
#include "array_erasure.hpp"


namespace md
{
    namespace detail
    {
        using type_hash_t = void const*;

        // type_hash assigns a unique value for each given type.
        template<typename Tag>
        inline type_hash_t type_hash() noexcept
        {
            static char dummy;
            return &dummy;
        }

        // attribute_table is a table of arrays (columns) of the same length.
        // Each column is keyed by a tag type.
        class attribute_table
        {
        public:
            // size returns the number of elements in the contained arrays.
            md::index size() const
            {
                return size_;
            }

            // resize resizes all the contained arrays to the given size.
            void resize(md::index size)
            {
                for (auto& node : arrays_) {
                    node.second->resize(size);
                }
                size_ = size;
            }

            // require creates a column for given key if it does not exist, or
            // does nothing otherwise.
            template<typename T, typename Tag>
            void require(md::attribute_key<T, Tag> key)
            {
                type_hash_t const tag_key = type_hash<Tag>();

                if (arrays_.find(tag_key) == arrays_.end()) {
                    arrays_.emplace(
                        tag_key,
                        detail::array_erasure::make<T>(size_, md::default_value(key))
                    );
                }
            }

            // view returns a mutable view into the array.
            template<typename T, typename Tag>
            md::array_view<T> view(md::attribute_key<T, Tag>)
            {
                return arrays_.at(type_hash<Tag>())->template recover<T>();
            }

            template<typename T, typename Tag>
            md::array_view<T const> view(md::attribute_key<T, Tag>) const
            {
                return arrays_.at(type_hash<Tag>())->template recover<T>();
            }

        private:
            md::index size_ = 0;
            std::unordered_map<type_hash_t, std::unique_ptr<detail::array_erasure>> arrays_;
        };
    }
}

#endif
