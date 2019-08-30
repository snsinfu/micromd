// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SYSTEM_DETAIL_ATTRIBUTE_TABLE_HPP
#define MD_SYSTEM_DETAIL_ATTRIBUTE_TABLE_HPP

// This module defines attribute_table class: A data structure like a columnar
// database (or a dataframe), keyed by attribute_keys.

#include <memory>
#include <vector>

#include "../../basic_types.hpp"
#include "../attribute.hpp"

#include "array_erasure.hpp"
#include "type_hash.hpp"


namespace md
{
    namespace detail
    {
        // attribute_table is a table of arrays (columns) of the same length.
        // Each column is keyed by a tag type.
        class attribute_table
        {
        public:
            attribute_table()
                : arrays_(detail::type_hash::size())
            {
            }

            // size returns the number of elements in the columns.
            md::index size() const
            {
                return size_;
            }

            // resize resizes all the columns to the given size.
            void resize(md::index size)
            {
                for (auto& array : arrays_) {
                    if (array) {
                        array->resize(size);
                    }
                }
                size_ = size;
            }

            // require creates a column for given key if it does not exist. It
            // does nothing otherwise.
            template<typename T, typename Tag>
            void require(md::attribute_key<T, Tag> key)
            {
                auto& entry = arrays_[detail::type_hash::hash<Tag>::value];
                if (!entry) {
                    entry = detail::array_erasure::make<T>(size_, md::default_value(key));
                }
            }

            // view returns a mutable view into the column with given key.
            template<typename T, typename Tag>
            md::array_view<T> view(md::attribute_key<T, Tag>)
            {
                return arrays_[detail::type_hash::hash<Tag>::value]->template recover<T>();
            }

            template<typename T, typename Tag>
            md::array_view<T const> view(md::attribute_key<T, Tag>) const
            {
                return arrays_[detail::type_hash::hash<Tag>::value]->template recover<T>();
            }

        private:
            md::index size_ = 0;
            std::vector<std::unique_ptr<detail::array_erasure>> arrays_;
        };
    }
}

#endif
