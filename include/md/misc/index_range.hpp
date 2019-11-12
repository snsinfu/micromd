// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_MISC_INDEX_RANGE_HPP
#define MD_MISC_INDEX_RANGE_HPP

// This module provides a lightweight iterator and range for indices.

#include <cstddef>
#include <iterator>

#include "../basic_types.hpp"


namespace md
{
    // `index_iterator` is a forward iterator for looping over a range of
    // `md::index` values.
    class index_iterator
    {
    public:
        using value_type = md::index;
        using reference = md::index;
        using pointer = void;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::forward_iterator_tag;

        // Default constructor creates an `index_iterator` that points to zero.
        index_iterator() = default;

        // Constructor creates an `index_iterator` that points to the given
        // index value.
        explicit index_iterator(md::index value)
            : value_{value}
        {
        }

        reference operator*() const
        {
            return value_;
        }

        index_iterator& operator++()
        {
            value_++;
            return *this;
        }

        index_iterator operator++(int)
        {
            index_iterator copy = *this;
            operator++();
            return copy;
        }

    private:
        md::index value_ = 0;
    };

    inline bool operator==(md::index_iterator const& a, md::index_iterator const& b)
    {
        return *a == *b;
    }

    inline bool operator!=(md::index_iterator const& a, md::index_iterator const& b)
    {
        return *a != *b;
    }

    // `index_range` is a contiguous range of `md::index` values.
    class index_range
    {
    public:
        // Default constructor creates an empty range.
        index_range() = default;

        // Constructor `index_range(b, e)` creates a range spanning from `b`
        // (inclusive) to `e` (exclusive).
        index_range(md::index b, md::index e)
            : begin_{b}, end_{e}
        {
        }

        // Constructor `index_range(n)` creates a range spanning from 0 to `n`
        // (exclusive).
        explicit index_range(md::index n)
            : begin_{0}, end_{n}
        {
        }

        // `size` returns the number of index values in the range.
        md::index size() const
        {
            return end_ - begin_;
        }

        md::index_iterator begin() const
        {
            return md::index_iterator{begin_};
        }

        md::index_iterator end() const
        {
            return md::index_iterator{end_};
        }

        // Indexing an `index_range` returns the index passed as the argument
        // offsetted by the start of the range.
        md::index operator[](md::index offset) const
        {
            return begin_ + offset;
        }

    private:
        md::index begin_ = 0;
        md::index end_ = 0;
    };
}


#endif
