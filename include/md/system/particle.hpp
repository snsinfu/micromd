// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SYSTEM_PARTICLE_HPP
#define MD_SYSTEM_PARTICLE_HPP

// This module defines convenience classes used to access particles in a system.

#include <cstddef>
#include <iterator>

#include "../basic_types.hpp"


namespace md
{
    class system;

    // particle_ref provides access to the attributes of a particle in a system.
    struct particle_ref
    {
        md::index index;

        md::scalar& mass;
        md::scalar& mobility;
        md::point& position;
        md::vector& velocity;

        particle_ref(md::system& system, md::index idx);
    };

    // particle_iterator is a forward iterator that scans particles in a system.
    class particle_iterator
    {
    public:
        using value_type = md::particle_ref;
        using reference = md::particle_ref;
        using pointer = void;
        using difference_type = std::ptrdiff_t;
        using iterator_category = std::forward_iterator_tag;

        particle_iterator() = default;

        // Constructor takes a system and the index of the particle pointed to
        // by the constructed iterator.
        particle_iterator(md::system& system, md::index idx)
            : system_{&system}, index_{idx}
        {
        }

        // oeprator== compares the index for equality.
        bool operator==(particle_iterator const& other) const
        {
            return index_ == other.index_;
        }

        // operator!= compares the index for inequality.
        bool operator!=(particle_iterator const& other) const
        {
            return !(*this == other);
        }

        // operator* returns a particle_ref for the particle pointed-to by this
        // iterator.
        md::particle_ref operator*() const
        {
            return md::particle_ref(*system_, index_);
        }

        // operator++ increments the index.
        particle_iterator operator++(int)
        {
            auto copy = *this;
            ++*this;
            return copy;
        }

        // operator++ increments the index.
        particle_iterator& operator++()
        {
            index_++;
            return *this;
        }

    private:
        md::system* system_ = nullptr;
        md::index index_ = 0;
    };
}

#endif
