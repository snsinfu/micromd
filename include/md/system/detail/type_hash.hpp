// Copyright snsinfu 2019.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SYSTEM_DETAIL_TYPE_HASH_HPP
#define MD_SYSTEM_DETAIL_TYPE_HASH_HPP

// This module provides a lightweight perfect hashing system for C++ static
// types.

#include <cstddef>


namespace md
{
    namespace detail
    {
        // type_hash hashes static type to a unique integer within [0,n) where
        // n is the number of hashed types.
        struct type_hash
        {
            // hash<T>::value is a unique integer associated to the type T. Do
            // not use this value before entering main.
            template<typename T>
            struct hash
            {
                static std::size_t const value;
            };

            // size returns the number of hashed types. Do not use this function
            // before entering main.
            static std::size_t size() noexcept
            {
                return get_counter();
            }

        private:
            static std::size_t& get_counter() noexcept
            {
                static std::size_t counter;
                return counter;
            }
        };

        template<typename T>
        std::size_t const type_hash::hash<T>::value = type_hash::get_counter()++;
    }
}

#endif
