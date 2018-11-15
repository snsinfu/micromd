// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SYSTEM_ATTRIBUTE_HPP
#define MD_SYSTEM_ATTRIBUTE_HPP

// This module defines the basics of the attribute mechanism implemented in the
// system class.
//
// Here the C++ type system is abused: An attribute is completely described by
// a function pointer:
//
//     T(* attribute)(Tag*).
//
// The return type T is used as the type of attribute values. And the dummy
// parameter type Tag, which is supposed to be a unique struct, is used to
// identify the attribute.

#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>


namespace md
{
    // attribute_key is a type used to statically index system attribute. T is
    // the type of attribute values and Tag is a unique tag type.
    template<typename T, typename Tag = T>
    using attribute_key = T(*)(Tag*);

    // default_value returns the default value associated with given key.
    template<typename T, typename Tag>
    constexpr T default_value(md::attribute_key<T, Tag> key)
    {
        return key ? key(nullptr) : T{};
    }
}

#endif
