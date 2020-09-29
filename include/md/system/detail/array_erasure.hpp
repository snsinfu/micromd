// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SYSTEM_DETAIL_ARRAY_ERASURE_HPP
#define MD_SYSTEM_DETAIL_ARRAY_ERASURE_HPP

// This module provides a type erasure for resizable arrays. The class is used
// by attribute_table to store differently typed arrays (columns).

#include <cstddef>
#include <memory>
#include <vector>


namespace md
{
    namespace detail
    {
        // array_erasure is a type erasure for resizable arrays.
        class array_erasure
        {
        private:
            array_erasure() = default;

        public:
            virtual ~array_erasure() = default;

            // clone creates a type-erased copy of this instance.
            virtual std::unique_ptr<array_erasure> clone() = 0;

            // size returns the length of the array.
            virtual std::size_t size() const = 0;

            // resize extends or shrinks the array to given length.
            virtual void resize(std::size_t size) = 0;

            // instance is the only allowed implementation of array_erasure.
            template<typename T>
            class instance;

            // make returns a unique_ptr of newly constructed instance<T>.
            template<typename T>
            static std::unique_ptr<instance<T>> make(std::size_t size, T def)
            {
                return std::unique_ptr<instance<T>>{new instance<T>(size, def)};
            }

            // recover returns a reference to the erased instance<T>. Behavior
            // is undefined if the actual type is not T.
            template<typename T>
            instance<T>& recover()
            {
                return static_cast<instance<T>&>(*this);
            }
        };

        template<typename T>
        class array_erasure::instance : public array_erasure
        {
        public:
            // Constructor takes the initial size of an array and the default
            // value used to create new elements.
            instance(std::size_t size, T def)
                : values_(size, def), default_{def}
            {
            }

            // data returns the pointer to the first element of the array.
            T* data()
            {
                return values_.data();
            }

            // size returns the length of the array.
            std::size_t size() const final override
            {
                return values_.size();
            }

            // clone creates a type-erased copy of this instance.
            std::unique_ptr<array_erasure> clone() final override
            {
                return std::unique_ptr<array_erasure>{new instance(*this)};
            }

            // resize extends or shrinks the array to given length. New elements
            // are initialized with the default value.
            void resize(std::size_t size) final override
            {
                values_.resize(size, default_);
            }

        private:
            std::vector<T> values_;
            T default_;
        };
    }
}

#endif
