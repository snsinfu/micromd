// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SYSTEM_DETAIL_ARRAY_ERASURE_HPP
#define MD_SYSTEM_DETAIL_ARRAY_ERASURE_HPP

#include <cassert>
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

            // recover returns a reference to the erased instance<T>. Assertion
            // fails if T is not the erased element type.
            template<typename T>
            instance<T>& recover()
            {
                assert(dynamic_cast<instance<T>*>(this) != nullptr);
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
            std::size_t size() const override
            {
                return values_.size();
            }

            // resize extends or shrinks the array to given length. New elements
            // are initialized with the default value.
            void resize(std::size_t size) override
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
