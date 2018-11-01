// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_DETAIL_DYNARRAY_HPP
#define MD_DETAIL_DYNARRAY_HPP

#include <vector>

#include "../typedef.hpp"


namespace md
{
    namespace detail
    {
        // dynarray_base is a type-erased base class of dynarray<T>.
        class dynarray_base
        {
        public:
            virtual ~dynarray_base() = default;

            // resize resizes the internal storage to given size.
            virtual void resize(md::index size) = 0;
        };

        // dynarray is a resizable array of Ts.
        template<typename T>
        class dynarray : public detail::dynarray_base
        {
        public:
            dynarray(md::index size, T def)
                : default_{def}, values_(size, def)
            {
            }

            // resize resizes the internal storage to given size. Newly created
            // elements are filled with the default value.
            void resize(md::index size) override
            {
                values_.resize(size, default_);
            }

            // view returns a view into the array.
            md::array_view<T> view()
            {
                return values_;
            }

        private:
            T default_;
            std::vector<T> values_;
        };
    }
}

#endif
