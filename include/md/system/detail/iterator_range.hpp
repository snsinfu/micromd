// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_SYSTEM_DETAIL_ITERATOR_RANGE_HPP
#define MD_SYSTEM_DETAIL_ITERATOR_RANGE_HPP

// This module provides a basic class for adapting iterator pair as a range.

namespace md
{
    namespace detail
    {
        template<typename Iterator>
        class iterator_range
        {
        public:
            iterator_range(Iterator begin, Iterator end)
                : begin_{begin}, end_{end}
            {
            }

            Iterator begin() const
            {
                return begin_;
            }

            Iterator end() const
            {
                return end_;
            }

        private:
            Iterator begin_;
            Iterator end_;
        };
    }
}

#endif
