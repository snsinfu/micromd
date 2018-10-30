#ifndef MD_HPP
#define MD_HPP

// class array_view
// class vector
// class point
//
// class neighbor_searcher
// class neighbor_list

//------------------------------------------------------------------------------
// vendor - array_view
//------------------------------------------------------------------------------

#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <type_traits>
#include <utility>


namespace md
{
    namespace array_view_detail
    {
        // data returns a pointer to the first element of a contiguous container
        // or an array.
        template<typename Cont>
        constexpr auto data(Cont& cont) noexcept -> decltype(cont.data())
        {
            return cont.data();
        }

        template<typename T, std::size_t N>
        constexpr auto data(T(& arr)[N]) noexcept -> T*
        {
            return arr;
        }

        // size returns the number of elements in a container or an array.
        template<typename Cont>
        constexpr auto size(Cont const& cont) noexcept -> decltype(cont.size())
        {
            return cont.size();
        }

        template<typename T, std::size_t N>
        constexpr auto size(T(&)[N]) noexcept -> std::size_t
        {
            return N;
        }
    }

    // array_view is an non-owning range view of an array. It is similar to
    // GSL's span but uses size_t for indexing, for better compatibility with
    // the standard library.
    template<typename T>
    class array_view
    {
    public:
        // value_type is the non-qualified type of the elements.
        using value_type = typename std::remove_cv<T>::type;

        // pointer is the type of a pointer to an element.
        using pointer = T*;

        // reference is the type of a reference to an element.
        using reference = T&;

        // size_type is the type of size and index values.
        using size_type = std::size_t;

        // iterator is the type of iterators. Guaranteed to be a random access
        // iterator. Currently this is implemented as an alias of T* but this
        // may change in the future.
        using iterator = T*;

        // reverse_iterator is the type of reverse iterators. Guaranteed to be
        // a random access iterator.
        using reverse_iterator = std::reverse_iterator<iterator>;

        // const_array_view is the type of read-only array_view with the same
        // value_type. This is the same as array_view<T> if T is already const
        // qualified.
        using const_array_view = array_view<T const>;

        // array_view is default constructible.
        array_view() = default;

        // array_view is copy constructible.
        array_view(array_view&) = default;

        // array_view is copy constructible.
        array_view(array_view const&) = default;

        // array_view is constructible from a raw memory span [p, p + n). The
        // behavior is undefined if the memory region is invalid.
        constexpr array_view(pointer p, size_type n)
            : data_{p}, size_{n}
        {
        }

        // array_view is implicitly convertible from a contiguous container like
        // std::vector.
        template<
            typename Container,
            typename P = decltype(array_view_detail::data(std::declval<Container&>())),
            typename S = decltype(array_view_detail::size(std::declval<Container&>())),
            typename = typename std::enable_if<
                std::is_convertible<typename std::remove_pointer<P>::type(*)[], T(*)[]>::value,
                int
            >::type,
            typename = typename std::enable_if<
                !std::is_same<Container, array_view<value_type>>::value,
                int
            >::type
        >
        constexpr array_view(Container& container) noexcept
            : array_view{array_view_detail::data(container), array_view_detail::size(container)}
        {
        }

        // array_view is implicitly convertible to a read-only view.
        constexpr operator const_array_view() const noexcept
        {
            return {data(), size()};
        }

        // empty returns true if the view is empty.
        constexpr bool empty() const noexcept
        {
            return size() == 0;
        }

        // size returns the number of elements.
        constexpr size_type size() const noexcept
        {
            return size_;
        }

        // data r eturns a pointer to the first element. The pointer may or may
        // not be null if the view is empty.
        constexpr pointer data() const noexcept
        {
            return data_;
        }

        // front returns a reference to the first element. The behavior is
        // undefined if the view is empty.
        constexpr reference front() const
        {
            return operator[](0);
        }

        // back returns a reference to the last element. The behavior is
        // undefined if the view is empty.
        constexpr reference back() const
        {
            return operator[](size() - 1);
        }

        // Indexing operator returns a reference to the idx-th element. The
        // behavior is undefined if the index is out of bounds.
        constexpr reference operator[](size_type idx) const
        {
            return data()[idx];
        }

        // at returns a reference to the idx-th element. It throws
        // std::out_of_range if the index is out of bounds.
        reference at(size_type idx) const
        {
            if (idx >= size()) {
                throw std::out_of_range("array_view access out-of-bounds");
            }
            return operator[](idx);
        }

        // begin returns an iterator to the beginning of the viewed array.
        constexpr iterator begin() const noexcept
        {
            return data();
        }

        // end returns an iterator to the end of the viewed array.
        constexpr iterator end() const noexcept
        {
            return data() + size();
        }

        // rbegin returns a reverse iterator to the reverse beginning.
        reverse_iterator rbegin() const noexcept
        {
            return reverse_iterator{end()};
        }

        // rend returns a reverse iterator to the reverse end.
        reverse_iterator rend() const noexcept
        {
            return reverse_iterator{begin()};
        }

        // as_const returns a read-only view of the same array.
        constexpr const_array_view as_const() const noexcept
        {
            return {data(), size()};
        }

        // swap swaps the viewed array with other.
        void swap(array_view& other) noexcept
        {
            auto const copy = *this;
            *this = other;
            other = copy;
        }

        // subview returns a view of the subarray with given region.
        constexpr array_view subview(size_type offset, size_type count) const
        {
            return {data() + offset, count};
        }

        // first returns a view of the first count elements.
        constexpr array_view first(size_type count) const
        {
            return subview(0, count);
        }

        /// last returns a view of the last count elements.
        constexpr array_view last(size_type count) const
        {
            return subview(size() - count, count);
        }

    private:
        pointer data_ = nullptr;
        size_type size_ = 0;
    };
}


//------------------------------------------------------------------------------
// vendor - point
//------------------------------------------------------------------------------

#include <cmath>
#include <cstddef>
#include <istream>
#include <ostream>


namespace md
{
    namespace point_detail
    {
        // coordinates hold three-dimensional cartesian coordinate values.
        struct coordinates
        {
            double x = 0;
            double y = 0;
            double z = 0;

            // Default constructor initializes coordinates to zero.
            coordinates() = default;

            // Constructor initializes coordinates to specified values.
            inline coordinates(double x_init, double y_init, double z_init)
                : x{x_init}, y{y_init}, z{z_init}
            {
            }

            // Index 0, 1 and 2 correspond to x, y and z, respectively.
            inline double& operator[](std::size_t idx)
            {
                double* coords[] = {
                    &x,
                    &y,
                    &z
                };
                return *coords[idx];
            }

            // Index 0, 1 and 2 correspond to x, y and z, respectively.
            inline double const& operator[](std::size_t idx) const
            {
                double const* coords[] = {
                    &x,
                    &y,
                    &z
                };
                return *coords[idx];
            }
        };

        // ostream operator formats coords as space-delimited floating-point
        // coordinate values.
        template<typename Char, typename Tr>
        std::basic_ostream<Char, Tr>& operator<<(
            std::basic_ostream<Char, Tr>& os,
            coordinates const& coords
        )
        {
            using sentry_type = typename std::basic_ostream<Char, Tr>::sentry;

            if (sentry_type sentry{os}) {
                Char const delim = os.widen(' ');
                os << coords.x << delim << coords.y << delim << coords.z;
            }

            return os;
        }

        // istream operator parses space-delimited floating-point coordinate
        // values into coords.
        template<typename Char, typename Tr>
        std::basic_istream<Char, Tr>& operator>>(
            std::basic_istream<Char, Tr>& is,
            coordinates& coords
        )
        {
            using sentry_type = typename std::basic_istream<Char, Tr>::sentry;

            if (sentry_type sentry{is}) {
                is >> coords.x >> coords.y >> coords.z;
            }

            return is;
        }
    }

    // vector is a three-dimensional vector with cartesian coordinates.
    struct vector : point_detail::coordinates
    {
        using point_detail::coordinates::coordinates;

        // Component-wise addition.
        inline vector& operator+=(vector const& other) noexcept
        {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }

        // Component-wise subtraction.
        inline vector& operator-=(vector const& other) noexcept
        {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }

        // Scalar multiplication is broadcast to all components.
        inline vector& operator*=(double mult) noexcept
        {
            x *= mult;
            y *= mult;
            z *= mult;
            return *this;
        }

        // Scalar division is broadcast to all components.
        inline vector& operator/=(double divisor)
        {
            return *this *= 1 / divisor;
        }

        // dot returns the dot product of this vector and other.
        inline double dot(vector const& other) const noexcept
        {
            return (
                x * other.x +
                y * other.y +
                z * other.z
            );
        }

        // cross returns the cross product of this vector and rhs.
        inline vector cross(vector const& rhs) const noexcept
        {
            return vector {
                y * rhs.z - z * rhs.y,
                z * rhs.x - x * rhs.z,
                x * rhs.y - y * rhs.x
            };
        }

        // norm returns the Euclidean norm of this vector.
        inline double norm() const noexcept
        {
            return std::sqrt(squared_norm());
        }

        // squared_norm returns the squared Euclidean norm of this vector.
        inline double squared_norm() const noexcept
        {
            return dot(*this);
        }

        // normalize returns the unit vector that is parallel to this vector.
        // Behavior is undefined if this vector is zero.
        inline vector normalize() const
        {
            return vector{*this} /= norm();
        }
    };

    // Unary no-op `+vec` returns a copy of vec.
    inline vector operator+(vector const& vec) noexcept
    {
        return vec;
    }

    // Unary negation `-vec` returns an anti-parallel copy of vec.
    inline vector operator-(vector const& vec) noexcept
    {
        return {-vec.x, -vec.y, -vec.z};
    }

    // Sum `lhs + rhs` returns a vector with component-wise sum of lhs and rhs.
    inline vector operator+(vector const& lhs, vector const& rhs) noexcept
    {
        return vector{lhs} += rhs;
    }

    // Difference `lhs - rhs` returns a vector with component-wise difference of
    // lhs and rhs.
    inline vector operator-(vector const& lhs, vector const& rhs) noexcept
    {
        return vector{lhs} -= rhs;
    }

    // Multiplication `vec * mult` returns vec with each component scaled by
    // mult.
    inline vector operator*(vector const& vec, double mult) noexcept
    {
        return vector{vec} *= mult;
    }

    // Multiplication `mult * vec` returns vec with each component scaled by
    // mult.
    inline vector operator*(double mult, vector const& vec) noexcept
    {
        return vector{vec} *= mult;
    }

    // Division `vec / divisor` returns vec with each component scaled by
    // 1/divisor.
    inline vector operator/(vector const& vec, double divisor)
    {
        return vector{vec} /= divisor;
    }

    // dot returns the dot product of lhs and rhs.
    inline double dot(vector const& lhs, vector const& rhs) noexcept
    {
        return lhs.dot(rhs);
    }

    // cross returns the cross product of lhs and rhs.
    inline vector cross(vector const& lhs, vector const& rhs) noexcept
    {
        return lhs.cross(rhs);
    }

    // norm returns the Euclidean norm of vec.
    inline double norm(vector const& vec) noexcept
    {
        return vec.norm();
    }

    // squared_norm returns the squares Euclidean norm of vec.
    inline double squared_norm(vector const& vec) noexcept
    {
        return vec.squared_norm();
    }

    // normalize returns the unit vector that is parallel to vec.
    inline vector normalize(vector const& vec) noexcept
    {
        return vec.normalize();
    }

    // point is a cartesian representation of a point in the three-dimensional
    // Euclidean space.
    struct point : point_detail::coordinates
    {
        using point_detail::coordinates::coordinates;
        using vector_t = struct vector;

        // vector returns the coordinate vector of this point.
        inline vector_t vector() const noexcept
        {
            return vector_t{x, y, z};
        }

        // Adding a vector translates the point.
        inline point& operator+=(vector_t const& disp) noexcept
        {
            x += disp.x;
            y += disp.y;
            z += disp.z;
            return *this;
        }

        // Subtracting a vector translates the point.
        inline point& operator-=(vector_t const& disp) noexcept
        {
            x -= disp.x;
            y -= disp.y;
            z -= disp.z;
            return *this;
        }

        // distance returns the Euclidean distance between this point and other.
        inline double distance(point const& other) const noexcept
        {
            return (vector() - other.vector()).norm();
        }

        // squared_distance returns the squared Euclidean distance between this
        // point and other.
        inline double squared_distance(point const& other) const noexcept
        {
            return (vector() - other.vector()).squared_norm();
        }
    };

    // Sum `pt + disp` returns the point pt translated by disp.
    inline point operator+(point const& pt, vector const& disp) noexcept
    {
        return point{pt} += disp;
    }

    // Difference `pt - disp` returns the point pt translated by -disp.
    inline point operator-(point const& pt, vector const& disp) noexcept
    {
        return point{pt} -= disp;
    }

    // Difference `pa - pb` returns the displacement vector from pb towards pa.
    inline vector operator-(point const& pa, point const& pb) noexcept
    {
        return pa.vector() - pb.vector();
    }

    // distance returns the Euclidean distance between pa and pb.
    inline double distance(point const& pa, point const& pb) noexcept
    {
        return pa.distance(pb);
    }

    // squared_distance returns the squared Euclidean distance between pa and
    // pb.
    inline double squared_distance(point const& pa, point const& pb) noexcept
    {
        return pa.squared_distance(pb);
    }
}


//------------------------------------------------------------------------------
// typedefs
//------------------------------------------------------------------------------

#include <cstddef>


namespace md
{
    // index is the integral type used to index arrays.
    using index = std::size_t;

    // scalar is the floating-point type of choice.
    using scalar = double;
}


//------------------------------------------------------------------------------
// neighbor_searcher
//------------------------------------------------------------------------------

#include <algorithm>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>


namespace md
{
    // linear_hash is a linear hash function for integral 3-vectors.
    struct linear_hash
    {
        using hash_t = std::uint32_t;

        hash_t x_coeff = 3929498747;
        hash_t y_coeff = 1008281837;
        hash_t z_coeff = 1832832077;
        hash_t modulus = 1021;

        inline hash_t operator()(hash_t x, hash_t y, hash_t z) const
        {
            // Avoid 32-bit wraparound (mod 2^32) with 64-bit extension.
            using hash2x_t = std::uint64_t;

            hash2x_t sum = 0;
            sum += hash2x_t(x_coeff) * x;
            sum += hash2x_t(y_coeff) * y;
            sum += hash2x_t(z_coeff) * z;

            return hash_t(sum % modulus);
        }
    };


    // neighbor_searcher
    class neighbor_searcher
    {
    public:
        // Constructor takes the cut-off distance and a hash function.
        neighbor_searcher(md::scalar dcut, md::linear_hash hash)
            : dcut_{dcut}, hash_{hash}, buckets_(hash.modulus)
        {
            // Pre-compute neighborhood of each bucket for faster query.
            using hash_t = linear_hash::hash_t;

            hash_t const coord_deltas[] = {
                hash.modulus - 1,
                hash.modulus,
                hash.modulus + 1
            };
            std::set<hash_t> hash_deltas;

            for (hash_t const dx : coord_deltas) {
                for (hash_t const dy : coord_deltas) {
                    for (hash_t const dz : coord_deltas) {
                        hash_deltas.insert(hash_(dx, dy, dz));
                    }
                }
            }

            for (md::index center = 0; center < hash.modulus; center++) {
                std::vector<md::index>& neighbors = buckets_[center].neighbors;

                for (hash_t const delta : hash_deltas) {
                    hash_t const neighbor = (center + delta) % hash.modulus;

                    // Leverage symmetry.
                    if (neighbor >= center) {
                        neighbors.push_back(neighbor);
                    }
                }
            }
        }

        // set_points fills hash table.
        void set_points(md::array_view<md::point const> points)
        {
            for (hash_bucket& bucket : buckets_) {
                bucket.members.clear();
            }

            for (md::index idx = 0; idx < points.size(); idx++) {
                md::point const pt = points[idx];
                hash_bucket& bucket = buckets_[locate_bucket(pt)];
                bucket.members.push_back({idx, pt});
            }
        }

        // search
        template<typename OutputIterator>
        void search(OutputIterator out)
        {
            for (hash_bucket const& center : buckets_) {
                for (md::index const neighbor : center.neighbors) {
                    search_among(center, buckets_[neighbor], out);
                }
            }
        }

    private:
        struct hash_bucket
        {
            struct member
            {
                md::index index;
                md::point point;
            };
            std::vector<member> members;
            std::vector<md::index> neighbors;
        };

        template<typename OutputIterator>
        inline void search_among(
            hash_bucket const& bucket_a,
            hash_bucket const& bucket_b,
            OutputIterator& out
        )
        {
            md::scalar const dcut2 = dcut_ * dcut_;

            for (hash_bucket::member member_j : bucket_b.members) {
                for (hash_bucket::member member_i : bucket_a.members) {
                    if (member_i.index == member_j.index) {
                        // Avoid double counting in intra-bucket search.
                        break;
                    }

                    if (md::squared_distance(member_i.point, member_j.point) > dcut2) {
                        continue;
                    }

                    *out++ = std::make_pair(
                        std::min(member_i.index, member_j.index),
                        std::max(member_i.index, member_j.index)
                    );
                }
            }
        }

        // locate_bin returns the hash index for a point.
        inline linear_hash::hash_t locate_bucket(md::point pt) const
        {
            using hash_t = linear_hash::hash_t;

            // Offset coordinates to avoid negative value.
            constexpr md::scalar offset = 1L << 30;

            md::scalar const freq = 1 / dcut_;
            hash_t const x = hash_t(offset + freq * pt.x);
            hash_t const y = hash_t(offset + freq * pt.y);
            hash_t const z = hash_t(offset + freq * pt.z);

            return hash_(x, y, z);
        }

    private:
        md::scalar dcut_;
        md::linear_hash hash_;
        std::vector<hash_bucket> buckets_;
    };

    namespace detail
    {
        // compute_variance returns the population variance of given points.
        inline md::vector compute_variance(md::array_view<md::point const> points)
        {
            if (points.empty()) {
                return {};
            }

            md::scalar sum_x = 0;
            md::scalar sum_y = 0;
            md::scalar sum_z = 0;

            md::scalar sum_x2 = 0;
            md::scalar sum_y2 = 0;
            md::scalar sum_z2 = 0;

            for (md::point const& pt : points) {
                sum_x += pt.x;
                sum_y += pt.y;
                sum_z += pt.z;

                sum_x2 += pt.x * pt.x;
                sum_y2 += pt.y * pt.y;
                sum_z2 += pt.z * pt.z;
            }

            md::scalar const n = md::scalar(points.size());

            return md::vector {
                sum_x2 / n - (sum_x / n) * (sum_x / n),
                sum_y2 / n - (sum_y / n) * (sum_y / n),
                sum_z2 / n - (sum_z / n) * (sum_z / n)
            };
        }

        // determine_hash returns a linear_hash object that is heuristically
        // parameterized to make neighbor_searcher perform good on given points.
        inline md::linear_hash determine_hash(md::array_view<md::point const> points)
        {
            // Heuristic: Squash lowest-variance axis.
            md::vector const var = detail::compute_variance(points);
            md::scalar const min = std::min({var.x, var.y, var.z});

            md::linear_hash hash;

            if (var.x == min) {
                hash.x_coeff = 0;
            } else if (var.y == min) {
                hash.y_coeff = 0;
            } else {
                hash.z_coeff = 0;
            }

            // Heuristic: Good bucket occupancy.
            constexpr md::index bucket_occupancy = 16;

            hash.modulus = md::linear_hash::hash_t(points.size() / bucket_occupancy);
            hash.modulus |= 1;

            return hash;
        }
    }
}


#endif
