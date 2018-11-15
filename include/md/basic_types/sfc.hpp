// Copyright snsinfu 2018.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef MD_BASIC_TYPES_SFC_HPP
#define MD_BASIC_TYPES_SFC_HPP

// This module provides a small fast random engine.
// Source: https://github.com/snsinfu/cxx-sfcrandom

#include <cstdint>
#include <istream>
#include <limits>
#include <ostream>


namespace md
{
    namespace sfc_detail
    {
        // rotate returns the left rotation of x.
        template<typename T>
        inline T rotate(T x, int n)
        {
            // Avoid UB.
            if (n == 0) {
                return x;
            } else {
                return T(x << n) | T(x >> (std::numeric_limits<T>::digits - n));
            }
        }

        // state holds sfc generator state variables.
        template<typename T>
        struct state
        {
            T a = 0;
            T b = 0;
            T c = 0;
            T counter = 1;

            // init initializes state variables using given 64-bit seed.
            void init(std::uint64_t seed_val) noexcept;

            // init_seq initializes state variables using words read from given
            // sequence. Lower bits are taken first.
            template<typename SeedSeq>
            void init_seq(SeedSeq& seed_seq);
        };

        template<>
        inline void state<std::uint16_t>::init(std::uint64_t seed_val) noexcept
        {
            a = std::uint16_t(seed_val);
            b = std::uint16_t(seed_val >> 16);
            c = std::uint16_t(seed_val >> 32);
            counter = std::uint16_t(seed_val >> 48);
        }

        template<>
        inline void state<std::uint32_t>::init(std::uint64_t seed_val) noexcept
        {
            a = 0;
            b = std::uint32_t(seed_val);
            c = std::uint32_t(seed_val >> 32);
            counter = 1;
        }

        template<>
        inline void state<std::uint64_t>::init(std::uint64_t seed_val) noexcept
        {
            a = seed_val;
            b = seed_val;
            c = seed_val;
            counter = 1;
        }

        template<>
        template<typename SeedSeq>
        inline void state<std::uint16_t>::init_seq(SeedSeq& seed_seq)
        {
            std::uint32_t uints[2];
            seed_seq.generate(uints, uints + 2);

            a = std::uint16_t(uints[0]);
            b = std::uint16_t(uints[0] >> 16);
            c = std::uint16_t(uints[1]);
            counter = 1;
        }

        template<>
        template<typename SeedSeq>
        inline void state<std::uint32_t>::init_seq(SeedSeq& seed_seq)
        {
            std::uint32_t uints[3];
            seed_seq.generate(uints, uints + 3);

            a = uints[0];
            b = uints[1];
            c = uints[2];
            counter = 1;
        }

        template<>
        template<typename SeedSeq>
        inline void state<std::uint64_t>::init_seq(SeedSeq& seed_seq)
        {
            std::uint32_t uints[6];
            seed_seq.generate(uints, uints + 6);

            a = std::uint64_t(uints[0]) << 32 | uints[1];
            b = std::uint64_t(uints[2]) << 32 | uints[3];
            c = std::uint64_t(uints[4]) << 32 | uints[5];
            counter = 1;
        }

        // Equality comparison operator compares the member-wise equality.
        template<typename T>
        bool operator==(state<T> const& s1, state<T> const& s2) noexcept
        {
            return s1.a == s2.a
                && s1.b == s2.b
                && s1.c == s2.c
                && s1.counter == s2.counter;
        }

        template<typename T>
        bool operator!=(state<T> const& s1, state<T> const& s2) noexcept
        {
            return !(s1 == s2);
        }

        // stream_restorer saves and restores stream flags.
        template<typename Char, typename Tr>
        class stream_restorer
        {
        public:
            using char_type = Char;
            using stream_type = std::basic_ios<Char, Tr>;
            using flags_type = typename stream_type::fmtflags;

            stream_restorer(stream_restorer const&) = delete;
            stream_restorer& operator=(stream_restorer const&) = delete;

            // Constructor saves the stream flags.
            explicit stream_restorer(stream_type& stream)
                : stream_(stream), flags_(stream.flags()), fill_(stream.fill())
            {
            }

            // Destructor restores the saved stream flags.
            ~stream_restorer()
            {
                stream_.flags(flags_);
                stream_.fill(fill_);
            }

        private:
            stream_type& stream_;
            flags_type flags_;
            char_type fill_;
        };

        // Stream output operator writes state variables to a stream in decimal
        // format with white-space delimiter.
        template<typename Char, typename Tr, typename T>
        std::basic_ostream<Char, Tr>& operator<<(
            std::basic_ostream<Char, Tr>& os,
            state<T> const& st
        )
        {
            using sentry_type = typename std::basic_ostream<Char, Tr>::sentry;

            if (sentry_type sentry{os}) {
                Char const space = os.widen(' ');

                stream_restorer<Char, Tr> restorer{os};
                os.flags(std::ios::dec | std::ios::left);
                os.fill(space);

                os << st.a << space
                   << st.b << space
                   << st.c << space
                   << st.counter;
            }

            return os;
        }

        // Stream input operator reads state variables from a stream.
        template<typename Char, typename Tr, typename T>
        std::basic_istream<Char, Tr>& operator>>(
            std::basic_istream<Char, Tr>& is,
            state<T>& st
        )
        {
            using sentry_type = typename std::basic_istream<Char, Tr>::sentry;

            if (sentry_type sentry{is}) {
                state<T> state;

                stream_restorer<Char, Tr> restorer{is};
                is.flags(std::ios::dec | std::ios::skipws);

                is >> state.a
                   >> state.b
                   >> state.c
                   >> state.counter;

                st = state;
            }

            return is;
        }

        // check_seed_seq is instanciatable if T satisfies the SeedSequence
        // concept.
        template<typename T>
        auto check_seed_seq(
            std::uint32_t* rb = nullptr,
            std::uint32_t* re = nullptr
        ) -> decltype(std::declval<T>().generate(rb, re));
    }

    // small_fast_counting_engine_v4 implements the sfc16, sfc32 and sfc64
    // random number generators found in PractRand. It implements the C++11
    // RandomNumberEngine concept.
    template<
        typename T,
        int Rot,
        int RSh,
        int LSh,
        unsigned long long Rnd1,
        unsigned long long Rnd3
    >
    class small_fast_counting_engine_v4
    {
        using state_type = sfc_detail::state<T>;

        static constexpr int rotation = Rot;
        static constexpr int right_shift = RSh;
        static constexpr int left_shift = LSh;

        // Initial round count for seed/1.
        static constexpr unsigned long long initial_round_1 = Rnd1;

        // Initial round count for seed/3.
        static constexpr unsigned long long initial_round_3 = Rnd3;

    public:
        // result_type is an alias of T.
        using result_type = T;

        // default_seed is the default seed value used by default constructor
        // and seed function.
        static constexpr result_type default_seed = 0;

        // min returns 0.
        static constexpr result_type min() noexcept
        {
            return 0;
        }

        // max returns a result_type with all bits set.
        static constexpr result_type max() noexcept
        {
            return result_type(-1);
        }

        // Constructor creates a seeded engine. This overload extends seed_val
        // to 64-bit and seeds the engine with the PractRand seed algorithm.
        explicit small_fast_counting_engine_v4(result_type seed_val = default_seed) noexcept
        {
            seed(seed_val);
        }

        // Constructor creates a seeded engine. This overload reads three words
        // from given sequence and seeds the engine with the PractRand seed/3
        // algorithm.
        template<
            typename SeedSeq,
            typename = decltype(sfc_detail::check_seed_seq<SeedSeq>())
        >
        explicit small_fast_counting_engine_v4(SeedSeq& seed_seq)
        {
            seed(seed_seq);
        }

        // seed re-seeds the engine. This overload extends seed_val to 64-bit
        // and seeds the engine with the PractRand seed/1 algorithm.
        void seed(result_type seed_val = default_seed) noexcept
        {
            seed64(std::uint64_t{seed_val});
        }

        // seed re-seeds the engine. This overload reads three words from given
        // sequence and seeds the engine with the PractRand seed/3 algorithm.
        template<
            typename SeedSeq,
            typename = decltype(sfc_detail::check_seed_seq<SeedSeq>())
        >
        void seed(SeedSeq& seed_seq)
        {
            state_.init_seq(seed_seq);
            discard(initial_round_3);
        }

        // seed64 seeds the engine using given 64-bit value using the PractRand
        // seed/1 algorithm. This function is a non-standard extension.
        void seed64(std::uint64_t seed_val) noexcept
        {
            state_.init(seed_val);
            discard(initial_round_1);
        }

        // Invocation operator advances the engine state and returns the
        // generated value.
        inline result_type operator()() noexcept
        {
            auto& a = state_.a;
            auto& b = state_.b;
            auto& c = state_.c;
            auto& counter = state_.counter;

            auto const output = result_type(a + b + counter++);
            a = result_type(b ^ result_type(b >> right_shift));
            b = result_type(c + result_type(c << left_shift));
            c = result_type(sfc_detail::rotate(c, rotation) + output);

            return output;
        }

        // discard advances the engine state by z. The time complexity is O(z).
        void discard(unsigned long long z) noexcept
        {
            for (; z != 0; z--) {
                (*this)();
            }
        }

        // Equality comparison operator compares two engines for their states.
        friend bool operator==(
            small_fast_counting_engine_v4 const& e1,
            small_fast_counting_engine_v4 const& e2
        ) noexcept
        {
            return e1.state_ == e2.state_;
        }

        friend bool operator!=(
            small_fast_counting_engine_v4 const& e1,
            small_fast_counting_engine_v4 const& e2
        ) noexcept
        {
            return !(e1 == e2);
        }

        // Stream output operator writes the engine state to stream.
        template<typename Char, typename Tr>
        friend std::basic_ostream<Char, Tr>& operator<<(
            std::basic_ostream<Char, Tr>& os,
            small_fast_counting_engine_v4 const& engine
        )
        {
            return os << engine.state_;
        }

        // Stream input operator loads engine state from stream.
        template<typename Char, typename Tr>
        friend std::basic_istream<Char, Tr>& operator>>(
            std::basic_istream<Char, Tr>& is,
            small_fast_counting_engine_v4& engine
        )
        {
            return is >> engine.state_;
        }

    private:
        state_type state_;
    };

    // sfc16 is the standard 16-bit variant of small_fast_counting_engine_v4.
    // The parameters are taken from PractRand 0.94.
    using sfc16 = small_fast_counting_engine_v4<std::uint16_t, 6, 5, 3, 10, 10>;

    // sfc32 is the standard 32-bit variant of small_fast_counting_engine_v4.
    // The parameters are taken from PractRand 0.94.
    using sfc32 = small_fast_counting_engine_v4<std::uint32_t, 21, 9, 3, 12, 15>;

    // sfc64 is the standard 64-bit variant of small_fast_counting_engine_v4.
    // The parameters are taken from PractRand 0.94.
    using sfc64 = small_fast_counting_engine_v4<std::uint64_t, 24, 11, 3, 12, 18>;
}

#endif
