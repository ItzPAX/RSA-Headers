#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <bit>

class big_int
{
public:
    using limb_t = std::uint32_t;
    using dlimb_t = std::uint64_t;

    std::vector<limb_t> _v;
    std::size_t         _bit_len = 0;

public:
    static int compare(const big_int& a, const big_int& b)
    {
        if (a._bit_len != b._bit_len)
            return a._bit_len > b._bit_len ? 1 : -1;

        std::size_t n = a._v.size();
        for (std::size_t i = 0; i < n; ++i) {
            limb_t ai = a._v[n - 1 - i];
            limb_t bi = b._v[n - 1 - i];
            if (ai != bi) return ai > bi ? 1 : -1;
        }
        return 0;
    }

    bool operator==(const big_int& rhs) const { return compare(*this, rhs) == 0; }
    bool operator!=(const big_int& rhs) const { return compare(*this, rhs) != 0; }

    bool operator< (const big_int& rhs) const { return compare(*this, rhs) < 0; }
    bool operator<=(const big_int& rhs) const { return compare(*this, rhs) <= 0; }
    bool operator> (const big_int& rhs) const { return compare(*this, rhs) > 0; }
    bool operator>=(const big_int& rhs) const { return compare(*this, rhs) >= 0; }

    std::uint32_t mod_uint32(std::uint32_t m) const
    {
        if (m == 0) throw std::invalid_argument("mod_uint32: divide by zero");
        if (_bit_len == 0) return 0;

        std::uint64_t r = 0;
        for (std::size_t i = _v.size(); i-- > 0; ) {
            r = (r << 32) + _v[i];                // build 64-bit word
            r %= m;                               // always < m afterwards
        }
        return static_cast<std::uint32_t>(r);
    }

    // O(N^2) IMPROVE IMPROVE IMPROVE
    std::string dec_to_bin(const std::string& dec) const
    {
        if (dec == "0") return "0";

        std::string num = dec;
        std::string bits;

        while (!(num.size() == 1 && num[0] == '0')) {
            int carry = 0;
            std::string next;

            for (char c : num) {
                int cur = carry * 10 + (c - '0');
                char q = static_cast<char>('0' + cur / 2);
                if (!next.empty() || q != '0')
                    next.push_back(q);
                carry = cur % 2;
            }
            bits.push_back(static_cast<char>('0' + carry));
            num = next.empty() ? "0" : next;
        }
        std::reverse(bits.begin(), bits.end());
        return bits;
    }

    void shl1()
    {
        if (_v.empty())
            _v.push_back(0);

        limb_t carry = 0;
        for (std::size_t i = 0; i < _v.size(); ++i) {
            limb_t new_carry = _v[i] >> 31;
            _v[i] = (_v[i] << 1) | carry;
            carry = new_carry;
        }
        if (carry) _v.push_back(carry);
        ++_bit_len;
    }

    void shr1()
    {
        if (_v.empty()) return;
        limb_t carry = 0;
        for (std::size_t i = _v.size(); i-- > 0;) {
            limb_t new_carry = _v[i] & 1;
            _v[i] = (_v[i] >> 1) | (carry << 31);
            carry = new_carry;
        }
        trim();
    }

    bool is_odd() const
    {
        return _bit_len && (_v[0] & 1u);
    }

    bool is_zero() const
    {
        return _bit_len == 0;
    }

    bool bit(std::size_t i) const
    {
        return (i < _bit_len) && ((_v[i >> 5] >> (i & 31)) & 1u);
    }

    void trim()
    {
        while (!_v.empty() && _v.back() == 0) _v.pop_back();
        _bit_len = _v.empty() ? 0
            : (_v.size() - 1) * 32 + (32 - std::countl_zero(_v.back()));
    }

public:
    big_int()
        : _bit_len(0)
    {

    }

    big_int(std::vector<limb_t> _)
    {
        _bit_len = _.empty() ? 0 :
            (_.size() - 1) * sizeof(limb_t) +
            std::bit_width(_.back());

        _v.resize(_.size());
        memcpy(_v.data(), _.data(), _.size() * sizeof(limb_t));
    }

    big_int(int dec)
    {
        if (dec < 0)
            throw std::invalid_argument("big_int(int): negative value not supported");

        if (dec == 0) {
            _bit_len = 0;
            return;
        }

        limb_t v = static_cast<limb_t>(dec);
        _v.push_back(v);

        _bit_len = (_v.size() - 1) * 32 +
            std::bit_width(_v.back());
    }

    big_int(std::string dec)
        : _bit_len(0)
    {
        std::string bin = dec_to_bin(dec);
        _bit_len = bin.size();

        _v.assign((_bit_len + 31) / 32, 0);

        std::size_t bit_pos = 0;
        for (auto it = bin.rbegin(); it != bin.rend(); ++it, ++bit_pos)
            if (*it == '1')
                _v[bit_pos >> 5]
                |= static_cast<limb_t>(1u) << (bit_pos & 31);
    }

public:
    std::string as_string() const
    {
        if (_bit_len == 0) return "0";

        std::string raw;
        raw.reserve(_bit_len);

        std::size_t full_limbs = _bit_len >> 5;
        std::size_t extra_bits = _bit_len & 31;

        if (extra_bits) {
            limb_t limb = _v[full_limbs];
            for (int b = static_cast<int>(extra_bits) - 1; b >= 0; --b)
                raw.push_back((limb & (1u << b)) ? '1' : '0');
        }

        for (std::size_t i = 0; i < full_limbs; ++i) {
            limb_t limb = _v[full_limbs - 1 - i];
            for (int b = 31; b >= 0; --b)
                raw.push_back((limb & (1u << b)) ? '1' : '0');
        }

        std::size_t pad = (4 - (raw.size() % 4)) & 3;
        std::string out;
        out.reserve(raw.size() + raw.size() / 4 + pad);

        out.append(pad, '0');

        std::size_t nibble_cnt = pad;
        for (char c : raw) {
            if (nibble_cnt == 4) {
                out.push_back(' ');
                nibble_cnt = 0;
            }
            out.push_back(c);
            ++nibble_cnt;
        }

        return out;
    }


public:
    big_int powmod(big_int exp, const big_int& mod)
    {
        if (mod._bit_len == 0)
            throw std::invalid_argument("powmod : modulus is zero");

        *this = *this % mod;
        big_int result("1");

        while (exp._bit_len) {
            if (exp.is_odd()) {
                result = (result * *this) % mod;
            }
            exp.shr1();
            if (exp._bit_len)
                *this = (*this * *this) % mod;
        }
        return result;
    }

    big_int operator+(const big_int& rhs) const
    {
        const std::size_t la = _v.size();
        const std::size_t lb = rhs._v.size();
        const std::size_t n = std::max(la, lb);

        std::vector<limb_t> sum;
        sum.resize(n);

        dlimb_t carry = 0;
        for (std::size_t i = 0; i < n; ++i) {
            dlimb_t a = (i < la) ? _v[i] : 0;
            dlimb_t b = (i < lb) ? rhs._v[i] : 0;

            dlimb_t t = a + b + carry;
            sum[i] = static_cast<limb_t>(t);
            carry = t >> 32;
        }
        if (carry) sum.push_back(static_cast<limb_t>(carry));

        big_int s(sum);
        s.trim();
        return s;
    }

    big_int operator-(const big_int& rhs) const
    {
        if (compare(*this, rhs) < 0)
            throw std::invalid_argument("big_int::operator- : negative result");

        big_int diff; diff._v.resize(_v.size());
        limb_t borrow = 0;

        for (std::size_t i = 0; i < _v.size(); ++i) {
            limb_t a = _v[i];
            limb_t b = (i < rhs._v.size()) ? rhs._v[i] : 0;

            limb_t d = a - b - borrow;
            borrow = (borrow ? a <= b : a < b);
            diff._v[i] = d;
        }
        diff.trim();
        return diff;
    }

    // O(n^2) OPTIMIZE ME
    big_int operator*(const big_int& rhs) const
    {
        std::size_t la = _v.size(), lb = rhs._v.size();
        std::vector<limb_t> prod(la + lb, 0);

        for (std::size_t i = 0; i < la; ++i) {
            dlimb_t carry = 0;
            for (std::size_t j = 0; j < lb; ++j) {
                dlimb_t cur = prod[i + j] +
                    static_cast<dlimb_t>(_v[i]) * rhs._v[j] + carry;
                prod[i + j] = static_cast<limb_t>(cur);
                carry = cur >> 32;
            }
            prod[i + lb] += static_cast<limb_t>(carry);
        }
        big_int r; r._v = std::move(prod); r.trim(); return r;
    }

    big_int operator/(const big_int& rhs) const
    {
        if (rhs._bit_len == 0)
            throw std::invalid_argument("big_int::operator/ : divide by zero");

        if (compare(*this, rhs) < 0)
            return big_int();

        big_int rem;
        big_int quo;

        quo._bit_len = 0;

        for (std::size_t i = _bit_len; i-- > 0; )
        {
            rem.shl1();

            std::size_t limb_idx = i >> 5;
            std::size_t bit_idx = i & 31;
            if ((_v[limb_idx] >> bit_idx) & 1u)
                rem._v[0] |= 1u;
            rem.trim();

            quo.shl1();

            if (compare(rem, rhs) >= 0) {
                rem = rem - rhs;
                rem.trim();
                quo._v[0] |= 1u;
            }
        }

        quo.trim();
        return quo;
    }

    big_int operator%(const big_int& rhs) const
    {
        if (rhs._bit_len == 0)
            throw std::invalid_argument("big_int::operator% : divide by zero");

        if (compare(*this, rhs) < 0)
            return *this;

        big_int rem;

        for (std::size_t i = _bit_len; i-- > 0; )
        {
            rem.shl1();

            std::size_t limb_idx = i >> 5;
            std::size_t bit_idx = i & 31;
            if ((_v[limb_idx] >> bit_idx) & 1u)
                rem._v[0] |= 1u;
            rem.trim();

            if (compare(rem, rhs) >= 0) {
                rem = rem - rhs;
                rem.trim();
            }
        }

        return rem;
    }

    big_int operator>>(std::size_t k) const
    {
        if (k == 0)               return *this;
        if (k >= _bit_len)        return big_int();

        std::size_t limb_shift = k >> 5;
        std::size_t bit_shift = k & 31;

        big_int r;
        std::size_t new_limbs = _v.size() - limb_shift;
        r._v.resize(new_limbs);

        if (bit_shift == 0) {
            for (std::size_t i = 0; i < new_limbs; ++i)
                r._v[i] = _v[i + limb_shift];
        }
        else {
            const std::size_t inv = 32 - bit_shift;
            for (std::size_t i = 0; i < new_limbs; ++i) {
                std::uint32_t low = _v[i + limb_shift] >> bit_shift;
                std::uint32_t high = (i + limb_shift + 1 < _v.size())
                    ? _v[i + limb_shift + 1] << inv
                    : 0;
                r._v[i] = low | high;
            }
        }

        r.trim();
        return r;
    }

    big_int operator<<(std::size_t k) const
    {
        if (k == 0 || _bit_len == 0) return *this;

        std::size_t limb_shift = k >> 5;
        std::size_t bit_shift = k & 31;

        big_int r;
        r._v.assign(_v.size() + limb_shift + 1, 0);

        if (bit_shift == 0) {
            for (std::size_t i = 0; i < _v.size(); ++i)
                r._v[i + limb_shift] = _v[i];
        }
        else {
            const std::size_t inv = 32 - bit_shift;
            std::uint32_t carry = 0;
            for (std::size_t i = 0; i < _v.size(); ++i) {
                std::uint32_t cur = _v[i];
                r._v[i + limb_shift] = (cur << bit_shift) | carry;
                carry = cur >> inv;
            }
            r._v[_v.size() + limb_shift] = carry;
        }

        r._bit_len = _bit_len + k;
        r.trim();
        return r;
    }
};

namespace barrett
{
    struct barrett_ctx
    {
        big_int n;
        big_int mu;
        std::size_t k;
    };

    barrett_ctx make_ctx(const big_int& n)
    {
        barrett_ctx c{ n };
        c.k = n._v.size();

        big_int b2k;
        b2k._v.assign(c.k * 2 + 1, 0);
        b2k._v.back() = 1;
        b2k._bit_len = (c.k * 2 + 1) * 32;

        c.mu = b2k / n;
        return c;
    }

    inline big_int barrett_red(big_int x, const barrett_ctx& c)
    {
        if (x < c.n) return x;

        if (c.k == 0) return x % c.n;

        const std::size_t shift1 = 32 * (c.k - 1);
        const std::size_t shift2 = 32 * (c.k + 1);

        big_int q1 = x >> shift1;
        big_int q2 = q1 * c.mu;
        big_int q3 = q2 >> shift2;

        big_int r = x - q3 * c.n;
        while (r >= c.n) r = r - c.n;      

        return r;
    }

    big_int f_powmod(big_int base,
        const big_int& exp,
        const barrett_ctx& ctx)
    {
        const int W = 4;
        const int TABLE = 1 << (W - 1);

        base = barrett_red(base, ctx);
        big_int one(1);

        big_int table[TABLE];
        table[0] = base;                        
        big_int base2 = barrett_red(base * base, ctx);
        for (int i = 1; i < TABLE; ++i)
            table[i] = barrett_red(table[i - 1] * base2, ctx);

        big_int acc = one;
        std::size_t bitpos = exp._bit_len;

        while (bitpos)
        {
            if (!exp.bit(bitpos - 1)) {
                acc = barrett_red(acc * acc, ctx);
                --bitpos;
                continue;
            }

            std::size_t width = std::min<std::size_t>(W, bitpos);
            std::size_t val = 0;
            for (std::size_t i = 0; i < width; ++i)
                val = (val << 1) | exp.bit(bitpos - 1 - i);

            while ((val & 1) == 0) {
                val >>= 1;
                --width;
            }

            for (std::size_t i = 0; i < width; ++i)
                acc = barrett_red(acc * acc, ctx);

            acc = barrett_red(acc * table[(val - 1) >> 1], ctx);

            bitpos -= width;
        }
        return acc;
    }
}

namespace big_int_tests
{
    bool run_tests()
    {
        std::cout << "Testing addition of 128 + 27...\n";
        big_int a1("128");
        big_int a2("27");
        std::cout << "Result: " << (a1 + a2).as_string() << " (expected: 1001 1011)\n\n";

        std::cout << "Testing subtraction of 4096 - 255...\n";
        big_int s1("4096");
        big_int s2("255");
        std::cout << "Result: " << (s1 - s2).as_string() << " (expected: 1111 0000 0001)\n\n";

        std::cout << "Testing multiplication of 510 * 258...\n";
        big_int m1("510");
        big_int m2("258");
        std::cout << "Result: " << (m1 * m2).as_string() << " (expected: 0010 0000 0001 1111 1100)\n\n";

        std::cout << "Testing division of 123456789012345678901234567890 / 98765432109876543210...\n";
        big_int d1("123456789012345678901234567890");
        big_int d2("98765432109876543210");
        std::cout << "Result: " << (d1 / d2).as_string() << " (expected: 0100 1010 1000 0001 0111 1100 0111 0100)\n\n";

        std::cout << "Testing modulo of 619 % 600...\n";
        big_int mod1("619");
        big_int mod2("600");
        std::cout << "Result: " << (mod1 % mod2).as_string() << " (expected: 0001 0011)\n\n";

        std::cout << "Testing powmod of (5113 ^ 252) % 21626...\n";
        big_int t_n("21626");
        big_int t_e("252");
        big_int t_msg("5113");
        std::cout << "Result: " << t_msg.powmod(t_e, t_n).as_string() << " (expected: 0100 1001 1100 1011)\n\n";

        std::cout << "Testing using powmod with big numbers...\n";
        big_int begin("123456734567890123456789012345678901678901234567890123456678901234567890123456789012345678901234567832109876546543210987654321098765432132109876546543210987654321098765432913579135791357913579135791357913579135791357913579135710987654321098709876543210987901234567890123456789012345678901234567890123456789012345678901234");
        big_int begin2("98765431032109876546543210987654321098765432109876543210987654321321098765491357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357654321098765432109876543210987654321098709876543210321098765465432109876543210987654321098765432109879876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876543210987654321098765432109876");
        big_int n = begin * begin2 * begin2;
        big_int e("65537");
        big_int m("135791125125357191357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913579135791357913");

        std::cout << "n bits : " << n._bit_len << '\n';
        std::cout << "e bits : " << e._bit_len << '\n';
        std::cout << "m bits : " << m._bit_len << '\n';
        std::cout << "Powmod (m^e mod n) bit length: " << m.powmod(e, n)._bit_len << "\n\n";

        std::cout << "Testing right shift (>>) of 128 >> 4...\n";
        big_int sh1("128");
        std::cout << "sh1: " << sh1.as_string() << " result: " << (sh1 >> 4).as_string() << "\n\n";

        std::cout << "Testing left shift (<<) of 128 << 12...\n";
        std::cout << "sh1: " << sh1.as_string() << " result: " << (sh1 << 12).as_string() << "\n\n";

        return true;
    }
}