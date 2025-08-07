#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <bit>

class big_int
{
public:
    std::vector<std::uint8_t> _value;   // little-endian: _value[0] = least-significant byte
    std::size_t               _bit_len; // total number of useful bits

public:
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

    static int compare(const big_int& a, const big_int& b)
    {
        if (a._bit_len != b._bit_len)
            return a._bit_len > b._bit_len ? 1 : -1;

        std::size_t n = (a._bit_len + 7) >> 3;
        for (std::size_t i = 0; i < n; ++i) {
            std::uint8_t ai = a._value[n - 1 - i];
            std::uint8_t bi = b._value[n - 1 - i];
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

    void shl1()
    {
        if (_value.empty())
            _value.push_back(0);

        std::uint16_t carry = 0;
        for (std::size_t i = 0; i < _value.size(); ++i) {
            std::uint16_t cur = static_cast<std::uint16_t>(_value[i]) << 1 | carry;
            _value[i] = static_cast<std::uint8_t>(cur & 0xFF);
            carry = cur >> 8;
        }
        if (carry) _value.push_back(static_cast<std::uint8_t>(carry));
        ++_bit_len;
    }

    void shr1()
    {
        if (_bit_len == 0) return;

        std::uint8_t carry = 0;
        for (std::size_t i = _value.size(); i-- > 0; ) {
            std::uint8_t cur = _value[i];
            std::uint8_t nextC = cur & 1u;
            _value[i] = (cur >> 1) | (carry << 7);
            carry = nextC;
        }
        trim();                         // updates _bit_len and pops MS zeros
    }

    bool is_odd() const
    {
        return _bit_len && (_value[0] & 1u);
    }

    void trim()
    {
        while (!_value.empty() && _value.back() == 0)
            _value.pop_back();

        if (_value.empty()) {
            _bit_len = 0;
            return;
        }
        _bit_len = (_value.size() - 1) * 8 +
            std::bit_width(_value.back());
    }

public:
    big_int()
        : _bit_len(0)
    {

    }

    big_int(std::vector<uint8_t> _, size_t bits)
        : _bit_len(bits)
    {
        _value.resize(_.size() * sizeof(uint8_t));
        memcpy(_value.data(), _.data(), _.size() * sizeof(uint8_t));
    }

    big_int(int dec)
    {
        if (dec < 0)
            throw std::invalid_argument("BigInt(int): negative value not supported");

        if (dec == 0) {                 // the canonical zero
            _bit_len = 0;
            return;
        }

        std::uint32_t v = static_cast<std::uint32_t>(dec);  // promote to unsigned
        while (v) {
            _value.push_back(static_cast<std::uint8_t>(v & 0xFF)); // LSB first
            v >>= 8;
        }

        _bit_len = (_value.size() - 1) * 8 +
            std::bit_width(_value.back());
    }

    big_int(std::string dec)
        : _bit_len(0)
    {
        std::string bin = dec_to_bin(dec);
        _bit_len = bin.size();
        _value.assign((_bit_len + 7) / 8, 0);

        std::size_t bit_pos = 0;
        for (auto it = bin.rbegin(); it != bin.rend(); ++it, ++bit_pos)
            if (*it == '1')
                _value[bit_pos >> 3] |= 1u << (bit_pos & 7);
    }

public:
    std::string as_string() const
    {
        if (_bit_len == 0) return "0";

        std::string raw;

        std::size_t full_bytes = _bit_len >> 3;
        std::size_t extra_bits = _bit_len & 7;

        if (extra_bits) {
            std::uint8_t msb = _value[full_bytes];
            for (int b = static_cast<int>(extra_bits) - 1; b >= 0; --b)
                raw.push_back((msb & (1u << b)) ? '1' : '0');
        }

        for (std::size_t i = 0; i < full_bytes; ++i) {
            std::uint8_t byte = _value[full_bytes - 1 - i];
            for (int b = 7; b >= 0; --b)
                raw.push_back((byte & (1u << b)) ? '1' : '0');
        }

        std::size_t first_group = raw.size() % 4;
        if (first_group == 0) first_group = 4;

        std::string out;
        out.reserve(raw.size() + raw.size() / 4 + raw.size() % 4);

        int pad = 4 - raw.size() % 4 == 4 ? 0 : 4 - raw.size() % 4;
        for (int i = 0; i < pad; i++)
        {
            out.push_back('0');
        }

        std::size_t count = 0;
        for (char c : raw) {
            if (count == first_group) {
                out.push_back(' ');
                count = 0;
                first_group = 4;
            }
            out.push_back(c);
            ++count;
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
        const std::size_t len_a = (_bit_len + 7) >> 3;
        const std::size_t len_b = (rhs._bit_len + 7) >> 3;
        const std::size_t max_len = std::max(len_a, len_b);

        std::vector<std::uint8_t> sum;
        sum.reserve(max_len + 1);

        std::uint16_t carry = 0;

        for (std::size_t i = 0; i < max_len; ++i)
        {
            std::uint16_t a = (i < len_a) ? _value[i] : 0;
            std::uint16_t b = (i < len_b) ? rhs._value[i] : 0;
            std::uint16_t s = a + b + carry;

            sum.push_back(static_cast<std::uint8_t>(s));

            carry = s >> 8;
        }
        if (carry) sum.push_back(static_cast<std::uint8_t>(carry));

        std::size_t bit_len = sum.empty() ? 0 :
            (sum.size() - 1) * 8 +
            std::bit_width(sum.back());

        return big_int(std::move(sum), bit_len);
    }

    big_int operator-(const big_int& rhs) const
    {
        if (compare(*this, rhs) < 0)
            throw std::invalid_argument("big_int::operator- : negative result");

        const std::size_t len_a = (_bit_len + 7) >> 3;
        const std::size_t len_b = (rhs._bit_len + 7) >> 3;

        std::vector<std::uint8_t> diff;
        diff.reserve(len_a);

        int borrow = 0;
        for (std::size_t i = 0; i < len_a; ++i) {
            int a = _value[i];
            int b = (i < len_b) ? rhs._value[i] : 0;

            int d = a - b - borrow;
            if (d < 0) { d += 256; borrow = 1; }
            else { borrow = 0; }

            diff.push_back(static_cast<std::uint8_t>(d));
        }

        while (!diff.empty() && diff.back() == 0)
            diff.pop_back();

        std::size_t bit_len = 0;
        if (!diff.empty())
            bit_len = (diff.size() - 1) * 8 +
            std::bit_width(diff.back());

        return big_int(std::move(diff), bit_len);
    }

    // O(n^2) OPTIMIZE ME
    big_int operator*(const big_int& rhs) const
    {
        if (_bit_len == 0 || rhs._bit_len == 0)
            return big_int();

        const std::size_t len_a = (_bit_len + 7) >> 3;
        const std::size_t len_b = (rhs._bit_len + 7) >> 3;

        std::vector<std::uint8_t> prod(len_a + len_b, 0);

        for (std::size_t i = 0; i < len_a; ++i) {
            std::uint16_t carry = 0;
            for (std::size_t j = 0; j < len_b; ++j) {
                std::uint32_t cur =
                    prod[i + j] +
                    static_cast<std::uint32_t>(_value[i]) *
                    static_cast<std::uint32_t>(rhs._value[j]) +
                    carry;

                prod[i + j] = static_cast<std::uint8_t>(cur & 0xFF);
                carry = cur >> 8;
            }
            std::size_t k = i + len_b;
            while (carry) {
                std::uint32_t cur = prod[k] + carry;
                prod[k] = static_cast<std::uint8_t>(cur & 0xFF);
                carry = cur >> 8;
                ++k;
            }
        }

        while (!prod.empty() && prod.back() == 0)
            prod.pop_back();

        std::size_t bit_len = 0;
        if (!prod.empty())
            bit_len = (prod.size() - 1) * 8 +
            std::bit_width(prod.back());

        return big_int(std::move(prod), bit_len);
    }

    big_int operator/(const big_int& rhs) const
    {
        if (rhs._bit_len == 0)
            throw std::invalid_argument("big_int::operator/ : divide by zero");

        if (compare(*this, rhs) < 0)
            return big_int();

        big_int rem;
        big_int quo;

        for (std::size_t i = _bit_len; i-- > 0; )
        {
            rem.shl1();
            std::uint8_t cur = _value[i >> 3];
            if ((cur >> (i & 7)) & 1)
                rem._value[0] |= 1;
            rem.trim();

            quo.shl1();

            if (compare(rem, rhs) >= 0) {             
                rem = rem - rhs;
                rem.trim();
                quo._value[0] |= 1;
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
            return *this;                          // |dividend| < |divisor|

        big_int rem;                                // running remainder (0)

        for (std::size_t i = _bit_len; i-- > 0; )
        {
            rem.shl1();
            std::uint8_t cur_byte = _value[i >> 3];
            if ((cur_byte >> (i & 7)) & 1)
                rem._value[0] |= 1;
            rem.trim();

            if (compare(rem, rhs) >= 0) {
                rem = rem - rhs;
                rem.trim();
            }
        }
        return rem;
    }
};

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
        std::cout << "Powmod (m^e mod n) bit length: " << m.powmod(e, n)._bit_len;

        return true;
    }
}