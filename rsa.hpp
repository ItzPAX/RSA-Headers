#pragma once

// https://github.com/czkz/base64/blob/master/base64.h
namespace base64
{
	[[nodiscard]] inline std::string to_base64(std::string_view data);
	[[nodiscard]] inline std::string to_base64(const void* data, size_t size_bytes);

	[[nodiscard]] inline std::string from_base64(std::string_view data);
	[[nodiscard]] inline std::string from_base64(const void* data, size_t size_bytes);

	inline void from_base64_inplace(std::string& data);
	inline size_t from_base64_inplace(void* data, size_t size_bytes);

	static constexpr uint8_t base64_lut[256][7] = {
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xf8\x03\xe0\x0f\x80\x3e",
		"\xff\xff\xff\xff\xff\xff", "\xf8\x03\xe0\x0f\x80\x3e", "\xff\xff\xff\xff\xff\xff", "\xfc\x03\xf0\x0f\xc0\x3f",
		"\xd0\x03\x40\x0d\x00\x34", "\xd4\x03\x50\x0d\x40\x35", "\xd8\x03\x60\x0d\x80\x36", "\xdc\x03\x70\x0d\xc0\x37",
		"\xe0\x03\x80\x0e\x00\x38", "\xe4\x03\x90\x0e\x40\x39", "\xe8\x03\xa0\x0e\x80\x3a", "\xec\x03\xb0\x0e\xc0\x3b",
		"\xf0\x03\xc0\x0f\x00\x3c", "\xf4\x03\xd0\x0f\x40\x3d", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\x00\x00\x00\x00\x00\x00", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\x00\x00\x00\x00\x00\x00", "\x04\x00\x10\x00\x40\x01", "\x08\x00\x20\x00\x80\x02",
		"\x0c\x00\x30\x00\xc0\x03", "\x10\x00\x40\x01\x00\x04", "\x14\x00\x50\x01\x40\x05", "\x18\x00\x60\x01\x80\x06",
		"\x1c\x00\x70\x01\xc0\x07", "\x20\x00\x80\x02\x00\x08", "\x24\x00\x90\x02\x40\x09", "\x28\x00\xa0\x02\x80\x0a",
		"\x2c\x00\xb0\x02\xc0\x0b", "\x30\x00\xc0\x03\x00\x0c", "\x34\x00\xd0\x03\x40\x0d", "\x38\x00\xe0\x03\x80\x0e",
		"\x3c\x00\xf0\x03\xc0\x0f", "\x40\x01\x00\x04\x00\x10", "\x44\x01\x10\x04\x40\x11", "\x48\x01\x20\x04\x80\x12",
		"\x4c\x01\x30\x04\xc0\x13", "\x50\x01\x40\x05\x00\x14", "\x54\x01\x50\x05\x40\x15", "\x58\x01\x60\x05\x80\x16",
		"\x5c\x01\x70\x05\xc0\x17", "\x60\x01\x80\x06\x00\x18", "\x64\x01\x90\x06\x40\x19", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xfc\x03\xf0\x0f\xc0\x3f",
		"\xff\xff\xff\xff\xff\xff", "\x68\x01\xa0\x06\x80\x1a", "\x6c\x01\xb0\x06\xc0\x1b", "\x70\x01\xc0\x07\x00\x1c",
		"\x74\x01\xd0\x07\x40\x1d", "\x78\x01\xe0\x07\x80\x1e", "\x7c\x01\xf0\x07\xc0\x1f", "\x80\x02\x00\x08\x00\x20",
		"\x84\x02\x10\x08\x40\x21", "\x88\x02\x20\x08\x80\x22", "\x8c\x02\x30\x08\xc0\x23", "\x90\x02\x40\x09\x00\x24",
		"\x94\x02\x50\x09\x40\x25", "\x98\x02\x60\x09\x80\x26", "\x9c\x02\x70\x09\xc0\x27", "\xa0\x02\x80\x0a\x00\x28",
		"\xa4\x02\x90\x0a\x40\x29", "\xa8\x02\xa0\x0a\x80\x2a", "\xac\x02\xb0\x0a\xc0\x2b", "\xb0\x02\xc0\x0b\x00\x2c",
		"\xb4\x02\xd0\x0b\x40\x2d", "\xb8\x02\xe0\x0b\x80\x2e", "\xbc\x02\xf0\x0b\xc0\x2f", "\xc0\x03\x00\x0c\x00\x30",
		"\xc4\x03\x10\x0c\x40\x31", "\xc8\x03\x20\x0c\x80\x32", "\xcc\x03\x30\x0c\xc0\x33", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
		"\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff", "\xff\xff\xff\xff\xff\xff",
	};

	[[nodiscard]] inline std::string to_base64(std::string_view data) {
		constexpr std::string_view a = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
		const std::basic_string_view<uint8_t> s{ reinterpret_cast<const uint8_t*>(data.data()), data.size() };
		std::string ret;
		const size_t origLen = s.size();
		const size_t newLen = (origLen / 3 + (origLen % 3 != 0)) * 4;
		ret.resize(newLen, '=');
		char* out = ret.data();

		size_t i; // Signed because (origLen - 2) can underflow
		for (i = 0; i < origLen - 2; i += 3) {
			out[i / 3 * 4 + 0] = a[(s[i] & 0b11111100) >> 2];
			out[i / 3 * 4 + 1] = a[(s[i + 0] & 0b00000011) << 4 | (s[i + 1] & 0b11110000) >> 4];
			out[i / 3 * 4 + 2] = a[(s[i + 1] & 0b00001111) << 2 | (s[i + 2] & 0b11000000) >> 6];
			out[i / 3 * 4 + 3] = a[(s[i + 2] & 0b00111111)];
		}
		switch (origLen % 3) {

		case 2: {                                         // TWO bytes remain
			out[i / 3 * 4 + 0] = a[(s[i] & 0xFC) >> 2];
			out[i / 3 * 4 + 1] = a[(s[i] & 0x03) << 4 | (s[i + 1] & 0xF0) >> 4];
			out[i / 3 * 4 + 2] = a[(s[i + 1] & 0x0F) << 2];     // last 2 bits are 0
			/* out[i/3*4+3] already '=' */
			break;
		}

		case 1: {                                         // ONE byte remains
			out[i / 3 * 4 + 0] = a[(s[i] & 0xFC) >> 2];
			out[i / 3 * 4 + 1] = a[(s[i] & 0x03) << 4];     // last 4 bits are 0
			/* out[i/3*4+2] and +3 stay '=' */
			break;
		}
		}
		return ret;
	}


	[[nodiscard]] inline std::string to_base64(const void* data, size_t size_bytes) {
		return to_base64({ reinterpret_cast<const char*>(data), size_bytes });
	}


	/// dstLen must be at least (3/4 * srcLen), even if there's less data encoded
	/// dst may be equal to src
	inline size_t from_base64(const void* src, size_t srcLen, void* dst, size_t dstLen) {
		std::basic_string_view<uint8_t> data{ reinterpret_cast<const uint8_t*>(src), srcLen };
		const size_t origLen = data.size();
		if (origLen == 0) {
			return 0;
		}
		if (origLen % 4 != 0) {
			return 0;
		}
		const size_t newLen = origLen / 4 * 3;
		if (dstLen < newLen) {
			return 0;
		}

		size_t outLen = newLen;
		uint8_t* out = reinterpret_cast<uint8_t*>(dst);

		for (size_t i = 0; i < origLen; i += 4) {
			out[i / 4 * 3 + 0] = base64_lut[data[i + 0]][0] | base64_lut[data[i + 1]][1];
			out[i / 4 * 3 + 1] = base64_lut[data[i + 1]][2] | base64_lut[data[i + 2]][3];
			out[i / 4 * 3 + 2] = base64_lut[data[i + 2]][4] | base64_lut[data[i + 3]][5];
		}
		// Min possible data.size() is 4, so rbegin() + 3 <= rend()
		// Min outLen is 3, so outLen-- will never underflow
		for (auto it = data.rbegin(); it != data.rbegin() + 3; ++it) {
			if (*it == '=') {
				outLen--;
			}
			else {
				break;
			}
		}
		return outLen;
	}

	inline size_t from_base64_inplace(void* data, size_t size_bytes) {
		return from_base64(data, size_bytes, data, size_bytes);
	}

	inline void from_base64_inplace(std::string& data) {
		auto newLen = from_base64(data.data(), data.size(), data.data(), data.size());
		data.resize(newLen);
	}

	[[nodiscard]] inline std::string from_base64(const void* data, size_t size_bytes) {
		std::string ret;
		ret.resize(size_bytes / 4 * 3);
		auto newLen = from_base64(data, size_bytes, ret.data(), ret.size());
		ret.resize(newLen);
		return ret;
	}

	[[nodiscard]] inline std::string from_base64(std::string_view data) {
		return from_base64(data.data(), data.size());
	}
}

#include <fstream>
#include <sstream>

#include "bigint.hpp"

#include <windows.h>
#include <bcrypt.h>
#pragma comment(lib, "bcrypt.lib")

namespace s_rng
{
	inline void secure_random_bytes(void* buf, std::size_t len)
	{
		if (BCryptGenRandom(nullptr,
			static_cast<PUCHAR>(buf),
			static_cast<ULONG>(len),
			BCRYPT_USE_SYSTEM_PREFERRED_RNG) != 0)
			throw std::system_error(GetLastError(),
				std::system_category(),
				"BCryptGenRandom failed");
	}
	class secure_rng
	{
	public:
		using result_type = std::uint64_t;

		static constexpr result_type _min() { return 0; }
		static constexpr result_type _max() { return ~result_type(0); }

		result_type operator()()
		{
			result_type r;
			secure_random_bytes(&r, sizeof(r));
			return r;
		}
	};
}

static const int MR_BASE[12] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };

enum KEY_LENGTH
{
	RSA32 = 32,    // for testing only
	RSA256 = 256,   // Unsafe
	RSA1024 = 1024, // Unsafe
	RSA2048 = 2048, // Recommended for short lived signatures
	RSA3072 = 3072, // Recommended for normal messages
	RSA4096 = 4096, // For high security to 2030+ (!!!)KEY GENERATION CAN TAKE A FEW SECONDS(!!!)
};

struct rsa_pub
{
	big_int n; // prime factor
	big_int e; // coprime to n
	barrett::barrett_ctx ctx;
};

struct rsa_priv
{
	big_int n; // prime factor
	big_int d; // inverse of e
	big_int p, q; // prime factors
	barrett::barrett_ctx ctx;
};

struct egcd_result { big_int g, x, y; };

class rsa
{
private:
	KEY_LENGTH _kl;
	s_rng::secure_rng rng;

	int small_primes[50000];

public:
	rsa(KEY_LENGTH kl)
		: _kl(kl)
	{
		fill_small_primes();
	}

private:
	void fill_small_primes()
	{
		constexpr std::size_t N_PRIMES = sizeof(small_primes) / sizeof(int);     // 50 000
		const std::size_t limit =
			static_cast<std::size_t>(N_PRIMES *
				(std::log(N_PRIMES) + std::log(std::log(N_PRIMES)))) + 50;

		std::vector<bool> is_composite(limit + 1, false);

		std::size_t count = 0;
		for (std::size_t p = 2; p <= limit && count < N_PRIMES; ++p) {
			if (!is_composite[p]) {
				small_primes[count++] = static_cast<int>(p);
				if (p * p <= limit)                 // avoid overflow
					for (std::size_t m = p * p; m <= limit; m += p)
						is_composite[m] = true;
			}
		}
	}

	big_int gcd(big_int a, big_int b)
	{
		while (b._bit_len) {
			big_int r = a % b;
			a = std::move(b);
			b = std::move(r);
		}
		return a;
	}

	egcd_result egcd(big_int a, big_int b)
	{
		big_int old_r = a, r = b;
		big_int old_s(1), s(0);
		big_int old_t(0), t(1);

		while (r._bit_len) {
			big_int q = old_r / r;
			big_int tmp;

			tmp = old_r - q * r;  old_r = std::move(r);  r = std::move(tmp);
			tmp = old_s - q * s;  old_s = std::move(s);  s = std::move(tmp);
			tmp = old_t - q * t;  old_t = std::move(t);  t = std::move(tmp);
		}
		return { old_r, old_s, old_t };
	}

	big_int choose_e(const big_int& phi)
	{
		static const big_int default_e(65537);
		if (gcd(default_e, phi) == big_int(1))
			return default_e;

		big_int e = default_e + big_int(2);
		while (gcd(e, phi) != big_int(1))
			e = e + big_int(2);
		return e;
	}

	big_int random_in_range(const big_int& min, const big_int& max)
	{
		if (big_int::compare(min, max) > 0)
			throw std::invalid_argument("random_in_range: min > max");

		if (min == max) return min;

		big_int span = (max - min) + big_int(1);
		std::size_t k = span._bit_len;

		big_int r;
		do {
			r = random_k_bits(k);
		} while (r >= span);

		return min + r;
	}

	bool divisible_by_small(const big_int& n)
	{
		for (std::uint32_t p : small_primes)
			if (n.mod_uint32(p) == 0)
				return true;
		return false;
	}

	bool probably_prime(const big_int& n)
	{
		if (n <= big_int(3))   return n == big_int(2) || n == big_int(3);
		if (!n.is_odd())       return false;
		if (divisible_by_small(n)) return false;

		big_int d = n - 1;
		std::size_t s = 0;
		while (!d.is_odd()) { d = d >> 1; ++s; }

		barrett::barrett_ctx ctx = barrett::make_ctx(n);

		for (int a_int : MR_BASE) {
			if (big_int(a_int) >= n) continue;

			big_int a(a_int);
			big_int x = barrett::f_powmod(a, d, ctx);

			if (x == 1 || x == n - 1) continue;

			bool passed = false;
			for (std::size_t r = 1; r < s; ++r) {
				x = barrett::f_powmod(x, big_int(2), ctx);
				if (x == n - 1) { passed = true; break; }
			}
			if (!passed) return false;
		}
		return true;
	}

	big_int random_k_bits(std::size_t k)
	{
		if (k == 0) return big_int();

		big_int r;
		r._v.resize((k + 31) / 32);
		for (auto& limb : r._v) limb = rng();

		std::size_t hi_bits = k & 31;

		if (hi_bits == 0) {
			r._v.back() |= 0x8000'0000u;
		}
		else {
			std::uint32_t mask = (1u << hi_bits) - 1u;
			r._v.back() &= mask;
			r._v.back() |= 1u << (hi_bits - 1);
		}

		r._v[0] |= 1u;
		r.trim();

		return r;
	}

public:
	std::pair<rsa_pub, rsa_priv> generate_key_pair()
	{
		std::size_t lp = (_kl + 1) / 2;
		std::size_t lq = _kl - lp;

		big_int p, q;
		do {
			p = random_k_bits(lp);
		} while (!probably_prime(p));
		do {
			q = random_k_bits(lp);
		} while (!probably_prime(q));

		big_int n = p * q;
		if (n._bit_len != _kl) return generate_key_pair();

		big_int phi = (p - 1) * (q - 1);
		big_int e = choose_e(phi);
		big_int d = modinv(e, phi);

		rsa_pub  pub{ n, e };
		rsa_priv prv{ n, d, p, q };
		pub.ctx = barrett::make_ctx(n);
		prv.ctx = barrett::make_ctx(n);

		return { pub, prv };
	}

	big_int rsa_encrypt(const big_int& m, const rsa_pub& pub)
	{
		if (m >= pub.n) throw std::invalid_argument("message too large");
		return barrett::f_powmod(m, pub.e, pub.ctx);
	}

	big_int rsa_decrypt(const big_int& c, const rsa_priv& prv)
	{
		if (c >= prv.n) throw std::invalid_argument("cipher > modulus");
		return barrett::f_powmod(c, prv.d, prv.ctx);
	}

	// returns a byte string not suitable for network traffic
	std::string rsa_encrypt(const std::string& str, const rsa_pub& pub)
	{
		std::vector<std::uint8_t> bytes(str.begin(), str.end());
		std::vector<std::uint8_t> encrypted;

		std::size_t modulus_bytes = (pub.n._bit_len + 7) / 8;
		std::size_t plain_bytes = (pub.n._bit_len - 1) / 8;

		for (std::size_t off = 0; off < bytes.size(); off += plain_bytes)
		{
			std::size_t len = min(plain_bytes, bytes.size() - off);

			std::vector<std::uint8_t> block(bytes.begin() + off,
				bytes.begin() + off + len);

			big_int m = bytes_to_bigint(block);
			big_int c = rsa_encrypt(m, pub);

			std::vector<std::uint8_t> out = bigint_to_bytes(c);
			out.resize(modulus_bytes, 0);
			encrypted.insert(encrypted.end(), out.begin(), out.end());
		}

		return std::string(encrypted.begin(), encrypted.end());
	}

	// can ONLY decrypt a byte string, to decrypt a base64 string call b64_decrypt
	std::string rsa_decrypt(const std::string& str, const rsa_priv& prv)
	{
		std::vector<std::uint8_t> bytes(str.begin(), str.end());
		std::vector<std::uint8_t> plain;

		std::size_t modulus_bytes = (prv.n._bit_len + 7) / 8;
		if (bytes.size() % modulus_bytes != 0)
			throw std::invalid_argument("ciphertext length not a multiple of block size");

		for (std::size_t off = 0; off < bytes.size(); off += modulus_bytes)
		{
			std::vector<std::uint8_t> block(bytes.begin() + off,
				bytes.begin() + off + modulus_bytes);

			big_int c = bytes_to_bigint(block);
			big_int m = rsa_decrypt(c, prv);

			std::vector<std::uint8_t> out = bigint_to_bytes(m);
			plain.insert(plain.end(), out.begin(), out.end());
		}

		return std::string(plain.begin(), plain.end());
	}

	// returns a encrypted base64 string suitable for network traffic
	std::string b64_encrypt(const std::string& str, const rsa_pub& pub)
	{
		std::string byte_string = rsa_encrypt(str, pub);
		return base64::to_base64(byte_string);
	}

	// can ONLY decrypt base64 strings, use rsa_decrypt for byte strings
	std::string b64_decrypt(const std::string& str, const rsa_priv& prv)
	{
		std::string byte_string = base64::from_base64(str);
		return rsa_decrypt(byte_string, prv);
	}

public:
	static big_int bytes_to_bigint(const std::vector<std::uint8_t>& bytes)
	{
		big_int x;
		x._v.resize((bytes.size() + 3) / 4, 0);
		for (std::size_t i = 0; i < bytes.size(); ++i)
			x._v[i >> 2] |= static_cast<std::uint32_t>(bytes[i]) << ((i & 3) * 8);
		x.trim();
		return x;
	}

	static std::vector<std::uint8_t> bigint_to_bytes(big_int x)
	{
		std::vector<std::uint8_t> out;
		out.reserve(x._v.size() * 4);
		for (std::size_t i = 0; i < x._v.size(); ++i) {
			std::uint32_t limb = x._v[i];
			out.push_back(limb & 0xFF);
			out.push_back((limb >> 8) & 0xFF);
			out.push_back((limb >> 16) & 0xFF);
			out.push_back((limb >> 24) & 0xFF);
		}
		while (!out.empty() && out.back() == 0) out.pop_back();  // trim
		return out;
	}

	static big_int modinv(big_int a, big_int m)
	{
		a = a % m;
		if (a == big_int(0) || m <= big_int(1))
			throw std::invalid_argument("modinv: invalid input");

		big_int t(0), new_t(1);
		big_int r = m, new_r = a;

		while (new_r._bit_len) {
			big_int q = r / new_r;

			big_int tmp = new_t * q;
			while (t < tmp)
				t = t + m;
			tmp = t - tmp;
			t = std::move(new_t);
			new_t = tmp;

			tmp = r - q * new_r;
			r = std::move(new_r);
			new_r = std::move(tmp);
		}

		if (r != big_int(1))
			throw std::invalid_argument("modinv: numbers are not coprime");

		return t;
	}
};

namespace key_export
{
	std::string der_len(std::size_t L)
	{
		if (L < 128) return std::string(1, static_cast<char>(L));
		std::string out;
		while (L) { out.insert(out.begin(), static_cast<char>(L & 0xFF)); L >>= 8; }
		out.insert(out.begin(), static_cast<char>(0x80 | out.size()));
		return out;
	}

	std::string asn1_integer(const big_int& x)
	{
		std::vector<std::uint8_t> v = rsa::bigint_to_bytes(x);
		if (v.empty()) v.push_back(0);
		if (v.back() & 0x80) v.push_back(0);
		std::string s(v.rbegin(), v.rend());
		return std::string(1, 0x02) + der_len(s.size()) + s;
	}

	// PKCS #1
	std::string export_public_pem(const rsa_pub& pub)
	{
		std::string seq = asn1_integer(pub.n) + asn1_integer(pub.e);
		std::string der = std::string(1, 0x30) + der_len(seq.size()) + seq;
		std::string b64 = base64::to_base64(der);

		/* 64-col line break */
		for (std::size_t i = 64; i < b64.size(); i += 65) b64.insert(i, "\n");

		return "-----BEGIN RSA PUBLIC KEY-----\n" + b64 +
			"\n-----END RSA PUBLIC KEY-----\n";
	}

	// PKCS #1
	std::string export_private_pem(const rsa_priv& prv)
	{
		big_int p = prv.p;
		big_int q = prv.q;
		big_int d = prv.d;
		big_int n = prv.n;
		big_int e = (prv.ctx.n == prv.n) ? prv.ctx.n : big_int(65537);

		big_int dP = d % (p - 1);
		big_int dQ = d % (q - 1);
		big_int qInv = rsa::modinv(q, p);

		std::string seq = asn1_integer(0) +
			asn1_integer(n) +
			asn1_integer(e) +
			asn1_integer(d) +
			asn1_integer(p) +
			asn1_integer(q) +
			asn1_integer(dP) +
			asn1_integer(dQ) +
			asn1_integer(qInv);

		std::string der = std::string(1, 0x30) + der_len(seq.size()) + seq;

		/* ---- 3. Base-64 with 64-column breaks --------------------------- */
		std::string b64 = base64::to_base64(der);
		for (std::size_t i = 64; i < b64.size(); i += 65) b64.insert(i, "\n");

		/* ---- 4. wrap in PEM armor --------------------------------------- */
		return "-----BEGIN RSA PRIVATE KEY-----\n" +
			b64 +
			"\n-----END RSA PRIVATE KEY-----\n";
	}

	void export_pub_to_file(const rsa_pub& pub, const std::string& filename)
	{
		std::ofstream out(filename);
		if (out.is_open())
		{
			out << export_public_pem(pub);
		}
	}

	void export_prv_to_file(const rsa_priv& prv, const std::string& filename)
	{
		std::ofstream out(filename);
		if (out.is_open())
		{
			out << export_private_pem(prv);
		}
	}
}

namespace key_import
{
	struct reader {
		const std::vector<std::uint8_t>& v; std::size_t pos{ 0 };
		std::uint8_t u8() { return v[pos++]; }
		std::vector<std::uint8_t> take(std::size_t n)
		{
			auto b = v.begin() + pos, e = b + n; pos += n; return { b,e };
		}
		std::size_t len() { return v.size() - pos; }
	};

	std::vector<std::uint8_t> base64_pem_body(const std::string& pem)
	{
		auto beg = pem.find('\n');
		auto end = pem.rfind("-----END");
		std::string b64 = pem.substr(beg + 1, end - beg - 1);
		b64.erase(std::remove_if(b64.begin(), b64.end(), isspace), b64.end());
		std::string bin = base64::from_base64(b64);
		return { bin.begin(), bin.end() };
	}

	std::size_t der_len(reader& R)
	{
		std::size_t len = R.u8();
		if ((len & 0x80) == 0) return len;
		std::size_t octets = len & 0x7F;  len = 0;
		while (octets--) len = (len << 8) | R.u8();
		return len;
	}

	big_int der_integer(reader& R)
	{
		if (R.u8() != 0x02) throw std::runtime_error("ASN.1: expected INTEGER");
		std::size_t L = der_len(R);
		auto bytes = R.take(L);
		while (bytes.size() && bytes[0] == 0) bytes.erase(bytes.begin());
		std::reverse(bytes.begin(), bytes.end());
		return rsa::bytes_to_bigint(bytes);
	}

	// PKCS #1
	rsa_pub import_public_pem(const std::string& pem)
	{
		auto der = base64_pem_body(pem);
		reader R{ der };

		if (R.u8() != 0x30) throw std::runtime_error("ASN.1: expected SEQUENCE");
		der_len(R);

		big_int n = der_integer(R);
		big_int e = der_integer(R);

		rsa_pub pub{ n, e };
		pub.ctx = barrett::make_ctx(pub.n);
		return pub;
	}

	// PKCS #1
	rsa_priv import_private_pem(const std::string& pem)
	{
		auto der = base64_pem_body(pem);
		reader R{ der };

		if (R.u8() != 0x30) throw std::runtime_error("ASN.1: expected SEQUENCE");
		der_len(R);

		big_int version = der_integer(R);
		big_int n = der_integer(R);
		big_int e = der_integer(R);
		big_int d = der_integer(R);
		big_int p = der_integer(R);
		big_int q = der_integer(R);

		rsa_priv prv{ n, d, p, q };
		prv.ctx = barrett::make_ctx(prv.n);
		return prv;
	}

	rsa_pub import_pub_from_file(const std::string& filename)
	{
		std::ifstream in(filename);
		if (in.is_open())
		{
			std::stringstream buffer;
			buffer << in.rdbuf();
			return import_public_pem(buffer.str());
		}
		std::cout << "No pem: " << filename << " found\n";
		return rsa_pub{};
	}

	rsa_priv import_prv_from_file(const std::string& filename)
	{
		std::ifstream in(filename);
		if (in.is_open())
		{
			std::stringstream buffer;
			buffer << in.rdbuf();
			return import_private_pem(buffer.str());
		}
		std::cout << "No pem: " << filename << " found\n";
		return rsa_priv{};
	}

	bool verify_keypair(const rsa_pub& pub, const rsa_priv& prv)
	{
		if (pub.n != prv.n) return false;

		big_int phi = (prv.p - 1) * (prv.q - 1);
		big_int ed = barrett::f_powmod(pub.e, big_int(1), prv.ctx);

		ed = (pub.e * prv.d) % phi;
		return ed == big_int(1);
	}
}