#pragma once
#include <random>
#include "bigint.hpp"

static const int MR_BASE[12] = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37 };

enum KEY_LENGTH
{
	RSA256 = 256,   // Unsafe
	RSA1024 = 1024, // Unsafe
	RSA2048 = 2048, // Recommended for short lived signatures
	RSA4096 = 4096, // For security to 2030+ (!!!)KEY GENERATION CAN TAKE A FEW SECONDS(!!!)
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
	barrett::barrett_ctx ctx;
};

struct egcd_result { big_int g, x, y; };

class rsa
{
private:
	KEY_LENGTH _kl;
	rsa_pub pub;
	rsa_priv priv;

	std::uniform_int_distribution<std::uint64_t> u64;

	int small_primes[512];

private:
	void fill_small_primes()
	{
		for (int i = 2, primes_found = 0; primes_found < sizeof(small_primes) / sizeof(int); i++)
		{
			bool prime = true;
			for (int j = 2; j < i; j++)
			{
				if (i % j == 0) { prime = false; break; }
			}

			if (prime)
			{
				small_primes[primes_found++] = i;
			}
		}
	}

public:
	rsa(KEY_LENGTH kl)
		: _kl(kl)
	{
		fill_small_primes();
		u64 = std::uniform_int_distribution<std::uint64_t>(0, ~0ULL);
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

	big_int modinv(big_int a, big_int m)
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

		return t;                            // already 0 … m-1
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

	big_int random_in_range(const big_int& min, const big_int& max, std::mt19937_64& gen)
	{
		if (big_int::compare(min, max) > 0)
			throw std::invalid_argument("random_in_range: min > max");

		if (min == max) return min;

		big_int span = (max - min) + big_int(1);
		std::size_t k = span._bit_len;

		big_int r;
		do {
			r = random_k_bits(k, gen);
		} while (r >= span);

		return min + r;
	}

	bool divisible_by_small(const big_int& n)
	{
		for (int p : small_primes)
			if ((n % big_int(p)).is_zero()) return true;
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

	big_int random_k_bits(std::size_t k, std::mt19937_64& gen)
	{
		if (k == 0) return big_int();

		static std::uniform_int_distribution<std::uint32_t> u32;

		big_int r;
		r._v.resize((k + 31) / 32);
		for (auto& limb : r._v) limb = u32(gen);

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

	std::pair<rsa_pub, rsa_priv> generate_key_pair()
	{
		std::cout << "generating key pair, this might take some time...\n";

		std::size_t lp = (_kl + 1) / 2;
		std::size_t lq = _kl - lp;

		std::random_device rd; std::mt19937_64 gen(rd());

		big_int p, q;
		do {
			p = random_k_bits(lp, gen);
		} while (!probably_prime(p));
		do {
			q = random_k_bits(lp, gen);
		} while (!probably_prime(q));

		big_int n = p * q;
		if (n._bit_len != _kl) return generate_key_pair();

		big_int phi = (p - 1) * (q - 1);
		big_int e = choose_e(phi);
		big_int d = modinv(e, phi);

		rsa_pub  pub{ n, e };
		rsa_priv prv{ n, d };

		std::cout << "N: " << n.as_string() << "\nd: " << d.as_string() << "\ne: " << e.as_string() << "\n";

		return { pub, prv };
	}

	big_int bytes_to_bigint(const std::vector<std::uint8_t>& bytes)
	{
		big_int x;
		x._v.resize((bytes.size() + 3) / 4, 0);
		for (std::size_t i = 0; i < bytes.size(); ++i)
			x._v[i >> 2] |= static_cast<std::uint32_t>(bytes[i]) << ((i & 3) * 8);
		x.trim();
		return x;
	}

	std::vector<std::uint8_t> bigint_to_bytes(big_int x)
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
};