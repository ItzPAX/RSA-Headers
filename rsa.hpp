#pragma once
#include <random>
#include "bigint.hpp"

enum KEY_LENGTH
{
	RSA256 = 256,
	RSA1024 = 1024,
	RSA2048 = 2048, // Recommended for short lived signatures
	RSA4096 = 4096, // For security to 2030+
};

struct rsa_pub
{
	big_int n; // prime factor
	big_int e; // coprime to n
};

struct rsa_priv
{
	big_int n; // prime factor
	big_int d; // inverse of e
};

class rsa
{
private:
	KEY_LENGTH _kl;
	rsa_pub pub;
	rsa_priv priv;

	std::uniform_int_distribution<std::uint64_t> u64;

public:
	rsa(KEY_LENGTH kl)
		: _kl(kl)
	{
		u64 = std::uniform_int_distribution<std::uint64_t>(0, ~0ULL);
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

	bool probably_prime(const big_int& n, std::mt19937_64& gen, int rounds = 64)
	{
		if (n <= big_int(3)) return n == big_int(2) || n == big_int(3);
		if (!n.is_odd())    return false;

		// write n-1 = d·2^s
		big_int d = n - big_int(1);
		std::size_t s = 0;
		while (!d.is_odd()) { d = d / big_int(2); ++s; }

		for (int i = 0; i < rounds; ++i) {
			big_int a = random_in_range(big_int(2), n - big_int(2), gen);
			big_int x = a.powmod(d, n);
			if (x == big_int(1) || x == n - big_int(1)) continue;

			bool continue_outer = false;
			for (std::size_t r = 1; r < s; ++r) {
				x = x.powmod(big_int(2), n);
				if (x == n - big_int(1)) { continue_outer = true; break; }
			}
			if (continue_outer) continue;
			return false;
		}
		return true;
	}

	big_int random_k_bits(std::size_t k, std::mt19937_64& gen)
	{
		big_int r;
		r._value.resize((k + 7) / 8);

		for (auto& b : r._value)
			b = static_cast<std::uint8_t>(u64(gen) & 0xFF);

		r._value.back() |= 1u << ((k - 1) & 7);
		r._value[0] |= 1u;
		r.trim();
		return r;
	}

	std::pair<rsa_pub, rsa_priv> generate_key_pair()
	{
		std::size_t lp = (_kl + 1) / 2;   // ceil
		std::size_t lq = _kl - lp;

		std::random_device rd; std::mt19937_64 gen(rd());

		big_int p, q;
		do { p = random_k_bits(lp, gen); } while (!probably_prime(p, gen, 16));
		do { q = random_k_bits(lq, gen); } while (q == p || !probably_prime(q, gen, 16));

		big_int n = p * q;
		if (n._bit_len != _kl) return generate_key_pair();   // very rare retry

		big_int phi = (p - 1) * (q - 1);
		big_int e("65537");
		//big_int d = modinv(e, phi);

		std::cout << "p: " << p.as_string() << "\nq: " << q.as_string() << std::endl;

		return std::make_pair<rsa_pub, rsa_priv>({}, {});
	}
};