# RSA-Headers
A (triple) header only way to implement RSA-OAEP and RSA-PSS functionality

## Disclaimer
Do **NOT** use this in production environments

While the key size may be large, this implementation lacks multiple critical protections found in hardened cryptographic libraries.

## Implemented and TODO
- [x] big_int class for supporting keys up to 4096 bits
- [x] OAEP RSA encryption
- [x] PSS RSA signing
- [x] Key export and import
- [x] Key pair verification
- [ ] RSA blinding on all private key operations
- [ ] fast randomized Miller-Rabin rounds
- [ ] Montgomery multiplication

## Requirements
C++20 and later is required for the bit length calculations

Currently only supports windows systems due to secure random number generation

## Installation
Just copy the three header files into your project and include "rsa.hpp" for examples see Usage or examples.cpp

## Usage
### Creating a key pair
The constructor parameter is the size of the key to be used, the current maximum is 4096 bits
```cpp
#include "rsa.hpp"
int main()
{
	rsa rsa(RSA4096);
	auto [pub_key, prv_key] = rsa.generate_key_pair();
	return 0;
}
```
### Encrypting / Decrypting a string
Allows for encryption of strings of any length through splitting into chunks
```cpp
#include "rsa.hpp"
int main()
{
	rsa rsa(RSA4096);
	auto [pub_key, prv_key] = rsa.generate_key_pair();
	std::string encrypted = rsa.b64_encrypt("Hello from RSA!", pub_key);
	std::string decrypted = rsa.b64_decrypt(encrypted, prv_key);
	return 0;
}
```
### Signing and verifying a message
Allows for signing and verification of messages of any length
```cpp
#include "rsa.hpp"
int main()
{
	rsa rsa(RSA4096);
	std::string signature = rsa.b64_sign("Sign me!", prv_key);
	bool is_signed = rsa.b64_verify("Sign me!", signature, pub_key);
	return 0;
}
```
