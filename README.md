# faursal
Finally a usable RSA library - a (triple) header only way to implement basic RSA features in your personal projects

## Disclaimer
Do **NOT** use this in production environments

## Implemented and TODO
- [x] big_int class for supporting keys up to 4096 bits
- [x] OAEP RSA encryption
- [x] PSS RSA signing
- [x] Key export and import
- [x] Key pair verification
- [ ] RSA blinding on all private key operations
- [ ] fast randomized Miller-Rabin rounds
- [ ] Montgomery multiplication

## Usage
### Creating a key pair
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
