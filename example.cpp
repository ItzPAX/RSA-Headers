#include <iostream>
#include "rsa.hpp"

int main()
{
	big_int_tests::run_tests();

	rsa rsa(RSA4096);
	auto [pub_key, prv_key] = rsa.generate_key_pair();
	
	std::cout << key_export::export_private_pem(prv_key) << std::endl;
	std::cout << key_export::export_public_pem(pub_key) << std::endl;
	std::cout << "VERIFY: " << (key_import::verify_keypair(pub_key, prv_key) ? "OK" : "FAIL") << "\n";
	
	std::string encrypted = rsa.b64_encrypt("Hello from RSA!", pub_key);
	std::cout << "ENC: " << encrypted << std::endl;
	std::string decrypted = rsa.b64_decrypt(encrypted, prv_key);
	std::cout << "DEC: " << decrypted << std::endl;
	
	std::string sign = rsa.b64_sign("Test Sign me!", prv_key);
	std::cout << "SIGNATURE: " << sign << std::endl;
	bool test = rsa.b64_verify("Test Sign me!", sign, pub_key);
	std::cout << (test ? "SIGN OK\n" : "SIGN FAIL\n");

	system("pause");

	return 0;
}