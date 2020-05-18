#include <iostream>
#include <sstream>

int main()
{
	std::stringstream ssl ("lm_0	1");
	std::string key;
	double val;

	ssl >> key;
	ssl >> val;

	std::cout << "key " << key << " and val " << val << std::endl;
	return 0;
}
