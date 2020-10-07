/* insensitive case comparison
 */


#ifndef ic_H
#define ic_H

#include <cstring>

// case insensitive comparison, available to anyone
struct ic { 
	bool operator()(const std::string& lhs, const std::string& rhs) const {
		return strcasecmp(lhs.c_str(), rhs.c_str()) < 0;
	}
};

#endif
