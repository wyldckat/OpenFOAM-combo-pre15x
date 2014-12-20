#ifndef IOEXCEPTION_H
#define IOEXCEPTION_H

#include <string>

class IOException : public std::exception
{
 public:
    IOException( const char *str ) : std::exception(), description(str) {};
    virtual ~IOException() throw() {};
    virtual const char* what() const throw()
	{ return(this->description.c_str()); };

 private:
    std::basic_string<char> description;
};

#endif // IOEXCEPTION_H
