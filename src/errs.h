#pragma once

#include <string>

using namespace std;

typedef const char* errorType;

class ErrorDescriptor {
public:
	ErrorDescriptor( const ErrorDescriptor& err )
		: type( err.type ), param( err.param )
	{
		lastError = this;
	}
	ErrorDescriptor( const errorType Type )
		: type(Type)
	{
		set_terminate();
		lastError = this;
	}
	ErrorDescriptor( const errorType Type, const char* Param )
		: type(Type), param(Param)
	{
		set_terminate();
		lastError = this;
	}
	ErrorDescriptor( const errorType Type, const string Param )
		: type(Type), param(Param)
	{
		set_terminate();
		lastError = this;
	}
	virtual ~ErrorDescriptor() {
		if( lastError == this ) lastError = 0;
	}

	virtual const char* Title() const { return "Application Error"; }
	virtual errorType   Type()  const { return type; }
	virtual const char* Param() const { return param.c_str(); }
	static const ErrorDescriptor* LastError();
	
private:
	ErrorDescriptor();
	ErrorDescriptor& operator=(ErrorDescriptor&);

	static void set_terminate();
	static ErrorDescriptor* lastError;

	errorType type;
	string param;
};

const char* LastErrorTitle();
const char* LastErrorText();
const char* LastErrorParam(); 

void Failure();
void Warning();

const char* BuildAssertMessage( const char* file, int line );

#define ErrorCode(code) extern const char* code;
#define DefineErrorMess(code,mess) const char* code = mess;

#define IF_NOT( cond ) ((cond)) ? ((void)0) : (void)
#define check( cond, error )  IF_NOT((cond)) ::Signal((error))
#define check3( cond, error, msg ) IF_NOT((cond)) ::Signal((error),(msg))
#define assert( cond )  IF_NOT((cond)) ::Signal( INTERNAL_ERROR, \
	BuildAssertMessage( __FILE__, __LINE__ ) )
#define checkS( cond, msg ) check3( (cond), INTERNAL_ERROR, (msg) )
#define _catch catch( ErrorDescriptor& )

inline void Signal( const char* error ) {
	ErrorDescriptor ed( error );
	throw( ed );
}

inline void Signal( const char* error, const char* msgParam ) {
	throw( ErrorDescriptor( error, msgParam ) );
}

inline void Signal( const char* error, const string msgParam ) {
	throw( ErrorDescriptor( error, msgParam ) );
}

inline errorType lastError() {
	return ErrorDescriptor::LastError()->Type();
}

ErrorCode(INTERNAL_ERROR)
ErrorCode(NOGMEM)
ErrorCode(NOLMEM)
ErrorCode(NOSMEM)
ErrorCode(NOMEM)

inline const ErrorDescriptor* ErrorDescriptor::LastError() { 
	assert( lastError );
	return lastError; 
}

