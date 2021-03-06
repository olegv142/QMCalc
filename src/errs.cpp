#include "errs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <eh.h>

#include "compat.h"

DefineErrorMess( INTERNAL_ERROR, "Internal error" )
DefineErrorMess( NOGMEM, "Not enougth global memory" )
DefineErrorMess( NOLMEM, "Not enougth local memory" )
DefineErrorMess( NOSMEM, "Not enougth stack memory" )
DefineErrorMess( NOMEM,  "Not enougth memory" )

#define ASSERT_BUFF_LEN 128
#define TEXT_BUFF_LEN   128

ErrorDescriptor ErrorDescriptor::lastError;

static void _terminate()
{
	Failure();
} 

static terminate_function old_handler = set_terminate( _terminate );

const char* LastErrorTitle() 
{ 
	if( const char* lastErrorTitle = ErrorDescriptor::LastError().Title() )
		return lastErrorTitle; 
	else
		return "Unregistered Error";
}

const char* LastErrorText()  
{ 
	const ErrorDescriptor& lastError = ErrorDescriptor::LastError();
	static char textBuffer[TEXT_BUFF_LEN+1];
	textBuffer[TEXT_BUFF_LEN] = 0;
	if( const char* lastErrorType = lastError.Type() ) {
		const char* lastErrorParam = lastError.Param();
		if( lastErrorParam && *lastErrorParam )
			snprintf( textBuffer, TEXT_BUFF_LEN, "%s : %s", lastErrorType, lastErrorParam );
		else
			snprintf( textBuffer, TEXT_BUFF_LEN, "%s", lastErrorType );
	}
	return textBuffer;
}

const char* LastErrorParam() 
{ 
	const ErrorDescriptor& lastError = ErrorDescriptor::LastError();
	if( const char* lastErrorParam = lastError.Param() ) 
		return lastErrorParam; 
	return "";
} 

const char* BuildAssertMessage( const char* file, int line )
{
	static char assertBuffer[ASSERT_BUFF_LEN+1];
	assertBuffer[ASSERT_BUFF_LEN] = 0;
	snprintf( assertBuffer, ASSERT_BUFF_LEN, "file %s, line %d", file, line );
	return assertBuffer;
}

void Failure()
{
	Warning();
	exit(1);
}

void Warning()
{
	puts( "\n" );
	puts( LastErrorTitle() );
	puts( LastErrorText() );
}


