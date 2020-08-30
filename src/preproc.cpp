#include "preproc.h"
#include "express.h"
#include <stdio.h>

DefineErrorMess( PP_SYNTAX, "Preprocessor syntax error" );
DefineErrorMess( PP_RPREQ, "Preprocessor missing )" );
DefineErrorMess( PP_ARGS, "Too many macro arguments" );
DefineErrorMess( PP_NOMACRO, "Macro not available" );
DefineErrorMess( PP_NOARG, "Macro argument expected" );
DefineErrorMess( PP_RECURSION, "Recursive macro" );
DefineErrorMess( PP_STRING, "Invalid string clipping" );

static char* next( char *p )
{
	assert( p );
	while( *p )
		switch( *(++p) ) {
			case ' ' :
			case '\t' :
			case '\r' :
			case '\n' :
				break;
			default :
				return p;
		}
	return p;
}

char* TPreprocessor::Run()
{
	char *ptr = 0;
	try {
		ptr = run();
	} _catch {
		if( lastError() == INTERNAL_ERROR )
			throw;
		else if( !*LastErrorParam() )
			Signal( lastError(), ErrorDump() );
		else
			throw;
	}
	return ptr;
}

string TPreprocessor::ErrorDump()
{
	if( right() > DUMP_LENGTH ) {
		char* p = Ptr() + DUMP_LENGTH;
		char  c = *p;
		*p = 0;
		string dump = string( Ptr() );
		*p = c;
		return dump;
	} else
		return string( Ptr() );
}

PPTokenValue TPreprocessor::lookToken()
{
	clear();
	char *p = Search( '#' );
	if( !p )
		return curToken = ENDM;
	p = next( p );
	if( *p == '#' )
		return curToken = TRUNC;
	if( *p == '(' )
		return curToken = MACRO;
	if( *p == '?' ) {
		p = next( p );
		check( IsDigit( *p ), PP_SYNTAX );
		tindex = *p - '0';
		return curToken = MACRO_COND;
	}
	if( IsDigit( *p ) ) {
		tindex = *p - '0';
		return curToken = MACRO_VAR;
	}
	if( IsAlpha( *p ) ) {
		tindex = ToUpper( *p ) - 'A';
		p = next( p );
		if( *p == '(' )
			return curToken = FOR;
		if( *p == '[' )
			return curToken = FOR_VAR;
		Signal( PP_SYNTAX );
	}
	Signal( PP_SYNTAX );
	return ENDM;
}

#define sscanf sscanf_s

void TPreprocessor::getToken()
{
	char *begin = Ptr();
	char *p = begin;
	int n = 0;
	assert( *p == '#' );
	p = next( p );
	switch( curToken ) {
		case FOR :
			p = next( p );
			p = argScan( p );
			check( argc == 2, PP_SYNTAX );
			check( argv[0] && argv[1], PP_SYNTAX );
			break;
		case MACRO :
			p = argScan( p );
			check( argv[0], PP_SYNTAX );
			break;
		case FOR_VAR :
			p = next( p );
			sscanf( p, "[ \t\r\n%n%d \t\r\n%n, \t\r\n%d \t\r\n%n",
						&n, &tleft, &n, &tright, &n );
			p += n;
			check( *p == ']', PP_SYNTAX );
			break;
		case MACRO_COND :
			p = next( p );
			p = next( p );
			check( *p == '(', PP_SYNTAX );
			p = argScan( p );
			check( argc <= 2, PP_SYNTAX );
			check( argv[0] || argv[1], PP_SYNTAX );
			break;
		case MACRO_VAR :
		case TRUNC :
			break;
		default :
			assert( false );
	}
	Delete( size_t( p - begin ) + 1 );
}

void TPreprocessor::skipToken()
{
	char *begin = Ptr();
	char *p = begin;
	assert( *p == '#' );
	p = next( p );
	if( *p == '#' )
		Skip( size_t( p - begin ) + 1 );
	else
		Skip( 1 );
}

void TPreprocessor::init()
{
	parent = 0;
	for( int i = 0 ; i < 10 ; i++ )
		argv[i] = 0;
	curToken = ENDM;
	tindex = 0;
	tleft  = 0;
	tright = -1;
	argc   = 0;
}

void TPreprocessor::clear()
{
	for( int i = 0 ; i < 10 ; i++ )
		if( argv[i] ) {
			delete argv[i];
			argv[i] = 0;
		}
	curToken = ENDM;
	tindex = 0;
	tleft  = 0;
	tright = -1;
	argc   = 0;
}

char* TPreprocessor::argScan( char *p )
{
	assert( !argc );
	assert( *p == '(' );
	for( ; *p != ')' ; argc++ ) {
		check( argc < 10, PP_ARGS );
		char *begin = p = next( p );
		char *end = 0;
		int c = 0;
		int d = 0;
		int e = 0;
		for( ; ( *p != ',' && *p != ')' ) || c || d || e ; p = next( p ) ) {
			check( *p, PP_RPREQ );
			end = p;
			if( *p == '\'' )
				e ^= 1;
			else if( !e )
				if( *p == '(' )
					c++;
				else if( *p == ')' )
					c--;
				else if( *p == '[' )
					d++;
				else if( *p == ']' )
					d--;
		}
		if( end )
			argv[argc] = new TCharBuffer( begin, size_t( end - begin ) + 1 );
	}
	return p;
}

TMacroProcessor::TMacroProcessor( TMacroProcessor* Parent )
{
	assert( Parent );
	parent = Parent;
	TCharBuffer* name = Parent->argv[0];
	assert( name );
	InsFile( name->Text() );
}

char* TMacroProcessor::run()
{
	static int recursion = 0;
	check( recursion++ < MAX_RECURSION, PP_RECURSION );
	for( Begin() ; lookToken() != ENDM ; )
		if( curToken == MACRO_VAR ) {
			check( parent, PP_NOMACRO );
			getToken();
			TCharBuffer* arg = parent->argv[tindex];
			check( arg, PP_NOARG );
			Insert( arg->Text() );
		} else if( curToken == MACRO_COND ) {
			check( parent, PP_NOMACRO );
			getToken();
			TCharBuffer* arg = parent->argv[tindex];
			if( arg && argv[0] )
				Insert( argv[0]->Text() );
			if( !arg && argv[1] )
				Insert( argv[1]->Text() );
		} else if( curToken == TRUNC ) {
			check( parent, PP_NOMACRO );
			Trunc();
		} else
			skipToken();
	for( Begin() ; lookToken() != ENDM ; )
		if( curToken == MACRO ) {
			getToken();
			TMacroProcessor macro( this );
			Insert( macro.Run() );
		} else
			skipToken();
	recursion--;
	return this->Text();
}

