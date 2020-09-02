#include "express.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

DefineErrorMess( EXPR_DIVZERO, "Dividing by zero" );
DefineErrorMess( EXPR_INVVAR, "Invalid variable name" );
DefineErrorMess( EXPR_NOTASS, "Variable not assigned" );
DefineErrorMess( EXPR_EXPREQ, "Required }" );
DefineErrorMess( EXPR_RPREQ, "Required )" );
DefineErrorMess( EXPR_SEMREQ, "Required ;" );
DefineErrorMess( EXPR_UNTN, "Unterminated name" );
DefineErrorMess( EXPR_EXPN, "Name expected" );
DefineErrorMess( EXPR_FUNC, "Unknown function" );
DefineErrorMess( EXPR_SYNTAX, "Expression syntax error" );
DefineErrorMess( EXPR_CHK, "Checking false" );
DefineErrorMess( EXPR_ERR, "User defined" );
DefineErrorMess( EXPR_STR, "String is empty" );
DefineErrorMess( EXPR_INVS, "Invalid string" );


string* TStringHandle::empty = new string();

TStringHandle::TStringHandle( double v )
{
	handle.str.id = empty;
	handle.val = v;
}

TStringHandle::TStringHandle( string s )
{
	handle.str.id = empty;
	if( !s.empty() ) {
		string* st = new string( s );
		handle.str.ptr = st;
	} else
		handle.val = 0;
}

string TStringHandle::Str()
{
	if( handle.val ) {
		check( handle.str.id == empty, EXPR_INVS );
		return *handle.str.ptr;
	} else
		return *empty;
}

static double chk( double val )
{
	check( val, EXPR_CHK );
	return val;
}

static double err( double str )
{
	check( str, EXPR_STR );
	TStringHandle mess = str;
	Signal( EXPR_ERR, mess.Str() );
	return str;
}

static function funList[FUNCTIONS] = {
{ "ABS",   4, abs },
{ "ACOS",  0, acos },
{ "ASIN",  0, asin },
{ "ATAN",  0, atan },
{ "COS",   4, cos },
{ "CH",    0, cosh },
{ "CHK",   0, chk },
{ "CEIL",  0, ceil },
{ "EXP",   2, exp },
{ "ERR",   0, err },
{ "FLOOR", 1, floor },
{ "LN",    2, log },
{ "LG",    0, log10 },
{ "SIN",   3, sin },
{ "SH",    0, sinh },
{ "SQRT",  0, sqrt },
{ "TAN",   2, tan },
{ "TH",    0, tanh },
{ "",      1, 0 }
};

void print( string name, double value )
{
	printf( DUMP_FORMAT, name.c_str(), value );
}

TVarNode::TVarNode()
{
	memset( this, 0, sizeof( TVarNode ) );
}

TVarNode::~TVarNode()
{
	for( int i = 0 ; i < LETTERS ; i++ )
		if( tab[i] )
			delete tab[i];
}

void TVarNode::Assign( const char* pname, const char* pend, double value )
{
	assert( pname && pend );
	if( pname < pend ) {
		TVarNode*& v = tab[index( *pname )];
		if( !v )
			v = new TVarNode;
		v->Assign( pname + 1, pend, value );
	} else {
		ass = true;
		val = value;
	}
}

bool TVarNode::Lookup( const char* pname, const char* pend, double & value ) const
{
	assert( pname && pend );
	if( pname < pend ) {
		TVarNode* v = tab[index( *pname )];
		if( !v )
			return false;
		else
			return v->Lookup( pname + 1, pend, value );
	} else {
		value = val;
		return ass;
	}
}

int TVarNode::index( char c )
{
	if( IsDigit( c ) )
		return c - '0';
	if( c == '_' )
		return 10;
	c = ToUpper( c );
	check( c >= 'A' && c <= 'Z', EXPR_INVVAR );
	return c - 'A' + 11;
}

char TVarNode::letter( int index )
{
	assert( index >= 0 && index < LETTERS );
	if( index < 10 )
		return index + '0';
	else if( index == 10 )
		return '_';
	else
		return index + 'A' - 11;
}

void TVarNode::Dump( string prefix ) const
{
	if( ass )
		print( prefix, val );
	for( int i = 0 ; i < LETTERS ; i++ )
		if( tab[i] )
			tab[i]->Dump( prefix + letter( i ) );
}

TExpression::TExpression( const char* expr )
{
	expression = expr;
	curToken = END;
	tvalue = 0;
	ptr = pname = pend = 0;
	Set( "PI", M_PI );
}

double TExpression::Eval()
{
	double e = 0;
	if( !( ptr = expression ) )
		return e;
	do {
		getToken();
		e = expr();
		check(  curToken == SEM ||
				curToken == END, EXPR_SEMREQ );
	} while( curToken != END );
	return e;
}

void   TExpression::assign( const char* pname, const char* pend, double value )
{
	assert( pname && pend );
	check( pname < pend, EXPR_INVVAR );
	check( IsFirstLetter( *pname ), EXPR_INVVAR );
	vars.Assign( pname, pend, value );
}

double TExpression::lookup( const char* pname, const char* pend ) const
{
	assert( pname && pend );
	check( pname < pend, EXPR_INVVAR );
	check( IsFirstLetter( *pname ), EXPR_INVVAR );
	double val = 0;
	if( !vars.Lookup( pname, pend, val ) )
		if( *pend )
			Signal( EXPR_NOTASS );
		else
			Signal( EXPR_NOTASS, pname );
	return val;
}

void   TExpression::Set( const char *name, double value )
{
	assign( name, strchr( name, 0 ), value );
}

void   TExpression::Set( const char name, double value )
{
	char s[] = { name, 0 };
	assign( s, s + 1, value );
}

double TExpression::Get( const char *name ) const
{
	return lookup( name, strchr( name, 0 ) );
}

double TExpression::Get( const char name ) const
{
	char s[] = { name, 0 };
	return lookup( s, s + 1 );
}

char TExpression::nextChar()
{
	const char* p = ptr;
	while( *p )
		switch( *p ) {
			case ' ' :
			case '\t' :
			case '\r' :
			case '\n' :
				p++;
				break;
			default :
				return *p;
		}
	return 0;
}

TokenValue TExpression::getToken()
{
	while( *ptr ) {
		int count;
		char next;
		char c = *ptr;
		pname = pend = ptr;
		switch( c ) {
			case ' ' :
			case '\t' :
			case '\r' :
			case '\n' :
				ptr++;
				break;
			case '=' :
				if( *++ptr == '=' ) {
					ptr++;
					return curToken = EQ;
				} else
					return curToken = ASSIGN;
			case ':' :
				return curToken = ELSE;
			case '!' :
				if( *++ptr == '=' ) {
					ptr++;
					return curToken = NOTEQ;
				} else
					return curToken = NOT;
			case '>' :
				if( *++ptr == '=' ) {
					ptr++;
					return curToken = MOREQ;
				} else
					return curToken = MOR;
			case '<' :
				if( *++ptr == '=' ) {
					ptr++;
					return curToken = LESEQ;
				} else
					return curToken = LES;
			case '{' :
				pname = pend = ++ptr;
				for( count = 1 ; *ptr ; ptr++ ) {
					if( *ptr == '{' )
						count++;
					else if( *ptr == '}' )
						count--;
					if( !count )
						break;
				}
				check( *ptr, EXPR_EXPREQ );
				pend = ptr++;
				return curToken = EXPRESS;
			case '\'' :
				pname = pend = ++ptr;
				for( ; *ptr ; ptr++ )
					if( *ptr == '\'' )
						break;
				check( *ptr, EXPR_UNTN );
				pend = ptr++;
				next = nextChar();
				if( next == '(' ) {
					check( pend > pname, EXPR_EXPN );
					return curToken = EFUNCTION;
				} else if( next == '{' ) {
					check( pend > pname, EXPR_EXPN );
					return curToken = CONTEXT;
				} else
					return curToken = NAME;
			case '\\' :
				tvalue = *++ptr;
				check( tvalue, EXPR_SYNTAX );
				pend = ++ptr;
				return curToken = NUMBER;
			default :
				if( c == '0' && *( ptr + 1 ) == 'x' ) {
					ptr += 2;
					tvalue = strtol( ptr, ( char** )&pend, 16 );
					check( pend > ptr, EXPR_SYNTAX );
					ptr = pend;
					return curToken = NUMBER;
				} else if( IsDigit( c ) || c == '.' ) {
					tvalue = strtod( ptr, ( char** )&pend );
					check( pend > ptr, EXPR_SYNTAX );
					ptr = pend;
					return curToken = NUMBER;
				} else if( IsFirstLetter( c ) ) {
					for( ; *ptr ; ptr++ )
						if( !IsLetter( *ptr ) )
							break;
					pend = ptr;
					if( nextChar() == '(' )
						return curToken = FUNCTION;
					else
						return curToken = VAR;
				} else {
					check( c > ' ', EXPR_SYNTAX );
					ptr++;
					return curToken = TokenValue( c );
				}
		}
	}
	return curToken = END;
}

double TExpression::expr()
{
	if( curToken == SEM  ||
		curToken == END )
		return 0;
	else
		return log();
}

double TExpression::log()
{
	double left = comp();
	for(;;) 
		if( curToken == CONCAT ) {
			getToken();
			TStringHandle l = left;
			TStringHandle r = comp();
			TStringHandle conc = string( l.Str() + r.Str() );
			left = conc;
		} else 
	                switch( curToken ) {
		                case AND :
			                getToken();
				        left = ( comp() && left );
					break;
	                        case OR :
		                        getToken();
			                left = ( comp() || left );
				        break;
	                        default :
		                        return left;
			}
}

double TExpression::comp()
{
	double left = add();
	for(;;)
		switch( curToken ) {
			case EQ :
				getToken();
				left = ( left == add() );
				break;
			case NOTEQ :
				getToken();
				left = ( left != add() );
				break;
			case MOR :
				getToken();
				left = ( left > add() );
				break;
			case LES :
				getToken();
				left = ( left < add() );
				break;
			case MOREQ :
				getToken();
				left = ( left >= add() );
				break;
			case LESEQ :
				getToken();
				left = ( left <= add() );
				break;
			default :
				return left;
		}
}

double TExpression::add()
{
	double left = mul();
	for(;;)
		switch( curToken ) {
			case PLUS :
				getToken();
				left += mul();
				break;
			case MINUS :
				getToken();
				left -= mul();
				break;
			default :
				return left;
		}
}

double TExpression::mul()
{
	double d;    
	double left = pow();
	for(;;)
		switch( curToken ) {
			case MUL :
				getToken();
				left *= pow();
				break;
			case DIV :
				getToken();
				d = pow();
				check( d, EXPR_DIVZERO );
				left /= d;
				break;
			case MOD :
				getToken();
				left = fmod( left, pow() );
				break;
			default :
				return left;
		}
}

double TExpression::pow()
{
	double left = prim();
	for(;;)
		switch( curToken ) {
			case POW :
				getToken();
				left = ::pow( left, prim() );
				break;
			default :
				return left;
	}
}

double TExpression::prim()
{
	if( curToken == NAME ) {
		TStringHandle ret( tokenName() );
		getToken();
		return ret;
	}        
	double val;
	TokenValue next;
	const char* name = pname;
	const char* end  = pend;
	switch( curToken ) {
		default :
			Signal( EXPR_SYNTAX );
		case NUMBER :
			val = tvalue;
			getToken();
			return val;
		case VAR :
			next = getToken();
			if( next == ASSIGN ) {
				getToken();
				val = log();
				assign( name, end, val );
				return val;
			} else
				return lookup( name, end );
		case MINUS :
			getToken();
			return -prim();
		case NOT :
			getToken();
			return !prim();
		case LP :
			getToken();
			val = log();
			check( curToken == RP, EXPR_RPREQ );
			getToken();
			return val;
		case FUNCTION :
			FunctionPointer fun = getFunction( pname, pend );
			getToken();
			val = prim();
			check( fun, EXPR_FUNC );
			return fun( val );
	}
}

FunctionPointer TExpression::getFunction( const char* pname, const char* pend )
{
	assert( pname );
	for( int i = 0 ; funList[i].ptr ; i += funList[i].shift )
		if( funList[i].name[0] == ToUpper( pname[0] ) )
			for(;;) {
				const char *f = funList[i].name;
				const char *n = pname;
				for( f++, n++ ; *f ; f++, n++ )
					if( *f != ToUpper( *n ) )
						break;
				if( !*f && n == pend )
					return funList[i].ptr;
				if( funList[++i].shift )
					return 0;
			}
	return 0;
}

string TExpression::ErrorDump() const
{
	const char *e = ( ptr - expression > DUMP_LENGTH ) ?
		ptr - DUMP_LENGTH : expression;
	char* p = (char*)ptr;
	char c = *p;
	*p = 0;
	string dump = string( e );
	*p = c;
	return dump;
}

string TExpression::tokenName()
{
	char* p = (char*)pend;
	char c = *p;
	*p = 0;
	string name = pname;
	*p = c;
	int pos = 0;
	while( ( pos = name.find( '\\', pos ) ) != string::npos )
		if( name[pos+1] == 'n' ) {
			name.insert( pos++, 1, 13 );
			name.insert( pos++, 1, 10 );
		} else
			pos++;
	return name;
}
