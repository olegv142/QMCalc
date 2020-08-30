#pragma once

#include "errs.h"

#define LETTERS  12+'Z'-'A'
#define FUNCTIONS       64
#define FUN_NAME_SIZE    9
#define DUMP_LENGTH     24
#define DUMP_FORMAT "%s = %.8g\n"

#define M_PI 3.14159265358979323846


ErrorCode( EXPR_DIVZERO );
ErrorCode( EXPR_INVVAR );
ErrorCode( EXPR_NOTASS );
ErrorCode( EXPR_EXPREQ );
ErrorCode( EXPR_RPREQ );
ErrorCode( EXPR_SEMREQ );
ErrorCode( EXPR_UNTN );
ErrorCode( EXPR_EXPN );
ErrorCode( EXPR_FUNC );
ErrorCode( EXPR_SYNTAX );
ErrorCode( EXPR_RET );
ErrorCode( EXPR_BRK );
ErrorCode( EXPR_CHK );
ErrorCode( EXPR_ERR );
ErrorCode( EXPR_STR );

union StringHandle {
        struct {
                string* ptr;
                string* id;
        } str;
        double   val;
};

class  TStringHandle {
public:
        TStringHandle( double v ) { handle.val = v; }
        TStringHandle( string s );
        operator double() { return handle.val; }
        string Str();
private:
        StringHandle handle;
		static string* empty;
};

typedef double ( *FunctionPointer )( double );

enum TokenValue {
        END, NUMBER, VAR, FUNCTION, EFUNCTION, EXPRESS, NAME, CONTEXT,
        EQ, NOT, MOR, LES, NOTEQ, MOREQ, LESEQ, STAT, AND = '&', OR = '|',
        PLUS = '+', MINUS = '-', MUL = '*', DIV = '/', MOD = '%', POW = '^',
        ASSIGN = '=', LP = '(', RP = ')', COL = ',', LOOP = '$',
        SEM = ';', IF = '?', ELSE = ':', WHILE = '@', CONCAT = '"'
};

struct function {
        char name[FUN_NAME_SIZE];
        int shift;
        FunctionPointer ptr;
};

struct TagVoid {
        char name[FUN_NAME_SIZE];
        int shift;
        void* ptr;
};

class  TVarNode {
public:
        TVarNode();
        ~TVarNode();
        void Stat( const char* pname, const char* pend, double value );
        void Assign( const char* pname, const char* pend, double value );
        bool Look( const char* pname, const char* pend, double & value );
        void Dump( string prefix );

private:
        int  index( char c );
        char letter( int index );

        double val;
        bool   ass;
        TVarNode*  tab[LETTERS];
};

class  TExpression {
friend double ret( double );
friend double brk( double );
public:
        TExpression( const char* expr = 0 );
	virtual ~TExpression() {}

        void SetExpression( const char* expr ) { expression = expr; }

        virtual double Calc();
        virtual void   Stat( const char* pname, const char* pend, double value );
        virtual void   Assign( const char* pname, const char* pend, double value );
        virtual double Look( const char* pname, const char* pend );

        void   stat( const char *name, double value );
        void   stat( const char name, double value );
        void   assign( const char *name, double value );
        void   assign( const char name, double value );
        double look( const char *name );
        double look( const char name );

        virtual void Dump( string prefix = "" )
                { var.Dump( prefix ); }

        string ErrorDump();

        static TExpression* CurExpression();

        static void RegisterFunction( char* name, FunctionPointer ptr );

protected:
        virtual FunctionPointer getFunction( const char* pname, const char* pend );
        virtual double expr();
        virtual double prim();

        double     log();
        TokenValue getToken();

        string tokenName();

        static void registerTag( TagVoid* array, char* name, void* ptr );

        const char*  expression;

        TokenValue   curToken;
        double       tvalue;
        const char*  ptr;
        const char*  pname;
        const char*  pend;

        TVarNode var;

        static double retVal;
        static TExpression* curExpression;

private:
        char   nextChar();

        double comp();
        double add();
        double mul();
        double pow();
};

inline bool IsDigit( unsigned char c )
{ return '0' <= c && c <= '9'; }

inline bool IsAlpha( unsigned char c )
{ return ( 'A' <= c && c <= 'Z') || ( 'a' <= c && c <= 'z'); }

inline bool IsFirstLetter( unsigned char c )
{ return c == '_' || IsAlpha( c ); }

inline bool IsLetter( unsigned char c )
{ return c == '_' || IsAlpha( c ) || IsDigit( c ); }

inline char ToUpper( unsigned char c )
{ if( 'a' <= c && c <= 'z') return c - 'a' + 'A'; return c; }


//--------------------------------------------------------------------------
//  Function   Operators             Prioritet
//
//  prim()     ! - = ()              high
//  pow()      ^
//  mul()      * / %
//  add()      + -
//  comp()     > >= < <= == !=
//  log()      & | "                 low
//--------------------------------------------------------------------------

//------ Gramma ------------------------------------------------------------
//  expression_list :
//      expression
//      expression_list ; expression
//
//  expression :
//      empty_string
//      log
//
//  log :
//      comp
//      log | comp
//      log & comp
//      log " comp
//
//  comp :
//      add
//      comp >  add
//      comp >= add
//      comp <  add
//      comp <= add
//      comp == add
//      comp != add
//
//  add :
//      mul
//      add + mul
//      add - mul
//
//  mul :
//      pow
//      mul * pow
//      mul / pow
//      mul % pow
//  pow :
//      prim
//      pow ^ prim
//
//  prim :
//      NUMBER
//      'string'
//      var
//      var = log
//      var := log
//      - prim
//      ( log )
//      function( log )
//
//  function :
//      sin         cos         tan
//      asin        acos        atan
//      sh          ch          th
//      exp         ln          lg
//      abs         floor       ceil
//      sqrt
//      err         - generate error
//      chk         - check condition
//      dump        - variables printout
//--------------------------------------------------------------------------
