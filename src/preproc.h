#pragma once

#include "chbuff.h"

#define DUMP_LENGTH   24
#define MAX_RECURSION 16

ErrorCode( PP_SYNTAX );
ErrorCode( PP_RPREQ );
ErrorCode( PP_ARGS );
ErrorCode( PP_NOMACRO );
ErrorCode( PP_NOARG );
ErrorCode( PP_RECURSION );
ErrorCode( PP_STRING );

enum PPTokenValue {
FOR, MACRO, FOR_VAR, MACRO_VAR, MACRO_COND, TRUNC, ENDM
};

//---------------------------------------------------------------------------
//   FOR         #Alpha( mask, expr )
//   MACRO       #( name, arg1, arg2, .. )
//   FOR_VAR     #Alpha[]
//               #Alpha[ left ]
//               #Alpha[ left, right ]
//   MACRO_VAR   #Num
//   MACRO_COND  #?Num( true )
//               #?Num( , false )
//               #?Num( true, false )
//   TRUNC       ##
//---------------------------------------------------------------------------

class  TPreprocessor : protected TCharBuffer {
friend class TMacroProcessor;
public:
        TPreprocessor() { init(); }
        TPreprocessor( const char *str ) : TCharBuffer( str ) { init(); }

        char* Run();

        TCharBuffer::Save;
        TCharBuffer::Text;
        string ErrorDump();

protected:
        virtual char* run() { return this->Text(); }

        PPTokenValue lookToken();
        void getToken();
        void skipToken();

        PPTokenValue curToken;
        int tindex;
        int tleft;
        int tright;
        int argc;
        TCharBuffer* argv[10];

        TPreprocessor* parent;

private:
        void  init();
        void  clear();
        char* argScan( char *p );
};

class  TMacroProcessor : public TPreprocessor {
public:
        TMacroProcessor( const char *str ) : TPreprocessor( str ) {}
        TMacroProcessor( TMacroProcessor* Parent );
protected:
        char* run();
};

