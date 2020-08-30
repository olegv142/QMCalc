#include "cexpress.h"

void   TCompoundExpression::Assign( const char* pname, const char* pend, double value )
{
        assert( pname );
        if( parent && *pname == '_' )
                parent->Assign( ++pname, pend, value );
        else
                TExpression::Assign( pname, pend, value );
}

double TCompoundExpression::Look( const char* pname, const char* pend )
{
        assert( pname );
        if( parent && *pname == '_' )
                return parent->Look( ++pname, pend );
        assert( pname && pend );
        check( pname < pend, EXPR_INVVAR );
        check( IsAlpha( *pname ), EXPR_INVVAR );
        double val = 0;
        if( var.Look( pname, pend, val ) )
                return val;
        else if( !parent )
                if( *pend )
                        Signal( EXPR_NOTASS );
                else
                        Signal( EXPR_NOTASS, pname );
        return parent->Look( pname, pend );
}

double TCompoundExpression::Calc()
{
        char   c = 0;
        double e = 0;
        if( endPtr ) {
                c = *endPtr;
                *(char*)endPtr = 0;
        }
        try {
                e = TExpression::Calc();
        } _catch {
                errorType err = lastError();
                if( endPtr )
                        *(char*)endPtr = c;
                if( err == INTERNAL_ERROR )
                        throw;
                if( err != EXPR_RET )
                        if( !*LastErrorParam() && err != EXPR_BRK )
                                Signal( err, ErrorDump() );
                        else
                                throw;
                else
                        e = retVal;
        }
        if( endPtr )
                *(char*)endPtr = c;
        curExpression = parent;
        return e;
}

double TCompoundExpression::expr()
{
        if( curToken == SEM  ||
                curToken == ELSE ||
                curToken == END )
                return 0;
        TokenValue  t = curToken;
        double      v = tvalue;
        const char *p = ptr;
        const char *n = pname;
        const char *e = pend;
        double val = log();
        if( curToken == IF )
                if( val ) {
                        getToken();
                        val = expr();
                        if( curToken == ELSE )
                                skipExpr();
                } else {
                        if( searchElse() == ELSE ) {
                                getToken();
                                val = expr();
                        }
                }
        else if( curToken == WHILE ) {
                try {
                        while( val ) {
                                getToken();
                                val = expr();
                                check(  curToken == SEM ||
                                                curToken == END, EXPR_SEMREQ );
                                curToken = t;
                                tvalue = v;
                                ptr = p;
                                pname = n;
                                pend = e;
                                val = log();
                        }
                } _catch {
                        if( lastError() != EXPR_BRK )
                                throw;
                        else
                                val = retVal;
                }
                skipExpr();
        }
        return val;
}

double TCompoundExpression::prim()
{
        double val;
        TCompoundExpression* expr;
        switch( curToken ) {
                case EFUNCTION :
                        Signal( EXPR_SYNTAX );
                case EXPRESS :
                        expr = new TCompoundExpression( pname, pend, this );
                        try {
                                val = expr->Calc();
                        } _catch {
                                delete expr;
                                throw;
                        }
                        delete expr;
                        getToken();
                        return val;
                default :
                        return TExpression::prim();
        }
}

void TCompoundExpression::Dump( string prefix )
{
        if( parent )
                parent->Dump( "_" + prefix );
        TExpression::Dump( prefix );
}

TokenValue TCompoundExpression::skipExpr()
{
        if( curToken == SEM || curToken == END )
                return curToken;
        int c = 0;
        int e = 0;
        for( ; *ptr ; ptr++ ) {
                char ch = *ptr;
                if( ch == '\'' )
                        e ^= 1;
                else if( !e ) {
                        if( ch == '{' )
                                c++;
                        else if( ch == '}' )
                                c--;
                        else if( ch == ';' && !c )
                                break;
                        if( c < 0 )
                                break;
                }
        }
        check( c <= 0 && !e, EXPR_SYNTAX );
        return getToken();
}

TokenValue TCompoundExpression::searchElse()
{
        int c = 0;
        int d = 0;
        int e = 0;
        for( ; *ptr ; ptr++ ) {
                char ch = *ptr;
                if( ch == '\'' )
                        e ^= 1;
                else if( !e ) {
                        if( ch == '{' )
                                c++;
                        else if( ch == '}' )
                                c--;
                        else if( !c ) {
                                if( ch == ';' )
                                        break;
                                else if( ch == '?' )
                                        d++;
                                else if( ch == ':' )
                                        d--;
                                if( d < 0 )
                                        break;
                        }
                        if( c < 0 )
                                break;
                }
        }
        check( c <= 0 && !e, EXPR_SYNTAX );
        return getToken();
}
