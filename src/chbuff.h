#pragma once

#include "errs.h"

typedef unsigned size_t;

ErrorCode( CHB_OUTMEM )

class  TCharBuffer {
public:
        TCharBuffer();
        TCharBuffer( const char *str );
        TCharBuffer( const char *str, size_t length );
        ~TCharBuffer() { delete pBuff; }

        char* Text() { return pBuff; }
        char* Ptr() { return ptr; }
        void  Begin() { ptr = pBuff; }
        char* Search( char c );
        void  Skip( size_t size );

        void Delete( size_t length );
        void Insert( const char *str );
        void InsFile( const char *name );
        void Trunc() { Delete( right() ); }

        void Save( const char *name );

protected:
        size_t length() { return len; }
        size_t left() { return size_t( ptr - pBuff ); }
        size_t right() { return len - left(); }

private:
        void   del( size_t length );
        void   ins( size_t length );

        char   *pBuff;
        char   *ptr;
        size_t len;
};
