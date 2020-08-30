#include "file.h"
#include "chbuff.h"
#include <string.h>
#include <malloc.h>

extern char **_argv;

DefineErrorMess( CHB_OUTMEM, "Char buffer out of memory" );

TCharBuffer::TCharBuffer()
{
	len = 0;
	ptr = pBuff = (char*)malloc( 1 );
	check( pBuff, NOMEM );
	*pBuff = 0;
}

TCharBuffer::TCharBuffer( const char *str )
{
		len = strlen( str );
	ptr = pBuff = (char*)malloc( len + 1 );
	check( pBuff, NOMEM );
	memcpy( pBuff, str, len + 1 );
}

TCharBuffer::TCharBuffer( const char *str, size_t length )
{
	len = length;
	ptr = pBuff = (char*)malloc( len + 1 );
	check( pBuff, NOMEM );
	memcpy( pBuff, str, len );
	pBuff[len] = 0;
}

char* TCharBuffer::Search( char c )
{
	ptr = strchr( ptr, c );
	if( ptr )
		return ptr;
	ptr = pBuff + len;
	return 0;
}

void  TCharBuffer::Skip( size_t size )
{
	ptr += size;
	assert( left() <= len );
}

void TCharBuffer::Delete( size_t length )
{
	if( length )
		del( length );
}

void TCharBuffer::Insert( const char *str )
{
	size_t length = strlen( str );
	if( length ) {
		ins( length );
		memcpy( ptr, str, length );
	}
}

void TCharBuffer::InsFile( const char *name )
{
	int handle = 0;
	handle = File::open( name, O_RDONLY | O_BINARY );
	long length = File::length( handle );
	check( length < UINT_MAX, CHB_OUTMEM );
	if( length ) {
		ins( (size_t)length );
		File::readRecord( handle, ptr, (uint16_t)length );
	}
	File::close( handle );
}

void TCharBuffer::Save( const char *name )
{
	int handle = File::open( name, O_WRONLY | O_CREAT |
					O_TRUNC | O_BINARY );
	File::write( handle, pBuff, len );
	File::close( handle );
}

void TCharBuffer::del( size_t length )
{
	long rest = (long)right() - length;
	assert( rest >= 0 );
	memmove( ptr, ptr + length, size_t( rest + 1 ) );
	len -= length;
	char *p = (char*)realloc( pBuff, len + 1 );
	check( p, NOMEM );
	ptr = p + left();
	pBuff = p;
}

void TCharBuffer::ins( size_t length )
{
	size_t rest = right();
	long size = (long)len + length;
	check( size < 0xff00l, CHB_OUTMEM );
	len = size_t( size );
	char *p = (char*)realloc( pBuff, len + 1 );
	check( p, NOMEM );
	ptr = p + left();
	pBuff = p;
	memmove( ptr + length, ptr, rest + 1 );
}
