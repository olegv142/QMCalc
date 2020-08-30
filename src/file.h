#pragma once

#include "errs.h"

#include <io.h>
#include <fcntl.h>
#include <sys\stat.h>
#include <limits.h>
#include "stdint.h"

ErrorCode( DE_NOFILE )
ErrorCode( DE_NOPATH )
ErrorCode( DE_NOOPEN )
ErrorCode( DE_NOREAD )
ErrorCode( DE_NOWRITE )

class File
{
public:
// Open \ close functions:
static  int  open( const char *fname, int mode, int perm = _S_IREAD | _S_IWRITE );
static  void close( int handle );

// Read functions:
static  int      read( int handle, void *buf, unsigned nBytes );
static  void     readRecord( int handle, void *buf, unsigned nBytes );
static  uint8_t  getByte( int handle );
static  uint16_t getWord( int handle );
static  uint32_t getDword( int handle );

// Write functions: 
static  void write( int handle, void *buf, unsigned nBytes );
static  void put( int handle, uint8_t v );
static  void put( int handle, uint16_t v );
static  void put( int handle, uint32_t v );

// Miscellanea:
static  long seek( int handle, long offset, int fromWhere );
static  long tell( int handle );
static  long length( int handle );
static  void chsize( int handle, long newSize );
static  bool eof( int handle );
};

inline  void File::close( int handle )
{
	_close( handle );
}

inline  long File::seek( int handle, long offset, int fromWhere )
{
	return _lseek( handle, offset, fromWhere );
}

inline  long File::tell( int handle )
{
	return _tell( handle );
}

inline  long File::length( int handle )
{
	return _filelength( handle );
}

inline  void File::chsize( int handle, long newSize )
{
	_chsize( handle, newSize );
}

inline  bool File::eof( int handle )
{
	return _eof( handle ) != 0;
}

// Read functions:
// These functions read Byte, Word, and DWORD accordingly and
// check that they have really read
//
inline uint8_t File::getByte( int handle )
{
	uint8_t ret;
	readRecord( handle, &ret, sizeof(ret) );
	return ret;
}

inline uint16_t File::getWord( int handle )
{
	uint16_t ret;
	readRecord( handle, &ret, sizeof(ret) );
	return ret;
}

inline uint32_t File::getDword( int handle )
{
	uint32_t ret;
	readRecord( handle, &ret, sizeof(ret) );
	return ret;
}

// These functions write Byte, Word, and DWORD accordingly
//
inline void File::put( int handle, uint8_t b )
{
	write( handle, &b, sizeof(b) );
}

inline void File::put( int handle, uint16_t w )
{
	write( handle, &w, sizeof(w) );
}

inline void File::put( int handle, uint32_t d )
{
	write( handle, &d, sizeof(d) );
}

