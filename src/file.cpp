#define _CRT_SECURE_NO_WARNINGS

#include "file.h"

DefineErrorMess( DE_NOFILE, "File not found" )
DefineErrorMess( DE_NOPATH, "Path not found" )
DefineErrorMess( DE_NOOPEN, "Can't open file" )
DefineErrorMess( DE_NOREAD, "Can't read the file" )
DefineErrorMess( DE_NOWRITE, "Can't write the file" )

int File::open( const char *fname, int mode, int perm )
{
	int handle = _open( fname, mode, perm );
	check3( handle != -1, DE_NOOPEN, fname );
	return handle;
}

// Read functions:
int     File::read( int handle, void *buf, unsigned nBytes )
{
	int size = _read( handle, buf, nBytes );
	check( size != -1, DE_NOREAD );
	return size;
}

void    File::readRecord( int handle, void *buf, unsigned nBytes )
{
	unsigned size = read( handle, buf, nBytes );
	check( size == nBytes, DE_NOREAD );
}

// Write function: 
void    File::write( int handle, void *buf, unsigned nBytes )
{
	unsigned size = _write( handle, buf, nBytes );
	check( size == nBytes, DE_NOWRITE );
}
