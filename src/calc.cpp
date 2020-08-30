#include "express.h"
#include "preproc.h"

int main(int argc, char* argv[])
{
	TExpression calc; 
	for( int i = 1 ; i < argc ; i++ ) {
		TMacroProcessor pp( argv[i] );
		calc.SetExpression( pp.Run() );
		calc.Calc();
	}
	calc.Dump();
}

