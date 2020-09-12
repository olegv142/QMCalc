#include "express.h"
#include "preproc.h"
#include "sqwe.h"

#include <iostream>

int main(int argc, char* argv[])
{
	if (argc < 2) {
		std::cout << "Single Quantum Well Electrons model. Usage:\n"
			<< "sqwe #(params file path) [param=val ...]\n";
		return 1;
	}
	TExpression x; 
	for( int i = 1 ; i < argc ; i++ ) {
		TMacroProcessor pp( argv[i] );
		x.SetExpression( pp.Run() );
		x.Eval();
	}
	struct SQWEParams params;
	GetSQWEParams(params, x);
	SQWESolver s(params);
	s.Solve();
	s.SaveResults("sqwe.dat");
	return 0;
}

