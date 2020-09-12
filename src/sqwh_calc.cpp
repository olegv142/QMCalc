#include "express.h"
#include "preproc.h"
#include "sqwh.h"

#include <iostream>

int main(int argc, char* argv[])
{
	if (argc < 2) {
		std::cout << "Single Quantum Well Holes model. Usage:\n"
			<< "sqwh #(params file path) [param=val ...]\n";
		return 1;
	}
	TExpression x; 
	for( int i = 1 ; i < argc ; i++ ) {
		TMacroProcessor pp( argv[i] );
		x.SetExpression( pp.Run() );
		x.Eval();
	}
	struct SQWHParams params;
	GetSQWHParams(params, x);
	SQWHSolver s(params);
	s.Solve();
	s.SaveResults("sqwh");
	return 0;
}

