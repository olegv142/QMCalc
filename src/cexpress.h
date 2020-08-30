#pragma once

#include "express.h"

class TCompoundExpression : public TExpression {
public:
	TCompoundExpression( const char* expr = 0, const char *end = 0, TCompoundExpression* par = 0 )
		: TExpression( expr ) { endPtr = end; parent = par; }
	~TCompoundExpression() { curExpression = parent; }

	void   Assign( const char* pname, const char* pend, double value );
	double Look( const char* pname, const char* pend );

	double Calc();

	void Dump( string prefix = "" );

protected:
	double  expr();
	double  prim();

	TCompoundExpression* parent;

private:
	TokenValue skipExpr();
	TokenValue searchElse();

	const char *endPtr;
};

//------ Gramma ------------------------------------------------------------
//  expression_list :
//      expression
//      expression_list ; expression
//
//  expression :
//      empty_string
//      log
//      log ? expression : expression    ( condition ? true : false )
//      log @ expression                 ( condition @ loop )
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
//      { expression_list }
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
//      ret         - return
//      brk         - break loop
//      cbrk        - break loop on Ctr-Break
//
//  external_function treated as data file, which define function
//--------------------------------------------------------------------------
