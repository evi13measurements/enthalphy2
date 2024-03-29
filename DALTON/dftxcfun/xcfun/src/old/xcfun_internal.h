#ifndef XCFUN_INTERNAL_H
#define XCFUN_INTERNAL_H
#define XCFUN_INTERNAL
#include "xcfun.h"
#include "config.h"
#include "functional.h"
//#include "array.h"
#include "parameters.h"

struct xc_functional_data
{
public:
  void initialize();
  void destroy();

  void set_mode(int mode);

  int get_type(void) const;
  int get_max_order(void) const;

  void eval(ireal_t *res, int order, const ireal_t *dens) const;

  int input_length(void) const;
  int output_length(int order) const;
  int derivative_index(const int derivative[]) const;

  void find_max_order(void);

  int mode; // One of XC_VARS_*
  int type; // LDA, GGA etc
  int max_order; // Maximum derivative order with current settings

  double parameters[XC_NR_PARAMS];
  //array<functional *> active_functionals;
};

void xcint_die(const char *message, int code);
int xcint_input_length(int mode, int type);
int xcint_output_length(int mode, int type, int order);

typedef void (*evaluator)(const xc_functional_data &fun, ireal_t *, const ireal_t *);
evaluator xc_evaluator_lookup(int mode, int type, int order);

void xcint_setup_functionals();

void xcint_assure_setup();

#endif
