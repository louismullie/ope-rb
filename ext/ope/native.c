#include <ruby.h>
#include <math.h>
#include <openssl/rand.h>
#include <gmp.h>
#include <mpfr.h>

/*
static VALUE ope_rb_afc(VALUE self, VALUE num) {
  
  mpq_t i; mp_frac_1; mpq_t frac_12; mpq_t frac_pi; mpq_t frac_360;
  VALUE num_str;
  mpq_t total;
  
  mp_rat_init(&i); mp_rat_init(&frac_1); mp_rat_init(&frac_12);
  mp_rat_init(&frac_360); mp_rat_init(&frac_pi);
  mp_rat_init(&total);
  
  num_str = rb_funcall(num, rb_intern("to_s"), 0);
  
  mp_rat_read_decimal(&i, 10, StringValuePtr(num_str));
  mp_rat_set_value(&frac_1, 1, 2);
  mp_rat_set_value(&frac_12, 1, 12);
  mp_rat_set_value(&frac_360, 1, 360);
  mp_rat_read_decimal(&frac_pi, 10, "0.9189385332");      // 0.5 * log(2*PI)
  
  mp_print_dec("frac_12", frac_12);
  
  mp_rat_copy(&i, &total);
  mp_rat_add(&total, &frac_1, &total);
  
  mp_rat_clear(&i); mp_rat_clear(&frac_12);
  mp_rat_clear(&frac_360); mp_rat_clear(&frac_pi);
  
  return self;
  
}*/

static int ope_afc(mpfr_t i, mpfr_t result) {
  
  mpfr_t frac_12, frac_360, frac_pi, tmp;
  mpfr_inits(frac_12, frac_360, frac_pi, tmp, NULL);
  
  mpfr_set_d(frac_12, 0.08333333333333, MPFR_RNDN);
  mpfr_set_d(frac_360, 0.00277777777777, MPFR_RNDN);
  mpfr_set_d(frac_pi, 0.9189385332, MPFR_RNDN);

  mpfr_add_d(result, i, 0.5, MPFR_RNDN);
  mpfr_log(tmp, i, MPFR_RNDN);
  mpfr_mul(result, result, tmp, MPFR_RNDN);
  mpfr_sub(result, result, i, MPFR_RNDN);
  mpfr_div(tmp, frac_12, i, MPFR_RNDN);
  mpfr_add(result, result, tmp, MPFR_RNDN);
  mpfr_div(tmp, frac_360, i, MPFR_RNDN);
  mpfr_div(tmp, tmp, i, MPFR_RNDN);
  mpfr_div(tmp, tmp, i, MPFR_RNDN);
  mpfr_add(result, result, tmp, MPFR_RNDN);
  mpfr_add(result, result, frac_pi, MPFR_RNDN);
  
  mpfr_clears(frac_12, frac_360, frac_pi, tmp, NULL);
  
  return 0;
  
}

static int num_to_mpfr(VALUE num, mpfr_t result) {
  
  VALUE num_str;
  
  num_str = rb_funcall(num, rb_intern("to_s"), 0);
  
  if (mpfr_set_str(result, StringValuePtr(num_str), 10, MPFR_RNDN) < 0) {
    rb_raise(rb_eRuntimeError, "Could not set string.");
  }
  
  return 0;
  
}

static VALUE afc_op_3(VALUE self, VALUE m_num, VALUE n1_num, VALUE n2_num, VALUE k_num) {

  mpfr_t m, n1, n2, k, tmp, tmp2, out;
  mpfr_inits(m, n1, n2, k, tmp, tmp2, out, NULL);
  
  char* str[32];
  
  num_to_mpfr(m_num, m); num_to_mpfr(n1_num, n1);
  num_to_mpfr(n2_num, n2); num_to_mpfr(k_num, k);

  ope_afc(m, out);
  mpfr_sub(tmp, n1, m, MPFR_RNDN);
  ope_afc(tmp, tmp2);
  mpfr_add(out, out, tmp2, MPFR_RNDN);
  mpfr_sub(tmp, k, m, MPFR_RNDN);
  ope_afc(tmp, tmp2);
  mpfr_add(out, out, tmp2, MPFR_RNDN);
  mpfr_sub(tmp, n2, k, MPFR_RNDN);
  mpfr_add(tmp2, tmp, m, MPFR_RNDN);
  ope_afc(tmp2, tmp);
  mpfr_add(out, out, tmp, MPFR_RNDN);
  
  mpfr_sprintf(str, "%.128RNf", out);
  return rb_str_new(str, 32);
  
}

// Math.exp(a - afc(xr-1) - afc(n1-xr+1) - afc(k-xr+1) - afc(n2-k+xr-1))

static VALUE afc_op_5(VALUE self, VALUE a_num, VALUE xr_num,
                      VALUE n1_num, VALUE n2_num, VALUE k_num) {

  mpfr_t a, xr, n1, n2, k, tmp, tmp2, out;
  mpfr_inits(a, xr, n1, n2, k, tmp, tmp2, out, NULL);
  
  char* str[32];
  
  num_to_mpfr(a_num, tmp);
  num_to_mpfr(xr_num, xr);
  num_to_mpfr(n1_num, n1);
  num_to_mpfr(n2_num, n2);
  num_to_mpfr(k_num, k);

  mpfr_sub(tmp, xr, m, MPFR_RNDN);
  ope_afc(tmp, tmp2);
  mpfr_add(out, out, tmp2, MPFR_RNDN);
  mpfr_sub(tmp, n2, k, MPFR_RNDN);
  mpfr_add(tmp2, tmp, m, MPFR_RNDN);
  ope_afc(tmp2, tmp);
  mpfr_add(out, out, tmp, MPFR_RNDN);
  
  mpfr_sprintf(str, "%.128RNf", out);
  return rb_str_new(str, 32);
  
}

static VALUE ope_rb_afc(VALUE self, VALUE num) {
  
  VALUE num_str;
  mpfr_t i, result;
  mpfr_inits(i, result, NULL);
  char* str[32];
  
  num_str = rb_funcall(num, rb_intern("to_s"), 0);
  
  if (mpfr_set_str(i, StringValuePtr(num_str), 10, MPFR_RNDN) < 0) {
    rb_raise(rb_eRuntimeError, "Could not set string.");
  }
  
  if (ope_afc(i, result) < 0) {
    rb_raise(rb_eRuntimeError, "Could not calculate AFC."); 
  }
  
  mpfr_sprintf (str, "%.128RNf", result);
  return rb_str_new(str, 32);
  
}

static VALUE ope_rb;

void Init_native(void) {
  
  static VALUE ope_rb_hgd;
  
	ope_rb = rb_define_module("OPE");
	ope_rb_hgd = rb_define_class_under(ope_rb, "HGD", rb_cObject);
	//rb_define_singleton_method(ope_rb_hgd, "rhyper_native", ope_rb_rhyper, 5);
	rb_define_singleton_method(ope_rb_hgd, "afc_native", ope_rb_afc, 1);
	rb_define_singleton_method(ope_rb_hgd, "afc_native_op_3", afc_op_3, 4);
	rb_define_singleton_method(ope_rb_hgd, "afc_native_op_5", afc_op_5, 5);

  return;
	
}