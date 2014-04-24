#include <ruby.h>
#include <math.h>
#include <openssl/rand.h>
#include <imath.h>
#include <imrat.h>

/*
 *  Helper methods from Mathlib : A C Library of Special Functions
 *  Copyright (C) 2005 The R Foundation, released under the GPL.
 */
static double afc(int i)
{
    double di, value;
/*
	di = i;
	value = (di + 0.5) * log(di) - di + 0.08333333333333 / di
	    - 0.00277777777777 / di / di / di + 0.9189385332;
    }*/
    return value;
}

static int imax2(int x, int y) {
  return (x < y) ? y : x;
}

static int imin2(int x, int y) {
  return (x < y) ? x : y;
}

static VALUE ope_rb;
static VALUE ope_rb_hgd;

static void mp_print_int(const char* msg, mpz_t a) {
  
  char  *buf;
  int    len;  
  
  len = mp_int_string_len(&a, 10);
  buf = calloc(len, sizeof(*buf));
  mp_int_to_string(&a, 10, buf, len);
  printf("(int) %s = %s\n", msg, buf);
  
}

static void mp_print_dec(const char* msg, mpq_t a) {
  
  char  *buf;
  int    len;  
  
  len = mp_rat_decimal_len(&a, 10, (mp_size) 10);
  buf = calloc(len, sizeof(*buf));
  mp_rat_to_decimal(&a, 10, (mp_size) 10, MP_ROUND_DOWN, buf, len);
  printf("(dec) %s = %s\n", msg, buf);
  
}

void int_to_rat(mpz_t a, mpq_t rat) {
  
  char *buf; int len;
  
  len = mp_int_string_len(&a, 10);
  buf = calloc(len, sizeof(*buf));
  mp_int_to_string(&a, 10, buf, len);
  
  mp_rat_init(&rat);
  mp_rat_read_string(&rat, 10, buf);
  
  return;
  
} 

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
  
}

/*(di + 0.5) * log(di) - di + 0.08333333333333 / di
    - 0.00277777777777 / di / di / di + 0.9189385332;*/

static mpz_t rhyper(mpz_t kk, mpz_t nn1, mpz_t nn2, long prec) {
  
  mpz_t result;
  
  // Static (constant) values
  static mpq_t con;
  static mpq_t deltal;
  static mpq_t deltau;
  static mpq_t scale;
  
  // Dynamic (computed) rationals
  mpq_t jx; mpq_t tn; mpq_t n1;
  mpq_t n2; mpq_t k; mpq_t p;
  mpq_t u; mpq_t v; mpq_t a;
  mpq_t ix; mpq_t xl; mpq_t xr;
  mpq_t m; mpq_t kl; mpq_t kr;
  mpq_t lamdl; mpq_t lamdr; mpq_t nk;
  mpq_t nm; mpq_t p1; mpq_t p2;
  mpq_t p3; mpq_t minjx; mpq_t maxjx;
  
  // Temporary variables
  mpq_t kk_dec; mpq_t kk2;
  
  // Initialize and read in constant values
  mp_rat_init(&con); mp_rat_init(&deltal);
  mp_rat_init(&deltau); mp_rat_init(&scale);
  
  mp_rat_read_decimal(&con, 10, "57.56462733");
  mp_rat_read_decimal(&deltal, 10, "0.0078");
  mp_rat_read_decimal(&deltau, 10, "0.0034");
  mp_rat_read_decimal(&scale, 10, "10000000000000000000000000");
  
  // Initialize computed rationals
  mp_rat_init(&jx); mp_rat_init(&tn); mp_rat_init(&n1);
  mp_rat_init(&n2); mp_rat_init(&k); mp_rat_init(&p);
  mp_rat_init(&u); mp_rat_init(&v); mp_rat_init(&a);
  mp_rat_init(&ix); mp_rat_init(&xl); mp_rat_init(&xr);
  mp_rat_init(&m); mp_rat_init(&kl); mp_rat_init(&kr);
  mp_rat_init(&lamdl); mp_rat_init(&lamdr); mp_rat_init(&nk);
  mp_rat_init(&nm); mp_rat_init(&p1); mp_rat_init(&p2);
  mp_rat_init(&p3); mp_rat_init(&minjx); mp_rat_init(&maxjx); 
  
  // Temporary vairables
  mp_rat_init(&kk_dec); mp_rat_init(&kk2);
  
  // Print input values
  mp_print_int("nn1", nn1);
  mp_print_int("nn2", nn2);
  mp_print_int("kk", kk);
  
  // Print constant values
  mp_print_dec("con", con);
  mp_print_dec("deltal", deltal);
  mp_print_dec("deltau", deltau);
  mp_print_dec("scale", scale);
  
  if (mp_int_compare(&nn1, &nn2) >= 0) {
    int_to_rat(nn2, n1);
    int_to_rat(nn1, n2);
  } else {
    int_to_rat(nn1, n1);
    int_to_rat(nn2, n2);
  }
  
  mp_rat_add(&n1, &n2, &tn);
  int_to_rat(kk, kk_dec);
  mp_rat_add(&kk_dec, &kk_dec, &kk2);
  
  if (mp_rat_compare(&kk2, &tn) >= 0) {
    mp_rat_sub(&tn, &kk_dec, &k);
  } else {
    k = kk_dec;
  }
  
  // Print n1, n2 values
  mp_print_dec("n1", n1);
  mp_print_dec("n2", n2);
  mp_print_dec("tn", tn);
  mp_print_dec("k", k);
  
  // Clear input values
  mp_int_clear(&nn1); mp_int_clear(&nn2); mp_int_clear(&kk);
  
  // Clear constant values
  mp_rat_clear(&con); mp_rat_clear(&deltal);
  mp_rat_clear(&deltau); mp_rat_clear(&scale);
  
  // Clear computed rationals.
  mp_rat_clear(&jx); mp_rat_clear(&tn); mp_rat_clear(&n1);
  mp_rat_clear(&n2); mp_rat_clear(&k); mp_rat_clear(&p);
  mp_rat_clear(&u); mp_rat_clear(&v); mp_rat_clear(&a);
  mp_rat_clear(&ix); mp_rat_clear(&xl); mp_rat_clear(&xr);
  mp_rat_clear(&m); mp_rat_clear(&kl); mp_rat_clear(&kr);
  mp_rat_clear(&lamdl); mp_rat_clear(&lamdr); mp_rat_clear(&nk);
  mp_rat_clear(&nm); mp_rat_clear(&p1); mp_rat_clear(&p2);
  mp_rat_clear(&p3); mp_rat_clear(&minjx); mp_rat_clear(&maxjx);
  
  mp_rat_clear(&kk_dec); mp_rat_clear(&kk2);
  
  // Get result
  mp_int_init(&result);
  mp_int_init_value(&result, 1);
  
  return result;
  
}

static VALUE ope_rb_rhyper(VALUE self, VALUE rb_nn1, VALUE rb_nn2, VALUE rb_kk, VALUE rb_seed, VALUE rb_prec) {
  
  // VALUE -> CSTRING -> MPZ_T
  char* nn1_str; char* nn2_str; char* kk_str;
  mpz_t nn1; mpz_t nn2; mpz_t kk; int prec;
  
  // Convert parameters to string
  nn1_str = StringValuePtr(rb_nn1);
  nn2_str = StringValuePtr(rb_nn2);
  kk_str = StringValuePtr(rb_kk);
  
  // Initialize the integers
  mp_int_init(&nn1); mp_int_read_string(&nn1, 10, nn1_str);
  mp_int_init(&nn2); mp_int_read_string(&nn2, 10, nn2_str);
  mp_int_init(&kk); mp_int_read_string(&kk, 10, kk_str);
  
  // Convert precision to integer
  prec = NUM2INT(rb_prec);
  
  // Call rhyper function
  rhyper(nn1, nn2, kk, prec);
  
  //
  
  return self;
  
/*

  int nn1, nn2, kk;
  int i, ix;
  int reject, setup1, setup2;

  double e, f, g, p, r, t, u, v, y;
  double de, dg, dr, ds, dt, gl, gu, nk, nm, ub;
  double xk, xm, xn, y1, ym, yn, yk, alv;

  static int ks = -1;
  static int n1s = -1, n2s = -1;

  static int k, m;
  static int minjx, maxjx, n1, n2;

  static double a, d, s, w;
  static double tn, xl, xr, kl, kr, lamdl, lamdr, p1, p2, p3;
  
  unsigned char* seed;
  int seedLen;
  
  if(!isfinite(nn1in) || !isfinite(nn2in) || !isfinite(kkin))
    return -1;

  nn1 = (int) nn1in;
  nn2 = (int) nn2in;
  kk	= (int) kkin;

  if (nn1 < 0 || nn2 < 0 || kk < 0 || kk > nn1 + nn2) {
    //rb_raise(rb_eRuntimeError, 'Incorrect arguments to HGD');
    return INT2NUM(-2);
  }

  reject = 1;
  if (nn1 != n1s || nn2 != n2s) {
  setup1 = 1;	setup2 = 1;
  } else if (kk != ks) {
  setup1 = 0;	setup2 = 1;
  } else {
  setup1 = 0;	setup2 = 0;
  }
  if (setup1) {
  n1s = nn1;
  n2s = nn2;
  tn = nn1 + nn2;
  if (nn1 <= nn2) {
    n1 = nn1;
    n2 = nn2;
  } else {
    n1 = nn2;
    n2 = nn1;
  }
  }
  if (setup2) {
  ks = kk;
  if (kk + kk >= tn) {
    k = (int)(tn - kk);
  } else {
    k = kk;
  }
  }
  if (setup1 || setup2) {
  m = (int) ((k + 1.0) * (n1 + 1.0) / (tn + 2.0));
  minjx = imax2(0, k - n2);
  maxjx = imin2(n1, k);
  }

  if (minjx == maxjx) {
  ix = maxjx;

  if (kk + kk >= tn) {
  if (nn1 > nn2) {
    ix = kk - nn2 + ix;
  } else {
    ix = nn1 - ix;
  }
  } else {
  if (nn1 > nn2)
    ix = kk - ix;
  }
  
  return DBL2NUM(ix);

  } else if (m - minjx < 10) {
  if (setup1 || setup2) {
    if (k < n2) {
  w = exp(con + afc(n2) + afc(n1 + n2 - k)
  	- afc(n2 - k) - afc(n1 + n2));
    } else {
  w = exp(con + afc(n1) + afc(k)
  	- afc(k - n2) - afc(n1 + n2));
    }
  }
    L10:
  p = w;
  ix = minjx;
  u = unif_rand(seed, seedLen) * scale;
    L20:
  if (u > p) {
    u -= p;
    p *= (n1 - ix) * (k - ix);
    ix++;
    p = p / ix / (n2 - k + ix);
    if (ix > maxjx)
  goto L10;
    goto L20;
  }
  } else {

  if (setup1 || setup2) {
    s = sqrt((tn - k) * k * n1 * n2 / (tn - 1) / tn / tn);

    d = (int) (1.5 * s) + .5;
    xl = m - d + .5;
    xr = m + d + .5;
    a = afc(m) + afc(n1 - m) + afc(k - m) + afc(n2 - k + m);
    kl = exp(a - afc((int) (xl)) - afc((int) (n1 - xl))
       - afc((int) (k - xl))
       - afc((int) (n2 - k + xl)));
    kr = exp(a - afc((int) (xr - 1))
       - afc((int) (n1 - xr + 1))
       - afc((int) (k - xr + 1))
       - afc((int) (n2 - k + xr - 1)));
    lamdl = -log(xl * (n2 - k + xl) / (n1 - xl + 1) / (k - xl + 1));
    lamdr = -log((n1 - xr + 1) * (k - xr + 1) / xr / (n2 - k + xr));
    p1 = d + d;
    p2 = p1 + kl / lamdl;
    p3 = p2 + kr / lamdr;
  }
    L30:
  u = unif_rand(seed, seedLen) * p3;
  v = unif_rand(seed, seedLen);
  if (u < p1) {
    ix = (int) (xl + u);
  } else if (u <= p2) {
    ix = (int) (xl + log(v) / lamdl);
    if (ix < minjx)
  goto L30;
    v = v * (u - p1) * lamdl;
  } else {
    ix = (int) (xr - log(v) / lamdr);
    if (ix > maxjx)
  goto L30;
    v = v * (u - p2) * lamdr;
  }

  if (m < 100 || ix <= 50) {
    f = 1.0;
    if (m < ix) {
  for (i = m + 1; i <= ix; i++)
      f = f * (n1 - i + 1) * (k - i + 1) / (n2 - k + i) / i;
    } else if (m > ix) {
  for (i = ix + 1; i <= m; i++)
      f = f * i * (n2 - k + i) / (n1 - i + 1) / (k - i + 1);
    }
    if (v <= f) {
  reject = 0;
    }
  } else {
    
    y = ix;
    y1 = y + 1.0;
    ym = y - m;
    yn = n1 - y + 1.0;
    yk = k - y + 1.0;
    nk = n2 - k + y1;
    r = -ym / y1;
    s = ym / yn;
    t = ym / yk;
    e = -ym / nk;
    g = yn * yk / (y1 * nk) - 1.0;
    dg = 1.0;
    if (g < 0.0)
  dg = 1.0 + g;
    gu = g * (1.0 + g * (-0.5 + g / 3.0));
    gl = gu - .25 * (g * g * g * g) / dg;
    xm = m + 0.5;
    xn = n1 - m + 0.5;
    xk = k - m + 0.5;
    nm = n2 - k + xm;
    ub = y * gu - m * gl + deltau
  + xm * r * (1. + r * (-0.5 + r / 3.0))
  + xn * s * (1. + s * (-0.5 + s / 3.0))
  + xk * t * (1. + t * (-0.5 + t / 3.0))
  + nm * e * (1. + e * (-0.5 + e / 3.0));
    alv = log(v);
    if (alv > ub) {
  reject = 1;
    } else {
  dr = xm * (r * r * r * r);
  if (r < 0.0)
      dr /= (1.0 + r);
  ds = xn * (s * s * s * s);
  if (s < 0.0)
      ds /= (1.0 + s);
  dt = xk * (t * t * t * t);
  if (t < 0.0)
      dt /= (1.0 + t);
  de = nm * (e * e * e * e);
  if (e < 0.0)
      de /= (1.0 + e);
  if (alv < ub - 0.25 * (dr + ds + dt + de)
      + (y + m) * (gl - gu) - deltal) {
      reject = 0;
  }
  else {
      if (alv <= (a - afc(ix) - afc(n1 - ix)
  		- afc(k - ix) - afc(n2 - k + ix))) {
  	reject = 0;
      } else {
  	reject = 1;
      }
  }
    }
  }
  if (reject)
    goto L30;
  }
  
  if (kk + kk >= tn) {
    if (nn1 > nn2) {
      ix = kk - nn2 + ix;
    } else {
      ix = nn1 - ix;
    }
  } else {
  if (nn1 > nn2)
    ix = kk - ix;
  }
  
  return DBL2NUM(ix);
*/
  
}

void Init_native(void) {
  
	ope_rb = rb_define_module("OPE");
	
	ope_rb_hgd = rb_define_class_under(ope_rb, "HGD", rb_cObject);
	rb_define_singleton_method(ope_rb_hgd, "rhyper_native", ope_rb_rhyper, 5);
	rb_define_singleton_method(ope_rb_hgd, "afc_native", ope_rb_afc, 1);

  return;
	
}