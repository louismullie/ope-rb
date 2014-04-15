#include <ruby.h>
#include <math.h>
#include <openssl/rand.h>

/*
 *  Helper methods from Mathlib : A C Library of Special Functions
 *  Copyright (C) 2005 The R Foundation, released under the GPL.
 */
static double afc(int i)
{
    static const double al[9] =
    {
	0.0,
	0.0,/*ln(0!)=ln(1)*/
	0.0,/*ln(1!)=ln(1)*/
	0.69314718055994530941723212145817,/*ln(2) */
	1.79175946922805500081247735838070,/*ln(6) */
	3.17805383034794561964694160129705,/*ln(24)*/
	4.78749174278204599424770093452324,
	6.57925121201010099506017829290394,
	8.52516136106541430016553103634712
	/*, 10.60460290274525022841722740072165*/
    };
    double di, value;

    if (i < 0) {
      printf(("rhyper.c: afc(i), i=%d < 0 -- SHOULD NOT HAPPEN!\n"),
		      i);
      return -1;/* unreached (Wall) */
    } else if (i <= 7) {
	value = al[i + 1];
    } else {
	di = i;
	value = (di + 0.5) * log(di) - di + 0.08333333333333 / di
	    - 0.00277777777777 / di / di / di + 0.9189385332;
    }
    return value;
}

static int imax2(int x, int y) {
  return (x < y) ? y : x;
}

static int imin2(int x, int y) {
  return (x < y) ? x : y;
}

static VALUE ope_rb;
static VALUE ope_rb_cipher;

/* 
 * Initialize an OPE::Cipher object with a key.
 */
static VALUE ope_rb_initialize(VALUE self, VALUE key) {
 
  int keyLen;
  
  // Replace key value with key.to_str
  StringValue(key);
  
  // Get the key length as an int.
  keyLen = RSTRING_LEN(key);
  
  // Make sure key is not empty
  if (keyLen == 0) {
    rb_raise(rb_eArgError, "Key must be non-empty.");
  }
  
  // Set key as instance variable
  rb_iv_set(self, "@key", key);
  
  return self;
  
}

static unsigned char* ope_rb_get_key(VALUE self) {
  
  VALUE rbKey;
  unsigned char* cKey;
  
  // Get the key instance variable
  rbKey = rb_iv_get(self, "@key");
  StringValue(rbKey);
  
  // Convert the key to a byte array
  cKey = (unsigned char*) RSTRING_PTR(rbKey);
  
  // Return 1 upon success
  return cKey;
  
}

double unif_rand(unsigned char* seed, int seedLen) {

  unsigned char* prBits;
  double result;
  
  RAND_bytes(prBits, 2*seedLen / 8);
  
  memcpy(&result, prBits, sizeof(double));
  
  return result;
  
}

static VALUE ope_rb_rhyper(VALUE self, VALUE rbNn1in, VALUE rbNn2in, VALUE rbKkin, VALUE rbSeed) {
  
  double nn1in, nn2in, kkin;

  static const double con = 57.56462733;
  static const double deltal = 0.0078;
  static const double deltau = 0.0034;
  static const double scale = 1e25;
  
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
  
  nn1in = NUM2DBL(rbNn1in);
  nn2in = NUM2DBL(rbNn2in);
  kkin = NUM2DBL(rbKkin);
  
  StringValue(rbSeed);
  seed = (unsigned char*) RSTRING_PTR(rbSeed);
  seedLen = RSTRING_LEN(rbSeed);
  RAND_seed(seed, seedLen);
  
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
    /* squeeze using upper and lower bounds */
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
  
}
/*
static VALUE ope_rb_domain_gap(VALUE self, double ndomain, double nrange, double rgap) {
  
  return ope_rb_hyper(rgap, ndomain, nrange-ndomain);
  
}

static VALUE ope_rb_lazy_sample(VALUE self, double, d_lo, double d_hi, double r_lo, double r_hi, go_low) {
  
  double ndomain = d_hi - d_lo + 1;
  double nrange = r_hi - r_lo + 1;
  
  if (nrange < ndomain) {
    rb_raise(rb_eRuntimeError, 'Range smaller than domain.');
  }
  
  if (ndomain == 1) {
    return ope_domain_range(d_lo, r_lo, r_hi);
  }
  
  dgap_cache 
  
}

static VALUE ope_rb_domain_range_initialize(VALUE self, double d_arg, double r_lo_arg, double r_hi_arg) {
  
  rb_iv_set(self, "@d_arg", DBL2NUM(d_arg));
  rb_iv_set(self, "@r_lo_arg", DBL2NUM(r_lo_arg));
  rb_iv_set(self, "@r_hi_arg", DBL2NUM(r_hi_arg));
  
  return self;
  
}
*/
/*
double sampleHGD(double lowD, double highD, double lowR, double highR, double y, double coins, unsigned int coinsLen) {
  
  double whiteBalls, blackBalls, ballsPicked;
  
  whiteBalls = highD - lowD + 1;
  blackBalls = highR - lowR + 1 - whiteBalls;
  
  ballsPicked = y - lowR;
  
  if ((whiteBalls > 0) && (blackBalls >= 0) && (y >= lowR)) {
    rb_raise(rb_eRuntimeError, "HGD sampling problem");
  } 
  
  unsigned int precision = (unsigned int) numBits(highR-lowR+1) + 10;
  
  return lowD + hgd(ballsPicked, whiteBalls, blackBalls, coins, coinsLen);
  
}

double tapeGen(double lowD, double highD, double y, unsigned int desiredNoBits) {

  unsigned int desiredBytes;
  
  unsigned int totalLen;
  totalLen = sizeof(lowD) + sizeof(highD) + sizeof(y);
  
  unsigned char lowDBytes[sizeof(lowD)];
  unsigned char highDBytes[sizeof(highD)];
  unsigned char yBytes[sizeof(y)];
  unsigned char concat[totalLen];
  unsigned char shaDigest[SHA256_DIGEST_LEN];
  unsigned char seed[AES_BLOCK_BYTES];
  
  int i;
  
  if (desiredNoBits % 8 == 0) {
    rb_raise(rb_eRuntimeError, "desiredNoBits is not a multiple of 8");
  }
  
  desiredBytes = desiredNoBits / bitsPerByte;
  
  memcpy(&lowD, lowDBytes, sizeof(lowD));
  memcpy(&highD, highDBytes, sizeof(highD));
  memcpy(&y, y, sizeof(y));
  
  for (i = 0; i < totalLen; i++) {
    concat[i] = i < sizeof(lowD) ? lowDBytes[i] :
                (i < sizeof(lowD) + sizeof(highD) ?
                highDBytes[i] : yBytes[i]);
  }
  
  unsigned char* buffer = malloc(SHA256_DIGEST_LEN);
  SHA256_CTX sha256;
  SHA256_Init(&sha256);
  SHA256_Update(&sha256, concat, totalLen);
  SHA256_Final(shaDigest, &sha256);
  
  AES_encrypt(shaDigest, seed, {0x00});
  
  if (AES_BLOCK_BYTES >= desiredBytes) {
    //...
  }
  
  setSeed(seed);
  
  return randomLen(desiredNoBits);
  
}*/

void Init_native(void) {
  
	ope_rb = rb_define_module("OPE");
	
	ope_rb_cipher = rb_define_class_under(ope_rb, "Cipher", rb_cObject);
	
	rb_define_method(ope_rb_cipher, "initialize", ope_rb_initialize, 1);
	//rb_define_method(ope_rb_cipher, "domain_gap", ope_rb_domain_gap, 3);
	rb_define_method(ope_rb_cipher, "rhyper", ope_rb_rhyper, 4);

	//ope_rb_domain_range = rb_define_class_under(ope_rb, "DomainRange", rb_cObject);
	//rb_define_method(ope_rb_domain_range, "initialize", ope_rb_domain_range_initialize, 3);
	
  return;
	
}