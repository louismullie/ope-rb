require 'mkmf'

dir_config('openssl')
have_header("openssl/ssl.h")
have_library("ssl", "SSLv23_method")
have_header('gmp.h')
have_library('gmp', '__gmpz_init')
have_header('mpfr.h')
have_header('mpf2mpfr.h')
have_library('mpfr', 'mpfr_init')
create_makefile('ope/native')