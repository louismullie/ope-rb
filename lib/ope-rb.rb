require 'ope/native'
require 'ope-rb/version'
require 'ope-rb/prng'
require 'ope-rb/hgd'
require 'ope-rb/crypto'
require 'ope-rb/errors'

module OPE
  
  class Cipher
    
    OPE_KEY_SIZE = 128
    
    include Crypto
    
    def initialize(key, pt_len = 8, ct_len = 16)
      
      check_key_size(key)
        
      @key, @pt_len, @ct_len =
      key, pt_len * 8, ct_len * 8
      
      create_cipher
      
    end
    
    def check_key_size(key)
      
      valid = key.bytesize * 8 == OPE_KEY_SIZE
      raise Errors::KeySizeError unless valid
      
    end
    
    def sample_hgd(lo_d, hi_d, lo_r, hi_r, y, coins)
      
      wb, bp = hi_d - lo_d + 1, y - lo_r
      bb = hi_r - lo_r + 1 - wb
      
      unless wb > 0 && bb >= 0 && y >= lo_r && y <= hi_r
        raise Errors::HGDParameterError
      end
      
      prec = (hi_r - lo_r + 1).size * 8 + 10

      lo_d + HGD.rhyper(bp, wb, bb, coins, prec)
      
    end
    
    def tape_gen(lo_d, hi_d, y, desired_no_bits)
      
      if desired_no_bits % 8 != 0
        raise Errors::NumBitsError
      elsif (desired_bytes = desired_no_bits / 8) > 16
        raise Errors::DesiredBytesError
      end
      
      digest = hash_concat(lo_d, hi_d, y)
      seed = get_seed(digest, desired_bytes)
      
    end
    
    def encrypt(m)
      
      lo_d, hi_d = 0, (1 << @pt_len) - 1
      lo_r, hi_r = 0, (1 << @ct_len) - 1
      
      encrypt_recurse(lo_d, hi_d, lo_r, hi_r, m)
      
    end
    
    def decrypt(m)

      lo_d, hi_d = 0, (1 << @pt_len) - 1
      lo_r, hi_r = 0, (1 << @ct_len) - 1

      decrypt_recurse(lo_d, hi_d, lo_r, hi_r, m)
      
    end
    
    private
    
    def get_seed(digest, desired_bytes)
      
      # Encrypt the hash using the cipher
      seed = cipher_encrypt(digest)
      
      # Take only the desired no of bytes
      seed = seed.byteslice(0, desired_bytes)
   
      # Convert the seed to a Bignum
      seed.unpack('H*')[0].hex
      
    end
    
    def encrypt_recurse(lo_d, hi_d, lo_r, hi_r, m)
      
      m2, n = hi_d - lo_d + 1, hi_r - lo_r + 1
      d, r = lo_d - 1, lo_r - 1; y = r + (n + 1) / 2
      
      raise Errors::IncorrectMValueError unless m2 > 0

      coins = nil
      
      if m2 == 1
        
        coins = tape_gen(lo_d, hi_d, m, @ct_len)
        return lo_r + (coins % n)
        
      end
      
      coins = tape_gen(lo_d, hi_d, y, @pt_len)
      
      x = sample_hgd(lo_d, hi_d, lo_r, hi_r, y, coins)
      
      lo_d, hi_d = *((m <= x) ? [d + 1, x] : [x + 1, d + m2])
      lo_r, hi_r = *((m <= x) ? [r + 1, y] : [y + 1, r + n])
      
      encrypt_recurse(lo_d, hi_d, lo_r, hi_r, m)
      
    end
    
    def decrypt_recurse(lo_d, hi_d, lo_r, hi_r, c)
      
      m2, n = hi_d - lo_d + 1, hi_r - lo_r + 1
      d, r = lo_d - 1, lo_r - 1; y = r + (n + 1) / 2
      
      raise Errors::IncorrectMValueError unless m2 > 0
      
      if m2 == 1

        m = lo_d.to_i
        
        coins = tape_gen(lo_d, hi_d, m, @ct_len)
        w = lo_r + coins % n
        
        return m.to_i if w == c
        
        raise Errors::BadDecryptError
        
      end
      
      coins = tape_gen(lo_d, hi_d, y, @pt_len)
      
      x = sample_hgd(lo_d, hi_d, lo_r, hi_r, y, coins)
      
      lo_d, hi_d = *((c <= y) ? [d + 1, x] : [x + 1, d + m2])
      lo_r, hi_r = *((c <= y) ? [r + 1, y] : [y + 1, r + n])

      decrypt_recurse(lo_d, hi_d, lo_r, hi_r, c)
      
    end
    
  end
  
end
